#!Python
# ===================================================================
# Name:		PACEtomo
# Purpose:	Runs parallel dose-symmetric tilt series on many targets with geometrical predictions using the first target as tracking tilt series.
# 		Make sure to run selectTargets script first to generate compatible Navigator settings and a target file.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2021/04/16
# Revision:	v1.3.2
# Last Change:	2022/07/12: added CTFfind option, added zero tilt exposure time option
#		2022/06/03: progress
#		2022/05/30: geoRefine
# ===================================================================

import serialem as sem
import os
import copy
import time
import numpy as np
from scipy import optimize

############ SETTINGS ############ 

startTilt	= 0		# starting tilt angle (should be divisible by step)
step		= 3		# tilt step
maxTilt		= 60		# maximum +/- tilt angle relative to startTilt
minDefocus	= -5		# minimum defocus of target range
maxDefocus	= -5		# maximum defocus of target range
stepDefocus	= 0.5		# step between target defoci (between TS)

focusSlope	= 0.0		# empirical linear focus correction [microns per degree] (obtained by linear regression of CTF fitted defoci over tilt series; microscope stage dependent)
delayIS		= 0.3 		# delay between applying image shift and Record
delayTilt	= 0 		# delay after stage tilt
zeroExpTime	= 0 		# set to exposure time [s] used for start tilt image, if 0: use same exposure time for all tilt images [NOT IN VIDEO TUTORIAL]

# Geometry settings
pretilt		= 0		# pretilt of sample in deg e.g. after FIB milling (if milling direction is not perpendicular to the tilt axis, estimate and add rotation)
rotation	= 0		# rotation of lamella vs tilt axis in deg (should be 0 deg if lamella is oriented perpendicular to tilt axis)

# Holey support settings
tgtPattern 	= False		# use same tgt pattern on different stage positions (useful for collection on holey support film)
alignToP	= False		# use generic image in buffer P to align to

# Persistent settings
beamTiltComp	= True		# use beam tilt compensation (uses coma vs image shift calibrations)
addAF		= False		# does autofocus at the start of every tilt group, increases exposure on tracking TS drastically
previewAli	= True		# adds initial dose, but makes sure start tilt image is on target (uses view image and aligns to buffer P if alignToP == True)
geoRefine	= False		# uses on-the-fly CtfFind results of first image to refine geometry before tilting (only use when CTF fits on your sample seem reliable) [NOT IN VIDEO TUTORIAL]

# Even more persistent settings
imageShiftLimit	= 15		# maximum image shift SerialEM is allowed to apply (this is a SerialEM property entry, default is 15 microns)
dataPoints	= 4		# number of recent specimen shift data points used for estimation of eucentric offset (default: 4)
alignLimit	= 0.5		# maximum shift in microns allowed for record tracking between tilts, should reduce loss of target in case of low contrast (not applied for tracking TS)
fitLimit	= 30		# geoRefine: minimum resolution in Angstroms needed for CTF fit to be considered for geoRefine
parabolTh	= 9		# geoRefine: minimum number of passable CtfFind values to fit paraboloid instead of plane 
minCounts	= 0 		# minimum mean counts per second of record image (if set > 0, tilt series branch will be aborted if mean counts are not sufficient)
splitLog	= False		# splits log file every 4 tilts (in case of >50 targets, writing to a large log file can become slow)
ignoreFirstNegShift = True	# ignore first shift on 2nd branch, which is usually very large on bad stages
slowTilt	= False		# do backlash step for all tilt angles, on bad stages large tilt steps are less accurate
taOffsetPos	= 0 		# additional tilt axis offset values applied to calculations for postitive and...
taOffsetNeg	= 0 		# ...negative branch of the tilt series (possibly useful for side-entry holder systems)

# Not changed unless problematic
extendedMdoc	= True		# saves additional info to .mdoc file
checkDewar	= True		# check if dewars are refilling before every acquisition
doCtfFind	= True		# set to False to skip CTFfind estimation (only necessary if it causes crashes) 

########## END SETTINGS ########## 

########### FUNCTIONS ###########

def checkFilling():
	filling = sem.AreDewarsFilling()
	while filling == 1:
		sem.Echo("Dewars are filling...")
		sem.Delay(60)
		filling = sem.AreDewarsFilling()

def parseTargets(targetFile):
	targets = []
	for line in targetFile:
		col = line.strip(os.linesep).split(" ")
		if col[0] == "": continue
		if line.startswith("_set") and len(col) == 4:
			if col[1] in globals():
				sem.Echo("WARNING: Read setting from tgts file and overwrite: " + col[1] + " = " + col[3])
				globals()[col[1]] = float(col[3])
			else:
				sem.Echo("WARNING: Attempted to overwrite " + col[1] + " but variable does not exist!")
		elif line.startswith("_tgt"):
			targets.append({})
		else:
			targets[-1][col[0]] = col[2]
	return targets

def geoPlane(x, a, b):
	return a * x[0] + b * x[1]

def geoPara(x, a, b, c, d, e):
	return a * x[0] + b * x[1] + c * (x[0]**2) + d * (x[1]**2) + e * x[0] * x[1]

def Tilt(tilt, sec):
	def calcSSChange(x, z0):			#x = array(tilt, n0) => needs to be one array for optimize.curve_fit()
		return x[1] * (np.cos(np.radians(x[0])) - np.cos(np.radians(x[0] - increment))) - z0 * (np.sin(np.radians(x[0])) - np.sin(np.radians(x[0] - increment)))

	def calcFocusChange(x, z0):			#x = array(tilt, n0) => needs to be one array for optimize.curve_fit()
		return z0 * (np.cos(np.radians(x[0])) - np.cos(np.radians(x[0] - increment))) + x[1] * (np.sin(np.radians(x[0])) - np.sin(np.radians(x[0] - increment)))

	sem.TiltTo(tilt)
	if tilt < startTilt:
		increment = -step
		sem.TiltBy(-step)
		sem.TiltTo(tilt)
		if (tilt - startTilt) % (2 * increment) == 0:		#figure out which section of the stack to align to
			ali = sec - 1
			newGroup = False
		else:
			ali = sec - 3
			newGroup = True
		pn = 2
	else:
		if slowTilt:						#on bad stages, better to do backlash as well to enhance accuracy
			sem.TiltBy(-step)
			sem.TiltTo(tilt)
		increment = step
		if (tilt - startTilt) % (2 * increment) == 0 or tilt - startTilt == increment:		#figure out which section of the stack to align to
			ali = sec - 1
			newGroup = False
		else:
			ali = sec - 3
			newGroup = True
		pn = 1

	sem.Delay(delayTilt)

	if zeroExpTime > 0 and tilt == startTilt:
		sem.SetExposure("R", zeroExpTime)

	for pos in range(0,len(position)):
		if pos != 0 and position[pos][pn]["skip"]: 
			sem.Echo("[" + str(pos) + "] was skipped on this branch.")
			continue
		if tilt != startTilt:
			sem.OpenOldFile(targets[pos]["tsfile"])
			sem.ReadFile(ali, "O")				#read last image of position for AlignTo
		else:
			if tgtPattern:
				targets[pos]["tsfile"] = targets[pos]["tsfile"].split(".mrc")[0] + "_" + str(round(time.time())%1000000) + ".mrc"
				sem.OpenNewFile(targets[pos]["tsfile"])
			else:
				sem.OpenNewFile(targets[pos]["tsfile"])
				sem.ReadOtherFile(0, "O", targets[pos]["tgtfile"])	#reads tgt file for first AlignTo instead

### Calculate and apply predicted shifts
		SSchange = 0 						#only apply changes if not startTilt
		focuschange = 0
		if tilt != startTilt:
			SSchange = calcSSChange([tilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])
			focuschange = calcFocusChange([tilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])

		SSYprev = position[pos][pn]["SSY"]
		SSYpred = position[pos][pn]["SSY"] + SSchange

		focuscorrection = focusSlope * (tilt - startTilt)
		position[pos][pn]["focus"] += focuscorrection
		position[pos][pn]["focus"] -= focuschange

		sem.SetDefocus(position[pos][pn]["focus"])
		sem.SetImageShift(position[pos][pn]["ISXset"], position[pos][pn]["ISYset"])
		sem.ImageShiftByMicrons(0, SSchange)

### Autofocus (optional)
		if pos == 0 and addAF and newGroup:
			sem.G(-1)
			(defocus, RepVal2) = sem.ReportAutoFocus()
			focuserror = float(defocus) - targetDefocus
			for i in range(0, len(position)):
				position[i][pn]["focus"] -= focuserror
			sem.SetDefocus(position[pos][pn]["focus"])

### Record
		if checkDewar: checkFilling()
		if beamTiltComp: 
			stig = sem.ReportObjectiveStigmator()
			sem.AdjustBeamTiltforIS()
		sem.Delay(delayIS)
		sem.R()
		sem.S()
		if beamTiltComp: 
			sem.RestoreBeamTilt()
			sem.SetObjectiveStigmator(stig[0], stig[1])

		if pos != 0: sem.LimitNextAutoAlign(alignLimit)		#gives maximum distance for AlignTo to avoid runaway tracking
		if tilt != startTilt or not tgtPattern: 
			sem.AlignTo("O")

		(bufISX, bufISY) = sem.ReportISforBufferShift()
		sem.ImageShiftByUnits(position[pos][pn]["ISXali"], position[pos][pn]["ISYali"])		#remove accumulated buffer shifts to calculate alignment to initial startTilt image

		position[pos][pn]["focus"] += focuscorrection						#remove correction or it accumulates

		(position[pos][pn]["ISXset"], position[pos][pn]["ISYset"], RepVal3, RepVal4, RepVal5, RepVal6) = sem.ReportImageShift()

		if pos == 0:						#apply measured shifts of first/tracking position to other positions
			for i in range(1, len(position)):
				position[i][pn]["ISXset"] += bufISX + position[pos][pn]["ISXali"]	#apply accumulated (stage dependent) buffer shifts of tracking TS to all targets
				position[i][pn]["ISYset"] += bufISY + position[pos][pn]["ISYali"]
				if tilt == startTilt:			#also save shifts from startTilt image for second branch since it will alignTo the startTilt image
					position[i][2]["ISXset"] += bufISX
					position[i][2]["ISYset"] += bufISY
					position[i][2]["ISXali"] += bufISX
					position[i][2]["ISYali"] += bufISY
			if tilt == startTilt:				#do not forget about 0 position
				position[0][2]["ISXset"] += bufISX
				position[0][2]["ISYset"] += bufISY

		position[pos][pn]["ISXali"] += bufISX
		position[pos][pn]["ISYali"] += bufISY
		if tilt == startTilt:					#save alignment of first tilt to tgt file for the second branch
			position[pos][2]["ISXali"] += bufISX
			position[pos][2]["ISYali"] += bufISY

		(position[pos][pn]["SSX"], position[pos][pn]["SSY"]) = sem.ReportSpecimenShift()

		sem.Echo("[" + str(pos) + "] Prediction: y = " + str(round(SSYpred, 3)) + " | z = " + str(round(position[pos][pn]["focus"], 3)) + " | z0 = " + str(round(position[pos][pn]["z0"], 3)))
		sem.Echo("[" + str(pos) + "] Reality: y = " + str(round(position[pos][pn]["SSY"], 3)))
		sem.Echo("[" + str(pos) + "] Focus change: " + str(round(focuschange, 3)) + " | Focus correction: " + str(round(focuscorrection, 3)))

### Calculate new z0

		ddy = position[pos][pn]["SSY"] - SSYprev
		if (ignoreFirstNegShift and pn == 2 and len(position[pos][pn]["shifts"]) == 0) or tilt == startTilt:		#ignore first shift to target image and ignore first shift of second branch, which is quite large on bad stages
			ddy = calcSSChange([tilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])

		position[pos][pn]["shifts"].append(ddy)
		position[pos][pn]["angles"].append(tilt)

		if len(position[pos][pn]["shifts"]) > dataPoints:
			position[pos][pn]["shifts"].pop(0)
			position[pos][pn]["angles"].pop(0)

		(position[pos][pn]["z0"], cov) = optimize.curve_fit(calcSSChange, np.vstack((position[pos][pn]["angles"], [position[pos][pn]["n0"] for i in range(0, len(position[pos][pn]["angles"]))])), position[pos][pn]["shifts"], p0=(position[pos][pn]["z0"]))
		position[pos][pn]["z0"] = position[pos][pn]["z0"][0]

		if doCtfFind:
			(cfind, RepVal2, RepVal3, RepVal4, RepVal5, cres) = sem.CtfFind("A", (targetDefocus - 2), (targetDefocus + 2))
			if geoRefine and tilt == startTilt:
				if cres < fitLimit:				# save vectors for geoRefine only if CTF fit has reasonable resolution
					geo[0].append(position[pos][pn]["SSX"])
					geo[1].append(position[pos][pn]["SSY"])
					geo[2].append(cfind)

		progress = sec * len(position) + pos + 1
		percent = round(100 * (progress / maxProgress), 1)
		bar = '#' * int(percent / 2) + '_' * (50 - int(percent / 2))
		remTime = (sem.ReportClock() - startTime) / percent * (100 - percent)
		sem.Echo("Progress: |" + bar + "| " + str(percent) + " % (" + str(int(remTime / 60)) + " min remaining)")

		if extendedMdoc:
			sem.AddToAutodoc("SpecimenShift", str(position[pos][pn]["SSX"]) + " " + str(position[pos][pn]["SSY"]))
			sem.AddToAutodoc("EucentricOffset", str(position[pos][pn]["z0"]))
			if doCtfFind:
				sem.AddToAutodoc("CtfFind", str(cfind))
			sem.WriteAutodoc()

		sem.CloseFile()

		if np.sqrt(np.array([position[pos][pn]["SSX"], position[pos][pn]["SSY"]], dtype=float).dot(np.array([position[pos][pn]["SSX"], position[pos][pn]["SSY"]], dtype=float))) > imageShiftLimit - alignLimit:
			position[pos][pn]["skip"] = True
			sem.Echo("WARNING: Target [" + str(pos) + "] is approaching the image shift limit. This branch will be aborted.")

		if minCounts > 0:
			meanCounts = sem.ReportMeanCounts()
			(expTime, RepVal2) = sem.ReportExposure("R")
			if meanCounts / expTime < minCounts:
				position[pos][pn]["skip"] = True
				sem.Echo("WARNING: Target [" + str(pos) + "] was too dark. This branch will be aborted.")	

	if zeroExpTime > 0 and tilt == startTilt:
		sem.RestoreCameraSet("R")			

######## END FUNCTIONS ########

if (maxTilt - startTilt) > 70 or (startTilt - maxTilt - step) < -70:
	print("ERROR: Maximal tilt angle too high! Stage limitations do not allow for symmetrical tilt series with these values!")
	sem.Return()

sem.SuppressReports()

sem.ReportNavItem()
navNote = sem.GetVariable("navNote")

tf = sem.DoesFileExist(navNote)
if tf == 0:
	sem.OKBox("Target file not found! Please choose the directory containing the target file! WARNING: All future target files will be searched here!")
	sem.UserSetDirectory("Please choose the directory containing the target file!")

curDir = sem.ReportDirectory()

with open(os.path.join(curDir, navNote)) as f:
	targetFile = f.readlines()

sem.SaveLogOpenNew(navNote.split("_tgts")[0])

sem.ResetClock()

targets = parseTargets(targetFile)

targetDefocus = maxDefocus						#use highest defocus for tracking TS
sem.SetTargetDefocus(targetDefocus)

sem.Echo("##### Starting new PACEtomo with parameters: #####")
sem.Echo("Start: " + str(startTilt) + " deg - Min/Max: " + str(startTilt - maxTilt) + "/" + str(startTilt + maxTilt) + "deg (" + str(step) + " deg increments)")
sem.Echo("Data points used: " + str(dataPoints))
sem.Echo("Target defocus range (min/max/step): " + str(minDefocus) + "/" + str(maxDefocus) + "/" + str(stepDefocus))
sem.Echo("Sample pretilt (rotation): " + str(pretilt) + " (" + str(rotation) + ")")
sem.Echo("Focus correction slope: " + str(focusSlope))

sem.SetProperty("ImageShiftLimit", imageShiftLimit)

sem.Echo("Realigning to target 1...")

sem.MoveToNavItem()
sem.Eucentricity(1)
sem.UpdateItemZ()

if alignToP:
	sem.V()
	sem.AlignTo("P")
else:
	sem.RealignToNavItem(1)

sem.Echo("Tilting to start tilt angle...")
#backlash correction
sem.V()
sem.Copy("A", "O")

sem.TiltBy(-2 * step)
if slowTilt: sem.TiltBy(step)
sem.TiltTo(startTilt)

sem.V()
sem.AlignTo("O")
sem.GoToLowDoseArea("R")

if not tgtPattern:
	sem.LoadNavMap("O")						#preview ali before first tilt image is taken
	sem.L()
	sem.AlignTo("O")

(ISX0, ISY0, RepVal3, RepVal4, RepVal5, RepVal6) = sem.ReportImageShift()
(SSX0, SSY0) = sem.ReportSpecimenShift()

sem.G()
focus0 = float(sem.ReportDefocus())
positionFocus = focus0 						#set maxDefocus as focus0 and add focus steps in loop
minFocus0 = focus0 - maxDefocus + minDefocus

sem.Echo("Setting up " + str(len(targets)) + " targets...")

position = []
for points in targets:
	if np.sqrt(np.array([points["SSX"], points["SSY"]], dtype=float).dot(np.array([points["SSX"], points["SSY"]], dtype=float))) > imageShiftLimit - alignLimit:
		sem.Echo("WARNING: Target [" + points["tgtfile"] + "] is too close to the image shift limit. This target will we skipped.")
		continue

	position.append([])
	position[-1].append({})

	tiltScaling = np.cos(np.radians(pretilt * np.cos(np.radians(rotation)) + startTilt)) / np.cos(np.radians(pretilt * np.cos(np.radians(rotation))))	#stretch shifts from 0 tilt to startTilt
	sem.ImageShiftByMicrons(float(points["SSX"]), float(points["SSY"]) * tiltScaling)	#apply relative shifts to find out absolute IS after realign to item
	if previewAli:							#adds initial dose, but makes sure start tilt image is on target
		if alignToP:
			sem.V()
			sem.AlignTo("P")
			sem.GoToLowDoseArea("R")
		else:
			sem.ReadOtherFile(0, "O", points["tgtfile"])	#reads tgt file for first AlignTo instead
			sem.L()
			sem.AlignTo("O")	
	(ISXset, ISYset, RepVal3, RepVal4, RepVal5, RepVal6) = sem.ReportImageShift()
	(SSX, SSY) = sem.ReportSpecimenShift()
	sem.SetImageShift(ISX0, ISY0)					#reset IS to center position	

	z0_ini = np.tan(np.radians(pretilt)) * (np.cos(np.radians(rotation)) * float(points["SSY"]) - np.sin(np.radians(rotation)) * float(points["SSX"]))
	correctedFocus = positionFocus - z0_ini * np.cos(np.radians(startTilt)) - float(points["SSY"]) * np.sin(np.radians(startTilt))

	position[-1][0]["SSX"] = float(SSX)
	position[-1][0]["SSY"] = float(SSY)
	position[-1][0]["focus"] = correctedFocus
	position[-1][0]["z0"] = z0_ini					#offset from eucentric height (will be refined during collection)
	position[-1][0]["n0"] = float(points["SSY"])			#offset from tilt axis
	position[-1][0]["shifts"] = []
	position[-1][0]["angles"] = []
	position[-1][0]["ISXset"] = float(ISXset)
	position[-1][0]["ISYset"] = float(ISYset)
	position[-1][0]["ISXali"] = 0
	position[-1][0]["ISYali"] = 0
	position[-1][0]["skip"] = False

	position[-1].append(copy.deepcopy(position[-1][0]))		#plus and minus branch start with same values
	position[-1].append(copy.deepcopy(position[-1][0]))

	position[-1][1]["n0"] -= taOffsetPos
	position[-1][2]["n0"] -= taOffsetNeg

	positionFocus += stepDefocus					#adds defocus step between targets and resets to initial defocus if minDefocus is surpassed
	if positionFocus > minFocus0: positionFocus = focus0

### TS
sem.Echo("Start tilt series...")
branchsteps = maxTilt / 2 / step
maxProgress = (maxTilt * 2 / step + 1) * len(position)
startTime = sem.ReportClock()

geo = [[], [], []]

plustilt = minustilt = startTilt
Tilt(startTilt, 0)

if geoRefine:
	if len(geo[2]) >= 3:							# if number of points > 3: fit z = a * x + b * y
		sem.Echo("Refining geometry...")
		sem.Echo(str(len(geo[2])) + " usable CtfFind results found.")
		if len(geo[2]) >= parabolTh:					# if number of points > 6: fit z = a * x + b * y + c * (x**2) + d * (y**2) + e * x * y
			sem.Echo("Fitting paraboloid...")
			geoF = geoPara
		else:							
			sem.Echo("Fitting plane...")
			geoF = geoPlane

		p, cov = optimize.curve_fit(geoF, [geo[0], geo[1]], [z - geo[2][0] for z in geo[2]])

		ss = 0
		for i in range(0, len(geo[2])):
			ss += (geo[2][i] - geo[2][0] - geoF([geo[0][i], geo[1][i]], *p))**2
		rmse = np.sqrt(ss / len(geo[2]))

		sem.Echo("Fit parameters: " + " # ".join(p.astype(str)))
		sem.Echo("RMSE: " + str(rmse))

		for pos in range(0, len(position)):				# calculate and adjust refined z0
			zs = geoF([position[pos][1]["SSX"], position[pos][1]["SSY"]], *p)
			z0_ref = position[pos][1]["z0"] + zs * np.cos(np.radians(startTilt)) + position[pos][1]["SSY"] * np.sin(np.radians(startTilt))

			position[pos][1]["z0"] = z0_ref
			position[pos][2]["z0"] = z0_ref
	else: 
		sem.Echo("WARNING: Not enough reliable CtfFind results (" + str(len(geo[2])) + ") to refine geometry. Continuing with initial geometry model.")

for i in range(0,int(branchsteps)):
	for j in range(0,2):
		plustilt += step
		sem.Echo("Tilt step " + str(i * 4 + j + 1) + " out of " + str(int(branchsteps) * 4 + 1) + " (" + str(plustilt) + " deg)...")
		Tilt(plustilt, i * 4 + j + 1)
	for j in range(0,2):
		minustilt -= step
		sem.Echo("Tilt step " + str(i * 4 + j + 3) + " out of " + str(int(branchsteps) * 4 + 1) + " (" + str(minustilt) + " deg)...")
		Tilt(minustilt, i * 4 + j + 3)
	if splitLog: sem.SaveLogOpenNew(navNote.split("_tgts")[0])

### Finish
sem.TiltTo(0)
sem.SetDefocus(focus0)
sem.ResetImageShift()
sem.CloseFile()

totalTime = round(sem.ReportClock() / 60, 1)
perTime = round(totalTime / len(position), 1)
print("##### All tilt series completed in " + str(totalTime) + " min (" + str(perTime) + " min per tilt series) #####")
sem.SaveLog()
sem.Exit()