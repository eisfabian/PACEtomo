#!Python
# ===================================================================
# Name:		PACEtomo
# Purpose:	Runs parallel dose-symmetric tilt series on many targets with geometrical predictions using the first target as tracking tilt series.
# 		Make sure to run selectTargets script first to generate compatible Navigator settings and a target file.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2021/04/16
# Revision:	v1.1
# Last Change:	2022/03/23
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

focusSlope	= 0.0		# empirical linear focus correction (obtained by linear regression of CTF fitted defoci over tilt series; microscope stage dependent)
delayIS		= 0.3 		# delay between applying image shift and Record
delayTilt	= 0 		# delay after stage tilt

# Lamella settings
pretilt		= 0		# pretilt of sample in deg e.g. after FIB milling (if milling direction is not perpendicular to the tilt axis, estimate and add rotation)
rotation	= 0		# rotation of lamella vs tilt axis in deg (should be 0 deg if lamella is oriented perpendicular to tilt axis)

# Holey support settings
tgtPattern 	= False		# use same tgt pattern on different stage positions (useful for collection on holey support film)
alignToP	= False		# use generic image in buffer P to align to

# Persistent settings
beamTiltComp	= True		# use beam tilt compensation (uses coma vs image shift calibrations)
addAF		= False		# does autofocus at the start of every tilt group, increases exposure on tracking TS drastically (did not improve defocus spread noticably during testing)
previewAli	= True		# adds initial dose, but makes sure start tilt image is on target (uses view image and aligns to buffer P if alignToP == True)

# Even more persistent settings
imageShiftLimit	= 15		# maximum image shift SerialEM is allowed to apply (this is a SerialEM property entry, default is 15 microns)
dataPoints	= 4		# number of recent specimen shift data points used for estimation of eucentric offset (default: 4)
alignLimit	= 0.5		# maximum shift in microns allowed for record tracking between tilts, should reduce loss of target in case of low contrast (not applied for tracking TS)
extendedMdoc	= True		# saves additional info to .mdoc file
splitLog	= False		# splits log file every 4 tilts (in case of >50 targets, writing to a large log file can become slow)
ignoreFirstNegShift = True	# ignore first shift on 2nd branch, which is usually very large on bad stages
slowTilt	= False		# do backlash step for all tilt angles, on bad stages large tilt steps are less accurate

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
		if line.startswith("_tgt"):
			targets.append({})
		else:
			targets[-1][col[0]] = col[2]
	return targets

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

	for pos in range(0,len(position)):

		if tilt != startTilt:
			sem.OpenOldFile(targets[pos]["tsfile"])
			sem.ReadFile(ali, "M")				#read last image of position for AlignTo
		else:
			if tgtPattern:
				targets[pos]["tsfile"] = targets[pos]["tsfile"].split(".mrc")[0] + "_" + str(round(time.time())%1000000) + ".mrc"
				sem.OpenNewFile(targets[pos]["tsfile"])
			else:
				sem.OpenNewFile(targets[pos]["tsfile"])
				sem.ReadOtherFile(0, "M", targets[pos]["tgtfile"])	#reads tgt file for first AlignTo instead

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
		checkFilling()
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
			sem.AlignTo("M")

		(bufISX, bufISY) = sem.ReportISforBufferShift()
		sem.ImageShiftByUnits(position[pos][pn]["ISXali"], position[pos][pn]["ISYali"])		#remove accumulated buffer shifts to calculate alignment to initial startTilt image

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

		sem.Echo("[" + str(pos) + "] Prediction: ### y = " + str(SSYpred) + " ### z = " + str(position[pos][pn]["focus"]) + " ### z0 = " + str(position[pos][pn]["z0"]))
		sem.Echo("[" + str(pos) + "] Reality: ###### y = " + str(position[pos][pn]["SSY"]) + " ###")
		sem.Echo("[" + str(pos) + "] Focus change: " + str(focuschange) + " ### Focus correction: " + str(focuscorrection))

		position[pos][pn]["focus"] += focuscorrection		#remove correction or it accumulates

### Calculate new z0

		ddy = position[pos][pn]["SSY"] - SSYprev
		if (ignoreFirstNegShift == True and pn == 2 and len(position[pos][pn]["shifts"]) == 0) or tilt == startTilt:		#ignore first shift to target image and ignore first shift of second branch, which is quite large on bad stages
			ddy = 0

		position[pos][pn]["shifts"].append(ddy)
		position[pos][pn]["angles"].append(tilt)

		if len(position[pos][pn]["shifts"]) > dataPoints:
			position[pos][pn]["shifts"].pop(0)
			position[pos][pn]["angles"].pop(0)


		(position[pos][pn]["z0"], cov) = optimize.curve_fit(calcSSChange, np.vstack((position[pos][pn]["angles"], [position[pos][pn]["n0"] for i in range(0, len(position[pos][pn]["angles"]))])), position[pos][pn]["shifts"], p0=(position[pos][pn]["z0"]))
		position[pos][pn]["z0"] = position[pos][pn]["z0"][0]

		if extendedMdoc == True:
			sem.AddToAutodoc("SpecimenShift", str(position[pos][pn]["SSX"]) + " " + str(position[pos][pn]["SSY"]))
			sem.AddToAutodoc("EucentricOffset", str(position[pos][pn]["z0"]))
			cfind = sem.CtfFind("A", (targetDefocus - 2), (targetDefocus + 2))[0]
			sem.AddToAutodoc("CtfFind", str(cfind))
			sem.WriteAutodoc()

		sem.CloseFile()

######## END FUNCTIONS ########

if (maxTilt - startTilt) > 70 or (startTilt - maxTilt - step) < -70:
	print("ERROR: Maximal tilt angle too high! Stage limitations do not allow for symmetrical tilt series with these values!")
	sem.Return()

sem.ReportNavItem()
navNote = sem.GetVariable("navNote")

with open(navNote) as f:
	targetFile = f.readlines()

sem.SaveLogOpenNew(navNote.split("_")[0])

sem.ResetClock()

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
	sem.LoadNavMap("M")						#preview ali before first tilt image is taken
	sem.L()
	sem.AlignTo("M")

(ISX0, ISY0, RepVal3, RepVal4, RepVal5, RepVal6) = sem.ReportImageShift()
(SSX0, SSY0) = sem.ReportSpecimenShift()

sem.G()
focus0 = float(sem.ReportDefocus())
positionFocus = focus0 						#set maxDefocus as focus0 and add focus steps in loop
minFocus0 = focus0 - maxDefocus + minDefocus

targets = parseTargets(targetFile)

sem.Echo("Setting up " + str(len(targets)) + " targets...")

position = []
for points in targets:
	position.append([])
	position[-1].append({})

	tiltScaling = np.cos(np.radians(pretilt + startTilt)) / np.cos(np.radians(pretilt))	#stretch shifts from 0 tilt to startTilt
	sem.ImageShiftByMicrons(float(points["SSX"]), float(points["SSY"]) * tiltScaling)	#apply relative shifts to find out absolute IS after realign to item
	if previewAli:							#adds initial dose, but makes sure start tilt image is on target
		if alignToP:
			sem.V()
			sem.AlignTo("P")
			sem.GoToLowDoseArea("R")
		else:
			sem.ReadOtherFile(0, "M", points["tgtfile"])	#reads tgt file for first AlignTo instead
			sem.L()
			sem.AlignTo("M")	
	(ISXset, ISYset, RepVal3, RepVal4, RepVal5, RepVal6) = sem.ReportImageShift()
	(SSX, SSY) = sem.ReportSpecimenShift()
	sem.SetImageShift(ISX0, ISY0)					#reset IS to center position	

	z0_ini = np.tan(np.radians(float(pretilt))) * (np.cos(np.radians(float(rotation))) * float(points["SSY"]) - np.sin(np.radians(float(rotation))) * float(points["SSX"]))

	position[-1][0]["SSX"] = float(SSX)
	position[-1][0]["SSY"] = float(SSY)
	position[-1][0]["focus"] = positionFocus
	position[-1][0]["z0"] = z0_ini					#offset from eucentric height (will be refined during collection)
	position[-1][0]["n0"] = float(points["SSY"])			#offset from tilt axis
	position[-1][0]["shifts"] = []
	position[-1][0]["angles"] = []
	position[-1][0]["ISXset"] = float(ISXset)
	position[-1][0]["ISYset"] = float(ISYset)
	position[-1][0]["ISXali"] = 0
	position[-1][0]["ISYali"] = 0

	position[-1].append(copy.deepcopy(position[-1][0]))		#plus and minus branch start with same values
	position[-1].append(copy.deepcopy(position[-1][0]))

	positionFocus += stepDefocus					#adds defocus step between targets and resets to initial defocus if minDefocus is surpassed
	if positionFocus > minFocus0: positionFocus = focus0

### TS
sem.Echo("Start tilt series...")
branchsteps = maxTilt / 2 / step

plustilt = minustilt = startTilt
Tilt(startTilt, 0)

for i in range(0,int(branchsteps)):
	for j in range(0,2):
		plustilt += step
		sem.Echo("Tilt step " + str(i * 4 + j + 1) + " out of " + str(int(branchsteps) * 4 + 1) + " (" + str(plustilt) + " deg)...")
		Tilt(plustilt, i * 4 + j + 1)
	for j in range(0,2):
		minustilt -= step
		sem.Echo("Tilt step " + str(i * 4 + j + 3) + " out of " + str(int(branchsteps) * 4 + 1) + " (" + str(minustilt) + " deg)...")
		Tilt(minustilt, i * 4 + j + 3)
	if splitLog: sem.SaveLogOpenNew(navNote.split("_")[0])

### Finish
sem.TiltTo(0)
sem.SetDefocus(focus0)
sem.ResetImageShift()
sem.CloseFile()

sem.ReportClock()
print("##### Tilt series completed #####")
sem.SaveLog()
sem.Exit()
