#!Python
# ===================================================================
# Name:		PACEtomo_selectTargets
# Purpose:	Selects targets, saves maps and writes text file with coords for PACEtomo.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2021/04/19
# Revision:	v1.4
# Last Change:	2022/05/25: beam diameter visualization
# ===================================================================

import serialem as sem
import os
import numpy as np
import scipy as sp
import scipy.optimize

# Please set the appropiate tilt axis offset to improve performance!

############ SETTINGS ############ 

delaytime 	= 10		# time in seconds to apply image shift before next image is taken
useSearch 	= False		# use Search mode instead of View mode to find targets by dragging [NOT IN VIDEO TUTORIAL]

targetByShift 	= False		# ask to enter image shifts instead of dragging manually
targetPattern 	= False		# regular pattern of targets (holey support film)
alignToP 	= False		# refine vectors by aligning to hole reference in buffer P
size 		= 1 		# size of collection pattern (1: 3x3, 2: 5x5, 3: 7x7, ...)

vecA 		= (2.4, 1.1)	# specimen shift in microns to neighbouring hole
vecB 		= (-vecA[1], vecA[0])

drawBeam 	= True		# draws navigator item representing beam diameter (only on Thermo Scientific microscopes, needs ReportIlluminatedArea) [NOT IN VIDEO TUTORIAL]
beamDiameter 	= 0 		# beam diameter in microns (if 0, ReportIlluminatedArea will be used, which is only available on some microscopes)
maxTilt 	= 60		# tilt angle to calculate stretching of beam perpendicular to tilt axis

########## END SETTINGS ########## 

########### FUNCTIONS ###########

##########
# Source: https://stackoverflow.com/a/52062369
def angles_in_ellipse(num, a, b):
	assert(num > 0)
	assert(a < b)
	angles = 2 * np.pi * np.arange(num) / num
	if a != b:
	        e2 = (1.0 - a ** 2.0 / b ** 2.0)
        	tot_size = sp.special.ellipeinc(2.0 * np.pi, e2)
        	arc_size = tot_size / num
        	arcs = np.arange(num) * arc_size
        	res = sp.optimize.root(
        		lambda x: (sp.special.ellipeinc(x, e2) - arcs), angles)
        	angles = res.x 
	return angles
##########

def parseNav(navFile):
	with open(navFile) as f:
		navContent = f.readlines()
	header = []
	items = []
	newItem = {}
	index = 0

	for line in navContent:
		col = line.rstrip().split(" ")
		if col[0] == "": 
			if "Item" in newItem.keys():
				items.append(newItem)
				newItem = {}
				continue
			else:
				continue
		if line.startswith("[Item"):
			index += 1
			newItem = {"index": index, "Item": col[2].strip("]")}
		elif "Item" in newItem.keys():
			newItem[col[0]] = [val for val in col[2:]]
		else:
			header.append(line)
	if "Item" in newItem.keys():	#append last target
		items.append(newItem)
	return header, items

def writeNav(header, items, filename):
	text = ""
	for line in header:
		text += line
	text += os.linesep
	for item in items:
		text += "[Item = " + item["Item"] + "]" + os.linesep
		item.pop("Item")
		item.pop("index")
		for key, attr in item.items():
			text += key + " = "
			for val in attr:
				text += str(val) + " "
			text += os.linesep
		text += os.linesep
	with open(filename, "w", newline="") as f:
		f.write(text)
	return

######## END FUNCTIONS ########

sem.Pause("Please make sure that 'Move stage for big mouse shifts' is unchecked!")

groupInfo = sem.ReportGroupStatus()
pointRefine = 0
if not targetPattern and not targetByShift and groupInfo[1] > 0:
	pointRefine = sem.YesNoBox("The selected navigator item is part of a group. Do you want to use the points of this group as initial target coordinates to be refined? (NOT IN VIDEO TUTORIAL)")
	if pointRefine == 1:
		groupID = groupInfo[1]
		groupPoints = groupInfo[2]
		firstPoint = sem.ReportNavItem()
		navLabel = sem.GetVariable("navLabel")
		sem.Pause("Please make sure that the first point of the group is selected! (Selected point: " + navLabel + ", points in group: " + str(int(groupPoints)) + ")")
		coordsRefine = []
		coordsRefine.append([0, 0])
		for i in range(int(firstPoint[0]) + 1, int(firstPoint[0] + groupPoints)):
			point = sem.ReportOtherItem(i)
			coordsRefine.append([-point[1] + firstPoint[1], -point[2] + firstPoint[2]])

sem.UserSetDirectory("Please choose a directory for saving targets and tilt series!")

sem.EnterString("userName","Please provide a rootname for the PACE-tomo collection area!")
userName = sem.GetVariable("userName")

sem.OpenTextFile("1", "W", 0, userName + "_tgts.txt")

imageShiftLimit = sem.ReportProperty("ImageShiftLimit")

#align center
sem.ResetImageShift()
if pointRefine == 1:
	sem.MoveToNavItem(int(firstPoint[0]))
elif alignToP:			#center hole for center of tgtPattern
	sem.V()
	sem.AlignTo("P")

userInput = 0
while userInput == 0:
	if useSearch: 
		sem.Search()
	else:
		sem.V()
	sem.OKBox("Please center your target by dragging the image using the right mouse button! (Delay: " + str(delaytime) + " s)")
	sem.Delay(delaytime)
	userConfirm = sem.YesNoBox("Do you want to take a preview image here?")
	if userConfirm == 1:
		sem.L()
		userInput = sem.YesNoBox("Do you want to use the current image and coordinates as target 1 (tracking target)? If you choose no, a view image is taken and you can drag to the target before another preview image is taken!")

#save map
sem.OpenNewFile(userName + "_tgt_001.mrc")
sem.S("A")
mapIndex = sem.NewMap(0, userName + "_tgts.txt")
sem.SetItemAcquire(int(mapIndex))

#save coords
(ISX0, ISY0, RepVal3, RepVal4, RepVal5, RepVal6) = sem.ReportImageShift()
(SSX0, SSY0) = sem.ReportSpecimenShift()
(stageX, stageY, stageZ) = sem.ReportStageXYZ()
stageX -= SSX0
stageY -= SSY0

sem.WriteLineToFile("1", "_tgt = 001")
sem.WriteLineToFile("1", "tgtfile = " + userName + "_tgt_001.mrc")
sem.WriteLineToFile("1", "tsfile = " + userName + "_ts_001.mrc")
sem.WriteLineToFile("1", "map = " + str(mapIndex))
sem.WriteLineToFile("1", "stageX = " + str(stageX))
sem.WriteLineToFile("1", "stageY = " + str(stageY))
sem.WriteLineToFile("1", "SSX = 0")
sem.WriteLineToFile("1", "SSY = 0")
sem.WriteLineToFile("1", "")

sem.CloseFile()

sem.Echo("Target 001 (" + userName + "_tgt_001.mrc) with image shifts 0, 0 was added.")

#create navigator entry to draw beam diameter
if drawBeam:
	if beamDiameter == 0:
		beamR = sem.ReportIlluminatedArea() * 100 / 2
	else:
		beamR = beamDiameter / 2
	a = 1 * beamR
	b = beamR / np.cos(np.radians(maxTilt))
	n = 32
	phi = angles_in_ellipse(n, a, b)

	sem.SaveNavigator()
	navFile = sem.ReportNavFile()
	navHeader, navItems = parseNav(navFile)

	ptsX = (a * np.cos(phi) + stageX).round(3)
	ptsX = np.append(ptsX, ptsX[0])
	ptsY = (b * np.sin(phi) + stageY).round(3)
	ptsY = np.append(ptsY, ptsY[0])

	navItems.append({'index': int(navItems[-1]["index"]) + 1, 'Item': str(int(navItems[-1]["Item"]) + 1), 'Color': ['3'], 'StageXYZ': [str(stageX), str(stageY), str(stageZ)], 'NumPts': [str(n+1)], 'Regis': ['1'], 'Type': ['1'], 'Note': ['beam_diameter'], 'PtsX': ptsX, 'PtsY': ptsY})

	writeNav(navHeader, navItems, navFile)
	sem.ReadNavFile(navFile)

#make view map tor realign to item
sem.V()
sem.OpenNewFile(userName + "_tgt_001_view.mrc")
sem.S("A")
sem.NewMap()
sem.CloseFile()

if targetPattern:
	if alignToP:			# refine grid vectors by aligning to hole reference in P
		sem.GoToLowDoseArea("R")
		sem.Echo("Vector A: " + str(vecA))

		shiftx = size * vecA[0]
		shifty = size * vecA[1]
		sem.ImageShiftByMicrons(shiftx, shifty)

		sem.V()
		sem.AlignTo("P")
		sem.GoToLowDoseArea("R")

		(SSX, SSY) = sem.ReportSpecimenShift()
		SSX -= SSX0
		SSY -= SSY0		

		vecA = (round(SSX / size, 4), round(SSY / size, 4))

		sem.Echo("Refined vector A: " + str(vecA))

		sem.ImageShiftByMicrons(-SSX, -SSY)		# reset IS to center position

		sem.Echo("Vector B: " + str(vecB))

		shiftx = size * vecB[0]
		shifty = size * vecB[1]
		sem.ImageShiftByMicrons(shiftx, shifty)

		sem.V()
		sem.AlignTo("P")
		sem.GoToLowDoseArea("R")

		(SSX, SSY) = sem.ReportSpecimenShift()
		SSX -= SSX0
		SSY -= SSY0		

		vecB = (round(SSX / size, 4), round(SSY / size, 4))

		sem.Echo("Refined vector B: " + str(vecB))


	targetNo = 1
	for i in range(-size,size+1):
		for j in range(-size,size+1):
			if i == j == 0: continue

			targetNo += 1

			SSX = i * vecA[0] + j * vecB[0]
			SSY = i * vecA[1] + j * vecB[1]

			sem.WriteLineToFile("1", "_tgt = " + str(targetNo).zfill(3))
			sem.WriteLineToFile("1", "tsfile = " + userName + "_ts_" + str(targetNo).zfill(3) + ".mrc")
			sem.WriteLineToFile("1", "SSX = " + str(SSX))
			sem.WriteLineToFile("1", "SSY = " + str(SSY))
			sem.WriteLineToFile("1", "")

			sem.CloseFile()

			sem.Echo("Target " + str(targetNo).zfill(3) + " (" + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc) with image shifts " + str(SSX) + ", " + str(SSY) + " was added.")

else:					#loop over other targets
	targetNo = 1
	pointNo = 1 			# counter for points in group
	addTargets = 1
	while addTargets == 1:
		userInput = 0
		userSkip = 1		# initialize as 1 to not apply shifts in case of redo
		while userInput == 0:
			if targetByShift or (pointRefine == 1 and userSkip == 1):
				sem.GoToLowDoseArea("R")
				if targetByShift:
					shiftx = sem.EnterDefaultedNumber(0, 1, "Enter X shift:")
					shifty = sem.EnterDefaultedNumber(0, 1, "Enter Y shift:")
				else:
					sem.SetImageShift(ISX0, ISY0)
					shiftx, shifty = coordsRefine[pointNo][0], coordsRefine[pointNo][1]
				sem.ImageShiftByMicrons(shiftx, shifty)
			if useSearch: 
				sem.Search()
			else:
				sem.V()
			sem.OKBox("Please center your target by dragging the image using the right mouse button! (Delay: " + str(delaytime) + " s)")
			sem.Delay(delaytime)
			userConfirm = sem.YesNoBox("Do you want to take a preview image here?")
			if userConfirm == 1:
				sem.L()
				userRefine = sem.YesNoBox("Do you want to refine the position at this mag?")
				if userRefine == 1:
					sem.OKBox("Please center your target by dragging the image using the right mouse button! (Delay: " + str(delaytime) + " s)")
					sem.Delay(delaytime)
					sem.L()
				userInput = sem.YesNoBox("Do you want to use the current image and coordinates as target position? If you choose no, a new view image is taken to align your target!")
			if pointRefine == 1 and userInput == 0:
				userSkip = sem.YesNoBox("Do you want to skip this point of the group?")
				if userSkip == 1:
					pointNo += 1
					if pointNo >= groupPoints:		# disable pointRefine when last point of a group is skipped
						pointRefine = 0
						addTargets = sem.YesNoBox("All points of the selected group have been viewed. Do you want to add another target manually?")
						if addTargets == 0: break
		if addTargets == 0: break

		(SSX, SSY) = sem.ReportSpecimenShift()
		SSX -= SSX0
		SSY -= SSY0

		if np.sqrt(np.array([SSX, SSY]).dot(np.array([SSX, SSY]))) > imageShiftLimit - 0.5:
			userAbort = sem.YesNoBox("WARNING: This target is close to the image shift limit of SerialEM (" + str(imageShiftLimit) + " microns). It will probably be dropped during acquisition. Do you still want to save this target?")
			if userAbort == 0: continue

		targetNo += 1
		pointNo += 1

#save target
		sem.OpenNewFile(userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc")
		sem.S("A")
		mapIndex = sem.NewMap()

		sem.WriteLineToFile("1", "_tgt = " + str(targetNo).zfill(3))
		sem.WriteLineToFile("1", "tgtfile = " + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc")
		sem.WriteLineToFile("1", "tsfile = " + userName + "_ts_" + str(targetNo).zfill(3) + ".mrc")
		sem.WriteLineToFile("1", "SSX = " + str(SSX))
		sem.WriteLineToFile("1", "SSY = " + str(SSY))
		sem.WriteLineToFile("1", "")

		sem.CloseFile()

		sem.Echo("Target " + str(targetNo).zfill(3) + " (" + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc) with image shifts " + str(round(SSX, 3)) + ", " + str(round(SSY, 3)) + " was added.")

		#create navigator entry to draw beam diameter
		if drawBeam:
			(RepVal1, RepVal2, RepVal3, RepVal4, ISX, ISY) = sem.ReportImageShift()
			sem.SaveNavigator()
			navHeader, navItems = parseNav(navFile)

			ptsX = (a * np.cos(phi) + stageX + ISX).round(3)
			ptsX = np.append(ptsX, ptsX[0])
			ptsY = (b * np.sin(phi) + stageY + ISY).round(3)
			ptsY = np.append(ptsY, ptsY[0])

			navItems.append({'index': int(navItems[-1]["index"]) + 1, 'Item': str(int(navItems[-1]["Item"]) + 1), 'Color': ['3'], 'StageXYZ': [str(stageX + ISX), str(stageY + ISY), str(stageZ)], 'NumPts': [str(n+1)], 'Regis': ['1'], 'Type': ['1'], 'Note': ['beam_diameter'], 'PtsX': ptsX, 'PtsY': ptsY})

			writeNav(navHeader, navItems, navFile)
			sem.ReadNavFile(navFile)

		if pointRefine == 1 and pointNo >= groupPoints:
			addTargets = sem.YesNoBox("All points of the selected group have been viewed. Do you want to add another target manually?")
			pointRefine = 0
		else:
			addTargets = sem.YesNoBox("Do you want to add another target?")


sem.SetImageShift(0,0)

if drawBeam:
	sem.SaveNavigator()
	navHeader, navItems = parseNav(navFile)
	cleanItems = []
	for item in navItems:
		if "Note" in item.keys() and item["Note"][0] == "beam_diameter": continue
		cleanItems.append(item)
	writeNav(navHeader, cleanItems, navFile)
	sem.ReadNavFile(navFile)

sem.OKBox("Target selection completed! " + str(targetNo) + " targets were selected.")
