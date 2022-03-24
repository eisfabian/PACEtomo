#!Python
# ===================================================================
# Name:		PACEtomo_selectTargets
# Purpose:	Selects targets, saves maps and writes text file with coords for PACEtomo.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2021/04/19
# Revision:	v1.1
# Last Change:	2022/03/23
# ===================================================================

import serialem as sem

# Please set the appropiate tilt axis offset to improve performance!

############ SETTINGS ############ 

delaytime = 5		# time in seconds to apply image shift before next image is taken
targetByShift = True	# ask to enter image shifts instead of dragging manually
targetPattern = True	# regular pattern of targets (holey support film)
alignToP = True		# refine vectors by aligning to hole reference in buffer P
size = 1 		# size of collection pattern (1: 3x3, 2: 5x5, 3: 7x7, ...)

vecA = (2.4, 1.1)	# specimen shift in microns to neighbouring hole
vecB = (-vecA[1], vecA[0])

########## END SETTINGS ########## 

sem.Pause("Please make sure that 'Move stage for big mouse shifts' is unchecked!")

sem.UserSetDirectory("Please choose a directory for saving targets and tilt series!")

sem.EnterString("userName","Please provide a rootname for the PACE-tomo collection area!")
userName = sem.GetVariable("userName")

sem.OpenTextFile("1", "W", 0, userName + "_tgts.txt")
#align center
sem.ResetImageShift()
if alignToP:			#center hole for center of tgtPattern
	sem.V()
	sem.AlignTo("P")

userInput = 0
while userInput == 0:
	sem.V()
	sem.Delay(delaytime)
	userConfirm = sem.YesNoBox("Do you want to take a preview image here?")
	if userConfirm == 1:
		sem.L()
		userInput = sem.YesNoBox("Do you want to use the current image and coordinates as target 1 (tracking target)? If you choose no, a view image is taken and you have several seconds to apply shifts, before another preview image is taken!")

#save map
sem.OpenNewFile(userName + "_tgt_001.mrc")
sem.S("A")
mapIndex = sem.NewMap(0, userName + "_tgts.txt")
sem.SetItemAcquire(int(mapIndex))

#save coords
(SSX0, SSY0) = sem.ReportSpecimenShift()
(stageX, stageY, stageZ) = sem.ReportStageXYZ()
(RepVal1, RepVal2, RepVal3, RepVal4, ISX, ISY) = sem.ReportImageShift()
stageX += ISX
stageY += ISY

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

#make view map tor realign to item
sem.V()
sem.OpenNewFile(userName + "_tgt_001_view.mrc")
sem.S("A")
sem.NewMap()
sem.CloseFile()

if targetPattern:
	if alignToP:			#refine grid vectors by aligning to hole reference in P
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

else:					#loop over other targets
	targetNo = 1
	addTargets = 1
	while addTargets == 1:
		userInput = 0
		while userInput == 0:
			if targetByShift:
				sem.GoToLowDoseArea("R")
				shiftx = sem.EnterDefaultedNumber(0, 1, "Enter X shift:")
				shifty = sem.EnterDefaultedNumber(0, 1, "Enter Y shift:")
				sem.ImageShiftByMicrons(shiftx, shifty)
			sem.V()
			sem.Delay(delaytime)
			userConfirm = sem.YesNoBox("Do you want to take a preview image here?")
			if userConfirm == 1:
				sem.L()
				userInput = sem.YesNoBox("Do you want to use the current image and coordinates as target position? If you choose no, a view image is taken and you have several seconds to apply shifts, before another preview image is taken!")
		targetNo += 1
#save map
		sem.OpenNewFile(userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc")
		sem.S("A")
		mapIndex = sem.NewMap()

#save coords
		(SSX, SSY) = sem.ReportSpecimenShift()
		SSX -= SSX0
		SSY -= SSY0

		sem.WriteLineToFile("1", "_tgt = " + str(targetNo).zfill(3))
		sem.WriteLineToFile("1", "tgtfile = " + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc")
		sem.WriteLineToFile("1", "tsfile = " + userName + "_ts_" + str(targetNo).zfill(3) + ".mrc")
		sem.WriteLineToFile("1", "stageX = " + str(stageX))
		sem.WriteLineToFile("1", "stageY = " + str(stageY))
		sem.WriteLineToFile("1", "SSX = " + str(SSX))
		sem.WriteLineToFile("1", "SSY = " + str(SSY))
		sem.WriteLineToFile("1", "")

		sem.CloseFile()

		addTargets = sem.YesNoBox("Do you want to add another target?")


sem.SetImageShift(0,0)
sem.OKBox("Target selection completed!")
