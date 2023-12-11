#!Python
# ===================================================================
#ScriptName	PACEtomo_targetsFromMontage
# Author:	Fabian Eisenstein
# Created:	2022/12/09
# Revision:	v0.13
# Last Change:	2023/06/21: added check for dummy property (requires >June2023), used SkipAcquiringNavItem, fixed navigor save popup
#		2023/06/07: added sanity check for template pixel sizes, limited binning to multiple of 2
#		2023/05/15: added new EndAcquireAtItems command
#		2023/02/28: automatically determine virt map binning
#		2023/02/22: added image shift limit warning
#		2023/02/20: added counter when tgts file exists, added binfactor for gatan (non-square) camera, added SerialEM error to skip maps during Acquire at Items
#		2023/02/13: added check for drawnID to enable Acquire at Items, added padding for points at edge of montage, fixed labels
#		2023/02/01: added check for DUMMY version
#		2023/01/31: added samePosID for view map, added view map for all targets
#		2023/01/23: added nav file manipulation
#		2023/01/20: fixed view map pixel size
#		2023/01/19: added ResetImageShift before stage move, added option to manually subtract tilt axis offset, added considerations for view pixel size
#		2023/01/17: added options to run in DUMMY version, added nav map in montage mag
#		2023/01/12: fixes after F200 test
#		2023/01/11: switched to mrcfile, added resize, fixes after Krios test

# Take a montage of your target area. 
# Select targets using the navigator with "Add Points" (points have to be part of the same group).
# Drag the tracking target position to the top of the group.
# Select one of the points and run the script.
# (Make sure drawing of points is not disabled in Navigator.)
# You can run the script in the DUMMY version. Just ignore the error messages. It does need a "template" map for Record and View mode in that case. These should have the Navigator Notes "Template Preview" and "Template View", respectively.
# You can also run the script using Acquire at Items. Just check "Skip stage move to item if possible" in the Task-related options (Acquire at Items window).

##### SETTINGS #####

noUI 	= False	# set to True to avoid folder/name selection (e.g. to run the script in a Acquire at Items routine), it will use the label of the montage as name template
prefix	= "pos"	# prefix for name when running noUI

### END SETTINGS ###

import serialem as sem
import os
import copy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import mrcfile
from skimage.transform import resize

### FUNCTIONS ###

def WriteTargets(targetFile, targets, savedRun=False, resume={"sec": 0, "pos": 0}, settings={}):
	output = ""
	if settings != {}:
		for key, val in settings.items():
			if val != "":
				output += "_set " + key + " = " + str(val) + "\n"
		output += "\n"
	if resume["sec"] > 0 or resume["pos"] > 0:
		output += "_spos = " + str(resume["sec"]) + "," + str(resume["pos"]) + "\n" * 2
	for pos in range(len(targets)):
		output += "_tgt = " + str(pos + 1).zfill(3) + "\n"
		for key in targets[pos].keys():
			output += key + " = " + str(targets[pos][key]) + "\n"
		if savedRun:
			output += "_pbr" + "\n"
			for key in savedRun[pos][0].keys():
				output += key + " = " + str(savedRun[pos][0][key]) + "\n"
			output += "_nbr" + "\n"
			for key in savedRun[pos][1].keys():
				output += key + " = " + str(savedRun[pos][1][key]) + "\n"		
		output += "\n"
	with open(targetFile, "w") as f:
		f.write(output)

def WriteMrc(outfilename, image, pixSize):
	with mrcfile.new(os.path.join(curDir, outfilename), overwrite=True) as mrc:
		mrc.set_data(image)
		mrc.voxel_size = (pixSize, pixSize, pixSize)
		mrc.update_header_from_data()
		mrc.update_header_stats()

def ParseNav(navFile):
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

def WriteNav(header, items, filename):
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

def CropImage(image, coords, fov):
	imageCrop = image[max(0, int(coords[0] - fov[0] / 2)):min(image.shape[0], int(coords[0] + fov[0] / 2)), max(0, int(coords[1] - fov[1] / 2)):min(image.shape[1], int(coords[1] + fov[1] / 2))]
	if imageCrop.shape != fov:
		mean = np.mean(imageCrop)
		if imageCrop.shape[0] < fov[0]:
			padding = np.full((int(fov[0] - imageCrop.shape[0]), imageCrop.shape[1]), mean)
			if int(coords[0] - fov[0] / 2) < 0:
				imageCrop = np.concatenate((padding, imageCrop), axis=0)
			if int(coords[0] + fov[0] / 2) > image.shape[0]:
				imageCrop = np.concatenate((imageCrop, padding), axis=0)
		if imageCrop.shape[1] < fov[1]:
			padding = np.full((imageCrop.shape[0], int(fov[1] - imageCrop.shape[1])), mean)
			if int(coords[1] - fov[1] / 2) < 0:
				imageCrop = np.concatenate((padding, imageCrop), axis=1)
			if int(coords[1] + fov[1] / 2) > image.shape[1]:
				imageCrop = np.concatenate((imageCrop, padding), axis=1)		
		sem.Echo("WARNING: Target position is close to the edge of the map and was padded.")
	return imageCrop

### END FUNCTIONS ###

s2ssMatrix = None
dummy = False
if sem.ReportProperty("DummyInstance") == 1:
	sem.Echo("WARNING: You are using SerialEM in DUMMY mode. Make sure a View image and Preview image template are available in the Navigator!")
	dummy = True

# Get nav item info

navID = int(sem.ReportNavItem()[0])
groupInfo = sem.ReportGroupStatus()
groupID = int(groupInfo[1])
drawnID = int(sem.NavIndexItemDrawnOn(navID))

if drawnID == 0:
	sem.SkipAcquiringNavItem()
	sem.EndAcquireAtItems()	# only available from 12-May-23
	sem.Echo("NOTE: This nav item is not a point on a map and will be ignored.")
	sem.Exit()

# Get user input

if not noUI:
	sem.UserSetDirectory("Please choose a directory for saving targets and tilt series!")

	sem.EnterString("userName","Please provide a rootname for the PACE-tomo collection area!")
	userName = sem.GetVariable("userName")
else:
	sem.ReportOtherItem(drawnID)
	navLabel = sem.GetVariable("navLabel")
	userName = prefix + navLabel.strip().replace(" ", "_")

curDir = sem.ReportDirectory()
tgtsFilePath = os.path.join(curDir, userName + "_tgts.txt")

# Make sure tgts file is unique
counter = 1
while os.path.exists(tgtsFilePath):
	counter += 1
	tgtsFilePath = os.path.join(curDir, userName + str(counter) + "_tgts.txt") 
if counter > 1:
	userName = userName + str(counter)

# Make sure nav contains template images

prevID = int(sem.NavIndexWithNote("Template Preview"))
viewID = int(sem.NavIndexWithNote("Template View"))
if prevID == 0 or viewID == 0:
	if not dummy:
		sem.SetColumnOrGunValve(0)
		sem.AllowFileOverwrite(1)
		if prevID == 0:
			sem.OpenNewFile("template_preview.mrc")
			sem.L()
			sem.S()
			prevID = int(sem.NewMap(0, "Template Preview"))
			sem.ChangeItemLabel(prevID, "TP")
			sem.CloseFile()
			if s2ssMatrix is None:
				s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
		if viewID == 0:
			sem.OpenNewFile("template_view.mrc")
			sem.V()
			sem.S()
			viewID = int(sem.NewMap(0, "Template View"))
			sem.ChangeItemLabel(viewID, "TV")
			sem.CloseFile()
		sem.SetColumnOrGunValve(1)
		sem.AllowFileOverwrite(0)
	else:
		sem.OKBox("ERROR: Template for View or Preview image was not found!")
		sem.Echo("ERROR: Template for View or Preview image was not found!")
		sem.Echo("If you have a View and Preview image in the Navigator that you want to use as template, please change the Navigator Note to 'Template View' and 'Template Preview', respectively.")
		sem.Exit()

# Load templates to get pixel sizes

sem.LoadOtherMap(prevID, "A")
imgProp = sem.ImageProperties("A")
camX = imgProp[0] * imgProp[2]
camY = imgProp[1] * imgProp[2]
recordPixSize = imgProp[4] * 10 / imgProp[2]

if camX != camY: 		# most non-gatan cameras are square
	binFactor = 2 		# gatan cameras need the extra binning factor because SR is bin 1 
else:
	binFactor = 1

sem.LoadOtherMap(viewID, "A")
imgProp = sem.ImageProperties("A")
viewPixSize = imgProp[4] * 10 / imgProp[2]

if recordPixSize == viewPixSize:
	sem.OKBox("ERROR: Template for View or Preview show the same pixel size! Make sure the maps could be loaded properly!")
	sem.Echo("ERROR: Template for View or Preview show the same pixel size! Make sure the maps could be loaded properly!")
	sem.Exit()	

# Make sure s2s matrix is defined

if not dummy and s2ssMatrix is None:
	sem.GoToLowDoseArea("R")
	s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))

# Load montage

sem.LoadOtherMap(drawnID, "A")

buffer, *_ = sem.ReportCurrentBuffer()

imgProp = sem.ImageProperties(buffer)
x = imgProp[0]
y = imgProp[1]
pixSize = imgProp[4] * 10

image = np.asarray(sem.bufferImage(buffer), dtype=float)

# Get point coords

sem.GetNavGroupStageCoords(groupID, "groupStageX", "groupStageY", "groupStageZ")
sem.GetNavGroupImageCoords(groupID, buffer, "groupImageX", "groupImageY")

groupStageX = np.array(sem.GetVariable("groupStageX").split(), dtype=float)
groupStageY = np.array(sem.GetVariable("groupStageY").split(), dtype=float)
stageZ = np.array(sem.GetVariable("groupStageZ").split(), dtype=float)[0]
groupImageX = np.array(sem.GetVariable("groupImageX").split(), dtype=float)
groupImageY = np.array(sem.GetVariable("groupImageY").split(), dtype=float)

# Calculate field of view

fov_recX = int(camX * recordPixSize / pixSize)
fov_recY = int(camY * recordPixSize / pixSize)
fov_viewX = int(camX * viewPixSize / pixSize)
fov_viewY = int(camY * viewPixSize / pixSize)

out_bin_rec = min(8, int(pixSize / recordPixSize)) 	# binning factor of the created virtual map for Record mode (needs to be int)
if out_bin_rec < 8:
	if out_bin_rec < 4:
		if out_bin_rec > 1:
			out_bin_rec = 2
	else:
		out_bin_rec = 4
out_bin_view = min(8, int(pixSize / viewPixSize))	# binning factor of the created virtual map for View mode (needs to be int)
if out_bin_view < 8:
	if out_bin_view < 4:
		if out_bin_view > 1:
			out_bin_view = 2
	else:
		out_bin_view = 4
out_recX = int(camX / out_bin_rec)
out_recY = int(camY / out_bin_rec)
out_viewX = int(camX / out_bin_view)
out_viewY = int(camY / out_bin_view)

# Crop out virt maps

tgtImages = []
targets = []
for i in range(groupStageX.size):
	sem.Echo("Creating virtual map of point " + str(i + 1) + "...")

	imageCrop = CropImage(image, (groupImageY[i], groupImageX[i]), (fov_recY, fov_recX))
	imageProc = np.flip(resize(imageCrop, (out_recY, out_recX), preserve_range=True, anti_aliasing=True).astype(np.float32), axis=0)
	tgtImages.append(imageProc)
	WriteMrc(userName + "_tgt_" + str(i + 1).zfill(3) + ".mrc", tgtImages[i], recordPixSize * out_bin_rec)

	imageCrop = CropImage(image, (groupImageY[i], groupImageX[i]), (fov_viewY, fov_viewX))
	viewImageProc = np.flip(resize(imageCrop, (out_viewY, out_viewX), preserve_range=True, anti_aliasing=True).astype(np.float32), axis=0)
	WriteMrc(userName + "_tgt_" + str(i + 1).zfill(3) + "_view.mrc", viewImageProc, viewPixSize * out_bin_view)

	if not dummy:
		SSX, SSY = s2ssMatrix @ np.array([groupStageX[i] - groupStageX[0], groupStageY[i] - groupStageY[0]])
		if SSX > 15 or SSY > 15:
			sem.Echo("WARNING: Point " + str(i + 1) + " requires image shifts (" + str(round(SSX, 1)) + "|" + str(round(SSY, 1)) + ") beyond the default image shift limit (15)!")
		targets.append({"tgtfile": userName + "_tgt_" + str(i + 1).zfill(3) + ".mrc", "tsfile": userName + "_ts_" + str(i + 1).zfill(3) + ".mrc", "viewfile": userName + "_tgt_" + str(i + 1).zfill(3) + "_view.mrc", "SSX": SSX, "SSY": SSY, "stageX": groupStageX[i], "stageY": groupStageY[i], "skip": "False"})
	else:	# leave out SS coords and calculate from stage coords on the fly in PACEtomo.py
		targets.append({"tgtfile": userName + "_tgt_" + str(i + 1).zfill(3) + ".mrc", "tsfile": userName + "_ts_" + str(i + 1).zfill(3) + ".mrc", "viewfile": userName + "_tgt_" + str(i + 1).zfill(3) + "_view.mrc", "stageX": groupStageX[i], "stageY": groupStageY[i], "skip": "False"})

# Load nav file

sem.SaveNavigator()
navFile = sem.ReportNavFile()
navHeader, navItems = ParseNav(navFile)

templatePrev = copy.deepcopy(navItems[prevID - 1])
templateView = copy.deepcopy(navItems[viewID - 1])
templatePrev.pop("RawStageXY")
templatePrev.pop("SamePosId")
templateView.pop("RawStageXY")
templateView.pop("SamePosId")

# Determine image dimenstions for polygons

ptsDX_rec = np.array(np.array(templatePrev["PtsX"], dtype=float) - float(templatePrev["StageXYZ"][0]))
ptsDY_rec = np.array(np.array(templatePrev["PtsY"], dtype=float) - float(templatePrev["StageXYZ"][1]))
ptsDX_view = np.array(np.array(templateView["PtsX"], dtype=float) - float(templateView["StageXYZ"][0]))
ptsDY_view = np.array(np.array(templateView["PtsY"], dtype=float) - float(templateView["StageXYZ"][1]))

# Add maps

newNavIDs = []
for i, tgt in enumerate(targets):
	tgtItem = copy.deepcopy(templatePrev)

	# Create view map first to get mapID
	viewItem = copy.deepcopy(templateView)
	viewItem["Item"] = "v" + str(i + 1).zfill(2)
	viewItem["StageXYZ"] = [tgt["stageX"], tgt["stageY"], stageZ]
	viewItem["PtsX"] = ptsDX_view + tgt["stageX"]
	viewItem["PtsY"] = ptsDY_view + tgt["stageY"]
	viewItem["MapFile"] = [os.path.join(curDir, os.path.splitext(tgt["tgtfile"])[0] + "_view.mrc")]
	viewItem["Note"] = [os.path.splitext(tgt["tgtfile"])[0] + "_view.mrc"]
	viewItem["MapBinning"] = [out_bin_view * binFactor]
	viewItem["MapMinMaxScale"] = [np.min(viewImageProc), np.max(viewImageProc)]

	# Create new map ID
	uniqueID = int(sem.GetUniqueNavID())
	while uniqueID in newNavIDs:
		uniqueID = int(sem.GetUniqueNavID())
	newNavIDs.append(uniqueID)

	newViewID = uniqueID

	viewItem["MapID"] = [newViewID]
	viewItem["RealignedID"] = [drawnID]
	viewItem["SamePosId"] = [newViewID]

	# Create target map
	tgtItem["Item"] = str(i + 1).zfill(3)
	tgtItem["StageXYZ"] = [tgt["stageX"], tgt["stageY"], stageZ]
	tgtItem["PtsX"] = ptsDX_rec + tgt["stageX"]
	tgtItem["PtsY"] = ptsDY_rec + tgt["stageY"]
	tgtItem["MapFile"] = [os.path.join(curDir, tgt["tgtfile"])]
	tgtItem["MapBinning"] = [out_bin_rec * binFactor]
	tgtItem["MapMinMaxScale"] = [np.min(tgtImages[i]), np.max(tgtImages[i])]

	# Create new map ID
	uniqueID = int(sem.GetUniqueNavID())
	while uniqueID in newNavIDs:
		uniqueID = int(sem.GetUniqueNavID())
	newNavIDs.append(uniqueID)

	tgtItem["MapID"] = [uniqueID]
	tgtItem["RealignedID"] = [drawnID]
	tgtItem["SamePosId"] = [newViewID]

	if i == 0:
		tgtItem["Note"] = [userName + "_tgts.txt"]
		tgtItem["Acquire"] = [1]
	else:
		tgtItem["Note"] = [tgt["tgtfile"]]
		tgtItem["Acquire"] = [0]

	navItems.append(tgtItem)
	navItems.append(viewItem)

# Write nav file

sem.SaveNavigator(navFile + "~")						# backup unaltered nav file
WriteNav(navHeader, navItems, navFile)
sem.ReadNavFile(navFile)

sem.Echo("Updated navigator file with target maps.")

# Write tgts file

WriteTargets(tgtsFilePath, targets)

sem.Echo("Target selection completed! " + str(len(targets)) + " targets were selected.")

