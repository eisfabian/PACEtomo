#!Python
# ===================================================================
#ScriptName	PACEtomo_selectTargets
# Purpose:	Selects targets, saves maps and writes text file with coords for PACEtomo.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2021/04/19
# Revision:	v1.6.1
# Last Change:	2023/05/08: fixed preTilt sign fix
# ===================================================================

############ SETTINGS ############ 

guidance	= True		# when set to False, pop-up windows with guidance are kept minimal and keyboard shortcuts are used instead

targetByShift 	= False		# ask to enter image shifts instead of dragging manually

targetPattern 	= False		# regular pattern of targets (holey support film)
alignToP 	= False		# refine vectors by aligning to hole reference in buffer P
size 		= 1 		# size of collection pattern (1: 3x3, 2: 5x5, 3: 7x7, ...)

drawBeam 	= True		# draws navigator item representing beam diameter
beamDiameter 	= 0 		# beam diameter [microns] (if 0, ReportIlluminatedArea will be used, which is only available on some Thermo Scientific microscopes)
maxTilt 	= 60		# tilt angle [degrees] to calculate stretching of beam perpendicular to tilt axis

# Advanced settings
sampleName	= ""		# optional prefix for all files created
useSearch 	= False		# use Search mode instead of View mode to find targets by dragging
vecA		= (0, 0)	# vectors for grid pattern [microns specimen shift] are determined automatically...
vecB		= (0, 0)	# ...only change if you want to setup pattern without alignToP reference
patternRot	= 0 		# rotation of pattern grid relative to tilt axis (used for filling a polygon with points)

########## END SETTINGS ########## 

import serialem as sem
import os
import copy
import glob
import numpy as np
import scipy as sp
import scipy.optimize
from scipy.signal import fftconvolve
import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import matplotlib.lines
import matplotlib.path
from matplotlib.backend_bases import MouseButton
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)

########### FUNCTIONS ###########

def parseTargets(targetFile):
	with open(targetFile) as f:
		content = f.readlines()
	targets = []
	savedRun = []
	branch = None
	resume = {"sec": 0, "pos": 0}
	settings = {}
	for line in content:
		col = line.strip(os.linesep).split(" ")
		if col[0] == "": continue
		if line.startswith("_set") and len(col) == 4:
			settings[col[1]] = col[3]
		elif line.startswith("_spos"):
			resume["sec"] = int(col[2].split(",")[0])
			resume["pos"] = int(col[2].split(",")[1])
		elif line.startswith("_tgt"):
			targets.append({})
			branch = None
		elif line.startswith("_pbr"):
			savedRun.append([{},{}])
			branch = 0
		elif line.startswith("_nbr"):
			branch = 1
		else:
			if branch is None:
				targets[-1][col[0]] = col[2]
			else:
				savedRun[-1][branch][col[0]] = col[2]
	for i in range(len(targets)):
		if "tgtfile" not in targets[i].keys(): targets[i]["tgtfile"] = False
		if "tsfile" not in targets[i].keys() or sem.DoesFileExist(targets[i]["tsfile"]) == 0: targets[i]["tsfile"] = False
		targets[i]["SSX"] = float(targets[i]["SSX"])
		targets[i]["SSY"] = float(targets[i]["SSY"])
		if "skip" not in targets[i].keys() or targets[i]["skip"] == "False": 
			targets[i]["skip"] = False 
		else: 
			targets[i]["skip"] = True

	if savedRun == []: savedRun = False
	sem.Echo("NOTE: Found " + str(len(targets)) + " targets in " + os.path.basename(targetFile) + ".")
	return targets, savedRun, resume, settings

def writeTargets(targetFile, targets, savedRun=False, resume={"sec": 0, "pos": 0}, settings={}):
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

def saveNewTarget(targetFile, targetNo, target):
	if targetPattern and usePolygon != 1:
		if targetNo == 1:
			mapIndex = int(sem.AddStagePosAsNavPoint(target["stageX"], target["stageY"], stageZ))
			sem.ChangeItemNote(mapIndex, userName + "_tgts.txt")
			sem.SetItemAcquire(mapIndex)
	else:
		# save map
		sem.OpenNewFile(userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc")
		sem.S("A")
		if targetNo == 1:
			mapIndex = int(sem.NewMap(0, userName + "_tgts.txt"))
			sem.SetItemAcquire(mapIndex)
		else:
			mapIndex = int(sem.NewMap(0, userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc"))
		sem.CloseFile()
	sem.ChangeItemLabel(mapIndex, str(targetNo).zfill(3))

	# write target file
	output = "_tgt = " + str(targetNo).zfill(3) + "\n"
	if not targetPattern:
		output += "tgtfile = " + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc" + "\n"
	output += "tsfile = " + userName + "_ts_" + str(targetNo).zfill(3) + ".mrc" + "\n"
	for key in target.keys():
		output += key + " = " + str(target[key]) + "\n"
	output += "skip = False" + 2 * "\n"
	with open(targetFile, "a") as f:
		f.write(output)
	sem.Echo("Target " + str(targetNo).zfill(3) + " (" + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc) with image shifts " + str(round(target["SSX"], 3)) + ", " + str(round(target["SSY"], 3)) + " was added.")

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

def drag():
	global userInput, targetNo
	if guidance:
		sem.OKBox("Please center your target by dragging the image using the right mouse button! Press the <b> key when finished!")
		sem.Echo("NOTE: Please center your target by dragging the image using the right mouse button! Press the <b> key when finished!")
	else:
		sem.Echo("NOTE: Please press the <b> key after centering your target! (To take another view image, immediately press the <v> key after <b>!)")
	while not sem.KeyBreak():
		sem.Delay(0.1)
	if guidance:
		userConfirm = sem.YesNoBox("\n".join(["PREVIEW?", "", "Do you want to take a preview image here?"]))
	else:
		userConfirm = 1 										# default to preview unless <v> key is pressed within one second after <b>
		for i in range(10):
			if sem.KeyBreak("v"):
				userConfirm = 0
				sem.Echo("Take new view image!")
				break
			sem.Delay(0.1)
	if userConfirm == 1:
		sem.L()
		if guidance: 
			userRefine = sem.YesNoBox("\n".join(["REFINE?", "", "Do you want to refine the position at this mag?"]))
		else:
			userRefine = 0 										# default to no refine unless <r> key is pressed within one second
			sem.Echo("NOTE: Please immediately press the <r> key to refine the position at this mag!")
			for i in range(10):
				if sem.KeyBreak("r"):
					userRefine = 1
					sem.Echo("Refine target!")
					break
				sem.Delay(0.1)			
		if userRefine == 1:
			if guidance:
				sem.OKBox("Please center your target by dragging the image using the right mouse button! Press the <b> key when finished!")
				sem.Echo("NOTE: Please center your target by dragging the image using the right mouse button! Press the <b> key when finished!")
			else:
				sem.Echo("NOTE: Please press the <b> key after centering your target!")
			while not sem.KeyBreak():
				sem.Delay(0.1)
			sem.L()
		userInput = sem.YesNoBox("\n".join(["SAVE?", "", "Do you want to use the current image and coordinates as target " + str(targetNo + 1) + "?","If you choose no, a new view image is taken to align your target!"]))

def loopAddTargets():
	global pointRefine, userInput, targetNo
	pointNo = 1 												# counter for points in group
	addTargets = 1
	while addTargets == 1:
		userInput = 0
		userSkip = 1											# initialize as 1 to not apply shifts in case of redo
		while userInput == 0:
			if targetByShift or (pointRefine == 1 and userSkip == 1):
				if targetByShift:
					sem.GoToLowDoseArea("R")
					shiftx = sem.EnterDefaultedNumber(0, 1, "Enter X shift:")
					shifty = sem.EnterDefaultedNumber(0, 1, "Enter Y shift:")
				else:
					sem.SetImageShift(0, 0)							# use 0,0 instead of ISX0,ISY0 to account for user shift of first target away from point
					shiftx, shifty = coordsRefine[pointNo]
				sem.ImageShiftByMicrons(shiftx, shifty)
			if useSearch: 
				sem.Search()
			else:
				sem.V()
			drag()
			if userInput == 0:
				if pointRefine == 1:
					userSkip = sem.YesNoBox("\n".join(["SKIP?", "", "Do you want to skip this point of the group?"]))
					if userSkip == 1:
						pointNo += 1
						if pointNo >= len(coordsRefine):				# disable pointRefine when last point of a group is skipped
							pointRefine = 0
							addTargets = sem.YesNoBox("All points of the selected group have been viewed. Do you want to add another target manually?")
				else:
					if guidance:
						addTargets = sem.YesNoBox("Do you want to keep looking for a target?")
					else:
						sem.Echo("NOTE: Please immediately press the <f> key to finish target search!")
						for i in range(10):
							if sem.KeyBreak("f"):
								addTargets = 0
								break
							sem.Delay(0.1)	
				if addTargets == 0: break
		if addTargets == 0: break

		target = {}
		target["SSX"], target["SSY"] = sem.ReportSpecimenShift()
		target["SSX"] -= SSX0
		target["SSY"] -= SSY0

		if np.linalg.norm(np.array([target["SSX"], target["SSY"]])) > imageShiftLimit - 0.5:
			userAbort = sem.YesNoBox("WARNING: This target is close to the image shift limit of SerialEM (" + str(imageShiftLimit) + " microns). It will probably be dropped during acquisition. Do you still want to save this target?")
			if userAbort == 0: continue

		targetNo += 1
		pointNo += 1

		ISX, ISY, *_ = sem.ReportImageShift()
		target["stageX"], target["stageY"] = sem.AdjustStagePosForNav(stageX, stageY, ISX, ISY)		
		saveNewTarget(tgtsFilePath, targetNo, target)

		if drawBeam:
			beamPolygons.append(drawBeamPolygon(target["stageX"], target["stageY"], stageZ, beamR, maxTilt))

		if pointRefine == 1 and pointNo >= len(coordsRefine):
			addTargets = sem.YesNoBox("All points of the selected group have been viewed. Do you want to add another target manually?")
			pointRefine = 0
		else:
			addTargets = sem.YesNoBox("\n".join(["ADD?", "", "Do you want to add another target?"]))

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

def drawBeamPolygon(stageX, stageY, stageZ, radius, angle):
	a = radius
	b = radius / np.cos(np.radians(angle))
	n = 32
	phi = angles_in_ellipse(n, a, b)

	ptsX = (a * np.cos(phi) + float(stageX)).round(3)
	ptsX = np.append(ptsX, ptsX[0])
	ptsY = (b * np.sin(phi) + float(stageY)).round(3)
	ptsY = np.append(ptsY, ptsY[0])

	sem.SetVariable('ptsX', listToSEMarray(ptsX))
	sem.SetVariable('ptsY', listToSEMarray(ptsY))

	polyID = int(sem.AddStagePointsAsPolygon("ptsX", "ptsY", stageZ))
	sem.ChangeItemColor(polyID, 3)

	return polyID

##########
# Source: https://github.com/Sabrewarrior/normxcorr2-python/blob/master/normxcorr2.py
def autoXcorr(image):
	image = image - np.mean(image)
	corr = fftconvolve(image, image[::-1, ::-1], mode="same")
	corr = corr / np.sum(np.square(image))
	return corr
##########

def vecByXcorr(diameter):
	# devSettings
	updatePlot 	= False		# show plot after finding each peak
	maxBinning 	= 10		# maximum binning of input image (depends on pixel size and diameter)
	pixPerHole 	= 20		# desired amount of pixels per hole diameter
	searchRangeFac	= 2 		# how many times the radius should be excluded around a found peak
	cropRangeFac 	= 5		# how many times should the searched image contain the search range
	maxPeaks 	= 6		# how many peaks should be searched
	accuracy 	= 0.1 		# dot product and normalized vector length difference needs to be less for vectors to be accepted

	buffer, *_ = sem.ReportCurrentBuffer()

	imgProp = sem.ImageProperties(buffer)
	pixSize = imgProp[4]
	binning = min(max(1, np.floor(1000 * diameter / pixPerHole / pixSize)), maxBinning)

	if binning > 1:
		sem.ReduceImage(buffer, binning)
		buffer = "A"
		imgProp = sem.ImageProperties(buffer)
		pixSize = imgProp[4]
	x = imgProp[0]
	y = imgProp[1]

	searchRadius = int(diameter * 1000 / pixSize / 2 * searchRangeFac)
	image = np.asarray(sem.bufferImage(buffer))

	imgStartX = int(max(0, x / 2 - cropRangeFac * searchRadius))
	imgStartY = int(max(0, y / 2 - cropRangeFac * searchRadius))
	image = image[imgStartX:(imgStartX + 2 * cropRangeFac * searchRadius), imgStartY:(imgStartY + 2 * cropRangeFac * searchRadius)] 

	output = autoXcorr(image)
	outsq = np.square(output)

	peaks = []
	while len(peaks) < maxPeaks:
		peaks.append(np.unravel_index(np.argmax(outsq), outsq.shape))
		rangeStartX = max(0, peaks[-1][0] - searchRadius)
		rangeStartY = max(0, peaks[-1][1] - searchRadius)
		outsq[rangeStartX:(rangeStartX + 2 * searchRadius), rangeStartY:(rangeStartY + 2 * searchRadius)] = np.zeros(outsq[rangeStartX:(rangeStartX + 2 * searchRadius), rangeStartY:(rangeStartY + 2 * searchRadius)].shape)
		if updatePlot:
			imgplot = plt.imshow(outsq)
			plt.show()

	peaks = np.array(peaks) + np.array([imgStartX, imgStartY])
	success = False

	ptID0 = int(sem.AddImagePosAsNavPoint(buffer, peaks[0][1], peaks[0][0], 0))

	vecA = np.zeros(2)
	i = 0
	while np.linalg.norm(vecA) < 2 * searchRadius:
		i += 1
		if i >= len(peaks):
			 sem.Echo("WARNING: Failed to find appropiately sized vector!")
			 sem.DeleteNavigatorItem(int(ptID0))
			 return success, np.zeros(2), np.zeros(2)
		vecA = peaks[i] - peaks[0]

	shortest = 1
	for i in range(1, len(peaks)):
		vecB = peaks[i] - peaks[0]
		if np.linalg.norm(vecB) < np.linalg.norm(vecA) and np.linalg.norm(vecB) > 2 * searchRadius:
			vecA = vecB
			shortest = i

	ptID1 = int(sem.AddImagePosAsNavPoint(buffer, peaks[shortest][1], peaks[shortest][0], 0))

	for i in range(1, len(peaks)):
		vecB = peaks[i] - peaks[0]
		dotprod = vecA.dot(vecB) / vecA.dot(vecA)
		if abs(dotprod) < accuracy and abs(np.linalg.norm(vecA) - np.linalg.norm(vecB)) / np.linalg.norm(vecA) < accuracy:
			ptID2 = int(sem.AddImagePosAsNavPoint(buffer, peaks[i][1], peaks[i][0], 0))
			success = True
			break

	if not success:
		sem.Echo("WARNING: Failed to find orthogonal vectors for pattern!")
	else:
		ptID0, x0, y0, *_ = sem.ReportOtherItem(ptID0)
		ptID1, x1, y1, *_ = sem.ReportOtherItem(ptID1)
		ptID2, x2, y2, *_ = sem.ReportOtherItem(ptID2)

		vecA = np.around((-x1 + x0, -y1 + y0), 3)
		vecB = np.around((-x2 + x0, -y2 + y0), 3)
		if np.linalg.norm(vecA) > diameter:
			if vecA[1] < 0: vecA = -vecA								# make vectors face towards y > 0
			if vecB[1] < 0: vecB = -vecB
			sem.Echo("NOTE: Orthogonal vectors for pattern were found!")
		else:
			sem.Echo("WARNING: Failed to find appropiately sized vectors for pattern!")
			success = False
	sem.DeleteNavigatorItem(int(ptID2))
	sem.DeleteNavigatorItem(int(ptID1))
	sem.DeleteNavigatorItem(int(ptID0))

	if vecB[1] > vecA[1]:											# ensure vecA has larger off tilt axis shift
		return success, vecB, vecA
	else:
		return success, vecA, vecB

def gui(targetFile):
	##########
	# Source: https://stackoverflow.com/a/25023944
	class DragDropListbox(tk.Listbox):									# listbox with drag'n'drop reordering of entries
		def __init__(self, master, **kw):
			kw['selectmode'] = tk.SINGLE
			tk.Listbox.__init__(self, master, kw)
			self.bind('<Button-1>', self.setCurrent)
			self.bind('<B1-Motion>', self.shiftSelection)
			self.curIndex = None

		def setCurrent(self, event):
			self.curIndex = self.nearest(event.y)
			if self.curIndex == 0: self.curIndex = None

		def shiftSelection(self, event):
			i = self.nearest(event.y)
			if i < self.curIndex and i != 0:
				x = self.get(i)
				self.delete(i)
				self.insert(i+1, x)
				self.curIndex = i
			elif i > self.curIndex and i != 0:
				x = self.get(i)
				self.delete(i)
				self.insert(i-1, x)
				self.curIndex = i	
	##########
	# Source: https://www.daniweb.com/programming/software-development/code/484591/a-tooltip-class-for-tkinter
	class CreateToolTip(object):
		def __init__(self, widget, text='widget info'):
			self.widget = widget
			self.text = text
			self.widget.bind("<Enter>", self.enter)
			self.widget.bind("<Leave>", self.close)
			self.widget.bind("<ButtonPress>", self.close)

		def enter(self, event=None):
			x = self.widget.winfo_rootx() + 25
			y = self.widget.winfo_rooty() + 30
			self.tw = tk.Toplevel(self.widget)							# creates a toplevel window
			self.tw.wm_overrideredirect(True)							# leaves only the label and removes the app window
			self.tw.wm_geometry("+%d+%d" % (x, y))
			labelTT = tk.Label(self.tw, text=self.text, justify='left', background='#ffffff', relief='solid', borderwidth=1, font=("Courier","10"))
			labelTT.pack(ipadx=1)

		def close(self, event=None):
			if self.tw:
				self.tw.destroy()
	##########

	def updateList():
		listbox.delete(0, tk.END)
		for i in range(len(targets)):
			listbox.insert(i, " ".join([str(i+1).zfill(3).ljust(5), str(targets[i]["tgtfile"]).ljust(fileNameLen + 2), str(targets[i]["tsfile"]).ljust(fileNameLen + 1), str(round(float(targets[i]["SSX"]), 2)).rjust(6), str(round(float(targets[i]["SSY"]), 2)).rjust(6), str(targets[i]["skip"])]))

	def plotTargets():
		colors = ["#5689bf" if not val["skip"] else '#aaaaaa' for val in targets]			# setup color array
		colors[0] = "#c92b27"										# color tracking TS
		legend_elements = [matplotlib.lines.Line2D([0], [0], color="#c92b27", label="tracking", lw="2"), matplotlib.lines.Line2D([0], [0], color="#5689bf", label="acquire", lw="2"), matplotlib.lines.Line2D([0], [0], color="#cccccc", label="skip", lw="2"), matplotlib.lines.Line2D([0], [0], color="#fab182", label="geo point", lw="2")]
		plt.legend(handles=legend_elements)
		if swapXY:											# set up axis to fit SerialEM view
			xval = np.array([float(val["SSY"]) for val in targets])
			yval = np.array([float(val["SSX"]) for val in targets])
			plt.axvline(x=0, color="#cccccc", ls="--")
			geoXval = np.array([val[1] for val in geoPoints])
			geoYval = np.array([val[0] for val in geoPoints])
		else:
			xval = np.array([float(val["SSX"]) for val in targets])
			yval = np.array([float(val["SSY"]) for val in targets])
			plt.axhline(y=0, color="#cccccc", ls="--")
			geoXval = np.array([val[0] for val in geoPoints])
			geoYval = np.array([val[1] for val in geoPoints])
		if invX:
			plt.gca().invert_xaxis()
		if invY:
			plt.gca().invert_yaxis()	
		plt.scatter(xval, yval, marker="s", facecolor="none", color=colors, s=1000, linewidths=2, picker=True, label="targets")
		if geoPoints != []:
			plt.scatter(geoXval, geoYval, marker="o", color="#fab182", s=100, picker=False)
		plt.margins(0.25, 0.25)
		plt.axis("equal")
		if showBeam:											# plot circles with plot size in microns (beamdiameter)
			xsize = abs(plt.gca().get_window_extent().width / (plt.gca().get_xlim()[1] - plt.gca().get_xlim()[0]))
			ysize = abs(plt.gca().get_window_extent().height / (plt.gca().get_ylim()[1] - plt.gca().get_ylim()[0]))
			size = (2 * beamR * min(xsize, ysize) * 72 / fig.dpi) **2
			plt.scatter(xval, yval, marker="o", facecolor="none", color="#ffd700", s=size, linewidths=2, picker=False)
			if geoPoints != []:
				plt.scatter(geoXval, geoYval, marker="o", facecolor="none", color="#fab182", s=size, linewidths=2, picker=False)
		for i in range(len(targets)):									# add target numbers to plot
			plt.annotate(str(i + 1).zfill(3), (xval[i], yval[i]), ha="center", va="center")		

	def realignTrack(rough=False, wiggle=0.5):
		global stageX, stageY, stageZ, SSX0, SSY0, ISX0, ISY0
		stageX, stageY, stageZ = sem.ReportStageXYZ()
		if "ISX0" in globals():
			sem.SetImageShift(ISX0, ISY0)
		if abs(stageX - float(targets[0]["stageX"])) > wiggle or abs(stageY - float(targets[0]["stageY"])) > wiggle or "ISX0" not in globals():	# test if stage was moved or tracking target was changed (with 0.5 micron wiggle room)
			sem.Echo("Realigning to tracking target...")
			if rough:
				if targets[0]["tgtfile"]:
					mapIndex = int(sem.NavIndexWithNote(targets[0]["tgtfile"].rsplit(".mrc", 1)[0] + "_view.mrc"))	# find map index of tracking tilt series view image from file name
				else:
					mapIndex = int(sem.NavIndexWithNote(userName + "_tgt_001_view.mrc"))	# use tgt_001_view in case tgtfile does not exist
			else:
				mapIndex = int(sem.NavIndexWithNote(userName + "_tgts.txt"))
			sem.RealignToOtherItem(mapIndex, 1)							# realign to tracking TS
			sem.GoToLowDoseArea("R")
			stageX, stageY, stageZ = sem.ReportStageXYZ()
			SSX0, SSY0 = sem.ReportSpecimenShift()
			ISX0, ISY0, *_ = sem.ReportImageShift()

	def onSelect(event):											# click on target in plot
		ind = event.ind[0]
		button = event.mouseevent.button
		dblclick = event.mouseevent.dblclick
		label = event.artist.get_label()
		if label == "targets":
			if button is MouseButton.RIGHT:								# right click to skip target
				skip(ind)

			if button is MouseButton.LEFT:								# left click to select target in list
				listbox.selection_clear(0, tk.END)
				listbox.selection_set(ind)
				listbox.see(ind)
				if dblclick:
					success = openTs(box=False)
					if not success:
						success = openTgt(box=False)
					if not success:
						tk.messagebox.showwarning(title="File not found", message="No tilt series or target file was found for this target!")

	def onClick(event):											# click in plot to add point
		button = event.button
		if button is MouseButton.MIDDLE:								# middle click to add geo point
			if swapXY:
				geoPoints.append([round(event.ydata, 3), round(event.xdata, 3)])
			else:
				geoPoints.append([round(event.xdata, 3), round(event.ydata, 3)])
			plt.clf()
			plotTargets()
			fig.canvas.draw()

	def onScroll(event):											# scroll in plot to zoom
		button = event.button
		xlim = plt.gca().get_xlim()
		ylim = plt.gca().get_ylim()
		if button == "up":
			scale_factor = 1 / 1.2
		else:
			scale_factor = 1.2
		plt.gca().set_xlim([xlim[0] * scale_factor, xlim[1] * scale_factor])
		plt.gca().set_ylim([ylim[0] * scale_factor, ylim[1] * scale_factor])
		fig.canvas.draw()

	def openTgt(box=True):											# open preview map with default application
		ind = listbox.curselection()[0]
		if targets[ind]["tgtfile"]:
			os.system("start " + os.path.join(curDir, targets[ind]["tgtfile"]))
			return True
		else:
			if box:
				tk.messagebox.showwarning(title="File not found", message="Target file was not found!")
			return False

	def openTs(box=True):											# open tilt series if it exists
		ind = listbox.curselection()[0]
		if targets[ind]["tsfile"]:
			os.system("start " + os.path.join(curDir, targets[ind]["tsfile"]))
			return True
		else:
			if box:
				tk.messagebox.showwarning(title="File not found", message="Tilt series file was not found!")
			return False

	def makeTrack():											# change the tracking target
		nonlocal targets, targetsOrig, geoPoints
		ind = listbox.curselection()[0]
		if not targets[ind]["tgtfile"]:
			tk.messagebox.showwarning(title="File not found", message="The chosen target does not have a tgt file for alignment!")
		else:
			confirm = tk.messagebox.askyesno(title="Confirmation", message="\n".join(["This will save all changes and cause additional exposures of the new tracking target.", "", " Do you want to proceed?"]))
			if not confirm:
				return
			realignTrack(rough=True)								# realign to view map of old tracking area
			mapIndex = int(sem.NavIndexWithNote(userName + "_tgts.txt"))				# find index of old tracking map and change acquire state and navNote
			sem.SetItemAcquire(mapIndex, 0)
			sem.ChangeItemNote(mapIndex, targets[0]["tgtfile"])

			SSXoffset = targets[ind]["SSX"]
			SSYoffset = targets[ind]["SSY"]
			for i in range(len(targets)):								# make new tracking target center of shifts
				targets[i]["SSX"] -= SSXoffset
				targets[i]["SSY"] -= SSYoffset
			targets.insert(0, targets.pop(ind))							# move target to start of list
			targets[0]["skip"] = False

			mapIndex = int(sem.NavIndexWithNote(targets[0]["tgtfile"]))				# find index of new tracking map and change acquire state and navNote
			sem.SetItemAcquire(mapIndex)
			sem.ChangeItemNote(mapIndex, userName + "_tgts.txt")

			viewIndex = int(sem.NavIndexWithNote(targets[0]["tgtfile"].rsplit(".mrc", 1)[0] + "_view.mrc"))
			if viewIndex > 0:									# if view map exists already, align to new tracking target
				realignTrack()
			else:											# if not, make view map
				sem.ImageShiftByMicrons(SSXoffset, SSYoffset)					# apply SSoffsets and take temp view map of new tracking target
				sem.V()
				sem.OpenNewFile(targets[0]["tgtfile"].rsplit(".mrc", 1)[0] + "_view_temp.mrc")
				sem.S("A")
				tempIndex = int(sem.NewMap(0, targets[0]["tgtfile"].rsplit(".mrc", 1)[0] + "_view_temp.mrc"))
				sem.CloseFile()
				realignTrack()									# move stage and realign to new tracking target
				sem.DeleteNavigatorItem(tempIndex)
				sem.V()										# take view map of new tracking target
				sem.OpenNewFile(targets[0]["tgtfile"].rsplit(".mrc", 1)[0] + "_view.mrc")
				sem.S("A")
				sem.NewMap(0, targets[0]["tgtfile"].rsplit(".mrc", 1)[0] + "_view.mrc")
				sem.CloseFile()
				sem.GoToLowDoseArea("R")

			if geoPoints != []:									# shift geoPoints accordingly
				for i in range(len(geoPoints)):
					geoPoints[i][0] -= SSXoffset
					geoPoints[i][1] -= SSYoffset

			targetsOrig = copy.deepcopy(targets)							# make new backup copy for reset
			saveFile(ask=False)									# save changes to file
			updateList()
			plt.clf()
			plotTargets()
			fig.canvas.draw()

	def skip(ind=None):											# skip selected target
		if ind == None:
			ind = listbox.curselection()[0]	
		if ind != 0:
			if savedRun:
				if savedRun[ind][0]["skip"] == "True" and savedRun[ind][1]["skip"] == "True":
					sem.Echo("WARNING: Targets previously skipped during acquisition cannot be unskipped!")
					return
			targets[ind]["skip"] = not targets[ind]["skip"]
			updateList()
			listbox.selection_set(ind)
			listbox.see(ind)
			plt.clf()
			plotTargets()
			fig.canvas.draw()
		else:
			sem.Echo("WARNING: Tracking target cannot be skipped!")

	def loadMap():
		global mapLabel
		if "mapLabel" not in globals():
			mapLabel = ""
		mapLabel = tk.simpledialog.askstring("Map Label", "\n".join(["Please enter the navigator label", "of the map you want to load!"]), initialvalue=mapLabel)
		mapIndex = int(sem.NavIndexWithLabel(mapLabel))
		if mapIndex == 0:
			tk.messagebox.showwarning(title="Map not found", message="Map with label '" + mapLabel + "' was not found!")
		else:
			sem.LoadOtherMap(mapIndex)

	def saveOrder():											# save order of listbox (after manual reordering by dragging)
		nonlocal targets
		order = [int(entry.split(" ", 1)[0]) - 1 for entry in listbox.get(0, tk.END)]
		if order[0] != 0:
			tk.messagebox.showwarning(title="Not allowed", message="You can't change the order of the tracking target!")
			return
		targets = [targets[i] for i in order]
		updateList()
		plt.clf()
		plotTargets()
		fig.canvas.draw()

	def moreTargets():											# append new targets (closes GUI)
		global userInput, pointRefine, reopen, SSX0, SSY0
		confirm = tk.messagebox.askyesno(title="Confirmation", message="\n".join(["This will save all changes and start the procedure to add extra targets by dragging.", "", " Do you want to proceed?"]))
		if not confirm:
			return
		saveFile(ask=False)										# save changes to file
		closeGUI()
		realignTrack()
		if drawBeam:											# make beam polygons visible
			for polyID in beamPolygons:
				sem.ChangeItemDraw(polyID, 1)
		sem.Echo("Adding targets...")		
		userInput = 0
		pointRefine = 0
		loopAddTargets()
		sem.SetImageShift(ISX0, ISY0)
		reopen = True											# reopen GUI after finishing
		return

	def saveViews():
		btnSaveViews.grid_forget()									# replace view button with progress bar
		progressBar.grid(column=1, row=4, sticky=tk.W, pady=6)
		disabledWidgets = [btnSave, btnReset, btnMakeTrack, btnLoadMap, btnAddTargets, btnReorder, btnMeasureGeometry]
		keepDisabled = []
		for widget in disabledWidgets:									# disable buttons that could interfere with measure geometry
			if widget["state"] == tk.DISABLED:							# keep already disabled buttons disabled
				keepDisabled.append(widget)
			else:
				widget["state"] = tk.DISABLED

		realignTrack(rough=True)
		sem.GoToLowDoseArea("R")
		for i in range(len(targets)):
			if targets[i]["tgtfile"]:
				viewName = targets[i]["tgtfile"].rsplit(".mrc", 1)[0] + "_view.mrc"		# if tgtfile is present, retain numbering
			else:
				viewName = userName + "_tgt_" + str(i + 1).zfill(3) + "_view.mrc"		# if not, renumber
			viewIndex = int(sem.NavIndexWithNote(viewName))
			if viewIndex == 0:									# take view image if it does not exist
				sem.ImageShiftByMicrons(targets[i]["SSX"], targets[i]["SSY"])
				sem.V()
				sem.OpenNewFile(viewName)
				sem.S("A")
				sem.NewMap(0, viewName)
				sem.CloseFile()
			else:
				sem.LoadOtherMap(viewIndex)
			targets[i]["viewfile"] = viewName
			sem.SnapshotToFile(0, 0, "0", "JPG", "JPG", viewName.rsplit(".mrc", 1)[0] + ".jpg")
			sem.GoToLowDoseArea("R")
			sem.SetImageShift(ISX0, ISY0)
			progressBar["value"] = (i + 1) / len(targets)						# update progress bar
			top.update()
		progressBar.grid_forget()									# restore button
		btnSaveViews.grid(column=1, row=4, sticky=tk.W, pady=3)
		for widget in [item for item in disabledWidgets if item not in keepDisabled]:			# reenable interfering buttons
			widget["state"] = tk.NORMAL

	def copyAcq():												# find nav points set to acquire and replace Note with tgts file
		copied = 0
		acqItems = []
		for i in range(int(sem.ReportNumTableItems())):
			if sem.ReportItemAcquire(i + 1) == 1:
				otherItem = sem.ReportOtherItem(i + 1)
				acqNote = sem.GetVariable("navNote")							
				if not acqNote.endswith(".txt"):						# check if selected nav item already has tgts file
					acqItems.append([i + 1, otherItem[1], otherItem[2], otherItem[3]])	# add to list of items to be changed (include stage xyz)
				elif acqNote.startswith(fileStem):						# check if tgts file is copy of current tgts file
					acqNo = acqNote.rsplit(".txt", 1)[0].rsplit("tgts_p", 1)		# figure out copy number
					if len(acqNo) > 1:
						acqNo = int(acqNo[-1])
						if acqNo > copied:
							copied = acqNo
		for item in acqItems:
			copied += 1
			targetFileCopy = os.path.join(curDir, fileStem + "_p" + str(copied).zfill(2) + ".txt")	# add stage position number to tgts file name
			sem.ChangeItemNote(item[0], os.path.basename(targetFileCopy))				# adjust nav note accordingly
			sem.ChangeItemLabel(item[0], str(1).zfill(3))
			targetsTemp = copy.deepcopy(targets)							# make temp deepcopy of targets
			targetsTemp[0]["stageX"] = item[1]
			targetsTemp[0]["stageY"] = item[2]
			groupIndex = int(sem.GetUniqueNavID())
			for i in range(len(targetsTemp)):							# add stage position number to ts file names
				if "tsfile" in targetsTemp[i].keys():
					targetsTemp[i]["tsfile"] = userName + "_p" + str(copied).zfill(2) + "_ts_" + str(i + 1).zfill(3) + ".mrc"
				if i > 0:									# if not tracking tgt, add tgt nav point
					stageShift = ss2sMatrix @ np.array([targetsTemp[i]["SSX"], targetsTemp[i]["SSY"]])
					ptIndex = int(sem.AddStagePosAsNavPoint(targetsTemp[0]["stageX"] + stageShift[0], targetsTemp[0]["stageY"] + stageShift[1], item[3], groupIndex))
					#ptIndex = int(sem.AddStagePosAsNavPoint(targetsTemp[0]["stageX"] + ss2sMatrix[0] * targetsTemp[i]["SSX"], targetsTemp[0]["stageY"] + ss2sMatrix[3] * targetsTemp[i]["SSY"], item[3], groupIndex))
					sem.ChangeItemLabel(ptIndex, str(i + 1).zfill(3))
			writeTargets(targetFileCopy, targetsTemp, settings=settings)				# write new tgts file for each stage position
		tk.messagebox.showinfo(title="Target file copied", message="The targets file was copied to " + str(copied) + " navigator points!")
		sem.Echo("NOTE: Copied tgts file to " + str(copied) + " navigator points!")

	def saveFile(ask=True):											# save targets to file
		if ask:
			confirm = tk.messagebox.askyesno(title="Confirmation", message="\n".join(["Save changes?", "", "Do you want to save changes and overwrite your targets file?"]))
			if not confirm:
				return
		targetsTemp = copy.deepcopy(targets)
		for i in range(len(targetsTemp)):
			if not targetsTemp[i]["tsfile"]:							# restore tsfile value in case file was not present
				if targetsTemp[i]["tgtfile"]:
					col = targetsTemp[i]["tgtfile"].split("_tgt_")				# if tgtfile is present, retain numbering
					targetsTemp[i]["tsfile"] = col[0] + "_ts_" + col[1]	
				else:
					targetsTemp[i]["tsfile"] = userName + "_ts_" + str(i + 1).zfill(3) + ".mrc"	# if not, renumber (there should not be a mixed case)
					targetsTemp[i].pop("tgtfile")						# don't save tgtfile to tgts file if False
		os.replace(targetFile, targetFile + "~")							# make backup
		writeTargets(targetFile, targetsTemp, savedRun, resume, settings)
		sem.Echo("NOTE: Changes to targets file were saved!")

	def resetOrder():											# restore targets array as read from file
		nonlocal targets, geoPoints
		targets = copy.deepcopy(targetsOrig)
		geoPoints = []											# also resets geoPoints since possible shifts cannot be recovered
		updateList()
		plt.clf()
		plotTargets()
		fig.canvas.draw()

	def readEntry(*args):											# read all entry fields upon change
		settings["startTilt"] = ecStartTilt.get()
		settings["minTilt"] = ecMinTilt.get()
		settings["maxTilt"] = ecMaxTilt.get()
		settings["step"] = ecStep.get()
		settings["pretilt"] = ecPretilt.get()
		settings["rotation"] = ecRotation.get()
		for key in settings.keys():
			try:
				settings[key] = float(settings[key])
			except ValueError:
				if settings[key] != "":
					settings[key] = ""
					sem.Echo("WARNING: " + key + " is not a number!")

	def resetGeo():
		nonlocal geoPoints
		geoPoints = []
		plt.clf()
		plotTargets()
		fig.canvas.draw()

	def measureGeo():											# measure geometry
		if len(geoPoints) < 3:
			tk.messagebox.showwarning(title="Not enough points", message="\n".join(["You need at least 3 points to measure the geometry!","","Add points by middle clicking in the plot. Points should not be on targets, dark areas or holes."]))
		else:
			btnMeasureGeometry.grid_forget()							# replace button with progress bar
			progressBar.grid(column=3, row=5, sticky=tk.E, pady=6, padx=5)
			disabledWidgets = [btnSave, btnReset, btnResetGeoPts, btnMakeTrack, btnSaveViews, btnLoadMap, btnAddTargets]
			keepDisabled = []
			for widget in disabledWidgets:								# disable buttons that could interfere with measure geometry
				if widget["state"] == tk.DISABLED:						# keep already disabled buttons disabled
					keepDisabled.append(widget)
				else:
					widget["state"] = tk.DISABLED

			sem.SetImageShift(0,0)
			realignTrack(rough=True)
			sem.Echo("Measuring geometry...")
			geoXYZ = [[], [], []]
			for i in range(len(geoPoints)):
				sem.ImageShiftByMicrons(geoPoints[i][0], geoPoints[i][1])
				sem.G(-1)
				defocus, _ = sem.ReportAutoFocus()
				if defocus != 0:
					geoXYZ[0].append(geoPoints[i][0])
					geoXYZ[1].append(geoPoints[i][1])
					geoXYZ[2].append(defocus)
				sem.SetImageShift(0,0)
				progressBar["value"] = (i + 1) / len(geoPoints)				# update progress bar
				top.update()
			##########
			# Source: https://math.stackexchange.com/q/99317
			# subtract out the centroid and take the SVD, extract the left singular vectors, the corresponding left singular vector is the normal vector of the best-fitting plane
			svd = np.linalg.svd(geoXYZ - np.mean(geoXYZ, axis=1, keepdims=True))
			left = svd[0]
			norm = left[:, -1]
			##########		
			sem.Echo("Fitted plane into cloud of " + str(len(geoPoints)) + " points.")
			sem.Echo("Normal vector: " + str(norm))
			sign = 1 if norm[1] <= 0 else -1
			tilty = sign * np.degrees(np.arccos(norm[2]))
			sem.Echo("Estimated pretilt: " + str(round(tilty, 1)) + " degrees")
			rotation = -np.degrees(np.arctan(norm[0]/norm[1]))
			sem.Echo("Estimated rotation: " + str(round(rotation, 1)) + " degrees")
			tk.messagebox.showinfo(title="measureGeometry", message="\n".join(["Geometry measurement:","","Estimated pretilt: " + str(round(tilty, 1)) + " degrees","Estimated rotation: " + str(round(rotation, 1)) + " degrees"]))
			entryPretilt.delete(0, tk.END)								# update entry fields
			entryPretilt.insert(0, str(round(tilty, 1)))
			entryRotation.delete(0, tk.END)
			entryRotation.insert(0, str(round(rotation, 1)))
			readEntry()
			progressBar.grid_forget()								# restore button
			btnMeasureGeometry.grid(column=3, row=5, sticky=tk.E, pady=3, padx=5)
			for widget in [item for item in disabledWidgets if item not in keepDisabled]:		# reenable interfering buttons
				widget["state"] = tk.NORMAL

	def toggleBeam():
		nonlocal showBeam
		showBeam = not showBeam
		plt.clf()
		plotTargets()
		fig.canvas.draw()
		if drawBeam:											# also toggle beam polygons in SerialEM
			for polyID in beamPolygons:
				sem.ChangeItemDraw(polyID)

	def swapAxes():
		nonlocal swapXY
		swapXY = not swapXY
		plt.clf()
		plotTargets()
		fig.canvas.draw()

	def invertX():
		nonlocal invX
		invX = not invX
		plt.clf()
		plotTargets()
		fig.canvas.draw()

	def invertY():
		nonlocal invY
		invY = not invY
		plt.clf()
		plotTargets()
		fig.canvas.draw()

	global reopen, targetNo
	reopen = False

	targets, savedRun, resume, settings = parseTargets(targetFile)
	targetNo = len(targets)	

	targetsOrig = copy.deepcopy(targets)									# make backup copy for reset
	geoPoints = []
	if targetPattern and "size" in settings.keys():
		vecA = (float(settings["vecA0"]), float(settings["vecA1"]))
		vecB = (float(settings["vecB0"]), float(settings["vecB1"]))
		size = int(float(settings["size"]))
		if size > 1:
			geoPoints.append([0.5 * (vecA[0] + vecB[0]), 0.5 * (vecA[1] + vecB[1])])
		geoPoints.append([(size - 0.5) * (vecA[0] + vecB[0]), (size - 0.5) * (vecA[1] + vecB[1])])
		geoPoints.append([(size - 0.5) * (vecA[0] - vecB[0]), (size - 0.5) * (vecA[1] - vecB[1])])
		geoPoints.append([(size - 0.5) * (-vecA[0] + vecB[0]), (size - 0.5) * (-vecA[1] + vecB[1])])
		geoPoints.append([(size - 0.5) * (-vecA[0] - vecB[0]), (size - 0.5) * (-vecA[1] - vecB[1])])

	showBeam = False
	if drawBeam:												# draw beam polygons but keep invisible until toggled
		if beamPolygons == []:
			for i in range(len(targets)):
				if "stageX" in targets[i].keys():
					beamPolygons.append(drawBeamPolygon(targets[i]["stageX"], targets[i]["stageY"], 0, beamR, maxTilt))
					sem.ChangeItemDraw(beamPolygons[-1], 0)
		else:
			for polyID in beamPolygons:
				sem.ChangeItemDraw(polyID, 0)

	s2cMatrix = sem.SpecimenToCameraMatrix(0)									# figure out axis of plot to match SerialEM view
	if(abs(s2cMatrix[0]) > abs(s2cMatrix[1])):
		swapXY = False
		invX = True if s2cMatrix[0] < 0 else False
		invY = True if s2cMatrix[3] < 0 else False		
	else:
		swapXY = True
		invX = True if s2cMatrix[1] < 0 else False
		invY = True if s2cMatrix[2] < 0 else False

	# create a root window.
	top = tk.Tk()
	top.option_add("*font", "Courier")
	top.title("Targets")
	top.columnconfigure(2, weight=2)
	top.rowconfigure(6, weight=2)
	pixel = tk.PhotoImage(width=1, height=1)

	# create target list
	fileNameLen = max(len(str(targets[0]["tgtfile"])), len(str(targets[0]["tsfile"])))
	labelList = tk.Label(top, text=" ".join(["Tgt".ljust(5), "Tgtfile".ljust(fileNameLen + 2), "TSfile".ljust(fileNameLen + 1), "SSX".rjust(6), "SSY".rjust(6), "Skip".ljust(70 - 24 - 2 * fileNameLen)])) 
	labelList.grid(column=0, row=0, sticky=tk.E)

	listbox = DragDropListbox(top, height=10, selectmode="SINGLE", width=70, activestyle='dotbox')
	updateList()
	listbox.grid(column=0, row=1, rowspan=6, padx=10, sticky=tk.NE)

	# create buttons
	btnWidth = 150
	btnHeight = 20
	
	btnOpenTgtfile = tk.Button(top, text="Open Tgtfile", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=openTgt)
	btnOpenTgtfile.grid(column=1, row=1, sticky=tk.W, pady=0)
	CreateToolTip(btnOpenTgtfile, "Opens .mrc file of selected target preview image.")

	btnOpenTSfile = tk.Button(top, text="Open TSfile", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=openTs)
	btnOpenTSfile.grid(column=2, row=1, sticky=tk.W, pady=0, padx=5)
	CreateToolTip(btnOpenTSfile, "Opens .mrc tilt series stack file of selected target.")

	btnSkip = tk.Button(top, text="Skip", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=skip)
	btnSkip.grid(column=1, row=2, sticky=tk.W, pady=3)
	CreateToolTip(btnSkip, "Toggles skipping of selected target during collection.")

	btnBorderMakeTrack = tk.Frame(top, highlightbackground="#ffd700", highlightthickness=2)		
	btnMakeTrack = tk.Button(btnBorderMakeTrack, text="Make Track", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=makeTrack)
	if savedRun: btnMakeTrack["state"] = tk.DISABLED
	btnMakeTrack.pack()
	btnBorderMakeTrack.grid(column=2, row=2, sticky=tk.W, pady=3, padx=5)
	CreateToolTip(btnMakeTrack, "\n".join(["Makes selected target the tracking target.", "This will cause additional exposures on the selected target for realignment."]))

	btnLoadMap = tk.Button(top, text="Load Map", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=loadMap)
	btnLoadMap.grid(column=1, row=3, sticky=tk.W, pady=0)
	CreateToolTip(btnLoadMap, "Loads map into SerialEM for cross referencing.")

	btnBorderAddTargets = tk.Frame(top, highlightbackground="#ffd700", highlightthickness=2)
	btnAddTargets = tk.Button(btnBorderAddTargets, text="Add Targets", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=moreTargets)
	if savedRun: btnAddTargets["state"] = tk.DISABLED
	btnAddTargets.pack()
	btnBorderAddTargets.grid(column=2, row=3, sticky=tk.W, pady=0, padx=5)
	CreateToolTip(btnAddTargets, "\n".join(["Continues the target selection process.", "This might cause additional exposures on the tracking target for realignment."]))

	btnSaveViews = tk.Button(top, text="Save Views", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=saveViews)
	btnSaveViews.grid(column=1, row=4, sticky=tk.W, pady=3)
	CreateToolTip(btnSaveViews, "Takes and saves view image for each target.")

	#xbutton = tk.Button(top, text="x", image=pixel, compound="center", height=btnHeight, width=btnWidth, command="")
	#if savedRun: xbutton["state"] = tk.DISABLED
	#xbutton.grid(column=2, row=4, sticky=tk.W, pady=3, padx=5)
	#CreateToolTip(xbutton, "Placeholder.")
	progressBar = tk.ttk.Progressbar(top, orient="horizontal", length=160, mode="determinate")
	progressBar["maximum"] = 1
	progressBar["value"] = 0

	btnReorder = tk.Button(top, text="Reorder", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=saveOrder)
	if savedRun: btnReorder["state"] = tk.DISABLED
	btnReorder.grid(column=1, row=5, sticky=tk.W, pady=0)
	CreateToolTip(btnReorder, "Applies current order as displayed in the list.")

	btnCopyToAcq = tk.Button(top, text="Copy to Acq", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=copyAcq)
	if savedRun or "size" not in settings.keys(): btnCopyToAcq["state"] = tk.DISABLED
	btnCopyToAcq.grid(column=2, row=5, sticky=tk.W, pady=0, padx=5)
	CreateToolTip(btnCopyToAcq, "Applies tgt pattern to all navigator points marked as Acquire.")

	btnSave = tk.Button(top, text="Save", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=saveFile)
	btnSave.grid(column=1, row=6, sticky=tk.SW, pady=3)
	CreateToolTip(btnSave, "Saves all changes to the tgts file.")

	btnReset = tk.Button(top, text="Reset", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=resetOrder)
	btnReset.grid(column=2, row=6, sticky=tk.SW, pady=3, padx=5)
	CreateToolTip(btnReset, "Resets changes made since window was opened.")
	  
	# create settings area
	labelSettings = tk.Label(top, text="Settings (optional)")
	labelSettings.grid(column=3, columnspan=5, row=0, )

	labelTilt = tk.Label(top, text = "Tilt angles [deg]")
	labelTilt.grid(column=3, row=2, sticky=tk.NE)
	labelStartTilt = tk.Label(top, text = "Start")
	labelStartTilt.grid(column=4, row=1, sticky=tk.SW)
	labelMinTilt = tk.Label(top, text = "Min")
	labelMinTilt.grid(column=5, row=1, sticky=tk.SW)
	labelMaxTilt = tk.Label(top, text = "Max")
	labelMaxTilt.grid(column=6, row=1, sticky=tk.SW)
	labelStepTilt = tk.Label(top, text = "Step")
	labelStepTilt.grid(column=7, row=1, sticky=tk.SW)

	ecStartTilt = tk.StringVar(top, settings["startTilt"] if "startTilt" in settings.keys() else "")
	ecStartTilt.trace_add("write", readEntry)
	entryStartTilt = tk.Entry(top, textvariable=ecStartTilt, width=6)
	if savedRun: entryStartTilt["state"] = tk.DISABLED
	entryStartTilt.grid(column=4, row=2, sticky=tk.NW)

	ecMinTilt = tk.StringVar(top, settings["minTilt"] if "minTilt" in settings.keys() else "")
	ecMinTilt.trace_add("write", readEntry)
	entryMinTilt = tk.Entry(top, textvariable=ecMinTilt, width=4)
	entryMinTilt.grid(column=5, row=2, sticky=tk.NW)

	ecMaxTilt = tk.StringVar(top, settings["maxTilt"] if "maxTilt" in settings.keys() else "")
	ecMaxTilt.trace_add("write", readEntry)
	entryMaxTilt = tk.Entry(top, textvariable=ecMaxTilt, width=4)
	entryMaxTilt.grid(column=6, row=2, sticky=tk.NW)

	ecStep = tk.StringVar(top, settings["step"] if "step" in settings.keys() else "")
	ecStep.trace_add("write", readEntry)
	entryStep = tk.Entry(top, textvariable=ecStep, width=5)
	if savedRun: entryStep["state"] = tk.DISABLED
	entryStep.grid(column=7, row=2, sticky=tk.NW)

	labelPretilt = tk.Label(top, text = "Pretilt [deg]")
	labelPretilt.grid(column=3, row=3, sticky=tk.NE)

	ecPretilt = tk.StringVar(top, settings["pretilt"] if "pretilt" in settings.keys() else "")
	ecPretilt.trace_add("write", readEntry)
	entryPretilt = tk.Entry(top, textvariable=ecPretilt, width=20)
	if savedRun: entryPretilt["state"] = tk.DISABLED
	entryPretilt.grid(column=4, columnspan=4, row=3, sticky=tk.NW)

	labelRotation = tk.Label(top, text = "Rotation [deg]")
	labelRotation.grid(column=3, row=4, sticky=tk.NE)

	ecRotation = tk.StringVar(top, settings["rotation"] if "rotation" in settings.keys() else "")
	ecRotation.trace_add("write", readEntry)
	entryRotation = tk.Entry(top, textvariable=ecRotation, width=20)
	if savedRun: entryRotation["state"] = tk.DISABLED
	entryRotation.grid(column=4, columnspan=4, row=4, sticky=tk.NW)

	# create geometry buttons
	btnMeasureGeometry = tk.Button(top, text="Measure Geometry", image=pixel, compound="center", height=btnHeight, command=measureGeo)
	if savedRun: btnMeasureGeometry["state"] = tk.DISABLED
	btnMeasureGeometry.grid(column=3, row=5, sticky=tk.E, pady=3, padx=5)
	CreateToolTip(btnMeasureGeometry, "\n".join(["Runs routine to measure pretilt and rotation.", "This will cause exposures on the chosen geometry points."]))

	btnResetGeoPts = tk.Button(top, text="Reset Geo Pts", image=pixel, compound="center", height=btnHeight, command=resetGeo)
	btnResetGeoPts.grid(column=4, columnspan=4, row=5, sticky=tk.W, pady=3, padx=5)
	CreateToolTip(btnResetGeoPts, "Deletes all geometry points.")

	# create target plot
	colors = ["#5689bf" if not val["skip"] else '#aaaaaa' for val in targets]
	colors[0] = "#c92b27"
	fig = plt.figure(tight_layout=True, figsize=(15, 8))

	plotTargets()

	canvas = FigureCanvasTkAgg(fig, master = top)  
	canvas.draw()
	canvas.get_tk_widget().grid(column=0, columnspan=8, row=7, padx=10, pady=10)				# placing the canvas on the Tkinter window
	fig.canvas.mpl_connect('pick_event', onSelect)
	fig.canvas.mpl_connect('button_press_event', onClick)
	fig.canvas.mpl_connect('scroll_event', onScroll)

	toolbar_frame = tk.Frame(top)										# creating the Matplotlib toolbar
	toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
	toolbar.update()
	toolbar_frame.grid(column=0, row=8, sticky=tk.W, padx=10)

	# create axis buttons
	btnToggleBeam = tk.Button(top, text="Toggle beam", image=pixel, compound="center", height=btnHeight, command=toggleBeam)
	btnToggleBeam.grid(column=1, row=8, sticky=tk.W, pady=5, padx=5)
	CreateToolTip(btnToggleBeam, "\n".join(["Toggles display of beam diameter.", "(Beam diameter is taken from IlluminatedArea or from value set in script.)"]))

	btnSwapXY = tk.Button(top, text="Swap XY", image=pixel, compound="center", height=btnHeight, command=swapAxes)
	btnSwapXY.grid(column=2, row=8, sticky=tk.E, pady=5, padx=5)
	CreateToolTip(btnSwapXY, "Swaps axes of target plot.")

	btnInvertX = tk.Button(top, text="Invert X", image=pixel, compound="center", height=btnHeight, command=invertX)
	btnInvertX.grid(column=3, row=8, sticky=tk.E, pady=5, padx=5)
	CreateToolTip(btnInvertX, "Inverts X axis of target plot.")

	btnInvertY = tk.Button(top, text="Invert Y", image=pixel, compound="center", height=btnHeight, command=invertY)
	btnInvertY.grid(column=4, columnspan=4, row=8, sticky=tk.W, pady=5, padx=5)
	CreateToolTip(btnInvertY, "Inverts Y axis of target plot.")

	def closeGUI():
		if targets != targetsOrig:
			saveFile()
		top.quit()
		top.destroy()

	top.protocol("WM_DELETE_WINDOW", closeGUI)
	top.lift()
	top.attributes("-topmost", True)
	top.after_idle(top.attributes, "-topmost", False)
	top.mainloop()

######## END FUNCTIONS ########

sem.SuppressReports()
sem.Pause("Please make sure that 'Move stage for big mouse shifts' is unchecked!")
if int(sem.ReportAxisPosition("F")[0]) != 0 and sem.IsVariableDefined("warningFocusArea") == 0:
	sem.Pause("WARNING: Position of Focus area is not 0! Please set it to 0 if you intend to use the Measure Geometry function!")
	sem.SetPersistentVar("warningFocusArea", "")
sem.UserSetDirectory("Please choose a directory for saving targets and tilt series!")

imageShiftLimit = sem.ReportProperty("ImageShiftLimit")

navSize = sem.ReportNumTableItems()
if navSize > 0:
	navInfo = sem.ReportNavItem()
	navID = int(navInfo[0])										# check if selected nav item already has tgts file
	navNote = sem.GetVariable("navNote")
else:
	navNote = ""
fileStem = navNote.rsplit(".txt", 1)[0]
curDir = sem.ReportDirectory()

if fileStem != "":
	userName = fileStem.split("_tgts")[0]
	tf = sorted(glob.glob(os.path.join(curDir, userName + "_tgts.txt")))					# find original tgts file
	tfp = sorted(glob.glob(os.path.join(curDir, userName + "_tgts_p??.txt")))				# find tgts files copied to other positions
	tfr = sorted(glob.glob(os.path.join(curDir, fileStem + "_run??.txt")))					# find run files of current nav item but not other copied tgts file
	tf.extend(tfr)												# only add run files to list of considered files
else:
	tf = []
editTgts = 0
if tf != []:
	editTgts = sem.YesNoBox("\n".join(["EDIT TARGETS?", "", "The selected navigator item already has a tgts file attached (" + os.path.basename(tf[-1]) + "). Do you want to edit these targets?"]))
	if editTgts == 1:
		tgtsFilePath = tf[-1]

sem.GoToLowDoseArea("R")											# need SS to stage matrix for conversion
ss2sMatrix = np.array(sem.SpecimenToStageMatrix(0)).reshape((2, 2))
s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
beamPolygons = []
beamR = beamDiameter / 2

if editTgts == 0:
	if navSize > 0:
		groupInfo = sem.ReportGroupStatus()
	else:
		groupInfo = [0, 0, 0]
	pointRefine = 0
	usePolygon = 0
	if not targetPattern and not targetByShift and groupInfo[1] > 0:
		pointRefine = sem.YesNoBox("\n".join(["GROUP OF POINTS?", "", "The selected navigator item is part of a group of " + str(groupInfo[2]) + "points. Do you want to use these points as initial target coordinates to be refined?"]))
		if pointRefine == 1:
			groupID = int(groupInfo[1])
			sem.GetNavGroupStageCoords(groupID, "groupStageX", "groupStageY", "groupStageZ")
			groupStage = np.column_stack((np.array(sem.GetVariable("groupStageX").split(), dtype=float), np.array(sem.GetVariable("groupStageY").split(), dtype=float)))
			groupStageZ = float(sem.GetVariable("groupStageZ").split()[0])
			coordsRefine = np.array([s2ssMatrix @ x for x in groupStage - groupStage[0]])
	elif targetPattern and not alignToP and navInfo[4] == 1:						# if targetPattern and nav item is polygon
		usePolygon = sem.YesNoBox("\n".join(["POLYGON?", "", "The selected navigator item is a polygon. Do you want to fill it with a grid of points based on the beam diameter?"]))
		if usePolygon == 1:
			sem.SaveNavigator()									# parse nav file to get polygon coords
			navFile = sem.ReportNavFile()
			navHeader, navItems = parseNav(navFile)

			vertices = np.vstack([navItems[navID - 1]["PtsX"], navItems[navID - 1]["PtsY"]]).transpose()
			polygon = matplotlib.path.Path(vertices)

	sem.EnterString("userName","Please provide a rootname for the PACE-tomo collection area!")
	userName = sampleName + sem.GetVariable("userName")

	tgtsFilePath = os.path.join(curDir, userName + "_tgts.txt")

	# Make sure tgts file is unique
	counter = 1
	while os.path.exists(tgtsFilePath):
		counter += 1
		tgtsFilePath = os.path.join(curDir, userName + str(counter) + "_tgts.txt") 
	if counter > 1:
		userName = userName + str(counter)

	# align center
	sem.TiltTo(0)
	sem.ResetImageShift()
	if pointRefine == 1:
		sem.MoveStageTo(*groupStage[0], groupStageZ)
	elif alignToP:												# center hole for center of tgtPattern
		sem.V()
		sem.AlignTo("P")
		if float(sem.ReportDefocus()) < -50:
			sem.Pause("WARNING: Large defocus offsets for View can cause image shift offsets when determining a target pattern by aligning to a hole reference!")

	targetNo = 0

	userInput = 0
	while userInput == 0 and (not targetPattern or usePolygon == 1):	# Only collect Preview image for non-target pattern or polygon setups
		if useSearch: 
			sem.Search()
		else:
			sem.V()
		drag()

	targetNo += 1
	target = {}

	# save coords
	sem.GoToLowDoseArea("R")
	ISX0, ISY0, *_ = sem.ReportImageShift()
	SSX0, SSY0 = sem.ReportSpecimenShift()
	stageX, stageY, stageZ = sem.ReportStageXYZ()
	target["stageX"], target["stageY"] = sem.AdjustStagePosForNav(stageX, stageY, ISX0, ISY0)
	target["SSX"], target["SSY"] = [0, 0]
	target["viewfile"] = userName + "_tgt_001_view.mrc"

	saveNewTarget(tgtsFilePath, targetNo, target)

	if drawBeam or usePolygon == 1:
		if beamR == 0:
			beamR = sem.ReportIlluminatedArea() * 100 / 2
		if drawBeam:
			beamPolygons.append(drawBeamPolygon(target["stageX"], target["stageY"], stageZ, beamR, maxTilt))

	# make view map tor realign to item
	sem.SetCameraArea("V", "F")
	sem.V()
	sem.OpenNewFile(userName + "_tgt_001_view.mrc")
	sem.S("A")
	sem.NewMap(0, userName + "_tgt_001_view.mrc")
	sem.CloseFile()
	sem.RestoreCameraSet("V")

	if targetPattern:
		if np.linalg.norm(vecA) == 0 and usePolygon == 0:
			holeDiameter = sem.EnterDefaultedNumber(1.2, 1, "Please enter the hole diameter in microns!")
			foundVecs = False
			if holeDiameter > 0:
				foundVecs, vecA, vecB = vecByXcorr(holeDiameter)

			if not foundVecs:
				sem.Echo("Please center the neighboring hole by dragging the image using the right mouse button and press the <b> key when finished!")
				sem.OKBox("\n".join(["Please center the neighboring hole [2] by dragging the image using the right mouse button!","","Press the <b> key when finished!","","Hole pattern:","0 0 0","0 1 2 <=","0 0 0"]))
				while not sem.KeyBreak():
					sem.Delay(0.1)
				sem.GoToLowDoseArea("R")
				SSX, SSY = sem.ReportSpecimenShift()
				SSX -= SSX0
				SSY -= SSY0
				vecA = (SSX, SSY)
				vecB = (-vecA[1], vecA[0])

		if np.linalg.norm(vecA) == 0 and usePolygon == 1:
			dist = [2 * beamR, 2 * beamR / np.cos(np.radians(maxTilt))]
			theta = np.arctan(np.tan(np.radians(patternRot)) * np.cos(np.radians(maxTilt)))
			rotM = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
			size = 10
			vecA = (rotM @ np.array([1, 0])) * dist
			vecB = (rotM @ np.array([0, 1])) * dist

		if alignToP:											# refine grid vectors by aligning to hole reference in P
			sem.SetImageShift(ISX0, ISY0)								# reset IS to central hole after dragging
			sem.Echo("Vector A: " + str(vecA))

			shiftx = size * vecA[0]
			shifty = size * vecA[1]
			sem.ImageShiftByMicrons(shiftx, shifty)

			sem.V()
			sem.AlignTo("P")
			sem.GoToLowDoseArea("R")

			SSX, SSY = sem.ReportSpecimenShift()
			SSX -= SSX0
			SSY -= SSY0		

			vecA = (round(SSX / size, 4), round(SSY / size, 4))

			sem.Echo("Refined vector A: " + str(vecA))

			sem.ImageShiftByMicrons(-SSX, -SSY)							# reset IS to center position

			sem.Echo("Vector B: " + str(vecB))

			shiftx = size * vecB[0]
			shifty = size * vecB[1]
			sem.ImageShiftByMicrons(shiftx, shifty)

			sem.V()
			sem.AlignTo("P")
			sem.GoToLowDoseArea("R")

			SSX, SSY = sem.ReportSpecimenShift()
			SSX -= SSX0
			SSY -= SSY0		

			vecB = (round(SSX / size, 4), round(SSY / size, 4))

			sem.Echo("Refined vector B: " + str(vecB))

		output = ""
		if usePolygon != 1:
			output += "_set size = " + str(size) + "\n"
		output += "_set vecA0 = " + str(vecA[0]) + "\n"
		output += "_set vecA1 = " + str(vecA[1]) + "\n"
		output += "_set vecB0 = " + str(vecB[0]) + "\n"
		output += "_set vecB1 = " + str(vecB[1]) + 2 * "\n"


		for i in range(-size,size+1):
			for j in range(-size,size+1):
				if i == j == 0: continue

				SSX = i * vecA[0] + j * vecB[0]
				SSY = i * vecA[1] + j * vecB[1]

				stageShift = ss2sMatrix @ np.array([SSX, SSY])

				if usePolygon == 1 and not polygon.contains_points([(target["stageX"] + stageShift[0], target["stageY"] + stageShift[1])])[0]:
					continue

				targetNo += 1

				output += "_tgt = " + str(targetNo).zfill(3) + "\n"
				output += "tsfile = " + userName + "_ts_" + str(targetNo).zfill(3) + ".mrc" + "\n"
				output += "SSX = " + str(SSX) + "\n"
				output += "SSY = " + str(SSY) + "\n"
				output += "stageX = " + str(target["stageX"] + stageShift[0]) + "\n"
				output += "stageY = " + str(target["stageY"] + stageShift[1]) + "\n"
				output += "skip = False" + 2 * "\n"

				ptIndex = int(sem.AddStagePosAsNavPoint(target["stageX"] + stageShift[0], target["stageY"] + stageShift[1], stageZ))
				#ptIndex = int(sem.AddStagePosAsNavPoint(target["stageX"] + ss2sMatrix[0, 0] * SSX + ss2sMatrix[0, 1] * SSY, target["stageY"] + ss2sMatrix[1, 0] * SSX + ss2sMatrix[1, 1] * SSY, stageZ))
				sem.ChangeItemLabel(ptIndex, str(targetNo).zfill(3))

				sem.Echo("Target " + str(targetNo).zfill(3) + " (" + userName + "_tgt_" + str(targetNo).zfill(3) + ".mrc) with image shifts " + str(SSX) + ", " + str(SSY) + " was added.")

		with open(tgtsFilePath, "a") as f:
			f.write(output)	

	else:													# loop over other targets
		loopAddTargets()
		sem.SetImageShift(ISX0, ISY0)

if beamR == 0:													# in case beamR was not defined during target selection
	sem.GoToLowDoseArea("R")
	beamR = sem.ReportIlluminatedArea() * 100 / 2

reopen = True
while reopen:
	gui(tgtsFilePath)											# open GUI after selection is done

sem.SetImageShift(0,0)
if len(beamPolygons) > 0:
	beamPolygons.reverse()
	for polyID in beamPolygons:											# needs to be reversed to keep IDs consistent during deletion
		sem.DeleteNavigatorItem(polyID)

sem.Echo("Target selection completed! " + str(targetNo) + " targets were selected.")
if guidance: 
	sem.OKBox("Target selection completed! " + str(targetNo) + " targets were selected.")
sem.Exit()
