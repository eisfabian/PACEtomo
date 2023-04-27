#!Python
# ===================================================================
#ScriptName	PACEtomo_selectTargets
# Purpose:	Selects targets, saves maps and writes text file with coords for PACEtomo.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2021/04/19
# Revision:	v1.5.1
# Last Change:	2023/04/27: fixed realign to tracking target
#		2022/09/27: bug fixes after Krios testing
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
useSearch 	= False		# use Search mode instead of View mode to find targets by dragging
vecA		= (0, 0)	# vectors for grid pattern [microns specimen shift] are determined automatically...
vecB		= (0, 0)	# ...only change if you want to setup pattern without alignToP reference

########## END SETTINGS ########## 

import serialem as sem
import os
import copy
import glob
import numpy as np
import scipy as sp
import scipy.optimize
import tkinter as tk
import matplotlib.pyplot as plt
import matplotlib.lines
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
	if targetPattern:
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
				sem.GoToLowDoseArea("R")
				if targetByShift:
					shiftx = sem.EnterDefaultedNumber(0, 1, "Enter X shift:")
					shifty = sem.EnterDefaultedNumber(0, 1, "Enter Y shift:")
				else:
					sem.SetImageShift(0, 0)							# use 0,0 instead of ISX0,ISY0 to account for user shift of first target away from point
					shiftx, shifty = coordsRefine[pointNo][0], coordsRefine[pointNo][1]
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
						if pointNo >= groupPoints:					# disable pointRefine when last point of a group is skipped
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
			polygons.append(drawBeamPolygon(target["stageX"], target["stageY"], stageZ, beamR, maxTilt))

		if pointRefine == 1 and pointNo >= groupPoints:
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
			label = tk.Label(self.tw, text=self.text, justify='left', background='#ffffff', relief='solid', borderwidth=1, font=("Courier","10"))
			label.pack(ipadx=1)

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
		plt.axis('equal')
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
				if ind != 0:
					targets[ind]["skip"] = not targets[ind]["skip"]
				updateList()
				listbox.selection_set(ind)
				listbox.see(ind)

				plt.clf()
				plotTargets()
				fig.canvas.draw()
				fig.canvas.flush_events()

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
			fig.canvas.flush_events()

	def openTgt(box=True):												# open preview map with default application
		ind = listbox.curselection()[0]
		if targets[ind]["tgtfile"]:
			os.system("start " + targets[ind]["tgtfile"])
			return True
		else:
			if box:
				tk.messagebox.showwarning(title="File not found", message="Target file was not found!")
			return False

	def openTs(box=True):												# open tilt series if it exists
		ind = listbox.curselection()[0]
		if targets[ind]["tsfile"]:
			os.system("start " + targets[ind]["tsfile"])
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

	def skip():												# skip selected target
		ind = listbox.curselection()[0]	
		if ind != 0:
			targets[ind]["skip"] = not targets[ind]["skip"]
			updateList()
			listbox.selection_set(ind)
			listbox.see(ind)
			plt.clf()
			plotTargets()
			fig.canvas.draw()

	def loadMap():
		mapLabel = tk.simpledialog.askstring("Map Label", "\n".join(["Please enter the navigator label", "of the map you want to load!"]))
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
			for polyID in polygons:
				sem.ChangeItemDraw(polyID, 1)
		sem.Echo("Adding targets...")		
		userInput = 0
		pointRefine = 0
		loopAddTargets()
		sem.SetImageShift(ISX0, ISY0)
		reopen = True											# reopen GUI after finishing
		return

	def saveViews():
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
			sem.SnapshotToFile(0, 0, "0", "JPG", "JPG", viewName.rsplit(".mrc", 1)[0] + ".jpg")
			sem.GoToLowDoseArea("R")
			sem.SetImageShift(ISX0, ISY0)

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
					ptIndex = int(sem.AddStagePosAsNavPoint(targetsTemp[0]["stageX"] - targetsTemp[i]["SSX"], targetsTemp[0]["stageY"] - targetsTemp[i]["SSY"], item[3], groupIndex))
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
			##########
			# Source: https://math.stackexchange.com/q/99317
			# subtract out the centroid and take the SVD, extract the left singular vectors, the corresponding left singular vector is the normal vector of the best-fitting plane
			svd = np.linalg.svd(geoXYZ - np.mean(geoXYZ, axis=1, keepdims=True))
			left = svd[0]
			norm = left[:, -1]
			##########		
			sem.Echo("Fitted plane into cloud of " + str(len(geoPoints)) + " points.")
			sem.Echo("Normal vector: " + str(norm))
			tilty = -np.degrees(np.arctan(np.linalg.norm(norm[0:2])))
			sem.Echo("Estimated pretilt: " + str(round(tilty, 1)) + " degrees")
			rotation = -np.degrees(np.arctan(norm[0]/norm[1]))
			sem.Echo("Estimated rotation: " + str(round(rotation, 1)) + " degrees")
			tk.messagebox.showinfo(title="measureGeometry", message="\n".join(["Geometry measurement:","","Estimated pretilt: " + str(round(tilty, 1)) + " degrees","Estimated rotation: " + str(round(rotation, 1)) + " degrees"]))
			entryPretilt.delete(0, tk.END)								# update entry fields
			entryPretilt.insert(0, str(round(tilty, 1)))
			entryRotation.delete(0, tk.END)
			entryRotation.insert(0, str(round(rotation, 1)))
			readEntry()

	def toggleBeam():
		nonlocal showBeam
		showBeam = not showBeam
		plt.clf()
		plotTargets()
		fig.canvas.draw()
		if drawBeam:											# also toggle beam polygons in SerialEM
			for polyID in polygons:
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
		if polygons == []:
			for i in range(len(targets)):
				if "stageX" in targets[i].keys():
					polygons.append(drawBeamPolygon(targets[i]["stageX"], targets[i]["stageY"], 0, beamR, maxTilt))
					sem.ChangeItemDraw(polygons[-1], 0)
		else:
			for polyID in polygons:
				sem.ChangeItemDraw(polyID, 0)

	cMatrix = sem.SpecimenToCameraMatrix(0)									# figure out axis of plot to match SerialEM view
	if(abs(cMatrix[0]) > abs(cMatrix[1])):
		swapXY = False
		invX = True if cMatrix[0] < 0 else False
		invY = True if cMatrix[3] < 0 else False		
	else:
		swapXY = True
		invX = True if cMatrix[1] < 0 else False
		invY = True if cMatrix[2] < 0 else False

	# create a root window.
	top = tk.Tk()
	top.option_add("*font", "Courier")
	top.title("Targets")
	top.columnconfigure(2, weight=2)
	top.rowconfigure(6, weight=2)
	pixel = tk.PhotoImage(width=1, height=1)

	# create target list
	fileNameLen = max(len(str(targets[0]["tgtfile"])), len(str(targets[0]["tsfile"])))
	label = tk.Label(top, text=" ".join(["Tgt".ljust(5), "Tgtfile".ljust(fileNameLen + 2), "TSfile".ljust(fileNameLen + 1), "SSX".rjust(6), "SSY".rjust(6), "Skip".ljust(70 - 24 - 2 * fileNameLen)])) 
	label.grid(column=0, row=0, sticky=tk.E)

	listbox = DragDropListbox(top, height=10, selectmode="SINGLE", width=70, activestyle='dotbox')
	updateList()
	listbox.grid(column=0, row=1, rowspan=6, padx=10, sticky=tk.NE)

	# create buttons
	btnWidth = 150
	btnHeight = 20
	
	tgtbutton = tk.Button(top, text="Open Tgtfile", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=openTgt)
	tgtbutton.grid(column=1, row=1, sticky=tk.W, pady=0)
	CreateToolTip(tgtbutton, "Opens .mrc file of selected target preview image.")

	tsbutton = tk.Button(top, text="Open TSfile", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=openTs)
	tsbutton.grid(column=2, row=1, sticky=tk.W, pady=0, padx=5)
	CreateToolTip(tsbutton, "Opens .mrc tilt series stack file of selected target.")

	skipbutton = tk.Button(top, text="Skip", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=skip)
	skipbutton.grid(column=1, row=2, sticky=tk.W, pady=3)
	CreateToolTip(skipbutton, "Toggles skipping of selected target during collection.")

	mtbuttonBorder = tk.Frame(top, highlightbackground="#ffd700", highlightthickness=2)		
	mtbutton = tk.Button(mtbuttonBorder, text="Make Track", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=makeTrack)
	if savedRun: mtbutton["state"] = tk.DISABLED
	mtbutton.pack()
	mtbuttonBorder.grid(column=2, row=2, sticky=tk.W, pady=3, padx=5)
	CreateToolTip(mtbutton, "\n".join(["Makes selected target the tracking target.", "This will cause additional exposures on the selected target for realignment."]))

	lmapbutton = tk.Button(top, text="Load Map", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=loadMap)
	lmapbutton.grid(column=1, row=3, sticky=tk.W, pady=0)
	CreateToolTip(lmapbutton, "Loads map into SerialEM for cross referencing.")

	addbuttonBorder = tk.Frame(top, highlightbackground="#ffd700", highlightthickness=2)
	addbutton = tk.Button(addbuttonBorder, text="Add Targets", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=moreTargets)
	if savedRun: addbutton["state"] = tk.DISABLED
	addbutton.pack()
	addbuttonBorder.grid(column=2, row=3, sticky=tk.W, pady=0, padx=5)
	CreateToolTip(addbutton, "\n".join(["Continues the target selection process.", "This might cause additional exposures on the tracking target for realignment."]))

	viewbutton = tk.Button(top, text="Save Views", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=saveViews)
	viewbutton.grid(column=1, row=4, sticky=tk.W, pady=3)
	CreateToolTip(viewbutton, "Takes and saves view image for each target.")

	#xbutton = tk.Button(top, text="x", image=pixel, compound="center", height=btnHeight, width=btnWidth, command="")
	#if savedRun: xbutton["state"] = tk.DISABLED
	#xbutton.grid(column=2, row=4, sticky=tk.W, pady=3, padx=5)
	#CreateToolTip(xbutton, "Placeholder.")

	rebutton = tk.Button(top, text="Reorder", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=saveOrder)
	if savedRun: rebutton["state"] = tk.DISABLED
	rebutton.grid(column=1, row=5, sticky=tk.W, pady=0)
	CreateToolTip(rebutton, "Applies current order as displayed in the list.")

	cpbutton = tk.Button(top, text="Copy to Acq", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=copyAcq)
	if savedRun or "size" not in settings.keys(): cpbutton["state"] = tk.DISABLED
	cpbutton.grid(column=2, row=5, sticky=tk.W, pady=0, padx=5)
	CreateToolTip(cpbutton, "Applies tgt pattern to all navigator points marked as Acquire.")

	sbutton = tk.Button(top, text="Save", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=saveFile)
	sbutton.grid(column=1, row=6, sticky=tk.SW, pady=3)
	CreateToolTip(sbutton, "Saves all changes to the tgts file.")

	rbutton = tk.Button(top, text="Reset", image=pixel, compound="center", height=btnHeight, width=btnWidth, command=resetOrder)
	rbutton.grid(column=2, row=6, sticky=tk.SW, pady=3, padx=5)
	CreateToolTip(rbutton, "Resets changes made since window was opened.")
	  
	# create settings area
	labelSet = tk.Label(top, text="Settings (optional)")
	labelSet.grid(column=3, columnspan=5, row=0, )

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
	labelPretilt.grid(column=3, row=4, sticky=tk.NE)

	ecPretilt = tk.StringVar(top, settings["pretilt"] if "pretilt" in settings.keys() else "")
	ecPretilt.trace_add("write", readEntry)
	entryPretilt = tk.Entry(top, textvariable=ecPretilt, width=20)
	if savedRun: entryPretilt["state"] = tk.DISABLED
	entryPretilt.grid(column=4, columnspan=4, row=4, sticky=tk.NW)

	labelRotation = tk.Label(top, text = "Rotation [deg]")
	labelRotation.grid(column=3, row=5, sticky=tk.NE)

	ecRotation = tk.StringVar(top, settings["rotation"] if "rotation" in settings.keys() else "")
	ecRotation.trace_add("write", readEntry)
	entryRotation = tk.Entry(top, textvariable=ecRotation, width=20)
	if savedRun: entryRotation["state"] = tk.DISABLED
	entryRotation.grid(column=4, columnspan=4, row=5, sticky=tk.NW)

	# create geometry buttons
	geombutton = tk.Button(top, text="Measure Geometry", image=pixel, compound="center", height=btnHeight, command=measureGeo)
	if savedRun: geombutton["state"] = tk.DISABLED
	geombutton.grid(column=3, row=6, sticky=tk.SE, pady=3, padx=5)
	CreateToolTip(geombutton, "\n".join(["Runs routine to measure pretilt and rotation.", "This will cause exposures on the chosen geometry points."]))

	georbutton = tk.Button(top, text="Reset Geo Pts", image=pixel, compound="center", height=btnHeight, command=resetGeo)
	georbutton.grid(column=4, columnspan=4, row=6, sticky=tk.SW, pady=3, padx=5)
	CreateToolTip(georbutton, "Deletes all geometry points.")

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

	toolbar_frame = tk.Frame(top)										# creating the Matplotlib toolbar
	toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
	toolbar.update()
	toolbar_frame.grid(column=0, row=8, sticky=tk.W, padx=10)

	# create axis buttons
	beambutton = tk.Button(top, text="Toggle beam", image=pixel, compound="center", height=btnHeight, command=toggleBeam)
	beambutton.grid(column=1, row=8, sticky=tk.W, pady=5, padx=5)
	CreateToolTip(beambutton, "\n".join(["Toggles display of beam diameter.", "(Beam diameter is taken from IlluminatedArea or from value set in script.)"]))

	swapbutton = tk.Button(top, text="Swap XY", image=pixel, compound="center", height=btnHeight, command=swapAxes)
	swapbutton.grid(column=2, row=8, sticky=tk.E, pady=5, padx=5)
	CreateToolTip(swapbutton, "Swaps axes of target plot.")

	invXbutton = tk.Button(top, text="Invert X", image=pixel, compound="center", height=btnHeight, command=invertX)
	invXbutton.grid(column=3, row=8, sticky=tk.E, pady=5, padx=5)
	CreateToolTip(invXbutton, "Inverts X axis of target plot.")

	invYbutton = tk.Button(top, text="Invert Y", image=pixel, compound="center", height=btnHeight, command=invertY)
	invYbutton.grid(column=4, columnspan=4, row=8, sticky=tk.W, pady=5, padx=5)
	CreateToolTip(invYbutton, "Inverts Y axis of target plot.")

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
sem.UserSetDirectory("Please choose a directory for saving targets and tilt series!")

imageShiftLimit = sem.ReportProperty("ImageShiftLimit")

navSize = sem.ReportNumTableItems()
if navSize > 0:
	sem.ReportNavItem()												# check if selected nav item already has tgts file
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

polygons = []
beamR = beamDiameter / 2

if editTgts == 0:
	if navSize > 0:
		groupInfo = sem.ReportGroupStatus()
	else:
		groupInfo = [0, 0, 0]
	pointRefine = 0
	if not targetPattern and not targetByShift and groupInfo[1] > 0:
		pointRefine = sem.YesNoBox("\n".join(["GROUP OF POINTS?", "", "The selected navigator item is part of a group. Do you want to use the points of this group as initial target coordinates to be refined?"]))
		if pointRefine == 1:
			groupID = groupInfo[1]
			groupPoints = groupInfo[2]
			firstPoint = sem.ReportNavItem()
			navLabel = sem.GetVariable("navLabel")
			if guidance:
				sem.Pause("Please make sure that the first point of the group is selected! (Selected point: " + navLabel + ", points in group: " + str(int(groupPoints)) + ")")
			sem.Echo("NOTE: Selected point: " + navLabel + ", points in group: " + str(int(groupPoints)))
			coordsRefine = []
			coordsRefine.append([0, 0])
			for i in range(int(firstPoint[0]) + 1, int(firstPoint[0] + groupPoints)):
				point = sem.ReportOtherItem(i)
				coordsRefine.append([-point[1] + firstPoint[1], -point[2] + firstPoint[2]])

	sem.EnterString("userName","Please provide a rootname for the PACE-tomo collection area!")
	userName = sem.GetVariable("userName")

	tgtsFilePath = os.path.join(curDir, userName + "_tgts.txt")

	# align center
	sem.TiltTo(0)
	sem.ResetImageShift()
	if pointRefine == 1:
		sem.MoveToNavItem(int(firstPoint[0]))
	elif alignToP:												# center hole for center of tgtPattern
		sem.V()
		sem.AlignTo("P")
		if float(sem.ReportDefocus()) < -50:
			sem.Pause("WARNING: Large defocus offsets for View can cause image shift offsets when determining a target pattern by aligning to a hole reference!")

	targetNo = 0

	userInput = 0
	while userInput == 0 and not targetPattern:
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

	saveNewTarget(tgtsFilePath, targetNo, target)

	if drawBeam:
		if beamR == 0:
			beamR = sem.ReportIlluminatedArea() * 100 / 2
		polygons.append(drawBeamPolygon(target["stageX"], target["stageY"], stageZ, beamR, maxTilt))

	# make view map tor realign to item
	sem.SetCameraArea("V", "F")
	sem.V()
	sem.OpenNewFile(userName + "_tgt_001_view.mrc")
	sem.S("A")
	sem.NewMap(0, userName + "_tgt_001_view.mrc")
	sem.CloseFile()
	sem.RestoreCameraSet("V")

	if targetPattern:
		if vecA == (0, 0):
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
		output += "_set vecA0 = " + str(vecA[0]) + "\n"
		output += "_set vecA1 = " + str(vecA[1]) + "\n"
		output += "_set vecB0 = " + str(vecB[0]) + "\n"
		output += "_set vecB1 = " + str(vecB[1]) + "\n"
		output += "_set size = " + str(size) + 2 * "\n"

		for i in range(-size,size+1):
			for j in range(-size,size+1):
				if i == j == 0: continue

				targetNo += 1

				SSX = i * vecA[0] + j * vecB[0]
				SSY = i * vecA[1] + j * vecB[1]

				output += "_tgt = " + str(targetNo).zfill(3) + "\n"
				output += "tsfile = " + userName + "_ts_" + str(targetNo).zfill(3) + ".mrc" + "\n"
				output += "SSX = " + str(SSX) + "\n"
				output += "SSY = " + str(SSY) + "\n"
				output += "skip = False" + 2 * "\n"

				ptIndex = int(sem.AddStagePosAsNavPoint(target["stageX"] - SSX, target["stageY"] - SSY, stageZ))
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
if len(polygons) > 0:
	polygons.reverse()
	for polyID in polygons:											# needs to be reversed to keep IDs consistent during deletion
		sem.DeleteNavigatorItem(polyID)

sem.Echo("Target selection completed! " + str(targetNo) + " targets were selected.")
if guidance: 
	sem.OKBox("Target selection completed! " + str(targetNo) + " targets were selected.")
sem.Exit()
