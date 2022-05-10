#!Python
# ===================================================================
# Name:		PACEtomo_measureOffset
# Purpose:	Estimates tilt axis offset for PACEtomo. (Thanks to Wim Hagen for the suggestion!)
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2022/05/10
# Revision:	v1.0
# Last Change:	2022/05/10
# ===================================================================

import serialem as sem
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

############ SETTINGS ############ 

increment	= 5		# tilt step
maxTilt		= 15		# maximum +/- tilt angle
offset 		= 5		# +/- offset for measured positions in microns from tilt axis

plot 		= True		# plot measurements

########## END SETTINGS ########## 

########### FUNCTIONS ###########

def dZ(alpha, y0):
	return y0 * np.tan(np.radians(-alpha))

def Tilt(tilt):
	sem.TiltTo(tilt)

	for i in range(0, len(offsets)):
		sem.ImageShiftByMicrons(0, offsets[i])
		sem.G(-1)
		(defocus, RepVal2) = sem.ReportAutoFocus()
		focus[i].append(float(defocus))
		sem.SetImageShift(0, 0)

	if tilt == 0:
		for j in range(0, len(offsets)): 
			focus0.append(focus[j][-1])

	angles.append(float(tilt))
	
###########################

sem.ResetClock()
sem.SuppressReports()

sem.Echo("##### Starting tilt axis offset estimation #####")
sem.Echo("Rough eucentricity...")
sem.Eucentricity(1)

sem.Echo("Autofocus...")
sem.G()

sem.Echo("Start tilt series...")
starttilt = -maxTilt
sem.TiltTo(starttilt)
sem.TiltBy(-increment)

offsets = [-offset, 0, offset]
angles = []
focus = [[], [], []]
focus0 = []

steps = 2 * maxTilt / increment + 1

tilt = starttilt
for i in range(0,int(steps)):
	sem.Echo("Tilt to " + str(tilt) + " deg")
	Tilt(tilt)
	tilt += increment

relFocus = focus
for i in range(0,len(angles)):
	for j in range(0, len(offsets)):
		relFocus[j][i] -= focus0[j]

y0 = [0, 0, 0]
for j in range(0, len(offsets)):
	y0[j], cov = optimize.curve_fit(dZ, angles, relFocus[j], p0=0)

sem.Echo("Tilt axis offsets:")
for i in range(0, len(offsets)):
	sem.Echo("[" + str(offsets[i]) + "]: " + str(round(y0[i][0] + offsets[i], 2)))
avgOffset = sum(y0) / len(offsets)
sem.Echo("##############################################")
sem.Echo("Average tilt axis offset: " + str(round(avgOffset[0], 2)))
sem.Echo("##############################################")

#sem.SetTiltAxisOffset(avgOffset)

sem.TiltTo(0)
sem.ResetImageShift()

sem.SuppressReports(0)
sem.ReportClock()

if plot:
	fig = plt.figure(figsize=(8, 6), tight_layout=True)
	plt.title('Z Shifts [microns]')
	for i in range(0, len(angles)):
		values = []
		for j in range(0, len(offsets)):
			values.append(relFocus[j][i])
		plt.plot(offsets, values, label=str(angles[i]) + " deg")

	plt.legend()
	plt.show()
