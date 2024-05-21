#!Python
# ===================================================================
#ScriptName     PACEtomo_measureOffset
# Purpose:      Estimates tilt axis offset for PACEtomo. (Thanks to Wim Hagen for the suggestion!)
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2022/05/10
# Revision:     v1.8
# Last Change:  2024/05/10: forced "Center Shift on Tilt axis" (1.8)
# ===================================================================

############ SETTINGS ############ 

increment   = 5     # tilt step
maxTilt     = 15    # maximum +/- tilt angle
offset      = 5     # +/- offset for measured positions in microns from tilt axis (also accepts lists e.g. [2, 4, 6])

plot        = True  # plot measurements

########## END SETTINGS ########## 

import serialem as sem
import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt

########### FUNCTIONS ###########

def dZ(alpha, y0):
    return y0 * np.tan(np.radians(-alpha))

def Tilt(tilt):
    sem.TiltTo(tilt)

    for i in range(len(offsets)):
        sem.ImageShiftByMicrons(0, offsets[i])
        sem.G(-1)
        defocus, *_ = sem.ReportAutoFocus()
        focus[i].append(float(defocus))
        sem.SetImageShift(0, 0)

    if tilt == 0:
        for j in range(len(offsets)): 
            focus0.append(focus[j][-1])

    angles.append(float(tilt))
    
###########################

sem.ResetClock()
sem.SuppressReports()
sem.SetUserSetting("ShiftToTiltAxis", 1, 1)

oldOffset = sem.ReportTiltAxisOffset()[0]

sem.Echo("Currently set tilt axis offset: " + str(oldOffset))

sem.Echo("##### Starting tilt axis offset estimation #####")
sem.Echo("Rough eucentricity...")
sem.Eucentricity(1)

sem.Echo("Autofocus...")
sem.G()

sem.Echo("Start tilt series...")
starttilt = -maxTilt
sem.TiltTo(starttilt)
sem.TiltBy(-increment)

offsets = [0]
if isinstance(offset, (list, tuple)):
    for val in offset:
        offsets.extend([-val, val])
else:
    offsets.extend([-offset, offset])

angles = []
focus = [[] for i in range(len(offsets))]
focus0 = []

steps = 2 * maxTilt / increment + 1

tilt = starttilt
for i in range(int(steps)):
    sem.Echo("Tilt to " + str(tilt) + " deg")
    Tilt(tilt)
    tilt += increment

relFocus = focus
for i in range(len(angles)):
    for j in range(len(offsets)):
        relFocus[j][i] -= focus0[j]

y0 = np.zeros(len(offsets))
y0_neg = np.zeros(len(offsets))
y0_pos = np.zeros(len(offsets))
for j in range(len(offsets)):
    y0[j], cov = optimize.curve_fit(dZ, angles, relFocus[j], p0=0)
    y0_neg[j], cov = optimize.curve_fit(dZ, [angle for angle in angles if angle <= 0], relFocus[j][:len([angle for angle in angles if angle <= 0])], p0=0)
    y0_pos[j], cov = optimize.curve_fit(dZ, [angle for angle in angles if angle >= 0], relFocus[j][len([angle for angle in angles if angle < 0]):], p0=0)

sem.Echo("Remaining tilt axis offsets:")
for i in range(0, len(offsets)):
    sem.Echo("[" + str(offsets[i]) + "]: " + str(round(y0[i] + offsets[i], 2)) + " (neg: " + str(round(y0_neg[i] + offsets[i], 2)) + ", pos: " + str(round(y0_pos[i] + offsets[i], 2)) + ")")
avgOffset = sum(y0) / len(offsets)
avgOffset_neg = sum(y0_neg) / len(offsets)
avgOffset_pos = sum(y0_pos) / len(offsets)
sem.Echo("Average remaining tilt axis offset: " + str(round(avgOffset, 2)) + " (neg: " + str(round(avgOffset_neg, 2)) + ", pos: " + str(round(avgOffset_pos, 2)) + ")")
totalOffset = round(avgOffset + oldOffset, 2)
sem.Echo("##############################################")
sem.Echo("Estimated total tilt axis offset: " + str(totalOffset))
sem.Echo("##############################################")

sem.TiltTo(0)
sem.ResetImageShift()

sem.SuppressReports(0)
sem.ReportClock()

if plot:
    offsets, relFocus = zip(*sorted(zip(offsets, relFocus)))    # ensure right order for plot points
    fig = plt.figure(figsize=(8, 6), tight_layout=True)
    plt.title('Z Shifts [microns]')
    for i in range(len(angles)):
        values = []
        for j in range(len(offsets)):
            values.append(relFocus[j][i])
        plt.plot(offsets, values, label=str(angles[i]) + " deg")

    plt.legend()
    plt.show()

userInput = sem.YesNoBox("The estimated total tilt axis offset is " + str(totalOffset) + ". Do you want to set the new tilt axis offset?")
if userInput == 1:
    sem.SetTiltAxisOffset(totalOffset)
    sem.Echo("The new tilt axis offset has been set!")
sem.Exit()