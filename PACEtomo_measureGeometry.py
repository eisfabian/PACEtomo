#!Python
# ===================================================================
#ScriptName     PACEtomo_measureGeometry
# Purpose:      Measures pretilt and rotation for PACEtomo.
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2022/04/22
# Revision:     v1.9
# Last Change:  2024/04/25: adopted fixes from selectTargets script (1.8)
# ===================================================================
    
import serialem as sem
import numpy as np

sem.OKBox("This script is deprecated. Please use the measure geometry function in the selectTargets script GUI!")

groupInfo = sem.ReportGroupStatus()
groupPoints = groupInfo[2]
if groupPoints < 3:
    sem.OKBox("You need at least 3 points in the group!")
    sem.Exit()
firstPoint = sem.ReportNavItem()
navLabel = sem.GetVariable("navLabel")
sem.Pause("Please make sure that the first point of the group is selected! (Selected point: " + navLabel + ", points in group: " + str(int(groupPoints)) + ")")

sem.GoToLowDoseArea("R")
s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
sem.SetImageShift(0,0)
sem.MoveToNavItem(int(firstPoint[0]))
points = [[], [], []]
for i in range(int(firstPoint[0]), int(firstPoint[0] + groupPoints)):
    point = sem.ReportOtherItem(i)
    #coords = [s2ssMatrix[0] * (point[1] - firstPoint[1]), s2ssMatrix[3] * (point[2] - firstPoint[2])]
    coords = s2ssMatrix @ np.array([point[1] - firstPoint[1], point[2] - firstPoint[2]])
    sem.ImageShiftByMicrons(coords[0], coords[1])
    sem.G(-1)
    d, *_ = sem.ReportAutoFocus()
    points[0].append(coords[0])
    points[1].append(coords[1])
    points[2].append(d)
    sem.SetImageShift(0,0)    

##########
# Source: https://math.stackexchange.com/q/99317
# subtract out the centroid and take the SVD
svd = np.linalg.svd(points - np.mean(points, axis=1, keepdims=True))
# Extract the left singular vectors
left = svd[0]
# the corresponding left singular vector is the normal vector of the best-fitting plane
norm = left[:, -1]
##########

print("Fitted plane into cloud of " + str(int(groupPoints)) + " points.")

print("Normal vector: " + str(norm))
sign = 1 if norm[1] <= 0 else -1
tilty = sign * np.degrees(np.arccos(norm[2]))
print("Estimated pretilt: " + str(round(tilty, 1)) + " degrees")
rotation = -np.degrees(np.arctan(norm[0]/norm[1]))
print("Estimated rotation: " + str(round(rotation, 1)) + " degrees")

sem.Exit()