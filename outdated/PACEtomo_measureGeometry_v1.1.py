#!Python
# ===================================================================
# Name:		PACEtomo_measureGeometry
# Purpose:	Measures pretilt and rotation for PACEtomo.
#		More information at http://github.com/eisfabian/PACEtomo
# Author:	Fabian Eisenstein
# Created:	2022/04/22
# Revision:	v1.1
# Last Change:	2022/05/20: fixed considering rotation in pretilt
# ===================================================================

import serialem as sem
import numpy as np

groupInfo = sem.ReportGroupStatus()
groupPoints = groupInfo[2]
if groupPoints < 3:
	sem.OKBox("You need at least 3 points in the group!")
	sem.Exit()
firstPoint = sem.ReportNavItem()
navLabel = sem.GetVariable("navLabel")
sem.Pause("Please make sure that the first point of the group is selected! (Selected point: " + navLabel + ", points in group: " + str(int(groupPoints)) + ")")

sem.SetImageShift(0,0)
sem.MoveToNavItem(int(firstPoint[0]))
points = [[], [], []]
for i in range(int(firstPoint[0]), int(firstPoint[0] + groupPoints)):
	point = sem.ReportOtherItem(i)
	coords = [-point[1] + firstPoint[1], -point[2] + firstPoint[2]]
	sem.ImageShiftByMicrons(coords[0], coords[1])
	sem.G(-1)
	(d, error) = sem.ReportAutoFocus()
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

#tiltx = -np.degrees(np.arctan(norm[0]))
#print(tiltx)

tilty = -np.degrees(np.arctan(np.sqrt(norm[0:2].dot(norm[0:2]))))
print("Estimated pretilt: " + str(round(tilty, 1)) + " degrees")

rotation = -np.degrees(np.arctan(norm[0]/norm[1]))
print("Estimated rotation: " + str(round(rotation, 1)) + " degrees")
