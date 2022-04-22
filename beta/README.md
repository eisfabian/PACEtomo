# beta
This folder contains scripts and versions that have not been extensively tested yet. Your feedback is appreciated!

### PACEtomo_selectTargets.py [v1.3]
- Added option to use SerialEM's Low Dose Search instead of View to allow for a different field of view during target selection. Search is not used for a *targetPattern* setup and the map for realignment of target 1 is still taken in View mode.
- Added option to run target selection on a group of points:
  - Use *Add Points* in the Navigator to select points on a montage or a low mag image.
  - Select the **first** point of the group and run the script in manual mode (*targetPattern* and *targetByShift* should be *False*).
  - You will be asked if you want to use the coordinates of the points in the group as targets to be refined by dragging.
  - The rest of the process is the same with the exception that the View images are taken at the coordinates of the next point after a target was selected.
  - You can skip points by saying no when asked to take a Preview image. 
  - You can add additional targets after going through all points.
