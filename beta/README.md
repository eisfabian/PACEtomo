# beta
This folder contains scripts and versions that have not been extensively tested yet. Your feedback is appreciated!

### PACEtomo_measureOffset.py [v1.0]
- This script should give you a tilt axis offset estimated by movement along the z-axis during tilting (Thanks to Wim Hagen for the suggestion!). Unfortunately, I could not yet test it on a Krios so I can't tell how accurate the estimation will be. Please let me know if you give it a try!
- How to use:
  - You can use the default values or adjust the maximum tilt and tilt increment for more or less datapoints used in the estimation.
  - Make sure that there are no dark images or holes where the script measures the defoci.
  - It will show you 3 tilt axis offsets (on tilt axis and with +/- the selected offset) and an average tilt axis offset.
  - Set the tilt axis offset in SerialEM (-> Tasks -> Eucentricity -> Set Tilt Axis Offset).
  - Make sure "Center image shift on tilt axis" is checked in the "Image Alignment & Focus" window.
  - Run the script again to see if there is a remaining offset.

### PACEtomo_selectTargets.py [v1.3]
- Added option to use SerialEM's Low Dose Search instead of View to allow for a different field of view during target selection. Search is not used for a *targetPattern* setup and the map for realignment of target 1 is still taken in View mode.
- Added option to run target selection on a group of points:
  - Use *Add Points* in the Navigator to select points on a montage or a low mag image.
  - Select the **first** point of the group and run the script in manual mode (*targetPattern* and *targetByShift* should be *False*).
  - You will be asked if you want to use the coordinates of the points in the group as targets to be refined by dragging.
  - The rest of the process is the same with the exception that the View images are taken at the coordinates of the next point after a target was selected.
  - You can skip points by saying no when asked to take a Preview image. 
  - You can add additional targets after going through all points.
- Minor text fixes.
