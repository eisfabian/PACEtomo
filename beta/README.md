# beta
This folder contains scripts and versions that have not been extensively tested yet. Your feedback is appreciated!

### PACEtomo [v1.2]
- New version of the PACEtomo acquisition script.
- Changes:
  - Independent abort of a target's tilt series branches in case of image shift approaching the limit or optionally in case of dark images (*minCounts* setting).
  - Numerical settings in the script can now be overwritten by settings in the target file (*rootname_tgts.txt*). This allows for varying settings during batch acquisition of several PACEtomo areas. To overwrite a setting add a line like ```_set varName = numericalValue``` at the beginning or end of the targets file.
  - Option to apply additional tilt axis offsets for positive and negative branches as used for the side-entry holder dataset in the manuscript (*taOffsetPos* and *taOffsetNeg* settings). These offsets are applied on top of the global offset set in SerialEM and are only used for the internal calculations.
  - Additional consideration of the *rotation* value should improve performance especially when estimating the geometry of a support film.
  - Changed alignment buffer from M to O to allow for more Roll buffers.
  - Minor text fixes.

### PACEtomo_selectTargets.py [v1.4]
- Small update to the target selection script, which still needs testing on a Krios.
- Changes:
  - Draws beam at tilted stage (ellipse) around targets during manual selection using the *beamDiameter* and the *maxTilt* settings. If the *beamDiameter* is set to 0 the script will attempt to read it from the microscope illuminated area value, which is only available on Thermo Scientific microscopes like the Krios. **WARNING:** To draw the beam diameter this script will add a polygon item to the navigator file and reload the navigator. Be careful and maybe make a backup of the navigator file to be safe!
  - New warning when you select a target close to the SerialEM image shift limit (default: 15 microns).

### PACEtomo_measureGeometry.py [v1.1]
- This script uses a group of navigator points, measures the defocus at these points and estimates pretilt and rotation of the sample. These values should give you an idea about the general geometry of your sample, but don't account for sample deformations. The more points you select, the better the fit should get. I usually use 5-9. Please let me know if you give it a try!
- How to use:
  - Use *Add Points* in the Navigator to select points (at least 3) on a montage or a low magnification image.
  - The first point will be used as the stage position so it should me somewhat centred.
  - Make sure to not surpass the beam shift limits of the microscope (usually within 15 μm).
  - Select the **first** point of the group and run the script.

### PACEtomo_measureOffset.py [v1.0]
- This script should give you a tilt axis offset estimated by movement along the z-axis during tilting (Thanks to Wim Hagen for the suggestion!). It usually gives results within 0.1-0.2 microns of the optimal position for PACEtomo and I adjust it depending on the focus slopes I observe during my PACEtomo runs. Please let me know if you give it a try!
- How to use:
  - You can use the default values or adjust the maximum tilt and tilt increment for more or less datapoints used in the estimation.
  - Make sure that there are no dark images or holes where the script measures the defoci.
  - It will show you 3 tilt axis offsets (on tilt axis and with +/- the selected offset) and an average tilt axis offset.
  - Set the tilt axis offset in SerialEM (-> Tasks -> Eucentricity -> Set Tilt Axis Offset).
  - Make sure "Center image shift on tilt axis" is checked in the "Image Alignment & Focus" window.
  - Run the script again to see if there is a remaining offset.
