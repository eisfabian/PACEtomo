# beta
This folder contains scripts and versions that have not been extensively tested yet. Your feedback is appreciated!

### PACEtomo_selectTargets.py [v1.4]
- Small update to the target selection script, which still needs testing on a Krios.
- Changes:
  - Draws beam at tilted stage (ellipse) around targets during manual selection using the *beamDiameter* and the *maxTilt* settings. If the *beamDiameter* is set to 0 the script will attempt to read it from the microscope illuminated area value, which is only available on Thermo Scientific microscopes like the Krios. **WARNING:** To draw the beam diameter this script will add a polygon item to the navigator file and reload the navigator. Be careful and maybe make a backup of the navigator file to be safe!
  - New warning when you select a target close to the SerialEM image shift limit (default: 15 microns).
