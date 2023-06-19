# beta
This folder contains scripts and versions that have not been extensively tested yet. Your feedback is appreciated!

### PACEtomo.py [v1.6]
This update includes mainly options for more robust tracking (e.g. for cryoARMs), CFEG functions and bug fixes.
- Notes:
  - When using lower mag tracking, make sure the record beam settings still cover the camera!
  - The CFEG flashing for TFS instruments was not tested!
- Changes:
  - Added mag offset for tracking tilt series.
  - Added trackTwice option to take a second tracking shot if the alignment was bad.
  - Output file containing all script settings.
  - Fixed progress not considering skipped targets.
  - Added check for existence of ComaVsIS calibration.
  - Added overwriting of existing tilt series file when restarting an acquisition (not for recovery).
  - Changed initial target alignment of targetPattern collection to only use View mag.
  - Added defocus offset for previewAli to boost contrast.
  - Added specimen coords obtainable from stage coords in tgts file. This allows for more flexible creation of tgts file in dummy mode.
  - Switched from targeted tilt angle to microscope tilt angle value when predicting specimen movement.
  - Added sanity check warning boxes.
  - Added functions for CFEG flashing.
  - Minor text fixes.

### PACEtomo_selectTargets.py [v1.6]
Lots of bug fixes, addition of target setup from polygons and initial grid vector estimation.
- Polygon setup: 
  - Draw the desired polygon using the navigator.
  - Move the stage to your desired tracking position.
  - Prepare the selectTargets script:
    - Set *targetPattern* to *True* and *alignToP* to *False*.
    - You can change the *beamDiameter* setting to spread out the targets.
    - The *maxTilt* setting is used to scale the distance between the targets perpendicular to the tilt axis.
    - You can rotate the grid vectors using the *patternRot* setting.
  - Run the script on the polygon item in the Navigator. 
  - You will be asked if you want to run the polygon setup if the script detected the polygon properly.
  - Select your tracking target as usual by dragging the image and confirming the target selection.
  - A grid of points will be created to fill the polygon starting from the tracking target as origin.
  - Note: No preview images will be saved for your targets and hence, the collection might be slightly off the targets shown. You can use the *Save Views* function in the GUI to save view images for every point and set *viewAli* to *True* in the PACEtomo script to use the view images for initial target alignment.
- Changes:
  - Added initial vector estimation for targetPattern using auto-correlation (needs at least 9 holes in the field of View).
  - Added targetPattern setup to fill a polygon using beamDiameter as distance between targets (needs some testing).
  - Fixed being able to unskip skipped targets and continue run.
  - Added progress bars for some GUI functions.
  - Switched to specimen to stage conversion matrix to avoid coordinate inversion on some systems.
  - Added counter in case tgts file already existed.
  - Added view files being saved to tgts file to use for initial target alignment (optional).
  - Fixed measureGeometry outputting only negative pretilts.
  - Added sampleName as additional prefix for all files.
  - Minor text fixes.

### PACEtomo_targetsFromMontage.py [v0.10]
This script can use a medium mag montage to crop out virtual maps as targets (similar to [py-EM](https://www.nature.com/articles/s41592-019-0396-9)). It is still experimental and works decently well on high-contrast features.
Huge thanks to Zhengyi Yang for doing a lot of testing and troubleshooting!
- How to use:
  - On top of the PACEtomo dependencies, this script requires the [mrcfile](https://pypi.org/project/mrcfile/) and the [scikit-image](https://pypi.org/project/scikit-image/) packages.
  - Collect a montage of your target area.
  - Add points using the navigator (as one group). The first point will be considered the tracking target.
  - Run the script on a point.
  - You can also run the script in the dummy version. Any error/warning messages can be ignored.
    - WARNING: There seems to be cross talk between the DUMMY instance and the real instance of SerialEM when running DUMMY on the microscope computer. I'm not sure if that's a new bug or if that's always been the case. I tested the script in a DUMMY version running on a separate computer.
  - You can run the script on several groups of points using Acquire at Items if you check "Skip stage move to item if possible". You can automate the file naming by setting *noUI* to *True*.

Please let me know if you run into any issues or have any feature requests. If there is enough interest, I will work on this further!
