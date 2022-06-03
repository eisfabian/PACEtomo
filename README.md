# PACE-tomo
Parallel cryo electron tomography (PACE-tomo) is a SerialEM script written in Python that allows for the collection of an arbitrary number of tilt series in parallel via beam image shift.
Please refer to the [publication](https://doi.org/10.1101/2022.04.07.487557) for more details.

## Hardware

PACE-tomo has been tested on a JEOL JEM-F200 and a Thermo Fisher Scientific Krios G3i and G4 transmission electron microscope using the dedicated computer for the Gatan K2 or K3 direct electron detector. The open-source microscope control software [SerialEM](https://bio3d.colorado.edu/SerialEM/) should be installed and calibrated.

## Requirements
PACE-tomo does not require the installation of any stand-alone software. However, it does require SerialEM 4.0 or higher capable of running Python scripts.

You can run the following lines of code in a SerialEM script window to test if Python is configured correctly:
```python
#!Python
import serialem as sem
sem.OKBox("Python works!")
```
If you get an error message, please consult the [SerialEM website](https://bio3d.colorado.edu/SerialEM/hlp/html/about_scripts.htm#Python) on how to setup Python for SerialEM.

Additionally, you require the numpy and scipy modules for Python. Depending on your Python installation the following commands in the Windows Command Prompt should take care of it (you will need network connection):

	pip install numpy
	pip install scipy

To check if SerialEM has access to the modules, run this script inside a SerialEM script window:
```python
#!Python
import serialem as sem
import numpy as np
from scipy import optimize
sem.OKBox("All necessary modules are installed!")
```
To use PACE-tomo just copy the content of *PACEtomo.py* and *PACEtomo_selectTargets.py* as well as any auxiliary scripts you want to use in an empty SerialEM script slot each.

## Usage
For optimal results, load the sample such that the lamellae are oriented with the milling direction perpendicular to the tilt axis ([schematic for a Krios autoloader](gridOrientation.png)). Set up SerialEM low dose mode like you would for conventional tilt series acquisition. Make sure buffer O is outside the range of Roll Buffers (Buffer controls window). Offsets for *Focus* and *Trial* areas can be set to 0 to preserve specimen area. Make sure to set the appropriate [tilt axis offset](https://bio3d.colorado.edu/SerialEM/hlp/html/menu_tasks.htm#hid_tasks_settiltaxisoffset) (more details [below](#pacetomo_measureoffsetpy)). It is recommended to do a [coma-free alignment](https://bio3d.colorado.edu/SerialEM/hlp/html/menu_focus.htm#hid_focus_coma_by_ctf) and a [coma vs image shift calibration](https://bio3d.colorado.edu/SerialEM/hlp/html/menu_calibration.htm#hid_focustuning_comavs) to minimise beam tilt for large image shifts (you might need a carbon film to get nice power spectra for CTF fitting). However, in most cases beam tilt will not be resolution limiting.

### Target selection

Before you run a PACE-tomo acquisition, you must define the targets using the *PACEtomo_selectTargets* script. 

**Caution:** To avoid problems with the SerialEM working directory, please choose the folder in which you save your frames during acquisition before running the *selectTargets* script.

There are 3 ways to define targets:

1. **Manually: Selecting targets by dragging the image and centring features of interest.**
	- Inside the script, choose a *delaytime* for dragging the image before taking the next image.
	- By default the script will use View mode to find targets. If you set *useSearch* to *True*, it will use Search mode instead.
	- Set all other settings to *False*.
	- If you want to choose your targets freely, just proceed with the next step. Alternatively, you can use a group of points as template for target selection (not in the video tutorial):
		- Use *Add Points* in the Navigator to select points on a montage or a low magnification image.
  		- Select the **first** point of the group and run the script as described below.
	- Move the stage to your first target (tracking target), which should have enough contrast to be tracked confidently.
	- In the SerialEM UI in the *Image Alignments & Focus* window, uncheck *Move stage for big mouse shifts*.
	- Run the script from the script window.
	- Choose the folder where all files related to PACE-tomo including the final tilt series are saved. This should be the same folder for all PACE-tomo acquisition areas if you plan to run them in batch via *Acquire at items*!
	- Choose a rootname for the current acquisition area. All files related with this acquisition area will be named accordingly.
	- The script will guide you step-by-step through the following process via text windows: <img width="400" src="selectTargets.png" alt="Target selection process" />
	- When dragging the image, make sure not to hit the "Shift" key as this will trigger stage movement.
	- The script finishes when you do not add any more targets.
 
2. **Fixed shifts: Selecting targets by specifying relative image shifts in specimen coordinates.**
	- The overall process is like 1., but instead of dragging to centre a target, you supply shifts in µm for X and Y that are applied from the current position to reach the next target. This is useful for (semi-)ordered patterns of targets.
	- Set only *targetByShift* to *True*. Set all other settings to *False*.
	- Move the stage to your first target (tracking target).
	- Run the script from the script window.
	- The script will guide you through the process.

3. **Pattern: Selecting a target pattern that can be applied to arbitrary stage positions.**
	- A target pattern is useful for the collection on regular holey support films and can be easily transferred to other stage positions.
	- Set *targetPattern* to *True*.
	- If you have a hole reference saved in buffer P and want to refine the manually entered grid vectors (*vecA* and *vecB*) according to hole positions, set *alignToP* to *True*.
	- Set *size* to the appropriate value n for your desired pattern (2n+1)x(2n+1), e.g. 2 for a 5x5 pattern.
	- Enter specimen shifts in x and y between neighbouring holes in *vecA*. *vecB* assumes a perpendicular pattern, but you can specify custom values as well.
	- One way to determine the specimen shift between holes:
		- Centre hole with stage (make sure image shift is 0)
		- Take view image and drag with the right mouse button to centre the neighbouring hole (make sure stage didn’t move).
		- Run the *ReportSpecimenShift* command in SerialEM and take the values output in the log window.
	- Run the script from the script window.

Once target selection is completed all targets are saved in the navigator and a *rootname_tgts.txt* file is created. Target 1 is set to *Acquire* and the name of the *rootname_tgts.txt* file is saved in its *Note* entry. 

**Caution:** Only target 1 of each PACE-tomo acquisition area should be set to *Acquire*!

### Acquisition

PACE-tomo runs a dose-symmetric tilt scheme with groups of 2 tilt images per branch in all cases. Before starting the PACE-tomo collection, please check the settings inside the *PACEtomo* script. Most settings are self-explanatory, but here is a more detailed description for some of them:

- The *startTilt* in degrees is usually 0 or, in case of a lamella, the compensating tilt for the *pretilt* (in case of a *pretilt* of -10 degrees, a *startTilt* of 10 degrees can be used). The *startTilt* has to be divisible by the tilt *step*.
- The tilt range is given by the *maxTilt* relative to the *startTilt* in degrees. This is usually ± 60 or less.
- If a defocus range is given, PACE-tomo will use different target defoci (separated by *stepDefocus*) for each target. If you want to use the same target defocus, keep *minDefocus* and *maxDefocus* the same.
- If your tilt axis offset is not appropriately set, there will be a pseudo-linear defocus slope throughout your tilt series. You can run PACE-tomo on a carbon film, estimate the defocus by CTF fitting and plot the change in µm per degree. Set this value as *focusSlope* to compensate in subsequent acquisitions. Alternatively, refine the tilt axis offset to minimise the slope. When using SerialEM’s fine eucentricity routine to obtain a tilt axis offset, a significant focus slope remains. You can use the [*PACEtomo_measureOffset.py*](#pacetomo_measureoffsetpy) script to get a PACE-tomo optimised estimate for the tilt axis offset.
- You can set delays to be applied after adjusting the image shift and after tilting. On modern state-of-the-art microscopes and resolutions typical of subtomogram averaging, such delays should not be necessary.
- The *pretilt* of the lamella (if applicable) is determined during the focused ion beam milling process and is usually between 8-12 degrees. The sign is important and depends on the orientation in which the grid was loaded into the microscope. For example in case of a FIB milling angle of 10 degrees: If the lamella appears thinner/brighter at +10 degrees stage tilt angle, the pretilt value should be -10 degrees and vice versa. (It is recommended to load lamella containing grids consistently in the same orientation.) You can use the [*PACEtomo_measureGeometry.py*](#pacetomo_measuregeometrypy) script to get a rough estimate of your lamella geometry.
- Lamellae should be oriented with the milling direction perpendicular to the tilt axis during sample loading. In this case the *rotation* should be 0 degrees. If there is a residual rotation you can estimate and enter it for the initial estimation of the eucentric offset. Again, the sign depends on the orientation of the lamella in the microscope.
- If you want to use PACE-tomo on a regular target pattern (e.g., holey support film), set *tgtPattern* to *True*. Additionally, set *alignToP* to *True* if you saved a hole template to buffer P to use for target alignment.
- *beamTiltComp* should be set to *True* if you did the [coma vs image shift calibration](https://bio3d.colorado.edu/SerialEM/hlp/html/menu_calibration.htm#hid_focustuning_comavs).
- By setting *addAF* to *True* you can add additional autofocus routines on target 1 at every branch switch of the dose-symmetric tilt series. This should help keeping the defocus spread low at the cost of overexposing the tracking tilt series.
- *previewAli* can be set to *True* if you want to align every target to its saved Preview image. This helps to keep your target centred if your *startTilt* is not 0. If your field of view is large and your feature of choice does not have to be centred precisely, *previewAli* can be set to *False* to reduce the initial dose on your targets. *previewAli* can also be useful when transferring a “targetPattern” to a different stage position, where the holey support film might have slightly different grid vectors. In this case you should also have *alignToP* set to *True*.
- SerialEM has a hard limit on applying image shifts, which is 15 µm by default. *imageShiftLimit* will overwrite this SerialEM property. The maximum amount of image shift is system dependent and the limit for Thermo Scientific TEM systems is always somewhere below 25 µm.
- The number of *dataPoints* used for the calculation of the eucentric offset of each target was kept at 4 throughout all experiments. Changes could be beneficial to performance.
- The *alignLimit* should keep the cross-correlation alignment in check in cases of low contrast. On good stages the alignment error should never be worse than 0.5 µm.
- If you set *minCounts* greater than 0, a branch of a target can be terminated independently if the image mean counts were below the threshold. The counts are considered per second of exposure and the *ReportExposure* command in SerialEM 4.0+ is used to obtain the exposure time.
- *ignoreFirstNegShift* usually improves the alignment of the first tilt images from the negative branch and should generally be set to *True*.
- *slowTilt* should only be set to *True* if you need additional tilt backlash corrections for the positive tilt branch, which should not be necessary for good stages.
- *taOffsetPos* and *taOffsetNeg* allow you to apply additional tilt axis offsets for positive and negative branches, respectively (as used for the side-entry holder dataset in the manuscript). These offsets are applied on top of the global offset set in SerialEM and are only used for the internal calculations.

Any numerical setting in the script can be overwritten by settings in the target file (*rootname_tgts.txt*). This allows for varying settings during batch acquisition of several PACE-tomo areas. To overwrite a setting add a line like ```_set varName = numericalValue``` at the beginning or end of the targets file.

You can run the PACE-tomo acquisition script either by selecting the entry of target 1 in the Navigator (its Note entry contains *rootname_tgts.txt*) and pressing *Run* in the script window or you can run it in batch via the *Acquire at Items...* dialogue. In case you want to run PACE-tomo on a regular grid *targetPattern*, you can copy the *Note* entry from target 1 it was defined on to any other point in the Navigator, allowing for batch PACE-tomo acquisition. 

### Output
All files are created in the folder you specified during target selection. Target images have the suffix *tgt_xxx* and collected tilt series have been saved with the suffix *ts_xxx* and their accompoanying *.mdoc* file.

### Auxiliary scripts [not shown in video tutorials]

#### [*PACEtomo_measureOffset.py*](https://github.com/eisfabian/PACEtomo/blob/main/PACEtomo_measureOffset.py)
SerialEM usually estimates the tilt axis offset using y-displacements measured during the fine eucentricity routine, which yields suboptimal results for PACE-tomo. This script will estimate the tilt axis offset optimized for movement along the z-axis during tilting <sub>(Thanks to Wim Hagen for the suggestion!)</sub>. It will use SerialEM's autofocus routine to measure the z-height at 3 positions (on tilt axis and ± the given *offset*) throughout a limited tilt series (given by *increment* and *maxTilt*). The results should be within 0.1-0.2 µm of the optimal position for PACE-tomo and you can adjust it depending on the focus slope you observe during a PACEtomo run.
- How to use:
  - You can use the default values or adjust the *maxTilt* and tilt *increment* for more or less datapoints used in the estimation.
  - Make sure that there are no dark images or holes where the script runs the autofocus routine.
  - It will show you 3 tilt axis offsets for the different positions and an average tilt axis offset. The results are relative to the offset already set in SerialEM.
  - Set the tilt axis offset in SerialEM (-> Tasks -> Eucentricity -> Set Tilt Axis Offset).
  - Make sure "Center image shift on tilt axis" is checked in the "Image Alignment & Focus" window.
  - Run the script again to see if there is a remaining offset.

#### [*PACEtomo_measureGeometry.py*](https://github.com/eisfabian/PACEtomo/blob/main/PACEtomo_measureGeometry.py)
This script uses a group of navigator points, measures the z-height using the autofocus routine at these points and estimates *pretilt* and *rotation* values for the sample assuming all points are on a plane. These values should give you an idea about the general geometry of your sample, but don't account for sample deformations. The more points you select, the better the fit should get (5-9 data points should be sufficient).
- How to use:
  - Use Add Points in the Navigator to select points (at least 3) on a montage or a low magnification image. **Caution:** Make sure not to expose areas you want to image later!
  - The first point will be used as the stage position so it should me somewhat centred.
  - Make sure to not surpass the beam shift limits of the microscope (usually within 20 μm).
  - Select the first point of the group and run the script.

### Video Tutorials

[![PACE-tomo setup and collection on a lamella](https://img.youtube.com/vi/NY3mjphVGfA/0.jpg)](https://www.youtube.com/watch?v=NY3mjphVGfA)
[![PACE-tomo setup on a holey film](https://img.youtube.com/vi/LtGu3t6dkfk/0.jpg)](https://www.youtube.com/watch?v=LtGu3t6dkfk)

### Troubleshooting
- The *rootname_tgts.txt* was not found: Run the *PACEtomo_selectTargets* script again to the point where you select the folder, then cancel it.
- Image shift limits exceeded after start tilt images were collected: Double check if you want to collect targets with such high image shifts. If yes, change the *imageShiftLimit* setting in the PACE-tomo acquisition script accordingly.
- to be continued...

### Recent changes
Please also check the [beta folder](https://github.com/eisfabian/PACEtomo/tree/main/beta) for the latest updates!

#### PACEtomo [v1.2]
- Added independent abort of a target's tilt series branches in case of image shift approaching the limit or optionally in case of dark images (*minCounts* setting).
- Numerical settings in the script can now be overwritten by settings in the target file. This allows for varying settings during batch acquisition of several PACE-tomo areas.
- Added option to apply additional tilt axis offsets for positive and negative branches as used for the side-entry holder dataset in the manuscript (*taOffsetPos* and *taOffsetNeg* settings). These offsets are applied on top of the global offset set in SerialEM and are only used for the internal calculations.
- Additional consideration of the *rotation* value should improve performance especially when estimating the geometry of a support film.
- Changed alignment buffer from M to O to allow for more Roll buffers.
- Minor text fixes.

#### PACEtomo_selectTargets.py [v1.3]
- Added option to use SerialEM's Low Dose Search instead of View to allow for a different field of view during target selection. Search is not used for a *targetPattern* setup and the map for realignment of target 1 is still taken in View mode.
- Added option to run target selection on a group of points.
- Minor text fixes.
