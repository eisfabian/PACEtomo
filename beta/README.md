# beta
This folder contains scripts and versions that have not been extensively tested yet. Your feedback is appreciated!

### PACEtomo_measureGeometry.py [v1.1]
- This script uses a group of navigator points, measures the defocus at these points and estimates pretilt and rotation of the sample. These values should give you an idea about the general geometry of your sample, but don't account for sample deformations. The more points you select, the better the fit should get. I usually use 5-9. Please let me know if you give it a try!
- How to use:
  - Use *Add Points* in the Navigator to select points (at least 3) on a montage or a low magnification image.
  - The first point will be used as the stage position so it should me somewhat centred.
  - Make sure to not surpass the beam shift limits of the microscope (usually within 15 μm).
  - Select the **first** point of the group and run the script.

### PACEtomo_measureOffset.py [v1.0]
- This script should give you a tilt axis offset estimated by movement along the z-axis during tilting (Thanks to Wim Hagen for the suggestion!). Unfortunately, I could not yet test it on a Krios so I can't tell how accurate the estimation will be. Please let me know if you give it a try!
- How to use:
  - You can use the default values or adjust the maximum tilt and tilt increment for more or less datapoints used in the estimation.
  - Make sure that there are no dark images or holes where the script measures the defoci.
  - It will show you 3 tilt axis offsets (on tilt axis and with +/- the selected offset) and an average tilt axis offset.
  - Set the tilt axis offset in SerialEM (-> Tasks -> Eucentricity -> Set Tilt Axis Offset).
  - Make sure "Center image shift on tilt axis" is checked in the "Image Alignment & Focus" window.
  - Run the script again to see if there is a remaining offset.
