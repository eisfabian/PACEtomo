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
