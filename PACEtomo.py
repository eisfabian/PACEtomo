#!Python
# ===================================================================
#ScriptName     PACEtomo
# Purpose:      Runs parallel dose-symmetric tilt series on many targets with geometrical predictions using the first target as tracking tilt series.
#               Make sure to run selectTargets script first to generate compatible Navigator settings and a target file.
#               More information at http://github.com/eisfabian/PACEtomo
# Author:       Fabian Eisenstein
# Created:      2021/04/16
# Revision:     v1.9.2
# Last Change:  2025/02/28: fixed tilt range issue, added simplified RealignToItem
# ===================================================================

############ SETTINGS ############ 

startTilt       = 0         # starting tilt angle [degrees] (should be divisible by step)
minTilt         = -60       # minimum absolute tilt angle [degrees]
maxTilt         = 60        # maximum absolute tilt angle [degrees]
step            = 3         # tilt step [degrees]
groupSize       = 2         # group size for grouped dose-symmetric scheme (number of contiguously acquired images on one side of the tilt series)
minDefocus      = -5        # minimum defocus [microns] of target range (low defocus)
maxDefocus      = -5        # maximum defocus [microns] of target range (high defocus)
stepDefocus     = 0.5       # step [microns] between target defoci (between TS)

focusSlope      = 0.0       # [DEPRECATED] empirical linear focus correction [microns per degree] (obtained by linear regression of CTF fitted defoci over tilt series; microscope stage dependent)
delayIS         = 0.3       # delay [s] between applying image shift and Record
delayTilt       = 0.3       # delay [s] after stage tilt
zeroExpTime     = 0         # set to exposure time [s] used for start tilt image, if 0: use same exposure time for all tilt images
zeroDefocus	    = 0 		# set to defocus [microns] used for start tilt image, if 0: use same defocus for all tilt images

# Track settings
trackExpTime    = 0         # set to exposure time [s] used for tracking tilt series, if 0: use same exposure time for all tilt series
trackDefocus    = 0         # set to defocus [microns] used for tracking tilt series, if 0: use same defocus range for all tilt series
trackMag        = 0         # set to nominal magnification for tracking tilt series (make sure detector is still covered under the same beam conditions), if 0: use same mag for all tilt series
trackTwice      = False     # track in 2 steps, useful when large tracking shifts cause inaccuracies in alignment and hence, in residual errors for all targets, but causes double exposure of tracking area

# Geometry settings
pretilt         = 0         # pretilt [degrees] of sample in deg e.g. after FIB milling (if milling direction is not perpendicular to the tilt axis, estimate and add rotation)
rotation        = 0         # rotation [degrees] of lamella vs tilt axis in deg (should be 0 deg if lamella is oriented perpendicular to tilt axis)
measureGeo      = False     # estimates pretilt and rotation values of sample by measuring defocus on geo points saved in target setup or automatically determined points within tgt pattern

# Holey support settings
tgtPattern      = False     # use same tgt pattern on different stage positions (useful for collection on holey support film)
alignToP        = False     # use generic View image in buffer P to align to
refineVec       = False     # refine tgt pattern for local stage position by aligning furthest targets along each axis to to buffer P
refineGeo       = False     # uses on-the-fly CtfFind results of first image to refine geometry before tilting (only use when CTF fits on your sample seem reliable)

# Session settings
beamTiltComp    = True      # use beam tilt compensation (uses coma vs image shift calibrations)
addAF           = False     # does autofocus at the start of every tilt group, increases exposure on tracking TS drastically
previewAli      = True      # adds initial dose, but makes sure start tilt image is on target (uses view image and aligns to buffer P if alignToP == True)
viewAli         = False     # adds an alignment step with a View image if it was saved during the target selection (only if previewAli is activated)

# Output settings
sortByTilt      = True      # sorts tilt series by tilt angle after acquisition is completed (takes additional time), requires mrcfile module
binFinalStack   = 1         # bin factor for final saved stack after acquisition (unbinned stack will be deleted to save storage space)
delFinalStack   = False     # delete final tilt series stacks (only keeps frames for reconstruction to save storage space) 
doCtfFind       = False     # set to False to skip CTFfind estimation (only necessary if it causes crashes => if it does crash, SerialEM will output some troubleshoot data that you should send to David!) 
doCtfPlotter    = True      # runs ctfplotter instead of CTFfind, needs standalone version of 3dmod on PATH
extendedMdoc    = True      # saves additional info to .mdoc file

# Hardware settings
slowTilt        = False     # do backlash step for all tilt angles, on bad stages large tilt steps are less accurate
taOffsetPos     = 0         # additional tilt axis offset values [microns] applied to calculations for positive and...
taOffsetNeg     = 0         # ...negative branch of the tilt series (possibly useful for side-entry holder systems)
checkDewar      = True      # check if dewars are refilling before every acquisition
cryoARM         = False     # if you use a JEOL cryoARM TEM, this will keep the dewar refilling in sync
coldFEG         = False     # if you use a cold FEG, this will flash the gun whenever the dewars are being refilled
flashInterval   = -1        # time in hours between cold FEG flashes, -1: flash only during dewar refill (interval is ignored on Krios, uses FlashingAdvised function instead)
slitInterval    = 0         # time in minutes between centering the energy filter slit using RefineZLP, ONLY works with tgtPattern (needs pattern vectors to find good position for alignment)

# Advanced settings
fitLimit        = 30        # refineGeo: minimum resolution [Angstroms] needed for CTF fit to be considered for refineGeo
parabolTh       = 9         # refineGeo: minimum number of passable CtfFind values to fit paraboloid instead of plane 
imageShiftLimit = 20        # maximum image shift [microns] SerialEM is allowed to apply (this is a SerialEM property entry, default is 15 microns)
dataPoints      = 4         # number of recent specimen shift data points used for estimation of eucentric offset (default: 4)
alignLimit      = 0.5       # maximum shift [microns] allowed for record tracking between tilts, should reduce loss of target in case of low contrast (not applied for tracking TS); also the threshold to take a second tracking image when using trackTwice
minCounts       = 0         # minimum mean counts per second of record image (if set > 0, tilt series branch will be aborted if mean counts are not sufficient)
ignoreNegStart  = True      # ignore first shift on 2nd branch, which is usually very large on bad stages
realignToItem   = False     # Use SerialEM's RealignToItem routine instead of simple image realignment (was default in PACEtomo <=v1.9.1)
refFromPreview  = False     # Makes temporary reference from Preview image collected during previewAli for use with first Record image
noZeroRecAli    = False     # Skip alignment of first tilt image to reference 
autoStartTilt   = False     # Uses measured pretilt to set compensating startTilt      
tiltTargets     = 0         # Stage tilt at which targets were selected (if not 0, it will be automatically used as startTilt!)

# Target montage settings
tgtMontage      = False     # collect montage for each target using the shorter camera dimension (e.g. for square aperture montage tomography)
tgtMntSize      = 1         # size of montage pattern (1: 3x3, 2: 5x5, 3: 7x7, ...)
tgtMntOverlap   = 0.05      # montage tile overlap as fraction of shorter camera dimension
tgtMntXOffset   = 0         # max offset [microns] applied along tilt axis throughout tilt series (+tgtMntXOffset is reached at maxTilt, -tgtMntXOffset at minTilt)
tgtMntFocusCor  = False     # do focus compensation for tiles of montage
tgtTrackMnt     = False     # set to True if you also want the tracking target to be a montage

debug           = False     # Enables additional output for a few processes (e.g. cross-correlation for all image alignments)
breakpoints     = False     # Waits at every debug output for user to press B key.

########## END SETTINGS ########## 

versionPACE = "1.9.2"

import serialem as sem
import os
import copy
import time
import struct
from datetime import datetime
import glob
from functools import wraps
import numpy as np
from scipy import optimize
if sortByTilt: import mrcfile

versionCheck = sem.IsVersionAtLeast("40200", "20240814")
if not versionCheck and sem.IsVariableDefined("warningVersion") == 0:
    runScript = sem.YesNoBox("\n".join(["WARNING: You are using a version of SerialEM that does not support all PACEtomo features. It is recommended to update to the latest SerialEM beta version!", "", "Do you want to run PACEtomo regardless?"]))
    if not runScript:
        sem.Exit()
    else:
        sem.SetPersistentVar("warningVersion", "")

########### FUNCTIONS ###########

def checkFilling():
    global dewarFillTime
    filling = sem.AreDewarsFilling()
    timerStart = 0
    if filling >= 1:
        timerStart = time.time()
        log(datetime.now().strftime("%d.%m.%Y %H:%M:%S") + ": Dewars are filling...", color=4)
        if cryoARM:                                                                             # make sure both tanks are being filled on cryoARM
            sem.LongOperation("RS", "0", "RT", "0")
        if coldFEG:                                                                             # flash gun while dewars refill
            sem.LongOperation("FF", "0")
    while filling >= 1:
        log("Dewars are still filling...")
        sem.Delay(60, "s")
        filling = sem.AreDewarsFilling()
    if timerStart > 0:
        log(f"Dewars finished filling after {(time.time() - timerStart) / 60} minutes.")
        dewarFillTime = time.time() - timerStart

def checkColdFEG():
    if not cryoARM:                                                                             # Routine for Krios CFEG with Advanced scripting >4
        flashLow = 0
        flashHigh = sem.IsFEGFlashingAdvised(1)
        if flashHigh == 1:
            sem.NextFEGFlashHighTemp(1)
        else:
            flashLow = sem.IsFEGFlashingAdvised(0)
        if flashLow == 1 or flashHigh ==1:
            sem.LongOperation("FF", "0")
    else:
            sem.LongOperation("FF", str(flashInterval))

def checkSlit(vec, size, tilt, pn):                                                             # check ZLP in hole outside of pattern along tilt axis
    global lastSlitCheck
    log("Refining ZLP...", style=1)
    sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
    shift = vec * (size + 1)
    shift[1] *= np.cos(np.radians(tilt))
    sem.ImageShiftByMicrons(*shift)
    sem.RefineZLP()
    sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
    lastSlitCheck = sem.ReportClock()

def checkValves():
    if not int(sem.ReportColumnOrGunValve()):
        sem.SetColumndOrGunValve(1)

def retryOpen(max_attempts=5, delay=5):
    """Decorator to retry function on permission exception."""

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            attempts = 0
            while attempts < max_attempts:
                try:
                    return func(*args, **kwargs)
                except PermissionError as e:
                    attempts += 1
                    if attempts == max_attempts:
                        raise e
                    log(f"WARNING: File [{args[0]}] could not be opened. Trying again... [{attempts}]")
                    time.sleep(delay)
            return None
        return wrapper
    return decorator

def openOldFile(file_path, max_attempts=5, delay=5):
    # Allow for multiple attempts to open file to prevent crash when file is being copied during access
    attempts = 0
    while attempts < max_attempts:
        try:
            sem.NoMessageBoxOnError(1)
            sem.OpenOldFile(file_path)
            # Break on success
            break
        except Exception as e:
            attempts += 1
            if attempts == max_attempts:
                raise e
            log(f"WARNING: File [{file_path}] could not be opened. Trying again... [{attempts}]")
            time.sleep(delay)
        finally:
            sem.NoMessageBoxOnError(0)

@retryOpen()
def parseTargets(file_path):
    """Reads targets file."""

    with open(file_path) as f:                                                                         # open last tgts or tgts_run file
        targetFile = f.readlines()

    targets = []
    geoPoints = []
    savedRun = []
    branch = None
    resume = {"sec": 0, "pos": 0}
    for line in targetFile:
        col = line.strip(os.linesep).split(" ")
        if col[0] == "": continue
        if line.startswith("_set") and len(col) == 4:
            if col[1] in globals():
                log(f"WARNING: Read setting from tgts file and overwrite: {col[1]} = {col[3]}")
                globals()[col[1]] = float(col[3])
            else:
                log(f"WARNING: Attempted to overwrite {col[1]} but variable does not exist!")
        elif line.startswith("_bset") and len(col) == 4:
            if col[1] in globals():
                val = True if col[3].lower() in ["true", "yes", "y", "on"] else False
                log(f"WARNING: Read setting from tgts file and overwrite: {col[1]} = {val}")
                globals()[col[1]] = val
            else:
                log(f"WARNING: Attempted to overwrite {col[1]} but variable does not exist!")                
        elif line.startswith("_spos"):
            resume["sec"] = int(col[2].split(",")[0])
            resume["pos"] = int(col[2].split(",")[1])
        elif line.startswith("_tgt"):
            targets.append({})
            branch = None
        elif line.startswith("_pbr"):
            savedRun.append([{},{}])
            branch = 0
        elif line.startswith("_nbr"):
            branch = 1
        elif line.startswith("_geo"):
            geoPoints.append({})
            branch = "geo"
        else:
            if branch is None:
                targets[-1][col[0]] = col[2]
            elif branch == "geo":
                geoPoints[-1][col[0]] = float(col[2])
            else:
                savedRun[-1][branch][col[0]] = col[2]
    if savedRun == []: savedRun = False
    return targets, savedRun, resume, geoPoints

@retryOpen()
def updateTargets(fileName, targets, position=[], sec=0, pos=0):
    output = ""
    if sec > 0 or pos > 0:
        output += "_set startTilt = " + str(startTilt) + "\n"
        output += "_set minTilt = " + str(minTilt) + "\n"
        output += "_set maxTilt = " + str(maxTilt) + "\n"
        output += "_set step = " + str(step) + "\n"
        output += "_set pretilt = " + str(pretilt) + "\n"
        output += "_set rotation = " + str(rotation) + "\n"
        output += "_spos = " + str(sec) + "," + str(pos) + "\n" * 2
    for pos in range(len(targets)):
        output += "_tgt = " + str(pos + 1).zfill(3) + "\n"
        for key in targets[pos].keys():
            output += key + " = " + targets[pos][key] + "\n"
        if position != []:
            output += "_pbr" + "\n"
            for key in position[pos][1].keys():
                output += key + " = " + str(position[pos][1][key]).strip("[]").replace("np.float64(", "").replace(")", "").replace(" ","") + "\n"   # temp fix for numpy>=2.0 changing str output for list of numbers to include np.float64()
            output += "_nbr" + "\n"
            for key in position[pos][2].keys():
                output += key + " = " + str(position[pos][2][key]).strip("[]").replace("np.float64(", "").replace(")", "").replace(" ","") + "\n"   # temp fix for numpy>=2.0 changing str output for list of numbers to include np.float64()     
        output += "\n"
    with open(fileName, "w") as f:
        f.write(output)

def geoPlane(x, a, b):
    return a * x[0] + b * x[1]

def geoPara(x, a, b, c, d, e):
    return a * x[0] + b * x[1] + c * (x[0]**2) + d * (x[1]**2) + e * x[0] * x[1]

def parseMdoc(mdocFile):
    with open(mdocFile, "r") as f:
        content = f.readlines()
    header = []
    items = []
    newItem = {}
    index = 0

    for line in content:
        col = line.strip().split()
        if len(col) == 0: 
            if "ZValue" in newItem.keys():
                items.append(newItem)
                newItem = {}
                continue
            else:
                continue
        if line.startswith("[ZValue"):
            index += 1
            newItem = {"index": index, "ZValue": col[2].strip("]")}
        elif "ZValue" in newItem.keys():
            newItem[col[0]] = [val for val in col[2:]]
        else:
            header.append(line.strip())
    if "ZValue" in newItem.keys(): # append last item
        items.append(newItem)
    return header, items

def writeMdoc(header, items, filename):
    content = ""
    for line in header:
        if line.startswith("[T"):
            content += os.linesep
        content += line + os.linesep
    content += os.linesep
    for i, item in enumerate(items):
        content += "[ZValue = " + str(i) + "]" + os.linesep
        item.pop("ZValue")
        item.pop("index")
        for key, attr in item.items():
            content += key + " = " + (" ".join(attr) if key != "DateTime" else "  ".join(attr)) # Date and time is separated by double space in SerialEM
            content += os.linesep
        if i < len(items) - 1: # Don't add additional linebreak at end of file
            content += os.linesep
    with open(filename, "w", newline="") as f:
        f.write(content)
    return

def getExtendedHeader(filename):
    with open(filename, "rb") as mrc_file:
        mrc = mrc_file.read(4096)

    # MRC header format SERI (SERI is not covered by mrcfile)
    format = ""
    for char in struct.iter_unpack("s", mrc[104:108]):
        format += char[0].decode("utf-8")
    log(f"DEBUG: MRC format: {format}")

    if format == "SERI":
        # Get number of sections
        section_number = struct.unpack("i", mrc[8: 12])[0]
        log(f"DEBUG: Number of sections: {section_number}")

        # Get bytes per section
        bytes_per_section = struct.unpack("h", mrc[128: 130])[0]
        log(f"DEBUG: Bytes per section: {bytes_per_section}")

        # Bitflags
        bitflags = struct.unpack("h", mrc[130: 132])[0]
        log(f"DEBUG: Bitflags: {bitflags}")
        """
        https://bio3d.colorado.edu/imod/doc/mrc_format.txt
        1 = Tilt angle in degrees * 100  (2 bytes)
        2 = X, Y, Z piece coordinates for montage (6 bytes)
        4 = X, Y stage position in microns * 25   (4 bytes)
        8 = Magnification / 100 (2 bytes)
        16 = Intensity * 25000  (2 bytes)
        32 = Exposure dose in e-/A2, a float in 4 bytes
        128, 512: Reserved for 4-byte items
        64, 256, 1024: Reserved for 2-byte items
        If the number of bytes implied by these flags does
        not add up to the value in nint, then nint and nreal
        are interpreted as ints and reals per section
        """

        section_data = []
        for i in range(1024, 1024 + bytes_per_section * section_number, bytes_per_section):     # extended header starts at byte 1024
            section_data.append(mrc[i:i + bytes_per_section])

        return section_data
    
    return None

def writeExtendedHeader(filename, section_data):
    with open(filename, "r+b") as mrc_file:
        mrc_file.seek(1024)
        mrc_file.write(b"".join(section_data))

@retryOpen()
def sortTS(ts_name):
    log(f"Sorting {ts_name} by tilt angle...")
    if os.path.exists(os.path.join(curDir, ts_name + ".mdoc")):
        # Read mdoc file
        header, tilts = parseMdoc(os.path.join(curDir, ts_name + ".mdoc"))
        # Make list of tilt angles
        tiltAngles = [float(tilt["TiltAngle"][0]) for tilt in tilts]
        log(f"DEBUG: Tilts before sorting: {tiltAngles}")
        # Get extended header data
        ext_header_sections = getExtendedHeader(os.path.join(curDir, ts_name))
        log(f"DEBUG: MRC header sections: {len(ext_header_sections)}")
        # Open mrc file
        with mrcfile.open(os.path.join(curDir, ts_name), "r+") as mrc:
            stack = mrc.data
            log(f"DEBUG: MRC data dims: {np.array(stack).shape}")
            # Sort tilts and stack according to tilt angle
            zippedTilts = sorted(zip(tiltAngles, tilts, stack, ext_header_sections), key=lambda x: x[0])
            tiltAngles, tilts, stack, ext_header_sections = zip(*zippedTilts)
            log(f"DEBUG: Tilts after sorting: {tiltAngles}")
            stack = np.array(stack)
            # Save new stack in same file
            mrc.set_data(stack)
        # Update extended header
        writeExtendedHeader(os.path.join(curDir, ts_name), ext_header_sections)
        # Rename original mdoc
        os.rename(os.path.join(curDir, ts_name + ".mdoc"), os.path.join(curDir, os.path.splitext(ts_name)[0] + "_unsorted.mrc.mdoc"))
        # Save sorted mdoc
        writeMdoc(header, tilts, os.path.join(curDir, ts_name + ".mdoc"))
        log(f"NOTE: {ts_name} was sorted by tilt angle!")
    else:
        log(f"WARNING: {ts_name}.mdoc file could not be found! Tilt series stack was not sorted!")

def bin2d(img, factor):
    # Add third dimension if img is not stack
    if img.ndim == 2:
        img = np.expand_dims(img, 0)
    # Only bin in 2D and keep stack size
    factors = (1, factor, factor)
    # Calculate the new dimensions after cropping to allow even binning
    new_shape = tuple((dim // factor) * factor for dim, factor in zip(img.shape, factors))
    # Center crop the array to the new dimensions
    slices = tuple(slice((dim - new_dim) // 2, (dim + new_dim) // 2) for dim, new_dim in zip(img.shape, new_shape))
    cropped_img = img[slices]
    # Determine the new shape for reshaping
    reshaped_shape = np.array([(dim // factor, factor) for dim, factor in zip(cropped_img.shape, factors)]).reshape(-1)
    # Reshape the array
    reshaped_img = cropped_img.reshape(reshaped_shape)
    # Calculate the mean along the new axes
    for i in range(-1, -cropped_img.ndim-1, -1):
        reshaped_img = reshaped_img.mean(axis=i)
    # Remove added dimension
    if reshaped_img.shape[0] == 1:
        reshaped_img = reshaped_img[0, :, :]
    return reshaped_img

@retryOpen()
def binStack(ts_name, factor):
    factor = int(factor)
    log(f"Binning {ts_name} by {factor}...")
    if checkFrames(ts_name):
        # Check if tilt stack file exists
        if os.path.exists(os.path.join(curDir, ts_name)):
            # Open mrc file
            with mrcfile.open(os.path.join(curDir, ts_name), "r+") as mrc:
                stack = mrc.data
                voxel_size = mrc.voxel_size.x
                # Bin stack by factor
                stack = bin2d(stack, factor)
                # Save new stack in same file
                mrc.set_data(stack.astype(np.int16))
                # Update header pixel size
                mrc.voxel_size = voxel_size * factor
            log(f"NOTE: {ts_name} was binned by {factor}. Please use saved frames to regenerate the unbinned tilt series.")
        else:
            log(f"WARNING: {ts_name} file could not be found!")

def checkFrames(ts_name):
    # Check if frame saving is available by checking if warning was accepted by user at start of script
    if sem.IsVariableDefined("warningFramePath") == 0:
        # Make sure frames were saved
        frame_file, frame_dir, frame_name = sem.ReportLastFrameFile()
        if checkFrames and len(glob.glob(os.path.join(frame_dir, os.path.splitext(ts_name)[0] + "*"))) > 0:
            return True

    log(f"WARNING: Frames for {ts_name} could not be found. Keeping tilt stack unprocessed.")
    return False

def alignTo(buffer, debug=False):
    sem.AlignTo(buffer, 0, 0, 0, int(debug))
    if debug:
        try:
            sem.AddBufToStackWindow("A", 0, 0, 0, 0, "CC") #M #S [#B] [#O] [title]
        except AttributeError:
            # Show CC briefly, then switch back to aligned buffer for buffer shift
            sem.Delay(1, "s")
        sem.Copy("B", "A")
        sem.AlignTo(buffer)

def realignTo(nav_id=None, target=None):
    if target is not None and not realignToItem:
        # Move stage to target position
        sem.MoveStageTo(float(target["stageX"]), float(target["stageY"]))
        if "viewfile" in target.keys():
            sem.ReadOtherFile(0, "O", target["viewfile"]) # reads view file for first AlignTo instead
            sem.V()
            alignTo("O", debug)
            ASX, ASY = sem.ReportAlignShift()[4:6]
            log(f"Alignment (View) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")    
        if "tgtfile" in target.keys():                
            sem.ReadOtherFile(0, "O", target["tgtfile"]) # reads tgt file for first AlignTo instead
            sem.L()
            alignTo("O", debug)
            AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
            log(f"Alignment (Prev) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
        elif "viewfile" in target.keys():
            # Use align between mags to align preview image to view image
            # If View image was already aligned, take new centered View image at startTilt and use as reference instead
            sem.V()
            sem.Copy("A", "O")
            # Check defocus offset
            sem.GoToLowDoseArea("R") # Switch to R before applying defocus offset to not mess with potential mP/nP offsets between View and Rec
            defocus_offset = max(-10, sem.ReportLDDefocusOffset("V"))
            if defocus_offset != 0:
                sem.ChangeFocus(defocus_offset) # Higher defocus for better correlation, but max at 10 to avoid major distortions
            sem.L()
            sem.AlignBetweenMags("O", -1, -1, -1)
            AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
            if defocus_offset != 0:
                sem.ChangeFocus(-defocus_offset) # Reset focus
            log(f"Alignment (Pv2V) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
        else:
            log(f"WARNING: No target file or view file found for realignment!")
    elif nav_id is not None:
        sem.RealignToOtherItem(nav_id, 1)
    else:
        log(f"WARNING: No target provided for realignment!")

def log(text, color=0, style=0):
    if text.startswith("DEBUG:") and not debug:
        return
    if text.startswith("NOTE:"):
        color = 4
    elif text.startswith("WARNING:"):
        color = 5
    elif text.startswith("ERROR:"):
        color = 2
        style = 1 
    elif text.startswith("DEBUG:"):
        color = 1
        if breakpoints:
            breakpoint()

    if sem.IsVersionAtLeast("40200", "20240205"):
        sem.SetNextLogOutputStyle(style, color)
    sem.EchoBreakLines(text)

def breakpoint():
    """Breakpoint for debugging in SerialEM."""

    while not sem.KeyBreak():
        sem.Delay(0.1, "s")
    for i in range(5):
        if sem.KeyBreak("d"):
            dumpVars()
            break
        sem.Delay(0.1, "s")

def Tilt(tilt):
    def calcSSChange(x, z0):                                                                    # x = array(tilt, n0) => needs to be one array for optimize.curve_fit()
        return x[1] * (np.cos(np.radians(x[0])) - np.cos(np.radians(x[0] - increment))) - z0 * (np.sin(np.radians(x[0])) - np.sin(np.radians(x[0] - increment)))

    def calcFocusChange(x, z0):                                                                 # x = array(tilt, n0) => needs to be one array for optimize.curve_fit()
        return z0 * (np.cos(np.radians(x[0])) - np.cos(np.radians(x[0] - increment))) + x[1] * (np.sin(np.radians(x[0])) - np.sin(np.radians(x[0] - increment)))

    def setTrack():
        global trackMag, origMag
        if trackDefocus < maxDefocus:
            sem.SetDefocus(position[0][pn]["focus"] + trackDefocus - targetDefocus)
        if trackExpTime > 0:
            if tilt == startTilt:
                sem.SetExposure("R", max(trackExpTime, zeroExpTime))
            else:
                sem.SetExposure("R", trackExpTime)
        if trackMag > 0:
            if tilt == startTilt:
                origMag, *_ = sem.ReportMag()
                sem.UpdateLowDoseParams("R")
            attempt = 0
            while sem.ReportMag()[0] == origMag:                                                # has to be checked, because Rec is sometimes not updated (JEOL)
                if attempt >= 10:
                    log("WARNING: Magnification could not be changed. Continuing with the same magnification for all tilt series.")
                    trackMag = 0
                    break
                sem.SetMag(trackMag)
                sem.GoToLowDoseArea("R")
                attempt += 1
            sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
            if not recover:
                sem.ImageShiftByMicrons(0, SSchange)    

    def resetTrack():
        if trackExpTime > 0:
            sem.RestoreCameraSet("R")
        if trackMag > 0:
            while sem.ReportMag()[0] != origMag:                                                # has to be checked, because Rec is sometimes not updated (JEOL)
                sem.SetMag(origMag)
                sem.GoToLowDoseArea("R")

    global recover

    # Tilt if within tilt range, skip branch if not
    if -tiltLimit <= tilt <= tiltLimit:
        if tilt != startTilt:
            sem.TiltTo(tilt)
        skip_branch = False
    else:
        log(f"WARNING: Tilt angle [{tilt} degrees] could not be reached. This branch of the tilt series will be aborted.")
        skip_branch = True

    if tilt < startTilt:
        increment = -step
        if tilt - step >= -tiltLimit:
            sem.TiltBy(-step)
            sem.TiltTo(tilt)
        pn = 2
    else:
        if slowTilt and startTilt < tilt <= tiltLimit:                                                       # on bad stages, better to do backlash as well to enhance accuracy
            sem.TiltBy(-step)
            sem.TiltTo(tilt)
        increment = step
        pn = 1

    # After branch was determined, set branch to be skipped and return
    if skip_branch:
        for pos in range(len(position)):
            position[pos][pn]["skip"] = True
        return

    sem.Delay(delayTilt, "s")
    realTilt = float(sem.ReportTiltAngle())

    if zeroExpTime > 0 and tilt == startTilt:
        sem.SetExposure("R", zeroExpTime)

    if recover:
        # preview align to last tracking TS
        openOldFile(targets[0]["tsfile"])
        sem.ReadFile(position[0][pn]["sec"], "O")                                               # read last image of position for AlignTo
        sem.SetDefocus(position[0][pn]["focus"])
        sem.SetImageShift(position[0][pn]["ISXset"], position[0][pn]["ISYset"])
        SSchange = 0                                                                            # needs to be defined for setTrack
        setTrack()
        if checkDewar: checkFilling()
        sem.L()
        alignTo("O", debug)
        bufISX, bufISY = sem.ReportISforBufferShift()
        sem.ImageShiftByUnits(position[0][pn]["ISXali"], position[0][pn]["ISYali"])             # remove accumulated buffer shifts to calculate alignment to initial startTilt image
        position[0][pn]["ISXset"], position[0][pn]["ISYset"], *_ = sem.ReportImageShift()
        for i in range(1, len(position)):
            position[i][pn]["ISXset"] += bufISX + position[0][pn]["ISXali"]                     # apply accumulated (stage dependent) buffer shifts of tracking TS to all targets
            position[i][pn]["ISYset"] += bufISY + position[0][pn]["ISYali"]
        resetTrack()
        sem.CloseFile()

        posStart = posResumed
    else:
        posStart = 0

    for pos in range(posStart, len(position)):
        log(f"\nTarget {pos + 1} / {len(position)}:", style=1)
        sem.SetStatusLine(2, "Target: " + str(pos + 1) + " / " + str(len(position)))
        if pos != 0 and position[pos][pn]["skip"]: 
            log(f"[{pos + 1}] was skipped on this branch.")
            continue
        if tilt != startTilt:
            openOldFile(targets[pos]["tsfile"])
            sem.ReadFile(position[pos][pn]["sec"], "O")                                         # read last image of position for AlignTo
        else:
            if os.path.exists(os.path.join(curDir, targets[pos]["tsfile"])):
                # Close all files incase file to be renamed is currently open
                while sem.ReportFileNumber() > 0:
                    sem.CloseFile()
                os.replace(os.path.join(curDir, targets[pos]["tsfile"]), os.path.join(curDir, targets[pos]["tsfile"]) + "~")
                log("WARNING: Tilt series file already exists. Existing file was renamed.")
            sem.OpenNewFile(targets[pos]["tsfile"])
            if not tgtPattern and "tgtfile" in targets[pos].keys():
                if refFromPreview:
                    temp_ref = os.path.splitext(targets[pos]["tgtfile"])[0] + "_tempref.mrc"
                    if os.path.exists(temp_ref) and "prevASX" in targets[pos].keys():
                        sem.ReadOtherFile(0, "O", temp_ref)                                     # reads temp reference from previewAli instead

                sem.ReadOtherFile(0, "O", targets[pos]["tgtfile"])                              # reads tgt file for first AlignTo instead

        sem.AreaForCumulRecordDose(pos + 1)                                                     # set area to accumulate record dose (counting from 1)

### Calculate and apply predicted shifts
        SSchange = 0                                                                            # only apply changes if not startTilt
        focuschange = 0
        if tilt != startTilt:
            SSchange = calcSSChange([realTilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])
            focuschange = calcFocusChange([realTilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])

        SSYprev = position[pos][pn]["SSY"]
        SSYpred = position[pos][pn]["SSY"] + SSchange

        focuscorrection = focusSlope * (tilt - startTilt)
        position[pos][pn]["focus"] += focuscorrection
        position[pos][pn]["focus"] -= focuschange

        sem.SetDefocus(position[pos][pn]["focus"])
        if zeroDefocus != 0 and tilt == startTilt:
            sem.ChangeFocus(zeroDefocus - maxDefocus)

        sem.SetImageShift(position[pos][pn]["ISXset"], position[pos][pn]["ISYset"])
        sem.ImageShiftByMicrons(0, SSchange)

        # Apply sliding offset along tilt axis throughout tilt series when collecting montage tilt series
        if tgtMontage and (tgtTrackMnt or pos != 0) and tgtMntXOffset > 0:
            mont_offset = np.array([tgtMntXOffset * 2 * (tilt - startTilt) / (maxTilt - minTilt), 0])
            sem.ImageShiftByMicrons(mont_offset)
        else:
            mont_offset = None

### Autofocus (optional) and tracking TS settings
        if pos == 0:
            if addAF and (tilt - startTilt) % (groupSize * increment) == step and abs(tilt - startTilt) > step:
                sem.G(-1)
                defocus, *_ = sem.ReportAutoFocus()
                focuserror = float(defocus) - targetDefocus
                for i in range(0, len(position)):
                    position[i][pn]["focus"] -= focuserror
                sem.SetDefocus(position[pos][pn]["focus"])

            setTrack()

### Record
        if checkDewar: checkFilling()
        checkValves()
        sem.SetFrameBaseName(0, 1, 0, os.path.splitext(targets[pos]["tsfile"])[0])              # change frame name in accordance with tilt series
        if beamTiltComp: 
            sem.AdjustBeamTiltforIS()
        sem.Delay(delayIS, "s")
        sem.R()
        sem.S()

        bufISXpre = 0                                                                           # only non 0 if two tracking images are taken
        bufISYpre = 0
        if tilt != startTilt or (not tgtPattern and "tgtfile" in targets[pos].keys() and not noZeroRecAli): # align to previous image if it exists 
            if pos != 0: 
                sem.LimitNextAutoAlign(alignLimit)                                              # gives maximum distance for AlignTo to avoid runaway tracking
            alignTo("O", debug)
            if trackTwice and pos == 0:                                                         # track twice if alignLimit for tracking area is surpassed
                ASX, ASY = sem.ReportAlignShift()[4:6]
                if abs(ASX) > alignLimit * 1000 or abs(ASY) > alignLimit * 1000:
                    bufISXpre, bufISYpre = sem.ReportISforBufferShift()                         # have to be added only to ISset but not ISali (since ali only considers the IS chain of ali images)
                    sem.R()
                    sem.S()
                    alignTo("O", debug)

        bufISX, bufISY = sem.ReportISforBufferShift()

        # Subtract montage offset if given
        if mont_offset is not None:
            bufISX, bufISY = np.array([bufISX, bufISY]) - ss2isMatrix @ mont_offset

        sem.ImageShiftByUnits(position[pos][pn]["ISXali"], position[pos][pn]["ISYali"])         # remove accumulated buffer shifts to calculate alignment to initial startTilt image

        if beamTiltComp: 
            sem.RestoreBeamTilt()

        position[pos][pn]["ISXset"], position[pos][pn]["ISYset"], *_ = sem.ReportImageShift()
        position[pos][pn]["SSX"], position[pos][pn]["SSY"] = sem.ReportSpecimenShift()

        # Collect surrounding tiles for montage tilt series
        if tgtMontage and (tgtTrackMnt or pos != 0):
            sem.ImageShiftByUnits(-bufISX - position[pos][pn]["ISXali"], -bufISY - position[pos][pn]["ISYali"]) # reset shifts to already taken center image
            for i in range(-tgtMntSize, tgtMntSize + 1):
                for j in range(-tgtMntSize, tgtMntSize + 1):
                    if i == j == 0: continue
                    if tilt != startTilt:
                        openOldFile(os.path.splitext(targets[pos]["tsfile"])[0] + "_" + str(i) + "_" + str(j) + ".mrc")
                    else:
                        sem.OpenNewFile(os.path.splitext(targets[pos]["tsfile"])[0] + "_" + str(i) + "_" + str(j) + ".mrc")

                    montX, montY = (i - i * tgtMntOverlap) * min([camX, camY]), (j - j * tgtMntOverlap) * min([camX, camY])
                    pixelShiftFromCenter = f"{montX} {montY}"

                    # Apply sliding offset along tilt axis throughout tilt series
                    if tgtMntXOffset > 0:
                        # Adjust x shift by tgtMntXOffset [microns at max tilt] by fraction of tilt series along SSX
                        montX, montY = np.array([montX, montY]) + ss2cMatrix @ np.array([tgtMntXOffset * 2 * (tilt - startTilt) / (maxTilt - minTilt), 0])

                    sem.ImageShiftByPixels(montX, montY)
                    if tgtMntFocusCor:
                        montSSX, montSSY = c2ssMatrix @ np.array([montX, montY])

                        # With sample geometry (needs to be tested)
                        correctedFocus = position[pos][pn]["focus"] - np.cos(np.radians(realTilt)) * np.tan(np.radians(pretilt)) * (np.cos(np.radians(rotation)) / np.cos(np.radians(realTilt)) * montSSY - np.sin(np.radians(rotation)) * montSSX) - np.tan(np.radians(realTilt)) * montSSY 
                        # Without sample geometry
                        #correctedFocus = position[pos][pn]["focus"] - np.tan(np.radians(realTilt)) * montSSY

                        sem.SetDefocus(correctedFocus)
                    if beamTiltComp: 
                        sem.AdjustBeamTiltforIS()
                    sem.Delay(delayIS, "s")
                    sem.R()
                    sem.S()

                    mont_SSX, mont_SSY = sem.ReportSpecimenShift()

                    sem.ImageShiftByPixels(-montX, -montY)
                    if beamTiltComp: 
                        sem.RestoreBeamTilt()

                    # Add shift to all montage tilt series mdoc files for auto stitching
                    if extendedMdoc:
                        sem.AddToAutodoc("PixelShiftFromCenter", pixelShiftFromCenter)
                        sem.WriteAutodoc()

                    sem.CloseFile()

        position[pos][pn]["focus"] -= focuscorrection                                           # remove correction or it accumulates

        dose = sem.ImageConditions("A")[0]
        if dose > 0:
            sem.AccumulateRecordDose(dose)
            position[pos][1]["dose"] += dose
            position[pos][2]["dose"] += dose

        if pos == 0:                                                                            # apply measured shifts of first/tracking position to other positions
            for i in range(1, len(position)):
                position[i][pn]["ISXset"] += bufISX + bufISXpre + position[pos][pn]["ISXali"]   # apply accumulated (stage dependent) buffer shifts of tracking TS to all targets
                position[i][pn]["ISYset"] += bufISY + bufISYpre + position[pos][pn]["ISYali"]
                if tilt == startTilt:                                                           # also save shifts from startTilt image for second branch since it will alignTo the startTilt image
                    position[i][2]["ISXset"] += bufISX + bufISXpre
                    position[i][2]["ISYset"] += bufISY + bufISYpre
            if tilt == startTilt:                                                               # do not forget about 0 position
                position[0][2]["ISXset"] += bufISX + bufISXpre
                position[0][2]["ISYset"] += bufISY + bufISYpre

            resetTrack()

        position[pos][pn]["ISXali"] += bufISX
        position[pos][pn]["ISYali"] += bufISY
        if tilt == startTilt:                                                                   # save alignment of first tilt to tgt file for the second branch
            position[pos][2]["ISXali"] += bufISX
            position[pos][2]["ISYali"] += bufISY

        aErrX, aErrY = is2ssMatrix @ np.array([position[pos][pn]["ISXali"], position[pos][pn]["ISYali"]])

        log(f"[{pos + 1}] Prediction: y = {round(SSYpred, 3)} microns | z = {round(position[pos][pn]['focus'], 3)} microns | z0 = {round(position[pos][pn]['z0'], 3)} microns")
        log(f"[{pos + 1}] Reality: y = {round(position[pos][pn]['SSY'], 3)} microns")
        log(f"[{pos + 1}] Focus change: {round(focuschange, 3)} microns | Focus correction: {round(focuscorrection, 3)} microns")
        log(f"[{pos + 1}] Alignment error: x = {round(aErrX * 1000)} nm | y = {round(aErrY * 1000)} nm")        

### Calculate new z0

        ddy = position[pos][pn]["SSY"] - SSYprev
        if (tilt == startTilt or
                (ignoreNegStart and pn == 2 and len(position[pos][pn]["shifts"]) == 0) or
                recover or
                (resumePN == 1 and tilt == resumePlus + step and pos < posResumed) or
                (resumePN == 1 and tilt == resumeMinus - step) or
                (resumePN == 2 and tilt == resumeMinus - step and pos < posResumed) or
                (resumePN == 2 and tilt == resumePlus + step)):        
                # ignore shift if first image or first shift of second branch or first image after resuming run (all possible conditions)
            ddy = calcSSChange([realTilt, position[pos][pn]["n0"]], position[pos][pn]["z0"])

        position[pos][pn]["shifts"].append(ddy)
        position[pos][pn]["angles"].append(realTilt)

        if len(position[pos][pn]["shifts"]) > dataPoints:
            position[pos][pn]["shifts"].pop(0)
            position[pos][pn]["angles"].pop(0)

        position[pos][pn]["z0"], cov = optimize.curve_fit(calcSSChange, np.vstack((position[pos][pn]["angles"], [position[pos][pn]["n0"] for i in range(0, len(position[pos][pn]["angles"]))])), position[pos][pn]["shifts"], p0=(position[pos][pn]["z0"]))
        position[pos][pn]["z0"] = position[pos][pn]["z0"][0]

        if doCtfFind:
            try:
                sem.NoMessageBoxOnError(1)
                cfind = sem.CtfFind("A", (min(maxDefocus, trackDefocus) - 2), min(-0.2, minDefocus + 2))
                log(f"[{pos + 1}] CtfFind: {round(cfind[0], 3)} microns ({round(cfind[-1], 2)} A)")
            except:
                log(f"WARNING: CtfFind crashed on {targets[pos]['tsfile']} section {int(sem.ReportFileZsize()) - 1}. Trying to continue...")
            finally:
                sem.NoMessageBoxOnError(0)

        if doCtfPlotter:
            try:
                sem.NoMessageBoxOnError(1)
                cplot = sem.Ctfplotter("A", (min(maxDefocus, trackDefocus) - 2), min(-0.2, minDefocus + 2), 1, 0, pretilt)
                log(f"[{pos + 1}] Ctfplotter: {round(cplot[0], 3)} microns")
            except:
                log(f"WARNING: Ctfplotter crashed on {targets[pos]['tsfile']} section {int(sem.ReportFileZsize()) - 1}. Trying to continue...")
            finally:
                sem.NoMessageBoxOnError(0)

        if refineGeo and tilt == startTilt:
            if doCtfPlotter:
                geo[0].append(position[pos][pn]["SSX"])
                geo[1].append(position[pos][pn]["SSY"])
                geo[2].append(cplot[0])
            elif doCtfFind and len(cfind) > 5:
                if cfind[5] < fitLimit:                                                         # save vectors for refineGeo only if CTF fit has reasonable resolution
                    geo[0].append(position[pos][pn]["SSX"])
                    geo[1].append(position[pos][pn]["SSY"])
                    geo[2].append(cfind[0])

        position[pos][pn]["sec"] = int(sem.ReportFileZsize()) - 1                               # save section number for next alignment

        # progress = collected images * (positions - skipped positions) + current position - skipped positions scaled assuming homogeneous distribution of skipped positions
        progress = position[pos][pn]["sec"] * (len(position) - skippedTgts) + pos - skippedTgts * pos / len(position) + 1
        percent = round(100 * (progress / maxProgress), 1)
        bar = '#' * int(percent / 2) + '_' * (50 - int(percent / 2))
        if percent - resumePercent > 0:
            remTime = int((sem.ReportClock() - startTime - dewarFillTime) / (percent - resumePercent) * (100 - percent) / 60)
        else:
            remTime = "?"
        log(f"Progress: |{bar}| {percent} % ({remTime} min remaining)")

        if extendedMdoc:
            sem.AddToAutodoc("SpecimenShift", str(position[pos][pn]["SSX"]) + " " + str(position[pos][pn]["SSY"]))
            sem.AddToAutodoc("EucentricOffset", str(position[pos][pn]["z0"]))
            if tgtMontage:
                sem.AddToAutodoc("PixelShiftFromCenter", "0 0")
            if doCtfFind:
                sem.AddToAutodoc("CtfFind", str(cfind[0]))
            if doCtfPlotter:
                sem.AddToAutodoc("Ctfplotter", str(cplot[0]))
            sem.WriteAutodoc()

        sem.CloseFile()

### Abort conditions
        if np.linalg.norm(np.array([position[pos][pn]["SSX"], position[pos][pn]["SSY"]], dtype=float)) > imageShiftLimit - alignLimit:
            position[pos][pn]["skip"] = True
            log(f"WARNING: Target [{pos + 1}] is approaching the image shift limit. This branch will be aborted.")

        if minCounts > 0:
            meanCounts = sem.ReportMeanCounts()
            expTime, *_ = sem.ReportExposure("R")
            if meanCounts / expTime < minCounts:
                position[pos][pn]["skip"] = True
                log(f"WARNING: Target [{pos + 1}] was too dark. This branch will be aborted.")

        if tilt >= maxTilt or tilt <= minTilt:
            position[pos][pn]["skip"] = True
            if maxTilt - startTilt != abs(minTilt - startTilt):
                log(f"WARNING: Target [{pos + 1}] has reached the final tilt angle. This branch will be aborted.")            

        updateTargets(runFileName, targets, position, position[pos][pn]["sec"], pos)    

### Refine energy filter slit if appropriate
    if tgtPattern and slitInterval > 0 and (lastSlitCheck - sem.ReportClock() / 60) > slitInterval:
        checkSlit(np.array([vecB0, vecB1]), size, realTilt, pn)

    if zeroExpTime > 0 and tilt == startTilt:
        sem.RestoreCameraSet("R")

    if recover:
        recover = False    

def dumpVars(filename):
    output = "# PACEtomo settings from " + datetime.now().strftime("%d.%m.%Y %H:%M:%S") + "\n"
    save = False
    for var in globals():                                                                       # globals() is ordered by creation, start and end points might have to be adjusted if script changes
        if var == "sem":                                                                        # first var after settings vars
            break
        if save:
            output += var + " = " + str(globals()[var]) + "\n"
        if var == "SEMflush":                                                                   # last var before settings vars
            save = True
    with open(filename + "_settings.txt", "w") as f:
        f.write(output)

######## END FUNCTIONS ########

# Adjust user settings
sem.SetProperty("ImageShiftLimit", imageShiftLimit)
tiltLimit = sem.ReportProperty("MaximumTiltAngle")

sem.SetUserSetting("DriftProtection", 1)
sem.SetUserSetting("ShiftToTiltAxis", 1)
sem.SetNewFileType(0)                                                                           # set file type to mrc in case user changed default file type
sem.SetFrameBaseName(0, 1, 0, "PACEtomo_setup")                                                 # change frame name at start to avoid overwriting in case sets other than Record save frames

# Warnings
log(f"DEBUG: Tilt limit is: {tiltLimit}")
if (maxTilt > tiltLimit or minTilt < -tiltLimit) and sem.IsVariableDefined("warningTiltAngle") == 0:
    sem.Pause("WARNING: Tilt angles go beyond +/- 70 degrees. Most stage limitations do not allow for symmetrical tilt series with these values!")
    sem.SetPersistentVar("warningTiltAngle", "")

if int(sem.ReportAxisPosition("F")[0]) != 0 and sem.IsVariableDefined("warningFocusArea") == 0:
    sem.Pause("WARNING: Position of Focus area is not 0! Please set it to 0 to autofocus on the tracking target!")
    sem.SetPersistentVar("warningFocusArea", "")

tiltAxisOffset = sem.ReportTiltAxisOffset()[0]
if float(tiltAxisOffset) == 0 and sem.IsVariableDefined("warningTAOffset") == 0:
    sem.Pause("WARNING: No tilt axis offset was set! Please run the PACEtomo_measureOffset script to determine appropiate tilt axis offset.")
    sem.SetPersistentVar("warningTAOffset", "")

# Saving SerialEM setup
sem.SaveSettings()
sem.SaveNavigator()

sem.SuppressReports()
if beamTiltComp:                                                                                # check if there is a calibration saved, throws error if not
    sem.ReportComaVsISmatrix()
if tgtPattern:                                                                                  # initialize in case tgts file contains values
    vecA0 = vecA1 = vecB0 = vecB1 = size = None

### Find target file
sem.ReportNavItem()
navID = int(sem.GetVariable("navIndex"))
navNote = sem.GetVariable("navNote")
fileStem, fileExt = os.path.splitext(navNote)
curDir = sem.ReportDirectory()

if fileStem != "" and fileExt == ".txt":
    tf = sorted(glob.glob(os.path.join(curDir, fileStem + ".txt")))                             # find  tgts file
    tfr = sorted(glob.glob(os.path.join(curDir, fileStem + "_run??.txt")))                      # find run files but not copied tgts file
    tf.extend(tfr)                                                                              # only add run files to list of considered files
    while tf == []:
        searchInput = sem.YesNoBox("\n".join(["Target file not found! Please choose the directory containing the target file!", "WARNING: All future target files will be searched here!"]))
        if searchInput == 0:
            sem.Exit()
        sem.UserSetDirectory("Please choose the directory containing the target file!")
        curDir = sem.ReportDirectory()
        tf = sorted(glob.glob(os.path.join(curDir, fileStem + ".txt")))                         # find  tgts file
        tfr = sorted(glob.glob(os.path.join(curDir, fileStem + "_run??.txt")))                  # find run files but not copied tgts file
        tf.extend(tfr)                                                                          # only add run files to list of considered files
else:
    sem.OKBox("The navigator item note does not contain a target file. Make sure to setup PACEtomo targets using the selectTargets script and select the Navigator item marked to be acquired!")
    sem.Exit()

# Check if frame folder is set reasonably
if sem.ReportCameraProperty(0, "K2Type") > 0:
    framePath = sem.ReportFrameSavingPath()
    framePar = os.path.abspath(os.path.join(framePath, os.pardir))
    curPar = os.path.abspath(os.path.join(curDir, os.pardir))
    if (framePath == "NONE" or (framePath not in curDir and curDir not in framePath and framePar not in curDir and curPar not in framePath)) and sem.IsVariableDefined("warningFramePath") == 0:
        sem.Pause("WARNING: Current frame path (" + framePath + ") does not seem plausible or camera does not save frames.")
        sem.SetPersistentVar("warningFramePath", "")
else:
    log("WARNING: Camera frame path could not be obtained for your camera.")

sem.SaveLogOpenNew(navNote.split("_tgts")[0])

log(f"PACEtomo Version {versionPACE}", color=5, style=1)
sem.ProgramTimeStamps()

# Open last tgts or tgts_run file and read contents
targets, savedRun, resume, geoPoints = parseTargets(tf[-1])

# Sanity check of settings
if maxDefocus > minDefocus:
    minDefocus, maxDefocus = maxDefocus, minDefocus
if maxTilt < minTilt:
    minTilt, maxTilt = maxTilt, minTilt
if (maxTilt - startTilt) % step != 0 or (startTilt - minTilt) % step != 0:
    maxTilt = round(int((maxTilt - startTilt) / step) * step + startTilt, 1)
    minTilt = round(int((minTilt - startTilt) / step) * step + startTilt, 1)
    print(f"WARNING: Tilt increment does not divide evenly into tilt range. Tilt range will be adjusted to: {minTilt}, {maxTilt}")

### Recovery data
recoverInput = 0
recover = False
realign = False
if savedRun != False and (resume["sec"] > 0 or resume["pos"] > 0):
    recoverInput = sem.YesNoBox("The target file contains recovery data. Do you want to attempt to continue the acquisition? Tracking accuracy might be impacted.")
    if recoverInput == 1:
        recover = True
        while sem.ReportFileNumber() > 0:
            sem.CloseFile()

        stageX, stageY, stageZ = sem.ReportStageXYZ()
        if abs(stageX - float(targets[0]["stageX"])) > 1.0 or abs(stageY - float(targets[0]["stageY"])) > 1.0: # test if stage was moved (with 1 micron wiggle room)
            userRealign = sem.YesNoBox("It seems that the stage was moved since stopping acquisition. Do you want to realign to the tracking target before resuming? This will also reset prediction parameters reducing tracking accuracy.")    
            realign = True if userRealign == 1 else False
    else:
        sem.AllowFileOverwrite(1)

### Start setup
dumpVars(os.path.splitext(os.path.basename(tf[-1]))[0])                                         # write settings vars to text file

sem.ResetClock()

targetDefocus = maxDefocus                                                                      # use highest defocus for tracking TS
sem.SetTargetDefocus(targetDefocus)

# Collect exposure settings
expTime = sem.ReportExposure("R")[0]

if recover:
    log("##### Recovery attempt of PACEtomo with parameters: #####", style=1)
else:
    log("##### Starting new PACEtomo with parameters: #####", style=1)
log(f"Start: {startTilt} deg - Min/Max: {minTilt}/{maxTilt} deg ({step} deg increments)")
log(f"Data points used: {dataPoints}")
log(f"Target defocus range (min/max/step): {minDefocus}/{maxDefocus}/{stepDefocus}")
log(f"Sample pretilt (rotation): {pretilt} ({rotation})")
log(f"Tilt axis offset: {round(tiltAxisOffset, 3)}")
log(f"Focus correction slope: {focusSlope}")
log(f"Exposure time per tilt: {round(expTime, 3)} s (total: {round(expTime * int((maxTilt - minTilt) / step), 3)} s)")

if trackMag > 0:
    log("WARNING: A magnification offset for the tracking target changes the Low Dose Record mode temporarily. Please double-check your Record mode in case the script is stopped prematurely or crashes!")

if startTilt * pretilt > 0:
    log("WARNING: Start tilt and pretilt have the same sign! If you want to compensate for the pretilt, the start tilt should have the opposite sign!")

branchsteps = max(maxTilt - startTilt, abs(minTilt - startTilt)) / groupSize / step

### Create run file
counter = 1
while os.path.exists(os.path.join(curDir, fileStem + "_run" + str(counter).zfill(2) + ".txt")):
    counter += 1
runFileName = os.path.join(curDir, fileStem + "_run" + str(counter).zfill(2) + ".txt")

### Initial actions
if not recover:
    log("Moving to target area...")

    sem.SetCameraArea("V", "F")                                                                 # set View to Full for Eucentricity
    sem.MoveToNavItem(navID)
    log("Refining eucentricity...")
    sem.Eucentricity(1)
    sem.UpdateItemZ()
    sem.RestoreCameraSet("V")

    log("Realigning to target 1...")
    if alignToP:
        x, y, binning, exp, *_ = sem.ImageProperties("P")
        sem.SetExposure("V", exp)
        sem.SetBinning("V", int(binning))
        sem.V()
        sem.CropCenterToSize("A", int(x), int(y))
        alignTo("P", debug)
        sem.RestoreCameraSet("V")
        if refineVec and tgtPattern and size is not None:
            if float(sem.ReportDefocus()) < -50:
                log("WARNING: Large defocus offsets for View can cause additional offsets in image shift upon mag change.")
            size = int(size)
            log("Refining target pattern...")
            sem.GoToLowDoseArea("R")
            ISX0, ISY0, *_ = sem.ReportImageShift()
            SSX0, SSY0 = sem.ReportSpecimenShift()
            log(f"Vector A: ({vecA0}, {vecA1})")
            shiftx = size * vecA0
            shifty = size * vecA1
            sem.ImageShiftByMicrons(shiftx, shifty)

            sem.V()
            alignTo("P", debug)
            sem.GoToLowDoseArea("R")

            SSX, SSY = sem.ReportSpecimenShift()
            SSX -= SSX0
            SSY -= SSY0        
            if np.linalg.norm([shiftx - SSX, shifty - SSY]) > 0.5:
                log("WARNING: Refined vector differs by more than 0.5 microns! Original vectors will be used.")
            else:
                vecA0, vecA1 = (round(SSX / size, 4), round(SSY / size, 4))
                log(f"Refined vector A: ({vecA0}, {vecA1})")

                sem.SetImageShift(ISX0, ISY0)                                                   # reset IS to center position
                log(f"Vector B: ({vecB0}, {vecB1})")
                shiftx = size * vecB0
                shifty = size * vecB1
                sem.ImageShiftByMicrons(shiftx, shifty)

                sem.V()
                alignTo("P", debug)
                sem.GoToLowDoseArea("R")

                SSX, SSY = sem.ReportSpecimenShift()
                SSX -= SSX0
                SSY -= SSY0
                if np.linalg.norm([shiftx - SSX, shifty - SSY]) > 0.5:
                    log("WARNING: Refined vector differs by more than 0.5 microns! Original vectors will be used.")
                else:
                    vecB0, vecB1 = (round(SSX / size, 4), round(SSY / size, 4))
                    log(f"Refined vector B: ({vecB0}, {vecB1})")

                    targetNo = 0
                    for i in range(-size,size+1):
                        for j in range(-size,size+1):
                            if i == j == 0: continue
                            targetNo += 1
                            SSX = i * vecA0 + j * vecB0
                            SSY = i * vecA1 + j * vecB1
                            targets[targetNo]["SSX"] = str(SSX)
                            targets[targetNo]["SSY"] = str(SSY)
                    log("NOTE: Target pattern was overwritten using refined vectors.")
            sem.SetImageShift(ISX0, ISY0)                                                       # reset IS to center position
    else:
        #sem.RealignToOtherItem(navID, 1) # <= sometimes unreliable
        realignTo(nav_id=navID, target=targets[0])

    if measureGeo:
        log("Measuring geometry...")
        if int(round(float(sem.ReportTiltAngle()))) != 0:
            sem.TiltTo(0)
        if len(geoPoints) > 0 and "SSX" in geoPoints[0].keys():                                 # if there are geo points in tgts file, adjust format from dict to list
            geoPoints = [[point["SSX"], point["SSY"]] for point in geoPoints]
        if len(geoPoints) < 3 and tgtPattern and size is not None:
            if size > 1:
                geoPoints.append([0.5 * (vecA0 + vecB0), 0.5 * (vecA1 + vecB1)])
            geoPoints.append([(size - 0.5) * (vecA0 + vecB0), (size - 0.5) * (vecA1 + vecB1)])
            geoPoints.append([(size - 0.5) * (vecA0 - vecB0), (size - 0.5) * (vecA1 - vecB1)])
            geoPoints.append([(size - 0.5) * (-vecA0 + vecB0), (size - 0.5) * (-vecA1 + vecB1)])
            geoPoints.append([(size - 0.5) * (-vecA0 - vecB0), (size - 0.5) * (-vecA1 - vecB1)])

        # Clean geo_points beyond image shift limit
        geoPoints = [point for point in geoPoints if np.linalg.norm(np.array([point[0], point[1]], dtype=float)) < imageShiftLimit]

        if len(geoPoints) >= 3:
            geoXYZ = [[], [], []]
            sem.GoToLowDoseArea("R")
            ISX0, ISY0, *_ = sem.ReportImageShift()
            for i in range(len(geoPoints)):
                sem.ImageShiftByMicrons(geoPoints[i][0], geoPoints[i][1])
                sem.G(-1)
                defocus, *_ = sem.ReportAutoFocus()
                drift = sem.ReportFocusDrift()
                if abs(defocus) >= 0.01 and np.linalg.norm(drift) >= 0.01:
                    geoXYZ[0].append(geoPoints[i][0])
                    geoXYZ[1].append(geoPoints[i][1])
                    geoXYZ[2].append(defocus)
                else:
                    log("WARNING: Measured defocus is 0. This geo point will not be considered.")
                sem.SetImageShift(ISX0, ISY0)                                                   # reset IS to center position
            if len(geoXYZ[0]) >= 3:
                ##########
                # Source: https://math.stackexchange.com/q/99317
                # subtract out the centroid and take the SVD, extract the left singular vectors, the corresponding left singular vector is the normal vector of the best-fitting plane
                svd = np.linalg.svd(geoXYZ - np.mean(geoXYZ, axis=1, keepdims=True))
                left = svd[0]
                norm = left[:, -1]
                ##########        
                log(f"Fitted plane into cloud of {len(geoXYZ[0])} points ({len(geoPoints) - len(geoXYZ[0])} discarded).")
                log(f"Normal vector: {norm}")

                # Errors
                errors = []
                for point in zip(*geoXYZ):
                    errors.append(np.dot(norm, point - np.mean(geoXYZ, axis=1)) ** 2)
                log(f"Fitting error: {np.mean(errors)}")

                if debug:
                    log("DEBUG:\nGeo points [x, y, z, err]:")
                    for point in zip(*geoXYZ, errors):
                        log(f"# {point}", color=1)

                # Calculate pretilt and rotation
                sign = 1 if norm[1] <= 0 else -1
                pretilt = round(sign * np.degrees(np.arccos(norm[2])), 1)
                log(f"Estimated pretilt: {pretilt} degrees", style=1)
                rotation = round(-np.degrees(np.arctan(norm[0]/norm[1])), 1)
                log(f"Estimated rotation: {rotation} degrees", style=1)

                if startTilt * pretilt > 0:
                    log("WARNING: Start tilt and pretilt have the same sign! If you want to compensate for the pretilt, the start tilt should have the opposite sign!")
            else:
                log("WARNING: Not enough geo points could be checked successfully. Geometry could not be measured.")
        else:
            log("WARNING: Not enough geo points were defined. Geometry could not be measured.")

    if autoStartTilt or tiltTargets != 0:
        startTiltOri = startTilt
        if tiltTargets != 0:
            # Use tilt at which targets were selected as startTilt
            startTilt = tiltTargets
        elif autoStartTilt:
            # Adjust start tilt to compensate for measured pretilt
            startTilt = -int(round(np.degrees(np.arctan(np.sin(np.radians(pretilt)) * np.cos(np.radians(rotation)) / np.cos(np.radians(pretilt))))))
        maxTilt = np.clip(maxTilt - startTiltOri + startTilt, -int(tiltLimit), int(tiltLimit))
        minTilt = np.clip(minTilt - startTiltOri + startTilt, -int(tiltLimit), int(tiltLimit))

        log("WARNING: Automatically adjusted tilt series parameters!")
        log(f"Start: {startTilt} deg - Min/Max: {minTilt}/{maxTilt} deg ({step} deg increments)", style=1)

        # Update branch steps
        branchsteps = max(maxTilt - startTilt, abs(minTilt - startTilt)) / groupSize / step

    log("Tilting to start tilt angle...")
    # backlash correction
    sem.V()
    sem.Copy("A", "O")

    curTilt = int(round(float(sem.ReportTiltAngle())))

    # Walk up if necessary
    while abs(startTilt - curTilt) > 10:
        log(f"DEBUG: Doing walkup to {curTilt + (10 if startTilt > 0 else -10)}...")
        sem.TiltTo(curTilt + (10 if startTilt > 0 else -10))
        sem.V()
        alignTo("O", debug)
        sem.V()
        sem.Copy("A", "O")
        curTilt = int(round(float(sem.ReportTiltAngle())))

    sem.TiltTo(startTilt - step)
    sem.TiltTo(startTilt)

    sem.V()
    alignTo("O", debug)
    sem.GoToLowDoseArea("R")

    if not tgtPattern and previewAli:
        sem.LoadOtherMap(navID, "O")                                                            # preview ali before first tilt image is taken
        sem.AcquireToMatchBuffer("O")                                                           # in case view image was saved for tracking target
        alignTo("O", debug)

    ISX0, ISY0, *_ = sem.ReportImageShift()
    SSX0, SSY0 = sem.ReportSpecimenShift()

    sem.G()
    focus0 = float(sem.ReportDefocus())
    positionFocus = focus0                                                                      # set maxDefocus as focus0 and add focus steps in loop
    minFocus0 = focus0 - maxDefocus + minDefocus

    sem.GoToLowDoseArea("R")
    s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
    is2ssMatrix = np.array(sem.ISToSpecimenMatrix(0)).reshape((2, 2))
    ss2isMatrix = np.array(sem.SpecimenToISMatrix(0)).reshape((2, 2))
    camX, camY, *_ = sem.CameraProperties()
    c2ssMatrix = np.array(sem.CameraToSpecimenMatrix(0)).reshape((2, 2))
    ss2cMatrix = np.array(sem.SpecimenToCameraMatrix(0)).reshape((2, 2))
    if debug:
        log("DEBUG: Conversion matrices:")
        log(f"    Stage to Specimen: {s2ssMatrix}", color=1)
        log(f"    IS to Specimen: {is2ssMatrix}", color=1)
        log(f"    Image to Specimen: {c2ssMatrix}", color=1)
        log(f"    Specimen to Camera: {ss2cMatrix}", color=1)

    if previewAli:
        sem.SetDefocus(min(focus0, focus0 - 5 - targetDefocus))                                 # set defocus for Preview to at least -5 micron
### Target setup
    log(f"Setting up {len(targets)} targets...")

    posTemplate = {"SSX": 0, "SSY": 0, "focus": 0, "z0": 0, "n0": 0, "shifts": [], "angles": [], "ISXset": 0, "ISYset": 0, "ISXali": 0, "ISYali": 0, "dose": 0, "sec": 0, "skip": False}
    position = []
    skippedTgts = 0
    for i, tgt in enumerate(targets):
        position.append([])
        position[-1].append(copy.deepcopy(posTemplate))

        log(f"Target {i + 1}...")
        skip = False
        if "skip" in tgt.keys() and tgt["skip"] == "True":
            log(f"WARNING: Target [{str(i + 1).zfill(3)}] was set to be skipped.")
            skip = True
        if "SSX" not in tgt.keys() and "stageX" in tgt.keys():                                  # if SS coords are missing but stage coords are present, calc SS coords
            tgt["SSX"], tgt["SSY"] = s2ssMatrix @ np.array([float(tgt["stageX"]) - float(targets[0]["stageX"]), float(tgt["stageY"]) - float(targets[0]["stageY"])])
        if np.linalg.norm(np.array([tgt["SSX"], tgt["SSY"]], dtype=float)) > imageShiftLimit - alignLimit:
            log(f"WARNING: Target [{str(i + 1).zfill(3)}] is too close to the image shift limit. This target will we skipped.")
            skip = True

        if skip: 
            position[-1][0]["skip"] = True
            position[-1].append(copy.deepcopy(position[-1][0]))
            position[-1].append(copy.deepcopy(position[-1][0]))
            skippedTgts += 1
            continue

        if tiltTargets == 0:
            tiltScaling = np.cos(np.radians(pretilt * np.cos(np.radians(rotation)) + startTilt)) / np.cos(np.radians(pretilt * np.cos(np.radians(rotation)))) # stretch shifts from 0 tilt to startTilt
        else:
            tiltScaling = 1
        log(f"DEBUG: Tilt scaling to start tilt [{startTilt}]: {tiltScaling}")

        sem.ImageShiftByMicrons(float(tgt["SSX"]), float(tgt["SSY"]) * tiltScaling)             # apply relative shifts to find out absolute IS after realign to item
        if (previewAli or viewAli):                                                             # adds initial dose, but makes sure start tilt image is on target
            if alignToP:
                x, y, binning, exp, *_ = sem.ImageProperties("P")
                sem.SetExposure("V", exp)
                sem.SetBinning("V", int(binning))
                sem.V()
                sem.CropCenterToSize("A", int(x), int(y))
                alignTo("P", debug)
                sem.RestoreCameraSet("V")
            else:
                if "viewfile" in tgt.keys() and viewAli and i != 0:                             # skip for tracking target since it was already aligned after tilt to startTilt   
                    sem.ReadOtherFile(0, "O", tgt["viewfile"])                                  # reads view file for first AlignTo instead
                    sem.V()
                    alignTo("O", debug)
                    ASX, ASY = sem.ReportAlignShift()[4:6]
                    log(f"Target alignment (View) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")    
                if "tgtfile" in tgt.keys() and previewAli:                
                    sem.ReadOtherFile(0, "O", tgt["tgtfile"])                                   # reads tgt file for first AlignTo instead
                    sem.L()
                    alignTo("O", debug)
                    AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
                    log(f"Target alignment (Prev) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")
                elif "viewfile" in tgt.keys() and previewAli:
                    # Use align between mags to align preview image to view image
                    if not viewAli:
                        #sem.GoToLowDoseArea("V")                                                # If ReadOtherFile while in Record, pixel size of Record is used and AlignBetweenMags fails (seems to be fixed in 4.2beta from 14.08.2024)
                        sem.ReadOtherFile(0, "O", tgt["viewfile"])                              # reads view file for first AlignTo instead
                    else:
                        # If View image was already aligned, take new centered View image at startTilt and use as reference instead
                        sem.V()
                        sem.Copy("A", "O")
                    # Check defocus offset
                    sem.GoToLowDoseArea("R")                                                    # Switch to R before applying defocus offset to not mess with potential mP/nP offsets between View and Rec
                    defocus_offset = max(-10, sem.ReportLDDefocusOffset("V"))
                    if defocus_offset != 0:
                        sem.ChangeFocus(defocus_offset)                                             # Higher defocus for better correlation, but max at 10 to avoid major distortions
                    sem.L()
                    sem.AlignBetweenMags("O", -1, -1, -1)
                    AISX, AISY, ASX, ASY = sem.ReportAlignShift()[2:6]
                    if defocus_offset != 0:
                        sem.ChangeFocus(-defocus_offset)                                            # Reset focus
                    log(f"Target alignment (Pv2V) error in X | Y: {round(ASX, 0)} nm | {round(ASY, 0)} nm")           

                # Save preview image as new reference
                if refFromPreview:
                    sem.OpenNewFile(os.path.splitext(tgt["tgtfile"])[0] + "_tempref.mrc")
                    sem.S()
                    sem.CloseFile()
                    position[-1][0]["ISXali"] = AISX                                            # Save shifts to real reference
                    position[-1][0]["ISYali"] = AISY 

            sem.GoToLowDoseArea("R")
        ISXset, ISYset, *_ = sem.ReportImageShift()
        SSX, SSY = sem.ReportSpecimenShift()
        sem.SetImageShift(ISX0, ISY0)                                                           # reset IS to center position    

        z0_ini = np.tan(np.radians(pretilt)) * (np.cos(np.radians(rotation)) * float(tgt["SSY"]) - np.sin(np.radians(rotation)) * float(tgt["SSX"]))
        correctedFocus = positionFocus - z0_ini * np.cos(np.radians(startTilt)) - float(tgt["SSY"]) * np.sin(np.radians(startTilt))

        position[-1][0]["SSX"] = float(SSX)
        position[-1][0]["SSY"] = float(SSY)
        position[-1][0]["focus"] = correctedFocus
        position[-1][0]["z0"] = z0_ini                                                          # offset from eucentric height (will be refined during collection)
        position[-1][0]["n0"] = float(tgt["SSY"])                                               # offset from tilt axis
        position[-1][0]["ISXset"] = float(ISXset)
        position[-1][0]["ISYset"] = float(ISYset)

        position[-1].append(copy.deepcopy(position[-1][0]))                                     # plus and minus branch start with same values
        position[-1].append(copy.deepcopy(position[-1][0]))

        position[-1][1]["n0"] -= taOffsetPos
        position[-1][2]["n0"] -= taOffsetNeg

        positionFocus += stepDefocus                                                            # adds defocus step between targets and resets to initial defocus if minDefocus is surpassed
        if positionFocus > minFocus0: positionFocus = focus0

### Start tilt
    log("Start tilt series...", style=1)
    log(f"Tilt step 1 out of {int((maxTilt - minTilt) / step + 1)} ({startTilt} deg)...")
    sem.SetStatusLine(1, f"Tilt step: 1 / {int((maxTilt - minTilt) / step + 1)}")

    dewarFillTime = 0
    maxProgress = ((maxTilt - minTilt) / step + 1) * (len(position) - skippedTgts)
    resumePercent = 0
    startTime = sem.ReportClock()
    lastSlitCheck = startTime

    geo = [[], [], []]

    plustilt = minustilt = startTilt
    Tilt(startTilt)

    if refineGeo:
        if len(geo[2]) >= 3:                                                                    # if number of points > 3: fit z = a * x + b * y
            log("Refining geometry...")
            log(f"{len(geo[2])} usable CtfFind results found.")
            if len(geo[2]) >= parabolTh:                                                        # if number of points > 6: fit z = a * x + b * y + c * (x**2) + d * (y**2) + e * x * y
                log("Fitting paraboloid...")
                geoF = geoPara
            else:                            
                log("Fitting plane...")
                geoF = geoPlane

            p, cov = optimize.curve_fit(geoF, [geo[0], geo[1]], [z - geo[2][0] for z in geo[2]])

            ss = 0
            for i in range(0, len(geo[2])):
                ss += (geo[2][i] - geo[2][0] - geoF([geo[0][i], geo[1][i]], *p))**2
            rmse = np.sqrt(ss / len(geo[2]))

            log("Fit parameters: " + " # ".join(p.astype(str)))
            log(f"RMSE: {round(rmse, 3)}")

            for pos in range(0, len(position)):                                                 # calculate and adjust refined z0
                zs = geoF([position[pos][1]["SSX"], position[pos][1]["SSY"]], *p)
                z0_ref = position[pos][1]["z0"] + zs * np.cos(np.radians(startTilt)) + position[pos][1]["SSY"] * np.sin(np.radians(startTilt))

                position[pos][1]["z0"] = z0_ref
                position[pos][2]["z0"] = z0_ref
        else: 
            log(f"WARNING: Not enough reliable CtfFind results ({len(geo[2])}) to refine geometry. Continuing with initial geometry model.")

    startstep = 0
    substep = [0, 0]
    posResumed = -1
    resumePN = 0

### Recovery attempt
else:
    if realign:
        sem.MoveToNavItem(navID)
        if alignToP:
            x, y, binning, exp, *_ = sem.ImageProperties("P")
            sem.SetExposure("V", exp)
            sem.SetBinning("V", int(binning))
            sem.V()
            sem.CropCenterToSize("A", int(x), int(y))
            alignTo("P", debug)
            sem.RestoreCameraSet("V")
        else:
            sem.RealignToOtherItem(navID, 1)
    position = []
    skippedTgts = 0
    for pos in range(len(targets)):
        position.append([{},{},{}])
        for i in range(2):
            position[-1][i+1]["SSX"] = float(savedRun[pos][i]["SSX"])
            position[-1][i+1]["SSY"] = float(savedRun[pos][i]["SSY"])
            position[-1][i+1]["focus"] = float(savedRun[pos][i]["focus"])
            position[-1][i+1]["z0"] = float(savedRun[pos][i]["z0"])
            position[-1][i+1]["n0"] = float(savedRun[pos][i]["n0"])
            if savedRun[pos][i]["shifts"] != "" and not realign:
                position[-1][i+1]["shifts"] = [float(shift) for shift in savedRun[pos][i]["shifts"].split(",")]
            else:
                position[-1][i+1]["shifts"] = []
            if savedRun[pos][i]["angles"] != "" and not realign:
                position[-1][i+1]["angles"] = [float(angle) for angle in savedRun[pos][i]["angles"].split(",")]
            else:
                position[-1][i+1]["angles"] = []
            position[-1][i+1]["ISXset"] = float(savedRun[pos][i]["ISXset"])
            position[-1][i+1]["ISYset"] = float(savedRun[pos][i]["ISYset"])
            position[-1][i+1]["ISXali"] = float(savedRun[pos][i]["ISXali"])
            position[-1][i+1]["ISYali"] = float(savedRun[pos][i]["ISYali"])
            position[-1][i+1]["dose"] = float(savedRun[pos][i]["dose"])
            position[-1][i+1]["sec"] = int(savedRun[pos][i]["sec"])
            position[-1][i+1]["skip"] = True if savedRun[pos][i]["skip"] == "True" or targets[pos]["skip"] == "True" else False

        sem.AreaForCumulRecordDose(pos + 1)                                                     # set dose accumulator to highest recorded prior dose
        sem.AccumulateRecordDose(max(position[-1][1]["dose"], position[-1][2]["dose"]))

        if targets[pos]["skip"] == "True":
            skippedTgts += 1

    startstep = (resume["sec"] - 1) // (2 * groupSize)                                          # figure out start values for branch loops
    substep = [min((resume["sec"] - 1) % (2 * groupSize), groupSize), (resume["sec"] - 1) % (2 * groupSize) // (groupSize + 1)]

    realTilt = float(savedRun[resume["pos"]][0]["angles"].split(",")[-1])
    if np.floor(realTilt) % step == 0:                                                          # necessary because angles array was switched to realTilt
        lastTilt = np.floor(realTilt)
    else:
        lastTilt = np.ceil(realTilt)
    plustilt = resumePlus = lastTilt                                                            # obtain last angle from savedRun in case position["angles"] was reset
    if substep[0] < groupSize:                                                                  # subtract step when stopped during positive branch
        plustilt -= step
        resumePN = 1                                                                            # indicator which branch was interrupted
        sem.TiltTo(plustilt)
    if savedRun[pos][1]["angles"] != "":
        realTilt = float(savedRun[resume["pos"]][1]["angles"].split(",")[-1])
        if np.floor(realTilt) % step == 0:                                                      # necessary because angles array was switched to realTilt
            lastTilt = np.floor(realTilt)
        else:
            lastTilt = np.ceil(realTilt)
        minustilt = resumeMinus = lastTilt
        if substep[0] == groupSize:                                                             # add step when stopped during negative branch
            minustilt += step
            resumePN = 2
            sem.TiltTo(minustilt)
    else:
        minustilt = resumeMinus = startTilt
        resumePN = 1

    posResumed = resume["pos"] + 1

    maxProgress = ((maxTilt - minTilt) / step + 1) * (len(position) - skippedTgts)
    # progress = collected images * (positions - skipped positions) + current position - skipped positions scaled assuming homogeneous distribution of skipped positions
    progress = resume["sec"] * (len(position) - skippedTgts) + resume["pos"] - skippedTgts * resume["pos"] / len(position)
    resumePercent = round(100 * (progress / maxProgress), 1)

    sem.GoToLowDoseArea("R")
    origMag, *_ = sem.ReportMag()
    s2ssMatrix = np.array(sem.StageToSpecimenMatrix(0)).reshape((2, 2))
    is2ssMatrix = np.array(sem.ISToSpecimenMatrix(0)).reshape((2, 2))
    camX, camY, *_ = sem.CameraProperties()
    c2ssMatrix = np.array(sem.CameraToSpecimenMatrix(0)).reshape((2, 2))
    ss2cMatrix = np.array(sem.SpecimenToCameraMatrix(0)).reshape((2, 2))
    if debug:
        log("DEBUG: Conversion matrices:")
        log(f"    Stage to Specimen: {s2ssMatrix}", color=1)
        log(f"    IS to Specimen: {is2ssMatrix}", color=1)
        log(f"    Image to Specimen: {c2ssMatrix}", color=1)
        log(f"    Specimen to Camera: {ss2cMatrix}", color=1)

    focus0 = (position[0][1]["focus"] + position[0][2]["focus"]) / 2                            # get estimate for original microscope focus value by taking average of both branches of tracking target

    startTime = sem.ReportClock()
    lastSlitCheck = startTime


### Tilt series
for i in range(startstep, int(np.ceil(branchsteps))):
    for j in range(substep[0], groupSize):
        plustilt += step
        if all([pos[1]["skip"] for pos in position]): continue
        log(f"\nTilt step {i * 2 * groupSize + j + 1 + 1} out of {int((maxTilt - minTilt) / step + 1)} ({plustilt} deg)...", style=1)
        sem.SetStatusLine(1, f"Tilt step: {i * 2 * groupSize + j + 1 + 1} / {int((maxTilt - minTilt) / step + 1)}")
        Tilt(plustilt)
    for j in range(substep[1], groupSize):
        minustilt -= step
        if all([pos[2]["skip"] for pos in position]): continue
        log(f"\nTilt step {i * 2 * groupSize + j + groupSize + 1 + 1} out of {int((maxTilt - minTilt) / step + 1)} ({minustilt} deg)...", style=1)
        sem.SetStatusLine(1, f"Tilt step: {i * 2 * groupSize + j + groupSize + 1 + 1} / {int((maxTilt - minTilt) / step + 1)}")
        Tilt(minustilt)
    substep = [0, 0]                                                                            # reset sub steps after recovery
    if coldFEG: checkColdFEG()                                                                  # check for flashing at the end of each step

### Finish
sem.ClearStatusLine(0)
if trackMag > 0:
    sem.RestoreLowDoseParams("R")                                                               # restore record mag before script just in case
sem.TiltTo(0)
sem.SetDefocus(focus0)
sem.SetImageShift(0, 0)
sem.CloseFile()
updateTargets(runFileName, targets)

# Format final tilt stacks
if delFinalStack:
    for target in targets:
        if checkFrames(target["tsfile"]):
            os.remove(target["tsfile"])
            log(f"NOTE: {target['tsfile']} was deleted. Please use saved frames to generate the tilt series.")
else:
    if sortByTilt or binFinalStack > 1:
        for target in targets:
            if sortByTilt:
                sortTS(target["tsfile"])
            if binFinalStack > 1:
                binStack(target["tsfile"], binFinalStack)

totalTime = round(sem.ReportClock() / 60, 1)
perTime = round(totalTime / len(position), 1)
if recoverInput == 1:
    perTime = f"since recovery: {perTime}"
log(datetime.now().strftime("%d.%m.%Y %H:%M:%S"))
log(f"##### All tilt series completed in {totalTime} min ({perTime} min per tilt series) #####", color=3, style=1)
sem.SaveLog()
sem.Exit()