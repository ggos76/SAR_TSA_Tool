#!/usr/bin/env python
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. May 2025                                      -
 * -------------------
 ---------------------------------------------------
'''

# ----------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# ----------------------------------------------------------------------------------------------
# A) Input Files
TSA_stacks_folder = r"E:\QC03M_7224_4504_CP01_A\7_TSA_stacks"
no_data_value = -32768.0000
prefix = "QC03M_7224_4504_CP01_A_"                              # leave blank if no prefix is needed

# B) spatial subset option # Accepted values are "yes" or "no"
subset_input_data = "yes"
AOI_vector_file = r"D:\RCMP_data_RCM\QC_03m\stacks\QC03M_7224_4504_CP01_A\AOI\AOI_QC03m_7224_4504_CP01_A_final_crop_tight_UTM18TD000.pix"    # mandatory when subset_input_data = "yes"
AOI_segment_number = 2

# C) Apply masking
# C.1) DEM file for slope calculation.
apply_slope_mask = "no"
DEM_file = r"D:\RCMP_data_aux\QC_prj\DEM\QC_Glo30DEM_LatLong.tif"
DEM_elevation_channel = 1
slope_max = 15.0                # Slopes above slope_max will be masked.

# C.2) CANUSA border buffer masking - area outise the buffer will be masked
apply_canusa_border_mask = "no"
canusa_border_file = r"D:\RCMP_data_aux\CANUSA_border_LatLongWGS84_buffer_2km.pix"

# C.3) Other exclusion masks
apply_water_mask = "yes"
water_mask_file = r"D:\RCMP_data_aux\QC_prj\Land_Cover\QC_20231004_Landsat\2023_LandCover_water_UTM19TD000.pix"  # This is a vector file with a single segment

#D) Other options
# Generate overviews - either yes or no,
generate_overviews = "yes"
delete_intermediary_files = "yes"   

# Behaviour when output file exists 
if_file_exists = "skip"     # Valid options are "skip" or "regenerate"

# -----------------------------------------------------------------------------------------------------
#  Part 2: Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
All input stack must: 
 1) Corresponds to the same geographical region: 
 2) have the same pojection
 3) have the same number of lines and columns

'''
# -----------------------------------------------------------------------------------------------------------------
#  Part 3: imports
# -----------------------------------------------------------------------------------------------------------------
import sys
import os
import time
import locale
import shutil
import fnmatch
from pathlib import Path

import pci
from pci.pcimod import pcimod
from pci.clip import clip
from pci.fexport import fexport
from pci.ras2poly import ras2poly
from pci.reproj import reproj
from pci.model import model
from pci.slasp import slasp
from pci.thr import thr
from pci.bit2poly import bit2poly
from pci.pyramid import pyramid
from pci.api import gobs
from pci.api import cts
from pci.api.cts import crs_to_mapunits
from pci.exceptions import PCIException
from pci.api import datasource as ds

from TSA_utilities.SAR_TSA_utilities_definitions import file_size_check
from TSA_utilities.SAR_TSA_utilities_definitions import ortho_raster_footprint
from TSA_utilities.SAR_TSA_utilities_definitions import stack_chans_masking
from TSA_utilities.SAR_TSA_utilities_definitions import get_folder_proctime_and_size
from TSA_utilities.SAR_TSA_version_control import version_control

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -----------------------------------------------------------------------------------------------------------------
#  Part 4: Parameters Validation
# -----------------------------------------------------------------------------------------------------------------
TSA_04_start = time.time()
val_start = time.time()

# Hardcoded parameters
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list + no_validation_list

# Version control
vs_catalyst = pci.version
vs_python = sys.version_info[:3]
version_control (vs_catalyst, vs_python)

# A) Check if input file exists and the quantity type (intensity or coherence)
if not os.path.exists(TSA_stacks_folder):
    print ("Error - The TSA_stacks_folder  does not exist or the path/filename is wrong")
    sys.exit()

# B) subset input data
subset_input_data = subset_input_data.lower()
if subset_input_data not in yes_no_validation_list:
    print("Error - The specified subset option is invalid")
    print('Accepted values are: "yes" or "no"')
    sys.exit()
elif subset_input_data in yes_validation_list:
    subset_input_data = True
    if not os.path.exists(AOI_vector_file):
        print ("Error - The AOI_vector_file does not exist or the path/filename is wrong")
        sys.exit()
else: 
    subset_input_data = False


# C) slope mask
apply_slope_mask = apply_slope_mask.lower()
if apply_slope_mask not in yes_no_validation_list:
    print("Error - The specified subset option is invalid")
    print('Accepted values are: "yes" or "no"')
    sys.exit()
elif apply_slope_mask in yes_validation_list:
    apply_slope_mask = True

    if not os.path.exists(DEM_file):
        print ("Error - The DEM_file does not exist or the path/filename is wrong")
        sys.exit()
    if slope_max < 1 or slope_max > 89: 
        print ("Error - The slope_max value must be be within [1, 89]")
        sys.exit()
else: 
    apply_slope_mask = False

apply_canusa_border_mask = apply_canusa_border_mask.lower()
if apply_canusa_border_mask not in yes_no_validation_list:
    print("Error - The specified subset option is invalid")
    print('Accepted values are: "yes" or "no"')
    sys.exit()
elif apply_canusa_border_mask in yes_validation_list:
    apply_canusa_border_mask = True

    if not os.path.exists(canusa_border_file):
        print ("Error - The canusa_border_file does not exist or the path/filename is wrong")
        sys.exit()
else: 
    apply_canusa_border_mask = False

# C.3) Other exclusion masks
# TBC
if apply_water_mask not in yes_no_validation_list:
    print("Error - The specified subset option is invalid")
    print('Accepted values are: "yes" or "no"')
    sys.exit()
elif apply_water_mask in yes_validation_list:
    apply_water_mask = True

    if not os.path.exists(water_mask_file):
        print ("Error - The water_mask_file does not exist or the path/filename is wrong")
        sys.exit()
else: 
    apply_water_mask = False

# D) Other options
generate_overviews = generate_overviews.lower()
if generate_overviews in yes_validation_list:
    generate_overviews = True
else:
    generate_overviews = False

if if_file_exists.lower() not in ["skip","regenerate"]:     
    print('Error - valid options for existing_data are "skip" or "regenerate"')
    sys.exit()
else: 
    info_message_skip = ("   output file already exists - skip (if_file_exists = skip)")
    info_message_regn =  ("   output file already exists - regenerate (if_file_exists = regenerate)")

intermediary_files_list = []
delete_intermediary_files = delete_intermediary_files.lower()
if delete_intermediary_files in yes_validation_list: 
    delete_intermediary_files = True
elif delete_intermediary_files in no_validation_list: 
    delete_intermediary_files = False
else: 
    print ('Error - The delete_intermediary_files parameter must be set with "yes" or "no"')
    sys.exit()


# E) Creating the output folder and the output file name.  
    # Note: In anticipation for full automatation, we will create automatically the output folder and 
    # the output file name.      
one_up = Path(TSA_stacks_folder).resolve().parents[0]
stack_folder = "8_TSA_stacks_postproc"

Fld_output_stacks = os.path.join (one_up, stack_folder)
if not os.path.exists(Fld_output_stacks):
    os.makedirs(Fld_output_stacks)


# All validations have suceeded, a time log file is open.
out_folder_script4 = one_up
script4_procTime = os.path.join(out_folder_script4, prefix + "TSA_part_04_b_stack_data_postproc_processingTime.txt")
time_log = open(script4_procTime, "w")

val_stop = time.time()
string_1 = "Validation step (secs): " + str(round((val_stop - val_start), 2)) 
time_log.write("%s\n" % string_1)


# -----------------------------------------------------------------------------------------------------------------   
# -----------------------------------------------------------------------------------------------------------------   
print("\t")
print("-----------------------------------------------------------------------------------------------------------")
print("                                          Main program                                                     ")
print("-----------------------------------------------------------------------------------------------------------")
print("\t")

# A) retrieving the stack to process: 
print(time.strftime("%H:%M:%S") + " Finding the input stack to process")
input_stack_files = []
for root, dirs, files in os.walk(TSA_stacks_folder):
    for filename in fnmatch.filter(files, "*.pix"):
        input_stack_files.append(os.path.join(root, filename))

if len(input_stack_files) == 0:
    print ("Error - No stack_statistics.txt file found in the stack_lists_folder")
    sys.exit()

print("   Number of stack to process: " + str(len(input_stack_files)))
for ii in input_stack_files:
    print ("   " + ii)


# B) Extra verification to be sure the input stack are in a projection with meters units (i.e. Not LatLong) 

for ii in input_stack_files: 
    with ds.open_dataset(ii, ds.eAM_READ) as ds5:
        raster_MapProjection = crs_to_mapunits(ds5.crs)

raster_MapProjection2 = raster_MapProjection.lower().replace(" ","")
prj_check =  ["utm", "lcc"]
found = any (s in raster_MapProjection2 for s in prj_check)
if not found:
    print("no match found")
    sys.exit()

# C) Checking for files size (columns and lines)
print ("\t")
print("Files size verification")
files_list = input_stack_files
file_size_check (files_list)
    
string_1 = ("   " + time.strftime("%H:%M:%S") + " input files size verification")
time_log.write("%s\n" % string_1)
    
proc_stop_time = time.time()

# -----------------------------------------------------------------------------------------------------------------
# C) Subsetting the input data
# It is not mandatory to have the subset vector file in the same projection as the input data but it is 
# recommended. The same input file subsetted twice with the same subset vector but in different projections may 
# yield slighly differet output file size (lines x columns) 

if subset_input_data is True:
    print ("\t")
    print('----------------------------------------------------------------------------------------------------------')
    print('                                 Subsetting the input stack                                               ')
    print('----------------------------------------------------------------------------------------------------------')
    print ("\t")

    count = 1
    nb_files = str(len(input_stack_files))
    output_files_sub = []
    for ii in input_stack_files: 
        print ("\t")
        print(time.strftime("%H:%M:%S") + " Subsetting input stack " + str(count) + " of " + nb_files)
        # retrieve the number of channels to process
        with ds.open_dataset(ii, ds.eAM_READ) as ds2:
            aux = ds2.aux_data
            num_channels = ds2.chan_count
            AOI_MapProjection = crs_to_mapunits(ds2.crs)
    
        fili = ii
        dbic = [1, -num_channels]
        dbsl = []
        sltype = ""
        base = os.path.basename(ii)
        output_file = os.path.join(Fld_output_stacks, "s" + base)
        filo = output_file
        ftype = "pix"
        foptions = ""
        clipmeth = "LAYERVEC"
        clipfil = AOI_vector_file
        cliplay = [2]
        laybnds = "extents"
        coordtyp = ""
        clipul = ""
        cliplr = ""
        clipwh = ""
        initvalu = [no_data_value]
        setnodat = "Y";
        oclipbdy = "Y"

        print ("   input stack-->" + fili)
        print ("   output stack-->" + filo)

        output_files_sub.append(filo)  # That will conver all cases of skip / regenerate and new files

        if os.path.exists (filo) and if_file_exists == "skip": 
            print (info_message_skip)
        else: 
            if os.path.exists (filo) and if_file_exists == "regenerate": 
                print (info_message_regn)
                os.remove (filo)
            try: 
                clip (fili, dbic, dbsl, sltype, filo, ftype, foptions, clipmeth, clipfil, cliplay, \
                laybnds, coordtyp, clipul, cliplr, clipwh, initvalu, setnodat, oclipbdy )
            except PCIException as e:
                print (e)
            except Exception as e:
                print (e)
        count = count + 1
    
    # We now use the subsetted version of the input stacks.
    prefix = "s" + prefix  
    input_stack_files = output_files_sub


#---------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------
# D ) Creating footprints for the input stacks
print ("\t")
print('----------------------------------------------------------------------------------------------------------')
print('                             Creating vector footprints for the input stacks                              ')
print('----------------------------------------------------------------------------------------------------------')
print ("\t")

# Create the stack vector footprints (original projection and LatLong WGS84)

input_raster = input_stack_files[0]
output_folder = Fld_output_stacks

stack_template, stack_footprint_utm, stack_footprint_LatLong = ortho_raster_footprint (prefix, input_raster, output_folder)

print (stack_template)
print (stack_footprint_utm)
print (stack_footprint_LatLong)

if apply_slope_mask is True:
    print ("\t")
    print('----------------------------------------------------------------------------------------------------------')
    print('                 Calculation of terrain derivatives for steep slopes masking                              ')
    print('----------------------------------------------------------------------------------------------------------')
    print ("\t")


    # Create a subset of the DEM in LatLong 
    print (time.strftime("%H:%M:%S") + " Subsetting the input DEM using the stack footprint in LatLong")
    
    fili = DEM_file
    dbic = [1]
    dbsl = []
    sltype = ""

    filo = os.path.join(Fld_output_stacks, prefix + "DEM_for_terrain_derivatives_LatLong.pix")
    out_dem_terrain_derivatives = filo
    ftype = "PIX"
    foptions = ""
    clipmeth = "LAYERVEC"
    clipfil = stack_footprint_LatLong
    cliplay = [2]
    laybnds = "extents"
    coordtyp = ""
    clipul = ""
    cliplr = ""
    clipwh = ""
    initvalu = [no_data_value]
    setnodat = "n";
    oclipbdy = "n"
    
    try:
        clip (fili, dbic, dbsl, sltype, filo, ftype, foptions, clipmeth, clipfil, cliplay,laybnds, coordtyp,\
            clipul, cliplr, clipwh, initvalu, setnodat, oclipbdy )
        intermediary_files_list.append(filo)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)
    
    print (time.strftime("%H:%M:%S") + " Reprojecting the subsetted DEM to the stack projection")
 
    fili = out_dem_terrain_derivatives
    dbic = [1] 
    dbsl = []
    sltype = ""
    out_proj = raster_MapProjection.replace(" ", "")
    filo = os.path.join (Fld_output_stacks, prefix + "DEM_for_terrain_derivatives_" + out_proj + ".pix")
    ftype   =   ""  # uses PIX format by default
    foptions    =   ""  # no file options are used
    repmeth =   "BR"  # uses bounds and resolution method
    dbsz    =   []  # not used for BR method
    pxsz    =   [30, 30]  
    maxbnds =   "YES"    
    mapunits = raster_MapProjection
    llbounds    =   "NO"
    ulx =   ""
    uly =   ""
    lrx =   ""
    lry =   ""
    resample = ""  
    proc =   ""  # uses AUTO by default
    tipostrn = "" # uses CORNER tile positioning transformation with 10 meter stride

    try: 
        reproj( fili, dbic, dbsl, sltype, filo, ftype, foptions, repmeth, dbsz, pxsz, maxbnds,\
            mapunits, llbounds, ulx, uly, lrx, lry, resample, proc, tipostrn )
        intermediary_files_list.append(filo)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)

    #------------------------------------------------------------------------------------------
    # Creating a model on the fly to remove the 0 values on the subsetted DEM edges due 
    # to UTM - LantLong reprojection(s) and to avoid artifical high slopes on the edges.  
    model_set_nodata = os.path.join(Fld_output_stacks, "DEM_model_set_nodata.txt")
    model_lines = []
    string_mod = "if %1<=0 then"
    model_lines.append(string_mod)
    string_mod = "%1 = " + str(no_data_value)
    model_lines.append(string_mod)
    string_mod = "endif"
    model_lines.append(string_mod)
    with open(model_set_nodata, "w") as f:
        f.write("\n".join(model_lines))

    file = filo
    source = model_set_nodata
    undefval = []
    model(file, source, undefval)
  
    out_dem_terrain_derivatives = filo

    #------------------------------------------------------------------------------------------
    # We create the terrain derivatives for the subsetted (and reprojected) DEM file.
    # Add two 32R channels to receive the slopes and aspects layers.
    print (time.strftime("%H:%M:%S") + " Calculating the terrain derivatives (slope, aspect)")

    file = out_dem_terrain_derivatives
    pciop = "add"
    pcival = [0,0,0,2] 
    pcimod (file, pciop, pcival)
 
    filedem	= out_dem_terrain_derivatives
    dbec = [1]      # input DEM channel
    filo = out_dem_terrain_derivatives
    dboc = [2,3]	# output slope, aspect channels
    elevunit = 'METER'	# default elevation units
    elfactor	=	[]	# default scale, offset factors
    backelev	=	[no_data_value]
    zeroslop	=	[0.0]	# default

    try:
        slasp (filedem, dbec, filo, dboc, elevunit, elfactor, backelev, zeroslop)
        intermediary_files_list.append(filo)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)
    # -----------------------------------------------------------------------
    # -----------------------------------------------------------------------
    print (time.strftime("%H:%M:%S") + " Creating the exclusion mask for the steep slopes") 

    file = filo
    dbic = [2]
    dbob = []	# create new bitmap
    #tval =	[slope_max,90]	  # threshold range (min,max)
    tval =	[slope_max, 90]	  # threshold range (min,max)
    comp=	'OFF'	
    dbsn=	'slopes'	# output segment name
    dbsd=	'slopes_mask'	# output segment description

    try: 
        thr (file, dbic, dbob, tval, comp, dbsn, dbsd )   
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)

    print (time.strftime("%H:%M:%S") + " Converting the steep slopes mask to a vector file") 
    fili = file 
    dbib = [2]                            
    filo = os.path.join(Fld_output_stacks, prefix + "steep_slopes_mask_" + str(slope_max) + ".pix")  
    steep_slope_mask = filo
    smoothv  = "NO"                          
    dbsd     = "slope"    
    ftype    = ""                              
    foptions = ""                              

    try:
        bit2poly( fili, dbib, filo, smoothv, dbsd, ftype, foptions )
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)

    #-----------------------------------------------------------------------------
    # Ready to apply the slope mask
    nb_files = str(len(input_stack_files))
    count = 1
    for input_stack in input_stack_files: 
        print("\t") 
        print (time.strftime("%H:%M:%S") + " Applying the steep slopes masks, stack " + str(count) + " of " + nb_files)
        print ("   input stack --> " + input_stack )    
       
        mask_type = "exclusion"  
        mask_file = steep_slope_mask
        mask_seg_number = [2]

        stack_chans_masking (input_stack, mask_type, mask_file, mask_seg_number,no_data_value)
        
        count = count + 1
    
if apply_canusa_border_mask is True:
    print ("\t")
    print('----------------------------------------------------------------------------------------------------------')
    print('                            Masking the areas outside the CANUSA border buffer                            ')
    print('----------------------------------------------------------------------------------------------------------')
    print ("\t")
    
    # projection check and automatic reprojection here?
    nb_files = str(len(input_stack_files))
    count = 1
    for input_stack in input_stack_files: 
        print("\t") 
        print (time.strftime("%H:%M:%S") + " Masking stack " + str(count) + " of " + nb_files)
        print ("   input stack --> " + input_stack )    
       
        mask_type = "inclusion"  
        mask_file = canusa_border_file
        mask_seg_number = [2]

        stack_chans_masking (input_stack, mask_type, mask_file, mask_seg_number,no_data_value)
        
        count = count + 1

if apply_water_mask is True:
    print ("\t")
    print('----------------------------------------------------------------------------------------------------------')
    print('                                     Masking the waterbodies                                              ')
    print('----------------------------------------------------------------------------------------------------------')
    print ("\t")
    
    # projection check and automatic reprojection here?
    nb_files = str(len(input_stack_files))
    count = 1
    for input_stack in input_stack_files: 
        print("\t") 
        print (time.strftime("%H:%M:%S") + " Masking stack " + str(count) + " of " + nb_files)
        print ("   input stack --> " + input_stack )    
       
        mask_type = "exclusion"  
        mask_file = water_mask_file
        mask_seg_number = [2]

        stack_chans_masking (input_stack, mask_type, mask_file, mask_seg_number,no_data_value)
        
        count = count + 1


print("\t")
print ('----------------------------------------------------------------------------------------------------------')
print ('----------------------------------------------------------------------------------------------------------')
print(time.strftime("%H:%M:%S") + "   All process completed")
print("\t")
TSA_04_stop = time.time()

ellapse_time_seconds = round((TSA_04_stop - TSA_04_start), 2)
ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

print("Processing time (seconds): " + str(ellapse_time_seconds))
print("Processing time (minutes): " + str(ellapse_time_minutes))
print("Processing time (hours): " + str(ellapse_time_hours))
 

string1 = ("\t"); 
time_log.write("%s\n" % string1)
string1 = ("------------------------------------------------------------------------------------------------------------")
time_log.write("%s\n" % string1)
string1 = "TSA_04 total processing time (secs):;" + str(ellapse_time_seconds)
time_log.write("%s\n" % string1)
string1 = "TSA_04 total processing time (mins):;" + str(ellapse_time_minutes)
time_log.write("%s\n" % string1)
string1 = "TSA_04 total processing time hours):;" + str(ellapse_time_hours)
time_log.write("%s\n" % string1)


folder = Fld_output_stacks
proc_start_time = TSA_04_start
proc_stop_time
out_folder_time_size, size_mb = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
string_1  = (time.strftime("%H:%M:%S") + " Total size of the output folder: " + size_mb + " MB")
time_log.write("%s\n" % string_1)
time_log.close()