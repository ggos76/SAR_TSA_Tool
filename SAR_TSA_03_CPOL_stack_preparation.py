#!/usr/bin/env python
'''-----------------------------------------------------------------------------------------------------------------
 * Gabriel Gosselin, CCMEO\ CCRS. 2024-2025                                                                        -
 * -----------------------------------------------------------------------------------------------------------------
'''
# -----------------------------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# -----------------------------------------------------------------------------------------------------------------

# A) Input/Output
Coregistered_Pairs_Report = r"C:\Users\ggosseli\Desktop\QC03M_7262_4500\2_Coregistered_Scenes\04_Coregistered_Pairs_Report.txt"
output_folder = r"C:\Users\ggosseli\Desktop\QC03M_7262_4500"
prefix = "QC03M_7262_4500"

# B) Filtering options / Filtering is mandatory for Compact Pol analysis
filter_type = "psboxcar"          # options are psboxcar or pspolfil
filter_size_X_Y = [9,9]           # PSPOLFIL is always square, only the X value will be used. 

# C) Time-series creation. Compact Pol discriminators to stack. Look at the Part 2: Notes section for more information. 
CPDIS_to_stack = [10, 11]
produce_CPDIS_RGB = "yes"           # valid options are "yes" or "no"

apply_masking = "no"       # Valid option are "yes" or "no"
mask_type = "exclusion"    # Valid option are "inclusion" or "exclusion"
mask_file = r"\\W-BSC-A157283\share_GG\CanUS_border_1m_UTM18T_D000_v4_2km_buffer.pix"
mask_seg_number = [2]

# D) Orthorectification options
# Ortho bounds options: 1 (from an AOI file)  or 2 (from the input file)
DEM_file = r"C:\Users\ggosseli\Desktop\DEM\Glo30DEM_CanUS_LatLong.tif"
DEM_elevation_channel = 1

ortho_bounds_option = 1
AOI_vector_file = r"C:\Users\ggosseli\Desktop\QC03M_7262_4500\AOI_stack_QC03M_7262_4500_ASC_UTM18TD000.pix"
AOI_segment_number = 2

ortho_resolution_X = "2"
ortho_resolution_Y = "2"

# E) General options
generate_overviews = "yes"         #either "yes" or "no"
delete_intermediary_files = "no"   #either "yes" or "no"

# Behaviour when output file exists 
if_file_exists = "skip"     # Valid options are "skip" or "regenerate"

# -----------------------------------------------------------------------------------------------------------------
#  Part 2: Notes
# -----------------------------------------------------------------------------------------------------------------
'''
1) The compact polarimetry analysis is based on the Catalyst PSCOMDIS algorithm. 
   The output consist of 11 CP discriminator (32R channels): 
   [1]  Degree of polarization 
   [2]  Degrees of circular polarization 
   [3]  Degrees of linear polarization 
   [4]  Circular polarization ratio
   [5]  linear polarization ratio
   [6]  Orientation 
   [7]  Ellipticity 
   [8]  Relative phase
   [9]  Coherence 
   [10] Entropy 
   [11] Alpha angle 

   2) Behaviour when output file exists
       Warning - if_file_exists = "skip"  is a usefull feature for debugging or to (re) generate a subsequent product 
            in the processing chain without having to regenerate good files, for example when a stack is to be
            regenerated with an exclusion mask or a different filter size. However, missing (or deleted) intermediary 
            files when the final product exists will cause the final product to be outdated regarding the newest 
            regenerated files. 
       
                 if_file_exists = "regenerate" is a more secure option. 
   
    3) Overviews won't be generated for intermediary files when delete_intermediary_files = "no"    
    
 '''
# -----------------------------------------------------------------------------------------------------------------
#  Part 3: Imports
# -----------------------------------------------------------------------------------------------------------------
import sys
import os
import fnmatch
import time
import shutil
import locale
import re
from tracemalloc import start

from numpy import stack
import pci
from pci.psboxcar import psboxcar
from pci.pssartsa import pssartsa
from pci.iia import iia
from pci.ortho import ortho
from pci.pspolfil import pspolfil
from pci.pscomdis import pscomdis
from pci.pyramid import pyramid
from pci.exceptions import PCIException
from pci.api import datasource as ds
from pci.api.cts import crs_to_mapunits

from TSA_utilities.SAR_TSA_utilities_definitions import ortho_run
from TSA_utilities.SAR_TSA_utilities_definitions import get_folder_proctime_and_size
from TSA_utilities.SAR_TSA_utilities_definitions import file_size_check
from TSA_utilities.SAR_TSA_utilities_definitions import stack_masking

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -----------------------------------------------------------------------------------------------------------------
#  Part 4: Parameters validation
# -----------------------------------------------------------------------------------------------------------------
TSA_03_start = time.time()
val_start = time.time()

# Hardcoded parameters
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list + no_validation_list
GB = 1073741824
no_data_value = -32768.0000

AOI_file = AOI_vector_file
AOI_file_segment_number = AOI_segment_number

#  Versions control
print("\t")
print(pci.version)

print("Installed python version: " + sys.version)
py1 = str(sys.version_info[0])
py2 = str(sys.version_info[1])
py3 = (py1 + "." + py2)
python_version = float(py3)
if python_version < 2.0:
    print("You are using Python v" + str(python_version))
    print("You need to update to Python 3.8 or newer versions")
    sys.exit()
print("\t")

# -----------------------------------------------------------------------------------------------------------------
# A) Input files and folder
if not os.path.exists(Coregistered_Pairs_Report):
    print ("Error - The SARINGESTAOI_Baselines_Report does not exist or the path/filename is wrong.")
    sys.exit()

# -----------------------------------------------------------------------------------------------------------------
# B) Filtering options
if filter_type.lower() not in ["psboxcar", "pspolfil"]:
    print ('Error - The filter type option must be "psboxcar" or "pspolfil"')
    sys.exit()

filter_size_X = filter_size_X_Y[0] 
filter_size_Y = filter_size_X_Y[1] 
    
# Check for odd and positive numbers
'''
        if filter_size < 5:
            print("Error - The Filter window size must be an odd integer greater or equal to 5")
            print("Specified value:" + str(filter_size))
            sys.exit()
        FWsize_modulo = filter_size % 2
        if FWsize_modulo != 1:
            print("Error - The Filter window size must be an odd integer")
            print("Specified value:" + str(filter_size))
            sys.exit()
    else:
        apply_insraw_filter = False

'''
if filter_type == "pspolfil": 
    if filter_size_X < 5: 
        print ("Error - The PSPOLFIL filter size (X) must be equal or greater than 5")
        sys.exit()

if filter_type == "psboxcar": 
    if (filter_size_X == 1) and (filter_size_Y == 1): 
        print ("Error - At least one of the PSBOXCAR filter size must be greater than 1")
        sys.exit()

# B.3) Apply an exclusion or an inclusion mask
if apply_masking.lower() in yes_validation_list:
    apply_masking = True
    mask_type = mask_type.lower()
    if mask_type not in ["exclusion", "inclusion"]:
        print("Error - the mask_type option is invalid")
        print('Valid options are "exclusion" or "inclusion"')
        sys.exit()
    if not os.path.exists(mask_file):
        print("Error - The input mask_file does not exist or the path/filename is invalid.")
        sys.exit()
else :
    apply_masking = False
# -----------------------------------------------------------------------------------------------------------------
# C) Compact Pol discriminators to stack. Look at the Part 2: Notes section for more information. 

CP_chans = [1,2,3,4,5,6,7,8,9,10,11]
for item in CPDIS_to_stack: 
    if item not in CP_chans: 
        print ("Error - The following channel in CPDIS_to_stack is invalid")
        print ("wrong value: " + str(item))
        sys.exit()

CP_chans_labels = ["degpol","degcirpol","deglinpol","cpr","lpr","orientation","ellipticity","relphase","coherency","entropy","alpha"]
produce_CPDIS_RGB = produce_CPDIS_RGB.lower()
if produce_CPDIS_RGB in yes_validation_list: 
    produce_CPDIS_RGB = True
elif produce_CPDIS_RGB in no_validation_list: 
    produce_CPDIS_RGB = False
else: 
    print ('Error - The produce_CPDIS_RGB parameter valid options are "yes" or "no"')
    sys.exit()

# -----------------------------------------------------------------------------------------------------------------
# D)  Orthorectification options
if ortho_bounds_option == 1:
    print ("The output ortho bounds will be taken from the specified AOI file")
    if not os.path.exists(AOI_vector_file):
        print("Error - The AOI file does not exist or the path/filename is wrong")
        sys.exit()
elif ortho_bounds_option == 2:
    print("The output ortho bounds will be taken from each input images")
else:
    print("Error - The ortho bounds option is not valid")
    print("Valid options are 1 or 2")
    sys.exit()

if isinstance(ortho_resolution_X, str) is False:
    print("Warning - The ortho_resolution_X must be a string")
    print("Converting to a string")
    ortho_resolution_X_temp = str(ortho_resolution_X)
    ortho_resolution_X = ortho_resolution_X_temp

if isinstance(ortho_resolution_Y, str) is False:
    print("Warning - The ortho_resolution_Y must be a string")
    print("Converting to a string")
    ortho_resolution_Y_temp = str(ortho_resolution_Y)
    ortho_resolution_Y = ortho_resolution_Y_temp
    sys.exit()

if not os.path.exists(DEM_file):
    print ("Error - The DEM_file does not exist or the path/filename is wrong.")
    sys.exit()

# -----------------------------------------------------------------------------------------------------------------
# E) General options
#    Check for scene size conformity (size and matrix type). 
if if_file_exists.lower() not in ["skip","regenerate"]:     
    print('Error - valid options for existing_data are "skip" or "regenerate"')
    sys.exit()
else: 
    info_message_skip = ("   output file already exists - skip (if_file_exists = skip)")
    info_message_regn =  ("   output file already exists - regenerate (if_file_exists = regenerate)")
 
Fld_Coregistration = os.path.dirname(Coregistered_Pairs_Report)
Coregistered_Pairs_Report

input_scenes_list = []
for root, dirs, files in os.walk(Fld_Coregistration):
    for filename in fnmatch.filter(files, "*.pix"):
        input_scenes_list.append(os.path.join(root, filename))

print ("\t")
print(time.strftime("%H:%M:%S") + " Checking the scenes dimensions (lines, columns)")
files_list = input_scenes_list
file_size_check (files_list)

print ("\t")
print(time.strftime("%H:%M:%S") + " Checking the scenes matrix type")
for ii in input_scenes_list: 
    print ("Checking the input scenes conformity--->" + ii)
    with ds.open_dataset(ii, ds.eAM_READ) as ds2:
        aux = ds2.aux_data
        Matrix_Type = aux.get_file_metadata_value("Matrix_Type")
        num_channels = ds2.chan_count
        print (Matrix_Type)

        if Matrix_Type != "s2c":
            print ("Errror - The input scene matrix type must be s2c for the compact pol analysis")
            sys.exit()
        if num_channels != 2: 
            print ("Errror - The input scene mush have exactly 2 channels")
            sys.exit()

generate_overviews = generate_overviews.lower()
if generate_overviews in yes_validation_list: 
    generate_overviews = True
elif generate_overviews in no_validation_list: 
    generate_overviews = False
else: 
    print ('Error - The generate_overviews parameter must be set with "yes" or "no"')
    sys.exit()

delete_intermediary_files = delete_intermediary_files.lower()
if delete_intermediary_files in yes_validation_list: 
    delete_intermediary_files = True
elif delete_intermediary_files in no_validation_list: 
    delete_intermediary_files = False
else: 
    print ('Error - The delete_intermediary_files parameter must be set with "yes" or "no"')
    sys.exit()

# Creating the ouput folder, we will add automatically a subfolder
Fld_CPOL = os.path.join (output_folder, "5_1_1_compactpol")
if not os.path.exists(Fld_CPOL):
    os.makedirs(Fld_CPOL)

Fld_CPOL_ortho = os.path.join (output_folder, "5_1_2_compactpol_ortho")
if not os.path.exists(Fld_CPOL_ortho):
    os.makedirs(Fld_CPOL_ortho)

# All validations have suceeded, a time log file is open.
script3_procTime = os.path.join(output_folder, prefix + "TSA_script3_CPOL_processingTime.txt")
time_log = open(script3_procTime, "w")

current_time =  time.localtime()
string_0 = time.strftime("%Y-%m-%d", current_time)
time_log.write("%s\n" % string_0)

string_0 = ("Process;proc.time (secs);Data size (MB);Number of files")
time_log.write("%s\n" % string_0)

val_stop = time.time()
string_1 = "Validation: ;" + str (round((val_stop - val_start), 2)) 
time_log.write("%s\n" % string_1)

# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
#  Part 5: Main program
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
proc_start_time = time.time()
print("\t")
print("----------------------------------------------------------------------------------------------------------")
print("                              Input CP scenes filtering  and CP parameters generation                     ")
print("----------------------------------------------------------------------------------------------------------")
print("\t")

outfiles_pscomdis = []
nb_files = str(len(input_scenes_list))
count = 1
for ii in input_scenes_list: 

    with ds.open_dataset(ii, ds.eAM_READ) as ds5:
        aux = ds5.aux_data
        Acquisition_DateTime = aux.get_file_metadata_value("Acquisition_DateTime")

    print("\t")
    print(time.strftime("%H:%M:%S") + "  Processing file " + str(count) + " of " + nb_files)
    print (" input file--> " + ii)
    print ("   filtering the input file")
    fili = ii

    # Remove the prefix  coreg* for cleaner ouput names
    base_out = os.path.basename(ii[:-4])
    index = base_out.find(prefix)
    if index == -1:
        print ("   prefix not found within the file name, will use the full name instead")
    else:
        base_out = base_out[index:]

    # Applying a filter 
    if filter_type == "psboxcar":  
        
        out_f = ("_psboxcar_" + str(filter_size_X) +"x"+ str(filter_size_Y))
        filo = os.path.join (Fld_CPOL, base_out + out_f + ".pix")
        flsz = filter_size_X_Y

        print ("   output file--> " + filo)
        if os.path.exists (filo) and if_file_exists == "skip": 
            print (info_message_skip)
        else: 
            if os.path.exists (filo) and if_file_exists == "regenerate": 
                print (info_message_regn)
                os.remove (filo)
            
            print ("   Filtering the file")
            try: 
                psboxcar (fili, filo, flsz)
            except PCIException as e:
                    print(e)
            except Exception as e:
                print(e)   

    if filter_type == "pspolfil": 
        out_f = ("_pspolfil_" + str(filter_size_X) +"x"+ str(filter_size_X))
        filo = os.path.join (Fld_CPOL, base_out + out_f + ".pix")
        flsz = filter_size_X_Y[0]
        nlook = [1]

        print (" output file--> " + filo)
        if os.path.exists (filo) and if_file_exists == "skip": 
            print (info_message_skip)
        else: 
            if os.path.exists (filo) and if_file_exists == "regenerate": 
                print (info_message_regn)
                os.remove (filo)

            print ("   Filtering the file")
            try: 
                pspolfil (fili, filo, flsz, nlook)
            except PCIException as e:
                 print(e)
            except Exception as e:
                print(e)  
    iii_input = filo

    # Use PSCOMDIS to generate the compact pol parameters
    print ("   " + time.strftime("%H:%M:%S") + " Generating the Compact Pol parameters")
    fili = filo
    base_out = os.path.basename(filo[:-4])
    filo = os.path.join (Fld_CPOL, base_out + "_pscomdis.pix")
    angletyp = "Degrees"
    
    print ("   output file--> " + filo)
    if os.path.exists (filo) and if_file_exists == "skip": 
        print (info_message_skip)
    else: 
        if os.path.exists (filo) and if_file_exists == "regenerate": 
            print (info_message_regn)
            os.remove (filo)
        try: 
           pscomdis (fili, filo, angletyp)
        except PCIException as e:
                print(e)
        except Exception as e:
            print(e)  
    
    outfiles_pscomdis.append(filo)
    iii_output = filo
    # Transfer the math model from the filtered file to the PSCOMDIS file
    fili = iii_input
    filo = iii_output
    dbsl = [2]         # Potential bug, assumes the math model is segment 2
    dbos = []

    try: 
        iia (fili, filo, dbsl, dbos)
    except PCIException as e:
        print(e)
    except Exception as e:
        print(e)  

    with ds.open_dataset(iii_output, ds.eAM_WRITE) as ds6:
        aux = ds6.aux_data
        metadata = aux.file_metadata
        metadata["Acquisition_DateTime"] = Acquisition_DateTime
        aux.file_metadata = metadata
        ds6.aux_data = aux

    count = count + 1

proc_stop_time = time.time()
folder = Fld_CPOL
out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
string_1 = ("CompactPol parameters, files filtering and parameters generation: " + out_folder_time_size) 
time_log.write("%s\n" % string_1)

print("\t")
print("-------------------------------------------------------------------------------------------------------")
print("                               Compact pol discriminators orthorectification                           ")
print("-------------------------------------------------------------------------------------------------------")
print("\t")
proc_start_time = time.time()
if ortho_bounds_option == 1:   # extends from the AOI file

    print(time.strftime("%H:%M:%S") + " Extracting the bounding box coordinate around the input AOI")
    print("   Input AOI: " + AOI_file)
    # X_coordinates = []
    # Y_coordinates = []
    # open the dataset in write mode
    with ds.open_dataset(AOI_file) as ds1:
        AOI_MapProjection = crs_to_mapunits(ds1.crs)

        print("\t")
        print("   AOI Map projection: " + AOI_MapProjection)
        ## Could DO A VALIDATION OF THE PROJECTION HERE,

        # get vector segment
        io = ds1.get_vector_io(AOI_file_segment_number)
        xs = []
        ys = []
        # iterate over shapes in the segment
        for index, shape in zip(io.shape_ids, io):
            # iterate over rings in the shape:
            xs.append(shape.extents[0])
            xs.append(shape.extents[2])
            ys.append(shape.extents[1])
            ys.append(shape.extents[3])

    # Upper left and lower right coordinates of the bounding box to
    # calculate the DBIW window for every image
    AOI_UL_X = min(xs)
    AOI_UL_Y = max(ys)
    AOI_LR_X = max(xs)
    AOI_LR_Y = min(ys)

    print("AOI Extends")
    print(AOI_UL_X, AOI_UL_Y, AOI_LR_X, AOI_LR_Y)
    print("\t")

    ulx_AOI = str(AOI_UL_X)
    uly_AOI = str(AOI_UL_Y)
    lrx_AOI = str(AOI_LR_X)
    lry_AOI = str(AOI_LR_Y)

# ----------------------------------------------------------------------------------------------------------
files_to_ortho_list = outfiles_pscomdis
CP_ortho_list = []
number_of_files = str(len(files_to_ortho_list))
print("\t")
count = 1
for input_scene in files_to_ortho_list:

    print(((time.strftime("%H:%M:%S")) + " Orthorectifying file " + str(count) + " of " + number_of_files))
    print("   Input file: " + input_scene)

    mfile = input_scene
    dbic = [1,2,3,4,5,6,7,8,9,10,11]
    mmseg = []
    dbiw = []
    srcbgd = "NONE"

    base = os.path.basename(input_scene)
    filo = os.path.join(Fld_CPOL_ortho, "o" + base)

    ftype = "PIX"
    foptions = "TILED256"
    outbgd = [-32768.00000]
    edgeclip = []

    if ortho_bounds_option == 1:  # extents from the AOI scenes
        mapunits = AOI_MapProjection
        ulx = ulx_AOI
        uly = uly_AOI
        lrx = lrx_AOI
        lry = lry_AOI
    if ortho_bounds_option == 2:  # extents from the input scenes
        mapunits = ""
        ulx = ""
        uly = ""
        lrx = ""
        lry = ""

    tipostrn = ""
    bxpxsz = ortho_resolution_X
    bypxsz = ortho_resolution_Y
    filedem = DEM_file
    dbec = [DEM_elevation_channel]
    backelev = []
    elevref = ""
    elevunit = ""
    elfactor = []
    proc = ""
    sampling = [1]
    resample = "near"

    print("   output orthorectified file: " + filo)
    if os.path.exists (filo) and if_file_exists == "skip": 
        print (info_message_skip)
    else: 
        if os.path.exists (filo) and if_file_exists == "regenerate": 
            print (info_message_regn)
            os.remove (filo)
        try:
            ortho(mfile, dbic, mmseg, dbiw, srcbgd, filo, ftype, foptions, outbgd,
                    ulx, uly, lrx, lry, edgeclip, tipostrn, mapunits, bxpxsz, bypxsz, filedem,
                    dbec, backelev, elevref, elevunit, elfactor, proc, sampling, resample)
            if generate_overviews is True and delete_intermediary_files is False:
                print ("   Generating file overviews")
                pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")

        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)  
    CP_ortho_list.append(filo)
    count = count + 1

proc_stop_time = time.time()
folder = Fld_CPOL_ortho
out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
string_1 = ("CompactPol scenes orthorectification: " + out_folder_time_size) 
time_log.write("%s\n" % string_1)


print("\t")
print("-------------------------------------------------------------------------------------------------------")
print("                       Time series creation for selected Compact Pol parameters                        ")
print("-------------------------------------------------------------------------------------------------------")
print("\t")
proc_stop_time = time.time()

if apply_masking is True:
    print("Stack data preprocessing - Applying the exclusion or inclusion vector mask")
    input_stack = CP_ortho_list
    stack_masking(input_stack, mask_type, mask_file, mask_seg_number,no_data_value, output_folder)
if apply_masking is False :
    print(time.strftime("%H:%M:%S") + " Stack preparation - Exclusion or inclusion mask not requested")

# Creating the ouput folder
Fld_TSA_stacks = os.path.join (output_folder, "7_TSA_stacks")
if not os.path.exists(Fld_TSA_stacks):
    os.makedirs(Fld_TSA_stacks)

# Create a temporary mfile
temp_mfile = os.path.join(Fld_TSA_stacks, "temp_mfile_for_CP_data_ingestion.txt")
mfile_input = open(temp_mfile, "w")
mfile_input.write('\n'.join(CP_ortho_list))
mfile_input.close()

nb_params = str(len(CPDIS_to_stack))
count = 1
for in_chan in CPDIS_to_stack:
    print ("\t")
    print(time.strftime("%H:%M:%S") + " processing CP parameter " + str (count) + " of " + nb_params)
    lab_index = in_chan - 1
    label_out = CP_chans_labels[lab_index]

    mfile = temp_mfile
    dbic = [in_chan]
    mask = []
    maskfile = "" 
    stack  = "yes"
    flsz = []
    filo = os.path.join (Fld_TSA_stacks, prefix + "_" + label_out + ".pix")

    print("   output stack-->" + filo)
    if os.path.exists (filo) and if_file_exists == "skip": 
        print (info_message_skip)
    else: 
        if os.path.exists (filo) and if_file_exists == "regenerate": 
            print (info_message_regn)
            os.remove (filo)
    try: 
        pssartsa(mfile, dbic, mask, maskfile, stack, flsz, filo)
        if generate_overviews is True:
            pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")
    except PCIException as e:
        print(e)
    except Exception as e:
        print(e)
    count = count + 1

os.remove(temp_mfile)

proc_stop_time = time.time()
folder = Fld_TSA_stacks
out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
string_1 = ("CompactPol stack generation: " + out_folder_time_size) 
time_log.write("%s\n" % string_1)

# -------------------------------------------------------------------------------------------------------------
# J) Deleting the intermediary files if requested
# -------------------------------------------------------------------------------------------------------------

if delete_intermediary_files is True:
    proc_start_time = time.time()
    print("\t")
    print("-------------------------------------------------------------------------------------------------------")
    print("                                Deleting the intermediary files                                        ")
    print("-------------------------------------------------------------------------------------------------------")
    print("\t")
    
    string1 = "Deleting intermediary files requested"
    time_log.write("%s\n" % string1)

    del_folders_list = []
    # List of all (possible) folders to delete
    del_folders_list.append(os.path.join(output_folder, "5_1_2_compactpol_ortho"))
    del_folders_list.append(os.path.join(output_folder, "5_1_1_compactpol"))

    for delete in del_folders_list:
        if os.path.isdir(delete):
            print(time.strftime("%H:%M:%S") + " Deleting " + delete)
            shutil.rmtree(delete)


print("\t")
print("-------------------------------------------------------------------------------------------------------")
print((time.strftime("%H:%M:%S")))
print("All processing completed")
print("\t")

TSA_03_stop = time.time()
ellapse_time_seconds = round((TSA_03_stop - TSA_03_start), 2)
ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

print("Processing time (seconds): " + str(ellapse_time_seconds))
print("Processing time (minutes): " + str(ellapse_time_minutes))
print("Processing time (hours): " + str(ellapse_time_hours))

string1 = "TSA_03 total processing time (secs):;" + str(ellapse_time_seconds)
time_log.write("%s\n" % string1)
string1 = "TSA_03 total processing time (mins):;" + str(ellapse_time_minutes)
time_log.write("%s\n" % string1)
string1 = "TSA_03 total processing time hours):;" + str(ellapse_time_hours)
time_log.write("%s\n" % string1)
time_log.close()
