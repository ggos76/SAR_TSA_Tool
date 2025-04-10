#!/usr/bin/env python
'''-----------------------------------------------------------------------------------------------------------------
 * Gabriel Gosselin, CCMEO\ CCRS. 2024-2025                                                                                -
 * -----------------------------------------------------------------------------------------------------------------
'''

# -----------------------------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# -----------------------------------------------------------------------------------------------------------------

# A) Input/Output
Coregistered_Pairs_Report = r"E:\RCMP_RCM\stack01_3MCP34_DESC\2_Coregistered_Scenes\04_Coregistered_Pairs_Report.txt"
output_folder = r"E:\RCMP_RCM\stack01_3MCP34_DESC"
prefix = "stk01_"

# B) Elevation source
DEM_file = r"D:\RCMP_border\aux_DEM\Glo30DEM_CanUS_LatLong.tif"
DEM_elevation_channel = 1

# C)  The channels to process for the time series (TSA) generation, can be a subset of the coregistered files channels
# C.1) These parameters must be common to all scenes to be processed. The TSA_channel_mapping list must include all
#      channels of the input file and not only the channels to be processed.
TSA_channel_mapping = [1,2]
TSA_channel_labels = ["RH", "RV"]

# C.2) Coherence  options
produce_coherence_layers = "yes"    # Only available for SLC data
TSA_coherence_layers = [1]      # Must be all subset of TSA_channel_mapping
apply_insraw_filter = "yes"
filter_size = 9

# C.3) Intensity options
produce_intensity_layers = "yes"
TSA_intensity_layers = [1, 2]    # Must be all subset of TSA_channel_mapping

# C.4) Other intensity layers
produce_incidence_angle_layer = "yes"

# C.5) Other math layers from the intensity files
produce_math_layers = "yes"               
TSA_math_models_definition = r"C:\Users\ggosseli\source\repos\SAR_TSA_Tool\TSA_utilities\model_00_band_ratio.txt"
TSA_math_xtra_channels = 2
TSA_math_xtra_labels = ["RH div RV", "math_2"]

# D) Orthorectification options
# Ortho bounds options: 1 (from an AOI file)  or 2 (from the input file)
ortho_bounds_option = 1
AOI_vector_file = r"E:\RCMP_RCM\stack01_3MCP34_DESC_aoi.pix"
AOI_segment_number = 2

ortho_resolution_X = "2"
ortho_resolution_Y = "2"

# E) General options
# Generate overviews - either yes or no,
generate_overviews = "yes"
# keep or delete intermediate files - either yes or no. No is recommended.
delete_intermediary_files = "no"

# -----------------------------------------------------------------------------------------------------------------
#  Part 2: Notes
# -----------------------------------------------------------------------------------------------------------------
'''
1) There is no verification if the specified DEM covers completely the spatial
   extents of the AOI. If not the script will run to completion but all
   interferograms will be blank after INSRAW.
'''
# -----------------------------------------------------------------------------------------------------------------
#  Part 3: Imports
# -----------------------------------------------------------------------------------------------------------------
import sys
import os
import fnmatch
import time
import locale
import re
from tracemalloc import start

import pci
from pci.insraw import insraw
from pci.pyramid import pyramid
from pci.exceptions import PCIException
from pci.api import datasource as ds

from TSA_utilities.SAR_TSA_utilities_definitions import ortho_run
from TSA_utilities.SAR_TSA_utilities_definitions import nan_replace
from TSA_utilities.SAR_TSA_utilities_definitions import psiqinterp_run
from TSA_utilities.SAR_TSA_utilities_definitions import create_list
from TSA_utilities.SAR_TSA_utilities_definitions import get_folder_proctime_and_size
from TSA_utilities.SAR_TSA_other_layers import inc_angle_layer
from TSA_utilities.SAR_TSA_other_layers import math_layers_prod
from TSA_utilities.SAR_TSA_other_layers import math_layers_split

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

Fld_Coregistration = os.path.join(output_folder, "2_Coregistered_Scenes")
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
if python_version < 3.8:
    print("You are using Python v" + str(python_version))
    print("You need to update to Python 3.8 or newer versions")
    sys.exit()
print("\t")

# -----------------------------------------------------------------------------------------------------------------
# A) Input files and folder
if not os.path.exists(Coregistered_Pairs_Report):
    print ("Error - The SARINGESTAOI_Baselines_Report does not exist or the path/filename is wrong.")
    sys.exit()
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

if not os.path.exists(Fld_Coregistration):
    print ("Error - The coregistered scenes folder does not exist")
    print (" Expected path: " + Fld_Coregistration)
    sys.exit()
# -----------------------------------------------------------------------------------------------------------------
# B) Elevation source
if not os.path.exists(DEM_file):
    print ("Error - The DEM_file does not exist or the path/filename is wrong.")
    sys.exit()

# -----------------------------------------------------------------------------------------------------------------
# C.0) Channel mapping.  We check for the input scenes conformity. All files must have the same number of 
#     lines and colums and the same matrix type. Matrix type must be complex if the coherence option is 
#     selected.  
produce_coherence_layers = produce_coherence_layers.lower()
if produce_coherence_layers in yes_validation_list:
    produce_coherence_layers = True
else:
    produce_coherence_layers = False

produce_intensity_layers = produce_intensity_layers.lower()
if produce_intensity_layers in yes_validation_list:
    produce_intensity_layers = True
else:
    produce_intensity_layers = False

produce_incidence_angle_layer = produce_incidence_angle_layer.lower()
if produce_incidence_angle_layer in yes_validation_list:
    produce_incidence_angle_layer = True
else:
    produce_incidence_angle_layer = False

produce_math_layers = produce_math_layers.lower()
if produce_math_layers in yes_validation_list:
    produce_math_layers = True
else:
    produce_math_layers = False

# Check if a least one layer type is selected:
if (produce_coherence_layers is False) and (produce_intensity_layers is False) :
    print ("Error - At least one type of output layer must be selected")
    sys.exit()

if (produce_intensity_layers is False) and (produce_math_layers is True): 
    print ("Error - produce_intensity_layers must be set to yes when produce_math_layers is set to yes")
    sys.exit()

if (produce_intensity_layers is False) and (produce_incidence_angle_layer is True): 
    print ("Error - produce_intensity_layers must be set to yes when produce_incidence_angle_layer is set to yes")
    sys.exit()
# ------------------------------------------------------------------------------
# C.1) Channel mapping.  We check for the input scenes conformity
if len (TSA_channel_mapping) != len (TSA_channel_labels):
    print ("Error - The length of TSA_channel_mapping and the length of TSA_channel_labels are not the same")
    sys.exit()

check_scenes_list = []
for root, dirs, files in os.walk(Fld_Coregistration):
    for filename in fnmatch.filter(files, "*.pix"):
        check_scenes_list.append(os.path.join(root, filename))

for ii in check_scenes_list:
    print ("Checking the input scenes conformity--->" + ii)
    with ds.open_dataset(ii, ds.eAM_READ) as ds2:
        aux = ds2.aux_data
        Matrix_Type = aux.get_file_metadata_value("Matrix_Type")
        SensorModelName = aux.get_file_metadata_value("SensorModelName")
        num_channels = ds2.chan_count
        print (Matrix_Type)

        # check size here!
        # zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        # check size here!
        if len(TSA_channel_mapping) != num_channels:
            print ("Error 11 - Mismatch between the length of TSA_channel_mapping  and the number of input channels")
            sys.exit()

        if produce_coherence_layers is True:
            if Matrix_Type not in ["s1c", "s2c", "s4c"]:

                print("Error 20 - The input scene channels must be complex when the coherence layer option is selected")
                sys.exit()
            for in_chan in TSA_coherence_layers:
                if in_chan not in TSA_channel_mapping:
                    print("Error - TSA_coherence_layers input channel " + str(in_chan) + " is not included in"
                          " TSA_channel_mapping")
                    sys.exit()

        if produce_intensity_layers is True:
            for in_chan in TSA_intensity_layers:
                if in_chan not in TSA_channel_mapping:
                    print("Error - TSA_intensity_layers input channel " + str(in_chan) + " is not included in"
                          " TSA_channel_mapping")
                    sys.exit()
# ------------------------------------------------------------------------------
# C.2) Coherence options - Validation for the Raw interferograms filter options
if produce_coherence_layers is True:
    apply_insraw_filter = apply_insraw_filter.lower()
    if apply_insraw_filter in yes_validation_list:
        apply_insraw_filter = True

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

    # The basic validation is done, we can create the empty folders
    Fld_Raw_Interferograms = os.path.join(output_folder, "3_1_1_RAW_interferograms")
    if not os.path.exists(Fld_Raw_Interferograms):
        os.makedirs(Fld_Raw_Interferograms)
    Fld_coherence = os.path.join(output_folder, "3_1_2_Coherence")
    if not os.path.exists(Fld_coherence):
        os.makedirs(Fld_coherence)
    Fld_coherence_ortho = os.path.join(output_folder, "3_1_3_Coherence_ortho")
    if not os.path.exists(Fld_coherence_ortho):
        os.makedirs(Fld_coherence_ortho)

if produce_intensity_layers is True:
    Fld_intensity = os.path.join (output_folder, "3_2_1_Intensity")
    if not os.path.exists(Fld_intensity):
        os.makedirs(Fld_intensity)
    Fld_intensity_ortho = os.path.join (output_folder, "3_2_2_Intensity_ortho")
    if not os.path.exists(Fld_intensity_ortho):
        os.makedirs(Fld_intensity_ortho)

#-------------------------------------------------------------------------------------------
# C.3) Parameters validation for the incidence angle and the other intensity layers options
# Check if there is an incidence angle segment for all input scenes.    
if produce_incidence_angle_layer is True:
    print("\t")
    print(time.strftime("%H:%M:%S") + " Checking for incidence angle arrays...")
    print("\t")
    for ii in check_scenes_list: 
        print ("Checking -->" + ii )
        with ds.open_dataset(ii, ds.eAM_READ) as ds3:
            arr_id = ds3.get_array_io_ids()  # returns a list of the incidence angle array segment numbers. 
            for jj in arr_id: 
                print ("   incidence angle array id -->"+ str(jj))
            ss2 = len(str(arr_id))
            if ss2 == 0:
                print ("Error - This file does not contain an incidence angle arrray")
                sys.exit()

# -----------------------------------------------------------------------------------------------
# C.5 Math Layers
if produce_math_layers is True:
    #A) check if the intensity layers are selected. 
    
    if len(TSA_math_xtra_labels) != TSA_math_xtra_channels: 
        print ("Error - The number of  TSA_math_xtra_channels must be equal to the number of TSA_math_xtra_labels")
        sys.exit()
    if not os.path.exists (TSA_math_models_definition): 
        print ("Error - The TSA_math_models_definition  file does't exist or the path/file is wrong")
        sys.exit()
    invalid_chars = r'[<>:"/\\|?*]'
    # Check for illegal characters in windows OS, this is somehow expected for math channels. 
    for in_label in TSA_math_xtra_labels:
        has_invalid = bool(re.search(invalid_chars, in_label))
        if has_invalid is True:
            print ("\t")
            print ("Error - The following TSA_math_xtra_labels contains an invalid character")
            print ('List of invalid characters: < > : " /  \\ | ? * "')
            print (in_label)
            sys.exit()
    Fld_math_layers = os.path.join (output_folder, "3_3_1_Math_layers")
    if not os.path.exists(Fld_math_layers):
        os.makedirs(Fld_math_layers)
    Fld_math_layers_ortho = os.path.join (output_folder, "3_3_2_Math_layers_ortho")
    if not os.path.exists(Fld_math_layers_ortho):
        os.makedirs(Fld_math_layers_ortho)
    Fld_math_layers_ortho_split = os.path.join (output_folder, "3_3_3_Math_layers_ortho")
    if not os.path.exists(Fld_math_layers_ortho_split):
        os.makedirs(Fld_math_layers_ortho_split)

#-------------------------------------------------------------------------------------------
# D) Orthorectification options
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

# F) Other options verification
generate_overviews = generate_overviews.lower()
if generate_overviews in yes_validation_list: 
    generate_overviews = True
elif generate_overviews in no_validation_list:
    generate_overviews = False
else: 
    print("\t")
    print("Error - Selected overviews option is not valid.")
    print('Accepted values are: "yes" or "no"')
    sys.exit()

delete_intermediary_files = delete_intermediary_files.lower()
if delete_intermediary_files in yes_validation_list: 
    delete_intermediary_files = True
elif delete_intermediary_files in no_validation_list:
    delete_intermediary_files = False
else: 
    print("\t")
    print("Error - Selected delete_intermediary_files option is not valid.")
    print('Accepted values are: "yes" or "no"')
    sys.exit()

Fld_Output_stack_lists = os.path.join(output_folder, "6_Output_stack_lists")
if not os.path.exists(Fld_Output_stack_lists):
    os.makedirs(Fld_Output_stack_lists)

# All validations have suceeded, a time log file is open.
script3_procTime = os.path.join(output_folder, prefix + "TSA_script3_processingTime.txt")
time_log = open(script3_procTime, "w")
string_0 = ("Process;proc.time (secs);Data size (MB);Number of files")
time_log.write("%s\n" % string_0)

val_stop = time.time()
string_1 = "Validation step:" + str (round((val_stop - val_start), 2)) 
time_log.write("%s\n" % string_1)

# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
#  Part 5: Main program
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------

if produce_coherence_layers is True:
    proc_start_time = time.time()
    print("\t")
    print("-------------------------------------------------------------------------------------------------------")
    print("                               RAW INTERFEROGRAMS GENERATION                                           ")
    print("-------------------------------------------------------------------------------------------------------")
    print("\t")

    Raw_Interferogram_list = []
    Input_file_lines = []
    pair_number = []
    Reference_file = []
    Reference_acquisition_date = []
    Dependent_file = []
    Dependent_acquisition_date = []
    print("\t")
    print(time.strftime("%H:%M:%S") + " Reading and parsing the Coregistered_Pairs_Report file")

    with open(Coregistered_Pairs_Report, "r") as input_file:
        for line in input_file:
            Input_file_lines.append(line.strip())
        for ii in Input_file_lines:
            temp = []
            temp = ii.split(';')
            pair_number.append(temp[0])
            Reference_file.append(temp[1])
            Reference_acquisition_date.append(temp[2])
            Dependent_file.append(temp[3])
            Dependent_acquisition_date.append(temp[4])

    print(time.strftime("%H:%M:%S") + " Genereating the raw interferograms")
    print ("\t")

    batch = 1
    batch_no = str(len(TSA_coherence_layers))
    for pol_chan  in TSA_coherence_layers:
        pol_labels = TSA_channel_labels [pol_chan -1]
        count = 1
        for ref_file, ref_date, dep_file, dep_date in \
            zip(Reference_file, Reference_acquisition_date, Dependent_file,
                Dependent_acquisition_date):

            rawinf_start = time.time()
            print(((time.strftime("%H:%M:%S")) + " RAW Interferogram generation, channel " + str (batch) +
                   " of " + batch_no + ", pair " + str(count) + " of " + str(len(Reference_file))))
            print("   Ref-->" + ref_file)
            print("   Dep-->" + dep_file)

            # FILREF=Reference
            filref = ref_file
            dbic_ref = [pol_chan]
            # fili=dependent
            fili = dep_file
            dbic = [pol_chan]

            if apply_insraw_filter is True:
                flsz = [filter_size]
                filo = os.path.join(Fld_Raw_Interferograms, "raw" +
                                    str(filter_size) + "_" + prefix + "ref" +
                                    ref_date + "_dep" + dep_date + "_" + pol_labels + ".pix")
            else:
                flsz = []
                filo = os.path.join(Fld_Raw_Interferograms, "raw_0_" + prefix +
                                    "ref" + ref_date + "_dep" + dep_date +"_" + pol_labels + ".pix")

            print("   Out-->" + filo)

            if not os.path.exists(filo):
                try:
                    insraw(filref, dbic_ref, fili, dbic, flsz, filo)
                    Raw_Interferogram_list.append(filo)

                    if generate_overviews is True:
                        pyramid(file=filo, force='yes', poption='aver',
                                dboc=[], olevels=[])
                except PCIException as e:
                        print(e)
                except Exception as e:
                    print(e)
            else:
                print ("Output file already exists - skip")
            print("\t")
            count = count + 1
        batch = batch + 1

    file = os.path.join(Fld_Raw_Interferograms, "MFILE_Insraw_Inteferograms_list.txt")
    with open(file, "w") as f:
        f.write("\n".join(Raw_Interferogram_list))

    proc_stop_time = time.time()
    folder = Fld_Raw_Interferograms
    out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
    string_1 = ("RAW interferogram generation: " + out_folder_time_size) 
    time_log.write("%s\n" % string_1)
    
    # ---------------------------------------------------------------------------------------------------
    # Replacings NANs (if any) by NoDATA in the output RAW interferograms
    proc_stop_time = time.time()
    step_nans = "INSRAW_"
    input_folder_nans = Fld_Raw_Interferograms
    nans_file = os.path.join (input_folder_nans, "INSRAW_replace_nans_info.txt")
    if os.path.exists(nans_file):
        print (" Checks for NANs already done - skip")
    else: 
        nan_replace(step_nans, input_folder_nans)
    
    proc_stop_time = time.time()
    folder = Fld_Raw_Interferograms
    out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
    string_1 = ("Checking for NANS:" + out_folder_time_size) 
    time_log.write("%s\n" % string_1)
   
    # -------------------------------------------------------------------------------------------------------------
    # Generating the coherence layers
    print("\t")
    print("-------------------------------------------------------------------------------------------------------")
    print("                              Generating the coherence layers                                          ")
    print("-------------------------------------------------------------------------------------------------------")
    print("\t")
    proc_start_time = time.time()
    print ("\t")
    print(time.strftime("%H:%M:%S") + " Converting the Raw Interferograms to phase coherence (0-1)")
    search_folder = Fld_Raw_Interferograms
    TSA_layers = TSA_coherence_layers
    TSA_labels = TSA_channel_labels
    unique_files = "no"
    keyword = "*.pix"
    interp_type = "amp"
    suffix = "_coh"
    ps_output_folder = Fld_coherence
    math_chans = "no"

    psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_layers,
                    TSA_labels, ps_output_folder, unique_files, prefix, math_chans)


    proc_stop_time = time.time()
    folder = Fld_coherence
    out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
    string_1 = ("Coherence layers generation:" + out_folder_time_size) 
    time_log.write("%s\n" % string_1)
    # --------------------------------------------------------------------------------------------------------------
    # Orthorectification of the coherences layers
    proc_start_time = time.time()
    print ("\t")
    print(time.strftime("%H:%M:%S") + " Orthorectification of the coherence layers ")

    input_folder_for_ortho = Fld_coherence
    output_folder_ortho = Fld_coherence_ortho

    ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel,
               ortho_bounds_option, AOI_file, AOI_file_segment_number, ortho_resolution_X,
               ortho_resolution_Y, generate_overviews, TSA_math_xtra_channels)

    proc_stop_time = time.time()
    folder = Fld_coherence_ortho
    out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
    string_1 = ("Coherence layers orthorectification:" + out_folder_time_size) 
    time_log.write("%s\n" % string_1)

if produce_intensity_layers is True:
    print("\t")
    print("-------------------------------------------------------------------------------------------------------")
    print("                              Generating the intensity layers                                          ")
    print("-------------------------------------------------------------------------------------------------------")
    print("\t")
    proc_start_time = time.time()

    print(time.strftime("%H:%M:%S") + " Creating the intensity layers")
    search_folder = Fld_Coregistration
    TSA_layers = TSA_intensity_layers
    TSA_labels = TSA_channel_labels
    unique_files = "yes"
    keyword = "*.pix"
    interp_type = "int"
    suffix = "int"
    ps_output_folder = Fld_intensity

    if produce_math_layers is True: 
        math_chans = "yes"
    else: 
        math_chans = "no"

    psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_layers,
                    TSA_labels, ps_output_folder, unique_files, prefix, math_chans)

    proc_stop_time = time.time()
    folder = Fld_intensity
    out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
    string_1 = ("Intensity layers generation:" + out_folder_time_size) 
    time_log.write("%s\n" % string_1)

    # --------------------------------------------------------------------------------
    if produce_incidence_angle_layer is True:
        proc_start_time = time.time()

        input_files_list = check_scenes_list
        output_folder_inc = Fld_intensity
        inc_angle_layer (input_files_list, output_folder_inc, prefix)
        
        proc_stop_time = time.time()
        folder = Fld_intensity
        out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
        string_1 = ("Incidence angle channel generation:" + out_folder_time_size) 
        time_log.write("%s\n" % string_1)


    # Orthorectification of the intensity layers
    proc_start_time = time.time()
    print(time.strftime("%H:%M:%S") + " Orthorectification of the intensity layers")
    input_folder_for_ortho = Fld_intensity
    output_folder_ortho = Fld_intensity_ortho
 
    ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel,
               ortho_bounds_option, AOI_file, AOI_file_segment_number, ortho_resolution_X,
               ortho_resolution_Y, generate_overviews, TSA_math_xtra_channels)
    
    proc_stop_time = time.time()
    folder = Fld_intensity_ortho
    out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
    string_1 = ("Intensity layers orthorectification:" + out_folder_time_size) 
    time_log.write("%s\n" % string_1)

    # --------------------------------------------------------------------------------
    if produce_math_layers is True:
        proc_start_time = time.time()

        math_layers_prod (TSA_math_models_definition, TSA_math_xtra_channels, TSA_math_xtra_labels,
                        Fld_math_layers, Fld_math_layers_ortho, prefix)
        
        proc_stop_time = time.time()
        folder = Fld_math_layers
        out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
        string_1 = ("Math layer(s) generation" + out_folder_time_size) 
        time_log.write("%s\n" % string_1)


        # Orthorectification of the Math layers
        proc_start_time = time.time()
        print("\t")
        print(time.strftime("%H:%M:%S") + " Orthorectification of the Math layers")
        input_folder_for_ortho = Fld_math_layers
        output_folder_ortho = Fld_math_layers_ortho

        ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel,
               ortho_bounds_option, AOI_file, AOI_file_segment_number, ortho_resolution_X,
               ortho_resolution_Y, generate_overviews, TSA_math_xtra_channels)
      
        # -------------------------------------------------------------------------------------------------------
        # Math layers splitting
        output_folder = Fld_math_layers_ortho_split
        suffix = " intensity"
        math_layers_split (Fld_math_layers_ortho, output_folder, prefix, TSA_math_xtra_channels, 
                           TSA_math_xtra_labels, suffix)

        proc_stop_time = time.time()
        folder = Fld_math_layers_ortho_split
        out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
        string_1 = ("Math layer(s) orthorectification:" + out_folder_time_size) 
        time_log.write("%s\n" % string_1)

# --------------------------------------------------------------------------------------------------------------
# E) Preparation for the Time Series analysis
# --------------------------------------------------------------------------------------------------------------
print("\t")
print(" ------------------------------------------------------------------------------------------------------")
print("                 Time Series Analysis - Files preparation  (Output lists creation)                     ")
print(" ------------------------------------------------------------------------------------------------------")
print("\t")
proc_start_time = time.time()

if produce_coherence_layers is True:
    # Creating output lists of files for the stacking.
    print(time.strftime("%H:%M:%S") + " Coherence - lists preparation for scene stacking")

    search_folder = Fld_coherence_ortho
    suffix = "_coh"
    TSA_channels = TSA_coherence_layers
    create_list (prefix, suffix, search_folder, Fld_Output_stack_lists, TSA_channels, TSA_channel_labels)


# Creating output lists of files for the stacking.
if produce_intensity_layers is True:
    print(time.strftime("%H:%M:%S") + " Intensity - lists preparation for scene stacking")

    search_folder = Fld_intensity_ortho
    suffix = "_int"
    TSA_channels = TSA_intensity_layers
    create_list (prefix, suffix, search_folder, Fld_Output_stack_lists, TSA_channels, TSA_channel_labels)

proc_stop_time = time.time()
folder = Fld_Output_stack_lists
out_folder_time_size = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
string_1 = ("Time series list preparation:" + out_folder_time_size) 
time_log.write("%s\n" % string_1)

'''
# ------------------------------------------------------
# J) Deleting intermediary files if requested
# ---------------------------------------------------------

if delete_intermediary_files in yes_validation_list:
    print("\t")
    print("--------------------------------------------------------------")
    print("\t")
    print((time.strftime("%H:%M:%S")) + "  Deleting intermediary files")

    del_folders_list = []
    # List of all (possible) folders to delete
    del_folders_list.append(os.path.join(output_folder, "1_INSINFO"))
    del_folders_list.append(os.path.join(output_folder, "2_SARINGESTAOI"))
    del_folders_list.append(os.path.join(output_folder, "4_RAW"))
    del_folders_list.append(os.path.join(output_folder, "5_DEFO"))

    for delete in del_folders_list:
        if os.path.isdir(delete):
            print("deleting " + delete)
            shutil.rmtree(delete)

'''
print("\t")
print(" ------------------------------------------------------------------------------------------------------")
print(" ------------------------------------------------------------------------------------------------------")
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

