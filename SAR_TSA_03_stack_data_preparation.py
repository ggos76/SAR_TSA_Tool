#!/usr/bin/env python
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. May 2024                                      -
 * ----------------------------------------------------------------------
'''

# ----------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# ----------------------------------------------------------------------------------------------
# A) Input/Output
Coregistered_Pairs_Report = r"E:\SAR_TSA_tests\RB2010_FQ22\2_Coregistered_Scenes\RB2010_04_Coregistered_Pairs_Report.txt"
output_folder = r"E:\SAR_TSA_tests\RB2010_FQ22"
prefix = "RB2010_"

# B) Elevation source
DEM_file = r"D:\HBL_S1A_mapping\sHBL_GLO30DEM_LongLat_D000_with_0_elev.pix"
DEM_elevation_channel = 1

# C)  The channels to process for the time series (TSA) generation, can be a subset of the coregistered files channels
# C.1) These parameters must be common to all scenes to be processed. The TSA_channel_mapping list must include all
#       channels of the input file and not only the channels to be processed.
TSA_channel_mapping = [1,2,3,4]
TSA_channel_labels = ["HH", "HV", "VH", "VV"]

# C.1) Coherence  options
produce_coherence_layers = "yes"    # Only available for SLC data
TSA_coherence_layers = [1, 4]      # Must be all subset of TSA_channel_mapping
apply_insraw_filter = "yes"
filter_size = 9

# C.2) Intensity options
produce_intensity_layers = "yes"
TSA_intensity_layers = [1, 2, 4]    # Must be all subset of TSA_channel_mapping

# C.3) Other intensity layers
produce_incidence_angle_layer = "yes"

produce_math_layers = "no"
TSA_math_models_definition = r""
TSA_math_labels = ["HH/VV", "HH/VV"]

# D) Orthorectification options
# Ortho bounds options: 1 (from an AOI file)  or 2 (from the input file)
ortho_bounds_option = 1
AOI_vector_file = r"D:\Wapusk_2010\Roberge_lake\AOI_Roberge_Lake_UTM15_D000.pix"
AOI_segment_number = 2

ortho_resolution_X = "7"
ortho_resolution_Y = "7"

# E) General options
# Generate overviews - either yes or no,
generate_overviews = "yes"
# keep or delete intermediate files - either yes or no. No is recommended.
delete_intermediary_files = "no"

# -----------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
1) There is no verification if the specified DEM covers completely the spatial
   extents of the AOI. If not the script will run to completion but all
   interferograms will be blank after INSRAW.

'''
# ---------------------------------------------------------------------------------------------
#  Main program
# ----------------------------------------------------------------------------------------------
import sys
import os
import fnmatch
import time
import locale

import pci
from pci.insraw import insraw
from pci.pyramid import pyramid
from pci.exceptions import PCIException
from pci.api import datasource as ds

from TSA_utilities.SAR_TSA_utilities_definitions import ortho_run
from TSA_utilities.SAR_TSA_utilities_definitions import nan_replace
from TSA_utilities.SAR_TSA_utilities_definitions import psiqinterp_run
from TSA_utilities.SAR_TSA_other_layers import inc_angle_layer
from TSA_utilities.SAR_TSA_utilities_definitions import create_list

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -------------------------------------------------------------------------
#  PART 2: Validation
# -------------------------------------------------------------------------

# Hardcoded parameters
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list + no_validation_list
GB = 1073741824

script2_procTime = os.path.join(output_folder, prefix + "script3_INSCOREG_ProcessingTime.txt")
time_log = open(script2_procTime, "w")
Fld_Coregistration = os.path.join(output_folder, "2_Coregistered_Scenes")
AOI_file = AOI_vector_file
AOI_file_segment_number = AOI_segment_number

#  Version control - do nothing for now.
print("\t")
print(pci.version)

print("Installed python version: " + sys.version)
py1 = str(sys.version_info[0])
py2 = str(sys.version_info[1])
py3 = (py1 + "." + py2)
python_version = float(py3)
if python_version < 3.6:
    print("You are using Python v" + str(python_version))
    print("You need to update to Python 3.6 or newer versions")
    sys.exit()
print("\t")

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

# B) Elevation source
if not os.path.exists(DEM_file):
    print ("Error - The DEM_file does not exist or the path/filename is wrong.")
    sys.exit()

#-------------------------------------------------------------------------------------------
# C.0) Channel mapping.  We check for the file input conformity

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
    produce_intensity_layers = False

produce_math_layers = produce_math_layers.lower()
if produce_math_layers in yes_validation_list:
    produce_math_layers = True
else:
    produce_math_layers = False

# Check if a least one layer type is selected:
if (produce_coherence_layers is False) and (produce_intensity_layers is False) and (produce_math_layers is False):
    print ("Error - At least one type of output layer must be selected")
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
    print ("Checking input scene conformity--->" + ii)
    with ds.open_dataset(ii, ds.eAM_READ) as ds2:
        aux = ds2.aux_data
        dataType = aux.get_file_metadata_value("dataType")
        SensorModelName = aux.get_file_metadata_value("SensorModelName")
        num_channels = ds2.chan_count

        # check size here!
        # zzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzzz
        # check size here!
        if len(TSA_channel_mapping) != num_channels:
            print ("Error 11 - Mismatch between the length of TSA_channel_mapping  and the number of input channels")
            sys.exit()

        if produce_coherence_layers is True:
            if dataType != "Complex":
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

# Validation for the Raw interferograms filter options
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

if produce_incidence_angle_layer is True:
    incidence_layer_scene = check_scenes_list[0]


# -------------------------------------------------------------------------------------------
# C.5 Math Layers



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
if generate_overviews not in yes_no_validation_list:
    print("\t")
    print("Error - Selected overviews option is not valid.")
    print('Accepted values are: "yes" or "no"')
    sys.exit()

delete_intermediary_files = delete_intermediary_files.lower()
if delete_intermediary_files not in yes_no_validation_list:
    print("\t")
    print("Error - Delete intermediary file option is not valid.")
    print('Accepted values are: "yes" or "no"')
    sys.exit()


Fld_Output_stack_lists = os.path.join(output_folder, "6_Output_stack_lists")
if not os.path.exists(Fld_Output_stack_lists):
    os.makedirs(Fld_Output_stack_lists)

start = time.time()
# ---------------------------------------------------------------------------------------------------------------
#                                                  Main program
# ---------------------------------------------------------------------------------------------------------------
if produce_coherence_layers is True:

    print("\t")
    print("-------------------------------------------------------------------------------------------------------")
    print("                               RAW INTERFEROGRAMS GENERATION                                           ")
    print("-------------------------------------------------------------------------------------------------------")
    print("\t")

    insraw_start = time.time()
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

                    if generate_overviews in yes_validation_list:
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

    insraw_stop = time.time()
    string1 = "3-INSRAW total proc time;"
    ellapse_time = str(round(insraw_stop - insraw_start, 2))
    string2 = string1 + ellapse_time
    time_log.write("%s\n" % string2)

    # Replacings NANs (if any) by NoDATA in the output RAW interferograms
    step_nans = "INSRAW_"
    input_folder_nans = Fld_Raw_Interferograms
    nan_replace(step_nans, input_folder_nans)

    # ----------------------------------------------------------------------------------------------------------------
    # Generating the coherence layers
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

    psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_layers,
                    TSA_labels, ps_output_folder, unique_files, prefix)

    # ----------------------------------------------------------------------------------------------------------------
    # Orthorectification
    print ("\t")
    print(time.strftime("%H:%M:%S") + " Orthorectification of the coherence layers ")

    input_folder_for_ortho = Fld_coherence
    #output_folder_ortho_raw = Fld_coherence_ortho
    output_folder_ortho = Fld_coherence_ortho

    ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel,
               ortho_bounds_option, AOI_file, AOI_file_segment_number, ortho_resolution_X,
               ortho_resolution_Y, generate_overviews)


if produce_intensity_layers is True:

    print("\t")
    print("-------------------------------------------------------------------------------------------------------")
    print("                              Generating the intensity layers                                          ")
    print("-------------------------------------------------------------------------------------------------------")
    print("\t")

    print(time.strftime("%H:%M:%S") + " Creating the intensity layers")
    search_folder = Fld_Coregistration
    TSA_layers = TSA_intensity_layers
    TSA_labels = TSA_channel_labels
    unique_files = "yes"
    keyword = "*.pix"
    interp_type = "int"
    suffix = "int"
    ps_output_folder = Fld_intensity

    psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_layers,
                    TSA_labels, ps_output_folder, unique_files, prefix)


if produce_incidence_angle_layer is True:

    input_file = check_scenes_list[0]
    output_folder_inc = Fld_intensity
    inc_angle_layer (input_file, output_folder_inc, prefix)


# ----------------------------------------------------------------------------------------------------------------
# Orthorectification
print(time.strftime("%H:%M:%S") + " Orthorectification of the intensity layers")
input_folder_for_ortho = Fld_intensity
#output_folder_ortho_coreg = Fld_intensity_ortho
output_folder_ortho = Fld_intensity_ortho

ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel,
           ortho_bounds_option, AOI_file, AOI_file_segment_number, ortho_resolution_X,
           ortho_resolution_Y, generate_overviews)





# -------------------------------------------------------------------------------------------------
# E) Preparation for the Time Series analysis
# -------------------------------------------------------------------------------------------------

'''
print("\t")
print(" ------------------------------------------------------------------------------------------")
print("                   Time Series Analysis - Files preparation                                ")
print(" ------------------------------------------------------------------------------------------")
print("\t")

# Orthorectified coregistered files, convert the complex channel(s) to intensity
print(time.strftime("%H:%M:%S") + " Converting the orthorectified scenes to intensity")
search_folder = output_folder_ortho_coreg
keyword = "o*.pix"
interp_type = "int"
suffix = "_int"

psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_channels, TSA_channels_label )


# Orthorectified Raw interferograms, convert the complex channel(s) to coherence
print ("\t")
print(time.strftime("%H:%M:%S") + " Converting the orthorectified RAW Inteferograms to coherence")
search_folder = output_folder_ortho_raw
keyword = "oraw*.pix"
interp_type = "amp"
suffix = "_coh"
psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_channels, TSA_channels_label)

'''

print("\t")
print(" ------------------------------------------------------------------------------------------")
print("                 Time Series Analysis - Files preparation  (Output lists creation)         ")
print(" ------------------------------------------------------------------------------------------")
print("\t")

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
print("--------------------------------------------------------------")
print((time.strftime("%H:%M:%S")))
print("All processing completed")
print("\t")
end = time.time()

ellapse_time_seconds = round((end - start), 2)
ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

print("Processing time (seconds): " + str(ellapse_time_seconds))
print("Processing time (minutes): " + str(ellapse_time_minutes))
print("Processing time (hours): " + str(ellapse_time_hours))
string1 = "10-Total proc time;" + str(ellapse_time_seconds)
time_log.write("%s\n" % string1)

time_log.close()

