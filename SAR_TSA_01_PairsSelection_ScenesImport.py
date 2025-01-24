#!/usr/bin/env python
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. May 2024                                      -
 * ----------------------------------------------------------------------
'''

# -----------------------------------------------------------------------------------------------------
#  User defined variables
# ----------------------------------------------------------------------------------------------
# The Input scenes can be retrieved from an input folder or an mfile (*.txt)
parent_folder_search = r"D:\Wapusk_2010\Roberge_lake\FQ22\unzip"
keyword = "product.xml"

# Applies only to Sentinel-1 data. Ingest a single swath or all swaths.
# Options are 1, 2, 3, 4 (all swaths)
sentinel_swath = 4

# Outputs
# Specify a prefix for the outputs - suggest ending it with '_' (optional)
prefix = "RB2010_"
output_folder = r"E:\SAR_TSA_tests\RB2010_FQ22"

# Elevation source
DEM_file = r"D:\HBL_S1A_mapping\sHBL_GLO30DEM_LongLat_D000_with_0_elev.pix"
DEM_elevation_channel = 1

# Start and Stop date        # Valid options are yes and no
# Must be in the following format YYYYMMDD
use_start_stop_date =  "no"
start_date_YYYYMMDD = 20180612
stop_date_YYYYMMDD = 20180925

# InSAR pairs selection mode
# 1: All pairs temporal             2: Single reference file
# 3: With baseline filters (SBAS)   4: Subsequent pairs
pairs_selection_mode = 4

# Only used when pairs_selection_mode=2
# Must be in the following format YYYYMMDD, for example   20190519
pairs_selection_mode_2_referencedate = 20140817

# Use only  with pairs_selection_mode=3; distances are in meters
absolute_perp_baseline_min = 0
absolute_perp_baseline_max = 250
# Temporal baseline in days,  min=0 ,  max=365
max_temporal_baseline = 70

# Subsequent pairs options.
# Number of pairs to form for every input scene. The files at the beginning and the end of the time series
# might have fewer pairs especially if number_subsequent_pairs is set with a large number.
# The candidate pair must have been acquired within the specified max_delta_days.
subsequent_pairs_no = 1     # Must be an integer
max_delta_days = 50

# Subset options # Accepted values are yes or no
subset_input_data = "no"

# Area of Interest settings
AOI_vector_file = r"D:\HBL_S1A_mapping\F186_clip_area_UTM_small.pix"
AOI_segment_number = 2

# How to fill areas beyond the AOI.
# Can be box or NoData, box is recommended for optimal pairs coregistration.
subset_filling_option = "box"

# Generate overviews - either yes or no,
generate_overviews = "yes"
# keep or delete intermediate files - either yes or no. No is recommended.
delete_intermediary_files = "no"

# ----------------------------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# ----------------------------------------------------------------------------------------------------------------------
'''
1) There is no verification if the specified DEM covers completely the spatial
   extents of the AOI. If not the script will run to completion but all
   interferograms will be blank after INSRAW.

'''
# ---------------------------------------------------------------------------------------------------------------------
#  Main program
# ---------------------------------------------------------------------------------------------------------------------
import sys
import os
import time
import locale
import fnmatch
from datetime import datetime

import pci
from pci.saringestaoi import saringestaoi
from pci.psinang import psinang
from pci.pcimod import pcimod
from pci.iii import iii
from pci.exceptions import PCIException
from pci.api import datasource as ds
from pci.api.inputsource import InputSourceFactoryBuilder

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")


# ---------------------------------------------------------------------------------------------------------------------
#  VALIDATION
# ---------------------------------------------------------------------------------------------------------------------
# A.1) Hardcoded parameters
Calibration_type = "sigma"
Output_filename_format = 2
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list+no_validation_list
GB = 1073741824

# PCI Version control - do nothing for now.
print("\t")
print(pci.version)

print("Installed python version: " + sys.version)
py1 = str(sys.version_info[0]); py2 = str(sys.version_info[1]); py3 = (py1 + "." + py2); python_version = float(py3)
if python_version < 3.6:
    print("You are using Python v" + str(python_version))
    print("You need to update to python 3.6 or newer versions")
    sys.exit()
print("\t")

if not os.path.exists(parent_folder_search):
    print("The specified parent_folder_search does not exists")
    print("Verify please the input scenes folder path: ")
    print(parent_folder_search)
    sys.exit()

# A.2) Creation of the INSINFO output folder and baseline report
fld_scenes_import = os.path.join(output_folder, "1_Scenes_Import")
if not os.path.exists(fld_scenes_import):
    os.makedirs(fld_scenes_import)
Output_InSARpairs_report_name = prefix + "02_Pair_selection_and_baselines_list.txt"

script1_procTime = os.path.join(output_folder, prefix + "script1_INSINFO_SARINGESTAOI_ProcessingTime.txt")
time_log = open(script1_procTime, "w")

# A.3) INSINFO EARLY VALIDATION
if absolute_perp_baseline_min < 0:
    print("Error - The perpendicular baseline minimum must be equal or superior to 0 (absolute value)")
    print("Program stopped")
    print("\t")
    sys.exit()

if absolute_perp_baseline_max <= absolute_perp_baseline_min:
    print("Error - The perpendicular baseline maximum value must be superior to the minimum value")
    print("Program stopped")
    print("\t")
    sys.exit()

if max_temporal_baseline <= 0 or max_temporal_baseline > 365:
    print("Error - The minimal temporal baseline is 0 days and the maximum is 365")
    print("Program stopped")
    print("\t")
    sys.exit()

if use_start_stop_date.lower() not in yes_no_validation_list:
    print("Error - Use Start and stop date option must be set with ""yes"" or ""no""")
    print("Program stopped")
    sys.exit()

if start_date_YYYYMMDD > stop_date_YYYYMMDD:
    print("Error - Start date is superior to Stop date")
    print("Program stopped")
    sys.exit()

# A.4)Sentinel-1 subswath selector
if sentinel_swath not in [1, 2, 3, 4]:
    print("\t")
    print("Error - sentinel_swath  option is not valid")
    print("Accepted values are: 1, 2, 3, 4")
    sys.exit()

# A.5) SARINGESTAOI EARLY VALIDATION
if Output_filename_format not in [1, 2, 3, 4]:
    print("\t")
    print("Error - Output file name format is not valid")
    print("Specified value is " + str(Output_filename_format))
    print("Accepted values are: 1, 2, 3 or 4")
    sys.exit()

if Calibration_type.lower() not in ["sigma", "gamma", "beta", "none"]:
    print("Error- Calibration type is invalid")
    print("Specified value is " + Calibration_type)
    print('"Accepted values are: "sigma", "beta", "gamma" or "none"')
    sys.exit()

if generate_overviews.lower() not in yes_no_validation_list:
    print("\t")
    print("Error - Selected overviews option is not valid.")
    print('Accepted values are: "yes" or "no"')
    sys.exit()

# A.6) Subset options validation
if subset_input_data.lower() not in yes_no_validation_list:
    print("Error- Specified subset option is invalid")
    print("Specified value is " + str(subset_input_data))
    print('Accepted values are: "yes" or "no"')
    sys.exit()

if subset_input_data.lower() in yes_validation_list:

    if not os.path.exists(AOI_vector_file):
        print("AOI vector file not found")
        print("Provide a valid path and filename")
        sys.exit()

    if not os.path.exists(DEM_file):
        print("DEM  file not found")
        print("Provide a valid path and filename")
        sys.exit()

generate_incidence_angle_layer = generate_incidence_angle_layer.lower()
if generate_incidence_angle_layer in yes_validation_list:
    generate_incidence_angle_layer = True
else:
    generate_incidence_angle_layer = False

#----------------------------------------------------------------------------------------------
#  B) Find input scenes from an MFILE or from an input folder with a specified keyword
#----------------------------------------------------------------------------------------------
start = time.time()

# Determine the input raw scenes - either from a folder or mfile textfile.
# in both cases a keyword can be used to restrict the input files

print ("\t")
print ("-------------------------------------------------------------------------------------------------------------")
print ("                      Search for input scenes and Metadata retrieval                                         ")
print ("-------------------------------------------------------------------------------------------------------------")

print((time.strftime("%H:%M:%S")) + " Search for input scenes")
print("   search folder: " + parent_folder_search)
print("   keyword: " + keyword)
factory = InputSourceFactoryBuilder().build_factory()
inputSource = factory.create_input_source(parent_folder_search, recursive=True,
                                          mask=keyword)
VendorInput = list(inputSource)
print ("   Number of scenes found: " + str (len(VendorInput)))
if len(VendorInput) < 2:
    print("A minimum of two images must be provided as input.")
    print("Program stopped")
    print("\t")
    sys.exit()


# write out an mfile of the discovered files for user info
file = os.path.join(fld_scenes_import, prefix + "01_Discovered_Scenes_list.txt")
with open(file, "w") as f:
    f.write("\n".join(VendorInput))

print ("\t")
print((time.strftime("%H:%M:%S")) + " Metadata Retrieval")
Acquisition_DateTime2 = []
Acquisition_DateTime3 = []
SensorModelName_list = []
count = 1

for ii in VendorInput:
    print("   " + ((time.strftime("%H:%M:%S")) + " Retrieving metadata, file " +
          str(count) + " of " + str(len(VendorInput)) + "..."))

    ii = ii.strip()
    with ds.open_dataset(ii, ds.eAM_READ) as ds2:
        aux = ds2.aux_data
        Acquisition_DateTime = aux.get_file_metadata_value(
            "Acquisition_DateTime")
        SensorModelName = aux.get_file_metadata_value("SensorModelName")
        SensorModelName_list.append(SensorModelName)
        # keep only the first 10 characters of the acquisition date
        Acquisition_DateTime2.append(Acquisition_DateTime[:10])
        count = count + 1

for ii in Acquisition_DateTime2:
    # Remove the hyphen for compact name.
    Acquisition_DateTime3.append(ii.replace("-", ""))

# A.4.2) Remove from the two lists all acquisitions that are not between the start and stop date.
if use_start_stop_date in yes_validation_list:
    Acquisition_DateTime3_b = []
    VendorInput_b = []
    for date, input_file in zip(Acquisition_DateTime3, VendorInput):
        date = int(date)

        if date < start_date_YYYYMMDD or date > stop_date_YYYYMMDD:
            print (str(date) + "--->" + input_file)
            print("   Acquisition date outside the range of start and stop date")
        else:
            Acquisition_DateTime3_b.append(date)
            VendorInput_b.append(input_file)

    Acquisition_DateTime3 = Acquisition_DateTime3_b
    VendorInput = VendorInput_b

# B.4) Sort the VendorInput list by acquisition date, important for the coregistration step. )
print ("\t")
print((time.strftime("%H:%M:%S")) + " Sorting the scenes by acquisition date")

pairs = list(zip(Acquisition_DateTime3, VendorInput))
pairs.sort()
# We split the two sorted list and crete two new sorted list.
Acquisition_DateTime3_sorted = []
VendorInputs_sorted = []

for ii in pairs:
    VendorInputs_sorted.append(ii[1].strip())
    Acquisition_DateTime3_sorted.append(ii[0])

print("   Number of sorted files: " + str(len(Acquisition_DateTime3_sorted)))
for ii, jj in zip(Acquisition_DateTime3_sorted, VendorInputs_sorted):
    print("   " + str(ii) + "-->" + jj)


print("\t")
print("--------------------------------------------------------------------------------------------------------------")
print("                                           Pairs selection                                                    ")
print("--------------------------------------------------------------------------------------------------------------")
print("\t")

insinfo_start = time.time()

maximg = len(VendorInputs_sorted)
insinfo_report_list = []
reference_file_list = []
dependent_file_list = []
date_YYYYMMDD_ref = []
date_YYYYMMDD_dep = []

#--------------------------------------------------------------------------------------------------
# 1- All pairs temporal
pair = 0
if pairs_selection_mode == 1:

    print("Selected coregistration option:1 - All pairs temporal")
    print("\t")

    for ii in range(0, maximg):
        for jj in range(0, maximg):
            if ii < jj:

                filref = VendorInputs_sorted[ii]
                mfile = VendorInputs_sorted[jj]

                ref = Acquisition_DateTime3_sorted[ii]
                dep = Acquisition_DateTime3_sorted[jj]
                tfile = os.path.join(fld_insinfo, prefix + "_" + str(pair)
                                     + "_INSINFO" + "_ref_" + ref
                                     + "_dep_" + dep + ".txt")

                print("\t")
                print(((time.strftime("%H:%M:%S")) + " Generating the INSINFO report"))
                print("Reference file: " + filref)
                print("Dependent file: " + mfile)

                try:
                    insinfo(filref, mfile, tfile)
                    insinfo_report_list.append(tfile)
                    reference_file_list.append(filref)
                    dependent_file_list.append(mfile)
                    date_YYYYMMDD_ref.append(Acquisition_DateTime3_sorted[ii])
                    date_YYYYMMDD_dep.append(Acquisition_DateTime3_sorted[jj])
                    pair = pair + 1

                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)
#--------------------------------------------------------------------------------------------------
# 2- Single reference file
pair = 0
if pairs_selection_mode == 2:

    print("Selected coregistration option: 2 - Single reference file")
    print("\t")

    ii_ref = -1000  # impossible value
    for c, value in enumerate(Acquisition_DateTime3_sorted):
        if value == pairs_selection_mode_2_referencedate:
            print(c, value)
            ii_ref=c
    # Validation in case where the specified date is not found
    if ii_ref == -1000:
        print("Error - The specified date is not valid ( " + pairs_selection_mode_2_referencedate + " )")
        print("Valid date are")
        for dates in Acquisition_DateTime3_sorted:
            print(dates)
        sys.exit()

    for jj in range(maximg):
        if ii_ref != jj:
            filref = VendorInputs_sorted[ii_ref]
            mfile = VendorInputs_sorted[jj]

            ref = Acquisition_DateTime3_sorted[ii_ref]
            dep = Acquisition_DateTime3_sorted[jj]
            tfile = os.path.join(fld_insinfo, prefix + str(pair)
                                 + "_INSINFO" + "_ref_" + ref + "_dep_"
                                 + dep + ".txt")

            print("\t")
            print(((time.strftime("%H:%M:%S")) +
                  " Generating the INSINFO report"))
            print("Reference file: " + filref)
            print("Dependent file: " + mfile)

            try:
                insinfo(filref, mfile, tfile)
                insinfo_report_list.append(tfile)
                reference_file_list.append(filref)
                dependent_file_list.append(mfile)
                date_YYYYMMDD_ref.append(
                    Acquisition_DateTime3_sorted[ii_ref])
                date_YYYYMMDD_dep.append(Acquisition_DateTime3_sorted[jj])
                pair = pair + 1

            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
#--------------------------------------------------------------------------------------------------
# 3- Baseline filters

pair = 0
if pairs_selection_mode == 3:
    from pci.insinfo import insinfo

    print("Selected coregistration option:3 - Baseline filters")
    print("\t")

    for ii in range(0, maximg):
        for jj in range(0, maximg):
            if ii < jj:

                filref = VendorInputs_sorted[ii]
                mfile = VendorInputs_sorted[jj]

                ref = Acquisition_DateTime3_sorted[ii]
                dep = Acquisition_DateTime3_sorted[jj]
                tfile = os.path.join(fld_insinfo, prefix + str(pair)
                                     + "_INSINFO" + "_ref_" + ref
                                     + "_dep_" + dep + ".txt")

                print("\t")
                print(((time.strftime("%H:%M:%S")) +
                        " Generating the INSINFO report"))
                print("Reference file: " + filref)
                print("Dependent file " + mfile)

                try:
                    insinfo(filref, mfile, tfile)
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)

                # Find the temporal baseline
                keyword = "Time Difference"
                with open(tfile) as myFile:
                    for num, line in enumerate(myFile, 1):
                        if keyword in line:
                            tbaseline = line
                            print(line)

                tbaseline2 = tbaseline.split(":", 4)
                tbaseline3 = str(tbaseline2[1]) + ":" + str(tbaseline2[2])
                tbaseline4 = tbaseline3.replace(":", ".")
                temp_base = float(tbaseline4)

                #Find the perpendicular baseline
                keyword = "Middle Line Perpendicular Baseline Distance"
                with open(tfile) as myFile:
                    for num, line in enumerate(myFile, 1):
                        if keyword in line:
                            pbaseline = line

                pbaseline2 = pbaseline.split(":")
                Pbase = abs(float(pbaseline2[1]))

                if (Pbase > absolute_perp_baseline_min and Pbase < absolute_perp_baseline_max and temp_base < max_temporal_baseline):

                    print("Baseline conditions satisfied")
                    insinfo_report_list.append(tfile)
                    reference_file_list.append(filref)
                    dependent_file_list.append(mfile)
                    date_YYYYMMDD_ref.append(
                        Acquisition_DateTime3_sorted[ii])
                    date_YYYYMMDD_dep.append(
                        Acquisition_DateTime3_sorted[jj])
                    pair = pair + 1
                else:
                    print("Baselines conditions not satisfied")
                    os.remove(tfile)
                # Break statement, since the scenes are time ordered, if the temporal baseline
                # is already exceeded, there is no need to continue to loop through the rest
                # of the dependent files list.
                if temp_base > max_temporal_baseline:
                    break
#--------------------------------------------------------------------------------------------------
# 4- Subsequent pairs
#    Modified version so we don't rely on INSINFO

if pairs_selection_mode == 4:

    date_format = "%Y%m%d"
    pair_selection_4_list = []

    # Hardcoded values since we don't have access to INSINFO
    Perpendicular_baseline = "1000"

    print("Selected coregistration option: 4- Subsequent pairs")
    print("\t")

    pair = 0
    attempt = 1
    for ii in range(0, maximg - 1):
        for kk in range (1, subsequent_pairs_no + 1):
            jj = ii + kk
            print("\t")
            print(((time.strftime("%H:%M:%S")) + " Evaluating the candidate pair"))
            print ("   subsequent pair number:" + str(subsequent_pairs_no) + " // Max delta days: "
                   + str (max_delta_days))

            if jj <= maximg - 1:
                filref = VendorInputs_sorted[ii]
                mfile = VendorInputs_sorted[jj]
                ref = Acquisition_DateTime3_sorted[ii]
                dep = Acquisition_DateTime3_sorted[jj]
                print("   Reference file: " + ref)
                print("   Dependent file: " + dep)

                a = datetime.strptime(str(ref), date_format)
                b = datetime.strptime(str(dep), date_format)
                delta = b - a
                if delta.days <= max_delta_days:
                    print ("   Time difference (days): " + str(delta.days) + "  ---> Valid pair")
                    print("   pair " + str(pair))

                    Temporal_baseline = str(delta.days)

                    string1 = (str(pair) + ";" + filref +";" + str(ref) + ";" + mfile + ";" + str(dep) + ";" +
                               Temporal_baseline+ ";"+ Perpendicular_baseline)
                    pair_selection_4_list.append(string1)

                    pair = pair + 1
                else:
                    print("   Pair rejected - Time difference exceeding the max_delta_days ")
                    print("   Temporal baseline: "+ str(delta.days) + " (max_delta_days: " +
                          str(max_delta_days) + ")")
            else:
                print("   Pair rejected - Out of list index")

            attempt = attempt + 1

    file=open(os.path.join(fld_scenes_import,Output_InSARpairs_report_name), "w")
    file.write('\n'.join(pair_selection_4_list))
    file.close()


#-----------------------------------------------------------------------------
# C) Data retrieval from the generated INSINFO report(s)
#-----------------------------------------------------------------------------
if pairs_selection_mode == 3:

    verbose_log = os.path.join(fld_scenes_import, Output_InSARpairs_report_name)
    with open(verbose_log, "w") as f:

        pair = 0
        for ii in insinfo_report_list:

            keyword = "Time Difference"
            with open(ii) as myFile:
                for num, line in enumerate(myFile, 1):
                    if keyword in line:
                        tbaseline = line

            tbaseline2 = tbaseline.split(":", 4)
            tbaseline3 = str(tbaseline2[1]) + ":" + str(tbaseline2[2])
            tbaseline4 = tbaseline3.replace(":", ".")
            Temporal_baseline = str(tbaseline4)

            keyword = "Middle Line Perpendicular Baseline Distance"
            with open(ii) as myFile:
                for num, line in enumerate(myFile, 1):
                    if keyword in line:
                        pbaseline = line

            pbaseline2 = pbaseline.split(":")
            print(pbaseline2[1])
            Perpendicular_baseline = str(pbaseline2[1])

            string1 = (str(pair) + ";" + reference_file_list[pair] + ";"
                       + date_YYYYMMDD_ref[pair] + ";" + dependent_file_list[pair]
                       + ";" + date_YYYYMMDD_dep[pair] + ";" + Temporal_baseline
                       + ";" + Perpendicular_baseline)

            f.write("%s\n" % string1)
            pair = pair + 1

    insinfo_stop = time.time()
    string1 = "1-INSTOPO proc time;"
    ellapse_time = str(round(insinfo_stop - insinfo_start, 2))
    string2 = string1 + ellapse_time
    time_log.write("%s\n" % string2)


print("\t")
print("--------------------------------------------------------------------------------------------------------------")
print("                                   Selected scenes ingestion                                                  ")
print("--------------------------------------------------------------------------------------------------------------")
print("\t")

metadata_start = time.time()
# B) Read and parse the baseline report
Input_file_lines = []
pair_number = []
Reference_file = []
Reference_acquisition_date = []
Dependent_file = []
Dependent_acquisition_date = []

print(time.strftime("%H:%M:%S") + " Reading and parsing the baseline report")
print("\t")

with open(os.path.join(fld_scenes_import, Output_InSARpairs_report_name), "r") as input_file:
    for line in input_file:
        if len(line.strip()) > 0:
            Input_file_lines.append(line)

for ii in Input_file_lines:
    temp = ii.split(";")
    pair_number.append(temp[0])
    Reference_file.append(temp[1])
    Reference_acquisition_date.append(temp[2])
    Dependent_file.append(temp[3])
    Dependent_acquisition_date.append(temp[4])

pairs = len(pair_number) + 1

# Merge temporarily the two lists it will contains some redondance as a
# file can be represented multiple times.
Reference_Dependent_file_list = Reference_file + Dependent_file

#-----------------------------------------------------------------------------------------
#C) Metadata acquisition

print(time.strftime("%H:%M:%S") + " Metadata acquistion")
print("\t")

Acquisition_DateTime2 = []
Acquisition_DateTime3 = []
SensorModelName2 = []
Acquisition_Type2 = []
SourceID2 = []
PlatformName2 = []
Orbit_Direction2 = []
BeamMode2 = []
ProductType2 = []

count = 1
for ii in Reference_Dependent_file_list:

    print(((time.strftime("%H:%M:%S")) + " Retrieving file metadata, file "
          + str(count) + " of " + str(len(Reference_Dependent_file_list))))
    print(ii)
    print("\t")

    with ds.open_dataset(ii, ds.eAM_READ) as ds2:
        aux = ds2.aux_data

        Acquisition_DateTime = aux.get_file_metadata_value(
            "Acquisition_DateTime")
        Acquisition_DateTime2.append(Acquisition_DateTime)

        Acquisition_Type = aux.get_file_metadata_value("Acquisition_Type")
        Acquisition_Type2.append(Acquisition_Type)

        SensorModelName = aux.get_file_metadata_value("SensorModelName")
        SensorModelName2.append(SensorModelName)

        SourceID = aux.get_file_metadata_value("SourceID")
        SourceID2.append(SourceID)

        PlatformName = aux.get_file_metadata_value("PlatformName")
        PlatformName2.append(PlatformName)

        Orbit_Direction = aux.get_file_metadata_value("OrbitDirection")
        Orbit_Direction2.append(Orbit_Direction)

        BeamMode = aux.get_file_metadata_value("BeamMode")
        BeamMode2.append(BeamMode)

        ProductType = aux.get_file_metadata_value("ProductType")
        ProductType2.append(ProductType)

    count = count + 1

for ii in Acquisition_DateTime2:
    a = ii.replace("-", "")
    b = a.replace(":", "")
    Acquisition_DateTime3.append(b[:8])

metadata_stop = time.time()
string1 = "1-Metadata retrieval proc time;"
ellapse_time = str(round(metadata_stop - metadata_start, 2))
string2 = string1 + ellapse_time
time_log.write("%s\n" % string2)

#-----------------------------------------------------------------------------------------
#D) Build output filename based on user choice
#-----------------------------------------------------------------------------------------

output_file_name = []

if Output_filename_format == 1:    # Source ID
    for ii in SourceID2:
        temp = prefix + ii + ".pix"
        output_file_name.append(temp)

elif Output_filename_format == 2:   # Compact date
    for date in Acquisition_DateTime3:
        temp = prefix + date + ".pix"
        output_file_name.append(temp)

elif Output_filename_format == 3:
    for type, date in zip(Acquisition_Type2, Acquisition_DateTime3):
        temp = prefix + type + "_" + date + ".pix"
        output_file_name.append(temp)

elif Output_filename_format == 4:
    for type, date, orbit in \
            zip(SensorModelName2, Acquisition_DateTime3, Orbit_Direction2):
        temp = prefix + type + "_" + date + "_" + orbit + ".pix"
        output_file_name.append(temp)

#-----------------------------------------------------------------------------------------
#E) SARINGESTAOI
#   Ingest the files stored in the MFILE , use the acquisition date to name the files.
#-----------------------------------------------------------------------------------------
SaringestAOI_start = time.time()
# Split the list here
#A = [1,2,3,4,5,6]
#B = A[:len(A)//2]
#C = A[len(A)//2:]

ReferenceFiles_output_FilenamesList=output_file_name[:len(
    output_file_name) // 2]   # First half of the list
DependentFiles_output_FilenamesList=output_file_name[len(
    output_file_name) // 2:]   # Second half of the list

print("\t")
print("---------------------------------------------------------------------------------------------------------------")
print("                                         SAR Scenes Ingestion                                                  ")
print("---------------------------------------------------------------------------------------------------------------")
print("\t")

# Subset options
if subset_input_data.lower() in no_validation_list:
    dbiw = []
    maskfile = r""
    mask = []
    filedem = r""
    dbec = []
if subset_input_data.lower() in yes_validation_list:
    dbiw = []
    maskfile = AOI_vector_file
    mask = [AOI_segment_number]
    filedem = DEM_file
    dbec = [DEM_elevation_channel]

calibtyp = Calibration_type
if generate_overviews.lower() in yes_validation_list:
    poption = 'AVER'
else:
    poption = 'OFF'

dblayout = 'TILED256'
fillop = subset_filling_option

IngestedFilesPath = []
Ref_count = 1
Dep_count = 1
for pair, ref_in, dep_in, ref_out, dep_out in zip(pair_number, Reference_file, Dependent_file,
                                                  ReferenceFiles_output_FilenamesList,
                                                  DependentFiles_output_FilenamesList):

    ref_filo = os.path.join(fld_scenes_import, ref_out)
    if os.path.exists(ref_filo):
        print("\t")
        print("\t")
        print((time.strftime("%H:%M:%S")) + " Ingesting reference file, pair " + str(Ref_count) + " of " +
              str(len(pair_number)) + "...")
        print("   Reference file already ingested")
        print("   " + ref_filo)
        print("\t")
        Ref_count = Ref_count + 1

    else:
        ingestref_start = time.time()
        print("\t")
        print((time.strftime("%H:%M:%S")) + " Ingesting reference file pair " + str(Ref_count) + " of " +
              str(len(pair_number)) + "...")
        print("   " + ref_in)

        # Sentinel-1 swath selection if requested.
        if SensorModelName2[0] == "SENTINEL_1" and sentinel_swath in [1, 2, 3]:
            subswath = "?t=" + str(sentinel_swath)
            a1 = ref_in
            if a1[1] == ":":
                a1 = "/" + a1
            a1 = a1.replace("\\", "/")
            a1 = parse.quote(a1)
            ref_in = "file:" + a1 + subswath

        fili = ref_in
        filo = ref_filo

        try:
            saringestaoi(fili, filo, calibtyp, dbiw, maskfile, mask,
                         fillop, filedem, dbec, poption, dblayout)
            IngestedFilesPath.append(filo)
            Ref_count = Ref_count + 1
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)

    # DEPENDENT FILES LIST
    dep_filo = os.path.join(fld_scenes_import, dep_out)

    if os.path.exists(dep_filo):
        print("\t")
        print((time.strftime("%H:%M:%S")) + " Ingesting dependent file pair " + str(Dep_count) + " of " +
              str(len(pair_number)) + "...")
        print(dep_in)
        print("   Dependent file already ingested")
        print(dep_filo)
        print("\t")
        Dep_count = Dep_count + 1

    else:
        ingestdep_start = time.time()
        print("\t")
        print((time.strftime("%H:%M:%S")) + " Ingesting dependent file pair " + str(Dep_count) + " of " +
              str(len(pair_number)) + "...")
        print("   " + dep_in)

        # Sentinel-1 swath selection if requested.
        if SensorModelName2[0] == "SENTINEL_1" and sentinel_swath in [1, 2, 3]:
            subswath = "?t=" + str(sentinel_swath)
            a1 = dep_in
            if a1[1] == ":":
                a1 = "/" + a1
            a1 = a1.replace("\\", "/")
            a1 = parse.quote(a1)
            dep_in = "file:" + a1 + subswath

        fili = dep_in
        filo = dep_filo

        try:
            saringestaoi(fili, filo, calibtyp, dbiw, maskfile, mask,
                         fillop, filedem, dbec, poption, dblayout)

            IngestedFilesPath.append(filo)

            Dep_count = Dep_count + 1
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)

# Creating a list of ingested files that can be reused later.
file = os.path.join(fld_scenes_import, prefix +
                    "03_Ingested_pairs_list_for_coregistration.txt")
with open(file, "w") as f:
    for pair, ref_ingested, ref_date, dep_ingested, dep_date in \
       zip(pair_number,
           ReferenceFiles_output_FilenamesList,
           Reference_acquisition_date,
           DependentFiles_output_FilenamesList,
           Dependent_acquisition_date):

        ref_out = os.path.join(fld_scenes_import, ref_ingested)
        dep_out = os.path.join(fld_scenes_import, dep_ingested)

        string1 = (pair + ";" + ref_out + ";" + ref_date + ";" + dep_out +
                   ";" + dep_date)

        f.write("%s\n" % string1)

SaringestAOI_stop = time.time()
string1 = "2- Data Ingestion proc time;"
ellapse_time = str(round(SaringestAOI_stop - SaringestAOI_start, 2))
string2 = string1 + ellapse_time
time_log.write("%s\n" % string2)


print("\t")
print("--------------------------------------------------------------------------------------------------------------")
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
string1 = "4-Total proc time (seconds);" + str(ellapse_time_seconds)
time_log.write("%s\n" % string1)
time_log.close()

