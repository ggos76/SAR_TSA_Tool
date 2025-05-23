#!/usr/bin/env python
'''-----------------------------------------------------------------------------------------------------------------
 * Gabriel Gosselin, CCRS.  2024-2025                                                                              -
 * -----------------------------------------------------------------------------------------------------------------
'''
# -----------------------------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# -----------------------------------------------------------------------------------------------------------------
# A) Input Files
Ingested_pairs_list_for_coregistration = r"E:\test_25\1_Scenes_Import\t25_03_Ingested_pairs_list_for_coregistration.txt"
output_folder = r"E:\test_25"
prefix = "t25_"

# B) Coregistration options.
# Select a co-registered file for line/column consistency for subsequent times
assign_specific_coreg_ref_file = "no"
coregistration_reference_file = r"\TSA_tests\3_COREG\Coregistered_Pairs_Report.txt"

reference_channel = 1
dep_output_channels = [1,2]
dep_output_label = ["RH", "RV"]
number_of_GCPs = 500
minimum_correlation_score = 0.72
search_radius_pixels = 150

# C) Other Options
generate_overviews = "yes"    # Valid optiopns are "yes" or "no"

# Behaviour when output file exists 
if_file_exists = "regenerate"     # Valid options are "skip" or "regenerate"
# -----------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
1) A fix needs to be done for assign_specific_coreg_ref_file, does not work when yes. 
   This concept needs to be rethought.

'''
# ---------------------------------------------------------------------------------------------
#  Main program
# ----------------------------------------------------------------------------------------------
import sys
import os
import fnmatch
import time
import locale
import glob
import shutil

import pci
from pci.inscoreg import inscoreg
from pci.pyramid import pyramid
from pci.exceptions import PCIException
from pci.api import datasource as ds

from TSA_utilities.SAR_TSA_utilities_definitions import nan_replace
from TSA_utilities.SAR_TSA_utilities_definitions import file_size_check
from TSA_utilities.SAR_TSA_utilities_definitions import get_folder_proctime_and_size
from TSA_utilities.SAR_TSA_version_control import version_control 

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -------------------------------------------------------------------------
#  PART 2: Input Parameters validation
# -------------------------------------------------------------------------
# Hardcoded parameters
TSA_02_start = time.time()
val_start = time.time()

yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list + no_validation_list
GB = 1073741824

# Version control
vs_catalyst = pci.version
vs_python = sys.version_info[:3]
version_control (vs_catalyst, vs_python)

# A) Input files
if not os.path.exists(Ingested_pairs_list_for_coregistration):
    print ("Error - The Ingested_pairs_list_for_coregistration file does not exist or the path/filename is wrong.")
    sys.exit()
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# B) Coregistration options (from INSCOREG)
if assign_specific_coreg_ref_file.lower() in no_validation_list:
    use_explicit_ref = False

   # coregistration_reference_file

'''
perform_coregistration = perform_coregistration.lower()
if perform_coregistration in yes_validation_list:
    use_explicit_ref = False
    perform_coregistration = True
    
elif perform_coregistration in no_validation_list:
    perform_coregistration = False
    if os.path.isfile(coregistered_pairs_report):
        print ("The coregistered file is found")
    else:
        print ("Error - The coregistered file report does not exist or the path/filename is wrong")
        sys.exit()
else:
    print ('Error - perform_coregistration must be set to "yes" or "no"')
    sys.exit()
'''

if minimum_correlation_score < 0.72 or minimum_correlation_score > 1:
    print("Error - The minimum correlation score must be between 0.72 and 1")
    print("Specified value:"+str(minimum_correlation_score))
    sys.exit()

if len(dep_output_channels) != len (dep_output_label):
    print ("Error - the number of dependent output channel is different than the number of labels")
    sys.exit()

# A.4)Creating the output folders
Fld_Coregistration = os.path.join(output_folder, "2_Coregistered_Scenes")
if not os.path.exists(Fld_Coregistration):
    os.makedirs(Fld_Coregistration)

# C) Other options
if generate_overviews in yes_validation_list: 
    generate_overviews = True

if if_file_exists.lower() not in ["skip","regenerate"]:     
    print('Error - valid options for existing_data are "skip" or "regenerate"')
    sys.exit()
else: 
    info_message_skip = ("   output file already exists - skip (if_file_exists = skip)")
    info_message_regn =  ("   output file already exists - regenerate (if_file_exists = regenerate)")


# All validations have suceeded, a time log file is open.
val_stop = time.time()
script2_procTime = os.path.join(output_folder, prefix + "TSA_part_02_Scenes_coregistration_processing_time.txt")
time_log = open(script2_procTime, "w")

current_time =  time.localtime()
string_0 = time.strftime("%Y-%m-%d", current_time)
time_log.write("%s\n" % string_0)

string_0 = ("Process;proc.time (secs);Data size (MB);Number of files")
time_log.write("%s\n" % string_0)

val_stop = time.time()
string_1 = "Validation: ;" + str (round((val_stop - val_start), 2)) 
time_log.write("%s\n" % string_1)


# -----------------------------------------------------------------------------------------------------------------
# Scenes coregistration
# -----------------------------------------------------------------------------------------------------------------

print("\t")
print("----------------------------------------------------------------------------------------------------------")
print("                                          Coregistration                                                  ")
print("----------------------------------------------------------------------------------------------------------")

proc_start_time = time.time()

# B.1) Read and parse the baseline report
Input_file_lines = []
pair_number = []
Reference_file = []
Reference_acquisition_date = []
Dependent_file = []
Dependent_acquisition_date = []
pair_number2 = []
Reference_acquisition_date2 = []
Dependent_acquisition_date2 = []

print(time.strftime("%H:%M:%S") + " Reading and parsing the saringest baseline report")
print("\t")

with open(Ingested_pairs_list_for_coregistration, "r") as input_file:
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

pairs = (len(pair_number)) + 1
print("\t")

# B.2) Coregistration
last_coregistered_list = []
dep_resampled_files_list = []
count = 1
new_Reference_file = []

# Restructure ref/dep lists to accomodate explicit reference, if using
if use_explicit_ref:
    # Add the original first ref file as the first dependent file
    Dependent_file.insert(0, Reference_file[0])
    Dependent_acquisition_date.insert(0, Dependent_acquisition_date[0])
    # Add the explicit reference file as first ref file
    Reference_file.insert(0, explicit_coreg_reference_file)
    explicit_coreg_date = explicit_coreg_reference_file[-12:-4]
    Reference_acquisition_date.insert(0, explicit_coreg_date)

for pair, input_ref, input_dep, date_ref, date_dep in zip(pair_number, Reference_file, Dependent_file,
                                                          Reference_acquisition_date, Dependent_acquisition_date):
    coreg_start = time.time()

    print("\t")
    print("---------------------------------------------------------------------------------------")
    print(((time.strftime("%H:%M:%S")) + " Processing pair " +
        str(count)+" of "+str(len(Reference_file))))
    print("   Original reference file:" + input_ref)
    print("   Original dependent file:" + input_dep)
    print("\t")

    # B.2.1) Take the reference file of the first pair (pair = 0) and create a copy of it the coreg folder. 
    # This file is not resampled but its size (pixels, columns, i.e. the geographical frame) will determine 
    # the size of all other coregistered files
    print("STEP 1")
    if pair == "0":
        print("STEP 2")
        outbasename = os.path.basename(input_ref)
        outputFilePath = os.path.join(Fld_Coregistration, "coreg_" + str(pair) + "_" + outbasename)
        shutil.copy2(input_ref, outputFilePath)
        print("   (copy):"+outputFilePath)

    # B.2.2) Check if a resampled version of the reference file exist
    # Any resampled dependent file that can be used as reference file
    # will be named coreg_pair#_prefix_YYYYMMDD.
    # Find in the output folder the file that corresponds to *date_dep
    # and use it as the reference file

    pathfile = os.path.join(Fld_Coregistration, "coreg" + "*" + date_ref + ".pix")
    resamp_ref_list = []  # If a file exist - Go to step 9
    print("STEP 3")
    for ref_resampled in glob.glob(pathfile):
        resamp_ref_list.append(ref_resampled)

    # The resampled version of the reference file list is empty We need to coregister the reference file to
    #  a previously coregistered file.
    if not resamp_ref_list:
        print("STEP 4")
        print("\t")
        print("   (A) Reference file coregistration")

        print("     (A.1) Resampled reference file does not exist")
        print("     (A.2) Creating resampled reference file")
        # This can happen if there is a GAP in the time series when
        # SBAS is used. For example A*B, B*C, C*D, E-F, F*G, G*H
        # In such case we need to connect D-E to temporally create
        # a coregistered file for D-E to ensure the same number of
        # lines and columns. We need to use the last resampled file.

        # retrieve the original reference file from SARINGESTAOI

        print("\t")
        for last_coregistered in reversed(last_coregistered_list):
            print("STEP 5")
            # The last_coregistered_list contains the sucessfully
            # coregistered pairs and others that failed to coregister
            # We check first if the coreg file exists to avoid creating
            # another problem
            print("  checking if " + last_coregistered + " exists")
            if os.path.isfile(last_coregistered):
                print("STEP 6")
                filref = last_coregistered
                dbic_ref = [reference_channel]

                # !!!!!!We temporarilly place the reference file to the dependent
                #        position to create a resampled version of it

                fili = input_ref
                dbic_dep = [reference_channel]
                dbic = dep_output_channels
                numpts = [number_of_GCPs]
                minscore = [minimum_correlation_score]
                searchr = [search_radius_pixels]
                filo = os.path.join(Fld_Coregistration, "coreg_" + str(pair) + "_" +
                                    os.path.basename(input_ref))
                print("    Temp reference file: " + filref)
                print("    Temp dependent file: " + fili)
                print("    Filo: " + filo)
                try:
                    inscoreg(filref=filref, dbic_ref=dbic_ref, fili=fili,
                                dbic_dep=dbic_dep, dbic=dbic, numpts=numpts,
                                minscore=minscore, searchr=searchr, filo=filo)
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)

            # Check if FILO exists, if yes it means that the coregistration worked.
            # if not we need to try again with the second last coregistered file.
                if os.path.isfile(filo):
                    print("STEP 8")
                    print("\t")
                    print("  Temporary reference file sucessfully coregistered ")
                    break
                else:
                    print("STEP 7")
                    print("ERROR 1! NO GCPs found for the temporary coregistered file")
                    print("Trying again with the previous coregistered file in the last_coregistered_list")
                    print("\t")

        print("\t")
        print("  New resampled reference file: " + filo)
        filref = filo    #

    else:
        print("STEP 9")
        filref = resamp_ref_list[0]
        print("   2-Using resampled reference file " + filref)

    # B.2.3) Main coregistration
    print("\t")
    print("   (B) Main coregistration for dependent file")

    print("    B.1) Reference file:"+filref)
    print("    B.2) Dependent file:"+input_dep)

    dbic_ref = [reference_channel]
    fili = input_dep
    dbic_dep = [reference_channel]
    dbic = dep_output_channels
    numpts = [number_of_GCPs]
    minscore = [minimum_correlation_score]
    searchr = [search_radius_pixels]

    filo = os.path.join(Fld_Coregistration, "coreg_" +
                str(pair) + "_" + os.path.basename(input_dep))
    print("    B.3) Output dependent resampled file: " + filo)
    print("\t")
    # to be used if there is a gap in the time series, the last
    # dependent file will be the closest in time to the next reference file

    print("STEP 10")
    try:
        inscoreg(filref=filref, dbic_ref=dbic_ref, fili=fili,
                    dbic_dep=dbic_dep, dbic=dbic, numpts=numpts,
                    minscore=minscore, searchr=searchr, filo=filo)

        if generate_overviews is True:
            pyramid(file=filo, force='no', poption='aver',
                    dboc=[], olevels=[])
    except PCIException as e:
        print(e)
    except Exception as e:
        print(e)
    count = count + 1

    if os.path.isfile(filo):
        print("STEP 12")
        print("--------------------------------------------------------")
        last_coregistered_list.append(filo)

        # Output coregistered pairs report. Record only sucesfully
        #  coregistered fikes
        new_Reference_file.append(filref)
        pair_number2.append(pair)
        Reference_acquisition_date2.append(date_ref)
        dep_resampled_files_list.append(filo)
        Dependent_acquisition_date2.append(date_dep)

        # Output text file.
        coreg_stop = time.time()
        string1 = " Coregistration pair " + str(pair) + ";"
        ellapse_time = str(round(coreg_stop - coreg_start, 2))
        fsize = round((os.path.getsize(filo) / GB), 3)
        string2 = string1 + ellapse_time + ";" + str(fsize)
        time_log.write("%s\n" % string2)

    else:
        print("STEP 11")
        print("ERROR 2! NO GCPs found for this INSCOREG file")
        print("SKIPPING TO THE NEXT PAIR")
        print("\t")
        print("\t")

# B.3) Update the list of coregistered pairs for INSRAW
if use_explicit_ref:
    # Remove the first pair that was added earlier
    pair_number2.pop(0)
    new_Reference_file.pop(0)
    Reference_acquisition_date2.pop(0)
    dep_resampled_files_list.pop(0)
    Dependent_acquisition_date2.pop(0)
print("Output files list size verification. ")
print(len(pair_number2))
print(len(new_Reference_file))
print(len(Reference_acquisition_date2))
print(len(dep_resampled_files_list))
print(len(Dependent_acquisition_date2))

coregistration_report = os.path.join(Fld_Coregistration, prefix + "04_Coregistered_Pairs_Report.txt")
coreg_report_when_true = coregistration_report
with open(coregistration_report, "w") as f:
    for pair, ref_file, ref_date, dep_file, dep_date in \
            zip(pair_number2, new_Reference_file, Reference_acquisition_date2,
                dep_resampled_files_list, Dependent_acquisition_date2):
        string1 = (pair + ";" + ref_file + ";" + ref_date +
                    ";" + dep_file + ";" + dep_date)
        f.write("%s\n" % string1)



proc_stop_time = time.time()
folder = Fld_Coregistration
out_folder_time_size, size_mb = get_folder_proctime_and_size (folder, proc_stop_time, proc_start_time)
string_1 = ("Input scenes coregistration:" + out_folder_time_size) 
time_log.write("%s\n" % string_1)


# B.4 - Safety check. All coregistered files are now supposed to have the same number of lines and columns as per design.
# This is essential for PSSARTSA.  If it's not the case something went wrong with SARINGESTAOI\INSCOREG.
files_list = []
for root, dirs, files in os.walk(Fld_Coregistration):
    for filename in fnmatch.filter(files, "*.pix"):
        files_list.append(os.path.join(root, filename))

file_size_check (files_list)


# Recreate a safery copy of the coregistered files. 

print("-----------------------------------------------------------------------------------------------------------------")
print("-----------------------------------------------------------------------------------------------------------------")
print((time.strftime("%H:%M:%S")))
print("All processing completed")
print("\t")

TSA_02_stop = time.time()
ellapse_time_seconds = round((TSA_02_stop - TSA_02_start), 2)
ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

print("Processing time (seconds): " + str(ellapse_time_seconds))
print("Processing time (minutes): " + str(ellapse_time_minutes))
print("Processing time (hours): " + str(ellapse_time_hours))

string1 = "TSA_02 total processing time (secs):;" + str(ellapse_time_seconds)
time_log.write("%s\n" % string1)
string1 = "TSA_02 total processing time (mins):;" + str(ellapse_time_minutes)
time_log.write("%s\n" % string1)
string1 = "TSA_02 total processing time hours):;" + str(ellapse_time_hours)
time_log.write("%s\n" % string1)
time_log.close()