#!/usr/bin/env python
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. May 2024                                      -
 * ----------------------------------------------------------------------
need to explore the following GPU wrapper for numpy
https://stackoverflow.com/questions/49605231/does-numpy-automatically-detect-and-use-gpu
look for cupy in particular
https://cupy.dev/
'''

# ----------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# ----------------------------------------------------------------------------------------------
# A) Input Files
stack_lists_folder = r"E:\QC03M_7224_4504_CP01_A\6_TSA_stack_lists"
prefix = "QC03M_7224_4504_CP01_A_"
no_data_value = -32768.0000

use_stack_subset = "yes"               # Valid option are "yes" or "no"
sub_stack_start_date = 20250308        # must be an interger in YYYYMMDD format
sub_stack_stop_date = 20250425
suffix = "_sl01"                       # optionnal- suffix for output stacks, leave blank for no suffix


#B.1) Time series preparation and data stacking
stacking_method = "catalyst"             # Valid option are "numpy or "catalyst"
output_type = 2                          # 1: analysis layers only    2: analysis and stack together   3: analysis and stack separate

# B.2) Filtering and mask options
apply_masking = "no"       # Valid option are "yes" or "no"
mask_type = "exclusion"    # Valid option are "inclusion" or "exclusion"
mask_file = r"\\W-BSC-A157283\share_GG\CanUS_border_1m_UTM18T_D000_v4_2km_buffer.pix"
mask_seg_number = [2]

apply_mean_filter = "yes"               # Valid option are "yes" or "no"
filter_size_X_Y = [5,5]                 # Tuple of odd integer

#C) Other options
# Generate overviews - either yes or no,
generate_overviews = "yes"

# Behaviour when output file exists 
if_file_exists = "skip"     # Valid options are "skip" or "regenerate"
# -----------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
2025-06-06 - The MUMPY solution needs updates to be working again. 

1) There is no verification if the specified DEM covers completely the spatial
   extents of the AOI. If not the script will run to completion but all
   interferograms will be blank after INSRAW.

   
 2) stacking_method
    2.1:  numpy. Doen't rely on Catalyst tools but could be slow and crash with very large datasets.
    2.2:  catalyst : Use PSSARTSA to stack intensity layers and inscohstats to stack coherence layers. 
                     These algotithms are not always available depending on the license type. 
   
    The two stacking method gives similar but not identical results. 
'''
# -----------------------------------------------------------------------------------------------------------------
#  Imports
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
from pci.fexport import fexport
from pci.fav import fav
from pci.iii import iii
from pci.pyramid import pyramid
from pci.api import gobs
from pci.api import cts
from pci.exceptions import PCIException
from pci.api import datasource as ds

from TSA_utilities.SAR_TSA_utilities_definitions import file_size_check
from TSA_utilities.SAR_TSA_utilities_definitions import stack_masking
from TSA_utilities.SAR_TSA_utilities_definitions import stack_min_max
from TSA_utilities.SAR_TSA_utilities_definitions import get_folder_proctime_and_size
from TSA_utilities.SAR_TSA_version_control import version_control

import numpy as np
locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -----------------------------------------------------------------------------------------------------------------
#  Parameters Validation
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
if not os.path.exists(stack_lists_folder):
    print ("Error - The stack_lists_folder does not exist or the path/filename is wrong.")
    sys.exit()
else: 
    Fld_stack_lists = stack_lists_folder

# Create a list of stack_list to proccess
input_stack_lists = []
for root, dirs, files in os.walk(stack_lists_folder):
    for filename in fnmatch.filter(files, "*stack_statistics.txt"):
        input_stack_lists.append(os.path.join(root, filename))

if len(input_stack_lists) == 0:
    print ("Error - No stack_statistics.txt file found in the stack_lists_folder")
    sys.exit()

use_stack_subset = use_stack_subset.lower()
if use_stack_subset not in yes_no_validation_list:
    print ('Error - the use_stack_subset option must be "yes" or "no"')
    sys.exit()
elif use_stack_subset in yes_validation_list: 
    use_stack_subset = True
    #-------------------------------------------------------------------------------------
    # Validation for the dates here. 
    #-----------------------------------------------------------------------------------------
else: 
    use_stack_subset = False

#B) Time series preparation
#B.1A) stacking_method
stacking_method = stacking_method.lower()
if stacking_method not in ["numpy", "catalyst"]: 
    print ('Error - the staking_method options are "numpy" or "catalyst"')
    sys.exit()
elif  stacking_method == "numpy":
    use_numpy = True
    use_catalyst = False
elif  stacking_method == "catalyst":
    from pci.pssartsa import pssartsa
    use_numpy = False
    use_catalyst = True
else:
    print ('Error - Undefine error')
    sys.exit()

#B.1B) Output type
if output_type not in [1,2,3]: 
    print ("Error - the output type must ne 1, 2 or 3")
    sys.exit()

# B.2) Mean filter
if apply_mean_filter.lower() in yes_validation_list:
    apply_mean_filter = True

    Fsize_X_modulo = filter_size_X_Y[0] % 2
    Fsize_Y_modulo = filter_size_X_Y[1] % 2
    if Fsize_X_modulo != 1 or Fsize_Y_modulo!= 1 :
        print ("Error - The filter_size_X_Y must only be compose of odd integers ")
        sys.exit()

    out_filter = ("_" + str(filter_size_X_Y[0]) + "x" + str(filter_size_X_Y[0]))
else:
    apply_mean_filter = False
    out_filter = ""

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

# C) Other options
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

# D) Creating the output folder and the output file name.  
    # Note: In anticipation for full automatation, we will create automatically the output folder and 
    # the output file name.      
two_up = Path(input_stack_lists[0]).resolve().parents[1]
stack_folder = "7_TSA_stacks"

Fld_output_stacks = os.path.join (two_up, "7_TSA_stacks")
if not os.path.exists(Fld_output_stacks):
    os.makedirs(Fld_output_stacks)

# All validations have suceeded, a time log file is open.
out_folder_script4 = two_up
script4_procTime = os.path.join(out_folder_script4, prefix + "TSA_part_04_stack_data_analysis_processingTime.txt")
time_log = open(script4_procTime, "w")

val_stop = time.time()
string_1 = "Validation step (secs): " + str (round((val_stop - val_start), 2)) 
time_log.write("%s\n" % string_1)


# -----------------------------------------------------------------------------------------------------------------   
# -----------------------------------------------------------------------------------------------------------------   
print("\t")
print("-----------------------------------------------------------------------------------------------------------")
print("                                          Main program                                                     ")
print("-----------------------------------------------------------------------------------------------------------")
print("\t")

print(time.strftime("%H:%M:%S") + " Finding the input stack lists to process")
input_stack_files = []
for root, dirs, files in os.walk(Fld_stack_lists):
    for filename in fnmatch.filter(files, "*stack_statistics.txt"):
        input_stack_files.append(os.path.join(root, filename))

print("   Number of stack to process: " + str(len(input_stack_files)))
for ii in input_stack_files:
    print ("   " + ii)

# -----------------------------------------------------------------------------------------------------------------
# Stack processing
# -----------------------------------------------------------------------------------------------------------------
print ("\t")
nb_files = str(len(input_stack_files))
count = 1
for in_stack_file in input_stack_files: 
    proc_start_time = time.time()
    print ("\t")
    time_log.write("%s\n" % "--------------------------------------------------------------------------------------")
    string_1 = ("Processing stack " + str (count) + " of " +  nb_files)
    time_log.write("%s\n" % string_1)
    time_log.write("%s\n" % in_stack_file)
    
    print(time.strftime("%H:%M:%S") + " " + string_1)
    print ("   " + in_stack_file)

    print ("   reading and parsinng the input stack file")
    Input_file_lines = []
    with open(in_stack_file, "r") as input_file:
        for line in input_file:
            Input_file_lines.append(line.strip())

    input_file_list = []
    mid_date_list = []
    ref_date_list = []
    dep_date_list = []

    for ii in Input_file_lines:
        temp = ii.split(';')
        input_file_list.append(temp[0])
        mid_date_list.append(temp[1])
        ref_date_list.append(temp[2])
        dep_date_list.append(temp[3])

    input_file_list.pop(0)
    mid_date_list.pop(0)
    ref_date_list.pop(0)
    dep_date_list.pop(0)

    nb_files = str(len(input_file_list))
    # Verification
    print("\t")
    print("   " + time.strftime("%H:%M:%S") + " files to process")

    print("   File -----> acquisition date")
    count2 = 1
    for ii, jj in zip(input_file_list, mid_date_list):
        print ("   file " + str(count2) + "|" + nb_files + " " + ii + " -----> " + jj)
        count2 = count2 + 1

    # Creating a subseted version of the stack
    if use_stack_subset is True:
        print ("\t")
        print("   " + time.strftime("%H:%M:%S") + " Extracting a subset of the input scenes")
        sub_input_file_list = []
        sub_mid_date_list = []
        sub_ref_date_list = []
        sub_dep_date_list = []

        print("   Start date: " + str(sub_stack_start_date))
        print("   Stop date: "+ str(sub_stack_stop_date))
        for ii, jj, kk, ll in zip(input_file_list, mid_date_list, ref_date_list, dep_date_list):
            
            if int(kk) >= sub_stack_start_date and int(ll)<= sub_stack_stop_date:
                sub_input_file_list.append(ii)
                sub_mid_date_list.append(jj)
                sub_ref_date_list.append(kk)
                sub_dep_date_list.append(ll)
                print ("   reference date: " + kk + " --> inside range")
            else: 
                print ("   reference date: " + kk + " --> outside range")
        
        # reset the lists of inoput files with their subset
        input_file_list = sub_input_file_list
        mid_date_list = sub_mid_date_list
        ref_date_list = sub_ref_date_list
        dep_date_list = sub_dep_date_list

    # -----------------------------------------------------------------------------------------------------------------
    # File size and projection verifications. All files must have the same number of lines and columns 
    # and the same projection
    # -----------------------------------------------------------------------------------------------------------------
 
    print("\t")
    print ("----------------------------------------------------------------------------------------------------------")
    print("File size verification")
    files_list = input_file_list
    file_size_check (files_list)
    
    string_1 = ("   " + time.strftime("%H:%M:%S") + " input files size verification")
    time_log.write("%s\n" % string_1)
    
    print("\t")
    print ("----------------------------------------------------------------------------------------------------------")
    proc_stop_time = time.time()
   
    # -----------------------------------------------------------------------------------------------------------------
    if apply_masking is True:
        string_1 = ("   " + time.strftime("%H:%M:%S") + " masking the input files (" + mask_type + ")")
        time_log.write("%s\n" % string_1)

        print("Stack data preprocessing - Applying the exclusion or inclusion vector mask")
        input_stack = input_file_list
        stack_masking(input_stack, mask_type, mask_file, mask_seg_number,no_data_value, Fld_output_stacks)
    if apply_masking is False :
        print(time.strftime("%H:%M:%S") + " Stack preparation - Exclusion or inclusion mask not requested")
    print("\t")


    proc_stop_time = time.time()
    # determine if the input data are coherence. 
    if "_coh_" in in_stack_file: 
        data_type_coh = True
    else: 
        data_type_coh = False

    #-----------------------------------------------------------------------------------------------------------------
    # STACKING METHOD USING NUMPY   --  NEED UPDATES, not working. 
    #------------------------------------------------------------------------------------------------------------------
    if use_numpy is True: 

        print(time.strftime("%H:%M:%S") + " Selected stacking method is NUMPY")

        # -------------------------------------------   Filtering   ---------------------------------------------------
        nb_files = str(len(input_file_list))    
        count  = 1
        if apply_mean_filter is True:
            print ("Stack data preprocessing - Average filtering option is selected. Average filter size: " + str (filter_size_X_Y[0]) +
                    "x" + str (filter_size_X_Y[1]))

            for ii in input_file_list: 
                print("   " + time.strftime("%H:%M:%S") + " Filtering file " + str(count) + " of " + str(nb_files))
                file = ii
                dbic = [1]  
                dboc = [1]  
                flsz = [filter_size_X_Y[0], filter_size_X_Y[1]]
                mask = []  
                bgrange = []  
                failvalu = []  # no failure value
                bgzero = ''  # default, set background to zero

                try:
                    fav(file, dbic, dboc, flsz, mask, bgrange, failvalu, bgzero)
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)
                count = count + 1
        # ---------------------------------------------------------------------------------------------------------------------
        # here we assume the list to be in the chronological order since it was created by TSA_02
        print("\t")
        print(time.strftime("%H:%M:%S") + " Copying and moving the first file of the stack")
        first_file = input_file_list[0]
        print ("   " + first_file)

        # Valid for all output type options
        copy_out_data =  os.path.join(outfile_tsa_data + out_filter + ".pix" )
        shutil.copy2(first_file, copy_out_data)

        print ("\t")
        print(time.strftime("%H:%M:%S") + " Transferring the other files in the data stack")

        file_to_transfer = input_file_list
        file_to_transfer.pop(0)
        nb_files = str(len(file_to_transfer))

        print(time.strftime("%H:%M:%S") + " Creating output stack empty channels")
        fsize = round((os.path.getsize(copy_out_data)/ 1073741824), 3)

        total_estimated_size = len(input_file_list) * fsize
        print ("   Estimated size of the output stack: " + str (total_estimated_size) + " GB")
        print ("\t")
        print(time.strftime("%H:%M:%S") + " Adding the empty channels to the stack")
        file = copy_out_data
        pciop = 'ADD'
        pcival = [0,0,0,len(file_to_transfer)]

        try:
            pcimod(file,pciop,pcival)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)

        print ("\t")
        count = 2
        for ii in file_to_transfer:

            print(time.strftime("%H:%M:%S") + " Transfering file " + str(count) + " of " + nb_files)
            print("    " + ii )

            fili = ii  #
            filo = copy_out_data
            dbic = [1]
            dboc = [count]
            print ("    output channel " + str(count))
            dbiw = []
            dbow = []
            options = ""
            try:
                iii(fili, filo, dbic, dboc, dbiw, dbow, options)
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
            count = count + 1

        print("\t")
        print('----------------------------------------------------------------------------------------------------------')
        print("                                      Time series preparation                                             ")
        print('----------------------------------------------------------------------------------------------------------')
        print("\t")

    
        # -------------------------------------------------------------------------------------------------------------
        # Apply_masking if needed and assign to do NoDataValue.
        #  if the mask is set with inclusion, we need to set the pixels outside the mask with No_data
        # if the mask is set with exclusion, we need to set the pixels inside the mask with No_data
   

        # -------------------------------------------------------------------------------------------------------------
        # Renaming the channels according to the mid date and creating a metadata of the same name.
        print("\t")
        print ("\t")
        print(time.strftime("%H:%M:%S") + " Stack preparation - Adding channels description and metadata")

        # open dataset in write mode
        channel = 1
        for mid_date, ref_date, dep_date in  zip (mid_date_list, ref_date_list, dep_date_list):

            # print (mid_date + "-->" + str(channel))
            with ds.open_dataset(copy_out_data, ds.eAM_WRITE) as ds1:
                # get the AuxiliaryData
                aux = ds1.aux_data

                # get the metadata for the first channel
                mdmap = aux.get_chan_metadata(channel)
                # set the NO_DATA_VALUE to -32768
                mdmap['MID_DATE'] = mid_date
                mdmap['REF_DATE'] = ref_date
                mdmap['DEP_DATE'] = dep_date
                mdmap['NO_DATA_VALUE'] = str(no_data_value)
                chan_desc = (mid_date + " (ref_" + ref_date + "_" + "dep_" + dep_date + ")")
                # set the metadata map back to the AuxiliaryData
                aux.set_chan_description(chan_desc, channel)
                aux.set_chan_metadata(mdmap, channel)

                # set the AuxiliaryData back to the dataset
                ds1.aux_data = aux
            channel = channel + 1


        print("\t")
        print('----------------------------------------------------------------------------------------------------------')
        print("                   Time series analysis (Absolute sum of changes and other statistics)                    ")
        print('----------------------------------------------------------------------------------------------------------')
        print("\t")

        with ds.open_dataset(copy_out_data) as ds9:

            print(time.strftime("%H:%M:%S") + " Stack to process: " + copy_out_data)
            reader = ds.BasicReader(ds9)
            # read the raster channels
            raster = reader.read_raster(0, 0, reader.width, reader.height)
            num_channels = ds9.chan_count
            #print (num_channels)
            ref_crs = ds9.crs                # coordinate system
            ref_geocoding = ds9.geocoding    # Geocoding


        # initialising the output raster file that will receive the data
        print("\t")
        print(time.strftime("%H:%M:%S") + " Initializing the output file ")
        output_ts_raster = np.zeros([reader.width, reader.height,8], dtype="float32")

        # [0] = sum of absolute difference (between subsequent pairs)
        # [1] = average value
        # [2] = standard-deviation
        # [3] = maximum value
        # [4] = minimum value
        # [5] = Span (max - min)
        # [6] = (Max value / sum) * 100
        # [7] = Index for the maximum value

        print("\t")
        print(time.strftime("%H:%M:%S") + " Computing the time series statistics")
        count = 1
        for x in range(reader.width):
            for y in range(reader.height):

                val = raster.data[x, y, :]
                if no_data_value in val:
                    output_ts_raster[x, y, 0] = no_data_value
                    output_ts_raster[x, y, 1] = no_data_value
                    output_ts_raster[x, y, 2] = no_data_value
                    output_ts_raster[x, y, 3] = no_data_value
                    output_ts_raster[x, y, 4] = no_data_value
                    output_ts_raster[x, y, 5] = no_data_value
                    output_ts_raster[x, y, 6] = no_data_value
                    output_ts_raster[x, y, 7] = no_data_value
                else:
                    #val = raster.data[x, y, :]
                    abs_diff = [abs(t - s) for s, t in zip(val, val[1:])]
                    max_t = np.max(val)
                    min_t = np.min(val)
                    sum_t = np.sum(val)
                    output_ts_raster[x, y, 0] = np.sum(abs_diff)
                    output_ts_raster[x, y, 1] = np.average(val)
                    output_ts_raster[x, y, 2] = np.std(val)
                    output_ts_raster[x, y, 3] = max_t
                    output_ts_raster[x, y, 4] = min_t
                    output_ts_raster[x, y, 5] = max_t - min_t
                    pct_contribution = ((max_t/sum_t) * 100)
                    output_ts_raster[x, y, 6] = pct_contribution
                    output_ts_raster[x, y, 7] = (np.argmax(val) + 1)

            sys.stdout.write("\r" +"Progress: " + str(round((count/reader.width) * 100,1)) +" %")
            sys.stdout.flush()
            count = count + 1

        print("\t")
        print("\t")
        print(time.strftime("%H:%M:%S") + " Exporting the result to a PCIDSK file")
    
        out_raster_analysis = gobs.copy_array_to_raster(output_ts_raster)
        out_tsa_analysis_pix =  os.path.join(outfile_tsa_analysis + out_filter + ".pix" )

        with ds.new_dataset(out_tsa_analysis_pix, 'PCIDSK', '') as write_ts:
            writer =  ds.BasicWriter(write_ts)     # create a writer to write the raster
            writer.create(out_raster_analysis)                # create the file on disk based on the size and datatype of raster
            writer.write_raster(out_raster_analysis)          # write the raster to the new file
            writer.crs = ref_crs                    # write the coordinate system
            writer.geocoding = ref_geocoding        # write the geocoding information / must be of same

        print("\t")
        print(time.strftime("%H:%M:%S") + " Writing  metadata to the output file")
        with ds.open_dataset(out_tsa_analysis_pix, ds.eAM_WRITE) as ds1:
            # get the AuxiliaryData
            aux = ds1.aux_data
            # Add the NO_DATA_VALUE at the file level.
            metadata = aux.file_metadata
            metadata['NO_DATA_VALUE'] = str(no_data_value)
            aux.file_metadata = metadata

            # set the metadata map back to the AuxiliaryData
            aux.set_chan_description("Absolute sum of changes", 1)
            aux.set_chan_description("Mean", 2)
            aux.set_chan_description("Standard deviation", 3)
            aux.set_chan_description("Maximum", 4)
            aux.set_chan_description("Minimum", 5)
            aux.set_chan_description("Span (max - min)", 6)
            aux.set_chan_description("Pct(%) contribution of max value", 7)
            aux.set_chan_description("Date of max value", 8)
            # set the AuxiliaryData back to the dataset (ds1)
            ds1.aux_data = aux

        #  Post processing depending on the stack output_type
        if output_type == 1:    # We need to delete the  *TSA_stack_data_.pix file
            print ("\t")
            print(time.strftime("%H:%M:%S") + " output_type 1 selected, deleting the stack data file")
            os.remove (copy_out_data)
        if output_type == 2:    # We need to transfert the *TSA_stack_data_.pix file into the *TSA_stack_analysis*.pix file 

            print(time.strftime("%H:%M:%S") + " output_type 2 selected, preparing the files.")
            with ds.open_dataset(copy_out_data) as ds10:

                num_channels = ds10.chan_count  # number of layers in the data stack
                file = out_tsa_analysis_pix
                pciop = 'ADD'
                pcival = [0,0,0,num_channels]   # Add the equivalent of empty channels to the statistic file. 

                try:
                    pcimod(file,pciop,pcival)
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)

            for ii in range (1, num_channels + 1): 

                fili = copy_out_data
                filo = out_tsa_analysis_pix
                dbic = [ii]
            
                out_chan =  8 + ii
                dboc = [out_chan]
                dbiw = []
                dbow = []
                options = ""
                try:
                    iii(fili, filo, dbic, dboc, dbiw, dbow, options)
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)
            
            new_tsa_analysis_name =  (outfile_tsa_analysis_stack + out_filter + ".pix")            
            os.rename(out_tsa_analysis_pix, new_tsa_analysis_name)
            os.remove (copy_out_data)
    #------------------------------------------------------------------------------------------------------------------------
    # STACKING METHOD USING CATALYST
    #------------------------------------------------------------------------------------------------------------------------
    if use_catalyst is True:
        
        string_1 = ("   " + time.strftime("%H:%M:%S") + " Selected stacking method is CATALYST")
        print (string_1)
        time_log.write("%s\n" % string_1)

       # 1) We need to create the output file name
        # hack to make coherence data works with PSSARTSA to ensure time ordered channels. 
        if data_type_coh is True:  
           # We create the "Acquisition DateTime" metadata and fili it with the  "Dep_Acquisition_DateTime"
           for ii in input_file_list:  
       
                with ds.open_dataset(ii,ds.eAM_WRITE) as ds12:   # Note .eAM_WRITE
                    aux = ds12.aux_data
                    Dep_DateTime = aux.get_file_metadata_value('Dep_Acquisition_DateTime')
                
                    mdmap = aux.file_metadata
                    mdmap['Acquisition_DateTime'] = Dep_DateTime
                    aux.file_metadata = mdmap
                    ds12.aux_data = aux

        # preparing the mfile
        temp_mfile = os.path.join(Fld_output_stacks, "temp_mfile_for_data_ingestion.txt")
        mfile_input = open(temp_mfile, "w")
        mfile_input.write('\n'.join(input_file_list))
        mfile_input.close() 

        if apply_mean_filter is True: 
            flsz = filter_size_X_Y
        if apply_mean_filter is False: 
            flsz = []

        mfile =  temp_mfile
        dbic = [1]	
        mask =	[]			 					
        maskfile =	''		

        base_stack = os.path.basename(in_stack_file[:-21])

        if output_type == 1: # analysis layers only
            print("   " + time.strftime("%H:%M:%S") + " Output type 1: analysis layers only")
            stack =  "no"							
            output_file =  os.path.join(Fld_output_stacks, base_stack + out_filter + "_TSA_stack_analysis" + suffix + ".pix")
        if output_type == 2:  # analysis and stack together
            print("   " + time.strftime("%H:%M:%S") + " Output type 2: analysis and data stack in the same file")
            stack =  "yes"
            output_file =  os.path.join(Fld_output_stacks, base_stack + out_filter + "_TSA_stack_analysis_and_data" + suffix  + ".pix")
        if output_type == 3:  # analysis and stack separate
            print("   " + time.strftime("%H:%M:%S") + " Output type 3: analysis and data stack in separate files")
            stack =  "yes"
            output_file = os.path.join(Fld_output_stacks, base_stack + out_filter + "_TSA_type3_TEMP.pix")
            out_type3_analysis =  os.path.join(Fld_output_stacks, base_stack + out_filter + "_TSA_stack_analysis" + suffix + "pix")
            out_type3_data =  os.path.join(Fld_output_stacks, base_stack + out_filter + "_TSA_stack_data" + suffix  + ".pix")
            
        print ("   Output file: " + output_file)
        filo = output_file
        
        string_1 = ("   " + time.strftime("%H:%M:%S") + " creating the TSA file")
        time_log.write("%s\n" % string_1)

        if os.path.exists (filo) and if_file_exists == "skip": 
            print (info_message_skip)
        else: 
            if os.path.exists (filo) and if_file_exists == "regenerate": 
                print (info_message_regn)
                os.remove (filo)
            try:
                pssartsa(mfile, dbic, mask, maskfile, stack, flsz, filo)
                if generate_overviews is True :
                    pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
        
        # -----------------------------------------------------------------------------------------------------------------
        # Split the TSA data analysis and the TSA data stack
        if output_type == 3: 
            
            with ds.open_dataset(output_file) as ds5:
                nb_chans = ds5.chan_count

            string_1 = ("   " + time.strftime("%H:%M:%S") + " Creating the TSA analysis file")
            time_log.write("%s\n" % string_1)
            print (string_1)
           
           # We export the fourteen first layers to create the idependent data analysis file. 
            fili =	output_file
            filo =	out_type3_analysis
            dbiw =	[]
            dbic =	[1,-14]
            dbib =	[]
            dbvs =	[]
            dblut =	[]
            dbpct =	[]
            ftype =	"PIX"
            foptions = ""

            if os.path.exists (filo) and if_file_exists == "skip": 
                print (info_message_skip)
            else: 
                if os.path.exists (filo) and if_file_exists == "regenerate": 
                    print (info_message_regn)
                    os.remove (filo)
                try:
                    fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions)
                    if generate_overviews is True :
                        pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)
        

            string_1 = ("   " + time.strftime("%H:%M:%S") + " Creating the TSA data file")
            time_log.write("%s\n" % string_1)
            print (string_1)
            # The we export the other stack data layers. 
            fili =	output_file
            filo =	out_type3_data
            dbiw =	[]
            dbic =	[15,-nb_chans]
            dbib =	[]
            dbvs =	[]
            dblut =	[]
            dbpct =	[]
            ftype =	"PIX"
            foptions = ""

            if os.path.exists (filo) and if_file_exists == "skip": 
                print (info_message_skip)
            else: 
                if os.path.exists (filo) and if_file_exists == "regenerate": 
                    print (info_message_regn)
                    os.remove (filo)
                try:
                    fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions)
                    if generate_overviews is True :
                        pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)

            # We delete the initial file. 
    

        proc_stop_time = time.time()
        ellapse_time_seconds = str(round(( proc_stop_time -  proc_start_time), 2))
        string_1 = ("   processing time (sec): " + ellapse_time_seconds)
        time_log.write("%s\n" % string_1)
    count = count + 1
# -----------------------------------------------------------------------------------------------------------------
# Deletintg the temporary file
del_files_list = []
for root, dirs, files in os.walk(stack_lists_folder):
    for filename in fnmatch.filter(files, "*_TEMP.pix"):
        del_files_list.append(os.path.join(root, filename))

if len (del_files_list) > 0:
    print("\t")

    string_1 = ("----------------- Deleting the temporary files ----------------------")
    time_log.write("%s\n" % string_1)
    print (string_1)
    for ii in del_files_list:
        os.remove(ii)
        print("   " + time.strftime("%H:%M:%S") + " Deleting temporary file: " + ii)


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