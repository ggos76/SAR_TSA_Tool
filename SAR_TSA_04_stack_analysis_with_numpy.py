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
input_scenes_list = r"E:\SAR_TSA_tests\2018_F186_S1A_sub\6_Output_stack_lists\VH_coh_stack_statistics.txt"
input_type = "coh"                       # Options are "int" or "coh"
output_folder = r"E:\SAR_TSA_tests\2018_F186_S1A_sub\7_TSA"
output_file_name = "2018_F186_S1A_VH_Coherence_0_1.pix"

no_data_value = -32768.0000

#B) Time series preparation
# B.1) Mean filter
apply_mean_filter = "yes"  # Valid option are "yes" or "no"
filter_size_X_Y = [5,5]

# B.2) Apply an exclusion or an inclusion mask
apply_masking = "no"       # Valid option are "yes" or "no"
mask_type = "exclusion"    # Valid option are "inclusion" or "exclusion"
mask_file = r"E:\TSA_test2\LSTP_clip_layer2.pix"
mask_seg_number = [2]

# B.3) Apply  Min / Max bounds
apply_min_max_bounds = "yes"
min_floor = 0
max_floor = 1
reassign_type = "to_no_data"            # Options are "to_no_data" or "to_min_max"


#C) Other options
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
#  Imports
# ----------------------------------------------------------------------------------------------
import sys
import os
import time
import locale
import shutil

import pci
from pci.pcimod import pcimod
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

import numpy as np
locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -----------------------------------------------------------------------------------------------------------------
#  Parameters Validation
# -----------------------------------------------------------------------------------------------------------------

# Hardcoded parameters
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list + no_validation_list
prefix = "x_"

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

# A) Input files
if not os.path.exists(input_scenes_list):
    print ("Error - The input_scenes_list does not exist or the path/filename is wrong.")
    sys.exit()

input_type = input_type.lower()
if input_type not in ["int", "coh"]:
    print("Error - The entered input_type is not valid.")
    print('Valid options are: "int" or "coh"')
    sys.exit()

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

#Check if the specified output file already exists. If yes error.

outfile_exists = os.path.join(output_folder, output_file_name)
if os.path.exists(outfile_exists):
    print("Warning - the output file already exists")
    print(outfile_exists)
    print("Delete the existing file or specify a different name")
    sys.exit()


#B) Time series preparation

# B.1) Mean filter
if apply_mean_filter.lower() in yes_validation_list:
    apply_mean_filter = True

    Fsize_X_modulo = filter_size_X_Y[0] % 2
    Fsize_Y_modulo = filter_size_X_Y[1] % 2
    if Fsize_X_modulo != 1 or Fsize_Y_modulo!= 1 :
        print ("Error - The filter_size_X_Y must only be compose of odd integers ")
        sys.exit()
else:
    apply_mean_filter = False

# B.2) Apply an exclusion or an inclusion mask
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

# B.3) Apply  Min / Max bounds
if apply_min_max_bounds.lower() in yes_validation_list:
    apply_min_max_bounds = True

    # check for other conditions
    reassign_type = reassign_type.lower()
    if reassign_type not in ["to_no_data", "to_min_max"]:
        print("Error - The reassign_type parameter is not valid.")
        print('Valid options are: "to_no_data" or "to_min_max"')
        sys.exit()

    min_floor = float(min_floor)
    max_floor = float(max_floor)
    if min_floor >= max_floor:
        print("Error - The min_floor must be inferior to the man_floor")
        sys.exit()
else:
    apply_min_max_bounds = False

# C) Other options
generate_overviews = generate_overviews.lower()
if generate_overviews in yes_validation_list:
    generate_overviews = True
else:
    generate_overviews = False


# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
#                                                           Main Program
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
start = time.time()

print("\t")
print("---------------------------------------------------------------------------------------------------")
print("                            Reading and parsing the input_scenes_list                              ")
print("---------------------------------------------------------------------------------------------------")
print("\t")

Input_file_lines = []
with open(input_scenes_list, "r") as input_file:
    for line in input_file:
        Input_file_lines.append(line.strip())

input_file_list = []
mid_date_list = []
ref_date_list = []
dep_date_list = []

for ii in Input_file_lines:
    # temp = []
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
print(time.strftime("%H:%M:%S") + " List of input files to process")

print("   File -----> acquisition date")
count = 1
for ii, jj in zip(input_file_list, mid_date_list):
    print ("   file " + str(count) + "|" + nb_files + " " + ii + " -----> " + jj)

    count = count + 1

# -----------------------------------------------------------------------------------------------------------------
# File size and projection verifications. All files must have the same number of lines and columns
# and the same projection
# -----------------------------------------------------------------------------------------------------------------
print("\t")
print("---------------------------------------------------------------------------------------------------")
print("                               File size (L X C) verification                                      ")
print("---------------------------------------------------------------------------------------------------")

print("\t")
files_list = input_file_list
file_size_check (files_list)

# -----------------------------------------------------------------------------------------------------------------
#                                              Creating the stack of data stack
# -----------------------------------------------------------------------------------------------------------------

print("\t")
print ('----------------------------------------------------------------------------------------------------------')
print("                                    Creating the data stack                                                ")
print ('----------------------------------------------------------------------------------------------------------')

# here we assume the list to be in the chronological order since it was created by TSA_02
print("\t")
print(time.strftime("%H:%M:%S") + " Copying and moving the first file of the stack")
first_file = input_file_list[0]
print ("   " + first_file)

base = os.path.basename(output_file_name[:-4])
base2 = base + "_stack.pix"
copy_out = os.path.join(output_folder, base2)
shutil.copy2(first_file, copy_out)

print ("\t")
print(time.strftime("%H:%M:%S") + " Transferring the other files in the data stack")

file_to_transfer = input_file_list
file_to_transfer.pop(0)
nb_files = str(len(file_to_transfer))

print(time.strftime("%H:%M:%S") + " Creating output stack empty channels")
fsize = round((os.path.getsize(copy_out)/ 1073741824), 3)

total_estimated_size = len(input_file_list) * fsize
print ("   Total estimated size of the output stack: " + str (total_estimated_size) + " GB")

print ("\t")
print(time.strftime("%H:%M:%S") + " Adding the empty channels to the stack")
file = copy_out
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
    filo = copy_out
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

if apply_mean_filter is True:
    print ("\t")
    print (time.strftime("%H:%M:%S") + " Stack preparation - Average filtering with a (X-Y) " +
           str (filter_size_X_Y[0]) + "x" + str (filter_size_X_Y[1]) + " filter")

    max_chan = len (mid_date_list)
    count = 1
    for ii in range (1, max_chan + 1):
        print("   " + time.strftime("%H:%M:%S") + " Filtering channel " + str(count) + " of " + str(max_chan))
        file = copy_out
        dbic = [ii]  # use elevation data
        dboc = [ii]  # overwrite input data
        flsz = [filter_size_X_Y[0], filter_size_X_Y[1]]
        mask = []  # process entire image
        bgrange = []  # no background values
        failvalu = []  # no failure value
        bgzero = ''  # default, set background to zero

        try:
            fav(file, dbic, dboc, flsz, mask, bgrange, failvalu, bgzero)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
        count = count + 1


# -------------------------------------------------------------------------------------------------------------
# Apply_masking if needed and assign to do data value.
#  if the mask is set with inclusion, we need to set the pixels outside the mask with No_data
 # if the mask is set with exclusion, we need to set the pixels inside the mask with No_data
print("\t")
print("\t")
if apply_masking is True:
    print(time.strftime("%H:%M:%S") + " Stack preparation - Applying the exclusion or inclusion vector mask")
    input_stack = copy_out
    stack_masking(input_stack, mask_type, mask_file, mask_seg_number,no_data_value, output_folder)
if apply_masking is False :
    print(time.strftime("%H:%M:%S") + " Stack preparation - Exclusion or inclusion mask not requested")

print("\t")
print("\t")
if apply_min_max_bounds is True:
    print(time.strftime("%H:%M:%S") + " Stack preparation - Minimum and maximum bounds not requested")
    input_stack = copy_out
    stack_min_max (input_stack,no_data_value, min_floor, max_floor, output_folder, reassign_type)
if apply_min_max_bounds is False:
    print("No min max requested")


# -------------------------------------------------------------------------------------------------------------
# Renaming the channels according to the mid date and creating a metadata of the same name.
print("\t")
print(time.strftime("%H:%M:%S") + "Stack preparation - Adding channels description and metadata")

# open dataset in write mode
channel = 1
for mid_date, ref_date, dep_date in  zip (mid_date_list, ref_date_list, dep_date_list):

    # print (mid_date + "-->" + str(channel))
    with ds.open_dataset(copy_out, ds.eAM_WRITE) as ds1:
        # get the AuxiliaryData
        aux = ds1.aux_data

        # get the metadata for the first channel
        mdmap = aux.get_chan_metadata(channel)
        # set the NO_DATA_VALUE to -32768
        mdmap['MID_DATE'] = mid_date
        mdmap['REF_DATE'] = ref_date
        mdmap['DEP_DATE'] = dep_date
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

with ds.open_dataset(copy_out) as ds9:

    print(time.strftime("%H:%M:%S") + " Stack to process: " + copy_out)
    reader=ds.BasicReader(ds9)
    # read the raster channels
    raster = reader.read_raster(0, 0, reader.width, reader.height)
    num_channels=ds9.chan_count
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
out_raster = gobs.copy_array_to_raster(output_ts_raster)
base = output_file_name[:-4]

out_file_ts = os.path.join(output_folder, base + "_TSA_file.pix")

with ds.new_dataset(out_file_ts, 'PCIDSK', '') as write_ts:
    writer =  ds.BasicWriter(write_ts)     # create a writer to write the raster
    writer.create(out_raster)                # create the file on disk based on the size and datatype of raster
    writer.write_raster(out_raster)          # write the raster to the new file
    writer.crs = ref_crs                    # write the coordinate system
    writer.geocoding = ref_geocoding        # write the geocoding information / must be of same

print("\t")
print(time.strftime("%H:%M:%S") + " Writing some metadata to the output file")
with ds.open_dataset(out_file_ts, ds.eAM_WRITE) as ds1:
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


print("\t")
print("---------------------------------------------------------------------------------------------")
print("---------------------------------------------------------------------------------------------")
print(time.strftime("%H:%M:%S") + "   All process completed")

print("\t")
end = time.time()

ellapse_time_seconds = round((end - start), 2)
ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

print("Processing time (seconds): " + str(ellapse_time_seconds))
print("Processing time (minutes): " + str(ellapse_time_minutes))
print("Processing time (hours): " + str(ellapse_time_hours))


