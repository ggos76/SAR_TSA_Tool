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
#input_stack_statistics_file = r"C:\Users\ggosseli\Desktop\S1_tests_Ottawa_2024\output\6_Output_stack_lists\VV_coh_stack_statistics.txt"
input_stack_statistics_file = r"D:\QC03M_7224_4504\6_Output_stack_lists\QC03M_7224_4504_RH_int_stack_statistics.txt"
prefix = ""                                # optionnal, leave blank for no prefix
no_data_value = -32768.0000

#B.1) Time series preparation and data stacking
stacking_method = "catalyst"             # Valid option are "numpy or "catalyst"
output_type = 2                          # 1: analysis layers only    2: analysis and stack together   3: analysis and stack separate

# B.2) Filtering and mask options
apply_masking = "yes"       # Valid option are "yes" or "no"
mask_type = "exclusion"    # Valid option are "inclusion" or "exclusion"
mask_file = r"\\W-BSC-A157283\share_GG\CanUS_border_1m_UTM18T_D000_v4_2km_buffer.pix"
mask_seg_number = [2]

apply_mean_filter = "yes"               # Valid option are "yes" or "no"
filter_size_X_Y = [3,3]                 # Tuple of odd integer

# B.3) Apply  Min / Max bounds
apply_min_max_bounds = "no"
min_floor = 0.0
max_floor = 1.5
reassign_type = "to_min_max"            # Options are "to_no_data" or "to_min_max"


#C) Other options
# Generate overviews - either yes or no,
generate_overviews = "yes"


# -----------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
1) There is no verification if the specified DEM covers completely the spatial
   extents of the AOI. If not the script will run to completion but all
   interferograms will be blank after INSRAW.

   
 2) stacking_method
    2.1:  numpy. Doen't rely on Catalyst tools but could be slow and crash with very large datasets.
    2.2:  catalyst : Use PSSARTSA to stack intensity layers and inscohstats to stack coherence layers. 
                     These algotithms are not always available depending on the license type. 
   
    The two stackling methond will give the same results. 
   
'''
# -----------------------------------------------------------------------------------------------------------------
#  Imports
# -----------------------------------------------------------------------------------------------------------------
import sys
import os
import time
import locale
import shutil
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
prefix = ""

#  Version control - do nothing for now.
print("\t")
print(pci.version)

print("Installed python version: " + sys.version)
py1 = str(sys.version_info[0])
py2 = str(sys.version_info[1])
py3 = (py1 + "." + py2)
python_version = float(py3)
if python_version < 2.1:
    print("You are using Python v" + str(python_version))
    print("You need to update to Python 3.6 or newer versions")
    sys.exit()
print("\t")

# A) Check if input file exists and the quantity type (intensity or coherence)
if not os.path.exists(input_stack_statistics_file):
    print ("Error - The input_stack_statistics_file does not exist or the path/filename is wrong.")
    sys.exit()

# determining the input type (int or coh)
check_type = os.path.basename (input_stack_statistics_file)

if "int_stack_statistics" in check_type: 
    data_type_int = True
    data_type_coh = False
    print ("The input layers type to stack is Intensity (int)")
elif "coh_stack_statistics" in check_type:    
    data_type_int = False
    data_type_coh = True
    print ("The input layers type to stack is Coherence (coh)")
else: 
    print ("Error - The input file must be a stack_statistics file  (*_stack_statistics.txt)")
    sys.exit()


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

# B.4) Apply  Min / Max bounds
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


# D) Creating the output folder and the output file name.  
    # Note: In anticipation for full automatation, we will create automatically the output folder and 
    # the output file name.      
two_up = Path(input_stack_statistics_file).resolve().parents[1]
stack_folder = "7_TSA_stacks"

output_folder = os.path.join (two_up, "7_TSA_stacks")
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

base = os.path.basename(input_stack_statistics_file)
base_out = base[:-20]
print (base_out)


# For output_type == 1  or output_type == 3
outfile_tsa_analysis = os.path.join(output_folder, (prefix + base_out + "TSA_stack_analysis"))
check_out = outfile_tsa_analysis
# For output_type == 2
outfile_tsa_analysis_stack = os.path.join (output_folder, (prefix + base_out + "TSA_stack_analysis_data"))
check_out = outfile_tsa_analysis_stack
# For output_type == 3
outfile_tsa_data = os.path.join(output_folder, (prefix + base_out + "TSA_stack_data"))

#Check if the specified output file already exists. If yes error.
if os.path.exists(check_out):
    print("Warning - the output file already exists")
    print("Output_file-->" + check_out )
    print("Delete the existing file or specify a different name")
    sys.exit()

# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------    
#                                           Main Program
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------    
start = time.time()

print("\t")
print("---------------------------------------------------------------------------------------------------------")
print("                             Reading and parsing the input_stack_statistics_file                         ")
print("---------------------------------------------------------------------------------------------------------")
print("\t")

Input_file_lines = []
with open(input_stack_statistics_file, "r") as input_file:
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
print ('----------------------------------------------------------------------------------------------------------')
print("                                   File size (L X C) verification                                          ")
print ('----------------------------------------------------------------------------------------------------------')

print("\t")
files_list = input_file_list
file_size_check (files_list)


print("\t")
print ('----------------------------------------------------------------------------------------------------------')
print("                                      Stack data preprocessing                                             ")
print ('----------------------------------------------------------------------------------------------------------')
print("\t")

# -----------------------------------------------------------------------------------------------------------------
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

# -----------------------------------------------------------------------------------------------------------------
if apply_min_max_bounds is True:
    print ("\t")
    print("Stack data preprocessing - Minimum and maximum bounds option is selected")
    input_stack = input_file_list
    stack_min_max (input_stack, no_data_value, min_floor, max_floor, reassign_type)
if apply_min_max_bounds is False:
    print("No min max requested")

# -----------------------------------------------------------------------------------------------------------------
if apply_masking is True:
    print("Stack data preprocessing - Applying the exclusion or inclusion vector mask")
    input_stack = input_file_list
    stack_masking(input_stack, mask_type, mask_file, mask_seg_number,no_data_value, output_folder)
if apply_masking is False :
    print(time.strftime("%H:%M:%S") + " Stack preparation - Exclusion or inclusion mask not requested")

print("\t")
  

#------------------------------------------------------------------------------------------------------------------
# STACKING METHOD USING NUMPY
#------------------------------------------------------------------------------------------------------------------
if use_numpy is True: 

    print(time.strftime("%H:%M:%S") + " Selected stacking method is NUMPY")
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
    print(time.strftime("%H:%M:%S") + " Selected stacking method is CATALYST")
    

    if data_type_coh is True:  # hack to make coherence data works with PSSARTSA to ensure time ordered channels. 
        # We create the "Acquisition DateTime" metadata and fiil it with the  "Dep_Acquisition_DateTime"
       for ii in input_file_list:  
       
            with ds.open_dataset(ii,ds.eAM_WRITE) as ds12:   # Note .eAM_WRITE
                aux = ds12.aux_data
                Dep_DateTime = aux.get_file_metadata_value('Dep_Acquisition_DateTime')
                
                mdmap = aux.file_metadata
                mdmap['Acquisition_DateTime'] = Dep_DateTime
                aux.file_metadata = mdmap
                ds12.aux_data = aux


    # preparing the mfile
    temp_mfile = os.path.join(output_folder, "temp_mfile_for_data_ingestion.txt")
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

    if output_type == 1: # analysis layers only
        print("   " + time.strftime("%H:%M:%S") + " Output type 1: analysis layers only")
        stack =  "no"							
        input_file =  os.path.join(output_folder, outfile_tsa_analysis + out_filter + ".pix" )
    if output_type == 2:  # analysis and stack together
        print("   " + time.strftime("%H:%M:%S") + " Output type 2: analysis and data stack in the same file")
        stack =  "yes"
        input_file =  os.path.join(output_folder, outfile_tsa_analysis_stack + out_filter + ".pix" )
    if output_type == 3:  # analysis and stack separate
        print("   " + time.strftime("%H:%M:%S") + " Output type 3: analysis and data stack in separate files")
        stack =  "yes"
        input_file =  os.path.join(output_folder, outfile_tsa_data + out_filter + ".pix" )
    print ("   Output file: " + input_file)

    filo = input_file

    try:
        pssartsa(mfile, dbic, mask, maskfile, stack, flsz, filo)
    except PCIException as e:
        print(e)
    except Exception as e:
        print(e)

    if generate_overviews is True: 
        print("   " + time.strftime("%H:%M:%S") + " generating the overviews for the output stack")
        file = filo
        dboc = []
        force = ""	
        olevels=[]
        poption = "AVER"
        pyramid( file, dboc, force, olevels, poption )


    if os.path.exists(temp_mfile):
        os.remove(temp_mfile)
    # -----------------------------------------------------------------------------------------------------------------
    # Split the TSA data analysis and the TSA data stack

    if output_type == 3: 
      
        # First we export the fourteen first layers to create the idependent data analysis file. 
        fili =	input_file
        filo =	os.path.join(output_folder, outfile_tsa_analysis + out_filter + ".pix" )
        dbiw =	[]
        dbic =	[1,-14]
        dbib =	[]
        dbvs =	[]
        dblut =	[]
        dbpct =	[]
        ftype =	"TIF"
        foptions = ""

        try:
            fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
        
        # Second we delete the same fourteen layers from the original pix file to only keep
        # the stack data.
        file = input_file
        pciop = "DEL"
        pcival =[1,-14]
        try:
            pcimod(file, pciop, pcival)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
 
# ----------------------------------------------------------------------------------------------------------------
print("\t")
print ('----------------------------------------------------------------------------------------------------------')
print ('----------------------------------------------------------------------------------------------------------')
print(time.strftime("%H:%M:%S") + "   All process completed")
print("\t")
end = time.time()

ellapse_time_seconds = round((end - start), 2)
ellapse_time_minutes = round((ellapse_time_seconds / 60), 2)
ellapse_time_hours = round((ellapse_time_seconds / 3600), 2)

print("Processing time (seconds): " + str(ellapse_time_seconds))
print("Processing time (minutes): " + str(ellapse_time_minutes))
print("Processing time (hours): " + str(ellapse_time_hours))
