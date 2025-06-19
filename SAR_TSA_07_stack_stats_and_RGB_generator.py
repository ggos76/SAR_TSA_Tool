#!/usr/bin/env python
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. May 2025                                      -
 * -------------------
 ---------------------------------------------------
'''

# ----------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# ----------------------------------------------------------------------------------------------
# A) Input Files -  should point toward 7_   or 8_ 
TSA_stacks_folder = r"E:\QC03M_7158_4495_CP05_A\8_TSA_stacks_postproc"
no_data_value = -32768.0000
prefix = ""                              # leave blank if no prefix is needed

# B) TSA stacks channels to process (muts be between 1 to 14 inclusively
#    See the notes section for the predefined channels
TSA_channels_to_process = [1,2,3,10]

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
The time series tool returns 14 channels named as follows:
01 - Absolute Sum of Changes
02 - Average Intensity
03 - Standard Deviation Intensity
04 - Coefficient of Variation (Intensity)
05 - Median
06 - Maximum Intensity
07 - Minimum Intensity
08 - Span: |Maximum - Minimum|
09 - Maximum (Absolute)  Change
10 - Maximum Change as a percentage of ASC
11 - Index of Maximum (Absolute)  Change
12 - Maximum 2nd (Absolute)  Change
13 - Maximum 2nd Change as a percentage of ASC
14 - Index of Maximum 2nd (Absolute)  Change


TSA_channels_labels contains a predefine set of labels. 

'''

# -----------------------------------------------------------------------------------------------------------------
#  Part 3: imports
# -----------------------------------------------------------------------------------------------------------------
import sys
import os
import time
import locale
import fnmatch
from pathlib import Path
import numpy as np
import csv
from openpyxl import Workbook

import pci
from pci.pcimod import pcimod
from pci.grey2rgb import grey2rgb   
from pci.fexport import fexport
from pci.thr import thr
from pci.bit2poly import bit2poly
from pci.pyramid import pyramid
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

TSA_channels_labels = ["ASC", "AVG", "STDEV", "CVAR", "MEDIAN", "MAX", "MIN", "SPAN", "1MAXCHG", "1MAXCHGPCT", 
                       "1MAXCHGPCTindex", "2MAXCHG", "2MAXCHGPCT", "2MAXCHGPCTindex"]

# Version control
vs_catalyst = pci.version
vs_python = sys.version_info[:3]
version_control (vs_catalyst, vs_python)

if if_file_exists.lower() not in ["skip","regenerate"]:     
    print('Error - valid options for existing_data are "skip" or "regenerate"')
    sys.exit()
else: 
    info_message_skip = ("   output file already exists - skip (if_file_exists = skip)")
    info_message_regn =  ("   output file already exists - regenerate (if_file_exists = regenerate)")

generate_stats_file = True
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

filtered_files = [f for f in input_stack_files if not os.path.basename(f).startswith("RGB_")]
input_stack_files = filtered_files


print("   Number of stack to process: " + str(len(input_stack_files)))
for ii in input_stack_files:
    print ("   " + ii)


outdirname = os.path.dirname(input_stack_files[0])
output_folder = os.path.join (outdirname, "stats_and_RGB_files")
if not os.path.exists(output_folder):
    os.makedirs(output_folder)



# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
# B) Create one statistic files per input stack for the selected channels to process. 
if generate_stats_file is True:

    nb_files = str(len(input_stack_files))
    count = 1
    for in_stack in input_stack_files:

        print("\t")
        print(time.strftime("%H:%M:%S") + " Creating the statistics files, stack " + str(count) + " of " + nb_files)
        print("   Reading the input stack --> " + in_stack) 

        base = os.path.basename (in_stack[:-4])
        output_stats_txt_file = os.path.join(output_folder, base + "_statsfile.txt")
        output_stats_xlsx_file = os.path.join(output_folder, base + "_statsfile.xlsx")

        if os.path.exists (output_stats_xlsx_file) and if_file_exists == "skip": 
            print (info_message_skip)
        else :
            if os.path.exists (output_stats_xlsx_file) and if_file_exists == "regenerate": 
                print (info_message_regn)
                os.remove (output_stats_xlsx_file)
            if os.path.exists (output_stats_xlsx_file) and if_file_exists == "regenerate": 
                os.remove (output_stats_txt_file)

            with ds.open_dataset(in_stack) as ds9:

                header_file = "channel;label;min;max;mean;stdev;count;--;min p1;max p99;count p1 p99"
                out_stats_write = []
                out_stats_write.append(header_file)
            
                reader = ds.BasicReader(ds9)
                # read the raster channels
                raster = reader.read_raster(0, 0, reader.width, reader.height)
                num_channels = ds9.chan_count
            
                nb_chans = str(len(TSA_channels_to_process))
                count_chan = 1
                for in_chan in TSA_channels_to_process: 
                # -----------------------------------------------------------------------------------
                    label_out = TSA_channels_labels[in_chan - 1]
                    stout1 = ("   " + time.strftime("%H:%M:%S") + " Computing the statistics, channel " + str(count_chan) + " of " + nb_chans)
                    stout2 = ("     input channel --> " + str (in_chan) + " - " + label_out)
                    print (stout1 + "   " + stout2)
                    in_chan_data = raster.data[:,:,in_chan - 1]

                    #in_chan_data = raster[in_chan - 1]      # in_chan is 1-based, raster is 0-based
                    in_chan_data = in_chan_data.reshape(-1)
                    # Compute some first order statistics
                
                    in_chan_data2 = np.delete(in_chan_data, np.where(in_chan_data == float(no_data_value)))
                    # raster3 = np.delete(raster2, np.where(raster2 == float(-32768.00000)))

                    arr_mean = str(np.mean(in_chan_data2))
                    arr_std = str(np.std(in_chan_data2))
                    arr_min = str(np.amin(in_chan_data2))
                    arr_max = str(np.amax(in_chan_data2))
                    arr_len = str(len(in_chan_data2))
                
                    p1 = np.percentile(in_chan_data2, 1)
                    p99 = np.percentile(in_chan_data2, 99)
                    filtered_data = in_chan_data2[(in_chan_data2 > p1) & (in_chan_data2 < p99)]
                    p1min = str(np.amin(filtered_data))
                    p99max = str(np.amax(filtered_data))
                    arr_len_p = str(len(filtered_data))

                    output_string1 = (str(in_chan) + ";" + label_out + ";" + arr_min + ";" + arr_max + ";" + arr_mean + ";" + arr_std + ";" + arr_len + ";")
                    output_string2 = (" " + ";" + p1min + ";" + p99max + ";" + arr_len_p)
                    output_string3 = output_string1 + output_string2
                    out_stats_write.append(output_string3)
                
                    count_chan =  count_chan + 1
                 # -----------------------------------------------------------------------------------
                 # Write the output stats file   
                with open(output_stats_txt_file, "w") as f:
                    for line in out_stats_write:
                        f.write(line + "\n")

                csv_filename = output_stats_txt_file
                xlsx_filename = output_stats_xlsx_file
                workbook = Workbook()
                sheet = workbook.active
                # Open the CSV file and read its content
                with open(csv_filename, 'r', newline='', encoding='utf-8') as csvfile:
                    reader = csv.reader(csvfile, delimiter=';')
                    for row in reader:
                        parsed_row = []
                        for value in row:
                            try:
                                parsed_row.append(float(value))
                            except ValueError:
                                parsed_row.append(value)
                        sheet.append(parsed_row)

                # Save the Excel file
                workbook.save(xlsx_filename)
                print(f"   Conversion complete: '{xlsx_filename}' created.")    
                os.remove (output_stats_txt_file)

        count = count + 1   # count for input_stack


# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
# C) CREATING THE RGB FILES ASSOCIATED TO THE INPUT STACKS

print("\t")
print("-----------------------------------------------------------------------------------------------------------")
print("                              Generating the RGB images for the selected parameters                        ")
print("-----------------------------------------------------------------------------------------------------------")
print("\t")

out_RGB_list = []
nb_files = str(len(input_stack_files))
nb_chans = str(len(TSA_channels_to_process))

count = 1
for in_stack in input_stack_files:
     print ("\t")   
     print (time.strftime("%H:%M:%S") + " Generating the RGB images, stack " + str (count) + " of " + nb_files )
     print ("   input stack --> " + in_stack)

     count_chan = 1
     for in_chan in TSA_channels_to_process:
        label_out = TSA_channels_labels[in_chan - 1]

        stout1 = ("   " + time.strftime("%H:%M:%S") + " Processing, channel " + str(count_chan) + " of " + nb_chans)
        stout2 = ("     input channel --> " + str (in_chan) + " - " + label_out)
        print (stout1 + "   " + stout2)

        
        filo = os.path.join(output_folder, "RGB_" + prefix + os.path.basename(in_stack)[:-4] + "_" + label_out + ".pix")
        out_RGB_list.append(filo)
        out_pix = filo

        if os.path.exists (out_pix) and if_file_exists == "skip": 
            print (info_message_skip)
        else :
            if os.path.exists (out_pix) and if_file_exists == "regenerate": 
                print (info_message_regn)
                os.remove (out_pix)

            mfile = in_stack
            dbic = [in_chan]
            colormap = "jet"
            inmin = []
            inmax = []
            clampmet = "clamp"
            numbin = [100]

            try:
                grey2rgb(mfile, dbic, colormap, inmin, inmax, clampmet, numbin, filo)
            except PCIException as e:
                print (e)
            except Exception as e:
                print (e)
           
        count_chan = count_chan + 1 

     count = count +1


print ("////////////////////////////////////////////////////////////////////////////////////")

nb_files = str(len(out_RGB_list))
print (nb_files)

for input_RGB in out_RGB_list: 

    print (" fili -->" + input_RGB )
    fili = input_RGB
    filo = os.path.join(output_folder, os.path.basename(input_RGB)[:-4] + ".jpg")
    print (" filo -->" + filo)
    dbiw	=	[]
    dbic	=	[1,2,3]
    dbib	=	[]
    dbvs	=	[]
    dblut	=	[]
    dbpct	=	[]
    ftype	=	"jpg"
    foptions	= ""

    fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions )




'''
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

'''
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
 
'''
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
'''