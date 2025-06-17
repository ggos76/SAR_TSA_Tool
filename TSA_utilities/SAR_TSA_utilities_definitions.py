
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS                                                    -
 * ----------------------------------------------------------------------
'''

# The following locale settings are used to ensure that Python is configured
# the same way as PCI's C/C++ code.  
import locale
locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

import os
import fnmatch
import time
import sys
import shutil
import zipfile
from datetime import datetime

import numpy as np

from pci.ortho import ortho
from pci.pyramid import pyramid
from pci.ras2poly import ras2poly
from pci.reproj import reproj
from pci.exceptions import PCIException
from pci.api import datasource as ds
from pci.api.cts import crs_to_mapunits
from pci.replacenans import replacenans
from pci.psiqinterp import psiqinterp
from pci.poly2bit import poly2bit
from pci.fexport import fexport
from pci.model import model
from pci import nspio


# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def nan_replace(step_nans, input_folder_nans):

    print ("\t")
    print ("--------------------------------------------------------------")
    print ("Checking for NANs")
    print ("\t")

    nan_report_txt = os.path.join(input_folder_nans, step_nans + "replace_nans_info.txt")
    nspio.enableDefaultReport(nan_report_txt)

    nans_files_list = []
    for root, dirs, files in os.walk(input_folder_nans):
        for filename in fnmatch.filter(files, "*.pix"):
            nans_files_list.append(os.path.join(root, filename))

    file_number = str(len(nans_files_list))
    count = 1
    for ii in nans_files_list:
        st1 = ((time.strftime("%H:%M:%S")))
        print ("\t")
        print (st1 + "...processing file " + str(count) + " of " +
               file_number)

        file = ii
        dbic = [1]
        newval = "-32768.00000"

        print ("--------------------------------------------------------------")
        nspio.Report.addInfo("fili-->" + file + "\n")

        try:
            replacenans(file, dbic, newval)
        except PCIException as e:
            print (e)
        except Exception as e:
            print (e)
        count = count + 1

    nspio.enableDefaultReport('term')
    return ()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#------------------------------------------------------------------------------------------------------------------
# Orthorectification
#------------------------------------------------------------------------------------------------------------------
def ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel, ortho_bounds_option, 
               AOI_file, AOI_file_segment_number, ortho_resolution_X, ortho_resolution_Y, generate_overviews, 
               TSA_math_xtra_channels, if_file_exists, info_message_skip, info_message_regn, delete_intermediary_files):


    ortho_details_dump = os.path.join (output_folder_ortho, "ortho_info_dump.txt")
    nspio.enableDefaultReport(ortho_details_dump)

    # A) Find the files to orthorectify
    files_to_ortho_list = []
    file_keyword = "*.pix"
    for root, dirs, files in os.walk(input_folder_for_ortho):
        for filename in fnmatch.filter(files, file_keyword):
            files_to_ortho_list.append(os.path.join(root, filename))

    # B) Find the bounds for the orthorectified files
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
            ## DO A VALIDATION OF THE PROJECTION HERE,

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


    # C) Find the channel(s) to orthorectify
    # Special case for orthorectifying the math layers. 
    # Check if the input folder contains "3_3_1_Math_layers"
    target = "3_3_1_Math_layers"

    if target in os.path.normpath(input_folder_for_ortho): 
        print ("Orthorectification of the match layers")
        with ds.open_dataset(files_to_ortho_list[0]) as ds4:
             chans_list = []
             chans = ds4.chan_count
             for ii in range (1,chans+1,1): 
                 chans_list.append(ii)

        print (chans_list)
        dbic_ortho = chans_list[-TSA_math_xtra_channels:] 
    else: 
        dbic_ortho = [] # Blank, will process all channels of the input file.
    
    print (dbic_ortho)
    # D) Orthorectification
    number_of_files = str(len(files_to_ortho_list))
    print("\t")
    count = 1
    for input_scene in files_to_ortho_list:
        print ("\t")
        print(" --------------------------------------------------------------------------------------------------")
        print(" --------------------------------------------------------------------------------------------------")
        print(((time.strftime("%H:%M:%S")) + " Orthorectifying file " + str(count) + " of " + number_of_files))
        print("   Input file: " + input_scene)

        mfile = input_scene
        dbic = dbic_ortho  
        mmseg = []
        dbiw = []
        srcbgd = "NONE"

        base = os.path.basename(input_scene)
        filo = os.path.join(output_folder_ortho, "o_" + base)

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

        print("   output file--> " + filo)
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
                    pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
            print ("\t")
        count = count + 1

    nspio.enableDefaultReport('term')
    return ()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#----------------------------------------------------------------------------------------------------------------------------
# Conversion from complex to intensity data
#----------------------------------------------------------------------------------------------------------------------------
def psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_layers, TSA_labels, ps_output_folder, unique_files, 
                    prefix, math_chans, if_file_exists, info_message_skip, info_message_regn):
    
    #-------------------------------------------------------------------------------------------------------------------------
    # A) Find the unique scenes (mainly for the coregistered files when SBAS mode is use)
    file_to_process = []
    for root, dirs, files in os.walk(search_folder):
        for filename in fnmatch.filter(files, keyword):
            file_to_process.append(os.path.join(root, filename))

        if unique_files.lower() == "yes":
            Acquisition_DateTime3 = []

            for input_scene in file_to_process:
                with ds.open_dataset(input_scene) as ds1:
                    aux = ds1.aux_data
                    Acquisition_DateTime = aux.get_file_metadata_value('Acquisition_DateTime')
                    Acquisition_DateTime2 = Acquisition_DateTime[:10]
                    Acquisition_DateTime3.append(Acquisition_DateTime2.replace("-", ""))

            Unique_coreg_list = []
            Redundant_date_list = []
            for file, date in zip(file_to_process, Acquisition_DateTime3):
                if date not in Redundant_date_list:
                    Redundant_date_list.append(date)
                    Unique_coreg_list.append(file)

            file_to_process = Unique_coreg_list

    print("\t")
    nb_files = str(len(file_to_process))
    for ii in file_to_process:
        print (ii)
    
    nb_chans = str(len(TSA_layers))

    #-------------------------------------------------------------------------------------------------------------------------
    # B) Converting the complex channels to intensity
    print ("\t")
    if interp_type == "int":
        count = 1
        for ii in file_to_process:
            print ("\t")
            print (time.strftime("%H:%M:%S") + " Processing file " + str (count) + " of " + nb_files)
            print ("   input file-->" + ii)
            # We need to  remove the "coreg: part for individual layers
            base1 = os.path.basename(ii[:-4])
            index = base1.find(prefix)
            base = base1[index:]
            chan = 1
            for in_chan in TSA_layers:
                nb_chan = str(len(TSA_layers))
                fili = ii
                dbic = [in_chan]
                cinterp = interp_type
                dboc = []
                in_label = TSA_labels[in_chan-1]
                filo = os.path.join(ps_output_folder, base + "_" + in_label +"_"+ suffix + ".pix")
                print("   Converting channel " + str(chan) + " of " + nb_chan + ": output file-->" + filo)
                ftype = "PIX"
                foptions = ""

                if os.path.exists (filo) and if_file_exists == "skip": 
                    print (info_message_skip)
                else: 
                    if os.path.exists (filo) and if_file_exists == "regenerate": 
                        print (info_message_regn)
                        os.remove (filo)
                    try:
                        psiqinterp(fili, dbic, cinterp, dboc, filo, ftype, foptions)
                    except PCIException as e:
                        print(e)
                    except Exception as e:
                        print(e)
                    chan = chan + 1
            
            if math_chans == "yes": 
            # The math layers are selected, we will  run PSIQINTERP  again for all channels 
            # Note: FEXPORT doesn't work since the math model cannot be exported.  
                print ("   Exporting the math layers")

                out_math = os.path.dirname (ps_output_folder)
                out_math2 = os.path.join(out_math, "3_3_1_Math_layers")

                fili = ii
                dbic = []
                cinterp = interp_type
                dboc = []
                filo = os.path.join(out_math2, base + "_" + suffix + ".pix")
                print ("   output file-->" + filo)
                ftype = "PIX"
                foptions = ""
             
                if os.path.exists (filo) and if_file_exists == "skip": 
                    print (info_message_skip)
                else: 
                    if os.path.exists (filo) and if_file_exists == "regenerate": 
                        print (info_message_regn)
                        os.remove (filo)
                    try:
                        psiqinterp(fili, dbic, cinterp, dboc, filo, ftype, foptions)
                    except PCIException as e:
                        print(e)
                    except Exception as e:
                        print(e)                    
            count = count + 1

    
    #-------------------------------------------------------------------------------------------------------------------------
    # C) Converting the complex channels to amplitude (coherence layers)
    elif interp_type == "amp":   # producing the coherence layer
        count = 1
        for ii in file_to_process:

            print("\t")
            print(time.strftime("%H:%M:%S") + " Converting file " + str(count) + " of " + nb_files)
            fili = ii
            dbic = []
            cinterp = interp_type
            dboc = []
            base = os.path.basename(ii[:-4])
            filo = os.path.join(ps_output_folder, base + suffix + ".pix")
            ftype = "PIX"
            foptions = ""

            print ("   output file-->" + filo)
            if os.path.exists (filo) and if_file_exists == "skip": 
                print (info_message_skip)
            else: 
                if os.path.exists (filo) and if_file_exists == "regenerate": 
                    print (info_message_regn)
                    os.remove (filo)

                try:
                    psiqinterp(fili, dbic, cinterp, dboc, filo, ftype, foptions)
                except PCIException as e:
                    print(e)
                except Exception as e:
                    print(e)
            count = count + 1
    else:
        print ('Error - The  cinterp parameter is invalid')
        sys.exit()

    return()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def create_list (prefix, suffix, search_folder, Fld_Output_stack_lists, TSA_channels, TSA_channel_labels):

    for in_chan in TSA_channels:
        print (in_chan)
        in_label = TSA_channel_labels[in_chan - 1]

        # A) Search for the file to process
        keyword = "*_" + in_label + suffix + ".pix"

        print (keyword)
        file_to_process = []
        for root, dirs, files in os.walk(search_folder):
            for filename in fnmatch.filter(files, keyword):
                file_to_process.append(os.path.join(root, filename))

        out_filename = (prefix + in_label + suffix + "_stack_list.txt")
        file = os.path.join(Fld_Output_stack_lists, out_filename)
        with open(file, "w") as f:
            f.write("\n".join(file_to_process))


        # B // Collect usefull statistics for all the files founds

        ref_date_list = []
        dep_date_list = []
        date_mid_point_list = []
        raster_lines_list = []
        raster_columns_list = []
        arr_mean_list = []
        arr_std_list = []
        arr_min_list = []
        arr_max_list = []
        arr_med_list = []
        arr_len_list = []
        date_format = "%Y%m%d"

        for input_scene in file_to_process:
            with ds.open_dataset(input_scene, ds.eAM_READ) as ds5:
                aux = ds5.aux_data

                # Extract the compact date and compute the midpoint in needed
                if suffix == "_int" or suffix =="":
                    Acquisition_DateTime = aux.get_file_metadata_value('Acquisition_DateTime')
                    date_ref1 = Acquisition_DateTime[:10]
                    date_ref2 = date_ref1.replace("-", "")
                    date_mid_point_list.append(date_ref2)
                    ref_date_list.append(date_ref2)
                    dep_date_list.append(date_ref2)


                if suffix == "_coh":

                    Ref_Acquisition_DateTime = aux.get_file_metadata_value('Ref_Acquisition_DateTime')
                    date_ref1 = Ref_Acquisition_DateTime[:10]
                    date_ref2 = date_ref1.replace("-", "")

                    Dep_Acquisition_DateTime = aux.get_file_metadata_value('Dep_Acquisition_DateTime')
                    date_dep1 = Dep_Acquisition_DateTime[:10]
                    date_dep2 = date_dep1 .replace("-", "")

                    a = datetime.strptime(str(date_ref2), date_format)
                    b = datetime.strptime(str(date_dep2), date_format)
                    mid = a + (b - a)/2
                    mid2 = str(mid)
                    mid3 = mid2[:10]
                    mid4 = mid3.replace("-", "")

                    ref_date_list.append(date_ref2)
                    dep_date_list.append(date_dep2)
                    date_mid_point_list.append(mid4)

                # Read raster data
                reader = ds.BasicReader(ds5)
                raster = reader.read_raster(0, 0, reader.width, reader.height)
                raster2 = raster.data[:,:,0]

                n_lines = str(reader.height)
                raster_lines_list.append(n_lines)
                n_columns = str(reader.width)
                raster_columns_list.append (n_columns)

                # Compute some first order statistics
                raster2 = raster2.reshape(-1)
                raster3 = np.delete(raster2, np.where(raster2 == float(-32768.00000)))
                if suffix == "_coh":
                    raster4 = np.delete(raster3, np.where(raster3 > 1))
                elif suffix == "_int":
                    raster4 = np.delete(raster3, np.where(raster3 > 1.5))
                else :
                    raster4 = np.delete(raster3, np.where(raster3 > 90))

                arr_mean = str(np.mean(raster4))
                arr_std = str(np.std(raster4))
                arr_min = str(np.amin(raster4))
                arr_max = str(np.amax(raster4))
                arr_med = str(np.median(raster4))
                arr_len = str(len(raster4))

                arr_mean_list.append(arr_mean)
                arr_std_list.append(arr_std)
                arr_min_list.append(arr_min)
                arr_max_list.append(arr_max)
                arr_med_list.append(arr_med)
                arr_len_list.append(arr_len)

        # C) Sorting the list base on the mid-date (coherence) and the ref-date (intensity)
        # Sorting with numpy and reconverting to a list is more secure (?) with different data type .
        date_mid_point_arr = np.array(date_mid_point_list)
        indices = date_mid_point_arr.argsort()                     # Here we short
        date_mid_point_arr_s = date_mid_point_arr[indices]
        date_mid_point_arr_s_out = date_mid_point_arr_s.tolist()

        file_to_process_arr = np.array(file_to_process)
        file_to_process_arr_s = file_to_process_arr[indices]
        file_to_process_s_out = file_to_process_arr_s.tolist()

        ref_date_list_arr = np.array(ref_date_list)
        ref_date_list_arr_s = ref_date_list_arr[indices]
        ref_date_list_s_out = ref_date_list_arr_s.tolist()

        dep_date_list_arr = np.array(dep_date_list)
        dep_date_list_arr_s = dep_date_list_arr[indices]
        dep_date_list_s_out = dep_date_list_arr_s.tolist()

        # ////////////// files dimensions
        raster_lines_list_arr = np.array(raster_lines_list)
        raster_lines_list_arr_s = raster_lines_list_arr[indices]
        raster_lines_list_s_out = raster_lines_list_arr_s.tolist()

        raster_columns_list_arr = np.array(raster_columns_list)
        raster_columns_list_arr_s = raster_columns_list_arr[indices]
        raster_columns_list_s_out = raster_columns_list_arr_s.tolist()

        # ///////////  First order statistics from numpy
        arr_mean_list_arr = np.array(arr_mean_list)
        arr_mean_list_arr_s = arr_mean_list_arr[indices]
        arr_mean_list_s_out = arr_mean_list_arr_s.tolist()

        arr_std_list_arr = np.array(arr_std_list)
        arr_std_list_arr_s = arr_std_list_arr[indices]
        arr_std_list_s_out = arr_std_list_arr_s.tolist()

        arr_min_list_arr = np.array(arr_min_list)
        arr_min_list_arr_s = arr_min_list_arr[indices]
        arr_min_list_s_out = arr_min_list_arr_s.tolist()

        arr_max_list_arr = np.array(arr_max_list)
        arr_max_list_arr_s = arr_max_list_arr[indices]
        arr_max_list_s_out = arr_max_list_arr_s.tolist()

        arr_med_list_arr = np.array(arr_med_list)
        arr_med_list_arr_s = arr_med_list_arr[indices]
        arr_med_list_s_out = arr_med_list_arr_s.tolist()

        arr_len_list_arr = np.array(arr_len_list)
        arr_len_list_arr_s = arr_len_list_arr[indices]
        arr_len_list_s_out = arr_len_list_arr_s.tolist()


        # ------------------------------------------------------------------------------------
        out_text_file = []

        string_header1 = ("file" + ";" + "mid_date"+ ";"+ "ref_date"+ ";"+ "dep_date"+ ";"+ "lines"+ ";"+ "columns"+ ";")
        string_header_2 = ("mean" + ";" + "std" + ";" + "min" + ";" + "max" + ";" + "med" + ";" + "length")
        string_header = string_header1 + string_header_2
        out_text_file.append (string_header)

        for (file_in, mid_in, ref_in, dep_in, lines, columns, mean_in,
             std_in, min_in, max_in, med_in, len_in) in zip (file_to_process_s_out, date_mid_point_arr_s_out,
                                                              ref_date_list_s_out,dep_date_list_s_out,raster_lines_list_s_out,
                                                              raster_columns_list_s_out, arr_mean_list_s_out, arr_std_list_s_out,
                                                              arr_min_list_s_out, arr_max_list_s_out, arr_med_list_s_out,arr_len_list_s_out):

            string1 = (file_in + ";" + mid_in + ";" + ref_in + ";" + dep_in + ";" + lines + ";" + columns + ";")
            string2 = (mean_in + ";" + std_in + ";" + min_in + ";" + max_in + ";" + med_in + ";"+ len_in)
            string_out = string1 + string2
            out_text_file.append(string_out)

        out_filename = (prefix + in_label + suffix + "_stack_statistics.txt")
        file = os.path.join(Fld_Output_stack_lists, out_filename)
        with open(file, "w") as f:
            f.write("\n".join(out_text_file))

    return()


# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def file_size_check (files_list):

    # Safety check, File size and projection verifications. All files must have the same number of lines and
    # columns. 
    
    # Using the first file of the list as a reference.
    with ds.open_dataset(files_list[0], ds.eAM_READ) as ds1:
        reference_width = ds1.width         # Number of columns
        reference_height = ds1.height       # Number of row

        count = 1
        error_count = 0
        error_file = []
        for input_coreg in files_list:

            with ds.open_dataset(input_coreg, ds.eAM_READ) as ds2:

                print(((time.strftime("%H:%M:%S")) + " Reading file dimensions, file " + str(count) +
                       " of " + str(len(files_list))))
                print("   " + input_coreg)

                width = ds2.width         # Number of columns
                height = ds2.height       # Number of row
                print("   " + str(width) +"C x "+str(height) +" L")

                if width != reference_width or height != reference_height:
                    print("**** Error **** Number of lines or columns is different than the reference file")
                    print(str(reference_width) + " C x " + str(reference_height) + " L")
                    error_count = error_count + 1
                    error_file.append(input_coreg)
            print("\t")
            count = count + 1

    if error_count != 0:
        print("\t")
        print("The following file(s) are of different size, reprocess the files or remove them from the input folder")
        for ii in error_file:
            print(ii)
        sys.exit()

    else:
        print("All input files have the same number of lines and columns")
    return()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def stack_masking (input_stack, mask_type, mask_file, mask_seg_number,no_data_value, output_folder):

    # First we create the EASI model that will be use for all the input_files. 
    output_model_file = []
    chan = '%1'
    no_data = str(no_data_value)

    if mask_type == "inclusion":
        string1 = ("if %%2 = 1 then")
        output_model_file.append(string1)
        string1 = (chan + " = " +  chan)
        output_model_file.append(string1)
        string1 = ("else")
        output_model_file.append(string1)
        string1 = (chan + " = " + no_data)
        output_model_file.append(string1)
        string1 = ("endif")
        output_model_file.append(string1)

        file_name = ("model_inclusion_mask_for_channel.txt")

    if mask_type == "exclusion":
        string1 = ("if %%2 = 1 then")
        output_model_file.append(string1)
        string1 = (chan + " = " + no_data)
        output_model_file.append(string1)
        string1 = ("else")
        output_model_file.append(string1)
        string1 = (chan + " = " +  chan)
        output_model_file.append(string1)
        string1 = ("endif")
        output_model_file.append(string1)

        file_name = ("model_exclusion_mask_for_channel.txt")

    output_folder_model = os.path.dirname(input_stack[0])
    file_model = os.path.join(output_folder_model, file_name)
    with open(file_model, "w") as f:
        f.write("\n".join(output_model_file))

    nb_files = str(len(input_stack))
    count = 1
    for ii in input_stack:
        print ("   " + time.strftime("%H:%M:%S") + " Masking file " + str(count) + " of " + nb_files) 
        
        fili = mask_file
        dbvs = mask_seg_number  # polygon layer
        filo = ii
        dbsd = ''
        pixres = []
        ftype = "PIX"  # output format type
        foptions = ""  # output options

        poly2bit(fili, dbvs, filo, dbsd, pixres, ftype, foptions)

        file = ii
        source = file_model
        undefval = []
        try:
            model(file, source, undefval)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)

        count = count + 1
    #os.remove(file_model)
    return ()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def stack_min_max (input_stack, no_data_value, min_floor, max_floor, reassign_type):

    # First we create the EASI model that will be use for all the input_files. 
    output_model_file = []
    chan = '%1'
    no_data = str(no_data_value)
    min_floor = str(min_floor)
    max_floor = str(max_floor)

    if reassign_type == "to_no_data":
        string1 = ("if " + chan + "<= " + min_floor + " then")
        output_model_file.append(string1)
        string1 = (chan + " = " + no_data)
        output_model_file.append(string1)
        string1 = "endif"
        output_model_file.append(string1)
        string1 = ("if " + chan + ">= " + str(max_floor) + " then")
        output_model_file.append(string1)
        string1 = (chan +  " = " + no_data)
        output_model_file.append(string1)
        string1 = "endif"
        output_model_file.append(string1)

        file_name = ("model_min_max_floor_to_nodata_channel.txt")

    if reassign_type == "to_min_max":
        string1 = ("if " + chan + "<>" + no_data + " and " + chan + " <=" + min_floor + " then")
        output_model_file.append(string1)
        string1 = (chan + " = " + min_floor)
        output_model_file.append(string1)
        string1 = "endif"
        output_model_file.append(string1)
        string1 = ("if " + chan + "<>" + no_data + " and " + chan + " >=" + max_floor + " then")
        output_model_file.append(string1)
        string1 = (chan + " = " + str(max_floor))
        output_model_file.append(string1)
        string1 = "endif"
        output_model_file.append(string1)

        file_name = ("model_min_max_floor_to_min_max_channel.txt")

    output_folder_model = os.path.dirname(input_stack[0])

    file_model = os.path.join(output_folder_model, file_name)
    with open(file_model, "w") as f:
        f.write("\n".join(output_model_file))
  
    nb_file = str(len(input_stack))
    count = 1
    for ii in input_stack: 
        print ("   " + time.strftime("%H:%M:%S") + " Applying min-max floor, file " + str (count)  + " of " + nb_file)
        file = ii
        source = file_model
        undefval = []
        
        try:
            model(file, source, undefval)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
        count = count + 1
    
    os.remove(file_model)
      
    return ()
    
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def stack_chans_masking (input_stack, mask_type, mask_file, mask_seg_number,no_data_value):

    output_folder = os.path.dirname(input_stack)

    #A) Convert the input vector mask to a bitmap
    fili = mask_file
    dbvs = mask_seg_number              
    filo = input_stack   
    dbsd     = "canusa_border_mask"                
    pixres   = []           
    ftype    = "pix"             
    foptions = ""                

    poly2bit( fili, dbvs, filo, dbsd, pixres, ftype, foptions )


    # Depending of the previous pre-procesing that were applied, we do not know the bitmap segment number
    # corresponding to the CANUSA mask. We retrieve here the segment number of the most recent bitmap.  
    with ds.open_dataset(input_stack, ds.eAM_READ) as ds3:
        bitmap_list = ds3.bitmap_ids

        last_bitmap = bitmap_list[-1]  # Get the last bitmap segment number
        chans_list = []
        chans = ds3.chan_count
        for ii in range (1,chans+1,1): 
            chans_list.append(ii)


        count = 1
        for in_chan in chans_list:
            aa = ("   masking channel " + str(count) + " of " + str(len(chans_list)))
            sys.stdout.write("\r" + aa)
            sys.stdout.flush()

            model_chans_mask = os.path.join(output_folder, "stack_chans_temp_model.txt")
            model_lines = []
            
            string_1 = "if %%" + str(last_bitmap) + " = 1 then"  # Bitmap is on
            model_lines.append(string_1)

            if mask_type == "inclusion": 
                string_1 = "%" + str(in_chan) + " = " + "%" + str(in_chan)
                model_lines.append(string_1)
                string_1 = "else"
                model_lines.append(string_1)
                string_1 = "%" + str(in_chan) + " = " + str(no_data_value)
                model_lines.append(string_1)
            
            if mask_type == "exclusion":    
                string_1 = "%" + str(in_chan) + " = " +  str(no_data_value)
                model_lines.append(string_1)
                string_1 = "else"
                model_lines.append(string_1)
                string_1 = "%" + str(in_chan) + " = " + "%" + str(in_chan)
                model_lines.append(string_1)

            string_1 = "endif"
            model_lines.append(string_1)

            with open(model_chans_mask, "w") as f:
                f.write("\n".join(model_lines))

            file = input_stack
            source = model_chans_mask
            undefval = []
            model(file, source, undefval)
        
            count = count + 1  

    return()        


# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def date_formater (input_date):

    gpsd = str(input_date)
    gps_MM = gpsd[4: 6]
    gps_DD = gpsd[6: 8]
    gps_YYYY = gpsd[0: 4]

    output_date_format = (gps_MM + "/" + gps_DD + "/" + gps_YYYY)
    return (output_date_format)

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def check_YYYYMMDD (date_YYYYMMDD, variable_name):

    # Check if the date is an integer
    if np.issubdtype(type(date_YYYYMMDD), int) is False:
        print("Error - The parameter <" + variable_name + "> number is not an integer")
        print ("Entered value: " + str(date_YYYYMMDD))
        sys.exit()
    # Check if the YYYYMMDD has 8 digits
    if len(str(date_YYYYMMDD)) != 8:
        print("Error - The parameter <" + variable_name + "> must be 8 digits")
        print("Entered value:" + str(date_YYYYMMDD))
        sys.exit()

    # check if the user entered date is valid
    check_year = int((str(date_YYYYMMDD)[0:4]))
    check_month = int((str(date_YYYYMMDD)[4:6]))
    check_day = int((str(date_YYYYMMDD)[6:8]))

    # Check the year validity
    if  check_year < 1972 or check_year > 2030:
        print("Error - The parameter <" + variable_name + "> year value seems not plausible")
        print(" Accepted values are >1971 and < 2030")
        print("Entered value: " + str(check_year))
        sys.exit()

    # check the month validity
    if check_month < 1 or check_month > 12:
        print("Error - The parameter <" + variable_name + "> month is not valid")
        print("Valid values are 01 to 12")
        print("Entered value: " + str(check_month))
        sys.exit()
    else:
        if check_day < 1 or check_day > 31:
            print("Error - The parameter <" + variable_name + "> day is not valid")
            print("Entered value: " + str(check_day))
            sys.exit()
        else:
            if check_month in [4, 6, 9, 11] and check_day > 30:
                print("Error - The parameter <" + variable_name +
                      "> day is not valid for the corresponding month (" + str (check_month) + "), must be <= 30")
                print("Entered value: " + str(check_day))
                sys.exit()
            elif check_month in [1, 3, 5, 7, 8, 10, 12] and check_day > 31:
                print("Error - The parameter <" + variable_name +
                      "> day is not valid for the corresponding month (" + str (check_month) + "), must be <=31")
                print("Entered value: " + str(check_day))
                sys.exit()
            elif check_month in [2] and check_day > 29:  # No check for leap years
                print("Error - The parameter <" + variable_name +
                      "> day is not valid for the corresponding month (" + str (check_month) + "), must be <=29")
                print("Entered value: " + str(check_day))
                sys.exit()
            else:
                print ("")  # day and month are ok.
    return ()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def get_folder_proctime_and_size (folder_path, proc_stop_time, proc_start_time):
    total_size_bytes = sum(
        os.path.getsize(os.path.join(dirpath, file))
        for dirpath, _, files in os.walk(folder_path)
        for file in files
        if os.path.isfile(os.path.join(dirpath, file))
    )
    size_mb = str(round(total_size_bytes / (1024 * 1024),4))
    time_sec = str (round((proc_stop_time - proc_start_time), 2))
    
   
    input_files_list = []
    for keyword in ["*.pix", "*.txt", "*.log"]:
        for root, dirs, files in os.walk(folder_path):
            for filename in fnmatch.filter(files, keyword):
                input_files_list.append(os.path.join(root, filename))
   
    total_files = str(len(input_files_list))

    out_folder_time_size= (";" + time_sec + ";" +size_mb + ";" + total_files)
    return out_folder_time_size, size_mb

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def unzip_batch (parent_folder_search, fld_unzip, if_file_exists, info_message_regn, info_message_skip):
      
    if not os.path.exists(fld_unzip):
        os.makedirs(fld_unzip)
  
    files_to_unzip = []
    for root, dirs, files in os.walk(parent_folder_search):
        for filename in fnmatch.filter(files, "*.zip"):
            files_to_unzip.append(os.path.join(root, filename))
    
    count = 1
    nb_files = str(len(files_to_unzip))
    unzip_folder = []
    for ii in files_to_unzip:
        print("   " + (time.strftime("%H:%M:%S")) + " unzipping file " + str (count) + " of " + nb_files)
        
        base1 = os.path.basename (ii[:-4])
        out_folder = os.path.join (fld_unzip, base1)

        if os.path.exists (out_folder) and if_file_exists == "skip":
            print (info_message_skip)
        else:  
            if os.path.exists (out_folder) and if_file_exists == "regenerate":
                print (info_message_regn)
                shutil.rmtree(out_folder)
            
            os.makedirs (out_folder)
            unzip_folder.append (out_folder) 
        
            with zipfile.ZipFile(ii, 'r') as zip_ref:
                zip_ref.extractall(out_folder)

        count = count + 1
    
    return ()
    
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////    

def ortho_raster_footprint (prefix, input_raster, output_folder): 

    # this function will output two vector files
    # 1) footprint in the input_rater native projection
    # 2) footprint in LatLong WGS84

    with ds.open_dataset(input_raster, ds.eAM_READ) as ds1:
        raster_MapProjection = crs_to_mapunits(ds1.crs)

    raster_projection = raster_MapProjection.strip()
    raster_projection = raster_projection.replace(" ","")
   

    #1) We export the first channel of the input stack to create a template for the footprints
    print (time.strftime("%H:%M:%S") + "  Exporting the template to create the stack footprint")
    fili = input_raster
    filo = os.path.join(output_folder, prefix + "footprint_template.pix")
    dbiw = []
    dbic = [1]
    dbib = []
    dbvs =[]
    dblut =	[]
    dbpct =	[]
    ftype =	"pix"
    foptions = ""
    try:
        fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions )
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)
    
    stack_template = filo
    
    #2) Setting the first channel to an arbitrary value (100) for the Raster to vector convertion
    file = filo
    source = "%1 = 100"
    undefval = []
    
    try:
        model(file, source, undefval)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)

    #3) Creating the stack vector footprint
    print (time.strftime("%H:%M:%S") + "  Creating the stacks vector footprint")
    fili = filo
    dbic = [1]
    filo = os.path.join(output_folder, prefix + "stack_footprint_" + raster_projection + ".pix")
    smoothv = "no"
    dbsd = "footprint"
    ftype = "pix"
    foptions = ""

    try:
        ras2poly (fili, dbic, filo, smoothv, dbsd, ftype, foptions)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)
    
    stack_footprint_utm = filo

    #4) Reprojecting the foorptint to LatLong WGS84
    print (time.strftime("%H:%M:%S") + "  Reprojecting the vector footprint to LatLong WSG84")
    fili = stack_footprint_utm
    dbic = [] 
    dbsl = [2]
    sltype = ""
    filo = os.path.join(output_folder, prefix + "stack_footprint_LatLong_WGS84.pix")   
    ftype   =   ""  # uses PIX format by default
    foptions    =   ""  # no file options are used
    repmeth =   "BR"  # uses bounds and resolution method
    dbsz    =   []  # not used for BR method
    pxsz    =   []  # uses 25 meters resolution
    maxbnds =   "YES"    # uses maximum bounds
    mapunits = "LONG/LAT D000"
    llbounds    =   "NO"
    ulx =   ""
    uly =   ""
    lrx =   ""
    lry =   ""
    resample = ""  # uses CUBIC resample
    proc =   ""  # uses AUTO by default
    tipostrn = "" # uses CORNER tile positioning transformation with 10 meter stride

    try:
        reproj( fili, dbic, dbsl, sltype, filo, ftype, foptions, repmeth, dbsz, pxsz, maxbnds,\
        mapunits, llbounds, ulx, uly, lrx, lry, resample, proc, tipostrn )
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)
     
    stack_footprint_LatLong = filo    
        
    return stack_template, stack_footprint_utm, stack_footprint_LatLong

    