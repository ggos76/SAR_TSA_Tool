
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS                                                    -
 * ----------------------------------------------------------------------
'''

# The following locale settings are used to ensure that Python is configured
# the same way as PCI's C/C++ code.  
import locale
locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

import os, fnmatch, time, sys
from datetime import datetime

import numpy as np

from pci.ortho import ortho
from pci.pyramid import pyramid
from pci.exceptions import PCIException
from pci.api import datasource as ds
from pci.api.cts import crs_to_mapunits
from pci.replacenans import replacenans
from pci.psiqinterp import psiqinterp
from pci.poly2bit import poly2bit
from pci.model import model
from pci import nspio

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def nan_replace(step_nans, input_folder_nans):

    print ("\t")
    print ("--------------------------------------------------------------")
    print ("Checking for NANs")
    print ("\t")

    calibration_report_txt = os.path.join(input_folder_nans, step_nans +
                                          "replace_nans_info.txt")
    nspio.enableDefaultReport(calibration_report_txt)

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

    return ()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def ortho_run (input_folder_for_ortho, output_folder_ortho, DEM_file, DEM_elevation_channel, ortho_bounds_option, AOI_file,
               AOI_file_segment_number, ortho_resolution_X, ortho_resolution_Y, generate_overviews):

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


    # C) Orthorectification
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
        dbic = []  # Blank, will process all channels of the input file.
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

        print("   Output file: " + filo)
        if os.path.exists(filo):
            print ("File already exist - skip")
        else:
            # pyramids options
            file = filo
            force = 'yes'
            poption = 'aver'
            dboc = []
            olevels = []

            try:
                ortho(mfile, dbic, mmseg, dbiw, srcbgd, filo, ftype, foptions, outbgd,
                      ulx, uly, lrx, lry, edgeclip, tipostrn, mapunits, bxpxsz, bypxsz, filedem,
                      dbec, backelev, elevref, elevunit, elfactor, proc, sampling, resample)
                if generate_overviews is True:
                    pyramid(file, dboc, force, olevels, poption)

            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
            print("   Output orthorectified file: " + filo)
        count = count + 1

    return ()

# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def psiqinterp_run (search_folder, keyword, interp_type, suffix, TSA_layers, TSA_labels,
                    ps_output_folder, unique_files, prefix):

    file_to_process = []
    for root, dirs, files in os.walk(search_folder):
        for filename in fnmatch.filter(files, keyword):
            file_to_process.append(os.path.join(root, filename))

        # A) Find the unique scenes (mainly for the coregistered files)
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

    print ("\t")

    if interp_type == "int":
         count = 1
         for ii in file_to_process:

            print((time.strftime("%H:%M:%S") + " Input Scene: " + ii))
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

                if os.path.exists(filo):
                    print ("File already exist - skip")
                else:
                    try:
                        psiqinterp(fili, dbic, cinterp, dboc, filo, ftype, foptions)
                    except PCIException as e:
                        print(e)
                    except Exception as e:
                        print(e)
                    chan = chan + 1
            print ("\t")
            count = count + 1

    elif interp_type == "amp":
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

            print("   output:" + filo)
            if os.path.exists (filo):
                print ("File already exist - skip")
            else:
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

def create_list (prefix, suffix, search_folder, Fld_Output_stack_lists, TSA_channels, TSA_channel_labels):

    for in_chan in TSA_channels:
        in_label = TSA_channel_labels[in_chan - 1]

        # A) Search for the file to process
        keyword = "*_" + in_label + suffix + ".pix"
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
                if suffix == "_int":
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
                raster4 = np.delete(raster3, np.where(raster3 > 1.0))


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

def file_size_check (files_list):


    # Safety check, File size and projection verifications. All files must have the same number of lines and
    # columns and the same projection

    # -------------------------------------------------------------------------------------------------------
    # using the first file of the list as a reference.
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

# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
def stack_masking (input_stack, mask_type, mask_file, mask_seg_number,no_data_value, output_folder):

    print ("   Mask type selected " + mask_type)
    fili = mask_file
    dbvs = mask_seg_number  # polygon layer
    filo = input_stack
    dbsd = ''
    pixres = []
    ftype = "PIX"  # output format type
    foptions = ""  # output options

    poly2bit(fili, dbvs, filo, dbsd, pixres, ftype, foptions)

    with ds.open_dataset(input_stack) as ds1:
        num_channels = ds1.chan_count
        print ("   There are " + str(num_channels) + " channels in the input stack")


    count = 1
    for input_chan in range (1,num_channels +1):

        out_print = (("   " + (time.strftime("%H:%M:%S")) + " Masking channel " + str(count) + " of " + str(num_channels)))
        sys.stdout.write("\r" + out_print)
        sys.stdout.flush()

        output_model_file = []

        if mask_type == "inclusion":
            string1 = ("if %%2 = 1 then")
            output_model_file.append(string1)
            string1 = ("%" + str(input_chan) +"= %" + str(input_chan))
            output_model_file.append(string1)
            string1 = ("else")
            output_model_file.append(string1)
            string1 = ("%" + str(input_chan) +"=" + str(no_data_value))
            output_model_file.append(string1)
            string1 = ("endif")
            output_model_file.append(string1)

            file_name = ("model_inclusion_mask_for_channel_" + str(input_chan) + ".txt")

        if mask_type == "exclusion":
            string1 = ("if %%2 = 1 then")
            output_model_file.append(string1)
            string1 = ("%" + str(input_chan) +"=" + str(no_data_value))
            output_model_file.append(string1)
            string1 = ("else")
            output_model_file.append(string1)
            string1 = ("%" + str(input_chan) +"= %" + str(input_chan))
            output_model_file.append(string1)
            string1 = ("endif")
            output_model_file.append(string1)

            file_name = ("model_exclusion_mask_for_channel_" + str(input_chan) + ".txt")

        file_model = os.path.join(output_folder, file_name)
        with open(file_model, "w") as f:
            f.write("\n".join(output_model_file))

        file = input_stack
        source = file_model
        undefval = []

        try:
            model (file, source, undefval)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
        count = count + 1
        os.remove(file_model)
    return ()
# ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def stack_min_max (input_stack,no_data_value, min_floor, max_floor, output_folder, reassign_type):

    with ds.open_dataset(input_stack) as ds1:
        num_channels = ds1.chan_count
        print ("   There are " + str(num_channels) + " channels in the input stack")

        '''
        # Could be the optimal solution with numpy, to be rewritten later
        arr = input raster channel de PCI
        np.place (arr, arr == 0, -10)
        np.place(arr, arr <= min_floor and arr != no_data_value, min_floor)
        np.place(arr, arr <= min_floor and arr != no_data_value, no_data_value)
        np.place(arr, min_floor == 0, no_data_value)
        '''
        count = 1
        for input_chan in range(1, num_channels + 1):

            print(("   " + (time.strftime("%H:%M:%S")) + " Applying min/max floor to channel " + str(count) +
                   " of " + str(num_channels)))
            output_model_file = []

            if reassign_type == "to_no_data":
                string1 = ("if %" + str(input_chan) + "<= " + str(min_floor) + " then")
                output_model_file.append(string1)
                string1 = ("%" + str(input_chan) + " = " + str(no_data_value))
                output_model_file.append(string1)
                string1 = "endif"
                output_model_file.append(string1)
                string1 = ("if %" + str(input_chan) + ">= " + str(max_floor) + " then")
                output_model_file.append(string1)
                string1 = ("%" + str(input_chan) + " = " + str(no_data_value))
                output_model_file.append(string1)
                string1 = "endif"
                output_model_file.append(string1)

                file_name = ("model_min_max_floor_to_nodata_channel_" + str(input_chan) + ".txt")

            if reassign_type == "to_min_max":
                string1 = ("if %" + str(input_chan) + "<= " + str(min_floor) + " then")
                output_model_file.append(string1)
                string1 = ("%" + str(input_chan) + " = " + str(min_floor))
                output_model_file.append(string1)
                string1 = "endif"
                output_model_file.append(string1)
                string1 = ("if %" + str(input_chan) + ">= " + str(max_floor) + " then")
                output_model_file.append(string1)
                string1 = ("%" + str(input_chan) + " = " + str(max_floor))
                output_model_file.append(string1)
                string1 = "endif"
                output_model_file.append(string1)

                file_name = ("model_min_max_floor_to_min_max_channel_" + str(input_chan) + ".txt")

            file_model = os.path.join(output_folder, file_name)
            with open(file_model, "w") as f:
                f.write("\n".join(output_model_file))

            file = input_stack
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
def date_formater (input_date):

    gpsd = str(input_date)
    gps_MM = gpsd[4: 6]
    gps_DD = gpsd[6: 8]
    gps_YYYY = gpsd[0: 4]

    output_date_format = (gps_MM + "/" + gps_DD + "/" + gps_YYYY)
    return (output_date_format)

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