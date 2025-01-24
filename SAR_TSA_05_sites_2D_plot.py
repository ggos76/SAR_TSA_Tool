#!/usr/bin/env python
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. July 2024                                      -
 * ----------------------------------------------------------------------
'''

# ----------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# ----------------------------------------------------------------------------------------------
# A) Input/output files
input_stack = r"E:\SAR_TSA_tests\RS22010_Roberge_lake\7_TSA\RB2010_TSA_VV_coherence_0_1.pix"
input_polygons_file = r"D:\Wapusk_2010\Roberge_lake\Polygons_for_stats_v2.pix"
polygons_file_vec_seg = 3
output_folder = r"E:\SAR_TSA_tests\RS22010_Roberge_lake\8_plots"

# B) Plot options
plot_type = "grp"   # options are "grp" (grouped) or "ind" (individual)
plot_x_axis_label = "Date"
plot_y_axis_label = "Intensity (VV sigma)"

use_x_axis_limit = "yes"
x_axis_min_max = [20100605, 20101015]    # Use the following date format: YYYYMMDD
use_y_axis_limit = "yes"
y_axis_min_max = [0, 1]

# C) Other options
# Delete previous output files if they already exist, useful for test runs
delete_file_if_existing = "yes"

# -----------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
TBD

'''

# ---------------------------------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------------------------------
import sys
import os
import time
import locale
import datetime as dt

import matplotlib.pyplot as plt
import matplotlib.dates as mdates

import pci
from pci.fexport import fexport
from pci.model import model
from pci.grdpol import grdpol
from pci.exceptions import PCIException
from pci.api import datasource as ds

from TSA_utilities.SAR_TSA_utilities_definitions import date_formater
from TSA_utilities.SAR_TSA_utilities_definitions import check_YYYYMMDD

import numpy as np
locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

# -----------------------------------------------------------------------------------------------------------------
#  Parameters Validation
# -----------------------------------------------------------------------------------------------------------------
#  Version control - do nothing for now.
print ("----------------------------------------------------------------------------------------------------------")
print(time.strftime("%H:%M:%S") + "   Input parameters validation")
print("\t")
print("PCI version: " + pci.version)

print("Installed python version: " + sys.version)
py1 = str(sys.version_info[0])
py2 = str(sys.version_info[1])
py3 = (py1 + "." + py2)
python_version = float(py3)
if python_version < 3.6:
    print("You are using Python v" + str(python_version))
    print("You need to update to Python 3.6 or newer versions")
    sys.exit()

# Hardcoded parameters
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list + no_validation_list

# A) Input/output files
if not os.path.exists(input_stack):
    print ("Error - The iinput_stack does not exist or the path/filename is wrong.")
    sys.exit()
if not os.path.exists(input_polygons_file):
    print ("Error - The input_polygons_file does not exist or the path/filename is wrong.")
    sys.exit()

# check if  polygons_file_vec_seg is integer and >=2
if np.issubdtype(type(polygons_file_vec_seg), int) is False:
    print('Error - The "polygons_file_vec_seg" number  is not an integer')
    sys.exit()
if polygons_file_vec_seg < 2:
    print('Error - The "polygons_file_vec_seg" number must be >=2')
    sys.exit()

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# B) Plot options
plot_type = plot_type.lower()   # options are "grp" (grouped) or "ind" (individual)
if plot_type not in ["grp", "ind"]:
    print ('Error -thr  parameter <plot_type> value must be "grp" or "ind"')
    sys.exit()

if plot_type == "ind":
    produce_individual_plots = True
elif plot_type  == "grp":
    produce_individual_plots = False
else:
    print ("Error - Undefined error")
    sys.exit()

use_x_axis_limit = use_x_axis_limit.lower()
if use_x_axis_limit in yes_validation_list:
    use_x_axis_limit = True

    # Check if the entered dates are valid
    variable_name = "x_axis_min_max"
    # Consistency check for most obvious error.
    date_YYYYMMDD = x_axis_min_max[0]
    check_YYYYMMDD(date_YYYYMMDD, variable_name)
    date_YYYYMMDD = x_axis_min_max[1]
    check_YYYYMMDD(date_YYYYMMDD, variable_name)

    # Check if the start date is inferior to the stop date
    if x_axis_min_max[0] > x_axis_min_max[1]:
        print ("Error -  the parameter <x_axis_min_max> start date must be anterior to the stop date")
        print ("   Entered start date : " + str (x_axis_min_max[0]))
        print ("   Entered stop date: " + str (x_axis_min_max[1]))
        sys.exit()
else:
    use_x_axis_limit = False

use_y_axis_limit = use_y_axis_limit.lower()
if use_y_axis_limit in yes_validation_list:
    use_y_axis_limit = True
else:
    use_y_axis_limit = False

# C) Other options
delete_file_if_existing = delete_file_if_existing.lower()
if delete_file_if_existing in yes_validation_list:
    delete_file_if_existing = True
else:
    delete_file_if_existing = False

print(time.strftime("%H:%M:%S") + "   Input parameters validation done")

# -----------------------------------------------------------------------------------------------------------------
#                                                   Main Program
# -----------------------------------------------------------------------------------------------------------------
print ("\t")
print ("----------------------------------------------------------------------------------------------------------")
print ("                                           Main program                                                   ")
print ("----------------------------------------------------------------------------------------------------------")
start = time.time()

print("\t")
print(time.strftime("%H:%M:%S") + " Reading the polygons layer")
with ds.open_dataset(input_polygons_file, ds.eAM_READ) as ds3:
    io = ds3.get_vector_io(polygons_file_vec_seg)
    poly_number = len(io.shape_ids)
    print ("   Number of polygons to sample: " + str(poly_number))

print("\t")
print(time.strftime("%H:%M:%S") + " Reading the input stack")

with ds.open_dataset(input_stack, ds.eAM_READ) as ds2:
    width = ds2.width  # Number of columns
    height = ds2.height  # Number of row
    chans_count = ds2.chan_count
    print ("   stack dimensions: " + str(width)+"C x "+ str(height) + "L")
    print ("   number of channels to sample: " + str (chans_count))

# We need to check if the projections are compatible and reproject the polygons file if needed.


# Create a raster channel (32R). This new raster channel will contain the shape_ids
print("\t")
print(time.strftime("%H:%M:%S") + " Creating the polygons file for data sampling")
print ("   Creating an empty channel")
fili = input_stack
base = os.path.basename (fili)
filo = os.path.join (output_folder, base[:-4] + "_polygons.pix")
dbiw =	[]
dbic =	[1]
dbib =	[]
dbvs =	[]
dblut =	[]
dbpct =	[]
ftype =	"PIX"
foptions = ""

if os.path.exists(filo) and delete_file_if_existing is True:
    os.remove(filo)

try:
    fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions)
except PCIException as e:
    print(e)
except Exception as e:
    print(e)

file = filo
source = "%1=0"
undefval = []

try:
    model(file, source, undefval)
except PCIException as e:
    print(e)
except Exception as e:
    print(e)
print ("   Burning the polygon layer into the empty channel")
print ("   The background value is -1")
filv = input_polygons_file
dbvs =	[polygons_file_vec_seg]
file =	filo
dboc =	[1]
fldnme = ""	# use default "ATTRIBUTE"
algo = "" # use default "POLYGON"

try:
    grdpol( filv, dbvs, file, dboc, fldnme, algo )
except PCIException as e:
    print(e)
except Exception as e:
    print(e)

raster_polygon = filo


# Creating the plots
print ("\t")
print(time.strftime("%H:%M:%S") + " Generating the plots")
print ("\t")
with ds.open_dataset(input_stack) as ds1:
    reader = ds.BasicReader(ds1)
    raster = reader.read_raster(0, 0, reader.width, reader.height)
    raster_stack = raster.data [:,:,:]
    aux_stack = ds1.aux_data

with ds.open_dataset(raster_polygon) as ds2:
    reader = ds.BasicReader(ds2)
    raster = reader.read_raster(0, 0, reader.width, reader.height)
    raster_polygon = raster.data
    raster_polygon2 = raster_polygon.reshape(-1)
# Do a double check to make sure the raster dimensions match.

count = 1
for in_poly in range (0, poly_number):

    print (time.strftime("%H:%M:%S") + "Processing polygon id " + str (in_poly) + " (" + str (count) +
           " of " + str(poly_number)+")")

    ras_poly_mean = []
    ras_poly_std = []
    chan_out2 = []
    chan_mid_date_YYYYMMDD = []

    for ii in range (0, chans_count):

        chan_aux = aux_stack.get_chan_metadata(ii+1)
        middate = chan_aux['MID_DATE']
        input_chan = raster_stack [:, :, ii]
        input_chan2 = input_chan.reshape(-1)

        data_subset = np.delete(input_chan2, np.where(raster_polygon2 != float(in_poly)))
        array_lenght = str(len(data_subset))
        if ii == 0 :
            print ("   Number of pixels in current polygon: " + array_lenght)

        data_subset_mean = np.mean(data_subset)
        data_subset_std = np.std(data_subset)
        stats_string = str (round (data_subset_mean, 4)) +" ± " + str (round (data_subset_std, 4))
        ras_poly_mean.append(data_subset_mean)
        ras_poly_std.append(data_subset_std)
        print ("   " + middate + "--> Channel " + str (ii + 1)  + " (mean ± std) : "  + stats_string )
        chan_out = ii + 1
        chan_out2.append (chan_out)
        chan_mid_date_YYYYMMDD.append(middate)

    # ////////////////////////////////////////////////////////////////////////////////////////
    # Reformating the dates from YYYYMMDD   to MM/DD/YYYY for Matplotlib
    chan_mid_date_date_MMDDYYYY = []
    for input_date in chan_mid_date_YYYYMMDD:

        (output_date_format) = date_formater (input_date)
        chan_mid_date_date_MMDDYYYY.append(output_date_format)

    plot_date = np.array (chan_mid_date_date_MMDDYYYY)
    # print (plot_date)
    date_x = [dt.datetime.strptime(d,'%m/%d/%Y').date() for d in plot_date]
    formatter = mdates.DateFormatter("%Y-%m-%d")

    ax = plt.gca()
    ax.xaxis.set_major_formatter(formatter)
    locator = mdates.MonthLocator()
    ax.xaxis.set_major_locator(locator)

    # Format the date to use for the X axis
    plot_mean = np.array(ras_poly_mean)
    plot_sd = np.array(data_subset_std)
    plot_chan_x = np.array(chan_out2)

    if use_x_axis_limit is True:
        (output_date_format) = date_formater(str(x_axis_min_max[0]))
        x_min = dt.datetime.strptime(output_date_format,'%m/%d/%Y').date()
        (output_date_format) = date_formater(str(x_axis_min_max[1]))
        x_max = dt.datetime.strptime(output_date_format,'%m/%d/%Y').date()
        plt.xlim (x_min, x_max)
    if use_y_axis_limit is True:
        plt.ylim(y_axis_min_max[0], y_axis_min_max[1])

    plt.grid(color='lightgray', linestyle='-', linewidth=0.25)

    base = os.path.basename (input_stack[:-4])
    if produce_individual_plots is True:
        plt.errorbar(date_x, plot_mean, plot_sd, label='_nolegend_', marker='o', linestyle="solid", linewidth=1.0, color=(0.9, 0.6, 0.6))
        fname = os.path.join(output_folder, base + "_polygon_" + str(count)+"_plot.png")
        if os.path.exists(fname) and delete_file_if_existing is True:
            os.remove(fname)
        plot_title = (base + " (polygon " + str(count - 1) + ": " + str (array_lenght)+ " pixels)")

    if produce_individual_plots is False:
        label_out = ("poly. " + str((count - 1)))

        plt.plot(date_x, plot_mean, label=label_out, marker='o', linestyle="solid", linewidth=1.0)
        fname = os.path.join(output_folder, base + "_all_polygons_grouped_plot.png")
        if os.path.exists(fname) and delete_file_if_existing is True:
            os.remove(fname)
        plot_title = base
        #plt.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
        #plt.tight_layout()
        plt.legend()


    plt.xlabel(plot_x_axis_label, fontsize=12)
    plt.xticks(rotation=45)
    plt.ylabel(plot_y_axis_label, fontsize=12)
    plt.title(plot_title, fontsize=12)

    plt.show(block=False)
    plt.pause(1)

    # Saving the plot to a PNG file
    plt.legend(bbox_to_anchor=(1.0, 1.0), loc='upper left')
    plt.tight_layout()
    figure = plt.gcf()
    figure.set_size_inches(10, 8)
    plt.savefig(fname, dpi=300, facecolor='w', edgecolor='w',
                orientation='landscape', format=None, transparent=False,
                bbox_inches='tight', pad_inches=0.05, metadata=None)

    if produce_individual_plots is True:
        plt.close()

    print ("\t")
    count = count + 1


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


