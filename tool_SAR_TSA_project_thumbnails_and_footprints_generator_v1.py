#!/usr/bin/env python
'''-----------------------------------------------------------------------------------------------------------------
 * Gabriel Gosselin, CCMEO\ CCRS. 2024-2025                                                                        -
 * -----------------------------------------------------------------------------------------------------------------
 This is a standalone and facultative script in the TSA workflow
'''

# -----------------------------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# -----------------------------------------------------------------------------------------------------------------
#A) Input parameters
parent_folder_search = r"D:\RCMP_data_RCM\QC_03m\other"
keyword = "manifest.safe"
unzip_files_first = "yes"                     # Valid options are "yes" or "no"

# B) Ortho options -  Use a coarse resolution to produce resonably sized thumbnails
Ortho_resolution_X = "25"
Ortho_resolution_Y = "25"
DEM_file = r"D:\RCMP_prj_overviews\aux_files\BC\DEM\BC_Glo30DEM_LatLong_D000.tif"
Elevation_channel = 1

#C) Output options
filename_output_format = 5           # Valid options are 1,2,3,4,5
prefix = "QC03m_"                    # Optional. Leave blank "" for no prefix.
Add_long_lat_suffix = "yes"          # Valid options are "yes" or "no" 
date_format = "unique"               # Valid options are "compact" or "unique"   

output_options = 2   
output_folder = r"E:\RCMP_QC_projet_overview"

# Behaviour when output file exists 
if_file_exists = "regenerate"     # Valid options are "skip" or "regenerate"
# -----------------------------------------------------------------------------------------------------------------
#  Part 2: Notes
# -----------------------------------------------------------------------------------------------------------------
'''
A) Ortho images extents options
    The ortho options are file defined and will correspond to its entire extents. The ortho resolution in X and Y are 
    in meters and are controled by the user specidied parameters Ortho_resolution_X and Ortho_resolution_Y. 
    These parameters are expected to be provided in the projection of the input scene. 
    
B) filename_output_format.  The following metadata are collected

            MATADATA       --> Example(s): 
    <Acquisition_DateTime> --> 2019-05-14T11:08:19.061659Z
        date_time_compact  --> 20190514
        date_time_unique   --> 20190514_110819061659

    <Acquisition_Type>  -->  Very High Resolution 3m, IW_SLC__1SDH

    <BeamMode>        --> 3MCP25, IW

    <PlatformName>         --> RCM, RASARSAT-2, SENTINEL-1
    
    <SourceID>   --> RCM2_OK3490411_PK3500804_1_3MCP25_20250226_224408_CH_CV_SLC_3500804_1
                 --> S1A_IW_SL1__1_DH_20250130T225222_20250130T225254_057680_071BE4_26B3.SAFE

    <OrbitDirection>      --> Ascending, Descending
      (OrbitDir_short)      --> ASC, DESC

    <ProductType>     --> SLC, MLC, GRD

    <Polarizations>     --> RH,RV, HH,VV
     (clean_version)    --> RH-RV, HH-VV

     Matrix_Type   --> s2c, s4c
--------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
filename_output_format == 1  # Basename
   RCM--> manifest.safe

filename_output_format == 2  # Source ID
    filename = prefix + <SourceID>  + ".pix"

filename_output_format == 3  # date_time_compact
    filename = prefix + <date_time_compact> + ".pix"

filename_output_format == 4  # date_time_unique
    filename = prefix + <Acquisition_Type2> + <date_time_unique> + ".pix"
    
filename_output_format == 5
    filename = prefix + <beam mode> + <Orbit_Direction2> + < date_time_unique> +  ".pix"

filename_output_format == 6

--------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
iutput_options
    1 =  Geotiff thumbnail  + vector footprint
    2 =  Geotiff thumbnail  + vector footprint + PIX degraded ortho
    3 =  Vector footprint
 
--------------------------------------------------------------------------------------------------------------
--------------------------------------------------------------------------------------------------------------
ortho_exist_behavior
1 = skip
2 = delete and regenerate the ortho
3 = go foward and add a unique id base on count, 

'''
# -----------------------------------------------------------------------------------------------------------------
#  Part 3: Imports
# -----------------------------------------------------------------------------------------------------------------

# The following locale settings are used to ensure that Python is configured
# the same way as PCI's C/C++ code.  
import locale
locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

import os
import fnmatch
import sys
import time
import math
import shutil

from pci.ortho import ortho
from pci.pyramid import pyramid
from pci.psiqinterp import psiqinterp
from pci.model import model
from pci.pcimod import pcimod
from pci.scale import scale
from pci.thr import thr
from pci.fav import fav
from pci.reproj import reproj
from pci.bit2poly import bit2poly
from pci.poly2pnt import poly2pnt
from pci.exceptions import PCIException
from pci.api import datasource as ds
from pci.api.cts import crs_to_mapunits
from pci import nspio

from TSA_utilities.SAR_TSA_utilities_definitions import unzip_batch

# -----------------------------------------------------------------------------------------------------------------
#  Part 4: Input parameters validation
# -----------------------------------------------------------------------------------------------------------------
print ("---------------------------------------------------------------------------------------------------------")
print (time.strftime("%H:%M:%S") + " Input parameters validation")
start = time.time()

# Hardcoded parameters / Discarded parameters
yes_validation_list = ["yes", "y", "yse", "ys"]
no_validation_list = ["no", "n", "nn"]
yes_no_validation_list = yes_validation_list+no_validation_list

#A) Input parameters
if not os.path.exists(parent_folder_search): 
    print ("Error - The input folder does not exists or the path is wrong")
    sys.exit()

unzip_files_first = unzip_files_first.lower()
if unzip_files_first in yes_validation_list: 
    unzip_files = True
elif unzip_files_first in no_validation_list: 
    unzip_files = False
else: 
    print ('Error - the unzip_files_first parameter must be set with "yes" or "no"')
    sys.exit()

# B) Ortho options 
if not os.path.exists(DEM_file): 
    print ("Error - The input DEM does not exists or the path/filename is wrong")
    sys.exit()

if isinstance(Ortho_resolution_X, str) is False:
    print ("Warning - The Ortho_resolution_X parameter must be in string format") 
    print ("Automatic conversion")
    Ortho_resolution_X = str(Ortho_resolution_X)

if isinstance(Ortho_resolution_X, str) is False:
    print ("Warning - The Ortho_resolution_Y parameter must be in string format") 
    print ("Automatic conversion")
    Ortho_resolution_X = str(Ortho_resolution_X)

#C) Output options
if filename_output_format not in [1,2,3,4,5]: 
    print ("Error - valid filename_output_format options are 1,2,3,4,or 5")
    sys.exit()

Add_long_lat_suffix = Add_long_lat_suffix.lower()
if Add_long_lat_suffix not in yes_no_validation_list: 
    print ('Error - The Add_long_lat_suffix parameter must be "yes" or "no"')
elif Add_long_lat_suffix in yes_validation_list: 
    Add_long_lat_suffix = True
else: 
    Add_long_lat_suffix = False

if date_format.lower() not in ["compact","unique"]:
    print("Error - specified date_format options is invalid")
    print('Accepted values are: "compact" or "unique"')
    sys.exit()

# Create the output file 
if not os.path.exists(output_folder):  
    os.makedirs(output_folder)

if output_options not in [1,2,3]: 
    print("Error - the output_options parameters must be set with 1,2 or 3")
    sys.exit()

if output_options in [1,2]: 
    output_thumbails = True
else: 
    output_thumbails = False

if if_file_exists.lower() not in ["skip","regenerate"]:     
    print('Error - valid options for existing_data are "skip" or "regenerate"')
    sys.exit()
else: 
    info_message_skip = ("   output file already exists - skip (if_file_exists = skip)")
    info_message_regn = ("   output file already exists - regenerate (if_file_exists = regenerate)")

# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------
#  Part 5: Main program
# -----------------------------------------------------------------------------------------------------------------
# -----------------------------------------------------------------------------------------------------------------

# A.1 ) Run the unzip fucntion if requested. Only *.zip files are supported. 
if unzip_files is True: 
    print ("-------------------------------------------------------------------------------------------------------")
    print ("Input scenes unzipping")

    fld_unzip = os.path.join(output_folder, "0_unzipped_scenes")
    
    unzip_batch (parent_folder_search, fld_unzip, if_file_exists, info_message_regn, info_message_skip)
    parent_folder_search = fld_unzip
  

# A.2 ) Find the files to process
files_to_ortho_list = []
for root, dirs, files in os.walk(parent_folder_search):
    for filename in fnmatch.filter(files, keyword):
        files_to_ortho_list.append(os.path.join(root,filename))                
                  
if len(files_to_ortho_list) == 0:
    print ("\t")
    print ("Error - no files found")
    print ("Specify another keyword or search folder")
    sys.exit()

total_files = str(len(files_to_ortho_list))
print("Number of inputs files: "+ total_files)
print("\t")

for input_files in files_to_ortho_list: 
    print (input_files)

# B ) Creatinng a list of output files name according to the selected format
# B.1) Extracting the metadata
print ("\t")
print ("-------------------------------------------------------------------------------------------------------")
print ("Input scenes metadata extraction")
print ("\t")

Acquisition_DateTime2 = []
Acquisition_Type2 = []
SourceID2 = []
PlatformName2 = []
Orbit_Direction2 = []
BeamMode2 = []
ProductType2 = []
projection_list = []
polarization_chan = []
Matrix_type2 = []
chan_count2 = []

nb_files = str(len(files_to_ortho_list))
count = 1
for ii in files_to_ortho_list:
    print ("\t")
    print ("   " + time.strftime("%H:%M:%S") + " Extracting metadata, file  " + str(count) + " of " + nb_files)
    print ("   Reading file -->" + ii)

    with ds.open_dataset(ii,ds.eAM_READ) as ds2:
        aux = ds2.aux_data

        num_channels = ds2.chan_count
        chan_count2.append(num_channels)
                                                                                   
        projection_list.append(crs_to_mapunits(ds2.crs))

        Acquisition_DateTime = aux.get_file_metadata_value('Acquisition_DateTime')
        Acquisition_DateTime2.append(Acquisition_DateTime) 

        Acquisition_Type=aux.get_file_metadata_value('Acquisition_Type')
        Acquisition_Type2.append(Acquisition_Type)

        matrix_type=aux.get_file_metadata_value('Matrix_Type')
        Matrix_type2.append(matrix_type)

        Polarizations = aux.get_file_metadata_value('Polarizations')
        temp_pol_c = Polarizations.replace(",", "-")
        polarization_chan.append (temp_pol_c)

        SourceID=aux.get_file_metadata_value('SourceID')
        SourceID2.append(SourceID)

        PlatformName=aux.get_file_metadata_value('PlatformName')
        if PlatformName == "RADARSAT-2":
            pname = "RS2"
        elif PlatformName == "SENTINEL-1": 
            pname = "S1"
        else:         
            pname = PlatformName
        PlatformName2.append(pname)
       
        Orbit_Direction = aux.get_file_metadata_value('OrbitDirection')
        if Orbit_Direction == "Descending":
            orbit_type = "DESC"
        elif Orbit_Direction == "Ascending":
            orbit_type = "ASC"
        else:
            orbit_type = Orbit_Direction      
        Orbit_Direction2.append(orbit_type)

        BeamMode = aux.get_file_metadata_value('BeamMode')
        BeamMode2.append(BeamMode)

        ProductType = aux.get_file_metadata_value('ProductType')
        ProductType2.append(ProductType)
        
    count = count + 1

# B.2 ) Formating the acquisition date
date_time_compact = []
date_time_unique = []  # necessary for scenes on the same swath.
string_replace = str.maketrans({'-': '', ':': '', 'T': '_', '.': '', 'Z': ''})

for ii in Acquisition_DateTime2:   
    output_string = ii.translate(string_replace)
    date_time_compact.append(output_string[:8])
    date_time_unique.append (output_string)

if date_format.lower() == "compact":
    date_out = date_time_compact
elif date_format.lower() == "unique":
    date_out = date_time_unique
else:
    sys.exit()

# B.3) Creating the final ouput file names according to the user selection
output_file_name=[]
if filename_output_format == 1:
    for input_file in files_to_ortho_list :
        basename = os.path.basename (input_file)
        output_file_name.append(basename)

elif filename_output_format == 2:    # Source ID
    for ii in SourceID2:
        temp = prefix + ii + ".pix"
        output_file_name.append(temp)

elif filename_output_format == 3:   # Compact date
    for date in  date_time_unique:
        temp = prefix + date + ".pix"
        output_file_name.append(temp)

elif filename_output_format == 4:
    for acq_type, date in zip (Acquisition_Type2, date_time_unique):
        temp=prefix + acq_type + "_" + date + ".pix"
        output_file_name.append(temp)

elif filename_output_format == 5:
    for acq_type, date, orbit in zip (BeamMode2, date_time_unique, Orbit_Direction2):
        temp = prefix + acq_type + "_" + orbit + "_" + date + ".pix"
        output_file_name.append(temp)

print("\t")
print ("Selected filename_output_format : " + str (filename_output_format))
count = 1
for ii, jj in zip (files_to_ortho_list, output_file_name): 
    print ("Output file name, file " + str(count) + " of " + nb_files)
    print ("   input-->" + ii)
    print ("   output-->" + jj)
    print ("\t")
    
    count = count + 1


# -----------------------------------------------------------------------------------------------------------------
# C) Orthorectification
# -----------------------------------------------------------------------------------------------------------------
print ("\t")
print ("-------------------------------------------------------------------------------------------------------")
print ("Input scenes Orthorectification")
print ("\t")

ortho_details_dump = os.path.join (output_folder, "ortho_info_dump.txt")
nspio.enableDefaultReport(ortho_details_dump)

Orthorectified_scene_list = []
nb_files = str(len(files_to_ortho_list))
print ("\t")
count = 1
for input_file, outname in zip(files_to_ortho_list, output_file_name):

    print ("   " + time.strftime("%H:%M:%S") + " orthorectifying file " + str(count) + " of " + nb_files)
    print ("   Input file: " + input_file)

    base = os.path.basename(input_file)
    mfile=input_file
    dbic = []          
    mmseg = []
    dbiw = []
    srcbgd = "NONE"

    filo = os.path.join(output_folder, outname)

    ftype = "PIX"
    foptions = "TILED256"
    outbgd = [-32768.00000]
    edgeclip = []
    tipostrn = ""
    filedem = DEM_file
    dbec = [Elevation_channel]
    backelev = []
    elevref =""
    elevunit =""
    elfactor = []
    proc = ""
    sampling = [1]
    resample = "near"
    mapunits=""
    ulx=""
    uly=""
    lrx=""
    lry=""
    bxpxsz = Ortho_resolution_X
    bypxsz = Ortho_resolution_Y

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
            pyramid(file = filo, dboc = [], force = "yes", olevels = [], poption= "aver")
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)
        print ("\t")

    Orthorectified_scene_list.append(filo) 
    count = count + 1
nspio.enableDefaultReport('term')

# -----------------------------------------------------------------------------------------------------------------
# D) Footprints and thumbnails creation
# -----------------------------------------------------------------------------------------------------------------
# D.1) Creating the ortho foorprints excluding the nodata values. Creating a bitmap using 

number_of_files = str(len(Orthorectified_scene_list))   
count = 1 
for input_file, m_type, n_chan in zip(Orthorectified_scene_list, Matrix_type2, chan_count2): 
    print ("\t")
    print (((time.strftime("%H:%M:%S")) + " Footprints creation, file " + str(count)+ " of " + number_of_files))
    print ("   input file-->" + ii)
    
    file = input_file
    dbic = [1]
    dbob = []
    tval = [-1000, 100000]  # threshold range (min,max)
    comp = 'off'
    dbsn = 'footprit'  # output segment name
    dbsd = 'bitmap to generate footprint'  # output segment description
   
    try:  
        thr(file, dbic, dbob, tval, comp, dbsn, dbsd)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)

    # ---------------------------------------------------------------------------------------------------------
    print ("   Creating the vector file")
    fili = input_file
    dbib = [2]
    outname = os.path.basename(input_file[:-4])
    filo = os.path.join(output_folder, outname + "_footprint.pix")
    smoothv = "no"
    dbsd = "ortho_footprint"  # output layer description
    ftype = ""  # defaults to PIX
    foptions = ""  # output format options
    
    try:  
        bit2poly(fili, dbib, filo, smoothv, dbsd, ftype, foptions)
    except PCIException as e:
        print (e)
    except Exception as e:
        print (e)
    
    if Add_long_lat_suffix is True: 
        # Converting the polygon to a point file (*)centroid)
        print ("   Add_long_lat_suffix selected, converting to Long\Lat (if necessary) and calculating the footoprint centroid")
        fili = filo
        dbvs = [2]                 
        filo =  os.path.join (output_folder, "temp_" + outname + "_centroid.pix")
        dbsd     = ""                   # use default
        ftype    = "pix"                # output format type
        foptions = ""                   # output options

        poly2pnt(fili, dbvs, filo, dbsd, ftype, foptions)

        with ds.open_dataset(filo) as ds1:
            centroid_map_projection = crs_to_mapunits(ds1.crs)
            centroid_map_projection = ''.join(centroid_map_projection.split())

            if centroid_map_projection != "LONG/LATD000": 
                print ("   reprojecting the centroid to lat/long") 

                fili    =  filo
                dbic    =   [] # use channel 1
                dbsl    =   [2]
                sltype  =   "VEC"

                centro_base = os.path.basename(filo)
                centro_out = os.path.basename(centro_base[:-4] + "_LatLong.pix")

                filo = os.path.join(output_folder, centro_out)
                ftype   =   ""  # uses PIX format by default
                foptions    =   ""  # no file options are used
                repmeth =   ""  # uses pixel/line bounds by default
                dbsz    =   []  # uses 512 x 512 by default
                pxsz    =   []  # uses 1.0, 1.0 by default
                maxbnds =   "NO"    # uses user-defined bounds
                mapunits    =   "LONG/LAT D000"
                llbounds    =   "YES"
                ulx =   ""
                uly =  ""
                lrx =   ""
                lry =   ""
                resample    =   ""  # uses NEAR by default
                proc    =   ""  # uses AUTO by default
                tipostrn    =   ""  # don't use tile positioning transformation

                reproj( fili, dbic, dbsl, sltype, filo, ftype, \
                foptions, repmeth, dbsz, pxsz, maxbnds,\
                mapunits, llbounds, ulx, uly, lrx, \
                lry, resample, proc, tipostrn )


                with ds.open_dataset(filo) as ds2:
                    io = ds2.get_vector_io(2)
                    xs = []
                    ys = []
                    for index, shape in zip(io.shape_ids, io):
                    # iterate over rings in the shape:
                        xs.append(shape.extents[0]) 
                        ys.append(shape.extents[1])

                    Long_pt = math.trunc(xs[0] * 100) /100.0
                    Lat_pt = math.trunc(ys[0] * 100) /100.0

                    if Long_pt < 0:
                        long_q = "W"
                        Long_pt = abs(Long_pt)
                    else: 
                        long_q = "E"
                    if Lat_pt < 0:
                        lat_q = "S"
                        Lat_pt = abs(Lat_pt)
                    else: 
                        lat_q = "N"
                
                    LongLat_suffix = (str(Long_pt) + long_q + "_" + str(Lat_pt) + lat_q)
                    print ("   output suffix-->" + LongLat_suffix)
    # ---------------------------------------------------------------------------------------------------------
   
    if output_thumbails is True: 
        print ("   thumbnail option is selected")

        if m_type in ["s1c","s2c","s4c"]:  #% Need to convert to real data
            
            print("    matrix type is " + m_type)
            fili = input_file
            dbic = [] 
            cinterp = "Int" 
            dboc = []
            outname = os.path.basename(input_file)
            filo = os.path.join(output_folder, "temp_iq_" + outname )

            ftype = "PIX" 
            foptions = "" 

            try:
                psiqinterp (fili, dbic, cinterp, dboc, filo, ftype, foptions)
            except PCIException as e:
                print (e)
            except Exception as e:
                print (e)   
            input_file = filo   
           
        if n_chan == 1:
            dbic_scale = [1]
        elif n_chan == 2:       # Dual pol and compact pol case
            dbic_scale = [1,2]
        elif n_chan == 4: 
            dbic_scale = [1,2,4] # Quand-pol data case
        else: 
            dbic_scale = [1]   

        # Filtering
        file = input_file
        dbic = dbic_scale
        dboc = dbic_scale
        flsz = [5,5] # use a 21 x 21 filter
        mask = [] # process entire image
        bgrange = [] # no background values
        failvalu = [] # no failure value
        bgzero = '' # default, set background to zero

        fav(file, dbic, dboc, flsz, mask, bgrange, failvalu, bgzero)

        if n_chan == 2:  # special case for dual pol data, we will create a RGB thumbnail by duplicating the first channel
            print ("   Dual pol data handling")
            file = input_file
            pciop = "add"
            pcival = [0,0,0,1]    # possible bug here, assume that the input file channels is always 32R
            pcimod (file, pciop, pcival)
        
            file = input_file
            source = "%3 = %1"
            undefval =  [0]

            model(file,source, undefval)
            dbic_scale = [1,2,3]

        # Thumbnail_creation
        print ("   Scaling to 8U and exporting the thumbnail to tif format")
        fili =  input_file
        if Add_long_lat_suffix is False: 
            filo = os.path.join(output_folder, outname[:-4] + "_thumbnail.tif")
        else: 
            filo = os.path.join(output_folder, outname[:-4] + "_thumbnail_" + LongLat_suffix + ".tif")
        dbic	=  dbic_scale		# 16-bit radar image
        dboc	=  []		# default
        dbiw	=  []		# default, process entire image
        dbow	=  []           # default, output entire image
        inrange	=  []		# default
        trim	=  [1]		# trim 1% of low and high values
        outrange = [1,254]		# default
        sfunct	=  "root"		# automatic normalized quantization
        datatype = "8U"		# default to 8U
        ftype	=  "tif"	# output to a TIFF
        foptions = ""		# output file options
        
        try:
            scale( fili, filo, dbic, dboc, dbiw, dbow, inrange, trim, outrange, sfunct, datatype, ftype, foptions)
        except PCIException as e:
            print (e)
        except Exception as e:
            print (e) 

        count = count +1

# -----------------------------------------------------------------------------------------------------------------
# E) Deleting the intermediary files
# -----------------------------------------------------------------------------------------------------------------
print("\t")
print (time.strftime("%H:%M:%S") + " Deleting the intermediary files")
print("\t")

if unzip_files is True: 
    shutil.rmtree(fld_unzip)

if output_options == 2: # Keep the degraded PIX orthos as well

    for ii in Orthorectified_scene_list: 
        outname = ii[:-4]
        outname2 = (outname + "_ortho.pix")
        os.rename(ii, outname2)

keywords = ["*thumbnail.tif","*thumbnail*.tif", "*footprint*", "*_ortho.pix"]  # The files we need to keep
try:
    for filename in os.listdir(output_folder):
        file_path = os.path.join(output_folder, filename)

        # Check if filename matches any pattern
        if os.path.isfile(file_path) and not any(fnmatch.fnmatch(filename, pattern) for pattern in keywords):
            os.remove(file_path)
            print(f"Deleted: {filename}")
        else:
            print(f"Kept: {filename}")
except Exception as e:
    print(f"Error: {e}")

 
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
