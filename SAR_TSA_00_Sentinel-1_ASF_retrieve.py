#!/usr/bin/env python
'''----------------------------------------------------------------------------------------------
 * Gabriel Gosselin, CCRS. July 2024                                                            -
 * ----------------------------------------------------------------------------------------------
'''

# ----------------------------------------------------------------------------------------------
#  Part 1: User defined variables
# ----------------------------------------------------------------------------------------------
# A) search_method options are "aoi"or "point"  (should in be in LatLong wgs84)/
search_method = "aoi"
AOI_vector_file = r"D:\Wapusk_2010\Roberge_lake\AOI_Roberge_Lake_UTM15_D000.pix"
AOI_segment_number = 2

point_LongLat_dec = "-93.122, 57.419"

start_date = '2016-01-20'
stop_date = '2022-01-20'

orbit_direction = "asc"    # options are "desc" or "asc"
product_type = "SLC"       # options are "SLC" or "GRD"

output_folder = r"E:\ASF_search1"


# -----------------------------------------------------------------------------------------------------
#  Scripts -  Notes.
# -----------------------------------------------------------------------------------------------------
'''
# https://github.com/asfadmin/Discovery-asf_search

# to install ASF
# pip install asf_search

'''
# ---------------------------------------------------------------------------------------------
#  Imports
# ----------------------------------------------------------------------------------------------
import os
import sys
import asf_search as asf
import numpy as np
import dateutil.parser

import zipfile
import fnmatch
import time

'''

'''
# -----------------------------------------------------------------------------------------------------------------
#  Parameters Validation
# -----------------------------------------------------------------------------------------------------------------

print ("---------------------------------------------------------------------------------------------------------")
print (((time.strftime("%H:%M:%S")) + "   Input parameters validation"))
print ("\t")

search_method = search_method.lower()
if search_method not in ["aoi", "point"]:
    print ("Error - The parameter <search_method> value is not valid")
    print ('   Valid options are "aoi" or "point"')
    sys.exit()


ASC_list = ["a", "asc", "ascending"]
DESC_list = ["d", "desc", "descending"]

orbit_direction = orbit_direction.lower()
if orbit_direction in ASC_list:
    orbit_direction = "ASCENDING"
elif orbit_direction in DESC_list:
    orbit_direction = "DESCENDING"
else:
    print ("Error - The parameter <"+ orbit_direction + "> value is not valid")
    sys.exit()

if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# -----------------------------------------------------------------------------------------------------------------
#                                                   Main Program
# -----------------------------------------------------------------------------------------------------------------
print ("\t")
print ("---------------------------------------------------------------------------------------------------------")
print ("                                ASF catalogue search                                                  ")
print ("---------------------------------------------------------------------------------------------------------")
print ("\t")
# A) # Formatting the input coordinates for the ASF search.
if search_method == "aoi":
    print(((time.strftime("%H:%M:%S")) + "   Selected search method: AOI"))
    print(((time.strftime("%H:%M:%S")) + "   Checking the input AOI conformity"))
    # import some PCI libraries
    from pci.exceptions import PCIException

    from pci.api import datasource as ds
    from pci.api.cts import crs_to_mapunits
    from pci.reproj import reproj
    from pci.fexport import fexport

    with ds.open_dataset(AOI_vector_file) as ds1:
        AOI_MapProjection = crs_to_mapunits(ds1.crs)
        AOI_MapProjection = AOI_MapProjection.replace(" ", "")
        print("\t")
        print("   AOI Map projection: " + AOI_MapProjection)

        # Get the file extension
        file_ext = os.path.basename (AOI_vector_file[-3:])
        file_ext = file_ext.lower()

        if AOI_MapProjection != "LONG/LATD000":
            # we need to reproject to LatLong WGS84
            print ("   The AOI file must be in LatLong WGS84 coordinates  - reprojecting...")

            fili = AOI_vector_file  # input image file
            dbic = []  # use channel 1,2,3
            dbsl = [AOI_segment_number]
            sltype = 'VEC'

            base = os.path.basename (AOI_vector_file)
            filo = os.path.join (output_folder, "LatLongWGS84_" + base)
            AOI_vector_file = filo

            ftype = ""
            foptions = ""  # no file options are used
            repmeth = "BR"  # uses bounds and resolution method
            dbsz = []  # not used for BR method
            pxsz = []  # uses 25 meters resolution
            maxbnds = "YES"  # uses maximum bounds
            mapunits = "LONG/LAT D000"
            llbounds = "NO"
            ulx = ""
            uly = ""
            lrx = ""
            lry = ""
            resample = "CUBIC"  # uses CUBIC resample
            proc = ""  # uses AUTO by default
            tipostrn = "CORNER,10"  # uses CORNER tile positioning transformation with 10 meter stride

            try:
                reproj(fili, dbic, dbsl, sltype, filo, ftype, \
                       foptions, repmeth, dbsz, pxsz, maxbnds, \
                       mapunits, llbounds, ulx, uly, lrx, \
                       lry, resample, proc, tipostrn)
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
        else:
            print("   The AOI file is already in LatLong WGS84 coordinates...")

        file_ext = os.path.basename (AOI_vector_file[-3:])
        file_ext = file_ext.lower()


        if file_ext != "shp":
            # We need to translate to the shapefile format
            print("   The AOI file must be in the shapefile (.shp) format - translating...")
            fili = AOI_vector_file

            base = os.path.basename (AOI_vector_file[:-4])
            filo = os.path.join (output_folder, base + ".shp")
            dbiw = []
            dbic = []
            dbib = []
            dbvs = [2]
            dblut = []
            dbpct = []
            ftype = "shp"
            foptions = ""

            try:
                fexport(fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions)
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
            AOI_vector_file = filo
        else:
            print("   The AOI file is already in the shapefile (.shp) format...")

    print ("\t")
    print(((time.strftime("%H:%M:%S")) + "   Computing the AOI bounding box coordinates"))
    # For the case of an input polygon that may have many vertices, we will simply define a tight bounding box.
    X_coordinates=[]
    Y_coordinates=[]
    # open the dataset in write mode
    with ds.open_dataset(AOI_vector_file, ds.eAM_WRITE) as dataset:
        
        ## DO A VALIDATION OF THE PROJECTION HERE, 
        
        # get vector segment 
        io = dataset.get_vector_io(1)
        xs = []
        ys = []
        # iterate over shapes in the segment
        for index, shape in zip(io.shape_ids, io):
            # iterate over rings in the shape:
            xs.append(round(shape.extents[0],3))
            xs.append(round(shape.extents[2],3))
            ys.append(round(shape.extents[1],3))
            ys.append(round(shape.extents[3],3)) 

    # Upper left and lower right coordinates of the bounding box to 
    # calculate the DBIW window for every image              
    min_X = str(min(xs))
    max_Y = str(max(ys))
    max_X = str(max(xs))
    min_Y = str(min(ys))  

    UL_XY = (min_X + " " + max_Y )
    UR_XY = (max_X + " " + max_Y )
    LL_XY = (min_X + " " + min_Y)
    LR_XY = (max_X + " " + min_Y)
    # Note: ASF  (mtk) needs anti-clockwise  coordinates stating from the upper_right
    coord_string = 'POLYGON((' + (UR_XY + "," + UL_XY + "," + LL_XY + "," + LR_XY + "," + UR_XY) + '))'
    print ("   AOI bounding box coordinates:" + coord_string)
    search_coordinates = coord_string

elif search_method == "point":
    print(((time.strftime("%H:%M:%S")) + "   Selected search method: point"))
    long_lat = point_LongLat_dec.split(",")
    long = str(long_lat[0])
    lat =  str(long_lat[1])
    search_coordinates = ("POINT("+long+" "+lat+")")
    print ("   Point coordinates: " + search_coordinates)
else: 
    print ("Error - Undefined error")
    sys.exit()


##########################################################
# SEARCHING THE ASF catalogue
search_results = asf.geo_search(
    platform = 'SENTINEL-1',
    processingLevel=product_type,
    beamMode = 'IW',
    flightDirection=orbit_direction,
    intersectsWith=search_coordinates,
    start = start_date,
    end = stop_date)

results_dic = search_results.geojson()

asf_scene_list = []
for ii in results_dic["features"]:
    asf_scene_list.append (ii["properties"])
asf_scene_list = list(sorted(asf_scene_list, reverse=True, key=lambda x:  dateutil.parser.isoparse((x['startTime']))))
print (list(x ['startTime'] for x in asf_scene_list ))


url_list2 = (list(x ["url"] for x in asf_scene_list ))
path_list2 = (list(x ["pathNumber"] for x in asf_scene_list ))
frame_list2 = (list(x ["frameNumber"] for x in asf_scene_list ))

count = 1
for out_url, out_path, out_frame in zip(url_list2, path_list2, frame_list2):
    print (str(count) + "---" + str(out_path) + "---" + str(out_frame) +"---" + out_url)
    count = count + 1

# ---------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------
# Find the number of uniques combination of path and frame
# Merge  the path and frame list into a single list
path_frame_list = []

for path_in, frame_in in zip (path_list2, frame_list2):
    out_path_frame = str(path_in) + "_" + str(frame_in)
    path_frame_list.append(out_path_frame)

values, counts = np.unique(path_frame_list, return_counts=True)
print ("\t")
print ("   Number of uniques path frame ")
print ("   Path_Frame ---> scenes count")
for  ii, jj in zip(values, counts):
    print ("   " + ii + "--->" + str(jj))

print ("\t")
print ("---------------------------------------------------------------------------------------------------------")
print ("                 Data downloading for the ASF catalogue ( asf.download_url)                              ")
print ("---------------------------------------------------------------------------------------------------------")
print ("\t")

# --------------------------------------------------------------------------------------------------
print ("\t")
print ("   Enter 1 for downloading all scenes  ")
print ("   Enter 2 for downloading scenes for a particular path_frame")
user_option =  input("   Enter 1 or 2: ")
max_try = 6
num_try = 1
while num_try < max_try:
    if user_option in ["1","2"]:
        break
    else:
        print ("   Error  - you must enter 1 or 2.... try " + str(num_try) + " of " + str(max_try-1))
        user_option = input("   Enter 1 or 2: ")

    num_try = num_try + 1

if user_option == "1":
    print ("   All scenes will be downloaded")

print("\t")
user_path_frame = input("   Enter a Path_Frame number (e.g: 91_24) : ")
num_try = 1
if user_option == "2":
    while num_try < max_try:
        if user_path_frame in values:
            sub_url_list = []
            comp_path_frame = user_path_frame.split("_")
            for out_url, out_path, out_frame in zip(url_list2, path_list2, frame_list2):
                if out_path == int(comp_path_frame[0]) and out_frame == int(comp_path_frame[1]):
                    sub_url_list.append(out_url)

            url_list2 = sub_url_list
            break
        else:
            print ("   Error  - you must enter a valid Path_Frame.... try " + str(num_try) + " of " + str(max_try-1))
            user_path_frame = input("   Enter a Path_Frame number (e.g: 91_24) : ")



print ("\t")
print ('\t')
print(((time.strftime("%H:%M:%S")) + "   Downloading the selected scenes (total = " + str (len(url_list2))+" )"))

fld_download_zip = os.path.join(output_folder, "zip_downloaded")
if not os.path.exists(fld_download_zip):
    os.makedirs(fld_download_zip)

session = asf.ASFSession().auth_with_creds('gosga76B', 'Mabe1994')
#asf.download_urls(urls=url_list2, path=fld_download_zip, session=session, processes=1)

count = 1
for ii in url_list2:
    print("\t")
    print(("   " + (time.strftime("%H:%M:%S")) + "   Downloading file " + str(count) + " of " +str (len(url_list2))))

    path_out = os.path.join(fld_download_zip, ii)
    if os.path.exists(output_folder):
        print ("   Skip Downloading - File already exists")
    else:
         asf.download_url(url=ii, path=fld_download_zip, session=session)

    count = count + 1

print ("\t")
print(((time.strftime("%H:%M:%S")) + "   Scenes download completed"))


print(((time.strftime("%H:%M:%S")) + " Files unzipping"))


file_unzip_list = []    
for root, dirs, files in os.walk(fld_download_zip):
    for filename in fnmatch.filter(files, "*.zip"):
        file_unzip_list.append(os.path.join(root, filename))    
        
for ii in file_unzip_list:        
    with zipfile.ZipFile(ii, 'r') as zip_ref:
        zip_ref.extractall(fld_download_zip)
               
print(((time.strftime("%H:%M:%S")) + " Unzipping completed"))