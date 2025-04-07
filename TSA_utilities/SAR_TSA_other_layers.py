
'''----------------------------------------------------------------------
 * Gabriel Gosselin, CCRS                                                    -
 * ----------------------------------------------------------------------
'''

# The following locale settings are used to ensure that Python is configured
# the same way as PCI's C/C++ code.  
import locale

import pci.psinang

locale.setlocale(locale.LC_ALL, "")
locale.setlocale(locale.LC_NUMERIC, "C")

import os, fnmatch, time, sys

from pci.iia import iia
from pci.psinang import psinang
from pci.pcimod import pcimod

from pci.exceptions import PCIException
from pci.api import datasource as ds


# /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

def inc_angle_layer (input_files_list, output_folder_inc, prefix):

    print("\t")
    print("----------------------------------------------------------------------------------------------------------")
    print("                      Generating the stack  incidence angle layer                                         ")
    print("----------------------------------------------------------------------------------------------------------")
    print("\t")
    
    '''
    Note 1: All input sar scene have been coregistered to a common frame, they all have the same number of lines and
    columns. We assume here:
      1) That all input scenes have the same viewing geometry
      2) Their footprint is reasonably the same
      3) The incidence angle range is very simmilar for all scenes
    Thus, we create a single incidence layer for all the scenes from the first scene in the stack
    
    Note 2: The list structure is kept for future use if needed.
    '''
    for input_file in input_files_list: 

        with ds.open_dataset(input_file, ds.eAM_READ) as ds2:
            num_channels = ds2.chan_count

        # Creating the incidence angle layer
        print("   Generating the incidence angle layer")
        fili = input_file
        base = os.path.basename (input_file[:-4])
        filo = os.path.join (output_folder_inc,prefix + "Incidence_angle_layer.pix")
        angletyp = "degree"

        if os.path.exists (filo): 
            print ("Output incidence angle file already exists - skip")
        else: 
            try:
                psinang(fili, filo, angletyp)
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)
            angle_map = filo

            fili = input_file
            filo = angle_map
            dbsl = [2, -6]  # input segments to transfer
            dbos = []  # overwrite existing segment 1

            try:
                iia(fili, filo, dbsl, dbos)
            except PCIException as e:
                print(e)
            except Exception as e:
                print(e)

    return ()