
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
from pci.model import model
from pci.fexport import fexport

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
    print (time.strftime("%H:%M:%S") + "  Generating the incidence angle layer")
    for input_file in input_files_list: 

        with ds.open_dataset(input_file, ds.eAM_READ) as ds2:
            num_channels = ds2.chan_count

        # Creating the incidence angle layer
        fili = input_file
        base = os.path.basename (input_file[:-4])
        filo = os.path.join (output_folder_inc,prefix + "Incidence_angle_layer.pix")
        angletyp = "degree"

        if os.path.exists (filo): 
            print ("Output incidence angle file already exists - skip")
        else: 
            print("   output file-->" + filo)
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


def math_layers_prod (TSA_math_models_definition, TSA_math_xtra_channels, TSA_math_xtra_labels, Fld_math_layers, 
                      Fld_math_layers_ortho, prefix):

    print("\t")
    print("----------------------------------------------------------------------------------------------------------")
    print("                                   Generating the math layers                                             ")
    print("----------------------------------------------------------------------------------------------------------")
    print("\t")

    # A)  We retrieve the scenes from the Fld_math_layers
    input_files_list = []
    for root, dirs, files in os.walk(Fld_math_layers):
        for filename in fnmatch.filter(files, "*_int.pix"):
            input_files_list.append(os.path.join(root, filename))

    # B) Adding empty channels to receive the model results
    nb_files = str(len(input_files_list))
    count = 1
    for input_file in input_files_list: 
        print ("\t")
        print (time.strftime("%H:%M:%S") + " Generating math layers, file " + str(count) + " of " + nb_files)
        print ("   input file-->" + input_file)
        print ("   Adding empty 32R channel(s)")
        file = input_file
        pciop = "ADD"
        pcival =[0,0,0,TSA_math_xtra_channels]

        try: 
            pcimod(file, pciop, pcival)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)


        print ("   Applying the model(s) to generate the math layer(s)")
        file = input_file
        source = TSA_math_models_definition
        undefval = [-32768.00]
        try: 
            model(file, source, undefval)
        except PCIException as e:
            print(e)
        except Exception as e:
            print(e)

        count = count +1
        # Code to add proper name to the generated channels
    return()

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def math_layers_split (Fld_math_layers_ortho, output_folder, prefix, TSA_math_xtra_channels,
                      TSA_math_xtra_labels, suffix):
    
    # A) Find the files to process
    input_files_list = []
    for root, dirs, files in os.walk(Fld_math_layers_ortho):
        for filename in fnmatch.filter(files, "*_int.pix"):
            input_files_list.append(os.path.join(root, filename))

    #chan_list = []
    #chan = 1
    #while chan <= TSA_math_xtra_channels: 
        #chan_list.append (chan)
       # chan = chan + 1

    # B) Rename the channels using mdmap
    chan_list = list(range(1, TSA_math_xtra_channels + 1))

    for in_file in input_files_list:     
        with ds.open_dataset(in_file, ds.eAM_WRITE) as ds5:
            # get the AuxiliaryData
            aux = ds5.aux_data

            for in_chan, in_label in zip (chan_list, TSA_math_xtra_labels): 

                # get the metadata for the first channel
                mdmap = aux.get_chan_metadata(in_chan)
                # set the NO_DATA_VALUE to -32768
                #mdmap['MID_DATE'] = mid_date
                #mdmap['REF_DATE'] = ref_date
                #mdmap['DEP_DATE'] = dep_date
                chan_desc = in_label
                # set the metadata map back to the AuxiliaryData
                aux.set_chan_description(chan_desc, in_chan)
                aux.set_chan_metadata(mdmap, in_chan)

                # set the AuxiliaryData back to the dataset
                ds5.aux_data = aux


    for in_file in input_files_list:
        for in_chan, in_label in zip (chan_list, TSA_math_xtra_labels):
            print (in_chan)
            print (in_label)
            # export the single channel
            fili =	in_file
            base = os.path.basename(in_file[:-7])
            outname = (base + in_label + ".pix")
            print ("   output_file_name-->" + outname)
            filo = os.path.join (output_folder, outname)
            dbiw =	[]
            dbic =	[in_chan]
            dbib =	[]
            dbvs =	[]
            dblut =	[]
            dbpct =	[]
            ftype =	"PIX"
            foptions	= ""
            try: 
                fexport( fili, filo, dbiw, dbic, dbib, dbvs, dblut, dbpct, ftype, foptions )
            except PCIException as e:
                print(e)
            except Exception as e:
                 print(e)
    return()
