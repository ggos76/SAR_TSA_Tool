
'''-----------------------------------------------------------------------------------------------
 * Gabriel Gosselin, CCRS   2024-2025                                                            -
 * -----------------------------------------------------------------------------------------------
'''
import sys

# ------------------------------------------------------------------------------------------------
def version_control (vs_catalyst, vs_python):
    
    min_python_vs = (3,10,0) 
    print ("Minimum Python version to run the TSA scripts: " + str(min_python_vs))
    print ("Current Python version: " + str(vs_python))
    if min_python_vs < vs_python: 
        print ("Python version requirement is meet")
    else: 
        print ("Error - The minimum python version requirement is not meet")
        sys.exit()
     
    min_catalyst_vs =(3,1,0)
    version_tuple = tuple(map(int, vs_catalyst.split(".")))
    vs_catalyst = version_tuple
    print ("Minimum Catalyst install to run the TSA scripts: " + str(min_catalyst_vs))
    print ("Current Catalyst install: " + str(vs_catalyst))
    if min_catalyst_vs < vs_catalyst:
        print ("Catalyst version requirement is meet")
    else:
        print ("Error - The minimum Catalyst version requirement is not meet")
        sys.exit()

    return ()