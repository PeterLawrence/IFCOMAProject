"""
Created on Fri March 16 09:59:58 2025

@author: P.J.Lawrence
"""
import sys
import os
import ifcopenshell
from ifcopenshell import geom
import exogetifcparam
import ifc_config


def scan_spaces(ifc_file, settings):
    print("-------------- IfcSpace ---------------")
    ifc_spaces = ifc_file.by_type('IfcSpace')
    for ifc_space in ifc_spaces:
        output_name = False
        space_loc = -1
        Elevation = -1
        Element_list = [] # used to identify repeated elements
        for boundary in ifc_space.BoundedBy:
            if boundary.is_a("IfcRelSpaceBoundary"):
                if space_loc == -1:
                    Elevation = exogetifcparam.get_storey_Elevation_spaceboundary_quiet(boundary)
                    floor_longname = exogetifcparam.get_storey_Longname_spaceboundary(boundary)
                valid_item = True
                #if boundary.InternalOrExternalBoundary!='INTERNAL':
                #    valid_item = False
                if valid_item:
                    if not output_name:
                        output_name = True
                        print('===============================')
                        print('ifcSpace: ', ifc_space.Name, ' is at floor:', ifc_space.Decomposes[0][4][2],end=",")
                        print('IfcRelSpaceBoundary: ',boundary.Name)
                    elem = boundary.RelatedBuildingElement
                    if elem:
                        item_repeated=False
                        if elem.Name in Element_list:
                            item_repeated=True
                        else:
                            Element_list.append(elem.Name)
                         
                        elem_story = exogetifcparam.get_entity_storey_long_name(elem)
                        OutputElem = True
                        if elem_story!=floor_longname:
                            if item_repeated:
                                print("XR ",end="")
                            else:
                                print("XX ",end="")
                            print(elem.is_a(),":",end="")
                            print(elem.Name," is on Wrong Level ",elem_story,end="")
                            if item_repeated:
                                print(" and repeated",end="")
                            print()
                        else:
                            if item_repeated:
                                print("+R ",end="")
                            else:
                                print("++ ",end="")
                            print(elem.is_a(),":",end="")
                            print(elem.Name,end="")
                            if item_repeated:
                                print(" is repeated",end="")
                            print()
                    

def generate_data(IFC_filename):
    """
     main routine to generate output file/data

    """
    global g_OMA_Class
    global g_plotter_class
    
    print("Loading IFC file ", IFC_filename)
    ifc_file = ifcopenshell.open(IFC_filename)

    settings = geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)

    scan_spaces(ifc_file, settings)

def main(argv):
    """
     main routine

    """
    if len(argv)>0:
        # use data entered on command line
        IFCFilename = argv[0]
    else:
        # use data defined in ifc_config.py
        IFCFilename = ifc_config.IFCFile

    generate_data(IFCFilename)

    

if __name__ == "__main__":
    main(sys.argv[1:])



