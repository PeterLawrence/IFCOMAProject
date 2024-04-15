# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence
"""

# the default input file names if no command line parameters  are provided

#FileName = "Office_Building-ArchitecturalModel.ifc"
#Folder = "E:\\Support\BIM\\2011-09-14-Office-IFC\\Model-files\\"

#FileName = "Clinic_Handover.ifc"
#Folder = "E:\\Support\BIM\\2011-09-14-Clinic-IFC\\Model files\\"

#Folder = "E:\\Support\\BIM\\Hotel\\v6\\"
#FileName = "hotel-v6-4x3-level2-SBs-3xRequiredDoorFlowRate.ifc"

#Folder = "E:\\Support\\BIM\\Hotel\\v7\\"
#FileName = "hotel-v7-4x3-level2-SBs-DoorParams.ifc"

Folder = "E:\\Support\\BIM\\Hotel\\"
FileName = "hotel-v12.ifc"

#Folder = "E:\\Support\\BIM\\TestModels\\"
#FileName = "Hospital-Original-withIntDoors-InitialDesign.ifc" # No spaces defined
#FileName = "StationNightclub-NoRoofVersion-final2.ifc" # No spaces defined
#FileName = "StationNightclub-final2.ifc"

IFCFile = Folder + FileName

# the default output file names if no command line parameters  are provided
EZFile = ".\\EZtest2.xml"
MTAFile = ".\\IFCMeshLevels2.mta"
GraphFile = "..\\BuildingGraph.html"

# the default output file type if no command line parameters  are provided
OutputType = "Graph"  # "DumpStairs" "Exodus" "EvacutionZ" "Graph"

# show a basic plot of the model
plot_model = True
