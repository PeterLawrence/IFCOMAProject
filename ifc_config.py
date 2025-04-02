# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence
"""
import os

# the default input file names if no command line parameters  are provided

Folder = "TestCases"
FileName = "hotel-v12.ifc"
#FileName = "HM-NoPass4.ifc"
#InfaFileName = "Infra-Road.ifc"
#FileName = InfaFileName

IFCFile = Folder + os.sep + FileName

# the default output file names if no command line parameters  are provided
OutputFolder = "Outputs"
EZMapFile = OutputFolder + os.sep + "map.xml"
EZPopFile = OutputFolder + os.sep + "population.xml"
MTAFile = OutputFolder + os.sep + "IFCMeshLevels2.mta"
GraphFile = OutputFolder + os.sep + "BuildingGraph.html"
CFastFile = OutputFolder + os.sep + "Building.in"

# the default output file type if no command line parameters  are provided
OutputType = "Graph"  # "DumpStairs" "Exodus" "EvacutionZ" "Graph" "CFAST"

# show a basic plot of the model
plot_model = True
