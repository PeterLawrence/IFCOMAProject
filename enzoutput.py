# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@brief Data output for EvacuationZ, see https://evacuationz.wordpress.com/
@note to install winget install Evacuationz 

@author: P.J.Lawrence
"""

import xml.etree.ElementTree as ET
import exoifcutils
from oma_class import OMAClass
import re


def output_node_enz(root, ifc_space):
    Node = ET.SubElement(root, "Node")
    ET.indent(root, space="\t", level=1)

    ET.SubElement(Node, "Name").text = ifc_space['Reference']
    ET.SubElement(Node, "Ref").text = str(ifc_space['enz_node_id'])
    ET.SubElement(Node, "Area").text = str(ifc_space['area'])
    if 'objects' in ifc_space: # also 'furniture'
        furniture = ET.SubElement(Node, "Furniture")
        ET.indent(furniture, space="\t", level=2)
        for object in ifc_space['objects']:
            ET.SubElement(furniture, "Object").text = object
        


def output_external_exit_enz(root, ifc_door):
    Node = ET.SubElement(root, "Node", type='enz_safe')
    ET.indent(root, space="\t", level=1)

    ET.SubElement(Node, "Name").text = ifc_door['Name']
    ET.SubElement(Node, "Ref").text = str(ifc_door['enz_node_id'])


def output_connection_link_enz(root, connectlist, width, name):
    Connection = ET.SubElement(root, "Connection")
    ET.indent(root, space="\t", level=1)

    ET.SubElement(Connection, "Name").text = name
    ET.SubElement(Connection, "NodeRef", refstyle="enz_ref").text = str(connectlist[0])
    ET.SubElement(Connection, "NodeRef", refstyle="enz_ref").text = str(connectlist[1])

    ConnectionType = ET.SubElement(Connection, "ConnectionType", type = 'enz_door')
    ET.indent(Connection, space="\t", level=2)

    if width>0:
        ET.SubElement(ConnectionType, "Width", units = 'm').text = str(width)


def output_connection_stair_enz(root, ifc_stair, space_list):
    SpaceIDs = []
    for SpaceGlobalID in ifc_stair['space_connecting']:
        index = exoifcutils.SpaceDefined(SpaceGlobalID,space_list)
        if index>-1:
            ifc_space = space_list[index]
            if 'enz_node_id' in ifc_space:
                SpaceIDs.append(ifc_space['enz_node_id'])

    if len(SpaceIDs)==2:
        Connection = ET.SubElement(root, "Connection")
        ET.indent(root, space="\t", level=1)

        ET.SubElement(Connection, "Name").text = ifc_stair['Name']

        ET.SubElement(Connection, "NodeRef").text = str(SpaceIDs[0])
        ET.SubElement(Connection, "NodeRef").text = str(SpaceIDs[1])
        ET.SubElement(Connection, "Length").text = str(ifc_stair['TravelDistance'])

        ConnectionType = ET.SubElement(Connection, "ConnectionType", type = 'enz_stairs')
        ET.indent(Connection, space="\t", level=2)
        
        ET.SubElement(ConnectionType, "Tread", units = 'm').text = str(ifc_stair['TreadLength'])
        ET.SubElement(ConnectionType, "Riser", units = 'm').text = str(ifc_stair['RiserHeight'])
        ET.SubElement(ConnectionType, "Width", units = 'm').text = str(ifc_stair['Width'])
    else:
        print("Rejected Door ",ifc_stair['Name'], ifc_stair['space_connecting'])


def output_connection_escalator_enz(root, ifc_escalator, space_list):
    SpaceIDs = []
    for SpaceGlobalID in ifc_escalator['space_connecting']:
        index = exoifcutils.SpaceDefined(SpaceGlobalID,space_list)
        if index>-1:
            ifc_space = space_list[index]
            if 'enz_node_id' in ifc_space:
                SpaceIDs.append(ifc_space['enz_node_id'])

    if len(SpaceIDs)==2:
        Connection = ET.SubElement(root, "Connection")
        ET.indent(root, space="\t", level=1)

        ET.SubElement(Connection, "Name").text = ifc_escalator['Name']

        ET.SubElement(Connection, "NodeRef").text = str(SpaceIDs[0])
        ET.SubElement(Connection, "NodeRef").text = str(SpaceIDs[1])
        if 'OverallLength' in ifc_escalator: # travel length
            ET.SubElement(Connection, "Length").text = str(ifc_escalator['OverallLength'])
        elif 'RunLength' in ifc_escalator: # horizonal length
            ET.SubElement(Connection, "RunLength").text = str(ifc_escalator['OverallLength'])

        ConnectionType = ET.SubElement(Connection, "ConnectionType", type = 'enz_stairs')
        ET.indent(Connection, space="\t", level=2)
        
        ET.SubElement(ConnectionType, "Tread", units = 'm').text = str(ifc_escalator['TreadLength'])
        ET.SubElement(ConnectionType, "Riser", units = 'm').text = str(ifc_escalator['RiserHeight'])
        ET.SubElement(ConnectionType, "Width", units = 'm').text = str(ifc_escalator['Width'])
    else:
        print("Rejected Door ",ifc_escalator['Name'], ifc_escalator['space_connecting'])


def output_node_pop_enz(root, ifc_space):
    people_allocated = 0
    if 'objects' in ifc_space: # also 'furniture'
        # count beds
        for object in ifc_space['objects']:
            if re.search('bed', object, re.IGNORECASE):
                if re.search('double', object, re.IGNORECASE): 
                    people_allocated+=2
                else:
                    people_allocated+=1

        # if no beds count chairs
        if people_allocated==0:
            for object in ifc_space['objects']:
                if re.search('chair', object, re.IGNORECASE): 
                    people_allocated+=1
        
    if people_allocated>0:
        PopulationDefinition = ET.SubElement(root, "PopulationDefinition")
        ET.indent(root, space="\t", level=1)
        ET.SubElement(PopulationDefinition, "Agents").text = str(people_allocated)
        ET.SubElement(PopulationDefinition, "NodeRef", refstyle="enz_ref").text = str(ifc_space['enz_node_id'])

    return people_allocated;
        

def enz_population(XMLFile, OMA_Class):
    total_people_count = 0
    root = ET.Element("EvacuatioNZ_Populate")
    for ifc_space in OMA_Class.m_space_list:
        if 'enz_node_id' in ifc_space:
            if ifc_space['enz_node_id']>-1:
                total_people_count+=output_node_pop_enz(root, ifc_space)
            

    print("Total people Count",total_people_count);

    tree = ET.ElementTree(root)
    tree.write(XMLFile)
            

def enz_output(XMLFileMap,XMLFilePop, OMA_Class):
    root = ET.Element("ENZ_Map")

    Description = ET.SubElement(root, "Description").text = "IFC Export file"
    ET.indent(root, space="\t", level=0)

    enz_node_id = 0
    for ifc_space in OMA_Class.m_space_list:
        output_node = False
        if 'elemIDs' in ifc_space:
            output_node = (len(ifc_space['elemIDs'])>0)
        if output_node == False:
            output_node = ('stairflightsUp' in ifc_space or 'stairflightsDown' in ifc_space)
        if output_node:
            enz_node_id+=1
            ifc_space['enz_node_id'] = enz_node_id
            output_node_enz(root, ifc_space)

    for ifc_door in OMA_Class.m_door_list:
        if ifc_door['IsExternal']:
            enz_node_id+=1
            ifc_door['enz_node_id'] = enz_node_id
            output_external_exit_enz(root, ifc_door)

    print("Max Enz nodes",enz_node_id);

    linkcount = 0
    for ifc_door in OMA_Class.m_door_list:
        if 'Spaces' in ifc_door:
            connectlist = []
            for spaceid in ifc_door['Spaces']:
                space_index = exoifcutils.SpaceDefined(spaceid,OMA_Class.m_space_list)
                if space_index>-1:
                    ifc_space = OMA_Class.m_space_list[space_index]
                    if 'enz_node_id' in ifc_space:
                        connectlist.append(ifc_space['enz_node_id'])
                    else:
                        print("rejected: door enz_node_id",space_index, ifc_space)
                else:
                    print("rejected: door index",space_index, spaceid)
                    
            OverallWidth=-1.0
            if 'OverallWidth' in ifc_door:
                OverallWidth = ifc_door['OverallWidth']
            if len(connectlist)==1:
                if ifc_door['IsExternal'] and 'enz_node_id' in ifc_door:
                    connectlist.append(ifc_door['enz_node_id'])
            if len(connectlist)==2:
                output_connection_link_enz(root,connectlist,OverallWidth,ifc_door['Name'])
                linkcount+=1
            else:
                print(f"Rejected: door connnect {ifc_door['Name']} ID:{ifc_door['GlobalId']} spaces {ifc_door['Spaces']} connections {connectlist}")


    for ifc_stair in OMA_Class.m_stair_list:
        if 'space_connecting' in ifc_stair:
            if len(ifc_stair['space_connecting'])==2:
                output_connection_stair_enz(root,ifc_stair,OMA_Class.m_space_list)
                linkcount+=1

    for ifc_escalator in OMA_Class.m_escalator_list:
        if 'space_connecting' in ifc_escalator:
            if len(ifc_escalator['space_connecting']) == 2:
                output_connection_escalator_enz(root,ifc_escalator,OMA_Class.m_space_list)
                linkcount+=1

    print("Max Enz links",linkcount);

    tree = ET.ElementTree(root)
    tree.write(XMLFileMap)

    if XMLFilePop is not None:
        enz_population(XMLFilePop, OMA_Class)
