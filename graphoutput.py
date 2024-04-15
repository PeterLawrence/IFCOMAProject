# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence
"""

from pyvis.network import Network
import exoifcutils
from oma_class import OMAClass


def graph_view(OMA_Class, ifc_file, filename):
    print("===================== output network ====================")
    space_free_links = []  # list of direct connections between spaces
    space_count = len(OMA_Class.m_space_list)
    for iSpace in range(space_count - 1):
        if 'boundarylist' in OMA_Class.m_space_list[iSpace]:
            for jSpace in range(iSpace + 1, space_count):
                if 'boundarylist' in OMA_Class.m_space_list[jSpace]:
                    # print('graph_view:',OMA_Class.m_space_list[iSpace]['Name'],OMA_Class.m_space_list[jSpace]['Name'])
                    if exoifcutils.are_adjacent(OMA_Class.m_space_list[iSpace]['boundarylist'], OMA_Class.m_space_list[jSpace]['boundarylist']):
                        print("Adjacent:", OMA_Class.m_space_list[iSpace]['Name'], OMA_Class.m_space_list[jSpace]['Name'])
                        OMA_Class.m_space_list[iSpace]['elemIDs'].append(["IfcSpace", OMA_Class.m_space_list[jSpace]['GlobalId']])
                        OMA_Class.m_space_list[jSpace]['elemIDs'].append(["IfcSpace", OMA_Class.m_space_list[iSpace]['GlobalId']])
                        space_free_links.append([OMA_Class.m_space_list[iSpace]['GlobalId'], OMA_Class.m_space_list[jSpace]['GlobalId']])

    net = Network()

    enz_node_id = 0
    for ifc_space in OMA_Class.m_space_list:
        output_node = False
        if 'elemIDs' in ifc_space:
            output_node = (len(ifc_space['elemIDs']) > 0)
        # removed 15/03/24 because ifc_stair is all that's required 
        # if output_node == False:
        #    output_node = ('stairflightsUp' in ifc_space or 'stairflightsDown' in ifc_space)
        if not output_node:
            output_node = ('stairsUp' in ifc_space or 'stairsDown' in ifc_space)
        if not output_node:
            output_node = ('escalatorsUp' in ifc_space or 'escalatorsDown' in ifc_space)
        if output_node:
            enz_node_id += 1
            ifc_space['enz_node_id'] = enz_node_id
            space_area = 1
            if 'area' in ifc_space:
                space_area = ifc_space['area']
            if 'Reference' in ifc_space:
                net.add_node(enz_node_id, label=ifc_space['Name'], title=ifc_space['Reference'], value=space_area,
                             shape="square")
            else:
                net.add_node(enz_node_id, label=ifc_space['Name'], value=space_area, shape="square")

    for ifc_door in OMA_Class.m_door_list:
        if ifc_door['IsExternal']:
            enz_node_id += 1
            ifc_door['enz_node_id'] = enz_node_id
            if 'RequiredDoorFlowrate' in ifc_door:
                net.add_node(enz_node_id, label=ifc_door['Name'], title=str(ifc_door['RequiredDoorFlowrate']),
                             color='#00ff1e')
            else:
                net.add_node(enz_node_id, label=ifc_door['Name'], color='#00ff1e')

    print("Max Enz nodes", enz_node_id)

    linkcount = 0
    for ifc_door in OMA_Class.m_door_list:
        if 'Spaces' in ifc_door:
            connectlist = []
            for spaceid in ifc_door['Spaces']:
                space_index = OMA_Class.SpaceDefined(spaceid)
                if space_index > - 1:
                    ifc_space = OMA_Class.m_space_list[space_index]
                    if 'enz_node_id' in ifc_space:
                        connectlist.append(ifc_space['enz_node_id'])
                    else:
                        print("rejected: door enz_node_id", space_index, ifc_space)
                else:
                    print("rejected: door index", space_index, spaceid)

            if len(connectlist) == 1:
                if ifc_door['IsExternal'] and 'enz_node_id' in ifc_door:
                    connectlist.append(ifc_door['enz_node_id'])
            if len(connectlist) == 2:
                net.add_edge(connectlist[0], connectlist[1], title=ifc_door['Name'])
                linkcount += 1
            else:
                if len(connectlist) == 4:  # temp fix for M_Single-Flush:Hotel Door:647825
                    net.add_edge(connectlist[2], connectlist[3], title=ifc_door['Name'])
                print(
                    f"Rejected: door connnect {ifc_door['Name']} ID:{ifc_door['GlobalId']} spaces {ifc_door['Spaces']} connections {connectlist}")
    
    for ifc_stair in OMA_Class.m_stair_list:
        if 'space_connecting' in ifc_stair:
            if len(ifc_stair['space_connecting']) == 2:
                SpaceIDs = []
                for SpaceGlobalID in ifc_stair['space_connecting']:
                    index = OMA_Class.SpaceDefined(SpaceGlobalID)
                    if index > -1:
                        ifc_space = OMA_Class.m_space_list[index]
                        if 'enz_node_id' in ifc_space:
                            SpaceIDs.append(ifc_space['enz_node_id'])
                if len(SpaceIDs) == 2:
                    print("Stair link", ifc_stair['Name'])
                    net.add_edge(SpaceIDs[0], SpaceIDs[1], title=ifc_stair['Name'], color='#FF0000', weight=20.0)
                    linkcount += 1

    for ifc_escalator in OMA_Class.m_escalator_list:
        if 'space_connecting' in ifc_escalator:
            if len(ifc_escalator['space_connecting']) == 2:
                SpaceIDs = []
                for SpaceGlobalID in ifc_escalator['space_connecting']:
                    index = OMA_Class.SpaceDefined(SpaceGlobalID)
                    if index > -1:
                        ifc_space = OMA_Class.m_space_list[index]
                        if 'enz_node_id' in ifc_space:
                            SpaceIDs.append(ifc_space['enz_node_id'])
                if len(SpaceIDs) == 2:
                    print("Escalator link", ifc_escalator['Name'])
                    net.add_edge(SpaceIDs[0], SpaceIDs[1], title=ifc_escalator['Name'], color='#FF00FF', weight=20.0)
                    linkcount += 1

    for space_link in space_free_links:
        index1 = OMA_Class.SpaceDefined(space_link[0])
        index2 = OMA_Class.SpaceDefined(space_link[1])
        print("space links", index1, index2, space_link[0], space_link[1])
        if index1 > -1 and index2 > -1:
            ifc_space1 = OMA_Class.m_space_list[index1]
            ifc_space2 = OMA_Class.m_space_list[index2]
            a_title = ifc_space1['Name'] + " to " + ifc_space2['Name']
            print(ifc_space1['enz_node_id'], ifc_space2['enz_node_id'])
            net.add_edge(ifc_space1['enz_node_id'], ifc_space2['enz_node_id'], title=a_title, color='#00FF00')
            linkcount += 1

    # virtual boundaries

    print("virtual boundaries")
    for ifc_space in OMA_Class.m_space_list:
        if 'elemIDs' in ifc_space:
            boundary_elems = ifc_space['elemIDs']
            for item in boundary_elems:
                if item[0] == "IfcRelSpaceBoundaryV":
                    print(ifc_space['Name'], item[0], item[1])
                    vb = ifc_file.by_id(item[1])
                    print("RelatedBuildingElement", vb.RelatedBuildingElement)
                elif item[0] == "IfcVirtualElement":
                    print("IfcVirtualElement", item[0], item[1])

    print("Max Enz links", linkcount)
    net.show_buttons(filter_=['physics'])
    net.show(filename, notebook=False)
