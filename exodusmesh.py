"""
Created on Fri Dec 9 10:59:58 2023

@author: P.J.Lawrence
"""

import ifcopenshell
from ifcopenshell import geom
import numpy as np
import math
import os
import exogetifcparam

import exooutput
import exoifcutils

from  plot_functions import PlotterClass
from oma_class import OMAClass

from shapely.geometry import MultiPoint, Polygon

# Used to position the geometry in EXODUS
g_rangeX = [0, -1]
g_rangeY = [0, -1]
g_rangeZ = [0, -1]

g_DoorID = 0  # Note external doors comes before g_NodeID
g_NodeStartID = 5  # allow space for 4 external doors
g_NodeID = g_NodeStartID

g_StairID = 1

g_space_node_list = []  # node lists for each room, is made up of a list of rows
g_subp_list = []

g_door_connections = []  # nodal connections between doors [ [NodeID, DoorId]... ]
g_stairflight_connections = []  # nodal connections between stairflight [ [NodeID, TransitNodeId]... ]

# boundary data
g_lines = []


def calc_door_box(shape):
    """ Calc the bounding box around a door for connecting doors
        Only works for horizonal and vertical doors in the x,y direction,
        see calc_door_region
    """
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    xp = []
    yp = []
    zp = []
    VertexCount = len(verts)
    iPoint = 0
    xmin = verts[iPoint]
    xmax = xmin
    iPoint += 1
    ymin = verts[iPoint]
    ymax = ymin
    iPoint += 1
    zmin = verts[iPoint]
    iPoint += 1

    while iPoint < VertexCount:
        xp.append(verts[iPoint])
        iPoint += 1
        yp.append(verts[iPoint])
        iPoint += 1
        zp.append(verts[iPoint])
        zmin = min(zmin, verts[iPoint])
        iPoint += 1
   
    myPoints = []
    for iPoint in range(0, len(xp)):
        if zp[iPoint] == zmin:
            aPoint = [xp[iPoint], yp[iPoint]]
            xmin = min(xmin, xp[iPoint])
            xmax = max(xmax, xp[iPoint])
            ymin = min(ymin, yp[iPoint])
            ymax = max(ymax, yp[iPoint])
            myPoints.append(aPoint)
        
    # since we know test case the doors are vert/horz. quick fix
    door_points = []
    door_box = []
    boxwidth = 0.60
    if ymax-ymin > xmax-xmin:
        door_points.append([(xmax+xmin)/2.0, ymin])
        door_points.append([(xmax+xmin)/2.0, ymax])
        door_box.append([(xmax+xmin)/2.0+boxwidth, ymax])
        door_box.append([(xmax+xmin)/2.0-boxwidth, ymax])
        door_box.append([(xmax+xmin)/2.0-boxwidth, ymin])
        door_box.append([(xmax+xmin)/2.0+boxwidth, ymin])
    else:
        door_points.append([xmin, (ymax+ymin)/2.0])
        door_points.append([xmax, (ymax+ymin)/2.0])
        door_box.append([xmax, (ymax+ymin)/2.0+boxwidth])
        door_box.append([xmin, (ymax+ymin)/2.0+boxwidth])
        door_box.append([xmin, (ymax+ymin)/2.0-boxwidth])
        door_box.append([xmax, (ymax+ymin)/2.0-boxwidth])
    print(door_points)
    print(door_box)
    #tri_idx = [[0, 1, 2], [2, 3, 0]]
    #xp = [door_box[0][0], door_box[1][0], door_box[2][0], door_box[3][0]]
    #yp = [door_box[0][1], door_box[1][1], door_box[2][1], door_box[3][1]]
    #zp = [zmin+boxwidth, zmin+boxwidth, zmin+boxwidth, zmin+boxwidth]
    print("===========================================================")
    print("door_box:", door_box)
    return door_box


def calc_door_region(shape):
    """ Calc the region around a door for connecting doors """
    verts = shape.verts # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    VertexCount = len(verts)
    iPoint = 0
    x = verts[iPoint]
    iPoint += 1
    y = verts[iPoint]
    iPoint += 1
    p1 = (x, y)
    p2 = p1
    distance = 0
    iPoint += 1

    # search for the two points furthest apart in 2D x,y space
    while iPoint < VertexCount:
        x = verts[iPoint]
        iPoint += 1
        y = verts[iPoint]
        iPoint += 1
        p3 = (x, y)
        iPoint += 1
        distancep1_3 = exoifcutils.calcDistance2D(p1, p3)
        distancep2_3 = exoifcutils.calcDistance2D(p2, p3)
        if distancep1_3 > distance or distancep2_3 > distance:
            if distancep1_3 > distance and distancep2_3 > distance:
                if distancep1_3 > distancep2_3:
                    p2 = p3
                    distance = distancep1_3
                else:
                    p1 = p3
                    distance = distancep2_3
            else:
                if distancep1_3 > distance:
                    p2 = p3
                    distance = distancep1_3
                else:
                    p1 = p3
                    distance = distancep2_3
                    
    extra_size = 0.5
    set1 = exoifcutils.CalcPositionSides(p1, p2, distance, extra_size)
    set2 = exoifcutils.CalcPositionSides(p2, p1, distance, extra_size)
    # TO DO - ensure points in an anti-clockwise directions
    DoorRegion = [set1[1], set1[0], set2[1], set2[0]]
    print("Points:", p1, p2, "distance", distance)
    print("Region:", DoorRegion)
    print("===========================================================")
    return DoorRegion
        

def calc_flightstair_box(ifc_stairflight):
    """ Calc the area box for connecting a stair case """
    if 'PathData' in ifc_stairflight:
        datalen = len(ifc_stairflight['PathData'])
        if datalen > 1:
            #print("PathData :",ifc_stairflight['PathData'])
            print("===========================================================")
            width = ifc_stairflight['Width']/2.0
            expand_size = ifc_stairflight['treadlength']/2.0
            extra_size = 0.5
            path_data = ifc_stairflight['PathData']
            level_data = ifc_stairflight['AStepData']
            # bottom
            BotPoint1 = [path_data[1][0], path_data[1][1]]
            BotPoint2 = [path_data[0][0], path_data[0][1]]
            BotElvation = level_data[0][0]
            BotLength = exoifcutils.calcDistance2D(BotPoint1, BotPoint2)

            BotlineAOut = exoifcutils.CalcPosition(BotPoint1, BotPoint2, BotLength, expand_size)
            BotlineBOut = exoifcutils.CalcPosition(BotPoint1, BotPoint2, BotLength, expand_size+extra_size)

            BotLength = exoifcutils.calcDistance2D(BotlineAOut, BotlineBOut)
            set1 = exoifcutils.CalcPositionSides(BotlineAOut, BotlineBOut, BotLength, width)
            set2 = exoifcutils.CalcPositionSides(BotlineBOut, BotlineAOut, BotLength, width)
            # TO DO - ensure points in an anti-clockwise directions
            LowBox = [set1[0], set1[1], set2[0], set2[1]]
            print("Points:", BotPoint1, BotPoint2)
            print("Extended Points:", BotlineAOut, BotlineBOut)
            print("Lower box", LowBox)
            # top
            TopPoint1 = [path_data[datalen-2][0], path_data[datalen-2][1]]
            TopPoint2 = [path_data[datalen-1][0], path_data[datalen-1][1]]
            TopElvation = level_data[len(level_data)-1][0]
            TopLength = exoifcutils.calcDistance2D(TopPoint1, TopPoint2)

            ToplineAOut = exoifcutils.CalcPosition(TopPoint1, TopPoint2, TopLength, expand_size)
            ToplineBOut = exoifcutils.CalcPosition(TopPoint1, TopPoint2, TopLength, expand_size+extra_size)
            
            TopLength = exoifcutils.calcDistance2D(ToplineAOut, ToplineBOut)
            set1 = exoifcutils.CalcPositionSides(ToplineAOut, ToplineBOut, TopLength, width)
            set2 = exoifcutils.CalcPositionSides(ToplineBOut, ToplineAOut, TopLength, width)
            # TO DO - ensure points in an anti-clockwise directions
            TopBox = [set1[0], set1[1], set2[0], set2[1]]
            print("Points:", TopPoint1, TopPoint2)
            print("Extended Points:", ToplineAOut, ToplineBOut)
            print("Upper box", TopBox)
            print("===========================================================")
            if BotElvation < TopElvation:
                return LowBox, TopBox
            else:
                return TopBox, LowBox


def calc_flightstair_direction(ifc_stairflight):
    """
    Calc the direction of a stair
    """
    print("Calc stair direction for ", ifc_stairflight['Name'])
    if 'PathData' in ifc_stairflight:
        if len(ifc_stairflight['PathData']) > 1:
            path_data = ifc_stairflight['PathData']
            lineA = ((path_data[0][0], path_data[0][1]), (path_data[0][0], path_data[0][1]+10.0))
            lineB = ((path_data[0][0], path_data[0][1]), (path_data[1][0], path_data[1][1]))
            anAngle = exoifcutils.calcAngle(lineA, lineB)
            if anAngle < 0:
                anAngle = 180 - anAngle              
            return anAngle
    return ifc_stairflight['Direction']


def calc_flightstair_physical_param(ifc_stairflight):
    travel_length = 0.0
    horizontal_length = 0.0
    if 'PathData' in ifc_stairflight:
        path_data = ifc_stairflight['PathData']
        if len(path_data) > 1:
            for i in range(1, len(path_data)):
                travel_length += exoifcutils.calcDistance3D(path_data[i-1], path_data[i])
                horizontal_length += exoifcutils.calcDistance2D(path_data[i-1], path_data[i])
    else:
        travel_length = ifc_stairflight['Length']
        horizontal_length = ifc_stairflight['HorizontalLength']
                    
    return horizontal_length, travel_length


def add_lines(shape, WinID):
    global g_lines
    edges = shape.edges
    verts = shape.verts
    xp = []
    yp = []
    zp = []
    VertexCount = len(verts)
    iPoint = 0
    while iPoint < VertexCount:
        xp.append(verts[iPoint])
        iPoint += 1
        yp.append(verts[iPoint])
        iPoint += 1
        zp.append(verts[iPoint])
        iPoint += 1
    edgecount = len(edges)
    iPoint = 0
    line_count = 0
    while iPoint < edgecount:
        p1 = edges[iPoint]
        iPoint += 1
        p2 = edges[iPoint]
        iPoint += 1
        line = [xp[p1], yp[p1], xp[p2], yp[p2], True, WinID]
        g_lines.append(line)
        line_count += 1
    return line_count


def add_element_tri_gen_exodus_nodes(shape, floorLevel):
    faces = shape.faces  # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    
    iPoint = 0
    xp = []
    yp = []
    zp = []
    VertexCount = len(verts)
    while iPoint < VertexCount:
        xp.append(verts[iPoint])
        iPoint += 1
        yp.append(verts[iPoint])
        iPoint += 1
        zp.append(verts[iPoint])
        iPoint += 1

    lminX = min(xp)
    lmaxX = max(xp)
    lminY = min(yp)
    lmaxY = max(yp)
    lminZ = min(zp)
    lmaxZ = max(zp)
        
    iPoint = 0
    tri_idx = []
    FaceCount = len(faces)
    while iPoint < FaceCount:
        ixp = int(faces[iPoint])
        iPoint += 1
        iyp = int(faces[iPoint])
        iPoint += 1
        izp = int(faces[iPoint])
        iPoint += 1
        tri_idx.append((ixp, iyp, izp))

    return gen_exodus_nodes(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ, None, floorLevel)


def add_element_tri_gen_exodus_nodes_h(shape, Elev, Name):
    global g_rangeX, g_rangeY, g_rangeZ
    faces = shape.faces  # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    
    iPoint = 0
    xp = []
    yp = []
    zp = []
    VertexCount = len(verts)
    while iPoint < VertexCount:
        xp.append(verts[iPoint])
        iPoint += 1
        yp.append(verts[iPoint])
        iPoint += 1
        zp.append(verts[iPoint])
        iPoint += 1

    lminX = min(xp)
    lmaxX = max(xp)
    lminY = min(yp)
    lmaxY = max(yp)
    lminZ = min(zp)
    lmaxZ = max(zp)
    
    if g_rangeX[1] < g_rangeX[0]:
        g_rangeX = [lminX, lmaxX]
        g_rangeY = [lminY, lmaxY]
        g_rangeZ = [lminZ, lmaxZ]
    else:
        g_rangeX = [min(lminX, g_rangeX[0]), max(lmaxX, g_rangeX[1])]
        g_rangeY = [min(lminY, g_rangeY[0]), max(lmaxY, g_rangeY[1])]
        g_rangeZ = [min(lminZ, g_rangeZ[0]), max(lmaxZ, g_rangeZ[1])]

    iPoint = 0
    tri_idx = []
    FaceCount = len(faces)
    while iPoint < FaceCount:
        ixp = int(faces[iPoint])
        iPoint += 1
        iyp = int(faces[iPoint])
        iPoint += 1
        izp = int(faces[iPoint])
        iPoint += 1
        if zp[izp] > Elev-0.1 and zp[izp] < Elev+0.1 and zp[iyp] > Elev-0.1 and zp[iyp] < Elev+0.1 and zp[ixp] > Elev-0.1 and zp[ixp] < Elev+0.1:
            tri_idx.append((ixp, iyp, izp))
        
    return gen_exodus_nodes(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ, Name, Elev)


def gen_exodus_nodes(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ,node_name, floorElevation):
    global g_NodeID
    global g_space_node_list
    global g_door_connections
    node_list = []
    
    lx = 0
    row_length = 0
    node_count = 0
    while lx < lmaxX:
        row_list = []
        if lx >= lminX:
            ly = 0
            while ly < lmaxY:
                if ly > lminY:
                    iFace = 0
                    AddPoint = False
                    for triface in tri_idx:
                        if exoifcutils.TriangleOrientation(xp[triface[0]], yp[triface[0]], xp[triface[1]], yp[triface[1]], xp[triface[2]], yp[triface[2]])==1:
                            AddPoint = exoifcutils.intri(xp[triface[0]], yp[triface[0]], xp[triface[1]], yp[triface[1]], xp[triface[2]], yp[triface[2]], lx, ly)
                        else:
                            AddPoint = exoifcutils.intri(xp[triface[2]], yp[triface[2]], xp[triface[1]], yp[triface[1]], xp[triface[0]], yp[triface[0]], lx, ly)
                        if AddPoint:
                            break
                        
                    if AddPoint:
                        NdeHeight = round(max(0.0, floorElevation-zp[triface[0]]), 4)
                        node_count += 1
                        if node_name:
                            row_list.append([g_NodeID, lx, ly, NdeHeight, node_name])
                        else:
                            row_list.append([g_NodeID, lx, ly, NdeHeight])
                        g_NodeID += 1
                    else:
                        row_list.append([-1, lx, ly, -1])
                ly += 0.5
        if len(row_list) > 0:
            node_list.append(row_list)
        lx += 0.5

    g_space_node_list.append(node_list)
    return node_count


def nodal_mesh_space_data(space_data, ifc_space, settings):
    global g_space_node_list
    NumNodes = 0
    space_data['subspaces'] = []
    
    if ifc_space:
        ifc_space_shape = ifcopenshell.geom.create_shape(settings, ifc_space)
        if ifc_space_shape:
            Elevation = space_data['Elevation']
            print("Adding nodes for space ", space_data['Name'], Elevation, space_data['FloorIndex'])
            NumNodes = add_element_tri_gen_exodus_nodes_h(ifc_space_shape.geometry, Elevation, space_data['Reference'])
            space_data['NodeCount'] = NumNodes
            if NumNodes > 0:
                space_data['NodeLocation'] = len(g_space_node_list)-1
            else:
                space_data['NodeLocation'] = -1
            
    return


def nodal_mesh_landing_data(landing_data, ifc_landing, settings):
    NumNodes = 0
    landing_data['subspaces'] = []
    global g_space_node_list
    
    if ifc_landing:
        ifc_landing_shape = ifcopenshell.geom.create_shape(settings, ifc_landing)
        if ifc_landing_shape:
            Elevation = landing_data['Height']
            print("Adding nodes for space ", landing_data['Name'], Elevation, landing_data['FloorIndex'])
            NumNodes = add_element_tri_gen_exodus_nodes_h(ifc_landing_shape.geometry, Elevation, landing_data['Name'])
            landing_data['NodeCount'] = NumNodes
            if NumNodes > 0:
                landing_data['NodeLocation'] = len(g_space_node_list)-1
            else:
                landing_data['NodeLocation'] = -1
            
    return


def add_floor_nodes(space_list, ifc_file, settings):
    global g_subp_list
    global g_lines
    print(' Meshing region')
    for space in space_list:
        ifc_space = ifc_file.by_id(space['GlobalId'])
        nodal_mesh_space_data(space, ifc_space, settings)
        for ifc_space_boundary in ifc_space.BoundedBy:
            if ifc_space_boundary.is_a("IfcRelSpaceBoundary"):
                if True:  #for ifc_space_boundary in ifc_space_boundarys:
                    ifc_relshape = ifc_space_boundary.ConnectionGeometry.SurfaceOnRelatingElement
                    
                    if ifc_relshape.is_a('IfcCurveBoundedPlane') and ifc_relshape.InnerBoundaries is None:
                        ifc_relshape.InnerBoundaries = ()
                    try:
                        entity_shape_geom = ifcopenshell.geom.create_shape(settings, ifc_relshape)
                    except:
                        entity_shape_geom=None
                            
                    if entity_shape_geom is not None:
                        point_count = len(entity_shape_geom.verts) // 3

                        plot_item = True
                        # just want to plot horizontal stuff
                        # which has a zero height, since Z value is relative to floor on
                        # This ensures only floor IfcRelSpaceBoundarys are used
                        for ip in range(0, point_count):
                            aPoint = ip*3+2
                            if entity_shape_geom.verts[aPoint] != 0:
                                plot_item = False
                                break
                            
                        if plot_item:
                            #print("Add Lines for ",space['Name'],space['Elevation'],space['FloorIndex'])
                            subspace_data = {}
                            exogetifcparam.get_basic_subspace_info(ifc_space_boundary, subspace_data)

                            subspace_data['LineIndex'] = len(g_lines)
                            line_count = add_lines(entity_shape_geom, space['FloorIndex'])
                            subspace_data['LineCount'] = line_count
                            SpaceIndex = exoifcutils.SpaceDefined(space['GlobalId'], space_list)
                            subspace_data['ParentIndex'] = SpaceIndex
                            
                            g_subp_list.append(subspace_data)
                            space['subspaces'].append(len(g_subp_list)-1)
    print("Done")


def add_landing_nodes(landings_list, ifc_file, settings):
    global g_subp_list
    global g_lines
    print(' Meshing landing')
    for landing in landings_list:
        ifc_landing = ifc_file.by_id(landing['GlobalId'])
        nodal_mesh_landing_data(landing, ifc_landing, settings)

    print("Done")


def connect_door(Space, NodeList, door_box, DoorID):
    global g_door_connections
    print("connect_door External Door Box ", door_box, " Space ", Space['Name'])
    DoorList = []
    for aRow in NodeList:
        row_length = len(aRow)
        aColID = 0
        for aNode in aRow:
            if aNode[0] > -1:
                if exoifcutils.inbox(door_box, aNode[1], aNode[2]):
                    DoorList.append(aNode)
    for Node in DoorList:
        g_door_connections.append([Node[0], DoorID])


def get_mid_elevation_stair_flight(stair_flight):
    Elevaction = stair_flight['Elevation']
    
    if 'PathData' in stair_flight:
        centre_data = int(len(stair_flight['PathData'])/2)
        Elevaction = stair_flight['PathData'][centre_data][2]  # z value

    return Elevaction


def order_stair_data(ifc_stair, space_list, stair_flights_list, landings_list):
    """
    Order the parts of the staircase from bottom floor to top
    """
    level_data = []
    if 'space_connecting' not in ifc_stair:
        return None
    
    if 'elemIDs' not in ifc_stair:
        return None
    
    for spaceId in ifc_stair['space_connecting']:
        spaceIndex = exoifcutils.SpaceDefined(spaceId, space_list)
        if spaceIndex > -1:
            space = space_list[spaceIndex]
            level_data.append([space['Elevation'], "IfcSpace", spaceIndex])

    if len(level_data) != 2:
        return None 
    
    for element in ifc_stair['elemIDs']:
        if element[0] == 'IfcStairFlight':
            flightIndex = exoifcutils.find_stairflight(element[1], stair_flights_list)
            if flightIndex > -1:
                stair_flight = stair_flights_list[flightIndex]
                elevation = get_mid_elevation_stair_flight(stair_flight)
                level_data.append([elevation, element[0], flightIndex])
        elif element[0] == 'IfcSlab':
            landingIndex = exoifcutils.find_landings(element[1], landings_list)
            if landingIndex > -1:
                landing = landings_list[landingIndex]
                level_data.append([landing['Elevation'], element[0], landingIndex])
        else:
            print("Unknown element in stair conectivity ", element)
        level_data = sorted(level_data)
    return level_data


def connect_virtual_spaces(space_list, ifc_file, settings):
    print("Adding Virtual elements================================")
    for space in space_list:
        if 'elemIDs' in space:
            for elems in space['elemIDs']:
                if elems[0] == 'IfcRelSpaceBoundaryV':  # note the V
                    print('IfcRelSpaceBoundaryV ', space['Reference'])
                    ifc_relspace_boundary = ifc_file.by_id(elems[1])
                    if ifc_relspace_boundary:
                        print(ifc_relspace_boundary)
                        ifc_relshape = ifc_relspace_boundary.ConnectionGeometry.SurfaceOnRelatingElement
                        if ifc_relshape:
                            if ifc_relshape.is_a('IfcCurveBoundedPlane') and ifc_relshape.InnerBoundaries is None:
                                ifc_relshape.InnerBoundaries = ()
                            entity_shape_geom = ifcopenshell.geom.create_shape(settings, ifc_relshape)
                            if entity_shape_geom:
                                print("Adding Virtual element ", ifc_relshape)
                            else:
                                print("No boundary")


'''
def connect_openings(ifc_file,settings):
    ifc_openings = ifc_file.by_type("IfcOpeningElement")
    for ifc_opening in ifc_openings:
        owner = ifc_opening.VoidsElements[0].RelatingBuildingElement
        if (not owner.is_a('IfcSlab')):
            a_story = ifcopenshell.util.element.get_container(owner)
            element_info = a_story.get_info()
            Elevation = exogetifcparam.get_property(element_info,'Elevation')
            plot_item =(Elevation>g_floorLevelMin and Elevation<g_floorLevelMax)
            if plot_item:
                print(ifc_opening)
                ifc_relshape = ifc_opening #.IsDefinedBy# ifc_opening.ConnectionGeometry.SurfaceOnRelatingElement
                if ifc_relshape:
                    if ifc_relshape.is_a('IfcCurveBoundedPlane') and ifc_relshape.InnerBoundaries is None:
                        ifc_relshape.InnerBoundaries = ()
                    entity_shape_geom = ifcopenshell.geom.create_shape(settings, ifc_relshape)
                    if entity_shape_geom:
                        print("Adding opening element ",ifc_relshape)
                    else:
                        print("No opening boundary")
'''


def connect_rooms(Space1, NodeList1, Space2, NodeList2, door_box):
    print("connect_room: DoorBox ", door_box, "Spaces ", Space1['Name'], Space2['Name'])
    DoorList1 = []
    for aRow in NodeList1:
        for aNode in aRow:
            if aNode[0] > -1:
                if exoifcutils.inbox(door_box, aNode[1], aNode[2]):
                    DoorList1.append(aNode)

    DoorList2 = []
    for aRow in NodeList2:
        for aNode in aRow:
            if aNode[0] > -1:
                if exoifcutils.inbox(door_box, aNode[1], aNode[2]):
                    DoorList2.append(aNode)

    print("Connect Nodes ", DoorList1, DoorList2)

    iConnects = min(len(DoorList1), len(DoorList2))
    
    for iConnect in range(iConnects):
        Node1 = DoorList1[iConnect][0]
        Node2 = DoorList2[iConnect][0]
        g_door_connections.append([Node1, Node2])


def connect_rooms_doors(door_list, space_list, ifc_file, settings):
    global g_space_node_list
    ifc_doors = ifc_file.by_type("IfcDoor")
    for ifc_door in ifc_doors:
        element_info = ifc_door.get_info()
        space_id = exogetifcparam.get_property(element_info, 'GlobalId')
        doorIndex = exoifcutils.DoorDefined(space_id, door_list)
        if doorIndex > -1:
            door_data = door_list[doorIndex]
            spaceIndexList = []
            for SpaceID in door_data['Spaces']:
                spaceIndex = exoifcutils.SpaceDefined(SpaceID, space_list)
                if spaceIndex > -1:
                    spaceIndexList.append(spaceIndex)
                                          
            if len(spaceIndexList) == 2:
                entity_shape = ifcopenshell.geom.create_shape(settings, ifc_door)
                door_box = calc_door_box(entity_shape.geometry)
                #door_region = calc_door_region(entity_shape.geometry)
                connect_rooms(space_list[spaceIndexList[0]], g_space_node_list[spaceIndexList[0]], space_list[spaceIndexList[1]], g_space_node_list[spaceIndexList[1]], door_box)
            else:
                if not door_data['IsExternal']:
                    door_data['IsExternal'] = True

      
def connect_external_doors(door_list, space_list, ifc_file,settings):
    global g_DoorID;
    for door_data in door_list:
        if door_data['IsExternal']:
            if len(door_data['Spaces']) == 1:
                print("External door ", door_data)
                ifc_door = ifc_file.by_id(door_data['GlobalId'])
                entity_shape = ifcopenshell.geom.create_shape(settings, ifc_door)
                door_box = calc_door_box(entity_shape.geometry)
                door_region = calc_door_region(entity_shape.geometry)
                door_data['NodeID']=g_DoorID
                g_DoorID += 1
                aSpaceIndex = exoifcutils.SpaceDefined(door_data['Spaces'][0], space_list)
                if aSpaceIndex > -1:
                    connect_door(space_list[aSpaceIndex], g_space_node_list[aSpaceIndex],door_box, door_data['NodeID'])
                    door_data['Location']=exoifcutils.midbox(door_box)
                else:
                    print("Can't find space external door ", door_data)
                print("")
            else:
                print("Odd external door ", door_data)
                
                
def connect_stairflight(Space1, NodeList1, stair_box):
    NodeList = []
    print("Area is :", stair_box)
    for aRow in NodeList1:
        for aNode in aRow:
            if aNode[0] > -1:
                if exoifcutils.PointInPolygon(aNode[1], aNode[2], stair_box) == 1:
                    NodeList.append(aNode)

    return NodeList


def plot_transit_boxes(transit_object, plotter_class):
    if 'TopBox' in transit_object:
        if transit_object['TopBox'] is not None:
            z = transit_object['Elevation']+transit_object['Height']-0.24
            if plotter_class is not None:
                plotter_class.add_tri_delaunay(transit_object['TopBox'], z)
    if 'LowBox' in transit_object:
        if transit_object['LowBox'] is not None:
            z = transit_object['Elevation']+0.24
            if plotter_class is not None:
                plotter_class.add_tri_delaunay(transit_object['LowBox'], z)

            
def assign_stairflight_id_direction(stairflight):
    """
    assign node id to stairflight and calc direction and distances
    """
    global g_NodeID
    stairflight['nodeid'] = g_NodeID
    g_NodeID+=1
    direction = calc_flightstair_direction(stairflight)
    print("Stairflight Direction Data ", stairflight['Direction'], direction)
    stairflight['Direction'] = direction
    horizontal_length, travel_length = calc_flightstair_physical_param(stairflight)
    print("Stairflight width, horizontal_length, travel_length  ", stairflight['Width'], horizontal_length, travel_length)
    print("Stairflight HorizontalLength", horizontal_length, stairflight['TotalRun'])


def connect_stairs(ifc_stair, space_list, landings_list, stair_flights_list, level_data, transit_node_offset):
    global g_stairflight_connections
    element1 = None
    for element2 in level_data:
        if element1 is not None:
            stairflightIndex = -1
            landingIndex = -1
            spaceIndex = -1
            top = False
            if element1[1] == 'IfcStairFlight':
                top = True
                stairflightIndex = element1[2]
                if element2[1] == 'IfcSlab':
                    landingIndex = element2[2]
                elif element2[1] == 'IfcSpace':
                    spaceIndex = element2[2]
            elif element2[1] == 'IfcStairFlight':
                stairflightIndex = element2[2]
                if element1[1] == 'IfcSlab':
                    landingIndex = element1[2]
                elif element1[1] == 'IfcSpace':
                    spaceIndex = element1[2]
            else:
                print("Badly defined staircase ", ifc_stair['Name'])

            stairflight = None
            if stairflightIndex > -1:
                stairflight = stair_flights_list[stairflightIndex]
                if 'LowBox' not in stairflight:
                    LowBox, TopBox = calc_flightstair_box(stairflight)
                    stairflight['LowBox'] = LowBox
                    stairflight['TopBox'] = TopBox
                    assign_stairflight_id_direction(stairflight)

            NodeList = None   
            if stairflightIndex > -1 and spaceIndex > -1:
                if top:
                    NodeList = connect_stairflight(space_list[spaceIndex], g_space_node_list[spaceIndex], stairflight['LowBox'])
                    stairflight['ConxTop'] = NodeList
                else:
                    NodeList = connect_stairflight(space_list[spaceIndex], g_space_node_list[spaceIndex], stairflight['TopBox'])
                    stairflight['ConxLower'] = NodeList
                print("Space Node List is ", NodeList)
            elif stairflightIndex is not None and landingIndex is not None:
                if top:
                    print("Stair connect top")
                else:
                    print("Stair connect bot")
            else:
                print("Badly defined staircase order ", ifc_stair['Name'])

            if NodeList is not None:
                for aNode in NodeList:
                    g_stairflight_connections.append([aNode[0], stairflight['nodeid']])
                    pass
                
        element1 = element2


def connect_escalator(escalator, landings_list, space_list):
    UpperPoints = escalator['GeomData']['UpperPoints']
    LowerPoints = escalator['GeomData']['LowerPoints']

    TopBox = MultiPoint(UpperPoints).minimum_rotated_rectangle
    LowBox = MultiPoint(LowerPoints).minimum_rotated_rectangle

    escalator['TopBox'] = list(TopBox.exterior.coords)
    escalator['LowBox'] = list(LowBox.exterior.coords)

    print("Escalator BOX DATA")
    print("Upper:", escalator['TopBox'])
    print("Lower:", escalator['LowBox'])

    # TO DO use boxes to connect up to spaces/landings
    
    
def exodus_ouput(OMA_Class, ifc_file, settings, plotter_class, MTAFile):
    global g_space_node_list
    global g_lines
    global g_DoorID
    global g_NodeID
    global g_StairID
    global g_subp_list
    global g_door_connections
    global g_stairfligh_connections
    
    add_floor_nodes(OMA_Class.m_space_list, ifc_file, settings)
    add_landing_nodes(OMA_Class.m_landings_list, ifc_file, settings)
    
    connect_rooms_doors(OMA_Class.m_door_list, OMA_Class.m_space_list, ifc_file, settings)
    connect_external_doors(OMA_Class.m_door_list, OMA_Class.m_space_list, ifc_file, settings)

    #connect_openings(ifc_file,settings)
    #connect_virtual_spaces(ifc_file,settings)
    
    transit_node_offset = g_NodeID;
    for stair in OMA_Class.m_stair_list:
        level_data = order_stair_data(stair, OMA_Class.m_space_list, OMA_Class.m_stair_flights_list, OMA_Class.m_landings_list)
        if level_data is not None:
            connect_stairs(stair, OMA_Class.m_space_list, OMA_Class.m_landings_list, OMA_Class.m_stair_flights_list, level_data, transit_node_offset)

    for escalator in OMA_Class.m_escalator_list:
        connect_escalator(escalator, OMA_Class.m_landings_list, OMA_Class.m_space_list)

    offsetX = 2.0
    offsetY = 2.0
    exooutput.AdjustLineLoc(g_rangeX, g_rangeY, offsetX, offsetY, g_lines)
    exooutput.AdjustNodeLoc(g_rangeX, g_rangeY, offsetX, offsetY, g_space_node_list)
    exooutput.AdjustDoorsLoc(g_rangeX, g_rangeY, offsetX, offsetY, OMA_Class.m_door_list)
    # AdjustTransitNodeLoc function also updates direction
    # however staircase direction is recalculated by exodusmesh.calc_flightstair_direction later
    exooutput.AdjustTransitNodeLoc(g_rangeX, g_rangeY, offsetX, offsetY, OMA_Class.m_stair_flights_list)
    exooutput.AdjustTransitNodeLoc(g_rangeX, g_rangeY, offsetX, offsetY, OMA_Class.m_escalator_list)
    exooutput.AdjustTransitNodeLoc(g_rangeX, g_rangeY, offsetX, offsetY, OMA_Class.m_elevator_list)
    exooutput.AdjustTransitNodeLoc(g_rangeX, g_rangeY, offsetX, offsetY, OMA_Class.m_movingwalkway_list)

    print("g_lines", len(g_lines))
    print("g_space_node_list", len(g_space_node_list))
    print("g_subp_list", len(g_subp_list))
    print("g_door_connections", len(g_door_connections))
 
    # Save EXODUS data
    exodusmtafile = open(MTAFile, "w")
    
    exooutput.OutputMTAFileHeader(exodusmtafile)
    floorIndex = 0
    for floor in OMA_Class.m_floor_list:
        exooutput.OutputWindowFloorData(floorIndex, floor['Name'], exodusmtafile)
        floorIndex += 1
    
    external_doors = 0
    for door in OMA_Class.m_door_list:
        if door['IsExternal']:
            if 'Location' in door:
                exooutput.output_external_door(door, door['FloorIndex'], g_door_connections, exodusmtafile)
                external_doors += 1

    print("External doors found:", external_doors)

    floorIndex = 0
    for floor in OMA_Class.m_floor_list:
        for space in OMA_Class.m_space_list:
            if space['FloorIndex'] == floorIndex:
                space_nodes = g_space_node_list[space['NodeLocation']]
                exooutput.output_room_nodes(space_nodes, floorIndex, g_door_connections, g_stairflight_connections, exodusmtafile)

        for landing in OMA_Class.m_landings_list:
            if landing['FloorIndex'] == floorIndex:
                landing_nodes = g_space_node_list[landing['NodeLocation']]
                height = landing['Height']-floor['Elevation']
                exooutput.output_landing_nodes(landing_nodes, height, floorIndex, exodusmtafile)
        
        floorIndex += 1
        
    for stairflight in OMA_Class.m_stair_flights_list:
        if 'nodedata' in stairflight:
            stair_node = stairflight['nodedata']
            #Elevation = stair_node[4]
            if stair_node[0]<transit_node_offset:
                assign_stairflight_id_direction(stairflight)
            stair_node[0] = stairflight['nodeid']
            
            aWinID = stairflight['FloorIndex']
            node_data = [stair_node[0], stair_node[1], stair_node[2], stair_node[3], stairflight['Name']]
            exooutput.OutputNodeLoc(node_data, 80, aWinID, exodusmtafile)
            exooutput.OutputTransitNodeConnections(stairflight, exodusmtafile)

    if OMA_Class.m_escalator_list is not None:
        for escalator in OMA_Class.m_escalator_list:
            if 'nodedata' in escalator:
                
                escalator_node_data = escalator['nodedata']
                if (escalator_node_data[0] < transit_node_offset):
                    escalator['nodeid'] = g_NodeID
                    g_NodeID += 1
                escalator_node_data[0] = escalator['nodeid']
                
                #Elevation = escalator_node_data[4]
                aWinID = escalator['FloorIndex']
                
                node_data = [escalator_node_data[0], escalator_node_data[1], escalator_node_data[2], escalator_node_data[3], escalator['Name']]
                exooutput.OutputNodeLoc(node_data, 80, aWinID, exodusmtafile)

    if OMA_Class.m_movingwalkway_list is not None:
        for movingwalkway in OMA_Class.m_movingwalkway_list:
            if 'nodedata' in movingwalkway:
                movingwalkway_node_data = movingwalkway['nodedata']
                if movingwalkway_node_data[0] < transit_node_offset:
                    movingwalkway['nodeid'] = g_NodeID
                    g_NodeID += 1
                movingwalkway_node_data[0] = movingwalkway['nodeid']
                
                #Elevation = movingwalkway_node_data[4]
                aWinID = movingwalkway['FloorIndex']
                
                node_data = [movingwalkway_node_data[0], movingwalkway_node_data[1], movingwalkway_node_data[2], movingwalkway_node_data[3], movingwalkway['Name']]
                exooutput.OutputNodeLoc(node_data, 80, aWinID, exodusmtafile)

    lift_count = 0
    if OMA_Class.m_elevator_list is not None:
        for elevator in OMA_Class.m_elevator_list:
            if 'nodedata' in elevator:
                elevator_node_data = elevator['nodedata']
                if elevator_node_data[0] < transit_node_offset:
                    elevator['nodeid'] = g_NodeID
                    g_NodeID += 1
                elevator_node_data[0] = elevator['nodeid']
                Count = 0
                LiftNodeList = []
                for floor_elevation in elevator['FloorListHeights']:
                    # need to generate a transit node per a floor - lift shaft opening
                    aWinID = exoifcutils.find_floor(floor_elevation, OMA_Class.m_floor_list)
                    if Count == 0:
                        NodeID = elevator_node_data[0]
                    else:
                        NodeID += 1
                        g_NodeID += 1
                    if aWinID > -1:
                        node_data = [NodeID, elevator_node_data[1], elevator_node_data[2], elevator_node_data[3], elevator['Name']]
                        exooutput.OutputNodeLoc(node_data, 80, aWinID, exodusmtafile)
                        LiftNodeList.append([NodeID, aWinID])
                        Count += 1
                if len(LiftNodeList) > 0:
                    lift_count += 1
                    elevator['LiftNodeList'] = LiftNodeList

    TransitID = 1
    for stairflight in OMA_Class.m_stair_flights_list:
        aWinID = stairflight['FloorIndex']
        print("Saving stair flight ", stairflight['Name'])
        plot_transit_boxes(stairflight, plotter_class)
        exooutput.OutputTransitNodeData(stairflight, TransitID, aWinID, 0, exodusmtafile)
        TransitID += 1

    if OMA_Class.m_escalator_list is not None:
        for escalator in OMA_Class.m_escalator_list:
            aWinID = escalator['FloorIndex']
            exooutput.OutputTransitNodeData(escalator, TransitID, aWinID, 2, exodusmtafile)
            TransitID += 1

    if OMA_Class.m_movingwalkway_list is not None:
        for movingwalkway in OMA_Class.m_movingwalkway_list:
            aWinID = movingwalkway['FloorIndex']
            exooutput.OutputTransitNodeData(movingwalkway, TransitID, aWinID, 3, exodusmtafile)
            TransitID += 1

    if OMA_Class.m_elevator_list is not None:
        for elevator in OMA_Class.m_elevator_list:
            if 'LiftNodeList' in elevator:
                for NodeData in elevator['LiftNodeList']:
                    elevator['nodeid'] = NodeData[0]
                    aWinID = NodeData[1]
                    exooutput.OutputTransitNodeData(elevator, TransitID, aWinID, 8, exodusmtafile)
                    TransitID += 1

    if lift_count > 0:
        LiftBankID = 0
        print("LiftsCount:", lift_count, file=exodusmtafile)
        for elevator in OMA_Class.m_elevator_list:
            if 'LiftNodeList' in elevator:
                LiftNodeList = elevator['LiftNodeList']
                exooutput.OutputLiftData(elevator, LiftBankID, LiftNodeList, exodusmtafile)
                LiftBankID += 1
        

    exodusmtafile.close()

    # Save EXODUS EGX
    pre, ext = os.path.splitext(MTAFile)
    EGXFile = pre+".egm"
    #exoifcutils.remove_virtual_lines(g_lines, OMA_Class.m_space_list, g_subp_list)

    exodusegxfile = open(EGXFile, "w")
    exooutput.OutputLineFileHeader(len(OMA_Class.m_floor_list), exodusegxfile)  # just 1 floor

    floorIndex = 0
    for floor in OMA_Class.m_floor_list:
        win_line_count = 0
        for line in g_lines:
            if line[5] == floorIndex:
                win_line_count += 1
        exooutput.OutputLineSectionStart(win_line_count, exodusegxfile)
        line_height = 0.0
        for line in g_lines:
            if line[5] == floorIndex:
                exooutput.OutputLineData(line, line_height, exodusegxfile)
        exooutput.OutputLineSectionEnd(exodusegxfile)
        floorIndex += 1

    exooutput.OutputLineTail(exodusegxfile)
    exodusegxfile.close()

    


