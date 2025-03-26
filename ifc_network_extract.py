# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence

packages used ifcopenshell, scipy, matplotlib,pyvis
"""

import sys
import os
import ifcopenshell
from ifcopenshell import geom
import ifcopenshell.util.element as Element
import numpy as np
import math
from operator import itemgetter

import exogetifcparam
import exoifcutils

import enzoutput
import graphoutput
import exodusmesh
import cfastoutput
import ifc_config

from plot_functions import PlotterClass
from oma_class import OMAClass


g_StairNodeID = 0  # used for numbering doors local IDs
g_OMA_Class = OMAClass()
g_plotter_class = None

def plot_geom(shape):
    global g_plotter_class
    faces = shape.faces  # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    VertexCount = len(verts)
    if VertexCount < 3:
        return None
    
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

    if g_plotter_class is not None:
        g_plotter_class.add_faces(xp, yp, zp, faces)


def sign_locations_xy(shape):
    """
    calculates location of sign
    assumes vertices 1 and 2 space the signs width
    """
    startloc = [0.0,0.0]
    endloc = [0.0,0.0]
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]

    VertexCount = len(verts)
    if VertexCount>13:
        startloc = [verts[9],verts[10]]
        endloc   = [verts[12],verts[13]]
        
    return startloc,endloc


def level_geom(shape, LevelElevation):
    """
    calculates Mid point (lx,ly) Height above floor level (NdeHeight )
    Direction and width of a object connecting between levels
    based on its 3D shape (i.e. staircase, ramp)
    """
    global g_plotter_class
    iPoint = 0
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    
    VertexCount = len(verts)
    xmin = verts[0]
    xmax = xmin
    ymin = verts[1]
    ymax = ymin
    zmin = verts[2]
    
    while iPoint < VertexCount:
        xmin = min(xmin, verts[iPoint])
        xmax = max(xmax, verts[iPoint])
        iPoint += 1
        ymin = min(ymin, verts[iPoint])
        ymax = max(ymax, verts[iPoint])
        iPoint += 1
        zmin = min(zmin, verts[iPoint])
        iPoint += 1
   
    stair_box = []

    width = -1.0
    flatLength = -1.0
    Direction = 0
    # Temp code to calc direction
    # Node for staircases direction is later recalculated using path data based on steps
    # see exodusmesh.calc_flightstair_direction
    if ymax-ymin > xmax-xmin:
        flatLength = ymax-ymin
        width = xmax-xmin
        Direction = 0
    else:
        width = ymax-ymin
        flatLength = xmax-xmin
        Direction = 270
    # end of temp code
    stair_box.append([xmin, ymax])
    stair_box.append([xmax, ymax])
    stair_box.append([xmax, ymin])
    stair_box.append([xmin, ymin])

    if LevelElevation!=None:
        NdeHeight = zmin - LevelElevation
    else:
        NdeHeight = zmin
    lx = (xmax+xmin)/2.0
    ly = (ymax+ymin)/2.0

    tri_idx = [[0, 1, 2], [2, 3, 0]]
    xp = [stair_box[0][0], stair_box[1][0], stair_box[2][0], stair_box[3][0]]
    yp = [stair_box[0][1], stair_box[1][1], stair_box[2][1], stair_box[3][1]]
    zp = [zmin, zmin, zmin, zmin]
    if g_plotter_class is not None:
        g_plotter_class.add_trisurf(xp, yp, zp, tri_idx)
    return lx, ly, NdeHeight, Direction, width, flatLength


def get_min_max_values(verts):
    VertexCount = len(verts)
    xmin = verts[0]
    xmax = xmin
    ymin = verts[1]
    ymax = ymin
    zmin = verts[2]
    zmax = zmin

    iPoint = 3
    while iPoint < VertexCount:
        xmin = min(xmin, verts[iPoint])
        xmax = max(xmax, verts[iPoint])
        iPoint += 1
        ymin = min(ymin, verts[iPoint])
        ymax = max(ymax, verts[iPoint])
        iPoint += 1
        zmin = min(zmin, verts[iPoint])
        zmax = max(zmax, verts[iPoint])
        iPoint += 1

    return xmin, xmax, ymin, ymax, zmin, zmax


def level_flight_geom(shape, LevelElevation, Tolerance):
    verts = shape.verts # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]

    VertexCount = len(verts)
    if VertexCount < 3:
        return None

    xmin, xmax, ymin, ymax, zmin, zmax = get_min_max_values(verts)
    
    upper_points = []
    lower_points = []
    iPoint = 0
    while iPoint < VertexCount:
        x = verts[iPoint]
        iPoint += 1
        y = verts[iPoint]
        iPoint += 1
        z = verts[iPoint]
        iPoint += 1
        if abs(z-zmax) < Tolerance:
            upper_points.append([x, y, z])
        elif abs(z-zmin) < Tolerance:
            lower_points.append([x, y, z])
   
    mid_x = (xmax+xmin)/2.0
    mid_y = (ymax+ymin)/2.0

    Level_data = {"MidX": mid_x, "MidY": mid_y}
    Level_data["UpperPoints"] = upper_points
    Level_data["LowerPoints"] = lower_points
    return Level_data


def add_step_level(levels, tri_idx, height, Tol):
    """
    Adds to the data entry for a list of points at a given level interval
    Used to store information about step heights in staircases
    """
    for entry in range(len(levels)):
        if abs(levels[entry][0]-height) < Tol:
            # include in this entry
            levels[entry][1].append(tri_idx)
            return True
        elif (height-levels[entry][0]) > 0:
            # insert here
            levels.insert(entry, (height, [tri_idx]))
            return True
    
    levels.append((height, [tri_idx]))
    return True


def cross2d(x, y):
    return x[..., 0] * y[..., 1] - x[..., 1] * y[..., 0]


def extract_flight_geom(shape, Tol):
    '''
    Extract path and step data for stairs
    :return: number of levels (steps) found, step level data and mid point path data
    '''
    global g_plotter_class
    faces = shape.faces  # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    VertexCount = len(verts)
    if VertexCount < 3:
        return None

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

    levels = []
    iPoint = 0
    FaceCount = len(faces)
    while iPoint < FaceCount:
        ip1 = int(faces[iPoint])
        iPoint += 1
        ip2 = int(faces[iPoint])
        iPoint += 1
        ip3 = int(faces[iPoint])
        iPoint += 1
        if abs(zp[ip1]-zp[ip2]) < Tol and abs(zp[ip1]-zp[ip3]) < Tol and abs(zp[ip2]-zp[ip3]) < Tol:
            add_step_level(levels, (ip1, ip2, ip3), zp[ip1], Tol)

    path_data = []
    for entry in levels:
        if g_plotter_class is not None:
            g_plotter_class.add_faces(xp, yp, zp, entry[1])
        if len(entry[1]) > 0:
            ax = 0
            ay = 0
            az = 0
            pc = len(entry[1])
            for ip in entry[1]:
                for i in range(3):
                    ax += xp[ip[i]]
                    ay += yp[ip[i]]
                    az += zp[ip[i]]
            ax = ax / (pc*3)
            ay = ay / (pc*3)
            az = az / (pc*3)
            path_data.append((ax, ay, az))

    width = 0.0
    if len(path_data) > 1:
        p1 = np.array([path_data[0][0], path_data[0][1]])
        p2 = np.array([path_data[1][0], path_data[1][1]])
        aStepData = levels[1][1]
        max_d = -1.0
        min_d = 1.0
        # calculate stair width based on maximum distance of step points from path line
        for ip in aStepData:
            for i in range(3):
                ax = xp[ip[i]]
                ay = yp[ip[i]]
                p3 = np.array([ax, ay])
                # calculate maximum distance from path line (p1,p2)
                # -ve or +ve depends on which side of the line the point is on
                d = cross2d(p2-p1, p3-p1)/np.linalg.norm(p2-p1)
                max_d = max(d, max_d)
                min_d = min(d, min_d)
        # width is the max distance minus the max -ve (i.e. min) distance
        width = max_d-min_d
        print("Step Max Data is ", min_d, max_d, width)

    return len(levels), levels, path_data, width


def get_element_top_bottom(shape):
    verts = shape.verts # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    
    iPoint = 0
    zp=[]
    VertexCount = len(verts)
    while iPoint < VertexCount:
        iPoint += 2
        zp.append(verts[iPoint])
        iPoint += 1

    lmaxZ = max(zp)
    lminZ = min(zp)
    return lmaxZ, lminZ


def get_floor_area(shape, Elev, aTol):
    faces = shape.faces  # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts  # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    global g_plotter_class
    
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

    if g_plotter_class is not None:
        g_plotter_class.update_ranges(lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ)
    
    # extract face area on floor
    iPoint = 0
    tri_idx = []
    FaceCount = len(faces)
    area = 0.0
    while iPoint < FaceCount:
        ip1 = int(faces[iPoint])
        iPoint += 1
        ip2 = int(faces[iPoint])
        iPoint += 1
        ip3 = int(faces[iPoint])
        iPoint += 1
        if zp[ip3] > Elev-aTol and zp[ip3] < Elev+aTol and zp[ip2] > Elev-aTol and zp[ip2] < Elev+aTol and zp[ip1] > Elev-aTol and zp[ip1] < Elev+aTol:
            area += abs(0.5*(xp[ip1]*(yp[ip2] - yp[ip3]) + xp[ip2]*(yp[ip3] - yp[ip1]) + xp[ip3]*(yp[ip1] - yp[ip2])))
            tri_idx.append((ip1, ip2, ip3))

    if g_plotter_class is not None:
        g_plotter_class.add_faces(xp, yp, zp, tri_idx)
    
    return area


def generate_boundary(shape, Elev, aTol):
    """
      Generate a list of boundary points from a shape at given elevation
      Note entries in edge_list_pt come in groups of 2, so 0,1 is an edge, 2,3 is an edge and so on
    """
    global g_plotter_class
    faces = shape.faces # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    
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

    if g_plotter_class is not None:
        g_plotter_class.update_ranges(lminX, lmaxX, lminY, lmaxY, lminZ, lmaxZ)
    
    # extract face area on floor
    iPoint = 0
    tri_idx = []
    FaceCount = len(faces)
    while iPoint < FaceCount:
        ip1 = int(faces[iPoint])
        iPoint += 1
        ip2 = int(faces[iPoint])
        iPoint += 1
        ip3 = int(faces[iPoint])
        iPoint += 1
        if zp[ip3] > Elev-aTol and zp[ip3] < Elev+aTol and zp[ip2] > Elev-aTol and zp[ip2] < Elev+aTol and zp[ip1] > Elev-aTol and zp[ip1] < Elev+aTol:
            tri_idx.append((ip1, ip2, ip3))
    
    edge_list = find_boundary(tri_idx)
    print("boundary edge list ", edge_list)
    edge_list_pt = None
    if len(edge_list) > 0:
        edge_list_pt = []
        for edge in edge_list:
            edge_list_pt.append([xp[edge[0]], yp[edge[0]], zp[edge[0]]])
            edge_list_pt.append([xp[edge[1]], yp[edge[1]], zp[edge[1]]])
    else:
        return "Invalid"
    
    return edge_list_pt


def find_boundary(tri_idx):
    """
      Generate a list edges from a set of triangles 
    """
    edge_list = []
    # Convert to a list of edges
    for itri in tri_idx:
        for i1 in range(3):
            i2 = (i1+1) % 3
            edge_list.append([itri[i1], itri[i2]])

    # Find edges which are defined twice, i.e. defined in two triangles
    edge_count = len(edge_list)
    removelist = []
    for ib in range(edge_count-1):
        if not(ib in removelist):
            for jb in range(ib+1, edge_count):
                if not(jb in removelist):
                    if ((edge_list[ib][0] == edge_list[jb][0] and edge_list[ib][1] == edge_list[jb][1]) or
                        (edge_list[ib][0] == edge_list[jb][1] and edge_list[ib][1] == edge_list[jb][0])):
                        removelist.append(ib)
                        removelist.append(jb)
                        break
						
    # The edges which are not defined twice must be on the boundary
    boundary_list = []
    for ib in range(edge_count):
        if not(ib in removelist):
            boundary_list.append(edge_list[ib])

    return boundary_list


def get_list_of_stair(stairflight_list):
    """
    :returns: a list if ifc_stairs from a list os ifc_stairflight, positions m_stair_list
    :param stairflight_list list of induces to m_stair_flights_list
    """
    global g_OMA_Class

    stair_list = []
    
    for stairflight_index in stairflight_list:
        astairflight = g_OMA_Class.m_stair_flights_list[stairflight_index]
        parent_GlobalID = astairflight['parent']
        stair_index = exoifcutils.StairsDefined(parent_GlobalID, g_OMA_Class.m_stair_list)
        if stair_index>-1:
            if stair_index not in stair_list:
                stair_list.append(stair_index)

    return stair_list


def find_stairflights_in_shape(shape, Elev):
    """
    find list of stairflights in on a given IfcSpace as specified by Elevation
    
    :return: two lists of locations in m_stair_flights_list, above and below
    """
    global g_OMA_Class
    
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
        
    stairflights_up = find_stairflights_in(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, Elev)
    stairflights_down = None
    below_Elev = exoifcutils.find_floor_below_elevation(Elev, g_OMA_Class.m_floor_list)
    if below_Elev is not None:
        stairflights_down = find_stairflights_in(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, below_Elev)
    return stairflights_down, stairflights_up


def find_escalators_in_shape(shape, Elev):
    """
    find list of escalators in or on a given IfcSpace as specified by Elevation
    
    :return: two lists of locations in m_escalator_list, above and beloe
    """
    global g_OMA_Class
    
    faces = shape.faces # Indices of vertices per triangle face e.g. [f1v1, f1v2, f1v3, f2v1, f2v2, f2v3,
    verts = shape.verts # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    
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
    #lminZ = min(zp)
    #lmaxZ = max(zp)
    
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
        
    escalators_up = find_escalators_in(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, Elev)
    escalators_down = None
    below_Elev = exoifcutils.find_floor_below_elevation(Elev, g_OMA_Class.m_floor_list)
    if below_Elev is not None:
        escalators_down = find_escalators_in(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, below_Elev)
    return escalators_down, escalators_up


def find_stairflights_in(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, Elevation):
    """
    find list of stairflights in or on a given IfcSpace as specified by Elevation
    
    :return: list of locations in m_stair_flights_list
    """
    global g_OMA_Class
    stairflightlist = []
    iPos = 0
    for a_stairflight in g_OMA_Class.m_stair_flights_list:
        if a_stairflight['Elevation'] == Elevation:
            node = a_stairflight['nodedata']
            lx = node[1]
            ly = node[2]+0.5  # had to shift stairs because there's a hole in the floor area in the hotel model at that location
            if lx > lminX and lx < lmaxX:
                if ly > lminY and ly < lmaxY:
                    print("Checking stair at ", lx, ly, lminX, lmaxX, lminY, lmaxY)
                    AddPoint = False
                    for triface  in tri_idx:
                        if exoifcutils.TriangleOrientation(xp[triface[0]], yp[triface[0]], xp[triface[1]], yp[triface[1]], xp[triface[2]], yp[triface[2]]) == 1:
                            AddPoint = exoifcutils.intri(xp[triface[0]], yp[triface[0]], xp[triface[1]], yp[triface[1]], xp[triface[2]], yp[triface[2]], lx, ly)
                        else:
                            AddPoint = exoifcutils.intri(xp[triface[2]], yp[triface[2]], xp[triface[1]], yp[triface[1]], xp[triface[0]], yp[triface[0]], lx, ly)
                        if AddPoint:
                            #print("Adding stair at ",lx,ly,lminX,lmaxX,lminY,lmaxY)
                            break
                    if AddPoint:
                        stairflightlist.append(iPos)
        iPos+=1
    return stairflightlist


def find_escalators_in(xp, yp, zp, tri_idx, lminX, lmaxX, lminY, lmaxY, Elevation):
    """
    find list of escalators in or on a given IfcSpace as specified by Elevation
    
    :return: list of locations in m_escalator_list
    """
    global g_OMA_Class
    
    escalatorlist = []
    for iPos in range(len(g_OMA_Class.m_escalator_list)):
        an_escalator = g_OMA_Class.m_escalator_list[iPos]
        if an_escalator['Elevation'] == Elevation:
            lx = an_escalator['xpos']
            ly = an_escalator['ypos']
            if lx > lminX and lx < lmaxX:
                if ly > lminY and ly < lmaxY:
                    print("Checking stair at ", lx, ly, lminX, lmaxX, lminY, lmaxY)
                    AddPoint = False
                    for triface  in tri_idx:
                        if exoifcutils.TriangleOrientation(xp[triface[0]], yp[triface[0]], xp[triface[1]], yp[triface[1]], xp[triface[2]], yp[triface[2]]) == 1:
                            AddPoint = exoifcutils.intri(xp[triface[0]], yp[triface[0]], xp[triface[1]], yp[triface[1]], xp[triface[2]], yp[triface[2]], lx, ly)
                        else:
                            AddPoint = exoifcutils.intri(xp[triface[2]], yp[triface[2]], xp[triface[1]], yp[triface[1]], xp[triface[0]], yp[triface[0]], lx, ly)
                        if AddPoint:
                            #print("Adding stair at ",lx,ly,lminX,lmaxX,lminY,lmaxY)
                            break
                    if AddPoint:
                        escalatorlist.append(iPos)
    return escalatorlist


def get_space_data(ifc_file, ifc_space, Elevation, floor_longname, settings):
    """
    Loads and adds ifc_space data to m_space_list
    
    :param ifc_space the ifc space data
    :param Elevation the space elevation
    :param floor_longname long floor name
    :param settings - ifc file settings

    :return location in g_space_data
    """
    global g_OMA_Class
    
    NumNodes = 0
    space_data = {}
    exogetifcparam.get_basic_space_info(ifc_space, space_data)
    # space relationship to doors etc
    boundary_elem = exogetifcparam.get_space_boundary_elem(ifc_space)
    space_data['elemIDs'] = boundary_elem
    ifc_space_shape = ifcopenshell.geom.create_shape(settings, ifc_space)
    #plot_geom(ifc_space_shape.geometry)
    # add door relationships to this space
    
    for item in boundary_elem:
        if item[0] == "IfcDoor":
            doorIndex = exoifcutils.DoorDefined(item[1], g_OMA_Class.m_door_list)
            
            if doorIndex > -1:
                # test if the door is in the correct ifc_space:
                if True:
                    # confirm the door is near the space
                    ifc_door = ifc_file.by_id(g_OMA_Class.m_door_list[doorIndex]['GlobalId'])
                    door_shape = ifcopenshell.geom.create_shape(settings, ifc_door)
                    aPoint = exoifcutils.get_centre(door_shape.geometry)
                    aDist = exoifcutils.CentreDistOfNearestFaceToPointTriangle(ifc_space_shape.geometry, aPoint)
                    if aDist>0.2:
                        # Not on the edge of the IfcSpace
                        print("Warning: Looks like this door ",ifc_door.Name, " is associated with the wrong space",ifc_space.Name)
                        doorIndex = -1 # reject door
               
            if doorIndex > -1:
                g_OMA_Class.m_door_list[doorIndex]['Spaces'].append(space_data['GlobalId'])

    space_data['subspaces'] = []
    g_OMA_Class.m_space_list.append(space_data)
    space_loc = len(g_OMA_Class.m_space_list)-1

                            
    if ifc_space_shape:
        #print("Checking IFC Space",ifc_space)
        space_data['Elevation'] = Elevation
        space_data['Floor'] = floor_longname
        space_data['FloorIndex'] = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
        area = get_floor_area(ifc_space_shape.geometry, Elevation, 0.1)
        space_data['area'] = area
        boundary_points = generate_boundary(ifc_space_shape.geometry, Elevation, 0.1)
        space_data['boundarylist'] = boundary_points

        # find stairs associated with this space
        stair_list_down, stair_list_up = find_stairflights_in_shape(ifc_space_shape.geometry, Elevation)
        if stair_list_up:
            # locations in m_stair_flights_list
            space_data['stairflightsUp'] = stair_list_up
            space_data['stairsUp'] = get_list_of_stair(stair_list_up)
        if stair_list_down:
            # locations in m_stair_flights_list
            space_data['stairflightsDown'] = stair_list_down
            space_data['stairsDown'] = get_list_of_stair(stair_list_down)

        # find escalators associated with this space
        escalator_list_down, escalator_list_up = find_escalators_in_shape(ifc_space_shape.geometry, Elevation)
        if escalator_list_up:
            # locations in m_escalator_list
            space_data['escalatorsUp'] = escalator_list_up 
        if escalator_list_down:
            # locations m_escalator_list m_stair_flights_list
            space_data['escalatorsDown'] = escalator_list_down 
        
    return space_loc


def stairs_data(ifc_file, settings):
    global g_OMA_Class
    global g_StairNodeID
    
    #global g_plotter_class
    
    ifc_stair = ifc_file.by_type("IfcStair")
    for stair in ifc_stair:
        if stair.IsDecomposedBy:
            Elevation = exogetifcparam.get_storey_Elevation(stair)
            floor_longname = exogetifcparam.get_storey_LongName(stair)
            floor_index = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
            if floor_index < 0:
                floor = {'Elevation': Elevation}
                floor['Name'] = floor_longname
                g_OMA_Class.m_floor_list.append(floor)
            plot_item = True
            if plot_item:
                print("Adding Staircase", Elevation)
                stair_dict = {'id': g_StairNodeID}
                exogetifcparam.get_basic_stair_info(stair, stair_dict)
                stair_dict['elemIDs'] = []
                g_OMA_Class.m_stair_list.append(stair_dict)
                # Retrieve all the components of IfcStair
                stair_components = stair.IsDecomposedBy[0].RelatedObjects
                for stair_element in stair_components:
                    if 'IfcStairFlight' == stair_element.is_a():
                        print('IfcStairFlight')
                        # Retrieved riser_number, tread_length, stair_width
                        stair_flight_dict = {'id': g_StairNodeID}
                        
                        exogetifcparam.get_stair_flight_info(stair_element, stair_flight_dict, stair_dict)
                        
                        entity_shape = ifcopenshell.geom.create_shape(settings, stair_element)
                        lx, ly, NdeHeight, Direction, width, CalcTotalRun = level_geom(entity_shape.geometry, Elevation)
                        geom_data = level_flight_geom(entity_shape.geometry, Elevation, 0.05)
                        stair_flight_dict['GeomData'] = geom_data
                        stair_flight_dict['ExoAngle'] = CalcAngle(geom_data['UpperPoints'],geom_data['LowerPoints'])
                        AStepCount, AStepData, PathData, calcwidth = extract_flight_geom(entity_shape.geometry, stair_flight_dict['RiserHeight']/2.0)
                        if calcwidth > 0.0:
                            width = calcwidth
                        stair_flight_dict['AStepCount'] = AStepCount
                        stair_flight_dict['AStepData'] = AStepData
                        stair_flight_dict['PathData'] = PathData
                        #stair_flight_dict['CalcWidth'] = calcwidth
                        stair_flight_dict['Elevation'] = Elevation
                        stair_flight_dict['Floor'] = floor_longname
                        #stair_flight_dict['FloorIndex'] = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
                        stair_flight_dict['Direction'] = Direction
                        stair_flight_dict['Width'] = width 
                        stair_flight_dict['CalcTotalRun'] = CalcTotalRun
                        stair_flight_dict['Lanes'] = math.floor(width/0.7)
                        stair_flight_dict['parent'] = stair_dict['GlobalId']
                        stair_dict['elemIDs'].append(['IfcStairFlight', stair_flight_dict['GlobalId']])
                        
                        stair_flight_dict['nodedata'] = [g_StairNodeID, lx, ly, NdeHeight, Elevation]
                        stair_flight_dict['nodeid'] = g_StairNodeID
                        g_StairNodeID += 1
                        
                        g_OMA_Class.m_stair_flights_list.append(stair_flight_dict)

                    elif 'IfcSlab' == stair_element.is_a():
                        print('IfcSlab Stairs===================')
                        element_info = stair_element.get_info()
                        GlobalId = exogetifcparam.get_property(element_info, 'GlobalId')
                        ifcslab_dict = {'GlobalId': GlobalId}
                        ifcslab_dict['parent'] = stair_dict['GlobalId']
                        ifcslab_dict['Name'] = exogetifcparam.get_property(element_info, 'Name')
                        
                        predefined_type = exogetifcparam.get_property(element_info, 'PredefinedType')
                        print(f"  Predefined Type: {predefined_type}")
                        stair_dict['elemIDs'].append(['IfcSlab', GlobalId])
                        
                        entity_shape = ifcopenshell.geom.create_shape(settings, stair_element)
                        #if g_plotter_class is not None:
                        #    g_plotter_class.add_element_tri(entity_shape.geometry)
                        maxheight, minheight = get_element_top_bottom(entity_shape.geometry)
                        ifcslab_dict['Elevation'] = Elevation
                        area = get_floor_area(entity_shape.geometry, maxheight, 0.02)
                        ifcslab_dict['Height'] = maxheight
                        ifcslab_dict['area'] = area
                        boundary_points = generate_boundary(entity_shape.geometry, maxheight, 0.02)
                        ifcslab_dict['boundarylist'] = boundary_points
                        
                        g_OMA_Class.m_landings_list.append(ifcslab_dict)


def calc_escalator_level_run(escalator_dict):
    if 'OverallLength' in escalator_dict and 'RunLength' in escalator_dict:
        # OverallLength and RunLength are property values for the escalator
        OverallLength = escalator_dict['OverallLength']
        RunLength = escalator_dict['RunLength']
        if OverallLength is not None and RunLength < OverallLength:
            LevelRun = (OverallLength - RunLength)/2.0
            return LevelRun
                
    # else calc Level run
    if 'RunLength' in escalator_dict:
        # OverallLength is a property value for the escalator
        RunLength = escalator_dict['RunLength']
    else:
        # HorizontalLength calculated from the escalator geometry representation
        RunLength = escalator_dict['HorizontalLength']
    
    if 'RunHeight' in escalator_dict:
        Height = escalator_dict['RunHeight']
    else:
        Height = escalator_dict['Height']
    
    riser_height = 0.2
    tread_depth = 0.3
    if 'RiserHeight' in escalator_dict:
        riser_height = escalator_dict['RiserHeight']
        
    if 'TreadLength' in escalator_dict and escalator_dict['TreadLength'] is not None:
        tread_depth = escalator_dict['TreadLength']
        
    LevelRun = 1.2
    if riser_height > 0.0 and riser_height > 0.0:
        RisersUp = math.ceil(Height/riser_height)
        if RisersUp > 1:
            LevelRun =  (RunLength - (RisersUp-1) * tread_depth)/2.0
    return LevelRun


def get_escalator_data(transport_info, settings, ifc_transport, Elevation):
    global g_OMA_Class
    global g_StairNodeID
    
    escalator_dict = {'id': g_StairNodeID}
    g_StairNodeID += 1
    escalator_dict['Elevation'] = Elevation
    exogetifcparam.get_basic_escalator_info(escalator_dict, ifc_transport)

    entity_shape = ifcopenshell.geom.create_shape(settings, ifc_transport)
    lx, ly, NdeHeight, Direction, width, HorizontalLength = level_geom(entity_shape.geometry, Elevation)
    geom_data = level_flight_geom(entity_shape.geometry, Elevation, 0.05)
   
    escalator_dict['GeomData'] = geom_data
    #escalator_dict['Floor'] = floor_longname
    
    escalator_dict['Direction'] = Direction
    escalator_dict['Width'] = width
    escalator_dict['HorizontalLength'] = HorizontalLength
    escalator_dict['Lanes'] = math.floor(width/0.7)
    escalator_dict['xpos'] = lx
    escalator_dict['ypos'] = ly
    if 'RunHeight' in escalator_dict:
        escalator_dict['Height'] = escalator_dict['RunHeight']
    else:
        escalator_dict['Height'] = 3.0

    LevelRun = calc_escalator_level_run(escalator_dict)
    escalator_dict['LevelRun'] = LevelRun
    
    escalator_dict['nodedata'] = [g_StairNodeID, lx, ly, NdeHeight, Elevation]
    escalator_dict['nodeid'] = g_StairNodeID
    g_StairNodeID += 1
    
    g_OMA_Class.m_escalator_list.append(escalator_dict)


def get_elevator_data(transport_info, settings, ifc_transport, Elevation):
    global g_OMA_Class
    global g_StairNodeID
    
    print('ELEVATOR ', transport_info['Name'])
    elevator_dict = {'id': g_StairNodeID }
    elevator_dict['Elevation'] = Elevation
    exogetifcparam.get_basic_elevator_info(elevator_dict, ifc_transport)
    
    entity_shape = ifcopenshell.geom.create_shape(settings, ifc_transport)
    lx, ly, NdeHeight, Direction, width, HorizontalLength = level_geom(entity_shape.geometry, Elevation)
    boundary_points = generate_boundary(entity_shape.geometry, Elevation, 0.1)
    maxheight, minheight = get_element_top_bottom(entity_shape.geometry)
     
    #elevator_dict['Direction'] = Direction
    elevator_dict['Width'] = width
    elevator_dict['Lanes'] = math.floor(width/0.5)
    elevator_dict['Direction'] = 180.0
    elevator_dict['HorizontalLength'] = HorizontalLength
    elevator_dict['xpos'] = lx
    elevator_dict['ypos'] = ly
    elevator_dict['Height'] = 3.0
    elevator_dict['boundary_points'] = boundary_points
    elevator_dict['maxheight'] = maxheight
    elevator_dict['minheight'] = minheight

    elevator_dict['FloorListHeights'] = exoifcutils.get_subfloor_heights(minheight, maxheight, g_OMA_Class.m_floor_list)

    elevator_dict['nodedata'] = [g_StairNodeID, lx, ly, NdeHeight, Elevation]
    elevator_dict['nodeid'] = g_StairNodeID
    g_StairNodeID += 1
    
    print("Elevator data:",elevator_dict['Name'],elevator_dict['GlobalId'])
    g_OMA_Class.m_elevator_list.append(elevator_dict)
    print("Elevator data: ", elevator_dict)


def get_movingwalkway_data(transport_info, settings, ifc_transport, Elevation):
    global g_OMA_Class
    global g_StairNodeID
    global g_plotter_class
    
    print('MOVINGWALKWAY ', transport_info['Name'])
    movingwalkway_dict = {'id': g_StairNodeID}
    
    movingwalkway_dict['Elevation'] = Elevation
    exogetifcparam.get_basic_movingwalkway_info(movingwalkway_dict, ifc_transport)
    
    entity_shape = ifcopenshell.geom.create_shape(settings, ifc_transport)
    if g_plotter_class is not None:
        g_plotter_class.add_element_tri(entity_shape.geometry)
    lx, ly, NdeHeight, Direction, width, HorizontalLength = level_geom(entity_shape.geometry, Elevation)
    #movingwalkway_dict['Floor'] = floor_longname
    #movingwalkway_dict['FloorIndex'] = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
    movingwalkway_dict['Direction'] = Direction
    movingwalkway_dict['Width'] = width
    movingwalkway_dict['HorizontalLength'] = HorizontalLength
    movingwalkway_dict['Lanes'] = math.floor(width/0.7)
    movingwalkway_dict['xpos'] = lx
    movingwalkway_dict['ypos'] = ly
    movingwalkway_dict['Height'] = 0.0
    
    movingwalkway_dict['nodedata'] = [g_StairNodeID, lx, ly, NdeHeight, Elevation]
    movingwalkway_dict['nodeid'] = g_StairNodeID
    g_StairNodeID += 1
    
    g_OMA_Class.m_movingwalkway_list.append(movingwalkway_dict)
    
    print("Movingwalkway data ", movingwalkway_dict)


def transport_data(ifc_file, settings):
    transports = ifc_file.by_type("IfcTransportElement")
    for transport in transports:
        Elevation = None
        aRelatingSpace = exogetifcparam.parent(transport)
        if aRelatingSpace:
            if aRelatingSpace.is_a("IfcSpace"):
                if aRelatingSpace.BoundedBy:
                    Elevation = exogetifcparam.get_storey_Elevation(aRelatingSpace.BoundedBy[0])
            elif aRelatingSpace.is_a("IfcBuildingStorey"):
                element_info = aRelatingSpace.get_info()
                Elevation = exogetifcparam.get_property(element_info, 'Elevation')
        
        transport_info = transport.get_info()
        if transport.PredefinedType == 'ESCALATOR':
            get_escalator_data(transport_info, settings, transport, Elevation)
        elif transport.PredefinedType == 'ELEVATOR':
            get_elevator_data(transport_info, settings, transport, Elevation)
        elif transport.PredefinedType == 'MOVINGWALKWAY':
            get_movingwalkway_data(transport_info, settings, transport, Elevation)


def extract_evelvator_data_from_group(IfcRelAssignsToGroups):
    global g_OMA_Class
    theElevatorPos = -1
    spaceIndexList = []
    doorIndexList = []
    if IfcRelAssignsToGroups is not None:
        for groups in IfcRelAssignsToGroups:
            for aRelatingSpace in groups.RelatedObjects:
                if aRelatingSpace.is_a("IfcSpace"):
                    print("IfcSpace ",aRelatingSpace.Name)
                    spaceIndex = g_OMA_Class.SpaceDefined(aRelatingSpace.GlobalId)
                    if spaceIndex > -1:
                        spaceIndexList.append(spaceIndex)
                        if aRelatingSpace.ContainsElements is not None:
                            for contained_element in aRelatingSpace.ContainsElements:
                                if contained_element.is_a("IfcRelContainedInSpatialStructure"):
                                    for element in contained_element.RelatedElements:
                                        if element.is_a("IfcTransportElement"):
                                            aElevatorPos = exoifcutils.ElevatorDefined(element.GlobalId, g_OMA_Class.m_elevator_list)
                                            if aElevatorPos>-1:
                                                theElevatorPos = aElevatorPos
                                
                elif aRelatingSpace.is_a("IfcTransportElement"):
                    print("IfcTransportElement ",aRelatingSpace.Name)
                    aElevatorPos = exoifcutils.ElevatorDefined(aRelatingSpace.GlobalId, g_OMA_Class.m_elevator_list)
                    if aElevatorPos>-1:
                        theElevatorPos = aElevatorPos
                elif aRelatingSpace.is_a("IfcDoor"):
                    doorIndex = exoifcutils.DoorDefined(aRelatingSpace.GlobalId, g_OMA_Class.m_door_list)
                    doorIndexList.append(doorIndex)
                            
    if theElevatorPos > -1:
        theElevator = g_OMA_Class.m_elevator_list[theElevatorPos]
        if theElevator is not None:
            if len(spaceIndexList)>0:
                theElevator['SpaceIndexList'] = spaceIndexList
                for space_index in spaceIndexList:
                    a_space = g_OMA_Class.m_space_list[space_index]
                    if a_space is not None:
                        a_space['elevator']=theElevatorPos
            if len(doorIndexList)>0:
                theElevator['DoorIndexList'] = doorIndexList
                        

def extract_elevator_data(ifc_file):
    zones = ifc_file.by_type("IfcZone")
    for zone in zones:
        if zone.ObjectType=='ElevatorShaft':
            print("elevator zone found:", zone.Name)
            IfcRelAssignsToGroups = zone.IsGroupedBy
            if IfcRelAssignsToGroups is not None:
                extract_evelvator_data_from_group(IfcRelAssignsToGroups)

    built_systems = ifc_file.by_type("IfcBuiltSystem")
    for built_system in built_systems:
        if built_system.PredefinedType=='TRANSPORT':
            print("elevator zone found:", built_system.Name)
            IfcRelAssignsToGroups = built_system.IsGroupedBy
            if IfcRelAssignsToGroups is not None:
                extract_evelvator_data_from_group(IfcRelAssignsToGroups)

    building_systems = ifc_file.by_type("IfcBuildingSystem")
    for building_system in building_systems:
        if building_system.PredefinedType=='TRANSPORT':
            print("elevator zone found:", building_system.Name)
            IfcRelAssignsToGroups = building_system.IsGroupedBy
            if IfcRelAssignsToGroups is not None:
                extract_evelvator_data_from_group(IfcRelAssignsToGroups)       
      

def CalcAngle(upper_points, lower_points):
    '''
    Calculates the direction of an object based on it's upper and lower points
    for example the lower and upper lines of a IfcStairFlight or IfcRampFlight
    '''
    UpperCount = len(upper_points)
    LowerCount = len(lower_points)
    if UpperCount==0 or LowerCount==0:
        return 0.0;
    Upper = [0,0]
    for point in upper_points:
        Upper[0] += point[0]
        Upper[1] += point[1]

    Upper[0] = Upper[0]/UpperCount
    Upper[1] = Upper[1]/UpperCount
    Lower = [0,0]
    for point in lower_points:
        Lower[0] += point[0]
        Lower[1] += point[1]

    Lower[0] = Lower[0]/LowerCount
    Lower[1] = Lower[1]/LowerCount

    lineA = [Lower,Upper]
    lineB = [Lower,[Upper[0],1]]
    return exoifcutils.calcAngle(lineA,lineB)


def ramp_data(ifc_file, settings):
    global g_OMA_Class
    global g_StairNodeID
    
    ramps = ifc_file.by_type("IfcRamp")
    for ramp in ramps:
        if ramp.IsDecomposedBy:
            Elevation = exogetifcparam.get_storey_Elevation(ramp)
            
            ramp_info = ramp.get_info()
            print(ramp.PredefinedType, ramp_info['Name'])
            ramp_dict = {'id': g_StairNodeID}
            g_StairNodeID += 1
            exogetifcparam.get_basic_ramp_info(ramp_dict, ramp)
            ramp_dict['Elevation'] = Elevation
            ramp_components = ramp.IsDecomposedBy[0].RelatedObjects
            for ramp_element in ramp_components:
                if 'IfcRampFlight' == ramp_element.is_a():
                    entity_shape = ifcopenshell.geom.create_shape(settings, ramp_element)
                    lx, ly, NdeHeight, Direction, width, HorizontalLength = level_geom(entity_shape.geometry, Elevation)
                    plot_geom(entity_shape.geometry)
                    geom_data = level_flight_geom(entity_shape.geometry, Elevation, 0.01)
                    ramp_dict['GeomData'] = geom_data
                    ramp_dict['ExoAngle'] = CalcAngle(geom_data['UpperPoints'],geom_data['LowerPoints'])
                    ramp_dict['Angle'] = Direction
                    #ramp_dict['Floor'] = floor_longname
                    #ramp_dict['FloorIndex'] = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
                    ramp_dict['Direction'] = Direction
                    ramp_dict['Width'] = width
                    ramp_dict['HorizontalLength'] = HorizontalLength
                    ramp_dict['Lanes'] = math.floor(width/0.7)
                    ramp_dict['xpos'] = lx
                    ramp_dict['ypos'] = ly
                    ramp_dict['Height'] = NdeHeight
                    
                    ramp_dict['nodedata'] = [g_StairNodeID, lx, ly, NdeHeight, Elevation]
                    ramp_dict['nodeid'] = g_StairNodeID
                    g_StairNodeID += 1
                    
            print(ramp_dict)
            g_OMA_Class.m_ramp_list.append(ramp_dict)


def sign_data(ifc_file, settings):
    global g_OMA_Class
    
    signs = ifc_file.by_type("IfcSign")
    sign_count = 0
    for ifcsign in signs:
        sign_count+=1
        sign_info = ifcsign.get_info()
        sign_dict = {'id': sign_count}
        print(ifcsign.PredefinedType, sign_info['Name'])
        exogetifcparam.get_basic_sign_info(sign_dict, ifcsign)
        
        Elevation = exogetifcparam.get_entity_storey_elevation(ifcsign)
        sign_dict['Elevation'] = Elevation
        entity_shape = ifcopenshell.geom.create_shape(settings, ifcsign)
        lx, ly, NdeHeight, Direction, width, HorizontalLength = level_geom(entity_shape.geometry, Elevation)
        startloc, endloc = sign_locations_xy(entity_shape.geometry)
        plot_geom(entity_shape.geometry)
        sign_dict['HeightFromTheGroundCalc'] = NdeHeight
        sign_dict['Direction'] = Direction
        sign_dict['StartLoc'] = startloc
        sign_dict['EndLoc'] = endloc

        g_OMA_Class.m_sign_list.append(sign_dict)


def road_data(ifc_file, settings):
    global g_OMA_Class
    
    roads = ifc_file.by_type("IfcRoad")
    road_count = 0
    for ifcroad in roads:
        road_count+=1
        road_info = ifcroad.get_info()
        road_dict = {'id': road_count}
        print(ifcroad.PredefinedType, road_info['Name'])
        
        Elevation = exogetifcparam.get_entity_storey_elevation(ifcroad)
        road_dict['Elevation'] = Elevation
        if ifcroad.Representation is not None:
            entity_shape = ifcopenshell.geom.create_shape(settings, ifcroad)
            plot_geom(entity_shape.geometry)
            
        if ifcroad.IsDecomposedBy:
            Elevation = exogetifcparam.get_storey_Elevation(ifcroad)
            
            road_components = ifcroad.IsDecomposedBy[0].RelatedObjects
            for road_element in road_components:
                if road_element.Representation is not None:
                    entity_shape = ifcopenshell.geom.create_shape(settings, road_element)
                    plot_geom(entity_shape.geometry)
                print(road_element.Name)
                if road_element.IsDecomposedBy:
                    roadpart_components = road_element.IsDecomposedBy[0].RelatedObjects
                    for roadpart_element in roadpart_components:
                        print(road_element.Name, "Subpart ", roadpart_element.Name)
                        if roadpart_element.Representation is not None:
                            entity_shape = ifcopenshell.geom.create_shape(settings, roadpart_element)
                            plot_geom(entity_shape.geometry)
                            startloc, endloc = sign_locations_xy(entity_shape.geometry)
                            print(startloc,endloc)
                        if roadpart_element.ContainsElements is not None:
                            for contained_element in roadpart_element.ContainsElements:
                                if contained_element.is_a("IfcRelContainedInSpatialStructure"):
                                    for element in contained_element.RelatedElements:
                                        if element.is_a("IfcCourse"):
                                            print(element.Name, "Subpart ", roadpart_element.Name," IfcCourse ",element.Name)
                                            entity_shape = ifcopenshell.geom.create_shape(settings, element)
                                            plot_geom(entity_shape.geometry)
                                            startloc, endloc = sign_locations_xy(entity_shape.geometry)
                                            print(startloc,endloc)  


def bridge_data(ifc_file, settings):
    global g_OMA_Class
    
    bridges = ifc_file.by_type("IfcBridge")
    bridge_count = 0
    for ifcbridge in bridges:
        bridge_count+=1
        bridge_info = ifcbridge.get_info()
        bridge_dict = {'id': bridge_count}
        print(ifcbridge.PredefinedType, bridge_info['Name'])
        
        Elevation = exogetifcparam.get_entity_storey_elevation(ifcbridge)
        bridge_dict['Elevation'] = Elevation
        if ifcbridge.Representation is not None:
            entity_shape = ifcopenshell.geom.create_shape(settings, ifcbridge)
            lx, ly, NdeHeight, Direction, width, HorizontalLength = level_geom(entity_shape.geometry, Elevation)
            startloc, endloc = sign_locations_xy(entity_shape.geometry)
            plot_geom(entity_shape.geometry)
            bridge_dict['HeightFromTheGroundCalc'] = NdeHeight
            bridge_dict['Direction'] = Direction
            bridge_dict['StartLoc'] = startloc
            bridge_dict['EndLoc'] = endloc

        #g_OMA_Class.m_sign_list.append(sign_dict)


def scan_spaces(ifc_file, settings):
    print("-------------- IfcSpace ---------------")
    ifc_spaces = ifc_file.by_type('IfcSpace')
    for ifc_space in ifc_spaces:
        output_name = False
        space_loc = -1
        for boundary in ifc_space.BoundedBy:
            if boundary.is_a("IfcRelSpaceBoundary"): # Should check name, which shoule be 2ndLevel
                if space_loc == -1:
                    Elevation = exogetifcparam.get_storey_Elevation_spaceboundary_quiet(boundary)
                    floor_longname = exogetifcparam.get_storey_Longname_spaceboundary(boundary)
                    # now load and add space data, returns location in g_space_data
                    space_loc = get_space_data(ifc_file, ifc_space, Elevation, floor_longname, settings)
                    
        if space_loc>-1:
            if ifc_space.ContainsElements and len(ifc_space.ContainsElements)>0:
                furnitureList = []
                objectList = []
                for ContainsElement in ifc_space.ContainsElements:
                    if ContainsElement.RelatedElements is not None:
                        for elem in ContainsElement.RelatedElements:
                            if elem:
                                if elem.is_a("IfcFurniture"):
                                    furnitureList.append(elem.Name)
                                    objectList.append(elem.ObjectType)
                                elif elem.is_a("IfcTransportElement"):
                                    print("ContainsElement: Transport ",elem.Name)
                if len(furnitureList)>0:
                    g_OMA_Class.m_space_list[space_loc]['furniture'] = furnitureList
                if len(objectList)>0:
                    g_OMA_Class.m_space_list[space_loc]['objects'] = objectList


def build_door_list(ifc_file):
    global g_OMA_Class
    
    ifc_doors = ifc_file.by_type("IfcDoor")
    for ifc_door in ifc_doors:
        a_story = Element.get_container(ifc_door)
        element_info = a_story.get_info()
        Elevation = exogetifcparam.get_property(element_info, 'Elevation')
        floor_longname = exogetifcparam.get_property(element_info, 'LongName')
        floor_index = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
        if floor_index < 0:
            floor = {'Elevation': Elevation}
            floor['Name'] = floor_longname
            g_OMA_Class.m_floor_list.append(floor)
            
        plot_item = True
        if plot_item:
            exogetifcparam.output_basic_door_info(ifc_door)
            door_data = exogetifcparam.get_basic_door_info(ifc_door)
            door_data['Elevation'] = Elevation
            door_data['Floor'] = floor_longname
            #door_data['FloorIndex'] = exoifcutils.find_floor(Elevation, g_OMA_Class.m_floor_list)
            door_data['Spaces'] = []
            g_OMA_Class.m_door_list.append(door_data)


def connect_rooms_doors(ifc_file):
    global g_OMA_Class
    
    ifc_doors = ifc_file.by_type("IfcDoor")
    for ifc_door in ifc_doors:
        element_info = ifc_door.get_info()
        space_id = exogetifcparam.get_property(element_info, 'GlobalId')
        doorIndex = exoifcutils.DoorDefined(space_id, g_OMA_Class.m_door_list)
        if doorIndex > -1:
            door_data = g_OMA_Class.m_door_list[doorIndex]
            spaceIndexList = []
            for SpaceID in door_data['Spaces']:
                spaceIndex = g_OMA_Class.SpaceDefined(SpaceID)
                if spaceIndex > -1:
                    spaceIndexList.append(spaceIndex)

            # check doors are assign correctly
            if len(spaceIndexList) % 2 == 0:
                pass
            else:
                if len(spaceIndexList)>2:
                    print("Warning: Door ", ifc_door.Name,", has odd number of connections,", len(spaceIndexList)," greater than 1")
                # door has an odd number of connections
                if door_data['IsExternal'] == False:
                    # not marked as external
                    print("Warning: Door ", ifc_door.Name,", has odd number of connections and not marked as external")
                    door_data['IsExternal'] = True


def build_conections_between_spaces():
    """
      assign 'space_connecting' for ifc_stair, ifc_stair_flights and ifc_escalators
    """
    global g_OMA_Class
    
    for space in g_OMA_Class.m_space_list:
        GlobalID = space['GlobalId']
        if 'stairsUp' in space or 'stairsDown' in space:
            if 'stairsUp' in space:
                for stairID in space['stairsUp']:
                    ifc_stair = g_OMA_Class.m_stair_list[stairID]
                    if 'space_connecting' in ifc_stair:
                        if not(GlobalID in ifc_stair['space_connecting']):
                            ifc_stair['space_connecting'].append(GlobalID)
                    else:
                        ifc_stair['space_connecting'] = [GlobalID]
            if 'stairsDown' in space:
                for stairID in space['stairsDown']:
                    ifc_stair = g_OMA_Class.m_stair_list[stairID]
                    if 'space_connecting' in ifc_stair:
                        if not(GlobalID in ifc_stair['space_connecting']):
                            ifc_stair['space_connecting'].append(GlobalID)
                    else:
                        ifc_stair['space_connecting'] = [GlobalID]
            
        if 'stairflightsUp' in space or 'stairflightsDown' in space:
            if 'stairflightsUp' in space:
                for stairID in space['stairflightsUp']:
                    ifc_stair_flight = g_OMA_Class.m_stair_flights_list[stairID]
                    if 'space_connecting' in ifc_stair_flight:
                        if not(GlobalID in ifc_stair_flight['space_connecting']):
                            ifc_stair_flight['space_connecting'].append(GlobalID)
                    else:
                        ifc_stair_flight['space_connecting'] = [GlobalID]
            if 'stairflightsDown' in space:
                for stairID in space['stairflightsDown']:
                    ifc_stair_flight = g_OMA_Class.m_stair_flights_list[stairID]
                    if 'space_connecting' in ifc_stair_flight:
                        if not(GlobalID in ifc_stair_flight['space_connecting']):
                            ifc_stair_flight['space_connecting'].append(GlobalID)
                    else:
                        ifc_stair_flight['space_connecting'] = [GlobalID]
                        
        if 'escalatorsUp' in space or 'escalatorsDown' in space:
            if 'escalatorsUp' in space:
                print(" escalatorsUp:", space['escalatorsUp'])
                for stairID in space['escalatorsUp']:
                    ifc_escalator = g_OMA_Class.m_escalator_list[stairID]
                    if 'space_connecting' in ifc_escalator:
                        if not(GlobalID in ifc_escalator['space_connecting']):
                            ifc_escalator['space_connecting'].append(GlobalID)
                    else:
                        ifc_escalator['space_connecting'] = [GlobalID]
            if 'escalatorsDown' in space:
                for stairID in space['escalatorsDown']:
                    ifc_escalator = g_OMA_Class.m_escalator_list[stairID]
                    if 'space_connecting' in ifc_escalator:
                        if not(GlobalID in ifc_escalator['space_connecting']):
                            ifc_escalator['space_connecting'].append(GlobalID)
                    else:
                        ifc_escalator['space_connecting'] = [GlobalID]


def calc_stair_travel_distance():
    """
    Calculates stair total length
    Note Routine not complete - the calculation of length/travel distanc is incomplete 15/03/24
    Note Width is based on maximum stairflight width - maybe should be minimum 15/03/24
    
    """
    global g_OMA_Class
    
    for stair in g_OMA_Class.m_stair_list:
        aWidth = 0.0
        aLength = 0.0
        for part in stair['elemIDs']:
            GlobalId = part[1]
            if part[0] == 'IfcStairFlight':
                stair_flight_index = exoifcutils.StairflightDefined(GlobalId, g_OMA_Class.m_stair_flights_list)
                if stair_flight_index > -1:
                    stair_flight = g_OMA_Class.m_stair_flights_list[stair_flight_index]
                    aWidth = max(aWidth, stair_flight['Width'])
                    travel_distance = math.sqrt(stair_flight['TotalRun']*stair_flight['TotalRun']+stair_flight['Height']*stair_flight['Height'])
                    aLength += travel_distance
            elif part[0] == 'IfcSlab':
                landing_index = exoifcutils.LandingDefined(GlobalId, g_OMA_Class.m_landings_list)
                if landing_index > -1:
                    landing = g_OMA_Class.m_landings_list[landing_index]
                    # TO DO
                    #aWidth=max(aWidth,landing['Width'])
                    #aLength+=landing['Length']
                
        stair['TravelDistance'] = aLength
        stair['Width'] = aWidth


def adjust_escalator_height(escalator):
    '''
    calculate escalator height based for now on the assumption
    it's the distance between floors
    '''
    global g_OMA_Class
    if 'Elevation' in escalator:
        elevation = escalator['Elevation']
        above_floor_elevaction = exoifcutils.find_floor_above_elevation(elevation, g_OMA_Class.m_floor_list)
        if above_floor_elevaction is not None:
           escalator['Height'] = above_floor_elevaction-elevation


def generate_data(IFC_filename, output_file, output_type):
    """
     main routine to generate output file/data

    """
    global g_OMA_Class
    global g_plotter_class
    
    print("Loading IFC file ", IFC_filename)
    ifc_file = ifcopenshell.open(IFC_filename)
    print("Load complete")

    settings = geom.settings()
    settings.set(settings.USE_WORLD_COORDS, True)

    build_door_list(ifc_file)

    stairs_data(ifc_file, settings)

    transport_data(ifc_file, settings)

    ramp_data(ifc_file, settings)

    sign_data(ifc_file, settings)

    road_data(ifc_file, settings)

    bridge_data(ifc_file, settings)

    # sorts floors before finding spaces - as we require an ordered list for linking levels
    ordered_floor_list = sorted(g_OMA_Class.m_floor_list, key=itemgetter('Elevation'))
    g_OMA_Class.m_floor_list = ordered_floor_list

    # Now need to correct FloorIndex values in doors, stairs, transport and ramps
    for door in g_OMA_Class.m_door_list:
        door['FloorIndex'] = exoifcutils.find_floor(door['Elevation'], g_OMA_Class.m_floor_list)
        
    for stair_flight in g_OMA_Class.m_stair_flights_list:
        stair_flight['FloorIndex'] = exoifcutils.find_floor(stair_flight['Elevation'], g_OMA_Class.m_floor_list)
        
    for landing in g_OMA_Class.m_landings_list:
        landing['FloorIndex'] = exoifcutils.find_floor(landing['Elevation'], g_OMA_Class.m_floor_list)
        
    for escalator in g_OMA_Class.m_escalator_list:
        escalator['FloorIndex'] = exoifcutils.find_floor(escalator['Elevation'], g_OMA_Class.m_floor_list)
        adjust_escalator_height(escalator)
        
    for elevator in g_OMA_Class.m_elevator_list:
        elevator['FloorIndex'] = exoifcutils.find_floor(elevator['Elevation'], g_OMA_Class.m_floor_list)
        
    for movingwalkway in g_OMA_Class.m_movingwalkway_list:
        movingwalkway['FloorIndex'] = exoifcutils.find_floor(elevator['Elevation'], g_OMA_Class.m_floor_list)

    for ramp in g_OMA_Class.m_ramp_list:
        ramp['FloorIndex'] = exoifcutils.find_floor(ramp['Elevation'], g_OMA_Class.m_floor_list)

    for sign in g_OMA_Class.m_sign_list:
        sign['FloorIndex'] = exoifcutils.find_floor(sign['Elevation'], g_OMA_Class.m_floor_list)
    
    scan_spaces(ifc_file, settings)
    extract_elevator_data(ifc_file)
    connect_rooms_doors(ifc_file)

    print("===================== door data ====================")
    #for door in g_OMA_Class.m_door_list:
    #    print(door)
    print(len(g_OMA_Class.m_door_list))

    print("===================== spaces data ====================")
    #for space in g_OMA_Class.m_space_list:
    #    print(space)
    print(len(g_OMA_Class.m_space_list))

    print("===================== floor data ====================")
    for floor in g_OMA_Class.m_floor_list:
        print(floor)
        
    print("===================== update stair connectionsdata ====================")
    build_conections_between_spaces()         
    
    # match g_OutputType: 
    if output_type == "Exodus":
        exodusmesh.exodus_ouput(g_OMA_Class, ifc_file, settings, g_plotter_class, output_file)
    elif output_type == "EvacutionZ":
        calc_stair_travel_distance()
        enzoutput.enz_output(output_file, ifc_config.EZPopFile, g_OMA_Class)
    elif output_type == "Graph":
        graphoutput.graph_view(g_OMA_Class, ifc_file, output_file)
    elif output_type == "CFAST":
        cfastoutput.cfast_output(g_OMA_Class, ifc_file, settings, output_file)
    elif output_type == "DumpStairs":
        g_OMA_Class.dump_data()
        g_OMA_Class.dump_ordered_space_boundaries()
    else:
        print("Unsupported format")

    if g_plotter_class is not None:
        g_plotter_class.show_geom()


def main(argv):
    """
     main routine

    """
    global g_plotter_class
    
    if len(argv)>0:
        # use data entered on command line
        IFCFilename = argv[0]
        
        if len(argv)>1:
            output_file = argv[1]
            pre, ext = os.path.splitext(output_file)
            ext = ext.lower()
            print(ext)
            if ext == '.mta':
                output_type = "Exodus"
            elif ext == '.xml':
                output_type = "EvacutionZ"
            elif ext == '.html':
                output_type = "Graph"
            elif ext == '.in':
                output_type = "CFAST"
            else:
              print("Unsupported format")
              return  
        else:
            output_type = "DumpStairs"
            output_file = None
    
    else:
        # use data defined in ifc_config.py
        IFCFilename = ifc_config.IFCFile

        output_type = ifc_config.OutputType
        if output_type == "Exodus":
            output_file = ifc_config.MTAFile
        elif output_type == "EvacutionZ":
             output_file = ifc_config.EZMapFile
        elif output_type == "Graph":
            output_file = ifc_config.GraphFile
        elif output_type == "CFAST":
            output_file = ifc_config.CFastFile
        elif output_type == "DumpStairs":
            output_file = None
        else:
            print("Unsupported format")
            return

    if ifc_config.plot_model:
        g_plotter_class = PlotterClass()
    
    generate_data(IFCFilename, output_file, output_type)


if __name__ == "__main__":
    main(sys.argv[1:])
