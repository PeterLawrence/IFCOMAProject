# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence
"""

import exoifcutils
from oma_class import OMAClass
import ifcopenshell


def calc_door_centre(shape):
    """ Calc the region around a door for connecting doors """
    verts = shape.verts # X Y Z of vertices in flattened list e.g. [v1x, v1y, v1z, v2x, v2y, v2z, ...]
    VertexCount = len(verts)
    x_sum = 0
    y_sum = 0
    iPoint = 0
    while iPoint < VertexCount:
        x = verts[iPoint]
        iPoint += 1
        y = verts[iPoint]
        iPoint += 1
        iPoint += 1
        x_sum+=x
        y_sum+=y
    if VertexCount > 3:
        x_sum = x_sum/(VertexCount/3)
        y_sum = y_sum/(VertexCount/3)
    return x_sum, y_sum


def calc_door_centre_ifc(GlobalID, ifc_file, settings):
    ifc_door = ifc_file.by_id(GlobalID)
    entity_shape = ifcopenshell.geom.create_shape(settings, ifc_door)
    x_centre, y_centre = calc_door_centre(entity_shape.geometry)
    return x_centre, y_centre
    

def get_width_height_location(boundary_points):
    if boundary_points is None:
        return False
    count = len(boundary_points)
    if count == 0:
        return False
    x = boundary_points[0][0]
    y = boundary_points[0][1]
    min_x = x
    min_y = y
    max_x = x
    max_y = y
    for point in range(1, count):
        x = boundary_points[point][0]
        y = boundary_points[point][1]
        min_x = min(min_x,x)
        min_y = min(min_y,y)
        max_x = max(max_x,x)
        max_y = max(max_y,y)

    width = max_x - min_x
    depth = max_y - min_y
    return width, depth, min_x, min_y
    

def output_cfast_header(cfast_file):
    print("&HEAD VERSION = 7500, TITLE = 'Users Guide Example Case' /", file=cfast_file)
    print("", file=cfast_file)


def output_cfast_scenario(cfast_file):
    print("!! Scenario Configuration", file=cfast_file)
    print("&TIME SIMULATION = 3600 PRINT = 60 SMOKEVIEW = 60 SPREADSHEET = 60 / ", file=cfast_file)
    print("&INIT PRESSURE = 101325 RELATIVE_HUMIDITY = 50 INTERIOR_TEMPERATURE = 20 EXTERIOR_TEMPERATURE = 20 /", file=cfast_file)
    print("&MISC  LOWER_OXYGEN_LIMIT = 0.1 / ", file=cfast_file)
    print("", file=cfast_file)


def output_cfast_material(cfast_file):
    print("!! Material Properties", file=cfast_file)
    print("&MATL ID = 'CONCRETE' MATERIAL = 'Concrete Normal Weight (6 in)', ", file=cfast_file)
    print("      CONDUCTIVITY = 1.75 DENSITY = 2200 SPECIFIC_HEAT = 1, THICKNESS = 0.15 EMISSIVITY = 0.94 /", file=cfast_file)
    print("", file=cfast_file)


def output_cfast_compartment(cfast_file, name, depth, height, width, x, y, z):
    print(f"&COMP ID = '{name}'", file=cfast_file)
    print(f"      DEPTH = {depth} HEIGHT = {height} WIDTH = {width} CEILING_MATL_ID = 'CONCRETE' WALL_MATL_ID = 'CONCRETE' FLOOR_MATL_ID = 'CONCRETE'", file=cfast_file)
    print(f"      ORIGIN = {x}, {y}, {z} GRID = 50, 50, 50 LEAK_AREA_RATIO = 0.00017, 5.2E-05 /", file=cfast_file)
    

def output_cfast_wall_vent(cfast_file, name, comp1, comp2, top, bottom, width, face, offset):
    print(f"&VENT TYPE = 'WALL' ID = '{name}' COMP_IDS = '{comp1}' '{comp2}'  TOP = {top}, BOTTOM = {bottom}, WIDTH = {width}", file=cfast_file)
    print(f"  FACE = '{face}' OFFSET = {offset} /", file=cfast_file)


def output_cfast_ceiling_vent(cfast_file, name, comp1, comp2, area, offset_x, offset_y):
    print(f"&VENT TYPE = 'CEILING' ID = '{name}' COMP_IDS = '{comp1}' '{comp2}'  AREA = {area}, SHAPE = 'SQUARE'", file=cfast_file)
    print(f"      OFFSETS = {offset_x}, {offset_y} /", file=cfast_file)


def output_cfast_floor_vent(cfast_file, name, comp1, comp2, area, offset_x, offset_y):
    print(f"&VENT TYPE = 'FLOOR' ID = '{name}' COMP_IDS = '{comp1}' '{comp2}'  AREA = {area}, SHAPE = 'SQUARE'", file=cfast_file)
    print(f"      OFFSETS = {offset_x}, {offset_y} /", file=cfast_file)


def find_wall_face(xmin, xmax, ymin, ymax, x_centre, y_centre):
    left_dist = abs(xmin-x_centre)
    right_dist = abs(xmax-x_centre)
    rear_dist = abs(ymax-y_centre)
    front_dist = abs(ymin-y_centre)
    if left_dist<right_dist:
        # LEFT side
        if rear_dist<front_dist:
            if rear_dist<left_dist:
                face = 'REAR'
            else:
                face = 'LEFT'
        else: # rear_dist>=front_dist
            if front_dist<left_dist:
                face = 'FRONT'
            else:
                face = 'LEFT'
    else: # left_dist>=right_dist
        # RIGHT side
        if rear_dist<front_dist:
            if rear_dist<right_dist:
                face = 'REAR'
            else:
                face = 'RIGHT'
        else: # rear_dist>=front_dist
            if front_dist<right_dist:
                face = 'FRONT'
            else:
                face = 'RIGHT'
    return face
      

def cfast_output(OMA_Class, ifc_file, settings, filename):
    # Save CFATS data
    cfast_file = open(filename, "w")
    output_cfast_header(cfast_file)
    output_cfast_scenario(cfast_file)
    output_cfast_material(cfast_file)
    
    print("===================== build network ====================")
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


    print("!! Compartments", file=cfast_file)
    compartment_id = 0
    for ifc_space in OMA_Class.m_space_list:
        output_node = False
        if 'elemIDs' in ifc_space:
            output_node = (len(ifc_space['elemIDs']) > 0)
        
        if output_node:
            compartment_id += 1

            name = ifc_space['Name']
            elevation = ifc_space['Elevation']
            above_floor_elevaction = exoifcutils.find_floor_above_elevation(elevation, OMA_Class.m_floor_list)
            if above_floor_elevaction is not None:
                height = above_floor_elevaction - elevation
            else:
                height = 3.0

            x = 0
            y = 0
            depth = 10.0
            width = 10.0
            if 'boundarylist' in ifc_space:
                width, depth, x, y = get_width_height_location(ifc_space['boundarylist'])

            print("!! ",x, y, x+width, y+depth , file=cfast_file)
            output_cfast_compartment(cfast_file, name, depth, height, width, x, y, elevation)
            ifc_space['compartment_id'] = compartment_id
            ifc_space['cfast_box_data'] = (width, depth, x, y)
            
    print("", file=cfast_file)

    print("!! Wall Vents", file=cfast_file)
    print("", file=cfast_file)

    linkcount = 0
    
    for ifc_door in OMA_Class.m_door_list:
        if 'Spaces' in ifc_door:
            connectlist = []
            # normally one or two connected spaces
            for spaceid in ifc_door['Spaces']:
                space_index = OMA_Class.SpaceDefined(spaceid)
                if space_index > - 1:
                    ifc_space = OMA_Class.m_space_list[space_index]
                    if 'compartment_id' in ifc_space:
                        connectlist.append(ifc_space)
                    else:
                        print("rejected: door compartment_id", space_index, ifc_space)
                else:
                    print("rejected: door index", space_index, spaceid)

            if len(connectlist)>0 and 'cfast_box_data' in connectlist[0]:
                the_compartment = connectlist[0]
                the_compartment_name = the_compartment['Name']
                if len(connectlist) == 1:
                    the_compartment_name = the_compartment['Name']
                    if ifc_door['IsExternal']:
                        conx_compartment_name = 'OUTSIDE'
                    else:
                        conx_compartment_name = 'OUTSIDE'
                elif len(connectlist) == 2:
                    conx_compartment = connectlist[1]
                    conx_compartment_name = conx_compartment['Name']
                    '''
                    if conx_compartment['area']<the_compartment['area']:
                        # use the smallest compartment as the 2nd compartment for wall (used for testing)
                        the_compartment = connectlist[1]
                        conx_compartment = connectlist[0]
                        the_compartment_name = the_compartment['Name']
                        conx_compartment_name = conx_compartment['Name']
                        
                    '''
                else:
                    conx_compartment_name = None

                if conx_compartment_name is not None:
                    linkcount += 1
                    width = ifc_door['OverallWidth']
                    x_centre, y_centre = calc_door_centre_ifc(ifc_door['GlobalId'], ifc_file, settings)
                    xmin = the_compartment['cfast_box_data'][2]
                    ymin = the_compartment['cfast_box_data'][3]
                    xmax = xmin + the_compartment['cfast_box_data'][0]
                    ymax = ymin + the_compartment['cfast_box_data'][1]
                    face = find_wall_face(xmin, xmax, ymin, ymax, x_centre, y_centre) # RIGHT, FRONT, LEFT, or REAR
                    
                    if face == 'FRONT' or face == 'REAR':
                        offset = x_centre - xmin -  width/2
                    else:
                        offset = y_centre - ymin -  width/2
                        
                                
                    bottom = 0.0
                    top = ifc_door['OverallHeight']
                    
                    print("!! ",xmin, xmax, ymin, ymax, x_centre, y_centre , file=cfast_file)
                    output_cfast_wall_vent(cfast_file, ifc_door['Name'], the_compartment_name, conx_compartment_name, top, bottom, width, face, offset)
                else:
                    print("Rejected: door connnect {ifc_door['Name']} ID:{ifc_door['GlobalId']} spaces {ifc_door['Spaces']} connections {connectlist}")

    print("!! Ceiling and Floor Vents", file=cfast_file)
    print("", file=cfast_file)
    
    for ifc_stair in OMA_Class.m_stair_list:
        if 'space_connecting' in ifc_stair:
            if len(ifc_stair['space_connecting']) == 2:
                Spaces = []
                for SpaceGlobalID in ifc_stair['space_connecting']:
                    index = OMA_Class.SpaceDefined(SpaceGlobalID)
                    if index > -1:
                        ifc_space = OMA_Class.m_space_list[index]
                        if 'compartment_id' in ifc_space:
                            Spaces.append(ifc_space)
                if len(Spaces)>0:
                    linkcount += 1
                    width = Spaces[0]['cfast_box_data'][0]
                    dept = Spaces[0]['cfast_box_data'][1]
                    gap = 0.2
                    x_offset = width-2.0*gap
                    y_offset = dept-2.0*gap
                    vent_area = x_offset*y_offset
                    x_offset = x_offset/2
                    y_offset = y_offset/2
                    if len(Spaces) == 2:
                        if Spaces[0]['Elevation']<Spaces[1]['Elevation']:
                            output_cfast_ceiling_vent(cfast_file, ifc_stair['Name'], Spaces[1]['Name'], Spaces[0]['Name'], vent_area, x_offset, y_offset)
                        else:
                            output_cfast_floor_vent(cfast_file, ifc_stair['Name'], Spaces[1]['Name'], Spaces[0]['Name'], vent_area, x_offset, y_offset)
                    else:
                        output_cfast_ceiling_vent(cfast_file, ifc_stair['Name'], "OUTSIDE", Spaces[0]['Name'], vent_area, x_offset, y_offset)

    print("", file=cfast_file)
    cfast_file.close()
    print("Max links", linkcount)
