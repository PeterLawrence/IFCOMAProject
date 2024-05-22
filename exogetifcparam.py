"""
Created on Fri Dec 9 09:59:58 2023

@author: P.J.Lawrence
"""
import ifcopenshell


def parent(instance):
    if instance.is_a("IfcOpeningElement"):
        return instance.VoidsElements[0].RelatingBuildingElement
    if instance.is_a("IfcElement"):
        fills = instance.FillsVoids
        if len(fills):
            return fills[0].RelatingOpeningElement
        contains = instance.ContainedInStructure
        if len(contains):
            return contains[0].RelatingStructure
    if instance.is_a("IfcObjectDefinition"):
        decompositions = instance.Decomposes
        if len(decompositions):
            return decompositions[0].RelatingObject


def get_property(element_info, tag):
    if tag in element_info:
        tag_property = element_info[tag]
        return tag_property
    return None


def outputElement(element_info, tag, spaces):
    for aSpace in range(0, spaces):
        print(" ", end='')
    if tag in element_info:
        tag_property = element_info[tag]
        print(f"{tag}: {tag_property}")
    else:
        print(f"{tag}: UNDEFINED")


def output_basic_door_info(door):
    print("=========== IFCDoor ===========")
    door_info = door.get_info()
    outputElement(door_info, "Name", 2)
    outputElement(door_info, "OverallWidth", 2)
    outputElement(door_info, "OverallHeight", 2)
    predefined_type = get_property(door_info, 'PredefinedType')
    print(f"  PredefinedType: {predefined_type}")
    psetdata = ifcopenshell.util.element.get_psets(door)
    if 'Pset_DoorCommon' in psetdata:
        print('  Pset_DoorCommon')
        entry = psetdata['Pset_DoorCommon']
        outputElement(entry, 'IsExternal', 3)
        outputElement(entry, 'Reference', 3)
        outputElement(entry, 'EffectiveWidth', 3)
        outputElement(entry, 'RequiredDoorFlowrate', 3)


def get_basic_door_info(door):
    print("=========== IFCDoor ===========")

    door_info = door.get_info()
    name = get_property(door_info, 'Name')
    OverallWidth = get_property(door_info, "OverallWidth")
    OverallHeight = get_property(door_info, "OverallHeight")
    psetdata = ifcopenshell.util.element.get_psets(door)
    door_data = {}
    door_data['GlobalId'] = get_property(door_info, 'GlobalId')
    if name is not None:
        door_data['Name'] = name
    if OverallWidth is not  None:
        door_data['OverallWidth'] = OverallWidth

    if OverallHeight is not  None:
        door_data['OverallHeight'] = OverallHeight

    if 'Pset_DoorCommon' in psetdata:
        entry = psetdata['Pset_DoorCommon']
        isExternal = get_property(entry, "IsExternal") # In IDS
        aReference = get_property(entry, 'Reference')
        EffectiveWidth = get_property(entry, 'EffectiveWidth')
        FireRating = get_property(entry, 'FireRating')
        FireExit = get_property(entry, 'FireExit') # In IDS
        ThermalTransmittance = get_property(entry, 'ThermalTransmittance')
        SelfClosing = get_property(entry, 'SelfClosing') # In IDS
        HandicapAccessible = get_property(entry, 'HandicapAccessible') # In IDS
        Direction = get_property(entry, 'Direction')
        OperationType = get_property(entry, 'OperationType')
        UsagePreferenceFactor = get_property(entry, 'UsagePreferenceFactor')
        if isExternal is not  None:
            door_data['IsExternal'] = isExternal
        if EffectiveWidth is not  None:
            door_data['EffectiveWidth'] = EffectiveWidth
        if FireRating is not  None:
            door_data['FireRating'] = FireRating
        if FireExit is not  None:
            door_data['FireExit'] = FireExit
        if ThermalTransmittance is not  None:
            door_data['ThermalTransmittance'] = ThermalTransmittance
        if SelfClosing is not  None:
            door_data['SelfClosing'] = SelfClosing
        if HandicapAccessible is not  None:
            door_data['HandicapAccessible'] = HandicapAccessible
        if Direction is not  None:
            door_data['Direction'] = Direction
        if OperationType is not  None:
            door_data['OperationType'] = OperationType
        if UsagePreferenceFactor is not  None:
            door_data['UsagePreferenceFactor'] = UsagePreferenceFactor

    if 'ePset_OMA_Door' in psetdata:
        entry = psetdata['ePset_OMA_Door']
        door_data['DesignPeopleFlowRate'] = get_property(entry, 'DesignPeopleFlowRate') # in IDS
        door_data['EffectiveWidth'] = get_property(entry, 'EffectiveWidth') # in IDS
        door_data['DirectionOfTravel'] = get_property(entry, 'DirectionOfTravel') # in IDS

    return door_data


def get_basic_stair_info(stair, stair_dict):
    print("=========== IFCStair ===========")
    stair_info = stair.get_info()

    stair_dict['GlobalId'] = get_property(stair_info, 'GlobalId')
    stair_dict['Name'] = get_property(stair_info, 'Name')

    psetdata = ifcopenshell.util.element.get_psets(stair)

    if 'Pset_StairCommon' in psetdata:
        entry = psetdata['Pset_StairCommon']
        outputElement(entry, 'IsExternal', 3)
        outputElement(entry, 'Reference', 3)
        stair_dict['treads'] = get_property(entry, 'NumberOfTreads')  # In IDS
        stair_dict['TreadLength'] = get_property(entry, 'TreadLength')  # In IDS
        stair_dict['risers'] = get_property(entry, 'NumberOfRiser')  # In IDS
        stair_dict['RiserHeight'] = get_property(entry, 'RiserHeight')  # In IDS
        stair_dict['NosingLength'] = get_property(entry, 'NosingLength')  # In IDS
        outputElement(entry, 'TreadLengthAtOffset', 3)
        outputElement(entry, 'TreadLengthAtInnerSide', 3)

    if 'ePset_OMA_Stair' in psetdata:
        entry = psetdata['ePset_OMA_Stair']
        stair_dict['EffectiveWidth'] = get_property(entry, 'EffectiveWidth')  # In IDS
        stair_dict['OccupancyNumberPeak'] = get_property(entry, 'OccupancyNumberPeak')  # In IDS
        stair_dict['RequiredStairFlowrate'] = get_property(entry, 'RequiredStairFlowrate')  # In IDS
        stair_dict['TravelDirection'] = get_property(entry, 'TravelDirection')  # In IDS


def get_stair_flight_info(stair_flight, stair_flight_dict, stair_dict):
    """
        Gets stair flight information from the IFC file
        :param stair_flight: the stair_flight IFC information
        :param dictionary stair_flight_dict: the stair_flight dictionary to collect information on
        :param dictionary stair_dict: the parent stair dictionary used to complete missing information in the IFC file
    """
    element_info = stair_flight.get_info()
    stair_flight_dict['Name'] = get_property(element_info, 'Name')
    riser_number = get_property(element_info, 'NumberOfRiser')
    riser_height = get_property(element_info, 'RiserHeight')
    tread_number = get_property(element_info, 'NumberOfTreads')
    tread_length = get_property(element_info, 'TreadLength')
    overall_width = get_property(element_info, 'OverallWidth')
    tread_length = get_property(element_info, 'TreadLengthAtInnerSide')

    stair_flight_dict['GlobalId'] = get_property(element_info, 'GlobalId')
    psetdata = ifcopenshell.util.element.get_psets(stair_flight)
    if 'Pset_StairFlightCommon' in psetdata:
        entry = psetdata['Pset_StairFlightCommon']
        riser_height = get_property(entry, 'RiserHeight')  # In IDS
        riser_number = get_property(entry, 'NumberOfRiser')  # In IDS
        tread_length = get_property(entry, 'TreadLength')  # In IDS
        stair_flight_dict['NosingLength'] = get_property(entry, 'NosingLength')  # In IDS
        tread_number = get_property(entry, 'NumberOfTreads')  # In IDS
        stair_flight_dict['EffectiveWidth'] = get_property(entry, 'EffectiveWidth')  # In IDS

    if 'ePset_OMA_Stairs' in psetdata:
        entry = psetdata['ePset_OMA_Stairs']
        stair_flight_dict['OccupancyNumberPeak'] = get_property(entry, 'OccupancyNumberPeak')  # In IDS
        stair_flight_dict['RequiredStairFlowrat'] = get_property(entry, 'RequiredStairFlowrat')  # In IDS

    if 'ePset_OMA_StairFlight' in psetdata:
        entry = psetdata['ePset_OMA_StairFlight']
        stair_flight_dict['TravelDirection'] = get_property(entry, 'TravelDirection')  # In IDS

    if overall_width is not  None:
        print(f"  OverallWidth:{overall_width}")
        stair_flight_dict['Width'] = overall_width
        stair_flight_dict['Lanes'] = math.floor(overall_width / 0.7)

    if riser_number == None:
        if 'risers' in stair_dict:
            riser_number = stair_dict['risers']

    if riser_number is not  None:
        print(f"  NumberOfRiser:{riser_number}")
        stair_flight_dict['risers'] = riser_number

    if riser_height is None:
        if 'RiserHeight' in stair_dict:
            riser_height = stair_dict['RiserHeight']

    if riser_height is not  None:
        print(f"  RiserHeight:{riser_height}")
        stair_flight_dict['RiserHeight'] = riser_height

    if riser_number is not  None and riser_height != None:
        stair_flight_dict['Height'] = riser_number * riser_height

    if tread_number is None:
        if 'treads' in stair_dict:
            tread_number = stair_dict['treads']

    if tread_number is not None:
        print(f"  NumberOfTreads:{tread_number}")
        stair_flight_dict['treads'] = tread_number

    if tread_length is None:
        if 'TreadLength' in stair_dict:
            tread_length = stair_dict['TreadLength']

    if tread_length is not  None:
        print(f"  TreadLength:{tread_length}")
        stair_flight_dict['TreadLength'] = tread_length

    if tread_length is not  None and tread_number != None:
        stair_flight_dict['TotalRun'] = tread_number * tread_length


def get_basic_escalator_info(escalator_dict, ifc_transport):
    transport_info = ifc_transport.get_info()
    escalator_dict['Name'] = get_property(transport_info, 'Name')
    print('ESCALATOR ', escalator_dict['Name'])

    psetdata = ifcopenshell.util.element.get_psets(ifc_transport)
    if 'Pset_ElementKinematics' in psetdata:
        entry = psetdata['Pset_ElementKinematics']  # In IDS
        escalator_dict['MaximumConstantSpeed'] = get_property(entry, "MaximumConstantSpeed")  # In IDS

    if 'ePset_OMA_TransportElementCommon' in psetdata:
        entry = psetdata['ePset_OMA_TransportElementCommon']
        escalator_dict['TreadLength'] = get_property(entry, "TreadLength")  # In IDS
        escalator_dict['ClearWidth'] = get_property(entry, "ClearWidth")  # In IDS
        escalator_dict['LoadingArea'] = get_property(entry, "LoadingArea")  # In IDS
        escalator_dict['DesignPeopleFlowRate'] = get_property(entry, "DesignPeopleFlowRate")  # In IDS
        escalator_dict['CapacityPeople'] = get_property(entry, "CapacityPeople")  # In IDS
        escalator_dict['RiserHeight'] = get_property(entry, "RiserHeight")  # In IDS
        escalator_dict['RunLength'] = get_property(entry, "RunLength")  # In IDS
        escalator_dict['RunHeight'] = get_property(entry, "RunHeight")  # In IDS
        escalator_dict['OverallLength'] = get_property(entry, "OverallLength")  # Not in IDS

    # also need level run and direction of travel


def get_basic_movingwalkway_info(movingwalkway_dict, ifc_transport):
    transport_info = ifc_transport.get_info()
    movingwalkway_dict['Name'] = get_property(transport_info, 'Name')
    print('MOVINGWALKWAY ', movingwalkway_dict['Name'])
    psetdata = ifcopenshell.util.element.get_psets(ifc_transport)
    if 'Pset_TransportElementCommon' in psetdata:
        entry = psetdata['Pset_TransportElementCommon']
        movingwalkway_dict['MaximumConstantSpeed'] = get_property(entry, "MaximumConstantSpeed") # In IDS

    if 'ePset_OMA_TransportElementCommon' in psetdata:
        entry = psetdata['ePset_OMA_TransportElementCommon']
        movingwalkway_dict['CapacityPeople'] = get_property(entry, "CapacityPeople") # In IDS
        movingwalkway_dict['ClearWidth'] = get_property(entry, "ClearWidth") # In IDS
        movingwalkway_dict['LoadingArea'] = get_property(entry, "LoadingArea")  # In IDS
        movingwalkway_dict['DesignPeopleFlowRate'] = get_property(entry, "DesignPeopleFlowRate")  # In IDS
        movingwalkway_dict['RunLength'] = get_property(entry, "RunLength")  # In IDS
        movingwalkway_dict['OverallLength'] = get_property(entry, "OverallLength")  # Not in IDS

    # also need ldirection of travel


def get_basic_elevator_info(elevator_dict, ifc_transport):
    transport_info = ifc_transport.get_info()
    elevator_dict['Name'] = get_property(transport_info, 'Name')
    print('ELEVATOR ', elevator_dict['Name'])
    psetdata = ifcopenshell.util.element.get_psets(ifc_transport)
    if "Pset_TransportElementCommon" in psetdata:
        entry = psetdata['Pset_TransportElementCommon']
        elevator_dict['Reference'] = get_property(entry, "Reference")  # Not in IDS
        elevator_dict['CapacityPeople'] = get_property(entry, "CapacityPeople")  # In IDS

    if "Pset_TransportElementElevator" in psetdata:
        entry = psetdata['Pset_TransportElementElevator']
        elevator_dict['ClearWidth'] = get_property(entry, "ClearWidth")  # In IDS
        elevator_dict['ClearDepth'] = get_property(entry, "ClearDepth")  # In IDS

    if 'Pset_ElementKinematics' in psetdata:
        entry = psetdata['Pset_ElementKinematics']
        elevator_dict['MaximumConstantSpeed'] = get_property(entry, "MaximumConstantSpeed")  # In IDS

    if 'ePset_OMA_TransportElementElevator' in psetdata:
        entry = psetdata['ePset_OMA_TransportElementElevator']
        elevator_dict['AccelerationRate'] = get_property(entry, "AccelerationRate")  # In IDS
        elevator_dict['DecelerationRate'] = get_property(entry, "DecelerationRate")  # In IDS
        elevator_dict['DoorOpeningTime'] = get_property(entry, "DoorOpeningTime")  # In IDS
        elevator_dict['DoorClosingTime'] = get_property(entry, "DoorClosingTime")  # In IDS
        elevator_dict['DoorDwellTime'] = get_property(entry, "DoorDwellTime")  # In IDS
        elevator_dict['MotorDelayTime'] = get_property(entry, "MotorDelayTime")  # In IDS
        elevator_dict['StoriesServed'] = get_property(entry, "StoriesServed")  # In IDS
        elevator_dict['Bank'] = get_property(entry, "Group")  # In IDS
        elevator_dict['JerkRate'] = get_property(entry, "JerkRate")  # In IDS
        elevator_dict['LevellingDelayTime'] = get_property(entry, "LevellingDelayTime")  # In IDS, but change to Pset_OMA_ElementKinematics maybe
        elevator_dict['AdvancedDoorOpeningTime'] = get_property(entry, "AdvancedDoorOpeningTime")  # In IDS, but change to Pset_OMA_ElementKinematics maybe
        elevator_dict['LoadingArea'] = get_property(entry, "LoadingArea")  # In IDS
        
        elevator_dict['Description'] = get_property(entry, "Description")  # Not in IDS


def get_basic_ramp_info(ramp_dict, ifc_ramp):
    ramp_info = ifc_ramp.get_info()
    ramp_dict['Name'] = ramp_info['Name']

    psetdata = ifcopenshell.util.element.get_psets(ifc_ramp)
    if 'Pset_RampCommon' in psetdata:
        entry = psetdata['Pset_RampCommon']
        ramp_dict['HandicapAccessible'] = get_property(entry, "HandicapAccessible")  # In IDS
        ramp_dict['RequiredSlope'] = get_property(entry, "RequiredSlope")  # In IDS
        ramp_dict['FireExit'] = get_property(entry, "FireExit")  # In IDS

    if 'ePset_OMA_RampCommon' in psetdata:
        entry = psetdata['ePset_OMA_RampCommon']
        ramp_dict['DirectionOfTravel'] = get_property(entry, "HandicapAccessible")  # In IDS
    

def get_basic_space_info(ifc_space, space_data):
    element_info = ifc_space.get_info()
    name = get_property(element_info, 'Name')
    space_id = get_property(element_info, 'GlobalId')
    psetdata = ifcopenshell.util.element.get_psets(ifc_space)
    reference = None
    if 'Pset_SpaceCommon' in psetdata:
        entry = psetdata['Pset_SpaceCommon']
        reference = get_property(entry, 'Reference')

    space_data['GlobalId'] = get_property(element_info, 'GlobalId')
    if name is not None:
        space_data['Name'] = name
    if reference is not None:
        space_data['Reference'] = reference


def get_basic_subspace_info(ifc_space_boundary, suspace_data):
    element_info = ifc_space_boundary.get_info()
    name = get_property(element_info, 'Name')
    psetdata = ifcopenshell.util.element.get_psets(ifc_space_boundary)
    reference = None
    parent_id = None
    if 'Pset_SpaceCommon' in psetdata:
        entry = psetdata['Pset_SpaceCommon']
        reference = get_property(entry, 'Reference')

    aRelatingSpace = ifc_space_boundary.RelatingSpace
    if aRelatingSpace:
        parent_info = aRelatingSpace.get_info()
        parent_id = get_property(parent_info, 'GlobalId')
        psetdata = ifcopenshell.util.element.get_psets(aRelatingSpace)
        if 'Pset_SpaceCommon' in psetdata:
            entry = psetdata['Pset_SpaceCommon']
            reference = get_property(entry, 'Reference')

    suspace_data['GlobalId'] = get_property(element_info, 'GlobalId')
    if name is not None:
        suspace_data['Name'] = name
    if reference is not None:
        suspace_data['Reference'] = reference
    if parent_id is not None:
        suspace_data['ParentId'] = parent_id


def get_storey_Elevation_spaceboundary_quiet(ifc_element):
    aRelatingSpace = ifc_element.RelatingSpace  # RelatedBuildingElement
    Elevation = 0
    if aRelatingSpace:
        a_storey = aRelatingSpace.Decomposes[0].RelatingObject
        element_info = a_storey.get_info()
        Elevation = get_property(element_info,
                                 'Elevation')  # moving to ElevationOfFFLRelative, under Pset_BuildingStoreyCommon
    return Elevation


def get_storey_Longname_spaceboundary(ifc_element):
    aRelatingSpace = ifc_element.RelatingSpace  # RelatedBuildingElement
    long_name = None
    if aRelatingSpace:
        a_storey = aRelatingSpace.Decomposes[0].RelatingObject
        element_info = a_storey.get_info()
        long_name = get_property(element_info, 'LongName')
    return long_name


def get_storey_Elevation_spaceboundary(ifc_element):
    aRelatingSpace = ifc_element.RelatingSpace  # RelatedBuildingElement
    Elevation = 0
    if aRelatingSpace:
        a_storey = aRelatingSpace.Decomposes[0].RelatingObject
        element_info = a_storey.get_info()
        print('  Level ', end='')
        outputElement(element_info, 'Name', 0)
        print('  Level ', end='')
        outputElement(element_info, 'Elevation', 0)
        print('  Level ', end='')
        outputElement(element_info, 'LongName', 0)
        print('  Level ', end='')
        outputElement(element_info, 'GlobalId', 0)
        Elevation = get_property(element_info, 'Elevation')

        element_info = aRelatingSpace.get_info()
        print('  Space ', end='')
        outputElement(element_info, 'GlobalId', 0)
        print('  Space ', end='')
        outputElement(element_info, 'Name', 0)
        print('  Space ', end='')
        outputElement(element_info, 'LongName', 0)
    return Elevation


def get_storey_GlobalId(ifc_element):
    a_relating_space = parent(ifc_element)
    global_id = None
    if a_relating_space:
        element_info = a_relating_space.get_info()
        global_id = get_property(element_info, 'GlobalId')
    return global_id


def get_storey_Elevation(ifc_element):
    aRelatingSpace = parent(ifc_element)
    Elevation = 0
    if aRelatingSpace:
        # print("Has Parent")
        element_info = aRelatingSpace.get_info()
        # outputElement(element_info,'GlobalId',0)
        # outputElement(element_info,'Name',0)
        Elevation = get_property(element_info, 'Elevation')

    return Elevation


def get_storey_LongName(ifc_element):
    a_relating_space = parent(ifc_element)
    LongName = None
    if a_relating_space:
        element_info = a_relating_space.get_info()
        LongName = get_property(element_info, 'LongName')
    return LongName


def get_storey_Name(ifc_element):
    a_relating_space = parent(ifc_element)
    Name = None
    if a_relating_space:
        element_info = a_relating_space.get_info()
        Name = get_property(element_info, 'Name')
    return Name


def get_storey_Description(ifc_element):
    a_relating_space = parent(ifc_element)
    Description = None
    if a_relating_space:
        element_info = a_relating_space.get_info()
        Description = get_property(element_info, 'Description')
    return Description


def get_space_boundary_elem(ifc_space):
    boundary_elem = []
    print("Scanning Boundaries================================")
    for boundary in ifc_space.BoundedBy:
        if boundary.is_a("IfcRelSpaceBoundary"):
            if boundary.PhysicalOrVirtualBoundary == "VIRTUAL":
                relspace_info = boundary.get_info()
                aGlobalId = get_property(relspace_info, 'GlobalId')
                boundary_elem.append(["IfcRelSpaceBoundaryV", aGlobalId])
            elem = boundary.RelatedBuildingElement
            if elem:
                if elem.is_a("IfcDoor"):
                    door_info = elem.get_info()
                    aGlobalId = get_property(door_info, 'GlobalId')
                    boundary_elem.append(["IfcDoor", aGlobalId])
                elif elem.is_a("IfcStair"):
                    stair_info = elem.get_info()
                    aGlobalId = get_property(stair_info, 'GlobalId')
                    boundary_elem.append(["IfcStair", aGlobalId])
                elif elem.is_a("IfcStairFlight"):
                    stairflight_info = elem.get_info()
                    aGlobalId = get_property(stairflight_info, 'GlobalId')
                    boundary_elem.append(["IfcStairFlight", aGlobalId])
                elif elem.is_a("IfcSpace"):
                    space_info = elem.get_info()
                    aGlobalId = get_property(space_info, 'GlobalId')
                    boundary_elem.append(["IfcSpace", aGlobalId])
                elif elem.is_a("IfcRelSpaceBoundary"):
                    relspace_info = elem.get_info()
                    aGlobalId = get_property(relspace_info, 'GlobalId')
                    boundary_elem.append(["IfcRelSpaceBoundary", aGlobalId])
                elif elem.is_a("IfcVirtualElement"):
                    relspace_info = elem.get_info()
                    aGlobalId = get_property(relspace_info, 'GlobalId')
                    boundary_elem.append(["IfcVirtualElement", aGlobalId])

    return boundary_elem


def get_relspace_boundary_elem(ifc_relspace):
    a_relating_space = ifc_relspace.RelatingSpace  # RelatedBuildingElement
    if a_relating_space:
        return get_space_boundary_elem(a_relating_space)
    return None
