import math


def OutputMTAFileHeader(exodusmtafile):
    print("ZoomFactor:1", file=exodusmtafile)  # global default zoom factor


def OutputWindowFloorData(aWinID, floor_name, exodusmtafile):
    print("ModelWindow", file=exodusmtafile)  # Data for a new model window
    print("WindowID:", aWinID, file=exodusmtafile)
    print("ZFactor:1", file=exodusmtafile)  # level's zoom factor
    print("WinName:", floor_name, file=exodusmtafile)
    print("WinFloorHeight:3", file=exodusmtafile)
    # print("FloorNodeSize:1145", file=exodusmtafile)
    # print("FloorArcSize:1858", file=exodusmtafile)
    print("DisplayWindow: 1", file=exodusmtafile)


def OutputTransitArcData(stairflight, exodusmtafile):
    TopNodeList = None
    if 'ConxTop' in stairflight:
        TopNodeList = stairflight['ConxTop']
    LowerNodeList = None
    if 'ConxLower' in stairflight:
        LowerNodeList = stairflight['ConxLower']
    if TopNodeList or LowerNodeList:
        count = 0
        top_count = 0
        lower_count = 0
        if TopNodeList:
            top_count = len(TopNodeList)
            count += top_count
        if LowerNodeList:
            lower_count = len(LowerNodeList)
            count += lower_count
        if count > 0:
            LaneCount = stairflight['Lanes']
            RiserCount = stairflight['risers'] - 2
            # print("ArcLanesSize:4 -5 0 874 5 0 875 5 0 876 15 0 877")
            print("ArcLanesSize: ", count, file=exodusmtafile, end=' ')

            if lower_count == LaneCount:
                start_pos = 0
            else:
                start_pos = -5
            for Node in LowerNodeList:
                print(start_pos, 0, Node[0], file=exodusmtafile, end=' ')
                start_pos += 5

            if top_count == LaneCount:
                start_pos = 0
            else:
                start_pos = -5
            for Node in TopNodeList:
                print(start_pos, RiserCount, Node[0], file=exodusmtafile, end=' ')
                start_pos += 5

            print("", file=exodusmtafile)


def OutputTransitNodeData(transitNode, TransitID, aWinID, terrain, exodusmtafile):
    """
    terrain 0 , 2 escalator, 8 lift/elevator
    """
    Width = transitNode['Width']
    if 'Lanes' in transitNode and transitNode['Lanes'] is not None:
        Lanes = max(1, transitNode['Lanes'])
    else:
        Lanes = math.floor(width/0.5)

    print("TransitNode:", file=exodusmtafile)
    print("Label:", TransitID, file=exodusmtafile)
    print("WinID:", aWinID, file=exodusmtafile)
    print("TransitNodeLabel:", transitNode['nodeid'], file=exodusmtafile)  # Associate node
    print("IconType:1", file=exodusmtafile)  # how's displayed in EXODUS
    print("Terrain:", terrain, file=exodusmtafile)  # 0 staircase
    if 'risers' in transitNode:
        print("NumberOfSubUnits:", transitNode['risers'] - 1, file=exodusmtafile)  # 'risers'-1
        UnitLength = math.sqrt(math.pow(transitNode['riserheight'], 2) + math.pow(transitNode['treadlength'], 2))
        print("UnitLengthSize:", UnitLength, file=exodusmtafile)  # diagonal tread distance 'riserheight' 'treadlength'
        print("TheDirectionAngle:", transitNode['Direction'], file=exodusmtafile)  # Direction of icon on screen
        print("TheMaxCapacity:", (transitNode['risers'] - 1) * 2, file=exodusmtafile)  # 'risers'-1
        print("LaneMaxCapacity:", transitNode['risers'] - 1, file=exodusmtafile)  # 'risers'-1
        print("TheDownMaxCapacity:", (transitNode['risers'] - 1) * 2,
              file=exodusmtafile)  # Lanes times (number of risers -1) 'risers'-1
        print("LaneDownMaxCapacity:", transitNode['risers'] - 1, file=exodusmtafile)  # 'risers'-1
    else:
        # for lifts we have transitNode['ClearWidth'] and transitNode['ClearDepth']
        # however at the moment not sure whether they are required
        UnitLengthSize = 0.5
        NumberOfUnits = 1
        TravelDist = 0.5
        if 'HorizontalLength' in transitNode:
            TravelDist = transitNode['HorizontalLength']
            NumberOfUnits = TravelDist / UnitLengthSize

        Capacity = 1
        if 'CapacityPeople' in transitNode and transitNode['CapacityPeople'] is not None:
            Capacity = transitNode['CapacityPeople']
        else:
            Capacity = NumberOfUnits * Lanes

        LaneCapacity = Capacity
        if Lanes > 1:
            LaneCapacity = LaneCapacity / Lanes

        print("NumberOfSubUnits:", NumberOfUnits, file=exodusmtafile)  # 'risers'-1
        print("UnitLengthSize:", UnitLengthSize,
              file=exodusmtafile)  # diagonal tread distance 'riserheight' 'treadlength'
        print("TheDirectionAngle:", transitNode['Direction'], file=exodusmtafile)  # Direction of icon on screen
        print("TheMaxCapacity:", Capacity, file=exodusmtafile)  # 'risers'-1
        print("LaneMaxCapacity:", LaneCapacity, file=exodusmtafile)  # 'risers'-1
        print("TheDownMaxCapacity:", Capacity, file=exodusmtafile)  # Lanes times (number of risers -1) 'risers'-1
        print("LaneDownMaxCapacity:", Capacity, file=exodusmtafile)  # 'risers'-1

    print("DefaultRates:1", file=exodusmtafile)
    if 'TotalRun' in transitNode:
        TravelDist = math.sqrt(math.pow(transitNode['TotalRun'], 2) + math.pow(transitNode['Height'], 2))
        print("Length:", TravelDist, file=exodusmtafile)  # Length Diagonal 'Length' 'Height'
    elif 'HorizontalLength' in transitNode:
        TravelDist = transitNode['HorizontalLength']
        print("Length:", TravelDist, file=exodusmtafile)  # Length Diagonal 'Length' 'Height'
    else:
        print("Length:", 0.5, file=exodusmtafile)

    print("Width:", Width, file=exodusmtafile)  # Width
    print("NumberOfLanes:", Lanes, file=exodusmtafile)  # Lanes
    print("TransitNodeAngle:", transitNode['Direction'],
          file=exodusmtafile)  # Direction related to movement (often same as TheDirectionAngle:)

    if 'riserheight' in transitNode and terrain == 0:
        # for stairs we remove top and bottom risers
        print("Height:", transitNode['Height'] - 2 * transitNode['riserheight'],
              file=exodusmtafile)  # 'Height' - 2*'riserheight' vertical Height, not including top and bottom risers
    else:
        print("Height:", transitNode['Height'], file=exodusmtafile)

    if 'noising' in transitNode:
        print("NosingDepth:", transitNode['noising'], file=exodusmtafile)  # 'noising'

    print("EWidth:0", file=exodusmtafile)
    print("UseEWidth:0", file=exodusmtafile)
    print("HDWSize:0", file=exodusmtafile)

    if 'MaximumConstantSpeed' in transitNode:
        print("HozSpeed", transitNode['MaximumConstantSpeed'], file=exodusmtafile)

    if 'LevelRun' in transitNode:
        print("LevelRun", transitNode['LevelRun'], file=exodusmtafile)

    OutputTransitArcData(transitNode, exodusmtafile)
    print("CoarseNodeType:80", file=exodusmtafile)


def OutputLiftData(LiftData, LiftID, LiftNodeList, exodusmtafile):
    print("LiftTitle:", LiftData["Name"], file=exodusmtafile)
    print("LiftBankID:", LiftID, file=exodusmtafile)  # LiftData['Bank']
    start_floor = 0
    if 'FloorList' in LiftData:
        start_floor = LiftData['FloorList'][0]
    print("LiftCurrentFloor:", start_floor, file=exodusmtafile)
    print("LiftCurrentCapacity:", LiftData['CapacityPeople'], file=exodusmtafile)
    print("LiftStatus:38", file=exodusmtafile)
    print("LiftAcceleration:", LiftData['AccelerationRate'], file=exodusmtafile)
    print("LiftDeceleration:", LiftData['DecelerationRate'], file=exodusmtafile)
    print("LiftConstantSpeed:", LiftData['MaximumConstantSpeed'], file=exodusmtafile)
    print("LiftOpeningTime:", LiftData['DoorOpeningTime'], file=exodusmtafile)
    print("LiftClosingTime:", LiftData['DoorClosingTime'], file=exodusmtafile)
    print("LiftDwellTime:", LiftData['DoorDwellTime'], file=exodusmtafile)
    print("LiftMotorDelay:", LiftData['MotorDelayTime'], file=exodusmtafile)
    print("LiftJerkRate:0", file=exodusmtafile)
    print("LiftSensorBreakDwellDelay:", LiftData['DoorDwellTime'], file=exodusmtafile)
    print("LiftStartDelay:0", file=exodusmtafile)
    # LiftData['LoadingArea']
    if LiftNodeList is not None:
        print("LiftNodeList:", file=exodusmtafile, end=' ')
        for NodeData in LiftNodeList:
            print(NodeData[0], file=exodusmtafile, end=' ')
        print('', file=exodusmtafile)


def OutputNodeLoc(aNode, aType, aWinID, exodusmtafile):
    print("Node", file=exodusmtafile)
    print("Label:", aNode[0], file=exodusmtafile)
    if len(aNode) > 4:
        title_str = aNode[4].replace(" ", "_")
        print("Title:", title_str, file=exodusmtafile)
    print("WinID:", aWinID, file=exodusmtafile)
    print("PosMetresX:", aNode[1], file=exodusmtafile)
    print("PosMetresY:", aNode[2], file=exodusmtafile)
    print("NodeHeight:", aNode[3], file=exodusmtafile)
    if aType == 1:
        print("Nodetype:FreeSpace", file=exodusmtafile)
    elif aType == 6:
        print("Nodetype:Landing", file=exodusmtafile)
    elif aType == 80:
        print("Nodetype:MegaNode", file=exodusmtafile)
    elif aType == 100:
        print("Nodetype:Door", file=exodusmtafile)
    print("Obstacle:0", file=exodusmtafile)


def OutputExternalDoor(aDoor, aWinID, exodusmtafile):
    print("Door", file=exodusmtafile)
    print("Label:", aDoor['NodeID'], file=exodusmtafile)
    if 'Name' in aDoor:
        title_str = aDoor['Name'].replace(" ", "_")
        print("Title:", title_str, file=exodusmtafile)
    print("WinID:", aWinID, file=exodusmtafile)
    print("PosMetresX:", aDoor['Location'][0], file=exodusmtafile)
    print("PosMetresY:", aDoor['Location'][1], file=exodusmtafile)
    print("NodeHeight:0", file=exodusmtafile)
    print("Obstacle:0", file=exodusmtafile)
    print("Open:1", file=exodusmtafile)
    print("Active:1", file=exodusmtafile)
    print("DoorType:0", file=exodusmtafile)
    print("DoorDelay:1.33 1.33", file=exodusmtafile)
    print("TimeModel:0", file=exodusmtafile)
    print("Attractiveness:100", file=exodusmtafile)
    print("CatchmentArea:7.5", file=exodusmtafile)
    print("Capacity:999999999", file=exodusmtafile)


def OutputLineFileHeader(NumberOfFloors, exodusegxfile):
    print("File Version:1", file=exodusegxfile)
    print("Floor Number:", NumberOfFloors, file=exodusegxfile)
    print("Fixed Flag:1", file=exodusegxfile)
    print("Line Option:1 1 1 1 1 1 1 1 1 ", file=exodusegxfile)
    print("Node Option:0 1 1 0 1 1 0 0 0 ", file=exodusegxfile)


def OutputLineSectionStart(LineCount, exodusegxfile):
    print("DXFFile Location:.", file=exodusegxfile)
    print("DXFFile Name:", file=exodusegxfile)
    print("Fixed Lines:0", file=exodusegxfile)
    print("User Lines:0", file=exodusegxfile)
    print("Additional Lines:", LineCount, file=exodusegxfile)


def OutputLineData(line, height, exodusegxfile):
    # X1,Y1,X2,Y2,Height,32bitColour
    if line[4]:
        print("Line:", line[0], " ", line[1], " ", line[2], " ", line[3], " ", height, " ", 16777216,
              file=exodusegxfile)


def OutputLineSectionEnd(exodusegxfile):
    print("Boundary Polygons:0", file=exodusegxfile)
    print("Convex Polygons:0", file=exodusegxfile)
    print("Screen Labels:0", file=exodusegxfile)


def OutputLineTail(exodusegxfile):
    print("NumberOfZones:0", file=exodusegxfile)


def AdjustLineLoc(RangeX, RangeY, offsetX, offsetY, lines):
    for line in lines:
        line[0] = offsetX + line[0]
        line[1] = offsetY + ((RangeY[1] - RangeY[0]) - line[1])
        line[2] = offsetX + line[2]
        line[3] = offsetY + ((RangeY[1] - RangeY[0]) - line[3])


def AdjustNodeLoc(RangeX, RangeY, offsetX, offsetY, room_nodes):
    for room in room_nodes:
        for aRow in room:
            for aNode in aRow:
                if aNode[0] > -1:
                    aNode[1] = offsetX + aNode[1]
                    aNode[2] = offsetY + ((RangeY[1] - RangeY[0]) - aNode[2])


def AdjustDoorsLoc(RangeX, RangeY, offsetX, offsetY, Doors):
    for aDoor in Doors:
        if 'Location' in aDoor:
            aDoor['Location'][0] = offsetX + aDoor['Location'][0]
            aDoor['Location'][1] = offsetY + ((RangeY[1] - RangeY[0]) - aDoor['Location'][1])


def AdjustTransitNodeLoc(RangeX, RangeY, offsetX, offsetY, transit_nodes):
    if transit_nodes is not None:
        for transit_node in transit_nodes:
            print("Transit Adjust:", transit_node['nodedata'])
            transit_node['nodedata'][1] = offsetX + transit_node['nodedata'][1]
            transit_node['nodedata'][2] = offsetY + ((RangeY[1] - RangeY[0]) - transit_node['nodedata'][2])
            # since changing Y direction need to update rotation attribute too
            if 'Direction' in transit_node:
                aDirection = transit_node['Direction']
                if aDirection > 0.0:
                    aDirection = 360.0 - aDirection
                elif aDirection < 0.0:
                    aDirection = aDirection + 360.0
                transit_node['Direction'] = aDirection


def find_door_connections(NodeID, nodal_door_connections):
    '''
    :return: a list of external doors this node is connected to
    '''
    door_connections = []
    if nodal_door_connections:
        for iConnection in nodal_door_connections:
            if iConnection[0] == NodeID:
                door_connections.append(iConnection[1])
            elif iConnection[1] == NodeID:
                door_connections.append(iConnection[0])
    return door_connections


def output_room_nodes(room_node_list, WinID, nodal_door_connections, node_stairflight_connections, exodusmtafile):
    """
    Outputs the nodes in a given room
    :param room_node_list: list of nodes in room
    :param nodal_door_connections: nodal connections between doors
    :param WinID: - floor ID
    :param exodusmtafile: MTA file
    """
    aRowID = 0
    node_count = len(room_node_list)
    # print("New node set ============================",node_count)
    for aRow in room_node_list:
        row_length = len(aRow)
        # print("New row ============================",row_length)
        aColID = 0
        for aNode in aRow:
            if aNode[0] > -1:
                OutputNodeLoc(aNode, 1, WinID, exodusmtafile)

                node_connections = find_door_connections(aNode[0], nodal_door_connections)
                stair_connections = find_door_connections(aNode[0], node_stairflight_connections)
                if stair_connections is not None and len(stair_connections) > 0:
                    node_connections.extend(stair_connections)

                Loc = aColID - 1
                leftNode = None
                if Loc > -1:
                    if aRow[Loc][0] > -1:
                        node_connections.append(aRow[Loc][0])

                Loc = aColID + 1
                rightNode = None
                if Loc < row_length:
                    if aRow[Loc][0] > -1:
                        node_connections.append(aRow[Loc][0])

                Loc = aRowID - 1
                aboveNode = None
                if Loc > -1:
                    if room_node_list[Loc][aColID][0] > -1:
                        node_connections.append(room_node_list[Loc][aColID][0])

                Loc = aRowID + 1
                belowNode = None
                if Loc < node_count:
                    if room_node_list[Loc][aColID][0] > -1:
                        node_connections.append(room_node_list[Loc][aColID][0])

                if len(node_connections) > 0:
                    print("adjacentNodes:", file=exodusmtafile, end='')
                    for NodeID in node_connections:
                        print(f" {NodeID} 0.5 0", file=exodusmtafile, end='')

                    print("", file=exodusmtafile)
            aColID += 1
        aRowID += 1


def OutputTransitNodeConnections(stairflight, exodusmtafile):
    TopNodeList = None
    if 'ConxTop' in stairflight:
        TopNodeList = stairflight['ConxTop']
    LowerNodeList = None
    if 'ConxLower' in stairflight:
        LowerNodeList = stairflight['ConxLower']
    if (TopNodeList or LowerNodeList):
        count = 0
        if TopNodeList:
            count += len(TopNodeList)
        if LowerNodeList:
            count += len(LowerNodeList)
        if count > 0:
            print("adjacentNodes: ", file=exodusmtafile, end=' ')

            for Node in LowerNodeList:
                print(f"{Node[0]} 0.5 0", file=exodusmtafile, end=' ')

            for Node in TopNodeList:
                print(f"{Node[0]} 0.5 0", file=exodusmtafile, end=' ')

            print("", file=exodusmtafile)


def output_landing_nodes(landing_node_list, height, WinID, exodusmtafile):
    """
    Outputs the nodes in a given landing
    param: landing_node_list - list of nodes on a landing
    param: height - height of landing from current level
    param: WinID - floor ID
    param: exodusmtafile  - MTA file
    """
    aRowID = 0
    node_count = len(landing_node_list)
    # print("New node set ============================",node_count)
    for aRow in landing_node_list:
        row_length = len(aRow)
        # print("New row ============================",row_length)
        aColID = 0
        for aNode in aRow:
            if aNode[0] > -1:
                # print ("Node", aColID, aRowID, aNode)
                if height > 0.0:
                    aNode[3] += height
                OutputNodeLoc(aNode, 6, WinID, exodusmtafile)

                node_connections = []

                Loc = aColID - 1
                leftNode = None
                if Loc > -1:
                    # print("left:",Loc, aRow[Loc])
                    if aRow[Loc][0] > -1:
                        node_connections.append(aRow[Loc][0])

                Loc = aColID + 1
                rightNode = None
                if Loc < row_length:
                    # print("right:", Loc, aRow[Loc])
                    if aRow[Loc][0] > -1:
                        node_connections.append(aRow[Loc][0])

                Loc = aRowID - 1
                aboveNode = None
                if Loc > -1:
                    # print("above:",Loc, aColID , landing_node_list[Loc][aColID])
                    if landing_node_list[Loc][aColID][0] > -1:
                        node_connections.append(landing_node_list[Loc][aColID][0])

                Loc = aRowID + 1
                belowNode = None
                if Loc < node_count:
                    # print("below:",Loc, aColID , landing_node_list[Loc][aColID])
                    if landing_node_list[Loc][aColID][0] > -1:
                        node_connections.append(landing_node_list[Loc][aColID][0])

                if len(node_connections) > 0:
                    print("adjacentNodes:", file=exodusmtafile, end='')
                    for NodeID in node_connections:
                        print(f" {NodeID} 0.5 0", file=exodusmtafile, end='')

                    print("", file=exodusmtafile)
            # else:
            #    print("Not valid ",aColID,aNode)
            aColID += 1
        aRowID += 1


def output_external_door(door_data, WinID, nodal_door_connections, exodusmtafile):
    """
    Outputs an external doors data
    """
    OutputExternalDoor(door_data, WinID, exodusmtafile)
    node_connections = find_door_connections(door_data['NodeID'], nodal_door_connections)
    if len(node_connections) > 0:
        print("adjacentNodes:", file=exodusmtafile, end='')
        for NodeID in node_connections:
            print(f" {NodeID} 0.0 0", file=exodusmtafile, end='')
        print("", file=exodusmtafile)
