# -*- coding: utf-8 -*-
"""
Created on Thu Dec 8 09:10:58 2023

@author: P.J.Lawrence
"""
import math
from typing import Any


def TriangleOrientation(Xa, Ya, Xb, Yb, Xc, Yc):
    ZERO_TOLERANCE2 = 0.0001
    Det = (Xb - Xa)*(Yc - Ya) - (Xc - Xa)*(Yb - Ya)
    if Det < -ZERO_TOLERANCE2:
        return -1
    if Det > ZERO_TOLERANCE2:
        return 1
    return 0


def intri(x1, y1, x2, y2, x3, y3, tpx, tpy):
    dist1 = PointSide(x1, y1, x2, y2, x3, y3)
    dist2 = PointSide(x1, y1, x2, y2, tpx, tpy)
    if dist1 * dist2 <= 0:
        return False

    dist1 = PointSide(x2, y2, x3, y3, x1, y1)
    dist2 = PointSide(x2, y2, x3, y3, tpx, tpy)
    if dist1 * dist2 <= 0:
        return False

    dist1 = PointSide(x3, y3, x1, y1, x2, y2)
    dist2 = PointSide(x3, y3, x1, y1, tpx, tpy)
    if dist1 * dist2 <= 0:
        return False
    return True


def PointInPolygon(X, Y, Region):
    '''
    Determine whether a point lies inside, outside, or on
    the boundary of a  polygon Region in 2 dimensional space.

    returns:-
      -2 fails
      -1 outside
       0 on edge
       1 inside
    '''
    DTOL = 0.00001
    '''
      Direction of ray is from X,Y through midpoint of first edge
      such that X,Y is not collinear with edge.
    
      DTol the distance tolerance for a point being on an edge of
      the polygon.
    
      Note NORMAL has 3 points because the direction of the polygon
      has to be established ie CW or CCW.
    '''
    
    INOUT = -2
    Direction = [0,0]
    Normal = [0,0,0]
    DE = [0,0]
    RHS = [0,0]
    
    i = 0
    LB = 0
    NumPoints = len(Region)
    while True:
        while True:
            i = i + 1
            if i>=NumPoints:
                return INOUT
            LA = LB
            LB = i
            aVA=Region[LA]
            aVB=Region[LB]
            DE[0] = aVB[0] - aVA[0]
            Direction[0] = X - aVA[0]
            DE[1] = aVB[1] - aVA[1]
            Direction[1] = Y - aVA[1]
            
            Len1 = DE[0]*DE[0]  + DE[1]*DE[1]
            Len2 = Direction[0]*Direction[0] + Direction[1]*Direction[1]
            if Len1 != 0.0:
                break
        
        if Len2 == 0.0:
            INOUT = 0
            return INOUT
        
        dotp = abs(DE[0]*Direction[0] + DE[1]*Direction[1])/math.sqrt(Len1*Len2)
        if dotp <= (1.0 - DTOL):
            break
    
    
    aVA=Region[LA]
    aVB=Region[LB]
    Direction[0] = 0.5*(aVA[0] + aVB[0]) - X
    Direction[1] = 0.5*(aVA[1] + aVB[1]) - Y
    
    Dist = math.sqrt(Direction[0]*Direction[0] + Direction[1]*Direction[1])
    Direction[0] = Direction[0]/Dist
    Direction[1] = Direction[1]/Dist
    
    Normal[0] = -Direction[1]
    Normal[1] =  Direction[0]
    Normal[2] = Normal[0]*X + Normal[1]*Y
    
    Dist = Normal[0]*aVB[0] + Normal[1]*aVB[1] - Normal[2]
    
    SB = 0
    if Dist > 0.0:
        SB = 1
    else:
        SB = -1
    
    M = 0
    if abs(Direction[1]) > abs(Direction[0]):
        M = 1
    
    K = 1
    '''
      For remaining edges of polygon, check whether ray intersects edge.
      Vertices or edges lying on ray are handled by looking at preceding
      and succeeding edges not lying on ray.
    '''
    
    N = i
    i = i + 1
    if i==NumPoints:
        i=0
    
    while True:
        LA = LB
        LB = i
        SA = SB
        aVB=Region[LB]
        Dist = Normal[0]*aVB[0] + Normal[1]*aVB[1] - Normal[2]
        
        if abs(Dist) <= DTOL:
            SB = 0
        elif Dist  > 0.0:
            SB = 1
        else:
            SB = -1
            
        S = SA*SB
        if S < 0:
            aVA=Region[LA]
            DE[0]  = aVA[0] - aVB[0]
            DE[1]  = aVA[1] - aVB[1]
            RHS[0] = aVA[0] - X
            RHS[1] = aVA[1] - Y
            
            T = (RHS[0]*DE[1] - RHS[1]*DE[0])/(Direction[0]*DE[1] - Direction[1]*DE[0])
            
            if T > DTOL:
                K = K + 1
            elif T>=-DTOL:
                INOUT = 0
                return INOUT
        elif S==0:
            L = LB
            while True:
                i = i + 1
                if i >= NumPoints:
                    i = 0
                if i == N:
                    return INOUT
                LA = LB
                LB = i
                aVB=Region[LB]
                Dist = Normal[0]*aVB[0] + Normal[1]*aVB[1] - Normal[2]
                if abs(Dist) >= DTOL:
                    break
            
            if Dist > 0.0:
                SB = 1
            else:
                SB = -1
            
            aL=Region[L]
            if M==0:
                T = (aL[0] - X)/Direction[M]
            else:
                T = (aL[1] - Y)/Direction[M]
            
            if abs(T) < DTOL:
                INOUT = 0
                return INOUT
            
            if LA!=L:
                aVA=Region[LA]
                if M==0:
                    TA = (aVA[0] - X)/Direction[M]
                else:
                    TA = (aVA[1] - Y)/Direction[M]
                if abs(TA) <= DTOL or T*TA < 0.0:
                    INOUT = 0
                    return INOUT
            
            if SA*SB < 0 and T > 0.0:
                K = K + 1
    
        i = i + 1
        if i >= NumPoints:
            i = 0
        if i == N:
            break
    '''
    Point lies inside polygon if number of intersections K is odd.
    '''
    if (K % 2) == 1:
        INOUT = 1
    else:
        INOUT = -1
    return INOUT


def PointSide(XS, YS, XE, YE, X, Y):
    Xlk = X - XE
    Ylk = Y - YE
    Xkj = XE - XS
    Ykj = YE - YS
    area = Xkj*Ylk - Xlk*Ykj
    
    if area > -0.001 and area < 0.001:
        return 0 # on edge
    if area < 0.0:
        return -1
    return 1


def same_point(x1,x2,y1,y2):
    thetol = 0.001
    if abs(x1-x2)<thetol:
        if abs(y1-y2)<thetol:
            return True
    return False

 
def LineIntersect( SL1,EL1, SL2, EL2  ):
    LX1 =  SL1[0]
    LY1 =  SL1[1]

    F1 = EL1[0] - LX1
    G1 = EL1[1] - LY1

    F2 = EL2[0] - SL2[0]
    G2 = EL2[1] - SL2[1]

    F1G2 = F1 * G2
    F2G1 = F2 * G1

    DET = F2G1 - F1G2
    if abs(DET) < 0.0001:
        return False

    X21 = SL2[0] - LX1
    Y21 = SL2[1] - LY1

    DET = 1.0 / DET

    S = (F2 * Y21 - G2 * X21) * DET
    T = (F1 * Y21 - G1 * X21) * DET
    Tolerance = 0.0001
    
    return T > -Tolerance and T < 1.0+Tolerance and S > -Tolerance and S < 1.0+Tolerance


def calcAngle(lineA, lineB):
    line1Y1 = lineA[0][1]
    line1X1 = lineA[0][0]
    line1Y2 = lineA[1][1]
    line1X2 = lineA[1][0]

    line2Y1 = lineB[0][1]
    line2X1 = lineB[0][0]
    line2Y2 = lineB[1][1]
    line2X2 = lineB[1][0]

    #calculate angle between pairs of lines
    angle1 = math.atan2(line1Y1-line1Y2,line1X1-line1X2)
    angle2 = math.atan2(line2Y1-line2Y2,line2X1-line2X2)
    angleDegrees = (angle1-angle2) * 360 / (2*math.pi)
    return angleDegrees


def CalcPosition(Point1, Point2, GapLength, expand_size):
    F = (Point2[0] - Point1[0])
    G = (Point2[1] - Point1[1])
    t = 1.0/GapLength*expand_size
    X = Point1[0] + t * F
    Y = Point1[1] + t * G
    return X, Y


def CalcPositionSides(Point1, Point2, GapLength, expand_size):
    F = (Point2[0] - Point1[0])
    G = (Point2[1] - Point1[1])
    print("GapLength, expand_size",GapLength,expand_size)
    t = 1.0/GapLength*expand_size
    X1 = Point1[0] + t * G
    Y1 = Point1[1] - t * F
    X2 = Point1[0] - t * G
    Y2 = Point1[1] + t * F
    return [X1,Y1],[X2,Y2]


def calcDistance3D(Point1, Point2):
    diffx = Point1[0]-Point2[0]
    diffy = Point1[1]-Point2[1]
    diffz = Point1[2]-Point2[2]
    length = math.sqrt(diffx*diffx + diffy*diffy + diffz*diffz)
    return length


def calcDistance2D(Point1, Point2):
    diffx = Point1[0]-Point2[0]
    diffy = Point1[1]-Point2[1]
    length = math.sqrt(diffx*diffx + diffy*diffy)
    return length


def SameValue(v1, v2):
    return abs(v1-v2)<0.06


def order_boundary_list(boundary_points):
    """
        orders a set of boundary_points edge lists
        for example {e1,e2,e2,e3,e4,e5,e3,e4,e5,e1} to get {e1,e2,e3,e4,e5}
        returns ordered list 
        :param boundary_points : list of point edges on boundary
    """
    ec1 = 0
    en1 = 1
    p1 = boundary_points[ec1]
    p2 = boundary_points[en1]
    boundary_ordered_poly = [p1,p2]
    max_loop = len(boundary_points)

    index_list = []
    for i in range(2, max_loop):
        index_list.append(i)

    #print('boundary_points ',boundary_points)
    p1 = boundary_ordered_poly[0]
    for i in range(0, max_loop):
        ec1 = len(boundary_ordered_poly)-2
        en1 = len(boundary_ordered_poly)-1
        #print('boundary_ordered_poly ',boundary_ordered_poly)
        if SameValue(boundary_ordered_poly[en1][0],p1[0]) and  SameValue(boundary_ordered_poly[en1][1],p1[1]):
            # loop found
            break;

        for pc2 in range(0, len(index_list), 2):
            pn2 = (pc2 + 1)
            pc2Index = index_list[pc2]
            pn2Index = index_list[pn2]
            #print("Test 1",boundary_ordered_poly[en1],boundary_points[pn2Index])
            #print("Test 2",boundary_ordered_poly[en1],boundary_points[pc2Index])
            if SameValue(boundary_ordered_poly[en1][0],boundary_points[pn2Index][0]) and  SameValue(boundary_ordered_poly[en1][1],boundary_points[pn2Index][1]):
                index_list.pop(pn2)
                index_list.pop(pc2)
                boundary_ordered_poly.append(boundary_points[pc2Index])
                break
            elif SameValue(boundary_ordered_poly[en1][0],boundary_points[pc2Index][0]) and  SameValue(boundary_ordered_poly[en1][1],boundary_points[pc2Index][1]):
                index_list.pop(pn2)
                index_list.pop(pc2)
                boundary_ordered_poly.append(boundary_points[pn2Index])
                break
        if len(index_list)==0:
            # loop found
            break;
    return boundary_ordered_poly
            

def are_adjacent(boundary_points1, boundary_points2):
    """
        returns true if to boundary lists are adjacent(edge groups)
        :param boundary_points1 : list of points on boundary 1
        :param boundary_points2 : list of points on boundary 2
    """
    if boundary_points1 is None or boundary_points2 is None:
        return False
    count1 = len(boundary_points1)
    count2 = len(boundary_points2)
    # note the steps size of two (each edge if in groups of two
    for point1 in range(0, count1, 2):
        # pn1, next point in edge (point1,pn1)
        pn1 = (point1 + 1) % count1
        # note the steps size of two
        for point2 in range(0, count2, 2):
            # pn2, next point in edge (point2,pn2)
            pn2 = (point2 + 1) % count2
            if boundary_points1[point1][2] != boundary_points2[point2][2]:
                return False
            if (SameValue(boundary_points1[point1][0],boundary_points2[point2][0]) and SameValue(boundary_points1[pn1][0],boundary_points2[pn2][0]) and
                SameValue(boundary_points1[point1][1],boundary_points2[point2][1]) and SameValue(boundary_points1[pn1][1],boundary_points2[pn2][1])):
                return True
            if (SameValue(boundary_points1[pn1][0],boundary_points2[point2][0]) and SameValue(boundary_points1[point1][0],boundary_points2[pn2][0]) and
                SameValue(boundary_points1[pn1][1],boundary_points2[point2][1]) and SameValue(boundary_points1[point1][1],boundary_points2[pn2][1])):
                return True


def LandingDefined(GlobalID, landings_list):
    iPos = 0
    for landing in landings_list:
        if landing['GlobalId'] == GlobalID:
            return iPos
        iPos+=1
    return (-1)


def ElevatorDefined(GlobalID, elevator_list):
    iPos = 0
    for elevator in elevator_list:
        if elevator['GlobalId'] == GlobalID:
            return iPos
        iPos+=1
    return (-1)


def StairflightDefined(GlobalID, stairflight_list):
    iPos = 0
    for stair_flight in stairflight_list:
        if stair_flight['GlobalId'] == GlobalID:
            return iPos
        iPos+=1
    return (-1)


def StairsDefined(GlobalID, stair_list):
    iPos = 0
    for stair in stair_list:
        if stair['GlobalId'] == GlobalID:
            return iPos
        iPos+=1
    return -1


def SpaceDefined(GlobalID, space_list):
    iPos = 0
    for space in space_list:
        if space['GlobalId'] == GlobalID:
            return iPos
        iPos+=1
    return -1


def DoorDefined(GlobalID, door_list):
    iPos = 0
    for door in door_list:
        if door['GlobalId'] == GlobalID:
            return iPos
        iPos += 1
    return -1


def find_floor_below_elevation(Elevation, floor_list):
    last_floor = None
    for floor in floor_list:
        if last_floor:
            if floor['Elevation'] == Elevation:
                return last_floor['Elevation']
        last_floor = floor           
    return None


def find_floor_above_elevation(Elevation, floor_list):
    last_floor = None
    for floor in floor_list:
        if last_floor:
            if last_floor['Elevation'] == Elevation:
                return floor['Elevation']
        last_floor = floor           
    return None


def find_floor(Elevation, floor_list):
    floor_index = 0
    for floor in floor_list:
        if floor['Elevation'] == Elevation:
            return floor_index
        floor_index+=1
                  
    return -1


def get_subfloor_list(min_elev, max_elev, floor_list):
    subfloorlist: list[Any] = []
    for floor in floor_list:
        anElevation = floor['Elevation']
        if anElevation >= min_elev and anElevation<max_elev:
            subfloorlist.append(floor)

    if len(subfloorlist)>0:
        return subfloorlist
    return None


def get_subfloor_indices(min_elev, max_elev, floor_list):
    subfloorlist = []
    iPos = 0
    for floor in floor_list:
        anElevation = floor['Elevation']
        if anElevation >= min_elev and anElevation<max_elev:
            subfloorlist.append(iPos)
        iPos+=1

    if len(subfloorlist)>0:
        return subfloorlist
    return None


def get_subfloor_heights(min_elev, max_elev, floor_list):
    subfloorlistheights = []
    iPos = 0
    for floor in floor_list:
        anElevation = floor['Elevation']
        if  anElevation >= min_elev and anElevation<max_elev:
            subfloorlistheights.append(anElevation)
        iPos+=1

    if len(subfloorlistheights)>0:
        return subfloorlistheights
    return None


def inbox(door_box,x,y):
    if x<door_box[2][0] or x>door_box[0][0]:
        return False
    if y<door_box[3][1] or y>door_box[0][1]:
        return False
    return True


def midbox(door_box):
    x = (door_box[2][0] + door_box[0][0])/2.0
    y = (door_box[3][1] + door_box[0][1])/2.0
    return [x,y]


def remove_virtual_lines(lines, spaces, subspaces):
    """
    remove virtual lines between spaces
    :param lines - list of lines on a given level
    :param spaces - spaces on floor
    :param subspaces - subspace on a given floot
    
    """
    for space in spaces:
        subspace_count = len(space['subspaces'])
        for iSpace in range(subspace_count-1):
            iSubSpace= space['subspaces'][iSpace]
            iIndex = subspaces[iSubSpace]['LineIndex']
            iLines = subspaces[iSubSpace]['LineCount']
            for jSpace in range(iSpace+1,subspace_count):
                jSubSpace= space['subspaces'][jSpace]
                jIndex = subspaces[jSubSpace]['LineIndex']
                jLines = subspaces[jSubSpace]['LineCount']
                if subspaces[iSubSpace]['ParentIndex']!=subspaces[jSubSpace]['ParentIndex']:
                    print("Space mismatch")
                for iLine in range(iIndex,iIndex+iLines):
                    for jLine in range(jIndex,jIndex+jLines):
                        if same_point(lines[iLine][0],lines[jLine][0],lines[iLine][1],lines[jLine][1]) and same_point(lines[iLine][2],lines[jLine][2],lines[iLine][3],lines[jLine][3]):
                            print("Remove Clockwise",iLine,jLine)
                            lines[iLine][4]=False
                            lines[jLine][4]=False
                        if same_point(lines[iLine][0],lines[jLine][2],lines[iLine][1],lines[jLine][3]) and same_point(lines[iLine][2],lines[jLine][0],lines[iLine][3],lines[jLine][1]):
                            lines[iLine][4]=False
                            lines[jLine][4]=False
                            print("Remove Anticlockwise",iLine,jLine)
