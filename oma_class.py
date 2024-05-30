# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 09:59:58 2023

@author: P.J.Lawrence

@brief OMA data storage class
"""

import ifcdataoutput
import exoifcutils

class OMAClass:
    m_floor_list = []  # floor data
    m_space_list = []  # list of spaces
    m_door_list = []   # list of doors

    m_stair_list = []  # stair lists
    m_stair_flights_list = []  # stair flights
    m_landings_list = []  # list of landings between stairs

    m_escalator_list = []  # escalator data
    m_elevator_list = []  # elevator data
    m_movingwalkway_list = []  # moving walkway data
    m_ramp_list = []  # ramp data
    m_sign_list = [] # signage data

    def dump_data(self):
        #ifcdataoutput.space_boundarys(self.m_space_list)
        ifcdataoutput.level_connection_data(self.m_space_list, self.m_stair_list, self.m_stair_flights_list, self.m_escalator_list)
        ifcdataoutput.stair_data(self.m_stair_list)
        ifcdataoutput.stair_flight_data(self.m_stair_flights_list)
        ifcdataoutput.landing_data(self.m_landings_list)
        ifcdataoutput.output_geom_data(self.m_stair_flights_list)
        ifcdataoutput.output_geom_data(self.m_escalator_list)
        ifcdataoutput.output_geom_data(self.m_ramp_list)
        ifcdataoutput.output_geom_data(self.m_sign_list)

    def SpaceDefined(self, GlobalID):
        iPos = 0
        for space in self.m_space_list:
            if space['GlobalId'] == GlobalID:
                return iPos
            iPos+=1
        return -1

    def dump_space_boundaries(self):
        '''
        Output boundary data
        Note this boundary data is form of edge groups
        i.e. points 0 and 1 are one edge, points 2 and 3 is another edge and so on
        '''
        print('Space Boundary Data')
        for space in self.m_space_list:
            print(space['Name'],space['Elevation'], space['boundarylist'])
        print('End of Space Boundary Data')

    def dump_ordered_space_boundaries(self):
        '''
        output boundary data as an ordered polygon boundary list
        '''
        print('Space Ordered Boundary Data')
        for space in self.m_space_list:
            boundary_list = exoifcutils.order_boundary_list(space['boundarylist'])
            print(space['Name'],space['Elevation'], boundary_list)
        print('End of Ordered Space Boundary Data')
