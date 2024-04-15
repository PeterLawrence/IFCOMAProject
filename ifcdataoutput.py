import json


def output_names(index_list,ifc_list):
    for index in index_list:
        print(ifc_list[index]['Name'])


def level_connection_data(space_list, stair_list, stair_flights_list, escalator_list):
    print("===================== Level Connections Data ====================")
    for space in space_list:
        GlobalID = space['GlobalId']
        if 'stairsUp' in space or 'stairsDown' in space or 'stairflightsUp' in space or 'stairflightsDown' in space or 'escalatorsUp' in space or 'escalatorsDown' in space: 
            print(f"Stairs in Space {space['Name']} Reference Stairs {space['Reference']}, {space['GlobalId']}")
            
        if 'stairsUp' in space or 'stairsDown' in space:
            
            if 'stairsUp' in space:
                print(" stairsUp:", space['stairsUp'])
                output_names(space['stairsUp'],stair_list)
                
            if 'stairsDown' in space:
                print(" stairsDown:",space['stairsDown'])
                output_names(space['stairsDown'],stair_list)
                
        if 'stairflightsUp' in space or 'stairflightsDown' in space:
                   
            if 'stairflightsUp' in space:
                print(" stairflightsUp:", space['stairflightsUp'])
                output_names(space['stairflightsUp'],stair_flights_list)
                
            if 'stairflightsDown' in space:
                print(" stairflightsDown:",space['stairflightsDown'])
                output_names(space['stairflightsDown'],stair_flights_list)
                      
        if 'escalatorsUp' in space or 'escalatorsDown' in space:
            
            if 'escalatorsUp' in space:
                print(" escalatorsUp:",space['escalatorsUp'])
                output_names(space['escalatorsUp'],escalator_list)
                    
            if 'escalatorsDown' in space:
                print(" escalatorsDown:",space['escalatorsDown'])
                output_names(space['escalatorsDown'],escalator_list)


def space_boundarys(space_list):
    for space in space_list:
        if 'boundarylist' in space:
            GlobalID = space['GlobalId']
            print(f"Space {space['Name']} {space['GlobalId']}")
            print("boundary points",len(space['boundarylist']))


def stair_data(stair_list):
    print("===================== Stair Data ====================")
    for stair in stair_list:
        print(json.dumps(stair,sort_keys=True, indent=4))


def stair_flight_data(stair_flights_list):
    print("===================== Stair Flight Data ====================")
    for stair_flight in stair_flights_list:
        print(json.dumps(stair_flight,sort_keys=True, indent=4))


def landing_data(landings_list):
    print("===================== Landing Data ====================")
    for landing in landings_list:
        print(json.dumps(landing,sort_keys=True, indent=4))


def output_attributes(an_object,attribute1, attribute2):
    if attribute1 in an_object and attribute2 in an_object:
            print(attribute1,":",an_object[attribute1],attribute2,":",an_object[attribute2])


def output_attribute(an_object,attribute):
    if attribute in an_object:
        print(attribute,":",an_object[attribute])


def output_geom_data(object_list):
    print("===================== Geometry Data ====================")
    for an_object in object_list:
        print ("Geometry Data for ",an_object['Name'])
        output_attributes(an_object,'risers','riserheight')
        output_attributes(an_object,'treads','treadlength')
        output_attribute(an_object,'Length')
        output_attribute(an_object,'Width')
        output_attribute(an_object,'HorizontalLength')
        output_attribute(an_object,'Direction')
        output_attribute(an_object,'Elevation')
        output_attribute(an_object,'AStepCount')
        #output_attribute(an_object,'AStepData')
        output_attribute(an_object,'PathData')
        if 'GeomData' in an_object:
            print(json.dumps(an_object['GeomData'],sort_keys=True, indent=4))
        else:
            print("No Geomerty Data")

