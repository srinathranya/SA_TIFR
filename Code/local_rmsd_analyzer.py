import kabsch
from math import sqrt

def distance_computer(max_value, min_value):
    distance = sqrt((max_value[0] - min_value[0])*(max_value[0] - min_value[0]) + (max_value[1] - min_value[1])*(max_value[1] - min_value[1]) + (max_value[2] - min_value[2])*(max_value[2] - min_value[2]))    
    return distance
    
def maximum_radius_computer(moved_alpha_carbons1):
    max_radius = 0
    a_list = []
    b_list = []
    c_list = []
    for i in range(0, len(moved_alpha_carbons1)):
        a_list.append(moved_alpha_carbons1[i][0])
        b_list.append(moved_alpha_carbons1[i][1])
        c_list.append(moved_alpha_carbons1[i][2])
    min_value = []
    max_value = []
    min_value.append(min(a_list))
    min_value.append(min(b_list))
    min_value.append(min(c_list))
    max_value.append(max(a_list))
    max_value.append(max(b_list))
    max_value.append(max(c_list))

    max_radius = distance_computer(max_value, min_value) 
    return int(max_radius)

def local_residue_list_computer(alpha_carbon_number, moved_alpha_carbon_list1, moved_alpha_carbon_list2, radius, res_list):
    local_residue_list1 = []
    local_residue_list2 = []

    for k in res_list:
        
        distance = 0
        distance = distance_computer(moved_alpha_carbon_list1[alpha_carbon_number], moved_alpha_carbon_list1[k])
        if (distance <= radius):
            local_residue_list1.append(k)
    for k in res_list:
        distance = 0
        distance = distance_computer(moved_alpha_carbon_list2[alpha_carbon_number], moved_alpha_carbon_list2[k])
        if (distance <= radius):
            local_residue_list2.append(k)    
    if (len(local_residue_list1) <= len(local_residue_list2)):
        local_residue_list = local_residue_list2
    else:
        local_residue_list = local_residue_list1
        
    return local_residue_list

    
def local_rmsd_analyzer(moved_alpha_carbons1, moved_alpha_carbons2, res_list):
    input_CA_1 = []
    for j in res_list:
        input_CA_1.append(moved_alpha_carbons1[j])

    radii_range = range(8, maximum_radius_computer(input_CA_1) + 5)
    radii_lrmsd_for_all_residues = {}

    if (res_list == []):
        print ("We have a problem")
        raw_input()

    for i in res_list:
        radii_lrmsd_for_all_residues[i] = {}
        for radius in radii_range:
            res_list1 = local_residue_list_computer(i, moved_alpha_carbons1, moved_alpha_carbons2, radius, res_list)
            CA_list_1 = []
            CA_list_2 = []
            for p in res_list1:
                CA_list_1.append(moved_alpha_carbons1[p])
                CA_list_2.append(moved_alpha_carbons2[p])
            t_RMSD, k_RMSD, moved_coords1, moved_coords2 = kabsch.kabsch_RMSD_calculator(CA_list_1, CA_list_2)
            radii_lrmsd_for_all_residues[i][radius] = k_RMSD
    return radii_lrmsd_for_all_residues, k_RMSD
