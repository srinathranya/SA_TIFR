from Bio.Seq import Seq
from Bio.PDB import * 

####################################################################################
### READING THE TWO STRUCTURES THAT HAVE BEEN PROVIDED FOR ANALYSIS BY THE USER ####
####################################################################################
def obtain_coordinates(file_input1, file_input2):
    parse1 = PDBParser()
    str1 = parse1.get_structure('str1','%s' %(file_input1))
    str2 = parse1.get_structure('str2','%s' %(file_input2))

    residue_list1 = []
    residue_list2 = []

    for residue in str1.get_residues():
        temp = str(residue).split()
        residue_list1.append(temp[1])
    for residue in str2.get_residues():
        temp = str(residue).split() 
        residue_list2.append(temp[1])
    
    alpha_carbon_coordinates1 = {}
    alpha_carbon_coordinates2 = {}
    
    all_coordinates1 = {}
    all_coordinates2 = {}
        
    counter1 = 1
    counter2 = 1
    
    CA_positions1 = []
    CA_positions2 = []
    
    if (len(residue_list1) != len(residue_list2)):
        print "The two proteins do not have the same number of amino acid residues"
    else: 
        print "The proteins have the same number of residues. We are continuing with the calculation."
        for atom in str1.get_atoms():
            s = atom.get_coord()
            all_coordinates1[counter1] = s
            if (atom.get_id() == 'CA'):
                CA_positions1.append(counter1)
            counter1 += 1                        
        for atom in str2.get_atoms():
            s = atom.get_coord()
            all_coordinates2[counter2] = s       
            if (atom.get_id() == 'CA'):
                CA_positions2.append(counter2)
            counter2 += 1         
        for k in range(0, len(CA_positions1)):
            alpha_carbon_coordinates1[k] = all_coordinates1[CA_positions1[k]] 
        for k in range(0, len(CA_positions2)):
            alpha_carbon_coordinates2[k] = all_coordinates2[CA_positions2[k]] 
    
    return all_coordinates1, all_coordinates2, alpha_carbon_coordinates1, alpha_carbon_coordinates2, CA_positions1, CA_positions2    
