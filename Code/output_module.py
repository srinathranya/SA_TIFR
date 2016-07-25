from Bio.PDB import *

def output_writer(CA_1, CA_2, clusters, residue_lrmsd_ALL, ref_rmsd_list, input_file1, input_file2):
     
    ### WRITING DOWN THE CLUSTER INFORMATION IN A FILE    
    f1 = open('clusters_%s_%s.txt'%(input_file1.split('.')[0], input_file2.split('.')[0]),'w')
    for i in range(0, len(clusters)):
        f1.write('ref-rmsd for loop %s is %s \n' %(i, ref_rmsd_list[i]))
        f1.write('cluster: %s %s \n\n'%(str(i), clusters[i]))
    f1.close()
    
    ### WRITING DOWN THE LOCAL RMSD INFORMATION IN A FILE - CAN BE PLOTTED IN MATHEMATICA AS A 2-D PlOT 
    for k in residue_lrmsd_ALL:
        f1 = open('all_residues_lrmsd_%s_%s_loop_%s.txt' %(input_file1.split('.')[0], input_file2.split('.')[0], k),'w')
        for i in residue_lrmsd_ALL[k]:
            new_str = ''
            for j in residue_lrmsd_ALL[k][i]:
                new_str = new_str + ' ' + str(residue_lrmsd_ALL[k][i][j])
            f1.write('%s \n'%(new_str))
        f1.close() 
        
   
    ### WRITING DOWN THE TWO ALIGNED STRUCTURES - USING THE MOST RIGID CLUSTER FOR ALIGNMENT (THAT IS OUR GLOBAL ALIGNMENT AXIS)
    
    parse1 = PDBParser()
    reference_structure = parse1.get_structure('str1','%s'%(input_file1))
    target_structure = parse1.get_structure('str2','%s'%(input_file2))    

    print CA_1
    print CA_2
    
    reference_atoms = []
    target_atoms = []

    residue_list = []

    for residues in reference_structure.get_residues():
        residue_list.append(residues.get_id()[1])

    for model in reference_structure:
        for chain in model:
            for residues in chain:
                if ((residues.get_id()[1] - residue_list[0]) in clusters[-1]):                
                    reference_atoms.append(residues['CA'])

    reference_for_alignment = reference_structure[0]

    for model in target_structure:
        for chain in model:
            for residues in chain:
                if ((residues.get_id()[1] - residue_list[0]) in clusters[-1]):                
                    target_atoms.append(residues['CA'])

    target_for_alignment = target_structure[0]

    align = Superimposer()
    align.set_atoms(reference_atoms, target_atoms)
    align.apply(target_for_alignment.get_atoms()) 
    
    output = PDBIO()
    output.set_structure(target_structure)
    output.save('Aligned_%s'%(input_file2))        

    return 0
