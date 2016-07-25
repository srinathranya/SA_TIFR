import local_rmsd_analyzer
from copy import deepcopy

def clusterer(moved_alpha_carbons1, moved_alpha_carbons2):

    residue_lrmsd_all_loop = {}

    counter = 0

    residue_list = range(0, len(moved_alpha_carbons1))
    tracker_residue_list = deepcopy(residue_list)
    ref_rmsd_list = []
    clusters = []
    counter = 1
    
    while (tracker_residue_list != []):
        input_alpha_carbons1 = []
        input_alpha_carbons2 = []
        
        #for i in tracker_residue_list:
            #input_alpha_carbons1.append(moved_alpha_carbons1[i])
            #input_alpha_carbons2.append(moved_alpha_carbons2[i])
    
        print ("This is loop %s" %(counter))
        #print ("tracker_residue_list initial - %s " %(tracker_residue_list))
        residue_lrmsd_all_loop[counter], reference_RMSD = local_rmsd_analyzer.local_rmsd_analyzer(moved_alpha_carbons1, moved_alpha_carbons2, tracker_residue_list)
        
        ref_rmsd_list.append(reference_RMSD)
        flexible_cluster = []
        rigid_cluster = []

        for residue in residue_lrmsd_all_loop[counter]:
            temp = []
            for radius in residue_lrmsd_all_loop[counter][residue]:
                temp.append(residue_lrmsd_all_loop[counter][residue][radius])
            
            if (max(temp) > reference_RMSD):
                flexible_cluster.append(residue)
            else:
                rigid_cluster.append(residue)
        
        residue_list = []
        
        for j in flexible_cluster:
            tracker_residue_list.remove(j)
            
        #print ("rigid_cluster - %s" %(rigid_cluster))
        #print ("tracker_residue_list new - %s " %(tracker_residue_list))

        if (flexible_cluster == []):
            clusters.append(rigid_cluster)
            break
        if (rigid_cluster == []):
            clusters.append(flexible_cluster)
            break

        clusters.append(flexible_cluster)
            
        counter += 1
        
    return residue_lrmsd_all_loop, ref_rmsd_list, clusters
