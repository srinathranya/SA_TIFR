#/*****************************************************************/#
# 									#
#             ALIGNMENT AXIS PROVIDER PROGRAM ( A-A-P )             #
# 									#
#          RAVINDRA VENKATRAMANI’S LAB (TIFR, MUMBAI, INDIA)        #
# 									#
#/*****************************************************************/#



####################################################################
##                                                                ## 
## IMPORTING ALL THE REQUIRED BIOPYTHON MODULES AND THE NECESSARY ##
## INPUT/OUTPUT AND LOCAL-RMSD MODULES FROM OUR SET OF MODULES    ## 
##                                                                ## 
####################################################################

import time
import input_analyzer    
import kabsch
import clustering_module

#/***********************************************************************************/

import output_module     /**** MODULE CAN BE IMPROVED DEPENDING ON FURTHER WORK  ****/

#/***********************************************************************************/

#import output_text_graphs    /** CAN BE CREATED OR CAN BE TREATED AS PART OF output_module **/        
#			       /** THE IDEA IS TO OUTPUT RELEVANT GRAPHS FROM THE DATA       **/
#/***********************************************************************************/




#/***********************************************************************************/
#import spectral_index_calculator      /*** THIS MODULE WILL BE WRITTEN BY SANJOY ***/
#/***********************************************************************************/


######################################################################
####                                                              ####
####    THIS COMMAND PRINTS THE MODULES THAT HAVE BEEN IMPORTED   ####
####                                                              ####
######################################################################
#                                            |  
#import sys                                  |
#print sys.modules.keys()            <=======|                     

#############################################################
##                                                      #####
## ASKING FOR THE FILENAMES FOR USER'S INPUT STRUCTURES #####
##                                                      #####
#############################################################

print ("Please make sure that the PDB file contains only the protein of interest and no water. If it does, please use softwares like VMD or PyMol to generate a new PDB which contains only the protein of interest")

file_input1 = raw_input("Please enter the file name for the first structure (with its extension) : ")
file_input2 = raw_input("Please enter the file name for the second structure (with its extension) : ")

print ("The structures that you have asked to be analyzed are %s and %s" %(file_input1, file_input2))

#######################################################
###                                                 ###
### OBTAINING THE COORDINATES OF THE ALPHA CARBONS  ###
### IN THE TWO INPUT STRUCTURES                     ###
###                                                 ###
#######################################################

start_time = time.time()

all_coordinates1, all_coordinates2, alpha_carbon_positions1, alpha_carbon_positions2, CA_positions1, CA_positions2 = input_analyzer.obtain_coordinates(file_input1, file_input2)

#/***************************************************/# --- TO BE USED TO DEBUG, IF NECESSARY
"""
print all_coordinates1
print all_coordinates2
print alpha_carbon_positions1
print alpha_carbon_positions2

print len(alpha_carbon_positions1)
print len(alpha_carbon_positions2)

print CA_positions1
print CA_positions2

raw_input()
"""
#/**************************************************/#

###########################################################
###                                                     ###
###                                                     ###
### COMPUTING THE OPTIMAL ROTATIONAL, THE TRANSLATIONAL ###
### MATRIX TO MOVE ONE SET OF ALPHA CARBONS ONTO THE    ###
### OTHER AND WRITING DOWN THE NEW COORDINATES.         ###
###                                                     ###
###                                                     ###
###########################################################

#/********************************************************/#
###  !!!!   DO NOT TOUCH THIS PART OF THE CODE  !!!!!  ####

coordinates1 = []
coordinates2 = []

for i in range(0, len(alpha_carbon_positions1)):
    coordinates1.append(alpha_carbon_positions1[i])
    coordinates2.append(alpha_carbon_positions2[i])                                                                                             

translated_rmsd, kabsch_RMSD, moved_alpha_carbons1, moved_alpha_carbons2 = kabsch.kabsch_RMSD_calculator(coordinates1, coordinates2)          

#### !!!!  DO NOT TOUCH THIS PART OF THE CODE  !!!!!! ######
#/********************************************************/#

############################################################
####                                                   #####
####                                                   #####
#### WE NOW INVOKE THE CLUSTERING MODULE AND BIN THE   #####
#### RESIDUES ACCORDING TO THEIR FLEXIBILITY. LOCAL    #####
#### RMSDS FOR EACH OF THE RESIDUES WILL BE STORED.    #####
####                                                   #####
####                                                   #####
############################################################

residue_lrmsd_all_residues_all_clusters, ref_RMSD_list_clusters, clusters_list = clustering_module.clusterer(moved_alpha_carbons1, moved_alpha_carbons2)

############################################################
####                                                   ##### 
####                                                   #####
#### HERE, WE MUST WRITE DOWN AN OUTPUT FUNCTION WHICH #####
#### WILL GIVE US THE ALIGNED PDBs.                    #####
####                                                   #####
####                                                   #####
############################################################

output_module.output_writer(CA_positions1, CA_positions2, clusters_list, residue_lrmsd_all_residues_all_clusters, ref_RMSD_list_clusters, file_input1, file_input2)

print ("This program ran successfully")
print ("Execution time = %s" %(time.time() - start_time))
