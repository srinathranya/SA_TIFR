import scipy.linalg
from math import sqrt
import numpy

# THIS FUNCTION COMPUTES THE CENTROID OF THE GIVEN SET OF COORDINATES

def centroid_computer(coordinates):
    centroid = []
    sum_1 = 0 
    sum_2 = 0 
    sum_3 = 0
    for i in range(0, len(coordinates)):
        sum_1 += coordinates[i][0]
        sum_2 += coordinates[i][1]
        sum_3 += coordinates[i][2]
    centroid.append(sum_1/len(coordinates))
    centroid.append(sum_2/len(coordinates))
    centroid.append(sum_3/len(coordinates))
    return centroid

# THIS FUNCTION COMPUTES THE TRANSLATED COORDINATES OF BOTH THE STRUCTURES
# WHEN THE CENTROIDS OF BOTH THE STRUCTURES HAVE BEEN MOVED TO THE ORIGIN

def structure_translator(coordinate, centroid):
    translated_coordinates = []
    for i in range(0, len(coordinate)):
        translated_coordinates.append([])
        translated_coordinates[i].append(coordinate[i][0] - centroid[0])
        translated_coordinates[i].append(coordinate[i][1] - centroid[1])
        translated_coordinates[i].append(coordinate[i][2] - centroid[2])
    return translated_coordinates

# THIS FUNCTION COMPUTES THE RMSD FOR A TWO SETS OF COORDINATES

def rmsd_calculator(coordinate_1, coordinate_2):
    rmsd = 0 
    rmsd_sum = 0 
    for i in range(0, len(coordinate_1)):
        for j in range(0, len(coordinate_1[i])):
            rmsd_sum = rmsd_sum + (float(coordinate_1[i][j]) - float(coordinate_2[i][j]))*(float(coordinate_1[i][j]) - float(coordinate_2[i][j]))
    rmsd = sqrt(rmsd_sum/(len(coordinate_1)))
    return rmsd

# THIS FUNCTION UTILIZES THE KABSCH ALGORITHM TO COMPUTE THE TRANSLATION
# AND ROTATIONS THAT ARE TO BE PERFORMED TO MINIMIZE THE RMSD BETWEEN THE 
# TWO STRUCTURES THAT HAVE BEEN PROVIDED AS INPUT

def kabsch_RMSD_calculator(coordinates_1, coordinates_2):
    kabsch_RMSD = 0
    centroid_1 = centroid_computer(coordinates_1)
    centroid_2 = centroid_computer(coordinates_2)
    translated_structure_1 = structure_translator(coordinates_1, centroid_1)
    translated_structure_2 = structure_translator(coordinates_2, centroid_2)
    translated_rmsd = rmsd_calculator(translated_structure_1, translated_structure_2)
    covariance_matrix = numpy.dot(numpy.transpose(translated_structure_1), translated_structure_2) 
    covariance_transpose = numpy.transpose(covariance_matrix)
    U, S, V = numpy.linalg.svd(covariance_matrix)
    d = numpy.linalg.det(numpy.dot(numpy.transpose(V),numpy.transpose(U)))
    optimal_rotation_matrix = numpy.dot(numpy.dot(numpy.transpose(V),numpy.array([[1,0,0],[0,1,0],[0,0,d]])),numpy.transpose(U))
    rotated_coordinates_2 = numpy.dot(translated_structure_2, optimal_rotation_matrix)
    kabsch_rmsd = rmsd_calculator(translated_structure_1, rotated_coordinates_2)    
    
    return translated_rmsd, kabsch_rmsd, translated_structure_1, rotated_coordinates_2
