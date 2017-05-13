#!/usr/bin/python
#
# Written by Sergio Hernandez Charpak
# May 2017
#---------------------------------------------------------
# Imports
#---------------------------------------------------------
#General Imports
#Reading the CIC
from ast import literal_eval
from struct import *
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
#Arguments
import sys
#Running bash
import subprocess
#Getting the filename
import ntpath
#---------------------------------------------------------
# Constants
#---------------------------------------------------------
# Usage
USAGE = "python pipeline_script_2.py file_search file_seeds file_FA_Trace_search fof_src_folder output_folder"
#---------------------------------------------------------
# Functions
#---------------------------------------------------------
#Function to extract the thresholds from a filename
def get_thresholds(filename):
    file_Dat = filename.strip('.dat')
    file_array = file_Dat.split('_')
    thresh_FA = float(file_array[2])
    thresh_Trace = float(file_array[4])
    return thresh_FA, thresh_Trace
#Function to determine is valid or not. If it contains a seed or not
#Each time it founds a seed in a group it deletes it from the seed array to make shorter the search.
# @param int id_group the id of the group to valid
def group_valid(id_group,seeds):
    index = where(fof_groups==id_group)
    index = index[0]
    x_group = particles[index,0]
    y_group = particles[index,1]
    z_group = particles[index,2]
    index_seeds_in_group = []
    for i in range(shape(seeds)[0]):
        seed_i_x = seeds[i,0]
        seed_i_y = seeds[i,1]
        seed_i_z = seeds[i,2]
        index_seed = np.where((x_group == seed_i_x) & (y_group == seed_i_y) & (z_group == seed_i_z) )
        index_seed = index_seed[0]
        if(size(index_seed) !=0):
            index_seeds_in_group.append(i)
    #Now it deletes the seeds from the seed array
    if((size(index_seeds_in_group) > 0)):
        return True, index_seeds_in_group
    else:
        return False, index_seeds_in_group
#Calculates the inertia tensor and eigenvalues given a group
def calcultate_inertia_tensor_overdensity_volume(group_number):
    index = where(fof_groups==group_number)
    index = index[0]
    x_group = particles[index,0]
    y_group = particles[index,1]
    z_group = particles[index,2]
    FA_group = FA_Trace[index,0]
    Trace_group = FA_Trace[index,1]
    #First we need to find the center of mass (CM)
    sum_mx = 0
    sum_my = 0
    sum_mz = 0
    sum_Trace = 0.0
    for i in range(len(x_group)):
        sum_mx +=  x_group[i]
        sum_my +=  y_group[i]
        sum_mz +=  z_group[i]
        sum_Trace += 1.0
    #We get now the x,y,z of the CM
    x_CM = sum_mx / sum_Trace
    y_CM = sum_my / sum_Trace
    z_CM = sum_mz / sum_Trace
    #We now form the Inertia Tensor for the group
    #http://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
    #The positions must be relative to their center of mass
    sum_xy = 0
    sum_yz = 0
    sum_zx = 0
    sum_x_square = 0
    sum_y_square = 0
    sum_z_square = 0
    for i in range(len(x_group)):
        sum_xy += ((x_group[i] - x_CM) * (y_group[i] - y_CM))
        sum_x_square+=((x_group[i] - x_CM)**2)
        sum_yz += ((y_group[i] - y_CM) * (z_group[i] - z_CM))
        sum_y_square+=((y_group[i] - y_CM)**2)
        sum_zx += ((z_group[i] - z_CM) * (x_group[i] - x_CM))
        sum_z_square+=((z_group[i] - z_CM)**2)
    I_xx = sum_y_square + sum_z_square
    I_yy = sum_x_square + sum_z_square
    I_zz = sum_x_square + sum_y_square
    I_xy = - sum_xy
    I_yz = - sum_yz
    I_xz = - sum_zx
    I_matrix = np.matrix([[I_xx, I_xy, I_xz], [I_xy, I_yy, I_yz], [I_xz, I_yz, I_zz]])
    I_a, I_b, I_c = np.linalg.eigvalsh(I_matrix)
    return I_a, I_b, I_c, x_CM, y_CM, z_CM, sum_Trace, len(x_group)
#Returns the name of the file given a path
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
#---------------------------------------------------------
# script
#---------------------------------------------------------
if(len(sys.argv)!=6):
    print ("Please use correctly")
    print USAGE
    sys.exit()
#Typical file_search/seeds: grid_FA_0.6_Trace_1.0.dat
#Typical file_values_FA_Trace: values_FA_0.6_Trace_1.0.dat
#Getting the arguments
file_search = sys.argv[1]
file_seeds = sys.argv[2]
file_values_FA_Trace = sys.argv[3]
fof_src_folder = sys.argv[4]
output_folder = sys.argv[5]
#Extracts the thresholds from the filename
filename_search = path_leaf(file_search)
filename_seeds  = path_leaf(file_seeds)
print ("Files (Search, Seeds):", filename_search, filename_seeds)
search_thresh_FA, search_thresh_Trace = get_thresholds(filename_search)
seeds_thresh_FA, seeds_thresh_Trace = get_thresholds(filename_seeds)
print ("Search file thresholds (FA, Trace): ", search_thresh_FA, search_thresh_Trace)
print ("Seeds file thresholds (FA, Trace): ", seeds_thresh_FA, seeds_thresh_Trace)
#Excecute the FoF script
# -e 0.4 means that the linking length is 0.5
# -m 100 means that it keeps groups of at least 100 particles
comando = fof_src_folder+'/fof -e 1.1 -m 20 < '+'\''+file_search+'\''
print (comando)
process =subprocess.Popen(comando,stdout=subprocess.PIPE, stderr=None, shell=True)
resultsString=process.communicate()
print (resultsString)
#loads the data from the results in the halo finder
fof_groups = loadtxt('fof.grp', skiprows=1)
#The search universe
particles = loadtxt(file_search, skiprows=6)
#The seeds to discriminate the groups
seeds = loadtxt(file_seeds, skiprows=6)
print (fof_groups.size)
print (fof_groups) #each entry correesponds to the ID of the fof group
id_groups = list(set(fof_groups))
print (size(id_groups), 'groups')#this is the number of groups found in the FOF
# Group validation
valid_id_groups = []
print (shape(seeds), 'shape seeds')
for id_group in id_groups:
    valid_condition, index_to_delete = group_valid(id_group,seeds)
    if(valid_condition):
        #It is a valid group
        valid_id_groups.append(id_group)
        seeds = np.delete(seeds, index_to_delete, axis=0)
    if(shape(seeds)[0]==0):
        #No more seeds to evaluate. Can break
        break
print ("There are ", size(valid_id_groups), "valid groups")
print ("Oposite to the ", size(id_groups), "groups at first")
#Now it proceeds to the inertia and group properties analysis
#It needs the values of the Trace
FA_Trace = loadtxt(file_values_FA_Trace, skiprows=6)
#Prepares the output file
output_name = output_folder+'/group_results_seeds_FA_'+str(seeds_thresh_FA)+'_Trace_'+str(seeds_thresh_Trace)+'_search_FA_'+str(search_thresh_FA)+'_Trace_'+str(search_thresh_Trace)+'.dat'
print ("Writing the group results values in file " + output_name)
fileout = open(output_name, 'w')
for id_valid_group in valid_id_groups:
    I_a, I_b, I_c, x_CM, y_CM, z_CM, sum_Trace, volume = calcultate_inertia_tensor_overdensity_volume(id_valid_group)
    fileout.write("%f %f %f %f %f %f %f %f %f \n"%(I_a, I_b, I_c, x_CM, y_CM, z_CM, sum_Trace, volume, id_valid_group ))
fileout.close()
print ("Finished writting the file")
#Now it prepares to write in a separate file the coordinates and the FA and the Trace
#Prepares the output file
output_name = output_folder+'/group_valid_positions_FA_Trace_seeds_FA_'+str(seeds_thresh_FA)+'_Trace_'+str(seeds_thresh_Trace)+'_search_FA_'+str(search_thresh_FA)+'_Trace_'+str(search_thresh_Trace)+'.dat'
print ("Writing the group positions and FA_Trace values in file " + output_name)
fileout = open(output_name, 'w')
for id_valid_group in valid_id_groups:
    index = where(fof_groups==id_valid_group)
    index = index[0]
    x_group = particles[index,0]
    y_group = particles[index,1]
    z_group = particles[index,2]
    FA_group = FA_Trace[index,0]
    Trace_group = FA_Trace[index,1]
    length_group = len(x_group)
    for j in range(length_group):
        fileout.write("%f %f %f %f %f %f \n"%(x_group[j], y_group[j], z_group[j], FA_group[j], Trace_group[j], id_valid_group ))
fileout.close()
print ("Finished pipeline 2")
