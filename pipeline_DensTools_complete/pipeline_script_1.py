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
#---------------------------------------------------------
# Constants
#---------------------------------------------------------
# Usage
USAGE = "python pipeline_script_2.py file_search file_seeds file_FA_Trace_search fof_src_folder output_folder"
#---------------------------------------------------------
# Functions
#---------------------------------------------------------
#Function to read the values from the DensTools pipeline
def read_grid(prefix,nfiles) :
    f=open(prefix+".0000","rb")
    num_grids,ngrid=np.fromfile(f,dtype=np.int32,count=2)
    f.close()
    print ("Will read %d fields"%num_grids+" with %d^3 nodes"%ngrid)
    grid_out=np.zeros([ngrid,ngrid,ngrid,num_grids])
    for ifil in np.arange(nfiles) :
        f=open(prefix+".%04d"%ifil,"rb")
        nug,ng=np.fromfile(f,dtype=np.int32,count=2)
        if (nug!=num_grids) or (ng!=ngrid) :
            print ("shit")
            sys.exit(1)
        nx_here=np.fromfile(f,dtype=np.int32,count=1)[0]
        print ("File #%d"%(ifil+1)+", %d slices found"%nx_here)
        for ix in np.arange(nx_here) :
            ix_this=np.fromfile(f,dtype=np.int32,count=1)[0]
            grid_out[ix_this,:,:,:]=np.fromfile(f,dtype=np.float32,count=ng*ng*nug).reshape([ng,ng,nug])
        f.close()
    if num_grids==1 :
        grid_out=grid_out[:,:,:,0]
    return grid_out
#---------------------------------------------------------
# script
#---------------------------------------------------------
if(len(sys.argv) < 3):
    print ("Please use correctly")
    print (USAGE)
    sys.exit()
#Number of files
n_files = int(sys.argv[1])
if(len(sys.argv) < (n_files+2)):
    print ("Please use correctly")
    print (USAGE)
    sys.exit()
i_files=2
#File list
prefix_list=[]
while(i_files<(n_files+2)):
    prefix_input_i = sys.argv[i_files]
    prefix_list.append(prefix_input_i)
    i_files+=1
print ("Input file list: ",prefix_list)
nfiles_input=1
#Output name formation
thresh_FA = float(sys.argv[n_files+2])
thresh_Trace = float(sys.argv[n_files+3])
print ("The Thresholds are: FA, Trace", str(thresh_FA), str(thresh_Trace))
output_folder = sys.argv[n_files+4]
print ("Output_folder: ", output_folder)
output_name = output_folder+'/grid_FA_'+str(thresh_FA)+'_Trace_'+str(thresh_Trace)+'.dat'
#First File - Template for the rest
prefix_input = prefix_list[0]
print ("Forming the FA and VDH original grid with #files: ", str(n_files))
print ("Reading from file with prefix: ", prefix_input)
vel1=read_grid(prefix_input+"_vel",nfiles_input)
FA = (vel1[:,:,:,0]-vel1[:,:,:,2])**2  + (vel1[:,:,:,1]-vel1[:,:,:,2])**2  + (vel1[:,:,:,0]-vel1[:,:,:,1])**2
FA = FA/(vel1[:,:,:,0]**2 + vel1[:,:,:,1]**2 + vel1[:,:,:,2]**2)
FA = np.sqrt(FA)/np.sqrt(3.0)
Trace = vel1[:,:,:,0] + vel1[:,:,:,1] + vel1[:,:,:,2]
print ("Shape FA: ", FA.shape)
print ("Shape VDH: ", Trace.shape)
print ("FA and VDH matrix formed")
n_x,n_y,n_z = shape(FA)
#Removing invalid values
value_FA_INV=-1.0
value_VDH_INV=-5000.0
where_are_NaNs = np.isnan(FA)
FA[where_are_NaNs] = value_FA_INV
Trace[where_are_NaNs] = value_VDH_INV
#The rest of files
i_files=1
while (i_files<n_files):
    prefix = prefix_list[i_files]
    print ("Reading from file with prefix: ", prefix)
    vel=read_grid(prefix+"_vel",nfiles_input)
    FA_new = (vel[:,:,:,0]-vel[:,:,:,2])**2  + (vel[:,:,:,1]-vel[:,:,:,2])**2  + (vel[:,:,:,0]-vel[:,:,:,1])**2
    FA_new = FA_new/(vel[:,:,:,0]**2 + vel[:,:,:,1]**2 + vel[:,:,:,2]**2)
    FA_new = np.sqrt(FA_new)/np.sqrt(3.0)
    Trace_new = vel[:,:,:,0] + vel[:,:,:,1] + vel[:,:,:,2]
    #Removing NaNs
    where_are_NaNs_new = np.isnan(FA_new)
    FA_new[where_are_NaNs_new] = value_FA_INV
    Trace_new[where_are_NaNs_new] = value_VDH_INV
    #Populating the FA and VDH grids
    i_files += 1
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                if(FA_new[i,j,k] != value_FA_INV):
                    if(FA[i,j,k] == value_FA_INV):
                        FA[i,j,k] = FA_new[i,j,k]
                        Trace[i,j,k] = Trace_new[i,j,k]
                    else:
                        FA[i,j,k] = np.mean([FA[i,j,k],FA_new[i,j,k]])
                        Trace[i,j,k] = np.mean([Trace[i,j,k],Trace_new[i,j,k]])
print ("Finished forming the input FA and VDH grids")
print ("Forming the resulting grid for Thresholds")
cell_list = []
cells_total = 0
for i in range (n_x):
    for j in range (n_y):
        for k in range (n_z):
            if((FA[i,j,k]<thresh_FA ) & (FA[i,j,k]>value_FA_INV ) & (Trace[i,j,k]>thresh_Trace) & (Trace[i,j,k]>value_VDH_INV)):
                cells_total +=1
                cell_list.append([i,j,k, FA[i,j,k], Trace[i,j,k]])
print ("Grid formed with number of cells:", cells_total)
print ("Writing in file " + output_name)
fileout = open(output_name, 'w')
fileout.write("%d\n"%cells_total) #points in total
fileout.write("%d\n"%cells_total) #points in 'DM'
fileout.write("0\n") #gas
fileout.write("0\n") #stars
fileout.write("0.01\n") # time
fileout.write("0\n") # nactive
for l in range(cells_total):
    fileout.write("%f %f %f\n"%(cell_list[l][0],cell_list[l][1],cell_list[l][2]))
fileout.close()
print ("Finished writting the file")
#It prints the values of the FA and Trace for the later analysis
output_name = output_folder+'/values_FA_'+str(thresh_FA)+'_Trace_'+str(thresh_Trace)+'.dat'
print ("Writing in FA and Trace values in file " + output_name)
fileout = open(output_name, 'w')
fileout.write("%d\n"%cells_total) #points in total
fileout.write("%d\n"%cells_total) #points in 'DM'
fileout.write("0\n") #gas
fileout.write("0\n") #stars
fileout.write("0.01\n") # time
fileout.write("0\n") # nactive
l = 0
for l in range(cells_total):
    fileout.write("%f %f\n"%(cell_list[l][3],cell_list[l][4]))
fileout.close()
print ("Finished pipeline 1")
