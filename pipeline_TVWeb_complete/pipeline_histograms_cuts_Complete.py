#Written by Sergio Daniel Hernandez Charpak
#May 2017
#Interact with the Fractional Anisotropy (FA) grid directly from the DensTools files.
#Here we read the results and plot them.
#---------------------------------------------
#Imports
#General Imports
from ast import literal_eval
from struct import *
import numpy as np
import matplotlib.pyplot as plt
import sys
from pylab import *
#---------------------------------------------
USAGE = "python pipeline_histograms_cuts_Complete.py file_eigen_vale_1 file_eigen_vale_2 file_eigen_vale_3 scale_factor_kpc grid cut_i output_folder"
#---------------------------------------------
#Functions
def writeFirstLine(filename,line):
    with open(filename,'r+') as f:
        content = f.read()
        f.seek(0,0)
        f.write(line.rstrip('\r\n') + '\n' + content)
def readFirstLine(filename):
    with open(filename, 'r') as f:
        first_line = f.readline()
        return first_line
#Function to read the EigenValues from the v_web CIC output
def read_CIC_scalar(filename):
    f = open(filename, "rb")
    dumb = f.read(38)

    dumb = f.read(4)
    n_x = f.read(4)
    n_y = f.read(4)
    n_z = f.read(4)
    nodes = f.read(8)
    x0 = f.read(4)
    y0 = f.read(4)
    z0 = f.read(4)
    dx = f.read(4)
    dy = f.read(4)
    dz = f.read(4)
    dumb = f.read(4)

    n_x = (unpack('i', n_x))[0]
    n_y = (unpack('i', n_y))[0]
    n_z = (unpack('i', n_z))[0]
    nodes = (unpack('q', nodes))[0]
    dx = (unpack('f', dx))[0]
    dy = (unpack('f', dy))[0]
    dz = (unpack('f', dz))[0]
    x0 = (unpack('f', x0))[0]
    y0 = (unpack('f', y0))[0]
    z0 = (unpack('f', z0))[0]
    print (n_x, n_y, n_z, nodes, dx, dy, dz)

    total_nodes = n_x * n_y *n_z
    dumb = f.read(4)
    array_data = f.read(total_nodes*4)
    dumb = f.read(4)
    format_s = str(total_nodes)+'f'
    array_data = unpack(format_s, array_data)
    f.close()
    array_data  = np.array(array_data)
    array_data.resize(n_z,n_y,n_x)
    array_data = array_data.transpose()
    return array_data
#---------------------------------------------------------
# script
#---------------------------------------------------------
if(len(sys.argv)!=8):
    print ("Please use correctly")
    print USAGE
    sys.exit()

inputfile_1 = sys.argv[1]
inputfile_2 = sys.argv[2]
inputfile_3 = sys.argv[3]
scale_factor_kpc = float(sys.argv[4])
scale_factor_Mpc = scale_factor_kpc/1000.0
scale_factor_Gpc = scale_factor_Mpc/1000.0
grid = int(sys.argv[5])
cut_i = int(sys.argv[6])
output_folder = sys.argv[7]
print ("Output_folder: ", output_folder)
#Reading Files
print ("Preparing to read the eigenvalues")
print (inputfile_1,inputfile_2,inputfile_3)
grid_1 = read_CIC_scalar(inputfile_1)
grid_2 = read_CIC_scalar(inputfile_2)
grid_3 = read_CIC_scalar(inputfile_3)
print ("Eigenvalues read")
print ("Forming the FA and the Trace grids")
FA = (grid_1-grid_3)**2  + (grid_2-grid_3)**2  + (grid_1-grid_2)**2
FA = FA/(grid_1**2 + grid_2**2 + grid_3**2)
FA = np.sqrt(FA)/np.sqrt(3.0)
Trace = grid_1 + grid_2 + grid_3
print ("Shape FA: ", FA.shape)
print ("Shape VDH: ", Trace.shape)
print ("FA and VDH matrix formed")

#Experiment - Removing invalid values
#Removing invalid values
value_FA_INV=-1.0
value_VDH_INV=-5000.0
where_are_NaNs = np.isnan(FA)
FA[where_are_NaNs] = value_FA_INV
Trace[where_are_NaNs] = value_VDH_INV
n_x,n_y,n_z = FA.shape

#---------------------------------------------
#Cuts
x_thicks = np.linspace(0, n_x, 5)
y_thicks = np.linspace(0, n_y, 5)
x_thicks_mpc = x_thicks*scale_factor_Mpc
xticklabels=["%.2f" % x for x in x_thicks_mpc]
y_thicks_mpc = y_thicks*scale_factor_Mpc
yticklabels=["%.2f" % y for y in y_thicks_mpc]
#FA

output_name = output_folder+'/FA_cut_i_'+str(cut_i)+'.png'
cut_i_mpc = cut_i*scale_factor_Mpc
cut = FA[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut.T,interpolation='none',origin='lower')
colorbar(label='FA value')
plt.xlabel("y [Mpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.ylabel("z [Mpc/h]", fontsize=16)
plt.title('FA cut for x = '+str(cut_i_mpc)+' Mpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," FA min: ", cut.min(),"FA max: ", cut.max())
plt.close(fig)
#DVH
output_name = output_folder+'/VDH_cut_i_'+str(cut_i)+'.png'
cut_trace = Trace[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut_trace.T, cmap = cm.Greys_r,interpolation='none',origin='lower')
plt.colorbar(label='VDH')
plt.xlabel("y [Mpc/h]", fontsize=16)
plt.ylabel("z [Mpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.title('VDH cut for \n x = '+str(cut_i_mpc)+' Mpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," VDH min: ", cut_trace.min()," max: ", cut_trace.max())
plt.close(fig)
#---------------------------------------------
#Histograms
#FA
#Preparation
FA_1D_orig = FA.ravel()
print("There are some invalid values (inf). Where the CIC wasn't succesfull (no particles in a cell.")
print ("Before: ",FA_1D_orig.shape)
FA_1D = FA_1D_orig[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",FA_1D.shape)
#Figure
output_name = output_folder+'/FA_hist_sim.png'
fig = plt.figure(figsize=(8,8))
binwidth=0.01
plt.hist(FA_1D, bins=np.arange(min(FA_1D), max(FA_1D) + binwidth, binwidth))
plt.gca().set_yscale("log")
plt.xlabel('FA',fontsize=18)
plt.xlim(0.0,1.0)
plt.ylabel('Number of cells (log scale)',fontsize=18)
plt.title('FA distribution in the simulation ', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)
#VDH
#Preparation
Trace_1D = Trace.ravel()
print ("VDH")
print ("Max, Min")
print (np.max(Trace_1D), np.min(Trace_1D))
print "There are some outliers as there where in the FA. These are the result of the CIC and v_web in cells with no particles at all."
print ("Before: ",Trace_1D.shape)
Trace_1D = Trace_1D[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Trace_1D.shape)
#Figure
#Big
output_name = output_folder+'/VDH_hist_sim_big.png'
fig = plt.figure(figsize=(8,8))
big_number = 10.0**(4.0)
x_min = -2.2*big_number
x_max = 2.2*big_number
binwidth=0.01*big_number
plt.hist(Trace_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Trace_1D)
plt.gca().set_yscale("log")
#plt.gca().set_xscale("log")
plt.xlabel('VDH',fontsize=18)
plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)')
plt.title('VDH distribution', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)
#Small
#Preparation
print ("There is still a strange behaviour due to missing data.")
print ("Before: ",Trace_1D.shape)
Trace_1D = Trace_1D[(Trace_1D <= 2100.0) & (Trace_1D>=-2100.0)]
print ("After: ",Trace_1D.shape)
#Figure
output_name = output_folder+'/VDH_hist_sim_small.png'
fig = plt.figure(figsize=(8,8))
x_min = -2200.0
x_max = 2200.0
binwidth=5.0
plt.hist(Trace_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Trace_1D)
plt.gca().set_yscale("log")
#plt.gca().set_xscale("log")
plt.xlabel('VDH',fontsize=18)
plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)',fontsize=18)
plt.title('VDH distribution ', fontsize=18)
plt.savefig(output_name,format = 'png')
plt.close(fig)
#Log scale
#Figure
output_name = output_folder+'/VDH_hist_sim_small_log.png'
fig = plt.figure(figsize=(8,8))
x_min = -2200.0
x_max = 2200.0
binwidth=5.0
plt.hist(Trace_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Trace_1D)
plt.gca().set_yscale("log")
plt.gca().set_xscale("log")
plt.xlabel('VDH',fontsize=18)
#plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)',fontsize=18)
plt.title('VDH distribution ', fontsize=18)
plt.savefig(output_name,format = 'png')
plt.close(fig)
#-------------------------------------------------
