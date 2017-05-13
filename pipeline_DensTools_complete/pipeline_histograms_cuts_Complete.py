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
USAGE = "python pipeline_histograms_cuts_Complete.py n_files prefix_files(n_files) scale_factor_kpc cut_i output_folder"
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
#Function to read the reasults from DensTools
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
#---------------------------------------------

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
scale_factor_kpc = float(sys.argv[n_files+2])
scale_factor_Mpc = scale_factor_kpc/1000.0
scale_factor_Gpc = scale_factor_Mpc/1000.0
cut_i = int(sys.argv[n_files+3])
#Output_folder
output_folder = sys.argv[n_files+4]
#Reading Files
print ("Output_folder: ", output_folder)
nfiles_input=1
#First file
prefix_input1 = prefix_list[0]
print ("Reading from file with prefix: ", prefix_input1)
#Reads the files with the velocity
vel1=read_grid(prefix_input1+"_vel",nfiles_input)

FA = (vel1[:,:,:,0]-vel1[:,:,:,2])**2  + (vel1[:,:,:,1]-vel1[:,:,:,2])**2  + (vel1[:,:,:,0]-vel1[:,:,:,1])**2
FA = FA/(vel1[:,:,:,0]**2 + vel1[:,:,:,1]**2 + vel1[:,:,:,2]**2)
FA = np.sqrt(FA)/np.sqrt(3.0)

print ("Shape FA: ", FA.shape)
Trace = vel1[:,:,:,0] + vel1[:,:,:,1] + vel1[:,:,:,2]
print ("Shape VDH: ", Trace.shape)

"""
dens=read_grid(prefix_input1+"_dens",nfiles_input)
dens_sm=read_grid(prefix_input1+"_dens_sm",nfiles_input)
lvel=read_grid(prefix_input1+"_linvel",nfiles_input)
tid=read_grid(prefix_input1+"_tidal",nfiles_input)
tidal = tid[:,:,:,0]+tid[:,:,:,1]+tid[:,:,:,2]
lineal_vel = np.sqrt(lvel[:,:,:,0]**2+lvel[:,:,:,1]**2+lvel[:,:,:,2]**2)/np.sqrt(vel1[:,:,:,0]**2+vel1[:,:,:,1]**2+vel1[:,:,:,2]**2) -1

print ("Raw Mins and Max for Density, DensitySM, Tidal and Linear Vel.")
print ("Density min: ", dens.min()," max: ", dens.max())
print ("Density SM min: ", dens_sm.min()," max: ", dens_sm.max())
print ("Tidal min: ", tidal.min()," max: ", tidal.max())
print ("Linear Vel min: ", lineal_vel.min()," max: ", lineal_vel.max())
"""
#Experiment - Removing invalid values
where_are_NaNs = np.isnan(FA)
FA[where_are_NaNs] = -0.1
Trace[where_are_NaNs] = -2500.0
#Getting only the valid values
"""
dens[where_are_NaNs] = -1.01
dens_sm[where_are_NaNs] = -1.01
tidal[where_are_NaNs] = -0.1
lineal_vel[where_are_NaNs] = -1.01

print ("No Nans Mins and Max for Density, DensitySM, Tidal and Linear Vel.")
print ("Density min: ", dens.min()," max: ", dens.max())
print ("Density SM min: ", dens_sm.min()," max: ", dens_sm.max())
print ("Tidal min: ", tidal.min()," max: ", tidal.max())
print ("Linear Vel min: ", lineal_vel.min()," max: ", lineal_vel.max())
"""
n_x,n_y,n_z = FA.shape

i_files=1
while (i_files<n_files):
    prefix = prefix_list[i_files]
    print ("Reading from file with prefix: ", prefix)
    vel=read_grid(prefix+"_vel",nfiles_input)

    FA_new = (vel[:,:,:,0]-vel[:,:,:,2])**2  + (vel[:,:,:,1]-vel[:,:,:,2])**2  + (vel[:,:,:,0]-vel[:,:,:,1])**2
    FA_new = FA_new/(vel[:,:,:,0]**2 + vel[:,:,:,1]**2 + vel[:,:,:,2]**2)
    FA_new = np.sqrt(FA_new)/np.sqrt(3.0)

    Trace_new = vel[:,:,:,0] + vel[:,:,:,1] + vel[:,:,:,2]
    """
    dens_new=read_grid(prefix+"_dens",nfiles_input)
    dens_sm_new=read_grid(prefix+"_dens_sm",nfiles_input)
    lvel_new=read_grid(prefix+"_linvel",nfiles_input)
    tid_new=read_grid(prefix+"_tidal",nfiles_input)
    tidal_new = tid_new[:,:,:,0]+tid_new[:,:,:,1]+tid_new[:,:,:,2]
    lineal_vel_new = np.sqrt(lvel[:,:,:,0]**2+lvel[:,:,:,1]**2+lvel[:,:,:,2]**2)/np.sqrt(vel1[:,:,:,0]**2+vel1[:,:,:,1]**2+vel1[:,:,:,2]**2) -1

    """
    where_are_NaNs_new = np.isnan(FA_new)
    FA_new[where_are_NaNs_new] = -0.1
    Trace_new[where_are_NaNs_new] = -2500.0
    """
    dens_new [where_are_NaNs_new] = -1.01
    dens_sm_new [where_are_NaNs_new] = -1.01
    tidal_new [where_are_NaNs_new] = -0.1
    lineal_vel_new [where_are_NaNs_new] = -1.01
    """
    i_files += 1
    for i in range(n_x):
        for j in range(n_y):
            for k in range(n_z):
                if(FA_new[i,j,k] != -0.1):
                    if(FA[i,j,k] == -0.1):
                        FA[i,j,k] = FA_new[i,j,k]
                        Trace[i,j,k] = Trace_new[i,j,k]
                        """
                        dens[i,j,k] = dens_new[i,j,k]
                        dens_sm[i,j,k] = dens_sm_new[i,j,k]
                        tidal[i,j,k] = tidal_new[i,j,k]
                        lineal_vel[i,j,k] = lineal_vel_new[i,j,k]
                        """
                    else:
                        FA[i,j,k] = np.mean([FA[i,j,k],FA_new[i,j,k]])
                        Trace[i,j,k] = np.mean([Trace[i,j,k],Trace_new[i,j,k]])
                        """
                        dens[i,j,k] = np.mean([dens[i,j,k],dens_new[i,j,k]])
                        dens_sm[i,j,k] = np.mean([dens_sm[i,j,k],dens_sm_new[i,j,k]])
                        tidal[i,j,k] = np.mean([tidal[i,j,k],tidal_new[i,j,k]])
                        lineal_vel[i,j,k] = np.mean([lineal_vel[i,j,k],lineal_vel_new[i,j,k]])
                        """


#---------------------------------------------
#Cuts
x_thicks = np.linspace(0, n_x, 5)
y_thicks = np.linspace(0, n_y, 5)
x_thicks_gpc = x_thicks*scale_factor_Gpc
xticklabels=["%.2f" % x for x in x_thicks_gpc]
y_thicks_gpc = y_thicks*scale_factor_Gpc
yticklabels=["%.2f" % y for y in y_thicks_gpc]
#FA

output_name = output_folder+'/FA_cut_i_'+str(cut_i)+'.png'
cut_i_gpc = cut_i*scale_factor_Gpc
cut = FA[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut.T,interpolation='none',origin='lower')
colorbar(label='FA value')
plt.xlabel("y [Gpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.ylabel("z [Gpc/h]", fontsize=16)
plt.title('FA cut for x = '+str(cut_i_gpc)+' Gpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," FA min: ", cut.min(),"FA max: ", cut.max())
plt.close(fig)
#DVH
output_name = output_folder+'/VDH_cut_i_'+str(cut_i)+'.png'
cut_trace = Trace[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut_trace.T, cmap = cm.Greys_r,interpolation='none',origin='lower')
plt.colorbar(label='VDH')
plt.xlabel("y [Gpc/h]", fontsize=16)
plt.ylabel("z [Gpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.title('VDH cut for \n x = '+str(cut_i_gpc)+' Gpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," VDH min: ", cut_trace.min()," max: ", cut_trace.max())
plt.close(fig)
#---------------------------------------------
#Histograms
#FA
#Preparation
FA_1D_orig = FA.ravel()
print ("There are some invalid values (inf). Where the CIC wasn't succesfull (no particles in a cell.")
print ("Before: ",FA_1D_orig.shape)
Hist_1D = FA_1D_orig[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Hist_1D.shape)
#Figure
output_name = output_folder+'/FA_hist_sim.png'
fig = plt.figure(figsize=(8,8))
binwidth=0.01
plt.hist(Hist_1D, bins=np.arange(min(Hist_1D), max(Hist_1D) + binwidth, binwidth))
plt.gca().set_yscale("log")
plt.xlabel('FA',fontsize=18)
plt.xlim(0.0,1.0)
plt.ylabel('Number of cells (log scale)',fontsize=18)
plt.title('FA distribution in the simulation ', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)
#VDH
#Preparation
Hist_1D = Trace.ravel()
print ("VDH")
print ("Max, Min")
print (np.max(Hist_1D), np.min(Hist_1D))
print ("There are some outliers as there where in the FA. These are the result of the CIC and v_web in cells with no particles at all.")
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Hist_1D.shape)
#Figure
#Big
output_name = output_folder+'/VDH_hist_sim_big.png'
fig = plt.figure(figsize=(8,8))
big_number = 10.0**(4.0)
x_min = -2.2*big_number
x_max = 2.2*big_number
binwidth=0.01*big_number
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
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
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(Hist_1D <= 2100.0) & (Hist_1D>=-2100.0)]
print ("After: ",Hist_1D.shape)
#Figure
output_name = output_folder+'/VDH_hist_sim_small.png'
fig = plt.figure(figsize=(8,8))
x_min = -2200.0
x_max = 2200.0
binwidth=5.0
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
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
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
plt.gca().set_yscale("log")
plt.gca().set_xscale("log")
plt.xlabel('VDH',fontsize=18)
#plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)',fontsize=18)
plt.title('VDH distribution ', fontsize=18)
plt.savefig(output_name,format = 'png')
plt.close(fig)
#-------------------------------------------------
"""
print ("Cuts for Density, DensitySM, Tidal and Linear Vel.")

#Density
output_name = output_folder+'/Density_cut_i_'+str(cut_i)+'.png'
cut_dens = dens[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut_dens.T,interpolation='none',origin='lower')
plt.colorbar(label='Density')
plt.xlabel("y [Gpc/h]", fontsize=16)
plt.ylabel("z [Gpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.title('Density cut for \n x = '+str(cut_i_gpc)+' Gpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," Density min: ", cut_dens.min()," max: ", cut_dens.max())
plt.close(fig)

#Density SM
output_name = output_folder+'/Density_SM_cut_i_'+str(cut_i)+'.png'
cut_dens_sm = dens_sm[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut_dens_sm.T, interpolation='none',origin='lower')
plt.colorbar(label='Density SM')
plt.xlabel("y [Gpc/h]", fontsize=16)
plt.ylabel("z [Gpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.title('Density SM cut for \n x = '+str(cut_i_gpc)+' Gpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," Density SM min: ", cut_dens_sm.min()," max: ", cut_dens_sm.max())
plt.close(fig)

#Tidal
output_name = output_folder+'/Tidal_cut_i_'+str(cut_i)+'.png'
cut_tidal = tidal[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut_tidal.T, interpolation='none',origin='lower')
plt.colorbar(label='Tidal')
plt.xlabel("y [Gpc/h]", fontsize=16)
plt.ylabel("z [Gpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.title('Tidal cut for \n x = '+str(cut_i_gpc)+' Gpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," Tidal min: ", cut_tidal.min()," max: ", cut_tidal.max())
plt.close(fig)

#Lineal Velocigy
output_name = output_folder+'/L_Vel_cut_i_'+str(cut_i)+'.png'
cut_lineal_vel = lineal_vel[cut_i,:,:]
fig = plt.figure(figsize=(8,8))
plt.imshow(cut_lineal_vel.T, interpolation='none',origin='lower')
plt.colorbar(label='Linear Velocity')
plt.xlabel("y [Gpc/h]", fontsize=16)
plt.ylabel("z [Gpc/h]", fontsize=16)
plt.xticks(x_thicks,xticklabels, fontsize=14)
plt.yticks(y_thicks,yticklabels, fontsize=14)
plt.title('Linear Velocity cut for \n x = '+str(cut_i_gpc)+' Gpc/h', fontsize=16)
plt.savefig(output_name,format = 'png')
print ("Cut ", cut_i," Linear Velocity min: ", cut_lineal_vel.min()," max: ", cut_lineal_vel.max())
plt.close(fig)

print ("Histograms for Density, DensitySM, Tidal and Linear Vel.")
#Density
Hist_1D = dens.ravel()
print ("Density")
print ("Max, Min")
print (np.max(Hist_1D), np.min(Hist_1D))
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Hist_1D.shape)
#Figure
#Big
output_name = output_folder+'/Density_hist_sim_big.png'
fig = plt.figure(figsize=(8,8))
x_min = np.min(Hist_1D)
x_max = np.max(Hist_1D)
binwidth=0.001
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
plt.gca().set_yscale("log")
#plt.gca().set_xscale("log")
plt.xlabel('Density',fontsize=18)
plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)')
plt.title('Density distribution', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)

#Density_SM
Hist_1D = dens_sm.ravel()
print ("Density SM")
print ("Max, Min")
print (np.max(Hist_1D), np.min(Hist_1D))
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Hist_1D.shape)
#Figure
#Big
output_name = output_folder+'/DensitySM_hist_sim_big.png'
fig = plt.figure(figsize=(8,8))
x_min = np.min(Hist_1D)
x_max = np.max(Hist_1D)
binwidth=0.001
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
plt.gca().set_yscale("log")
#plt.gca().set_xscale("log")
plt.xlabel('Density SM',fontsize=18)
plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)')
plt.title('Density SM distribution', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)

#Tidal
Hist_1D = tidal.ravel()
print ("Tidal")
print ("Max, Min")
print (np.max(Hist_1D), np.min(Hist_1D))
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Hist_1D.shape)
#Figure
#Big
output_name = output_folder+'/Tidal_hist_sim_big.png'
fig = plt.figure(figsize=(8,8))
x_min = np.min(Hist_1D)
x_max = np.max(Hist_1D)
binwidth=0.001
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
plt.gca().set_yscale("log")
#plt.gca().set_xscale("log")
plt.xlabel('Tidal',fontsize=18)
plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)')
plt.title('Tidal distribution', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)

#Lineal Velocity
Hist_1D = lineal_vel.ravel()
print ("Lineal Velocity")
print ("Max, Min")
print (np.max(Hist_1D), np.min(Hist_1D))
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(FA_1D_orig <= 1.0) & (FA_1D_orig>=0.0)]
print ("After: ",Hist_1D.shape)
#Again
print ("Before: ",Hist_1D.shape)
Hist_1D = Hist_1D[(Hist_1D<=1000.0)]
print ("After: ",Hist_1D.shape)
#Figure
#Big
output_name = output_folder+'/LinealVel_hist_sim_big.png'
fig = plt.figure(figsize=(8,8))
x_min = np.min(Hist_1D)
x_max = np.max(Hist_1D)
binwidth=1.0
plt.hist(Hist_1D, bins=np.arange(x_min, x_max + binwidth, binwidth))
#plt.hist(Hist_1D)
plt.gca().set_yscale("log")
#plt.gca().set_xscale("log")
plt.xlabel('Lineal Velocity',fontsize=18)
plt.xlim(x_min,x_max)
plt.ylabel('Number of cells (log scale)')
plt.title('Lineal Velocity distribution', fontsize=16)
plt.savefig(output_name,format = 'png')
plt.close(fig)
"""
