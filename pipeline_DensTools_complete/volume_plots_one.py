#!/usr/bin/python

#Plotting the volumes of the detected groups vs Laniakea-s volume
#By Sergio Daniel Hernandez Charpak

#-------------------------------------------
#Imports
import numpy as np
import matplotlib.pyplot as plt
import sys
#Getting the filename
import ntpath
#-------------------------------------------
# Constants
#-------------------------------------------
USAGE = "python pipeline_histograms.py input_folder thresh_VDH scale_factor n_grid output_folder"
#-------------------------------------------
#Functions
#-------------------------------------------
#Function to extract the thresholds from a group filename
def get_thresholds(filename):
    file_Dat = filename.strip(".dat")
    file_array = file_Dat.split("_")
    thresh_search_FA = float(file_array[9])
    thresh_search_Trace = float(file_array[11])
    thresh_seeds_FA = float(file_array[4])
    thresh_seeds_Trace = float(file_array[6])
    return thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace
#Function to get the length of a file
def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1
#Returns the name of the file given a path
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)
#---------------------------------------------------------
# script
#---------------------------------------------------------
if(len(sys.argv)!=6):
    print ("Please use correctly")
    print (USAGE)
    sys.exit()

input_folder = sys.argv[1]
thresh_VDH = float(sys.argv[2])
SCALE_FACTOR = float(sys.argv[3])
SCALE_FACTOR_Mpc = SCALE_FACTOR/1000.0
SCALE_FACTOR_Gpc = SCALE_FACTOR_Mpc/1000.0
GRID = int(sys.argv[4])
VOLUME = GRID ** 3.0
VOLUME_Kpc = VOLUME * (SCALE_FACTOR ** 3.0)
VOLUME_Mpc = VOLUME * (SCALE_FACTOR_Mpc ** 3.0)
VOLUME_Gpc = VOLUME * (SCALE_FACTOR_Gpc ** 3.0)
output_folder = sys.argv[5]

#Seed 5, 6, 7, 8, 9 Growth 9
array_seeds_FA = np.array([0.5,0.6,0.7,0.8,0.9])
#group_results_seeds_FA_0.6_Trace_1.0_search_FA_0.9_Trace_0.0.dat


#with Laniakea
#Radius Laniakea: 80 Mpc/h
#Volume Laniakea: 2.14E6(Mpc/h)^3
radius_laniakea_kpc = 80000.0
volume_laniakea_kpc = 2.14 * (10.0**15.0)
log_10_vol_laniakea_kpc = np.log10(volume_laniakea_kpc)
radius_laniakea_Mpc = 80.0
volume_laniakea_Mpc = 2.14 * (10.0**6.0)
log_10_vol_laniakea_Mpc = np.log10(volume_laniakea_Mpc)
radius_laniakea_Gpc = 80.0/1000.0
volume_laniakea_Gpc = 2.14 * (10.0**(-3.0))
log_10_vol_laniakea_gpc = np.log10(volume_laniakea_Gpc)

volume_sim_mpc = VOLUME_Mpc
volume_sim_gpc = VOLUME_Gpc
log_10_volume_sim_mpc = np.log10(volume_sim_mpc)
log_10_volume_sim_gpc = np.log10(volume_sim_gpc)

# ---------------------Mpc---------------------
print ("MPC scale")
fig = plt.figure(figsize = (14,14))
plt.axvline(log_10_vol_laniakea_Mpc, linewidth=4, color='m', label='Laniakea')
plt.axvline(log_10_volume_sim_mpc, linewidth=4, color='chocolate', label='Volume Simulation')
max_vol = 0
for thresh_seeds_FA in array_seeds_FA:

    inputfile = input_folder+"/group_results_seeds_FA_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_VDH)+"_search_FA_0.9_Trace_0.0.dat"
    print ("Processing the file: " + inputfile)
    filename = path_leaf(inputfile)
    thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace = get_thresholds(filename)
    group_data = np.loadtxt(inputfile)
    n_groups = file_len(inputfile)
    print ("There are", n_groups,"groups")
    #---------------------------Volume distribution----------------
    if(n_groups==1):
        volumes = group_data[7]
    else:
        volumes = group_data[:,7]
    volumes = volumes * (SCALE_FACTOR_Mpc**3.0)
    bins=100
    #plt.hist(np.log10(volumes), bins=bins, histtype='step',linewidth=3,label = "FA = "+str(thresh_seeds_FA))
    hist_vol, bins = np.histogram(np.log10(volumes),bins=bins)
    print ("Bins")
    print (bins)
    n_bins = np.zeros(len(hist_vol))
    for i in range(len(hist_vol)):
        n_bins[i] = (bins[i]+bins[i+1])/2.0
    #hist_array = np.array(hist_vol)
    print ("hist_vol")
    print (hist_vol)
    log_10_hist_vol = np.log10(hist_vol)
    print ("log 10 hist_vol")
    print (log_10_hist_vol)
    where_are_NaNs = np.isinf(log_10_hist_vol)
    log_10_hist_vol[where_are_NaNs] = 0
    print ("log 10 hist_vol no Nans")
    print (log_10_hist_vol)
    plt.step( n_bins,log_10_hist_vol,linewidth=4, label = "FA = "+str(thresh_seeds_FA) )
    if(np.max(log_10_hist_vol)>max_vol):
        max_vol=np.max(log_10_hist_vol)
#plt.yticks(range(0, int(max_vol)+1, 10))
plt.xlabel("$log_{10}$ of Volume $(Mpc/h)^3$", fontsize=25)
plt.ylabel("$log_{10}$ of Number of Groups", fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=25)
plt.grid()
plt.legend(loc='upper right', fontsize='25')
plt.title("Distribution of volume $(Mpc/h)^3$", fontsize=28)
plt.savefig(output_folder+"/volumes_distr_Mpc_laniakea_all_plot.png",format = "png")
#plt.show()
plt.close(fig)
# ---------------------Gpc---------------------
print ("GPC scale")
fig = plt.figure(figsize = (14,14))
plt.axvline(log_10_vol_laniakea_gpc, linewidth=4, color='m', label='Laniakea')
plt.axvline(log_10_volume_sim_gpc, linewidth=4, color='chocolate', label='Volume Simulation')
max_vol = 0
for thresh_seeds_FA in array_seeds_FA:
    inputfile = input_folder+"/group_results_seeds_FA_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_VDH)+"_search_FA_0.9_Trace_0.0.dat"
    print ("Processing the file: " + inputfile)
    filename = path_leaf(inputfile)
    thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace = get_thresholds(filename)
    group_data = np.loadtxt(inputfile)
    n_groups = file_len(inputfile)
    print ("There are", n_groups,"groups")
    if(n_groups==1):
        volumes = group_data[7]
    else:
        volumes = group_data[:,7]
    volumes = volumes * (SCALE_FACTOR_Gpc**3.0)
    bins=50
    hist_vol, bins = np.histogram(np.log10(volumes),bins=bins)
    print ("Bins")
    print (bins)
    n_bins = np.zeros(len(hist_vol))
    for i in range(len(hist_vol)):
        n_bins[i] = (bins[i]+bins[i+1])/2.0
    print ("hist_vol")
    print (hist_vol)
    log_10_hist_vol = np.log10(hist_vol)
    print ("log 10 hist_vol")
    print (log_10_hist_vol)
    where_are_NaNs = np.isinf(log_10_hist_vol)
    log_10_hist_vol[where_are_NaNs] = 0
    print ("log 10 hist_vol no Nans")
    print (log_10_hist_vol)
    plt.step( n_bins,log_10_hist_vol,linewidth=3, label = "FA = "+str(thresh_seeds_FA) )
    if(np.max(log_10_hist_vol)>max_vol):
        max_vol=np.max(log_10_hist_vol)
#y_thicks = range(0, int(max_vol)+1, 5)
#plt.yticks(y_thicks)
plt.xlabel("$log_{10}$ of Volume $(Gpc/h)^3$", fontsize=25)
plt.ylabel("$log_{10}$ of Number of Groups", fontsize=25)
plt.tick_params(axis='both', which='major', labelsize=25)
plt.grid()
plt.legend(loc='upper right', fontsize='25')
plt.title("Distribution of volume $(Gpc/h)^3$", fontsize=25)
plt.savefig(output_folder+"/volumes_distr_Gpc_laniakea_all_plot.png",format = "png")
#plt.show()
plt.close(fig)
