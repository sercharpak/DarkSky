#!/usr/bin/python

#Plotting the general properties of the detected groups
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
USAGE = "python pipeline_histograms.py file_groups scale_factor n_grid output_folder"
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

def file_len(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
    return i+1
#Returns the name of the file given a path
def path_leaf(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#Function to analyze the Inertia Tensor
#Based in http://scipython.com/book/chapter-6-numpy/problems/p65/the-moment-of-inertia-tensor/
def inertia_tensor_case (I_a, I_b, I_c):
    case = 0
    I_sum = I_a+I_b+I_c
    I_b_a = (I_b - I_a)/I_sum
    I_c_b = (I_c - I_b)/I_sum
    #Case 1 Ia=Ib=Ic: spherical top;
    if((I_a==I_b==I_c) or ((I_b_a<0.05) and (I_c_b<0.05))):
        case = 1
        return case
    #Case 2 Ia=Ib<Ic: oblate symmetric top;
    elif((I_a==I_b < I_c) or ((I_b_a<0.05) and (I_c_b>0.05))):
        case = 2
        return case
    #Case 3 Ia<Ib=Ic: prolate symmetric top;
    elif((I_a < I_b == I_c) or ((I_b_a>0.05) and (I_c_b<0.05))):
        case = 3
        return case
    #Case 4 Ia<Ib<Ic: asymmetric top
    elif((I_a < I_b < I_c) or ((I_b_a>0.05) and (I_c_b>0.05))):
        case = 4
        return case
    #Case 0 - ?
    else:
        case = 0
        return case
#---------------------------------------------------------
# script
#---------------------------------------------------------
if(len(sys.argv)!=5):
    print ("Please use correctly")
    print (USAGE)
    sys.exit()
inputfile = sys.argv[1]
SCALE_FACTOR = float(sys.argv[2])
SCALE_FACTOR_Mpc = SCALE_FACTOR/1000.0
SCALE_FACTOR_Gpc = SCALE_FACTOR_Mpc/1000.0
GRID = int(sys.argv[3])
VOLUME = GRID ** 3.0
VOLUME_Kpc = VOLUME * (SCALE_FACTOR ** 3.0)
VOLUME_Mpc = VOLUME * (SCALE_FACTOR_Mpc ** 3.0)
VOLUME_Gpc = VOLUME * (SCALE_FACTOR_Gpc ** 3.0)
output_folder = sys.argv[4]
print ("Processing the file: " + inputfile)
filename_groups = path_leaf(inputfile)
thresh_search_FA, thresh_search_Trace, thresh_seeds_FA, thresh_seeds_Trace = get_thresholds(filename_groups)
group_data = np.loadtxt(inputfile)
n_groups = file_len(inputfile)
print ("There are", n_groups,"groups")
print ("The grid size is: ", GRID)
print ("The scale factor is: ", SCALE_FACTOR)
print ("The real boxsize is [kpc/h]:", GRID*SCALE_FACTOR)
print (" or in [Mpc/h]: ", (GRID*SCALE_FACTOR_Mpc))
print (" or in [Gpc/h]: ", (GRID*SCALE_FACTOR_Gpc))
#-------------------------------------------
#Volume distribution
if(n_groups==1):
    volumes = group_data[7]
else:
    volumes = group_data[:,7]
print ("Volume distribution for \n growth: FA: " + str(thresh_search_FA) + " Trace: " + str(thresh_search_Trace) +"\n seeds: FA: " + str(thresh_seeds_FA) + " Trace: " + str(thresh_seeds_Trace))
#In kpc^3/h
volumes_kpc = volumes * (SCALE_FACTOR**3.0)
fig = plt.figure(figsize = (8,8))
binwidth=0.05
if(n_groups==1):
    width = 0.35
    plt.bar(np.log10(volumes_kpc)+width, width)
    plt.xlim(np.min(np.log10(volumes_kpc))-(binwidth/10.0), np.max(np.log10(volumes_kpc))+(binwidth/10.0))
else:
    plt.xlim(np.min(np.log10(volumes_kpc))-(binwidth/10.0), np.max(np.log10(volumes_kpc))+(binwidth/10.0))
    plt.hist(np.log10(volumes_kpc), bins=np.arange(min(np.log10(volumes_kpc)), max(np.log10(volumes_kpc)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume $(kpc/h)^3$", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.title("Distribution of volume $(kpc/h)^3$ with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/volumes_distr_kpc_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

#In Mpc^3/h
volumes_mpc = volumes_kpc * ((1.0/1000.0)**3.0)
fig = plt.figure(figsize = (8,8))
binwidth=0.05
if(n_groups==1):
    width = 0.35
    plt.bar(np.log10(volumes_mpc)+width, width)
    plt.xlim(np.min(np.log10(volumes_mpc))-(binwidth/10.0), np.max(np.log10(volumes_mpc))+(binwidth/10.0))
else:
    plt.xlim(np.min(np.log10(volumes_mpc))-(binwidth/10.0), np.max(np.log10(volumes_mpc))+(binwidth/10.0))
    plt.hist(np.log10(volumes_mpc), bins=np.arange(min(np.log10(volumes_mpc)), max(np.log10(volumes_mpc)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume $(Mpc/h)^3$", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.title("Distribution of volume $(Mpc/h)^3$ with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/volumes_distr_Mpc_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

#In Gpc^3/h
volumes_gpc = volumes_mpc * ((1.0/1000.0)**3.0)
fig = plt.figure(figsize = (8,8))
binwidth=0.05
if(n_groups==1):
    width = 0.35
    plt.bar(np.log10(volumes_gpc)+width, width)
    plt.xlim(np.min(np.log10(volumes_gpc))-(binwidth/10.0), np.max(np.log10(volumes_gpc))+(binwidth/10.0))
else:
    plt.xlim(np.min(np.log10(volumes_gpc))-(binwidth/10.0), np.max(np.log10(volumes_gpc))+(binwidth/10.0))
    plt.hist(np.log10(volumes_gpc), bins=np.arange(min(np.log10(volumes_gpc)), max(np.log10(volumes_gpc)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume $(Gpc/h)^3$", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.title("Distribution of volume $(Gpc/h)^3$ with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/volumes_distr_Gpc_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)
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
log_10_vol_laniakea_Gpc = np.log10(volume_laniakea_Gpc)
#In kpc^3/h
fig = plt.figure(figsize = (8,8))
plt.axvline(log_10_vol_laniakea_kpc, linewidth=2, color='r', label='Laniakea')
binwidth=0.05
if(n_groups==1):
    width = 0.35
    plt.bar(np.log10(volumes_kpc)+width, width)
    #plt.xlim(np.min((np.min(np.log10(volumes_kpc))-(binwidth/10.0)), (log_10_vol_laniakea_kpc-(binwidth/10.0))), np.max((np.max(np.log10(volumes_kpc))+(binwidth/10.0)), (log_10_vol_laniakea_kpc+(binwidth/10.0))))
else:
    #plt.xlim(np.min((np.min(np.log10(volumes_kpc))-(binwidth/10.0)), (log_10_vol_laniakea_kpc-(binwidth/10.0))), np.max((np.max(np.log10(volumes_kpc))+(binwidth/10.0)), (log_10_vol_laniakea_kpc+(binwidth/10.0))))
    plt.hist(np.log10(volumes_kpc), bins=np.arange(min(np.log10(volumes_kpc)), max(np.log10(volumes_kpc)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume $(kpc/h)^3$", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.legend(loc='upper right', fontsize='20')
plt.title("Distribution of volume $(kpc/h)^3$ with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/volumes_distr_kpc_laniakea_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

#In Mpc^3/h
fig = plt.figure(figsize = (8,8))
plt.axvline(log_10_vol_laniakea_Mpc, linewidth=2, color='r', label='Laniakea')
binwidth=0.05
if(n_groups==1):
    width = 0.35
    plt.bar(np.log10(volumes_mpc)+width, width)
    #plt.xlim(np.min((np.min(np.log10(volumes_mpc))-(binwidth/10.0)), (log_10_vol_laniakea_Mpc-(binwidth/10.0))), np.max((np.max(np.log10(volumes_mpc))+(binwidth/10.0)), (log_10_vol_laniakea_Mpc+(binwidth/10.0))))
else:
    #plt.xlim(np.min((np.min(np.log10(volumes_mpc))-(binwidth/10.0)), (log_10_vol_laniakea_Mpc-(binwidth/10.0))), np.max((np.max(np.log10(volumes_mpc))+(binwidth/10.0)), (log_10_vol_laniakea_Mpc+(binwidth/10.0))))
    plt.hist(np.log10(volumes_mpc), bins=np.arange(min(np.log10(volumes_mpc)), max(np.log10(volumes_mpc)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume $(Mpc/h)^3$", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.legend(loc='upper right', fontsize='20')
plt.title("Distribution of volume $(Mpc/h)^3$ with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/volumes_distr_Mpc_laniakea_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

#In Gpc^3/h
fig = plt.figure(figsize = (8,8))
plt.axvline(log_10_vol_laniakea_Gpc, linewidth=2, color='r', label='Laniakea')
binwidth=0.05
if(n_groups==1):
    width = 0.35
    plt.bar(np.log10(volumes_gpc)+width, width)
    #plt.xlim(np.min((np.min(np.log10(volumes_gpc))-(binwidth/10.0)), (log_10_vol_laniakea_Mpc-(binwidth/10.0))), np.max((np.max(np.log10(volumes_gpc))+(binwidth/10.0)), (log_10_vol_laniakea_Mpc+(binwidth/10.0))))
else:
    #plt.xlim(np.min((np.min(np.log10(volumes_gpc))-(binwidth/10.0)), (log_10_vol_laniakea_Mpc-(binwidth/10.0))), np.max((np.max(np.log10(volumes_gpc))+(binwidth/10.0)), (log_10_vol_laniakea_Mpc+(binwidth/10.0))))
    plt.hist(np.log10(volumes_gpc), bins=np.arange(min(np.log10(volumes_gpc)), max(np.log10(volumes_gpc)) + binwidth, binwidth))
plt.xlabel("$log_{10}$ of Volume $(Gpc/h)^3$", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.grid()
plt.legend(loc='upper right', fontsize='20')
plt.title("Distribution of volume $(Gpc/h)^3$ with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/volumes_distr_Gpc_laniakea_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

#-------------------------------------------
#Inertia Tensor Analysis
#Constants for writting the result of the shape of the respective detected group
CASES = []
CASE_0 = "none of the above"
CASES.append(CASE_0)
CASE_1 = "spherical top"
CASES.append(CASE_1)
CASE_2 = "oblate symmetric top"
CASES.append(CASE_2)
CASE_3 = "prolate symmetric top"
CASES.append(CASE_3)
CASE_4 = "asymmetric top"
CASES.append(CASE_4)
if(n_groups==1):
    lambdas_a = group_data[0]
    lambdas_b = group_data[1]
    lambdas_c = group_data[2]
else:
    lambdas_a = group_data[:,0]
    lambdas_b = group_data[:,1]
    lambdas_c = group_data[:,2]

#This counting could be done way better
inertia_cases = np.zeros(5)

if(n_groups==1):
    cases = inertia_tensor_case(lambdas_a, lambdas_b, lambdas_c)
    inertia_cases[int(cases)] += 1
else:
    cases = np.zeros(n_groups)
    for i in range(n_groups):
        cases[i] = inertia_tensor_case(lambdas_a[i], lambdas_b[i], lambdas_c[i])
        inertia_cases[int(cases[i])] += 1

print ("Inertia Analysis for \n growth: FA: " + str(thresh_search_FA) + " Trace: " + str(thresh_search_Trace) +"\n seeds: FA: " + str(thresh_seeds_FA) + " Trace: " + str(thresh_seeds_Trace))
fig = plt.figure(figsize = (8,8))
if(n_groups==1):
    width = 0.35
    plt.bar(cases+width, width)
else:
    plt.hist(cases)
plt.grid()
plt.xlim(0,4.1)
plt.xticks([0.0,1.0,2.0,3.0,4.0])
plt.xlabel("Inertia Case", fontsize=20)
plt.ylabel("Number of Groups", fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.title("Inertia Cases with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/inertia_cases_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)

print ("-------------------------------------------")
print ("Case", "Case Number", "Number of occurences")
for i in range (len(inertia_cases)):
    print (str(CASES[i]), str(i), str(inertia_cases[i]))

#EigenValues differences
#la - lb vs. lb - lc
#normalized by (la+lb+lc)

diff_b_a = lambdas_b - lambdas_a
diff_b_a = diff_b_a /(lambdas_a +lambdas_b + lambdas_c)

diff_c_b = lambdas_c - lambdas_b
diff_c_b = diff_c_b /(lambdas_a +lambdas_b + lambdas_c)

fig = plt.figure(figsize = (10,10))
plt.axvline(0, linewidth=2, color='r')
plt.axhline(0, linewidth=2, color='r')

plt.scatter(diff_c_b, diff_b_a)
plt.xlabel("$\lambda_b$ - $\lambda_a$ / ($\lambda_a$ + $\lambda_b$ + $\lambda_c$)", fontsize=20)
plt.ylabel("$\lambda_c$ - $\lambda_b$ / ($\lambda_a$ + $\lambda_b$ + $\lambda_c$)", fontsize=20)
plt.grid()
plt.tick_params(axis='both', which='major', labelsize=15)
plt.title("Inertia Tensor Eigenvalues differences with "+str(n_groups) + " groups \n Growth: FA: " + str(thresh_search_FA) + " Divergence: " + str(thresh_search_Trace) +"\n Seeds: FA: " + str(thresh_seeds_FA) + " Divergence: " + str(thresh_seeds_Trace), fontsize=16)
plt.savefig(output_folder+"/inertia_diff_"+str(thresh_seeds_FA)+"_Trace_"+str(thresh_seeds_Trace)+"_search_FA_"+str(thresh_search_FA)+"_Trace_"+str(thresh_search_Trace)+".png",format = "png")
plt.close(fig)
