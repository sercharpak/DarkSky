#!/bin/bash

#Script to Excecute the FOF.
#Written by Sergio Hernandez Charpak

#Example: ./script_awk_bash.sh /home/laptop/Documentos/COLAB/skydata/grids_large_halos/256_15B/FoF_Results/grid_FA_0.9_Trace_0.0.dat /home/laptop/Documentos/COLAB/skydata/grids_large_halos/FoF/src /home/laptop/Documentos/COLAB/skydata/grids_large_halos/256_15B/FoF_Results
if [ $# -ne 3 ]
then
    echo "Usage: $(basename $0) file_search fof_src_folder output_folder"
    exit 1
fi

file_search=$1
fof_src_folder=$2
output_folder=$3

fof_exe="$fof_src_folder/fof"

$fof_exe -e 1.1 -m 20 < "$file_search"

mv fof.grp "$output_folder"

#Excecutes the rest. Should be a Makefile
cd 0509
make
cd ../0609
make
cd ../0709
make
cd ../0809
make
cd ../0909
make
