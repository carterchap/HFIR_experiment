#!/bin/bash
usage="$0 directory \nArranges simulation output into a directory tree of volume/isotope/files"
 
if [ $# -lt 1 ]; then
		echo -e $usage;
		exit
fi
#gets absolute path to directory
directory=$1
cd $directory
if [ $? -ne 0 ]; then
		echo "Directory not found, exiting..."
		exit 1;
		echo "Error"
fi
files=(`ls *.root`)
vols=(`ls *.root | grep -Eo ".*PV"`)
unique_vols=(`echo "${vols[@]}" | tr ' ' '\n' | sort -u`)
for volume in "${unique_vols[@]}"; do
		mkdir $volume
		mv ${volume}*.root $volume
		cd $volume
		if [ $? != 0 ]; then
				echo "error encountered, exiting..."
				exit
		fi
		isos=(`ls *.root | grep -Eo "[0-9]{1,3}a[0-9]{1,3}z"`)
		unique_isos=(`echo "${isos[@]}" | tr ' ' '\n' | sort -u`)
		for isotope in "${unique_isos[@]}"; do
				ls *${isotope}*.root >> /dev/null
				if [ $? -eq 0 ]; then
						mkdir $isotope
						mv *${isotope}*.root $isotope
				fi
		done
		cd ../
done
