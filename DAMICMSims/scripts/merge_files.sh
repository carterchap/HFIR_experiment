#!/bin/bash

if [ $# -lt 2 ]; then
		echo "Usage:"
		echo "$0 volume_directory output_directory [log_directory]"
		echo "If log_directory is given, will run ToTree on the merged files"
		exit
fi
doToTree=1
if [ $# -gt 2 ]; then
		doToTree=0
		log_dir=`readlink -f $3`
fi
working_dir=$PWD
TOTREE="ToTree.C"

volume_dir=$1
volume=`echo "${volume_dir%/}"| rev|cut -d"/" -f-1 | rev`
if [ ! -d $2 ]; then
		echo "Output directory not found, exiting..."
		exit
fi
out_dir=`readlink -f $2`

echo "{" > compile.C
echo "gROOT->ProcessLine(\".L $TOTREE+\");" >> compile.C
echo "}" >> compile.C
rm -f compile.C
cd $volume_dir
if [ $? -ne 0 ]; then
		echo "Volume not found, exiting..."
		exit
fi
isotopes=(`ls -d */`)

for iso in "${isotopes[@]}"; do
		isotope=${iso%?} #Strips the trailing /
		cd $isotope

		if [ -d "processed" ]; then	
				ofname="${out_dir}/${volume}${isotope}.root"
				hadd -f ${ofname} processed/temp/*.root > /dev/null
				if [ $doToTree -eq 0 ]; then
						TIME=$(date +"%s")
						script_file="${log_dir}/pb_ps_$TIME.sh"
 						echo "#PBS -j oe" > $script_file
						echo "#PBS -o $log_dir/pbs_ps_$TIME.log" >> $script_file
						echo "export PATH=$PATH" >> $script_file
						echo "export PBS_O_PATH=$working_dir" >> $script_file
						echo "export ROOT_HIST=0" >> $script_file
						echo "cd $working_dir" >> $script_file
						echo "root -l -b -q '$TOTREE+(\"$ofname\",\"default.cfg\")'" >> $script_file
						chmod u+x $script_file
						qsub $script_file						
						sleep 1
				fi
		fi
		cd ../
done
