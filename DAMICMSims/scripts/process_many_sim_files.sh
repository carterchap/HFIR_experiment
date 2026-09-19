#!/bin/bash

if [ $# -lt 2 ]; then
		echo -e "Usage:"
		echo -e "$0 volume_directory log_directory (max_tasks) (cfg_file)\n Note that if max_tasks is set, damic1k_sim.cfg file will be used by default for configuration"
		exit
fi

volume=`readlink -f $1`
if [ ! -d $volume ]; then
		echo "Volume not found, exiting..."
		exit
fi

log_dir=`readlink -f $2`
if [ ! -d $log_dir ]; then
		echo "Log directory not valid, trying to find directory automatically..."
		log_dir=$volume/../../logs
		if [ ! -d $log_dir ]; then
				echo "Cannot find log directory, exiting..."
				exit
		fi
fi

cfg_file="default.cfg"
if [ $# -gt 2 ]; then
		max_tasks=$3
		cfg_file="damic1k_sim.cfg"
else 
		max_tasks=100
fi
if [ $# -gt 3 ]; then
		cfg_file=$4
fi


#Code for midway/other clusters
declare -a hn
hn=$(hostname)
#Check if we're on Midway
if [ -z "${hn##*midway*}" ]; then
		source /etc/profile
		module load ROOT
fi

#Create arrays to hold names of jobs
declare -a ofname_array
declare -a ifname_array
declare -a jobid_array

#Scripts to load
SCRIPT="ReadSimFast.c"
TOTREE="ToTree.C"
working_dir=$PWD

cd $volume
isotopes=(`ls -d */`)

for iso in "${isotopes[@]}"; do
		isotope=${iso%?} #Strips the trailing /
#		echo "Isotope: $isotope"
		cd $isotope
		mkdir processed
		mkdir processed/temp
#		shopt -s nullglob
		#note: if we do more than 10 files, this will need to be changed
		#maybe to *z[1-9]*.root or something
		base_files=(`ls *z[1-9].root`)
		shopt -u nullglob
		for file in "${base_files[@]}"; do
				ifname="$volume/$isotope/${file%.*}*.root"
#				echo $ifname
				ofname="$volume/$isotope/processed/temp/$file"
				ifname_array+=("$ifname")
				ofname_array+=("$ofname")
#				echo $ofname
		done
		cd ../
done

cd $working_dir
#compile our script
echo "{" > compile.C
echo "gROOT->ProcessLine(\".L $SCRIPT+\");" >> compile.C
echo "gROOT->ProcessLine(\".L $TOTREE+\");" >> compile.C
echo "}" >> compile.C

root -l -b -q compile.C
if [ -$? -ne 0 ]; then
		echo "Script compilation error, exiting..."
		exit
fi
sleep 1
rm -f compile.C


for i in "${!ofname_array[@]}"; do
		ifile=${ifname_array[$i]}
		ofile=${ofname_array[$i]}

		TIME=$(date +"%s%N")
		script_file="${log_dir}/pb_ps_$TIME.sh"
		echo "#!/bin/sh" > $script_file

		if [ -z "${hn##*midway*}" ]; then
				echo "#SBATCH --workdir=$log_dir" >> $script_file
#				echo "#SBATCH --array=1-$filenum" >> $script_file
				echo "#SBATCH --time=10:00:00" >> $script_file
				echo "#SBATCH --partition=sandyb" >> $script_file
#				echo "#SBATCH --ntasks=1" >> $script_file
#				echo "#SBATCH --cpus-per-task=16" >> $script_file
				echo "#SBATCH --output=${log_dir}/pbs_ps_$TIME.log" >> $script_file
				echo "#SBATCH --mem=5G" >> $script_file
				echo "source /etc/profile" >> $script_file
				echo "module load ROOT" >> $script_file
		else				
 				echo "#PBS -j oe" >> $script_file
				echo "#PBS -o $log_dir/pbs_ps_$TIME.log" >> $script_file
				echo "#PBS -l cput=10:00:00,walltime=10:00:00,mem=10gb" >> $script_file
		fi
		echo "export PATH=$PATH" >> $script_file
		#		echo "source ~/.profile" >> $script_file
		echo "export PBS_O_PATH=$PWD" >> $script_file
		echo "export ROOT_HIST=0" >> $script_file
		echo "cd $PWD" >> $script_file
		echo "root -l -b -q '$SCRIPT+(\"$cfg_file\",\"$ifile\",\"$ofile\")'" >> $script_file
		#Note that totree can fail in some cases, so this may need to be commented out
		echo "root -l -b -q '$TOTREE+(\"$ofile\",\"$cfg_file\")'" >> $script_file
		chmod u+x $script_file
		if [ "${#jobid_array[@]}" -lt $max_tasks ]; then
				jobid_array+=(`qsub $script_file`)
				echo ${jobid_array[$i]}
		else
				jobid_array+=(`qsub -W depend=afterany:"${jobid_array[$i-$max_tasks]%%.*}" $script_file`)
				echo ${jobid_array[$i]}
		fi		
		sleep .1
done
