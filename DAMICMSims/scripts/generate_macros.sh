#!/bin/bash
usage="USAGE: $0 <in_file> <simulation_build_folder>"

if [ $# -lt 2 ]; then
		echo $usage
		exit
fi
in_file=$1
macname=${in_file%.*}
macname=${macname##*/}
out_folder=$2
echo $out_folder
if [ ! -d $out_folder ]; then
		echo "Error, build folder doesn't exist, exiting..."
		exit
fi
declare -a vols #Volume
declare -a atomic_num #Atomic number
declare -a charge #Particle charge
declare -a do_surface #0=dovolume, 1=dosurface (uniform), 2=dosurface (diffusion)
declare -a min_depth
declare -a max_depth
declare -a temp
unset numfiles
while IFS='' read -r line || [[ -n "$line" ]];do
		#Check for comments first, any line that contains a # (not ideal, but easy)
		echo $line | grep "#"  >> /dev/null
		if [ $? -eq 0 ]; then 
				continue
		fi
		nums=()
		for val in $line; do
				nums+=("$val")
		done
		echo $line | grep -E "{[0-9]+:[0-9]+} {[0-9]+:[0-9]+}"
		if [ $? -eq 0 ] && [ -z ${numfiles+x} ]; then
				numfiles=${nums[0]}
				numdecays=${nums[1]}
				continue
		fi
				

		# Loop for isotope lines
		echo $line | grep -E "[0-9]{1,3} [0-9]{1,3}" >> /dev/null
		if [ $? -eq 0 ]; then
				if [ -z ${numfiles+x} ]; then
						numfiles=${nums[0]}
						numdecays=${nums[1]}
				else
						atomic_num+=(${nums[0]})
						charge+=(${nums[1]})
						if [ ${#nums[@]} -gt 3 ]; then
								if [ "${nums[2]}" = "uni" ];then
										do_surface+=(1)
								elif [ "${nums[2]}" = "dif" ]; then
										do_surface+=(2)
								else 
										echo "Warning, surface decay parameter not recognized, doing bulk decays"
										do_surface+=(0)
								fi
								min_depth+=(${nums[3]})
								max_depth+=(${nums[4]})
						else 
								do_surface+=(0)
								min_depth+=(0)
								max_depth+=(0)
						fi
				fi
				continue
		fi
		echo $line | grep PV >> /dev/null
		if [ $? -eq 0 ]; then
				vols+=($line)
				continue
		fi		
done < $in_file
macdir="$out_folder/macros"
outdir="$out_folder/outmac"
mkdir $macdir &> /dev/null
mkdir $macdir/done &> /dev/null
mkdir $outdir &> /dev/null
mkdir $out_folder/logs &> /dev/null

let "totalfiles=$numfiles * ${#atomic_num[@]} * ${#vols[@]}"
# echo $totalfiles
# echo ${atomic_num[@]}
filenum=0
for vol in "${!vols[@]}"; do
		for atomic in "${!atomic_num[@]}"; do
				for j in `seq 1 $numfiles`; do

						#Certain special decays equire setting the excitation energy of the state
						if [ ${charge[atomic]} -eq "91" ] && [ ${atomic_num[atomic]} -eq "234" ]; then
								excite_energy=73.92
						elif [ ${charge[atomic]} -eq "47" ] && [ ${atomic_num[atomic]} -eq "110" ]; then
                excite_energy=117.59
						else
								excite_energy=0
						fi

						#Surface decays using advanced diffusion model
						if [ "${do_surface[atomic]}" -eq "2" ]; then
								vol_command="dosurface \n /damic/gun/position/setadvdiffmodel true \n 
/damic/gun/position/setmaxembeddist ${max_depth[atomic]} nm"
								#surface decays using standard uniform embedding
						elif [ "${do_surface[atomic]}" -eq "1" ]; then
								vol_command="dosurface \n /damic/gun/position/setminembeddist ${min_depth[atomic]} nm \n 
/damic/gun/position/setmaxembeddist ${max_depth[atomic]} nm"
						else
								#Otherwise do a bulk decay
								vol_command="dovolume"
						fi

						let "filenum=$filenum+1"
						output="/control/verbose 0 \n
/run/verbose 0 \n
/run/initialize \n
/tracking/verbose 0 \n
/event/verbose 0 \n
/random/setSeeds $RANDOM $RANDOM \n
/damic/gun/particle ion \n 
/damic/gun/ion ${charge[atomic]} ${atomic_num[atomic]} 0 ${excite_energy} \n
/damic/gun/energy/mono 0 eV \n 
/damic/gun/direction/oned \n 
/damic/gun/direction/onedX 1\n
/damic/gun/direction/onedY 0\n
/damic/gun/direction/onedZ 0\n
/damic/gun/position/${vol_command}\n
/damic/gun/position/addvolume ${vols[vol]} 1\n
/analysis/setFileName $outdir/${vols[vol]}${atomic_num[atomic]}a${charge[atomic]}z$j\n
/run/beamOn $numdecays"
	
			 			echo -e $output > "$macdir/${macname}$filenum.mac"
				done
		done
done

#now write our pbs script


#check what system we're on
declare -a hn
hn=$(hostname)

if [ -z "${hn##*zev*}" ]; then
		fname="$out_folder/g4_pbs_${macname}.sh"
		echo "Creating pbs script for zev"
		echo "#!/bin/sh" > $fname
		echo "#PBS -o $out_folder/logs/g4_pbs_${macname}.log" >> $fname
		echo "#PBS -j oe" >> $fname
		echo "#PBS -l cput=600:00:00,walltime=100:00:00,ncpus=55" >> $fname
		echo "#PBS -d $out_folder" >> $fname
		echo "#PBS -t 1-$filenum%1" >> $fname
		echo "source $out_folder/geant4_env_10.2" >> $fname
		echo "$out_folder/DAMICG4 macros/${macname}\${PBS_ARRAYID}.mac 55" >> $fname
		echo "mv macros/${macname}\${PBS_ARRAYID}.mac macros/done/" >> $fname

		chgrp damic $fname
		chmod g+w $fname

elif [ -z "${hn##*midway*}" ]; then
		fname="$out_folder/g4_slurm_${macname}.sh"
		echo "Creating Slurm script for midway"
		echo "#!/bin/sh" > $fname
		echo "#SBATCH --workdir=$out_folder" >> $fname
		echo "#SBATCH --array=1-$filenum" >> $fname
		echo "#SBATCH --time=06:00:00" >> $fname
		echo "#SBATCH --partition=sandyb" >> $fname
		echo "#SBATCH --ntasks=1" >> $fname
		echo "#SBATCH --cpus-per-task=16" >> $fname
		echo "#SBATCH --output=$out_folder/logs/g4_slurm_${macname}_%A_%a.log" >> $fname
		echo "source /etc/profile" >> $fname
		echo "source $out_folder/geant4_env_10.2" >> $fname
		echo "$out_folder/DAMICG4 macros/${macname}\${SLURM_ARRAY_TASK_ID}.mac 16" >> $fname
		echo "mv macros/${macname}\${SLURM_ARRAY_TASK_ID}.mac macros/done/" >> $fname		

elif [ -z "${hn##*sdu*}" ]; then
		fname="$out_folder/g4_slurm_${macname}.sh"
		echo "Creating Slurm script for Abacus"
		echo "#!/bin/sh" > $fname
		echo "#SBATCH --workdir=$out_folder" >> $fname
		echo "#SBATCH --array=1-$filenum" >> $fname
		echo "#SBATCH --time=06:00:00" >> $fname
		echo "#SBATCH --ntasks=1" >> $fname
		echo "#SBATCH --cpus-per-task=24" >> $fname
		echo "#SBATCH --output=$out_folder/logs/g4_slurm_${macname}_%A_%a.log" >> $fname
		echo "source /etc/profile" >> $fname
		echo "source $out_folder/geant4_env_10.2" >> $fname
		echo "$out_folder/DAMICG4 macros/${macname}\${SLURM_ARRAY_TASK_ID}.mac 24" >> $fname
		echo "mv macros/${macname}\${SLURM_ARRAY_TASK_ID}.mac macros/done/" >> $fname		

		
elif [ -z "${hn##*altamira*}" ]; then
		fname="$out_folder/g4_slurm_${macname}.sh"
		echo "Creating Slurm script for altamira"
		echo "#!/bin/sh" > $fname		
		echo "#@ job_name = $macname " >> $fname
		echo "$@ initialdir = $out_folder " >> $fname
		echo "$@ output =$out_folder/logs/g4_slurm_${macname}_%j.log" >> $fname
		echo "$@ wall_clock_limit=06:00:00" >> $fname
		echo "$@ total_tasks = $filenum" >> $fname
		echo "#@ cpus_per_task = 16" >> $fname
		echo "source /etc/profile" >> $fname
		echo "source $out_folder/geant4_env_10.2" >> $fname
		echo "$out_folder/DAMICG4 macros/${macname}\${SLURM_ARRAY_TASK_ID}.mac 16" >> $fname
		echo "mv macros/${macname}\${SLURM_ARRAY_TASK_ID}.mac macros/done/" >> $fname		

else
		echo "Not sure what system we're on, not creating run script..."
fi
chmod a+x $fname
