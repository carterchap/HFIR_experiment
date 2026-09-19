#!/bin/bash

usage="$0 simulation_output\nOutputs a moderately pretty table of which isotopes in which volumes have been processed\nNote: only lists isotopes hard coded into the script"
isotopes=(208a81z 210a82z 210a83z 212a82z 212a83z 214a82z 214a83z 228a88z 228a89z 234a90z 234a91z 40a19z 54a25z 56a27z 57a27z 58a27z 59a26z 60a27z 87a37z)

if [ -$# -eq 0 ]; then 
		echo -e $usage
		exit
fi

echo > temp.txt
volume_dir=`readlink -f $1`
volumes=(`ls $volume_dir | grep -E "*PV"`)
for volume in "${volumes[@]}"; do
		status=()
		declare -a status
		for i in "${!isotopes[@]}"; do
				ls $volume_dir/$volume/*${isotopes[$i]}* >> /dev/null 2>&1
				if [ $? -ne 0 ]; then
						#						color='\033[0;31m'
						color='\e[0;31m'
				else
						color='\e[0;32m'
				fi
				status+=("${color}${isotopes[$i]}\033[0m")
		done
		echo -e "${volume}: ${status[@]}" >> temp.txt
		#printf "%15s %b" ${volume} "${status[@]}"
done

cat temp.txt | column -t
rm temp.txt
