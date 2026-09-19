#!/bin/bash
tagged=1329

#Script to pull a DAMIC100 simulation from svn on gev.uchicago.edu and compile it
#If you encounter any issues, contact Ryan Thomas (ryant@uchicago.edu) or Karthik Ramanathan (ramanathan@uchicago.edu)

usage="Usage: $0 working_directory [revision] \nCreates build and source directories in working_directory, checks out latest DAMICSIM code from SVN, and compiles simulation\n 'revision' specifies SVN revision"

if [ "$#" -lt 1 ]; then
		echo -e $usage
		exit
fi

if [ "$#" -gt 1 ]; then
		revision=$2;
else
		revision=$tagged
fi

directory=`readlink -f $1`
if [ -f /data/simulation/geant4_env_10.2 ]; then #For zev (Chicago) cluster
		geant4_env=/data/simulation/geant4_env_10.2
elif [ -f /project2/priviter/damic/simulation/builds/geant4_env_10.2 ]; then #for Midway (Chicago) cluster
		geant4_env=/project2/priviter/damic/simulation/builds/geant4_env_10.2
		source /etc/profile
		#TODO: move this to a better place
elif [ -f /gpfs/res_home/damic/rthomas/geant4_env_10.2 ]; then #For Altamira (Spanish) cluster
		geant4_env=/gpfs/res_home/damic/rthomas/geant4_env_10.2
		source /etc/profile
elif [ -f /gpfs/gss1/work/sdudamic/geant4_env_10.2 ]; then #For Abacus (Danish) cluster
		geant4_env=/gpfs/gss1/work/sdudamic/geant4_env_10.2
		source /etc/profile
else
		echo "Error, geant4_env not found"
		echo "Either add your server or contact Ryan Thomas for help"
		exit
fi
source $geant4_env
echo "Using source from SVN to build simulation working directory in: $directory"
echo -e "Enter y to proceed:"
read line
if [ "$line" = y ]; then
		mkdir $directory
		if [ $? != 0 ]; then 
				echo "Error, directory already exists, exiting..."
				exit
		fi;
		cd $directory
		mkdir build 
		svn co -r $revision https://gev.uchicago.edu/svn/damic-analysis/trunk/sim/DamicSimu/source		
		#Try to figure out our svn revision number
		version=`sh -c "cd source; svnversion"`
		echo "Using svn revision number: $version."
		cd build
		echo $version > SVN_VERSION
		cmake ../source/
		mkdir logs macros
		make -j 6
		cp $geant4_env ./
fi
