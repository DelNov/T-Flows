#!/bin/bash

# Specify the number of processors here (for parallel mode): 
NUMBEROFCPUS=4

# Input from user on the phase to be run 
if [ $1 = "fluid" ]; then 

	# Cleaning any old files 
	rm -f out* *.faces *.monit readme
	rm -f *.cfn *.dim *.pvtu *.vtu
	rm -rf Sub*
	rm -f Process
	rm -f Generate
	rm -f Divide
	rm -f Convert
	
	# Remove any old output file for fluid phase (can be refined) 
        rm -f out_fluidResults_parallelMode.dat

	# Creating soft links from Binaries/ 
	ln -i -s ../../../Binaries/* . 
	
	# Generating mesh
	./Generate < generate.scr

	# Creating a soft link for the fluid flow control file 
	ln -i -s control_a_flow_development control 

	if [ $2 = "serial" ]
	then
	    echo "#-------------------------------------------------#"
	    echo "# Running T-Flows (fluid phase) in serial mode... #"
	    echo "#-------------------------------------------------#"
	    ./Process >| out_fluidResults_serialMode.dat &
	    tail -100f out_fluidResults_serialMode.dat
	elif [ $2 = "parallel" ]
	then
	    # Decomposing the domain
	    ./Divide < divide.scr
	    echo "#---------------------------------------------------#"
	    echo "# Running T-Flows (fluid phase) in parallel mode... #"
	    echo "#---------------------------------------------------#"
	    mpirun -np $NUMBEROFCPUS ./Process >| out_fluidResults_parallelMode.dat & 
	    tail -100f out_fluidResults_parallelMode.dat
	else 
	    echo "ERROR! Please specify one of the available run modes (i.e., 'serial' or 'parallel')"
	    exit 1
	fi

elif [ $1 = "particles" ]; then 

	# Saving backup files - just to be on the safe side (need to be done in a more elegant way!)
	mkdir BackupFiles
	cp *.backup BackupFiles 

	# Removing the link to the fluid control file 
	rm -f control

	# Remove any old output file for particles (can be refined) 
        rm -f out_particleResults_parallelMode.dat

	# Creating a soft link for the particle control file 
	ln -i -s control_b_particles control 
	
 	# Assign a backup file to retrieve the solution from 
	sed -i 's/#LOAD_BACKUP_NAME/LOAD_BACKUP_NAME    chan-ts000400.backup #/' control 

	if [ $2 = "serial" ]
	then
	    echo "#-----------------------------------------------------#"
	    echo "# Running T-Flows (dispersed phase) in serial mode... #"
	    echo "#-----------------------------------------------------#"
	    ./Process >| out_particleResults_serialMode.dat &
	    tail -100f out_particleResults_serialMode.dat
	elif [ $2 = "parallel" ]
	then
	    # Decomposing the domain !Careful: here the domain should be divided already from the fluid phase
	    #./Divide < divide.scr
	    echo "#-------------------------------------------------------#"
	    echo "# Running T-Flows (dispersed phase) in parallel mode... #"
	    echo "#-------------------------------------------------------#"
	    mpirun -np $NUMBEROFCPUS ./Process >| out_particleResults_parallelMode.dat & 
	    tail -100f out_particleResults_parallelMode.dat
	else 
	    echo "ERROR! Please specify one of the available run modes (i.e., 'serial' or 'parallel')"
	    exit 1
	fi
fi
