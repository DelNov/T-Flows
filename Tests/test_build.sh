#!/bin/bash

# Description: This script is made for debug purposes.
# It automatically builds and runs most cases in T-Flows

# Requires:
# mpi, gfortran, git to launch tests
# python-matplotlib + texlive-base to plot results

# Exit when any command fails
set -e

# Compilation flags used in makefiles
FCOMP="gnu"
# Conduct tests with DEBUG=yes
DEBUG="no"
# Conduct tests with CGNS or not
CGNS="no"
CGNS_MPI="openmpi"

# A small reminder how to set up alternatives if you have mpich and openmpi:
#update-alternatives --install /usr/bin/mpif90 mpif90 /usr/bin/mpif90.openmpi 20
#update-alternatives --install /usr/bin/mpif90 mpif90 /usr/bin/mpif90.mpich   80

#update-alternatives --install /usr/bin/mpirun mpirun /usr/bin/mpirun.openmpi 20
#update-alternatives --install /usr/bin/mpirun mpirun /usr/bin/mpirun.mpich   80

#-------------------------
# Folders with test cases
#-------------------------
LAMINAR_BACKSTEP_ORTH_DIR=Laminar/Backstep/Orthogonal
LAMINAR_BACKSTEP_NON_ORTH_DIR=Laminar/Backstep/Nonorthogonal
LAMINAR_CAVITY_LID_DRIVEN_DIR=Laminar/Cavity/Lid_Driven/Re_1000
LAMINAR_CAVITY_THERM_DRIVEN_106_DIR=Laminar/Cavity/Thermally_Driven/Ra_10e6
LAMINAR_CAVITY_THERM_DRIVEN_108_DIR=Laminar/Cavity/Thermally_Driven/Ra_10e8
LAMINAR_T_JUNCTION_DIR=Swarm/T_Junction_Square
LAMINAR_CHANNEL=Laminar/Accuracy_Test/Channel_Re_2000

MULTDOM_MEMBRANE_DIR=Rans/Membrane
SWARM_PERIODIC_CYL_DIR=Swarm/Cylinders_Periodic

# Add compressed meshes for these:
# MULTDOM_SINGLE_ROD_DIR=Rans/Single_Rod
# MULTDOM_COPY_INLET_DIR=Laminar/Copy_Inlet
# MULTDOM_HEAT_EXCHANGER_2_DIR=Laminar/Heat_Exchanger/2_Domains
# MULTDOM_HEAT_EXCHANGER_3_DIR=Laminar/Heat_Exchanger/3_Domains

LES_CHANNEL_180_PD_DIR=Les/Channel_Re_Tau_180/Periodic_Domain

# Generate takes ages for this one
LES_CHANNEL_180_LD_DIR=Les/Channel_Re_Tau_180/Long_Domain

LES_PIPE_DIR=Les/Pipe_Re_Tau_180
LES_RB_109_DIR=Les/Rayleigh_Benard_Convection_Ra_10e09

RANS_BACKSTEP_05100_DIR=Rans/Backstep_Re_05100
RANS_BACKSTEP_28000_DIR=Rans/Backstep_Re_28000
RANS_CHANNEL_LR_LONG_DIR=Rans/Channel_Re_Tau_590/Long_Domain
RANS_CHANNEL_LR_RSM_DIR=Rans/Channel_Re_Tau_590/Rsm
RANS_CHANNEL_LR_STRETCHED_DIR=Rans/Channel_Re_Tau_590/Stretched_Mesh
RANS_CHANNEL_LR_UNIFORM_DIR=Rans/Channel_Re_Tau_590/Uniform_Mesh
RANS_FUEL_BUNDLE_DIR=Rans/Fuel_Bundle
RANS_IMPINGING_JET_DIR=Rans/Impinging_Jet_2d_Distant_Re_23000

HYB_CHANNEL_HR_STRETCHED_DIR=Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh
HYB_CHANNEL_HR_UNIFORM_DIR=Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh

#----------------------------------------------------------------------------
# All compilation tests including those with User_Mod/
#----------------------------------------------------------------------------
ALL_COMPILE_TESTS=("$LAMINAR_CAVITY_LID_DRIVEN_DIR" \
                   "$LAMINAR_CAVITY_THERM_DRIVEN_106_DIR" \
                   "$LAMINAR_CAVITY_THERM_DRIVEN_108_DIR" \
                   "$LAMINAR_T_JUNCTION_DIR" \
                   "$LES_CHANNEL_180_LD_DIR" \
                   "$LES_CHANNEL_180_PD_DIR" \
                   "$LES_PIPE_DIR" \
                   "$LES_RB_109_DIR" \
                   "$RANS_BACKSTEP_28000_DIR" \
                   "$RANS_CHANNEL_LR_UNIFORM_DIR" \
                   "$RANS_CHANNEL_LR_STRETCHED_DIR" \
                   "$RANS_CHANNEL_LR_RSM_DIR" \
                   "$HYB_CHANNEL_HR_UNIFORM_DIR" \
                   "$HYB_CHANNEL_HR_STRETCHED_DIR" \
                   "$MULTDOM_MEMBRANE_DIR" \
                   "$SWARM_PERIODIC_CYL_DIR" \
                   )
DONE_COMPILE_TESTS=0

#----------------------------------
# All directories to test Generate
#----------------------------------
ALL_GENERATE_TESTS=("$LAMINAR_BACKSTEP_ORTH_DIR" \
                    "$LAMINAR_BACKSTEP_NON_ORTH_DIR" \
                    "$LAMINAR_CAVITY_LID_DRIVEN_DIR" \
                    "$LAMINAR_CAVITY_THERM_DRIVEN_106_DIR" \
                    "$LAMINAR_CAVITY_THERM_DRIVEN_108_DIR" \
                    "$LES_CHANNEL_180_PD_DIR" \
                    "$RANS_BACKSTEP_05100_DIR" \
                    "$RANS_BACKSTEP_28000_DIR" \
                    "$RANS_CHANNEL_LR_LONG_DIR" \
                    "$RANS_CHANNEL_LR_RSM_DIR" \
                    "$RANS_CHANNEL_LR_STRETCHED_DIR" \
                    "$RANS_CHANNEL_LR_UNIFORM_DIR" \
                    "$HYB_CHANNEL_HR_UNIFORM_DIR" \
                    "$HYB_CHANNEL_HR_STRETCHED_DIR" \
                    "$MULTDOM_COPY_INLET" \
                    )
DONE_GENERATE_TESTS=0

#---------------------------------
# All directories to test Convert
#---------------------------------
ALL_CONVERT_TESTS=("$RANS_IMPINGING_JET_DIR" \
                   "$RANS_FUEL_BUNDLE_DIR" \
                   "$LES_PIPE_DIR" \
                   "$MULTDOM_MEMBRANE_DIR" \
                   "$SWARM_PERIODIC_CYL_DIR" \
                   )
DONE_CONVERT_TESTS=0

#--------------------------------
# All directories to test Divide
#--------------------------------
ALL_DIVIDE_TESTS=("$LAMINAR_BACKSTEP_ORTH_DIR" \
                  "$LAMINAR_BACKSTEP_NON_ORTH_DIR" \
                  "$LAMINAR_CAVITY_LID_DRIVEN_DIR" \
                  "$LAMINAR_CAVITY_THERM_DRIVEN_106_DIR" \
                  "$LAMINAR_CAVITY_THERM_DRIVEN_108_DIR" \
                  "$LES_CHANNEL_180_PD_DIR" \
                  "$RANS_BACKSTEP_05100_DIR" \
                  "$RANS_BACKSTEP_28000_DIR" \
                  "$RANS_CHANNEL_LR_LONG_DIR" \
                  "$RANS_CHANNEL_LR_RSM_DIR" \
                  "$RANS_CHANNEL_LR_STRETCHED_DIR" \
                  "$RANS_CHANNEL_LR_UNIFORM_DIR" \
                  "$HYB_CHANNEL_HR_STRETCHED_DIR" \
                  "$HYB_CHANNEL_HR_UNIFORM_DIR" \
                  "$RANS_IMPINGING_JET_DIR" \
                  "$RANS_FUEL_BUNDLE_DIR"
                  "$LES_PIPE_DIR" \
                  "$MULTDOM_MEMBRANE_DIR" \
                  "$SWARM_PERIODIC_CYL_DIR" \
                  )
DONE_DIVIDE_TESTS=0

#-----------------------------------------------
# All directories to test save_now and exit_now
#-----------------------------------------------
ALL_SAVE_EXIT_NOW_TESTS=("$LAMINAR_BACKSTEP_ORTH_DIR")

#----------------------------------------------------------------
# All directories to test Process, followed by turbulence models
#----------------------------------------------------------------
ALL_PROCESS_TESTS=("$LAMINAR_CAVITY_LID_DRIVEN_DIR" \
                   "$LAMINAR_CAVITY_THERM_DRIVEN_106_DIR" \
                   "$LAMINAR_CAVITY_THERM_DRIVEN_108_DIR" \
                   "$RANS_CHANNEL_LR_UNIFORM_DIR" \
                   "$RANS_CHANNEL_LR_STRETCHED_DIR" \
                   "$RANS_CHANNEL_LR_RSM_DIR" \
                   "$RANS_CHANNEL_LR_RSM_DIR" \
                   "$HYB_CHANNEL_HR_STRETCHED_DIR" \
                   "$HYB_CHANNEL_HR_UNIFORM_DIR" \
                   "$LES_PIPE_DIR" \
                   "$SWARM_PERIODIC_CYL_DIR" \
                   )
ALL_PROCESS_MODELS=("none" \
                    "none" \
                    "none" \
                    "k_eps" \
                    "k_eps_zeta_f" \
                    "rsm_manceau_hanjalic" \
                    "rsm_hanjalic_jakirlic" \
                    "hybrid_les_rans" \
                    "hybrid_les_rans" \
                    "les_dynamic" \
                    "lagrangian_particle_tracking")
DONE_PROCESS_TESTS=0

#----------------------------------------------------------------------------
# All directories to Process accuracy test
#----------------------------------------------------------------------------
ALL_PROCESS_ACCURACY_TESTS=("$LAMINAR_CHANNEL")
DONE_PROCESS_ACCURACY_TESTS=0

# Folder structure
TEST_DIR=$PWD                      # Dir with tests
GENE_DIR=$PWD/../Sources/Generate  # Generate src folder
CONV_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIVI_DIR=$PWD/../Sources/Divide    # Divide   src folder
PROC_DIR=$PWD/../Sources/Process   # Process  src folder
BINA_DIR=$PWD/../Binaries          # Binaries folder

# Executables
GENE_EXE=$BINA_DIR/Generate        # Generate executable
CONV_EXE=$BINA_DIR/Convert         # Convert  executable
DIVI_EXE=$BINA_DIR/Divide          # Divide   executable
PROC_EXE=$BINA_DIR/Process         # Process  executable

# Start time measurements from this moment
current_time=$(date +%s)

# Script logs
FULL_LOG=$TEST_DIR/test_build.log # logs of current script
if [ -f $FULL_LOG ]; then cp /dev/null $FULL_LOG; fi
# echo "Full log is being written in file" "$FULL_LOG"

# Keep track of the last executed command
# trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
# trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#------------------------------------------------------------------------------#
# time in seconds
#------------------------------------------------------------------------------#
function time_in_seconds {
  previous_time=$current_time
  current_time=$(date +%s)
  echo "time elapsed:" \
  "$(echo ""$current_time" - "$previous_time"" | bc -l)" "seconds"
}

#------------------------------------------------------------------------------#
# clean_compile
#
# return success
#------------------------------------------------------------------------------#
function clean_compile {
  # $1 = dir
  # $2 = MPI = yes/no
  # $3 = DIR_CASE path

  if [ -z "${1+xxx}" ]; then 
    echo "directory with sources is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "MPI flag is not set at all"
    exit 1
  fi

  cd $1
  echo "clean compile in:" "$1"
  make clean >> $FULL_LOG 2>&1

  if [ -z "${3+xxx}" ]; then
    echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS=$CGNS CGNS_MPI=$CGNS_MPI MPI=$2"
    make \
      FORTRAN=$FCOMP \
      DEBUG=$DEBUG \
      CGNS=$CGNS \
      CGNS_MPI=$CGNS_MPI \
      MPI=$2 >> $FULL_LOG 2>&1
    success=$?
  else
    echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS=$CGNS CGNS_MPI=$CGNS_MPI MPI=$2 DIR_CASE=$3"
    make \
      FORTRAN=$FCOMP \
      DEBUG=$DEBUG \
      CGNS=$CGNS \
      CGNS_MPI=$CGNS_MPI \
      MPI=$2 \
      DIR_CASE=$3 >> $FULL_LOG 2>&1
    success=$?
  fi

  time_in_seconds

  cd - > /dev/null
  return $success
  if [ $success -eq 0 ]; then
    echo "Clean compile passed."
  else
    echo "Clean compile in " $1 " failed!"
  fi
}
#------------------------------------------------------------------------------#
# make links
#------------------------------------------------------------------------------#
function make_links {
  # $1 = dir
  cd $1
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
}

#------------------------------------------------------------------------------#
# launch Generate
#
# return success
#------------------------------------------------------------------------------#
function launch_generate {
  # $1 = relative dir
  # $2 = quiet

  make_links $TEST_DIR/$1
  if [ "$2" != "quiet" ]; then
    echo ""
    echo "#=================================================================="
    echo "#   Generate test:" $1
    echo "#------------------------------------------------------------------"
    echo "generate.scr: " >> $FULL_LOG 2>&1
  fi

  if [ $(ls -1 generate*.scr 2>/dev/null | wc -l) -gt 0 ]; then
    if [ -L gen.scr ]; then unlink gen.scr; fi
    for gen in $(ls -1 generate*.scr); do
      ln -s "$gen" gen.scr
      cat gen.scr >> $FULL_LOG 2>&1
      $GENE_EXE < gen.scr >> $FULL_LOG 2>&1
      unlink gen.scr
    done
  fi

  success=$?
  time_in_seconds
  return $success
}

#------------------------------------------------------------------------------#
# launch Divide
#
# return success
#------------------------------------------------------------------------------#
function launch_divide {
  # $1 = relative dir
  # $2 = quiet

  make_links $TEST_DIR/$1
  if [ "$2" != "quiet" ]; then
    echo ""
    echo "#=================================================================="
    echo "#   Divide test:" $1
    echo "#------------------------------------------------------------------"
    echo "divide.scr: " >> $FULL_LOG 2>&1
  fi

  if [ $(ls -1 divide*.scr 2>/dev/null | wc -l) -gt 0 ]; then
    if [ -L div.scr ]; then unlink div.scr; fi
    for div in $(ls -1 divide*.scr); do
      ln -s "$div" div.scr
      cat div.scr >> $FULL_LOG 2>&1
      $DIVI_EXE < div.scr >> $FULL_LOG 2>&1
      unlink div.scr
    done
  fi

  success=$?
  time_in_seconds
  return $success
}

#------------------------------------------------------------------------------#
# launch Convert
#
# returns success
#------------------------------------------------------------------------------#
function launch_convert {
  # $1 = relative dir
  echo ""
  echo "#=================================================================="
  echo "#   Convert test:" $1
  echo "#------------------------------------------------------------------"

  make_links $TEST_DIR/$1
  unpack_mesh $TEST_DIR/$1

  echo "convert.scr: " >> $FULL_LOG 2>&1

  if [ $(ls -1 convert*.scr 2>/dev/null | wc -l) -gt 0 ]; then
    if [ -L conv.scr ]; then unlink conv.scr; fi
    for conv in $(ls -1 convert*.scr); do
      ln -s "$conv" conv.scr
      cat conv.scr >> $FULL_LOG 2>&1
      $CONV_EXE < conv.scr >> $FULL_LOG 2>&1
      unlink conv.scr
    done
  fi

  success=$?
  time_in_seconds
  return $success
}

#------------------------------------------------------------------------------#
# launch Process
#
# return success
#------------------------------------------------------------------------------#
function launch_process {
  # $1 - seq/par
  # $2 - number of threads

  if [ -z "${1+xxx}" ]; then 
    echo "seq/par is not specified"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "threads count is not set"
    exit 1
  fi

  if [ "$1" == "seq" ]; then
    echo "launching Process"
    $PROC_EXE >> $FULL_LOG 2>&1
    success=$?
  elif [ "$1" == "par" ]; then
    echo "launching mpirun -np "$2" Process"
    mpirun -np $2 $PROC_EXE >> $FULL_LOG 2>&1
    success=$?
  else
    echo "wrong argument in launch_process"
    exit 1
  fi

  echo "control: " >> $FULL_LOG 2>&1
  cat control >> $FULL_LOG 2>&1
  time_in_seconds
  return $success
}

#------------------------------------------------------------------------------#
# Generate tests
#------------------------------------------------------------------------------#
function generate_tests {

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Generate tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  clean_compile $GENE_DIR no # dir MPI

  for CASE_DIR in ${ALL_GENERATE_TESTS[@]}; do
    launch_generate $CASE_DIR
  done

  DONE_GENERATE_TESTS=1
}

#------------------------------------------------------------------------------#
# unpack mesh (.tar & tar.xz, .tgz & tar.gz, .gz, 7z)
#------------------------------------------------------------------------------#
function unpack_mesh {
  # $1 = directory with test case

  if [ -z "${1+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  cd "$1"

  # .gz, *.tgz
  if [ $(ls -1 {*.gz,*.tgz} 2>/dev/null | wc -l) -gt 0 ]; then
    for archive in $(ls -1 {*.gz,*.tgz} 2>/dev/null); do
      echo "uncompressing archive:" "$archive"
      gunzip -kf "$archive" >> $FULL_LOG 2>&1
    done
  fi

  # .tar, .tar.xz
  if [ $(ls -1 {*.tar,*.tar.xz} 2>/dev/null | wc -l) -gt 0 ]; then
    for archive in $(ls -1 {*.tar,*.tar.xz} 2>/dev/null); do
      echo "uncompressing archive:" "$archive"
      tar -xf "$archive" >> $FULL_LOG 2>&1
    done
  fi

  # .7z
  if [ $(ls -1 *.7z 2>/dev/null | wc -l) -gt 0 ]; then
    for archive in $(ls -1 *.7z 2>/dev/null); do
      echo "uncompressing archive:" "$archive"
      7z x -y "$archive" >> $FULL_LOG 2>&1
    done
  fi

  cd - > /dev/null
}

#------------------------------------------------------------------------------#
# convert tests
#------------------------------------------------------------------------------#
function convert_tests {

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Convert tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  clean_compile $CONV_DIR no # dir MPI

  for CASE_DIR in ${ALL_CONVERT_TESTS[@]}; do
    launch_convert $CASE_DIR
  done

  DONE_CONVERT_TESTS=1
}

#------------------------------------------------------------------------------#
# Divide tests
#------------------------------------------------------------------------------#
function divide_tests {

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Divide tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  clean_compile $DIVI_DIR no # dir MPI

  for CASE_DIR in ${ALL_DIVIDE_TESTS[@]}; do
    launch_divide $CASE_DIR
  done

  DONE_DIVIDE_TESTS=1
}

#------------------------------------------------------------------------------#
# Replace line with first occurence of string to arg in file
#------------------------------------------------------------------------------#
function replace_line_with_first_occurence_in_file {
  #$1 string to search in file
  #$2 new line
  #$3 file

  new_line=$(echo "$2")
  line_to_replace="$(grep -ni "$1" $3 | cut -d: -f1 | head -n1)"
  if [ -z "$line_to_replace" ]; then
     # if no such string is found-> add it
    echo $new_line >> $3
    echo "Warning: added "$2" to the end of $3"
  else
    sed --follow-symlinks -i ""$line_to_replace"s%.*%$new_line%" $3
  fi

  # previous approach was with
  # awk -v var="$line" 'NR==var {$2=3}1'  control > .control.tmp
}

#------------------------------------------------------------------------------#
# Get value next to the keyword
#------------------------------------------------------------------------------#
function get_value_next_to_keyword {
  #$1 string to search in file
  #$2 file

  # This function searches for a first line with occurence of $1 in $2
  #   and gets a value right next to it regardless of 
  #   preceeding or following tabs & spaces

  # Thus, from a line
  #   \tab   \tab  MASS_DENSITY  \tab \tab     1.0 \tab   \tab
  # it returns 1.0

  if [ -z "${1+xxx}" ]; then 
    echo "Keyword to search is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "File to search in is not set at all"
    exit 1
  fi

  echo "$(grep -i $1 $2 | tr -s '\t' | sed 's/^[ \t]*//' | \
    sed 's/\t/ /' | tr -s ' ' | cut -d' ' -f2 | head -n1)"
}

#------------------------------------------------------------------------------#
# processor: backup test
#------------------------------------------------------------------------------#
function process_backup_test {
  # 1 = test_dir

  if [ -z "${1+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  echo "Process backup tests.."
  cd "$1"

  n1=$(printf "%06d" 1)
  n2=$(printf "%06d" 2)

  # Unlink all possible links, which are used below
  if [ -L control ]; then unlink control; fi
  for i in {0..9}; do
    if [ -L control."$i" ]; then unlink control."$i"; fi
  done
  if [ -L divide.1.scr ]; then unlink divide.1.scr; fi

  n_dom=$(ls -1 control.? 2>/dev/null | wc -l)
  if [ $n_dom -eq 0 ]; then
    n_dom=1
    ln -s control control.0
    ln -s control control.1
    ln -s divide.scr divide.1.scr
  elif [ $n_dom -gt 1 ]; then
    n_dom=$(($n_dom-1))
    ln -s control.0 control
  fi

  nproc_in_div=$(head -n2 divide.1.scr | tail -n1)

  # BEGIN:---------------------------------------#
  echo "np=1, MPI=no, start from 0, make a backup"
  clean_compile $PROC_DIR no # dir MPI

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    # comment line with LOAD_BACKUP_NAME
    replace_line_with_first_occurence_in_file \
      "LOAD_BACKUP_NAME" \
      "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
      control."$i"

    # change number of timesteps to 3
    replace_line_with_first_occurence_in_file \
      "NUMBER_OF_TIME_STEPS" \
      "NUMBER_OF_TIME_STEPS 3" \
      control.0

    # change backup interval to 1 ts
    replace_line_with_first_occurence_in_file \
      "BACKUP_SAVE_INTERVAL" \
      "BACKUP_SAVE_INTERVAL 1" \
      control.0
  done

  launch_process seq 1
  #-----------------------------------------:END #


  # BEGIN:---------------------------------------------#
  echo "np=1, MPI=no, load from backup(produced by seq)"

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    # uncomment line with LOAD_BACKUP_NAME
    replace_line_with_first_occurence_in_file \
      "LOAD_BACKUP_NAME" \
      "LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
      control."$i"
  done

  launch_process seq 1
  #-----------------------------------------------:END #


  # BEGIN:----------------------------------------------#
  echo "np=1, MPI=yes, load from backup(produced by seq)"
  clean_compile $PROC_DIR yes # dir MPI

  launch_process par 1
  #------------------------------------------------:END #


  # BEGIN:----------------------------------------------#
  echo "np=2, MPI=yes, load from backup(produced by seq)"

  launch_process par $nproc_in_div
  #------------------------------------------------:END #


  # BEGIN:---------------------------------------------------#
  echo "np=2, MPI=yes, load from backup(produced by par.np=1)"

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    replace_line_with_first_occurence_in_file \
      "LOAD_BACKUP_NAME" \
      "LOAD_BACKUP_NAME "$name_in_div"-ts"$n2".backup" \
      control."$i"
  done

  launch_process par $nproc_in_div
  #-----------------------------------------------------:END #


  # BEGIN:----------------------------------------#
  echo "np=2, MPI=yes, start from 0, make a backup"

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    replace_line_with_first_occurence_in_file \
      "LOAD_BACKUP_NAME" \
      "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
      control."$i"
  done

  launch_process par $nproc_in_div
  #------------------------------------------:END #


  # BEGIN:----------------------------------------------------#
  echo "np=2, MPI=yes, load from backup(produced by par.np=2)"

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    replace_line_with_first_occurence_in_file \
      "LOAD_BACKUP_NAME" \
      "LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
      control."$i"
  done

  launch_process par $nproc_in_div
  #------------------------------------------------------:END #


  # BEGIN:------------------------------------------#
  echo "np=1, MPI=yes, backup=(produced by par.np=2)"
  clean_compile $PROC_DIR no # dir MPI
  launch_process par 1
  #--------------------------------------------:END #


  # Restore control
  if [ $n_dom -eq 1 ]; then
    unlink control.0
    unlink control.1
    unlink divide.1.scr

    git checkout control
  else
    unlink control
    git checkout control.?
  fi

  cd - > /dev/null

}

#------------------------------------------------------------------------------#
# process backup tests
#------------------------------------------------------------------------------#
function process_backup_tests {

  # Grasp/embrace as many different model combinations as you can
  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor backup tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  echo ""
  echo "#======================================================================"
  echo "#   Test 1: "$RANS_CHANNEL_LR_UNIFORM_DIR" [k_eps model + T]"
  echo "#----------------------------------------------------------------------"
  replace_line_with_first_occurence_in_file \
    "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL k_eps" \
    $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR/control

  process_backup_test \
    $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR

  echo ""
  echo "#======================================================================"
  echo "#   Test 2: "$RANS_CHANNEL_LR_UNIFORM_DIR" [k_eps_zeta_f model + T]"
  echo "#----------------------------------------------------------------------"
  replace_line_with_first_occurence_in_file \
    "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL k_eps_zeta_f" \
    $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR/control

  process_backup_test \
    $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR

  echo ""
  echo "#======================================================================"
  echo "#   Test 3: "$RANS_CHANNEL_LR_RSM_DIR" [rsm_hanjalic_jakirlic model + T]"
  echo "#----------------------------------------------------------------------"
  replace_line_with_first_occurence_in_file \
    "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL rsm_hanjalic_jakirlic" \
    $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR/control

  process_backup_test \
    $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR

  echo ""
  echo "#======================================================================"
  echo "#   Test 4: "$RANS_CHANNEL_LR_RSM_DIR" [rsm_manceau_hanjalic model + T]"
  echo "#----------------------------------------------------------------------"
  replace_line_with_first_occurence_in_file \
    "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL rsm_manceau_hanjalic" \
    $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR/control

  process_backup_test \
    $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR

  echo ""
  echo "#======================================================================"
  echo "#   Test 5: "$LES_PIPE_DIR" les_dynamic"
  echo "#----------------------------------------------------------------------"
  process_backup_test \
    $TEST_DIR/$LES_PIPE_DIR

  echo ""
  echo "#======================================================================"
  echo "#   Test 6: "$LAMINAR_CAVITY_LID_DRIVEN_DIR" none"
  echo "#----------------------------------------------------------------------"
  process_backup_test \
    $TEST_DIR/$LAMINAR_CAVITY_LID_DRIVEN_DIR

  echo ""
  echo "#======================================================================"
  echo "#   Test 7: "$MULTDOM_MEMBRANE_DIR" multidom: k_eps_zeta_f model + T"
  echo "#----------------------------------------------------------------------"
  process_backup_test \
    $TEST_DIR/$MULTDOM_MEMBRANE_DIR

}
#------------------------------------------------------------------------------#
# process: save_now / exit_now test
#------------------------------------------------------------------------------#
function process_save_exit_now_test {
  # $1 = relative dir with test

  if [ -z "${1+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test save_now & exit_now on:" $1
  echo "#----------------------------------------------------------------------"

  cd "$TEST_DIR/$1"

  # get rid of "save_now" and "exit_now" files if they happen to exist
  if [ -f "save_now" ]; then rm save_now; fi
  if [ -f "exit_now" ]; then rm exit_now; fi

  n1=$(printf "%06d" 1)

  # Unlink all possible links, which are used below
  if [ -L control ]; then unlink control; fi
  for i in {0..9}; do
    if [ -L control."$i" ]; then unlink control."$i"; fi
  done
  if [ -L divide.1.scr ]; then unlink divide.1.scr; fi

  n_dom=$(ls -1 control.? 2>/dev/null | wc -l)
  if [ $n_dom -eq 0 ]; then
    n_dom=1
    ln -s control control.0
    ln -s control control.1
    ln -s divide.scr divide.1.scr
  elif [ $n_dom -gt 1 ]; then
    n_dom=$(($n_dom-1))
    ln -s control.0 control
  fi

  nproc_in_div=$(head -n2 divide.1.scr | tail -n1)

  # BEGIN:---------------------------------------#
  echo "np=1, MPI=no, start from 0, make a backup"
  clean_compile $PROC_DIR no # dir MPI

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    # change number of timesteps to 3
    replace_line_with_first_occurence_in_file \
      "NUMBER_OF_TIME_STEPS" \
      "NUMBER_OF_TIME_STEPS 3" \
      control.0

    # change backup interval to 1 ts
    replace_line_with_first_occurence_in_file \
      "BACKUP_SAVE_INTERVAL" \
      "BACKUP_SAVE_INTERVAL 10" \
      control.0
  done

  for i in {1..3}; do
    echo ""
    echo "#===================================================================="
    if [ "$i" = 1 ]; then
      echo "#   Test np=1, MPI=no"
    fi
    if [ "$i" = 2 ]; then
      echo "#   Test np=1, MPI=yes"
    fi
    if [ "$i" = 3 ]; then
      echo "#   Test np=2, MPI=yes"
    fi
    echo "#--------------------------------------------------------------------"

    # BEGIN:---------------------------------------------#
    echo "np=1, MPI=no, load from backup(produced by seq)"

    for (( i=1; i<=$n_dom; i++ )); do
      name_in_div=$(head -n1 divide."$i".scr)

      # comment line with LOAD_BACKUP_NAME
      replace_line_with_first_occurence_in_file \
        "LOAD_BACKUP_NAME" \
        "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
        control."$i"
    done

    launch_process seq 1
    #-----------------------------------------------:END #

    if [ "$i" = 1 ]; then
      clean_compile $PROC_DIR no
    elif [ "$i" = 2 ]; then
      clean_compile $PROC_DIR yes
    elif [ "$i" = 3 ]; then
      clean_compile $PROC_DIR yes
    fi

    echo "#   Forcing to save: save_now"
    touch save_now

    # get current line count where search starts
    n_start="$(echo "$(wc -l $FULL_LOG | cut -d" " -f1) + 1" | bc -l)"

    # start from scratch
    if [ "$i" = 1 ]; then
      launch_process seq 1
    elif [ "$i" = 2 ]; then
      launch_process par 1
    elif [ "$i" = 3 ]; then
      launch_process par $nproc_in_div
    fi

    # find if save was made in the range [n_start: n_finish]
    if tail -n+$n_start $FULL_LOG | \
      grep -q "# Creating the file: "$name_in_div"-ts"$n1"\|\
               # Creating the file with fields: "$name_in_div"-ts"$n1""; then

      echo "save_now was successfull"

      echo "Forcing to exit: exit_now"
      touch exit_now

      for (( i=1; i<=$n_dom; i++ )); do
        name_in_div=$(head -n1 divide."$i".scr)

        # uncomment line with LOAD_BACKUP_NAME
        replace_line_with_first_occurence_in_file \
          "LOAD_BACKUP_NAME" \
          "LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
          control."$i"
      done

      # start from ts=1
      n_start="$(echo "$(wc -l $FULL_LOG | cut -d" " -f1) + 1" | bc -l)"

      if [ "$i" = 1 ]; then
        launch_process seq 1
      elif [ "$i" = 2 ]; then
        launch_process par 1
      elif [ "$i" = 3 ]; then
        launch_process par $nproc_in_div
      fi

      if tail -n+$n_start $FULL_LOG | \
        tr -s " " | \
        grep -q "Time step : 3"; then

          echo "exit_now was NOT successfull"
      else
        echo "exit_now was successfull"
      fi
    else
      echo "save_now was NOT successfull"
    fi

  done

  # Restore control
  if [ $n_dom -eq 1 ]; then
    unlink control.0
    unlink control.1
    unlink divide.1.scr

    git checkout control
  else
    unlink control
    git checkout control.?
  fi
}

#------------------------------------------------------------------------------#
# process save_now / exit_now tests
#------------------------------------------------------------------------------#
function process_save_exit_now_tests {

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor save_now and exit_now tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  for CASE_DIR in ${ALL_SAVE_EXIT_NOW_TESTS[@]}; do
    process_save_exit_now_test \
      $CASE_DIR
  done
}

#------------------------------------------------------------------------------#
# Launch matplotlib script which plots and prints results in .png format
#------------------------------------------------------------------------------#
function launch_matplotlib {
  # $1 = dir with matplotlib script
  # $2 = matplotlib script name
  # $3 = input file  [or "" if none]
  # $4 = output file [or "" if none]

  if [ -z "${1+xxx}" ]; then 
    echo "directory with matplotlib script is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "matplotlib script name is not set at all"
    exit 1
  elif [ -z "${3+xxx}" ]; then 
    echo "input file name is not set at all"
    exit 1
  elif [ -z "${4+xxx}" ]; then 
    echo "output file name is not set at all"
    exit 1
  fi

  cd "$1"
  if [ ! -f "$2" ]; then
    echo "Warning: $1/$2 does not exist"
    return
  fi
  sed "s%DAT_FILE_WITH_RESULTS_MACRO%$3%" "$2" > ./tmp
  sed -i "s%PNG_FILE_WITH_RESULTS_MACRO%$4%"     ./tmp
  # Launch script:
  $(sed -e '/#/d' ./tmp) >> $FULL_LOG 2>&1
  if [ "$4" == "" ]; then
    echo "new figure was created:" "$1"/results.png
  else
    echo "new figure was created:" "$1"/"$4".png
  fi
  time_in_seconds
}
#------------------------------------------------------------------------------#
# Individual process tests for compilation
#------------------------------------------------------------------------------#
function process_compilation_test {
  # $1 = relative dir with test

  echo ""
  echo "#=================================================================="
  echo "#  Compilation test in:" $1
  echo "#------------------------------------------------------------------"
  if [ -z "$TEST_DIR/${1+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  cd "$TEST_DIR/$1"

  # Unlink all possible links, which are used below
  if [ -L control ]; then unlink control; fi
  for i in {0..9}; do
    if [ -L control."$i" ]; then unlink control."$i"; fi
  done
  if [ -L divide.1.scr ]; then unlink divide.1.scr; fi

  n_dom=$(ls -1 control.? 2>/dev/null | wc -l)
  if [ $n_dom -eq 0 ]; then
    n_dom=1
    ln -s control control.1
    ln -s divide.scr divide.1.scr
  elif [ $n_dom -gt 1 ]; then
    n_dom=$(($n_dom-1))
    ln -s control.0 control
  fi

  nproc_in_div=$(head -n2 divide.1.scr | tail -n1)
  name_in_div=$(head -n1 divide.1.scr)

  # rel_dir to User_Mod/ from Process/
  rel_dir=$(realpath --relative-to="$PROC_DIR" "$TEST_DIR/$1")

  clean_compile $PROC_DIR yes $rel_dir # dir MPI DIR_CASE

  # Restore control
  if [ $n_dom -eq 1 ]; then
    unlink control.1
    unlink divide.1.scr

    git checkout control
  else
    unlink control
    git checkout control.?
  fi
}

#------------------------------------------------------------------------------#
# All process compilation tests
#------------------------------------------------------------------------------#
function process_compilation_tests {
  # $1 = dir with test

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor compilation tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  for CASE_DIR in ${ALL_COMPILE_TESTS[@]}; do
    process_compilation_test $CASE_DIR
  done
}

#------------------------------------------------------------------------------#
# Individual process tests for full length
#------------------------------------------------------------------------------#
function process_full_length_test {
  # $1 = dir with test
  # $2 = model
  # $3 = dir with results

  if [ -z "${1+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "model is not set at all"
    exit 1
  elif [ -z "${3+xxx}" ]; then 
    echo "directory with results is not set at all"
    exit 1
  fi

  cd "$1"

  n1=$(printf "%06d" 1)

  # Unlink all possible links, which are used below
  if [ -L control ]; then unlink control; fi
  for i in {0..9}; do
    if [ -L control."$i" ]; then unlink control."$i"; fi
  done
  if [ -L divide.1.scr ]; then unlink divide.1.scr; fi

  n_dom=$(ls -1 control.? 2>/dev/null | wc -l)
  if [ $n_dom -eq 0 ]; then
    n_dom=1
    ln -s control control.1
    ln -s divide.scr divide.1.scr
  elif [ $n_dom -gt 1 ]; then
    n_dom=$(($n_dom-1))
    ln -s control.0 control
  fi

  nproc_in_div=$(head -n2 divide.1.scr | tail -n1)

  # BEGIN:-------------------------#
  echo "np="$nproc_in_div", MPI=yes"
  # rel_dir to User_Mod/ from Process/
  rel_dir=$(realpath --relative-to="$PROC_DIR" "$1")

  clean_compile $PROC_DIR yes $rel_dir # dir MPI DIR_CASE

  for (( i=1; i<=$n_dom; i++ )); do
    name_in_div=$(head -n1 divide."$i".scr)

    # comment line with LOAD_BACKUP_NAME
    replace_line_with_first_occurence_in_file \
      "LOAD_BACKUP_NAME" \
      "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" \
      control."$i"

    # change model to $2
    replace_line_with_first_occurence_in_file \
      "TURBULENCE_MODEL" \
      "TURBULENCE_MODEL "$2"" \
      control."$i"
  done

  launch_process par $nproc_in_div
  #---------------------------:END #

  # Restore control
  if [ $n_dom -eq 1 ]; then
    unlink control.1
    unlink divide.1.scr

    git checkout control
  else
    unlink control
    git checkout control.?
  fi

  # If results are present [produced by User_Mod_Save_Results]
  if [ $(ls -1 "$name_in_div"-res-plus-ts??????.dat 2>/dev/null | wc -l)\
    -gt 0 ]; then

    # extract essential data from produced .dat files
    last_results_plus_dat_file=$(realpath --relative-to="$3" \
      $(ls -tr1 "$name_in_div"-res-plus-ts??????.dat | tail -n1))

    echo "results are:"
    echo "$(head -n8 $(ls -tr1 "$name_in_div"-res-plus-ts??????.dat | \
      tail -n1))"

    launch_matplotlib \
      "$3" \
      readme_python_matplotlib_script \
      "$last_results_plus_dat_file" \
      "result_plus_"$2""
  else
      echo "Warning: file "$name_in_div"-res-plus-ts??????.dat does not exist"
  fi
}
#------------------------------------------------------------------------------#
# All process tests
#------------------------------------------------------------------------------#
function process_full_length_tests {
  # $1 = dir with test
  # $2 = model
  # $3 = dir with results
  # it requires a template file readme_python_matplotlib_script.sh in Results/

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor full simulation tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  for i in ${!ALL_PROCESS_TESTS[@]}; do
    CASE_DIR="${ALL_PROCESS_TESTS[$i]}"
    CASE_MOD="${ALL_PROCESS_MODELS[$i]}"
    #MASS_DENSITY=$(get_value_next_to_keyword "MASS_DENSITY" \
    #  "${ALL_PROCESS_TESTS[$i]}""/control")
    #DYNAMIC_VISCOSITY=$(get_value_next_to_keyword "DYNAMIC_VISCOSITY" \
    #  "${ALL_PROCESS_TESTS[$i]}""/control")

    echo ""
    echo "#===================================================================="
    echo "#   Process test: "     $CASE_DIR
    echo "#   with physical model: " $CASE_MOD
    echo "#--------------------------------------------------------------------"

    process_full_length_test \
      "$TEST_DIR/$CASE_DIR" \
      "$CASE_MOD" \
      "$TEST_DIR/$CASE_DIR/Results"
  done
}

#------------------------------------------------------------------------------#
# Individual process accuracy test
#------------------------------------------------------------------------------#
function process_accuracy_test {
  # $1 = case dir

  # Full path to test case
  path="$TEST_DIR"/"$1"

  if [ -z "${1+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  for i in {3..7..1}
  do
    Ny="$(echo "1 + 2^"$i"" | bc)"
    echo "Ny: "$Ny""

    cd $path

    # change Ny in chan.dom
    replace_line_with_first_occurence_in_file \
      "1  65  33" \
      "  1  65  33  "$Ny" # Nx Nz Ny" \
      chan.dom

    if [ ! -f $GENE_EXE ]; then
      clean_compile $GENE_DIR no # dir MPI
    fi
    launch_generate "$1" "quiet"

    if [ ! -f $DIVI_EXE ]; then
      clean_compile $DIVI_DIR no # dir MPI
    fi
    launch_divide   "$1" "quiet"

    name_in_div=$(head  -n1 divide.scr)
    nproc_in_div=$(head -n2 divide.scr | tail -n1)

    echo "np="$nproc_in_div", MPI=yes"

    # rel_dir to User_Mod/ from Process/
    rel_dir=$(realpath --relative-to="$PROC_DIR" "$path")

    clean_compile $PROC_DIR yes $rel_dir # dir MPI DIR_CASE

    launch_process par $nproc_in_div

    # Restore control
    git checkout control

    # Process chan-ts??????-res.dat
    if ls "$name_in_div"-res-ts??????.dat 1> /dev/null 2>&1; then

      # extract essential data from produced .dat files
      last_results_dat_file=\
        $(ls -tr1 "$name_in_div"-res-ts??????.dat | tail -n1)

      # Store this file with result
      nNy=$(printf "%06d" $((Ny-1)))
      cp "$last_results_dat_file" Results/"$nNy".dat

    else
        echo "Warning: file "$name_in_div"-res-ts??????.dat does not exist"
    fi
  done # for loop

  launch_matplotlib \
    "$path"/Results \
    readme_python_matplotlib_script "" ""
}

#------------------------------------------------------------------------------#
# Process accuracy tests (simple tests with given analytical solution)
#------------------------------------------------------------------------------#
function process_accuracy_tests {
  # $1 = test dir
  # it requires a template file readme_python_matplotlib_script.sh in Results/

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor accuracy test check"
  echo "#"
  echo "#----------------------------------------------------------------------"

  for i in ${!ALL_PROCESS_ACCURACY_TESTS[@]}; do
    CASE_DIR="${ALL_PROCESS_ACCURACY_TESTS[$i]}"

    echo ""
    echo "#===================================================================="
    echo "#   Process test: "     $CASE_DIR
    echo "#   with tolerance for u, p, SIMPLE : 1.e-8"
    echo "#--------------------------------------------------------------------"

    process_accuracy_test "$CASE_DIR"
  done

}

#------------------------------------------------------------------------------#
# chose test
#------------------------------------------------------------------------------#
function chose_test {

  # First argument is the option
  option=$1

  if [ $option -eq 0 ]; then exit 1; fi
  if [ $option -eq 1 ]; then
    generate_tests
  fi
  if [ $option -eq 2 ]; then
    convert_tests
  fi
  if [ $option -eq 3 ]; then
    if [ $DONE_GENERATE_TESTS -eq 0 ]; then generate_tests; fi
    if [ $DONE_CONVERT_TESTS  -eq 0 ]; then convert_tests;  fi
    divide_tests
  fi
  if [ $option -eq 4 ]; then process_compilation_tests;    fi
  if [ $option -eq 5 ]; then
    if [ $DONE_GENERATE_TESTS -eq 0 ]; then generate_tests; fi
    if [ $DONE_CONVERT_TESTS  -eq 0 ]; then convert_tests;  fi
    if [ $DONE_DIVIDE_TESTS   -eq 0 ]; then divide_tests;   fi
    process_backup_tests
  fi
  if [ $option -eq 6 ]; then
    if [ $DONE_GENERATE_TESTS -eq 0 ]; then generate_tests; fi
    if [ $DONE_CONVERT_TESTS  -eq 0 ]; then convert_tests;  fi
    if [ $DONE_DIVIDE_TESTS   -eq 0 ]; then divide_tests;   fi
    process_save_exit_now_tests
  fi
  if [ $option -eq 7 ]; then
    if [ $DONE_GENERATE_TESTS -eq 0 ]; then generate_tests; fi
    if [ $DONE_CONVERT_TESTS  -eq 0 ]; then convert_tests;  fi
    if [ $DONE_DIVIDE_TESTS   -eq 0 ]; then divide_tests;   fi
    process_full_length_tests
  fi
  if [ $option -eq 8 ]; then
    process_accuracy_tests
  fi
  if [ $option -eq 9 ]; then
    generate_tests
    convert_tests
    divide_tests
    process_compilation_tests
    process_backup_tests
    process_save_exit_now_tests
    process_full_length_tests
    process_accuracy_tests
  fi
  if [ $option -eq 10 ]; then
    git clean -dfx $TEST_DIR/..
  fi

}

#------------------------------------------------------------------------------#
# Actual script
#------------------------------------------------------------------------------#

# Desired option was provided by command line, execute it and end
if [ $1 ]; then
  echo ""
  echo "#===================================================================="
  echo "#"
  echo "#   T-Flows testing in non-interactive mode"
  echo "#"
  echo "#--------------------------------------------------------------------"
  echo ""
  chose_test $1

# Desired option was not provided, enter interactive mode
else
  while [ 0 -eq 0 ]; do
    echo ""
    echo "#===================================================================="
    echo "#"
    echo "#   T-Flows testing"
    echo "#"
    echo "#--------------------------------------------------------------------"
    echo ""
    echo "  Chose the type of test you want to perform:"
    echo ""
    echo "  0. Exit"
    echo "  1. Generate tests"
    echo "  2. Convert tests"
    echo "  3. Divide tests"
    echo "  4. Processor compilation tests with User_Mod"
    echo "  5. Processor backup tests"
    echo "  6. Processor save_now/exit_now tests"
    echo "  7. Processor full lenght tests"
    echo "  8. Process accuracy test"
    echo "  9. Perform all tests"
    echo " 10. Clean all test directories"
    echo ""
    read -p "  Enter the desired type of test: " option
    chose_test $option
  done
fi;
