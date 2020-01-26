#!/bin/bash

# Description: This script is made for debug purposes.
# It automatically builds and runs most cases in T-Flows

# Requires: gnuplot, texlive-base, mpi, gfortran

# exit when any command fails
set -e

# Compilation flags used in makefiles
FCOMP="gnu"
# Conduct tests with DEBUG=yes
DEBUG="yes"
# Repeat tests with CGNS=yes
CGNS="yes"

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
LAMINAR_T_JUNCTION_DIR=Laminar/T_Junction

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
                   "$HYB_CHANNEL_HR_STRETCHED_DIR")
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
                    "$HYB_CHANNEL_HR_STRETCHED_DIR")
DONE_GENERATE_TESTS=0

#---------------------------------
# All directories to test Convert
#---------------------------------
ALL_CONVERT_TESTS=("$RANS_IMPINGING_JET_DIR" \
                   "$RANS_FUEL_BUNDLE_DIR" \
                   "$LES_PIPE_DIR")
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
                  "$LES_PIPE_DIR")
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
                   "$HYB_CHANNEL_HR_UNIFORM_DIR")
ALL_PROCESS_MODELS=("none" \
                    "none" \
                    "none" \
                    "k_eps" \
                    "k_eps_zeta_f" \
                    "rsm_manceau_hanjalic" \
                    "rsm_hanjalic_jakirlic" \
                    "hybrid_les_rans" \
                    "hybrid_les_rans")
DONE_PROCESS_TESTS=0

# Folder structure
TEST_DIR=$PWD                      # dir with tests
GENE_DIR=$PWD/../Sources/Generate  # Generate src folder
CONV_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIVI_DIR=$PWD/../Sources/Divide    # Divide   src folder
PROC_DIR=$PWD/../Sources/Process   # Process  src folder
BINA_DIR=$PWD/../Binaries/         # binaries folder

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
echo "Full log is being written in file" "$FULL_LOG"

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
  # $2 = CGNS = yes/no
  # $3 = MPI = yes/no
  # $4 = DIR_CASE path

  if [ -z "${1+xxx}" ]; then 
    echo "directory with sources is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "CGNS flag is not set at all"
    exit 1
  elif [ -z "${3+xxx}" ]; then 
    echo "MPI flag is not set at all"
    exit 1
  fi

  cd $1
  echo "clean compile in:" "$1"
  make clean >> $FULL_LOG 2>&1

  if [ -z "${4+xxx}" ]; then 
    echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS=$2 MPI=$3"
              make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS=$2 MPI=$3 \
              >> $FULL_LOG 2>&1
              success=$?
  else
    echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS=$2 MPI=$3 DIR_CASE=$4"
              make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS=$2 MPI=$3 DIR_CASE=$4 \
              >> $FULL_LOG 2>&1
              success=$?
  fi

  time_in_seconds

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
  make_links $TEST_DIR/$1
  echo ""
  echo "#=================================================================="
  echo "#   Generate test:" $1
  echo "#------------------------------------------------------------------"
  echo "generate.scr: " >> $FULL_LOG 2>&1
  cat generate.scr >> $FULL_LOG 2>&1
  $GENE_EXE < generate.scr >> $FULL_LOG 2>&1
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
  make_links $TEST_DIR/$1
  echo ""
  echo "#=================================================================="
  echo "#   Divide test:" $1
  echo "#------------------------------------------------------------------"
  echo "divide.scr: " >> $FULL_LOG 2>&1
  cat divide.scr >> $FULL_LOG 2>&1
  $DIVI_EXE < divide.scr >> $FULL_LOG 2>&1
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
  make_links $TEST_DIR/$1
  echo ""
  echo "#=================================================================="
  echo "#   Convert test:" $1
  echo "#------------------------------------------------------------------"
  echo "convert.scr: " >> $FULL_LOG 2>&1
  cat convert.scr >> $FULL_LOG 2>&1
  $CONV_EXE < convert.scr >> $FULL_LOG 2>&1
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

  #-- seq, no cgns
  clean_compile $GENE_DIR  no  no  # dir CGNS MPI

  for CASE_DIR in ${ALL_GENERATE_TESTS[@]}; do
    launch_generate $CASE_DIR
  done

  #-- seq, cgns(hdf5)
  if [ "$CGNS" = "yes" ]; then
    clean_compile $GENE_DIR  yes  no  # dir CGNS MPI

    for CASE_DIR in ${ALL_GENERATE_TESTS[@]}; do
      launch_generate $CASE_DIR
    done
  fi

  DONE_GENERATE_TESTS=1
}

#------------------------------------------------------------------------------#
# unpack .neu mesh
#------------------------------------------------------------------------------#
function unpack_neu_mesh {
  # $1 = archive name

   if [ -z "${1+xxx}" ]; then 
    echo "archive name is not set at all"
    exit 1
  fi

  echo "unpacking geometry" "$1"

  filename=$(basename -- "$1") # filename without path
  fname="${filename%%.*}"      # filename without extension
  ext="${filename#*.}"         # extension

  if [ -f "$fname"."neu" ]; then
    return
  elif [ -f "$fname"."tar.gz" ]; then
    tar -zxvf "$fname"."tar.gz"
  elif [ -f "$fname"."tgz" ]; then
    tar -zxvf "$fname"."tgz"
  elif [ -f "$fname"."neu.tgz" ]; then
    tar -zxvf "$fname"."neu.tgz"
  elif [ -f "$fname"."gz" ]; then
    gunzip -dv "$fname"."gz"
  elif [ -f "$fname"."neu.gz" ]; then
    gunzip -dv "$fname"."neu.gz"
  else
    echo "could not extract" "$fname"."$ext"
    exit 1
  fi
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

  cd $TEST_DIR/$RANS_IMPINGING_JET_DIR;  unpack_neu_mesh  jet.neu
  cd $TEST_DIR/$RANS_FUEL_BUNDLE_DIR;    unpack_neu_mesh  subflow.neu.gz
  cd $TEST_DIR/$LES_PIPE_DIR;            unpack_neu_mesh  pipe.neu.gz

  #-- seq, no cgns
  clean_compile $CONV_DIR no no # dir CGNS MPI

  for CASE_DIR in ${ALL_CONVERT_TESTS[@]}; do
    launch_convert $CASE_DIR
  done

  #-- seq, cgns
  if [ "$CGNS" = "yes" ]; then
    clean_compile $CONV_DIR yes no # dir CGNS MPI

    for CASE_DIR in ${ALL_CONVERT_TESTS[@]}; do
      launch_convert $CASE_DIR
    done
  fi

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

  #-- seq
  clean_compile $DIVI_DIR no no no # dir CGNS MPI

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
    sed -i ""$line_to_replace"s%.*%$new_line%" $3
  fi

  # previous approach was with
  # awk -v var="$line" 'NR==var {$2=3}1'  control > .control.tmp
}

#------------------------------------------------------------------------------#
# processor: backup test
#------------------------------------------------------------------------------#
function process_backup_test {
  # $1 = CGNS = yes
  # $2 = test_dir

  if [ -z "${1+xxx}" ]; then 
    echo "CGNS flag is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  echo "Process backup tests.."

  # save original contol file
  cp $2/control $2/control.backup

  name_in_div=$(head -n1 "$2"/divide.scr)
  nproc_in_div=$(head -n2  "$2"/divide.scr | tail -n1)

  n1=$(printf "%06d" 1)
  n2=$(printf "%06d" 2)

  #----------------------------------------#
  echo "np=1, MPI=no, start from 0, make a backup"
  clean_compile $PROC_DIR $1 no # dir CGNS MPI
  cd $2

  # comment line with LOAD_BACKUP_NAME
  replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
    "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

  # change number of timesteps to 3
  replace_line_with_first_occurence_in_file "NUMBER_OF_TIME_STEPS" \
    "NUMBER_OF_TIME_STEPS 3" control

  # change backup interval to 1 ts
  replace_line_with_first_occurence_in_file "BACKUP_SAVE_INTERVAL" \
    "BACKUP_SAVE_INTERVAL 1" control

  launch_process seq 1
  #----------------------------------------#
  echo "np=1, MPI=no, load from backup(produced by seq)"

  # uncomment line with LOAD_BACKUP_NAME
  replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
    "LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

  launch_process seq 1
  #----------------------------------------#
  echo "np=1, MPI=yes, load from backup(produced by seq)"
  clean_compile $PROC_DIR $1 yes # dir CGNS MPI
  cd $2

  launch_process par 1
  #----------------------------------------#
  echo "np=2, MPI=yes, load from backup(produced by seq)"

  launch_process par $nproc_in_div
  #----------------------------------------#
  echo "np=2, MPI=yes, load from backup(produced by par.np=1)"

  replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
    "LOAD_BACKUP_NAME "$name_in_div"-ts"$n2".backup" control

  launch_process par $nproc_in_div
  #----------------------------------------#
  echo "np=2, MPI=yes, start from 0, make a backup"

  replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
    "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

  launch_process par $nproc_in_div
  #----------------------------------------#
  echo "np=2, MPI=yes, load from backup(produced by par.np=2) "

  replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
    "LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

  launch_process par $nproc_in_div
  #----------------------------------------#
  echo "np=1, MPI=yes, backup=(produced by par.np=2)"
  clean_compile $PROC_DIR $1 no # dir CGNS MPI
  cd $2
  launch_process par 1
  #----------------------------------------#

  cp $2/control.backup $2/control
}

#------------------------------------------------------------------------------#
# process backup tests
#------------------------------------------------------------------------------#
function process_backup_tests {

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor backup tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  # Grasp/embrace as many different model combinations as you can

  echo ""
  echo "#======================================================================"
  echo "#   Test 1: "$RANS_CHANNEL_LR_UNIFORM_DIR" [k_eps model + T]"
  echo "#----------------------------------------------------------------------"
  #-- Channel_Re_Tau_590 [k_eps model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL k_eps" $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR/control
  process_backup_test no  $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test 2: "$RANS_CHANNEL_LR_UNIFORM_DIR" [k_eps_zeta_f model + T]"
  echo "#----------------------------------------------------------------------"
  #-- Channel_Re_Tau_590 [k_eps_zeta_f model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL k_eps_zeta_f" $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR/control
  process_backup_test no  $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $TEST_DIR/$RANS_CHANNEL_LR_UNIFORM_DIR
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test 3: "$RANS_CHANNEL_LR_RSM_DIR" [rsm_hanjalic_jakirlic model + T]"
  echo "#----------------------------------------------------------------------"
  #-- Channel_Re_Tau_590_Rsm [rsm_hanjalic_jakirlic model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL rsm_hanjalic_jakirlic" $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR/control
  process_backup_test no  $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test 4: "$RANS_CHANNEL_LR_RSM_DIR" [rsm_manceau_hanjalic model + T]"
  echo "#----------------------------------------------------------------------"
  #-- Channel_Re_Tau_590_Rsm [rsm_manceau_hanjalic model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL rsm_manceau_hanjalic" $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR/control
  process_backup_test no  $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $TEST_DIR/$RANS_CHANNEL_LR_RSM_DIR
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test 5: "$LES_PIPE_DIR" les_dynamic"
  echo "#----------------------------------------------------------------------"
  #-- Pipe_Re_Tau_180 [les_dynamic]
  process_backup_test no  $TEST_DIR/$LES_PIPE_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $TEST_DIR/$LES_PIPE_DIR
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test 6: "$LAMINAR_CAVITY_LID_DRIVEN_DIR" none"
  echo "#----------------------------------------------------------------------"
  #-- Cavity_Lid_Driven_Re_1000 [none]
  process_backup_test no  $TEST_DIR/$LAMINAR_CAVITY_LID_DRIVEN_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $TEST_DIR/$LAMINAR_CAVITY_LID_DRIVEN_DIR
  fi

}
#------------------------------------------------------------------------------#
# process: save_now / exit_now test
#------------------------------------------------------------------------------#
function process_save_exit_now_test {
  # $1 = CGNS = yes
  # $2 = relative dir with test

  if [ -z "${1+xxx}" ]; then 
    echo "CGNS flag is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  echo ""
  echo "#======================================================================"
  echo "#   Test save_now & exit_now on:" $2
  echo "#----------------------------------------------------------------------"

  cd "$TEST_DIR/$2"
  name_in_div=$(head -n1 divide.scr)
  nproc_in_div=$(head -n2  divide.scr | tail -n1)

  # get rid of "save_now" and "exit_now" files if they happen to exist
  if [ -f "save_now" ]; then rm save_now; fi
  if [ -f "exit_now" ]; then rm exit_now; fi

  # change number of timesteps to 3
  replace_line_with_first_occurence_in_file "NUMBER_OF_TIME_STEPS" \
    "NUMBER_OF_TIME_STEPS 3" control

  # change backup interval to 1 ts
  replace_line_with_first_occurence_in_file "BACKUP_SAVE_INTERVAL" \
    "BACKUP_SAVE_INTERVAL 10" control

  for i in {1..3}
  do
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

    # comment line with LOAD_BACKUP_NAME
    n1=$(printf "%06d" 1)
    replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
      "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

    if [ "$i" = 1 ]; then
      clean_compile $PROC_DIR $1 no
    fi
    if [ "$i" = 2 ]; then
      clean_compile $PROC_DIR $1 yes
    fi
    if [ "$i" = 3 ]; then
      clean_compile $PROC_DIR $1 yes
    fi

    cd $TEST_DIR/$2

    echo "#   Forcing to save: save_now"
    touch save_now

    # get current line count where search starts
    n_start="$(echo "$(wc -l $FULL_LOG | cut -d" " -f1) + 1" | bc -l)"

    # start from scratch
    if [ "$i" = 1 ]; then
      launch_process seq 1
    fi
    if [ "$i" = 2 ]; then
      launch_process par 1
    fi
    if [ "$i" = 3 ]; then
      launch_process par $nproc_in_div
    fi

    # find if save was made in the range [n_start: n_finish]
    if tail -n+$n_start $FULL_LOG | \
      grep -q "# Creating the file: "$name_in_div"-ts"$n1"\|# Creating the file with fields: "$name_in_div"-ts"$n1""; then

      echo "save_now was successfull"

      echo "Forcing to exit: exit_now"
      touch exit_now

      # uncomment line with LOAD_BACKUP_NAME
      replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
        "LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

      # start from ts=1
      n_start="$(echo "$(wc -l $FULL_LOG | cut -d" " -f1) + 1" | bc -l)"

      if [ "$i" = 1 ]; then
        launch_process seq 1
      fi
      if [ "$i" = 2 ]; then
        launch_process par 1
      fi
      if [ "$i" = 3 ]; then
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
    process_save_exit_now_test no $CASE_DIR
  done

  if [ "$CGNS" = "yes" ]; then
    for CASE_DIR in ${ALL_SAVE_EXIT_NOW_TESTS[@]}; do
      process_save_exit_now_test yes $CASE_DIR
    done
  fi
}

#------------------------------------------------------------------------------#
# launch execute script and convert results to .png format
#------------------------------------------------------------------------------#
function launch_gnuplot {
  # $1 = dir with gnuplot script
  # $2 = gnuplot script name
  # $3 = input file
  # $4 = output file

  if [ -z "${1+xxx}" ]; then 
    echo "directory with gnuplot script is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "gnuplot script name is not set at all"
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
  sed -i "s%TEX_FILE_WITH_PLOT_MACRO%$4%" ./tmp
  gnuplot ./tmp >> $FULL_LOG 2>&1
  latex "$4".tex >> $FULL_LOG 2>&1
  dvips -o "$4".ps "$4".dvi >> $FULL_LOG 2>&1
  ps2pdf "$4".ps "$4".pdf >> $FULL_LOG 2>&1
  convert -density 300 "$4".pdf -quality 100 -flatten "$4".png >> $FULL_LOG 2>&1
  echo "new graph was created:" "$1"/"$4".png
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
  name_in_div=$(head -n1 divide.scr)
  nproc_in_div=$(head -n2 divide.scr | tail -n1)

  # rel_dir to User_Mod/ from Process/
  rel_dir=$(realpath --relative-to="$PROC_DIR" "$TEST_DIR/$1")

  if [ "$CGNS" = "yes" ]; then
    clean_compile $PROC_DIR yes yes $rel_dir # dir CGNS MPI DIR_CASE
  else
    clean_compile $PROC_DIR no yes $rel_dir # dir CGNS MPI DIR_CASE
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
  name_in_div=$(head -n1 divide.scr)
  nproc_in_div=$(head -n2 divide.scr | tail -n1)

  echo "np="$nproc_in_div", MPI=yes"
  # rel_dir to User_Mod/ from Process/
  rel_dir=$(realpath --relative-to="$PROC_DIR" "$1")

  if [ "$CGNS" = "yes" ]; then
    clean_compile $PROC_DIR yes yes $rel_dir # dir CGNS MPI DIR_CASE
  else
    clean_compile $PROC_DIR no yes $rel_dir # dir CGNS MPI DIR_CASE
  fi

  cd "$1"

  n1=$(printf "%06d" 1)

  # comment line with LOAD_BACKUP_NAME
  replace_line_with_first_occurence_in_file "LOAD_BACKUP_NAME" \
    "#LOAD_BACKUP_NAME "$name_in_div"-ts"$n1".backup" control

#  # change number of timesteps to 2000
#  replace_line_with_first_occurence_in_file "NUMBER_OF_TIME_STEPS" \
#    "NUMBER_OF_TIME_STEPS 2000" control
#
#  # change backup interval to 100 ts
#  replace_line_with_first_occurence_in_file "BACKUP_SAVE_INTERVAL" \
#    "BACKUP_SAVE_INTERVAL 1000" control

  # change model to $2
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL "$2"" control

  launch_process par $nproc_in_div

  if ls *-res-plus.dat 1> /dev/null 2>&1; then # case-ts??????-res-plus.dat

    # extract essential data from produced .dat files
    last_results_plus_dat_file=$(realpath --relative-to="$3" \
      $(ls -tr1 *-res-plus.dat | tail -n1))

    echo "results are:"
    echo "$(head -n9 $(ls -tr1 *-res-plus.dat | tail -n1))"

    launch_gnuplot "$3" gnuplot_script_template.sh \
      "$last_results_plus_dat_file" "result_plus_"$2""
  else
      echo "Warning: file *-res-plus.dat does not exist"
  fi
}
#------------------------------------------------------------------------------#
# All process tests
#------------------------------------------------------------------------------#
function process_full_length_tests {
  # $1 = dir with test
  # $2 = model
  # $3 = dir with results
  # it requires a new file in Xmgrace/ dir called gnuplot_script_template.sh

  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   Running Processor full simulation tests"
  echo "#"
  echo "#----------------------------------------------------------------------"

  for i in ${!ALL_PROCESS_TESTS[@]}; do
    CASE_DIR="${ALL_PROCESS_TESTS[$i]}"
    CASE_TUR="${ALL_PROCESS_MODELS[$i]}"
    echo ""
    echo "#===================================================================="
    echo "#   Process test: "     $CASE_DIR
    echo "#   with turbulence model: " $CASE_TUR
    echo "#--------------------------------------------------------------------"

    process_full_length_test \
      "$TEST_DIR/$CASE_DIR" \
      "$CASE_TUR" \
      "$TEST_DIR/$CASE_DIR/Xmgrace"
  done

#  # Issue: pipe does not pass process_backup_tests
#  process_full_length_test \
#    "$TEST_DIR/$LES_PIPE_DIR" \
#    "les_dynamic" \
#    "$TEST_DIR/$LES_PIPE_DIR/Xmgrace"
}

#------------------------------------------------------------------------------#
# actual script
#------------------------------------------------------------------------------#
while [ 0 -eq 0 ]; do
  echo ""
  echo "#======================================================================"
  echo "#"
  echo "#   T-Flows testing"
  echo "#"
  echo "#----------------------------------------------------------------------"
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
  echo "  8. Perform all tests"
  echo "  9. Clean all test directories"
  echo ""
  read -p "  Enter the desired type of test: " option
  if [ $option -eq 0 ]; then exit 1;                       fi
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
    process_backup_tests;
  fi
  if [ $option -eq 6 ]; then
    if [ $DONE_GENERATE_TESTS -eq 0 ]; then generate_tests; fi
    if [ $DONE_CONVERT_TESTS  -eq 0 ]; then convert_tests;  fi
    if [ $DONE_DIVIDE_TESTS   -eq 0 ]; then divide_tests;   fi
    process_save_exit_now_tests;
  fi
  if [ $option -eq 7 ]; then
    if [ $DONE_GENERATE_TESTS -eq 0 ]; then generate_tests; fi
    if [ $DONE_CONVERT_TESTS  -eq 0 ]; then convert_tests;  fi
    if [ $DONE_DIVIDE_TESTS   -eq 0 ]; then divide_tests;   fi
    process_full_length_tests;
  fi
  if [ $option -eq 8 ]; then
    generate_tests
    convert_tests
    divide_tests
    process_compilation_tests
    process_backup_tests
    process_save_exit_now_tests
    process_full_length_tests
  fi
  if [ $option -eq 9 ]; then
    git clean -dfx ./
  fi
done

