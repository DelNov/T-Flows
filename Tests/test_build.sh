#!/bin/bash

# Description: This script is made for debug purposes.
# It automatically builds and runs most cases in T-Flows

# Requires: gnuplot, texlive-base, mpi, gfortran

# Compilation flags used in makefiles
FCOMP="gnu"
# Conduct tests with DEBUG=yes
DEBUG="no"
# Repeat tests with CGNS_HDF5=yes
CGNS="no"

# A small reminder how to set up alternatives if you have mpich and openmpi:
#update-alternatives --install /usr/bin/mpif90 mpif90 /usr/bin/mpif90.openmpi 20
#update-alternatives --install /usr/bin/mpif90 mpif90 /usr/bin/mpif90.mpich   80

#update-alternatives --install /usr/bin/mpirun mpirun /usr/bin/mpirun.openmpi 20
#update-alternatives --install /usr/bin/mpirun mpirun /usr/bin/mpirun.mpich   80

# Folder structure
TEST_DIR=$PWD                      # dir with tests
GENE_DIR=$PWD/../Sources/Generate  # Generate src folder
CONV_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIVI_DIR=$PWD/../Sources/Divide    # Divide   src folder
PROC_DIR=$PWD/../Sources/Process   # Process  src folder
BINA_DIR=$PWD/../Binaries/         # binaries folder

# Executables
GENE_EXE=$BINA_DIR/Generate        # Generate ex
CONV_EXE=$BINA_DIR/Convert         # Convert  ex
DIVI_EXE=$BINA_DIR/Divide          # Divide   ex
PROC_EXE=$BINA_DIR/Process         # Process  ex

# Folders with geometry
# Generator file (.dom)
LAMINAR_BACKSTEP_ORTH_DIR=$TEST_DIR/Laminar/Backstep/Orthogonal
LAMINAR_BACKSTEP_NON_ORTH_DIR=$TEST_DIR/Laminar/Backstep/Nonorthogonal

LES_CAVITY_LID_DRIVEN_DIR=$TEST_DIR/Laminar/Cavity/Lid_Driven/Re_1000
LES_CAVITY_THERM_DRIVEN_DIR_106=$TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e6
LES_CAVITY_THERM_DRIVEN_DIR_108=$TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e8

RANS_BACKSTEP_5100_DIR=$TEST_DIR/Rans/Backstep_Re_05100
RANS_BACKSTEP_28000_DIR=$TEST_DIR/Rans/Backstep_Re_28000

RANS_CHANNEL_LR_LONG_DIR=$TEST_DIR/Rans/Channel_Re_Tau_590/Long_Domain
RANS_CHANNEL_LR_RSM_DIR=$TEST_DIR/Rans/Channel_Re_Tau_590/Rsm
RANS_CHANNEL_LR_STRETCHED_DIR=$TEST_DIR/Rans/Channel_Re_Tau_590/Stretched_Mesh
RANS_CHANNEL_LR_UNIFORM_DIR=$TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh

HYB_CHANNEL_HR_UNIFORM_DIR=\
$TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh
HYB_CHANNEL_HR_STRETCHED_DIR=\
$TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh

# Mesh file (.cgns/.neu)
LES_PIPE_DIR=$TEST_DIR/Les/Pipe_Re_Tau_180

RANS_IMPINGING_JET_DIR=$TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000
RANS_FUEL_BUNDLE_DIR=$TEST_DIR/Rans/Fuel_Bundle
#RANS_PIPE_DIR=$TEST_DIR/Rans/Pipe_Re_Tau_2000
#------------------------------------------------------------------------------#

# Start time measurements from this moment
current_time=$(date +%s)

# Script logs
FULL_LOG=$TEST_DIR/test_build.log # logs of current script
if [ -f $FULL_LOG ]; then cp /dev/null $FULL_LOG; fi
echo "Full log is being written in file" "$FULL_LOG"

# exit when any command fails
set -e

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
#------------------------------------------------------------------------------#
function clean_compile {
  # $1 = dir 
  # $2 = CGNS_HDF5 = yes/no
  # $3 = MPI = yes/no
  # $4 = DIR_CASE path

  if [ -z "${1+xxx}" ]; then 
    echo "directory with sources is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "CGNS_HDF5 flag is not set at all"
    exit 1
  elif [ -z "${3+xxx}" ]; then 
    echo "MPI flag is not set at all"
    exit 1
  fi

  cd $1
  echo "clean compile in:" "$1"
  make clean >> $FULL_LOG 2>&1

  if [ -z "${4+xxx}" ]; then 
    echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3"
          make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 >> $FULL_LOG 2>&1
  else
    echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 DIR_CASE=$4"
          make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 DIR_CASE=$4 \
          >> $FULL_LOG 2>&1
  fi

  time_in_seconds
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
#------------------------------------------------------------------------------#
function launch_generate {
  echo "launching Generate in: " $PWD
  echo "generate.scr: " >> $FULL_LOG 2>&1
  cat generate.scr >> $FULL_LOG 2>&1
  $GENE_EXE < generate.scr >> $FULL_LOG 2>&1
  time_in_seconds
}
#------------------------------------------------------------------------------#
# launch Divide
#------------------------------------------------------------------------------#
function launch_divide {
  echo "launching Divide in: " $PWD
  echo "divide.scr: " >> $FULL_LOG 2>&1
  cat divide.scr >> $FULL_LOG 2>&1
  $DIVI_EXE < divide.scr >> $FULL_LOG 2>&1
  time_in_seconds
}
#------------------------------------------------------------------------------#
# launch Convert
#------------------------------------------------------------------------------#
function launch_convert {
  echo "launching Convert in: " $PWD
  echo "convert.scr: " >> $FULL_LOG 2>&1
  cat convert.scr >> $FULL_LOG 2>&1
  $CONV_EXE < convert.scr >> $FULL_LOG 2>&1
  time_in_seconds
}
#------------------------------------------------------------------------------#
# launch Process
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
    echo "launching Process in: " $PWD
    $PROC_EXE >> $FULL_LOG 2>&1
  elif [ "$1" == "par" ]; then
    echo "launching mpirun -np "$2" Process in: " $PWD
    mpirun -np $2 $PROC_EXE >> $FULL_LOG 2>&1
  else
    echo "wrong argument in launch_process"
    exit 1
  fi

  echo "control: " >> $FULL_LOG 2>&1
  cat control >> $FULL_LOG 2>&1
  time_in_seconds
}
#------------------------------------------------------------------------------#
# generator tests
#------------------------------------------------------------------------------#
function generator_tests {
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!"
  echo "!!"
  echo "!!    Running Generate Tests"
  echo "!!"
  echo "!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  #-- seq, no cgns
  clean_compile $GENE_DIR  no  no  # dir CGNS_HDF5 MPI

  make_links $LAMINAR_BACKSTEP_ORTH_DIR;        launch_generate
  make_links $LAMINAR_BACKSTEP_NON_ORTH_DIR;    launch_generate

  make_links $LES_CAVITY_LID_DRIVEN_DIR;        launch_generate
  make_links $LES_CAVITY_THERM_DRIVEN_DIR_106;  launch_generate
  make_links $LES_CAVITY_THERM_DRIVEN_DIR_108;  launch_generate

  make_links $RANS_BACKSTEP_5100_DIR;           launch_generate
  make_links $RANS_BACKSTEP_28000_DIR;          launch_generate

  make_links $RANS_CHANNEL_LR_LONG_DIR;         launch_generate
  make_links $RANS_CHANNEL_LR_RSM_DIR;          launch_generate
  make_links $RANS_CHANNEL_LR_STRETCHED_DIR;    launch_generate
  make_links $RANS_CHANNEL_LR_UNIFORM_DIR;      launch_generate

  make_links $HYB_CHANNEL_HR_STRETCHED_DIR;     launch_generate
  make_links $HYB_CHANNEL_HR_UNIFORM_DIR;       launch_generate

  #-- seq, cgns(hdf5)
  if [ "$CGNS" = "yes" ]; then
    clean_compile $GENE_DIR  yes  no  # dir CGNS_HDF5 MPI

    make_links $LAMINAR_BACKSTEP_ORTH_DIR;        launch_generate
    make_links $LAMINAR_BACKSTEP_NON_ORTH_DIR;    launch_generate

    make_links $LES_CAVITY_LID_DRIVEN_DIR;        launch_generate
    make_links $LES_CAVITY_THERM_DRIVEN_DIR_106;  launch_generate
    make_links $LES_CAVITY_THERM_DRIVEN_DIR_108;  launch_generate

    make_links $RANS_BACKSTEP_5100_DIR;           launch_generate
    make_links $RANS_BACKSTEP_28000_DIR;          launch_generate

    make_links $RANS_CHANNEL_LR_LONG_DIR;         launch_generate
    make_links $RANS_CHANNEL_LR_RSM_DIR;          launch_generate
    make_links $RANS_CHANNEL_LR_STRETCHED_DIR;    launch_generate
    make_links $RANS_CHANNEL_LR_UNIFORM_DIR;      launch_generate

    make_links $HYB_CHANNEL_HR_STRETCHED_DIR;     launch_generate
    make_links $HYB_CHANNEL_HR_UNIFORM_DIR;       launch_generate
  fi
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
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!"
  echo "!!"
  echo "!!    Running Convert Tests"
  echo "!!"
  echo "!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  cd $RANS_IMPINGING_JET_DIR;  unpack_neu_mesh  jet.neu
  cd $RANS_FUEL_BUNDLE_DIR;    unpack_neu_mesh  subflow.neu.gz
  cd $LES_PIPE_DIR;            unpack_neu_mesh  pipe.neu.gz

  #cd $RANS_PIPE_DIR;           unpack_neu_mesh  pipe.neu.gz

  #-- seq, no cgns
  clean_compile $CONV_DIR no no # dir CGNS_HDF5 MPI

  make_links $RANS_IMPINGING_JET_DIR;  launch_convert
  make_links $RANS_FUEL_BUNDLE_DIR;    launch_convert
  make_links $LES_PIPE_DIR;            launch_convert

  #-- seq, cgns_hdf5
  if [ "$CGNS" = "yes" ]; then
    clean_compile $CONV_DIR yes no # dir CGNS_HDF5 MPI

    make_links $RANS_IMPINGING_JET_DIR;  launch_convert
    make_links $RANS_FUEL_BUNDLE_DIR;    launch_convert
    make_links $LES_PIPE_DIR;            launch_convert
  fi
}
#------------------------------------------------------------------------------#
# Divide tests
#------------------------------------------------------------------------------#
function divide_tests {
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!"
  echo "!!"
  echo "!!    Running Divide Tests"
  echo "!!"
  echo "!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  #-- seq
  clean_compile $DIVI_DIR no no no # dir CGNS_HDF5 MPI
  cd $LAMINAR_BACKSTEP_ORTH_DIR;        launch_divide
  cd $LAMINAR_BACKSTEP_NON_ORTH_DIR;    launch_divide
  cd $LES_CAVITY_LID_DRIVEN_DIR;        launch_divide
  cd $LES_CAVITY_THERM_DRIVEN_DIR_106;  launch_divide
  cd $LES_CAVITY_THERM_DRIVEN_DIR_108;  launch_divide
  cd $RANS_BACKSTEP_5100_DIR;           launch_divide
  cd $RANS_BACKSTEP_28000_DIR;          launch_divide
  cd $RANS_CHANNEL_LR_LONG_DIR;         launch_divide
  cd $RANS_CHANNEL_LR_RSM_DIR;          launch_divide
  cd $RANS_CHANNEL_LR_STRETCHED_DIR;    launch_divide
  cd $RANS_CHANNEL_LR_UNIFORM_DIR;      launch_divide
  cd $RANS_IMPINGING_JET_DIR;           launch_divide
  cd $RANS_FUEL_BUNDLE_DIR;             launch_divide
  cd $HYB_CHANNEL_HR_STRETCHED_DIR;     launch_divide
  cd $HYB_CHANNEL_HR_UNIFORM_DIR;       launch_divide
  cd $LES_PIPE_DIR;                     launch_divide

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
  # $1 = CGNS_HDF5 = yes
  # $2 = test_dir

  if [ -z "${1+xxx}" ]; then 
    echo "CGNS_HDF5 flag is not set at all"
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
  clean_compile $PROC_DIR $1 no # dir CGNS_HDF5 MPI
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
  clean_compile $PROC_DIR $1 yes # dir CGNS_HDF5 MPI
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
  clean_compile $PROC_DIR $1 no # dir CGNS_HDF5 MPI
  cd $2
  launch_process par 1
  #----------------------------------------#

  cp $2/control.backup $2/control
}
#------------------------------------------------------------------------------#
# processor backup tests
#------------------------------------------------------------------------------#
function processor_backup_tests {

  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "  !!"
  echo "  !!    Running Processor Backup Tests"
  echo "  !!"
  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  # Grasp/embrace as many different model combinations as you can

  echo "================================= TEST 1 =============================="
  #-- Channel_Re_Tau_590 [k_eps model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL k_eps" $RANS_CHANNEL_LR_UNIFORM_DIR/control
  process_backup_test no  $RANS_CHANNEL_LR_UNIFORM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $RANS_CHANNEL_LR_UNIFORM_DIR
  fi

  echo "================================= TEST 2 =============================="
  #-- Channel_Re_Tau_590 [k_eps_zeta_f model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL k_eps_zeta_f" $RANS_CHANNEL_LR_UNIFORM_DIR/control
  process_backup_test no  $RANS_CHANNEL_LR_UNIFORM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $RANS_CHANNEL_LR_UNIFORM_DIR
  fi

  echo "================================= TEST 3 =============================="
  #-- Channel_Re_Tau_590_Rsm [rsm_hanjalic_jakirlic model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL rsm_hanjalic_jakirlic" $RANS_CHANNEL_LR_RSM_DIR/control
  process_backup_test no  $RANS_CHANNEL_LR_RSM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $RANS_CHANNEL_LR_RSM_DIR
  fi

  echo "================================= TEST 4 =============================="
  #-- Channel_Re_Tau_590_Rsm [rsm_manceau_hanjalic model + T]
  replace_line_with_first_occurence_in_file "TURBULENCE_MODEL" \
    "TURBULENCE_MODEL rsm_manceau_hanjalic" $RANS_CHANNEL_LR_RSM_DIR/control
  process_backup_test no  $RANS_CHANNEL_LR_RSM_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $RANS_CHANNEL_LR_RSM_DIR
  fi

  echo "================================= TEST 5 =============================="
  #-- Pipe_Re_Tau_180 [les_dynamic]
  process_backup_test no  $LES_PIPE_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $LES_PIPE_DIR
  fi

  echo "================================= TEST 6 =============================="
  #-- Cavity_Lid_Driven_Re_1000 [none]
  process_backup_test no  $LES_CAVITY_LID_DRIVEN_DIR
  if [ "$CGNS" = "yes" ]; then
    process_backup_test yes $LES_CAVITY_LID_DRIVEN_DIR
  fi

}
#------------------------------------------------------------------------------#
# processor: save_now / exit_now test
#------------------------------------------------------------------------------#
function process_save_exit_now_test {
  # $1 = CGNS_HDF5 = yes
  # $2 = dir with test

  if [ -z "${1+xxx}" ]; then 
    echo "CGNS_HDF5 flag is not set at all"
    exit 1
  elif [ -z "${2+xxx}" ]; then 
    echo "directory is not set at all"
    exit 1
  fi

  echo "Test: save_now & exit now on " $2

  cd "$2"
  name_in_div=$(head -n1 divide.scr)
  nproc_in_div=$(head -n2  divide.scr | tail -n1)

  # change number of timesteps to 3
  replace_line_with_first_occurence_in_file "NUMBER_OF_TIME_STEPS" \
    "NUMBER_OF_TIME_STEPS 3" control

  # change backup interval to 1 ts
  replace_line_with_first_occurence_in_file "BACKUP_SAVE_INTERVAL" \
    "BACKUP_SAVE_INTERVAL 1" control

  echo "================================= TEST 1 =============================="
  echo "np=1, MPI=no"
  clean_compile $PROC_DIR $1 no # dir CGNS_HDF5 MPI
  cd $2

  echo "save_now"
  touch save_now

  if launch_process seq 1 | grep -q ""$name_in_div"-ts000001"; then
    echo "exit_now"
    touch exit_now
    launch_process seq 1 | grep -q "# Exiting !"
    echo "save_exit_now_test was successfull"
  fi

  echo "================================= TEST 2 =============================="
  echo "np=1, MPI=yes"
  clean_compile $PROC_DIR $1 yes # dir CGNS_HDF5 MPI
  cd $2

  echo "save_now"
  touch save_now

  if launch_process par 1 | grep -q ""$name_in_div"-ts000001"; then
    echo "exit_now"
    touch exit_now
    launch_process par 1 | grep -q "# Exiting !"
    echo "save_exit_now_test was successfull"
  fi

  echo "================================= TEST 3 =============================="
  echo "np=2, MPI=yes"

  echo "save_now"
  touch save_now

  if launch_process par $nproc_in_div | grep -q ""$name_in_div"-ts000001"; then
    echo "exit_now"
    touch exit_now
    launch_process par $nproc_in_div | grep -q "# Exiting !"
    echo "save_exit_now_test was successfull"
  fi
}
#------------------------------------------------------------------------------#
# processor save_now / exit_now tests
#------------------------------------------------------------------------------#
function process_save_exit_now_tests {

  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "  !!"
  echo "  !!    Running Processor save_now and exit_now Tests"
  echo "  !!"
  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  process_save_exit_now_test no  $LAMINAR_BACKSTEP_ORTH_DIR
  if [ "$CGNS" = "yes" ]; then
    process_save_exit_now_test yes $LAMINAR_BACKSTEP_ORTH_DIR
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
# Individual processor tests for compilation
#------------------------------------------------------------------------------#
function processor_compilation_test {
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

  echo "Test:  compilation in " "$1" with "$2" model

  echo "np="$nproc_in_div", MPI=yes"
  # rel_dir to User_Mod/ from Process/
  rel_dir=$(realpath --relative-to="$PROC_DIR" "$1")

  if [ "$CGNS" = "yes" ]; then
    clean_compile $PROC_DIR yes yes $rel_dir # dir CGNS_HDF5 MPI DIR_CASE
  else
    clean_compile $PROC_DIR no yes $rel_dir # dir CGNS_HDF5 MPI DIR_CASE
  fi
}
#------------------------------------------------------------------------------#
# All processor compilation tests
#------------------------------------------------------------------------------#
function processor_compilation_tests {
  # $1 = dir with test
  # $2 = model
  # $3 = dir with results
  # it requires a new file in Xmgrace/ dir called gnuplot_script_template.sh

  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "  !!"
  echo "  !!    Running Processor compilation tests"
  echo "  !!"
  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  echo "================================= TEST 1 =============================="
  # no User_Mod/ dir !!!
  processor_compilation_test \
    "$LES_CAVITY_LID_DRIVEN_DIR" \
    "none" \
    "$LES_CAVITY_LID_DRIVEN_DIR/Xmgrace"

  echo "================================= TEST 2 =============================="
  # no User_Mod/ dir !!!
  processor_compilation_test \
    "$LES_CAVITY_THERM_DRIVEN_DIR_106" \
    "none" \
    "$LES_CAVITY_THERM_DRIVEN_DIR_106/Xmgrace"

  echo "================================= TEST 3 =============================="
  # no User_Mod/ dir !!!
  processor_compilation_test \
    "$LES_CAVITY_THERM_DRIVEN_DIR_108" \
    "none" \
    "$LES_CAVITY_THERM_DRIVEN_DIR_108/Xmgrace"

  echo "================================= TEST 4 =============================="
  # User_Mod/ dir exists
  processor_compilation_test \
    "$RANS_CHANNEL_LR_UNIFORM_DIR" \
    "k_eps" \
    "$RANS_CHANNEL_LR_UNIFORM_DIR/Xmgrace"

  echo "================================= TEST 5 =============================="
  # User_Mod/ dir exists
  processor_compilation_test \
    "$RANS_CHANNEL_LR_STRETCHED_DIR" \
    "k_eps_zeta_f" \
    "$RANS_CHANNEL_LR_STRETCHED_DIR/Xmgrace"

  echo "================================= TEST 6 =============================="
  # User_Mod/ dir exists
  processor_compilation_test \
    "$RANS_CHANNEL_LR_RSM_DIR" \
    "rsm_manceau_hanjalic" \
    "$RANS_CHANNEL_LR_RSM_DIR/Xmgrace"

  echo "================================= TEST 7 =============================="
  # User_Mod/ dir exists
  processor_compilation_test \
    "$RANS_CHANNEL_LR_RSM_DIR" \
    "rsm_hanjalic_jakirlic" \
    "$RANS_CHANNEL_LR_RSM_DIR/Xmgrace"

  echo "================================= TEST 8 =============================="
  # User_Mod/ dir exists
  processor_compilation_test \
    "$HYB_CHANNEL_HR_UNIFORM_DIR" \
    "hybrid_les_rans" \
    "$HYB_CHANNEL_HR_UNIFORM_DIR/Xmgrace"

  echo "================================= TEST 9 =============================="
  # User_Mod/ dir exists
  processor_compilation_test \
    "$HYB_CHANNEL_HR_STRETCHED_DIR" \
    "hybrid_les_rans" \
    "$HYB_CHANNEL_HR_STRETCHED_DIR/Xmgrace"

#  # Issue: pipe does not pass processor_backup_tests
#  processor_compilation_test \
#    "$LES_PIPE_DIR" \
#    "les_dynamic" \
#    "$LES_PIPE_DIR/Xmgrace"
}
#------------------------------------------------------------------------------#
# Individual processor tests for full length
#------------------------------------------------------------------------------#
function processor_full_length_test {
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

  echo "Test: full length simulation in " "$1" with "$2" model

  echo "np="$nproc_in_div", MPI=yes"
  # rel_dir to User_Mod/ from Process/
  rel_dir=$(realpath --relative-to="$PROC_DIR" "$1")

  if [ "$CGNS" = "yes" ]; then
    clean_compile $PROC_DIR yes yes $rel_dir # dir CGNS_HDF5 MPI DIR_CASE
  else
    clean_compile $PROC_DIR no yes $rel_dir # dir CGNS_HDF5 MPI DIR_CASE
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
# All processor tests
#------------------------------------------------------------------------------#
function processor_full_length_tests {
  # $1 = dir with test
  # $2 = model
  # $3 = dir with results
  # it requires a new file in Xmgrace/ dir called gnuplot_script_template.sh

  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "  !!"
  echo "  !!    Running Processor full simulation tests"
  echo "  !!"
  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  echo "================================= TEST 1 =============================="
  # no User_Mod/ dir !!!
  processor_full_length_test \
    "$LES_CAVITY_LID_DRIVEN_DIR" \
    "none" \
    "$LES_CAVITY_LID_DRIVEN_DIR/Xmgrace"

  echo "================================= TEST 2 =============================="
  # no User_Mod/ dir !!!
  processor_full_length_test \
    "$LES_CAVITY_THERM_DRIVEN_DIR_106" \
    "none" \
    "$LES_CAVITY_THERM_DRIVEN_DIR_106/Xmgrace"

  echo "================================= TEST 3 =============================="
  # no User_Mod/ dir !!!
  processor_full_length_test \
    "$LES_CAVITY_THERM_DRIVEN_DIR_108" \
    "none" \
    "$LES_CAVITY_THERM_DRIVEN_DIR_108/Xmgrace"

  echo "================================= TEST 4 =============================="
  # [~2 min test]
  processor_full_length_test \
    "$RANS_CHANNEL_LR_UNIFORM_DIR" \
    "k_eps" \
    "$RANS_CHANNEL_LR_UNIFORM_DIR/Xmgrace"

  echo "================================= TEST 5 =============================="
  # [~2 min test]
  processor_full_length_test \
    "$RANS_CHANNEL_LR_STRETCHED_DIR" \
    "k_eps_zeta_f" \
    "$RANS_CHANNEL_LR_STRETCHED_DIR/Xmgrace"

  echo "================================= TEST 6 =============================="
  # [~5 min test]
  processor_full_length_test \
    "$RANS_CHANNEL_LR_RSM_DIR" \
    "rsm_manceau_hanjalic" \
    "$RANS_CHANNEL_LR_RSM_DIR/Xmgrace"

  echo "================================= TEST 7 =============================="
  # [~5 min test]
  processor_full_length_test \
    "$RANS_CHANNEL_LR_RSM_DIR" \
    "rsm_hanjalic_jakirlic" \
    "$RANS_CHANNEL_LR_RSM_DIR/Xmgrace"

  echo "================================= TEST 8 =============================="
  # [~4.5 HOURS test]
  processor_full_length_test \
    "$HYB_CHANNEL_HR_UNIFORM_DIR" \
    "hybrid_les_rans" \
    "$HYB_CHANNEL_HR_UNIFORM_DIR/Xmgrace"

  echo "================================= TEST 9 =============================="
  # [~4 HOURS test]
  processor_full_length_test \
    "$HYB_CHANNEL_HR_STRETCHED_DIR" \
    "hybrid_les_rans" \
    "$HYB_CHANNEL_HR_STRETCHED_DIR/Xmgrace"

#  # Issue: pipe does not pass processor_backup_tests
#  processor_full_length_test \
#    "$LES_PIPE_DIR" \
#    "les_dynamic" \
#    "$LES_PIPE_DIR/Xmgrace"
}
#------------------------------------------------------------------------------#
# actual script
#------------------------------------------------------------------------------#
# generator_tests
# convert_tests
# divide_tests
processor_compilation_tests
# processor_backup_tests
# process_save_exit_now_tests
# processor_full_length_tests
