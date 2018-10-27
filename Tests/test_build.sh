#!/bin/bash

# it is an useful script in debug perposus to automatically build and run most 
# cases in T-Flows

# put your compilers here
FCOMP="gnu"; # or ifort/gfortran/mpif90/mpifort/mpiifort
DEBUG="no";  # run tests in debug mode

# Folder structure
TEST_DIR=$PWD                      # dir with tests
GEN_DIR=$PWD/../Sources/Generate  # Generate src folder
CON_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIV_DIR=$PWD/../Sources/Divide    # Divide   src folder
PRO_DIR=$PWD/../Sources/Process   # Process  src folder
BIN_DIR=$PWD/../Binaries/         # binaries folder

# Executables
GEN_EXE=$BIN_DIR/Generate        # Generate ex
CON_EXE=$BIN_DIR/Convert         # Convert  ex
DIV_EXE=$BIN_DIR/Divide          # Divide   ex
PRO_EXE=$BIN_DIR/Process         # Process  ex

#-------------------------------------------------#
#---------   READ ABOVE UP TO THIS ROW   ---------#
#-------------------------------------------------#

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#------------------------------------------------------------------------------#
# chd
#------------------------------------------------------------------------------#
function chd {
  # $1 = dir

  dir=$1
  dir=${dir//$TEST_DIR\//}

  echo "    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "    !!  Current directory: " $dir
  echo "    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  cd $1
}

#------------------------------------------------------------------------------#
# clean_compile
#------------------------------------------------------------------------------#
function clean_compile {
  # $1 = dir 
  # $2 = CGNS_HDF5 = yes/no
  # $3 = MPI = yes/no

  chd $1
  make clean
  echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3"
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3
}

#------------------------------------------------------------------------------#
# clean_compile_user
#------------------------------------------------------------------------------#
function clean_compile_case {
  # $1 = dir
  # $2 = CGNS_HDF5 = yes/no
  # $3 = MPI = yes/no
  # $4 = case_dir

  chd $1
  make clean
  echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 CASE_DIR=$4"
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 CASE_DIR=$4
}

#------------------------------------------------------------------------------#
# Make links
#------------------------------------------------------------------------------#
function make_links {

  chd $TEST_DIR/Laminar/Backstep/Orthogonal/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Laminar/Backstep/Nonorthogonal/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Laminar/Cavity/Lid_Driven/Re_1000/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e6/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Backstep_Re_26000_Rsm/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Backstep_Re_28000/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Long_Domain/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Rsm/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Stretched_Mesh/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Rans/Fuel_Bundle/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
  chd $TEST_DIR/Les/Pipe_Re_Tau_180/;
  ln -rsf $GEN_EXE $CON_EXE $DIV_EXE $PRO_EXE .
}

#------------------------------------------------------------------------------#
# Generator tests
#------------------------------------------------------------------------------#
function generator_tests {

  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!"
  echo "!!"
  echo "!!    Running Generator Tests"
  echo "!!"
  echo "!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  # Seq, no CGNS
  clean_compile $GEN_DIR no no  # dir CGNS_HDF5 MPI
  chd $TEST_DIR/Laminar/Backstep/Orthogonal;             $GEN_EXE < generate.scr
  chd $TEST_DIR/Laminar/Backstep/Nonorthogonal;          $GEN_EXE < generate.scr
  chd $TEST_DIR/Laminar/Cavity/Lid_Driven/Re_1000;       $GEN_EXE < generate.scr
  chd $TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e6; $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Backstep_Re_28000;                  $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Long_Domain;     $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Rsm;             $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Stretched_Mesh;  $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh;    $GEN_EXE < generate.scr
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh/;
  $GEN_EXE < generate.scr
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh/;
  $GEN_EXE < generate.scr

  # Seq, CGNS(HDF5)
  clean_compile $GEN_DIR yes no  # dir CGNS_HDF5 MPI
  chd $TEST_DIR/Laminar/Backstep/Orthogonal;             $GEN_EXE < generate.scr
  chd $TEST_DIR/Laminar/Backstep/Nonorthogonal;          $GEN_EXE < generate.scr
  chd $TEST_DIR/Laminar/Cavity/Lid_Driven/Re_1000;       $GEN_EXE < generate.scr
  chd $TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e6; $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Backstep_Re_28000;                  $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Long_Domain;     $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Rsm;             $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Stretched_Mesh;  $GEN_EXE < generate.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh;    $GEN_EXE < generate.scr
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh/;
  $GEN_EXE < generate.scr
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh/;
  $GEN_EXE < generate.scr
}

#------------------------------------------------------------------------------#
# Convert tests
#------------------------------------------------------------------------------#
function convert_tests {

  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!"
  echo "!!"
  echo "!!    Running Convert Tests"
  echo "!!"
  echo "!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  # Unpacking geometry
  chd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;
  if [ -f jet.neu.gz ]; then
    gunzip jet.neu.gz 
  fi
  chd $TEST_DIR/Rans/Fuel_Bundle;
  if [ -f subflow.neu.gz ]; then
    gunzip subflow.neu.gz
  fi
  chd $TEST_DIR/Les/Pipe_Re_Tau_180;
  if [ -f pipe.neu.gz ]; then
    gunzip pipe.neu.gz
  fi

  # Seq, no CGNS
  clean_compile $CON_DIR no no  # dir CGNS_HDF5 MPI
  chd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CON_EXE < convert.scr
  chd $TEST_DIR/Rans/Fuel_Bundle;                         $CON_EXE < convert.scr
  chd $TEST_DIR/Les/Pipe_Re_Tau_180;                      $CON_EXE < convert.scr

  # Seq, CGNS(HDF5)
  clean_compile $CON_DIR yes no  # dir CGNS_HDF5 MPI
  chd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CON_EXE < convert.scr
  chd $TEST_DIR/Rans/Fuel_Bundle;                         $CON_EXE < convert.scr
  chd $TEST_DIR/Les/Pipe_Re_Tau_180;                      $CON_EXE < convert.scr
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

  # Seq
  clean_compile $DIV_DIR no no  # dir CGNS_HDF5 MPI
  chd $TEST_DIR/Laminar/Backstep/Orthogonal;               $DIV_EXE < divide.scr
  chd $TEST_DIR/Laminar/Backstep/Nonorthogonal;            $DIV_EXE < divide.scr
  chd $TEST_DIR/Laminar/Cavity/Lid_Driven/Re_1000;         $DIV_EXE < divide.scr
  chd $TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e6;   $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;                $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Backstep_Re_28000;                    $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Long_Domain;       $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Rsm;               $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Stretched_Mesh;    $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh;      $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;    $DIV_EXE < divide.scr
  chd $TEST_DIR/Rans/Fuel_Bundle;                          $DIV_EXE < divide.scr
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh/;
  $DIV_EXE < divide.scr
  chd $TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh/;
  $DIV_EXE < divide.scr
  chd $TEST_DIR/Les/Pipe_Re_Tau_180/;                      $DIV_EXE < divide.scr
}

#------------------------------------------------------------------------------#
# Processor: backup test
#------------------------------------------------------------------------------#
function backup_test {
  # $1 = CGNS_HDF5 = yes
  # $2 = test_dir

  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "  !!"
  echo "  !!    Running Processor Backup Tests"
  echo "  !!"
  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  #----------------------------------#
  # np=1, MPI=no, backup=no
  #----------------------------------#
  echo "np=1, MPI=no, backup=no"
  clean_compile $PRO_DIR $1 no  # dir CGNS_HDF5 MPI
  chd $2;
  cp control control.backup
  line="$(grep -ni "NUMBER_OF_TIME_STEPS"       control | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=3}1'     control > .control.tmp
  line="$(grep -ni "BACKUP_SAVE_INTERVAL"  .control.tmp | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=1}1'     .control.tmp > control
  $PRO_EXE

  #----------------------------------#
  # np=1, MPI=no, backup=(from np=1)
  #----------------------------------#
  echo "np=1, MPI=no, backup=(from np=1)"
  sed -e  's/#LOAD_BACKUP_NAME/LOAD_BACKUP_NAME/g' .control.tmp > control
  $PRO_EXE

  #----------------------------------#
  # np=1, MPI=yes, backup=(from np=1)
  #----------------------------------#
  echo "np=1, MPI=yes, backup=(from np=1)"
  clean_compile $PRO_DIR $1 yes  # dir CGNS_HDF5 MPI
  chd $2;
  mpirun -np 1 $PRO_EXE

  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=1)
  #----------------------------------#
  nproc_in_div=$(tail -1  divide.scr)
  echo "np=2, MPI=yes, backup=(from np=1)"
  mpirun -np $nproc_in_div $PRO_EXE

  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=2)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=2)"
  line="$(grep -ni "BACKUP_SAVE_INTERVAL"  .control.tmp | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=2}1'    .control.tmp > control
  mpirun -np $nproc_in_div $PRO_EXE

  #----------------------------------#
  # np=2, MPI=yes, backup=no
  #----------------------------------#
  echo "np=2, MPI=yes, backup=no"
  sed -e  's/LOAD_BACKUP_NAME/#LOAD_BACKUP_NAME/g' .control.tmp > control
  mpirun -np $nproc_in_div $PRO_EXE

  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=2)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=2)"
  sed -e  's/#LOAD_BACKUP_NAME/LOAD_BACKUP_NAME/g' .control.tmp > control
  mpirun -np $nproc_in_div $PRO_EXE

  #----------------------------------#
  # np=1, MPI=no, backup=(from np=2)
  #----------------------------------#
  echo "np=1, MPI=yes, backup=(from np=2)"
  clean_compile $PRO_DIR $1 no  # dir CGNS_HDF5 MPI
  chd $2;
  $PRO_EXE

  cp control.backup control
  rm control.backup .control.tmp
}

#------------------------------------------------------------------------------#
# Processor: save_now / exit_now test
#------------------------------------------------------------------------------#
function save_exit_now_test {
  # $1 = CGNS_HDF5 = yes
  # $2 = test_dir

  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "  !!"
  echo "  !!    Running Processor \"save_now\" and \"exit_now\" Tests"
  echo "  !!"
  echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  #---------------#
  # np=1, MPI=no
  #---------------#
  echo "np=1, MPI=no"
  clean_compile $PRO_DIR $1 no  # dir CGNS_HDF5 MPI CGNS_ADF5
  chd $2;
  echo  save_now
  touch save_now
  if $PRO_EXE | grep -q "chan-ts000001"; then # exit if match is found
    echo  exit_now
    touch exit_now
    $PRO_EXE | grep -q "# Exiting !"
    echo save_exit_now_test was successfull
  fi

  #---------------#
  # np=1, MPI=yes
  #---------------#
  echo "np=1, MPI=yes"
  clean_compile $PRO_DIR $1 yes  # dir CGNS_HDF5 MPI CGNS_ADF5
  chd $2;
  echo  save_now
  touch save_now
  if mpirun -np 1 $PRO_EXE | grep -q "chan-ts000001"; then
    echo  exit_now
    touch exit_now
    mpirun -np 1 $PRO_EXE | grep -q "# Exiting !"
    echo save_exit_now_test was successfull
  fi

  #---------------#
  # np=2, MPI=yes
  #---------------#
  nproc_in_div=$(tail -1 divide.scr)
  echo $nproc_in_div
  echo "np=2, MPI=yes"
  echo  save_now
  touch save_now
  if mpirun -np $nproc_in_div $PRO_EXE | grep -q "chan-ts000001"; then
    echo  exit_now
    touch exit_now
    mpirun -np $nproc_in_div $PRO_EXE | grep -q "# Exiting !"
    echo save_exit_now_test was successfull
  fi
}

#------------------------------------------------------------------------------#
# Processor tests
#------------------------------------------------------------------------------#
function processor_tests {

  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
  echo "!!"
  echo "!!"
  echo "!!    Running Processor Tests"
  echo "!!"
  echo "!!"
  echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  #--------------#
  # Backup tests
  #--------------#

  # Backstep/Orthogonal
  backup_test no  $TEST_DIR/Laminar/Backstep/Orthogonal
# Issue: CGNS doesn't work with OpenMPI
# backup_test yes $TEST_DIR/Laminar/Backstep/Orthogonal

  # Channel_Re_Tau_590/Uniform_Mesh
  backup_test no  $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh
# Issue: CGNS doesn't work with OpenMPI
# backup_test yes $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh

  #-----------------------------#
  # Tests for save_now/exit_now
  #-----------------------------#
  save_exit_now_test no  $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh
# Issue: CGNS doesn't work with OpenMPI
# save_exit_now_test yes $TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh

  #------#
  # Runs
  #------#
  for test_case in {1..10}
  do
    case $test_case in
      "1") case_dir=$TEST_DIR/Laminar/Cavity/Lid_Driven/Re_1000
           clean_compile $PRO_DIR no yes
      ;;
      "2") case_dir=$TEST_DIR/Laminar/Cavity/Thermally_Driven/Ra_10e6
           clean_compile $PRO_DIR no yes
      ;;
      "3") case_dir=$TEST_DIR/Rans/Channel_Re_Tau_590/Uniform_Mesh
           clean_compile_case $PRO_DIR no yes $case_dir
      ;;
      "4") case_dir=$TEST_DIR/Rans/Channel_Re_Tau_590/Stretched_Mesh
           clean_compile_case $PRO_DIR no yes $case_dir
      ;;
      "5") case_dir=$TEST_DIR/Rans/Channel_Re_Tau_590/Rsm
           clean_compile_case $PRO_DIR no yes $case_dir
      ;;
      "6") case_dir=$TEST_DIR/Rans/Channel_Re_Tau_590/Long_Domain
           clean_compile $PRO_DIR no yes
      ;;
      "7") case_dir=$TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000
           clean_compile_case $PRO_DIR no yes $case_dir
      ;;
      "8") case_dir=$TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Uniform_Mesh
           clean_compile_case $PRO_DIR no yes $case_dir
      ;;
      "9") case_dir=$TEST_DIR/Hybrid_Les_Rans/Channel_Re_Tau_2000/Stretched_Mesh
           clean_compile_case $PRO_DIR no yes $case_dir
      ;;
      "10") case_dir=$TEST_DIR/Les/Pipe_Re_Tau_180
            clean_compile_case $PRO_DIR no yes $case_dir
      ;;
    esac

    chd $case_dir
    $DIV_EXE < divide.scr
    nproc_in_div=$(tail -1 divide.scr)
    mpirun -np $nproc_in_div $PRO_EXE > out_test
  done
}

#------------------------------------------------------------------------------#
# Actual script
#------------------------------------------------------------------------------#
make_links
generator_tests
convert_tests
divide_tests
processor_tests
