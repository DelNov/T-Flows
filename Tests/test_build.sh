#!/bin/bash

# it is an useful script to automatically build and run most cases in T-Flows.

# put your compilers here
FCOMP="gnu"; # or ifort/gfortran/mpif90/mpifort/mpiifort
DEBUG="yes"; # run tests in debug mode

# folder structure
TEST_DIR=$PWD                      # dir with tests
GENE_DIR=$PWD/../Sources/Generate  # Generate src folder
CONV_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIVI_DIR=$PWD/../Sources/Divide    # Divide   src folder
PROC_DIR=$PWD/../Sources/Process   # Process  src folder
BINA_DIR=$PWD/../Binaries/         # binaries folder

# executables
GENE_EXE=$BINA_DIR/Generate        # Generate ex
CONV_EXE=$BINA_DIR/Convert         # Convert  ex
DIVI_EXE=$BINA_DIR/Divide          # Divide   ex
PROC_EXE=$BINA_DIR/Process         # Process  ex

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
# make links
#------------------------------------------------------------------------------#
function make_links {
  ln -rsf $BINA_DIR/* $TEST_DIR/Laminar/Backstep_Orthogonal/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Laminar/Backstep_Nonorthogonal/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Backstep_Re_26000_Rsm/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Backstep_Re_28000/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Channel_Re_Tau_590/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000/;
  ln -rsf $BINA_DIR/* $TEST_DIR/Rans/Fuel_Bundle/;
}
#------------------------------------------------------------------------------#
# generator tests
#------------------------------------------------------------------------------#
function generator_tests {
  #-- seq, no cgns
  cd $GENE_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr

  #-- seq, cgns(adf5)
  cd $GENE_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_ADF5=yes
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr

  #-- seq, cgns(hdf5)
  cd $GENE_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr
}
#------------------------------------------------------------------------------#
# convert tests
#------------------------------------------------------------------------------#
function convert_tests {
  # unpacking geometry
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                 tar -zxvf    chan.tar.gz
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;  tar -zxvf     jet.tar.gz
  cd $TEST_DIR/Rans/Fuel_Bundle;                        tar -zxvf subflow.tar.gz

  #-- seq, no cgns
  cd $CONV_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr

  #-- seq, cgns (adf5)
  cd $CONV_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_ADF5=yes
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr

  #-- seq, cgns (hdf5)
  cd $CONV_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=yes
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr
}
#------------------------------------------------------------------------------#
# Divide tests
#------------------------------------------------------------------------------#
function divide_tests {
  #-- seq
  cd $DIVI_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;               $DIVI_EXE < divide.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;            $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;                $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                    $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;     $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                   $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;    $DIVI_EXE < divide.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                          $DIVI_EXE < divide.scr
}
#------------------------------------------------------------------------------#
# backup test
#------------------------------------------------------------------------------#
function backup_test {
  # $1 = CGNS_HDF5 = yes
  # $2 = test_dir

  #----------------------------------#
  # np=1, MPI=no, backup=no
  #----------------------------------#
  echo "np=1, MPI=no, backup=no"
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$1 MPI=no
  cd $2;
  cp control control.backup
  line="$(grep -ni "NUMBER_OF_TIME_STEPS"       control | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=3}1'     control > .control.tmp
  line="$(grep -ni "BACKUP_SAVE_INTERVAL"  .control.tmp | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=1}1'     .control.tmp > control
  $PROC_EXE
  #----------------------------------#
  # np=1, MPI=no, backup=(from np=1)
  #----------------------------------#
  echo "np=1, MPI=no, backup=(from np=1)"
  sed -e  's/#LOAD_BACKUP_NAME/LOAD_BACKUP_NAME/g' .control.tmp > control
  $PROC_EXE
  #----------------------------------#
  # np=1, MPI=yes, backup=(from np=1)
  #----------------------------------#
  echo "np=1, MPI=yes, backup=(from np=1)"
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$1 MPI=yes
  cd $2;
  mpirun -np 1 $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=1)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=1)"
  mpirun -np 2 $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=2)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=2)"
  line="$(grep -ni "BACKUP_SAVE_INTERVAL"  .control.tmp | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=2}1'    .control.tmp > control
  mpirun -np 2 $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=no
  #----------------------------------#
  echo "np=2, MPI=yes, backup=no"
  sed -e  's/LOAD_BACKUP_NAME/#LOAD_BACKUP_NAME/g' .control.tmp > control
  mpirun -np 2 $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=2)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=2)"
  sed -e  's/#LOAD_BACKUP_NAME/LOAD_BACKUP_NAME/g' .control.tmp > control
  mpirun -np 2 $PROC_EXE
  #----------------------------------#
  # np=1, MPI=no, backup=(from np=2)
  #----------------------------------#
  echo "np=1, MPI=yes, backup=(from np=2)"
  cd $PROC_DIR; make clean
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$1 MPI=no
  cd $2;
  $PROC_EXE

  cp control.backup control
  rm control.backup .control.tmp
}
#------------------------------------------------------------------------------#
# backup tests
#------------------------------------------------------------------------------#
function processor_tests {
  #-- Backstep_Orthogonal
  backup_test no  $TEST_DIR/Laminar/Backstep_Orthogonal
  backup_test yes $TEST_DIR/Laminar/Backstep_Orthogonal
  #-- Backstep_Orthogonal
  backup_test no  $TEST_DIR/Rans/Channel_Re_Tau_590
  backup_test yes $TEST_DIR/Rans/Channel_Re_Tau_590
}
#------------------------------------------------------------------------------#
# actual script
#------------------------------------------------------------------------------#
make_links
generator_tests
convert_tests
divide_tests
processor_tests
