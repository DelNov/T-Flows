#!/bin/bash

# it is an useful script in debug perposus to automatically build and run most 
# cases in T-Flows

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
# clean_compile
#------------------------------------------------------------------------------#
function clean_compile {
  # $1 = dir 
  # $2 = CGNS_HDF5 = yes/no
  # $3 = MPI = yes/no
  # $4 = CGNS_ADF5 = yes/no

  cd $1
  make clean
  echo "make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 CGNS_ADF5=$4"
  make FORTRAN=$FCOMP DEBUG=$DEBUG CGNS_HDF5=$2 MPI=$3 CGNS_ADF5=$4
}
#------------------------------------------------------------------------------#
# make links
#------------------------------------------------------------------------------#
function make_links {
  cd $TEST_DIR/Laminar/Backstep_Orthogonal/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Laminar/Backstep_Orthogonal/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Rans/Backstep_Re_28000/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Rans/Channel_Re_Tau_590/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
  cd $TEST_DIR/Rans/Fuel_Bundle/;
  ln -rsf $GENE_EXE $CONV_EXE $DIVI_EXE $PROC_EXE .
}
#------------------------------------------------------------------------------#
# generator tests
#------------------------------------------------------------------------------#
function generator_tests {
  #-- seq, no cgns
  clean_compile $GENE_DIR no no yes # dir CGNS_HDF5 MPI CGNS_ADF5

  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr

  #-- seq, cgns(adf5)
  clean_compile $GENE_DIR no no yes # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $TEST_DIR/Laminar/Backstep_Orthogonal;             $GENE_EXE < generate.scr
  cd $TEST_DIR/Laminar/Backstep_Nonorthogonal;          $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_26000_Rsm;              $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Backstep_Re_28000;                  $GENE_EXE < generate.scr
  cd $TEST_DIR/Rans/Channel_Re_Tau_590_Wall_Function;   $GENE_EXE < generate.scr

  #-- seq, cgns(hdf5)
  clean_compile $GENE_DIR yes no no # dir CGNS_HDF5 MPI CGNS_ADF5
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
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                tar -zxvf    chan.neu.tgz
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000; tar -zxvf     jet.neu.tgz
  cd $TEST_DIR/Rans/Fuel_Bundle;                       tar -zxvf subflow.neu.tgz

  #-- seq, no cgns
  clean_compile $CONV_DIR no no yes # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr

  #-- seq, cgns (adf5)
  clean_compile $CONV_DIR no no yes # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr

  #-- seq, cgns (hdf5)
  clean_compile $CONV_DIR yes no no # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $TEST_DIR/Rans/Channel_Re_Tau_590;                  $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Impinging_Jet_2d_Distant_Re_23000;   $CONV_EXE < convert.scr
  cd $TEST_DIR/Rans/Fuel_Bundle;                         $CONV_EXE < convert.scr
}
#------------------------------------------------------------------------------#
# Divide tests
#------------------------------------------------------------------------------#
function divide_tests {
  #-- seq
  clean_compile $DIVI_DIR no no no # dir CGNS_HDF5 MPI CGNS_ADF5
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
# processor: backup test
#------------------------------------------------------------------------------#
function backup_test {
  # $1 = CGNS_HDF5 = yes
  # $2 = test_dir

  #----------------------------------#
  # np=1, MPI=no, backup=no
  #----------------------------------#
  echo "np=1, MPI=no, backup=no"
  clean_compile $PROC_DIR $1 no no # dir CGNS_HDF5 MPI CGNS_ADF5
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
  clean_compile $PROC_DIR $1 yes no # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $2;
  mpirun -np 1 $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=1)
  #----------------------------------#
  nproc_in_div=$(tail -n2  divide.scr | head -n1)
  echo "np=2, MPI=yes, backup=(from np=1)"
  mpirun -np $nproc_in_div $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=2)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=2)"
  line="$(grep -ni "BACKUP_SAVE_INTERVAL"  .control.tmp | cut -d: -f1)"
  awk -v var="$line" 'NR==var {$2=2}1'    .control.tmp > control
  mpirun -np $nproc_in_div $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=no
  #----------------------------------#
  echo "np=2, MPI=yes, backup=no"
  sed -e  's/LOAD_BACKUP_NAME/#LOAD_BACKUP_NAME/g' .control.tmp > control
  mpirun -np $nproc_in_div $PROC_EXE
  #----------------------------------#
  # np=2, MPI=yes, backup=(from np=2)
  #----------------------------------#
  echo "np=2, MPI=yes, backup=(from np=2)"
  sed -e  's/#LOAD_BACKUP_NAME/LOAD_BACKUP_NAME/g' .control.tmp > control
  mpirun -np $nproc_in_div $PROC_EXE
  #----------------------------------#
  # np=1, MPI=no, backup=(from np=2)
  #----------------------------------#
  echo "np=1, MPI=yes, backup=(from np=2)"
  clean_compile $PROC_DIR $1 no no # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $2;
  $PROC_EXE

  cp control.backup control
  rm control.backup .control.tmp
}
#------------------------------------------------------------------------------#
# processor: save_now / exit_now test
#------------------------------------------------------------------------------#
function save_exit_now_test {
  # $1 = CGNS_HDF5 = yes
  # $2 = test_dir

  #---------------#
  # np=1, MPI=no
  #---------------#
  echo "np=1, MPI=no"
  clean_compile $PROC_DIR $1 no no # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $2;
  echo  save_now
  touch save_now
  if $PROC_EXE | grep -q "chan-ts000001"; then # exit if match is found
    echo  exit_now
    touch exit_now
    $PROC_EXE | grep -q "# Exiting !"
    echo save_exit_now_test was successfull
  fi
  #---------------#
  # np=1, MPI=yes
  #---------------#
  echo "np=1, MPI=yes"
  clean_compile $PROC_DIR $1 yes no # dir CGNS_HDF5 MPI CGNS_ADF5
  cd $2;
  echo  save_now
  touch save_now
  if mpirun -np 1 $PROC_EXE | grep -q "chan-ts000001"; then
    echo  exit_now
    touch exit_now
    mpirun -np 1 $PROC_EXE | grep -q "# Exiting !"
    echo save_exit_now_test was successfull
  fi
  #---------------#
  # np=2, MPI=yes
  #---------------#
  nproc_in_div=$(tail -n2  divide.scr | head -n1)
  echo $nproc_in_div
  echo "np=2, MPI=yes"
  echo  save_now
  touch save_now
  if mpirun -np $nproc_in_div $PROC_EXE | grep -q "chan-ts000001"; then
    echo  exit_now
    touch exit_now
    mpirun -np $nproc_in_div $PROC_EXE | grep -q "# Exiting !"
    echo save_exit_now_test was successfull
  fi
}
#------------------------------------------------------------------------------#
# processor tests
#------------------------------------------------------------------------------#
function processor_tests {
#-- backup tests
  #-- Backstep_Orthogonal
  backup_test no  $TEST_DIR/Laminar/Backstep_Orthogonal
  backup_test yes $TEST_DIR/Laminar/Backstep_Orthogonal
  #-- Channel_Re_Tau_590
  backup_test no  $TEST_DIR/Rans/Channel_Re_Tau_590
  backup_test yes $TEST_DIR/Rans/Channel_Re_Tau_590
  ##-- Fuel_Bundle
  #backup_test no  $TEST_DIR/Rans/Fuel_Bundle
  #backup_test yes $TEST_DIR/Rans/Fuel_Bundle
#-- save_now/exit_now tests
  save_exit_now_test no  $TEST_DIR/Rans/Channel_Re_Tau_590
  save_exit_now_test yes $TEST_DIR/Rans/Channel_Re_Tau_590

}
#------------------------------------------------------------------------------#
# actual script
#------------------------------------------------------------------------------#
make_links
generator_tests
convert_tests
divide_tests
processor_tests