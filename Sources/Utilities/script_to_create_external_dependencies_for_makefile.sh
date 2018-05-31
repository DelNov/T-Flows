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
# make file explicit dependencies constructor
#------------------------------------------------------------------------------#
function make_file_explicit_dependencies_constructor {

  # find directories
  cd $PROC_DIR  in Shared
  cd ../Shared
  #shared_mods=$(find . -type d -print | sed -e 's%\.\/%%' | grep -v '^\.' | cut -d'/' -f1 | sort -u)
  

  # find directories in Shared
  cd $PROC_DIR
  proc_mods=$(find . -type f -print | grep -iI "_Mod" | sed -e 's%^\.\/%%' |grep -v '^\.' |  sed -e 's%_Mod/.*%_Mod/*%g' | sed -e 's%.f90$%.o%' | sort -ur)



      while read -r line; do
        if [ -d "$deps"/ ]; then
          echo ---------------------
          echo "$line"/
          echo ---------------------
        else
          echo "$line"
        fi
      done <<< "$proc_mods"




  for f90_file in $(find . -print | grep -i .f90); do
    # find "use " and "Mod", but not with srarting with "!", remove everything after ",", cut 3 column

    deps=$(grep -Iie "use .*Mod" $f90_file | grep -v "^!" | sed '/,/,+1 d' | tr -s ' ' | cut -d' ' -f3)

#    # if non-empty
#    if [ ! -z "$deps" ]; then
#      echo "\$(DIR_OBJECT)/"$(basename -- "${f90_file%.*}").o: \\
#      while read -r line; do
#        #if [ "${line,,}" == "const_mod" ] || \
#        #   [ "${line,,}" == "cgns_mod" ]  || \
#        #   [ "${line,,}" == "comm_mod" ]  || \
#        #   [ "${line,,}" == "work_mod" ]  || \
#        #   [ "${line,,}" == "grid_mod" ]; then
#        #  echo "\$(DIR_SHARED)/""$line".o
#        #else
#        #  echo "\$(DIR_OBJECT)/""$line".o \\
#        #fi
#
#        if [ -d "$deps"/ ]; then
#          echo ---------------------
#          echo "$line"/
#          echo ---------------------
#        else
#          echo "$line"
#        fi
#
#
#      done <<< "$deps"
#      echo ""
#      
#    fi
#
#    #while read line; do
#    #done < <
#  done

  # (DIR_SHARED)
  #for f90_file in ../Shared/*.f90; do
  #  #statements
  #  #echo "\$(DIR_OBJECT)/""${f90_file%.*}".o: \\
  #  echo "${f90_file%.*}".o:
  #  grep -Ii "use " $f90_file | grep -Ii Mod || true
  #done
}