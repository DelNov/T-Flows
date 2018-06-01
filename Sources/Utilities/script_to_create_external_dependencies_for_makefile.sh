#!/bin/bash

# it is an useful script in debug perposus to automatically build and run most 
# cases in T-Flows

# put your compilers here
FCOMP="gnu"; # or ifort/gfortran/mpif90/mpifort/mpiifort
DEBUG="yes"; # run tests in debug mode

# folder structure
TEST_DIR=$PWD                      # dir with tests
GENE_DIR=$PWD/../Generate  # Generate src folder
CONV_DIR=$PWD/../Convert   # Convert  src folder
DIVI_DIR=$PWD/../Divide    # Divide   src folder
PROC_DIR=$PWD/../Process   # Process  src folder
SHAR_DIR=$PWD/../Shared   # Process  src folder
BINA_DIR=$PWD/../../Binaries/         # binaries folder

# executables
GENE_EXE=$BINA_DIR/Generate        # Generate ex
CONV_EXE=$BINA_DIR/Convert         # Convert  ex
DIVI_EXE=$BINA_DIR/Divide          # Divide   ex
PROC_EXE=$BINA_DIR/Process         # Process  ex

cd $PROC_DIR; cp /dev/null  tmp

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
# produces correct module structure
#------------------------------------------------------------------------------#
function module_list() {
  # $1 - dir

  cd $1
  # find file |
  # with _Mod |
  # remove ^./ |
  # not ^. |
  # turn _Mod.f90$ to _Mod |
  # turn /abc/ to #/ |
  # turn /#abc.f90$ to #/.f90 |
  # turn /abc.f90$ to /*.f90 |
  # turn # to /* |
  # turn Mod.f90 to Mod.o |
  # turn _Mod$ to _Mod.f90 |
  # sort reverse unique
  local result=$(find . -type f -print \
               | grep -i '_Mod' \
               | sed  -e 's%^\.\/%%' \
               | grep -v '^\.' \
               | sed  -e 's%_Mod.f90$%_Mod%g' \
               | sed  -e 's%/.*/%#/%g' \
               | sed  -e 's%\#/.*f90$%#/*.f90%g' \
               | sed  -e 's%/.*.f90$%/*.f90%g' \
               | sed  -e 's%\#%/*%g' \
               | sed  -e 's%_Mod.f90%_Mod.o%g' \
               | sed -e 's%_Mod$%_Mod.o%g' \
               | sort -ur)
  echo "$result"
}

#------------------------------------------------------------------------------#
# produces correct module structure
#------------------------------------------------------------------------------#
function search_string_in_list() {
  # $1 - string
  # $2 - line_of_list
  # $3 - text to append in front

  while read -r line_of_list; do # read $2 line by line
    if [[ $line_of_list == *"$1"* ]]; then # string 1 contain string 2
      if [[ $line_of_list == *"\*"* ]]; then
        echo "$line_of_list" \\
      else
        echo "$3""$line_of_list" \\
      fi
    fi
  done <<< "$2"
}
#------------------------------------------------------------------------------#
# make file explicit dependencies constructor
#------------------------------------------------------------------------------#
function make_file_explicit_dependencies_constructor() {

  #--------------------
  # search for Modules
  #--------------------

  proc_mods=$(module_list $PROC_DIR)
  echo "$proc_mods"

  shared_mods=$(module_list $SHAR_DIR)
  echo "$shared_mods"
#  sleep 10

  cd $PROC_DIR
  for f90_file in $(find . -print | grep -i .f90); do
  # with _Mod and use |
  # not ^! |
  # remove ,* |
  # usite delimiter
  # print 3 column
    echo '#-- '$f90_file >> tmp

    deps=$(grep -ie "use .*Mod" $f90_file \
         | grep -v "^!" \
         | sed 's/\,.*$//' \
         | tr -s ' ' \
         | cut -d' ' -f3)

    if [ ! -z "$deps" ]; then #if non-empty

      echo "\$(DIR_OBJECT)/"$(basename -- "${f90_file%.*}").o: \\>> tmp

      while read -r line; do # line of current.f90

        #--------------------
        res=$(search_string_in_list "$line" "$proc_mods" "\$(DIR_OBJECT)/")
        if [ ! -z "$res" ]; then echo "$res" >> tmp; fi

        #--------------------
        res=$(search_string_in_list "$line" "$shared_mods" "\$(DIR_SHARED)/")
        if [ ! -z "$res" ]; then echo "$res" >> tmp; fi

      done <<< "$deps"
      sed -i '$ s/.$//' tmp
      echo '' >> tmp

     fi
  done #for
}

make_file_explicit_dependencies_constructor
mv tmp makefile_explicit_dependencies


#    res=''
#    if [ ! -z "$deps" ]; then #if non-empty
#
#      res=$(echo "\$(DIR_OBJECT)/"$(basename -- "${f90_file%.*}").o: \\)
#
#      while read -r line; do # line of current.f90
#
#        #--------------------
#        res=$res$(search_str_in_list "$line" "$proc_mods" "\$(DIR_OBJECT)/")
#        if [ ! -z "$res" ]; then echo "$res"; fi
#
#        #--------------------
#        res=$res$(search_str_in_list "$line" "$shared_mods" "\$(DIR_SHARED)/")
#        if [ ! -z "$res" ]; then echo "$res"; fi
#
#      done <<< "$deps"
#      echo "$res"
#      echo '----------------'
#
#     fi
#  done #for