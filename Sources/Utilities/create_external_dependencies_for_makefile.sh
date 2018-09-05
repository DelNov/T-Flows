#!/bin/bash

# it is an useful script to create "makefile_explicit_dependencies"
# reason why you may wnat to use this:
# when you edit some file in Sources/, all dependent functions/subroutines must
# be updated to link program later.
# program make by default does not track dependencies.
# Therefore, you need to type "make clean; make ..." everytime you make changes
# It is allowed to specify such dependencies yourself.
# This script does it for you.

# folder structure
TEST_DIR=$PWD                      # Sources/Utilities
GENE_DIR=$PWD/../Generate          # Generate src folder
CONV_DIR=$PWD/../Convert           # Convert  src folder
DIVI_DIR=$PWD/../Divide            # Divide   src folder
PROC_DIR=$PWD/../Process           # Process  src folder
SHAR_DIR=$PWD/../Shared            # Process  src folder
BINA_DIR=$PWD/../../Binaries/      # binaries folder

# tmp file name and location
tmp_file=$BINA_DIR/tmp

#-------------------------------------------------#
#---------   READ ABOVE UP TO THIS ROW   ---------#
#-------------------------------------------------#

set -e # exit when any command fails
#set -v #Prints shell input lines as they are read.
#set -x #Print command traces before executing command.
#set +x #disable previous line

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
  # find file
  # with _Mod
  # remove ^./
  # not ^.
  # turn _Mod.f90$ to _Mod.o
  # turn _Mod/ to #
  # turn #abc/ to ?
  # turn ?abc.f90$ to ?*.f90
  # remove ^abc/
  # turn #abc.f90$ to #*.f90$
  # turn # to _Mod
  # turn ? to _Mod/*/
  # sort reverse unique
  # delete empty line
  local result=$(find . -type f -print \
               | grep -i '_Mod' \
               | sed  -e 's%^\.\/%%' \
               | grep -v '^\.' \
               | sed  -e 's%_Mod.f90$%_Mod.o%g' \
               | sed  -e 's%_Mod/%#%' \
               | sed  -e 's%#.*/%?%g' \
               | sed  -e 's%\?.*f90$%?*.f90%g' \
               | sed  -e 's%^.*/%%g' \
               | sed  -e 's%\#.*.f90$%#*.f90%g' \
               | sed  -e 's%\#%_Mod/%g' \
               | sed  -e 's%\?%_Mod/*/%g'  \
               | sort -ur \
               | sed '/^\s*$/d')


  echo "$result"
}
#------------------------------------------------------------------------------#
# produces correct module structure
#------------------------------------------------------------------------------#
function search_string_in_list() {
  # $1 - string
  # $2 - list
  # $3 - f90_file name without extension
  # $4 - text to append in front

  # to deal with "Comm_Mod_Par" and "Cgns_Mod_Par"
  str=$(echo "$1" | sed 's%Mod_.*$%Mod%g')
  while read -r line_of_list; do # read $2 line by line
    if [[ $line_of_list == *"$str"* ]]; then # if string 1 contains string 2
      if [[ $line_of_list == *"/*"* ]]; then
        echo "$4""$line_of_list" \\
      elif [ "${line_of_list:(-5)}" == "Mod.o" ]; then
        if [ "$str" != "$3" ]; then
          echo "\$(DIR_OBJECT)/""$line_of_list" \\
        fi
      else
        echo "$4""$line_of_list" \\
      fi
    fi
  done <<< "$2"
}
#------------------------------------------------------------------------------#
# make file explicit dependencies constructor
#------------------------------------------------------------------------------#
function make_file_explicit_dependencies() {
  # $1 - folder to create makefile list (Generate, Divide, Process, ..)

  # create empty file or remove all content
  cd $1; cp /dev/null $tmp_file

  #--------------------
  # search for Modules
  #--------------------

  proc_mods=$(module_list $1)
  #echo -e "$proc_mods"

  shared_mods=$(module_list $SHAR_DIR)
  #echo -e "$shared_mods"

  for folder in "$1" "$SHAR_DIR"; do
    cd "$folder"

    for f90_file in $(find . -print | grep -i .f90); do

      base_name=$(basename -- "${f90_file%.*}") #f90_file name without extension

      #---------------------------
      # determine deps of f90 file
      #---------------------------

      dependencies='' # all deps for $f90_file are store in this var

      #---------------------------------------------------------------------
      # dependencies of modules with "module_name/subroutines.f90" structure
      #---------------------------------------------------------------------
      mod_included_this=$(grep -ie "include " $f90_file \
           | cut -d"!" -f1 \
           | sed "s%'%%g" \
           | sed 's%"%%g' \
           | grep -viI '.h"$\|.h$' \
           | tr -s ' ' \
           | cut -d' ' -f3- \
           | sort -ur \
           | sed '/^\s*$/d')
      # with include
      # ignore commented
      # remove '
      # remove "
      # not .h" or .h'
      # unite delimiter
      # print from 3rd column
      # unique
      # delete empty line


      if [ ! -z "$mod_included_this" ]; then
        while read -r mod_included_this_line; do # line of mod_included_this
          dependencies=$(echo -e "$dependencies\n$(\
              grep -ie "use .*Mod" $mod_included_this_line \
            | cut -d"!" -f1 \
            | sed 's/\,.*$//' \
            | tr -s ' ' \
            | cut -d' ' -f3 \
            | sort -ur \
            | sed '/^\s*$/d')")
          # with use and Mod
          # ignore commented
          # remove anything after ,
          # unite delimiter
          # print 3rd column
          # unique
          # delete empty line

        done <<< "$mod_included_this"
        #echo "$f90_file"
        #echo "$f90_file" | sed 's%.f90%%g' | sed 's%./%%g'
        dependencies=$(echo -e "$dependencies\n$(echo "$f90_file" | sed 's%.f90%%g' | sed 's%./%%g')")
      fi
      #------------------------------------------
      # dependencies of functions and subroutines
      #------------------------------------------
      file_included_this=$(grep -ie "use .*Mod" $f90_file \
        | cut -d"!" -f1 \
        | sed 's/\,.*$//' \
        | sed 's%\/.*$%%' \
        | sed "s%'%%g" \
        | tr -s ' ' \
        | cut -d' ' -f3 \
        | sort -ur \
        | sed '/^\s*$/d')
        # with use and Mod
        # ignore commented
        # remove anything after ,
        # remove anything after /
        # unite delimiter
        # print 3rd column
        # unique
        # delete empty line


      if [ ! -z "$file_included_this" ]; then #if non-empty
        if [ ! -z "$dependencies" ]; then #if non-empty
          dependencies=$(echo -e "$dependencies\n$file_included_this")
        else
          dependencies="$file_included_this"
        fi
      fi

      if [ ! -z "$dependencies" ]; then
        dependencies=$(echo "$dependencies" | sort -ur \
        | sed '/^\s*$/d')
        #echo -e "$(basename -- "${f90_file%.*}")"" depends on\n""$dependencies"
      fi

      #-------------------------------------------------
      # search for deps in list proc_mods and shar_mods
      # if found : add them to $tmp_file
      #-------------------------------------------------
      echo '#-- '$f90_file >> $tmp_file

      if [ ! -z "$dependencies" ]; then #if non-empty

        echo "\$(DIR_OBJECT)/""$base_name".o:\\ >> $tmp_file

        while read -r line_of_deps_list; do # line of current.f90

          # $1 deps
          dep_list=$(search_string_in_list \
            "$line_of_deps_list" "$proc_mods" "$base_name" "")
          if [ ! -z "$dep_list" ]; then
            echo "$dep_list" >> $tmp_file
          fi

          # Shared deps
          dep_list=$(search_string_in_list \
            "$line_of_deps_list" "$shared_mods" "$base_name" "\$(DIR_SHARED)/")
          if [ ! -z "$dep_list" ]; then
            echo "$dep_list" >> $tmp_file
          fi

        done <<< "$dependencies"
        sed -i '$ s/.$//' $tmp_file
        echo '' >> $tmp_file

       fi
    done #for f90_file

  done #for folder

  cd $1; mv $tmp_file makefile_explicit_dependencies
}
#------------------------------------------------------------------------------#
# actual script
#------------------------------------------------------------------------------#
echo creating makefile_explicit_dependencies in Generete
make_file_explicit_dependencies $GENE_DIR
echo creating makefile_explicit_dependencies in Convert
make_file_explicit_dependencies $CONV_DIR
echo creating makefile_explicit_dependencies in Divide
make_file_explicit_dependencies $DIVI_DIR
echo creating makefile_explicit_dependencies in Processor
make_file_explicit_dependencies $PROC_DIR
echo done