#!/bin/bash

# it is an useful script to create "makefile_explicit_dependencies"
# reason why you may want to use this:
# when you edit some file in Sources/, all dependent functions/subroutines must
# be updated to link program later.
# program make by default does not track dependencies.
# Therefore, you need to type "make clean; make ..." everytime you make changes
# It is allowed to specify such dependencies yourself.
# This script does it for you.

# folder structure
CURR_DIR=$PWD                      # Generate/, Convert/, Divide/, Process_?pu/
SHAR_DIR=$PWD/../Shared            # Process_?pu  src folder

# tmp file name and location
tmp_file=$CURR_DIR/tmp_makefile_explicit_dependencies

#--------------------------#
#   Read above this line   #
#--------------------------#

#set -e # exit when any command fails
#set -v #Prints shell input lines as they are read.
#set -x #Print command traces before executing command.
#set +x #disable previous line

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'if [ $? -ne 0 ]; then echo "\"${last_command}\" command failed with exit code $?."; fi' EXIT

#---------------------------------------#
#   Produces correct module structure   #
#---------------------------------------#
function module_list() {
  # $1 - dir

  cd $1
  # find file including soft links
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
  local result=$(find -L . -type f -print \
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
#---------------------------------------#
#   Produces correct module structure   #
#---------------------------------------#
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
#-------------------------------------------------#
#   Make file explicit dependencies constructor   #
#-------------------------------------------------#
function make_file_explicit_dependencies() {
  # $1 - folder to create makefile list (Generate, Divide, Process_?pu, ..)

  # create empty file or remove all content
  cd $1; cp /dev/null $tmp_file

  #------------------------#
  #   Search for modules   #
  #------------------------#

  proc_mods=$(module_list $1)
  #echo -e "$proc_mods"

  shared_mods=$(module_list $SHAR_DIR)
  #echo -e "$shared_mods"

  for folder in "$1" "$SHAR_DIR"; do
    cd "$folder"

    for f90_file in $(find . -print | grep -i .f90); do

      base_name=$(basename -- "${f90_file%.*}") #f90_file name without extension

      #--------------------------------#
      #   Determine deps of f90 file   #
      #--------------------------------#

      dependencies='' # all deps for $f90_file are store in this var

      #-----------------------------------------------------------------------#
      #   Dependencies of modules with "mod_name/subroutines.f90" structure   #
      #-----------------------------------------------------------------------#
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
      if [ ! -z "$dependencies" ]; then #if non-empty

        echo '#---------------------------------------------------' >> $tmp_file
        echo '# Dependencies for: '$f90_file >> $tmp_file
        echo '#---------------------------------------------------' >> $tmp_file

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
#------------------------------------#
#   Check T-Flows folder structure   #
#------------------------------------#
function check_if_lauched_in_correct_folder() {
  stop_execution=false

  case "${CURR_DIR##*/}" in
    "Process_Cpu"|"Process_Gpu"|"Divide"|"Convert"|"Generate")
    parent_dir="$(dirname "$CURR_DIR")"
    if [ ! "${parent_dir##*/}" == "Sources" ]; then
      echo "$parent_dir"
      stop_execution=true
    fi
    ;;
    *)
      stop_execution=true
    ;;
  esac

  if [ $stop_execution == true ]; then
    echo This script can only be launched from Generate/, Convert/,
    echo Divide/, Process_Cpu/ and Process_Gpu/ folders of T-Flows
    echo Launch it this way: bash ../Utilities/create_makefile_dependencies.sh
    exit 3
  fi
}
#-------------------#
#   Actual script   #
#-------------------#
check_if_lauched_in_correct_folder
echo creating makefile_explicit_dependencies in $CURR_DIR
make_file_explicit_dependencies $CURR_DIR
echo done
