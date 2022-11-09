#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
LIGHT_CYAN='\033[1;36m'
LIGHT_PURPLE='\033[1;35m'
NC='\033[0m' # No Color

tabs 60

#------------------------------------------------------------------------------#
# Print_usage
#------------------------------------------------------------------------------#
print_usage() {
  echo "#======================================================================"
  echo "# Utility for extraction of module hierarchy/dependencies"
  echo "#----------------------------------------------------------------------"
  echo "# Proper usage: "
  echo "#"
  echo "# ./Utilities/seek_procedure_hierarcy.sh <Target_Mod> [-e <Exclude_Dir>]"
  echo "#"
  echo "# where Target_Mod is the module name for which you want to perform"
  echo "# the analysis, such as: Grid_Mod, Convert_Mod, Generate_Mod, hence"
  echo "# case sensitive, with the _Mod suffix, without the .f90 extension."
  echo "#"
  echo "# In cases where the same module name is used in more than one direc-"
  echo "# tory, you can use the second (always -e) and the third argument to "
  echo "# exclude some directories from the search.  At the time of writing "
  echo "# this, only Point_Mod is defined in Generate and in Process."
  echo "#"
  echo "# NOTE: The script is supposed to be executed from: T-Flows/Sources!"
  echo "#----------------------------------------------------------------------"
}

#------------------------------------------------------------------------------#
# Browse through all directories looking for module dependencies
#------------------------------------------------------------------------------#
extract_call_graph() {

  #-----------------------
  #   Handle parameters
  #-----------------------
  echo "ENTERED EXTRACT WITH PARAMETERS: " $1 $2

  # First parameter is the procedure name you seek
  procedure_name_you_seek="$1"
  procedure_file_you_seek="$procedure_name_you_seek"".f90"

  # Second parameter is the module in which the procedure resides
  module_in_which_you_seek="$2"

  local next_level=`expr $next_level + 1`
  local this_level=`expr $next_level - 1`

  #----------------------------------------------#
  #   Get the full path of the module you seek   #
  #----------------------------------------------#
  if [ $module_in_which_you_seek ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep $module_in_which_you_seek)
  else
    local full_path_you_seek=$(find . -name $procedure_file_you_seek)
  fi

  # If the definition of this procedure is found, carry on
  echo $full_path_you_seek
  if [ $full_path_you_seek ]; then
    echo "Found something!"

    # This command counts number of occurrences of modules name in the result of
    # command find. If it is more than one, the same file is in more directories
    n=$(echo "$full_path_you_seek" | tr " " "\n" | grep -c "$procedure_name_you_seek")
    if [ $n -gt 1 ]; then
      echo "Ambiguity: procedure "$procedure_name_you_seek" found in more than one directory, here is the list:"
      for path in ${full_path_you_seek[*]}; do
        echo $path
      done
      echo "Exclude all but one directory with the command line argument -e <directory>"
      exit
    fi

    #-----------------------------------------------------
    #   Storing results of the grep command in an array
    #-----------------------------------------------------
    local called_procedures=($(grep '\ \ call' $full_path_you_seek | awk '{print $2$3$4}' | tr -d ,))
    local called_modules=$called_procedures   # just to declare

    # At this point, $called procedures has a form like: Profiler%Start('Main')
    # From this mess, extract the module name and the procedure element wise
    for proc in "${!called_procedures[@]}"; do
      if [[ ${called_procedures[proc]} == *"%"* ]]; then
#       echo "MODERN PROCEDURE " ${called_procedures[$proc]}
        called_modules[$proc]=$(cut -d %  -f 1 <<< ${called_procedures[$proc]})
        called_modules[$proc]=$(cut -d \( -f 1 <<< ${called_modules[$proc]})
        called_procedures[$proc]=$(cut -d %  -f 2 <<< ${called_procedures[$proc]})
        called_procedures[$proc]=$(cut -d \( -f 1 <<< ${called_procedures[$proc]})
      else
#       echo "OLD FASHIONED OR GLOBAL FUNCITON " ${called_procedures[$proc]}
        called_modules[$proc]=""  # empty string (maybe Shared?)
        called_procedures[$proc]=$(cut -d \( -f 1 <<< ${called_procedures[$proc]})
      fi
    done

    #------------------------------------------------------------------
    #   Print out the name of the module you are currently analysing
    #------------------------------------------------------------------
    echo "-----------------------------------------------------------------------------------------------------------------------"
    if [ ! -z "$called_procedures" ]; then
      echo -e ${YELLOW}"${indent}"+ $procedure_name_you_seek "("$this_level")"${NC}
    else
      echo -e ${GREEN}"${indent}"⨯ $procedure_name_you_seek "("$this_level")"${NC}
    fi

    # Increase indend for the next level by appending 6 spaces to it
    local indent="${indent}"'      '

    #-------------------------------------------------------------
    #   If the list of called procedures in not empty, carry on
    #-------------------------------------------------------------
    if [ ! -z "$called_procedures" ]; then

      # Print the procedures which are called from here
      for proc in "${!called_procedures[@]}"; do
        if [ ${called_modules[proc]} ]; then
#         echo ${called_procedures[proc]} " @ " ${called_modules[proc]}
          echo -e "${indent}"• ${called_procedures[proc]}" \t (from: "${called_modules[proc]}")"
        else
          echo -e "${indent}"• ${LIGHT_CYAN}${called_procedures[proc]}${NC}" \t (global or external)"
#         echo ${called_procedures[proc]}
        fi
      done

      exit

      # Print the modules you have found here
      for procedure in ${called_procedures[*]}; do
        echo "${indent}"• $procedure
      done

      # Print the files you will want to seek next
      for procedure in ${called_procedures[*]}; do

        echo "ABOUT TO MAKE RECURSIVE CALL FOR " $procedure " " $called_module

        # Avoid standard T-Flows modules: extension _Mod
        if [[ ! "$procedure" == *"_Mod"* ]]; then
          extract_call_graph "$procedure" "$called_module"
        fi
      done
    fi
  else
    echo "Found nothing"
  fi
}

#---------------------------------------------------------------------
#  One command line argument is sent - must be the procedure name
#  Use the names without extension - say Grid_Mod, Convert_Mod ...
#---------------------------------------------------------------------
if [ $1 ]; then
  extract_call_graph $1

#-----------------------------------------------------------------------
#  Wrong number of command line argument is sent - describe the usage
#-----------------------------------------------------------------------
else
  print_usage
fi

