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
  echo "# Utility for extract call graphs from T-Flows"
  echo "#----------------------------------------------------------------------"
  echo "# Proper usage: "
  echo "#"
  echo "# ./Utilities/extract_call_graph.sh <Source> [-f <Force_Dir>]"
  echo "#"
  echo "# where Source  is the procedure name for which you want to perform"
  echo "# the analysis, such as: Main_Con, Compute_Energy, Main_Div, hence"
  echo "# case sensitive, without the .f90 extension."
  echo "#"
  echo "# In cases where the same module name is used in more than one direc-"
  echo "# tory, you can use the second (always -e) and the third argument to "
  echo "# exclude some directories from the search.  At the time of writing "
  echo "# this, Allocate_Cells is defined in Grid_Mod and Refines_Mod, but ."
  echo "# there could be more such examples."
  echo "#"
  echo "# NOTE: The script is supposed to be executed from: T-Flows/Sources!"
  echo "#----------------------------------------------------------------------"
}

glo_procedure=""
glo_module=""

#------------------------------------------------------------------------------#
# Sets $glo_procedure and $glo_module
#------------------------------------------------------------------------------#
extract_procedure_and_module() {

  full_name=$1

  if [[ "$full_name" == *"%"* ]]; then

    # The following four lines would work for:
    # Profiler%Start()
    # Grid(d)%Calculate() ...
    glo_module=$(cut -d %  -f 1 <<< $full_name)
    glo_module=$(cut -d \( -f 1 <<< $glo_module)
    glo_procedure=$(cut -d %  -f 2 <<< $full_name)
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

    # ... but get messed up for this
    # system_clock(Profiler % i_time_curr)
    if [[ "$glo_procedure" == *")"* ]]; then
      glo_procedure=$glo_module
      glo_module=""  # empty string (maybe none or Shared?)
    fi

  elif [[ "$full_name" == *"_Mod_"* ]]; then

    # This contraption is for friend functions like Comm_Mod_End
    glo_module=$(echo $full_name | awk -F '_Mod_' '{print $1}')
    glo_procedure=$full_name

  else

    # This is for global functions and external functions like Probe_1d
    glo_module=""  # empty string (maybe none or Shared?)
    glo_procedure=$(cut -d \( -f 1 <<< $full_name)

  fi
}

#------------------------------------------------------------------------------#
# Browse through all directories looking for a call graph
#------------------------------------------------------------------------------#
extract_call_graph() {

  #-----------------------
  #   Handle parameters
  #-----------------------

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
  if [ $module_in_which_you_seek ] && [ $exclude_dir ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep $module_in_which_you_seek | grep -v $exclude_dir)
  elif [ $exclude_dir ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep $exclude_dir)
  elif [ $module_in_which_you_seek ]; then
    local full_path_you_seek=$(find . -name $procedure_file_you_seek | grep $module_in_which_you_seek)
  else
    local full_path_you_seek=$(find . -name $procedure_file_you_seek)
  fi

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

  # If the definition of this procedure is found, carry on
  if [ $full_path_you_seek ]; then

    #-----------------------------------------------------
    #   Storing results of the grep command in an array
    #-----------------------------------------------------
    local called_procedures=($(grep '\ \ call' $full_path_you_seek | awk '{print $2$3$4}' | tr -d ,))
    local called_modules=$called_procedures   # just to declare

    echo "FETHED PROCEDURES:"
    for proc in "${!called_procedures[@]}"; do
      echo ${called_procedures[proc]}
    done

    # At this point, $called procedures has a form like: Profiler%Start('Main')
    # From this mess, extract the module name and the procedure element wise
    for proc in "${!called_procedures[@]}"; do
      extract_procedure_and_module ${called_procedures[proc]}  # sets $glo_procedure and $glo_module
      called_modules[$proc]=$glo_module
      called_procedures[$proc]=$glo_procedure
    done

    #------------------------------------------------------------------
    #   Print out the name of the module you are currently analysing
    #------------------------------------------------------------------
    echo "-----------------------------------------------------------------------------------------------------------------------"
    if [ ! -z "$called_procedures" ]; then
      echo -e ${YELLOW}"${indent}"+ $procedure_name_you_seek "("$this_level")"${NC}" calls: ""\t $module_in_which_you_seek"
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
          echo -e "${indent}"• ${called_procedures[proc]}" \t (from: "${called_modules[proc]}")"
        else
          echo -e "${indent}"• ${LIGHT_CYAN}${called_procedures[proc]}${NC}" \t (global or external)"
        fi

        #-------------------------------------------------------#
        #                                                       #
        #   The very important recursive call to its own self   #
        #                                                       #
        #-------------------------------------------------------#
        extract_call_graph ${called_procedures[proc]} ${called_modules[proc]}
      done
    fi
# else
#   echo "Found nothing"
  fi
}

#------------------------------------------------------------------------------
#  Three command line arguments are sent - process the second and pass them on
#------------------------------------------------------------------------------
if [ $3 ]; then
  if [ $2 == "-e" ]; then
    exclude_dir=""
    if [ "$3" ]; then
      exclude_dir="$3"
    fi
    extract_call_graph $1
  else
    print_usage
  fi

#---------------------------------------------------------------------
#  One command line argument is sent - must be the procedure name
#  Use the names without extension - say Grid_Mod, Convert_Mod ...
#---------------------------------------------------------------------
elif [ $1 ]; then
  extract_call_graph $1

#-----------------------------------------------------------------------
#  Wrong number of command line argument is sent - describe the usage
#-----------------------------------------------------------------------
else
  print_usage
fi

