#!/bin/bash

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
PURPLE='\033[0;35m'
LIGHT_CYAN='\033[1;36m'
LIGHT_PURPLE='\033[1;35m'
NC='\033[0m' # No Color

#------------------------------------------------------------------------------#
# Print_usage
#------------------------------------------------------------------------------#
print_usage() {
  echo "#======================================================================"
  echo "# Utility for extraction of module hierarchy/dependencies from T-Flows"
  echo "#----------------------------------------------------------------------"
  echo "# Proper usage: "
  echo "#"
  echo "# ./Utilities/extract_module_hierarchy.sh <Target_Mod> [-e <Exclude_Dir>]"
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
extract_hierarchy() {

  #-----------------------
  #   Handle parameters
  #-----------------------

  # First parameter is the module name you seek
  module_name_you_seek="$1"
  module_file_you_seek="$module_name_you_seek"".f90"

  # Second parameter is the level at which you currently are
  local next_level=`expr $next_level + 1`
  local this_level=`expr $next_level - 1`

  #----------------------------------------------#
  #   Get the full path of the module you seek   #
  #----------------------------------------------#
  if [ $exclude_dir ]; then
    local full_path_you_seek=$(find . -name $module_file_you_seek | grep -v $exclude_dir)
  else
    local full_path_you_seek=$(find . -name $module_file_you_seek)
  fi

  # This command counts number of occurrences of modules name in the result of
  # command find. If it is more than one, the same file is in more directories
  n=$(echo "$full_path_you_seek" | tr " " "\n" | grep -c "$module_name_you_seek")
  if [ $n -gt 1 ]; then
    echo "Ambiguity: module "$module_name_you_seek" found in more than one directory, here is the list:"
    for path in ${full_path_you_seek[*]}; do
      echo $path
    done
    echo "Exclude all but one directory with the command line argument -e <directory>"
    exit
  fi

  #-----------------------------------------------------
  #   Storing results of the grep command in an array
  #-----------------------------------------------------
  local used_modules=($(grep '\ \ use' $full_path_you_seek | awk '{print $2}' | tr -d ,))

  #------------------------------------------------------------------
  #   Print out the name of the module you are currently analysing
  #------------------------------------------------------------------
  if [ ! -z "$used_modules" ]; then
    if [ $this_level -lt 2 ]; then
      echo "-----------------------------------------------------------------"
      echo -e "${RED}${indent}"• $module_name_you_seek "("$this_level")${NC}"
    else
      echo -e "${indent}"• $module_name_you_seek "("$this_level")"
    fi
  else
    echo -e "${LIGHT_CYAN}${indent}"⨯ $module_name_you_seek "("$this_level")${NC}"
  fi

  # Increase indend for the next level by appending spaces to it
  local indent="${indent}"'     '

  #--------------------------------------------------------
  #   If the list of used modules in not empty, carry on
  #--------------------------------------------------------
  if [ ! -z "$used_modules" ]; then

    # Print the modules you have found here
    for mod in "${!used_modules[@]}"; do
#->      echo "${indent}"• ${used_modules[mod]}
#->    done
#->
#->    # Print the files you will want to seek next
#->    for mod in "${!used_modules[@]}"; do

      #-------------------------------------------------------#
      #                                                       #
      #   The very important recursive call to its own self   #
      #                                                       #
      #-------------------------------------------------------#
      if [[ "${used_modules[mod]}" == *"_Mod"* ]]; then   # only standard T-Flows modules
        extract_hierarchy "${used_modules[mod]}" $2 $3
      fi
    done
  fi
}

#----------------------------------------------------------------------------
#  Three command line arguments are sent - process the second and carry on
#----------------------------------------------------------------------------
if [ $3 ]; then
  if [ $2 == "-e" ]; then
    exclude_dir=""
    if [ "$3" ]; then
      exclude_dir="$3"
    fi
    extract_hierarchy $1
  else
    print_usage
  fi

#---------------------------------------------------------------------
#  One command line argument is sent - must be the module name
#  Use the names without extension - say Grid_Mod, Convert_Mod ...
#---------------------------------------------------------------------
elif [ $1 ]; then
  extract_hierarchy $1

#-----------------------------------------------------------------------
#  Wrong number of command line argument is sent - describe the usage
#-----------------------------------------------------------------------
else
  print_usage
fi

