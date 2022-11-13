#!/bin/bash

BLACK='\U001B[30m'
RED='\U001B[31m'
GREEN='\U001B[32m'
YELLOW='\U001B[33m'
BLUE='\U001B[34m'
MAGENTA='\U001B[35m'
CYAN='\U001B[36m'
WHITE='\U001B[37m'
RESET='\U001B[0m'

LIGHT_BLACK='\U001B[30;1m'
LIGHT_RED='\U001B[31;1m'
LIGHT_GREEN='\U001B[32;1m'
LIGHT_YELLOW='\U001B[33;1m'
LIGHT_BLUE='\U001B[34;1m'
LIGHT_MAGENTA='\U001B[35;1m'
LIGHT_CYAN='\U001B[36;1m'
LIGHT_WHITE='\U001B[37;1m'
RESET='\U001B[0m'

#------------------------------------------------------------------------------#
#   Global variables affecting the looks of the output
#------------------------------------------------------------------------------#

# The following four affect the width of the output
glo_indent="    "      # four characters wide
glo_separate="----"    # four characters wide, should be the same as glo_indent
glo_out_width=72       # should be multiple of indent and separator widhts

# The following lines desribe the color scheme
glo_color_mu=$LIGHT_CYAN       # module users other
glo_color_mn=$GREEN            # module not using others

# Global list of directories to be excluded from the search
glo_exclude_dir=""

#==============================================================================#
#   Print the separator line
#------------------------------------------------------------------------------#
print_separator() {

  ind=$1  # current indentation
  lev=$2  # current level

  printf "%s" "$ind"
  end=`expr $glo_out_width / ${#glo_separate} - $lev`
  for (( c=1; c<=$end; c++ ))
  do
    echo -n $glo_separate
  done
  echo ""
}

#==============================================================================#
#   Print a line
#------------------------------------------------------------------------------#
print_line() {

  ind=$1     # indentation
  color=$2   # color
  bullet=$3  # shape of the bullet
  proced=$4  # procedure
  lev=$5     # level
  module=$6  # module

  echo -e "$ind"${color}"$bullet""$proced"" ($lev)"${RESET}"\t ""$module"
}

#==============================================================================#
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
  exit
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

  if [ "$this_level" -eq 0 ]; then
    echo "#======================================================================="
    echo "# Extracting module hierarchy for "$module_name_you_seek
    echo "#"
    echo "# Modules are designated as follows:"
    echo "#"
    echo -n "# - modules using other modules:    "
    echo -e "${glo_color_mu}${indent}"• Using_Others_Mod "(level)${RESET}"
    echo -n "# - modules not using any other:    "
    echo -e "${glo_color_mn}${indent}"⨯ Not_Using_Any_Mod "(level)${RESET}"
    echo "#-----------------------------------------------------------------------"
  fi

  #---------------------------------------------------------------------
  #   Is the module you are currently analyzing on the excluded list?
  #---------------------------------------------------------------------
  if [[ ! "$glo_ignore_mod" == *"$module_name_you_seek"* ]]; then

    #----------------------------------------------
    #   Get the full path of the module you seek
    #----------------------------------------------
    if [ $glo_exclude_dir ]; then
      local full_path_you_seek=$(find . -name $module_file_you_seek | grep -v $glo_exclude_dir)
    else
      local full_path_you_seek=$(find . -name $module_file_you_seek)
    fi

    #--------------------------------
    #   If there is anything to do
    #--------------------------------
    if [ "$full_path_you_seek" ]; then

      # This command counts number of occurrences of modules name in
      # the result of command find. If it is more than one, the same
      # file is in more directories
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
      if [ "$used_modules" ]; then
        print_separator "$indent" $this_level
        print_line "$indent"                 \
                   $glo_color_mu             \
                   "• "                      \
                   $module_name_you_seek     \
                   $this_level
      else
        print_line "$indent"                 \
                   $glo_color_mn             \
                   "⨯ "                      \
                   $module_name_you_seek     \
                   $this_level
      fi

      # Increase indend for the next level by appending spaces to it
      local indent="${indent}"$glo_indent

      #--------------------------------------------------------
      #   If the list of used modules in not empty, carry on
      #--------------------------------------------------------
      if [ "$used_modules" ]; then

        # Print the modules you have found here
        for mod in "${!used_modules[@]}"; do

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

    fi  # if the paths is found (was not excluded by grep -v ...

  fi  # if this module is not ignored
}

#------------------------------------------------------------------------
#   Wrong number of command line argument is sent - describe the usage
#------------------------------------------------------------------------
if [ ! "$1" ]; then
  print_usage
fi

#-----------------------------------------------
#   Parse command line options like a pro :-)
#-----------------------------------------------

# Fetch the name and shift on
name=$1
shift

current_opt=""

while [[ $# -gt 0 ]]; do
  case $1 in
    -i)
      current_opt=$1
      glo_ignore_mod=$glo_ignore_mod" $2"
      shift # past argument
      shift # past value
      ;;    # part of the case construct
    -e)
      current_opt=$1
      glo_exclude_dir=$2
      shift # past argument
      shift # past value
      ;;    # part of the case construct
    --default)
      current_opt=$1
      default=yes
      shift # past argument
      ;;    # part of the case construct
    -*)
      echo "Unknown option $1"
      exit 1
      ;;    # part of the case construct
    *)
      if [ "$current_opt" == "-i" ]; then
        glo_ignore_mod=$glo_ignore_mod" $1"
      else
        echo "Unknown option $1"
        exit 1
      fi
      shift # past argument
      ;;    # part of the case construct
  esac
done

extract_hierarchy $name
# echo "default = ${default}"
