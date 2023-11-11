#==============================================================================#
#   Script for extracting call graphs from T-Flows sources
#------------------------------------------------------------------------------#
#   A brief description of the scripts main components and functionalities:
#
#   * Color Definitions: The script starts by defining color variables for
#     output formatting
#
#   * Global Settings: It sets global variables for the root source directory,
#     output formatting, and color schemes for different types of procedures
#     (e.g., member callers, global callers).
#
#   * Utility Functions: Functions like print_separator, print_line, and
#     print_usage are defined for formatting and providing usage information,
#     which is consistent with your module hierarchy script.
#
#   * Main Functionality - extract_call_graph:
#     - This function recursively explores procedure calls within your code.
#     - It takes a procedure name as input and analyses its calls by reading
#       .f90 files.
#     - The script distinguishes between member callers, global callers, and
#       external procedures.
#     - It handles various cases like procedures calling others, procedures not
#       calling any others, and excluded or ignored modules.
#     - The script uses grep and awk for text processing within Fortran source
#       files.
#
#   * Command-Line Argument Parsing: The script processes command-line arguments
#     for options like expanding all calls (-a), excluding certain directories
#     (-e), and ignoring specific modules (-i).
#
#   * Procedure and Module Extraction: An additional function
#     extract_procedure_and_module is used to parse and determine the procedure
#     and module names from the call lines.
#
#   * Execution Logic: The script executes the extract_call_graph function based
#    on the input parameters and command-line arguments.
#------------------------------------------------------------------------------#

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
#   Settings and global variables affecting the looks of the output
#------------------------------------------------------------------------------#

# Find the "root" of all sources
cur=$(pwd)
src=$(echo "$cur" | awk -F'/Sources/' '{print ($2 == "" ? $0 : $1 "/Sources")}')

# The following four affect the width of the output
tabs 60
glo_indent="    "    # four characters wide
glo_separate="----"  # four characters wide, should be the same as glo_indent
glo_out_width=72     # should be multiple of indent and separator widhts

# The following lines desribe the color scheme
glo_color_mc=$LIGHT_GREEN      # member caller
glo_color_mm=$GREEN            # member mute
glo_color_gc=$LIGHT_YELLOW     # global caller
glo_color_gm=$YELLOW           # global mute
glo_color_ex=$BLUE             # external

glo_exclude_dir=""   # directory to be excluded from the search
glo_expand_all="no"  # are you expanding them all?
glo_ignore_mod=""    # global list of modules to ignore

#------------------------------------------------------------------------------#
#   Some other global variables needed for functionality
#------------------------------------------------------------------------------#
analyzed_units=""   # list of analyzed units, to avoid duplication
glo_procedure=""    # result from extract_procedure_and_module
glo_module=""       # result from extract_procedure_and_module

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

}  # print_separator

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

  if [[ $lev ]]; then
    echo -e "$ind"${color}"$bullet""$proced"" ($lev)"${RESET}"\t ""$module"
  else
    echo -e "$ind"${color}"$bullet""$proced"${RESET}"\t ""$module"
  fi

}  # print_line

#==============================================================================#
# Print_usage
#------------------------------------------------------------------------------#
print_usage() {
  echo "#======================================================================"
  echo "# Utility for extract call graphs from T-Flows"
  echo "#----------------------------------------------------------------------"
  echo "# Proper usage: "
  echo "#"
  echo "# "$0" <Source> [options]"
  echo "#"
  echo "# where Source is the procedure name for which you want to perform"
  echo "# the analysis, such as: Main_Con, Compute_Energy, Main_Div, hence"
  echo -e "# case sensitive, with or without the .f90 extension."
  echo "#"
  echo "# Valid options are:"
  echo "#"
  echo "# -a"
  echo "#    Expand all. Don't contract units which have been expanded above."
  echo "#"
  echo "# -e <list of directories to exclude>"
  echo "#    In cases where the same module name is used in more than one"
  echo "#    directory, use this option to exclude one from the search."
  echo "#"
  echo "# -i <list of modules to ignore>"
  echo "#    You may want to exclude some of the smaller modules, such as"
  echo "#    Comm_Mod, Message_Mod, Work_Mod, Profiler_Mod, String_Mod,"
  echo "#    Tokenizer_Mod to reduce the amoun of information printed."
  echo "#"
  echo "# Note that:"
  echo -e "#   ${LIGHT_RED} The script is supposed to be executed from:" \
          "T-Flows/Sources,"${RESET}
  echo -e "#   ${LIGHT_RED} or any of its sub-directories. The <Target_Mod>" \
          "doesnot have"${RESET}
  echo -e "#   ${LIGHT_RED} to reside in the directory from which you launch" \
          "the script!"${RESET}
  echo "#----------------------------------------------------------------------"
  exit

}  # print_usage

#==============================================================================#
# Sets $glo_procedure and $glo_module
#------------------------------------------------------------------------------#
extract_procedure_and_module() {

  full_name=$1

  # It rarely gets cooler than this: extract the calling function pattern :-)
  pattern="${full_name//[^%()]}"

  # Take care of the lines in which call is not the first statement
  if [[ ${#pattern} == 1 ]]; then
    if [[ ${pattern:0:1} == "%" || ${pattern:0:1} == ")" ]]; then
      full_name=$(awk -F'call ?' '{print $2}' <<< ${called_procedures[proc]})
      pattern="${full_name//[^%()]}"
    fi
  fi

  # But it gets even better after all - patterns
  # can be characterized by four first characters
  first_one=${pattern:0:1}
  first_two=${pattern:0:2}
  first_four=${pattern:0:4}

  # This should take care of Comm_Mod_Friendly_Function
  if [[                                 \
        $pattern   == ""            ||  \
        $full_name == *"_Mod_"*         \
     ]]; then
    glo_module=$(echo $full_name | awk -F '_Mod_' '{print $1}')"_Mod"
    glo_procedure=$(echo $full_name | awk -F '_Mod_' '{print $2}')
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  # Take care of cases such as Grid%Comm%Sendrecv_Real_Arrays
  # Front%Elem(e)%Initialize_Elem()
  elif [[                               \
          $first_two  == "%%"    ||     \
          $first_four == "%()%"         \
       ]]; then
    glo_module=$(cut -d % -f 2 <<< $full_name)
    glo_module=$(cut -d \( -f 1 <<< $glo_module)
    glo_procedure=$(cut -d % -f 3 <<< $full_name)
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  # Likes of Profiler%Start() and Grid(d)%Calculate() ...
  # Profiler%Update_By_Rank(Profiler%previously_running)
  elif [[                               \
          $first_one  == "%"     ||     \
          $first_two  == "%("    ||     \
          $first_two  == ")%"    ||     \
          $first_four == "()%("         \
       ]]; then
    glo_module=$(cut -d %  -f 1 <<< $full_name)
    glo_module=$(cut -d \( -f 1 <<< $glo_module)
    glo_procedure=$(cut -d %  -f 2 <<< $full_name)
    glo_procedure=$(cut -d \( -f 1 <<< $glo_procedure)

  # This is for global functions and external functions like Probe_1d()
  # or system_clock(count_rate=Profiler%sys_count_rate)
  elif [[                               \
          $first_one == "("      ||     \
          $first_one == ""              \
       ]]; then
    glo_module="" # empty string (maybe none or Shared?)
    glo_procedure=$(cut -d \( -f 1 <<< $full_name)

  elif [[ $pattern == ")%()" ]]; then
    echo "Warning: complex pattern 'call' is not in the beginning"
    echo "Full name is: " $full_name
    exit

  else
    echo "Unknown pattern: " $pattern
    echo "in             : " $full_name
    exit
  fi

}  # extract_procedure_and_module

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
  local module_in_which_you_seek="$2"

  #   Level counters - these are used to indent the tree
  local next_level=`expr $next_level + 1`
  local this_level=`expr $next_level - 1`

  if [[ $this_level == 0 ]]; then
    echo "#======================================================================="
    echo "# Extracting module hierarchy for "$procedure_name_you_seek
    echo "#"
    echo "# Procedures are designated as follows:"
    echo "#"
    echo -n "# - members calling others:    "
    echo -e "${glo_color_mc}"• Member_Caller" (level)${RESET}       Module_Mod"
    echo -n "# - members not calling any:   "
    echo -e "${glo_color_mm}"⨯ Member_Mute" (level)${RESET}         Module_Mod"
    echo -n "# - global calling others:     "
    echo -e "${glo_color_gc}"• Global_Caller"${RESET}"
    echo -n "# - global not calling any:    "
    echo -e "${glo_color_gm}"⨯ Global_Mute  "${RESET}"
    echo -n "# - external procedure:        "
    echo -e "${glo_color_ex}"⨯ External"${RESET}"
    echo "#-----------------------------------------------------------------------"
  fi

  #---------------------------------------------------------------------
  #   Is the module you are currently analyzing on the excluded list?
  #---------------------------------------------------------------------
  if [[ $glo_ignore_mod != *$module_in_which_you_seek* ]] ||  \
     [[ ! $module_in_which_you_seek ]]; then

    #-----------------------------------------------------------------------------
    #   Get the full path of the module you seek
    #
    #   Some typical directories are excluded from here.  For example, sources
    #   in "Seqential" part of the Comm_Mod are just empty hooks, no one in
    #   the sane mind will be interested in analyzing them.
    #-----------------------------------------------------------------------------

    # Start building the find command as a string
    find_command="find $src -name $procedure_file_you_seek"

    # Add a condition to grep for the module if it's specified
    if [[ $module_in_which_you_seek ]]; then
      find_command+=" | grep $module_in_which_you_seek"
    fi

    # Exclude standard directories
    find_command+=" | grep -v No_Checking"
    find_command+=" | grep -v Sequential"
    find_command+=" | grep -v Unused"
    find_command+=" | grep -v Fake"

    # Convert glo_exclude_dir to an array and add a grep -v for each directory
    if [[ $glo_exclude_dir ]]; then
      read -r -a exclude_dirs <<< "$glo_exclude_dir"
      for exclude_dir in "${exclude_dirs[@]}"; do
        if [[ $exclude_dir ]]; then
          find_command+=" | grep -v $exclude_dir"
        fi
      done
    fi

    # Execute the find command
    local full_path_you_seek=$(eval $find_command)

    # This command counts number of occurrences of modules name in the result of
    # command find. If it is more than one, the same file is in more directories
    n=$(echo "$full_path_you_seek" | tr " " "\n" | grep -c "$procedure_name_you_seek")
    if [[ $n > 1 ]]; then
      echo "Ambiguity: procedure "$procedure_name_you_seek
           "found in more than one directory, here is the list:"
      for path in ${full_path_you_seek[*]}; do
        echo $path
      done
      echo "Exclude all but one directory with"\
           "the command line argument -e <directory>"
      exit
    fi

    # If the definition of this procedure is found, carry on
    if [[ $full_path_you_seek ]]; then

      #-----------------------------------------------------
      #   Storing results of the grep command in an array
      #-----------------------------------------------------
      local called_procedures=($(grep '  call ' $full_path_you_seek  \
                                 | awk '{print $2$3$4$5$6$7$8$9}' | tr -d ,))
      local called_modules=$called_procedures   # just to declare

      # Correct the lines in which call was not the first comman, something
      # like: if(some_condition_is_met) call The_Function_You_Seek
      for proc in "${!called_procedures[@]}"; do
        if [[ ${called_procedures[proc]} == *"call"* ]]; then  # call exists
          called_procedures[proc]=$(awk -F'call ?' '{print $2}' <<< ${called_procedures[proc]})
        fi
      done

      # At this point, $called procedures has a form like: Profiler%Start('Main')
      # From this mess, extract the module name and the procedure element wise
      for proc in "${!called_procedures[@]}"; do

        # Set $glo_procedure and $glo_module
        extract_procedure_and_module ${called_procedures[proc]}

        # Typical substitues:
        if [[ $glo_module != "" ]]; then
          if [[ $glo_module == "Backup" ]];       then glo_module=Backup_Mod;       fi
          if [[ $glo_module == "Comm" ]];         then glo_module=Comm_Mod;         fi
          if [[ $glo_module == "Control" ]];      then glo_module=Control_Mod;      fi
          if [[ $glo_module == "Convert" ]];      then glo_module=Convert_Mod;      fi
          if [[ $glo_module == "Divide" ]];       then glo_module=Divide_Mod;       fi
          if [[ $glo_module == "Dom" ]];          then glo_module=Domain_Mod;       fi
          if [[ $glo_module == "Elem" ]];         then glo_module=Elem_Mod;         fi
          if [[ $glo_module == "File" ]];         then glo_module=File_Mod;         fi
          if [[ $glo_module == "Flow" ]];         then glo_module=Field_Mod;        fi
          if [[ $glo_module == "Front" ]];        then glo_module=Front_Mod;        fi
          if [[ $glo_module == "Generate" ]];     then glo_module=Generate_Mod;     fi
          if [[ $glo_module == "Grid" ]];         then glo_module=Grid_Mod;         fi
          if [[ $glo_module == "Isoap" ]];        then glo_module=Isoap_Mod;        fi
          if [[ $glo_module == "Iso_Polygons" ]]; then glo_module=Iso_Polygons_Mod; fi
          if [[ $glo_module == "Info" ]];         then glo_module=Info_Mod;         fi
          if [[ $glo_module == "Iter" ]];         then glo_module=Iter_Mod;         fi
          if [[ $glo_module == "Line" ]];         then glo_module=Tokenizer_Mod;    fi
          if [[ $glo_module == "Math" ]];         then glo_module=Math_Mod;         fi
          if [[ $glo_module == "Message" ]];      then glo_module=Message_Mod;      fi
          if [[ $glo_module == "Metis" ]];        then glo_module=Metis_Mod;        fi
          if [[ $glo_module == "Monitor" ]];      then glo_module=Monitor_Mod;      fi
          if [[ $glo_module == "Msg" ]];          then glo_module=Message_Mod;      fi
          if [[ $glo_module == "Nat" ]];          then glo_module=Native_Mod;       fi
          if [[ $glo_module == "Particle" ]];     then glo_module=Particle_Mod;     fi
          if [[ $glo_module == "Pet" ]];          then glo_module=Petsc_Mod;        fi
          if [[ $glo_module == "Pol" ]];          then glo_module=Polyhedron_Mod;   fi
          if [[ $glo_module == "Polyhedron" ]];   then glo_module=Polyhedron_Mod;   fi
          if [[ $glo_module == "Por" ]];          then glo_module=Porosity_Mod;     fi
          if [[ $glo_module == "Process" ]];      then glo_module=Process_Mod;      fi
          if [[ $glo_module == "Prof" ]];         then glo_module=Profiler_Mod;     fi
          if [[ $glo_module == "Profiler" ]];     then glo_module=Profiler_Mod;     fi
          if [[ $glo_module == "Read_Control" ]]; then glo_module=Read_Control_Mod; fi
          if [[ $glo_module == "Results" ]];      then glo_module=Results_Mod;      fi
          if [[ $glo_module == "Sol" ]];          then glo_module=Solver_Mod;       fi
          if [[ $glo_module == "Sort" ]];         then glo_module=Sort_Mod;         fi
          if [[ $glo_module == "Stl" ]];          then glo_module=Stl_Mod;          fi
          if [[ $glo_module == "String" ]];       then glo_module=String_Mod;       fi
          if [[ $glo_module == "Surf" ]];         then glo_module=Surf_Mod;         fi
          if [[ $glo_module == "Swarm" ]];        then glo_module=Swarm_Mod;        fi
          if [[ $glo_module == "Tok" ]];          then glo_module=Tokenizer_Mod;    fi
          if [[ $glo_module == "Time" ]];         then glo_module=Time_Mod;         fi
          if [[ $glo_module == "Turb" ]];         then glo_module=Turb_Mod;         fi
          if [[ $glo_module == "Vof" ]];          then glo_module=Vof_Mod;          fi
          if [[ $glo_module == "Work" ]];         then glo_module=Work_Mod;         fi
        fi
        called_modules[$proc]=$glo_module
        called_procedures[$proc]=$glo_procedure

      done

      #------------------------------------------------------------------
      #   Form the current procedure / module combination for the unit
      #------------------------------------------------------------------
      if [[ $module_in_which_you_seek ]]; then
        local current_unit="$procedure_name_you_seek@$module_in_which_you_seek "
      else
        local current_unit="$procedure_name_you_seek@global "
      fi

      #------------------------------------------------------------------
      #   Print out the name of the module you are currently analysing
      #------------------------------------------------------------------
      if [[ $called_procedures ]]; then

        if [[ $module_in_which_you_seek ]]; then
          if [[ $analyzed_units != *$current_unit* ]]; then
            print_separator "$indent" $this_level
          fi
          print_line "$indent"                 \
                     $glo_color_mc             \
                     "• "                      \
                     $procedure_name_you_seek  \
                     $this_level               \
                     $module_in_which_you_seek
        else
          print_line "$indent"                 \
                     $glo_color_gc             \
                     "• "                      \
                     $procedure_name_you_seek  \
                     $this_level               \
                     "global"
        fi
      else
        if [[ $module_in_which_you_seek ]]; then
          print_line "$indent"                 \
                     $glo_color_mm             \
                     "⨯ "                      \
                     $procedure_name_you_seek  \
                     $this_level               \
                     $module_in_which_you_seek
        else
          print_line "$indent"                 \
                     $glo_color_gm             \
                     "⨯ "                      \
                     $procedure_name_you_seek  \
                     $this_level               \
                     "global"
        fi
      fi

      #-----------------------------------------------------
      #   If current units was analyzed, get out of here.
      #   Oterwise, update the list of units and continue
      #-----------------------------------------------------
      if [[ $analyzed_units == *$current_unit* ]]; then
        return
      fi
      if [[ $analyzed_units != *$current_unit* ]]; then
        if [[ $glo_expand_all == "no" ]]; then
          analyzed_units=$analyzed_units" $current_unit"
        fi
      fi

      # Increase indend for the next level by appending 6 spaces to it
      local indent="${indent}"$glo_indent

      #-------------------------------------------------------------
      #   If the list of called procedures in not empty, carry on
      #-------------------------------------------------------------
      if [[ $called_procedures ]]; then

        # Print the procedures which are called from here
        for proc in "${!called_procedures[@]}"; do

          # This seems to be the only place from which you can print
          # non-member / external functions which don't make any calls.  Period!
          # Yet, at this point you can't tell one from another, and that's why
          # another call to find is done here
          if [[ ! ${called_modules[proc]} ]]; then
            internal=$(find $src -type f -name ${called_procedures[proc]}".f90")
            if [[ ! $internal ]]; then
              print_line "$indent"                   \
                         $glo_color_ex               \
                         "⨯ "                        \
                         ${called_procedures[proc]}  \
                         ""                          \
                         "external"
            fi
          fi

          #-------------------------------------------------------#
          #                                                       #
          #   The very important recursive call to its own self   #
          #                                                       #
          #-------------------------------------------------------#
          extract_call_graph ${called_procedures[proc]} ${called_modules[proc]}
        done
      fi

    fi  # if the paths is found (was not excluded by grep -v ...

  fi  # not on the ignore list

}  # extract_call_graph

#------------------------------------------------------------------------
#   Wrong number of command line argument is sent - describe the usage
#------------------------------------------------------------------------
if [[ ! $1 ]]; then
  print_usage
fi

#-----------------------------------------------
#   Parse command line options like a pro :-)
#-----------------------------------------------

# Fetch the name
name=$1

# Remove the .f90 extension if it exists
name="${name%.f90}"

# Shift on
shift

current_opt=""

while [[ $# > 0 ]]; do
  case $1 in
    # All - expand all, don't contract already analyzed units
    -a)
      current_opt=$1
      glo_expand_all="yes"
      shift  # past argument
      ;;     # part of the case construct

    # Exclude - accumulate arguments
    -e)
      current_opt=$1
      glo_exclude_dir=$glo_exclude_dir" $2"
      shift  # past argument
      shift  # past value
      ;;     # part of the case construct

    # Ignore - accumulates arguments
    -i)
      current_opt=$1
      glo_ignore_mod=$glo_ignore_mod" $2"
      shift  # past argument
      shift  # past value
      ;;     # part of the case construct

    # Still don't have use for this, but is not asking for food, so keep it
    --default)
      current_opt=$1
      default=yes
      shift  # past argument
      ;;     # part of the case construct

    -*)
      echo "Unknown option $1"
      exit 1
      ;;     # part of the case construct

    # Accumulates additonal strings to glo_ignore
    *)
      if [[ $current_opt == -e ]]; then
        glo_exclude_dir=$glo_exclude_dir" $1"
      elif [[ $current_opt == -i ]]; then
        glo_ignore_mod=$glo_ignore_mod" $1"
      else
        echo "Unknown option $1"
        exit 1
      fi
      shift  # past argument
      ;;     # part of the case construct
  esac
done

extract_call_graph $name
