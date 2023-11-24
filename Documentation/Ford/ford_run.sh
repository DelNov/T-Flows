#==============================================================================#
#   Script Description: Generate FORD Documentation for T-Flows Subprojects
#------------------------------------------------------------------------------#
#   Overview
#
#   This script automates the process of generating FORD (FORtran Documentation)
#   for specific subprojects within the T-Flows project. It ensures that only
#   relevant source files are included in the documentation and standardizes
#   preprocessor #include statements within these files.
#
#   Key Steps
#
#   1. Extract Source File Names:
#      - The script reads the specified makefile (makefile_path) to extract
#        names of Fortran source files (.f90).
#      - It ignores placeholders and patterns, such as % and *, ensuring only
#        actual file names are extracted.
#
#   2. Copy Relevant Files to Temporary Directory:
#      - A temporary directory (src_temp) is created to hold the sources for
#        FORD processing.
#      - The script searches for and copies each source file and its correspon-
#        ding directory from the root directory (src_root) to this temporary
#        directory.
#      - This step includes both the local source files and the shared modules
#         specified in the makefile.
#
#   3. Standardize Include Statements:
#      - The script then modifies all .f90 files in the temporary directory:
#        - Transforms all #include statements for .f90 files into standard
#          include statements, ensuring FORD can process them correctly.
#        - Removes any include statements that refer to .h90 files, which are
#          not compatible with FORD's processing.
#
#   4. Generate Documentation with FORD:
#      - Runs the FORD tool to generate documentation based on the sources in
#        the temporary directory and the FORD configuration file (ford_con.md).
#
#   5. Cleanup:
#      - After documentation generation, the temporary directory is removed,
#        ensuring no residual files are left behind.
#
#   Usage
#   - This script is designed to be run from specific subdirectories within the
#     T-Flows project structure (Convert, Generate, Divide, Process).
#   - It reads the project's makefile to determine the relevant source files
#     and modules for documentation.
#
#   Dependencies
#   - Requires the FORD tool to be installed and configured for the T-Flows
#     project.
#   - Assumes a specific project structure and makefile format for T-Flows.
#==============================================================================#

#!/bin/bash

# Paths to the makefile and source directories
makefile_path="$1"/makefile

# Array of root directories for sources
src_roots=("$1" "$1"../Shared)

# Temporary directory for FORD
src_temp="./Temporary"

# Create temporary directory
mkdir -p "$src_temp"

# Array to hold all source files without the .f90 extension
declare -a src_files

rm -fR ./Html

#-------------------------------------------------
# Read the makefile and extract source file names
#-------------------------------------------------
while IFS= read -r line; do
  if [[ "$line" =~ \.f90 ]]; then
    # Split the line into words and check each word
    for word in $line; do
      # Add to src_files if the word ends with .f90 and does not contain % or *
      if [[ "$word" == *.f90 && "$word" != *%* && "$word" != *\** ]]; then
        # Remove the .f90 extension and add to array
        src_file="${word%.f90}"
        if [[ ! " ${src_files[*]} " =~ " ${src_file} " ]]; then
          src_files+=("$src_file")
        fi
      fi
    done
  fi
done < "$makefile_path"

#----------------------------------------------------------------
# Search for, copy the paths of the source files and directories
#----------------------------------------------------------------
for file in "${src_files[@]}"; do
  for src_root in "${src_roots[@]}"; do
    found_files=( $(find "$src_root" -name "$file.f90" -o -name "$file" -type d) )
    for found_file in "${found_files[@]}"; do
      echo "Copying $found_file to $src_temp"
      cp -R "$found_file" "$src_temp"
    done
  done
done

#-------------------------------------------------------------------
# Apply transformation to the include statement in the copied files
#-------------------------------------------------------------------
# find "$src_temp" -type f -name "*.f90" -exec sed -i -e '/#.*include.*\.f90/ s/#\s*include/include/g' {} \;
find "$src_temp" -type f -name "*.f90" -exec sed -i -e 's/#\s\+include/ include/g' {} \;
find "$src_temp" -type f -name "*.f90" -exec sed -i -e '/include.*\.h90/d' {} \;

echo "Include statements have been standardized in the temporary source directory."

#----------
# Run FORD
#----------
ford ford_project.md

#----------
# Clean up
#----------
rm -fR $src_temp
