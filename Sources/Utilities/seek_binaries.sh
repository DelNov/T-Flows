#!/bin/bash

one_lev="../"
cur_lev=""
dir_nam="Binaries"

#---------------------------------------------------------------------
#
#   Browse through upper levels until you find directory "Binaries"
#
#---------------------------------------------------------------------
for i in {1..16}
do

  #---------------------------------------
  # Set the name of the current directory
  #---------------------------------------
  cur_lev="$cur_lev$one_lev"
  cur_dir="$cur_lev$dir_nam"

  #---------------------------------------------
  # If this directory exists, create soft links
  #---------------------------------------------
  if [ -d "$cur_dir" ]; then
    echo "Found" $dir_nam "in:" $cur_dir
    echo "Creating soft links to executables here!"

    # Link Convert
    if [ -f "${cur_dir}"/Conv* ]; then
      ln -s -f ${cur_dir}/Conv* .
    fi

    # Link Divide
    if [ -f "${cur_dir}"/Divi* ]; then
      ln -s -f ${cur_dir}/Divi* .
    fi

    # Link Generate
    if [ -f "${cur_dir}"/Gene* ]; then
      ln -s -f ${cur_dir}/Gene* .
    fi

    # Link Process
    if [ -f "${cur_dir}"/Proc* ]; then
      ln -s -f ${cur_dir}/Proc* .
    fi

  fi

done

