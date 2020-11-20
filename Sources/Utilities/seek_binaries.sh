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
    if [ -f "${cur_dir}"/C* ]; then
      ln -s -f ${cur_dir}/C* .
    fi

    # Link Divide
    if [ -f "${cur_dir}"/D* ]; then
      ln -s -f ${cur_dir}/D* .
    fi

    # Link Generate
    if [ -f "${cur_dir}"/G* ]; then
      ln -s -f ${cur_dir}/G* .
    fi

    # Link Process
    if [ -f "${cur_dir}"/P* ]; then
      ln -s -f ${cur_dir}/P* .
    fi

  fi

done

