#!/bin/bash
# This script automatically builds 'openmpi' and 'mpich' libs
# on current C, C++, Fortran compilers
# and compiles T-Flows in parallel mode on each of them
#
# mpich lib   major versions:
# 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 3.0, 3.1, 3.2, 3.3
#
# openmpi lib major versions:
# 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.10, 2.0, 2.1, 3.0, 3.1, 4.0
#
# Prerequisites: C, C++, Fortran compilers
# makefile for Process must have only one row with 'mpif90'
#
# This script should use only compiled mpif90 and mpi.h,
# however I recommend to purge all mpi libs installed before launching this


# Repeat test for CGNS=yes
TEST_CGNS=false

#-----------------------------
# URLs with all major releases
#-----------------------------

MPICH_LIB="https://www.mpich.org/static/downloads"
OPENMPI_LIB="https://download.open-mpi.org/release/open-mpi"

MPICH_TESTS=("$MPICH_LIB""/1.0/mpich2-1.0.tar.gz" \
             "$MPICH_LIB""/1.1/mpich2-1.1.tar.gz" \
             "$MPICH_LIB""/1.2/mpich2-1.2.tar.gz" \
             "$MPICH_LIB""/1.3/mpich2-1.3.tar.gz" \
             "$MPICH_LIB""/1.4/mpich2-1.4.tar.gz" \
             "$MPICH_LIB""/1.5/mpich2-1.5.tar.gz" \
             "$MPICH_LIB""/3.0/mpich-3.0.tar.gz"  \
             "$MPICH_LIB""/3.1/mpich-3.1.tar.gz"  \
             "$MPICH_LIB""/3.2/mpich-3.2.tar.gz"  \
             "$MPICH_LIB""/3.3/mpich-3.3.tar.gz"  )

OPENMPI_TESTS=("$OPENMPI_LIB""/v1.0/openmpi-1.0.tar.gz"     \
               "$OPENMPI_LIB""/v1.1/openmpi-1.1.tar.gz"     \
               "$OPENMPI_LIB""/v1.2/openmpi-1.2.tar.gz"     \
               "$OPENMPI_LIB""/v1.3/openmpi-1.3.tar.gz"     \
               "$OPENMPI_LIB""/v1.4/openmpi-1.4.tar.gz"     \
               "$OPENMPI_LIB""/v1.5/openmpi-1.5.tar.gz"     \
               "$OPENMPI_LIB""/v1.6/openmpi-1.6.tar.gz"     \
               "$OPENMPI_LIB""/v1.7/openmpi-1.7.tar.gz"     \
               "$OPENMPI_LIB""/v1.8/openmpi-1.8.tar.gz"     \
               "$OPENMPI_LIB""/v1.10/openmpi-1.10.0.tar.gz" \
               "$OPENMPI_LIB""/v2.0/openmpi-2.0.0.tar.gz"   \
               "$OPENMPI_LIB""/v2.1/openmpi-2.1.0.tar.gz"   \
               "$OPENMPI_LIB""/v3.0/openmpi-3.0.0.tar.gz"   \
               "$OPENMPI_LIB""/v3.1/openmpi-3.1.0.tar.gz"   \
               "$OPENMPI_LIB""/v4.0/openmpi-4.0.0.tar.gz"   )

#-----------------
# Folder structure
#-----------------

TOP_DIR=$PWD                 # this is top dir for this test
INS_DIR=$TOP_DIR/install_dir # this dir contains compiled mpi libs
SRC_DIR=$TOP_DIR/src_dir     # this dir contains downloaded sources
rm -rf $INS_DIR $SRC_DIR
mkdir -p $SRC_DIR/
mkdir -p $INS_DIR/

TEST_DIR="$PWD"
GENE_DIR=$PWD/../Sources/Generate  # Generate src folder
CONV_DIR=$PWD/../Sources/Convert   # Convert  src folder
DIVI_DIR=$PWD/../Sources/Divide    # Divide   src folder
PROC_DIR=$PWD/../Sources/Process   # Process  src folder
BINA_DIR=$PWD/../Binaries/         # binaries folder

# Current environment vars
PATH_CLEAN="$PATH"
LD_LIBRARY_PATH_CLEAN="$LD_LIBRARY_PATH"

# Script temp dir
TEMP_DIR=`mktemp -d 2>/dev/null || mktemp -d -t 'TEMP_DIR'`

#------------------------------------------------------------------------------#

# C compiler
CCOMP='gcc'
# C++ compiler
CPPCOMP='g++'
# F90 (or later) compiler
FCOMP='gfortran'

# Start time measurements from this moment
current_time=$(date +%s)

#------------------------------------------------------------------------------#
# time in seconds
#------------------------------------------------------------------------------#
function time_in_seconds {
  previous_time=$current_time
  current_time=$(date +%s)
  echo "time elapsed:" \
  "$(echo ""$current_time" - "$previous_time"" | bc -l)" "seconds"
}
#------------------------------------------------------------------------------#
#   Sets compilers, optimization flags and environment                         #
#------------------------------------------------------------------------------#
function switch_compilers_and_environment_to {
  #$1 - "sequential", "parallel"
  #$2 - path

  # renew environment
  unset PATH
  export PATH="$PATH_CLEAN"
  unset LD_LIBRARY_PATH
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH_CLEAN"

  export CCLAGS='-O2 -fPIC'
  export CPPFLAGS='-O2 -fPIC'
  export CXXFLAGS='-O2 -fPIC'
  export CFLAGS='-O2 -fPIC'
  export FFLAGS='-O2 -fPIC'
  export FCFLAGS='-O2 -fPIC'

  if [ $1 = "parallel" ]; then
    export CC=$2/bin/mpicc
    export FC=$2/bin/mpif90
    export CXX=$2/bin/mpicxx
    export PATH=$2/bin:$PATH
    export LD_LIBRARY_PATH=$2/lib/:$LD_LIBRARY_PATH
  elif [ $1 = "sequential" ]; then
    export CC=$CCOMP
    export CXX=$CPPCOMP
    export FC=$FCOMP
  else
    exit 1
    echo "Error in switch_compilers_and_environment_to: wrong argument"
  fi
}
#------------------------------------------------------------------------------#
# Replace line with first occurence of string to arg in file
#------------------------------------------------------------------------------#
function replace_line_with_first_occurence_in_file {
  #$1 string to search in file
  #$2 new line
  #$3 file

  new_line=$(echo "$2")
  line_to_replace="$(grep -ni "$1" $3 | cut -d: -f1 | head -n1)"
  if [ -z "$line_to_replace" ]; then
     # if no such string is found-> add it
    echo $new_line >> $3
    echo "Warning: added "$2" to the end of $3"
  else
    sed -i ""$line_to_replace"s%.*%$new_line%" $3
  fi

  # previous approach was with
  # awk -v var="$line" 'NR==var {$2=3}1'  control > .control.tmp
}
#------------------------------------------------------------------------------#
#   MPICH TESTS                                                                #
#------------------------------------------------------------------------------#
#   Builds mpich lib and mpicc, mpicxx, mpif90 compilers                       #
#------------------------------------------------------------------------------#
function mpich_full_test {
  for i in ${!MPICH_TESTS[@]}; do
    ARCHIVE_TO_DOWNLOAD="${MPICH_TESTS[$i]}"
    FILENAME="$(echo "${ARCHIVE_TO_DOWNLOAD##*/}")"
    FILENAME_NO_EXT="$(echo "$FILENAME" | sed 's/.tar.gz//g')"

    # Log for MPICH_TESTS[$i]
    LOG_FILE="$SRC_DIR""/$FILENAME_NO_EXT"".log"
    if [ -f $LOG_FILE ]; then
      cp /dev/null $LOG_FILE
    fi

    echo ""
    echo "#===================================================================="
    echo "#   MPICH TEST: ""$FILENAME"
    echo "#--------------------------------------------------------------------"
    echo "#  Downloading sources: ""$ARCHIVE_TO_DOWNLOAD"
    cd $SRC_DIR
    wget -N "$ARCHIVE_TO_DOWNLOAD" >> $LOG_FILE 2>&1

    tar -zxf "$FILENAME" >> $LOG_FILE 2>&1

    MPICH_SRC_DIR="$SRC_DIR""/$FILENAME_NO_EXT" # directory with sources
    MPICH_INS_DIR="$INS_DIR""/$FILENAME_NO_EXT" # directory with binaries
    cd $MPICH_SRC_DIR

    echo '   Setting compilers'
    switch_compilers_and_environment_to "sequential"

    echo '   Configuring installation'

    ./configure \
    --prefix=$MPICH_INS_DIR \
    --enable-fast=all,O2 \
    --enable-fortran=all \
    >> $LOG_FILE 2>&1

    echo '   Building'
    make >> $LOG_FILE 2>&1

    echo '   Installing'
    make install >> $LOG_FILE 2>&1

    echo '   Testing build'
    switch_compilers_and_environment_to "parallel" "$MPICH_INS_DIR"
    #make check-recursive >> $LOG_FILE 2>&1

    echo '   Compiling T-Flows with this build'
    cd "$PROC_DIR"
    make clean >> $LOG_FILE 2>&1

    replace_line_with_first_occurence_in_file "mpif90" \
      "FC = ""$MPICH_INS_DIR""/bin/mpif90"" -I""$MPICH_INS_DIR""/include/" \
      makefile

    make   MPI=yes CGNS_HDF5=no         MPI_TYPE=mpich >> $LOG_FILE 2>&1

    if [ $TEST_CGNS == true ]; then
      make MPI=yes CGNS_HDF5=$TEST_CGNS MPI_TYPE=mpich >> $LOG_FILE 2>&1
    fi

    echo "   mpich lib ""$FILENAME"" has finished"

  done
}
#------------------------------------------------------------------------------#
#   OPENMPI TESTS                                                                #
#------------------------------------------------------------------------------#
#   Builds openmpi lib and mpicc, mpicxx, mpif90 compilers                       #
#------------------------------------------------------------------------------#
function openmpi_full_test {
  for i in ${!OPENMPI_TESTS[@]}; do
    ARCHIVE_TO_DOWNLOAD="${OPENMPI_TESTS[$i]}"
    FILENAME="$(echo "${ARCHIVE_TO_DOWNLOAD##*/}")"
    FILENAME_NO_EXT="$(echo "$FILENAME" | sed 's/.tar.gz//g')"

    # Log for OPENMPI_TESTS[$i]
    LOG_FILE="$SRC_DIR""/$FILENAME_NO_EXT"".log"
    if [ -f $LOG_FILE ]; then
      cp /dev/null $LOG_FILE
    fi

    echo ""
    echo "#===================================================================="
    echo "#   OPENMPI TEST: ""$FILENAME"
    echo "#--------------------------------------------------------------------"
    echo "#  Downloading sources: ""$ARCHIVE_TO_DOWNLOAD"
    cd $SRC_DIR
    wget -N "$ARCHIVE_TO_DOWNLOAD" >> $LOG_FILE 2>&1

    tar -zxf "$FILENAME" >> $LOG_FILE 2>&1

    OPENMPI_SRC_DIR="$SRC_DIR""/$FILENAME_NO_EXT" # directory with sources
    OPENMPI_INS_DIR="$INS_DIR""/$FILENAME_NO_EXT" # directory with binaries
    cd $OPENMPI_SRC_DIR

    echo '   Setting compilers'
    switch_compilers_and_environment_to "sequential"

    echo '   Configuring installation'

    ./configure \
    --prefix=$OPENMPI_INS_DIR \
    >> $LOG_FILE 2>&1

    echo '   Building'
    make >> $LOG_FILE 2>&1

    echo '   Installing'
    make install >> $LOG_FILE 2>&1

    echo '   Testing build'
    switch_compilers_and_environment_to "parallel" "$OPENMPI_INS_DIR"
    #make check-recursive >> $LOG_FILE 2>&1

    echo '   Compiling T-Flows with this build'
    cd "$PROC_DIR"
    make clean >> $LOG_FILE 2>&1

    replace_line_with_first_occurence_in_file "mpif90" \
      "FC = ""$OPENMPI_INS_DIR""/bin/mpif90"" -I""$OPENMPI_INS_DIR""/include/" \
      makefile

    make   MPI=yes CGNS_HDF5=no         MPI_TYPE=openmpi >> $LOG_FILE 2>&1

    if [ $TEST_CGNS == true ]; then
      make MPI=yes CGNS_HDF5=$TEST_CGNS MPI_TYPE=openmpi >> $LOG_FILE 2>&1
    fi

    echo "   openmpi lib ""$FILENAME"" has finished"

  done
}
#------------------------------------------------------------------------------#
#   Remove temp                                                                #
#------------------------------------------------------------------------------#
function remove_temp_dir {
  # remove tmp dir
  rm -rf $TEMP_DIR
}
#------------------------------------------------------------------------------#
#   Actual script with functions defined above                                 #
#------------------------------------------------------------------------------#

mpich_full_test
openmpi_full_test
remove_temp_dir

echo done
