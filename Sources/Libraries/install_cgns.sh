#!/bin/bash
# Requires: C, C++, Fortran compilers
# If you want cgnstools(graphical tools to view .cgns files) then install:
# libx11-dev libxmu-headers libxmu-dev packages for debian/ubuntu,
# libX11-devel for openSUSE,
# xorg-x11-devel for Fedora

# This script will automatically build cgns lib in several releases:
# 1) with hdf5 and parallel   access (optional)
# 2) with hdf5 and sequential access + cgnstools (optional)
# (parallel cgns lib can only be based on hdf5, not adf5)

# Hdf5_Par lib most likely can not be built with your current mpif90
# therefore script builds mpif90 first to compile Hdf5_Par later

# You can optionally build cgnstools to view/edit .cgns files
# https://cgns.github.io/CGNS_docs_current/cgnstools/cgnsview/index.html
# if this case script will download and compile dependencies: tcl, tk libs

# T-Flows supports currently only 3.2.1 cgns lib, but
# you can build with this script the latest cgns lib if you wish
# since 3.3 version fortran 77 and 90 support was dropped,
# moreover user is forced to use cgsize_t type instead of integer in cgns API

# CGNS lib does not use compression:
# https://cgns.github.io/FAQs.html

# build mpich, otherwise Hdf5_Par can not be built
BUILD_MPI=false

# build cgnstools, which contains cgnsview - small program to view .cgns files
# https://cgns.github.io/CGNS_docs_current/cgnstools/cgnsview/index.html
BUILD_CGNS_TOOLS=true

# build latest cgns lib, otherwise build 3.2.1 version
BUILD_LATEST_CGNS=false

# folder structure
CGNS_DIR=$PWD                     # this is top dir for this lib
INSTALL_DIR=$CGNS_DIR/install_dir # this dir contains Cgns/ Hdf5/ Mpich/
SRC_DIR=$CGNS_DIR/src_dir         # this dir contains downloaded sources
rm -rf $INSTALL_DIR $SRC_DIR

# script temp dir
TEMP_DIR=`mktemp -d 2>/dev/null || mktemp -d -t 'TEMP_DIR'`
# script logs
LOG_FILE=$SRC_DIR/build_cgns.log # logs of current script
if [ -f $LOG_FILE ]; then
  cp /dev/null $LOG_FILE
fi

# C compiler
CCOMP='gcc'
# C++ compiler
CPPCOMP='g++'
# F90 (or later) compiler
FCOMP='gfortran'

# Use "update-alternatives" to switch compiler version.
# For example, if you have gcc-4.8 and gcc-6:
# update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 80
# update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-6   20
# update-alternatives --config gcc

#-------------------------------#
#   Read above up to this row   #
#-------------------------------#

# mpich folder structure
MPICH_DIR=$INSTALL_DIR/Mpich

# openmpi folder structure
OPENMPI_DIR=$INSTALL_DIR/Openmpi
# hdf5 seq folder structure
# hdf5 par based on mpich folder structure
HDF5_PAR_MPICH_DIR=$INSTALL_DIR/Hdf5_Par_Mpich

# hdf5 par based on mpich folder structure
HDF5_PAR_OPENMPI_DIR=$INSTALL_DIR/Hdf5_Par_Openmpi

# hdf5 par based on mpich folder structure
HDF5_SEQ_DIR=$INSTALL_DIR/Hdf5_Seq

# tcl folder structure
TCL_DIR=$INSTALL_DIR/Tcl

# tk folder structure
TK_DIR=$INSTALL_DIR/Tk

#current environment vars
PATH_CLEAN="$PATH"
LD_LIBRARY_PATH_CLEAN="$LD_LIBRARY_PATH"

if [ $BUILD_LATEST_CGNS == false ]; then
  # cgns_3.2.1+hdf5+cgnstools folder structure
  CGNS_HDF5_SEQ_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Seq

  # cgns_3.2.1+ parallel hdf5 based on mpich folder structure
  CGNS_HDF5_PAR_MPICH_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Mpich

  # cgns_3.2.1+ parallel hdf5 based on openmpi folder structure
  CGNS_HDF5_PAR_OPENMPI_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Par_Openmpi

else
  # cgns_Latest+hdf5+mpi folder structure
  CGNS_HDF5_PAR_MPICH_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Par

  # cgns_Latest+hdf5+cgnstools folder structure
  CGNS_HDF5_SEQ_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Seq
fi

mkdir -p $SRC_DIR/
mkdir -p $INSTALL_DIR/

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT

#------------------------------------------------------------------------------#
#   switch_compilers_and_environment_to                                        #
#------------------------------------------------------------------------------#
#   Sets compilers, optimization flags and environment                         #
#------------------------------------------------------------------------------#
function switch_compilers_and_environment_to {
  #$1 - "mpich", "openmpi", "sequential", seq_openmpi_only

  # renew environment
  unset PATH
  export PATH="$PATH_CLEAN"
  unset LD_LIBRARY_PATH
  export LD_LIBRARY_PATH="$LD_LIBRARY_PATH_CLEAN"

  if [ $1 = "mpich" ]; then
    export CCLAGS='-O2 -fPIC'
    export CPPFLAGS='-O2 -fPIC'
    export CXXFLAGS='-O2 -fPIC'
    export CFLAGS='-O2 -fPIC'
    export FFLAGS='-O2 -fPIC'
    export FCFLAGS='-O2 -fPIC'

    export CC=$MPICH_DIR/bin/mpicc
    export FC=$MPICH_DIR/bin/mpif90
    export CXX=$MPICH_DIR/bin/mpicxx
    export PATH=$MPICH_DIR/bin:$PATH
    export LD_LIBRARY_PATH=$MPICH_DIR/lib/:$LD_LIBRARY_PATH
  elif [ $1 = "openmpi" ]; then
    export CCLAGS='-O2 -fPIC'
    export CPPFLAGS='-O2 -fPIC'
    export CXXFLAGS='-O2 -fPIC'
    export CFLAGS='-O2 -fPIC'
    export FFLAGS='-O2 -fPIC'
    export FCFLAGS='-O2 -fPIC'

    export CC=$OPENMPI_DIR/bin/mpicc
    export FC=$OPENMPI_DIR/bin/mpif90
    export CXX=$OPENMPI_DIR/bin/mpicxx
    export PATH=$OPENMPI_DIR/bin:$PATH
    export LD_LIBRARY_PATH=$OPENMPI_DIR/lib/:$LD_LIBRARY_PATH
  elif [ $1 = "sequential" ]; then
    export CCLAGS='-O2'
    export CPPFLAGS='-O2'
    export CXXFLAGS='-O2'
    export CFLAGS='-O2'
    export FFLAGS='-O2'
    export FCFLAGS='-O2'

    export CC=$CCOMP
    export CXX=$CPPCOMP
    export FC=$FCOMP
  else
    exit 1
    echo "Error in switch_compilers_and_environment_to: wrong argument"
  fi
}

#------------------------------------------------------------------------------#
#   MPICH 3.2.1                                                                #
#------------------------------------------------------------------------------#
#   Builds mpich lib and mpicc, mpicxx, mpif90 compilers                       #
#------------------------------------------------------------------------------#
function build_mpich_lib {
  echo 'Mpich library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  wget -N http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz \
    >> $LOG_FILE 2>&1
  tar -zxf mpich-3.2.1.tar.gz
  mkdir -p $SRC_DIR/Mpich/
  rsync -azh mpich-3.2.1/* $SRC_DIR/Mpich/
  rm -r mpich-3.2.1/

  echo '  Setting compilers'
  switch_compilers_and_environment_to sequential

  echo '  Configuring installation'
  cd $SRC_DIR/Mpich/

  ./configure \
  --prefix=$MPICH_DIR \
  --enable-fast=all,O3 \
  --enable-fortran=all \
  >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  make install >> $LOG_FILE 2>&1

  echo '  Testing build'
  switch_compilers_and_environment_to mpich
  #make check-recursive >> $LOG_FILE 2>&1

  cd $CGNS_DIR
  echo '  Mpich lib was successfully built [dynamically]'
}
#------------------------------------------------------------------------------#
#   OpenMPI 2.0.2                                                              #
#------------------------------------------------------------------------------#
#   Builds openmpi lib and mpicc, mpicxx, mpif90 compilers                     #
#------------------------------------------------------------------------------#
function build_openmpi_lib {
  echo 'OpenMPI library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  wget -N https://download.open-mpi.org/release/open-mpi/v2.0/openmpi-2.0.2.tar.gz \
    >> $LOG_FILE 2>&1
  tar -zxf openmpi-2.0.2.tar.gz
  
  mkdir -p $SRC_DIR/Openmpi/
  rsync -azh openmpi-2.0.2/* $SRC_DIR/Openmpi/
  rm -r openmpi-2.0.2/

  echo '  Setting compilers'
  switch_compilers_and_environment_to sequential

  echo '  Configuring installation'
  cd $SRC_DIR/Openmpi/

  ./configure \
  --prefix=$OPENMPI_DIR \
  >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  make install >> $LOG_FILE 2>&1

  echo '  Testing build'
  switch_compilers_and_environment_to openmpi
  #make check-recursive >> $LOG_FILE 2>&1

  cd $CGNS_DIR
  echo '  Openmpi lib was successfully built [dynamically]'
}
#------------------------------------------------------------------------------#
#   HDF5 1.8.21                                                                #
#------------------------------------------------------------------------------#
#   Builds hdf5 lib in parallel() and 
#   (Paraview 5.4.1 & Visit 2.12.3 work with HDF5 5.1.8 and not with 1.10.*)   #
#------------------------------------------------------------------------------#
function build_hdf5_lib {
  echo 'HDF5 library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  mkdir -p $SRC_DIR/Hdf5/
  git clone --depth=1 https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git \
  --branch hdf5_1_8_21 ./hdf5 >> $LOG_FILE 2>&1
  rm -rf ./hdf5/.git
  rsync -azh hdf5/* $SRC_DIR/Hdf5/
  rm -rf ./hdf5

  if [ $BUILD_MPI == true ]; then
    # --------------------------------------------------------------------------
    echo '  Working with parallel version based on mpich:'
    echo '    Setting compilers'
    switch_compilers_and_environment_to mpich

    echo '    Configuring installation'
    cd $SRC_DIR/Hdf5/

    ./configure \
    --prefix=$HDF5_PAR_MPICH_DIR \
    --enable-fortran=no \
    --disable-hl \
    --enable-shared=no \
    --enable-production \
    --enable-parallel \
    --with-zlib=no \
    >> $LOG_FILE 2>&1

    echo '    Building'
    make lib >> $LOG_FILE 2>&1

    echo '    Installing'
    mkdir  -p $HDF5_PAR_MPICH_DIR/
    make install >> $LOG_FILE 2>&1

    echo '    Testing build'
    cd testpar/
    make test >> $LOG_FILE 2>&1
    # --------------------------------------------------------------------------
    echo '  Working with parallel version based on openmpi:'
    echo '    Setting compilers'
    switch_compilers_and_environment_to openmpi

    echo '    Configuring installation'
    cd $SRC_DIR/Hdf5/
    make clean >> $LOG_FILE 2>&1

    ./configure \
    --prefix=$HDF5_PAR_OPENMPI_DIR \
    --enable-fortran=no \
    --disable-hl \
    --enable-shared=no \
    --enable-production \
    --enable-parallel \
    --with-zlib=no \
    >> $LOG_FILE 2>&1

    echo '    Building'
    make lib >> $LOG_FILE 2>&1

    echo '    Installing'
    mkdir  -p $HDF5_PAR_OPENMPI_DIR/
    make install >> $LOG_FILE 2>&1

    echo '    Testing build [fails naturally on openmpi]'
    cd testpar/
    #make test >> $LOG_FILE 2>&1
    # --------------------------------------------------------------------------
  fi

  # ----------------------------------------------------------------------------
  echo '  Working with sequential version:'

  echo '    Setting compilers'
  switch_compilers_and_environment_to sequential

  echo '    Configuring installation'
  cd $SRC_DIR/Hdf5/; rm -rf .git
  if [ $BUILD_MPI == true ]; then
    make clean >> $LOG_FILE 2>&1
  fi

  ./configure \
  --prefix=$HDF5_SEQ_DIR \
  --enable-fortran=no \
  --disable-hl \
  --enable-shared=no \
  --enable-production \
  --with-zlib=no \
  >> $LOG_FILE 2>&1

  echo '    Building'
  make lib >> $LOG_FILE 2>&1

  echo '    Installing'
  mkdir  -p $HDF5_SEQ_DIR/
  make install >> $LOG_FILE 2>&1

  echo '    Testing build'
  cd test/
  make test >> $LOG_FILE 2>&1
  # ----------------------------------------------------------------------------

  cd $CGNS_DIR
  echo '  Hdf5 lib was successfully built [statically]'
}
#------------------------------------------------------------------------------#
#   TCL 8.6.8                                                                  #
#------------------------------------------------------------------------------#
function build_tcl_lib {
  echo 'Tcl library:'

  cd $SRC_DIR

  echo '  Downloading sources'
  cd $SRC_DIR
  wget -N https://prdownloads.sourceforge.net/tcl/tcl8.6.8-src.tar.gz \
    >> $LOG_FILE 2>&1
  tar -zxf tcl8.6.8-src.tar.gz;
  mkdir -p $SRC_DIR/Tcl/
  rsync -azh tcl8.6.8/* $SRC_DIR/Tcl/
  rm -r tcl8.6.8/

  echo '  Setting compilers'
  switch_compilers_and_environment_to sequential

  echo '  Configuring installation'
  cd Tcl/unix/
  ./configure \
  --prefix=$TCL_DIR \
  >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  mkdir  -p $TCL_DIR/
  find $TCL_DIR/ -maxdepth 1 -type l -exec unlink {} \;
  make install >> $LOG_FILE 2>&1

  echo '  Making links and Lib_Config'
  ln -srf $TCL_DIR/lib            $TCL_DIR/unix
  ln -srf $TCL_DIR/include        $TCL_DIR/generic

  cd $CGNS_DIR
  echo '  Tcl lib was successfully built'
}
#------------------------------------------------------------------------------#
#   TK 8.6.8                                                                   #
#------------------------------------------------------------------------------#
#   Requires libx11-dev from repo                                              #
#------------------------------------------------------------------------------#
function build_tk_lib {
  echo 'Tk library:'

  cd $SRC_DIR

  echo '  Downloading sources'
  wget -N https://prdownloads.sourceforge.net/tcl/tk8.6.8-src.tar.gz \
    >> $LOG_FILE 2>&1
  tar -zxf tk8.6.8-src.tar.gz;
  mkdir -p $SRC_DIR/Tk/
  rsync -azh tk8.6.8/* $SRC_DIR/Tk/
  rm -r tk8.6.8/

  echo '  Setting compilers'
  switch_compilers_and_environment_to sequential

  echo '  Configuring installation'
  cd Tk/unix/
  ./configure \
  --prefix=$TK_DIR \
  --with-tcl=$TCL_DIR/lib \
    >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  mkdir -p $TK_DIR/
  find $TK_DIR/ -maxdepth 1 -type l -exec unlink {} \;
  make install >> $LOG_FILE 2>&1

  echo '  Making links and Lib_Config'
  ln -srf $TK_DIR/lib            $TK_DIR/unix
  ln -srf $TK_DIR/include        $TK_DIR/generic
  ln -srf $TK_DIR/lib/tk8.6/     $TK_DIR/library

  cd $CGNS_DIR
  echo '  Tk lib was successfully built'
}
#------------------------------------------------------------------------------#
#   CGNS latest                                                                #
#------------------------------------------------------------------------------#
#   Does not support 'include ' module fortran syntax, .f77 and f90.           #
#------------------------------------------------------------------------------#
function build_cgns_lib_latest {

  echo 'Cgns latest library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  mkdir -p $SRC_DIR/Cgns_Latest/
  git clone --depth=1 https://github.com/CGNS/CGNS.git ./cgns >> $LOG_FILE 2>&1
  rsync -azh cgns/* $SRC_DIR/Cgns_Latest/
  rm -rf cgns
  cd $SRC_DIR/Cgns_Latest/; rm -rf .git; cd src/

  if [ $BUILD_MPI == true ]; then
    # --------------------------------------------------------------------------
    echo '  Working with parallel version + HDF5:'

    echo '    Setting compilers'
    switch_compilers_and_environment_to mpich

    echo '    Configuring installation'
    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_PAR_MPICH_DIR \
    --with-hdf5=$HDF5_PAR_MPICH_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-parallel \
    --with-mpi=$MPICH_DIR/bin \
    >> $LOG_FILE 2>&1

    echo '    Building'
    make >> $LOG_FILE 2>&1

    echo '    Installing'
    make install >> $LOG_FILE 2>&1

  fi

  # ----------------------------------------------------------------------------
  echo '    Setting compilers'
  switch_compilers_and_environment_to sequential

  if [ $BUILD_CGNS_TOOLS == true ]; then
    echo '  Working with sequential version + HDF5 + cgnstools:'

    echo '    Configuring installation'
    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_DIR \
    --with-hdf5=$HDF5_SEQ_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-cgnstools \
    --with-tcl=$TCL_DIR \
    --with-tk=$TK_DIR \
    >> $LOG_FILE 2>&1

  else # no cgns GUI tools
    echo '  Working with sequential version + HDF5:'

    echo '    Configuring installation'
    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_DIR \
    --with-hdf5=$HDF5_SEQ_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    >> $LOG_FILE 2>&1
  fi

  echo '    Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  cd $CGNS_DIR
}
#------------------------------------------------------------------------------#
#   CGNS 3.2.1                                                                 #
#------------------------------------------------------------------------------#
#   Supports 'include ' module fortran syntax, .f77 and f90.                   #
#------------------------------------------------------------------------------#
function build_cgns_lib_3.2.1 {
  echo 'Cgns 3.2.1 library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  wget -N https://github.com/CGNS/CGNS/archive/v3.2.1.tar.gz >> $LOG_FILE 2>&1
  tar -zxf v3.2.1.tar.gz
  mkdir -p $SRC_DIR/Cgns_3.2.1/
  rsync -azh CGNS-3.2.1/* $SRC_DIR/Cgns_3.2.1/
  rm -r CGNS-3.2.1/
  cd $SRC_DIR/Cgns_3.2.1/; rm -rf .git; cd src/

  if [ $BUILD_MPI == true ]; then
    # --------------------------------------------------------------------------
    echo '  Working with parallel version based on mpich + HDF5:'

    echo '    Setting compilers'
    switch_compilers_and_environment_to mpich

    echo '    Configuring installation'

    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_PAR_MPICH_DIR \
    --with-hdf5=$HDF5_PAR_MPICH_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-parallel \
    --with-mpi=$MPICH_DIR/bin \
    >> $LOG_FILE 2>&1

    # fix a bug in 3.2.1 makefile[missing -l mpi]
    sed -i -e "s%-L$MPICH_DIR/lib\ -l%-L$MPICH_DIR/lib\ -l mpi%g" make.defs
    sed -i -e "s%$MPICH_DIR/bin/mpiexec%$MPICH_DIR/bin/mpiexec -n \$\${NPROCS:=4}%g" make.defs

    echo '    Building'
    make lib >> $LOG_FILE 2>&1

    echo '    Installing'
    make install >> $LOG_FILE 2>&1

    echo '  Testing build'
    switch_compilers_and_environment_to mpich
    #make ptests >> $LOG_FILE 2>&1

    # --------------------------------------------------------------------------
    echo '  Working with parallel version based on openmpi + HDF5:'

    echo '    Setting compilers'
    switch_compilers_and_environment_to openmpi

    echo '    Configuring installation'
    make clean >> $LOG_FILE 2>&1

    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_PAR_OPENMPI_DIR \
    --with-hdf5=$HDF5_PAR_OPENMPI_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-parallel \
    --with-mpi=$OPENMPI_DIR/bin \
    >> $LOG_FILE 2>&1

    # fix a bug in 3.2.1 makefile[missing -l mpi]
    sed -i -e "s%-L$OPENMPI_DIR/lib\ -l%-L$OPENMPI_DIR/lib\ -l mpi%g" make.defs
    sed -i -e "s%$OPENMPI_DIR/bin/mpiexec%$OPENMPI_DIR/bin/mpiexec -n \$\${NPROCS:=4}%g" make.defs

    echo '    Building'
    make lib >> $LOG_FILE 2>&1

    echo '    Installing'
    make install >> $LOG_FILE 2>&1

    echo '  Testing build'
    switch_compilers_and_environment_to openmpi
    #make ptests >> $LOG_FILE 2>&1

  fi

  # ----------------------------------------------------------------------------
  echo '    Setting compilers'
  switch_compilers_and_environment_to sequential

  if [ $BUILD_CGNS_TOOLS == true ]; then
    echo '  Working with sequential version + HDF5 + cgnstools:'

    echo '    Configuring installation'

    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_DIR \
    --with-hdf5=$HDF5_SEQ_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-cgnstools \
    --with-tcl=$TCL_DIR \
    --with-tk=$TK_DIR \
    >> $LOG_FILE 2>&1

  else # no cgns GUI tools
    # --------------------------------------------------------------------------
    echo '  Working with sequential version + HDF5:'

    echo '    Configuring installation'
    if [ $BUILD_MPI == true ]; then
      make clean >> $LOG_FILE 2>&1
    fi

    FLIBS="-Wl,--no-as-needed -ldl" \
    LIBS="-Wl,--no-as-needed -ldl" \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_DIR \
    --with-hdf5=$HDF5_SEQ_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    >> $LOG_FILE 2>&1
  fi

  echo '    Building'
  make lib >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  echo '  Testing build'
  switch_compilers_and_environment_to sequential
  #make tests >> $LOG_FILE 2>&1

  cd $CGNS_DIR
  echo '  Cgns 3.2.1 library was successfully built [statically]'
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

if [ $BUILD_MPI == true ]; then
  build_mpich_lib
  build_openmpi_lib
fi

build_hdf5_lib

if [ $BUILD_CGNS_TOOLS == true ]; then
  build_tcl_lib
  build_tk_lib
fi

if [ $BUILD_LATEST_CGNS == true ]; then
  build_cgns_lib_latest
else
  build_cgns_lib_3.2.1
fi

remove_temp_dir
switch_compilers_and_environment_to sequential
echo done
