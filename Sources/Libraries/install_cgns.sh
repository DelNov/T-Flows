#!/bin/bash
# Requires: C, C++, Fortran compilers
# Requires: libx11-dev from debian/ubuntu repository if BUILD_CGNS_TOOLS=true

# This script will automatically build cgns lib in several releases:
# 1) with adf5 and sequential access
# 2) with hdf5 and parallel   access (optional)
# 3) with hdf5 and sequential access + cgnstools (optional)
# (parallel cgns lib can only be based on hdf5, not adf5)
#
# Some grid meshing software(like pointwise) use hdf5, therefore it is not
# possible to use Convert on mesh with adf5

# Hdf5_Par lib most likely can not be built with your current mpif90
# therefore script builds mpif90 first to compile Hdf5_Par later

# you can optionally build cgnstools to view/edit .cgns files
# https://cgns.github.io/CGNS_docs_current/cgnstools/cgnsview/index.html
# if this case script will download and compile dependencies: tcl, tk libs

# T-Flows supports currently only 3.2.1 cgns lib, but
# you can build with this script the latest cgns lib if you wish
# since 3.3 version fortran 77 and 90 support was dropped,
# moreover user is forced to use cgsize_t type instead of integer in cgns API

# CGNS lib does not use compression, but script still builds zlib:
# https://cgns.github.io/FAQs.html

# build mpich, otherwise Hdf5_Par can not be built
BUILD_MPI=true

# build cgnstools, which contains cgnsview - small program to view .cgns files
# https://cgns.github.io/CGNS_docs_current/cgnstools/cgnsview/index.html
BUILD_CGNS_TOOLS=false

# build latest cgns lib, otherwise build 3.2.1 version
BUILD_LATEST_CGNS=false

# folder structure
CGNS_DIR=$PWD                     # this is top dir for this lib
INSTALL_DIR=$CGNS_DIR/install_dir # this dir contains Cgns/ Hdf5/ Mpich/
SRC_DIR=$CGNS_DIR/src_dir         # this dir contains downloaded sources

# script temp dir
TEMP_DIR=`mktemp -d 2>/dev/null || mktemp -d -t 'TEMP_DIR'`
# script logs
LOG_FILE=$SRC_DIR/build_cgns.log # logs of current script
cp /dev/null $LOG_FILE

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

# optimization flags
export CCLAGS='-O3'
export FFLAGS='-O3'
export FCFLAGS='-O3'
export CFLAGS='-O3'
export CXXFLAGS='-O3'
export CPPFLAGS='-O3'

# zlib folder structure
ZLIB_ROOT_DIR=$INSTALL_DIR/Zlib
ZLIB_LIB_DIR=$INSTALL_DIR/Zlib/Lib
ZLIB_INCLUDE_DIR=$INSTALL_DIR/Zlib/Include

# mpich folder structure
MPICH_ROOT_DIR=$INSTALL_DIR/Mpich
MPICH_BIN_DIR=$INSTALL_DIR/Mpich/Bin
MPICH_LIB_DIR=$INSTALL_DIR/Mpich/Lib
MPICH_INCLUDE_DIR=$INSTALL_DIR/Mpich/Include

# hdf5 seq folder structure
HDF5_SEQ_ROOT_DIR=$INSTALL_DIR/Hdf5_Seq
HDF5_SEQ_BIN_DIR=$INSTALL_DIR/Hdf5_Seq/Bin
HDF5_SEQ_LIB_DIR=$INSTALL_DIR/Hdf5_Seq/Lib
HDF5_SEQ_INCLUDE_DIR=$INSTALL_DIR/Hdf5_Seq/Include

# hdf5 par folder structure
HDF5_PAR_ROOT_DIR=$INSTALL_DIR/Hdf5_Par
HDF5_PAR_BIN_DIR=$INSTALL_DIR/Hdf5_Par/Bin
HDF5_PAR_LIB_DIR=$INSTALL_DIR/Hdf5_Par/Lib
HDF5_PAR_INCLUDE_DIR=$INSTALL_DIR/Hdf5_Par/Include

# tcl folder structure
TCL_ROOT_DIR=$INSTALL_DIR/Tcl
TCL_BIN_DIR=$INSTALL_DIR/Tcl/Bin
TCL_LIB_DIR=$INSTALL_DIR/Tcl/Lib
TCL_INCLUDE_DIR=$INSTALL_DIR/Tcl/Include

# tk folder structure
TK_ROOT_DIR=$INSTALL_DIR/Tk
TK_BIN_DIR=$INSTALL_DIR/Tk/Bin
TK_LIB_DIR=$INSTALL_DIR/Tk/Lib
TK_INCLUDE_DIR=$INSTALL_DIR/Tk/Include


if [ $BUILD_LATEST_CGNS == false ]; then

  # cgns_3.2.1+hdf5+mpi folder structure
  CGNS_HDF5_PAR_ROOT_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Par
  CGNS_HDF5_PAR_BIN_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Bin
  CGNS_HDF5_PAR_LIB_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Lib
  CGNS_HDF5_PAR_INCLUDE_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Par/Include

  # cgns_3.2.1+hdf5+cgnstools folder structure
  CGNS_HDF5_SEQ_ROOT_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Seq
  CGNS_HDF5_SEQ_BIN_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Bin
  CGNS_HDF5_SEQ_LIB_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Lib
  CGNS_HDF5_SEQ_INCLUDE_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Hdf5_Seq/Include

  # cgns_3.2.1+adf5 folder structure
  CGNS_ADF5_SEQ_ROOT_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Adf5_Seq
  CGNS_ADF5_SEQ_BIN_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Bin
  CGNS_ADF5_SEQ_LIB_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Lib
  CGNS_ADF5_SEQ_INCLUDE_DIR=$INSTALL_DIR/Cgnslib_3.2.1_Linux_64_Adf5_Seq/Include
else
  # cgns_Latest+hdf5+mpi folder structure
  CGNS_HDF5_PAR_ROOT_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Par
  CGNS_HDF5_PAR_BIN_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Par/Bin
  CGNS_HDF5_PAR_LIB_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Par/Lib
  CGNS_HDF5_PAR_INCLUDE_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Par/Include

  # cgns_Latest+hdf5+cgnstools folder structure
  CGNS_HDF5_SEQ_ROOT_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Seq
  CGNS_HDF5_SEQ_BIN_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Seq/Bin
  CGNS_HDF5_SEQ_LIB_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Seq/Lib
  CGNS_HDF5_SEQ_INCLUDE_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Hdf5_Seq/Include

  # cgns_Latest+adf5 folder structure
  CGNS_ADF5_SEQ_ROOT_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Adf5_Seq
  CGNS_ADF5_SEQ_BIN_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Adf5_Seq/Bin
  CGNS_ADF5_SEQ_LIB_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Adf5_Seq/Lib
  CGNS_ADF5_SEQ_INCLUDE_DIR=$INSTALL_DIR/Cgnslib_Latest_Linux_64_Adf5_Seq/Include
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
#   ZLIB 1.2.11                                                                #
#------------------------------------------------------------------------------#
#   Builds zlib                                                                #
#------------------------------------------------------------------------------#
function build_z_lib {
  echo 'Zlib library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  wget -N http://zlib.net/zlib-1.2.11.tar.gz >> $LOG_FILE 2>&1
  tar -zxf zlib-1.2.11.tar.gz
  mkdir -p $SRC_DIR/Zlib/
  rsync -azh zlib-1.2.11/* $SRC_DIR/Zlib/
  rm -r zlib-1.2.11/

  echo '    Setting compilers'
  export CC=$CCOMP
  export CXX=$CPPCOMP
  export FC=$FCOMP

  echo '  Configuring installation'
  cd $SRC_DIR/Zlib/
  ./configure \
  --prefix=$ZLIB_ROOT_DIR \
  --libdir=$ZLIB_LIB_DIR \
  --includedir=$ZLIB_INCLUDE_DIR \
  --sharedlibdir=$TEMP_DIR \
  --64 \
  --static \
  >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  mkdir  -p $ZLIB_ROOT_DIR/
  find $ZLIB_ROOT_DIR/ -maxdepth 1 -type l -exec unlink {} \;
  make install >> $LOG_FILE 2>&1

  echo '  Making links and Lib_Config'
  mv $ZLIB_LIB_DIR/pkgconfig/zlib.pc $ZLIB_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g"        $ZLIB_LIB_DIR/Lib_Config
  ln -srf $ZLIB_INCLUDE_DIR/zlib.h   $ZLIB_ROOT_DIR/zlib.h
  ln -srf $ZLIB_LIB_DIR/libz.a       $ZLIB_ROOT_DIR/libz.a
  ln -srf $ZLIB_INCLUDE_DIR          $ZLIB_ROOT_DIR/include
  ln -srf $ZLIB_LIB_DIR              $ZLIB_ROOT_DIR/lib
  rm -r $ZLIB_LIB_DIR/pkgconfig/     $ZLIB_ROOT_DIR/share/

  cd $CGNS_DIR
  echo '  Zlib lib was successfully built'  
}
#------------------------------------------------------------------------------#
#   MPICH 3.2.1                                                                #
#------------------------------------------------------------------------------#
#   Builds mpich lib and mpicc, mpicxx, mpif90 compilers                       #
#------------------------------------------------------------------------------#
function build_mpi_lib {
  echo 'Mpich library:'

  echo '  Downloading sources'
  cd $SRC_DIR
  wget -N http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz \
    >> $LOG_FILE 2>&1
  tar -zxf mpich-3.2.1.tar.gz
  mkdir -p $SRC_DIR/Mpich/
  rsync -azh mpich-3.2.1/* $SRC_DIR/Mpich/
  rm -r mpich-3.2.1/

  echo '  Configuring installation'
  cd $SRC_DIR/Mpich/
  ./configure \
  --prefix=$MPICH_ROOT_DIR \
  --bindir=$MPICH_BIN_DIR \
  --libdir=$MPICH_LIB_DIR \
  --includedir=$MPICH_INCLUDE_DIR \
  --datarootdir=$TEMP_DIR \
  CC=$CCOMP CXX=$CPPCOMP F77=$FCOMP FC=$FCOMP \
  --enable-fast=all,O3 \
  --enable-fortran=all \
  --disable-shared \
  >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  echo '  Making links and Lib_Config'
  mv $MPICH_LIB_DIR/pkgconfig/mpich.pc $MPICH_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g"          $MPICH_LIB_DIR/Lib_Config
  rm -r $MPICH_LIB_DIR/pkgconfig/

  cd $CGNS_DIR
  echo '  Mpich lib was successfully built'  
}
#------------------------------------------------------------------------------#
#   HDF5 1.8.20                                                                #
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
  --branch hdf5_1_8_20 ./hdf5 >> $LOG_FILE 2>&1
  rm -rf ./hdf5/.git
  rsync -azh hdf5/* $SRC_DIR/Hdf5/
  rm -rf ./hdf5

  if [ $BUILD_MPI == true ]; then
    echo '  Working with parallel version:'
    echo '    Setting compilers'
    export CC=$MPICH_BIN_DIR/mpicc
    export FC=$MPICH_BIN_DIR/mpif90
    export CXX=$MPICH_BIN_DIR/mpicxx

    echo '    Configuring installation'
    cd $SRC_DIR/Hdf5/
    ./configure \
    --prefix=$HDF5_PAR_ROOT_DIR \
    --bindir=$HDF5_PAR_BIN_DIR \
    --libdir=$HDF5_PAR_LIB_DIR \
    --includedir=$HDF5_PAR_INCLUDE_DIR \
    --datarootdir=$TEMP_DIR \
    --enable-fortran=yes \
    --enable-shared=no \
    --enable-production \
    --with-zlib=$ZLIB_LIB_DIR \
    --enable-parallel \
    >> $LOG_FILE 2>&1

    echo '    Building'
    make >> $LOG_FILE 2>&1

    echo '    Installing'
    mkdir  -p $HDF5_PAR_ROOT_DIR/
    find $HDF5_PAR_ROOT_DIR/ -maxdepth 1 -type l -exec unlink {} \;
    make install >> $LOG_FILE 2>&1

    echo '    Making links and Lib_Config'
    rm -r $HDF5_PAR_ROOT_DIR/share/
    mv $HDF5_PAR_LIB_DIR/libhdf5.settings $HDF5_PAR_LIB_DIR/Lib_Config
    sed -i -e "s%$CGNS_DIR/%%g"           $HDF5_PAR_LIB_DIR/Lib_Config
    ln -srf $HDF5_PAR_LIB_DIR             $HDF5_PAR_ROOT_DIR/lib
    ln -srf $HDF5_PAR_INCLUDE_DIR         $HDF5_PAR_ROOT_DIR/include
  fi

  echo '  Working with sequential version:'

  echo '    Setting compilers'
  export CC=$CCOMP
  export CXX=$CPPCOMP
  export FC=$FCOMP

  echo '    Configuring installation'
  cd $SRC_DIR/Hdf5/; rm -rf .git
  ./configure \
  --prefix=$HDF5_SEQ_ROOT_DIR \
  --bindir=$HDF5_SEQ_BIN_DIR \
  --libdir=$HDF5_SEQ_LIB_DIR \
  --includedir=$HDF5_SEQ_INCLUDE_DIR \
  --datarootdir=$TEMP_DIR \
  --enable-fortran=yes \
  --enable-shared=no \
  --enable-production \
  --with-zlib=$ZLIB_LIB_DIR \
  >> $LOG_FILE 2>&1

  echo '    Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  mkdir  -p $HDF5_SEQ_ROOT_DIR/
  find $HDF5_SEQ_ROOT_DIR/ -maxdepth 1 -type l -exec unlink {} \;
  make install >> $LOG_FILE 2>&1

  echo '    Making links and Lib_Config'
  rm -r $HDF5_SEQ_ROOT_DIR/share/
  mv $HDF5_SEQ_LIB_DIR/libhdf5.settings $HDF5_SEQ_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g"           $HDF5_SEQ_LIB_DIR/Lib_Config
  ln -srf $HDF5_SEQ_LIB_DIR             $HDF5_SEQ_ROOT_DIR/lib
  ln -srf $HDF5_SEQ_INCLUDE_DIR         $HDF5_SEQ_ROOT_DIR/include

  cd $CGNS_DIR
  echo '  Hdf5 lib was successfully built'
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
  export CC=$CCOMP
  export CXX=$CPPCOMP
  export FC=$FCOMP

  echo '  Configuring installation'
  cd Tcl/unix/
  ./configure \
  --prefix=$TCL_ROOT_DIR \
  --bindir=$TCL_BIN_DIR \
  --includedir=$TCL_INCLUDE_DIR \
  --datadir=$TEMP_DIR \
  --mandir=$TEMP_DIR \
  >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  mkdir  -p $TCL_ROOT_DIR/
  find $TCL_ROOT_DIR/ -maxdepth 1 -type l -exec unlink {} \;
  make install >> $LOG_FILE 2>&1

  echo '  Making links and Lib_Config'
  if [ -d "$TCL_LIB_DIR" ]; then
    rm -rf $TCL_LIB_DIR
  fi
  mv      $TCL_ROOT_DIR/lib/       $TCL_LIB_DIR
  ln -srf $TCL_LIB_DIR             $TCL_ROOT_DIR/lib
  ln -srf $TCL_LIB_DIR             $TCL_ROOT_DIR/unix
  ln -srf $TCL_INCLUDE_DIR         $TCL_ROOT_DIR/generic

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
  export CC=$CCOMP
  export CXX=$CPPCOMP
  export FC=$FCOMP

  echo '  Configuring installation'
  cd Tk/unix/
  ./configure \
  --prefix=$TK_ROOT_DIR \
  --bindir=$TK_BIN_DIR \
  --includedir=$TK_INCLUDE_DIR \
  --datadir=$TEMP_DIR \
  --mandir=$TEMP_DIR \
  --with-tcl=$TCL_LIB_DIR \
    >> $LOG_FILE 2>&1

  echo '  Building'
  make >> $LOG_FILE 2>&1

  echo '  Installing'
  mkdir -p $TK_ROOT_DIR/
  find $TK_ROOT_DIR/ -maxdepth 1 -type l -exec unlink {} \;
  make install >> $LOG_FILE 2>&1

  echo '  Making links and Lib_Config'
  if [ -d "$TK_LIB_DIR" ]; then
    rm -rf $TK_LIB_DIR
  fi
  mv      $TK_ROOT_DIR/lib/      $TK_LIB_DIR
  ln -srf $TK_LIB_DIR            $TK_ROOT_DIR/lib
  ln -srf $TK_LIB_DIR            $TK_ROOT_DIR/unix
  ln -srf $TK_INCLUDE_DIR        $TK_ROOT_DIR/generic
  ln -srf $TK_LIB_DIR/tk8.6/     $TK_ROOT_DIR/library

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

  # fix a bug in 3.3.1 configure
  sed -i -e "s%unset HAVE_ZLIB%echo intervened%g" configure

  if [ $BUILD_MPI == true ]; then
    echo '  Working with parallel version + HDF5:'

    echo '    Setting compilers'
    export CC=$MPICH_BIN_DIR/mpicc
    export FC=$MPICH_BIN_DIR/mpif90
    export CXX=$MPICH_BIN_DIR/mpicxx

    echo '    Configuring installation'
    FLIBS=-Wl,--no-as-needed\ -ldl \
    LIBS=-Wl,--no-as-needed\ -ldl \
    ./configure \
    --prefix=$CGNS_HDF5_PAR_ROOT_DIR \
    --bindir=$CGNS_HDF5_PAR_BIN_DIR \
    --libdir=$CGNS_HDF5_PAR_LIB_DIR \
    --includedir=$CGNS_HDF5_PAR_INCLUDE_DIR \
    --datadir=$TEMP_DIR \
    --mandir=$TEMP_DIR \
    --with-hdf5=$HDF5_PAR_ROOT_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-parallel \
    --with-mpi=$MPICH_BIN_DIR \
    --with-zlib=$ZLIB_ROOT_DIR \
    >> $LOG_FILE 2>&1

    echo '    Building'
    make >> $LOG_FILE 2>&1

    echo '    Installing'
    make install >> $LOG_FILE 2>&1

    echo '    Making links and Lib_Config'
    mv $CGNS_HDF5_PAR_INCLUDE_DIR/cgnsBuild.defs \
      $CGNS_HDF5_PAR_LIB_DIR/Lib_Config
    sed -i -e "s%$CGNS_DIR/%%g" $CGNS_HDF5_PAR_LIB_DIR/Lib_Config
    rm $CGNS_HDF5_PAR_INCLUDE_DIR/cgnsconfig.h
  fi

  echo '    Setting compilers'
  export CC=$CCOMP
  export CXX=$CPPCOMP
  export FC=$FCOMP

  if [ $BUILD_CGNS_TOOLS == true ]; then

    echo '  Working with sequential version + HDF5 + cgnstools:'

    echo '    Configuring installation'
    FLIBS=-Wl,--no-as-needed\ -ldl \
    LIBS=-Wl,--no-as-needed\ -ldl \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_ROOT_DIR \
    --bindir=$CGNS_HDF5_SEQ_BIN_DIR \
    --libdir=$CGNS_HDF5_SEQ_LIB_DIR \
    --includedir=$CGNS_HDF5_SEQ_INCLUDE_DIR \
    --datadir=$TCL_LIB_DIR/Tcl_Scripts \
    --mandir=$TEMP_DIR \
    --with-hdf5=$HDF5_SEQ_ROOT_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-cgnstools \
    --with-tcl=$TCL_ROOT_DIR \
    --with-tk=$TK_ROOT_DIR \
    --with-zlib=$ZLIB_ROOT_DIR \
    >> $LOG_FILE 2>&1

  else # no cgns GUI tools
    echo '  Working with sequential version + HDF5:'

    echo '    Configuring installation'
    FLIBS=-Wl,--no-as-needed\ -ldl \
    LIBS=-Wl,--no-as-needed\ -ldl \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_ROOT_DIR \
    --bindir=$CGNS_HDF5_SEQ_BIN_DIR \
    --libdir=$CGNS_HDF5_SEQ_LIB_DIR \
    --includedir=$CGNS_HDF5_SEQ_INCLUDE_DIR \
    --datadir=$TEMP_DIR \
    --mandir=$TEMP_DIR \
    --with-hdf5=$HDF5_SEQ_ROOT_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --with-zlib=$ZLIB_ROOT_DIR \
    >> $LOG_FILE 2>&1
  fi

  echo '    Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  echo '    Making links and Lib_Config'
  mv $CGNS_HDF5_SEQ_INCLUDE_DIR/cgnsBuild.defs \
    $CGNS_HDF5_SEQ_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g" $CGNS_HDF5_SEQ_LIB_DIR/Lib_Config
  rm $CGNS_HDF5_SEQ_INCLUDE_DIR/cgnsconfig.h

  echo '  Working with sequential version + ADF5:'
  echo '    Configuring installation'
  FLIBS=-Wl,--no-as-needed\ -ldl \
  LIBS=-Wl,--no-as-needed\ -ldl \
  ./configure \
  --prefix=$CGNS_ADF5_SEQ_ROOT_DIR \
  --bindir=$CGNS_ADF5_SEQ_BIN_DIR \
  --libdir=$CGNS_ADF5_SEQ_LIB_DIR \
  --includedir=$CGNS_ADF5_SEQ_INCLUDE_DIR \
  --datadir=$TEMP_DIR \
  --mandir=$TEMP_DIR \
  --with-fortran \
  --enable-lfs \
  --enable-64bit \
  --disable-shared \
  --with-zlib=$ZLIB_ROOT_DIR \
  >> $LOG_FILE 2>&1

  echo '    Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  echo '    Making links and Lib_Config'
  mv $CGNS_ADF5_SEQ_INCLUDE_DIR/cgnsBuild.defs \
    $CGNS_ADF5_SEQ_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g" $CGNS_ADF5_SEQ_LIB_DIR/Lib_Config
  rm $CGNS_ADF5_SEQ_INCLUDE_DIR/cgnsconfig.h

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
    echo '  Working with parallel version + HDF5:'

    echo '    Setting compilers'
    export CC=$MPICH_BIN_DIR/mpicc
    export FC=$MPICH_BIN_DIR/mpif90
    export CXX=$MPICH_BIN_DIR/mpicxx

    echo '    Configuring installation'
    FLIBS=-Wl,--no-as-needed\ -ldl \
    LIBS=-Wl,--no-as-needed\ -ldl \
    ./configure \
    --prefix=$CGNS_HDF5_PAR_ROOT_DIR \
    --bindir=$CGNS_HDF5_PAR_BIN_DIR \
    --libdir=$CGNS_HDF5_PAR_LIB_DIR \
    --includedir=$CGNS_HDF5_PAR_INCLUDE_DIR \
    --datadir=$TEMP_DIR \
    --mandir=$TEMP_DIR \
    --with-hdf5=$HDF5_PAR_ROOT_DIR \
    --with-fortran \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-parallel \
    --with-mpi=$MPICH_BIN_DIR \
    --with-zlib=$ZLIB_LIB_DIR \
    >> $LOG_FILE 2>&1

    # fix a bug in 3.2.1 makefile
    sed -i -e "s%-L$MPICH_LIB_DIR\ -l%-L$MPICH_LIB_DIR\ -l mpi%g" make.defs
    sed -i -e "s%$MPICH_BIN_DIR/mpiexec%$MPICH_BIN_DIR/mpiexec  -n \$\${NPROCS:=4}%g" make.defs

    echo '    Building'
    make >> $LOG_FILE 2>&1

    echo '    Installing'
    make install >> $LOG_FILE 2>&1

    echo '    Making links and Lib_Config'
    mv $CGNS_HDF5_PAR_INCLUDE_DIR/cgnsBuild.defs \
      $CGNS_HDF5_PAR_LIB_DIR/Lib_Config
    sed -i -e "s%$CGNS_DIR/%%g" $CGNS_HDF5_PAR_LIB_DIR/Lib_Config
    rm $CGNS_HDF5_PAR_INCLUDE_DIR/cgnsconfig.h
  fi

  echo '    Setting compilers'
  export CC=$CCOMP
  export CXX=$CPPCOMP
  export FC=$FCOMP

  if [ $BUILD_CGNS_TOOLS == true ]; then
    echo '  Working with sequential version + HDF5 + cgnstools:'

    echo '    Configuring installation'
    FLIBS=-Wl,--no-as-needed\ -ldl \
    LIBS=-Wl,--no-as-needed\ -ldl \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_ROOT_DIR \
    --bindir=$CGNS_HDF5_SEQ_BIN_DIR \
    --libdir=$CGNS_HDF5_SEQ_LIB_DIR \
    --includedir=$CGNS_HDF5_SEQ_INCLUDE_DIR \
    --datadir=$TCL_LIB_DIR/Tcl_Scripts \
    --mandir=$TEMP_DIR \
    --with-hdf5=$HDF5_SEQ_ROOT_DIR \
    --with-fortran \
    --enable-gcc \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --enable-cgnstools \
    --with-tcl=$TCL_ROOT_DIR \
    --with-tk=$TK_ROOT_DIR \
    --with-zlib=$ZLIB_LIB_DIR \
    >> $LOG_FILE 2>&1

  else # no cgns GUI tools
    echo '  Working with sequential version + HDF5:'

    echo '    Configuring installation'
    FLIBS=-Wl,--no-as-needed\ -ldl \
    LIBS=-Wl,--no-as-needed\ -ldl \
    ./configure \
    --prefix=$CGNS_HDF5_SEQ_ROOT_DIR \
    --bindir=$CGNS_HDF5_SEQ_BIN_DIR \
    --libdir=$CGNS_HDF5_SEQ_LIB_DIR \
    --includedir=$CGNS_HDF5_SEQ_INCLUDE_DIR \
    --datadir=$TEMP_DIR \
    --mandir=$TEMP_DIR \
    --with-hdf5=$HDF5_SEQ_ROOT_DIR \
    --with-fortran \
    --enable-gcc \
    --enable-lfs \
    --enable-64bit \
    --disable-shared \
    --with-zlib=$ZLIB_LIB_DIR \
    >> $LOG_FILE 2>&1
  fi

  echo '    Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  echo '    Making links and Lib_Config'
  mv $CGNS_HDF5_SEQ_INCLUDE_DIR/cgnsBuild.defs \
    $CGNS_HDF5_SEQ_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g" $CGNS_HDF5_SEQ_LIB_DIR/Lib_Config
  rm $CGNS_HDF5_SEQ_INCLUDE_DIR/cgnsconfig.h

  echo '    Working with sequential version + ADF5:'
  echo '    Configuring installation'
  FLIBS=-Wl,--no-as-needed\ -ldl \
  LIBS=-Wl,--no-as-needed\ -ldl \
  ./configure \
  --prefix=$CGNS_ADF5_SEQ_ROOT_DIR \
  --bindir=$CGNS_ADF5_SEQ_BIN_DIR \
  --libdir=$CGNS_ADF5_SEQ_LIB_DIR \
  --includedir=$CGNS_ADF5_SEQ_INCLUDE_DIR \
  --datadir=$TEMP_DIR \
  --mandir=$TEMP_DIR \
  --with-fortran \
  --enable-gcc \
  --enable-lfs \
  --enable-64bit \
  --disable-shared \
  --with-zlib=$ZLIB_LIB_DIR \
  >> $LOG_FILE 2>&1

  echo '    Building'
  make >> $LOG_FILE 2>&1

  echo '    Installing'
  make install >> $LOG_FILE 2>&1

  echo '    Making links and Lib_Config'
  mv $CGNS_ADF5_SEQ_INCLUDE_DIR/cgnsBuild.defs \
    $CGNS_ADF5_SEQ_LIB_DIR/Lib_Config
  sed -i -e "s%$CGNS_DIR/%%g" $CGNS_ADF5_SEQ_LIB_DIR/Lib_Config
  rm $CGNS_ADF5_SEQ_INCLUDE_DIR/cgnsconfig.h

  cd $CGNS_DIR
  echo '  Cgns 3.2.1 library was successfully built'
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

build_z_lib

if [ $BUILD_MPI == true ]; then
  build_mpi_lib
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
echo done