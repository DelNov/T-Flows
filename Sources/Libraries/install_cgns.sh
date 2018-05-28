#!/bin/bash

# it is an useful script to automatically build cgns lib in several releases:
# 1) with adf5 and sequential access + cgnstools (optional)
# 2) with hdf5 and parallel   access
# 3) with hdf5 and sequential access + cgnstools (optional)
# (parallel cgns lib can only be based on hdf5, not adf5)

# some grid meshing software(like pointwise) use hdf5, therefore it is not possible
# to convert such mesh with adf5

# hdf5_par lib most likely can not be built with your current mpif90
# therefore script builds mpif90 first to compile hdf5_par later

# you can optionally build cgnstools (requires tcl and tk libs) to view/edit .cgns files
# https://cgns.github.io/CGNS_docs_current/cgnstools/cgnsview/index.html

# T-FlowS supports currently only 3.2.1 cgns lib, but
# you can decide also if you want to use latest cgns lib
# since 3.3 version fortran 77 and 90 support was dropped, moreover user is forced
# to use cgsize_t type instead of integer in cgns API

# CGNS lib does not use compression:
# https://cgns.github.io/FAQs.html

# build mpich? [-> hdf5_par]
BUILD_MPI=true

# build cgnstools? [-> cgnsview]
CGNS_TOOLS=true

# build latest cgns lib, overwise build 3.2.1 version
BUILD_LATEST_CGNS=true

# folder structure
CGNS_DIR=$PWD                     # this is top dir for this lib
INSTALL_DIR=$CGNS_DIR/install_dir # this dir contains CNGS/ HDF5/ MPICH/
SRC_DIR=$CGNS_DIR/src_dir         # this dir contains downloaded and built sources

# put your compilers here (gcc, gfortran are allowed if BUILD_MPI=true)
CCOMP="gcc";
FCOMP="gfortran"; # or ifort/gfortran/mpif90/mpifort/mpiifort

# optimization flags
export CCLAGS="-O3"
export FFLAGS="-O3"
export FCFLAGS="-O3"
export CFLAGS="-O3"

#-------------------------------------------------#
#---------   READ ABOVE UP TO THIS ROW   ---------#
#-------------------------------------------------#

mkdir -p $SRC_DIR/
mkdir -p $INSTALL_DIR/

# exit when any command fails
set -e

# keep track of the last executed command
trap 'last_command=$current_command; current_command=$BASH_COMMAND' DEBUG
# echo an error message before exiting
trap 'echo "\"${last_command}\" command filed with exit code $?."' EXIT


#--------- functions definitions

#------ MPICH 3.2.1
function build_mpi_lib {

	# download sources
	cd $SRC_DIR/
	wget -N http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
	tar -zxvf mpich-3.2.1.tar.gz
	mkdir -p $SRC_DIR/MPICH/
	rsync -azvh mpich-3.2.1/* $SRC_DIR/MPICH/
	rm -r mpich-3.2.1/

	# compilers
	export CC=$CCOMP
	export FC=$FCOMP

	# configure
	cd $SRC_DIR/MPICH/
	./configure \
	--prefix=$INSTALL_DIR/MPICH \
	--enable-fast=all,O3 \
	--enable-fortran=all \
	--disable-shared \
	--disable-dependency-tracking

	# build
	make
	# install
	make install

	# return
	cd $CGNS_DIR
}
#------------------------------------------

#------ HDF5 5.1.8 (Paraview 5.4.1 & Visit 2.12.3 work with HDF5 5.1.8 and not with 5.1.10)
function build_hdf5_lib {

#------ parallel version

	# download sources
	cd $SRC_DIR/
	mkdir -p $SRC_DIR/HDF5/
	git clone --depth=1 https://bitbucket.hdfgroup.org/scm/hdffv/hdf5.git --branch hdf5_1_8 ./hdf5
	rsync -azvh hdf5/* $SRC_DIR/HDF5/
	rm -rf hdf5

	# compilers
	export CC=$INSTALL_DIR/MPICH/bin/mpicc
	export FC=$INSTALL_DIR/MPICH/bin/mpif90

	# configure
	cd $SRC_DIR/HDF5/; rm -rf .git
	./configure \
	--prefix=$INSTALL_DIR/HDF5_Par \
	--enable-fortran \
	--enable-parallel \
	--disable-shared \
	--enable-production

	# build
	make
	# install
	make install

#------ sequential version
	# compilers
	export CC=$CCOMP
	export FC=$FCOMP

	# configure
	cd $SRC_DIR/HDF5/; rm -rf .git
	./configure \
	--prefix=$INSTALL_DIR/HDF5_Seq \
	--enable-fortran \
	--disable-shared \
	--enable-production

	# build
	make
	# install
	make install

	cd $CGNS_DIR
}
#------------------------------------------

#------ TCL
function build_tcl_lib {

	cd $SRC_DIR/

	# compilers
	export CC=$CCOMP
	export FC=$FCOMP

	# download sources
	cd $SRC_DIR/
	wget -N https://prdownloads.sourceforge.net/tcl/tcl8.6.8-src.tar.gz
	tar -zxvf tcl8.6.8-src.tar.gz;
	mkdir -p $SRC_DIR/TCL/
	rsync -azvh tcl8.6.8/* $SRC_DIR/TCL/
	rm -r tcl8.6.8/

	# configure
	cd TCL/unix/
	./configure \
	--prefix=$INSTALL_DIR/TCL

	# build
	make
	# install
	make install

	cd $CGNS_DIR
}
#------------------------------------------

#------ TK (requires libx11-dev from repo)
function build_tk_lib {

	cd $SRC_DIR/

	# compilers
	export CC=$CCOMP
	export FC=$FCOMP

	# download sources
	cd $SRC_DIR/
	wget -N https://prdownloads.sourceforge.net/tcl/tk8.6.8-src.tar.gz
	tar -zxvf tk8.6.8-src.tar.gz;
	mkdir -p $SRC_DIR/TCL/
	rsync -azvh tk8.6.8/* $SRC_DIR/TK/
	rm -r tk8.6.8/

	# configure
	cd TK/unix/
	./configure \
	--prefix=$INSTALL_DIR/TK \
	--with-tcl=$INSTALL_DIR/TCL/lib/

	# build
	make
	# install
	make install

	cd $CGNS_DIR
}
#------------------------------------------

function build_cgns_lib_latest {

	# download sources (latest version)
	cd $SRC_DIR/
	mkdir -p $SRC_DIR/CGNS/
	git clone --depth=1 https://github.com/CGNS/CGNS.git ./cgns
	rsync -azvh cgns/* $SRC_DIR/CGNS/
	rm -rf cgns

#------ parallel CGNS with HDF5

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	# compilers
	export CC=$INSTALL_DIR/MPICH/bin/mpicc
	export FC=$INSTALL_DIR/MPICH/bin/mpif90

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/cgnslib_latest_linux_64_hdf5_par \
	--with-hdf5=$INSTALL_DIR/HDF5_Par \
	--with-fortran \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug

	# build
	make
	# install
	make install

#------ sequential CGNS with HDF5
	# compilers
	export CC=$CCOMP
	export FC=$FCOMP

if [ $CGNS_TOOLS == true ]; then

	# make links to TCL and TK where cgns searches for them
	cd $INSTALL_DIR
	rm -rf TCL/unix TCL/generic TK/unix TK/generic TK/library

	ln -s -r -f TCL/lib        TCL/tmp; mv TCL/tmp TCL/unix
	ln -s -r -f TCL/include    TCL/tmp; mv TCL/tmp TCL/generic
	ln -s -r -f TK/lib         TK/tmp;  mv TK/tmp  TK/unix
	ln -s -r -f TK/include     TK/tmp;  mv TK/tmp  TK/generic
    ln -s -r -f TK/lib/tk8.6/  TK/tmp;  mv TK/tmp  TK/library

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/cgnslib_latest_linux_64_hdf5_seq \
	--with-hdf5=$INSTALL_DIR/HDF5_Seq \
	--with-fortran \
	--enable-gcc \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--enable-cgnstools \
	--with-tcl=$INSTALL_DIR/TCL \
	--with-tk=$INSTALL_DIR/TK \
	--datarootdir=$INSTALL_DIR/cgnslib_latest_linux_64_hdf5_seq/tcl_scripts

else # no cgns GUI tools

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/cgnslib_latest_linux_64_hdf5_seq \
	--with-hdf5=$INSTALL_DIR/HDF5_Seq \
	--with-fortran \
	--enable-gcc \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug
fi

# build
make
# install
make install

#------ sequential CGNS without HDF5
# configure
cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
./configure \
--prefix=$INSTALL_DIR/cgnslib_latest_linux_64_adf5_seq \
--with-fortran \
--enable-gcc \
--enable-lfs \
--enable-64bit \
--disable-shared \
--disable-debug

# build
make
# install
make install

cd $CGNS_DIR

echo ------------------------------------------------------------------------------------------
echo parallel   CNGS with HDF5 is now installed in $INSTALL_DIR/cgnslib_latest_linux_64_hdf5_par
echo sequential CNGS with HDF5 is now installed in $INSTALL_DIR/cgnslib_latest_linux_64_hff5_seq
echo sequential CNGS with ADF5 is now installed in $INSTALL_DIR/cgnslib_latest_linux_64_adf5_seq
if [ $CGNS_TOOLS == true ]; then
echo CNGS tools are installed in $INSTALL_DIR/cgnslib_latest_linux_64_hdf5_seq/bin/
echo You can make relative links to them for convinience
fi
echo you can safely remove $SRC_DIR/ folder with its content
echo ------------------------------------------------------------------------------------------
}
#------------------------------------------

function build_cgns_lib_3.2.1 {

	# download sources (v3.2.1)
	cd $SRC_DIR/
	wget -N https://github.com/CGNS/CGNS/archive/v3.2.1.tar.gz
	tar -zxvf v3.2.1.tar.gz;
	mkdir -p $SRC_DIR/CGNS/
	rsync -azvh CGNS-3.2.1/* $SRC_DIR/CGNS/
	rm -r CGNS-3.2.1/

#------ parallel CGNS with HDF5

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	# compilers
	export CC=$INSTALL_DIR/MPICH/bin/mpicc
	export FC=$INSTALL_DIR/MPICH/bin/mpif90

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/cgnslib_3.2.1_linux_64_hdf5_par \
	--with-hdf5=$INSTALL_DIR/HDF5_Par \
	--with-fortran \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--enable-parallel \
	--with-mpi=$INSTALL_DIR/MPICH/bin/

    # fix a bug in 3.2.1 makefile
	sed -i -e "s%-L$INSTALL_DIR/MPICH/lib\ -l%-L$INSTALL_DIR/MPICH/lib\ -l mpi%g" make.defs
	sed -i -e "s%$INSTALL_DIR/MPICH/bin/mpiexec%$INSTALL_DIR/MPICH/bin/mpiexec  -n \$\${NPROCS:=4}%g" make.defs

	# build
	make
	# install
	make install

#------ sequential CGNS with HDF5
	# compilers
	export CC=$CCOMP
	export FC=$FCOMP

if [ $CGNS_TOOLS == true ]; then

	# make links to TCL and TK where cgns searches for them
	cd $INSTALL_DIR
	rm -rf TCL/unix TCL/generic TK/unix TK/generic TK/library

	ln -s -r -f TCL/lib        TCL/tmp; mv TCL/tmp TCL/unix
	ln -s -r -f TCL/include    TCL/tmp; mv TCL/tmp TCL/generic
	ln -s -r -f TK/lib         TK/tmp;  mv TK/tmp  TK/unix
	ln -s -r -f TK/include     TK/tmp;  mv TK/tmp  TK/generic
    ln -s -r -f TK/lib/tk8.6/  TK/tmp;  mv TK/tmp  TK/library

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/cgnslib_3.2.1_linux_64_hdf5_seq \
	--with-hdf5=$INSTALL_DIR/HDF5_Seq \
	--with-fortran \
	--enable-gcc \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug \
	--enable-cgnstools \
	--with-tcl=$INSTALL_DIR/TCL \
	--with-tk=$INSTALL_DIR/TK

else # no cgns GUI tools

	# configure
	cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

	FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
	./configure \
	--prefix=$INSTALL_DIR/cgnslib_3.2.1_linux_64_hdf5_seq \
	--with-hdf5=$INSTALL_DIR/HDF5_Seq \
	--with-fortran \
	--enable-gcc \
	--enable-lfs \
	--enable-64bit \
	--disable-shared \
	--disable-debug
fi

# build
make
# install
make install

#------ sequential CGNS without HDF5
# configure
cd $SRC_DIR/CGNS/; rm -rf .git; cd src/

FLIBS=-Wl,--no-as-needed\ -ldl\ -lz \
LIBS=-Wl,--no-as-needed\ -ldl\ -lz \
./configure \
--prefix=$INSTALL_DIR/cgnslib_3.2.1_linux_64_adf5_seq \
--with-fortran \
--enable-gcc \
--enable-lfs \
--enable-64bit \
--disable-shared \
--disable-debug
# build
make
# install
make install

cd $CGNS_DIR

echo ---------------------------------------------------------------------------
echo parallel CNGS with HDF5 is now installed in $INSTALL_DIR/cgnslib_3.2.1_linux_64_hdf5_par
echo sequential CNGS with HDF5 is now installed in $INSTALL_DIR/cgnslib_3.2.1_linux_64_hff5_seq
echo sequential CNGS with ADF5 is now installed in $INSTALL_DIR/cgnslib_3.2.1_linux_64_adf5_seq
if [ $CGNS_TOOLS == true ]; then
echo CNGS tools are installed in $INSTALL_DIR/cgnslib_3.2.1_linux_64_hdf5_seq/bin/
echo You can make relative links to them for convinience
fi
echo you can safely remove $SRC_DIR/ folder with its content
echo ---------------------------------------------------------------------------
}

#--------- script with functions defined above

if [ $BUILD_MPI == true ]; then
	build_mpi_lib
fi

build_hdf5_lib

if [ $CGNS_TOOLS == true ]; then
	build_tcl_lib
	build_tk_lib
fi

if [ $BUILD_LATEST_CGNS == true ]; then
	build_cgns_lib_latest
else
	build_cgns_lib_3.2.1
fi

echo done

# build and run tests for self-confidence in cgsn_src dirs

# cd tests
# make
# make test
# cd ../examples/fortran
# make
# make test
# cd ../../Test_UserGuideCode/Fortran_code
# make
# make test
# cd ../C_code
# make
# make test
# build and run tests for self-confidence
