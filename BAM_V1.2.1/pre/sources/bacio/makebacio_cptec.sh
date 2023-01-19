#! /bin/bash -x

COMP=${1}

if [ ${COMP} == "linux_intel_egeon" ]
then
  CC=/opt/intel/oneapi/mpi/2021.4.0/bin/mpiicc
elif [ ${COMP} == "linux_gnu_egeon" ]
then
  CC=/opt/ohpc/pub/mpi/openmpi4-gnu9/4.1.1/bin/mpicc
else
  CC=mpicc
fi

mkdir -p NCEPLIBS-bacio-2.5.0/build
cd NCEPLIBS-bacio-2.5.0/build
CC=${CC} cmake -DCMAKE_INSTALL_PREFIX=../../NCEPLIBS-bacio-2.5.0/rls ..
make -j2
make install
