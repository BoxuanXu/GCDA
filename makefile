
INS_IBM=-I/share/soft/rhel6.2/netcdf-4.1.1Parallel/include
LD_IBM=-L/share/soft/rhel6.2/netcdf-4.1.1Parallel/lib

INS_HDF5=-I/share/soft/rhel6.2/hdf5-1.8.8Parallel/include/
LD_HDF5=-L/share/soft/rhel6.2/hdf5-1.8.8Parallel/lib

INS_EIGEN=-I/data1/xubx/develop_tools_az/eigen-eigen-36fd1ba04c12/eigen-eigen-36fd1ba04c12


NC_LIB=-lnetcdf_c++ -lnetcdf
HDF_LIB=-lhdf5_hl -lhdf5
FLAG=-std=c++0x
CXX=mpic++
CXX_IBM=mpic++ 


all:main


run:main
	bsub -I -n 20 -R "span[ptile=1]" mpirun -np 50 ./main
    
debug:main_d
	mpirun -np 2 ./main_d

main:main.cpp *.h
	${CXX_IBM} ${FLAG} -o $@ $< ${INS_EIGEN} ${INS_IBM} ${INS_HDF5} -lz ${LD_IBM} ${NC_LIB} ${LD_HDF5} ${HDF_LIB}

main_d:main.cpp *.h
	${CXX} -lnetcdf ${FLAG} -DDEBUG -g ${INS_local} ${LD_local} ${LIB} -o $@ $<


