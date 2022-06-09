#!/bin/bash
src=implicit.c
out=implicit
/localhome1/home/u12132414/petsc-build/bin/mpicc -wd1572 -Wno-unknown-pragmas -g -O0  -wd1572 -Wno-unknown-pragmas -g -O0    -I/localhome1/home/u12132414/petsc/include -I/localhome1/home/u12132414/petsc/test/include -I/localhome1/home/u12132414/petsc-build/include     ${src}  -Wl,-rpath,/localhome1/home/u12132414/petsc/test/lib -L/localhome1/home/u12132414/petsc/test/lib -Wl,-rpath,/localhome1/home/u12132414/petsc-build/lib -L/localhome1/home/u12132414/petsc-build/lib -lpetsc -lf2clapack -lf2cblas -lpthread -lX11 -ldl -o ${out}
/localhome1/home/u12132414/petsc-build/bin/mpirun -n 1 ./${out}
