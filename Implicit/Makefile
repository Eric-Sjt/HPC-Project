
# -*- mode: makefile -*-
PETSC_DIR  := /work/mae-songjt/lib/petsc-3.17.2/
PETSC_ARCH := arch-darwin-c-debug

# This sample (GNU) Makefile can be used to compile PETSc applications with a single
# source file and can be easily modified to compile multi-file applications.
# It relies on pkg_config tool, and PETSC_DIR and PETSC_ARCH variables.
# Copy this file to your source directory as "Makefile" and modify as needed.
#
# For example - a single source file can be compiled with:
#
#  $ cd src/snes/tutorials/
#  $ make -f $PETSC_DIR/share/petsc/Makefile.user ex17
#
# The following variable must either be a path to petsc.pc or just "petsc" if petsc.pc
# has been installed to a system location or can be found in PKG_CONFIG_PATH.
petsc.pc := $(PETSC_DIR)/$(PETSC_ARCH)/lib/pkgconfig/petsc.pc

# Additional libraries that support pkg-config can be added to the list of PACKAGES below.
PACKAGES := $(petsc.pc)

CC := $(shell pkg-config --variable=ccompiler $(PACKAGES))
CXX := $(shell pkg-config --variable=cxxcompiler $(PACKAGES))
FC := $(shell pkg-config --variable=fcompiler $(PACKAGES))
CFLAGS_OTHER := $(shell pkg-config --cflags-only-other $(PACKAGES))
CFLAGS := $(shell pkg-config --variable=cflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
CXXFLAGS := $(shell pkg-config --variable=cxxflags_extra $(PACKAGES)) $(CFLAGS_OTHER)
FFLAGS := $(shell pkg-config --variable=fflags_extra $(PACKAGES))
CPPFLAGS := $(shell pkg-config --cflags-only-I $(PACKAGES))
LDFLAGS := $(shell pkg-config --libs-only-L --libs-only-other $(PACKAGES))
LDFLAGS += $(patsubst -L%, $(shell pkg-config --variable=ldflag_rpath $(PACKAGES))%, $(shell pkg-config --libs-only-L $(PACKAGES)))
LDLIBS := $(shell pkg-config --libs-only-l $(PACKAGES)) -lm
CUDAC := $(shell pkg-config --variable=cudacompiler $(PACKAGES))
CUDAC_FLAGS := $(shell pkg-config --variable=cudaflags_extra $(PACKAGES))
CUDA_LIB := $(shell pkg-config --variable=cudalib $(PACKAGES))
CUDA_INCLUDE := $(shell pkg-config --variable=cudainclude $(PACKAGES))

print:
	@echo CC=$(CC)
	@echo CXX=$(CXX)
	@echo FC=$(FC)
	@echo CFLAGS=$(CFLAGS)
	@echo CXXFLAGS=$(CXXFLAGS)
	@echo FFLAGS=$(FFLAGS)
	@echo CPPFLAGS=$(CPPFLAGS)
	@echo LDFLAGS=$(LDFLAGS)
	@echo LDLIBS=$(LDLIBS)
	@echo CUDAC=$(CUDAC)
	@echo CUDAC_FLAGS=$(CUDAC_FLAGS)
	@echo CUDA_LIB=$(CUDA_LIB)
	@echo CUDA_INCLUDE=$(CUDA_INCLUDE)

# Files
EXEC := implicit.out
SRC := implicit.c

#Rules
explicit.out : explicit.c
	/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/bin/mpiicc -fPIC -wd1572 -Wno-unknown-pragmas -g -O0  -fPIC -wd1572 -Wno-unknown-pragmas -g -O0    -I/work/mae-songjt/lib/petsc-3.17.2/include -I/work/mae-songjt/lib/petsc-3.17.2/arch-darwin-c-debug/include -I/work/mae-songjt/lib/petsc-3.17.2-build/include -I/share/intel/2018u4/compilers_and_libraries_2018.5.274/linux/mpi/intel64/include     $<  -Wl,-rpath,/work/mae-songjt/lib/petsc-3.17.2/arch-darwin-c-debug/lib -L/work/mae-songjt/lib/petsc-3.17.2/arch-darwin-c-debug/lib -Wl,-rpath,/work/mae-songjt/lib/petsc-3.17.2-build/lib -L/work/mae-songjt/lib/petsc-3.17.2-build/lib -lpetsc -lf2clapack -lf2cblas -lhdf5_hl -lhdf5 -lX11 -lquadmath -ldl -o $@
.PHONY : clean
clean:
	$(RM) $(EXEC) $(OBJ)
