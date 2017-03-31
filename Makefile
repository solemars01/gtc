#key parameters:
OPENMP=y
DOUBLE_PRECISION=n
ADIOS=n
BIGDATA=n
PETSc=y
DEBUG=n
GPU=n
GPU_UM=n
COMPILER=default
ICONFIG=TOKAMAK
#ICONFIG=FRC
#ICONFIG=CYLINDER
#ICONFIG=TOROIDAL3D
#COORDSYS=CYLINDRICAL
COORDSYS=BOOZER
############################################################################
#             Makefile to build the GTC code
#           ==================================
#
# You only need to type "gmake" to build the code on the platforms
# defined below. The makefile runs the "uname -s" command to detect
# the operating system automatically. Other options are:
#
#  % gmake OPENMP=y       Builds the code with OpenMP support
#  % gmake OPENMP=n       Builds the code WITHOUT OpenMP support
#                         The default is with OpenMP support	
#
#  % gmake DOUBLE_PRECISION=y  Builds with 8-byte floating point precision
#                              The default is single precision
#
#  % gmake PETSc=y        Builds the code using the PETSc parallel matrix
#                         solving library.  This is necessary to run the
#                         electromagnetic version of GTC.
#
#  % gmake DEBUG=y        Compiles the files with debugging flags.
#                         The default is no debug option.
#
#  If COMPILER is set to 'default' the Makefile will select the appropriate
#  compiler for the target.  If the user wished to override the target
#  machine default, COMPILER may be set to 'portland', 'intel', or 'gnu'.
#
# Special targets:
#
#  % gmake clean      Removes the executable and all object files (*.o)
#
#  % gmake cleanomp   Removes the executable and the object files
#                     containing OpenMP directives
#
#  % gmake doc        Rebuilds the documentation.
#
#############################################################################
# Default executable name
CMD := gtc
ALL:${CMD}
LIB :=
LD_LIB :=

##### inquire the hostname to judge the system name, to determine what target we're compiling on
ifneq (,${HOST})
  SYSTEMS := ${HOST}
else
  SYSTEMS := $(shell hostname)
endif

##### target machine determines the default compiler, and paths to libraries #####
# edison
ifneq (,$(findstring edison,$(SYSTEMS)))
   ifeq ($(COMPILER),default)
      COMPILER := intel
   endif
   CMP := ftn
#   PETSC_OPT = -D_PETSc35BEFORE ##only needed for kink/tearing
   PETSC_LIB = -I${PETSC_DIR}/include -I${PETSC_DIR}/include/finclude
   NETCDF_LIB = -I${NETCDF_DIR}/include
   NETCDF_LD_LIB = -L${NETCDF_DIR}/lib -lnetcdf
   LIBD = /global/u2/x/xiao/AdioFranklin/lib
   INCD = /global/u2/x/xiao/AdioFranklin/include
   ADIOS_LIB = -I$(INCD)
   ADIOS_LD_LIB = -L$(LIBD) -ladios -lmxml
endif

# cori
ifneq (,$(findstring cori,$(SYSTEMS)))
   ifeq ($(COMPILER),default)
      COMPILER := intel
   endif
   CMP := ftn
#   PETSC_OPT = -D_PETSc35BEFORE ##only needed for kink/tearing
   PETSC_LIB = -I${PETSC_DIR}/include -I${PETSC_DIR}/include/finclude
   NETCDF_LIB = -I${NETCDF_DIR}/include
   NETCDF_LD_LIB = -L${NETCDF_DIR}/lib -lnetcdf
   LIBD = /global/u2/x/xiao/AdioFranklin/lib
   INCD = /global/u2/x/xiao/AdioFranklin/include
   ADIOS_LIB = -I$(INCD)
   ADIOS_LD_LIB = -L$(LIBD) -ladios -lmxml
endif

# hopper
ifneq (,$(findstring hopper,$(SYSTEMS)))
   ifeq ($(COMPILER),default)
      COMPILER := portland
   endif
   CMP := ftn
   PETSC_LIB = -I${PETSC_DIR}/include -I${PETSC_DIR}/include/finclude
   NETCDF_LIB = -I${NETCDF_DIR}/include
   NETCDF_LD_LIB = -L${NETCDF_DIR}/lib -lnetcdf
   LIBD = /global/u2/x/xiao/AdioFranklin/lib
   INCD = /global/u2/x/xiao/AdioFranklin/include
   ADIOS_LIB = -I$(INCD)
   ADIOS_LD_LIB = -L$(LIBD) -ladios -lmxml
endif

# titan
ifneq (,$(findstring titan,$(SYSTEMS)))
   ifeq ($(COMPILER),default)
      ifeq ($(GPU),y)
	COMPILER := portland
      else
	COMPILER := intel
      endif
   endif
   CMP :=ftn
#   PETSC_OPT = -D_PETSc35BEFORE##only needed for kink/tearing
   PETSC_LIB = -I${PETSC_DIR}/include/petsc -I${PETSC_DIR}/include/petsc/finclude
   NETCDF_LIB = -I${NETCDF_DIR}/include
   NETCDF_LD_LIB = -L${NETCDF_DIR}/lib
   LIBD = /ccs/home/hardes/gtc_adios/jaguarlib
   INCD = /ccs/home/hardes/gtc_adios/jaguarinclude
   ADIOS_LIB = -I$(INCD)
   ADIOS_LD_LIB = -L$(LIBD) -ladios -lmxml
endif

# Tianhe-1A
ifneq (,$(findstring ln,$(SYSTEMS)))
   ifeq ($(COMPILER),default)
      ifeq ($(GPU),y)
	COMPILER := portland
      else
	COMPILER := intel
      endif
   endif
   CMP := mpif90

   PETSC_solver = petsc.o
   PETSC_OPT = -D_PETSc31P8ANDBEFORE
	 #PETSC_OPT += -D_PETSc30ANDBEFORE
   ifneq (,${PETSC_DIR})
      include ${PETSC_DIR}/conf/variables
   endif

   PETSC_LIB := -I${PETSC_DIR}/include -I${PETSC_DIR}/include/finclude 
   PETSC_LD_LIB := ${PETSC_KSP_LIB}

   NETCDF_LIB := -I${NETCDF_DIR}/include
   NETCDF_LD_LIB := -L${NETCDF_DIR}/lib -lnetcdf -lnetcdff -lnetcdf_c++\
		 -Wl,--rpath -Wl,${NETCDF_DIR}/lib
endif

# PKU cluster
ifneq (,$(findstring mgt,$(SYSTEMS)))
   PETSC_LIB = -Wl,-rpath,${PETSC_LIB_DIR} -L${PETSC_LIB_DIR} -lpetsc -lX11\
							 -L/share/soft/lapack-3.3.1 -llapack_LINUX -lblas_LINUX
endif


#iop
ifneq (,$(findstring indac,$(SYSTEMS)))
   ifeq ($(COMPILER),default)
      COMPILER := intel
   endif
   CMP := mpif90

   #PETSC_solver = petsc.o
   PETSC_OPT = -D_PETSc31P8ANDBEFORE
         #PETSC_OPT += -D_PETSc30ANDBEFORE
   ifneq (,${PETSC_DIR})
      include ${PETSC_DIR}/conf/variables
   endif

   PETSC_LIB := -I${PETSC_DIR}/include -I${PETSC_DIR}/include/finclude
   PETSC_LD_LIB := ${PETSC_KSP_LIB}

   NETCDF_LIB := -I${NETCDF_DIR}/include
   NETCDF_LD_LIB := -L${NETCDF_DIR}/lib -lnetcdf -lnetcdff -lnetcdf_c++\
                 -Wl,--rpath -Wl,${NETCDF_DIR}/lib
endif

#### compiler determines debugging, openmp, and double precision flags #####
# intel fortran
ifeq ($(COMPILER),intel)
   OPTIMOPT := -O
   #OPTIMOPT := -O3
   #OPTIMOPT += -fp-model fast=2
   #OPTIMOPT += -align array64byte
   #OPTIMOPT += -no-vec

   DEBUGOPT := -g -free -pc64 -check all -warn nounused -warn all -debug all -debug-parameters all -fp-stack-check
   OMPOPT := -qopenmp
   DPOPT := -fpconstant
endif

# portland group fortran
ifeq ($(COMPILER),portland)
   OPTIMOPT := -O3 -Mfree -Kieee
   DEBUGOPT := -g -C -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf3 -Melf -Mpgicoff -traceback -Minform=inform
   OMPOPT := -mp
   DPOPT := -DDOUBLE_PRECISION
endif

# GNU gfortran
ifeq ($(COMPILER),gnu)
   OPTIMOPT := 
   DEBUGOPT :=
   OMPOPT :=
   DPOPT :=
endif


### assemble all the options
OPT :=
LIB :=
ifeq ($(ICONFIG),FRC)
  ICONFIG_OPT = -D_FRC
  OPT += $(ICONFIG_OPT)
else ifeq ($(ICONFIG),TOKAMAK)
  ICONFIG_OPT = -D_TOKAMAK
  OPT += $(ICONFIG_OPT)
else ifeq ($(ICONFIG),CYLINDER)
  ICONFIG_OPT = -D_CYLINDER
  OPT += $(ICONFIG_OPT)
else ifeq ($(ICONFIG),TOROIDAL3D)
  ICONFIG_OPT = -D_TOROIDAL3D
  OPT += $(ICONFIG_OPT)
endif

ifeq ($(COORDSYS),CYLINDRICAL)
  COORDSYS_OPT = -D_CYLINDRICALCOORDSYS
  OPT += $(COORDSYS_OPT)
else ifeq ($(COORDSYS),BOOZER)
  COORDSYS_OPT = -D_BOOZERCOORDSYS
  OPT += $(COORDSYS_OPT)
endif

ifeq ($(BIGDATA),y)
  DATAOUT3D ?= dataout3d.o
  BIGDATA_OPT = -D_BIGDATA
  OPT += $(BIGDATA_OPT)
  ifeq ($(ADIOS),y)
    OPT := -DADIOS=1
    LIB := $(ADIOS_LIB)
    LD_LIB := $(ADIOS_LD_LIB)
  else
    LIB += $(NETCDF_LIB)
    LD_LIB += $(NETCDF_LD_LIB)
  endif
endif
ifeq ($(PETSc),y)
   PETSC_solver ?= petsc.o
   PETSC_OPT += -D_PETSc
   OPT += $(PETSC_OPT)
   LIB += $(PETSC_LIB)
   LD_LIB += $(PETSC_LD_LIB)
else
   PETSC_solver :=
endif
ifeq ($(OPENMP),y)
   OPT += $(OMPOPT)
endif
ifeq ($(DOUBLE_PRECISION),y)
   OPT += -DDOUBLE_PRECISION $(DPOPT)
endif
ifeq ($(DEBUG),y)
   OPT += $(DEBUGOPT)
else
   OPT += $(OPTIMOPT)
endif

ifeq ($(GPU),y)
  ifeq ($(GPU_UM),y)
    OPT := -DGPU_UM -Mcuda=cuda7.5 -Minfo=accel -acc -ta=nvidia:cc35,ptxinfo,maxregcount:64,managed $(OPT)
   else
    OPT := -Mcuda=cuda7.5 -Minfo=accel -acc -ta=nvidia:cc35,ptxinfo,maxregcount:64 $(OPT)
  endif
  #OPT := -Mcuda=cuda7.0 -Minfo=all -acc -ta=nvidia:cc35,ptxinfo,maxregcount:64,time,host $(OPT)
  NVCC := nvcc
  CUDA_OPT := -arch=sm_35 -O3\
	-I/opt/cray/mpt/7.4.0/gni/mpich-pgi/15.3/include
  CUDA_OBJ := shift_cuda.o
endif


##################################################################
# We add ".F90" to the list of suffixes to allow source files on which the
# co-processor will be run automatically.
.SUFFIXES: .o .F90 .F
.PHONY: clean

# List of all the object files needed to build the code
OBJ:=module.o function.o M3DC1_interface.o main.o setup.o restart.o diagnosis.o\
 snapshot.o poisson.o smooth.o field.o pushifk.o pushifk_boris.o\
 tracking.o $(DATAOUT3D) pushfield.o taehdf5.o eqdata.o eqplot.o\
 collision.o $(PETSC_solver) fft_gl.o SDP_interface.o\
 load.o push.o shift.o charge.o $(CUDA_OBJ)
HYBRID_OBJ:=hybridChargeParticle.o hybridPushParticle.o
KINETIC_OBJ:=gkChargeParticle.o gkPushParticle.o

OBJ += $(HYBRID_OBJ) $(KINETIC_OBJ)

ifeq ($(COORDSYS),CYLINDRICAL)
OBJ +=loadeRZ.o pusheRZ.o
endif

$(CMD): $(OBJ)
	$(CMP) $(OPT) -o $(CMD) $(OBJ) $(LD_LIB)

$(filter-out module.o,$(OBJ)): module.o

module.o : module.F90
	$(CMP) $(OPT) $(LIB) -c module.F90

M3DC1_interface.o: M3DC1_interface.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c M3DC1_interface.F90

eqdata.o: eqdata.F90 module.o M3DC1_interface.o function.o
	$(CMP) $(OPT) $(LIB) -c eqdata.F90

eqplot.o: eqplot.F90 module.o SDP_interface.o function.o
	$(CMP) $(OPT) $(LIB) -c eqplot.F90

snapshot.o: snapshot.F90 module.o SDP_interface.o function.o
	$(CMP) $(OPT) $(LIB) -c snapshot.F90

SDP_interface.o: SDP_interface.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c SDP_interface.F90

poisson.o: poisson.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c poisson.F90

setup.o: setup.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c setup.F90

load.o: load.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c load.F90

ifeq ($(COORDSYS),CYLINDRICAL)
loadeRZ.o: loadeRZ.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c loadeRZ.F90
pusheRZ.o: pusheRZ.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c pusheRZ.F90
endif

pushifk_boris.o: pushifk_boris.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c pushifk_boris.F90 

function.o: function.F90 module.o
	$(CMP) $(OPT) $(LIB) -c function.F90

#json_module.o: json_module.F90
#	$(CMP) $(OPT) $(LIB) -c json_module.F90

hybrid%.F90 : module.o PushParticle.F90 ChargeParticle.F90
	$(CMP) $(OPT) $(LIB) -D_hybrid -E $*.F90 > $@

gk%.F90 : module.o PushParticle.F90 ChargeParticle.F90
	$(CMP) $(OPT) $(LIB) -E $*.F90 > $@

.F90.o : module.o
	$(CMP) $(OPT) $(LIB) -c $<

shift_cuda.o : shift_cuda.cu
	$(NVCC) $(CUDA_OPT) -c shift_cuda.cu

diagnosis.o: diagnosis.F90 module.o function.o
	$(CMP) $(OPT) $(LIB) -c diagnosis.F90

# The following tag is meant to "clean" the directory by removing the
# executable along with all the object files created by the compilation 
# One only has to run:  gmake clean

clean::
	rm -f $(CMD) $(OBJ) *.mod *.o *genmod.f90 *.optrpt
