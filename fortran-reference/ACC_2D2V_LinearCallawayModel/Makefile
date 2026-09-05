PROG = DGACC
PREFIX = ./bin/
# IFROOT=/opt/intel/oneapi/compiler/2021.2.0/linux/compiler/lib/intel64_lin
# MKLROOT=/opt/intel/mkl/lib/intel64

FCOMP = ifort
# FLAG = -module $(PREFIX) -I$(PREFIX) -g -qopenmp -warn #-check all
FLAG = -O3 -march=broadwell -module $(PREFIX) -I$(PREFIX) -qopenmp
LIB_LIST = -qmkl

# -lpardiso500-INTEL1301-X86-64

# FCOMP = gfortran
# FLAG =  -O3 -J$(PREFIX) -I$(PREFIX) -fopenmp
# FLAG =  -g -J$(PREFIX) -I$(PREFIX) -Wall -Wextra -fcheck=all -fbacktrace -fopenmp
# FLAG =  -g -J$(PREFIX) -I$(PREFIX) -fcheck=all -fbacktrace -fopenmp
# LIB_LIST = -L${MKLROOT} -Wl,--no-as-needed -lmkl_sequential -lmkl_core -lpthread -fopenmp

OBJ = $(PREFIX)Global_Parameter.o $(PREFIX)USD_Math.o $(PREFIX)Boundary_Conditions.o $(PREFIX)Spatial_Mesh.o \
	 $(PREFIX)Velocity_Mesh.o $(PREFIX)Matrix.o $(PREFIX)Basis_Function.o $(PREFIX)Integration.o $(PREFIX)Velocity_Distribution.o \
	 $(PREFIX)Solvers.o $(PREFIX)Synthetic_Acceleration.o $(PREFIX)Out_Put_Result.o $(PREFIX)Callaway_2D_2V_DG.o


all : $(PREFIX) $(PROG)

$(PREFIX) :
	mkdir $(PREFIX)

$(PROG) : $(OBJ)
	$(FCOMP) $(OBJ) $(LIB_LIST) -o $(PROG)

$(PREFIX)%.o : %.f90
	$(FCOMP) $(FLAG) -c $< -o $@

#$(PREFIX)%.mod : %.mod
#	mv %.mod -f -t $(PREFIX)

clean :
	rm $(PREFIX)*.o $(PREFIX)*.mod $(PROG)


