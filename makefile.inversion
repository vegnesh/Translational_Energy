DIR=$(PWD)
SRCDIR=$(DIR)
OBJ=$(DIR)/objects

# Fortran Compiler
FC = gfortran
LD = gfortran

FFLAGS= -g  -fdefault-real-8 
LDFLAGS= -g  -fdefault-real-8 

SRCF90=  special_functions.f90 global_parameters.f90 bin_basic.f90   bin_parameters.f90 jacobian_energy.f90  inversion_onevar.f90
OBJF90=$(SRCF90:.f90=.obj)
OBJS=$(patsubst %,$(OBJ)/%,$(OBJF90))

EXEC=out_inversion.out

default: mytest


mytest: $(OBJS) $(EXEC)

$(EXEC): $(OBJS)
	@echo ""
	@echo "Building Test code"
	@echo ""
	$(LD) $(LDFLAGS) $(OBJS) -o $@ 

$(OBJ)/%.obj: $(SRCDIR)/%.f90
	$(FC) $(FFLAGS) -c $< -o $@

clean:  
	rm *.mod *.out $(OBJ)/*.obj
