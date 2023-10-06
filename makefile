FC := gfortran
FCFLAGS = -O3

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

PROGRAM = get_field_calc
OBJECTS = field_module.f90 fieldcode.f90

$(PROGRAM): $(OBJECTS)
	$(FC) $(FCFLAGS) $(OBJECTS) -o  $(PROGRAM)