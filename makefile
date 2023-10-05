FC := gfortran
FCFLAGS = -O3

%.o: %.f90
	$(FC) $(FCFLAGS) -c $<

PROGRAM = get_field_calc
OBJECTS = field_module.f90 calculate_field.f90 read_field_input.f90 fieldcode.f90

$(PROGRAM): $(OBJECTS)
	$(FC) $(FCFLAGS) $(OBJECTS) -o  $(PROGRAM)