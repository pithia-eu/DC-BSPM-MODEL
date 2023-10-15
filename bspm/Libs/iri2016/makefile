RM      = rm -f

.phony: all
all: irilib64.so

irilib64.so : proton_temp_7jun2016.f90 electron_temp_7jun2016.f90  oxygen_ratio_7jun2016.f90 helium_ratio_7jun2016.f90 proton_ratio_wrap.f90 proton_ratio_7jun2016.f90 elecdens7jun2016.f90 irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for
	gfortran -shared -fPIC -fdefault-integer-8 -fdefault-real-8 -o irilib64.so proton_temp_7jun2016.f90 electron_temp_7jun2016.f90 oxygen_ratio_7jun2016.f90  helium_ratio_7jun2016.f90 proton_ratio_7jun2016.f90 elecdens7jun2016.f90 irisub.for irifun.for iritec.for iridreg.for igrf.for cira.for iriflip.for

.phony: clean
clean:
	$(RM) *.o irilib64.so
