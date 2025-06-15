# Makefile for compiling the Serial Spectral Code Version 2.1 
# This uses the Intel Fortran Compiler, with FFTW 
# ---------------------------------------------------------------------
#
compile =gfortran -g -ffree-line-length-none -fcheck=bounds -fopenmp
link = -lrfftw_threads -lfftw_threads -lrfftw -lfftw  -L/home/ritwik.mukherjee/Software/FFTW2/lib

#link = -lrfftw_threads -lfftw_threads -lrfftw -lfftw  -L/home/prateek.anand/softwares/fftw-2.1.5/lib
#-----------------------------------------------------------------------
OBJ_FILE=spectral_serial.o  global_array.o  \
energy_serial.o rnkt_serial.o force_serial.o adbsh_serial.o nonlin_serial.o \
energy_injection.o dermat_serial.o  omplib.o 
 
MOD_FILE=mod_serial_fluid.o
MOD_FILE1=mod_initial.o
MOD_FILE2=mod_particle.o
test : 
	echo $(OBJ_FILE) 
# ---------------------------------------------------------------------
part.exe:$(OBJ_FILE) $(MOD_FILE) $(MOD_FILE1) $(MOD_FILE2)  
	$(compile) $(OBJ_FILE) $(MOD_FILE) $(MOD_FILE1) $(MOD_FILE2) -o part.exe $(link)
#                 ------------------------------------               
spectral_serial.o: spectral_serial.f90 $(MOD_FILE) $(MOD_FILE1) $(MOD_FILE2) 
	$(compile) $(link) -c -o spectral_serial.o \
	spectral_serial.f90 $(MOD_FILE) $(MOD_FILE1) $(MOD_FILE2)
#         
global_array.o: global_array.f90 $(MOD_FILE) 
	$(compile) -c -o global_array.o \
	global_array.f90 $(MOD_FILE)
#         
energy_serial.o: energy_serial.f90 $(MOD_FILE) 
	$(compile) -c -o energy_serial.o \
	energy_serial.f90 $(MOD_FILE)
#         
adbsh_serial.o: adbsh_serial.f90 $(MOD_FILE) 
	$(compile) -c -o adbsh_serial.o \
	adbsh_serial.f90 $(MOD_FILE)
#         
nonlin_serial.o: nonlin_serial.f90 $(MOD_FILE) $(MOD_FILE2)
	$(compile) -c -o nonlin_serial.o  \
	nonlin_serial.f90 $(MOD_FILE) $(MOD_FILE2)
#         
force_serial.o: force_serial.f90 $(MOD_FILE) 
	$(compile) -c -o force_serial.o \
	force_serial.f90 $(MOD_FILE)
#
energy_injection.o: energy_injection.f90 
	$(compile) -c -o energy_injection.o \
	energy_injection.f90 $(MOD_FILE)

#         
rnkt_serial.o: rnkt_serial.f90 $(MOD_FILE) 
	$(compile) -c -o rnkt_serial.o \
	rnkt_serial.f90 $(MOD_FILE)
#         
dermat_serial.o: dermat_serial.f90 $(MOD_FILE) 
	$(compile) -c -o dermat_serial.o \
	dermat_serial.f90 $(MOD_FILE)
#        
mod_serial_fluid.o:  mod_serial_fluid.f90 
	$(compile) -c -o mod_serial_fluid.o \
	mod_serial_fluid.f90
# ------------------------------------------------------------------------ 
mod_initial.o:  mod_initial.f90 $(MOD_FILE) omplib.o
	$(compile) -c -o mod_initial.o \
	mod_initial.f90 $(MOD_FILE) omplib.o
# ------------------------------------------------------------------------ 
mod_particle.o:  mod_particle.f90 $(MOD_FILE)
	$(compile) -c -o mod_particle.o \
	mod_particle.f90 $(MOD_FILE)
# ------------------------------------------------------------------------ 
omplib.o:  omplib.f90 $(MOD_FILE) 
	$(compile) -c -o omplib.o \
	omplib.f90 $(MOD_FILE)
# ------------------------------------------------------------------------ 
decimation.o: decimation.f90 $(MOD_FILE) 
	$(compile) -c -o decimation.o \
	decimation.f90 $(MOD_FILE)
#-------------------------------------------         
clean:
	rm -f *.o *.mod core *.cmp.f90 *.kmo work.* *.d
total_clean:
	rm -f *.o *.mod core *.exe *.out
clean_out:
	rm -f *.out
## --------------------------------------------------------------------##  
