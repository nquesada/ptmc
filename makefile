####################
# Makefile for GSL #
# 25-V-08          #
# Nicolas Quesada  #
####################

CC=/usr/bin/mpicc
CFLAGS=-O3 -Wall -I/usr/include -I.
LFLAGS=  -L/usr/lib/openmpi/lib -lgsl -lgslcblas -lm -L/usr/lib


%.out:%.o
	$(CC) $^ $(LFLAGS) -o $@

example_grad.out:example_grad.o lj_grad.o lj_params.o
	$(CC) $^ $(LFLAGS) -o $@


example_hess.out:example_hess.o lj_hess.o lj_params.o
	$(CC) $^ $(LFLAGS) -o $@

lj_hessian.out:lj_hessian.o lj_params.o
	$(CC) $^ $(LFLAGS) -o $@

ptmc.out:ptmc.o energy.o lj_params.o ptmc.h
	$(CC) $^ $(LFLAGS) -o $@

ptmc_data.out:ptmc_data.o ptmc.h
	$(CC) $^ $(LFLAGS) -o $@


clean:
	rm -rf *~
	rm -rf *.o
	rm -rf *.out
