#Makefile for MUCG program.

fc = gfortran
fflags = -J mods -I mods -O3 -finit-integer=2147483648 -finit-real=nan -finit-character=35 -fmax-stack-var-size=512

run_aio: mucg_aio | rund
	cd rund && ../mucg_aio

run: mucg | rund
	cd rund && ../mucg


mucg_aio: aiosrc.f90 | mods
	$(fc) $(fflags)  $^ -o $@

aiosrc.f90: src/functions.f90 src/mucg.f90
# 	sed s:'ALL\.':'src/ALL.':g $^ > $@
	cat $^ > $@

mucg:  src/functions.f90 src/mucg.f90 | mods
	$(fc) $(fflags)  $^ -o $@


rund mods:
	mkdir $@

clean:
	rm -rf rund mucg mucg_aio mods aiosrc.f90
