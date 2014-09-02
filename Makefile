
f2py = f2py 
f2pyflags = -c

strategy2opsimout:
	$(f2py) -c fortran/lsstdfacg.f -m strategy2opsimin

opsim2simlib:
	$(f2py) -c fortran/opsim2simlib.f90 -m opsim2simlib
default: strategy2opsimout opsim2simlib 

clean:
	rm -f *.so
