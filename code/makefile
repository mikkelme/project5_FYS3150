# Unfinished, copied from project 4

CPPflags = c++ -O3 -std=c++11
LIB =  -larmadillo #-llapack -lblas
PROG = main

${PROG}: 						${PROG}.o solver.o
										${CPPflags} ${PROG}.o solver.o ${LIB} -o ${PROG}.exe


${PROG}.o: 					${PROG}.cpp
										${CPPflags} -c ${PROG}.cpp

solver.o: 				solver.cpp solver.h
										${CPPflags} -c solver.cpp


clean:
	rm -f *.o *~
