# Unfinished, copied from project 4

CPPflags = c++ -O3 -std=c++11
LIB =  -larmadillo #-llapack -lblas
PROG = Gmain

${PROG}: 						${PROG}.o Gsolver.o
										${CPPflags} ${PROG}.o Gsolver.o ${LIB} -o ${PROG}.exe


${PROG}.o: 					${PROG}.cpp
										${CPPflags} -c ${PROG}.cpp

Gsolver.o: 				Gsolver.cpp Gsolver.h
										${CPPflags} -c Gsolver.cpp


clean:
	rm -f *.o *~
