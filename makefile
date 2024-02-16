include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules


testPetsc: simplePetsc.o  
	-${CLINKER} -o $@ $^  ${PETSC_LIB}

testPetsc2: simplePetsc2.o
	-${CLINKER} -o $@ $^ ${PETSC_LIB}

testPetsc3: simplePetsc3.o
	-${CLINKER} -o $@ $^ ${PETSC_LIB}

firstPde: petscEx.o
	-${CLINKER} -o $@ $^ ${PETSC_LIB}
test: test.o
	-${CLINKER} -o $@ $^ ${PETSC_LIB}