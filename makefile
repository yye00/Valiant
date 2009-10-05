include ${PETSC_DIR}/conf/base

all:		Valiant

VALIANT_OBJS = Valiant.o PEnKF.o PEnOpt.o

Valiant:  ${VALIANT_OBJS} chkopts
	-${CLINKER} -g3 -O0 -o Valiant.exe ${VALIANT_OBJS} ${PETSC_SNES_LIB}
	${RM} Defiant.o

	
