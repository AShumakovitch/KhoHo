CC = gcc
CFLAGS = -g -O3 -Wall -fPIC

# uncomment the following line if PARI/GP is not installed in one of the
# standard places; edit the path to PARI's header files appropriately
# PARI_INPUT = -I/path/to/my/copy/of/pari

UNAME := ${shell uname}
ifeq (${UNAME}, Darwin)  # Mac OS X
	LDFLAGS = -flat_namespace -bundle -undefined suppress
	STRIP = true   # Does nothing 
else   # Linux 
	LDFLAGS = -shared
	STRIP = strip -p ${SH_OBJ}
endif

SH_OBJ = print_ranks.so nicematr.so sparreduce.so

SPARSE_MAT_LIB = sparmat.o
sparreduce_EXTRA_LIBS = ${SPARSE_MAT_LIB}

%.o: %.c
	${CC} ${CFLAGS} ${PARI_INPUT} -c $< -o $@

%.so: %.o
	${CC} ${LDFLAGS} $< ${$*_EXTRA_LIBS} -o $@

all: binary strip

binary: ${SH_OBJ}

strip:
	${STRIP} ${SH_OBJ}

sparreduce.so: sparmat.c sparmat.h ${SPARSE_MAT_LIB} 

sparmat.o: sparmat.h

clean:
	rm -f ${SH_OBJ} ${SH_OBJ:.so=.o} ${SPARSE_MAT_LIB}

.PHONY: all binary strip clean
