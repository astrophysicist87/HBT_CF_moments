SHELL=/bin/sh

SRCS= \
main.cpp \
gauss/gauss_quadrature.cpp \
lib.cpp \
HBT_at_K.cpp

HDRS= \
gauss/gauss_quadrature.h \
HBT_at_K.h \
lib.h \
main.h \
parameters.h

MAKEFILE=makefile

COMMAND=run_HBT_CF_moments

OBJS= $(addsuffix .o, $(basename $(SRCS)))

CC= g++
CFLAGS= -g -O3 -fopenmp

WARNFLAGS=
LDFLAGS= -lgsl -lgslcblas
LIBS= -L/sw/lib -I/sw/include

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(COMMAND): $(OBJS) $(HDRS) $(MAKEFILE) 
	$(CC) -o $(COMMAND) $(OBJS) $(LDFLAGS) $(CFLAGS) $(LIBS)

clean:
	rm -f $(OBJS) $(COMMAND)
 
tarz:
	tar zcf - $(MAKEFILE) $(SRCS) $(HDRS) > $(COMMAND).tar.gz
 
