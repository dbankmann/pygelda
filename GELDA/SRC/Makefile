#
#	SRC Makefile for DGELDA
#	Version 1.1
#       July 19, 1995
# 
include ../make.inc

OBJS = dbdstp.o dbdtrp.o dbinom.o dgeerm.o dgeini.o dgelda.o dgenrm.o \
       dnfind.o dnfini.o dnfix1.o dnfotr.o dnfpen.o dnfred.o \
       dnfval.o drkstp.o drkrls.o drkcxs.o 

all: ../dgelda.a

../dgelda.a: $(OBJS)
	$(ARCH) $(ARCHFLAGS) ../dgelda.a $(OBJS)
	$(RANLIB) ../dgelda.a

.f.o:
	$(FORTRAN) $(FCOPTS) -c $<
	
clean:
	rm -f *.o 
