#
#       INSTALL Makefile for DGELDA
#       Version 1.1.1
#       October 12, 1995
#
include ../make.inc

all:  testlsame testdlamch testdsecnd 

slamch.o: slamch.f ; $(FORTRAN) $(NOOPT) -c $<
dlamch.o: dlamch.f ; $(FORTRAN) $(NOOPT) -c $<
dsecnd.o: $(TIMESRC) ; $(TIMECOMP) $(TIMEOPTS) -c $(TIMESRC)

testlsame:  lsame.o lsametst.o
	$(LOADER) -o testlsame lsame.o lsametst.o

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(LOADER) -o testdlamch dlamch.o lsame.o dlamchtst.o 

testdsecnd: dsecnd.o dsecndtst.o
	$(LOADER) -o testdsecnd dsecnd.o dsecndtst.o $(TIMELIB)

.f.o: ; $(FORTRAN) $(OPTS) -c $<

clean:
	rm -f *.o
	rm -f *.log *.aux *.toc

cleanup:
	$(MAKE) clean
	rm -f testdlamch testdsecnd testlsame lawn81.dvi
	cp dsecnd.f.org dsecnd.f
	cp dsecnd.f.org ../LAPACK/dsecnd.f


