FFLAGS = -O2 -I$(REG_STEER_HOME)/include 
LDFLAG = -L$(REG_STEER_HOME)/lib32 -L/usr/lib -lReG_Steer -lxml2 -lm -lz
COMPILER = g95
alnpar : alnpar.o npcoff.o ntcoff.o npstep.o ntstep.o npenergy.o ntenergy.o \
         eig.o nparity.o 
	$(COMPILER) $(FFLAGS) -o alnpar alnpar.o npcoff.o ntcoff.o ntstep.o npstep.o \
                  npenergy.o ntenergy.o eig.o nparity.o $(LDFLAG)
nalnpar : nalnpar.o npcoff.o ntcoff.o npstep.o ntstep.o npenergy.o ntenergy.o \
         eig.o nparity.o 
	$(COMPILER) $(FFLAGS) -o nalnpar nalnpar.o npcoff.o ntcoff.o ntstep.o npstep.o \
                  npenergy.o ntenergy.o eig.o nparity.o $(LDFLAG)
surface : surface.f
	$(COMPILER) $(FFLAGS) -o surface surface.f -lnag

tsurface : tsurface.f
	$(COMPILER) $(FFLAGS) -o tsurface tsurface.f -lnag

vacuum : vacuum.f conv.f
	$(COMPILER) -o vacuum vacuum.f conv.f -lnag

frmvac : frmvac.f conv.f
	$(COMPILER) -o frmvac frmvac.f conv.f -lnag

vpxvac : vpxvac.f conv.f
	$(COMPILER) -o vpxvac vpxvac.f conv.f -lnag

parset : parset.f npenergy.o ntenergy.o
	$(COMPILER) $(FFLAGS) -o parset parset.f npenergy.o ntenergy.o

maxmin : dmax_find.o max_find.o
	$(COMPILER) -o maxmin dmax_find.o max_find.o

andysymm : andysymm.f nturn_find.f
	$(COMPILER) -o andysymm +E6 andysymm.f nturn_find.f
alquen.o : alquen.f
	$(COMPILER) -c $(FFLAGS) alquen.f
alnpar.o : alnpar.f
	$(COMPILER) -c $(FFLAGS) alnpar.f
pcoeff.o : pcoeff.f
	$(COMPILER) -c $(FFLAGS) pcoeff.f
npcoff.o : npcoff.f
	$(COMPILER) -c $(FFLAGS) npcoff.f
tcoeff.o : tcoeff.f
	$(COMPILER) -c $(FFLAGS) tcoeff.f
ntcoff.o : ntcoff.f
	$(COMPILER) -c $(FFLAGS) ntcoff.f
pstep.o : pstep.f
	$(COMPILER) -c $(FFLAGS) pstep.f
npstep.o : npstep.f
	$(COMPILER) -c $(FFLAGS) npstep.f
tstep.o : tstep.f
	$(COMPILER) -c $(FFLAGS) tstep.f
ntstep.o : ntstep.f
	$(COMPILER) -c $(FFLAGS) ntstep.f
tstep1.o : tstep1.f
	$(COMPILER) -c $(FFLAGS) tstep1.f
penergy.o : penergy.f
	$(COMPILER) -c $(FFLAGS) penergy.f
npenergy.o : npenergy.f
	$(COMPILER) -c $(FFLAGS) npenergy.f
tenergy.o : tenergy.f
	$(COMPILER) -c $(FFLAGS) tenergy.f
ntenergy.o : ntenergy.f
	$(COMPILER) -c $(FFLAGS) ntenergy.f
eig.o : eig.f
	$(COMPILER) -c $(FFLAGS) eig.f

nparity.o : nparity.f
	$(COMPILER) -c $(FFLAGS) nparity.f

clean :
	rm alnpar  *gr *log ftn* *.upi *.o

