c**** subroutine to calculate the symmetric and antisymmetric components
c**** of the energy by decomposing the field into these symmetries.
      SUBROUTINE PARITY(PPLUS,TPLUS,APPLUS,SPPLUS,ATPLUS,STPLUS,
     +ETOT,RP2T,PAR,TIME,
     +NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS,ESYMM,EANTI) 
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PPLUS(NI,NJ),TPLUS(NI,NJ),APPLUS(NI,NJ),SPPLUS(NI,NJ),
     +ATPLUS(NI,NJ),STPLUS(NI,NJ)
      EXTERNAL PEN,TEN
c***** start of executable statements **************************
c***** a refers to anti, s to symm ******************************
c***** returns rpst - ratio of poloidal to toroidal energies ***
c***** etot - total energy **************************************
c***** parity = (esymm - eanti)/(esymm + eanti) *****************
c***** first of all we decompose the field into symmetric and antisymmetric
c***** components
      NJHALF = NJ/2
c***** first of all we do the interior mesh points, over half the theta range
      DO I=2,NI
        DO J=2,NJHALF
          APPLUS(I,J) = 0.5D0*(PPLUS(I,J) + PPLUS(I,NJ-J+2))
          ATPLUS(I,J) = 0.5D0*(TPLUS(I,J) - TPLUS(I,NJ-J+2))
          SPPLUS(I,J) = 0.5D0*(PPLUS(I,J) - PPLUS(I,NJ-J+2))
          STPLUS(I,J) = 0.5D0*(TPLUS(I,J) + TPLUS(I,NJ-J+2))
        END DO
        DO J=2,NJHALF
          APPLUS(I,NJ-J+2) = APPLUS(I,J)
          ATPLUS(I,NJ-J+2) = -ATPLUS(I,J)
          SPPLUS(I,NJ-J+2) = -SPPLUS(I,J)
          STPLUS(I,NJ-J+2) = +STPLUS(I,J)
        END DO
c***** now we do the points at the boundaries of the theta sweep treating
c***** the poles and the equator as special cases.
        APPLUS(I,1) = PPLUS(I,1)
        ATPLUS(I,1) = 0.0D0
        APPLUS(I,NJHALF+1) = PPLUS(I,NJHALF+1)
        ATPLUS(I,NJHALF+1) = 0.0D0
        SPPLUS(I,1) = 0.0D0
        STPLUS(I,1) = TPLUS(I,1)
        SPPLUS(I,NJHALF+1) = 0.0D0
        STPLUS(I,NJHALF+1) = TPLUS(I,NJHALF+1)
      END DO
c**** now we do the values at the centre of the torus.      
      DO J=1,NJ
        APPLUS(1,J) = PPLUS(1,NJ)
        ATPLUS(1,J) = 0.0D0
        SPPLUS(1,J) = 0.0D0
        STPLUS(1,J) = TPLUS(1,NJ)
      END DO
 
c***** PEN, TEN are functions to calulate the energy of the scalar field, 
c***** they are in npenergy.f and ntenergy.f
      EAP = PEN(APPLUS,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS) 
      ESP = PEN(SPPLUS,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS) 
      EAT = TEN(ATPLUS,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS)     
      EST = TEN(STPLUS,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS)     

c**** calculate symmetric and antisymmetric fields also poloidal and toroidal
c**** fields. 

      ESYMM = ESP + EST
      EANTI = EAP + EAT 
      EPOL = EAP + ESP
      ETOR = EAT + EST

      IF (ETOR .GE. 1D-60) THEN
        RP2T = EPOL/ETOR
      ENDIF

      ETOT = EPOL + ETOR
      
      ES = ESP + EST
      EA = EAP + EAT

      IF ((ES + EA) .GE. 1D-60) THEN
        PAR = (ES - EA)/(ES + EA)
      ELSE 
        PAR = 2.0D0
      ENDIF

      RETURN
      END
        
