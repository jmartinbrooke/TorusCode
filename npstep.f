
c***** solves the equations for the next time step ********
c***** subroutine to time-step the poloidal equation ***     
    
      SUBROUTINE NPSTEP(PCOEF1,PCOEF2,PCOEF3,PCOEF4,PCOEF5,
     +                 PTCOF5,PMINUS,PZERO,PPLUS,TZERO,
     +                 DIFFA,UNDER,
     +                 ALPHA0,ALPHA,ALAS2,BOUND,
     +                 DX,DTH,DTIME,XSTART,NI,NJ,INLIN,
     +                 HSN,HCS,SN,CS) 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PMINUS(NI,NJ),PZERO(NI,NJ),PPLUS(NI,NJ),
     +          TZERO(NI,NJ),
     +          PCOEF1(NI,NJ),PCOEF2(NI,NJ),PCOEF3(NI,NJ),
     +          PCOEF4(NI,NJ),PCOEF5(NI,NJ),PTCOF5(NI,NJ),
     +          BOUND(NJ,2*NJ),ALAS2(2*NJ),ALPHA0(NJ),ALPHA(NI,NJ),
     +          DIFFA(NI,NJ),UNDER(NI,NJ),
     +          HSN(NI),HCS(NI),SN(NJ),CS(NJ)
      COMMON /PAR/PI,BD,CALPHA,COMEGA,A 

      DO 20 I=2,NI-1
        X=DX*(I-1)
        DO 30 J=2,NJ-1
          TH = (J-1)*DTH
          C = HCS(I) - CS(J)
          BN = -PZERO(I,J)*SN(J) + C*(PZERO(I,J+1)-PZERO(I,J-1))
     +          *0.5/DTH
          BTH = -PZERO(I,J)*(1-HCS(I)*CS(J))/HSN(I)  + 
     +          X*C*(PZERO(I+1,J)-PZERO(I-1,J))*0.5/DX 
          BSQ = BN*BN + BTH*BTH + TZERO(I,J)**2 
          ALPHA(I,J) = ALPHA0(J)/(1+INLIN*BSQ)
30      CONTINUE
      ALPHA(I,1) = 0.0D0
        J=NJ
        TH = (J-1)*DTH
        C = HCS(I) - CS(J)
        BN = -PZERO(I,J)*SN(J) + C*(PZERO(I,1)-PZERO(I,J-1))
     +        *0.5/DTH
        BTH = -PZERO(I,J)*(1-HCS(I)*CS(J))/HSN(I)  + 
     +        X*C*(PZERO(I+1,J)-PZERO(I-1,J))*0.5/DX 
        BSQ = BN*BN + BTH*BTH + TZERO(I,J)**2 
        ALPHA(I,J) = ALPHA0(J)/(1+INLIN*BSQ)
20    CONTINUE

      AV = 0.0D0
      DO 1 J=1,NJ
        AV = AV + ALPHA(2,J) 
1     CONTINUE
      AV = AV/NJ

      DO 2 J=1,NJ
        ALPHA(1,J) = AV 
2     CONTINUE
 
      DO 75 I=2,NI-1
        DO 85 J=2,NJ-1
           PPLUS(I,J) = PCOEF1(I,J)*PZERO(I-1,J)
     +                  + PCOEF2(I,J)*PZERO(I,J-1) 
     +                  + PCOEF3(I,J)*PZERO(I,J+1)
     +                  + PCOEF4(I,J)*PZERO(I+1,J) 
     +                  + PCOEF5(I,J)*PMINUS(I,J)
     +                  + PTCOF5(I,J)*TZERO(I,J)*ALPHA(I,J)
c***** to calculate the diffusion operator on a *****
        DIFFA(I,J) = (PPLUS(I,J) - PMINUS(I,J) -
     +                ALPHA(I,J)*TZERO(I,J)*2.0D0*DTIME)
     +               *UNDER(I,J)
c***** end of j-loop **********
   85   CONTINUE
        J=1
        PPLUS(I,J) = PCOEF1(I,J)*PZERO(I-1,J)
     +                  + PCOEF2(I,J)*PZERO(I,NJ) 
     +                  + PCOEF3(I,J)*PZERO(I,J+1)
     +                  + PCOEF4(I,J)*PZERO(I+1,J) 
     +                  + PCOEF5(I,J)*PMINUS(I,J)
     +                  + PTCOF5(I,J)*TZERO(I,J)*ALPHA(I,J)
c***** to calculate the diffusion operator on a *****
        DIFFA(I,J) = (PPLUS(I,J) - PMINUS(I,J) -
     +                ALPHA(I,J)*TZERO(I,J)*2.0D0*DTIME)
     +               *UNDER(I,J)
        J=NJ
         PPLUS(I,J) = PCOEF1(I,J)*PZERO(I-1,J)
     +                + PCOEF2(I,J)*PZERO(I,J-1) 
     +                + PCOEF3(I,J)*PZERO(I,1)
     +                + PCOEF4(I,J)*PZERO(I+1,J) 
     +                + PCOEF5(I,J)*PMINUS(I,J)
     +                + PTCOF5(I,J)*TZERO(I,J)*ALPHA(I,J)
c***** to calculate the diffusion operator on a *****
        DIFFA(I,J) = (PPLUS(I,J) - PMINUS(I,J) -
     +                ALPHA(I,J)*TZERO(I,J)*2.0D0*DTIME)
     +               *UNDER(I,J)


c***** end of i loop ******
   75 CONTINUE

C******* surface boundary condition **********************

C******* surface boundary condition **********************
C******  poloidal field read in jepps matrix ************* 
c******* second row in *************************************
      DO 21 J = 1, NJ
        ALAS2(J) = PPLUS(NI-2,J)
   21 CONTINUE
C******* first row in **************************************
      DO 22 J = 1, NJ
        ALAS2(J+NJ) = PPLUS(NI-1,J)
   22 CONTINUE
c******* calculation of boundary row ***********************
       DO 24 I=1,NJ
         PPLUS(NI,I) = 0
         DO 25 J=1,2*NJ
           PPLUS(NI,I) = PPLUS(NI,I) +BOUND(I,J)*ALAS2(J)
   25    CONTINUE
   24  CONTINUE
c******* end of surface boundary IF **************************
c*******
c******* calculate boundary condition at centre *********************
c*******
      AV = 0.0D0
      DO 10 J = 1,NJ
        AV = AV + PPLUS(2,J)
   10 CONTINUE
      AV = AV/DFLOAT(NJ)
      DO 11 J = 1,NJ
        PPLUS(1,J) = AV
   11 CONTINUE

c***** return to main *****
      RETURN

      END
c***** end of PSTEP *****************


