

c***** call routine to calculate coeffs of the toroidal eqn ***
    
      SUBROUTINE NTCOFF(TCOEF1,TCOEF2,TCOEF3,TCOEF4,TCOEF5,
     +                  TPCOF1,TPCOF2,TPCOF3,TPCOF4,TPCOF5,
     +                  UNDER,ALPHA0,
     +                  NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +                  HSN,HCS,SN,CS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TCOEF1(NI,NJ),TCOEF2(NI,NJ),TCOEF3(NI,NJ),
     +          TCOEF4(NI,NJ),TCOEF5(NI,NJ),
     +          TPCOF1(NI,NJ),TPCOF2(NI,NJ),TPCOF3(NI,NJ),
     +          TPCOF4(NI,NJ),TPCOF5(NI,NJ),
     +          UNDER(NI,NJ),ALPHA0(NJ),
     +          HSN(NI),HCS(NI),SN(NJ),CS(NJ) 
      COMMON /PAR/PI,BD,CALPHA,COMEGA,A

      Q1 = DTIME/ (DX**2)
      Q2 = DTIME/ (DTH**2)
      Q3 = DTIME/(2.0*DX)
      Q4 = DTIME/(2.0*DTH)
  
c***** start of i-loop on interior points ********
 
      DO 20 I = 2,NI - 1
          X = (I-1)*DX + XSTART
          HSN3=HSN(I)**3
          HSN32=DSQRT(HSN3)
          DO 30 J = 1,NJ
              TH = (J-1)*DTH 
              C = HCS(I) - CS(J)
              C2=C**2
              C3=C*C2
              C32=DSQRT(C3)
              HSN3=HSN(I)**3
              HSN32=DSQRT(HSN3)
              MULT=-1.5D0*COMEGA*C32/HSN32
C             MULT=-1.0D0*COMEGA*HSN/C
              AX=C2*X**2
              BX=X*C2-X*C2*HCS(I)/HSN(I)+X*C*HSN(I)
              PBX=-X*MULT*SN(J)*HSN(I)
              ATH=C2
              BTH=-C*SN(J)
              PBTH=MULT*(-HCS(I)*CS(J) + 1.0D0)
              CXTH=-C2/HSN(I)**2
  
c***** determination of decay coefficients ******
  
           DENOM = 1.0/ (1+2.0*Q1*AX+2.0*Q2*ATH-0.5*CXTH*DTIME)
           UNDER(I,J) = DENOM
           TCOEF1(I,J) = 2.0*(Q1*AX-BX*Q3)*DENOM
           TCOEF2(I,J) = 2.0* (Q2*ATH-BTH*Q4)*DENOM
           TCOEF3(I,J) = 2.0* (Q2*ATH+BTH*Q4)*DENOM
           TCOEF4(I,J) = 2.0*(Q1*AX + BX*Q3)*DENOM
           TCOEF5(I,J) = (1.0D0-2.0*Q1*AX-2.0*Q2*ATH+0.5*CXTH*DTIME)
     +                   *DENOM
           TPCOF1(I,J) = -2.0*DENOM*Q3*PBX
           TPCOF2(I,J) = -2.0*DENOM*Q4*(PBTH)
           TPCOF3(I,J) = 2.0*DENOM*Q4*(PBTH)
           TPCOF4(I,J) = 2.0*DENOM*Q3*PBX
           TPCOF5(I,J) =  +2.0*DTIME*DENOM 
           IF (I .EQ. 2) THEN
             ALPHA0(J) = CALPHA*SN(J)
           ENDIF

c***** end of j loop *****
   30    CONTINUE
   20    CONTINUE

      
      RETURN

      END
c***** end of subroutine TCOEFF ****
    
