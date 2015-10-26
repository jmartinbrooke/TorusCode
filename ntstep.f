

c***** subroutine to time-step the toroidal equation ****    

      SUBROUTINE NTSTEP(TCOEF1,TCOEF2,TCOEF3,TCOEF4,TCOEF5,
     +                 TPCOF1,TPCOF2,TPCOF3,TPCOF4,TPCOF5,
     +                 TMINUS,TZERO,TPLUS,PZERO,
     +                 DX,DTH,DTIME,NI,NJ,XSTART,
     +                 INLIN,DIFFA,UNDER,AL0,ALPHA,
     +                 HSN,HCS,SN,CS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TMINUS(NI,NJ),TZERO(NI,NJ),TPLUS(NI,NJ),
     +          TCOEF1(NI,NJ),TCOEF2(NI,NJ),TCOEF3(NI,NJ),
     +          TCOEF4(NI,NJ),TCOEF5(NI,NJ),
     +          TPCOF1(NI,NJ),TPCOF2(NI,NJ),TPCOF3(NI,NJ),
     +          TPCOF4(NI,NJ),TPCOF5(NI,NJ),PZERO(NI,NJ),
     +          DIFFA(NI,NJ),UNDER(NI,NJ),
     +          ALPHA(NI,NJ),AL0(NI,NJ),
     +          HCS(NI),HSN(NI),SN(NJ),CS(NJ)
  
      COMMON /PAR/PI,BD,CALPHA,COMEGA,A

 
      DO 29 I=2,NI-1
        X = (I-1)*DX
        DO 39 J=2,NJ-1
          TH = (J-1)*DTH
          C = HCS(I) - CS(J)
          BN = -PZERO(I,J)*SN(J) + C*(PZERO(I,J+1)-PZERO(I,J-1))
     +          *0.5/DTH
          BTH = -PZERO(I,J)*(1-HCS(I)*CS(J))/HSN(I)  + 
     +          X*C*(PZERO(I+1,J)-PZERO(I-1,J))*0.5/DX 
          BSQ = BN*BN + BTH*BTH + TZERO(I,J)**2 
          TPLUS(I,J) = TCOEF1(I,J)*TZERO(I-1,J) 
     +                 + TCOEF2(I,J)*TZERO(I,J-1) 
     +                 + TCOEF3(I,J)*TZERO(I,J+1) 
     +                 + TCOEF4(I,J)*TZERO(I+1,J) 
     +                 + TCOEF5(I,J)*TMINUS(I,J)
     +                 + TPCOF1(I,J)*PZERO(I-1,J) 
     +                 + TPCOF2(I,J)*PZERO(I,J-1) 
     +                 + TPCOF3(I,J)*PZERO(I,J+1) 
     +                 + TPCOF4(I,J)*PZERO(I+1,J) 
     +                 - ALPHA(I,J)*DIFFA(I,J)
c***** non-linear term ****************************************
            ADETA = -X*(ALPHA(I+1,J)-ALPHA(I-1,J))*0.5/DX
            ADTH = (ALPHA(I,J+1)-ALPHA(I,J-1))*0.5/DTH
            TPLUS(I,J) = TPLUS(I,J) + (C*BTH*ADETA - C*BN*ADTH)
     +                  *UNDER(I,J)*2.0D0*DTIME
   39   CONTINUE

c**** j=1 ***************************************
         J = 1
          TH = (J-1)*DTH
          C = HCS(I) - CS(J)
          BN = -PZERO(I,J)*SN(J) + C*(PZERO(I,J+1)-PZERO(I,NJ))
     +          *0.5/DTH
          BTH = -PZERO(I,J)*(1-HCS(I)*CS(J))/HSN(I)  + 
     +          X*C*(PZERO(I+1,J)-PZERO(I-1,J))*0.5/DX 
          BSQ = BN*BN + BTH*BTH + TZERO(I,J)**2 
          TPLUS(I,J) = TCOEF1(I,J)*TZERO(I-1,J) 
     +                 + TCOEF2(I,J)*TZERO(I,NJ) 
     +                 + TCOEF3(I,J)*TZERO(I,J+1) 
     +                 + TCOEF4(I,J)*TZERO(I+1,J) 
     +                 + TCOEF5(I,J)*TMINUS(I,J)
     +                 + TPCOF1(I,J)*PZERO(I-1,J) 
     +                 + TPCOF2(I,J)*PZERO(I,NJ) 
     +                 + TPCOF3(I,J)*PZERO(I,J+1) 
     +                 + TPCOF4(I,J)*PZERO(I+1,J) 
     +                 - ALPHA(I,J)*DIFFA(I,J)
c***** non-linear term ****************************************
            ADETA = -X*(ALPHA(I+1,J)-ALPHA(I-1,J))*0.5/DX
            ADTH = (ALPHA(I,J+1)-ALPHA(I,NJ))*0.5/DTH
            TPLUS(I,J) = TPLUS(I,J) + (C*BTH*ADETA - C*BN*ADTH)
     +                  *UNDER(I,J)*2.0D0*DTIME

c***** j = nj ***************************************
          J = NJ
          TH = (J-1)*DTH
          C = HCS(I) - CS(J)
          BN = -PZERO(I,J)*SN(J) + C*(PZERO(I,1)-PZERO(I,J-1))
     +          *0.5/DTH
          BTH = -PZERO(I,J)*(1-HCS(I)*CS(J))/HSN(I)  + 
     +          X*C*(PZERO(I+1,J)-PZERO(I-1,J))*0.5/DX 
          BSQ = BN*BN + BTH*BTH + TZERO(I,J)**2 
          TPLUS(I,J) = TCOEF1(I,J)*TZERO(I-1,J) 
     +                 + TCOEF2(I,J)*TZERO(I,J-1) 
     +                 + TCOEF3(I,J)*TZERO(I,1) 
     +                 + TCOEF4(I,J)*TZERO(I+1,J) 
     +                 + TCOEF5(I,J)*TMINUS(I,J)
     +                 + TPCOF1(I,J)*PZERO(I-1,J) 
     +                 + TPCOF2(I,J)*PZERO(I,J-1) 
     +                 + TPCOF3(I,J)*PZERO(I,1) 
     +                 + TPCOF4(I,J)*PZERO(I+1,J) 
     +                 - ALPHA(I,J)*DIFFA(I,J)
c***** non-linear term ****************************************
            ADETA = -X*(ALPHA(I+1,J)-ALPHA(I-1,J))*0.5/DX
            ADTH = (ALPHA(I,1)-ALPHA(I,J-1))*0.5/DTH
            TPLUS(I,J) = TPLUS(I,J) + (C*BTH*ADETA - C*BN*ADTH)
     +                  *UNDER(I,J)*2.0D0*DTIME

c***** end of j-loop **********

c***** end of i-loop *********************************
   29 CONTINUE

C***** surface boundary condition **********************
C*****   toroidal field is zero at surface ***         
      DO 26 J=1,NJ
        TPLUS(NI,J) = 0
 26   CONTINUE

C*****
C***** calculate boundary condition at centre *********************
C*****
      AV = 0.0D0
      DO 10 J = 1,NJ
        AV = AV + TPLUS(2,J)
   10 CONTINUE
      AV = AV/DFLOAT(NJ)
      DO 11 J = 1,NJ
        TPLUS(1,J) = AV
   11 CONTINUE


c***** return to main *******
      RETURN

      END
c***** end of TSTEP *****

