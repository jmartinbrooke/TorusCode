
 
c***** subroutine to calculate coeffs of the poloidal eqn ***
    
      SUBROUTINE NPCOFF(PCOEF1,PCOEF2,PCOEF3,PCOEF4,PCOEF5,
     +               PTCOF5,NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +               IEQ,SN,CS,HSN,HCS)
    
      DOUBLE PRECISION PCOEF1(NI,NJ),PCOEF2(NI,NJ),PCOEF3(NI,NJ),
     +          PCOEF4(NI,NJ),PCOEF5(NI,NJ),PTCOF5(NI,NJ),
     +          SN(NJ),CS(NJ),HSN(NI),HCS(NI)

      DOUBLE PRECISION Q1,Q2,Q3,Q4,X,TH,DX,DTH,C,C2,AX,BX,ATH,BTH,CXTH,
     +DENOM,DTIME,XSTART,THSTAR
      COMMON /PAR/PI,BD,CALPHA,COMEGA,A

      Q1 = DTIME/ (DX**2)
      Q2 = DTIME/ (DTH**2)
      Q3 = DTIME/(2.0*DX)
      Q4 = DTIME/(2.0*DTH)
  
c***** start of i-loop on interior points ********
 
      DO 20 I = 2,NI - 1
          X = (I-1)*DX + XSTART
          DO 30 J = 1,NJ 
              TH = (J-1)*DTH 
C*****
C******** setting up coefficients  ***********
C******
               C = HCS(I) - CS(J)
               C2=C**2                                       
               AX=C2*X**2
               BX=X*C2-X*C2*HCS(I)/HSN(I)+X*C*HSN(I)
               ATH=C2
               BTH=-C*SN(J)
               CXTH=-C2/HSN(I)**2
  
c***** determination of coefficients for dufort-fraenkel solver **
  
              DENOM = 1.0/ (1+2.0*Q1*AX+2.0*Q2*ATH-0.5*CXTH*DTIME)
              PCOEF1(I,J) = 2.0*(Q1*AX-BX*Q3)*DENOM
              PCOEF2(I,J) = 2.0* (Q2*ATH-BTH*Q4)*DENOM
              PCOEF3(I,J) = 2.0* (Q2*ATH+BTH*Q4)*DENOM
              PCOEF4(I,J) = 2.0*(Q1*AX + BX*Q3)*DENOM
              PCOEF5(I,J) = (1.0D0-2.0*Q1*AX-2.0*Q2*ATH+0.5*CXTH*DTIME)
     +                 *DENOM
              PTCOF5(I,J) = (2.0D0*DTIME)*DENOM
 
c***** 30 marks the end of the main j-sweep *************
   30     CONTINUE
   20 CONTINUE
c***** return to main ****
      RETURN
      END
c****    end of subroutine PCOEFF *******

