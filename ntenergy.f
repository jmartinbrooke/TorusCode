
      REAL*8 FUNCTION TEN(TPLUS,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS)     
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 TPLUS(NI,NJ),HSN(NI),HCS(NI),SN(NJ),CS(NJ)
      energ = 0.0d0
      ebd = dexp(bd)
      DO 30 I= 2,NI-1
        X1 = (I-1)*DX
        X2 = X1 + DX
        HSN1 = HSN(I)
        HSN2 = HSN(I+1)
        HCS1 = HCS(I)
        HCS2 = HCS(I+1)
        DO 40 J=1,NJ-1
          TH1 = J*DTH
          TH2 = TH1 + DTH 
          cs1 = CS(J)
          cs2 = CS(J+1)
          sn1 = SN(J)
          sn2 = SN(J+1)
          c11 = hcs1 - cs1
          c12 = hcs1 - cs2
          c21 = hcs2 - cs1
          c22 = hcs2 - cs2
          bsq11 = tplus(i,j)*tplus(i,j) 
          bsq12 = tplus(i,j+1)*tplus(i,j+1) 
          bsq21 = tplus(i+1,j)*tplus(i+1,j) 
          bsq22 = tplus(i+1,j+1)*tplus(i+1,j+1) 
          deng = hsn1*bsq11/c11**3/x1 + hsn1*bsq21/c21**3/x2 +
     *           hsn1*bsq12/c12**3/x1 + hsn2*bsq22/c22**3/x2
          deng = deng*pi*dx*dth
          energ = energ + deng
c*****
c***** end of main j-loop ************************************
40      continue
        J=NJ
          cs1 = CS(J)
          cs2 = CS(1)
          sn1 = SN(J)
          sn2 = SN(1)
          c11 = hcs1 - cs1
          c12 = hcs1 - cs2
          c21 = hcs2 - cs1
          c22 = hcs2 - cs2
          bsq11 = tplus(i,j)*tplus(i,j) 
          bsq12 = tplus(i,1)*tplus(i,1) 
          bsq21 = tplus(i+1,j)*tplus(i+1,j) 
          bsq22 = tplus(i+1,1)*tplus(i+1,1) 
          deng = hsn1*bsq11/c11**3/x1 + hsn1*bsq21/c21**3/x2 +
     *           hsn1*bsq12/c12**3/x1 + hsn2*bsq22/c22**3/x2
          deng = deng*pi*dx*dth
          energ = energ + deng
c*****
c***** end of i loop ******************************************
30    continue
      ten = 0.5d0*energ
      return
      end
