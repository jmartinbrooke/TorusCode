      REAL*8 FUNCTION PEN(PPLUS,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS) 
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL *8 PPLUS(NI,NJ),HSN(NI),HCS(NI),SN(NJ),CS(NJ)
      ENERG = 0.0D0
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
          CS1 = CS(J)
          cs2 = CS(J+1)
          sn1 = SN(J)
          sn2 = SN(J+1)
          c11 = hcs1 - cs1
          c12 = hcs1 - cs2
          c21 = hcs2 - cs1
          c22 = hcs2 - cs2
          c21 = hcs2 -cs1
          if (i .eq. 2) then
            dads1 = (-1.5d0*pplus(i,j) + 2.0d0*pplus(i+1,j)
     *               -0.5D0*pplus(i+2,j)) / dx
            dads2 = (pplus(i+2,j) - pplus(i,j)) / (2.0*dx)
          else if (i .eq. ni-1) then
            dads1 = (pplus(i+1,j) - pplus(i-1,j)) / (2.0*dx)
            dads2 = (1.5d0*pplus(i+1,j) - 2.0d0*pplus(i,j)
     *               +0.5D0*pplus(i-1,j)) / dx
          else
            dads1 = (pplus(i+1,j) - pplus(i-1,j)) / (2.0*dx)
            dads2 = (pplus(i+2,j) - pplus(i,j)) / (2.0*dx)
          endif
          if (j .eq. 1) then
            dadth1 = (pplus(i,j+1) - pplus(i,nj)) / (2.0*dth)
            dadth2 = (pplus(i,j+2) - pplus(i,j)) / (2.0*dth)
          else if (j .eq. nj-1) then
            dadth1 = (pplus(i,j+1) - pplus(i,j-1)) / (2.0*dth)
            dadth2 = (pplus(i,1) - pplus(i,j)) / (2.0*dth)
          else
            dadth1 = (pplus(i,j+1) - pplus(i,j-1)) / (2.0*dth)
            dadth2 = (pplus(i,j+2) - pplus(i,j)) / (2.0*dth)
          endif
          bn11 = -pplus(i,j)*sn1 + c11*dadth1
          bn12 = -pplus(i,j+1)*sn2 + c12*dadth2
          bn21 = -pplus(i+1,j)*sn1 + c21*dadth1
          bn22 = -pplus(i+1,j+1)*sn2 + c22*dadth2
          bth11 = pplus(i,j)*(1.0d0 - hcs1*cs1)/hsn1 - x1*c11*dads1
          bth12 = pplus(i,j+1)*(1.0d0 - hcs1*cs2)/hsn1 - x1*c12*dads1
          bth21 = pplus(i+1,j)*(1.0d0 - hcs2*cs1)/hsn2 - x2*c21*dads2
          bth22 = pplus(i+1,j+1)*(1.0d0 - hcs2*cs2)/hsn2-x2*c22*dads2
          bsq11 = bn11*bn11 + bth11*bth11
          bsq12 = bn12*bn12 + bth12*bth12
          bsq21 = bn21*bn21 + bth21*bth21
          bsq22 = bn22*bn22 + bth22*bth22
          deng = hsn1*bsq11/c11**3/x1 + hsn2*bsq21/c21**3/x2 +
     *           hsn1*bsq12/c12**3/x1 + hsn2*bsq22/c22**3/x2
          deng = deng*pi*dx*dth
          energ = energ + deng
40      continue
30    continue
      pen = energ
      return
      end
