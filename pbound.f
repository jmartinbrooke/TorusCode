C**** program pbound to read pzero and produce     ** 
c**** pbound.dat ************************************
      PROGRAM PBOUND
      PARAMETER (NI=41,NJ=160)
      REAL*8 MATRI1(NI,NJ),bound(2*nj)
      open (1,file='pzero.bin',form='unformatted',status='unknown')
      open  (10,file='pbound.dat')
C     open  (20,file='outer.dat')

      READ (1) MATRI1
      close(1)

       do 1 j=1,nj
1         write(10,*) MATRI1(NI-2,J)
       do 2 j=1,nj
2         write(10,*) MATRI1(NI-1,J)
      close(10)
C     open  (10,file='pbound.dat')
C     do 3 j=1,nj
C3     READ(10,1000) BOUND(J)
C     do 4 j=1,nj
C4     READ(10,1000) BOUND(J+NJ)
          

      STOP
99999 FORMAT(D10.3)
1000  FORMAT(d20.13)
      END

