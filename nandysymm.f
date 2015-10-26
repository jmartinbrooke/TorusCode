      PROGRAM dturn
      INTEGER ntime , nmax, ilag, iturn, npts, nturn  
      PARAMETER (nmax = 10000)
      DOUBLE PRECISION a_n(nmax), a_s(nmax), b_n(nmax),b_s(nmax),
     +b_nturn(nmax/2), b_sturn(nmax/2), time(nmax),
     +a_nturn(nmax/2),tturn(nmax/2),a_sturn(nmax/2),
     +period(2:nmax/2)
      DOUBLE PRECISION pi, bnabs(2), bsabs(2), anabs(2), asabs(2),dtime
      DOUBLE PRECISION avper
      CHARACTER*20 filename


c**** the value of dtime for Andrew's stuff ***
      dtime = 0.0005
      ilag = 4
      iphase = 1
      PRINT *, 'Enter file name to be read'
      READ (*,'(A20)') filename
      
      open (1,FILE=filename)
      open (2,FILE="nandysymm.dat")
      open (3,FILE="nandyval.dat")
C     READ (1,*,END=10) time(i), a_n(i), a_s(i), b_n(i), b_s(i)
      npts = 0
      DO i = 1, nmax
      READ (1,*,END=10) time(i), a_n(i), a_s(i), b_n(i), b_s(i)
      npts=npts+1
      time(i) = i*dtime
      END DO
10    PRINT *,'End of reading file ',npts
      CALL turn(b_n,time,dtime,b_nturn,tturn,period,npts,nturn)
      CALL interp(a_n,time,dtime,a_nturn,tturn,npts,nturn)
      CALL interp(b_s,time,dtime,b_sturn,tturn,npts,nturn)
      CALL interp(a_s,time,dtime,a_sturn,tturn,npts,nturn)
      DO i = 1,nturn
      write (3,'(5e15.8)') b_nturn(i),b_sturn(i),a_nturn(i),a_sturn(i), 
     *         tturn(i)
      END DO
      bnabs(1) = DABS(b_nturn(iphase))
      bsabs(1) = DABS(b_sturn(iphase))
      anabs(1) = DABS(a_nturn(iphase))
      asabs(1) = DABS(a_sturn(iphase))
      avper = 0.
**** in this loop we write the same-time symmetries as n,
c**** the forward symmetries as s, thus we overwrite the original
c**** values
      DO i = ilag+iphase, nturn, ilag 
      bnabs(2) = DABS(b_nturn(i))
      bsabs(2) = DABS(b_sturn(i))
      anabs(2) = DABS(a_nturn(i))
      asabs(2) = DABS(a_sturn(i))
      b_n(i) = 0.5*(b_nturn(i) - b_sturn(i))
     * /DMAX1(bnabs(2),bsabs(2))
      b_s(i) = 0.5*(b_nturn(i-ilag) - b_sturn(i))
     * /  DMAX1(bnabs(1),bsabs(2))
      a_n(i) = 0.5*(a_nturn(i) - a_sturn(i))
     * / DMAX1(anabs(2),asabs(2))
      a_s(i) = 0.5*(a_nturn(i-ilag) - a_sturn(i))
     * /DMAX1(anabs(1),asabs(2))
      avper = avper + period(i)
      bnabs(1) = bnabs(2)
      bsabs(1) = bsabs(2)
      anabs(1) = anabs(2)
      asabs(1) = asabs(2)
      END DO 

      avper = ilag*avper/nturn
      write(*,*) 'avper = ',avper

      WRITE(2,'("    b_inst         b_forw         a_inst",
     *"        a_forw         time")')
      DO i = ilag + iphase, nturn -ilag, ilag
      WRITE(2,'(5e15.8)') b_n(i),b_s(i),a_n(i),a_s(i),tturn(i)
      END DO 
      STOP
      END
