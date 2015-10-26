c***** Actual arguments are: 
c***** energy(nmax) array of values whose maxima and minima we need to find
c***** time(nmax) array of times for the values of energy
c***** eturn(nmax/2) the array of turning values, interpolated if neccessary
c***** tturn(nmax/2) the array of turning times, ...
c***** iturn counts the number of maxima and minima
c***** IMPORTANT Please note that we record maxima AND minima and do not
c***** distinguish between them
      SUBROUTINE turn(energy,time,dtime,eturn,tturn,period,ncount,nturn)

      PARAMETER (nmax = 10000)
      DOUBLE PRECISION energy(nmax),time(nmax)
      DOUBLE PRECISION eturn(nmax/2),tturn(nmax/2),period(2:nmax/2)
      DOUBLE PRECISION  oldslp, newslp
      INTEGER ncount,iturn,nturn    


      oldslp = (energy(3) - energy(1))*0.5/dtime
      iturn = 0
c***** loop over the array supplied, from 2 to nmax
c***** we need at least three values to determine a slope
      DO i = 3, ncount
       newslp = (energy(i)  - energy(i-2))*0.5/dtime
      
c***** if turn is <= 0 we have a turning point 
       turnpt = oldslp*newslp
c***** have we found an exact turning point?
c     IF (DABS(newslp) .LT. 1e-12) THEN 
      IF (DABS(newslp) .GT. 1e+12) THEN 
        iturn = iturn + 1
c***** the point to which the derivative belongs is midway between
c***** the points used to calculate newslp, ie i and i-2
        tturn(iturn) = time(i-1)
        eturn(iturn) = energy(i-1)
        IF (iturn .NE. 1) period(iturn) = tturn(iturn) - tturn(iturn-1)
        PRINT *, 'zero slope encountered'
c*****  next line is to ensure that we record sign of lhs slope
        newslp = -oldslp
c***** we have ended the code needed for an exact turning point
c***** we now move onto the code for a crossing of the zero line
c***** interpolation needed here to find the turning point
      ELSE IF (turnpt .LT. 0) THEN
c***** its a turning point
        iturn = iturn + 1
        dabso = DABS(oldslp)
        dabsn = DABS(newslp)
        IF ((dabso + dabsn) .EQ. 0) THEN
          write (*,*) 'error at turn', iturn
          RETURN
        END IF
        eturn(iturn) = (energy(i-2)*dabsn + energy(i-1)*dabso)
     +  /(dabso + dabsn)
          
        tturn(iturn) = (time(i-2)*dabsn + time(i-1)*dabso)
     +  /(dabso + dabsn)
        IF (iturn .NE. 1) period(iturn) = tturn(iturn) 
     +   - tturn(iturn -1)
c**** end of outer if block 
        PRINT*,'not an exact turning point'
      END IF  

c**** reset oldslp ****
      oldslp = newslp
      
c***** end of the do loop
      END DO
      nturn = iturn
    
      RETURN
      END 
       
       
c***** this subroutine takes the times of turning values as 
c***** given by maxget in the array tturn(nmax/2)
c**** and a timeseries given by val and returns the values of 
c***** val in the array valturn
      SUBROUTINE interp(val,time,dtime,valturn,tturn,npts,nturn)
      INTEGER npts, nturn
      PARAMETER (nmax = 10000)
      DOUBLE PRECISION val(nmax),time(nmax),tturn(nmax/2),
     +valturn(nmax/2)
      DOUBLE PRECISION dtime
      INTEGER i, j

      j = 1
      DO i = 1, npts
        IF (time(i) .GT. tturn(j)) THEN
         valturn(j) = (val(i-1)*(time(i)-tturn(j)) + 
     +   val(i)*(tturn(j) - time(i-1))) / dtime
         j = j + 1
        END IF
        IF (j .EQ. nturn) GO TO 100
      END DO
100   RETURN
      END
