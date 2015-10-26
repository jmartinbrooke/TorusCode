c***** Actual arguments are: 
c***** energy(nmax) array of values whose maxima and minima we need to find
c***** time(nmax) array of times for the values of energy
c***** emax(nmax/2) the array of maximum values, interpolated if neccessary
c***** tmax(nmax/2) the array of maximimun times, ...
c***** emin(nmax/2), tmin(nmax/2) arrays for minimum values
c***** imax counts the number of maxima, imin counts minima
      SUBROUTINE maxget(energy,time,dtime,emax,tmax,emin,tmin,period,
     + imax,imin,ncount)

      PARAMETER (nmax = 10000)
      DOUBLE PRECISION energy(nmax),time(nmax)
      DOUBLE PRECISION emax(nmax/2),tmax(nmax/2),period(2:nmax/2)
      DOUBLE PRECISION emin(nmax/2),tmin(nmax/2)
      DOUBLE PRECISION  oldslp, newslp
      LOGICAL lmax        
      INTEGER ncount,imax, imin


      oldslp = (energy(3) - energy(1))*0.5/dtime
      imax = 0
      IF (oldslp .LT. 0) THEN
        lmax = .FALSE.
      ELSE
        lmax = .TRUE.
      END IF
        
c***** loop over the array supplied, from 2 to nmax
c***** we need at least three values to determine a slope
      DO i = 3, ncount
       newslp = (energy(i)  - energy(i-2))*0.5/dtime
      
c***** if turn is <= 0 we have a turning point 
       turn = oldslp*newslp
c***** have we found an exact turning point?
       IF (DABS(newslp) .LT. 1e-7) THEN 
C***** if lmax is true a maximum has been achieved
        IF (lmax) THEN
          imax = imax + 1
          tmax(imax) = time(i)
          emax(imax) = energy(i)
          IF (imax .NE. 1) period(imax) = tmax(imax) - tmax(imax-1)
c***** next turning point must be a minimum
          lmax = .FALSE.
          oldslp = -1
c***** we have a minimum, no special case necessary since period
c***** is calculated from the maxima
        ELSE
          imin = imin + 1
          emin(imin) = energy(i)
          tmin(imin) = time(i)
c***** next turning point must be a minimum
          lmax = .TRUE.
          oldslp = +1
        END IF
        PRINT *, 'zero slope encountered'
c***** we have ended the code needed for an exact turning point
c***** we now move onto the code for a crossing of the zero line
c***** interpolation needed here to find the turning point
        ELSE IF (turn .LT. 0) THEN
c***** its a maximum
        IF (lmax) THEN
          imax = imax + 1
          emax(imax) = (energy(i-2)*DABS(newslp) + 
     +    energy(i-1)*DABS(oldslp))
     +    /(DABS(oldslp) + DABS(newslp))
          
          tmax(imax) = (time(i-2)*DABS(newslp) + 
     +    time(i-1)*DABS(oldslp))
     +    /(DABS(oldslp) + DABS(newslp))
          IF (imax .NE. 1) period(imax) = tmax(imax) - tmax(imax -1)
          lmax = .FALSE.
        
c***** its a minimum
        ELSE
          imin = imin + 1
          emin(i) = (energy(i-2)*DABS(newslp) + 
     +    energy(i-1)*DABS(oldslp))
     +    /(DABS(oldslp) + DABS(newslp))
          
          tmin(i) = (time(i-2)*DABS(newslp) + 
     +    time(i-1)*DABS(oldslp))
     +    /(DABS(oldslp) + DABS(newslp))
        
          lmax = .TRUE.
        END IF
      END IF

c***** reset the slope
      oldslp = newslp  
c***** end of the do loop
      END DO
    
      RETURN
      END 
       
       
