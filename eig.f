      subroutine eig(enew,elast,eold,knew,kold,islop,
     *                dtime,period,grow,icon)
      IMPLICIT REAL*8 (A-H,O-Z)
      real*8 enew,eold,dtime,period,size,enewl,eoldl,elast
      integer knew,kold,islop,islopn,ichange
      intrinsic dabs,dnint
      size = dabs(enew - elast)
      islopn = (enew - elast)/size
      ichange = islop*islopn
      if (ichange .lt. 1) then
        period = (knew - kold)*dtime
        enewl = dlog(dabs(enew))
        eoldl = dlog(dabs(eold))
        grow = (enewl - eoldl)/period
        eold = enew
        kold = knew
        icon = 1
        write (20,*) 'knew = ',KNEW,' period = ',PERIOD
      endif
      islop = islopn
      elast = enew
      return
      end

