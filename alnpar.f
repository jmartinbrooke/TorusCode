      PROGRAM ALNPAR
C***** This program was written by John Brooke of Manchester Computing 
C***** University of Manchester, UK. May be used for academic use provided
C***** program's authorship is acknowledged. Any use for commerical 
C***** purposes must be negotiated with the author.
C*****
C***** dufort-fraenkel alpha-squared dynamo program**************
C***** rotation law- Keplerian
C*****
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

c     header file for steering function 
     
      include 'reg_steer_f77.inc'
c     declaration for steering components
      CHARACTER*(REG_MAX_NUM_STR_PARAMS)app_name
      INTEGER num_cmds,status,icmd,finished
      INTEGER commands(REG_INITIAL_NUM_CMDS)
      CHARACTER*(40) param_label_calpha,param_label_comega
      INTEGER param_steerable
      INTEGER param_ptrs
      INTEGER param_type
      CHARACTER*(40) param_minima
      CHARACTER*(40) param_maxima
      INTEGER num_recvd_cmds
      INTEGER num_params_changed
      INTEGER recvd_cmds(REG_MAX_NUM_STR_CMDS)
      CHARACTER*(REG_MAX_STRING_LENGTH)
     &   recvd_cmd_params(REG_MAX_NUM_STR_CMDS)
      CHARACTER*(REG_MAX_STRING_LENGTH) 
     &   changed_param_labels(REG_MAX_NUM_STR_PARAMS)


c for IO type  
      INTEGER num_types   
      
      CHARACTER*(40) io_labels(REG_INITIAL_NUM_IOTYPES)
      INTEGER  iotype_handles(REG_INITIAL_NUM_IOTYPES)
      INTEGER io_dirn(REG_INITIAL_NUM_IOTYPES)
      INTEGER io_freqs(REG_INITIAL_NUM_IOTYPES)

      INTEGER iohandle
      INTEGER data_type
      INTEGER data_count

      DOUBLE PRECISION array_data(4)


      
c***** tests for poloidal and toroidal convergence and tests for 
c***** restart using logical variable LSTART
c**** 12-04-97 latest version writes north and south strobed values
c**** also calculates length of icicles
      LOGICAL LSTART
c**** NI and NJ give the dimensions of the mesh, NTIME gives size of time
c**** steps in units of a diffusion time, NEFOLD gives the number of 
c**** magnetic diffusion times for which we must run the program

      PARAMETER (NI=41,NJ=160,NTIME=6000,NEFOLD=2400,IEQ=2)
      PARAMETER (N=NI*NJ)
c**** nimid, niqtr, njegt, nj3egt, njsxn, nj7egt, nj15sx, nj5egt are 
c**** used to place the control and data points for giving Poincare
c**** maps
      PARAMETER (nimid=ni/2,niqtr=ni/4,njegt=nj/8+1,nj3egt=3*nj/8+1,
     +njsxn=nj/16+1,nj7egt=nj-njegt+2,nj15sx=nj-njsxn+2,
     +nj5egt=nj-nj3egt+2)
c**** This is old style Fortran from the days before f90 compilers were
c**** widespread
c**** PMINUS, PZERO, PPLUS give the poloidal mesh values at three time levels
c**** TMINUS, TZERO, TPLUS give the toroidal mesh values at three time levels.
c**** PCOEF1 to PCOEF5 give the coefficients of the poloidal equation
c**** TCOEF1 to TCOEF5 give the coefficients of the toroidal equation
      DOUBLE PRECISION  PMINUS(NI,NJ),PZERO(NI,NJ),PPLUS(NI,NJ),
     +          TMINUS(NI,NJ),TZERO(NI,NJ),TPLUS(NI,NJ),
     +          PCOEF1(NI,NJ),PCOEF2(NI,NJ),PCOEF3(NI,NJ),
     +          PCOEF4(NI,NJ),PCOEF5(NI,NJ),TCOEF1(NI,NJ),
     +          TCOEF2(NI,NJ),TCOEF3(NI,NJ),TCOEF4(NI,NJ),
     +          TCOEF5(NI,NJ),
     +          PTCOF5(NI,NJ),
     +          TPCOF1(NI,NJ),TPCOF2(NI,NJ),TPCOF3(NI,NJ),
     +          TPCOF4(NI,NJ),TPCOF5(NI,NJ)
c**** DIFFA, UNDER are used in the timestepping routines npstep, ntstep
c**** BOUND AND ALAS2 are used in the jepps matrix boundary calculation 
c**** for the poloidal time step npstep.
c**** APPLUS, ATPLUS, SPPLUS, STPLUS are used in routine nparity that 
c**** separates the field into symmetric and antisymmetric components.
      DOUBLE PRECISION DIFFA(NI,NJ),UNDER(NI,NJ),BOUND(NJ,2*NJ),
     + ALAS2(2*NJ),ALPHA0(NJ),ALPHA(NI,NJ),SN(NJ),CS(NJ),HCS(NI),
     + HSN(NI),APPLUS(NI,NJ),ATPLUS(NI,NJ),SPPLUS(NI,NJ),STPLUS(NI,NJ)
c**** PI is the calculated value of pi, BD gives the thickness of the torus
c**** CALPHA and COMEGA are the dynamo parameters 
c**** TOL is a tolerance parameter, DX is the radial mesh separation, DTH is
c**** the angular mesh separation (angle goes 0 to pi), DTIME is the timestep
c**** XSTART and THSTAR give starting values, ESYMM and EANTI are the 
c**** symmetric and antisymmetric energies, ETOT is total energy, RP2T is the
c**** ratio of poloidal to toroidal energy, PAR is parity, EBD is 
      DOUBLE PRECISION PI,BD,CALPHA,COMEGA,A,TOL,P,HSN0, 
     +XSTART,THSTAR,DX,DTH,DTIME,PELAST,PEOLD,ESYMM,EANTI, 
     +ETOT,RP2T,PAR,EBD,X,TH,TIME,TVAL,DTVAL
c**** common block, in f90 this should go into a module
      COMMON /PAR/PI,BD,CALPHA,COMEGA,A
 

c***** start of executable statements *********************************
c***** first we read in previous saved values of two time steps from 
c***** end of previous run for both poloidal and toroidal field, this
c***** is pminus, pzero, tminus, tzero
c***** alnparn.px contains the poincare section data points for the 
c***** northern hemisphere, alnpars.px contains these points for the
c***** southern hemisphere. pbound.dat 
c***** esymmn.dat gives symmetric field

      

c     enable steering function   
      num_cmds = 2
      app_name = "magnetic simulation"
      commands(1) = REG_STR_STOP
      commands(2) = REG_STR_PAUSE      
      
      CALL steering_enable_f(reg_true)

      CALL steering_initialize_f(app_name,num_cmds,commands, status)
      IF(status .ne. REG_SUCCESS)THEN
         STOP 'Call to steering_initialize_f failed!'
      END IF      
      


      OPEN (1,FILE='pminus.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +      position='append')
      OPEN (2,FILE='pzero.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +      position='append')
C     OPEN (3,FILE='alnpar.par')
C     OPEN (8,FILE='jepps.bin',FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN (8,FILE='jepps.txt')
      OPEN (11,FILE='tminus.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +      position='append')
      OPEN (12,FILE='tzero.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +      position='append')
      OPEN (14,FILE='alnparn.px')
      OPEN (15,FILE='alnpars.px')
      OPEN (16,FILE='pbound.dat')
      OPEN (17,FILE='esymmn.dat', status='replace')
      OPEN (18,FILE='esymms.dat', status='replace')

      close(17)
      close(18)

C*****
C* Start of executable statements ***********
C*****
      PI  = 2.0D0*DASIN(1.0D0)
c***** next line means that we read the restart files from previous run
      LSTART = .FALSE.
c***** this tells us the thickness of the torus    
      BD=1.00D0
      TOL = 1.0D-2
c***** used for calculating how far out to produce poloidal streamlines 
c***** for plotting
      NIEXTR = NI + NI/2
c***** gives starting values to avoid crashes
      PENOLD = 1.0
      TENOLD = 1.0
c***** if INLIN = 0 we have a linear run, if its 1 we have a nonlinear run
      INLIN  = 1
c***** sets the dynamo parameters
      CALPHA = 5.57
      COMEGA= -1.0d+3
      P = CALPHA*COMEGA
      WRITE(*,*) CALPHA,COMEGA,P
      A=1.0D0
      PRINT *,PI
c**** calculates radial coordinate at the boundary of the torus
      HSN0 = DSINH(BD)
c**** if these are zero we use the whole domain of the torus
      XSTART = 0.00D0
      THSTAR = 0.0D0
c**** the values of the mesh and time distances for derivative calculation
      DX = (1.0D0-XSTART)/DFLOAT((NI-1))
      DTH = 2.0D0*PI* (1.0D0-THSTAR)/DFLOAT(NJ)
      DTIME = 1.0d0/(6.0D0*NTIME*HSN0**2)
c**** how often do we want to print values, higher NPRINT the less frequently
      NPRINT = 20
c**** now we set some default values to avoid crashes on runs with LSTART
c**** set to FALSE
      KOLD = 1
      JOLD = -1
      PELAST = +1D1 
      PEOLD = 1.0
      ISLOP = 1
      ESYMM = 0
      EANTI = 0
      EPS=0.0
      ETOT = 0.0D0
      R2PT = 0.0D0
      PAR  = 0.0D0

C****** read in the jepps boundary matrix, in text format ****

      READ(8,110) BOUND
110   FORMAT (4D20.13)

C     do  i=1,nj
C       READ (8,*) (BOUND(I,J),J=1,2*NJ)
C     end do  
c****** set up the circular and hyperbolic functions for all mesh points
      EBD = DEXP(BD)
      HCS(1) = 0.0D0
      HSN(1) = 0.0D0
      DO 310 I=2,NI
        X=(I-1)*DX
        HCS(I) =0.5D0*(X/EBD+EBD/X)                                    
        HSN(I) =0.5D0*(EBD/X-X/EBD)                                     
310   CONTINUE  
 
c**** calculate trig functions for all mesh points
      CS(1) = 1.0d0
      SN(1) = 0.0D0
      DO 320 J=2,NJ
        TH = (J-1)*DTH
        CS(J) = DCOS(TH)
        SN(J) = DSIN(TH)
320   CONTINUE  
C****** read in from restart file if LSTART is true******
  
      IF (LSTART) THEN 
          READ (1) TIME0
          READ (1) PMINUS
          READ (2) PZERO
          READ (11) TMINUS
          READ (12) TZERO

c***** create values if LSTART is false 
      ELSE
        DO 13 I=1,NI
          X = DFLOAT(I-1)*DX + 0.0
          DO 14 J = 1,NJ
              TH = DFLOAT(J-1)*DTH + 0.0
              SN(J) = DSIN(TH)
              CS(J) = DCOS(TH)
              PMINUS(I,J) = X*(1.0-X)*(SN(J) + CS(J))*1.0d-3  
              TMINUS(I,J) = X*(1.0-X)*(SN(J) + CS(J))*1.0d-3 +
     +          EPS*X*X*(1.0 - X)*SN(J)*CS(J)
              PZERO(I,J) = 0.9995*PMINUS(I,J)
              PPLUS(I,J) = 1.0D-3
              TZERO(I,J) = 0.9995*TMINUS(I,J)
              TPLUS(I,J) = 1.0D-3
   14      CONTINUE
   13    CONTINUE
         TIME0=0.0d0
C***** end of IF on LSTART **********
        ENDIF
        close(1)
        close(2)
        close(11)
        close(12)
c***** first call to parity function to print the parity of the field.
        TIME = TIME0

        ETOR = TEN(TZERO,NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS)     
        PRINT *,'ETOR = ',ETOR,' TPLUS(10,10)',TZERO(10,10)
        PRINT *,'TMINUS(10,10)= ', TMINUS(10,10)
        PRINT *,'first call to parity'
        CALL PARITY(PZERO,TZERO,APPLUS,SPPLUS,ATPLUS,STPLUS,
     +  ETOT,RP2T,PAR,TIME,
     +  NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS,ESYMM,EANTI) 
      
        WRITE(*,*) 'Initial parity',PAR,'Initial rp2t',RP2T

        WRITE(*,*) 'Initial energy',ETOT,'Initial time',TIME
        PRINT *,'end of first call to parity'
  
c***** call routine to calculate coeffs of the poloidal eqn ***
    
    
       CALL NPCOFF(PCOEF1,PCOEF2,PCOEF3,PCOEF4,PCOEF5,
     +             PTCOF5,NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +             IEQ,SN,CS,HSN,HCS)
c***** call routine to calculate coeffs of the toroidal eqn ***
    
      CALL NTCOFF(TCOEF1,TCOEF2,TCOEF3,TCOEF4,TCOEF5,
     +                  TPCOF1,TPCOF2,TPCOF3,TPCOF4,TPCOF5,
     +                  UNDER,ALPHA0,
     +                  NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +                  HSN,HCS,SN,CS)














c     Register some IO channels

      io_labels(1) = "POINTS_OUTPUT"
      io_dirn(1)   = REG_IO_OUT
      io_freqs(1)  = 1
  
      num_types = 1

      CALL register_iotypes_f(num_types, io_labels, io_dirn, 
     +            io_freqs, iotype_handles, status)

      IF(status .ne. REG_SUCCESS)THEN

        CALL steering_finalize_f(status)
        STOP 'Failed to register IO types'
      END IF



      WRITE(*,*) 'Returned IOtype = ', iotype_handles(1)

c     register steering parameters
      finished = reg_false
      param_label_comega = "comega"
      param_steerable = reg_true

      param_type = REG_DBL
      param_minima = ""
      param_maxima = ""
      CALL register_param_f(param_label_comega, param_steerable, +
     &       COMEGA, param_type,param_minima,param_maxima,status)      
  
      IF(status .ne. REG_SUCCESS)THEN
       CALL steering_finalize_f(status)
       STOP 'Failed to register IO params'
      END IF            
      
      param_label_calpha = "calpha"
      param_steerable = reg_true

      param_type = REG_DBL
      param_minima = ""
      param_maxima = ""
      CALL register_param_f(param_label_calpha, param_steerable, +
     &       CALPHA, param_type,param_minima,param_maxima,status)      
  
      IF(status .ne. REG_SUCCESS)THEN
       CALL steering_finalize_f(status)
       STOP 'Failed to register IO params'
      END IF          
      











c***** start of time-loop ******
     
c***** NTIME gives number of time steps per diffusion time, NEFOLD gives
c***** the length of the calculation in e-folding times.
c      DO 20 K = 1,NTIME*NEFOLD   
        

      k=1
      do  WHILE(k <= NTIME*NEFOLD .AND. (finished .ne. 1))
      
      TIME=K*1.0d0/NTIME + TIME0

      CALL steering_control_f(K, num_params_changed,
     &             changed_param_labels, num_recvd_cmds, 
     &             recvd_cmds, recvd_cmd_params, status)

     
      IF(status .eq. REG_SUCCESS)THEN
         IF(num_params_changed > 0)THEN
           DO iparam=1, num_params_changed, 1

              WRITE(*,*) 'Changed param no. ', iparam,' = ', 
     +                     changed_param_labels(iparam)

             IF(changed_param_labels(iparam) == param_label_comega)THEN
c............... comega parameter changed, update initial value.
               P = CALPHA*COMEGA
       CALL NPCOFF(PCOEF1,PCOEF2,PCOEF3,PCOEF4,PCOEF5,
     +             PTCOF5,NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +             IEQ,SN,CS,HSN,HCS)
c***** call routine to calculate coeffs of the toroidal eqn ***
    
      CALL NTCOFF(TCOEF1,TCOEF2,TCOEF3,TCOEF4,TCOEF5,
     +                  TPCOF1,TPCOF2,TPCOF3,TPCOF4,TPCOF5,
     +                  UNDER,ALPHA0,
     +                  NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +                  HSN,HCS,SN,CS)

             END IF
             IF(changed_param_labels(iparam) == param_label_calpha)THEN
c............... calpha parameter changed, update initial value.
               P = CALPHA*COMEGA
       CALL NPCOFF(PCOEF1,PCOEF2,PCOEF3,PCOEF4,PCOEF5,
     +             PTCOF5,NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +             IEQ,SN,CS,HSN,HCS)
c***** call routine to calculate coeffs of the toroidal eqn ***
    
      CALL NTCOFF(TCOEF1,TCOEF2,TCOEF3,TCOEF4,TCOEF5,
     +                  TPCOF1,TPCOF2,TPCOF3,TPCOF4,TPCOF5,
     +                  UNDER,ALPHA0,
     +                  NI,NJ,DX,DTH,DTIME,XSTART,THSTAR,
     +                  HSN,HCS,SN,CS)

             END IF
           END DO
         END IF







         IF(num_recvd_cmds > 0)THEN           
	   icmd = 1
           DO
c	      write(*,*) 'icmd = ', icmd
	      SELECT CASE (recvd_cmds(icmd))
                 CASE(REG_STR_PAUSE)
                   CALL steering_pause_f(num_params_changed, 
     &	              changed_param_labels, num_recvd_cmds, recvd_cmds,
     &                         recvd_cmd_params, status)

                IF(status .ne. REG_SUCCESS)THEN
                  WRITE(*,*) 'steering_pause_f returned error'
                  EXIT
                ELSE
                   icmd = 0
                END IF

               CASE(REG_STR_STOP)
                  WRITE (*,*) 'Received stop command from steerer'
                  finished = reg_true 
                  EXIT

               CASE DEFAULT
c                  write(*,*) 'default icmd = ', icmd


                   IF(recvd_cmds(icmd) .EQ. iotype_handles(1))THEN

c                    WRITE(*,*) 'Emitting data...'
c                     WRITE(6,*)k,TAU,EM,QM,gmx,pm2,amom
c                     WRITE(6,*)array_data(5),array_data(6),array_data(7),array_data(4),array_data(5),array_data(6),array_data(7)

            	    IF(status .eq. REG_SUCCESS)THEN
                      CALL emit_start_f(iotype_handles(1), i, 
     +                   iohandle, status)
                    END IF


                    IF( status .eq. REG_SUCCESS )THEN


                       data_count = 4
                       data_type  = REG_DBL
                       CALL emit_data_slice_f(iohandle, data_type, 
     +                       data_count,array_data, status)

                       CALL emit_stop_f(iohandle, status)

                    END IF
c                     WRITE (*,*) '...done'

                  END IF


             END SELECT
	     icmd = icmd + 1	   
	      
	       IF ((icmd .gt. num_recvd_cmds) 
     &	         .OR. (finished .eq. reg_true)) EXIT
	   END DO
	  

	  
	 END IF
       END IF     






   
c***** call routine to time-step the poloidal equation ***     
    
      CALL NPSTEP(PCOEF1,PCOEF2,PCOEF3,PCOEF4,PCOEF5,
     +                 PTCOF5,PMINUS,PZERO,PPLUS,TZERO,
     +                 DIFFA,UNDER,
     +                 ALPHA0,ALPHA,ALAS2,BOUND,
     +                 DX,DTH,DTIME,XSTART,NI,NJ,INLIN,
     +                 HSN,HCS,SN,CS) 
     

c***** call routine to time-step the toroidal equation ****    

C     PRINT *,'(before tstep)'
      CALL NTSTEP(TCOEF1,TCOEF2,TCOEF3,TCOEF4,TCOEF5,
     +                 TPCOF1,TPCOF2,TPCOF3,TPCOF4,TPCOF5,
     +                 TMINUS,TZERO,TPLUS,PZERO,
     +                 DX,DTH,DTIME,NI,NJ,XSTART,
     +                 INLIN,DIFFA,UNDER,ALPHA0,ALPHA,
     +                 HSN,HCS,SN,CS)



c***** write values for Poincare map (north hemisphere) 
c***** this routine monitors a control point to see when the field changes
c***** sign. Therefore we must set up some old values before the first
c***** call to the monitor routine
      tval = tplus(nimid,nj3egt) 
      dtval = DABS(tval)
      jnewn = IDNINT(tval/dtval)
      itestv = jnewn*joldn
      IF (itestv .LT. 0) THEN 
C     CALL poin(pminus(niqtr,njsxn),pzero(niqtr,njsxn),
C    +pplus(niqtr,njsxn),
C    +tplus(nimid,njegt),time,dtime,k,"north",0)
c***** call routine to print the values of the Poincare crossing
      CALL poin(pplus(nimid,njegt),
     +tplus(nimid,njegt),time,k,"north",0)
c**** call routine to calculate the parity each time we cross the Poincare
c**** hyperplane
      CALL PARITY(PPLUS,TPLUS,APPLUS,SPPLUS,ATPLUS,STPLUS,
     +  ETOT,RP2T,PAR,TIME,
     +  NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS,ESYMM,EANTI) 

c***** old way is to write the parity, new way is to write symmetric and 
c***** antisymmetric field energies      
C       WRITE(3,99995) PAR,RP2T,ETOT,TIME
        OPEN (17,FILE='esymmn.dat', STATUS='UNKNOWN', position='append')
        WRITE(17,99997) ESYMM, EANTI, TIME, '0'
        CLOSE(17)
        array_data(1) = ESYMM
        array_data(2) = EANTI
        array_data(3) = TIME
        array_data(4) = 0
      ENDIF
c***** put new value of control point to old value
      joldn = jnewn
         
c***** write values for Poincare map (south hemisphere), all details are as 
c***** for the northern

      tval = tplus(nimid,nj5egt) 
      dtval = DABS(tval)
      jnews = IDNINT(tval/dtval)
      itestv = jnews*jolds
      IF (itestv .LT. 0) THEN
C     CALL poin(pminus(niqtr,nj15sx),pzero(niqtr,nj15sx),
C    +pplus(niqtr,nj15sx),
C    +tplus(nimid,nj7egt),time,dtime,k,"south",1)
      CALL poin(pplus(nimid,nj7egt),
     +tplus(nimid,nj7egt),time,k,"south",1)
        CALL PARITY(PPLUS,TPLUS,APPLUS,SPPLUS,ATPLUS,STPLUS,
     +  ETOT,RP2T,PAR,TIME,
     +  NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS,ESYMM,EANTI) 
      
C       WRITE(3,99995) PAR,RP2T,ETOT,TIME
        OPEN (18,FILE='esymms.dat', STATUS='UNKNOWN', position='append')
        WRITE(18,99997) ESYMM, EANTI, TIME, '1'
        CLOSE(18)
        array_data(1) = ESYMM
        array_data(2) = EANTI
        array_data(3) = TIME
        array_data(4) = 1
      ENDIF
      jolds = jnews
         
c***** print observables every 1/1200 difftime ********

C     KPRINT = NTIME/200 
  
C     IF (MOD(K,KPRINT) .EQ. 0)  THEN
  
c***** call to parity routine, returns symmetric and antisymmetric    
c***** fields
C       CALL PARITY(PPLUS,TPLUS,APPLUS,SPPLUS,ATPLUS,STPLUS,
C    +  ETOT,RP2T,PAR,TIME,
C    +  NI,NJ,BD,DX,DTH,PI,HSN,HCS,SN,CS,ESYMM,EANTI) 
      
C       WRITE(3,99995) PAR,RP2T,ETOT,TIME
C       WRITE(17,99997) ESYMM, EANTI, TIME
C       WRITE(*,*) PAR,RP2T,ETOT,TIME
          
   
C     ENDIF 
  
c***** swap time levels *********************** *********
 
          DO 30 I = 1,NI
              DO 31 J = 1,NJ
                  PMINUS(I,J) = PZERO(I,J)
                  PZERO(I,J) = PPLUS(I,J)
                  PPLUS(I,J) = PMINUS(I,J)
                  TMINUS(I,J) = TZERO(I,J)
                  TZERO(I,J) = TPLUS(I,J)
                  TPLUS(I,J) = TMINUS(I,J)
   31         CONTINUE
   30     CONTINUE
    
c***** end of time loop, program stops if convergence not reached**
   
c   20 CONTINUE
      k = k+1
      end do       
c***** no convergence, therefore *****
c***** write results to channels 1 and 2 (poloidal) ***
c***** write results to channes 11 and 12 (toroidal) ***
       PRINT *,'end of time-loop'
       OPEN (1,FILE='pminus.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +       position='append')
       OPEN (2,FILE='pzero.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +       position='append')
       OPEN (11,FILE='tminus.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +       position='append')
       OPEN (12,FILE='tzero.bin',FORM='UNFORMATTED',STATUS='UNKNOWN', 
     +       position='append')
       WRITE(1) TIME
       WRITE(1) PMINUS
       WRITE(2) PZERO
       WRITE(11) TMINUS
       WRITE(12) TZERO
       WRITE (16,*) (PZERO(NI-2,J),J=1,NJ)
       WRITE (16,*) (PZERO(NI-1,J),J=1,NJ)
       PRINT *,'end of write-statements'
   
999   CLOSE(1)
      CLOSE(2)
      CLOSE(3)
C     CLOSE(4)
      CLOSE(8)
      CLOSE(10)
      CLOSE(11)
      CLOSE(12)
      CLOSE(13)
      CLOSE(14)
      CLOSE(15)
      CLOSE(16)
       PRINT *,'end of program'
1000  STOP
 
150   FORMAT('CALPHA = ',1P,D10.3,'   COMEGA = ',1P,D10.3,'    P = ',
     *       1P,D10.3)
151   FORMAT('BD = ',D10.3,/,'NUMBER OF ITERATIONS
     * = ',I10,/,'DTIME = ',1P,D10.3,/,'NI = ',I5,'    NJ = ',I5,
     * /,
     * 'IEQ = ',I5)
152   FORMAT(1X,'EQUATION DID NOT CONVERGE')
99999 FORMAT (1P,8D10.3)
99998 FORMAT (/,'RATE AT CENTRE =',1P,D10.3,/,'RATE AT SURFACE =',1P,
     +       D10.3,/,'RATE AT RSIDE =',1P,D10.3,/,'RATE AT LSIDE =',
     +       1P,D10.3)
99997 FORMAT (3D20.13,1X,A)
99996 FORMAT (4D10.3,D20.13)
99995 FORMAT (4D20.13)
99994 FORMAT (2D10.3)
90000 FORMAT (1P,4D20.13)
89999 FORMAT  (1P,d20.13)
      END
 
c***** this subroutine writes the values of the crossings of the 
c***** Poincare section
      SUBROUTINE poin(pplus,tplus,time,k,hemisp,funit)
      DOUBLE PRECISION pplus,tplus,time   
      INTEGER k
      CHARACTER*6 hemisp
      INTEGER funit

      WRITE (*,89998) hemisp,pplus,tplus,k
      IF (funit .EQ. 0) THEN
        WRITE (UNIT=14,FMT='(3e21.13)') pplus,tplus,time
      ELSE 
        WRITE (UNIT=15,FMT='(3e21.13)') pplus,tplus,time
      ENDIF
        
89998 FORMAT(A10,2e10.3,I8)
         
      RETURN
      END
 
 
C     SUBROUTINE poin(pmin,pzer,pplus,tplus,time,dtime,k,hemisp,funit)
C     DOUBLE PRECISION pmin,pzer,pplus,tplus,time   
C     INTEGER k
C     CHARACTER*6 hemisp
C     INTEGER funit
C     DOUBLE PRECISION dadt 

C     dadt = (1.5d0*pplus - 2.0d0*pzer
C    +      + 0.5d0*pmin)/dtime
C     WRITE (*,89998) hemisp,dadt,tplus,k
C     IF (funit .EQ. 0) THEN
C       WRITE (UNIT=14,FMT='(3e20.13)') dadt,tplus,time
C     ELSE 
C       WRITE (UNIT=15,FMT='(3e20.13)') dadt,tplus,time
C     ENDIF
C       
C     WRITE (UNIT=funit,FMT='(3e20.13)') dadt,tplus,time
C89998 FORMAT(A10,2e10.3,I8)
C        
C     RETURN
C     END
