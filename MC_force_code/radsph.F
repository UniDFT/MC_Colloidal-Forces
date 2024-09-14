      program  main
      implicit none
#include<mpif.h>
       integer myid, nprocs, ierror, iseeds(10), iseed
       data iseeds/11225641,8762346761,765123981,9787125,8712678901
     &   ,9873456723,97234627,2345691,9876123,651236919/
c initilize parallel mode
      call mpi_init(ierror)
      call mpi_comm_size(mpi_comm_world, nprocs, ierror)
      call mpi_comm_rank(mpi_comm_world, myid, ierror)
c initilize random number generator for multiprocessors
       iseed = iseeds(myid+1)
       iseed = ranset(iseed)
c  generate a new seed each time program is run and initialize iseed
       call pforce(myid)
c close parallel mode
c       call mpi_finalize(ierror)
       end


!	***********************************************************************
	subroutine pforce(myid)
!          NVT_MC FOR A DISCRETELY CHARGED MACROION IN 
!                     AN ELECTROLYTE SOLUTION
!	***********************************************************************
      IMPLICIT NONE
      INTEGER equil, prod, ii, icycl,npmax, ndispl
      INTEGER  attempt,nacc, ncycl, nmoves, imove, IL,ibeg 
      integer NPART,NCHARGE1,nstart,myid
      parameter (npmax=5000)
      real x(0:npmax),y(0:npmax),z(0:npmax)
      real sig(0:npmax),q(0:npmax)
      REAL     Ener, Dr, GW,angle
      integer atmp,naccpv
      common/conf1/x,y,z,npart,NCHARGE1
	  common/pot1/sig,q
!   ---initialize sysem
      CALL READDAT(myid,ibeg,equil, prod, ndispl, dr, angle,GW)
	  IF (ibeg.LT.2) THEN
		nstart=1
		ELSE
		nstart=2
	  ENDIF
      nmoves = ndispl
c_---total energy of the system
	  
	 
	  CALL TOTERG(GW, Ener)
      WRITE (*, 99001) Ener
!   ---start MC-cycle
      DO 100 ii = nstart, 2
!   --- ii=1 equilibration
!   --- ii=2 production
         IF (ii.EQ.1) THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (*, *) ' Start equilibration '
         ELSE
            ncycl = prod
            IF (ncycl.NE.0) WRITE (*, *) ' Start production '
         ENDIF
!   --- initilize sampling variables 
		CALL SAMPLE(myid,ibeg,IL, 0, Ener,Dr)
	 attempt = 0
         nacc = 0
	 atmp=0
	 naccpv=0
!   ---intialize the subroutine that adjust the maximum displacement
         CALL ADJUST(attempt, nacc, dr)
	 CALL ADJUSTANGLE(ATMP, NACCPV, ANGLE)
         DO 50 icycl = 1, ncycl
	   DO imove = 1, nmoves
!   ---attempt to displace a particle
         CALL MCMOVE(myid,IL, Ener, attempt, nacc, dr,atmp,naccpv,angle)
	 ENDDO
	       IF (ii.EQ.2) THEN
			 CALL SAMPLE(myid,ibeg,IL, 1, Ener,Dr)
               END IF

            IF (MOD(icycl,ncycl/5).EQ.0) THEN
         WRITE(*, *) '======>> Done ', icycl, ' out of ', ncycl
!   ---adjust maximum displacements	and maximum angular displacement
           CALL ADJUST(attempt, nacc, Dr)
           call adjustangle(atmp,naccpv,angle)						
           END IF
50       END DO

         IF (ncycl.NE.0) THEN
         IF (attempt.NE.0) then 
	write(6,99003) attempt,nacc,100.*FLOAT(nacc)/FLOAT(attempt)
	 endif
	 IF (ATMP.NE.0) then 
        write(6,99004) atmp, naccpv,100*FLOAT(NACCPV)/FLOAT(ATMP)
	 endif
	ENDIF
100     END DO

	CALL SAMPLE(myid,ibeg,IL, 2, Ener,Dr)
	WRITE(*,*) 'FINAL AVERAGE ENERGY/NkT = ', Ener      
 
99001  FORMAT (' Total energy at initial configuration: ', f12.5, /)
99002  FORMAT (' Total energy at the end of simulation: ', f12.5, /,
     &	    '		 running energy              : ', f12.5, /,    
     &       '		 difference                  :  ', e12.5, /) 
99003  FORMAT (' Number of att. to displ. a part.  : ', i10, /,  
     &	      ' success: ', i10, '(= ', f5.2, '%)')
99004  FORMAT (' Number of att. to displ. a charge.  : ', i10, /, 
     &	      ' success: ', i10, '(= ', f5.2, '%)')

       END
	
!	****************************************************************
!	THIS IS THE FIRST SUBROUTINE CALLED FROM THE PFORCE UNIT.
!	****************************************************************	
	SUBROUTINE READDAT(myid,ibeg,Equil,Prod,Ndispl,Dr,angle,GW)
!     ---input solution condition and model parameters

!     ---input parameters: file: (3), input.dat
!    ibeg  =  0 : initialize from a lattice
!             1 : read configuration from disk
!             2 : restart sampling    
!    Equil      : number of Monte Carlo cycles during equilibration
!    Prod       : number of Monte Carlo cycles during production
!    Dr         : maximum displacement
!    Ndispl     : number of attemps to displace a particle per MC cycle
!    NPART      : total number of particles
!    NPARTI(I)  : number of particles for species i
!    NION   = number of ion species
!    NIONmax= maxmum number of ion species
!    sig(i) = diameters of ions and macroions/Bjerrum length
!    sigi(i)= diameters of ions and macroions species, in Angstroms
!    q(i)   = valence for an ion
!    qi(i)  = valence for an ion, for input
!    CM(i)  = molar concentration of ions
!    lB0    = Bjerrum length for water at 298.15 K, in Angstroms
!    lB     = Bjerrum length at system temperature, in Angstroms

!    BOX    = simulation box length in Angstroms (FOR INPUT)
!     ---input parameters: file: iniconf.dat (restart file
!                to continue a simulation from disk)
!    npmax  = maximum number of particles
!    Dr     = optimized maximum displacement old configurations
!    GW     = Gaussian Width parameter, adjust for efficiency of Ewald Sum
!    X(1),Y(1),Z(1)            : position first particle 1
!        ...
!    X(NPART),Y(NPART),Z(NPART): position particle last particle
 
      IMPLICIT NONE
      REAL    pi,sig12, box,  hbox,con1, con2
      real Dr,angle, lB0, lB, GW, SHBOX
      INTEGER ibeg,Equil,Prod,i,j,k,Ndispl,NION,ngrid
      integer NIONmax,npmax,npart,NHIsmax
      integer NCHARGE1
      PARAMETER (npmax=5000, NIONmax=10, lB0=7.144464E0,ngrid=450)
      parameter(con1=2.338874e4, con2=6.022137e-4,NHIsmax=200)
      REAL sig(0:npmax), q(0:npmax), ZERO
      REAL IS, sigi(NIONmax), qi(NIONmax), CM(NIONmax)
      REAL     x(0:npmax), y(0:npmax), z(0:npmax)
      real Enert0
      integer ns0,nnode
      INTEGER  NPARTI(NIONmax),myid
	real ngcounter0(1:ngrid),ngco0(1:ngrid)
      COMMON /sys1/ box,hbox
      COMMON /pot1/sig,q
      COMMON /conf1/x,y,z,NPART,NCHARGE1
      COMMON /samm/NION,NPARTI
      COMMON /sam/SHBOX,lB
      COMMON /restart/ns0,Enert0,ngcounter0,ngco0
      character filename(20)*5,filename2(20)*5
      character OUTPUT1*5,INICONF*5 ,INPUT*5
      data filename/'OUT1','OUT2','OUT3','OUT4','OUT5','OUT6','OUT7',
     & 		'OUT8','OUT9','OUT10','OUT11','OUT12','OUT13','OUT14',
     &	 	'OUT15','OUT16','OUT17','OUT18','OUT19','OUT20'/
      data filename2/'INP1','INP2','INP3','INP4','INP5','INP6','INP7',
     &	 	'INP8','INP9','INP10','INP11','INP12','INP13','INP14',
     &		'INP15','INP16','INP17','INP18','INP19','INP20'/

	  INICONF=filename2(myid+1)
	  OUTPUT1=filename(myid+1)
          OPEN (65, FILE='input')
!     ---read simulation parameters
      READ (65, *)
      READ (65, *) nnode
      READ (65, *)
      READ (65, *) ibeg, Equil, Prod
      READ (65, *)
      READ (65, *) Dr,angle, GW
      READ (65, *)
      READ (65, *) Ndispl
      READ (65, *)
      READ (65, *) BOX
      READ (65, *)
!     ---read solution condition and model parameters
      READ (65, *) NION
      READ(65, *)
      DO i=1,NION
      READ (65, *) sigi(i), qi(i), NPARTI(i)
      ENDDO
      CLOSE(65)
! check charge neutrality at initial condition
	ZERO=0.0E0
        do i=1,NION
 	ZERO=ZERO+qi(i)*NPARTI(I)
	enddo
	IF (ZERO.GT.1.e-8) THEN
	WRITE(*,*)'Charge neutrality is unsatisfied at initial condition!' 
	PAUSE
	STOP
	ENDIF

	NPART=0
	DO i=1,NION-1
	NPART=NPART+NPARTI(i)
	ENDDO
	NCHARGE1=1
	  IF (NPART.GT.NPMax) THEN
         WRITE (*, *) ' ERROR: number of particles too large'
         STOP
        ENDIF

!     ---calculate parameters:
      PI = 3.14159265
	lB=lB0

!     --- normalized with Bjerrum length
	BOX = BOX/lB
        HBOX = BOX/2.E0
	SHBOX = HBOX*HBOX

	k=0
	DO 10 i = 1,NION-1
	DO j = 1, NPARTI(i)
	k=k+1
	sig(k) = sigi(i)/lB
	q(k) = qi(i)
	ENDDO
10	ENDDO

	DO J=1,NCHARGE1
	K=K+1
	sig(k)=0.4/lB
	Q(K)=-20
	ENDDO
	
	
	SIG(0)=SIGI(3)/lB
!	here define the charge of the macroion within the program. 	
	q(0)=0

!   --- assign LOCATION of the macroions
							 
	X(0)=0.0 
	Y(0)=0.0 
	Z(0)=0.0
	
  	
!     ---read/generate initial configuration
      IF (ibeg.EQ.0) THEN
         WRITE (*, *) 'generate randomly initial configuration'
         CALL CONFIG
         OPEN(64,FILE=INICONF)
         WRITE (64, *) NPART,Dr
         ns0=0
         Enert0=0.0
	do i=1,ngrid
	ngcounter0(i)=0
	ngco0(i)=0
	enddo

         WRITE (64, *) ns0, Enert0
	do i=1,ngrid
	write(64, *)ngcounter0(i),ngco0(i)
	enddo


         DO i = 1, NPART+NCHARGE1
         WRITE (64, * ) sig(i), q(i), X(i), Y(i), Z(i)
         END DO
         CLOSE(64)

      ELSEIF (ibeg.EQ.1) THEN
         WRITE (*, *) ' read conf from disk '
	 OPEN(64,FILE=INICONF)
	 READ (64, *) NPART, Dr
         READ(64, *) ns0, Enert0
	do i=1,ngrid
	READ(64, *)ngcounter0(i),ngco0(i)
	enddo
        do i=1,ngrid
        ngcounter0(i)=0
        ngco0(i)=0
        enddo

         ns0=0
         Enert0=0.0
         DO i = 1, NPART+NCHARGE1
         READ (64, *) sig(i), q(i), X(i), Y(i), Z(i)
	 ENDDO
	 CLOSE(64)

      ELSEIF (ibeg.EQ.2) THEN
	 WRITE(*,*) 'restart sampling'
         OPEN(64,FILE=INICONF)
         READ (64, *) NPART, Dr
         READ(64,*) ns0, Enert0
	do i=1,ngrid
	READ(64, *)ngcounter0(i),ngco0(i)
	enddo

         DO i = 1, NPART+NCHARGE1
         READ (64, *) sig(i), q(i), X(i), Y(i), Z(i)
         ENDDO
         CLOSE(64)

      ENDIF

      END


!	************************************************************************
      SUBROUTINE CONFIG
!	***********************************************************************
!     place `npart' particles on a simple cubic lattice
!
      IMPLICIT NONE
      INTEGER I, K, npart, npmax,N,FACTOR2,NCHARGE1,s
      PARAMETER (npmax=5000)
      REAL x(0:npmax), y(0:npmax), z(0:npmax),FACTOR1
      REAL box,hbox
	real phi,theta
      real costheta, cosphi, sintheta, sinphi
      REAL sig(0:npmax), q(0:npmax), RADIUS
      REAL SQXYZ,XD, YD, ZD, SRIK,RAN
      real*8 ranf
      
      COMMON /conf1/ x, y, z, npart,NCHARGE1
      COMMON /sys1/ box,hbox
      COMMON /pot1/ sig,q

	  DO 2 I=1,NPART
1       RAN=ranf()
	  X(I)=HBOX*(RAN-0.5)*2.
	  RAN=ranf()
	  Y(I)=HBOX*(RAN-0.5)*2.
        RAN=ranf()
	  Z(I)=HBOX*(RAN-0.5)*2.

      DO 3 K=0,i-1
      XD=X(I)-X(K)
      IF(XD.GT.HBOX) XD=XD-BOX
      IF(XD.LT.-HBOX) XD=XD+BOX
      YD=Y(I)-Y(K)
      IF(YD.GT.HBOX) YD=YD-BOX
      IF(YD.LT.-HBOX) YD=YD+BOX
      ZD=Z(I)-Z(K)
      IF(ZD.GT.HBOX) ZD=ZD-BOX
      IF(ZD.LT.-HBOX) ZD=ZD+BOX
      SRIK=(XD*XD+YD*YD+ZD*ZD)
      IF ((SRIK).LT.((SIG(I)+SIG(K))/2.)**2) GOTO 1

3     CONTINUE
2     CONTINUE
	
!	RANDOMLY ARRAGE NCHARGE1 CHARGES AROUND A SPHERES		
		
	DO 112 N=NPART+1, NPART+NCHARGE1
	radius= 0.0 		! (sig(0)-sig(n))/2	
	phi=2*(ranf())*3.1415926535
	sinphi=sin(phi)
	cosphi=cos(phi)
114	theta=2*(ranf())*3.1415926535
	costheta=cos(theta)
	sintheta=sin(theta)
	X(N)=radius*costheta*sinphi+X(0) 
	y(N)=radius*sinphi*sintheta+Y(0)
	z(N)=radius*cosphi+Z(0)
	do 113 i=npart+1,n-1
c       Checking overlap of point charges
        XD=X(I)-X(N)
        YD=Y(I)-Y(N)
        ZD=Z(I)-Z(N)
        sqxyz=(xd**2+yd**2+zd**2)
        if((sqxyz).lt.(sig(n)**2)) then
        goto 114
        endif
113     enddo
	
112	enddo

  	RETURN
	END
!	*****************************************************************
        SUBROUTINE SAMPLE(myid,ibeg,IL,Switch,Ener,Dr)

!    Switch (input) = 0: initialize varibales 
!                   = 1: sample averages
!                   = 2: return average results
!
!  Ener (input) : total energy/kT
!  Enert        : average total energy/NkT
!  ns           : number of sample points
!  F            : force between macroions
c  SN           : Shell Number
      INTEGER IL, npmax, ns, Switch, NIONmax, NHIsmax,ngrid
      PARAMETER (npmax=5000, NIONmax=10, NHIsmax=200, ngrid=450)
      INTEGER  NION, NPARTI(NIONmax),ibeg,sn,myid,npart
      INTEGER  I,nhcc,n(ngrid)
      real ngcounter(1:ngrid),ngco(1:ngrid)
	real ngcounter0(1:ngrid),ngco0(1:ngrid)
      real ngcounter1(1:ngrid),ngco1(1:ngrid)
	real rhoco(1:ngrid),rhocounter(1:ngrid)
      REAL     Ener, Enert,delr,delri
      REAL     Enert1,deltar,add
      REAL     DNEW(0:npmax),DOLD(0:npmax),XIL,YIL,ZIL
      REAL     sig(0:npmax),q(0:npmax),lB
      REAL     box, hbox,SHBOX
      REAL     x(0:npmax), y(0:npmax), z(0:npmax),sqrt3
      REAL     Dr,rxij,ryij,rzij,rij,rijsq
	real Enert0
      integer ns0
      INTEGER ncx(4),ncoll,NCHARGE1

      PARAMETER (sqrt3=1.7320508076)
	character filename(20)*5,filename2(20)*5
	character OUTPUT1*5,INICONF*5

      data filename/'OUT1','OUT2','OUT3','OUT4','OUT5','OUT6',
     &		'OUT7','OUT8','OUT9','OUT10','OUT11','OUT12',
     &		'OUT13','OUT14','OUT15','OUT16','OUT17','OUT18',
     &		'OUT19','OUT20'/
      data filename2/'INP1','INP2','INP3','INP4','INP5','INP6',
     &		'INP7','INP8','INP9','INP10','INP11','INP12',
     &		'INP13','INP14','INP15','INP16','INP17','INP18',
     &		'INP19','INP20'/

      COMMON /sys1/  box,hbox
      COMMON /pot1/  sig,q
      COMMON /DISTAN/ DOLD, DNEW,XIL,YIL,ZIL
      COMMON /samm/ NION, NPARTI
      COMMON /sam/ SHBOX, lB
      COMMON /conf1/ x, y, z, npart,NCHARGE1
      COMMON /restart/ ns0,Enert0,ngcounter0,ngco0
      SAVE ns, Enert,ngcounter,ngco

	OUTPUT1=filename(myid+1)
	INICONF=filename2(myid+1)
           delr =hbox/250
           delri = 1.0/delr
c 	write(*, *)'sig(0)=',sig(0)
c	write(*, *)'sig(npart_1)=',sig(npart-1)
c	write(*, *)'sig(npart+1)=',sig(npart+1)
c	write(*, *)'hbox=',hbox
c	write(*, *)'delr=',delr
c       write(*, *) '1st ns=',INT((sig(0)+sig(npart-1))*0.5/delr) + 1
c	write(*, *)'Avg.Den=',70/(box*box*box)
c	write(*, *)'hbox/delr=',hbox/delr
c	write(*, *)'sqrt3*hbox/delr',sqrt3*hbox/delr
c	write(*, *) 'hbox/sig(si)=',hbox/sig(npart-1)
c	write(*, *)'sig(0)*0.5/delr=',sig(0)*0.5/delr
c	write(*, *)'sig(npart-1)/delr=',sig(npart-1)/delr
c	write(*, *)'(hbox-sig(0)*0.5)/sig(npart-1)='
c     &			,(hbox-sig(0)*0.5)/sig(npart-1)
c	write(*, *)'sig(0)/sig(npart-1)=',sig(0)/sig(npart-1)
c	write(*, *)'sig(npart-1)*10+sig(0)=',10*sig(npart-1)+sig(0)
c	read(*, *)

!        ---Initialize
      IF (Switch.EQ.0) THEN
		ns=ns0
		Enert=Enert0*ns
	do i=1,ngrid
	ngcounter(i)=ngcounter0(i)*ns
	ngco(i)=ngco0(i)*ns
	enddo

	ELSE IF (Switch.EQ.1) THEN
		ns=ns+1
	do i=1,nparti(1)
	rxij=x(i)
	ryij=y(i)
	rzij=z(i)
	rijsq=(rxij**2+ryij**2+rzij**2)
	rij=sqrt(rijsq)
	sn=INT((rij/delr))+1
	ngcounter(sn)=ngcounter(sn)+1
	enddo

	do i=nparti(1)+1,npart
	rxij=x(i)
	ryij=y(i)
	rzij=z(i)
	rijsq=(rxij**2+ryij**2+rzij**2)
	rij=sqrt(rijsq)
	sn=INT((rij/delr))+1
	ngco(sn)=ngco(sn)+1
	enddo

		Enert = Enert + Ener

         IF (MOD(ns,10**2).EQ.0) THEN

	Enert1=Enert/Float(ns)

	do i=1,ngrid
	ngcounter1(i)=ngcounter(i)/float(ns)
	ngco1(i)=ngco(i)/float(ns)
	enddo

         OPEN(67,FILE=OUTPUT1)
         WRITE(67,*)ns,Enert1
	do i=1,ngrid
	write(67, *) ngcounter1(i),ngco1(i)
	enddo
        CLOSE(67)

      OPEN(69,file=INICONF)
   	  WRITE (69, *) NPART,Dr
        WRITE(69,*) ns, Enert1
	do i=1,ngrid
	write(69, *) ngcounter1(i),ngco1(i)
	enddo

      	DO i = 1, NPART+ncharge1
      	WRITE (69, 10) sig(i), q(i), X(i), Y(i), Z(i)
      	END DO
      	CLOSE(69)
10    	FORMAT(2X,6(F11.6,1X))
	ENDIF

	ELSE
	Enert1=Enert/Float(ns)

	do i=1,ngrid
	ngcounter1(i)=ngcounter(i)/float(ns)
	ngco1(i)=ngco(i)/float(ns)
	enddo

      OPEN(69,file=INICONF)
      WRITE (69, *) NPART, Dr
      WRITE(69,*) ns, Enert1
	do i=1,ngrid
	write(69, *) ngcounter1(i),ngco1(i)
	enddo

      DO i = 1, NPART+ncharge1
      WRITE (69, *) sig(i), q(i), X(i), Y(i), Z(i)
      END DO
      CLOSE(69)

      OPEN(77,FILE='results')
	write(77,*) Enert1
	do i=1,250
	rhocounter(i)=ngcounter1(i)/((4)*3.141592*delr*delr*delr*i**2)
	rhoco(i)=ngco1(i)/((4)*3.141592*delr*delr*delr*i**2)
	write(77, *)ns,i,delr*i/sig(npart-1),rhocounter(i),rhoco(i)
	enddo
      CLOSE(77)
30      FORMAT(1x,F4.2,1x,4E14.6)

      END IF
      RETURN
      END

c	**********************************************************************
      SUBROUTINE ADJUST(Attemp, Nacc, Dr)
!
!     adjusts maximum displacement such that 50% of the
!     moves will be accepted
!
!  Attemp (input) number of attemps that have been performed to displace a particle
!  Nacc   (input) number of successful attemps to displace a particle
!  Dr     (output) new maximum displacement
!

      IMPLICIT NONE
      INTEGER Attemp, Nacc, attempp, naccp
      REAL     box,hbox
      REAL     dro, frac, Dr
      SAVE naccp, attempp
      COMMON /sys1/ box,hbox

	  IF (Attemp.EQ.0.OR.attempp.GE.Attemp) THEN
         naccp = Nacc
         attempp = Attemp
      ELSE
         frac = FLOAT(Nacc-naccp)/FLOAT(Attemp-attempp)
         dro = Dr
         Dr = Dr*ABS(frac/0.5E0)
!        ---limit the change:
         IF (Dr/dro.GT.1.5E0) Dr = dro*1.5E0
         IF (Dr/dro.LT.0.5E0) Dr = dro*0.5E0
         IF (Dr.GT.HBOX/2.E0) Dr = HBOX/2.E0
         WRITE (*, 99001) Dr, dro, frac, Attemp - attempp, Nacc - naccp
!        ---store nacc and attemp for next use
         naccp = Nacc
         attempp = Attemp
	  END IF
99001 FORMAT (' Max. displ. set to : ', f6.3, ' (old : ', f6.3, ')', /,
     &		 ' Frac. acc.: ', f4.2, ' attempts: ', i7, ' succes: ', i7)
      END

!	********************************************************************
      subroutine ADJUSTANGLE(atmp,naccpv,angle)
!	********************************************************************
!     this subroutine adjusts the maximum rotation of a macroion so that
!     the acceptance ratio is around 50%

	integer atmp, naccpv, attemppiv, naccpiv
	real angleo, fracangle, angle

	save naccpiv, attemppiv

!     the variables are initialized at the beginning of the program (left case) 
!     or when it is restarted (right case)

	if ((atmp.eq.0).or.(attemppiv.ge.atmp)) then
		naccpiv=naccpv
		attemppiv=atmp
	else
		fracangle=float(naccpv-naccpiv)/ float(atmp-attemppiv)
		angleo=angle
		angle=angle*abs(fracangle/0.50)
		if (angle/angleo.gt.1.50) angle=angleo*1.5
		if (angle/angleo.lt.0.50) angle=angleo*0.5 
		write (*,2) angle,angleo,fracangle,atmp-attemppiv,naccpv-naccpiv
		naccpiv=naccpv
		attemppiv=atmp
	end if
2     format('max. angle displ. set to:',f6.3,'(old:',f6.3,')',/,
     &		' frac. acc.:',f4.2,' attempts:',i7,' success:',i7)
	end

!	***************************************************************
       SUBROUTINE MCMOVE(myid,IL,Ener,attempt,nacc,dr,atmp,naccpv,angle)
!    *******************************************************************
!     attempts to displace a randomly selected particle
!  Ener   (input/output) : total energy
!  Attemp (input/output) number of attemps that have been
!                  performed to displace a particle
!  Nacc   (input/output) number of successful attemps
!                  to displace a particle
!  Dr     (input) maximum displacement
!  atmp   (input/output) number of attempts that have been
!					to move a charge on a sphere
!	naccpv   (input/output)	number of successful attempts
!					to move a charge on a sphere
!	angle   (input) maximum angular displacement

!    *******************************************************************
      INTEGER IL, Attempt, Nacc,  npmax, npart
	integer MAXK,MYID,I,NCHARGE1,CHECK
	integer atmp,naccpv 
	PARAMETER (npmax=5000, MAXK = 1000)
      REAL     box,hbox, DNEW(0:npmax), DOLD(0:npmax)
      REAL     x(0:npmax), y(0:npmax), z(0:npmax),XIL,YIL,ZIL
      REAL     Ener, xn, yn, zn, DER, DEK, Dr, ANGLE
	  REAL     VD, VDN,GW
	  real*8 ranf
      REAL     sig(0:npmax), q(0:npmax)
      COMPLEX   SUMI(MAXK),DSUMI(MAXK)
      COMMON /pot1/sig,q
      COMMON / sys1 / box, hbox
      COMMON / conf1 /x,y,z,npart,NCHARGE1
	  COMMON / MOVE /SUMI,VD,DSUMI,VDN
	  COMMON /DISTAN/DOLD,DNEW,XIL,YIL,ZIL
!    *******************************************************************

!     ---select a particle at random
      IL = INT(float(NPART)*ranf())+1
	IF(IL.LT.(NPART+1)) THEN 	
! 	---give the particle a random displacement
	Attempt = Attempt + 1
		xn = X(IL) + (ranf()-0.5E0)*Dr
		yn = Y(IL) + (ranf()-0.5E0)*Dr
		zn = Z(IL) + (ranf()-0.5E0)*Dr
	ELSE
		CALL MOVECHARGE(IL,ATMP,XN,YN,ZN,angle)
	  ENDIF
	GW=1.0
!     ---calculate energy change and assing the check variable as true, i.e as 0.
	CHECK=0
	CALL DELTAE( IL, XN, YN, ZN, DER, DEK,CHECK)
c	write(*, *) 'DER+DEK=',der+dek
c	read(*, *)

c	If the DELTAE subroutine returns a CHECK value of 1 then the move is 
c	considered invalid and we get out of the MCMOVE subroutine.

	IF(CHECK.EQ.1) then
	RETURN
	endif

!     ---acceptance test
      IF (ranf().LT.EXP(-(DER+DEK))) THEN
!        --accepted
c	write(*, *)'ACCEPTED move =',il
	if(il.lt.(npart+1)) then
	Nacc = Nacc + 1
	else
	naccpv=naccpv + 1
	endif

	DO I=0,NPART+ncharge1
	DOLD(I)=DNEW(I)
	ENDDO
		
	VD=VDN
	DO I=1,MAXK
	SUMI(I) = DSUMI(I)
	ENDDO
	Ener = Ener + DER + DEK
!       ---put particle in simulation box

	IF(IL.LT.(NPART+1)) THEN
	IF (xn.LT.-hbox) xn = xn + BOX
	IF (xn.GT.hbox) xn = xn - BOX
	IF (yn.LT.-hbox) yn = yn + BOX
	IF (yn.GT.hbox) yn = yn - BOX
	IF (zn.LT.-hbox) zn = zn + BOX
	IF (zn.GT.hbox) zn = zn - BOX
	
	x(il)=xn
	y(il)=yn
	z(il)=zn

	ELSE
		X(IL) = xn
		Y(IL) = yn
		Z(IL) = zn
	ENDIF
c        write(*, *)'Ener-after mcmove=',Ener
c        call toterg(gw,engr)
c        read(*, *)
	ELSE
!       --- rejected
c	write(*, *)'REJECTED move =',il 
	ENDIF
 5    FORMAT(2X,6F12.6)
      RETURN
      END
!	*******************************************************************
	SUBROUTINE MOVECHARGE(N,ATMP,XN,YN,ZN,anglE)
	
!	DISPLACES ONE CHARGE ON THE SURFACE OF THE MACROION.

	PARAMETER(NPMAX=5000)
	INTEGER N,ATMP
	REAL SIG(0:npmax),X(0:npmax)
	real Y(0:npmax),Z(0:npmax),q(0:npmax),XN,YN,ZN
	real*8 RANF
        real alpha, bet, gamma, cosa, cosb, cosg, sena, senb, seng
        real delx, dely, delz,cas,dx,dy,dz
	COMMON/CONF1/X,Y,Z,NPART,NCHARGE1
        COMMON /pot1/sig,q
	
	atmp=atmp+1

        ang=angle*2*3.141592653589793238

c       Euler angles, alpha, bet, gamma (randomly distributed around 0, 2pi)
        cas=ranf()
        alpha=ang*cas
        cas=ranf()
        bet=ang*cas
        cas=ranf()
        gamma=ang*cas

        cosa=cos(alpha)
        cosb=cos(bet)
        cosg=cos(gamma)

        sena=sin(alpha)
        senb=sin(bet)
        seng=sin(gamma)

c     rotate the macroion and check for overlaps
        dx=x(n)-x(0)
        dy=y(n)-y(0)
        dz=z(n)-z(0)

        delx=(cosg*cosb*cosa-seng*sena)*dx
        delx=delx+(cosg*cosb*sena+seng*cosa)*dy+(-cosg*senb)*dz

        dely=(-seng*cosb*cosa-cosg*sena)*dx+(seng*senb)*dz
        dely=dely+(-seng*cosb*sena+cosg*cosa)*dy
        delz=(senb*cosa)*dx+(senb*sena)*dy+cosb*dz

        xn=x(0)+delx
        yn=y(0)+dely
        zn=z(0)+delz

	RETURN
	END
!	******************************************************************

	  SUBROUTINE DELTAE( IL, XN, YN, ZN, DER, DEK,CHECK)
!	*********************************************************************
!	**	Calculates energy change using Ewald method                  **
!	**	IL index of particle to be changed			 **
!	**	XN new x coordinate for particle IL 			 **
!	**	YN new y coordinate for particle IL 			 **
!	**	ZN		new x coordinate for particle IL 	 **
!	**	DER		energy change in R-space after the move	 **
!	**	DEK		energy change in K-space after the move	    **
!	**	dold	old distance between particle IL and other particles **
!	**	DNEW	new distance between particle IL and other particles **
!	*********************************************************************
        IMPLICIT NONE
        INTEGER  MAXK,TOTK,KMAX,KX,KY,KZ,I,KSQMAX,KSQ,IL,npmax,npart
		INTEGER NCHARGE1,CHECK 
        PARAMETER ( MAXK = 1000, KMAX = 5, KSQMAX = 27, npmax=5000)
        REAL         TWOPI, TWOPL, FACTOR, VD, VDN, RSQPI
        REAL       KVEC(MAXK), KAPPA, DEK,DER
        REAL       box,  hbox, DNEW(0:npmax), DOLD(0:npmax)
        REAL       sig(0:npmax), q(0:npmax)
        REAL       x(0:npmax), y(0:npmax), z(0:npmax), XN, YN, ZN
        REAL       E1, E2, XIL, YIL, ZIL, RXIJ, RYIJ, RZIJ
	    REAL       RIJSQ, RIJ, KRIJ
	    REAL       ERFC
        PARAMETER   ( TWOPI = 6.2831853, RSQPI = 0.5641896 )

        COMPLEX   EIKX(0:KMAX), EIKY(-KMAX:KMAX), EIKZ(-KMAX:KMAX)
        COMPLEX   DEIKX(0:KMAX),DEIKY(-KMAX:KMAX),DEIKZ(-KMAX:KMAX)
        COMPLEX   DEIKR, SUMI(MAXK), DSUMI(MAXK)

        COMMON / BLOCK2 / KVEC, KAPPA,twopl
	    COMMON / MOVE /   SUMI, VD, DSUMI, VDN
        COMMON / conf1 /x, y, z, npart,NCHARGE1
        COMMON /pot1/ sig,q
        COMMON / sys1 / box,hbox
	    COMMON /DISTAN/ DOLD, DNEW,XIL,YIL,ZIL


!    ** CALCULATE ENERGY CHANGE IN R-SPACE... HERE IL IS THE INDEX OF THE MICROION
!		THAT IS GOING TO BE DISPLACED.
		E1=0.
		E2=0.

		XIL=X(IL)
		YIL=Y(IL)
		ZIL=Z(IL)

		DO 5 I=0,NPART+ncharge1		
		IF(I.EQ.IL) GOTO 5	

! IL-th PARTICLE AT OLD POSITION

			RXIJ = X(I) - XIL
			RYIJ = Y(I) - YIL
			RZIJ = Z(I) - ZIL

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

	        RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
			RIJ=SQRT(RIJSQ)
              DOLD(I)=RIJ

		KRIJ=KAPPA*RIJ/BOX
	        E1=E1 + Q(I)*Q(IL)*ERFC(REAL(KRIJ))/RIJ
5		CONTINUE

!	NOW WE NEED TO FIND THE ENERGY AT THE NEW POSITION. SO THE NEW ENERGY IS 
!	THE ENERGY DUE TO THE DISPLACEMENT OF ONE MICROION INDEXED AS IL AND 
!	REARRANGING THE CHARGES ON THE SPHERE RANDOMLY. INODRDER TO GET A 
!	REARRANGED DISTRIBUTION OF CHARGES ON THE SPHERE WE CALL THE ARRANGECHARGES 
!	SUBROUTINE.

	if(il.lt.(npart+1)) then
		DO 6 I=0,NPART+ncharge1
		IF(I.EQ.IL) GOTO 6

!	 IL-th PARTICLE AW POSITION
			RXIJ = X(I) - XN
			RYIJ = Y(I) - YN
			RZIJ = Z(I) - ZN

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

	        RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
			RIJ=SQRT(RIJSQ)
	        DNEW(I)=RIJ

		      IF (RIJ .LT. (SIG(IL)+SIG(I))/2. ) THEN
		      DER = 300.
			  DEK=0.0
c	write(*, *)'Small ion-il,i=',il,i
c	read(*, *)
			CHECK=1
			  RETURN
		      ELSE
			  ENDIF

			KRIJ=KAPPA*RIJ/BOX
	        E2=E2 + Q(I)*Q(IL)*ERFC(REAL(KRIJ))/RIJ

6		CONTINUE
	else
                DO 77 I=1,NPART+ncharge1
                IF(I.EQ.IL) GOTO 77 

!        IL-th PARTICLE AW POSITION
                        RXIJ = X(I) - XN
                        RYIJ = Y(I) - YN
                        RZIJ = Z(I) - ZN

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

                RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
                        RIJ=SQRT(RIJSQ)
                DNEW(I)=RIJ

                      IF (RIJ .LT. (SIG(IL)+SIG(I))/2. ) THEN
                      DER = 300.
                          DEK=0.0
c	write(*, *)'Charge-il,i=',il,i
c	read(*, *)
			CHECK=1
                          RETURN
                      ELSE
                          ENDIF

                        KRIJ=KAPPA*RIJ/BOX
                E2=E2 + Q(I)*Q(IL)*ERFC(REAL(KRIJ))/RIJ

77               CONTINUE
	endif

      DER = E2 - E1

!    ** CALCULATE ENERGY CHANGE IN K-SPACE
!    ** CONSTRUCT EXP(IK.R) FOR IL-th ION AND K-VECTORS **
!    **   AT NEW AND OLD POSITIONS                      **
!    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY  **
       TWOPL=TWOPI / BOX

       EIKX(0)=(1.0,0.0)
       EIKY(0)=(1.0,0.0)
       EIKZ(0)=(1.0,0.0)
       EIKX(1)=CMPLX(COS(TWOPL*X(IL)),SIN(TWOPL*X(IL)))
       EIKY(1)=CMPLX(COS(TWOPL*Y(IL)),SIN(TWOPL*Y(IL)))
       EIKZ(1)=CMPLX(COS(TWOPL*Z(IL)),SIN(TWOPL*Z(IL)))
       EIKY(-1)=CONJG(EIKY(1))
       EIKZ(-1)=CONJG(EIKZ(1))

       DEIKX(0)=(1.0,0.0)
       DEIKY(0)=(1.0,0.0)
       DEIKZ(0)=(1.0,0.0)
       DEIKX(1)=CMPLX(COS(TWOPL*XN),SIN(TWOPL*XN))
       DEIKY(1)=CMPLX(COS(TWOPL*YN),SIN(TWOPL*YN))
       DEIKZ(1)=CMPLX(COS(TWOPL*ZN),SIN(TWOPL*ZN))
       DEIKY(-1)=CONJG(DEIKY(1))
       DEIKZ(-1)=CONJG(DEIKZ(1))

!    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **

		DO KX=2,KMAX
		DEIKX(KX)=DEIKX(KX-1)*DEIKX(1)
		EIKX(KX)=EIKX(KX-1)*EIKX(1)
		ENDDO

		DO KY=2,KMAX
		DEIKY(KY)=DEIKY(KY-1)*DEIKY(1)
		DEIKY(-KY)=CONJG(DEIKY(KY))
		EIKY(KY)=EIKY(KY-1)*EIKY(1)
		EIKY(-KY)=CONJG(EIKY(KY))
		ENDDO

		DO KZ=2,KMAX
		DEIKZ(KZ)=DEIKZ(KZ-1)*DEIKZ(1)
		DEIKZ(-KZ)=CONJG(DEIKZ(KZ))
		EIKZ(KZ)=EIKZ(KZ-1)*EIKZ(1)
		EIKZ(-KZ)=CONJG(EIKZ(KZ))
		ENDDO

!    ** SUM OVER ALL VECTORS **

      VDN = 0.0
      TOTK=0
      DO 30 KX=0,KMAX
      IF(KX.EQ.0) THEN
       FACTOR=1.0
      ELSE
       FACTOR=2.0
      ENDIF
       DO 20 KY=-KMAX,KMAX
        DO 10 KZ=-KMAX,KMAX
         KSQ=KX*KX+KY*KY+KZ*KZ
         IF((KSQ.LT.KSQMAX).AND.(KSQ.NE.0)) THEN
          TOTK=TOTK+1
        DEIKR=DEIKX(KX)*DEIKY(KY)*DEIKZ(KZ)-EIKX(KX)*EIKY(KY)*EIKZ(KZ)
        DSUMI(TOTK)=SUMI(TOTK)+Q(IL)*DEIKR
        VDN=VDN+FACTOR*KVEC(TOTK)*CONJG(DSUMI(TOTK))*DSUMI(TOTK)
         ENDIF
10      CONTINUE
20     CONTINUE
30    CONTINUE

		VDN = VDN / BOX
		DEK = VDN - VD
      RETURN
      END


!	****************************************************************

      subroutine toterg(gw, ener)

! ener (output) : total energy/kt
! gw: gaussian width parameter, adjust for efficiency of ewald sum


      IMPLICIT NONE
      INTEGER npmax,npart,NCHARGE1
      PARAMETER (npmax=5000)

      REAL     x(0:npmax),y(0:npmax),z(0:npmax)
      REAL     box,hbox
      REAL     sig(0:npmax), q(0:npmax)
      REAL     VR, VK, Ener, GW

      COMMON /sys1/  box,hbox
      COMMON /conf1/ x, y, z, npart,NCHARGE1
      COMMON /pot1/  sig,q

	  call setup(gw) 
      call rwald ( npart+ncharge1, vr, x, y, z, q, box)
      call kwald ( npart+ncharge1, vk, x, y, z, q, box )

      ener = vr + vk
c	write(*, *)'Ener=',ener
c	write(*, *)'vr+vk=',vr+vk
      end

!   *******************************************************************
!   ** real-space and reciprocal-space parts of ewald sum for ions.  **
!   **                                                               **
!   ** references:                                                   **
!   **                                                               *
!   ** woodcock and singer, trans. faraday soc. 67, 12, 1971.        **
!   ** de leeuw et al., proc. roy. soc. a 373, 27, 1980.             **
!   ** heyes, j. chem. phys. 74, 1924, 1981.                         **
!   ** see also fincham, mdions, ccp5 program library.               **
!   **                                                               **
!   ** routines supplied:                                            **
!   **                                                               **
!   ** subroutine setup                                              **
!   **    sets up the wavevectors for use in the ewald sum           **
!   ** subroutine rwald ( vr )                                       **
!   **    calculates the r-space part of the sum                     **
!   ** subroutine kwald ( vk )                                       **
!   **    calculates the k-space part of the sum                     **
!   ** real function erfc ( x )                                      **
!   **    returns the complementary error function                   **
!   **                                                               **
!   ** principal variables:                                          **
!   **                                                               **
!   **  totk         the total number of k-vectors stored            **
!   **  maxk         maximum possible number of k-vectors            **
!   **  kmax         max integer component of the k-vector           **
!   **  ksqmax       max square mod of the k-vector required         **
!   **  vr           energy from r-space sum/kt                      **
!   **  vk           energy from k-space sum/kt                      **
!   **  kvec(maxk)   array used to store k-vectors                   **
!   **  kappa        width of cancelling distribution                **
!   **                                                               **
!   **  (kappa is a input value that should be determined by 	    **
!   **   computational efficiency, kappa=5/box for initial try       **
!   **   kappa=(alpha)**0.5 in frenkel & smit's notation)            **
!   **                                                               **
!   ** usage:                                                        **
!   **                                                               **
!   ** setup is called once at the beginning of the simulation       **
!   ** to calculate all the k-vectors required in the ewald sum.     **
!   ** these vectors are used throughout the simulation in the       **
!   ** subroutine kwald to calculate the k-space contribution to the **
!   ** potential energy at each configuration. the self XTERM IS

!   ** subtracted from the k-space contribution in kwald.            **
!   ** the surface term for simulations in vacuum is not included.   **
!   ** routine rwald returns the r-space contribution to the ewald   **
!   ** sum and is called for each configuration in the simulation.   **
!   ** a cubic box and unit box length are assumed throughout.       **
!   *******************************************************************



        subroutine setup(gw) 

!   *******************************************************************
!   ** routine to set up the wave-vectors for the ewald sum.         **
!   ** gw: gaussian width parameter                                  *
!   ** in this example we allow a maximum of 1000 wavevectors.       **
!   *******************************************************************
        integer     maxk
        parameter ( maxk = 1000 )
        real      kvec(maxk), kappa, gw
        integer     kmax, ksqmax, ksq, kx, ky, kz, totk
        real     twopi, b, rkx, rky, rkz, rksq,twopl
        parameter ( kmax = 5, ksqmax = 27 , twopi = 6.2831853 )

        common / block2 / kvec, kappa,twopl
!   *******************************************************************
	  kappa = gw
        b = 1.0 / 4.0 / kappa / kappa

!   ** loop over k-vectors. note kx is non-negative **

        totk = 0

        do 100 kx = 0, kmax

           rkx = twopi * float ( kx )

           do 99 ky = -kmax, kmax

              rky = twopi * float ( ky )

              do 98 kz = -kmax, kmax

                 rkz = twopi * float ( kz )

                 ksq = kx * kx + ky * ky + kz * kz

                 if ( ( ksq .lt. ksqmax ) .and. ( ksq .ne. 0 ) ) then

                    totk = totk + 1

                    if ( totk .gt. maxk ) stop 'kvec is too small'

                    rksq = rkx * rkx + rky * rky + rkz * rkz
                    kvec(totk) = twopi * exp ( -b * rksq ) / rksq

                 endif

98            continue

99         continue

100     continue

        write( *, ' ( '' ewald sum setup complete ''     ) ' )
        write( *, ' ( '' number of wavevectors is '', i5 ) ' ) totk

        return
        end



        subroutine rwald ( n, vr ,rx, ry, rz, q, box)

!   *******************************************************************
!   ** calculates r-space part of potential energy by ewald method.  **
!   **                                                               **
!   ** principal variables:                                          **
!   **                                                               **
!   **         n                     number of ions (=npart)         **
!   **         rx(n),ry(n),rz(n)     positions of ions/lb            **
!   **         q(n)                  ionic valence                   **
!   **         vr                    r-space potential energy/kt     **
!   **         box                   simulation box length/lb        **
!   **                                                               **
!   ** routine referenced:                                           **
!   **                                                               **
!   ** function erfc ( x )                          **
!   **    returns the complementary error function                   **
!   *******************************************************************
        integer     maxk, npmax
        parameter ( maxk = 1000, npmax=5000 )
        real      kvec(maxk), kappa, twopl

        integer              n, i, j
        real      rx(0:npmax), ry(0:npmax), rz(0:npmax), q(0:npmax)
        real         vr, box
        real         rxi, ryi, rzi, zi, rxij, ryij, rzij
        real         rijsq, rij, krij, vij
    	real erfc

        common / block2 / kvec, kappa,twopl
!   *******************************************************************
	  vr = 0.0
        do 100 i =1,n-1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           zi  = q(i)

           do 99 j = i+1,n

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)

              rxij = rxij - box * anint ( rxij /box )
              ryij = ryij - box * anint ( ryij /box )
              rzij = rzij - box * anint ( rzij /box )
              rijsq = rxij * rxij + ryij * ryij + rzij * rzij
              rij   = sqrt ( rijsq )
              krij  = kappa * rij/ box
              vij   = zi*q(j)*erfc(real(krij))/rij

              vr    = vr + vij
99      continue
100     continue
        return
        end



        subroutine kwald ( n, vk, rx, ry, rz, q, box )

!   *******************************************************************
!   ** calculates k-space part of potential energy by ewald method.  **
!   **                                                               **
!   ** the self term is subtracted.                                  **
!   ** in one coordinate direction (x), symmetry is used to reduce   **
!   ** the sum to include only positive k-vectors.                   **
!   ** the negative vectors in this direction are included by use    **
!   ** of the multiplicative variable 'factor'.                      **
!   **                                                               **
!   ** principal variables:                                          **
!   **                                                               **
!   **  n                   number of ions                           **
!   **  rx(n),ry(n),rz(n)   positions of ions, in dimension of lb    **
!   **  q(n)                ionic valences                           **
!   **  vk                  k-space potential energy/kt              **
!   **  vks                 self part of k-space sum/kt              **
!   **  box                 simulation box length/lb                 **
!   *******************************************************************
        integer  maxk, n, totk, kmax, kx, ky, kz, i, ksqmax, ksq, npmax
        parameter   ( maxk = 1000, kmax = 5, ksqmax = 27, npmax=5000 )
        real       kvec(maxk)
        real      rx(0:npmax), ry(0:npmax), rz(0:npmax), q(0:npmax)
        real       kappa, vk, box
        real         twopi, twopl, factor, vd, vs, rsqpi, vdn
        parameter   ( twopi = 6.2831853, rsqpi = 0.5641896 )

        complex       eikx(1:npmax, 0:kmax)
        complex       eiky(1:npmax, -kmax:kmax)
        complex       eikz(1:npmax, -kmax:kmax)
        complex       eikr(1:npmax), sumi(maxk), dsumi(maxk)

        common / block2 / kvec, kappa,twopl
	 common / move /   sumi, vd, dsumi, vdn
!   *******************************************************************

!   ** construct dexp(ik.r) for all ions and k-vectors **

!   ** calculate kx, ky, kz = 0 , -1 and 1 dexplicitly **
        TWOPL=TWOPI / BOX

        do 10 i = 1,n

           eikx(i, 0) = (1.0, 0.0)
           eiky(i, 0) = (1.0, 0.0)
           eikz(i, 0) = (1.0, 0.0)

!           eikx(i, 1) = cmplx(cos(twopl*rx(i)),sin(twopl*rx(i)))
!           eiky(i, 1) = cmplx(cos(twopl*ry(i)),sin(twopl*ry(i)))
!           eikz(i, 1) = cmplx(cos(twopl*rz(i)),sin(twopl*rz(i)))
           EIKX(I, 1) = CMPLX ( COS ( TWOPL * RX(I) ) ,
     :                           SIN ( TWOPL * RX(I) ) )
           EIKY(I, 1) = CMPLX ( COS ( TWOPL * RY(I) ) ,
     :                           SIN ( TWOPL * RY(I) ) )
           EIKZ(I, 1) = CMPLX ( COS ( TWOPL * RZ(I) ) ,
     :                           SIN ( TWOPL * RZ(I) ) )
           eiky(i, -1) = conjg ( eiky(i, 1) )
           eikz(i, -1) = conjg ( eikz(i, 1) )

10      continue

!   ** calculate remaining kx, ky and kz by recurrence **

        do 12 kx = 2, kmax

           do 11 i = 1,n

              eikx(i, kx) = eikx(i, kx-1) * eikx(i, 1)

11         continue

12      continue

        do 14 ky = 2, kmax

           do 13 i =1,n 

              eiky(i,  ky) = eiky(i, ky-1) * eiky(i, 1)
              eiky(i, -ky) = conjg ( eiky(i, ky) )

13         continue

14      continue

        do 16 kz = 2, kmax

           do 15 i = 1,n

              eikz(i,  kz) = eikz(i, kz-1) * eikz(i, 1)
              eikz(i, -kz) = conjg ( eikz(i, kz) )

15         continue

16      continue

!   ** sum over all vectors **

        vd   = 0.0
        totk = 0

        do 24 kx = 0, kmax

        if ( kx .eq. 0 ) then

        factor = 1.0

        else

        factor = 2.0

        endif

        do 23 ky = -kmax, kmax

        do 22 kz = -kmax, kmax

        ksq = kx * kx + ky * ky + kz * kz

      if ( ( ksq .lt. ksqmax ) .and. ( ksq .ne. 0 ) ) then

       totk = totk + 1
        sumi(totk)  = (0.0, 0.0)

        do 21 i = 1,n

         eikr(i) = eikx(i, kx) * eiky(i, ky) * eikz(i, kz)
         sumi(totk) = sumi(totk) + q(i) * eikr(i)

21       continue

         vd=vd+factor*kvec(totk)*conjg(sumi(totk))*sumi(totk)

         endif

22       continue

23       continue

24       continue

	 vd=vd / box

!   ** calculates self part of k-space sum **

        vs = 0.0

        do 25 i =1,n

           vs = vs + q(i) * q(i)

25      continue

        vs = rsqpi * kappa * vs / box

!   ** calculate the total k-space potential **

        vk = vd - vs
        return
        end

!       ************************************
