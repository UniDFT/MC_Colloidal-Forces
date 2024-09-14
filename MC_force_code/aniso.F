       program  main
       implicit none

#include <mpif.h>

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
       call mpi_finalize(ierror)
       end

	subroutine pforce (myid) 
      implicit none
      integer equil, prod, ii, icycl, ndispl, attempt,nstart,
     *   attemptpivot,nacc,naccpivot,ncycl,nmoves, imove, ibeg,
     *   npart, ncharges,nmax,myid
      parameter (nmax=300) 
      real     ener, enercheck, dr, gw,deltar,angle,cas,ahelp,ranf
	real a,b,c,x(-8:nmax),y(-8:nmax),z(-8:nmax)
	logical alfa(nmax),beta(nmax)
      common /conf1/ x,y,z,npart,ncharges,a,b,c,alfa,beta

      call readdat(myid,ibeg,equil, prod, ndispl, dr,
     *	angle,gw,deltar)
	if (ibeg.lt.2) then
		nstart=1
	else
		nstart=2
	end if
	 
      nmoves = ndispl
c     number of trial moves before sampling the force

!    ---total energy of the system
      call toterg(gw, ener)
      write (*, 99001) ener

!    ---start mc-cycle
      do  ii = nstart, 2
!       --- ii=1 equilibration
!       --- ii=2 production
         if (ii.eq.1) then
            ncycl = equil
		  write (*, *) ' start equilibration '
	   else
            ncycl = prod
		  write (*, *) ' start production '
         endif
c     the procedure "sample" is called now to initialize the sampling variables
c     as well as to save the initial configuration of the system
	   call sample(myid,0, ener,dr,angle,deltar)
c        initialize running average quantities
	   attempt = 0
	   nacc = 0
	   attemptpivot=0
	   naccpivot=0
c        intialize the subroutine that adjust the maximum displacement and the maximum 
c        rotation of a macroion 
         call adjust(attempt, nacc, dr)
	   call adjustangle(attemptpivot,naccpivot,angle)
	   ahelp=float(ncharges)/float(ncharges+npart)	
         do icycl = 1, ncycl
              do imove = 1, nmoves
	            cas=ranf()
				if (cas.le.ahelp) then
				  call pivot(ener,attemptpivot,naccpivot,angle)
				else
				  call mcmove(ener, attempt, nacc, dr)
				end if
			end do
				if (ii.eq.2) then
c					the force is sampled once every nmoves times a mc move is tryed
c					in this procedure there's the instruction to
c					sample the force interacting between the two macroions  
	                 call sample(myid,1, ener,dr,angle,deltar)
c					this procedure also stores the intermediate configuration
				end if
             if (mod(icycl,ncycl/10).eq.0) then
			write (*, *) '======>> done ', icycl,' out of ', ncycl
c				adjust displacements
				call adjust(attempt, nacc, dr)
				call adjustangle(attemptpivot,naccpivot,angle)
             end if
		end do
            write (*, 99003) attempt, nacc, 
     &                               100.*float(nacc)/float(attempt)
            write (*, 99004) attemptpivot,  
     &              naccpivot, 100.*float(naccpivot)/float(attemptpivot)
	end do
c     now, the procedure "sample" is called to export the results

	call sample(myid,2, ener,dr,angle,deltar)
	write(*,*) 'final average energy = ', ener   
	call toterg (gw,enercheck)
	write (*,*) 'final configurational energy=', enercheck

      stop
 
99001 format (' total energy at initial configuration: ', f12.5, /)
99002 format (' total energy at the end of simulation: ', f12.5, /, 
     &        '       running energy              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5, /) 
99003 format (' number of att. to displ. a part.  : ', i10, /, 
     &        ' success: ', i10, '(= ', f6.2, '%)')
99004 format (' number of att. to rotate a macroi.: ', i10, /, 
     &        ' success: ', i10, '(= ', f6.2, '%)')
      end

      subroutine readdat(myid,ibeg,equil,prod,ndispl,dr,angle,gw,deltar)

c     ndispl is the number of mc moves before a force sampling
!    ---input solution condition and model parameters
!
!    ---input parameters: file: (3), input.dat
!   ibeg  =  0 : initialize 
!            1 : read configuration from disk	
!            1 : restart sampling
!   equil      : number of monte carlo cycles during equilibration
!   prod       : number of monte carlo cycles during production
!   dr         : maximum displacement
c   angle      : maximum rotation angle 
!   npart      : total number of particles
!   nparti(i)  : number of particles for species i
!
!   nion   = number of ion species
!   nionmax= maxmum number of ion species
!   siga,b,c(i) = square collision distance between microion and macroion/bjerrum length
!   sigi(i)= diameters of ions and macroions species, in angstroms
!   q(i)   = valence for a microion
!   qi(i)  = valence for an ion, for input
!   lb0    = bjerrum length for water at 298.15 k, in angstroms
!   lb     = bjerrum length at system temperature, in angstroms
!
!   box    = simulation box length in angstroms (for input)
!   rab    = separation between two macroions/main diameter
!    ---input parameters: file: iniconf.dat (restart file
!               to continue a simulation from disk)
!   npmax  = maximum number of particles
!   dr     = optimized maximum displacement old configurations
!   gw     = gaussian width parameter, adjust for efficiency of ewald sum
!   x(1),y(1),z(1)            : position microion 1
!       ...
!   x(npart),y(npart),z(npart): position particle last particle
!   x(-1),y(-1),z(-1)         : center of the first macroion
!   x,y,z(-2,-3,-4)			: position of three convex points of macroion 1
!   x(-5),y(-5),z(-5)         : center of the second macroion
!   x,y,z(-6,-7,-8)			: position of three convex points of macroion 2

      implicit none
      integer ibeg,equil,prod,i,j,k,ndispl,maxk,
     &        nion,nionmax,npmax,npart,zero,ns0,myid,nnode 
      parameter (npmax=300, nionmax=10,maxk=1000)
      integer  nparti(nionmax),ncharges, qi(nionmax)
      real box, hbox,con2, dr, lb0, lb, gw, shbox,deltar
	real x(-8:npmax),	y(-8:npmax), z(-8:npmax)
      real sigb(npmax), sigc(npmax),q(npmax),angle,siga(npmax)
      real sig(npmax),sigi(nionmax) 
      real rab,a,b,c,pi,a0,b0,c0,rab0(20)
	logical alfa(npmax),beta(npmax)
      real Enert0, MF0, MFHS0, MFCC0,phcc0
	parameter (con2=0.0006022137,lb0=7.144464) 

      common /rest1/ ns0, Enert0, MF0, MFHS0, MFCC0,phcc0,rab0
      common /sys1/  box,hbox,rab
      common /pot1/  sig,siga,sigb,sigc,q
      common /conf1/ x,y,z,npart,ncharges,a,b,c,alfa,beta
      character filename(20)*5,filename2(20)*5
      character output1*5,iniconf*5
      data filename/'out1','out2','out3','out4','out5','out6'
     &      ,'out7','out8','out9','out10','out11','out12'
     &      ,'out13','out14','out15','out16','out17','out18'
     &      ,'out19','out20'/
      data filename2/'inp1','inp2','inp3','inp4','inp5','inp6'
     &      ,'inp7','inp8','inp9','inp10','inp11','inp12'
     &      ,'inp13','inp14','inp15','inp16','inp17','inp18'
     &      ,'inp19','inp20'/

        iniconf=filename2(myid+1)
        output1=filename(myid+1) 
      open(65,file='input')
!    ---read simulation parameters
      read (65, *)
      read (65, *) nnode       
      read (65, *)
      read (65, *) ibeg, equil, prod
      read (65, *)
      read (65, *) dr, angle, gw, deltar
      read (65, *)
      read (65, *) ndispl
      read (65, *)
      read (65, *) box
      read (65, *)
      read (65, *) (rab0(i),i=1,nnode)          
!    ---read solution condition and model parameters
c    ---nion represents only microions
      read (65, *)
      read (65, *) nion
      read (65, *)
      do i=1,nion
	read (65, *) sigi(i), qi(i), nparti(i)
      enddo
c    --- main radii of ellipsoid a0,b0,c0
	read(65,*)
	read(65,*) a0,b0,c0
      close(65)


c assume that each macroion has 20 unit charges
	ncharges=20
c     check charge neutrality at initial condition
	zero=-40
	do i=1,nion
			zero=zero+qi(i)*nparti(i)
	enddo
	if (zero.ne.0) then
			write(*,*)'charge neutrality is unsatisfied 
     *		at initial condition!' 
			stop
	endif

!    ---calculate parameters:
c     pi,bjerrum length
      pi = 4.*atan(1.0e0)
	lb=lb0

!    --- normalized all lenghtes (Angstrom)with bjerrum length
	box = box/lb
        hbox = box/2.0
	shbox = hbox*hbox
	a=a0/lb
	b=b0/lb
	c=c0/lb

c diameter,charge, and square collision diameter for each microion
c total number of small ions

	npart=0
	k=0
	do i = 1,nion
		npart=npart+nparti(i)
		do j = 1, nparti(i)
			k=k+1
			sig(k)=sigi(i)/lb
			q(k) = float(qi(i))
			siga(k) = (a + sig(k)/2.0)**2
			sigb(k) = (b + sig(k)/2.0)**2
			sigc(k) = (c + sig(k)/2.0)**2
		end do
      end do

      if (npart.gt.npmax) then
         write (*, *) ' error: number of particles too large'
         stop
      end if

c separation between two macroions
	rab=rab0(myid+1)*2.0*max(a,b,c)

!    ---read/generate initial configuration
      if (ibeg.eq.0) then
         call config 
         open(64,file=iniconf)
         write (64, *) npart, dr, angle
	 ns0=0
         Enert0=0.0
         MF0=0.0
         MFHS0=0.
         MFCC0=0.
         phcc0=0.
         write (64,*) ns0, Enert0, MF0, MFHS0, MFCC0,phcc0

         do i = 1, npart
			write (64,*) alfa(i),beta(i),x(i),y(i),z(i)
         end do 
	   do i=-1,-8,-1
		   write(64,*) x(i), y(i), z(i)
	   end do
	   do i=npart+1,npart+2*ncharges
		   write(64,*) q(i),x(i),y(i),z(i)
	   end do
         close(64)

      elseif (ibeg.eq.1) then
         write (*, *) ' read conf from disk '
	   open(64,file=iniconf)
	   read (64, *) npart, dr, angle
         read (64,*) ns0, Enert0, MF0, MFHS0, MFCC0,phcc0
	   ns0=0
         Enert0=0.0
         MF0=0
         MFHS0=0
         MFCC0=0
         phcc0=0
         do i = 1, npart
			read (64, *) alfa(i),beta(i),x(i),y(i),z(i)
         end do
		do i=-1,-8,-1
			read(64,*) x(i),y(i),z(i)
		end do
		do i=npart+1,npart+2*ncharges
			read(64,*) q(i),x(i),y(i),z(i)
		end do
	    close(64)

	else if (ibeg.eq.2) then

		write(*,*) 'restart sampling'
	    open(64,file=iniconf)
	    read(64,*) npart, dr, angle
          read (64,*) ns0, Enert0, MF0, MFHS0, MFCC0,phcc0
          do i = 1, npart
			read (64, *) alfa(i),beta(i),x(i),y(i),z(i)
          end do
		do i=-1,-8,-1
			read(64,*) x(i),y(i),z(i)
		end do
		do i=npart+1,npart+2*ncharges
			read(64,*) q(i),x(i),y(i),z(i)
		end do
	    close(64)
	end if

!    ---write input data
      write (*,  1) equil, prod
      write (*,  2) ndispl, dr
      write (*,  3) npart, box, lb
      write (*, *) '#ions    diameter(A)   valence' 
	do i=1,nion
		write (*,  4) nparti(i), sigi(i), qi(i)
	enddo 
	write(*,*) 'main radii of ellipsoid' ,a,b,c
	write(*,*) 'macroions charges',ncharges
 1    format ('  number of equilibration cycles             :', i10, /, 
     &        '  number of production cycles                :', i10, /) 
 2    format ('  number of att. to displ. a part. per cycle :', i10, /, 
     &        '  maximum displacement                       :', f10.3, 
     &        //)
 3    format ('  total number of particles                  :', i10, /, 
     &        '  box length/bjerrum length                  :', f10.3/,
     &        '  bjerrum length                             :', f10.6/)

 4    format (2x,i5,2x,f8.5,2x,i6)
      end

	subroutine config 
c     initial configurations of macroions and microions
c     assume that macroions have ellipsoid shape with semiaxes a, b, c.
c     for spherical shape, a=b=c=radius
c     the body coordinate systemorigins at macroion center with the x, y and z in the semiaxe directions.
c     x,y,z(-1,-8): centers and main axes of macroions
c     x,y,z(1,npart): positions of small ions
c     x,y,z(npart+1,npart+20): locations of the unit charges
c     alfa, beta(i) labels whether or not microion i is within the collision diameter of macroions 
c This subroutine returns positions of ions and charges 
      implicit none
	integer i,k,npart,npmax,ncharges,ii,jj
	parameter(npmax=300)
	real x(-8:npmax), y(-8:npmax), z(-8:npmax), box, hbox,
     * sig(npmax),siga(npmax),sigb(npmax),sigc(npmax),q(npmax)
     * ,xx, xd, yd, zd, rik,yy,zz,check,
     * a, b, c, rab,ranf,xc(4,5,3)
     	logical alfa(npmax),beta(npmax)
      common /conf1/ x,y,z,npart,ncharges,a,b,c,alfa,beta
	common /sys1/ box, hbox, rab
	common /pot1/ sig, siga, sigb, sigc, q


c    set body coordinates
c   --- assign LOCATIONs of two macroions centers (macroion -1 on the top)
		x(-1)= rab/sqrt(12.e0)
		y(-1)= rab/sqrt(12.e0)
		z(-1)= rab/sqrt(12.e0)

		x(-5)=-rab/sqrt(12.e0)
		y(-5)=-rab/sqrt(12.e0)
		z(-5)=-rab/sqrt(12.e0)

c	locations of 6 main axes:
		x(-2)= x(-1) + a
		y(-2)= y(-1)	
		z(-2)= z(-1)
		x(-3)= x(-1) 
		y(-3)= y(-1) + b	
		z(-3)= z(-1)
		x(-4)= x(-1) 
		y(-4)= y(-1)	
		z(-4)= z(-1) + c

		x(-6)= x(-5) + a
		y(-6)= y(-5)	
		z(-6)= z(-5)
		x(-7)= x(-5) 
		y(-7)= y(-5) + b	
		z(-7)= z(-5)
		x(-8)= x(-5) 
		y(-8)= y(-5)	
		z(-8)= z(-5) + c

c     charge distribution of macroions (20 charges per macroion)
		ncharges = 20
		call dodeca(xc,2.e0)
		do i=1,20
                q(npart+20+i)=-1.e0
                q(npart+i)=-1.e0
		enddo
		i=0
		do 1 ii=1,4
		do 1 jj=1,5
		i=i+1
		x(npart+i)    = x(-1) + xc(ii,jj,1)*a
		y(npart+i)    = y(-1) + xc(ii,jj,2)*b
		z(npart+i)    = z(-1) + xc(ii,jj,3)*c
                x(npart+20+i) = x(-5) + xc(ii,jj,1)*a
                y(npart+20+i) = y(-5) + xc(ii,jj,2)*b
                z(npart+20+i) = z(-5) + xc(ii,jj,3)*c
1		continue

c	create a random configuration of microions
          do i=npart,1, -1
3		xx=ranf()
		x(i)=hbox*(xx-0.5)*2.0
		xx=ranf()
		y(i)=hbox*(xx-0.5)*2.0
		xx=ranf()
		z(i)=hbox*(xx-0.5)*2.0

c check collision with the macroions 
		do k=-1,-5,-4
			xx=(x(i)-x(k))**2
			yy=(y(i)-y(k))**2
			zz=(z(i)-z(k))**2
             if ((xx+yy+zz).lt.max(siga(i),sigb(i),sigc(i))) then
				if (k.eq.-1) then
					alfa(i)=.true.
				else
					beta(i)=.true.
				end if
				check=xx/siga(i)+yy/sigb(i)+zz/sigc(i)-1.e0
				if (check.le.0) goto 3
			else
				if (k.eq.-1) then
					alfa(i)=.false.
				else
					beta(i)=.false.
				end if
			end if													
		end do

		do k=i+1, npart
			xd=x(i)-x(k)
			if(xd.gt.hbox)  xd=xd-box
			if(xd.lt.-hbox) xd=xd+box
			yd=y(i)-y(k)
			if(yd.gt.hbox)  yd=yd-box
			if(yd.lt.-hbox) yd=yd+box
			zd=z(i)-z(k)
			if(zd.gt.hbox)  zd=zd-box
			if(zd.lt.-hbox) zd=zd+box
			rik=sqrt(xd*xd+yd*yd+zd*zd)
			if (rik .lt. (sig(i)+sig(k))) goto 3
		end do
	end do
      end

      subroutine mcmove(ener, attempt, nacc, dr)
!   *******************************************************************
!    attempts to displace a randomly selected particle
! ener   (input/output) : total energy
! attemp (input/output) number of attemps that have been
!                 performed to displace a particle
! nacc   (input/output) number of successful attemps
!                 to displace a particle
! dr     (input) maximum displacement
!   *******************************************************************
      implicit none
      integer il, attempt, nacc,  npmax, npart, maxk, i,ncharges
      parameter (npmax=300, maxk = 1000)
      real     box,hbox
      real     x(-8:npmax), y(-8:npmax), z(-8:npmax)
      real     ener, xn, yn, zn, der, dek, dr,rab,ranf
      real     vd, vdn,xept,a,b,c, q(npmax), cas 
      real     sig(npmax),siga(npmax),sigb(npmax),sigc(npmax)
	logical    alfa(npmax), beta(npmax), alfan, betan
      complex   sumi(maxk),dsumi(maxk)

      common /pot1/  sig,siga,sigb,sigc,q
      common / sys1 / box, hbox, rab
      common / conf1 /  x, y, z, npart,ncharges,a,b,c,alfa,beta
	common / move /   sumi, vd, dsumi, vdn

	attempt = attempt + 1
!    ---select a microion at random
	cas=ranf()
      il = int(npart*cas) + 1
!    ---give the particle a random displacement
	cas=ranf()
      xn = x(il) + (cas-0.5)*dr
	cas=ranf()
      yn = y(il) + (cas-0.5)*dr
	cas=ranf()
      zn = z(il) + (cas-0.5)*dr

!    ---calculate energy change
	call deltae( il, xn, yn, zn, alfan, betan, der, dek )
	
!    ---acceptance test
	if(der.gt.299.0) goto 10
	xept=exp(-der-dek)
	cas=ranf()
      if (cas .lt. xept) then
!       --accepted
		nacc = nacc + 1
		vd=vdn
		do i=1,maxk
			sumi(i) = dsumi(i)
		enddo
		ener = ener + der + dek

!       ---put particle in simulation box
		if (xn.lt.-hbox) xn = xn + box
		if (xn.gt.hbox) xn = xn - box
		if (yn.lt.-hbox) yn = yn + box
		if (yn.gt.hbox) yn = yn - box
		if (zn.lt.-hbox) zn = zn + box
		if (zn.gt.hbox) zn = zn - box

		x(il) = xn
		y(il) = yn
		z(il) = zn
		alfa(il)=alfan
		beta(il)=betan
	else
!		 --- rejected
      end if
10	continue
      end

	  subroutine deltae( il, xn, yn, zn, alfan, betan, der, dek )
c     this procedure is to be used when a microion is displaced
!*********************************************************************
!**	calculates energy change using ewald method                  **
!**	il		index of particle to be changed						 **
!**	xn		new x coordinate for particle il 					 **
!**	yn		new y coordinate for particle il 					 **
!**	zn		new x coordinate for particle il 					 **
!**	der		energy change in r-space after the move	    		 **
!**	dek		energy change in k-space after the move	 
c     alfa(i) indicates when a microparticle is close to the 1st macroion   	
c     beta(i) indicates when a microparticle is close to the 2nd macroion
c     alfan and betan have the same meaning referred to the new position of 
c     the microparticle being moved
!*********************************************************************
        implicit none
        integer     maxk, totk, kmax, kx, ky, kz, i, ksqmax, ksq, 
     &              il, npmax, npart,ncharges,part

c     ncharge is the number of charges in a macroion  
        parameter   ( maxk = 1000, kmax = 5, ksqmax = 27, npmax=300)
        real         twopi, twopl, factor, vd, vdn
        real       kvec(maxk), kappa, dek,der
        real       box, hbox
        real       sig(npmax),siga(npmax),sigb(npmax),
     *	   sigc(npmax),q(npmax), a,b,c
        real       x(-8:npmax), y(-8:npmax), z(-8:npmax), xn, yn, zn
        real       e1, e2, xil, yil, zil, rxij, ryij, rzij,xm1,xm2,xm3
	  real       rijsq, rij, krij,ym1,ym2,ym3,zm1,zm2,zm3,xx,yy,zz
	  real       erfuc,rab
	  logical 	   alfa(npmax), alfan, beta(npmax), betan
        parameter   ( twopi = 6.2831853)

        complex eikx(0:kmax), eiky(-kmax:kmax), eikz(-kmax:kmax)
        complex deikx(0:kmax),deiky(-kmax:kmax),deikz(-kmax:kmax)
        complex deikr, sumi(maxk), dsumi(maxk)

        common / block2 / kvec, kappa,twopl
	  common / move /   sumi, vd, dsumi, vdn
        common / conf1 /  x, y, z, npart, ncharges, a,b,c, alfa, beta
        common / pot1 /   sig,siga, sigb, sigc,q
        common / sys1 /   box, hbox, rab

		xil=x(il)
		yil=y(il)
		zil=z(il)
c     check overlaps with the macroions 
         do part=1,2
			xx=xn-x(-1+((part-1)*(-4)))
			yy=yn-y(-1+((part-1)*(-4)))
			zz=zn-z(-1+((part-1)*(-4)))
	        if((xx**2+yy**2+zz**2).le.
     &             (max(siga(il),sigb(il),sigc(il)))) then
			       if (part.eq.1) then
					alfan=.true.
				else if (part.eq.2) then
					betan=.true.
				end if
c     check collisions
				xm1=x(-2+((part-1)*(-4)))-x(-1+((part-1)*(-4)))
				ym1=y(-2+((part-1)*(-4)))-y(-1+((part-1)*(-4)))
				zm1=z(-2+((part-1)*(-4)))-z(-1+((part-1)*(-4)))

				xm2=x(-3+((part-1)*(-4)))-x(-1+((part-1)*(-4)))
				ym2=y(-3+((part-1)*(-4)))-y(-1+((part-1)*(-4)))
				zm2=z(-3+((part-1)*(-4)))-z(-1+((part-1)*(-4)))

				xm3=x(-4+((part-1)*(-4)))-x(-1+((part-1)*(-4)))
				ym3=y(-4+((part-1)*(-4)))-y(-1+((part-1)*(-4)))
				zm3=z(-4+((part-1)*(-4)))-z(-1+((part-1)*(-4)))
c     now the program is computing the distance of a microion from the macroion surface
				rxij = (xx*xm1 + yy*ym1 + zz*zm1)/a
				ryij = (xx*xm2 + yy*ym2 + zz*zm2)/b
				rzij = (xx*xm3 + yy*ym3 + zz*zm3)/c
				rijsq=(rxij*rxij / siga(il) + ryij*ryij / 
     *				sigb(il) + rzij*rzij / sigc(il))-1.e0
				if (rijsq .lt. 0.e0 )  then
					der=300.
					dek=0.
					return
				end if
			else
	   				if (part.eq.1) then
					alfan=.false.
					else if (part.eq.2) then
					betan=.false.
					end if
			end if
		end do

c     computation of the variation in energy due to a microion displacement
		e1=0.
		e2=0.
c     computation of electrostatic energy in the new position
		do  i=1,npart
		if(i.eq.il) then
	        e2=e2
	    else
			rxij = x(i) - xn
			ryij = y(i) - yn
			rzij = z(i) - zn

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

	        rijsq=rxij*rxij + ryij*ryij + rzij*rzij
			rij=sqrt(rijsq)

		    if (rij .lt. (sig(il)+sig(i))) then
				der = 300.
				dek=0.0
				return
			endif
			krij=kappa*rij/box
	        e2=e2 + q(i)*q(il)*erfuc(krij)/rij
		endif
		end do


		do  i=npart+1,npart+2*ncharges
			rxij = x(i) - xn
			ryij = y(i) - yn
			rzij = z(i) - zn

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

	        rijsq=rxij*rxij + ryij*ryij + rzij*rzij
			rij=sqrt(rijsq)

			krij=kappa*rij/box
	        e2=e2 + q(i)*q(il)*erfuc(krij)/rij
		end do

c     computation of electrostatic energy in the old position
		do i=1,npart+2*ncharges
		if(i.eq.il) then
			e1=e1
		else
			rxij = x(i) - xil
			ryij = y(i) - yil
			rzij = z(i) - zil

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

	        rijsq=rxij*rxij + ryij*ryij + rzij*rzij
			rij=sqrt(rijsq)

			krij=kappa*rij/box
	        e1=e1 + q(i)*q(il)*erfuc(krij)/rij
		endif
		end do

      der = e2 - e1

!    ** calculate the energy change in k-space
	 twopl=twopi/box
       eikx(0)=(1.0,0.0)
       eiky(0)=(1.0,0.0)
       eikz(0)=(1.0,0.0)
       eikx(1)=cmplx(cos(twopl*x(il)),sin(twopl*x(il)))
       eiky(1)=cmplx(cos(twopl*y(il)),sin(twopl*y(il)))
       eikz(1)=cmplx(cos(twopl*z(il)),sin(twopl*z(il)))
       eiky(-1)=conjg(eiky(1))
       eikz(-1)=conjg(eikz(1))

       deikx(0)=(1.0,0.0)
       deiky(0)=(1.0,0.0)
       deikz(0)=(1.0,0.0)
       deikx(1)=cmplx(cos(twopl*xn),sin(twopl*xn))
       deiky(1)=cmplx(cos(twopl*yn),sin(twopl*yn))
       deikz(1)=cmplx(cos(twopl*zn),sin(twopl*zn))
       deiky(-1)=conjg(deiky(1))
       deikz(-1)=conjg(deikz(1))

!    ** calculate remaining kx, ky and kz by recurrence **
!
	do kx=2,kmax
		deikx(kx)=deikx(kx-1)*deikx(1)
		eikx(kx)=eikx(kx-1)*eikx(1)
	enddo

	do ky=2,kmax
		deiky(ky)=deiky(ky-1)*deiky(1)
		deiky(-ky)=conjg(deiky(ky))
		eiky(ky)=eiky(ky-1)*eiky(1)
		eiky(-ky)=conjg(eiky(ky))
	enddo

	do kz=2,kmax
		deikz(kz)=deikz(kz-1)*deikz(1)
		deikz(-kz)=conjg(deikz(kz))
		eikz(kz)=eikz(kz-1)*eikz(1)
		eikz(-kz)=conjg(eikz(kz))
	enddo
!
!    ** sum over all vectors **

      vdn = 0.0
      totk=0
      do kx=0,kmax
		if(kx.eq.0) then
			factor=1.0
		else
			factor=2.0
		endif
		do ky=-kmax,kmax
			do  kz=-kmax,kmax
				ksq=kx*kx+ky*ky+kz*kz
				if((ksq.lt.ksqmax).and.(ksq.ne.0)) then
					totk=totk+1
					deikr=deikx(kx)*deiky(ky)*deikz(kz)
     &					-eikx(kx)*eiky(ky)*eikz(kz)
					dsumi(totk)=sumi(totk)+q(il)*deikr
					vdn=vdn+factor*kvec(totk)*
     *					conjg(dsumi(totk))*dsumi(totk)
				endif
			end do
		end do
	end do

	vdn = vdn / box
	dek = vdn - vd
      end

	subroutine pivot (ener,attemptpivot,naccpivot,angle)
c	this subroutine creates a new conformation by generating three
c     random euler angles and than rotating the macroion around its center
      implicit none
	integer maxpart,part,nhelp,help,npart, ncharges	,maxk
	parameter(maxpart=300,maxk=1000)
	real pi,angle,ang,box,hbox,siga(maxpart),sigb(maxpart)
      real	sigc(maxpart), der, dek,vd,vdn,ranf
	parameter (pi=3.141459265359796) 
      real x(-8:maxpart), y(-8:maxpart), z(-8:maxpart)
	real xs(-8:maxpart), ys(-8:maxpart), zs(-8:maxpart)
	integer  i,attemptpivot,naccpivot
	logical check1,check2,alfa(maxpart),beta(maxpart)
	real alpha, bet, gamma, cosa, cosb, cosg, sena, senb, seng
	real dx, dy, dz, a, b, c, xm1,xm2,xm3,ym1,ym2,ym3,xx,yy,zz
	real delx, dely, delz,cas,zm1,zm2,zm3,rxij,ryij,rzij,rijsq
	real xept,deltaer,ener,q(maxpart),sig(maxpart),rab 
	complex sumi(maxk), dsumi(maxk)
      common /conf1/ x,y,z,npart,ncharges,a,b,c,alfa,beta
	common /pot1/  sig,siga,sigb,sigc, q
	common / move / sumi, vd, dsumi, vdn
	common /sys1/  box, hbox,rab

      attemptpivot=attemptpivot+1
c pick a macroion
c     the variable part defines which macroion is to be rotated
	cas=ranf()
	if (cas.le.0.5) then
		part=1
	else
		part=2
	end if

c nhelp: starting number of the index of the charges to be affected
	nhelp=npart+(part-1)*ncharges
	help=(part-1)*(-4)	 
	check1=part.eq.1
	check2=part.eq.2

c     old conformation storage
c     xs, ys,zs are the variables over which the move is made:
c     only if the move is accepted, the program adjorns x, y and z.
	do i=nhelp+1,nhelp+ncharges
		xs(i)=x(i)
		ys(i)=y(i)
		zs(i)=z(i)
	end do

	do i=help-4,help-1
		xs(i)=x(i)
		ys(i)=y(i)
		zs(i)=z(i)
	end do

c	maximum angular displacement, as computed in the subroutine 
c     adjustangle, 0<angle<1
	ang=angle*2*pi

c	Euler angles, alpha, bet, gamma	(randomly distributed around 0, 2pi)
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
	do i=help-4,help-2
			dx=xs(i)-xs(-1+help)
			dy=ys(i)-xs(-1+help)
			dz=zs(i)-xs(-1+help)

			delx=(cosg*cosb*cosa-seng*sena)*dx
			delx=delx+(cosg*cosb*sena+seng*cosa)*dy+(-cosg*senb)*dz
	
			dely=(-seng*cosb*cosa-cosg*sena)*dx+(seng*senb)*dz
			dely=dely+(-seng*cosb*sena+cosg*cosa)*dy
	
			delz=(senb*cosa)*dx+(senb*sena)*dy+cosb*dz
	
			xs(i)=xs(-1+help)+delx
			ys(i)=ys(-1+help)+dely
			zs(i)=zs(-1+help)+delz 
	end do

c     check over all microions:
	if (a.eq.b.and.b.eq.c) goto 100
	xm1=xs(-2+help)-x(-1+help)
	ym1=ys(-2+help)-y(-1+help)
	zm1=zs(-2+help)-z(-1+help)

	xm2=xs(-3+help)-x(-1+help)
	ym2=ys(-3+help)-y(-1+help)
	zm2=zs(-3+help)-z(-1+help)

	xm3=xs(-4+help)-x(-1+help)
	ym3=ys(-4+help)-y(-1+help)
	zm3=zs(-4+help)-z(-1+help)

	do i=1,npart
			if ((check1.and.alfa(i)).or.
     *			(check2.and.beta(i))) then

				xx=x(i)-xs(-1+help)
				yy=y(i)-ys(-1+help)
				zz=z(i)-zs(-1+help)
c     now the program is computing the distance of a microion from the macroion surface

				rxij = (xx*xm1 + yy*ym1 + zz*zm1)/a
				ryij = (xx*xm2 + yy*ym2 + zz*zm2)/b
				rzij = (xx*xm3 + yy*ym3 + zz*zm3)/c
c     periodic boundary conditions
				rxij = rxij - box * anint ( rxij / box )
				ryij = ryij - box * anint ( ryij / box )
				rzij = rzij - box * anint ( rzij / box )

				rijsq=rxij*rxij / siga(i) + ryij*ryij / 
     *				sigb(i) + rzij*rzij	/ sigc(i) -1.0e0

				if (rijsq .lt.0.e0 )  return
			end if
	end do
100	continue
c     rotating the charges in the macroion
		do i=nhelp+1,nhelp+ncharges
			dx=xs(i)-xs(-1+help)
			dy=ys(i)-xs(-1+help)
			dz=zs(i)-xs(-1+help)

			delx=(cosg*cosb*cosa-seng*sena)*dx
			delx=delx+(cosg*cosb*sena+seng*cosa)*dy+(-cosg*senb)*dz
		
			dely=(-seng*cosb*cosa-cosg*sena)*dx+(seng*senb)*dz
			dely=dely+(-seng*cosb*sena+cosg*cosa)*dy
		
			delz=(senb*cosa)*dx+(senb*sena)*dy+cosb*dz
		
			xs(i)=xs(-1+help)+delx
			ys(i)=ys(-1+help)+dely
			zs(i)=zs(-1+help)+delz 
		end do
		call deltaerbig( part, xs, ys, zs, der, dek ) 
		if (der.lt.299.e0) then
			deltaer=-der-dek
			xept=exp(deltaer)
			cas=ranf()
			if (cas.lt.xept) then
				vd=vdn
				do i=1,maxk
				sumi(i)=dsumi(i)
				end do
				ener=ener+der+dek
				do i=nhelp+1,nhelp+ncharges
				x(i)=xs(i)
				y(i)=ys(i)
				z(i)=zs(i)
				end do

				do i=help-4,help-2
				x(i)=xs(i)
				y(i)=ys(i)
				z(i)=zs(i)
				end do
				naccpivot=naccpivot+1
			else
			end if
		else
		endif
	end
	subroutine deltaerbig( part, xs, ys, zs, der, dek )

c     this procedure is to be used when a macroion is rotated
!*********************************************************************
!**	calculates energy change using ewald method                  **
!**	il		index of particle to be changed						 **
!**	xn		new x coordinate for particle il 					 **
!**	yn		new y coordinate for particle il 					 **
!**	zn		new x coordinate for particle il 					 **
!**	der		energy change in r-space after the move	    		 **
!**	dek		energy change in k-space after the move	    		 **
!**	dold	old distance between particle il and other particles **
!**	dnew	new distance between particle il and other particles **
!*********************************************************************
 
        implicit none
        integer     maxk, totk, kmax, kx, ky, kz, i, ksqmax, ksq, 
     &              part, npmax, npart,ncharges,nhelp,k

        parameter   ( maxk = 1000, kmax = 5, ksqmax = 27, npmax=300)
        real         twopi, twopl, factor, vd, vdn, rsqpi
        real       kvec(maxk), kappa, dek,der
        real    box, hbox ,siga(npmax), sigb(npmax), sigc(npmax)
        real    sig(npmax), q(npmax), a,b,c
        real       x(-8:npmax), y(-8:npmax), z(-8:npmax) 
	  real        xs(-8:npmax), ys(-8:npmax), zs(-8:npmax)
        real       e1, e2, rxij, ryij, rzij
	  real       rijsq, rij, krij
	  real       erfuc,rab
	  logical 	   alfa(npmax), beta(npmax)
        parameter   ( twopi = 6.2831853, rsqpi = 0.5641896 )

        complex eikx(0:kmax), eiky(-kmax:kmax), eikz(-kmax:kmax)
        complex deikx(0:kmax),deiky(-kmax:kmax),deikz(-kmax:kmax)
        complex deikr, sumi(maxk), dsumi(maxk)

        common / block2 / kvec, kappa, twopl
	  common / move /   sumi, vd, dsumi, vdn
        common / conf1 /  x, y, z, npart, ncharges, a,b,c, alfa, beta
        common / pot1 /   sig,siga,sigb,sigc,q
        common / sys1 /   box,hbox,rab

        nhelp=npart+(part-1)*ncharges
	 
c     the procedure is calculating the variation in the real part of the ewald sum
c     e1 corresponds to the old configuration, e2 to the new

		e1=0.
		e2=0.
c     energy at old position

c     now considering the interaction between macroions

		do  i=npart+1,npart+ncharges
			do k=npart+ncharges+1,npart+2*ncharges

				rxij = x(i) - x(k)
				ryij = y(i) - y(k)
				rzij = z(i) - z(k)

				rijsq=rxij*rxij + ryij*ryij + rzij*rzij
				rij=sqrt(rijsq)

				krij=kappa*rij/box
				e1=e1 + q(i)*q(k)*erfuc(krij)/rij
			end do
		end do

c     now considering the interactions with microions

		do i=1,npart
			do k=1,ncharges

				rxij = x(i) - x(nhelp+k)
				ryij = y(i) - y(nhelp+k)
				rzij = z(i) - z(nhelp+k)
c     periodic boundary conditions 
				rxij = rxij - box * anint ( rxij / box )
				ryij = ryij - box * anint ( ryij / box )
				rzij = rzij - box * anint ( rzij / box )

				rijsq=rxij*rxij + ryij*ryij + rzij*rzij
				rij=sqrt(rijsq)

				krij=kappa*rij/box
				e1=e1+q(i)*q(nhelp+k)*erfuc(krij)/rij

			end do

		end do
c     computation of electrostatic energy in the new position

c     now considering the interactions with microions

		do i=1,npart
			do k=1,ncharges

				rxij = x(i) - xs(nhelp+k)
				ryij = y(i) - ys(nhelp+k)
				rzij = z(i) - zs(nhelp+k)
c     periodic boundary conditions
				rxij = rxij - box * anint ( rxij / box )
				ryij = ryij - box * anint ( ryij / box )
				rzij = rzij - box * anint ( rzij / box )

				rijsq=rxij*rxij + ryij*ryij + rzij*rzij
				rij=sqrt(rijsq)

				krij=kappa*rij/box
				e2=e2+q(i)*q(nhelp+k)*erfuc(krij)/rij
			end do

		end do
c     now considering the interaction between macroions

		do  i=npart+1,npart+ncharges
			do k=npart+ncharges+1,npart+2*ncharges

				if (part.eq.1) then
					rxij = xs(i) - x(k)
					ryij = ys(i) - y(k)
					rzij = zs(i) - z(k)

				else if (part.eq.2) then
					rxij = x(i) - xs(k)
					ryij = y(i) - ys(k)
					rzij = z(i) - zs(k)

				end if

				rijsq=rxij*rxij + ryij*ryij + rzij*rzij
				rij=sqrt(rijsq)

				krij=kappa*rij/box
				e2=e2 + q(i)*q(k)*erfuc(krij)/rij
			end do
		end do


      der = e2 - e1



!    ** calculate energy change in k-space
!    ** construct dexp(ik.r) for il-th ion and k-vectors **
!    **   at new and old positions                      **
!    ** calculate kx, ky, kz = 0 , -1 and 1 dexplicitly  **

c     the following cycle is repeated over all charges in a macroion

c      initializing variables

       vdn = 0.0

	 do i=1,maxk
		dsumi(i)=sumi(i)
	 end do

       eikx(0)=(1.0,0.0)
       eiky(0)=(1.0,0.0)
       eikz(0)=(1.0,0.0)

       deikx(0)=(1.0,0.0)
       deiky(0)=(1.0,0.0)
       deikz(0)=(1.0,0.0)

c     considering all the moved charges 
 
	do i=1,ncharges
		eikx(1)=cmplx(cos(twopl*x(nhelp+i))
     *	 ,sin(twopl*x(nhelp+i)))
		eiky(1)=cmplx(cos(twopl*y(nhelp+i))
     *	 ,sin(twopl*y(nhelp+i)))
		eikz(1)=cmplx(cos(twopl*z(nhelp+i))
     *	 ,sin(twopl*z(nhelp+i)))


		 eiky(-1)=conjg(eiky(1))
		 eikz(-1)=conjg(eikz(1))


		deikx(1)=cmplx(cos(twopl*xs(nhelp+i))
     *    ,sin(twopl*x(nhelp+i)))
		deiky(1)=cmplx(cos(twopl*ys(nhelp+i))
     *    ,sin(twopl*y(nhelp+i)))
		deikz(1)=cmplx(cos(twopl*zs(nhelp+i))
     *    ,sin(twopl*z(nhelp+i)))


		 deiky(-1)=conjg(deiky(1))
		 deikz(-1)=conjg(deikz(1))


!    ** calculate remaining kx, ky and kz by recurrence **
!


		do kx=2,kmax
			deikx(kx)=deikx(kx-1)*deikx(1)
			eikx(kx)=eikx(kx-1)*eikx(1)
		enddo

		do ky=2,kmax
			deiky(ky)=deiky(ky-1)*deiky(1)
			deiky(-ky)=conjg(deiky(ky))
			eiky(ky)=eiky(ky-1)*eiky(1)
			eiky(-ky)=conjg(eiky(ky))
		enddo

		do kz=2,kmax
			deikz(kz)=deikz(kz-1)*deikz(1)
			deikz(-kz)=conjg(deikz(kz))
			eikz(kz)=eikz(kz-1)*eikz(1)
			eikz(-kz)=conjg(eikz(kz))
		enddo

!
!    ** sum over all vectors **


c     computation of the variation in density

		totk=0
		do kx=0,kmax
			if(kx.eq.0) then
				factor=1.0
			else
				factor=2.0
			endif
			do ky=-kmax,kmax
				do  kz=-kmax,kmax
					ksq=kx*kx+ky*ky+kz*kz
					if((ksq.lt.ksqmax).and.(ksq.ne.0)) then
						totk=totk+1
						deikr=deikx(kx)*deiky(ky)*deikz(kz)
     *					    -eikx(kx)*eiky(ky)*eikz(kz)
						dsumi(totk)=dsumi(totk)+q(nhelp+i)*deikr
					endif
				end do
			end do
		end do
	end do
c     actual computation of the variation in fourier part of the electrostatic energy
	totk=0
	do kx=0,kmax
		if(kx.eq.0) then
				factor=1.0
		else
				factor=2.0
		endif
		do ky=-kmax,kmax
				do  kz=-kmax,kmax
					ksq=kx*kx+ky*ky+kz*kz
					if((ksq.lt.ksqmax).and.(ksq.ne.0)) then
						totk=totk+1
						vdn=vdn+factor*kvec(totk)*
     *					conjg(dsumi(totk))*dsumi(totk)
					endif
				end do
		end do
	end do
	vdn = vdn / box
	dek = vdn - vd
 	end

      subroutine sample(myid,switch, ener,dr,angle,deltar)
!
!   switch (input) = 0: initialize varibales 
!                  = 1: sample averages
!                  = 2: return average results
!
! ener (input) : total energy/kt
! enert        : average total energy/nkt
! ns           : number of sample points
! f            : force between macroions

      implicit none
      integer  npmax, ns, switch, npart, ncharges,nhcc,myid
      parameter (npmax=300)
      integer  i, nsf,ns0
      real     ener, enert, f, mf, fhs, fcc, mfhs, mfcc
      real     enert1, mf1, mfhs1, mfcc1,phcc,deltar
      real     siga(npmax),sigb(npmax),sigc(npmax),sig(npmax)
      real     box, hbox, q(npmax)
      real     x(-8:npmax), y(-8:npmax), z(-8:npmax)
      real     dr,a,b,c,rab,angle,rab0(20) 

	logical alfa(npmax), beta(npmax)
      real Enert0, MF0, MFHS0, MFCC0,phcc0
      common /rest1/ ns0, Enert0, MF0, MFHS0, MFCC0,phcc0,rab0
      common /sys1/  box, hbox,rab
      common /pot1/  sig, siga,sigb,sigc, q
      common /conf1/ x, y, z, npart,ncharges,a,b,c,alfa,beta

      character filename(20)*5,filename2(20)*5
      character output1*5,iniconf*5
      data filename/'out1','out2','out3','out4','out5','out6'
     &      ,'out7','out8','out9','out10','out11','out12'
     &      ,'out13','out14','out15','out16','out17','out18'
     &      ,'out19','out20'/
      data filename2/'inp1','inp2','inp3','inp4','inp5','inp6'
     &      ,'inp7','inp8','inp9','inp10','inp11','inp12'
     &      ,'inp13','inp14','inp15','inp16','inp17','inp18'
     &      ,'inp19','inp20'/

      save ns,enert, mf, mfcc,mfhs,nhcc,phcc 
        iniconf=filename2(myid+1)
        output1=filename(myid+1) 

!       ---initialize
      if (switch.eq.0) then

		ns=ns0*10**5
		nhcc=INT(float(ns*4)*phcc0 )
c        ---Total energy
		Enert=float(ns)*Enert0
c        ---Mean force between macroions
		MF = float(ns*4)*MF0
		MFHS= float(ns*4)*MFHS0
		MFCC= float(ns*4)*MFCC0

	else if (switch.eq.1) then
			ns=ns+1
! ---sample excess energy 
			enert = enert + ener

! ---sample mean force between macroions
c     with this instructions, the force is sampled 4 times at the given conditions
				call force(-1,deltar,nhcc,f, fhs, fcc)
				mf = mf + f
				mfhs=mfhs+fhs
				mfcc=mfcc+fcc
				call force(-1,-deltar,nhcc,f, fhs, fcc)
				mf = mf + f
				mfhs=mfhs+fhs
				mfcc=mfcc+fcc
				call force(-5,deltar,nhcc,f, fhs, fcc)
				mf = mf + f
				mfhs=mfhs+fhs
				mfcc=mfcc+fcc
				call force(-5,-deltar,nhcc,f, fhs, fcc)
				mf = mf + f
				mfhs=mfhs+fhs
				mfcc=mfcc+fcc
               if (mod(ns,10**5).eq.0) then
                enert1 = enert/ns
			  nsf=4*ns
	          phcc=float(nhcc)/float(nsf)
                mf1    = mf/nsf
                mfhs1  = mfhs/nsf
                mfcc1  = mfcc/nsf
			  open (69,file=output1,position='append')
				write(69,223) rab0(myid+1),ns/10**5,enert1,
     &				mf1,mfhs1,mfcc1,phcc,deltar
			  close (69)
  
			  open(64,file=iniconf)
   			  write (64, *) npart, dr, angle
			  write (64,*)  ns/10**5, Enert1, MF1, MFHS1,MFCC1,phcc

			  do i = 1, npart
			  write (64,*) alfa(i),beta(i),x(i),y(i),z(i)
			  end do 
			  do i=-1,-8,-1
      			write (64, *)  x(i), y(i), z(i)
			  end do
			  do i=npart+1,npart+2*ncharges
       			write (64, *)  q(i), x(i), y(i), z(i)
			  end do
      		  close(64)
		endif


      else if (switch.eq.2) then
	    nsf=4*ns
		open( 66, file='results',position='append')
		enert = enert/ns
	    phcc=float(nhcc)/float(nsf)
		mf    = mf/nsf
		mfhs  = mfhs/nsf
		mfcc  = mfcc/nsf
	write(66,223) rab0(myid+1),ns/10**5,enert,mf,mfhs,mfcc,phcc,deltar
		close(66)
      end if

223   format(1x,f5.2,i5,1x,4f10.2,1x,f8.5,1x,f6.4)   
      end

	  subroutine force(il,deltar,nhcc,f,fhs,fcc)
*********************************************************************
**	calculates energy change using ewald method                  **
**	xn		new x coordinate for particle il 					 **
**	yn		new y coordinate for particle il 					 **
**	zn		new x coordinate for particle il 					 **
**	der		energy change in r-space after the move	    		 **
**	dek		energy change in k-space after the move	    		 **
**	f	    force between macroions                              **
*********************************************************************
        implicit none
        integer     maxk, totk, kmax, kx, ky, kz, i, ksqmax, ksq,k, 
     &              npmax, npart, ncoll,il,nhelp,ncharges,nhcc
	  
        parameter   ( maxk = 1000, kmax = 5, ksqmax = 27, npmax=300)
        real   twopi, twopl, factor, vd, vdn,kvec(maxk), 
     *     kappa, dek,der, box, hbox, rab, xx, deltar, f, fcc,yy,zz, 
     *     fhs, siga(npmax),sigb(npmax),sigc(npmax), q(npmax), 
     *     x(-8:npmax), y(-8:npmax),xm3,ym3,zm3,a,b,c, 
     *     z(-8:npmax), xn(-8:npmax), yn(-8:npmax), zn(-8:npmax), e1,
     *     e2, rxij, ryij,xm1,ym1,zm1,xm2,ym2,zm2,sig(npmax), 
     *     rzij,rijsq, rij, krij,erfuc,delx,dely,delz

	  logical check1,check2,alfa(npmax),beta(npmax)

        parameter   ( twopi = 6.2831853)

        complex eikx(0:kmax), eiky(-kmax:kmax), eikz(-kmax:kmax)
        complex deikx(0:kmax),deiky(-kmax:kmax),deikz(-kmax:kmax)
        complex deikr, sumi(maxk), dsumi(maxk)

        common / block2 / kvec, kappa,twopl
	common / move /   sumi, vd, dsumi, vdn
        common / conf1 /  x, y, z, npart,ncharges,a,b,c, alfa, beta
        common / pot1 /   sig,siga,sigb,sigc,q
        common / sys1 /   box,hbox,rab

c     --- make a small displacement for one of the macroions 
c     the variable il define which macroion is moved, while ia denotes
c     the one which is fixed

c     computes new coordinates for the chosen macroion
	check1=(il.eq.(-1))
	check2=(il.eq.(-5))
	nhelp=npart+ncharges*(-il-1)/4
	delx=deltar/1.732050808
	dely=deltar/1.732050808
	delz=deltar/1.732050808

c     the macroion is translated along the main diagonal
c     so every coordinated has to be shifted of the same length
	do i=-3+il,il	
c             geometry: position of the main radii of inertia
			xn(i)  = x(i) + delx
			yn(i)  = y(i) + dely 
			zn(i)  = z(i) + delz
	end do
	do i=nhelp+1,nhelp+ncharges
c             charges: position of every charge in the macroion
			xn(i)  = x(i) + delx 
			yn(i)  = y(i) + dely 
			zn(i)  = z(i) + delz
	end do

c     energy change in real space: e1 is the old energy, e2 the new
		e1=0.
		e2=0.
c     particle at old position
		do i=1,npart
				do k=nhelp+1,nhelp+ncharges
					rxij = x(i) - x(k)
					ryij = y(i) - y(k)
					rzij = z(i) - z(k)
c                     periodic boundary conditions
					rxij = rxij - box * anint ( rxij / box )
					ryij = ryij - box * anint ( ryij / box )
					rzij = rzij - box * anint ( rzij / box )

					rijsq=rxij*rxij + ryij*ryij + rzij*rzij
					rij=sqrt(rijsq)

					krij=kappa*rij/box
					e1=e1 + q(i)*q(k)*erfuc(krij)/rij
				end do
		end do
c     interactions between macroions
	    do i=1,ncharges
			do k=1,ncharges
			rxij = x(npart+i) - x(npart+ncharges+k)
			ryij = y(npart+i) - y(npart+ncharges+k)
			rzij = z(npart+i) - z(npart+ncharges+k)
	        rijsq= rxij*rxij + ryij*ryij + rzij*rzij
			rij=sqrt(rijsq)

			krij=kappa*rij/box
	        e1=e1 + q(npart+i)*q(npart+ncharges+k)*erfuc(krij)/rij
			end do
		end do

c     particle at new position

c     initializing variables to count the number of collisions following from
c     a random displacement of a macroion
		ncoll=0
c     the next are the three position vectors which describe
c     the position in space of a macroion
		xm1=xn(-1+il)-xn(il)
		ym1=yn(-1+il)-yn(il)
		zm1=zn(-1+il)-zn(il)

		xm2=xn(-2+il)-xn(il)
		ym2=yn(-2+il)-yn(il)
		zm2=zn(-2+il)-zn(il)

		xm3=xn(-3+il)-xn(il)
		ym3=yn(-3+il)-yn(il)
		zm3=zn(-3+il)-zn(il)

c     the procedure checks for collisions only if the microion is labeled 
c     within the macrosphere defined in the beginning procedure
		do i=1,npart	
c			if ((check1.and.alfa(i)).or.
c     *			(check2.and.beta(i))) then
					xx=x(i)-xn(il)
					yy=y(i)-yn(il)
					zz=z(i)-zn(il)
c     now the program is computing the distance of a microion from the macroion surface
c	rx(y,z)ij are three coordinates of the micorion in the macroion coordinate system
					rxij = (xx*xm1 + yy*ym1 + zz*zm1)/a
					ryij = (xx*xm2 + yy*ym2 + zz*zm2)/b
					rzij = (xx*xm3 + yy*ym3 + zz*zm3)/c
c     siga,sigb,sigc are squares of collision radii
					rijsq=(rxij*rxij / siga(i) + ryij*ryij / 
     *					sigb(i) + rzij*rzij	/ sigc(i)) -1.e0
					if (rijsq .lt. 0.e0 )  then
						ncoll=ncoll+1
					end if
c			end if
		end do
 
c     electrostatic energy due to the presence of microions around a macroion 
c     in the new position

		do i=1,npart
				do k=nhelp+1,nhelp+ncharges
					rxij = x(i) - xn(k)
					ryij = y(i) - yn(k)
					rzij = z(i) - zn(k)
					rxij = rxij - box * anint ( rxij / box )
					ryij = ryij - box * anint ( ryij / box )
					rzij = rzij - box * anint ( rzij / box )

					rijsq=rxij*rxij + ryij*ryij + rzij*rzij
					rij=sqrt(rijsq)

					krij=kappa*rij/box
					e2=e2 + q(i)*q(k)*erfuc(krij)/rij
				end do
		end do
c     computation of the electrostatic energy due to the interaction between macroions

 	do i=1,ncharges
		do k=1,ncharges
			if (check1) then
c                 the first macroion is moved
				rxij = xn(npart+i) - x(npart+ncharges+k)
				ryij = yn(npart+i) - y(npart+ncharges+k)
				rzij = zn(npart+i) - z(npart+ncharges+k)
			end if
			if (check2) then
c                 the second macroion is moved
				rxij = x(npart+i) - xn(npart+ncharges+k)
				ryij = y(npart+i) - yn(npart+ncharges+k)
				rzij = z(npart+i) - zn(npart+ncharges+k)
			end if

			rijsq=rxij*rxij + ryij*ryij + rzij*rzij
			rij=sqrt(rijsq)

			krij=kappa*rij/box
			e2=e2 + q(npart+i)*q(npart+ncharges+k)*erfuc(krij)/rij
		end do
      end do

      der = e2 - e1

c     now for each charge in the macroion surface the program is computing
c     the fourier part of the ewald sum for both the old and the new positions
c     initializing variable
	twopl=twopi/box
	do i=1,maxk
		dsumi(i)=sumi(i)
	end do
	vdn = 0.0
c     start cycle over all charges
	do i=1,ncharges
				eikx(0)=(1.0,0.0)
				eiky(0)=(1.0,0.0)
				eikz(0)=(1.0,0.0)

				deikx(0)=(1.0,0.0)
				deiky(0)=(1.0,0.0)
				deikz(0)=(1.0,0.0)

				eikx(1)=
     *				dcmplx(cos(twopl*x(nhelp+i)),
     *				sin(twopl*x(nhelp+i)))
				eiky(1)=
     *				dcmplx(cos(twopl*y(nhelp+i)),
     *				sin(twopl*y(nhelp+i)))
				eikz(1)=
     *				dcmplx(cos(twopl*z(nhelp+i)),
     *				sin(twopl*z(nhelp+i)))

				deikx(1)=
     *				dcmplx(cos(twopl*xn(nhelp+i)),
     *				sin(twopl*xn(nhelp+i)))
				deiky(1)=
     *				dcmplx(cos(twopl*yn(nhelp+i)),
     *				sin(twopl*yn(nhelp+i)))
				deikz(1)=
     *				dcmplx(cos(twopl*zn(nhelp+i)),
     *				sin(twopl*zn(nhelp+i)))

				eiky(-1)=conjg(eiky(1))
				eikz(-1)=conjg(eikz(1))

				deiky(-1)=conjg(deiky(1))
				deikz(-1)=conjg(deikz(1))

c    				* calculate remaining kx, ky and kz by recurrence **
c
				do kx=2,kmax
					deikx(kx)=deikx(kx-1)*deikx(1)
					eikx(kx)=eikx(kx-1)*eikx(1)
				enddo

				do ky=2,kmax
					deiky(ky)=deiky(ky-1)*deiky(1)
					deiky(-ky)=conjg(deiky(ky))
					eiky(ky)=eiky(ky-1)*eiky(1)
					eiky(-ky)=conjg(eiky(ky))
				enddo

				do kz=2,kmax
					deikz(kz)=deikz(kz-1)*deikz(1)
					deikz(-kz)=conjg(deikz(kz))
					eikz(kz)=eikz(kz-1)*eikz(1)
					eikz(-kz)=conjg(eikz(kz))
				enddo

c    				* sum over all vectors **

c             computation of the variation in density:
			
			totk=0
			do kx=0,kmax
				if(kx.eq.0) then
				 factor=1.0
				else
				 factor=2.0
				endif
				do  ky=-kmax,kmax
					  do  kz=-kmax,kmax
						   ksq=kx*kx+ky*ky+kz*kz
						   if((ksq.lt.ksqmax).and.(ksq.ne.0)) then
								totk=totk+1
								deikr=deikx(kx)*deiky(ky)*deikz(kz)
     *							-eikx(kx)*eiky(ky)*eikz(kz)
								dsumi(totk)=dsumi(totk)+q(nhelp+i)
     *							*deikr
						   endif
					  end do
				end do
			end do
	end do

c             computation of the effective variation in energy

	totk=0
	do kx=0,kmax
		if(kx.eq.0) then
				 factor=1.0
		else
				 factor=2.0
		endif
		do  ky=-kmax,kmax
				do  kz=-kmax,kmax
						   ksq=kx*kx+ky*ky+kz*kz
						   if((ksq.lt.ksqmax).and.(ksq.ne.0)) then
								totk=totk+1
								vdn=vdn+factor*kvec(totk)*conjg(dsumi
     *							(totk))*dsumi(totk)
						   endif
				end do
		end do
	end do

c     concluding:

      nhcc=nhcc+ncoll
	vdn = vdn / box
	dek = vdn - vd
	if(check1) then
      fcc = -(dek+der)/deltar   
      fhs = -2.0e0*float(ncoll)/deltar
	else
      fcc = (dek+der)/deltar   
      fhs = 2.0e0*float(ncoll)/deltar
	endif
      f = fcc + fhs
      end

      subroutine adjust(attempdisp,naccdisp,dr)
c     this subroutine adjusts the maximum displacement of a microion so that
c     the acceptance ratio is around 50%
	implicit none
	integer attempdisp, naccdisp, attempp, naccp
	real box, hbox, dro, frac, dr,rab
	save naccp, attempp
	common /sys1/ box, hbox, rab

c     the variables are initialized at the beginning of the program (left case) 
c     or when it is restarted (right case)
	if ((attempdisp.eq.0).or.(attempp.ge.attempdisp)) then
		naccp=naccdisp
		attempp=attempdisp
	else
		frac=float(naccdisp-naccp)/ float(attempdisp-attempp)
		dro=dr
		dr=dr*abs(frac/0.50)
		if (dr/dro.gt.1.50) dr=dro*1.5
		if (dr/dro.lt.0.50) dr=dro*0.5 
		if (dr.gt.hbox/2.0) dr=hbox/2.0
		naccp=naccdisp
		attempp=attempdisp
	end if
1     format('max. displ. set to:',f6.3,'(old:',f6.3,')',/,
     *   'frac. acc.:',f4.2,'attempts:',i7,'success:',i7)

	end

	
      subroutine adjustangle(attemppivot,naccpivot,angle)
c     this subroutine adjusts the maximum rotation of a macroion so that
c     the acceptance ratio is around 50%
	implicit none
	integer attemppivot, naccpivot, attemppiv, naccpiv
	real angleo, frac, angle

	save naccpiv, attemppiv

c     the variables are initialized at the beginning of the program (left case) 
c     or when it is restarted (right case)

	if ((attemppivot.eq.0).or.(attemppiv.ge.attemppivot)) then
		naccpiv=naccpivot
		attemppiv=attemppivot
	else
		frac=float(naccpivot-naccpiv)/ float(attemppivot-attemppiv)
		angleo=angle
		angle=angle*abs(frac/0.50)
		if (angle/angleo.gt.1.50) angle=angleo*1.5
		if (angle/angleo.lt.0.50) angle=angleo*0.5 
		write (*,2) angle,angleo,frac,
     *		attemppivot-attemppiv,naccpivot-naccpiv
		naccpiv=naccpivot
		attemppiv=attemppivot
	end if
2     format('max. angle displ. set to:',f6.3,'(old:',f6.3,')',/,
     *   'frac. acc.:',f4.2,'attempts:',i7,'success:',i7)

	end
     
	subroutine dodeca(A,sigm)
c Assign unit charges on the vertexes of a regular dodecahedron
c R1: separation between pentagon center and any vertex
c edge2: separation between next nearest neighbors in a pentagon
c h1: height of the first middle pentagon
c Circumdiameter of dodecahedron
	implicit none
	integer i,j
	real pi, edge, R1, R2, h1, A(4,5,3),edge2,sigm
	pi=4.e0*atan(1.e0)
c Circumdiameter of a dodecahedron = 2.802517e0 * edge
	edge=sigm/2.802517e0
	R1=sqrt(edge**2/(2.e0*(1.e0-cos(72.e0/180.e0*pi))))
	edge2=2.e0*edge*sin(54.e0/180.e0*pi)
	R2=sqrt(edge2**2/(2.e0*(1.e0-cos(72.e0/180.e0*pi))))

c Dihedral angle of dodecahedron is 116.565051o, Inradius=1.113516*edge
	h1=edge*sin(72.e0/180.e0*pi)*sin(116.565051e0/180.e0*pi)

c Bottom pentagon
	A(1,1,1)=0.0e0
	A(1,1,2)=R1

	A(1,2,1)=-R1*cos(18.e0/180.e0*pi)
	A(1,5,1)= R1*cos(18.e0/180.e0*pi)

	A(1,2,2)= R1*sin(18.e0/180.e0*pi)
	A(1,5,2)= R1*sin(18.e0/180.e0*pi)

	A(1,3,1)=-R1*sin(36.e0/180.e0*pi)
	A(1,4,1)= R1*sin(36.e0/180.e0*pi)

	A(1,3,2)=-R1*cos(36.e0/180.e0*pi)
	A(1,4,2)=-R1*cos(36.e0/180.e0*pi)


c Pentagon in the middle
	A(2,1,1)=0.0e0
	A(2,1,2)=R2

	A(2,2,1)=-R2*cos(18.e0/180.e0*pi)
	A(2,5,1)= R2*cos(18.e0/180.e0*pi)

	A(2,2,2)= R2*sin(18.e0/180.e0*pi)
	A(2,5,2)= R2*sin(18.e0/180.e0*pi)

	A(2,3,1)=-R2*sin(36.e0/180.e0*pi)
	A(2,4,1)= R2*sin(36.e0/180.e0*pi)

	A(2,3,2)=-R2*cos(36.e0/180.e0*pi)
	A(2,4,2)=-R2*cos(36.e0/180.e0*pi)

c Next pentagon in the middle
	A(3,1,1)=0.0e0
	A(3,1,2)=-R2

	A(3,2,1)=-R2*cos(18.e0/180.e0*pi)
	A(3,5,1)= R2*cos(18.e0/180.e0*pi)

	A(3,2,2)=-R2*sin(18.e0/180.e0*pi)
	A(3,5,2)=-R2*sin(18.e0/180.e0*pi)

	A(3,3,1)=-R2*sin(36.e0/180.e0*pi)
	A(3,4,1)= R2*sin(36.e0/180.e0*pi)

	A(3,3,2)= R2*cos(36.e0/180.e0*pi)
	A(3,4,2)= R2*cos(36.e0/180.e0*pi)

c Top pentagon
	A(4,1,1)=0.0e0
	A(4,1,2)=-R1

	A(4,2,1)=-R1*cos(18.e0/180.e0*pi)
	A(4,5,1)= R1*cos(18.e0/180.e0*pi)

	A(4,2,2)=-R1*sin(18.e0/180.e0*pi)
	A(4,5,2)=-R1*sin(18.e0/180.e0*pi)

	A(4,3,1)=-R1*sin(36.e0/180.e0*pi)
	A(4,4,1)= R1*sin(36.e0/180.e0*pi)

	A(4,3,2)= R1*cos(36.e0/180.e0*pi)
	A(4,4,2)= R1*cos(36.e0/180.e0*pi)

c assign z coordinates
	do i=1,5
	A(1,i,3)=-1.113516*edge
	A(2,i,3)=-1.113516*edge+h1
	A(3,i,3)= 1.113516*edge-h1 
	A(4,i,3)= 1.113516*edge
	enddo
	end

      function erfuc ( x )

!   *******************************************************************
!   ** approximation to the complementary error function             **
!   **                                                               **
!   ** reference:                                                    **
!   **                                                               **
!   ** abramowitz and stegun, handbook of mathematical functions,    **
!   **    national bureau of standards, formula 7.1.26               **
!   *******************************************************************
        implicit none

        real       a1, a2, a3, a4, a5, p

        parameter ( a1 = 0.254829592, a2 = -0.284496736 )
        parameter ( a3 = 1.421413741, a4 = -1.453152027 )
        parameter ( a5 = 1.061405429, p  =  0.3275911   )

        real        t, xsq, tp
		real        erfuc, x

!   *******************************************************************
        t  = 1.0 / ( 1.0 + p * x )
        xsq = x * x
        tp = t * ( a1 + t * ( a2 + t * ( a3 + t * ( a4 + t * a5 ) ) ) )
        erfuc = (tp * exp ( -xsq ) )


      return
      end

      subroutine toterg(gw, ener)
!*******************************************************************
!
! ener (output) : total energy/kt
! gw: gaussian width parameter, adjust for efficiency of ewald sum
!
!*******************************************************************

      implicit none
      integer npmax,npart,ncharges
      parameter (npmax=300)

      real     x(-8:npmax),y(-8:npmax),z(-8:npmax)
      real     box,hbox,rab,sig(1:npmax),siga(1:npmax)
      real     q(1:npmax),sigb(1:npmax),sigc(1:npmax) 
      real     vr, vk, ener, gw,a,b,c

      logical alfa(1:npmax),beta(1:npmax)

      common /sys1/  box,hbox,rab
      common /pot1/  sig,siga,sigb,sigc,q
      common /conf1/ x, y, z, npart,ncharges,a,b,c,alfa,beta

      call setup(gw) 
      call rwald ( npart+2*ncharges, vr, x, y, z, q, box)
      call kwald ( npart+2*ncharges, vk, x, y, z, q, box )

      ener = vr + vk
      end

!************************************************************************************************************************
!   *******************************************************************
!   ** real-space and reciprocal-space parts of ewald sum for ions.  **
!   **                                                               **
!   ** references:                                                   **
!   **                                                               **
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
!   ** real function erfuc ( x )                                      **
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
!   ** potential energy at each configuration. the self term is      **
!   ** subtracted from the k-space contribution in kwald.            **
!   ** the surface term for simulations in vacuum is not included.   **
!   ** routine rwald returns the r-space contribution to the ewald   **
!   ** sum and is called for each configuration in the simulation.   **
!   ** a cubic box and unit box length are assumed throughout.       **
!   *******************************************************************



        subroutine setup(gw) 

!   *******************************************************************
!   ** routine to set up the wave-vectors for the ewald sum.         **
!   ** gw: gaussian width parameter                                  **
!   ** in this example we allow a maximum of 1000 wavevectors.       **
!   *******************************************************************
        implicit none
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
!   ** function erfuc ( x )                          **
!   **    returns the complementary error function                   **
!   *******************************************************************
        implicit none
        integer     maxk, npmax
        parameter ( maxk = 1000, npmax=300 )
        real      kvec(maxk), kappa, twopl

        integer              n, i, j
        real      rx(-8:npmax), ry(-8:npmax), rz(-8:npmax), q(1:npmax)
        real         vr, box
        real         rxi, ryi, rzi, zi, rxij, ryij, rzij
        real         rijsq, rij, krij, vij
    	real erfuc

        common / block2 / kvec, kappa, twopl
!   *******************************************************************

        vr = 0.0
        do 100 i = 1, n - 1

           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           zi  = q(i)

           do 99 j = i + 1, n

              rxij = rxi - rx(j)
              ryij = ryi - ry(j)
              rzij = rzi - rz(j)

              rxij = rxij - box * anint ( rxij /box )
              ryij = ryij - box * anint ( ryij /box )
              rzij = rzij - box * anint ( rzij /box )

              rijsq = rxij * rxij + ryij * ryij + rzij * rzij
              rij   = sqrt ( rijsq )
              krij  = kappa * rij	/ box
              vij   = zi * q(j) * erfuc( krij ) / rij

              vr    = vr + vij


99         continue

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
        implicit none
        integer  maxk, n, totk, kmax, kx, ky, kz, i, ksqmax, ksq, npmax
        parameter   ( maxk = 1000, kmax = 5, ksqmax = 27, npmax=5000 )
        real       kvec(maxk)
        real      rx(-8:npmax), ry(-8:npmax), rz(-8:npmax), q(1:npmax)
        real       kappa, vk, box
        real         twopi, twopl, factor, vd, vs, rsqpi, vdn
        parameter   ( twopi = 6.2831853, rsqpi = 0.5641896 )

        complex       eikx(1:npmax, 0:kmax)
        complex       eiky(1:npmax, -kmax:kmax)
        complex       eikz(1:npmax, -kmax:kmax)
        complex       eikr(npmax), sumi(maxk), dsumi(maxk)

        common / block2 / kvec, kappa, twopl
	  common / move /   sumi, vd, dsumi, vdn
!   *******************************************************************

!   ** construct dexp(ik.r) for all ions and k-vectors **

!   ** calculate kx, ky, kz = 0 , -1 and 1 dexplicitly **


        do 10 i = 1, n

           eikx(i, 0) = (1.0, 0.0)
           eiky(i, 0) = (1.0, 0.0)
           eikz(i, 0) = (1.0, 0.0)

           eikx(i, 1) = cmplx ( cos ( twopl * rx(i) ) ,
     :                           sin ( twopl * rx(i) ) )
           eiky(i, 1) = cmplx ( cos ( twopl * ry(i) ) ,
     :                           sin ( twopl * ry(i) ) )
           eikz(i, 1) = cmplx ( cos ( twopl * rz(i) ) ,
     :                           sin ( twopl * rz(i) ) )

           eiky(i, -1) = conjg ( eiky(i, 1) )
           eikz(i, -1) = conjg ( eikz(i, 1) )

10      continue

!   ** calculate remaining kx, ky and kz by recurrence **

        do 12 kx = 2, kmax

           do 11 i = 1, n

              eikx(i, kx) = eikx(i, kx-1) * eikx(i, 1)

11         continue

12      continue

        do 14 ky = 2, kmax

           do 13 i = 1, n

              eiky(i,  ky) = eiky(i, ky-1) * eiky(i, 1)
              eiky(i, -ky) = conjg ( eiky(i, ky) )

13         continue

14      continue

        do 16 kz = 2, kmax

           do 15 i = 1, n

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

                    do 21 i = 1, n

                       eikr(i) = eikx(i, kx) * eiky(i, ky) * eikz(i, kz)
                       sumi(totk) = sumi(totk) + q(i) * eikr(i)

21                  continue

                    vd = vd + factor * kvec(totk) 
     &	         		  * conjg ( sumi(totk) ) * sumi(totk)

                 endif

22            continue

23         continue

24      continue

	 vd=vd / box

!   ** calculates self part of k-space sum **

        vs = 0.0

        do 25 i = 1, n

           vs = vs + q(i) * q(i)

25      continue

        vs = rsqpi * kappa * vs / box

!   ** calculate the total k-space potential **

        vk = vd - vs

        return
        end





