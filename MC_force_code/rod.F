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




      subroutine adjust(attemp, nacc, dr)
c
c     adjusts maximum displacement such that 50% of the
c     movels will be accepted
c
c  attemp (input) number of attemps that have been performed to displace a particle
c  nacc   (input) number of successful attemps to displace a particle
c  dr     (output) new maximum displacement
c
 
      implicit none
      integer attemp, nacc, attempp, naccp, zdna
      real     box,temp,hbox
      real     dro, frac, dr, rab
      save naccp, attempp
      common /sys1/ box,hbox,temp,rab,zdna
 
      if (attemp.eq.0.or.attempp.ge.attemp) then
         naccp = nacc
         attempp = attemp
      else
         frac = float(nacc-naccp)/float(attemp-attempp)
         dro = dr
         dr = dr*abs(frac/0.5e0)
c        ---limit the change:
         if (dr/dro.gt.1.5e0) dr = dro*1.5e0
         if (dr/dro.lt.0.5e0) dr = dro*0.5e0
         if (dr.gt.hbox/2.e0) dr = hbox/2.e0
         write (*, 99001) dr, dro, frac, attemp - attempp, nacc - naccp
c        ---store nacc and attemp for next use
         naccp = nacc
         attempp = attemp
      end if
      return
99001 format (' max. displ. set to : ', f6.3, ' (old : ', f6.3, ')', /, 
     &        ' frac. acc.: ', f4.2, ' attempts: ', i7, ' succes: ', i7)
      end


	  subroutine deltae( il, xn, yn, zn, der )
*********************************************************************
**	calculates energy change using ewald method                  **
**	il		index of particle to be changed						 **
**	xn		new x coordinate for particle il 					 **
**	yn		new y coordinate for particle il 					 **
**	zn		new x coordinate for particle il 					 **
**	der		energy change in r-space after the move	    		 **
*********************************************************************
        implicit none
        integer     i,il, npmax, npart ,zdna 
        parameter   ( npmax=500)
        real       twopi, rsqpi
        real       der
        real       box, temp, hbox
        real       sig(npmax), q(npmax), pi, dcol(npmax)
        real       x(npmax), y(npmax), z(npmax), xn, yn, zn
        real       e1, e2, xil, yil, zil, rxij, ryij, rzij
	  real       rijsq, rij
	  real       rab,a1,a2
        parameter   ( twopi = 6.2831853, rsqpi = 0.5641896 )

        common / conf1 /  x, y, z, npart
        common / pot1 /   pi,sig,q,dcol
        common / sys1 /   box,hbox,temp,rab,zdna
	  common /distan/   xil,yil,zil

	  a1=(xn-rab)**2+(yn-rab)**2
	  a2=(xn+rab)**2+(yn+rab)**2
	  if((a1.lt.dcol(il)).or.(a2.lt.dcol(il))) then
	  der=300.e0
	  return
	  else
	  endif

		e1=0.
		e2=0.

		xil=x(il)
		yil=y(il)
		zil=z(il)

		do 5 i=1,npart+2*zdna
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
	    e1=e1 + q(i)*q(il)/rij
		endif
5		continue

		do 6 i=1,npart
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
		if (rij .lt. (sig(il)+sig(i))/2. ) then
		der = 300.e0
		return
		else
		endif
	    e2=e2 + q(i)*q(il)/rij
		endif
6		continue

          do 7 i=npart+1,2*zdna+npart
          rxij = x(i) - xn
          ryij = y(i) - yn
          rzij = z(i) - zn
          rxij = rxij - box * anint ( rxij / box )
          ryij = ryij - box * anint ( ryij / box )
          rzij = rzij - box * anint ( rzij / box )
          rijsq= rxij*rxij + ryij*ryij + rzij*rzij
          rij  = sqrt(rijsq)
          e2   = e2 + q(i)*q(il)/rij
7         continue
      der = e2 - e1
      end



        subroutine rwald ( n, vr ,rx, ry, rz, q, box)

        implicit none
        integer     npmax
        parameter ( npmax=500 )

        integer      n, i, j
        real         rx(npmax), ry(npmax), rz(npmax), q(npmax)
        real         vr, box
        real         rxi, ryi, rzi, zi, rxij, ryij, rzij
        real         rijsq, rij, vij

c    *******************************************************************

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
              vij   = zi * q(j) / rij
              vr    = vr + vij
99         continue
100     continue
        end




	 subroutine force(deltar,idna,ncoll,f, fhs, fcc)
*********************************************************************
**	calculates energy change using ewald method              **
**	idna	index of the particle to be changed		 **
**	xn		new x coordinate for particle il         **
**	yn		new y coordinate for particle il 	 **
**	zn		new x coordinate for particle il 	 **
**	der		energy change in r-space after the move	 **
**	f	    force between macroions                      **
*********************************************************************
        implicit none
        integer     i, ii,jj,npmax, npart,  ncoll ,il,ia,zdna,idna  
        parameter   (npmax=500)
        real       twopi, rsqpi, der
	real       box, temp, hbox, rab, deltar, f, fcc, fhs
	real       sig(npmax), q(npmax), pi,dcol(npmax)
        real       x(npmax), y(npmax), z(npmax), xn, yn, zn
        real       e1, e2, xil, yil, zil, rxij, ryij, rzij
	real       rijsq, rij,rij2,hccterm
        parameter   ( twopi = 6.2831853, rsqpi = 0.5641896 )

        common / conf1 /  x, y, z, npart
        common / pot1 /   pi,sig,q,dcol
        common / sys1 /   box,hbox,temp,rab,zdna

        ncoll=0
        e1=0.
        e2=0.

c collision
      do 3 i=1,npart
      rxij = x(i) - ((-1)**idna*rab+deltar*0.70710681e0)
      ryij = y(i) - ((-1)**idna*rab+deltar*0.70710681e0)
      rij2=rxij*rxij + ryij*ryij
      if (rij2 .lt. dcol(i) ) ncoll=ncoll+1
3     continue
      hccterm=float(ncoll)/deltar
      fhs = (-1)**(idna+1)*2.*hccterm

      do 22 ii=1,zdna
      il=npart+(idna-1)*zdna+ii
      xn  = x(il) + deltar*0.707106781e0
      yn  = y(il) + deltar*0.707106781e0
      zn  = z(il) 

      xil = x(il)
      yil = y(il)
      zil = z(il)

c    ** calculate energy change in r-space

c il_th particle at old position
              do 5 i=1,npart
              rxij = x(i) - xil
              ryij = y(i) - yil
              rzij = z(i) - zil

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

              rijsq=rxij*rxij + ryij*ryij + rzij*rzij
              rij=sqrt(rijsq)
              e1=e1 + q(i)*q(il)/rij
5             continue

              do 11 jj=1,zdna
              ia=npart+zdna*(2-idna)+jj
              rxij = x(ia) - xil
              ryij = y(ia) - yil
              rzij = z(ia) - zil

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

              rijsq=rxij*rxij + ryij*ryij + rzij*rzij
              rij=sqrt(rijsq)
              e1=e1 + q(ia)*q(il)/rij
11            continue

c  at new position  
              do 6 i=1,npart
              rxij = x(i) - xn
              ryij = y(i) - yn
              rzij = z(i) - zn

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

              rijsq=rxij*rxij + ryij*ryij + rzij*rzij
	        rij=sqrt(rijsq)

              e2=e2 + q(i)*q(il)/rij
6             continue

	        do 12 jj=1,zdna
              ia=npart+zdna*(2-idna)+jj
              rxij = x(ia) - xn
              ryij = y(ia) - yn
              rzij = z(ia) - zn

              rxij = rxij - box * anint ( rxij / box )
              ryij = ryij - box * anint ( ryij / box )
              rzij = rzij - box * anint ( rzij / box )

              rijsq=rxij*rxij + ryij*ryij + rzij*rzij
              rij=sqrt(rijsq)
              e2=e2 + q(ia)*q(il)/rij
12            continue
22	continue
      
      der = e2 - e1
      fcc = (-1)**(idna+1)*der/deltar
      f=fcc+fhs
      end

 
      subroutine pforce(myid)
c________________________________________________________________________
c          nvt_mc for interaction between two particles in 
c                     an electrolyte solution
c__________________________________________________________________________
 
      implicit none
      integer equil, prod, ii, icycl, ndispl, attempt, 
     &  nacc, ncycl, nmoves, imove, il,ibeg,nstart,myid
      real     ener, dr, deltar

c     ---initialize sysem
      call readdat(myid,ibeg,equil, prod, ndispl, dr, deltar)

	if (ibeg.lt.2) then
	nstart=1
	else
	nstart=2
	endif

      nmoves = ndispl
c     ---total energy of the system
      call toterg(ener)
      write (*, 99001) ener

c     ---start mc-cycle
      do 100 ii = nstart, 2
c        --- ii=1 equilibration
c        --- ii=2 production
         if (ii.eq.1) then
            ncycl = equil
            if (ncycl.ne.0) write (*, *) ' start equilibration '
         else
            ncycl = prod
            if (ncycl.ne.0) write (*, *) ' start production '
c       --- initilize sampling variables 
	      call sample(myid,0,ener,dr,deltar)
         endif


c   ---intialize the subroutine that adjust the maximum displacement
         attempt = 0
         nacc = 0
         call adjust(attempt, nacc, dr)


         do 50 icycl = 1, ncycl
               do imove = 1, nmoves
c              ---attempt to displace a particle
               call mcmove(il, ener, attempt, nacc, dr)
               enddo

	       if (ii.eq.2) then
               call sample(myid,1, ener,dr,deltar)
               end if

	 if (mod(icycl,ncycl/10).eq.0) then
		write (*, *) '======>> done ', icycl,' out of ', ncycl
c              ---adjust maximum displacements
               call adjust(attempt, nacc, dr)					
             end if
50       end do

      write (6, 99003) attempt,nacc,float(nacc)/float(attempt) 

100   end do

	call sample(myid, 2, ener,dr,deltar)

	write(*,*) 'final average energy/nkt = ', ener      
 
99001 format (' total energy at initial configuration: ', f12.5, /)
99002 format (' total energy at the end of simulation: ', f12.5, /, 
     &        '       running energy              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5, /) 
99003 format (' number of att. to displ. a part.  : ', i10, /, 
     &        ' success: ', i10, '(= ', f5.2, '%)')
      end


      subroutine mcmove(il, ener, attempt, nacc, dr)
c    *******************************************************************
c
c     attempts to displace a randomly selected particle
c
c
c  ener   (input/output) : total energy
c  attemp (input/output) number of attemps that have been
c                  performed to displace a particle
c  nacc   (input/output) number of successful attemps
c                  to displace a particle
c  dr     (input) maximum displacement
c    *******************************************************************
      implicit none
      integer il, attempt, nacc,  npmax, npart, maxk
	integer zdna
      parameter (npmax=500, maxk = 1000)
      real     box,temp,hbox
      real     x(npmax), y(npmax), z(npmax),xil,yil,zil
      real     ener, xn, yn, zn, der, dr
      real     ranf,rab
      real     sig(npmax), q(npmax), pi,dcol(npmax)
      common /pot1/  pi,sig,q,dcol
      common / sys1 / box, hbox, temp, rab,zdna
      common / conf1 /  x, y, z, npart
	common /distan/   xil,yil,zil
c    *******************************************************************
	attempt = attempt + 1
c     ---select a particle at random
      il = int( float(npart)*ranf() ) + 1

c     ---give the particle a random displacement
      xn = x(il) + (ranf() -0.5e0)*dr
      yn = y(il) + (ranf() -0.5e0)*dr
      zn = z(il) + (ranf() -0.5e0)*dr

c     ---calculate energy change
	call deltae( il, xn, yn, zn, der )
c     ---acceptance test
      if ( (der.lt.299.e0).and.(ranf().lt.exp(-der)) ) then
c        --accepted
		nacc = nacc + 1
		ener = ener + der 
c        ---put particle in simulation box
         if (xn.lt.-hbox) xn = xn + box
         if (xn.gt.hbox)  xn = xn - box
         if (yn.lt.-hbox) yn = yn + box
         if (yn.gt.hbox)  yn = yn - box
         if (zn.lt.-hbox) zn = zn + box
         if (zn.gt.hbox)  zn = zn - box
         x(il) = xn
         y(il) = yn
         z(il) = zn
	else
c       --- rejected
      end if
      end

 
      subroutine readdat(myid,ibeg,equil,prod,ndispl,dr,deltar)
c     ---input solution condition and model parameters
c
c     ---input parameters: file: (3), input.dat
c    ibeg  =  0 : initialize from a lattice
c             1 : read configuration from disk
c             2 : restart sampling    
c    equil      : number of monte carlo cycles during equilibration
c    prod       : number of monte carlo cycles during production
c    dr         : maximum displacement
c    ndispl     : number of attemps to displace a particle per mc cycle
c    npart      : total number of particles
c    nparti(i)  : number of particles for species i
c    temp       : temperature, kelvin
c
c    nion   = number of ion species
c    nionmax= maxmum number of ion species
c    sig(i) = diameters of ions and macroions/bjerrum length
c    sigi(i)= diameters of ions and macroions species, in angstroms
c    q(i)   = valence for an ion
c    qi(i)  = valence for an ion, for input
c    cm(i)  = molar concentration of ions
c    lb0    = bjerrum length for water at 298.15 k, in angstroms
c    lb     = bjerrum length at system temperature, in angstroms
c    debye_l= debye length in units of a 
c    is     = ionic strength mol/l
c
c    box    = simulation box length in angstroms (for input)
c    rab    = separation between two macroions/diameter of the first macroion 
c     ---input parameters: file: iniconf.dat (restart file
c                to continue a simulation from disk)
c    npmax  = maximum number of particles
c    dr     = optimized maximum displacement old configurations
c    x(1),y(1),z(1)            : position first particle 1
c        ...
c    x(npart),y(npart),z(npart): position particle last particle
 
      implicit none
      real  pi, box, temp, hbox,deltar,con1, con2, dr, lb0, lb, shbox
      integer ibeg,equil,prod,i,j,k,ndispl,
     &        zero,nion,nionmax,npmax,npart,nhismax
      parameter (npmax=500, nionmax=10, lb0=7.144464e0,
     &	  con1=2.338874e4, con2=6.022137e-4,nhismax=20)
      real sig(npmax),q(npmax),debye_l,ddna,dcol(npmax)
      real is, sigi(nionmax), qi(nionmax), cm(nionmax)
      real x(npmax), y(npmax), z(npmax), rab,rab0(20)
      real enert0,mf0,mfhs0,mfcc0,phcc0
      integer ns0,nnode
      integer  zdna,nparti(nionmax),myid
      common /sys1/  box,hbox,temp, rab,zdna
      common /pot1/  pi,sig,q,dcol
      common /conf1/ x, y, z, npart
      common /samm/ nion, nparti
      common /sam/ shbox, lb
      common /restart/ ns0,enert0,mf0,mfhs0,mfcc0,phcc0,rab0
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
      open (65, file='input')
c     ---read simulation parameters
      read (65, *)
      read (65, *) nnode
      read (65, *)
      read (65, *) ibeg, equil, prod
      read (65, *)
      read (65, *) dr, deltar
      read (65, *)
      read (65, *) ndispl
      read (65, *)
      read (65, *) box
      read (65, *)
      read (65, *) (rab0(i),i=1,nnode)
c     ---read solution condition and model parameters
      read (65, *)
      read (65, *) temp
      read (65, *)
      read (65, *) nion
      read (65, *)
      do i=1,nion
      read (65, *) sigi(i), qi(i), nparti(i)
      enddo
      close(65)

c check charge neutrality at initial condition
	zdna=int(box/1.7e0)
c	write(*,*) zdna
	zero=0
        do i=1,nion
	zero=zero+int(qi(i))*nparti(i)
	enddo
	zero=zero-2*zdna
	if (zero.ne.0) then
	write(*,*)'charge neutrality is unsatisfied!',zero 
	stop
	endif

c total number of small ions
	npart=0
	do i=1,nion
	npart=npart+nparti(i)
	enddo
      if (npart.gt.npmax) then
      write (*, *) ' error: number of particles too large'
      stop
      end if

c     ---calculate parameters:
      pi = 3.14159265
      lb=lb0
c cross-section diameter of dna chain
	ddna=20.e0
c concentration, ionic strength, M, Debye length, 1/A
	do i=1,nion
	cm(i)=float(nparti(i))/(box**3*con2)
	enddo

	is=0.0
	do i=1,nion
	is=is+qi(i)**2*cm(i)/2.
        enddo
	debye_l=1./sqrt(8.*pi*lb*is*con2)

c     --- normalized with bjerrum length
	box = box/lb
        hbox = box/2.e0
	shbox = hbox*hbox

	k=0
	do 10 i = 1,nion
	do j = 1, nparti(i)
	k=k+1
	sig(k) = sigi(i)/lb
	dcol(k)=((ddna+sigi(i))/(2.*lb))**2
	q(k) = qi(i)
	enddo
10	continue

c   --- assign locations of two dna chains (20 a-dna diameter)
	rab=rab0(myid+1)*ddna*sqrt(2.e0)/4.e0/lb 
c  --- assign positions of charges
	do i=npart+1,npart+zdna
	x(i)=-rab
	y(i)=-rab
	z(i)=(-float(zdna)/2.0e0+float(i)-0.5)*1.7e0/lb
	x(i+zdna)=rab
	y(i+zdna)=rab
	z(i+zdna)=z(i)
	q(i)=-1.e0
	q(i+zdna)=-1.e0
	enddo

c     ---read/generate initial configuration
      if (ibeg.eq.0) then
         write (*, *) 'generate randomly initial configuration'
         call config
         open(64,file=iniconf)
         write (64, *) npart, dr
	 ns0=0
         enert0=0.0
         mf0=0
         mfhs0=0
         mfcc0=0
         phcc0=0
         write (64,*) ns0, enert0, mf0, mfhs0, mfcc0,phcc0
         do i = 1, npart
         write (64,  5) x(i), y(i), z(i)
         end do
         close(64)

      elseif (ibeg.eq.1) then
         write (*, *) ' read conf from disk '
	 open(64,file=iniconf)
	 read (64, *) npart, dr
         read(64,*) ns0, enert0, mf0, mfhs0, mfcc0,phcc0
         ns0=0
         enert0=0.0
         mf0=0
         mfhs0=0
         mfcc0=0
         phcc0=0
         do i = 1, npart
         read (64, *) x(i), y(i), z(i)
         end do
	 close(64)
      elseif (ibeg.eq.2) then
	 write(*,*) 'restart sampling'
         open(64,file=iniconf)
         read (64, *) npart, dr
         read(64,*) ns0, enert0, mf0, mfhs0, mfcc0,phcc0
         do i = 1, npart
         read (64, *) x(i), y(i), z(i)
         end do
         close(64)

      end if

c     ---write input data
      write (*,  1) equil, prod
      write (*,  2) ndispl, dr,rab0(myid+1)
      write (*,  3) npart, box*lb, debye_l, lb,is
      write (*, *) '#ions    diameter   charge  concentration(m)' 
      do i=1,nion
      write (*,  4) nparti(i), sigi(i), qi(i), cm(i)
      enddo 
      return
 1    format ('  number of equilibration cycles             :', i10, /, 
     &        '  number of production cycles                :', i10, /) 
 2    format ('  number of att. to displ. a part. per cycle :', i10, /, 
     &        '  maximum displacement                       :',f10.3,/, 
     &        '  macroion separation distance:              :',f10.2
     &        //)
 3    format ('  total number of particles                  :', i10, /, 
     &        '  box length,a                 :', f10.3/,
     &        '  debye length,a               :', f10.3/,
     &        '  bjerrum length,a             :', f10.6/,
     &        '  ionic strength, m            :', f10.6/)

 4    format (2x,i5,2x,f8.1,2x,f6.1,2x,e12.6)
 5    format(2x,3f12.6)
      end

      subroutine config
c
c     randomly place `npart' particles
c
      implicit none
      integer i, k, npart, npmax,zdna
      parameter (npmax=500)
      real x(npmax), y(npmax), z(npmax)
      real box, temp, hbox, rab
      real sig(npmax), q(npmax), dcol(npmax),pi
      real xx,xd, yd, zd, rik
      real ranf,a1,a2
      common /conf1/ x, y, z, npart
      common /sys1/ box,hbox,temp, rab,zdna
      common /pot1/  pi,sig,q,dcol

      do 1 i=npart,1,-1
3     xx=ranf()
      x(i)=hbox*(xx-0.5e0)*2.e0
      xx=ranf()
      y(i)=hbox*(xx-0.5e0)*2.e0
      a1=(x(i)+rab)**2+(y(i)+rab)**2
      a2=(x(i)-rab)**2+(y(i)-rab)**2
      if((a1.lt.dcol(i)).or.(a2.lt.dcol(i))) goto 3
      xx=ranf()
      z(i)=hbox*(xx-0.5e0)*2.e0
      do 2 k=i+1,npart
      xd=x(i)-x(k)
      if(xd.gt.hbox) xd=xd-box
      if(xd.lt.-hbox) xd=xd+box
      yd=y(i)-y(k)
      if(yd.gt.hbox) yd=yd-box
      if(yd.lt.-hbox) yd=yd+box
      zd=z(i)-z(k)
      if(zd.gt.hbox) zd=zd-box
      if(zd.lt.-hbox) zd=zd+box
      rik=sqrt(xd*xd+yd*yd+zd*zd)
      if (rik .lt. (sig(i)+sig(k))/2. ) go to 3
2     continue
1     continue
      end


 
      subroutine sample(myid,switch, ener,dr,deltar)
c
c    switch (input) = 0: initialize varibales 
c                   = 1: sample averages
c                   = 2: return average results
c
c  ener (input) : total energy/kt
c  enert        : average total energy/nkt
c  ns           : number of sample points
c  f            : force between macroions

      implicit none
      integer npmax, ns, switch, nionmax, nhismax
      parameter (npmax=500, nionmax=10, nhismax=20)
      integer  nion, nparti(nionmax)
      integer  i, nhcc
      integer  npart,myid 
      real     ener, enert, f, mf, fhs, fcc, mfhs, mfcc
      real     enert1, mf1, mfhs1, mfcc1,phcc,deltar
      real     xil,yil,zil
      real     sig(npmax), q(npmax), pi, lb,dcol(npmax)
      real     box, temp, hbox, shbox
      real     x(npmax), y(npmax), z(npmax), rab,sqrt3,rab0(20)
      real     dr
      real     enert0,mf0,mfhs0,mfcc0,phcc0
      real     f2(4),fhs2(4),fcc2(4)
      integer  ns0, ncx(4),ncoll,zdna

      parameter (sqrt3=1.7320508076)
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

c
      common /sys1/  box, hbox, temp, rab,zdna
      common /pot1/  pi, sig, q,dcol
      common /distan/ xil,yil,zil
      common /samm/ nion, nparti
      common /sam/ shbox, lb
      common /conf1/ x, y, z, npart
      common /restart/ ns0,enert0,mf0,mfhs0,mfcc0,phcc0,rab0
      save ns, enert, mf, mfcc,mfhs,nhcc,phcc
      output1=filename(myid+1)
      iniconf=filename2(myid+1)
c        ---initialize
      if (switch.eq.0) then
	  ns=ns0*10**5
	  nhcc=int( float(ns*4)*phcc0 )
c        ---total energy
	  enert=float(ns)*enert0
c        ---mean force between macroions
	  mf = float(ns*4)*mf0
	  mfhs= float(ns*4)*mfhs0
	  mfcc= float(ns*4)*mfcc0

	else if (switch.eq.1) then
	  ns=ns+1
c  ---sample excess energy 
	  enert = enert + ener
c  ---sample mean force between macroions
	call force( deltar,1,ncx(1),f2(1),fhs2(1),fcc2(1)) 
        call force(-deltar,1,ncx(2),f2(2),fhs2(2),fcc2(2))
        call force( deltar,2,ncx(3),f2(3),fhs2(3),fcc2(3))
        call force(-deltar,2,ncx(4),f2(4),fhs2(4),fcc2(4))
        f=0.
        fhs=0.
        fcc=0.
        ncoll=0
        do i=1,4
        f=f+f2(i)
        fhs=fhs+fhs2(i)
        fcc=fcc+fcc2(i)
        ncoll=ncoll+ncx(i)
        enddo
	mf = mf + f
	mfhs=mfhs+fhs
	mfcc=mfcc+fcc
        nhcc=nhcc+ncoll

        if (mod(ns,10**5).eq.0) then
        phcc=float(nhcc)/float(ns*4)
        enert1 = enert/float(ns)
        mf1    = mf/float(ns*4)
        mfhs1  = mfhs/float(ns*4)
        mfcc1  = mfcc/float(ns*4)
        open(67,file=output1,position='append')
        write(67,223) rab0(myid+1), ns/10**5, enert1,
     &  mf1, mfhs1, mfcc1,phcc,deltar
        close(67)

	open(69,file=iniconf)
   	write (69, *) npart, dr
        write(69,222) ns/10**5, enert1, mf1, mfhs1, mfcc1,phcc
        do i = 1, npart
        write (69, 10) x(i), y(i), z(i)
        end do
        close(69)
	endif
      else
        phcc=float(nhcc)/float(ns*4)
        enert1 = enert/float(ns)
        mf1    = mf/float(ns*4)
        mfhs1  = mfhs/float(ns*4)
        mfcc1  = mfcc/float(ns*4)

        open(69,file=iniconf)
        write (69, *) npart, dr
        write(69,222) ns/10**5, enert1, mf1, mfhs1, mfcc1,phcc
        do i = 1, npart
        write (69, 10) x(i), y(i), z(i)
        end do
        close(69)

        open(77,file='results',position='append')
	write(77,30) rab0(myid+1),enert1,mf1,mfhs1,mfcc1
	close(77)
      end if

10    format(2x,3(f11.6,1x))
30    format(1x,f4.2,1x,4e14.6)
222   format(1x,i5,1x,4e14.6,1x,f6.3)
223   format(1x,f5.2,i5,1x,4f10.2,1x,f5.2,1x,f6.4)
      end


      subroutine toterg(ener)
c *******************************************************************
c
c  ener (output) : total energy/kt
c
c *******************************************************************

      implicit none
      integer npmax,npart,zdna,i
      parameter (npmax=500)

      real     x(npmax),y(npmax),z(npmax)
      real     box,temp,hbox
      real     pi, sig(npmax),dcol(npmax), q(npmax) 
      real     vr, ener, rab,shbox,lb,eself

      common /sys1/  box,hbox,temp,rab,zdna
      common /pot1/  pi,sig,q,dcol
      common /conf1/ x, y, z, npart
      common /sam/ shbox, lb
        call rwald ( npart+2*zdna, vr, x, y, z, q, box)
	eself=0.e0
	do i=1,zdna/2
	eself=eself+2.e0*lb/(1.7e0*float(i))
	enddo
	eself=eself*float(zdna)*2.0e0
	ener=vr-eself
      end

