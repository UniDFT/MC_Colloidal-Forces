	Program  Main
       implicit none

#include <mpif.h>
#include "sprng_f.h"

       integer streamnum, nstreams, seed
       SPRNG_POINTER stream
       integer i, myid, nprocs, ierror, junk

C Initilize parallel mode
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)

c Initilize random number generator for multiprocessors
       streamnum = myid
       nstreams = nprocs

C  generate a new seed each time program is run and initialize stream
       seed = make_sprng_seed()
       stream = init_sprng(streamnum,nstreams,seed,SPRNG_DEFAULT)
c       print *, 'Printing information about new stream'
       junk = print_sprng(stream)

        call pforce(myid,stream)

c close parallel mode
       junk = free_sprng(stream)
       call MPI_FINALIZE(ierror)
       end



      SUBROUTINE ADJUST(Attemp, Nacc, Dr)
c
c     adjusts maximum displacement such that 50% of the
c     movels will be accepted
c
c  Attemp (input) number of attemps that have been performed to displace a particle
c  Nacc   (input) number of successful attemps to displace a particle
c  Dr     (output) new maximum displacement
c
 
      IMPLICIT NONE
      INTEGER Attemp, Nacc, attempp, naccp
      REAL     box,hbox
      REAL     dro, frac, Dr, RAB
      SAVE naccp, attempp
      COMMON /sys1/ box,hbox, RAB
 
      IF (Attemp.EQ.0.OR.attempp.GE.Attemp) THEN
         naccp = Nacc
         attempp = Attemp
      ELSE
         frac = FLOAT(Nacc-naccp)/FLOAT(Attemp-attempp)
         dro = Dr
         Dr = Dr*ABS(frac/0.5E0)
c        ---limit the change:
         IF (Dr/dro.GT.1.5E0) Dr = dro*1.5E0
         IF (Dr/dro.LT.0.5E0) Dr = dro*0.5E0
         IF (Dr.GT.HBOX/2.E0) Dr = HBOX/2.E0
         WRITE (*, 99001) Dr, dro, frac, Attemp - attempp, Nacc - naccp
c        ---store nacc and attemp for next use
         naccp = Nacc
         attempp = Attemp
      END IF
      RETURN
99001 FORMAT (' Max. displ. set to : ', f6.3, ' (old : ', f6.3, ')', /, 
     &        ' Frac. acc.: ', f4.2, ' attempts: ', i7, ' succes: ', i7)
      END


	  SUBROUTINE DELTAE( IL, XN, YN, ZN, DER, DEK )
*********************************************************************
**	Calculates energy change using Ewald method                  **
**	IL		index of particle to be changed						 **
**	XN		new x coordinate for particle IL 					 **
**	YN		new y coordinate for particle IL 					 **
**	ZN		new x coordinate for particle IL 					 **
**	DER		energy change in R-space after the move	    		 **
**	DEK		energy change in K-space after the move	    		 **
**	DOLD	old distance between particle IL and other particles **
**	DNEW	new distance between particle IL and other particles **
*********************************************************************
        IMPLICIT NONE
        INTEGER     MAXK, TOTK, KMAX, KX, KY, KZ, I, KSQMAX, KSQ, 
     &              IL, npmax, npart 
        PARAMETER   ( MAXK = 1000, KMAX = 5, KSQMAX = 27, npmax=5000)
        REAL         TWOPI, TWOPL, FACTOR, VD, VDN, RSQPI
        REAL       KVEC(MAXK), KAPPA, DEK,DER
        REAL       box, hbox, DNEW(npmax), DOLD(npmax)
        REAL       sig(npmax), q(npmax), pi
        REAL       x(npmax), y(npmax), z(npmax), XN, YN, ZN
        REAL       E1, E2, XIL, YIL, ZIL, RXIJ, RYIJ, RZIJ
	  REAL       RIJSQ, RIJ, KRIJ
	  REAL       ERFC, RAB
        PARAMETER   ( TWOPI = 6.2831853, RSQPI = 0.5641896 )

        COMPLEX   EIKX(0:KMAX), EIKY(-KMAX:KMAX), EIKZ(-KMAX:KMAX)
        COMPLEX   DEIKX(0:KMAX),DEIKY(-KMAX:KMAX),DEIKZ(-KMAX:KMAX)
        COMPLEX   DEIKR, SUMI(MAXK), DSUMI(MAXK)

        COMMON / BLOCK2 / KVEC, KAPPA
	  COMMON / MOVE /   SUMI, VD, DSUMI, VDN
        COMMON / conf1 /  x, y, z, npart
        COMMON / pot1 /   pi,sig,q
        COMMON / sys1 /   box, hbox, RAB
	  COMMON /DISTAN/   DOLD, DNEW,XIL,YIL,ZIL
C    ** CALCULATE ENERGY CHANGE IN R-SPACE
		E1=0.
		E2=0.

		XIL=X(IL)
		YIL=Y(IL)
		ZIL=Z(IL)

		DO 5 I=1,NPART+2
		IF(I.EQ.IL) GOTO 5

C IL-th PARTICLE AT OLD POSITION

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

		DO 6 I=1,NPART+2
		IF(I.EQ.IL) GOTO 6
C IL-th PARTICLE AT NEW POSITION
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
			  RETURN
		      ELSE
			  ENDIF

			KRIJ=KAPPA*RIJ/BOX
	        E2=E2 + Q(I)*Q(IL)*ERFC(REAL(KRIJ))/RIJ
6		CONTINUE


      DER = E2 - E1

C    ** CALCULATE ENERGY CHANGE IN K-SPACE
C    ** CONSTRUCT EXP(IK.R) FOR IL-th ION AND K-VECTORS **
C    **   AT NEW AND OLD POSITIONS                      **
C    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY  **
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

C    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **
C
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
C
C    ** SUM OVER ALL VECTORS **

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
            DEIKR=DEIKX(KX)*DEIKY(KY)*DEIKZ(KZ)
     &			    -EIKX(KX)*EIKY(KY)*EIKZ(KZ)
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

************************************************************************************************************************
C    *******************************************************************
C    ** REAL-SPACE AND RECIPROCAL-SPACE PARTS OF EWALD SUM FOR IONS.  **
C    **                                                               **
C    ** REFERENCES:                                                   **
C    **                                                               **
C    ** WOODCOCK AND SINGER, TRANS. FARADAY SOC. 67, 12, 1971.        **
C    ** DE LEEUW ET AL., PROC. ROY. SOC. A 373, 27, 1980.             **
C    ** HEYES, J. CHEM. PHYS. 74, 1924, 1981.                         **
C    ** SEE ALSO FINCHAM, MDIONS, CCP5 PROGRAM LIBRARY.               **
C    **                                                               **
C    ** ROUTINES SUPPLIED:                                            **
C    **                                                               **
C    ** SUBROUTINE SETUP                                              **
C    **    SETS UP THE WAVEVECTORS FOR USE IN THE EWALD SUM           **
C    ** SUBROUTINE RWALD ( VR )                                       **
C    **    CALCULATES THE R-SPACE PART OF THE SUM                     **
C    ** SUBROUTINE KWALD ( VK )                                       **
C    **    CALCULATES THE K-SPACE PART OF THE SUM                     **
C    ** REAL FUNCTION ERFC ( X )                                      **
C    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    **  TOTK         THE TOTAL NUMBER OF K-VECTORS STORED            **
C    **  MAXK         MAXIMUM POSSIBLE NUMBER OF K-VECTORS            **
C    **  KMAX         MAX INTEGER COMPONENT OF THE K-VECTOR           **
C    **  KSQMAX       MAX SQUARE MOD OF THE K-VECTOR REQUIRED         **
C    **  VR           ENERGY FROM R-SPACE SUM/kT                      **
C    **  VK           ENERGY FROM K-SPACE SUM/kT                      **
C    **  KVEC(MAXK)   ARRAY USED TO STORE K-VECTORS                   **
C    **  KAPPA        WIDTH OF CANCELLING DISTRIBUTION                **
C    **                                                               **
C    **  (KAPPA IS A INPUT VALUE THAT SHOULD BE DETERMINED BY 	    **
C    **   COMPUTATIONAL EFFICIENCY, KAPPA=5/BOX FOR INITIAL TRY       **
C    **   KAPPA=(ALPHA)**0.5 IN FRENKEL & SMIT'S NOTATION)            **
C    **                                                               **
C    ** USAGE:                                                        **
C    **                                                               **
C    ** SETUP IS CALLED ONCE AT THE BEGINNING OF THE SIMULATION       **
C    ** TO CALCULATE ALL THE K-VECTORS REQUIRED IN THE EWALD SUM.     **
C    ** THESE VECTORS ARE USED THROUGHOUT THE SIMULATION IN THE       **
C    ** SUBROUTINE KWALD TO CALCULATE THE K-SPACE CONTRIBUTION TO THE **
C    ** POTENTIAL ENERGY AT EACH CONFIGURATION. THE SELF TERM IS      **
C    ** SUBTRACTED FROM THE K-SPACE CONTRIBUTION IN KWALD.            **
C    ** THE SURFACE TERM FOR SIMULATIONS IN VACUUM IS NOT INCLUDED.   **
C    ** ROUTINE RWALD RETURNS THE R-SPACE CONTRIBUTION TO THE EWALD   **
C    ** SUM AND IS CALLED FOR EACH CONFIGURATION IN THE SIMULATION.   **
C    ** A CUBIC BOX AND UNIT BOX LENGTH ARE ASSUMED THROUGHOUT.       **
C    *******************************************************************



        SUBROUTINE SETUP(GW) 

C    *******************************************************************
C    ** ROUTINE TO SET UP THE WAVE-VECTORS FOR THE EWALD SUM.         **
C    ** GW: Gaussian Width Parameter                                  **
C    ** IN THIS EXAMPLE WE ALLOW A MAXIMUM OF 1000 WAVEVECTORS.       **
C    *******************************************************************
        IMPLICIT NONE
        INTEGER     MAXK
        PARAMETER ( MAXK = 1000 )
        REAL      KVEC(MAXK), KAPPA, GW
        INTEGER     KMAX, KSQMAX, KSQ, KX, KY, KZ, TOTK
        REAL     TWOPI, B, RKX, RKY, RKZ, RKSQ
        PARAMETER ( KMAX = 5, KSQMAX = 27 , TWOPI = 6.2831853 )

        COMMON / BLOCK2 / KVEC, KAPPA
C    *******************************************************************
C
	  KAPPA = GW
C
        B = 1.0 / 4.0 / KAPPA / KAPPA

C    ** LOOP OVER K-VECTORS. NOTE KX IS NON-NEGATIVE **

        TOTK = 0

        DO 100 KX = 0, KMAX

           RKX = TWOPI * FLOAT ( KX )

           DO 99 KY = -KMAX, KMAX

              RKY = TWOPI * FLOAT ( KY )

              DO 98 KZ = -KMAX, KMAX

                 RKZ = TWOPI * FLOAT ( KZ )

                 KSQ = KX * KX + KY * KY + KZ * KZ

                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN

                    TOTK = TOTK + 1

                    IF ( TOTK .GT. MAXK ) STOP 'KVEC IS TOO SMALL'

                    RKSQ = RKX * RKX + RKY * RKY + RKZ * RKZ
                    KVEC(TOTK) = TWOPI * EXP ( -B * RKSQ ) / RKSQ

                 ENDIF

98            CONTINUE

99         CONTINUE

100     CONTINUE

        WRITE( *, ' ( '' EWALD SUM SETUP COMPLETE ''     ) ' )
        WRITE( *, ' ( '' NUMBER OF WAVEVECTORS IS '', I5 ) ' ) TOTK

        RETURN
        END



        SUBROUTINE RWALD ( N, VR ,RX, RY, RZ, Q, BOX)

C    *******************************************************************
C    ** CALCULATES R-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    **         N                     NUMBER OF IONS (=NPART)         **
C    **         RX(N),RY(N),RZ(N)     POSITIONS OF IONS/lB            **
C    **         Q(N)                  IONIC VALENCE                   **
C    **         VR                    R-SPACE POTENTIAL ENERGY/kT     **
C    **         BOX                   SIMULATION BOX LENGTH/lB        **
C    **                                                               **
C    ** ROUTINE REFERENCED:                                           **
C    **                                                               **
C    ** FUNCTION ERFC ( X )                          **
C    **    RETURNS THE COMPLEMENTARY ERROR FUNCTION                   **
C    *******************************************************************
        IMPLICIT NONE
        INTEGER     MAXK, npmax
        PARAMETER ( MAXK = 1000, npmax=5000 )
        REAL      KVEC(MAXK), KAPPA

        INTEGER              N, I, J
        REAL         RX(npmax), RY(npmax), RZ(npmax), Q(npmax)
        REAL         VR, BOX
        REAL         RXI, RYI, RZI, ZI, RXIJ, RYIJ, RZIJ
        REAL         RIJSQ, RIJ, KRIJ, VIJ
	REAL ERFC

        COMMON / BLOCK2 / KVEC, KAPPA
C    *******************************************************************

        VR = 0.0

        DO 100 I = 1, N - 1

           RXI = RX(I)
           RYI = RY(I)
           RZI = RZ(I)
           ZI  = Q(I)

           DO 99 J = I + 1, N

              RXIJ = RXI - RX(J)
              RYIJ = RYI - RY(J)
              RZIJ = RZI - RZ(J)

              RXIJ = RXIJ - BOX * ANINT ( RXIJ /BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ /BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ /BOX )

              RIJSQ = RXIJ * RXIJ + RYIJ * RYIJ + RZIJ * RZIJ
              RIJ   = SQRT ( RIJSQ )
              KRIJ  = KAPPA * RIJ	/ BOX
              VIJ   = ZI * Q(J) * ERFC(real( KRIJ )) / RIJ

              VR    = VR + VIJ

99         CONTINUE

100     CONTINUE

        RETURN
        END



        SUBROUTINE KWALD ( N, VK, RX, RY, RZ, Q, BOX )

C    *******************************************************************
C    ** CALCULATES K-SPACE PART OF POTENTIAL ENERGY BY EWALD METHOD.  **
C    **                                                               **
C    ** THE SELF TERM IS SUBTRACTED.                                  **
C    ** IN ONE COORDINATE DIRECTION (X), SYMMETRY IS USED TO REDUCE   **
C    ** THE SUM TO INCLUDE ONLY POSITIVE K-VECTORS.                   **
C    ** THE NEGATIVE VECTORS IN THIS DIRECTION ARE INCLUDED BY USE    **
C    ** OF THE MULTIPLICATIVE VARIABLE 'FACTOR'.                      **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    **  N                   NUMBER OF IONS                           **
C    **  RX(N),RY(N),RZ(N)   POSITIONS OF IONS, IN DIMENSION OF lB    **
C    **  Q(N)                IONIC VALENCES                           **
C    **  VK                  K-SPACE POTENTIAL ENERGY/kT              **
C    **  VKS                 SELF PART OF K-SPACE SUM/kT              **
C    **  BOX                 SIMULATION BOX LENGTH/lB                 **
C    *******************************************************************
        IMPLICIT NONE
        INTEGER  MAXK, N, TOTK, KMAX, KX, KY, KZ, I, KSQMAX, KSQ, npmax
        PARAMETER   ( MAXK = 1000, KMAX = 5, KSQMAX = 27, npmax=5000 )
        REAL       KVEC(MAXK)
        REAL       RX(npmax), RY(npmax), RZ(npmax), Q(npmax)
        REAL       KAPPA, VK, BOX
        REAL         TWOPI, TWOPL, FACTOR, VD, VS, RSQPI, VDN
        PARAMETER   ( TWOPI = 6.2831853, RSQPI = 0.5641896 )

        COMPLEX       EIKX(1:npmax, 0:KMAX)
        COMPLEX       EIKY(1:npmax, -KMAX:KMAX)
        COMPLEX       EIKZ(1:npmax, -KMAX:KMAX)
        COMPLEX       EIKR(npmax), SUMI(MAXK), DSUMI(MAXK)

        COMMON / BLOCK2 / KVEC, KAPPA
	  COMMON / MOVE /   SUMI, VD, DSUMI, VDN
C    *******************************************************************

C    ** CONSTRUCT EXP(IK.R) FOR ALL IONS AND K-VECTORS **

C    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY **

        TWOPL=TWOPI / BOX

        DO 10 I = 1, N

           EIKX(I, 0) = (1.0, 0.0)
           EIKY(I, 0) = (1.0, 0.0)
           EIKZ(I, 0) = (1.0, 0.0)

           EIKX(I, 1) = CMPLX ( COS ( TWOPL * RX(I) ) ,
     :                           SIN ( TWOPL * RX(I) ) )
           EIKY(I, 1) = CMPLX ( COS ( TWOPL * RY(I) ) ,
     :                           SIN ( TWOPL * RY(I) ) )
           EIKZ(I, 1) = CMPLX ( COS ( TWOPL * RZ(I) ) ,
     :                           SIN ( TWOPL * RZ(I) ) )

           EIKY(I, -1) = CONJG ( EIKY(I, 1) )
           EIKZ(I, -1) = CONJG ( EIKZ(I, 1) )

10      CONTINUE

C    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **

        DO 12 KX = 2, KMAX

           DO 11 I = 1, N

              EIKX(I, KX) = EIKX(I, KX-1) * EIKX(I, 1)

11         CONTINUE

12      CONTINUE

        DO 14 KY = 2, KMAX

           DO 13 I = 1, N

              EIKY(I,  KY) = EIKY(I, KY-1) * EIKY(I, 1)
              EIKY(I, -KY) = CONJG ( EIKY(I, KY) )

13         CONTINUE

14      CONTINUE

        DO 16 KZ = 2, KMAX

           DO 15 I = 1, N

              EIKZ(I,  KZ) = EIKZ(I, KZ-1) * EIKZ(I, 1)
              EIKZ(I, -KZ) = CONJG ( EIKZ(I, KZ) )

15         CONTINUE

16      CONTINUE

C    ** SUM OVER ALL VECTORS **

        VD   = 0.0
        TOTK = 0

        DO 24 KX = 0, KMAX

           IF ( KX .EQ. 0 ) THEN

              FACTOR = 1.0

           ELSE

              FACTOR = 2.0

           ENDIF

           DO 23 KY = -KMAX, KMAX

              DO 22 KZ = -KMAX, KMAX

                 KSQ = KX * KX + KY * KY + KZ * KZ

                 IF ( ( KSQ .LT. KSQMAX ) .AND. ( KSQ .NE. 0 ) ) THEN

                    TOTK = TOTK + 1
                    SUMI(TOTK)  = (0.0, 0.0)

                    DO 21 I = 1, N

                       EIKR(I) = EIKX(I, KX) * EIKY(I, KY) * EIKZ(I, KZ)
                       SUMI(TOTK) = SUMI(TOTK) + Q(I) * EIKR(I)

21                  CONTINUE

                    VD = VD + FACTOR * KVEC(TOTK) 
     &	         		  * CONJG ( SUMI(TOTK) ) * SUMI(TOTK)

                 ENDIF

22            CONTINUE

23         CONTINUE

24      CONTINUE

	 VD=VD / BOX

C    ** CALCULATES SELF PART OF K-SPACE SUM **

        VS = 0.0

        DO 25 I = 1, N

           VS = VS + Q(I) * Q(I)

25      CONTINUE

        VS = RSQPI * KAPPA * VS / BOX

C    ** CALCULATE THE TOTAL K-SPACE POTENTIAL **

        VK = VD - VS

        RETURN
        END


	  SUBROUTINE FORCE(DELTAR,IL,IA,ncoll,F, FHS, FCC)
*********************************************************************
**	Calculates energy change using Ewald method                  **
**	IL		index of particle to be changed						 **
**	XN		new x coordinate for particle IL 					 **
**	YN		new y coordinate for particle IL 					 **
**	ZN		new x coordinate for particle IL 					 **
**	DER		energy change in R-space after the move	    		 **
**	DEK		energy change in K-space after the move	    		 **
**	F	    force between macroions                              **
*********************************************************************
        IMPLICIT NONE
        INTEGER     MAXK, TOTK, KMAX, KX, KY, KZ, I, KSQMAX, KSQ, 
     &               npmax, npart,  ncoll ,IL,IA,stream
        PARAMETER   ( MAXK = 1000, KMAX = 5, KSQMAX = 27, npmax=5000)
        REAL       TWOPI, TWOPL, FACTOR, VD, VDN, RSQPI
        REAL       KVEC(MAXK), KAPPA, DEK,DER
	  REAL       box, hbox, RAB, XX, DELTAR, F, FCC, FHS
	  REAL       sig(npmax), q(npmax), pi
        REAL       x(npmax), y(npmax), z(npmax), XN, YN, ZN
        REAL       E1, E2, XIL, YIL, ZIL, RXIJ, RYIJ, RZIJ
	  REAL       RIJSQ, RIJ, KRIJ,hccterm
	  REAL       ERFC
        PARAMETER   ( TWOPI = 6.2831853, RSQPI = 0.5641896 )

        COMPLEX   EIKX(0:KMAX), EIKY(-KMAX:KMAX), EIKZ(-KMAX:KMAX)
        COMPLEX   DEIKX(0:KMAX),DEIKY(-KMAX:KMAX),DEIKZ(-KMAX:KMAX)
        COMPLEX   DEIKR, SUMI(MAXK), DSUMI(MAXK)

        COMMON / BLOCK2 / KVEC, KAPPA
	COMMON / MOVE /   SUMI, VD, DSUMI, VDN
        COMMON / conf1 /  x, y, z, npart
        COMMON / pot1 /   pi,sig,q
        COMMON / sys1 /   box, hbox, RAB

      XN  = X(IL) + DELTAR/RAB*(X(IL)-X(IA))
      YN  = Y(IL) + DELTAR/RAB*(Y(IL)-Y(IA))
      ZN  = Z(IL) + DELTAR/RAB*(Z(IL)-Z(IA))

      XIL = X(IL)
      YIL = Y(IL)
      ZIL = Z(IL)

C    ** CALCULATE ENERGY CHANGE IN R-SPACE
                E1=0.
                E2=0.

C IL_th PARTICLE AT OLD POSITION
                DO 5 I=1,NPART
                        RXIJ = X(I) - XIL
                        RYIJ = Y(I) - YIL
                        RZIJ = Z(I) - ZIL

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

                RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
                        RIJ=SQRT(RIJSQ)
                        KRIJ=KAPPA*RIJ/BOX
                E1=E1 + Q(I)*Q(IL)*erfc(KRIJ)/RIJ
5               CONTINUE

                        RXIJ = X(IA) - XIL
                        RYIJ = Y(IA) - YIL
                        RZIJ = Z(IA) - ZIL

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

                RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
                        RIJ=SQRT(RIJSQ)
                        KRIJ=KAPPA*RIJ/BOX
                E1=E1 + Q(IA)*Q(IL)*erfc(KRIJ)/RIJ

C IL_th PARTICLE AT NEW POSITION
                hccterm=0.0
                ncoll=0
                DO 6 I=1,NPART
                        RXIJ = X(I) - XN
                        RYIJ = Y(I) - YN
                        RZIJ = Z(I) - ZN

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

                RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
                        RIJ=SQRT(RIJSQ)
                        IF (RIJ .LT. (SIG(IL)+SIG(I))/2.0 )
     &         ncoll=ncoll+1

                        KRIJ=KAPPA*RIJ/BOX
                E2=E2 + Q(I)*Q(IL)*erfc(KRIJ)/RIJ
6               CONTINUE

                        RXIJ = X(IA) - XN
                        RYIJ = Y(IA) - YN
                        RZIJ = Z(IA) - ZN

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

                RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
                        RIJ=SQRT(RIJSQ)
                        KRIJ=KAPPA*RIJ/BOX
                E2=E2 + Q(IA)*Q(IL)*erfc(KRIJ)/RIJ
      hccterm=float(ncoll)/deltar
      DER = E2 - E1
C    ** CALCULAT
C    ** CONSTRUCT EXP(IK.R) FOR IL-th ION AND K-VECTORS **
C    **   AT NEW AND OLD POSITIONS                      **
C    ** CALCULATE KX, KY, KZ = 0 , -1 AND 1 EXPLICITLY  **
       TWOPL=TWOPI / BOX
       EIKX(0)=(1.0,0.0)
       EIKY(0)=(1.0,0.0)
       EIKZ(0)=(1.0,0.0)
       EIKX(1)=CMPLX(COS(TWOPL*XIL),SIN(TWOPL*XIL))
       EIKY(1)=CMPLX(COS(TWOPL*YIL),SIN(TWOPL*YIL))
       EIKZ(1)=CMPLX(COS(TWOPL*ZIL),SIN(TWOPL*ZIL))
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

C    ** CALCULATE REMAINING KX, KY AND KZ BY RECURRENCE **
C
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
C
C    ** SUM OVER ALL VECTORS **

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
            DEIKR=DEIKX(KX)*DEIKY(KY)*DEIKZ(KZ)
     &			    -EIKX(KX)*EIKY(KY)*EIKZ(KZ)
            DSUMI(TOTK)=SUMI(TOTK)+Q(IL)*DEIKR
          VDN=VDN+FACTOR*KVEC(TOTK)*CONJG(DSUMI(TOTK))*DSUMI(TOTK)
         ENDIF
10      CONTINUE
20     CONTINUE
30    CONTINUE

	VDN = VDN / BOX
	DEK = VDN - VD

	FCC = -(DEK+DER) / DELTAR
	FHS = -2.*hccterm
	F=FCC+FHS
      RETURN
      END

 
	subroutine pforce(myid,stream)
c________________________________________________________________________
c
c          NVT_MC FOR INTERACTION BETWEEN TWO PARTICLES IN 
c                     AN ELECTROLYTE SOLUTION
c__________________________________________________________________________
 
      IMPLICIT NONE
      INTEGER equil, prod, nsamp, ii, icycl, ndispl, attempt, 
     &  stream,nacc, ncycl, nmoves, imove, IL,ibeg,nstart,myid
      REAL     Ener, Enert, Dr, GW,deltar
c     ---initialize sysem
      CALL READDAT(myid,ibeg,equil, prod, nsamp, ndispl, dr, GW,deltar)
	IF (ibeg.LT.2) THEN
	nstart=1
	ELSE
	nstart=2
	ENDIF
      nmoves = ndispl
c     ---total energy of the system
      CALL TOTERG(GW, Ener)
      WRITE (*, 99001) Ener

c     ---start MC-cycle
      DO 100 ii = nstart, 2
c        --- ii=1 equilibration
c        --- ii=2 production
         IF (ii.EQ.1) THEN
            ncycl = equil
            IF (ncycl.NE.0) WRITE (*, *) ' Start equilibration '
         ELSE
            ncycl = prod
            IF (ncycl.NE.0) WRITE (*, *) ' Start production '
         ENDIF
c       --- initilize sampling variables 
		CALL SAMPLE(myid,ibeg,IL, 0, Ener,Dr,deltar,stream)
         attempt = 0
         nacc = 0

c   ---intialize the subroutine that adjust the maximum displacement
         CALL ADJUST(attempt, nacc, dr)

         DO 50 icycl = 1, ncycl
               DO imove = 1, nmoves
c              ---attempt to displace a particle
               CALL MCMOVE(myid,IL, Ener, attempt, nacc, dr,stream)
               ENDDO

	       IF (ii.EQ.2) THEN
               CALL SAMPLE(myid,ibeg,IL, 1, Ener,Dr,deltar,stream)
               END IF

            IF (MOD(icycl,ncycl/5).EQ.0) THEN
               WRITE (*, *) '======>> Done ', icycl, ' out of ', ncycl
c              ---adjust maximum displacements
               CALL ADJUST(attempt, nacc, Dr)						
            END IF
50       END DO

         IF (ncycl.NE.0) THEN
            IF (attempt.NE.0) WRITE (6, 99003) attempt, nacc, 
     &                               100.*FLOAT(nacc)/FLOAT(attempt)
	ENDIF
100   END DO
	CALL SAMPLE(myid,ibeg,IL, 2, Ener,Dr,deltar,stream)
	WRITE(*,*) 'FINAL AVERAGE ENERGY/NkT = ', Ener      
 
99001 FORMAT (' Total energy at initial configuration: ', f12.5, /)
99002 FORMAT (' Total energy at the end of simulation: ', f12.5, /, 
     &        '       running energy              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5, /) 
99003 FORMAT (' Number of att. to displ. a part.  : ', i10, /, 
     &        ' success: ', i10, '(= ', f5.2, '%)')
      END

      SUBROUTINE MCMOVE(myid,IL, Ener, Attempt, Nacc, Dr,stream)
C    *******************************************************************
c
c     attempts to displace a randomly selected particle
c
c
c  Ener   (input/output) : total energy
c  Attemp (input/output) number of attemps that have been
c                  performed to displace a particle
c  Nacc   (input/output) number of successful attemps
c                  to displace a particle
c  Dr     (input) maximum displacement
C    *******************************************************************
      IMPLICIT NONE
      INTEGER IL, Attempt, Nacc,  npmax, npart, MAXK, I,stream,myid
      PARAMETER (npmax=5000, MAXK = 1000)
      REAL     box,hbox, DNEW(npmax), DOLD(npmax)
      REAL     x(npmax), y(npmax), z(npmax),XIL,YIL,ZIL
      REAL     Ener, xn, yn, zn, DER, DEK, Dr
      real sprng
	REAL     VD, VDN, RAB 
        REAL     sig(npmax), q(npmax), pi
      COMMON /pot1/  pi,sig,q
      COMPLEX   SUMI(MAXK),DSUMI(MAXK)
        COMMON / sys1 / box, hbox, RAB
      COMMON / conf1 /  x, y, z, npart
	COMMON / MOVE /   SUMI, VD, DSUMI, VDN
	COMMON /DISTAN/   DOLD, DNEW,XIL,YIL,ZIL
C    *******************************************************************
	Attempt = Attempt + 1
c     ---select a particle at random
      IL = INT( float(NPART)*sprng(stream) ) + 1

c     ---give the particle a random displacement
      xn = X(IL) + (sprng(stream) -0.5E0)*Dr
      yn = Y(IL) + (sprng(stream) -0.5E0)*Dr
      zn = Z(IL) + (sprng(stream) -0.5E0)*Dr

c     ---calculate energy change
	CALL DELTAE( IL, XN, YN, ZN, DER, DEK )
c     ---acceptance test
      IF (sprng(stream) .LT.EXP(-(DER+DEK))) THEN
c        --accepted
		Nacc = Nacc + 1

		DO I=1,NPART
		DOLD(I)=DNEW(I)
		ENDDO
			
		VD=VDN
		DO I=1,MAXK
		SUMI(I) = DSUMI(I)
		ENDDO

		Ener = Ener + DER + DEK
c	write(*,*) Ener, Ener - der - dek, der,dek
c        ---put particle in simulation box
         IF (xn.LT.-hbox) xn = xn + BOX
         IF (xn.GT.hbox) xn = xn - BOX
         IF (yn.LT.-hbox) yn = yn + BOX
         IF (yn.GT.hbox) yn = yn - BOX
         IF (zn.LT.-hbox) zn = zn + BOX
         IF (zn.GT.hbox) zn = zn - BOX
         X(IL) = xn
         Y(IL) = yn
         Z(IL) = zn
	ELSE
c       --- rejected
      END IF
 5    FORMAT(2X,6F12.6)
      RETURN
      END

 
      SUBROUTINE
     &  READDAT(myid,ibeg,Equil, Prod, Nsamp, Ndispl, Dr, GW,deltar)
C     ---input solution condition and model parameters
c
c     ---input parameters: file: (3), input.dat
c    ibeg  =  0 : initialize from a lattice
c             1 : read configuration from disk
c             2 : restart sampling    
c    Equil      : number of Monte Carlo cycles during equilibration
c    Prod       : number of Monte Carlo cycles during production
c    Nsamp      : number of Monte Carlo cycles between two sampling periods
c    Dr         : maximum displacement
c    Ndispl     : number of attemps to displace a particle per MC cycle
c    NPART      : total number of particles
C    NPARTI(I)  : number of particles for species i
c
c    NION   = number of ion species
c    NIONmax= maxmum number of ion species
c    sig(i) = diameters of ions and macroions/Bjerrum length
c    sigi(i)= diameters of ions and macroions species, in Angstroms
c    q(i)   = valence for an ion
c    qi(i)  = valence for an ion, for input
c    CM(i)  = molar concentration of ions
c    lB0    = Bjerrum length, input
c    lB     = Bjerrum length, in Angstroms
c    Debye_L= Debye length in units of A 
c    IS     = Ionic Strength mol/l
c
c    BOX    = simulation box length in Angstroms (FOR INPUT)
c    RAB    = separation between two macroions/diameter of the first macroion 
c     ---input parameters: file: iniconf.dat (restart file
c                to continue a simulation from disk)
c    npmax  = maximum number of particles
c    Dr     = optimized maximum displacement old configurations
c    GW     = Gaussian Width parameter, adjust for efficiency of Ewald Sum
c    X(1),Y(1),Z(1)            : position first particle 1
c        ...
c    X(NPART),Y(NPART),Z(NPART): position particle last particle
 
      IMPLICIT NONE
      REAL      pi, box, hbox,deltar,
     &                 con1, con2, Dr, lB, GW, SHBOX
      INTEGER ibeg,Equil,Prod,i,j,k,Ndispl,Nsamp,
     &        NION,NIONmax,npmax,npart,NHIsmax
      PARAMETER (npmax=5000, NIONmax=10,
     &	  con1=2.338874e4, con2=6.022137e-4,NHIsmax=200)
      REAL     sig(npmax), q(npmax), ZERO, Debye_L
      REAL     IS, sigi(NIONmax), qi(NIONmax), CM(NIONmax)
      REAL     x(npmax), y(npmax), z(npmax), RAB, RAB0, lB0(10)
        real Enert0,MF0,MFHS0,MFCC0,phcc0
        integer ns0,nnode
      INTEGER  NPARTI(NIONmax),myid
      COMMON /sys1/  box, hbox, RAB
      COMMON /pot1/  pi,sig,q
      COMMON /conf1/ x, y, z, npart
      COMMON /samm/ NION, NPARTI
      COMMON /sam/ SHBOX, lB
        COMMON /restart/ ns0,Enert0,MF0,MFHS0,MFCC0,phcc0
        character filename(20)*5,filename2(20)*5
        character OUTPUT1*5,INICONF*5,INPUT*5
        data filename/'OUT1','OUT2','OUT3','OUT4','OUT5','OUT6'
     &      ,'OUT7','OUT8','OUT9','OUT10','OUT11','OUT12'
     &      ,'OUT13','OUT14','OUT15','OUT16','OUT17','OUT18'
     &      ,'OUT19','OUT20'/
        data filename2/'INP1','INP2','INP3','INP4','INP5','INP6'
     &      ,'INP7','INP8','INP9','INP10','INP11','INP12'
     &      ,'INP13','INP14','INP15','INP16','INP17','INP18'
     &      ,'INP19','INP20'/

	INICONF=filename2(myid+1)
	OUTPUT1=filename(myid+1)
      OPEN (65, FILE='INPUT')
c     ---read simulation parameters
      READ (65, *)
      READ (65, *) nnode
      READ (65, *)
      READ (65, *) ibeg, Equil, Prod, Nsamp
      READ (65, *)
      READ (65, *) Dr, GW, deltar
      READ (65, *)
      READ (65, *) Ndispl
      READ (65, *)
      READ (65, *) BOX
      READ (65, *)
      READ (65, *) RAB0
c     ---read solution condition and model parameters
      READ (65, *)
      READ (65, *) (lB0(i), i=1,nnode)
      READ (65, *)
      READ (65, *) NION
      READ (65, *)
      DO i=1,NION
      READ (65, *) sigi(i), qi(i), NPARTI(i)
      ENDDO
      CLOSE(65)

	lB=lB0(myid+1)

	NPART=0
	DO i=1,NION-1
	NPART=NPART+NPARTI(i)
	ENDDO

c     ---calculate parameters:
      PI = 3.14159265
c     ---ion concentrations, M
	DO i=1,NION
	CM(I)=FLOAT(NPARTI(I))/(BOX**3*CON2)
	ENDDO

c     ---ionic strength and Debye screening length (Angstrom), not include macroions
	IS=0.0
	do i=1,NION
	IS=IS+qi(i)**2*CM(I)/2.
        enddo
	Debye_L=1./sqrt(8.*PI*lB*IS*CON2)

c     --- normalized with Bjerrum length which should have unit of Angstrom
	BOX = BOX/lB
      HBOX = BOX/2.E0
	SHBOX = HBOX*HBOX

	k=0
	DO 10 i = 1,NION
	DO j = 1, NPARTI(i)
	k=k+1
	sig(k) = sigi(i)/lB
	q(k) = qi(i)
	ENDDO
10	CONTINUE

c   --- assign LOCATION of two macroions ( NPART+1, NPART+2 particles)
	RAB = RAB0*SIG(NPART + 1) 
	X(NPART+1)= RAB/3.464102 
	Y(NPART+1)= X(NPART+1) 
	Z(NPART+1)= X(NPART+1) 
	X(NPART+2)=-X(NPART+1) 
	Y(NPART+2)=-X(NPART+1) 
	Z(NPART+2)=-X(NPART+1)

c     ---read/generate initial configuration
      IF (ibeg.EQ.0) THEN
         WRITE (*, *) 'generate randomly initial configuration'
         CALL CONFIG
         OPEN(64,FILE=INICONF)
         WRITE (64, *) NPART, Dr
	 ns0=0
         Enert0=0.0
         MF0=0
         MFHS0=0
         MFCC0=0
         phcc0=0
         WRITE (64,*) ns0, Enert0, MF0, MFHS0, MFCC0,phcc0
         DO i = 1, NPART
         WRITE (64,  5) sig(i), q(i), X(i), Y(i), Z(i)
         END DO
         CLOSE(64)

      ELSEIF (ibeg.EQ.1) THEN
         WRITE (*, *) ' read conf from disk '
	 OPEN(64,FILE=INICONF)
	 READ (64, *) NPART, Dr
         READ(64,*) ns0, Enert0, MF0, MFHS0, MFCC0,phcc0
         ns0=0
         Enert0=0.0
         MF0=0
         MFHS0=0
         MFCC0=0
         phcc0=0
         DO i = 1, NPART
         READ (64, *) sig(i), q(i), X(i), Y(i), Z(i)
         END DO
	 CLOSE(64)
      ELSEIF (ibeg.EQ.2) THEN
	 WRITE(*,*) 'restart sampling'
         OPEN(64,FILE=INICONF)
         READ (64, *) NPART, Dr
         READ(64,*) ns0, Enert0, MF0, MFHS0, MFCC0,phcc0
         DO i = 1, NPART
         READ (64, *) sig(i), q(i), X(i), Y(i), Z(i)
         END DO
         CLOSE(64)

      END IF

c     ---write input data
      WRITE (*,  1) Equil, Prod, Nsamp
      WRITE (*,  2) Ndispl, Dr,RAB0
      WRITE (*,  3) NPART, BOX*lB, Debye_L, LB,IS
      WRITE (*, *) '#IONS    DIAMETER   CHARGE  CONCENTRATION(M)' 
	DO i=1,NION
      WRITE (*,  4) NPARTI(i), sigi(i), qi(i), CM(i)
	ENDDO 
      RETURN
 1    FORMAT ('  Number of equilibration cycles             :', i10, /, 
     &        '  Number of production cycles                :', i10, /, 
     &        '  Sample frequency                           :'i10, /)
 2    FORMAT ('  Number of att. to displ. a part. per cycle :', i10, /, 
     &        '  Maximum displacement                       :',f10.3,/, 
     &        '  Macroion separation distance:              :',f10.2
     &        //)
 3    FORMAT ('  Total Number of particles                  :', i10, /, 
     &        '  Box length,A                 :', f10.3/,
     &        '  Debye length,A               :', f10.3/,
     &        '  Bjerrum length,A             :', f10.6/,
     &        '  Ionic Strength, M            :', f10.6/)

 4    FORMAT (2x,I5,2x,f8.1,2x,f6.1,2x,e12.6)
 5    FORMAT(2X,6F12.6)
      END

      SUBROUTINE CONFIG
c
c     place `npart' particles on a simple cubic lattice
C
      IMPLICIT NONE
      INTEGER I, K, npart, npmax
      PARAMETER (npmax=5000)
      REAL x(npmax), y(npmax), z(npmax)
      REAL box, hbox, RAB
      REAL sig(npmax), q(npmax), pi
      REAL XX,XD, YD, ZD, RIK
      real ranf
      COMMON /conf1/ x, y, z, npart
      COMMON /sys1/ box, hbox, RAB
      COMMON /pot1/  pi,sig,q

      DO 1 I=NPART,1,-1
3     XX=ranf()
      X(I)=HBOX*(XX-0.5)*2.
      XX=ranf()
      Y(I)=HBOX*(XX-0.5)*2.
      XX=ranf()
      Z(I)=HBOX*(XX-0.5)*2.

      DO 2 K=I+1,NPART+2
      XD=X(I)-X(K)
      IF(XD.GT.HBOX) XD=XD-BOX
      IF(XD.LT.-HBOX) XD=XD+BOX
      YD=Y(I)-Y(K)
      IF(YD.GT.HBOX) YD=YD-BOX
      IF(YD.LT.-HBOX) YD=YD+BOX
      ZD=Z(I)-Z(K)
      IF(ZD.GT.HBOX) ZD=ZD-BOX
      IF(ZD.LT.-HBOX) ZD=ZD+BOX
      RIK=SQRT(XD*XD+YD*YD+ZD*ZD)
      IF (RIK .LT. (SIG(I)+SIG(K))/2. ) GO TO 3
2     CONTINUE
1     CONTINUE
      RETURN
      END


 
      SUBROUTINE SAMPLE(myid,ibeg,IL, Switch, Ener,Dr,deltar,stream)
c
c    Switch (input) = 0: initialize varibales 
c                   = 1: sample averages
c                   = 2: return average results
c
c  Ener (input) : total energy/kT
c  Enert        : average total energy/NkT
c  ns           : number of sample points
c  F            : force between macroions

      IMPLICIT NONE
      INTEGER IL, npmax, ns, Switch, NIONmax, NHIsmax
      PARAMETER (npmax=5000, NIONmax=10, NHIsmax=200)
      INTEGER  NION, NPARTI(NIONmax),ibeg,stream
      INTEGER  I, II, J, K, M, LOCAT, LOCATJ, nhcc
      INTEGER  NPART,myid 
      REAL     Ener, Enert, F, MF, FHS, FCC, MFHS, MFCC
      REAL     Enert1, MF1, MFHS1, MFCC1,phcc,deltar
      REAL     DNEW(npmax), DOLD(npmax),XIL,YIL,ZIL
      REAL     sig(npmax), q(npmax), pi, lB
      REAL     box, hbox, SHBOX, delr, delri, dV, rhoj
      REAL     x(npmax), y(npmax), z(npmax), RAB,sqrt3,RAB0
      REAL     signi,cosi,dist,Dr
      real Enert0,MF0,MFHS0,MFCC0,phcc0
      integer ns0
        REAL F2(4),FHS2(4),FCC2(4)
        INTEGER ncx(4),ncoll

      PARAMETER (sqrt3=1.7320508076)
	character filename(20)*5,filename2(20)*5
	character OUTPUT1*5,INICONF*5
        data filename/'OUT1','OUT2','OUT3','OUT4','OUT5','OUT6'
     &      ,'OUT7','OUT8','OUT9','OUT10','OUT11','OUT12'
     &      ,'OUT13','OUT14','OUT15','OUT16','OUT17','OUT18'
     &      ,'OUT19','OUT20'/
        data filename2/'INP1','INP2','INP3','INP4','INP5','INP6'
     &      ,'INP7','INP8','INP9','INP10','INP11','INP12'
     &      ,'INP13','INP14','INP15','INP16','INP17','INP18'
     &      ,'INP19','INP20'/

c
      COMMON /sys1/  box, hbox, RAB
      COMMON /pot1/  pi, sig, q
      COMMON /DISTAN/ DOLD, DNEW,XIL,YIL,ZIL
      COMMON /samm/ NION, NPARTI
      COMMON /sam/ SHBOX, lB
      COMMON /conf1/ x, y, z, npart
      COMMON /restart/ ns0,Enert0,MF0,MFHS0,MFCC0,phcc0
      SAVE ns, Enert, MF, MFCC,MFHS,nhcc,phcc,RAB0
	OUTPUT1=filename(myid+1)
	INICONF=filename2(myid+1)
c        ---Initialize
      IF (Switch.EQ.0) THEN
		ns=ns0*10**6
		nhcc=INT( float(ns*4)*phcc0 )
c        ---Total energy
		Enert=float(ns)*Enert0
c        ---Mean force between macroions
		MF = float(ns*4)*MF0
		MFHS= float(ns*4)*MFHS0
		MFCC= float(ns*4)*MFCC0

        RAB0= RAB/SIG(NPART + 1)

	ELSE IF (Switch.EQ.1) THEN
		ns=ns+1
c  ---sample excess energy 
		Enert = Enert + Ener

c  ---sample mean force between macroions
	CALL FORCE( DELTAR,NPART+1,NPART+2,ncx(1),F2(1),FHS2(1),FCC2(1)) 
        CALL FORCE(-DELTAR,NPART+1,NPART+2,ncx(2),F2(2),FHS2(2),FCC2(2))
        CALL FORCE( DELTAR,NPART+2,NPART+1,ncx(3),F2(3),FHS2(3),FCC2(3))
        CALL FORCE(-DELTAR,NPART+2,NPART+1,ncx(4),F2(4),FHS2(4),FCC2(4))
        F=0.
        FHS=0.
        FCC=0.
        ncoll=0
        do i=1,4
        F=F+F2(i)
        FHS=FHS+FHS2(i)
        FCC=FCC+FCC2(i)
        ncoll=ncoll+ncx(i)
        enddo

		MF = MF + F
		MFHS=MFHS+FHS
		MFCC=MFCC+FCC
                nhcc=nhcc+ncoll

               IF (MOD(ns,100000).EQ.0) THEN
                phcc=FLOAT(nhcc)/FLOAT(ns*4)
                Enert1 = Enert/FLOAT(ns)
                MF1    = MF/FLOAT(ns*4)
                MFHS1  = MFHS/FLOAT(ns*4)
                MFCC1  = MFCC/FLOAT(ns*4)
                OPEN(67,FILE=OUTPUT1,POSITION='APPEND')
                WRITE(67,223) RAB0, ns/10**6, Enert1,
     &                        MF1, MFHS1, MFCC1,phcc,deltar
                CLOSE(67)
223             FORMAT(1x,f5.2,I5,1x,4f10.2,1x,f5.2,1x,f6.4)


	        OPEN(69,file=INICONF)
   		WRITE (69, *) NPART, Dr
                WRITE(69,222) ns/10**6, Enert1, MF1, MFHS1, MFCC1,phcc
222             FORMAT(1x,I5,1x,4E14.6,1x,f6.3)
      		DO i = 1, NPART+2
      		WRITE (69, 10) sig(i), q(i), X(i), Y(i), Z(i)
      		END DO
      		CLOSE(69)
10    		FORMAT(2X,6(F11.6,1X))
		ENDIF

      ELSE
                phcc=FLOAT(nhcc)/FLOAT(ns*4)
                Enert1 = Enert/FLOAT(ns)
                MF1    = MF/FLOAT(ns*4)
                MFHS1  = MFHS/FLOAT(ns*4)
                MFCC1  = MFCC/FLOAT(ns*4)


      OPEN(69,file=INICONF)
      WRITE (69, *) NPART, Dr
      WRITE(69,222) ns/10**6, Enert1, MF1, MFHS1, MFCC1,phcc
      DO i = 1, NPART+2
      WRITE (69, 20) sig(i), q(i), X(i), Y(i), Z(i)
      END DO
      CLOSE(69)
20    FORMAT(2X,6(F11.6,1X))

                OPEN(77,FILE='RESULTS',POSITION='APPEND')
		write(77,30) RAB0,Enert1,MF1,MFHS1,MFCC1
		CLOSE(77)
30             FORMAT(1x,F4.2,1x,4E14.6)

      END IF
      RETURN
      END


      SUBROUTINE TOTERG(GW, Ener)
C *******************************************************************
c
c  Ener (output) : total energy/kT
c  GW: Gaussian Width parameter, adjust for efficiency of Ewald Sum
c
C *******************************************************************

      IMPLICIT NONE
      INTEGER npmax,npart
      PARAMETER (npmax=5000)

      REAL     x(npmax),y(npmax),z(npmax)
      REAL     box,hbox
      REAL     pi, sig(npmax), q(npmax) 
	REAL     VR, VK, Ener, GW, RAB

      COMMON /sys1/  box, hbox, RAB
      COMMON /pot1/  pi,sig,q
      COMMON /conf1/ x, y, z, npart

	CALL SETUP(GW) 
	CALL RWALD ( npart+2, VR, X, Y, Z, Q, BOX)
      call KWALD ( npart+2, VK, X, Y, Z, Q, BOX )

	Ener = VR + VK
      RETURN
      END


