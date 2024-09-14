       Program  Main
       implicit none
#include <mpif.h>
       integer i, myid, nprocs, ierror, junk

C Initilize parallel mode
       call MPI_INIT(ierror)
       call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
       call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
       call pforce(myid)
c close parallel mode
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
      REAL     box,temp,hbox
      REAL     dro, frac, Dr
      SAVE naccp, attempp
      COMMON /sys1/ box,hbox,temp
 
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
*********************************************************************
        IMPLICIT NONE
        INTEGER     MAXK, TOTK, KMAX, KX, KY, KZ, I, KSQMAX, KSQ, 
     &              IL, npmax, npart 
        PARAMETER   ( MAXK = 1000, KMAX = 5, KSQMAX = 27, npmax=5000)
        REAL         TWOPI, TWOPL, FACTOR, VD, VDN, RSQPI
        REAL       KVEC(MAXK), KAPPA, DEK,DER
        REAL       box, temp, hbox
        REAL       sig(npmax), q(npmax), pi
        REAL       x(npmax), y(npmax), z(npmax), XN, YN, ZN
        REAL       E1, E2, XIL, YIL, ZIL, RXIJ, RYIJ, RZIJ
	REAL       RIJSQ, RIJ, KRIJ
	REAL       ERFC
        PARAMETER   ( TWOPI = 6.2831853, RSQPI = 0.5641896 )

        COMPLEX   EIKX(0:KMAX), EIKY(-KMAX:KMAX), EIKZ(-KMAX:KMAX)
        COMPLEX   DEIKX(0:KMAX),DEIKY(-KMAX:KMAX),DEIKZ(-KMAX:KMAX)
        COMPLEX   DEIKR, SUMI(MAXK), DSUMI(MAXK)

        COMMON / BLOCK2 / KVEC, KAPPA
	COMMON / MOVE /   SUMI, VD, DSUMI, VDN
        COMMON / conf1 /  x, y, z, npart
        COMMON / pot1 /   pi,sig,q
        COMMON / sys1 /   box,hbox,temp
C    ** CALCULATE ENERGY CHANGE IN R-SPACE
		E1=0.
		E2=0.

		XIL=X(IL)
		YIL=Y(IL)
		ZIL=Z(IL)

		DO 5 I=1,NPART+3
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
		KRIJ=KAPPA*RIJ/BOX
	        E1=E1 + Q(I)*Q(IL)*ERFC(REAL(KRIJ))/RIJ
5		CONTINUE

C IL-th PARTICLE AT NEW POSITION
		DO 6 I=1,NPART+3
		IF(I.EQ.IL) GOTO 6
		RXIJ = X(I) - XN
		RYIJ = Y(I) - YN
		RZIJ = Z(I) - ZN
                RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
                RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
                RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

	        RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
		RIJ=SQRT(RIJSQ)

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

	FUNCTION Dielec(T)
c This function calculates the dielectric constant of pure water at different
c temperatures (T in Kelvin)
	REAL     Dielec, T, T0,Tr,tt
	REAL     A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,denw0,denw,denr
      PARAMETER (A1=7.62571E0,A2=2.44003E2,A3=-1.40569E2,A4=2.77841E1,
     &  A5=-9.62805E1,A6=4.17909E1,A7=-1.02099E1,A8=-4.52059E1,
     &  A9=8.46395E1,A10=-3.58644E1,T0=298.15,denw0=1000.0E0)   

	tt=T-273.15E0
      denw=999.65E0+2.0438e-1*tt-6.174e-2*tt**1.5
      Tr=T/T0
      denr=denw/denw0
      Dielec=1.0E0+(A1/Tr)*denr+(A2/Tr+A3+A4*Tr)*denr**2
     &       +(A5/Tr+A6*Tr+A7*Tr**2)*denr**3
     &       +(A8/Tr**2+A9/Tr+A10)*denr**4  
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


	  SUBROUTINE FORCE(itag,IL,DELTAR,ncoll,F, FHS, FCC)
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
     &               npmax, npart, nhcc, ncoll ,IL,IA,itag
        PARAMETER   ( MAXK = 1000, KMAX = 5, KSQMAX = 27, npmax=5000)
        REAL       TWOPI, TWOPL, FACTOR, VD, VDN, RSQPI
        REAL       KVEC(MAXK), KAPPA, DEK,DER
        real       ranf0
	REAL       box, temp, hbox, RAB, XX, DELTAR, F, FCC, FHS
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
        COMMON / sys1 /   box,hbox,temp
c     ---make a small displacement for one of the macroions 

c New position
c the sign of DELTAR determines which direction moves
	if(IL.eq.NPART+1) then
		if(itag.gt.0) then
         	XN  = X(IL) + DELTAR
         	YN  = Y(IL)
         	ZN  = Z(IL)  
		else
	 	XN  = X(IL) + DELTAR*0.5
         	YN  = Y(IL) - DELTAR*0.8660254
         	ZN  = Z(IL)
		endif
	elseif (IL.eq.NPART+2) then
                if(itag.gt.0) then
                XN  = X(IL) - DELTAR
                YN  = Y(IL)
                ZN  = Z(IL) 
                else
                XN  = X(IL) - DELTAR*0.5
                YN  = Y(IL) - DELTAR*0.8660254
                ZN  = Z(IL)
                endif
        elseif (IL.eq.NPART+3) then
                if(itag.gt.0) then
                XN  = X(IL) + DELTAR*0.5
                YN  = Y(IL) + DELTAR*0.8660254
                ZN  = Z(IL) 
                else
                XN  = X(IL) - DELTAR*0.5
                YN  = Y(IL) + DELTAR*0.8660254
                ZN  = Z(IL)
                endif
	endif

      XIL = X(IL)
      YIL = Y(IL)
      ZIL = Z(IL)

C    ** CALCULATE ENERGY CHANGE IN R-SPACE
                E1=0.
                E2=0.

C IL_th PARTICLE AT OLD POSITION
                DO 5 I=1,NPART+3
		if(IL.eq.I) goto 5
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

C IL_th PARTICLE AT NEW POSITION
                ncoll=0
                DO 6 I=1,NPART+3
		if(IL.eq.I) goto 6
                        RXIJ = X(I) - XN
                        RYIJ = Y(I) - YN
                        RZIJ = Z(I) - ZN

              RXIJ = RXIJ - BOX * ANINT ( RXIJ / BOX )
              RYIJ = RYIJ - BOX * ANINT ( RYIJ / BOX )
              RZIJ = RZIJ - BOX * ANINT ( RZIJ / BOX )

                RIJSQ=RXIJ*RXIJ + RYIJ*RYIJ + RZIJ*RZIJ
                        RIJ=SQRT(RIJSQ)

                        IF (RIJ .LT. (SIG(IL)+SIG(I))/2.0 )
     &                           ncoll=ncoll+1

                        KRIJ=KAPPA*RIJ/BOX
                E2=E2 + Q(I)*Q(IL)*erfc(KRIJ)/RIJ
6               CONTINUE

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

 
	subroutine pforce(myid)
c________________________________________________________________________
c
c          NVT_MC FOR INTERACTION BETWEEN TWO PARTICLES IN 
c                     AN ELECTROLYTE SOLUTION
c__________________________________________________________________________
 
      IMPLICIT NONE
      INTEGER equil, prod, nsamp, ii, icycl, ndispl, attempt, 
     &  iseed,nacc, ncycl, nmoves, imove, IL,ibeg,nstart,myid,iseeds(10)
      REAL Ener, Enert, Dr, GW,deltar,ranf0
c     ---initialize sysem
        data iseeds/89743,76234,89734,876234,98721,63283,87324,39734,
     &        763425,653245/
        EXTERNAL RANSET,ranf0
        iseed=iseeds(myid)
        CALL RANSET(Iseed)
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
		CALL SAMPLE(myid,ibeg,IL, 0, Ener,Dr,deltar,iseed)
         attempt = 0
         nacc = 0

c   ---intialize the subroutine that adjust the maximum displacement
         CALL ADJUST(attempt, nacc, dr)

         DO 50 icycl = 1, ncycl

c              ---attempt to displace a particle
             DO imove = 1, nmoves
             CALL MCMOVE(myid,IL, Ener, attempt, nacc, dr,iseed)
             ENDDO

c              ---sample averages
               IF (ii.EQ.2) THEN
               CALL SAMPLE(myid,ibeg,IL, 1, Ener,Dr,deltar,iseed)
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
	CALL SAMPLE(myid,ibeg,IL, 2, Ener,Dr,deltar,iseed)
	WRITE(*,*) 'FINAL AVERAGE ENERGY/NkT = ', Ener      
      STOP
 
99001 FORMAT (' Total energy at initial configuration: ', f12.5, /)
99002 FORMAT (' Total energy at the end of simulation: ', f12.5, /, 
     &        '       running energy              : ', f12.5, /, 
     &        '       difference                  :  ', e12.5, /) 
99003 FORMAT (' Number of att. to displ. a part.  : ', i10, /, 
     &        ' success: ', i10, '(= ', f5.2, '%)')
      END

      SUBROUTINE MCMOVE(myid,IL, Ener, Attempt, Nacc, Dr,iseed)
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
      INTEGER IL, Attempt, Nacc,  npmax, npart, MAXK, I,iseed,myid
      PARAMETER (npmax=5000, MAXK = 1000)
      REAL     box,temp,hbox
      REAL     x(npmax), y(npmax), z(npmax)
      REAL     Ener, xn, yn, zn, DER, DEK, Dr
      real ranf0,xtest
	REAL     VD, VDN 
        REAL     sig(npmax), q(npmax), pi
      COMMON /pot1/  pi,sig,q
      COMPLEX   SUMI(MAXK),DSUMI(MAXK)
        COMMON / sys1 / box, hbox, temp
      COMMON / conf1 /  x, y, z, npart
	COMMON / MOVE /   SUMI, VD, DSUMI, VDN
C    *******************************************************************
	Attempt = Attempt + 1
c     ---select a particle at random
      IL = INT( float(NPART)*ranf0(iseed) ) + 1

c     ---give the particle a random displacement
      xn = X(IL) + (ranf0(iseed) -0.5E0)*Dr
      yn = Y(IL) + (ranf0(iseed) -0.5E0)*Dr
      zn = Z(IL) + (ranf0(iseed) -0.5E0)*Dr

c     ---calculate energy change
	CALL DELTAE( IL, XN, YN, ZN, DER, DEK )
c     ---acceptance test
	xtest=DER+DEK
	if (xtest.gt.20.) then
	xtest=0.0
	else
	xtest=EXP(-(DER+DEK))
	endif

      IF (ranf0(iseed) .LT. xtest) THEN
c        --accepted
		Nacc = Nacc + 1
			
		VD=VDN
		DO I=1,MAXK
		SUMI(I) = DSUMI(I)
		ENDDO
		Ener = Ener + DER + DEK

c        ---put the particle in the simulation box
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
c    TEMP       : temperature, Kelvin
c
c    NION   = number of ion species
c    NIONmax= maxmum number of ion species
c    sig(i) = diameters of ions and macroions/Bjerrum length
c    sigi(i)= diameters of ions and macroions species, in Angstroms
c    q(i)   = valence for an ion
c    qi(i)  = valence for an ion, for input
c    Dielec(T) = dielectric constant of the medium 
c                at temperature T (external function)
c    CM(i)  = molar concentration of ions
c    lB0    = Bjerrum length for water at 298.15 K, in Angstroms
c    lB     = Bjerrum length at system temperature, in Angstroms
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
      REAL      pi, box, temp, hbox,deltar,
     &                 con1, con2, Dr, lB0, lB, GW, SHBOX
      INTEGER ibeg,Equil,Prod,i,j,k,Ndispl,Nsamp,
     &        NION,NIONmax,npmax,npart,NHIsmax
      PARAMETER (npmax=5000, NIONmax=10, lB0=7.144464E0,
     &	  con1=2.338874e4, con2=6.022137e-4,NHIsmax=200)
      REAL     sig(npmax), q(npmax), ZERO, Dielec, Debye_L
      REAL     IS, sigi(NIONmax), qi(NIONmax), CM(NIONmax)
      REAL     x(npmax), y(npmax), z(npmax), RAB,RAB0(10)
        real Enert0,MF0,MFHS0,MFCC0,phcc0
        integer ns0,nnode
      INTEGER  NPARTI(NIONmax),myid
      COMMON /sys1/  box,hbox,temp
      COMMON /pot1/  pi,sig,q
      COMMON /conf1/ x, y, z, npart
      COMMON /samm/ NION, NPARTI
      COMMON /sam/ SHBOX, lB
      EXTERNAL Dielec
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
      READ (65, *) (RAB0(i),i=1,nnode)
c     ---read solution condition and model parameters
      READ (65, *)
      READ (65, *) TEMP
      READ (65, *)
      READ (65, *) NION
      READ (65, *)
      DO i=1,NION
      READ (65, *) sigi(i), qi(i), NPARTI(i)
      ENDDO
      CLOSE(65)

	PI = 3.14159265
c check charge neutrality at initial condition
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
      IF (NPART.GT.NPMax) THEN
         WRITE (*, *) ' ERROR: number of particles too large'
         STOP
      END IF

c     ---calculate parameters:
	lB=lB0*con1/(Dielec(TEMP)*TEMP)

c     --- normalized with Bjerrum length
	BOX = BOX/lB
        HBOX = BOX/2.E0
	SHBOX = HBOX*HBOX

	DO i=1,NION
	CM(I)=FLOAT(NPARTI(I))/((BOX*lB)**3*CON2)
	ENDDO

	IS=0.0
	do i=1,NION-1
	IS=IS+qi(i)**2*CM(I)/2.
        enddo
	Debye_L=1./sqrt(8.*PI*lB*IS*CON2)

	k=0
	DO 10 i = 1,NION
	DO j = 1, NPARTI(i)
	k=k+1
	sig(k) = sigi(i)/lB
	q(k) = qi(i)
	ENDDO
10	CONTINUE

c   --- assign LOCATION of three macroions ( NPART+1, NPART+2,NPART+3)
	RAB=RAB0(myid+1)*SIG(NPART+1)/2.
	X(NPART+1)= RAB 
	Y(NPART+1)=-RAB*0.57735027 
	Z(NPART+1)= 0.0
	X(NPART+2)=-RAB 
	Y(NPART+2)=-RAB*0.57735027 
	Z(NPART+2)= 0.0
        X(NPART+3)= 0.0
        Y(NPART+3)= RAB*1.15470054
        Z(NPART+3)=-0.0

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
      WRITE (*,  2) Ndispl, Dr,RAB0(myid+1)
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
      REAL box, temp, hbox
      REAL sig(npmax), q(npmax), pi
      REAL XX,XD, YD, ZD, RIK
      real ranf
      COMMON /conf1/ x, y, z, npart
      COMMON /sys1/ box,hbox,temp
      COMMON /pot1/  pi,sig,q

      DO 1 I=NPART,1,-1
3     XX=ranf()
      X(I)=HBOX*(XX-0.5)*2.
      XX=ranf()
      Y(I)=HBOX*(XX-0.5)*2.
      XX=ranf()
      Z(I)=HBOX*(XX-0.5)*2.

      DO 2 K=I+1,NPART+3
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


 
      SUBROUTINE SAMPLE(myid,ibeg,IL, Switch, Ener,Dr,deltar,iseed)
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
      INTEGER  NION, NPARTI(NIONmax),ibeg,iseed
      INTEGER  I, II, J, K, M, LOCAT, LOCATJ, nhcc
      INTEGER  NPART,myid 
      REAL     Ener, Enert, F, MF, FHS, FCC, MFHS, MFCC
      REAL     Enert1, MF1, MFHS1, MFCC1,phcc,deltar
      REAL     sig(npmax), q(npmax), pi, lB
      REAL     box, temp, hbox, SHBOX, delr, delri, dV, rhoj
      REAL     x(npmax), y(npmax), z(npmax), RAB,sqrt3,RAB0(10)
      REAL     signi,cosi,dist,Dr
      real Enert0,MF0,MFHS0,MFCC0,phcc0
      real F2(12),FHS2(12),FCC2(12)
      integer ns0, ncx(12), ncoll
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
      COMMON /sys1/  box, hbox, temp
      COMMON /pot1/  pi, sig, q
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
		nhcc=INT( float(ns)*phcc0 )
c        ---Total energy
		Enert=float(ns)*Enert0
c        ---Mean force between macroions
		MF = float(ns)*MF0
		MFHS= float(ns)*MFHS0
		MFCC= float(ns)*MFCC0

        RAB=(X(NPART+1)-X(NPART+2))**2 + (Y(NPART+1)-Y(NPART+2))**2
        RAB0(myid+1) = sqrt(RAB)/SIG(NPART + 1)

	ELSE IF (Switch.EQ.1) THEN
		ns=ns+1
c  ---sample excess energy 
		Enert = Enert + Ener

c  ---sample mean force between macroions
	CALL FORCE( 1,NPART+1, deltar,ncx(1),F2(1),FHS2(1),FCC2(1))
        CALL FORCE( 1,NPART+1,-deltar,ncx(2),F2(2),FHS2(2),FCC2(2))
        CALL FORCE( 1,NPART+2, deltar,ncx(3),F2(3),FHS2(3),FCC2(3))
        CALL FORCE( 1,NPART+2,-deltar,ncx(4),F2(4),FHS2(4),FCC2(4))
        CALL FORCE( 1,NPART+3, deltar,ncx(5),F2(5),FHS2(5),FCC2(5))
        CALL FORCE( 1,NPART+3,-deltar,ncx(6),F2(6),FHS2(6),FCC2(6))

        CALL FORCE(-1,NPART+1, deltar,ncx(7),F2(7),FHS2(7),FCC2(7))
        CALL FORCE(-1,NPART+1,-deltar,ncx(8),F2(8),FHS2(8),FCC2(8))
        CALL FORCE(-1,NPART+2, deltar,ncx(9),F2(9),FHS2(9),FCC2(9))
        CALL FORCE(-1,NPART+2,-deltar,ncx(10),F2(10),FHS2(10),FCC2(10))
        CALL FORCE(-1,NPART+3, deltar,ncx(11),F2(11),FHS2(11),FCC2(11))
        CALL FORCE(-1,NPART+3,-deltar,ncx(12),F2(12),FHS2(12),FCC2(12))

	F=0.
	FHS=0.
	FCC=0.
	ncoll=0
	do i=1,12
	F=F+F2(i)
	FHS=FHS+FHS2(i)
	FCC=FCC+FCC2(i)
	ncoll=ncoll+ncx(i)
	enddo
	F=F/12.
	FHS=FHS/12.
	FCC=FCC/12.
	ncoll=ncoll/12.

		MF = MF + F
		MFHS=MFHS+FHS
		MFCC=MFCC+FCC
		nhcc=nhcc+ncoll

               IF (MOD(ns,100000).EQ.0) THEN
                phcc=FLOAT(nhcc)/FLOAT(ns)
                Enert1 = Enert/FLOAT(ns)
                MF1    = MF/FLOAT(ns)
                MFHS1  = MFHS/FLOAT(ns)
                MFCC1  = MFCC/FLOAT(ns)
                OPEN(67,FILE=OUTPUT1,POSITION='APPEND')
                WRITE(67,222) RAB0(myid+1), ns/10**6, Enert1, 
     &                        MF1, MFHS1, MFCC1,phcc,deltar
                CLOSE(67)
222             FORMAT(1x,f5.2,2x,I5,2x,4f10.2,2x,f5.2,2x,f6.4)

	        OPEN(69,file=INICONF)
   		WRITE (69, *) NPART, Dr
                WRITE(69,223) ns/10**6, Enert1, MF1, MFHS1, MFCC1,phcc
223             FORMAT(1x,I5,1x,4E14.6,1x,f6.3)
      		DO i = 1, NPART
      		WRITE (69, 10) sig(i), q(i), X(i), Y(i), Z(i)
      		END DO
      		CLOSE(69)
10    		FORMAT(2X,6(F11.6,1X))
		ENDIF

      ELSE
                phcc=FLOAT(nhcc)/FLOAT(ns)
                Enert1 = Enert/FLOAT(ns)
                MF1    = MF/FLOAT(ns)
                MFHS1  = MFHS/FLOAT(ns)
                MFCC1  = MFCC/FLOAT(ns)


      OPEN(69,file=INICONF)
      WRITE (69, *) NPART, Dr
      WRITE(69,223) ns/10**6, Enert1, MF1, MFHS1, MFCC1,phcc
      DO i = 1, NPART
      WRITE (69, 20) sig(i), q(i), X(i), Y(i), Z(i)
      END DO
      CLOSE(69)
20    FORMAT(2X,6(F11.6,1X))

                OPEN(77,FILE='RESULTS',POSITION='APPEND')
		write(77,30) RAB0(myid+1),Enert1,MF1,MFHS1,MFCC1
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
      REAL     box,temp,hbox
      REAL     pi, sig(npmax), q(npmax) 
	REAL     VR, VK, Ener, GW

      COMMON /sys1/  box,hbox,temp
      COMMON /pot1/  pi,sig,q
      COMMON /conf1/ x, y, z, npart

	CALL SETUP(GW) 
	CALL RWALD ( npart+3, VR, X, Y, Z, Q, BOX)
        call KWALD ( npart+3, VK, X, Y, Z, Q, BOX )
	Ener = VR + VK
      RETURN
      END


      SUBROUTINE RANSET(Iseed)
C  --- initializes random number generator
      IMPLICIT NONE
      INTEGER Iseed
        EXTERNAL RSTART
      CALL RSTART(Iseed)
      RETURN
      END


      SUBROUTINE RSTART(Iseeda)
C----------------------------------------------------------------------C
C       Initialize Marsaglia list of 24 random numbers.
C----------------------------------------------------------------------C
      IMPLICIT NONE
      real CARRY, ran, RANDX, SEED
      INTEGER i, I24, ISEED, Iseeda, J24
      COMMON /RANDOM/ SEED(24), CARRY, I24, J24, ISEED

      I24 = 24
      J24 = 10
      CARRY = 0.E+0
      ISEED = Iseeda
c
c       get rid of initial correlations in rand by throwing
c       away the first 100 random numbers generated.
c
      DO i = 1, 100
         ran = RANDX(ISEED)
      END DO
c
c       initialize the 24 elements of seed
c

      DO i = 1, 24
         SEED(i) = RANDX(ISEED)
      END DO

      RETURN
      END


      FUNCTION RANDX(Iseed)
C----------------------------------------------------------------------C
C  Random number generator, fast and rough, machine independent.
C  Returns an uniformly distributed deviate in the 0 to 1 interval.
C  This random number generator is portable, machine-independent and
C  reproducible, for any machine with at least 32 bits / real number.
C  REF: Press, Flannery, Teukolsky, Vetterling, Numerical Recipes (1986)
C----------------------------------------------------------------------C
      IMPLICIT NONE
      INTEGER IA, IC, Iseed, M1
      real RANDX, RM
      PARAMETER (M1=714025, IA=1366, IC=150889, RM=1.E+0/M1)
c
      Iseed = MOD(IA*Iseed+IC, M1)
      RANDX = Iseed*RM
      IF (RANDX.LT.0.E+0) THEN
         STOP '*** Random number is negative ***'
      END IF
      RETURN
      END

      FUNCTION ranf0(Idum)
c
c     random number generator
c     Idum (input): can be used as seed (not used in present
c                   random number generator.

      IMPLICIT NONE
      INTEGER Idum
      real ranf0, RCARRY
      ranf0 = RCARRY()
        Idum=Idum
      RETURN
C ----------------------------------------------------C
      END

      FUNCTION RCARRY()
C----------------------------------------------------------------------C
C       Random number generator from Marsaglia.
C----------------------------------------------------------------------C
      IMPLICIT NONE
      real CARRY, RCARRY, SEED, TWOm24, TWOp24, uni
      INTEGER I24, ISEED, J24
      PARAMETER (TWOp24=16777216.E+0, TWOm24=1.E+0/TWOp24)
      COMMON /RANDOM/ SEED(24), CARRY, I24, J24, ISEED
c
c       F. James Comp. Phys. Comm. 60, 329  (1990)
c       algorithm by G. Marsaglia and A. Zaman
c       base b = 2**24  lags r=24 and s=10
c
      uni = SEED(I24) - SEED(J24) - CARRY
      IF (uni.LT.0.E+0) THEN
         uni = uni + 1.E+0
         CARRY = TWOm24
      ELSE
         CARRY = 0.E+0
      END IF
      SEED(I24) = uni
      I24 = I24 - 1
      IF (I24.EQ.0) I24 = 24
      J24 = J24 - 1
      IF (J24.EQ.0) J24 = 24
      RCARRY = uni
      RETURN
      END

