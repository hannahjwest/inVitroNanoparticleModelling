*     MODIFIED 9/8/01 TO CHANGE PROLIFERATION TERM
*     (FROM fullflux.f)
*     OPTIMAL PROLIF WAS AT N=1 SHOULD BE AT N=N_inf
*
*     Markus Owen, 11/6/01
*     To simulate model for tumour growth and macrophage infiltration
*     Substitution of velocity gives system of parabolic PDEs
*
*     Modified from MPHIPHF to include multiple macrophage population
*
*     Modified from:
*     D03PHF Example Program Text
*     Mark 16 Revised. NAG Copyright 1993.

*     .. Parameters ..
      INTEGER NOUT,NOUT2
      PARAMETER (NIN=4,NOUT=10,NOUT1=11,NOUT2=12)
      INTEGER NPDE, NPTS, NCODE, M, NXI, NEQN, NIW, NWKRES,
     +     LENODE, NW
      PARAMETER (NPDE=3,NPTS=301,NCODE=2,M=2,NXI=2,
     +     NEQN=NPDE*NPTS+NCODE,NIW=24,
     +     NWKRES=NPDE*(NPTS+6*NXI+3*NPDE+15)
     +     +NCODE+NXI+7*NPTS+2,LENODE=11*NEQN+50,
     +     NW=NEQN*NEQN+NEQN+NWKRES+LENODE)

*..   Local Scalars ..
      DOUBLE PRECISION TOUT, NEC, NRR, RCDM, RCP, LAMBDA
      DOUBLE PRECISION DUX1,DUX3,DUX4,DUX5,VELC,DUX6
*      DOUBLE PRECISION KBIT,VBIT
      INTEGER I, IFAIL, IND, IT, ITASK, ITOL, ITRACE, ICDM, ICP
      LOGICAL THETA, NONEC, NORCDM, NORCP
      LOGICAL POP1,POP2,POPR1,POPR2,POPK1,POPK2
      CHARACTER LAOPT, NORM

*..   Local Arrays ..
      DOUBLE PRECISION ALGOPT(30), ATOL(1), RTOL(1), U(NEQN),
     +     W(NW), X(NPTS), XI(2)
      INTEGER IW(NIW)

*..   External Subroutines ..
      EXTERNAL BNDARY, D03PHF, ODEDEF, PDEDEF, UVINIT, SOURCE, D01GAF
*..   time
      DOUBLE PRECISION TS
      COMMON /TAXIS/ TS

*..   parameters used in groupings only in main program
      DOUBLE PRECISION BCL1,BCL2,BCM,CAL1,CAL2,CAM,CCL1,CCL2,CCM,CDL1
     $     ,CDL2,CDM,CK1,CK2,CP,CSPAT,K1TMP,K2TMP,L1STMP,L2STMP
     $     ,MINAL1,MINAL2,MINAM,MINDL1,MINDL2,MINDM,MINK1,MINK2
     $     ,USQRTPRT1

*..   parameters for functions returned by SOURCE 
      DOUBLE PRECISION BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG
      COMMON /SRC/ BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG

*..   parameters for spatial terms returned by SOURCE
      DOUBLE PRECISION CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      COMMON /SPATIAL/ CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      
*..   parameters for boundary conditions
      DOUBLE PRECISION HA,HL1,HL2,VINF,LINF2,P0,QP,HV
      COMMON /BCS/ HA,HL1,HL2,VINF,LINF2,P0,QP,HV

*..   parameters for initial data
      DOUBLE PRECISION RINIT
      COMMON /INIT/ RINIT

*..   parameters governing plots, time, time of m-phi introduction...
      INTEGER NPLOTS
      DOUBLE PRECISION TPLOT,TL1,TL2,TKL1,TKL2,TRL1,TRL2


c$$$*..   binding params
c$$$      DOUBLE PRECISION UKA,UKD,UKI,UL,UF0
c$$$      COMMON /BINDING/ UKA,UKD,UKI,UL,UF0

*..   binding params
      DOUBLE PRECISION UKA,UKD,UKI,UL,UF0,URHOL,URHOF
      COMMON /BINDING/ UKA,UKD,UKI,UL,UF0,URHOL,URHOF

*..   arrays of variables for integration
      DOUBLE PRECISION UM(NPTS),UL1(NPTS),UL2(NPTS),FOURPI,
     &     UMNV,UL1NV,UL2NV,UMQV,UL1QV,UL2QV,UMPV,UL1PV,UL2PV,
     &     UMINT,UL1INT,UL2INT,UBB1(NPTS),UBB2(NPTS),USQPRT1(NPTS)
      COMMON /FPI/ FOURPI

      CHARACTER*40 FSTEM, IDATA, ODATA, ODATA1, CAPLIN, ODATA2
      INTEGER ISTEM

      WRITE(6,*) 'Enter file stem:'
      READ(5,'(A40)') FSTEM
      ISTEM = INDEX(FSTEM,' ') - 1
      IDATA = FSTEM(1:ISTEM) // '.pars'

*..   Open parameter file
      OPEN (UNIT=NIN,FILE=IDATA,STATUS='old')
      Write(6,*) 'Opening data file ',IDATA

      FOURPI=16.0d0*atan(1.0d0)
      WRITE(6,*) FOURPI

C     .. Read in Parameters ..   
      READ(NIN,*) MM
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) CP
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) BL1
      READ(NIN,*) BL2
      READ(NIN,*) BM
      READ(NIN,*) MINDL1
      READ(NIN,*) MINDL2
      READ(NIN,*) MINDM
      READ(NIN,*) CDL1
      READ(NIN,*) CDL2
      READ(NIN,*) CDM
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) K1TMP
      READ(NIN,*) K2TMP
      READ(NIN,*) MINK1
      READ(NIN,*) MINK2
      READ(NIN,*) CK1
      READ(NIN,*) CK2
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) BCL1
      READ(NIN,*) BCL2
      READ(NIN,*) BCM
      READ(NIN,*) BCP
      READ(NIN,*) CCL1
      READ(NIN,*) CCL2
      READ(NIN,*) CCM
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) L1A
      READ(NIN,*) L2A
      READ(NIN,*) MA
      READ(NIN,*) MINAL1
      READ(NIN,*) MINAL2
      READ(NIN,*) MINAM
      READ(NIN,*) CAL1
      READ(NIN,*) CAL2
      READ(NIN,*) CAM
      READ(NIN,*) L1AD
      READ(NIN,*) L2AD
      READ(NIN,*) MAD
      READ(NIN,*) DF
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) CHI1
      READ(NIN,*) CHI2
      READ(NIN,*) DMPHI1
      READ(NIN,*) DMPHI2
      READ(NIN,*) DTUM
      READ(NIN,*) MINSP
      READ(NIN,*) CSPAT
      READ(NIN,*) DP
      READ(NIN,*) DV
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) HV
      READ(NIN,*) HL2
      READ(NIN,*) L1STMP
      READ(NIN,*) L2STMP
      READ(NIN,*) P0
      READ(NIN,*) PINF
      READ(NIN,*) PCAP
      READ(NIN,*) QP
      READ(NIN,*) HA
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) LAMBDA
      READ(NIN,*) NPLOTS
      READ(NIN,*) RINIT
      READ(NIN,*) TL1
      READ(NIN,*) TL2
      READ(NIN,*) TKL1
      READ(NIN,*) TKL2
      READ(NIN,*) TRL1
      READ(NIN,*) TRL2
      READ(NIN,*) TPLOT
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*)
      READ(NIN,*) UKA
      READ(NIN,*) UKD
      READ(NIN,*) UKI
      READ(NIN,*) UF0
      READ(NIN,*) UL
      READ(NIN,*) UG
      READ(NIN,*) URHOL
      READ(NIN,*) URHOF


      CLOSE(NIN) 

*..   Open output files
      WRITE(6,*) 'Capped (C, default) or Linear (L):'
      READ(5,'(A1)') CAPLIN
      IF (CAPLIN .EQ. 'L' .OR. CAPLIN .EQ. 'l') THEN
         ODATA = FSTEM(1:ISTEM) // '.splin'
         ODATA1 = FSTEM(1:ISTEM) // '.tlin' 
         ODATA2 = FSTEM(1:ISTEM) // '.vlin'
         PCAP = 1.0d0
      ELSE
         ODATA = FSTEM(1:ISTEM) // '.spcap'
         ODATA1 = FSTEM(1:ISTEM) // '.tcap'
         ODATA2 = FSTEM(1:ISTEM) // '.vcap'
      END IF

      OPEN (UNIT=NOUT,FILE=ODATA,STATUS='replace')
      OPEN (UNIT=NOUT1,FILE=ODATA1,STATUS='replace')
      OPEN (UNIT=NOUT2,FILE=ODATA2,STATUS='replace')
      WRITE(6,*) 'Opening output files ',ODATA,ODATA1
      WRITE(6,*) 'and ',ODATA2

*.. Executable Statements ..
      
*..   scalings for proliferation
      CCPC = CP**MM

*.. scalings for death terms
      CCDM = CDM**MM
      SCALDM = (1.0e0-MINDM/BM)*(CCDM+1.0e0)
      CCDL1 = CDL1**MM
      SCADL1 = (1.0e0-MINDL1/BL1)*(CCDL1+1.0e0)
      CCDL2 = CDL2**MM
      SCADL2 = (1.0e0-MINDL2/BL2)*(CCDL2+1.0e0)
      CCK1 = CK1**MM
      SCALK1 = (1.0e0-MINK1)*(CCK1+1.0e0)
      CCK2 = CK2**MM
      SCALK2 = (1.0e0-MINK2)*(CCK2+1.0e0)

*..   scalings for nutrient consumption
      CCCM = CCM**MM
      SCALCM = BCM*(CCCM+1.0e0)
      CCCL1 = CCL1**MM
      SCACL1 = BCL1*(CCCL1+1.0e0)
      CCCL2 = CCL2**MM
      SCACL2 = BCL2*(CCCL2+1.0e0)

*..   scalings for chemoattractant production
      CCAM = CAM**MM
      SCALAM = (1.0e0-MINAM)*(CCAM+1.0e0)
      CCAL1 = CAL1**MM
      SCAAL1 = (1.0e0-MINAL1)*(CCAL1+1.0e0)
      CCAL2 = CAL2**MM
      SCAAL2 = (1.0e0-MINAL2)*(CCAL2+1.0e0)

*..   scaling for spatial terms
      CCSP = CSPAT**MM
      SCALSP = (1.0e0-MINSP)*(CCSP+1.0e0)

      ITRACE = 0
      ITOL = 1
c      ATOL(1) = 1.0e-4
      ATOL(1) = 1.0e-7
      RTOL(1) = ATOL(1)

*
*     Set break-points
*
      X(1) = 0.0e0
      X(NPTS) = 1.0e0
      DO 20 I = 2, NPTS-1
         IF (LAMBDA .LT. 1.0e0) THEN 
            X(I) = (1-LAMBDA**(I-1))/(1-LAMBDA**(NPTS-1))
         ELSE
            X(I) = (I-1.0e0)/(NPTS-1.0e0)
         END IF
 20   CONTINUE

*
      XI(1) = 0.0e0
      XI(2) = 1.0e0
      NORM = 'A'
      LAOPT = 'F'
      IND = 0
      ITASK = 1

*
*Set THETA to .TRUE. if the Theta integrator is required
*
      THETA = .FALSE.
      DO 40 I = 1, 30
         ALGOPT(I) = 0.0e0
 40   CONTINUE
      IF (THETA) THEN
         ALGOPT(1) = 2.0e0
      ELSE
         ALGOPT(1) = 0.0e0
      END IF

*
*Loop over output value of t
*
      TS = 0.0e0
      TOUT = 0.0e0

      CALL UVINIT(NPDE,NPTS,X,U,NCODE,NEQN)
 
      WRITE (NOUT1,99998) TS,0.0,0.0,0.0,U(NEQN-3),
     &     0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
     &     0.0,0.0,0.0,U(NEQN-2)


c$$$         URHOFN = URHOL*UL-1.0D0
c$$$         DO 35 I = 1, NPTS
c$$$
c$$$         USQPRT1(I)=
c$$$     $  SQRT((UKD**2.0D0)*(UF0)*(U(NPDE*(I-1)+1)**2.0D0)*(URHOFN**2.0D0)
c$$$     $     + 2.0d0*U(NPDE*(I-1)+3)*UF0*U(NPDE*(I-1)+1)*UKD*URHOFN*(UKI  
c$$$     $     + 2.0D0*UKD))
c$$$
c$$$         UBB1(I) = 
c$$$     $        0.5d0*UF0/UL
c$$$     $     +0.5D0*USQPRT1(I)/(U(NPDE*(I-1)+1)*UL*URHOFN*UKD)
c$$$         UBB2(I) = 
c$$$     $      URHOFN*UKA*UBB1(I)*U(1)*(UF0
c$$$     $        -UL*UBB1(I))/(2.0D0*UKD+UKI)
c$$$
c$$$ 35   CONTINUE



         DO 50 I = 1, NPTS

         WRITE (NOUT,99999) TS, X(I), 
     $           U(NPDE*(I-1)+1), 
     $           U(NPDE*(I-1)+2), 
     $           U(NPDE*(I-1)+3),
c     $           UBB1(I),
c     $           UBB2(I),
     $           U(NEQN-1)
 50      CONTINUE

      WRITE (NOUT,*)
      WRITE (NOUT,*)

      WRITE(6,*) 'Initial conditions set up and written'

      POP1 = .TRUE.
      POP2 = .TRUE.
      POPR1 = .TRUE.
      POPR2 = .TRUE.
      POPK1 = .TRUE.
      POPK2 = .TRUE.


      DO 60 IT = 1, NPLOTS
         TOUT = TPLOT*IT/NPLOTS
*     time to include vesicles on the surface?
*     time to remove vesicles from surface?
*     time to switch on lysis?
         IF (TS .LT. TL1) THEN
            VINF = 0.0e0
         ELSE
            VINF = L1STMP
            IF (POP1) THEN
               WRITE(6,*) 'Vesicles applied: T =',TS,', VINF =',VINF  
            END IF
            POP1 = .FALSE.
         END IF
         IF (TS .GE. TRL1) THEN
            VINF = 0.0e0
            IF (POPR1) THEN
               WRITE(6,*) 'Vesicles off: T =',TS,', VINF =',VINF  
            END IF
            POPR1 = .FALSE.
         END IF
c$$$         IF (TS .LT. TKL1) THEN
c$$$            KILL1 = 0.0e0
c$$$         ELSE
c$$$            KILL1 = K1TMP
c$$$            IF (POPK1) THEN
c$$$               WRITE(6,*) 'Mphi1 Lysis: T =',TS,', K1 =',KILL1  
c$$$            END IF
c$$$            POPK1 = .FALSE.
c$$$         END IF
c$$$         IF (TS .LT. TL2) THEN
c$$$            LINF2 = 0.0e0
c$$$         ELSE
c$$$            LINF2 = L2STMP
c$$$            IF (POP2) THEN
c$$$               WRITE(6,*) 'Mphi2: T =',TS,', VINF =',LINF2  
c$$$            END IF
c$$$            POP2 = .FALSE.
c$$$         END IF
c$$$         IF (TS .GE. TRL2) THEN
c$$$            LINF2 = 0.0e0
c$$$            VINF = L1STMP
c$$$            IF (POPR2) THEN
c$$$               WRITE(6,*) 'Mphi2 off: T =',TS,', LINF2 =',LINF2  
c$$$               WRITE(6,*) 'Mphi1: T =',TS,', VINF =',VINF  
c$$$            END IF
c$$$            POPR2 = .FALSE.
c$$$         END IF
c$$$         IF (TS .LT. TKL2) THEN
c$$$            KILL2 = 0.0e0
c$$$         ELSE
c$$$            KILL2 = K2TMP
c$$$            IF (POPK2) THEN
c$$$               WRITE(6,*) 'Mphi2 Lysis: T =',TS,', K2 =',KILL2  
c$$$            END IF
c$$$            POPK2 = .FALSE.
c$$$         END IF

         IFAIL = -1

    

*     integrate system from ts to tout
         CALL D03PHF(NPDE,M,TS,TOUT,PDEDEF,BNDARY,U,NPTS,X,NCODE,ODEDEF,
     +        NXI,XI,NEQN,RTOL,ATOL,ITOL,NORM,LAOPT,ALGOPT,W,NW,
     +        IW,NIW,ITASK,ITRACE,IND,IFAIL)

   
         WRITE(6,*) 'Time= ',TOUT,', TEND=',TPLOT
         WRITE(6,*) 'm(',TOUT,',end) = ',U(901)

         DO 90 I = 1, NPTS

            WRITE (NOUT,99999) TS, X(I), 
     $           U(NPDE*(I-1)+1), 
     $           U(NPDE*(I-1)+2), 
     $           U(NPDE*(I-1)+3),
c     $           UBB1(I),
c     $           UBB2(I),
     $           U(NEQN-1)
 90      CONTINUE
         WRITE (NOUT,*)
         WRITE (NOUT,*)
*..
*..	  Calculate necrotic/prolif/death radii
         NRR = 0.1d0
         NONEC = .TRUE.
         NORCDM = .TRUE.
         NORCP = .TRUE.
         ICDM = 1
         ICP = 1
         IF (U(1) .GT. NRR) THEN
            NEC = 0.0d0
         ELSE
            DO 70 I = 2, NPTS
               IF (NONEC) THEN
        	  IF (U(NPDE*(I-1)+1) .GT. NRR) THEN
                     NEC = U(NEQN-3)*(X(I) - (U(NPDE*(I-1)+1) - NRR)
     $   *(X(I)-X(I-1))/(U(NPDE*(I-1)+1) - U(NPDE*(I-2)+1))) 
                     NONEC = .FALSE.
                  END IF
               END IF
 70         CONTINUE
         END IF
*..   calculate radius for death threshold
         IF (U(2) .GT. CDM) THEN
            RCDM = 0.0d0
         ELSE
            DO 71 I = 2, NPTS
               IF (NORCDM) THEN
        	  IF (U(NPDE*(I-1)+2) .GT. CDM) THEN
                     RCDM = U(NEQN-3)*(X(I) - (U(NPDE*(I-1)+2) - CDM)
     $   *(X(I)-X(I-1))/(U(NPDE*(I-1)+2) - U(NPDE*(I-2)+2))) 
                     RCDM = U(NEQN-3)*X(I)
                     ICDM = I
                     NORCDM = .FALSE.
        	  END IF
               END IF
 71         CONTINUE
         END IF
*..   calculate radius for prolif threshold
         IF (U(2) .GT. CP) THEN
            RCP = 0.0d0
         ELSE
            DO 72 I = 2, NPTS
               IF (NORCP) THEN
        	  IF (U(NPDE*(I-1)+2) .GT. CP) THEN
                     RCP = U(NEQN-3)*(X(I) - (U(NPDE*(I-1)+2) - CP)
     $   *(X(I)-X(I-1))/(U(NPDE*(I-1)+2) - U(NPDE*(I-2)+2))) 
                     RCP = U(NEQN-3)*X(I)
                     ICP = I
                     NORCP = .FALSE.
        	  END IF
               END IF
 72         CONTINUE
         END IF

*     calculate integrands for total volumes
         DO 80 I = 1, NPTS
            UM(I)  = U(NPDE*(I-1)+1)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $           *FOURPI*U(NEQN-3)
            UL1(I) = U(NPDE*(I-1)+4)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $           *FOURPI*U(NEQN-3)
            UL2(I) = U(NPDE*(I-1)+5)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $           *FOURPI*U(NEQN-3)
 80      CONTINUE
*..   Integrate for total volume
         CALL D01GAF(X,UM,NPTS,UMINT,ERROR,IFAIL)
         CALL D01GAF(X,UL1,NPTS,UL1INT,ERROR,IFAIL)
         CALL D01GAF(X,UL2,NPTS,UL2INT,ERROR,IFAIL)

*     calculate integrands for necrotic volumes
*     set integrand to zero outside region of interest
         IF (ICDM .GT. 1) THEN
            DO 81 I = 1, ICDM-1
               UM(I)  = U(NPDE*(I-1)+1)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $              *FOURPI*U(NEQN-3)
               UL1(I) = U(NPDE*(I-1)+4)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $              *FOURPI*U(NEQN-3)
               UL2(I) = U(NPDE*(I-1)+5)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $              *FOURPI*U(NEQN-3)
 81         CONTINUE
            DO 82 I = ICDM, NPTS
               UM(I)  = 0.0d0
               UL1(I) = 0.0d0
               UL2(I) = 0.0d0
 82         CONTINUE
*..   Integrate for necrotic volume
            CALL D01GAF(X,UM,NPTS,UMNV,ERROR,IFAIL)
            CALL D01GAF(X,UL1,NPTS,UL1NV,ERROR,IFAIL)
            CALL D01GAF(X,UL2,NPTS,UL2NV,ERROR,IFAIL)
         ELSE
            UMNV = 0.0d0
            UL1NV = 0.0d0
            UL2NV = 0.0d0
         END IF

*     calculate integrands for quiescent volumes
*     set integrand to zero outside region of interest
         IF (ICP .GT. 1) THEN
            DO 83 I = 1, ICDM-1
               UM(I)  = 0.0d0
               UL1(I) = 0.0d0
               UL2(I) = 0.0d0
 83         CONTINUE
            DO 84 I = ICDM, ICP-1
               UM(I)  = U(NPDE*(I-1)+1)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $              *FOURPI*U(NEQN-3)
               UL1(I) = U(NPDE*(I-1)+4)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $              *FOURPI*U(NEQN-3)
               UL2(I) = U(NPDE*(I-1)+5)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $              *FOURPI*U(NEQN-3)
 84         CONTINUE
            DO 85 I = ICP, NPTS
               UM(I)  = 0.0d0
               UL1(I) = 0.0d0
               UL2(I) = 0.0d0
 85         CONTINUE
*..   Integrate for quiescent volume
            CALL D01GAF(X,UM,NPTS,UMQV,ERROR,IFAIL)
            CALL D01GAF(X,UL1,NPTS,UL1QV,ERROR,IFAIL)
            CALL D01GAF(X,UL2,NPTS,UL2Qv,ERROR,IFAIL)
         ELSE
            UMQV = 0.0d0
            UL1QV = 0.0d0
            UL2QV = 0.0d0
         END IF

*     calculate integrands for prolif volumes
*     set integrand to zero outside region of interest
         DO 86 I = 1, ICP-1
            UM(I)  = 0.0d0
            UL1(I) = 0.0d0
            UL2(I) = 0.0d0
 86      CONTINUE
         DO 87 I = ICP, NPTS
            UM(I)  = U(NPDE*(I-1)+1)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $           *FOURPI*U(NEQN-3)
            UL1(I) = U(NPDE*(I-1)+4)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $           *FOURPI*U(NEQN-3)
            UL2(I) = U(NPDE*(I-1)+5)*X(I)*X(I)*U(NEQN-3)*U(NEQN-3)
     $           *FOURPI*U(NEQN-3)
 87      CONTINUE
*..   Integrate for prolif volume
         CALL D01GAF(X,UM,NPTS,UMPV,ERROR,IFAIL)
         CALL D01GAF(X,UL1,NPTS,UL1PV,ERROR,IFAIL)
         CALL D01GAF(X,UL2,NPTS,UL2PV,ERROR,IFAIL)
         
*..   Multiply by spheroid radius as we've integrated from 0 to 1
*         UMINT = UMINT*U(NEQN-3)
*         UL1INT = UL1INT*U(NEQN-3)
*         UL2INT = UL2INT*U(NEQN-3)
* should be no need for this as I've taken the factor in the integrand
* may actually be more efficient to multiply by R^3 outside for all these...

         WRITE (NOUT1,99998) TS,NEC,RCDM,RCP,U(NEQN-3),
     &        UMNV,UL1NV,UL2NV,UMQV,UL1QV,UL2QV,UMPV,UL1PV,UL2PV,
     &        UMINT,UL1INT,UL2INT,U(NEQN-2)

*..   calculate advection velocity
         DO 88 I=1,NPTS
            IF (I.EQ.1) THEN
            DUX1=(U(NPDE*I+1)-U(NPDE*(I-1)+1))/(X(I+1)-X(I))
            DUX3=(U(NPDE*I+3)-U(NPDE*(I-1)+3))/(X(I+1)-X(I))
            DUX4=(U(NPDE*I+4)-U(NPDE*(I-1)+4))/(X(I+1)-X(I))
            DUX5=(U(NPDE*I+5)-U(NPDE*(I-1)+5))/(X(I+1)-X(I))
            DUX6=(U(NPDE*I+6)-U(NPDE*(I-1)+6))/(X(I+1)-X(I))
            ENDIF
            IF (I.EQ.NPTS) THEN
            DUX1=(U(NPDE*(I-1)+1)-U(NPDE*(I-2)+1))/(X(I)-X(I-1))
            DUX3=(U(NPDE*(I-1)+3)-U(NPDE*(I-2)+3))/(X(I)-X(I-1))
            DUX4=(U(NPDE*(I-1)+4)-U(NPDE*(I-2)+4))/(X(I)-X(I-1))
            DUX5=(U(NPDE*(I-1)+5)-U(NPDE*(I-2)+5))/(X(I)-X(I-1))
            DUX6=(U(NPDE*(I-1)+6)-U(NPDE*(I-2)+6))/(X(I)-X(I-1))
            ENDIF
            IF ((I.NE.1).AND.(I.NE.NPTS)) THEN
            DUX1=(U(NPDE*I+1)-U(NPDE*(I-2)+1))/(X(I+1)-X(I-1))
            DUX3=(U(NPDE*I+3)-U(NPDE*(I-2)+3))/(X(I+1)-X(I-1))
            DUX4=(U(NPDE*I+4)-U(NPDE*(I-2)+4))/(X(I+1)-X(I-1))
            DUX5=(U(NPDE*I+5)-U(NPDE*(I-2)+5))/(X(I+1)-X(I-1))
            DUX6=(U(NPDE*I+6)-U(NPDE*(I-2)+6))/(X(I+1)-X(I-1))
            ENDIF
            VELC=DTUMFN*DUX1  + DMPFN1*DUX4 
     $     - CHIFN1*U(NPDE*(I-1)+4)*DUX3
     $     + DMPFN2*DUX5 - CHIFN2*U(NPDE*(I-1)+5)*DUX3  
     $     - DP*(DUX1+DUX4+DUX5)
            VELC=DTUMFN*DUX1  + DMPFN1*DUX4 
     $     - CHIFN1*U(NPDE*(I-1)+4)*DUX6
     $     + DMPFN2*DUX5 - CHIFN2*U(NPDE*(I-1)+5)*DUX3  
     $     - DP*(DUX1+DUX4+DUX5)
      	 VELC=(DMPFN1*DUX4/(U(NPDE*I+4)+0.001e0)-CHIFN1*DUX6-VELC)*DUX6
          VBIT= (DMPFN1*DUX4/(U(NPDE*I+4)+0.001e0)-CHIFN1*DUX6
     $           -VELC)*DUX6/(U(NEQN-3)*U(NEQN-3))
            KBIT=BKA/BBA*U(NPDE*I+3)*(1-U(NPDE*I+6))+DL1*U(NPDE*I+6)
     $    -(BKD+BKI)*U(NPDE*I+6)
            WRITE(NOUT2,99997) TS,X(I),VELC,U(NEQN-3)
            WRITE(NOUT2,99997) TS,X(I),VBIT,KBIT,U(NEQN-3)
 88      CONTINUE
         WRITE (NOUT2,*)
         WRITE (NOUT2,*)



*..
*..   continue main loop
 60   CONTINUE
      STOP
*
99997 FORMAT (4(1X,F15.5))
*99997 FORMAT (5(1X,F15.5))
99998 FORMAT (18(1X,F15.5))
99999 FORMAT (7(1X,F15.5))
      End
*

***********************************************************************
***********************************************************************
      SUBROUTINE UVINIT(NPDE,NPTS,X,U,NCODE,NEQN)
***********************************************************************
***********************************************************************
*     Routine for PDE initial values
*..   Scalar Arguments ..
      INTEGER NCODE, NEQN, NPDE, NPTS
*..   Array Arguments ..
      DOUBLE PRECISION U(NEQN), X(NPTS)
*..   Scalars in Common ..
      DOUBLE PRECISION TS
*..   Local Scalars ..
      INTEGER I
*..   Intrinsic Functions ..
      INTRINSIC EXP
*..   Common blocks ..
      COMMON /TAXIS/TS

      DOUBLE PRECISION RINIT
      COMMON /INIT/ RINIT

*.. Executable Statements ..
      DO 20 I = 1, NPTS
*     live tumour cells
         U(NPDE*(I-1)+1) = 0.8e0
*     Nutrient
         U(NPDE*(I-1)+2) = 0.2e0
*     free vesicles
         U(NPDE*(I-1)+3) = 0.0e0
 20   CONTINUE
      U(NEQN-1) = RINIT
      U(NEQN) = 0.0e0
     
      
      RETURN
      END
 
      
*
***********************************************************************
***********************************************************************
      SUBROUTINE ODEDEF(NPDE,T,NCODE,V,VDOT,NXI,XI,UCP,UCPX,RCP,UCPT,
     +     UCPTX,F,IRES)
***********************************************************************
***********************************************************************
*..   Scalar Arguments ..
      DOUBLE PRECISION T
      INTEGER IRES, NCODE, NPDE, NXI
*..   Array Arguments ..
      DOUBLE PRECISION F(*), RCP(NPDE,*), UCP(NPDE,*), UCPT(NPDE,*),
     +     UCPTX(NPDE,*), UCPX(NPDE,*), V(*), VDOT(*),
     +     XI(*)

*..   parameters for spatial terms
      DOUBLE PRECISION CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      COMMON /SPATIAL/ CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      
*..   parameters for boundary conditions
      DOUBLE PRECISION HA,HL1,HL2,VINF,LINF2,P0,QP,HV
      COMMON /BCS/ HA,HL1,HL2,VINF,LINF2,P0,QP,HV

*..   parameters for functions returned by SOURCE 
      DOUBLE PRECISION BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG
      COMMON /SRC/ BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG

 
*..   Executable Statements ..
      IF (IRES.EQ.1) THEN
         F(1) = VDOT(1) - QP*(P0-1.0e0+UCP(1,2))
         F(2) = VDOT(2)
      ELSE IF (IRES.EQ.-1) THEN
         F(1) = VDOT(1)
         F(2) = VDOT(2)
      END IF

      RETURN
      END

*
***********************************************************************
***********************************************************************
      SUBROUTINE PDEDEF(NPDE,T,X,U,UX,NCODE,V,VDOT,P,Q,R,IRES)
***********************************************************************
***********************************************************************
*     .. Parameters ..
c      INTEGER NOUT2
c      PARAMETER (NOUT2=12)
*..   Scalar Arguments ..
      DOUBLE PRECISION T, X
      INTEGER IRES, NCODE, NPDE
*..   Array Arguments ..
      DOUBLE PRECISION P(NPDE,NPDE), Q(NPDE), R(NPDE), U(NPDE),
     +     UX(NPDE), V(*), VDOT(*)

*..   External Subroutines ..
      EXTERNAL SOURCE

*..   Functions returned by SOURCE
      DOUBLE PRECISION DA,DC,DL1,DL2,DM,PA,PM
     $     ,CHIFN1,DMPFN1,CHIFN2,DMPFN2,DTUMFN
     $     ,DVFN,GM

*..   parameters for functions returned by SOURCE 
      DOUBLE PRECISION BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG
      COMMON /SRC/ BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG

*..   parameters for spatial terms returned by SOURCE
      DOUBLE PRECISION CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      COMMON /SPATIAL/ CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV

c$$$
c$$$*..   binding params
c$$$      DOUBLE PRECISION UKA,UKD,UKI,UL,UF0
c$$$      COMMON /BINDING/ UKA,UKD,UKI,UL,UF0

*..   binding params
      DOUBLE PRECISION UKA,UKD,UKI,UL,UF0,URHOL,URHOF
      COMMON /BINDING/ UKA,UKD,UKI,UL,UF0,URHOL,URHOF

*..   internal definition for common terms
      DOUBLE PRECISION VEL,UB1,UB2,UB,UGAMMA1,UBETA,UALPHA,ULL,UGAMMA2
     $     ,UFF

*..   Executable Statements ..
      CALL SOURCE(U,DA,DC,DL1,DL2,DM,PA,PM
     $     ,CHIFN1,DMPFN1,CHIFN2,DMPFN2,DTUMFN,DVFN,NPDE)

*..   advection velocity
      VEL = (DTUMFN-DP)*UX(1)
*..   expression for bound receptors
      UB = UKA*UL*U(3)*UF0/(UKA*UL*U(3)+UKD+UKI)
      
      ULL = UL*U(3)
      UGAMMA1 = UKA*URHOF*(URHOL*UL-1.0D0)
c      UGAMMA1 = 0.001d0
      UALPHA = (2.0D0*UKD+UKI)*UKA*ULL + 2.0D0*(UKD**2.0D0)  
     $     + (UKI+UGAMMA1*UF0*U(1)+3.0D0*UKD)*UKI
      
      UBETA = ((4.0D0*UKD*UKI+4.0D0*UKD**2.0D0+(UKI**2.0D0))*UKA*ULL
     $     + ((6.0D0*UKI**2.0D0)+20.0D0*UKD*UKI
     $                   + 16.0D0*(UKD**2.0D0))*UF0*U(1)*UGAMMA1
     $     + UGAMMA1*(UKD**3.0D0)
     $     + ((2.0D0*UKI**2.0D0)+10.0D0*UKD*UKI
     $                   + 16.0D0*(UKD**2.0D0))*UKI)*UKA*ULL
     $     + ((4.0D0*(UKD**2.0D0)+6.0D0*UKD*UKI+2.0D0*UKI**2.0D0)*UKI
     $     + UF0*U(1)*UGAMMA1*(UKI**2.0D0))*UGAMMA1*UF0*U(1)
     $     + (12.0D0*(UKD**3.0D0)+13.0d0*UKI*(UKD**2.0D0)
     $     +6.0D0*UKD*UKI**2.0D0 + UKI**3.0D0)*UKI + 4.0D0*UKD**4.0D0


         UB1=-1.0D0*(UALPHA-SQRT(UBETA))/(2.0D0*UGAMMA1*(2.0D0*UKD+UKI))

         UB2 = UGAMMA1*(UF0*U(1)-UB1)*UB1/(2.0D0*UGAMMA1*UB1
     $        +2.0D0*UKD+UKI)

         UFF = UF0*U(1)-UB1-2.0D0*UB2

c         write(6,*) UB,UALPHA,SQRT(UBETA),UB1,UB2
c         write(6,*) UGAMMA2,


*..   tumour lysis rate      
c      GM = UG*UKI*UB
      GM = UG*UKI*(UB1+UB2)

*..   Live tumour cell equation
      P(1,1) = 1.0e0
      P(1,2) = 0.0e0
      P(1,3) = 0.0e0
      R(1) = (DTUMFN*UX(1) - VEL*U(1) )/(V(1)*V(1))
      Q(1) = -X*VDOT(1)*UX(1)/V(1) 
     $        - PM*U(1) 
     $        + DM*U(1)
     $        + GM*U(1)
      
*..   Nutrient equation:
      P(2,1) = 0.0e0
      P(2,2) = 1e-3
      P(2,3) = 0.0e0
      R(2) = UX(2)
      Q(2) = V(1)*V(1)*DC

*..   Vesicle equation:
      P(3,1) = 0.0e0
      P(3,2) = 0.0e0
      P(3,3) = 1.0e0
      R(3) = (DVFN*UX(3)-VEL*U(3))/(V(1)*V(1))
      Q(3) = -X*VDOT(1)*UX(3)/V(1)
     $     +UKA*U(3)*UFF - UKD*UB1/UL
c     $        + UKI*UB/UL
      
 
      RETURN
      END

*
***********************************************************************
***********************************************************************
      SUBROUTINE BNDARY(NPDE,T,U,UX,NCODE,V,VDOT,IBND,BETA,GAMMA,IRES)
***********************************************************************
***********************************************************************
*..   Scalar Arguments ..
      DOUBLE PRECISION T
      INTEGER IBND, IRES, NCODE, NPDE
*..   Array Arguments ..
      DOUBLE PRECISION BETA(NPDE), GAMMA(NPDE), U(NPDE), UX(NPDE),
     +     V(*), VDOT(*)
*..   Intrinsic Functions ..
      INTRINSIC EXP

*..   parameters for boundary conditions
      DOUBLE PRECISION HA,HL1,HL2,VINF,LINF2,P0,QP,HV
      COMMON /BCS/ HA,HL1,HL2,VINF,LINF2,P0,QP,HV

*..   parameters for spatial terms
      DOUBLE PRECISION CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      COMMON /SPATIAL/ CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV

*..   Executable Statements ..
      IF (IBND.EQ.0) THEN
         BETA(1) = 0.0e0
         GAMMA(1) = UX(1)
         BETA(2) = 1.0e0
         GAMMA(2) = 0.0e0
         BETA(3) = 0.0e0
         GAMMA(3) = UX(3)
      ELSE 
         BETA(1) = 0.0e0
         GAMMA(1) = DP*UX(1)*U(1)
     $        + V(1)*QP*(P0-1.0e0+U(1))*U(1)
     $        + DTUM*(1.0e0-U(1))*UX(1)
         BETA(2) = 0.0e0
         GAMMA(2) = U(2)-1
         BETA(3) = 0.0e0
         GAMMA(3) = DV*UX(3)*U(1) - DTUM*UX(1)*U(3)
     $        -V(1)*HV*(VINF-U(3))*U(1) 
      END IF
      RETURN
      END

***********************************************************************
***********************************************************************
      SUBROUTINE SOURCE(U,DA,DC,DL1,DL2,DM,PA,PM
     $     ,CHIFN1,DMPFN1,CHIFN2,DMPFN2,DTUMFN,DVFN,NPDE)
***********************************************************************
***********************************************************************
*..   Functional forms returned by SOURCE ..
      DOUBLE PRECISION DA,DC,DL1,DL2,DM,PA,PM
     $     ,CHIFN1,DMPFN1,CHIFN2,DMPFN2,DTUMFN
     $     ,DVFN
      INTEGER NPDE
*..   Array Arguments ..
      DOUBLE PRECISION U(NPDE)

*..   internal definition for common terms
      DOUBLE PRECISION SCAL,UCAL1,UCAL2,UCAM,UCCL1,UCCL2,UCCM,UCDL1
     $     ,UCDL2,UCDM,UCK1,UCK2,UCMMSP,UCPC,UP

*..   parameters for functions returned by SOURCE 
      DOUBLE PRECISION BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG
      COMMON /SRC/ BCP,BL1,BL2,BM,CCAL1,CCAL2,CCAM,CCCL1,CCCL2,CCCM
     $     ,CCDL1,CCDL2,CCDM,CCK1,CCK2,CCPC,DF,KILL1,KILL2,L1A,L1AD,L2A
     $     ,L2AD,MA,MAD,MM,PCAP,PINF,SCAAL1,SCAAL2,SCACL1,SCACL2
     $     ,SCADL1,SCADL2,SCALAM,SCALCM,SCALDM,SCALK1,SCALK2,SCALNC,UG

*..   parameters for spatial terms returned by SOURCE
      DOUBLE PRECISION CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV
      COMMON /SPATIAL/ CCSP,CHI1,CHI2,DMPHI1,DMPHI2,DP,DTUM,MINSP
     $     ,SCALSP,DV

*..   powers for proliferation
      UCPC = U(2)**MM
*..   powers for death
      UCDM = U(2)**MM
      UCDL1 = U(2)**MM
      UCDL2 = U(2)**MM
      UCK1 = U(2)**MM
      UCK2 = U(2)**MM
*..   powers for nutrient consumption
      UCCM = U(2)**MM
      UCCL1 = U(2)**MM 
      UCCL2 = U(2)**MM
*..   powers for chemoattractant production
      UCAM = U(2)**MM
      UCAL1 = U(2)**MM
      UCAL2 = U(2)**MM
*..   power for movement scaling
      UCMMS = U(2)**MM
*.. proportion of cellular material
      UP = 1.0e0-U(1)

*.. capped prolif - no increase when N>PCAP
*.. normally PCAP=PINF, i.e. capping when 
*.. material is optimal
*.. to get simple linear dependence (no capping) just have pcap=1
      IF (UP .GT. PCAP) THEN
        SCALNC = (CCPC+1.0e0)*PCAP/PINF
      ELSE
        SCALNC = (CCPC+1.0e0)*UP/PINF
      END IF
      PM = SCALNC*UCPC/(CCPC+UCPC)
*..   cell death rate
      DM = BM*(1.0e0 - SCALDM*UCDM/(CCDM+UCDM))
*..   nutrient consumption rate
      DC = SCALCM*U(1)*UCCM/(CCCM+UCCM)+ BCP*PM*U(1)


*..   scale movement rates?
      SCAL = MINSP + SCALSP*(UCMMSP)/(CCSP+UCMMSP)
      CHIFN1 = CHI1*SCAL
      DMPFN1 = DMPHI1*SCAL
      CHIFN2 = CHI2*SCAL
      DMPFN2 = DMPHI2*SCAL
      DTUMFN = DTUM*SCAL
      DVFN = DV*SCAL

      RETURN
      END


***********************************************************************
***********************************************************************









