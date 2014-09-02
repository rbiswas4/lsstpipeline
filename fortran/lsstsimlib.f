      PROGRAM LSSTSIMLIB(OPSIMF, SIMLIB) 
C
C Take LSST opsim output and produce sdss simlib input
C
      IMPLICIT NONE
C
      CHARACTER LINE*500
C
      INTEGER OID,SID,FID,PID,SNO,EDATE,ETIME,FINT,FFINT,
     &        PAIRN,NIGHT
      REAL SLEWT,SLEWD,ROTS,ROTT,PRANK,FRANK,MSEE,RSEE,SEE,
     &     XPAR,CSEE,AMASS,SKYBV,VSKYB,RA,DEC,LST,ALT,AZI,
     &     MOOND,MOONRA,MOONDEC,MOONALT,MOONPHA,
     &     SUNA,SUNAZ,PHANG,RSCAT,MIESCAT,MOONIL,MOONB,DARKB,
     &     F5SIG,PSKYB,P5SIG,MSKYB,M5SIG,HEXDIRA,HEXDIDEC,FVIS,
     &     VERTEX,SOLEL
      REAL*8 MJD,MJDL
      CHARACTER FILT*1,SUBQ*7
C These in OpsimObs, but not Obsim
      REAL MOONAZ,MOONPHAO
      INTEGER DOWN
      CHARACTER MAKM*1
C
C List of all LSST fields that are covered in the opsim file
C
      INTEGER MXFLD,MXEP,MSURV,MXOFLD,MXOEP
      PARAMETER (MXFLD=4000,MXEP=2000,MSURV=1,MXOFLD=100,MXOEP=100000)
      INTEGER NFLD,IFLD(MXFLD),IEP,IFL,IFLC,IOFLD(MXFLD)
      INTEGER NFLDW,NFLTW
      LOGICAL OVFL(MXFLD)
      INTEGER IOVF
C
C List of data for each field
C
      INTEGER NEP(MXFLD),FIDID(MXFLD,MXEP)
      REAL*8 FMJD(MXFLD,MXEP)
      REAL FSKYB(MXFLD,MXEP),
     &     PSF(MXFLD,MXEP),FSIG5(MXFLD,MXEP)
      REAL FRA(MXFLD),FDEC(MXFLD),
     &     FMWEBV(MXFLD),FPIXS(MXFLD)
      CHARACTER FFLT(MXFLD,MXEP)
      INTEGER IFLT(MXFLD,MXEP),IORD(MXFLD,MXEP)
C
C Overflow fields
C
      INTEGER NEPO(MXOFLD),FIDIDO(MXOFLD,MXOEP)
      REAL*8 FMJDO(MXOFLD,MXOEP)
      REAL FSKYBO(MXOFLD,MXOEP),
     &     PSFO(MXOFLD,MXOEP),FSIG5O(MXOFLD,MXOEP)
      REAL FRAO(MXOFLD),FDECO(MXOFLD),
     &     FMWEBVO(MXOFLD),FPIXSO(MXOFLD)
      CHARACTER FFLTO(MXOFLD,MXOEP)
      INTEGER IFLTO(MXOFLD,MXOEP),IORDO(MXOFLD,MXOEP)
C
      INTEGER NLSSTF,IFIL,IFILT
      PARAMETER(NLSSTF=6)
      CHARACTER FILTS(NLSSTF),FILTR(NLSSTF)
      DATA FILTS/'u','g','r','i','z','y'/
      DATA FILTR/'u','g','r','i','z','Y'/
C
C Input file
C
      CHARACTER OPSIMF*100
C
C Output file
C
      CHARACTER SIMLIB*100
C
C Type of input
C
      LOGICAL OPSIM
      LOGICAL OPSIMOBS
C
C On the fly derived
C
      REAL SKYS,PSF1,RAD,DECD
      INTEGER COUNT
C
C Limitation on how many epoch to write
C
      INTEGER MXEPWR,NEPW
      PARAMETER (MXEPWR=MXEP+MXOEP)
      INTEGER NTGDW
C
C Zeropoint
C
      REAL ZPT5,F5,D51,D52,D5,N5
C
C Sky noise
C
      REAL SKYM,SKYF
C
      REAL TSL
      INTEGER NNIT,JEP,INIT,IWEP,JEPN,JEPX
C
C Paramters for simlib
C
      REAL GAIN,SNRLM
      PARAMETER (GAIN=1.0,SNRLM=5.0) 
C
C Stuff to parse line from OpsimObs
      INTEGER IH,ILOC,ISP,ISTART,ISTOP,IVAR
      CHARACTER FORM*5
C
C------------------------------------------------------------------
C
C Zero out the field data
C
      COUNT = 0
C
C Open the LSST opsim file
C
C     WRITE (6,'(A)') 'What Opsim or OpsimObs to process?'
C     READ (5,'(A)') OPSIMF
      OPEN(UNIT=30,FILE=OPSIMF,FORM='FORMATTED',STATUS='OLD')
C
C Initialize some stuff
C
      MJDL = 0.0
C
C Read the header line
C
      READ (30,'(A)') LINE
C
      IF (LINE(1:1).EQ.'o') THEN
C This is an Opsim output
        OPSIM = .TRUE.
        OPSIMOBS = .FALSE.
      ELSEIF (LINE(1:1).EQ.'#') THEN
C This is an OpsimObs output
        OPSIM = .FALSE.
        OPSIMOBS = .TRUE.
      ELSE
        WRITE (6,*) 'Unrecognized type of input.'
        STOP
      ENDIF
C
      NFLD = 0
      IOVF = 0
      DO IFL = 1,MXFLD
        NEP(IFL) = 0
        OVFL(IFL) = .FALSE.
        IOFLD(IFL) = 0
      ENDDO
      DO IFL = 1,MXOFLD
        NEPO(IFL) = 0
      ENDDO
C
C Read the opsim data
C
   20 CONTINUE
        READ (30,'(A)',END=30) LINE
C       WRITE (6,'(A)') 'What Opsim or OpsimObs to process?'
C Skip comments
        IF (LINE(1:1).EQ.'#') GOTO 20
        COUNT = COUNT + 1
C Debugging
C          IF (COUNT.GT.10) GOTO 30
        IF (OPSIM) THEN
C Parse the line, tab separated data
          READ (LINE,*) OID,SID,PID,FID,FILT,SNO,SUBQ,PAIRN,EDATE,MJD,
     &                NIGHT,ETIME,SLEWT,SLEWD,ROTS,ROTT,FVIS,FINT,FFINT,
     &                PRANK,FRANK,MSEE,RSEE,SEE,XPAR,CSEE,AMASS,
     &                SKYBV,VSKYB,RA,DEC,LST,ALT,AZI,MOOND,MOONRA,
     &                MOONDEC,MOONALT,MOONPHA,SUNA,SUNAZ,PHANG,RSCAT,
     &                MIESCAT,MOONIL,MOONB,DARKB,SOLEL,F5SIG,PSKYB,
     &                P5SIG,MSKYB,M5SIG,HEXDIRA,HEXDIDEC,VERTEX
          DOWN = 0
        ELSEIF (OPSIMOBS) THEN
          READ (LINE,*) RA,DEC,MJD,FILT,ALT,AZI,AMASS,XPAR,RSEE,SEE,
     &                  SUNA,SUNAZ,MOONALT,MOONAZ,MOONPHA,MOONPHAO,
     &                  MSKYB,M5SIG,NIGHT,DOWN,FID,MAKM
C No proposal ID.  Set it to 1.
          PID = 1
C Different format for RA and DEC than Opsim (deg -> radians)
           RA  =  RA*3.1416/180.0
           DEC = DEC*3.1416/180.0
        ENDIF
C Down?
        IF (DOWN.EQ.1) GOTO 20
        IF (MJD.LT.MJDL) THEN
          WRITE (6,*) 'Opsim data needs to be sorted by time.'
          STOP
        ELSEIF (MJD.EQ.MJDL) THEN
C Skip repeat visit.  Simply on to the next.
          GOTO 20
        ELSE
          MJDL = MJD
        ENDIF
C
C Weather
C
        IF (XPAR.GT.(5.0/8.0)) THEN
C Weather too bad for observing.  Simply on to the next.
          GOTO 20
        ENDIF
C
        IF (MOD(COUNT,10000).EQ.0) write (6,*) mjd,count
C
C OK.  Index this view
C
        IF (NFLD.GT.0) THEN
          DO IFL = 1,NFLD
            IF (FID.EQ.IFLD(IFL)) THEN
              IFLC = IFL
              GOTO 11
            ENDIF
          ENDDO
        ENDIF
C This must be a new field
        NFLD = NFLD + 1
        IF (NFLD.GT.MXFLD) THEN
          WRITE (6,*) 'Too many fields.  Increase MXFLD'
          STOP
        ENDIF
        IFLC = NFLD
        IFLD(IFLC) = FID
        FRA(IFLC) = RA
        FDEC(IFLC) = DEC
C for the moment set this to zero
        FMWEBV(IFLC) = 0.0
C LSST pixels are 0.2 arc sec
        FPIXS(IFLC) = 0.2
 11     CONTINUE
        IF (.NOT.OVFL(IFLC)) THEN
          NEP(IFLC) = NEP(IFLC) + 1
          FMJD(IFLC,NEP(IFLC)) = MJD
C These should use the modified values
          FSKYB(IFLC,NEP(IFLC)) = MSKYB
          FSIG5(IFLC,NEP(IFLC)) = M5SIG
          PSF(IFLC,NEP(IFLC)) = SEE
          DO IFIL = 1,NLSSTF
            IF (FILT.EQ.FILTS(IFIL)) THEN 
              IFLT(IFLC,NEP(IFLC)) = IFIL
              IFILT = IFIL
            ENDIF
          ENDDO
          FFLT(IFLC,NEP(IFLC)) = FILTR(IFILT)
C Make an ID for this from the propID and field ID
          FIDID(IFLC,NEP(IFLC)) = PID*10000 + FID
          IF (NEP(IFLC).EQ.MXEP) THEN
            WRITE (6,*)'Normal Epoch limit of ',MXEP,' for FIELD ',FID
            WRITE (6,*) 'Overflowing'
            OVFL(IFLC) = .TRUE.
            IOVF = IOVF + 1
            IF (IOVF.GT.MXOFLD) THEN
              WRITE (6,*) 'Too many overflow fields.  Increase MXOFLD'
              STOP
            ENDIF
            IOFLD(IFLC) = IOVF
            NEPO(IOVF) = NEP(IFLC)
C And copy in everything from the regular field into the overflow
            DO IEP = 1,NEP(IFLC)
              FMJDO(IOVF,IEP) = FMJD(IFLC,IEP)
              FSKYBO(IOVF,IEP) = FSKYB(IFLC,IEP)
              FSIG5O(IOVF,IEP) = FSIG5(IFLC,IEP)
              PSFO(IOVF,IEP) = PSF(IFLC,IEP)
              IFLTO(IOVF,IEP) = IFLT(IFLC,IEP)
              FFLTO(IOVF,IEP) = FFLT(IFLC,IEP)
              FIDIDO(IOVF,IEP) = FIDID(IFLC,IEP)
            ENDDO
          ENDIF
        ELSE
          IFL = IOFLD(IFLC)
          NEPO(IFL) = NEPO(IFL) + 1
          FMJDO(IFL,NEPO(IFL)) = MJD
C These should use the modified values
          FSKYBO(IFL,NEPO(IFL)) = MSKYB
          FSIG5O(IFL,NEPO(IFL)) = M5SIG
          PSFO(IFL,NEPO(IFL)) = SEE
          DO IFIL = 1,NLSSTF
            IF (FILT.EQ.FILTS(IFIL)) THEN 
              IFLTO(IFL,NEPO(IFL)) = IFIL
              IFILT = IFIL
            ENDIF
          ENDDO
          FFLTO(IFL,NEPO(IFL)) = FILTR(IFILT)
C Make an ID for this from the propID and field ID
          FIDIDO(IFL,NEPO(IFL)) = PID*10000 + FID
          IF (NEPO(IFL).EQ.MXOEP) THEN
            WRITE (6,*) 'Epoch limit of ',MXOEP,
     &                  ' for overflow FIELD ',FID
            STOP
          ENDIF
        ENDIF
C New line in opsim file
        GOTO 20
C Done with opsim file
   30 CONTINUE
C Reorder by filter for one night.
      DO IFLC = 1,NFLD
        IF (.NOT.OVFL(IFLC)) THEN
          IF (NEP(IFLC).GT.1) THEN
            NNIT = 0
            JEPN = 0
            DO IEP = 1,NEP(IFLC)
              NNIT = NNIT + 1
              IF (IEP.EQ.1) THEN
C first Epoch is special.
                TSL = FMJD(IFLC,2) - FMJD(IFLC,1)
              ELSE
                TSL = FMJD(IFLC,IEP) - FMJD(IFLC,IEP-1)
              ENDIF
              IF (TSL.GT.0.6) THEN
                JEPN = IEP - NNIT
                IF (JEPN.EQ.0) THEN
                  JEPN = 1
                  NNIT = NNIT - 1
                ENDIF
                JEPX = IEP - 1
                IF (JEPX.EQ.0) JEPN = 0
              ENDIF
              IF (IEP.EQ.NEP(IFLC)) THEN
                JEPN = IEP - NNIT
                IF (JEPN.EQ.0) THEN
                  JEPN = 1
                  NNIT = NNIT - 1
                ENDIF
                JEPX = IEP
              ENDIF
              IF (JEPN.NE.0) THEN
                INIT = 0
                DO IFIL = 1,NLSSTF
                  DO JEP = JEPN,JEPX
                    IF (IFLT(IFLC,JEP).EQ.IFIL) THEN
                      INIT = INIT + 1
                      IORD(IFLC,JEP) = IEP - NNIT + INIT - 1
                    ENDIF
                  ENDDO
                ENDDO
                NNIT = 0
                JEPN = 0
              ENDIF
            ENDDO
          ELSE
            IORD(IFLC,1) = 1
          ENDIF
C And the same for the overflow
        ELSE
          IFL = IOFLD(IFLC)
          IF (NEP(IFLC).GT.1) THEN
            NNIT = 0
            JEPN = 0
            DO IEP = 1,NEPO(IFL)
              NNIT = NNIT + 1
              IF (IEP.EQ.1) THEN
                TSL = FMJDO(IFL,2) - FMJDO(IFL,1)
              ELSE
                TSL = FMJDO(IFL,IEP) - FMJDO(IFL,IEP-1)
              ENDIF
              IF (TSL.GT.0.6) THEN
                JEPN = IEP - NNIT
                IF (JEPN.EQ.0) THEN
                  JEPN = 1
                  NNIT = NNIT - 1
                ENDIF
                JEPX = IEP - 1
                IF (JEPX.EQ.0) JEPN = 0
              ENDIF
              IF (IEP.EQ.NEPO(IFL)) THEN
                JEPN = IEP - NNIT
                IF (JEPN.EQ.0) THEN
                  JEPN = 1
                  NNIT = NNIT - 1
                ENDIF
                JEPX = IEP
              ENDIF
              IF (JEPN.NE.0) THEN
                INIT = 0
                DO IFIL = 1,NLSSTF
                  DO JEP = JEPN,JEPX
                    IF (IFLTO(IFL,JEP).EQ.IFIL) THEN
                      INIT = INIT + 1
                      IORDO(IFL,JEP) = IEP - NNIT + INIT - 1
                    ENDIF
                  ENDDO
                ENDDO
                NNIT = 0
                JEPN = 0
              ENDIF
            ENDDO
          ELSE
            IORDO(IFL,1) = 1
          ENDIF
        ENDIF
      ENDDO

C Open the Simlib file
C
C     WRITE (6,'(A)') 'What simlib to write?'
C     READ (5,'(A)') SIMLIB
      OPEN(UNIT=32,FILE=SIMLIB,FORM='FORMATTED',STATUS='NEW')
C
C Write the header
C
      WRITE (32,'(A)') 
     & 'SURVEY: LSST    FILTERS: ugrizY  TELESCOPE: LSST'
      WRITE (32,'(A)') 'USER: cinabro     HOST: motor1'
      WRITE (32,'(A,A)') 'COMMENT: LSST ',OPSIM
      WRITE (32,'(A)') 'BEGIN LIBGEN'
      WRITE (32,'(A)') ''
      WRITE (32,'(A)') ''
C
C Loop over the fields
C
      NFLDW = 0
      DO IFLC = 1,NFLD
        NFLDW = NFLDW + 1
C
        RAD = FRA(IFLC)*180.0/3.1416
        DECD = FDEC(IFLC)*180.0/3.1416
C Write the view header
        IF (OVFL(IFLC)) THEN
          NEPW = MAX(NEP(IFLC),NEPO(IOFLD(IFLC)))
        ELSE
          NEPW = NEP(IFLC)
        ENDIF
        WRITE (32,'(A)') '# ----------------------------------------'
        WRITE (32,'(A,I4)') 'LIBID: ',IFLD(IFLC)
        WRITE (32,'(A,F10.6,A,F10.6,A,I5,A,F6.3,A,F6.3,A,I5)')
     & 'RA: ',RAD,'    DECL: ',DECD,'   NOBS: ',
     & MIN(NEPW,MXEPWR),
     & '    MWEBV:',FMWEBV(IFLC),'   PIXSIZE:',FPIXS(IFLC),
     & '    FIELD:',IFLD(IFLC)
        WRITE (32,'(A)') ''
        WRITE (32,'(A1,27X,A)') '#','CCD  CCD        PSF1 PSF2 PSF2/1'
        WRITE (32,'(A,A)') '#     MJD        IDEXPT  FLT GAIN NOISE',
     &    ' SKYSIG (pixels)  RATIO  ZPTAVG ZPTSIG  MAG'
        NTGDW = 0
        IF (.NOT.OVFL(IFLC)) THEN
          DO IEP = 1,NEP(IFLC)
            DO JEP = 1,NEP(IFLC)
              IF (IORD(IFLC,JEP).EQ.IEP) IWEP = JEP
            ENDDO 
            NTGDW = NTGDW + 1
C Use Kessler routine to compute psf, zpt and skysig
            CALL OPSIM2SIMLIB(FSIG5(IFLC,IWEP),SNRLM,
     &                          FSKYB(IFLC,IWEP),
     &                          PSF(IFLC,IWEP),ZPT5,SKYS,PSF1)
            IF (NTGDW.LE.MXEPWR) THEN
              WRITE (32,
     &    '(A,F11.5,2X,I9,1X,A,2X,A,F6.2,1X,F5.2,A,F5.2,A)')
     &    'S: ',FMJD(IFLC,IWEP),FIDID(IFLC,IWEP),
     &    FFLT(IFLC,IWEP),
     &    '1.00  0.25 ',SKYS,PSF1,
     &    ' 0.00 0.000  ',ZPT5,'  0.005 -99.000'
            ENDIF
          ENDDO
C And over flows
        ELSE
          IFL = IOFLD(IFLC)
          DO IEP = 1,NEPO(IFL)
            DO JEP = 1,NEPO(IFL)
              IF (IORDO(IFL,JEP).EQ.IEP) IWEP = JEP
            ENDDO 
            NTGDW = NTGDW + 1
C Use Kessler routine to compute psf, zpt and skysig
            CALL OPSIM2SIMLIB(FSIG5O(IFL,IWEP),SNRLM,
     &                            FSKYBO(IFL,IWEP),
     &                          PSFO(IFL,IWEP),ZPT5,SKYS,PSF1)
            IF (NTGDW.LE.MXEPWR) THEN
              WRITE (32,
     &    '(A,F11.5,2X,I9,1X,A,2X,A,F6.2,1X,F5.2,A,F5.2,A)')
     &    'S: ',FMJDO(IFL,IWEP),FIDIDO(IFL,IWEP),
     &    FFLTO(IFL,IWEP),
     &    '1.00  0.25 ',SKYS,PSF1,
     &    ' 0.00 0.000  ',ZPT5,'  0.005 -99.000'
            ENDIF
          ENDDO
        ENDIF
C And the field footer
        WRITE (32,'(A,I4)') 'END_LIBID: ',IFLD(IFLC)
      ENDDO
C
C and the simlib footer
C
      WRITE (32,'(A)') ''
      WRITE (32,'(A,I4,A)') 'END_OF_SIMLIB: ',NFLDW,' ENTRIES'
C
C and that is it
C
      CLOSE (32)      
C
      END

      SUBROUTINE OPSIM2SIMLIB(
     &    OPSIM_MAGLIM,  OPSIM_SNR, OPSIM_MAGSKY, OPSIM_SEEING  ! (I)
     &   ,SIMLIB_ZPTAVG, SIMLIB_SKYSIG, SIMLIB_PSF1    ) ! (O)

      IMPLICIT NONE

c Created Jun 10, 2009 by R.Kessler
c convert LSST OPSIM args into SIMLIB args for the SNANA sim.

c subroutine args

      REAL
     &    OPSIM_MAGLIM   ! (I) 5sigma_ps (limiting mag)
     &  , OPSIM_SNR      ! (I) SNR for limiting mag
     &  , OPSIM_MAGSKY   ! (I) perry_skybrightness (mag/arcsec^2)
     &  , OPSIM_SEEING   ! (I) FWHMA, arcsec
c
     &  , SIMLIB_ZPTAVG  ! (O) Flux = 10^[ -0.4*(mag - ZPTAVG)]
     &  , SIMLIB_SKYSIG  ! (O) sky noise (p.e.) per pixel
     &  , SIMLIB_PSF1    ! (O) Gaussian sigma, pixels


      REAL SNR_MAGLIM, PIXSIZE, NPIX_ASEC, TEN
      PARAMETER ( PIXSIZE    = 0.2 )  ! pixel size, arcsec
      PARAMETER ( NPIX_ASEC = 1./(PIXSIZE*PIXSIZE) )  ! Npix per arcsec
      PARAMETER ( TEN = 10.0 )
c local variables.

      REAL
     &   AREA, ARG, TMP
     &  ,ZPT_APPROX, ZPT_COR
     &  ,FSKY_ASEC, FSKY_PIX

C -------------- BEGIN ----------
 
      SIMLIB_PSF1 = (OPSIM_SEEING/2.35)/PIXSIZE

      AREA = (1.51*OPSIM_SEEING)**2
   
      arg = AREA * OPSIM_SNR * OPSIM_SNR
      zpt_approx = 2.*OPSIM_MAGLIM - OPSIM_MAGSKY + 2.5*log10(ARG)

c tack  on source/noise correction term
      ARG = -0.4*( OPSIM_MAGSKY - OPSIM_MAGLIM )  ! usually this is positive
      TMP = TEN**(ARG)
      ZPT_COR = 2.5*log10( 1. + 1./(AREA*TMP) )

      SIMLIB_ZPTAVG = ZPT_APPROX + ZPT_COR

c get sky noise in photoelectrons per pixel

      ARG       = -0.4*( OPSIM_MAGSKY - SIMLIB_ZPTAVG )
      FSKY_ASEC = TEN**(ARG)               ! sky Npe/arcsec
      FSKY_PIX  = FSKY_ASEC / NPIX_ASEC    ! sky Npe per pixel

      SIMLIB_SKYSIG = sqrt(FSKY_PIX)

      RETURN
      END

