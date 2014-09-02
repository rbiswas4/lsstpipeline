      SUBROUTINE LSSTDFACG(DFACGF, OUTFILE , ACT) 
C
C D.Cinabro, version 1.0, 1 April 2013
C
C Take a description of the cadence for a set of LSST deep 
C fields and produce the input file for Lynne Jones' OpsimObs
C program.
C Input looks like:
C # Skip Moon Bright Time, default is 75.0
C MOONL: 75.0 (moon illumination level (0.0-100.0) above which not to observe)
C # Make up for missed observing epochs, default is 1, try to make up
C MAKEUP: 1 (1 = make up for missed obs, 0 = not)
C # Do not allow observations at high air mass, default cut above 2.0
C AIRMASS: 2.0
C # Field ID RA DEC (degrees) Start(Years) End(Years)
C FIELDID: 1 
C FIELDRA: 1.0 
C FIELDDEC: 1.0 
C FEILDSTART: 0.0 
C FEILDEND: 10.0
C # Filter Cadence(integer days) Exposure (integer seconds, 34 = 1 LSST standard exposure)
C FIELDFILTER: u 
C FIELDCAD: 4
C FIELDEXP: 34
C . . .
C
      IMPLICIT NONE
C
      CHARACTER LINE*500
      CHARACTER ACT*5
C
C Cut for weather.  No data if cloud.gt.5/8.
C
      REAL CLOUDM
      PARAMETER (CLOUDM=5.0/8.0)
C
C Default maximum moon phase to take data
C
      REAL MOONL,MOONLD
      PARAMETER (MOONLD=0.75)
C
C Make up or not, default is to make up
C
      INTEGER MAKEUP,MAKEUPD
      PARAMETER (MAKEUPD=1)
      LOGICAL MAKEU
C
C Default Airmass cut.
C
      REAL AIRMC,AIRM
      PARAMETER (AIRMC=2.0)
C
C Maximum number of fields
C
      INTEGER MXFLD
      PARAMETER (MXFLD=100)
      INTEGER FID(MXFLD)
      REAL RA(MXFLD),DEC(MXFLD),STA(MXFLD),STO(MXFLD)
      INTEGER FIL(MXFLD,6)
      INTEGER EXPT(MXFLD,6),CAD(MXFLD,6)
      CHARACTER*10 CFIL(MXFLD,6)
      CHARACTER*100 DFACGF
      LOGICAL MAKEO(MXFLD,6)
      CHARACTER*100 OPSOBO
      CHARACTER*100 OUTFILE
C
      INTEGER IFID,NFID
C
C Time of start of survey
C
      REAL*8 MJD0,MJD,MJDS
      PARAMETER (MJD0=49353.032079)
C
C Seconds per day
C
      REAL*8 SPD
      PARAMETER (SPD=24.0*60.0*60.0)
C
C Time to read out (readout + shutter = 3 seconds)
C
      REAL*8 RDT
      PARAMETER (RDT=3.0/SPD)
C
C Time to slew (5 seconds)
C
      REAL*8 SLT
      PARAMETER (SLT=5.0/SPD)
C
C Time to change filters (120 seconds)
C
      REAL*8 FCT
      PARAMETER (FCT=120.0/SPD)
C
C Filters as characters
C
      CHARACTER*1 FILTER(6)
      DATA FILTER/'u','g','r','i','z','y'/
      INTEGER NFIL(MXFLD)
C
C Standard exposure
C
      REAL STET
      PARAMETER (STET=34.0)
C
C Keeping track of the time
C
      INTEGER SMTM
      INTEGER INEXT,IFIL,IFN,ILAST,ILFN
      INTEGER NIGHT(MXFLD,6),LNIGHT
      REAL STOP(MXFLD,6)
C
      LOGICAL CLEAN
C
C OpsimObs output, labeled with L for Lynne
C
      REAL RAL,DECL,ALTL,AZL,AIRL,CLOUDL,RAWL,SEEL,SUNALL,SUNAZL,MOALL,
     &     MOAZL,MOPHL,MOPOL,SKYL,MAGLL
      INTEGER SNL,DOWL,FIDL
      CHARACTER*1 FILL,MAKL
      REAL*8 MJDL
      LOGICAL OKTMU(6)
C
C Long trailing blanks are a drag in output.
C
      CHARACTER FORM*6
      INTEGER ILL
C
C Keep track so that we never violate the 5 filters/night rule
C
      INTEGER FPN
      LOGICAL FTN(6)
C
C------------------------------------------------------------------
C
C Defaults
C
      MOONL = MOONLD
      MAKEUP = MAKEUPD
      AIRM = AIRMC
      NFID = 0
C
C Initial or cleanup
C
C     WRITE (6,'(A)') 'What are we doing?  Initial(I) or cleanup(C)?'
C     READ (5,'(A)') LINE
C     IF (LINE(1:1).EQ.'I') THEN
      IF (ACT(1:1).EQ.'I') THEN
        CLEAN = .FALSE.
      ELSEIF (ACT(1:1).EQ.'C') THEN
        CLEAN = .TRUE.
      ELSE
        WRITE (6,*) 'LSSTDFACG:  Invalid response: ',ACT
        STOP
      ENDIF

C Read input
C  Get file name
C
C      WRITE (6,'(A)') 'What Arbitrary Cadence file to process?'
C      READ (5,'(A)') DFACGF
      OPEN(UNIT=30,FILE=DFACGF,FORM='FORMATTED',STATUS='OLD')
C
C Read file
C
 10   CONTINUE
        READ (30,'(A)',END=11) LINE
        IF (LINE(1:1).EQ.'#') THEN
C Comment
          GOTO 10
        ELSEIF (LINE.EQ.' ') THEN
C Blank line
          GOTO 10
        ELSEIF (LINE(1:6).EQ.'MOONL:') THEN
C Moon Level to skip
          LINE = LINE(7:200)
          READ (LINE,*) MOONL
        ELSEIF (LINE(1:7).EQ.'MAKEUP:') THEN
C Makeup or not
          LINE = LINE(8:200)
          READ (LINE,*) MAKEUP
        ELSEIF (LINE(1:8).EQ.'AIRMASS:') THEN
C Makeup or not
          LINE = LINE(9:200)
          READ (LINE,*) AIRM
        ELSEIF (LINE(1:8).EQ.'FIELDID:') THEN
C A new field
          NFID = NFID + 1
          NFIL(NFID) = 0 
          IF (NFID.GT.MXFLD) THEN
            WRITE (6,*) 'LSSTDFACG: Too many fields.  Increase MXFLD'
            STOP
          ENDIF
          LINE = LINE(10:200)
          READ (LINE,*) FID(NFID)
        ELSEIF (LINE(1:8).EQ.'FIELDRA:') THEN
          LINE = LINE(10:200)
          READ (LINE,*) RA(NFID)
        ELSEIF (LINE(1:9).EQ.'FIELDDEC:') THEN
          LINE = LINE(11:200)
          READ (LINE,*) DEC(NFID)
        ELSEIF (LINE(1:11).EQ.'FIELDSTART:') THEN
          LINE = LINE(12:200)
          READ (LINE,*) STA(NFID)
        ELSEIF(LINE(1:10).EQ.'FIELDEND:') THEN
          LINE = LINE(11:200)
          READ (LINE,*) STO(NFID)
        ELSEIF (LINE(1:12).EQ.'FIELDFILTER:') THEN
          NFIL(NFID) = NFIL(NFID) + 1
          IF (NFIL(NFID) .GT.6) THEN
            WRITE (6,*) 'LSSTDFACG: There are only 6 filters?'
            STOP
          ENDIF
          LINE = LINE(13:200)
          READ (LINE,*) CFIL(NFID,NFIL(NFID))
        ELSEIF (LINE(1:9).EQ.'FIELDCAD:') THEN
          LINE = LINE(10:200)
          READ (LINE,*) CAD(NFID,NFIL(NFID))
        ELSEIF (LINE(1:9).EQ.'FIELDEXP:') THEN
          LINE = LINE(10:200)
          READ (LINE,*) EXPT(NFID,NFIL(NFID))
        ELSE
          WRITE (6,*) 'LSSTDFACG: I do not understand '
          WRITE (6,*) 'this in the input file?'
          WRITE (6,*) LINE
          STOP
        ENDIF
        GOTO 10
 11   CONTINUE
C End of input file
      CLOSE (30)

C Compute the initial nights, nexts, and stop times for the fields
C Also set makeup observation to false
      DO IFID = 1,NFID
        DO IFIL = 1,NFIL(IFID)
          NIGHT(IFID,IFIL) = NINT(MJD0)
          STOP(IFID,IFIL) = MJD0 + STA(IFID)*365.0 + STO(IFID)*365.0
          MAKEO(IFID,IFIL) = .FALSE.
        ENDDO
      ENDDO
C Check on make up
      MAKEU = .FALSE.
      IF (MAKEUP.EQ.1) THEN
        MAKEU = .TRUE.
      ENDIF
C Initial pass is to write out for input into OpsimObs
      IF (.NOT.CLEAN) THEN
C Define the output file
C        WRITE (6,'(A)') 'What should I call the output?'
C        READ (5,'(A)') OUTFILE
C Open output file
       OPEN(UNIT=31,FILE=OUTFILE,FORM='FORMATTED',STATUS='NEW')
C write a header
        WRITE (31,'(A,A)') '# RA(deg) DEC(deg) mjd ',
     &                     'filter exptime(sec) Fieldid'
        FPN = 0
        DO IFIL = 1,6
          FTN(IFIL) = .FALSE.
        ENDDO
C Start at time + 2 hours to avoid twilight with the first field
        MJDS = MJD0 - DFLOAT(NINT(MJD0)) + 2.0/24.0
C Initial observation at:
        LNIGHT = 0.0
        MJD = MJD0 + 2.0/24.0
        ILAST = 1
        ILFN = 1
 30     CONTINUE
C Find the smallest next
          INEXT = 0
          IFN = 0
          SMTM = 100000000
          DO IFID = 1,NFID
            DO IFIL = 1,NFIL(IFID)
              IF (NIGHT(IFID,IFIL).LE.SMTM) THEN
                INEXT = IFID
                IFN = IFIL
                SMTM = NIGHT(IFID,IFIL)
              ENDIF
            ENDDO
          ENDDO
          IF (INEXT.EQ.0) THEN
      WRITE (6,*) 'Something has gone wrong.  Cannot find next field'
            STOP
          ENDIF
C Check if over filter/night limit
          IF (FPN.GE.5) THEN
C Can we do this?, that is are we using a filter already in the camera?
            DO IFIL = 1,6
              IF (CFIL(INEXT,IFN).EQ.FILTER(IFIL)) THEN
                IF (FTN(IFIL)) GOTO 33
              ENDIF
            ENDDO
C Cannot do this observation.  Shift it forward one day.
            NIGHT(INEXT,IFN) = NIGHT(INEXT,IFN) + 1
C Still have to check if we are done
            GOTO 32
          ENDIF
 33       CONTINUE
C Time for this observation
          IF (NIGHT(INEXT,IFN).NE.LNIGHT) THEN
C New night
            MJD = DFLOAT(NIGHT(INEXT,IFN)) + MJDS
C Reset filter/night
            FPN = 0
            DO IFIL = 1,6
              FTN(IFIL) = .FALSE.
            ENDDO
          ELSE
            IF (IFN.NE.ILFN.AND.ILAST.EQ.INEXT) THEN
C Same night, same field, new filter
              MJD = RDT + FCT + 
     &              (DFLOAT(EXPT(ILAST,ILFN))/SPD) +
     &              MJD
            ELSE
C Same night, new field
              MJD = RDT + FCT + SLT +
     &              (DFLOAT(EXPT(ILAST,ILFN))/SPD) +
     &              MJD
            ENDIF
          ENDIF
C OK to write out? + an ordinary observation mark
          IF ((MJD.GE.STA(INEXT)).AND.(MJD.LE.STOP(INEXT,IFN))) THEN
            IF (.NOT.MAKEO(INEXT,IFN)) THEN
       WRITE (31,'(F15.10,1X,F15.10,1X,F20.10,1X,1A,1X,I4,1X,I5,1X,1A)')
     &        RA(INEXT),DEC(INEXT),MJD,CFIL(INEXT,IFN),
     &        EXPT(INEXT,IFN),FID(INEXT),'O'
            ELSE
C Make up observation with a make up mark
       WRITE (31,'(F15.10,1X,F15.10,1X,F20.10,1X,1A,1X,I4,1X,I5,1X,1A)') 
     &        RA(INEXT),DEC(INEXT),MJD,CFIL(INEXT,IFN),
     &        EXPT(INEXT,IFN),FID(INEXT),'M'
            ENDIF
C This observation
            LNIGHT = NIGHT(INEXT,IFN)
            ILAST = INEXT
            ILFN = IFN
C Increment filter/night
            DO IFIL = 1,6
              IF (CFIL(INEXT,IFN).EQ.FILTER(IFIL)) FTN(IFIL) = .TRUE.
            ENDDO
            FPN = 0
            DO IFIL = 1,6
              IF (FTN(IFIL)) FPN = FPN + 1
            ENDDO
          ENDIF
C Update night for this
          IF (.NOT.MAKEU) THEN
C No making up, just do observations on cadence
            NIGHT(INEXT,IFN) = NIGHT(INEXT,IFN) + CAD(INEXT,IFN)
          ELSE
            IF (CAD(INEXT,IFN).EQ.1) THEN
C Just do regular observations every night
              NIGHT(INEXT,IFN) = NIGHT(INEXT,IFN) + 1
            ELSE
              IF (.NOT.MAKEO(INEXT,IFN)) THEN
C Next observation is a make up
                NIGHT(INEXT,IFN) = NIGHT(INEXT,IFN) + 1 
                MAKEO(INEXT,IFN) = .TRUE.
              ELSE
C Next observation is a back on regular cadence, but this was a make up
                NIGHT(INEXT,IFN) = NIGHT(INEXT,IFN) + CAD(INEXT,IFN) - 1
                MAKEO(INEXT,IFN) = .FALSE.
              ENDIF
            ENDIF
          ENDIF
C Check if we are done.  Are all NIGHT greater than all STOP?
 32       CONTINUE
          DO IFID = 1,NFID
            DO IFIL = 1,NFIL(IFID)
              IF (DFLOAT(NIGHT(IFID,IFIL)).LT.STOP(IFID,IFIL)) THEN
                GOTO 30
              ENDIF
            ENDDO
          ENDDO
C Get here are and we are done
 31     CONTINUE
C Close the output file
        CLOSE(31)
      ELSE
C Clean up output of Opsim Ops.
C Open OpsimObs output.  Ask the name.
C        WRITE (6,'(A)') 'What OpsimObs output file to process?'
        READ (5,'(A)') OPSOBO
        OPEN(UNIT=30,FILE=OPSOBO,FORM='FORMATTED',STATUS='OLD')
C Define the output file
        WRITE (6,'(A)') 'What should I call the output?'
        READ (5,'(A)') OUTFILE
C Output processed ObsimObs file
        OPEN(UNIT=31,FILE=OUTFILE,FORM='FORMATTED',
     &       STATUS='NEW')
C Reading the OpsimObs file
 40     CONTINUE
          READ (30,'(A)',END=41) LINE
          DO 42 ILL = 200,1,-1
            IF (LINE(ILL:ILL).NE.' ') THEN
              WRITE (FORM,'(A1,I3,A2)') '(',ILL,'A)'
              GOTO 43
            ENDIF
 42       CONTINUE
 43       CONTINUE
          IF (LINE(1:1).EQ.'#') THEN
C Lynne's comment.  Pass to processed output
            WRITE (31,FORM) LINE(1:ILL)
            GOTO 40
          ELSE
C Read the data.
            READ (LINE,*) RAL,DECL,MJDL,FILL,ALTL,AZL,AIRL,CLOUDL,RAWL,
     &                    SEEL,SUNALL,SUNAZL,MOALL,MOAZL,MOPHL,MOPOL,
     &                    SKYL,MAGLL,SNL,DOWL,FIDL,MAKL
            IFIL = 0
            DO IFN = 1,6
              IF (FILL.EQ.FILTER(IFN)) IFIL = IFN
            ENDDO
            IF (MAKL.EQ.'O') OKTMU(IFIL) = .FALSE.
          ENDIF
C Do not write if down, OK to make up
          IF (DOWL.EQ.1) THEN
            WRITE (6,*) 'Skipping view at ',MJDL,'.  LSST down.'
            OKTMU(IFIL) = .TRUE.
            GOTO 40
          ENDIF
C Do not write if weathered out, OK to make up
          IF (CLOUDL.GT.CLOUDM) THEN
            WRITE (6,*) 'Skipping view at ',MJDL,'.  Clouds.'
            OKTMU(IFIL) = .TRUE.
            GOTO 40
          ENDIF
C Do not write if Moon is too large, OK to make up
          IF (MOPOL.GT.MOONL) THEN
            WRITE (6,*) 'Skipping view at ',MJDL,'.  Moon.'
            OKTMU(IFIL) = .TRUE.
            GOTO 40
          ENDIF
C Do not write if Airmass is too large, OK to make up
          IF (AIRL.GT.AIRM) THEN
            WRITE (6,*) 'Skipping view at ',MJDL,'.  Airmass.'
            OKTMU(IFIL) = .TRUE.
            GOTO 40
          ENDIF
C Handle make ups.  Are we doing that?
          IF (MAKEU) THEN
C Is this a make up observtion?
            IF (MAKL.EQ.'M') THEN
C If so did the last regular observation in this filter fail?  If so make it up.
C Otherwise do not write it out and go on.
              IF (.NOT.OKTMU(IFIL)) THEN
                WRITE (6,*) 'Making up view at ',MJDL,'.'
                GOTO 40
              ENDIF
            ENDIF
          ENDIF
C OK to write
          WRITE (31,FORM) LINE(1:ILL)
          GOTO 40
C Done
 41     CONTINUE
C Close the files
        CLOSE(30)
        CLOSE(31)
      ENDIF
C
      END

