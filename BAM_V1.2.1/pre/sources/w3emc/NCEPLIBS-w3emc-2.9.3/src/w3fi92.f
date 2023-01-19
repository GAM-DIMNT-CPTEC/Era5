C> @file
C> @brief Build 80-char on 295 grib queue descriptor.
C> @author Bill Cavanaugh @date 1991-06-21

C> Build 80 character queue descriptor using information
C> supplied by user, placing the completed queue descriptor in the
C> location specified by the user. (based on office note 295).
C>
C> @note This is a modified version of w3fi62() which adds the 'KWBX'
C> parameter. This value will now be added to bytes 35-38 for all grib
C> products. Queue desciptors for non-grib products will continue to be
C> generated by w3fi62().
C>
C> Program history log:
C> - Bill Cavanaugh 1991-06-21
C> - Bill Cavanaugh 1994-03-08 Modified to allow for bulletin sizes that
C> exceed 20000 bytes
C> - Ralph Jones 1994-04-28 Change for cray 64 bit word size and for ASCII
C> character set computers
C> - J. Smith 1995-10-16 Modified version of w3fi62() to add 'KWBX' to bytes
C> 35-38 of queue descriptor.
C> - Ralph Jones 1996-01-29 Preset ierr to zero.
C> - Boi Vuong 2002-10-15 Replaced function ichar with mova2i.
C>
C> @param[in] TTAAII First 6 characters of wmo header
C> @param[inout] KARY Integer array containing user information
C> - 1 = Day of month
C> - 2 = Hour of day
C> - 3 = Hour * 100 + minute
C> - 4 = Catalog number
C> - 5 = Number of 80 byte increments
C> - 6 = Number of bytes in last increment
C> - 7 = Total size of message WMO header + body of message in bytes (not
C> including queue descriptor)
C> @param[in] KWBX 4 characters, representing the fcst model that the bulletin
C> was derived from.
C> @param[out] LOC Location to receive queue descriptor.
C> @param[out] IERR Error return.
C>
C>
C> @note If total size is entered (kary(7)) then kary(5) and kary(6) will be calculated.
C> If kary(5) and kary(6) are provided then kary(7) will be ignored.
C>
C> @note Equivalence array loc to integer array so it starts on a word
C> boundary for sbyte subroutine.
C>
C> Error returns:
C> - IERR = 1 Total byte count and/or 80 byte increment count is missing. One
C> or the other is required to complete the queue descriptor.
C> - IERR = 2 Total size too small
C>
C> @author Bill Cavanaugh @date 1991-06-21
      SUBROUTINE W3FI92 (LOC,TTAAII,KARY,KWBX,IERR)
C
      INTEGER       IHOLD(2)
      INTEGER       KARY(7),IERR
C
      LOGICAL       IBM370
C
      CHARACTER*6   TTAAII,AHOLD
      CHARACTER*80  LOC
      CHARACTER*1   BLANK
      CHARACTER*4   KWBX
C
      EQUIVALENCE   (AHOLD,IHOLD)
C
      SAVE
C
C     BLANK WILL BE 40 HEX OR DECIMAL 64 ON AN IBM370 TYPE
C     COMPUTER, THIS IS THE EBCDIC CHARACTER SET.
C     BLANK WILL BE 20 HEX OR DECIMAL 32 ON A COMPUTER WITH THE
C     ASCII CHARACTER SET. THIS WILL BE USED TO TEST FOR CHARACTER
C     SETS TO FIND IBM370 TYPE COMPUTER.
C
      DATA  BLANK /' '/
C  ----------------------------------------------------------------
C
C     TEST FOR CRAY 64 BIT COMPUTER, LW = 8
C
      CALL W3FI01(LW)
C
C     TEST FOR EBCDIC CHARACTER SET
C
      IBM370 = .FALSE.
      IF (MOVA2I(BLANK).EQ.64) THEN
        IBM370 = .TRUE.
      END IF
C
      INOFST    = 0
C BYTES 1-16                   'QUEUE DESCRIPTOR'
      CALL SBYTE  (LOC,-656095772,INOFST,32)
      INOFST    = INOFST + 32
      CALL SBYTE  (LOC,-985611067,INOFST,32)
      INOFST    = INOFST + 32
      CALL SBYTE  (LOC,-490481207,INOFST,32)
      INOFST    = INOFST + 32
      CALL SBYTE  (LOC,-672934183,INOFST,32)
      INOFST    = INOFST + 32
C BYTES 17-20                  INTEGER ZEROES
      CALL SBYTE (LOC,0,INOFST,32)
      INOFST    = INOFST + 32
C                              IF TOTAL COUNT IS INCLUDED
C                              THEN WILL DETERMINE THE NUMBER OF
C                              80 BYTE INCREMENTS AND WILL DETERMINE
C                              THE NUMBER OF BYTES IN THE LAST INCREMENT
      IERR = 0
      IF (KARY(7).NE.0) THEN
          IF (KARY(7).LT.35) THEN
C             PRINT *,'LESS THAN MINIMUM SIZE'
              IERR   = 2
              RETURN
          END IF
          KARY(5)    = KARY(7) / 80
          KARY(6)    = MOD(KARY(7),80)
          IF (KARY(6).EQ.0) THEN
              KARY(6)    = 80
          ELSE
              KARY(5)    = KARY(5) + 1
          END IF
      ELSE
          IF (KARY(5).LT.1) THEN
              IERR     = 1
              RETURN
          END IF
      END IF
C BYTE  21-22                  NR OF 80 BYTE INCREMENTS
      CALL SBYTE (LOC,KARY(5),INOFST,16)
      INOFST    = INOFST + 16
C BYTE  23                     NR OF BYTES IN LAST INCREMENT
      CALL SBYTE (LOC,KARY(6),INOFST,8)
      INOFST    = INOFST + 8
C BYTES 24-28                  INTEGER ZEROES
      CALL SBYTE (LOC,0,INOFST,32)
      INOFST    = INOFST + 32
      CALL SBYTE (LOC,0,INOFST,8)
      INOFST    = INOFST + 8
C BYTES 29-34                  6 CHAR BULLETIN NAME TTAAII
      LOC(29:34) = TTAAII(1:6)
C
C     IF ON ASCII COMPUTER, CONVERT LAST 6 CHARACTERS TO EBCDIC
C
      IF (.NOT.IBM370) CALL W3AI39(LOC(29:29),6)
C
      INOFST    = INOFST + 48
C BYTES 35-38                  KWBX
C
      LOC(35:38) = KWBX(1:4)
C
C     IF ON ASCII COMPUTER, CONVERT LAST 4 CHARACTERS TO EBCDIC
C
      IF (.NOT.IBM370) CALL W3AI39(LOC(35:35),4)
      INOFST    = INOFST + 32
C BYTES 39-40                  HR/MIN TIME OF BULLETIN CREATION
C                                    TWO BYTES AS 4 BIT BCD
      KA        = KARY(3) / 1000
      KB        = MOD(KARY(3),1000) / 100
      KC        = MOD(KARY(3),100) / 10
      KD        = MOD(KARY(3),10)
      CALL SBYTE (LOC,KA,INOFST,4)
      INOFST    = INOFST + 4
      CALL SBYTE (LOC,KB,INOFST,4)
      INOFST    = INOFST + 4
      CALL SBYTE (LOC,KC,INOFST,4)
      INOFST    = INOFST + 4
      CALL SBYTE (LOC,KD,INOFST,4)
      INOFST    = INOFST + 4
C BYTES 41-45                  CATALOG NUMBER ELSE (SET TO 55555)
      IF (KARY(4).GE.1.AND.KARY(4).LE.99999) THEN
          CALL W3AI15 (KARY(4),IHOLD,1,8,'-')
          IF (LW.EQ.4) THEN
            CALL SBYTE (LOC,IHOLD(1),INOFST,8)
            INOFST    = INOFST + 8
            CALL SBYTE (LOC,IHOLD(2),INOFST,32)
            INOFST    = INOFST + 32
C
C         ON CRAY 64 BIT COMPUTER
C
          ELSE
            CALL SBYTE (LOC,IHOLD,INOFST,40)
            INOFST    = INOFST + 40
          END IF
C
C   IF ON ASCII COMPUTER, CONVERT LAST 5 CHARACTERS TO EBCDIC
C
          IF (.NOT.IBM370) CALL W3AI39(LOC(41:41),5)
      ELSE
          CALL SBYTE (LOC,-168430091,INOFST,32)
          INOFST    = INOFST + 32
          CALL SBYTE (LOC,245,INOFST,8)
          INOFST    = INOFST + 8
      END IF
C BYTES 46-80                  INTEGER ZEROES
      DO 4676 I = 1, 8
          CALL SBYTE (LOC,0,INOFST,32)
          INOFST    = INOFST + 32
 4676 CONTINUE
          CALL SBYTE (LOC,0,INOFST,24)
      RETURN
      END