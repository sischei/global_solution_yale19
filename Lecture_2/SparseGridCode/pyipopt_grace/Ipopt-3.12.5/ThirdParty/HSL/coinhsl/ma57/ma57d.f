C COPYRIGHT (c) 1999 CCLRC Council for the Central Laboratory
*                    of the Research Councils
C
C Version 3.7.0 (5th May 2011)
C See ChangeLog for version history
C

      SUBROUTINE MA57ID(CNTL, ICNTL)
C
C****************************************************************
C
C  Purpose
C  =======
C
C  The entries of the arrays CNTL and ICNTL control the action of
C  MA57. Default values for the entries are set in this routine.
C
      DOUBLE PRECISION    CNTL(5)
      INTEGER             ICNTL(20)
C
C  Parameters
C  ==========
C
C  The entries of the arrays CNTL and ICNTL control the action of
C     MA57. Default values are set by MA57ID.
C
C  CNTL(1) has default value 0.01 and is used for threshold pivoting.
C     Values less than 0.0 will be treated as 0.0 and values greater
C     than 0.5 as 0.5. Values near 0.0 may perhaps give faster
C     factorization times and less entries in the factors but may result
C     in a less stable factorization.  This parameter is only accessed
C     if ICNTL(7) is equal to 1.
C
C  CNTL(2)  has default value 1.0D-20. MA57B/BD will treat any pivot
C     whose modulus is less than CNTL(2) as zero. If ICNTL(16) = 1,
C     then blocks of entries less than CNTL(2) can be discarded during
C     the factorization and the corresponding pivots are placed at the
C     end of the ordering. In this case, a normal value for CNTL(2)
C     could be 1.0D-12.
C
C CNTL(3) has default value 0.5.  It is used by MA57D/DD to monitor
C     the convergence of the iterative refinement.  If the norm of
C     the scaled residuals does not decrease by a factor of at least
C     CNTL(3), convergence is deemed to be too slow and MA57D/DD
C     terminates with INFO(1) set to -8.

C CNTL(4) has default value 0.0. It is used by MA57B/BD to control
C     the static pivoting. If CNTL(4) is greater than
C     zero, then small pivots may be replaced by entries of value
C     CNTL(4) so that the factorization could be inaccurate. If CNTL(5)
C     is also greater than zero, then CNTL(4) is treated as zero (that
C     is uneliminated variables are delayed) until CNTL(5)*N 
C     fully summed variables have been delayed.  If static pivots are
C     used, it is recommended that iterative refinement is used when
C     computing the solution.
C
C  CNTL(5) has default value 0.0.  Static pivoting is invoked if CNTL(4)
C     is greater than 0.0 and the accumulated number of delayed pivots
C     exceeds CNTL(5)*N.
C
C  ICNTL(1) has default value 6.
C     It is the output stream for error messages. If it is negative,
C     these messages are suppressed.
C
C  ICNTL(2) has default value 6.
C     It is the output stream for warning messages. If it is negative,
C     these messages are suppressed.
C
C  ICNTL(3) has default value 6. It is the output stream for monitoring
C     printing. If it is negative, these messages are suppressed.
C
C  ICNTL(4) has default value -1.
C     It is the output stream for the printing of statistics.
C     If it is negative, the statistics are not printed.
C
C  ICNTL(5) is used by MA57 to control printing of error,
C     warning, and monitoring messages. It has default value 2.
C     Possible values are:
C
C    <1       No messages output.
C     1       Only error messages printed.
C     2       Errors and warnings printed.
C     3       Errors and warnings and terse monitoring
C             (only first ten entries of arrays printed).
C    >3       Errors and warnings and all information
C             on input and output parameters printed.
C
C  ICNTL(6) has default value 5 and must be set by the user to 1
C     if the pivot order in KEEP is to be used by MA57AD. For any
C     other value of ICNTL(6), a suitable pivot order will be chosen
C     automatically.  The choices available so far are:
C  ICNTL(6) = 0     AMD using MC47 (with dense row detection disabled)
C  ICNTL(6) = 2     AMD using MC47 (this was previously the MC50 option)
C  ICNTL(6) = 3     MA27 minimum degree ordering
C  ICNTL(6) = 4     METIS_NODEND ordering from MeTiS package.
C  ICNTL(6) = 5     Ordering chosen depending on matrix characteristics.
C                   At the moment choices are MC47 or METIS.
C                   INFO(36) is set to ordering used.
C  ICNTL(6) > 5     At the moment this is treated as 5 (the default).
C
C  ICNTL(7) is used by MA57BD to control numerical pivoting.  It
C     has default value 1.  Values out of range cause an error return
C     with INFO(1) equal to -10.
C     Possible values are:
C
C  1   Numerical pivoting is performed using the threshold value in
C      CNTL(1).
C
C  2   No pivoting will be performed and an error exit will occur
C      immediately a sign change or a zero is detected among the pivots.
C      This is suitable for cases when A is thought to be definite
C      and is likely to decrease the factorization time while still
C      providing a stable decomposition.
C
C  3   No pivoting will be performed and an error exit will occur if a
C      zero pivot is detected. This is likely to decrease the
C      factorization time, but may be unstable if there is a sign
C      change among the pivots.
C
C  4   No pivoting will be performed but the matrix will be altered
C      so that all pivots are of the same sign.
C
C  ICNTL(8) has default value 0. If MA57BD is called with ICNTL(8) NE 0,
C      then the factorization will discard factors and try
C      to continue the factorization to determine the amount of space
C      needed for a successful factorization.  In this case, a
C      factorization will not have been produced.  If the default value
C      of 0 is used and the factorization stops because of lack of
C      space, the user should reallocate the real or integer space for
C      FACT or IFACT, respectively and reset LFACT or LIFACT
C      appropriately, using MA57ED before recalling MA57BD.
C
C  ICNTL(9) has default value 10.  It corresponds to the maximum number
C      of steps of iterative refinement.
C
C  ICNTL(10) has default value 0. A positive value will return the
c      infinity norm of the input matrix, the computed solution, and
C      the scaled residual in RINFO(5) to RINFO(7), respectively,
C      a backward error estimate in RINFO(8) and RINFO(9), and an
C      estimate of the forward error in RINFO(10).  If ICNTL(10) is
C      negative or zero no estimates are returned.
C
C  ICNTL(11) The block size to be used by the Level 3 BLAS.
C
C  ICNTL(12) Two nodes of the assembly tree are merged only if both
C      involve less than ICNTL(12) eliminations.
C
C  ICNTL(13) Threshold on number of rows in a block for using BLAS2 in
C      MA57CD.
C
C  ICNTL(14) Threshold on number of entries in a row for declaring row
C      full on a call to MA57H/HD. Set as percentage of N.
C      So 100 means row must be full.
C
C  ICNTL(15) should be set to 1 (the default) if MC64 scaling is
C     requested.
C
C  ICNTL(16) should be set to 1 if "small" entries are to be removed
C     removed from the frontal matrices.  The default is 0.
C
C  ICNTL(17) to ICNTL(20) are set to zero by MA57ID but are not
C     currently used by MA57.
C
C Local variables
      INTEGER I
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C===============================================
C Default values for variables in control arrays
C===============================================
C     Threshold value for pivoting
      CNTL(1)   = 0.01D0
C Test for zero pivot
      CNTL(2)   = 1.0D-20
C Iterative refinement convergence
      CNTL(3)   = 0.5D0
C Static pivoting control
      CNTL(4) = ZERO
C Control to allow some delayed pivots
      CNTL(5) = ZERO
C     Printing control
      ICNTL(1)  = 6
      ICNTL(2)  = 6
      ICNTL(3)  = 6
      ICNTL(4)  = -1
      ICNTL(5)  = 2
C     Provide pivot order (1=NO)
C     Set to make automatic choice between METIS and MC47
      ICNTL(6)  = 5
C     Pivoting control
      ICNTL(7)  = 1
C     Restart facility
      ICNTL(8)  = 0
C     IR steps
      ICNTL(9)  = 10
C     Error estimates
      ICNTL(10) = 0
C     Blocking for Level 3 BLAS
      ICNTL(11) = 16
C     Node amalgamation parameter (NEMIN)
      ICNTL(12) = 16
C     Switch for use of Level 2 BLAS in solve
      ICNTL(13) = 10
C Flag to indicate threshold will be set to N
      ICNTL(14) = 100
C Flag to invoke MC64 scaling (0 off, 1 on)
      ICNTL(15) = 1
C Flag to invoke dropping small entries from front
C Default is not to drop (set to 1 to drop)
      ICNTL(16) = 0

C Set unused parameters
      DO 110 I=17,20
        ICNTL(I) = 0
  110 CONTINUE

      RETURN
      END


      SUBROUTINE MA57AD(N,NE,IRN,JCN,LKEEP,KEEP,IWORK,ICNTL,INFO,RINFO)
C This subroutine is a user-callable driver for the analysis phase of
C     MA57. It performs an ordering, a symbolic factorization and
C     computes information for the numerical factorization by MA57B/BD.
      INTEGER N,NE,IRN(NE),JCN(NE),IWORK(5*N),LKEEP,KEEP(LKEEP),
     *        ICNTL(20),INFO(40)
      DOUBLE PRECISION RINFO(20)
C
C  N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A.  It is not altered by the subroutine.
C    Restriction: N > 0.
C
C  NE is an INTEGER variable which must be set by the user to the
C    number of entries being input.  It is not altered by the
C    subroutine. Restriction: NE >= 0.
C
C  IRN and JCN are INTEGER  arrays of length NE. The user
C    must set them so that each off-diagonal nonzero $a_{ij}$ is
C    represented by IRN(k)=i and JCN(k)=j or by IRN(k)=j
C    and JCN(k)=i.  Multiple entries are allowed and any with
C    IRN(k) or JCN(k) out of range are ignored. These arrays will
C    be unaltered by the subroutine.
C
C  IWORK is an INTEGER  array of length 5*N. This need not be set
C    by the user and is used as workspace by MA57AD.
C
C  LKEEP is an INTEGER that must be set to length of array KEEP.
C    Restriction: LKEEP >= 5*N+NE+MAX(N,NE)+42
C
C  KEEP is an INTEGER  array of length LKEEP. It need not be set
C    by the user and must be preserved between a call to MA57AD
C    and subsequent calls to MA57BD and MA57CD.
C    If the user wishes to input
C    the pivot sequence, the position of variable {i} in the pivot
C    order should be placed in KEEP(i), i = 1, 2,..., n and ICNTL(6)
C    should be set to 1.  The subroutine may replace the given order
C    by another that gives the same fill-in pattern and
C    virtually identical numerical results.
C
C ICNTL is an INTEGER array of length 10
C    that contains control parameters and must be set by the user.
C    Default values for the components may be set by a call to
C    MA57I/ID. Details of the control parameters are given in MA57I/ID.
C
C INFO is an INTEGER array of length 40 that need not be set by the
C    user.  On return from MA57AD, a value of zero for INFO(1)
C    indicates that the subroutine has performed successfully.

C RINFO is a REAL (DOUBLE_PRECISION in the D version) array of length 20
C    that need not be set by the
C    user.  This array supplies information on the execution of MA57AD.
C
C**** Still to be updated
C    INFO(1):
C       0  Successful entry.
C      -1 N < 1
C      -2 NE < 0.
C      -9 Invalid permutation supplied in KEEP.
C     -15 LKEEP < 5*N+NE+MAX(N,NE)+42
C      +1 One or more indices out of range.
C      +2 One or more duplicate entries found.
C      +3. Combination of warnings +1 and +2.
C    INFO(2):
C        if INFO(1) = -1, the value input for N.
C        if INFO(1) = -2, the value input for NE.
C        if INFO(1) = -9, index at which error first detected.
C        if INFO(1) = -15, the value input for LKEEP.
C    INFO(3) Number of entries with out-of-range indices.
C    INFO(4) Number of off-diagonal duplicate entries.
C    INFO(5) Forecast number of reals to hold the factorization.
C    INFO(6) Forecast number of integers to hold the factorization.
C    INFO(7) Forecast maximum front size.
C    INFO(8) Number of nodes in the assembly tree.
C    INFO(9) Minimum size for LA of MA57BD (without compress).
C    INFO(10) Minimum size for LIFACT of MA57BD (without compress).
C    INFO(11) Minimum size for LA of MA57BD (with compress).
C    INFO(12) Minimum size for LIFACT of MA57BD (with compress).
C    INFO(13) Number of compresses.
C    INFO(14:40) Not used.

C Procedures
      INTRINSIC MIN
      EXTERNAL MA57GD,MC47ID,MC47BD,MA57VD,MA57HD,MA57JD,MA57KD,
     *         MA57LD,MA57MD,MA57ND
C MA57GD Expand representation to whole matrix and sort.
C MC47BD Approximate Minimum Degree avoiding problems with dense rows.
C MA57VD Is same as MA27GD. Sort for using MA27HD/MA57HD.
C MA57HD Is same as MA27HD. Minimum degree ordering from MA27.
C MA57JD Sort to upper triangular form (pivot sequence given).
C MA57KD Construct tree pointers given output from MA57JD.
C MA57LD Depth-first search of tree.
C MA57MD Construct map.
C MA57ND Calculate storage and operation counts.

C Local variables
      INTEGER I,IL,IN,IPE,IRNPRM,COUNT,FILS,FRERE,HOLD,IFCT,INVP,IPS,
     +        IW,IWFR,K,LDIAG,LP,LW,LROW,MAP,EXPNE,
     +        MP,NCMPA,NEMIN,NODE,NST,NSTEPS,NV,PERM,
     +        IW1,IW2,IW3,IW4,IW5,NSTK,ND,NELIM,NZE,ALENB,
     +        J,JJ,J1,J2,SIZE22,OXO
C Local variables for MA27 ordering call
      INTEGER IF27H
C Local variables for MeTiS call
      INTEGER METOPT(8),METFTN,ICNTL6,INF47(10),ICNT47(10)
      DOUBLE PRECISION ZERO,THRESH,AVNUM,MC47FI,RINF47(10)
      PARAMETER (ZERO=0.0D0)

C I        Temporary DO index.
C IPE      Displacement in array KEEP for array IPE of MA57GD,
C          MA57HD, MA57JD, and MA57KD.
C COUNT    Displacement in array KEEP for array COUNT of MA57GD,
C          MA57HD, MA57JD, MA57KD, MA57LD, and MA57ND.
C IWFR     First unused location in IFCT(1:LW).
C K        Temporary variable.
C LDIAG    Control for amount of information output.
C LP       Stream number for error messages.
C LW       Length of IFACT when four arrays of length N+1 are excluded.
C LROW     Subscript in array KEEP for array LAST of MA57KD and
C          MA57LD and LROW of MA57MD.
C MAP      Subscript in array KEEP for MAP array.
C MP       Stream number for monitoring.
C NODE     Subscript in array KEEP for array FLAG of MA57KD, and
C          NODE of MA57LD.
C NSTEPS   Number of nodes in the assembly tree.
C NV       Displacement in array IFACT for array NV of MA57GD and
C          MA57HD, PERM of MA57JD, IPR of MA57KD, and NE of MA57LD.
C PERM     Subscript in array KEEP for array PERM of MA57LD/ND.
C SP       Stream number for statistics.
C SIZES    subscript in array KEEP.  On exit,
C          KEEP(SIZES) is number of faulty entries.
C          KEEP(SIZES+1) is the number of nodes in the assembly tree.


C Set local print variables
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)

C Initialize information array.
      DO 10 I = 1,40
        INFO(I) = 0
   10 CONTINUE
      DO 11 I = 1,20
        RINFO(I) = ZERO
   11 CONTINUE

C Check N, NE, and LKEEP for obvious errors
      IF (N.LT.1)  GO TO 20
      IF (NE.LT.0) GO TO 30
      IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40

      IF (ICNTL(6).EQ.1) THEN
C Check permutation array
        DO 12 I = 1,N
          IWORK(I) = 0
   12   CONTINUE
        DO 14 I=1,N
          K = KEEP(I)
          IF (K.LE.0 .OR. K.GT.N) GO TO 80
          IF (IWORK(K).NE.0) GO TO 80
          IWORK(K) = I
   14   CONTINUE
      ENDIF

C If requested, print input variables.
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE(MP,99980) N,NE,(ICNTL(I),I=1,7),ICNTL(12),ICNTL(15)
99980 FORMAT (//'Entering analysis phase (MA57AD) with ...'/
     1      'N         Order of matrix                     =',I12/
     2      'NE        Number of entries                   =',I12/
     6      'ICNTL(1)  Stream for errors                   =',I12/
     7      ' --- (2)  Stream for warnings                 =',I12/
     8      ' --- (3)  Stream for monitoring               =',I12/
     9      ' --- (4)  Stream for statistics               =',I12/
     1      ' --- (5)  Level of diagnostic printing        =',I12/
     2      ' --- (6)  Flag for input pivot order          =',I12/
     2      ' --- (7)  Numerical pivoting control (st est) =',I12/
     2      ' --- (12) Node amalgamation parameter         =',I12/
     2      ' --- (15) Scaling control (storage estimate)  =',I12)
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(3(I6,A,2I8,A)))') ' Matrix entries:',
     +        (I,': (',IRN(I),JCN(I),')',I=1,K)
        IF (K.LT.NE) WRITE (MP,'(A)') '     . . .'

        IF (ICNTL(6).EQ.1) THEN
C Print out permutation array.
          K = MIN(10,N)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(A,10I6:/(7X,10I6))') ' KEEP =', (KEEP(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(7X,A)') '     . . .'
        END IF

      END IF

C Partition IWORK
      IW1 = 1
      IW2 = IW1 + N
      IW3 = IW2 + N
      IW4 = IW3 + N
      IW5 = IW4 + N
      FILS  = IW1
      FRERE = IW2
      ND    = IW3
      NELIM = IW4
      NV    = IW5

C Partition KEEP
      PERM = 1
      NSTEPS = PERM + N
      EXPNE  = NSTEPS + 1
      HOLD   = EXPNE + 1
      LROW = HOLD + 40
      NODE = LROW + N
      NSTK = NODE + N
      MAP  = NSTK + N
      IRNPRM = MAP + MAX(N,NE)
      INVP  = NODE
      IW    = NODE
      IPE   = LROW
      IFCT  = MAP
C This is set for MA57VD/MA57HD (MA27 orderings) only to allow more
C     space for expanded factors
      IF27H = NODE
      IPS   = MAP
      COUNT = NSTK

C Set HOLD(1)
      KEEP(HOLD) = 0

C Sort and order ... generate tree

C Set local value for ICNTL(6)
      ICNTL6 = ICNTL(6)
      IF (ICNTL(6).GT.5) ICNTL6 = 5

      IF (ICNTL6.EQ.4 .OR. ICNTL6.EQ.5) THEN
C MeTiS ordering requested.  Use dummy call to see if it has been
C     installed.
C Set flag for Fortran-style numbering of arrays
        METFTN    = 1
C Set default values for parameters.
        METOPT(1) = 0
C Dummy call with 1 by 1 matrix
        KEEP(IPE)   = 1
        KEEP(IPE+1) = 2
        KEEP(IFCT)  = 1
        CALL METIS_NODEND(1,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                    KEEP(NSTK),KEEP(PERM))
C Flag set if dummy code for METIS_NODEND has been used
        IF (KEEP(PERM).EQ.-1) THEN
          IF (ICNTL6 .EQ. 4) GO TO 90
C Reset ICNTL6 to use MC47
          ICNTL6 = 2
        ENDIF
      ENDIF
      IF (ICNTL6.NE.1) THEN
C Ordering is to be calculated by program

        CALL MC47ID(ICNT47)

        IF (ICNTL6 .NE. 3) THEN
C ELSE clause (ICNTL6.EQ.3) for MA27 ordering
C MC47 (with dense row detection disabled) used if ICNTL6 equal to 0.
C MC47 used if ICNTL6 equal to 2.
C METIS used if ICNTL6 equal to 4.
C Automatic choice of METIS or MC47 if ICNTL6 equal to 5.

C Sort matrix to obtain complete pattern (upper and lower triangle)
C     but omitting diagonals, duplicates, and out-of-range entries.
C     On exit, sorted matrix is given by IPE(pointers), IFCT (indices),
C     and COUNT (lengths).  IWFR is position after last index in IFCT.
          CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +                KEEP(IW),IWFR,ICNTL,INFO)

          IF (ICNTL6.EQ.5) THEN
C Calculate matrix statistics to determine ordering
            IF (ICNTL(7).EQ.2) THEN
C Action if positive definite option has been chosen.
              AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
              IF (N.GE.50000) THEN
                ICNTL6 = 4
C               IF (AVNUM.LE.6.0) ICNTL6 = 2
                GO TO 97
              ENDIF
              IF (N.LE.30000) THEN
                ICNTL6 = 2
                IF (AVNUM.GT.100.0) ICNTL6 = 4
                GO TO 97
              ENDIF
              IF (N.GT.30000 .AND. N.LT.50000) THEN
                IF (AVNUM.GT.46.0) THEN
                  ICNTL6 = 4
                ELSE
                  ICNTL6 = 2
                ENDIF
                GO TO 97
              ENDIF
            ELSE
C Matrix has not been declared positive definite.
              AVNUM = FLOAT(IWFR+N-1)/FLOAT(N)
C Set flag for detection of OXO matrix
              OXO = 0
C Check for KKT and OXO.
C Calculate size of possible  trailing block
              J2 = IWFR - 1
              SIZE22 = 0
              DO 100 J = N,1,-1
                J1 = KEEP(IPE+J-1)
C Note that MA57GD does not sort within order in columns
                DO  99 JJ = J1,J2
                  IF (KEEP(IFCT+JJ-1).GT.J) GO TO 101
   99           CONTINUE
                SIZE22 = SIZE22 + 1
                J2 = J1-1
  100         CONTINUE
  101         IF (SIZE22 .GT. 0) THEN
C Check to see if there are no entries in (1,1) block.
                DO 98 I = 1,NE
                  IF (IRN(I) .LE. N-SIZE22
     *          .AND. JCN(I) .LE. N-SIZE22) THEN
                      AVNUM = FLOAT(IWFR+N-SIZE22-1)/FLOAT(N)
                      GO TO 96
                  ENDIF
   98           CONTINUE
C The (1,1) block is zero.
                OXO = 1
                AVNUM = FLOAT(IWFR-1)/FLOAT(N)
              ENDIF
   96         IF (N .GE. 100000) THEN
                IF (AVNUM.GT.5.42) THEN
                  ICNTL6 = 4
                ELSE
                  ICNTL6 = 2
                ENDIF
                GO TO 97
              ENDIF
C Logic for OXO matrices
              IF (OXO.EQ.1) THEN
                IF (FLOAT(N-SIZE22)/FLOAT(SIZE22) .GT .1.8D0) THEN
                  ICNTL6 = 2
                ELSE
                  ICNTL6 = 4
                ENDIF
                GO TO 97
              ENDIF
C We can try further simple logic here ... then ...
C Call MC47 to test whether fill-in projected is large
              LW = LKEEP-IFCT+1
              CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
     +                    KEEP(IFCT),IWORK(NV),
     +                    KEEP(INVP),KEEP(PERM),IWORK(IW1),
     +                    IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +                    ICNT47,INF47,RINF47)
              INFO(13) = INF47(2)
              ICNTL6 = 2
              NEMIN    = ICNTL(12)
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
              NST = KEEP(NSTEPS)
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))
              KEEP(EXPNE) = IWORK(IW5)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)
C Check relative fill-in of MC47
              IF (FLOAT(INFO(5))/FLOAT(NE) .LT. 10.0) THEN
C We will run with the MC47 ordering
                GO TO 93
              ELSE
C Save value of relative fill-in for testing against METIS
                MC47FI = FLOAT(INFO(5))/FLOAT(NE)
C Must test METIS ordering now ... ugh
        CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +              KEEP(IW),IWFR,ICNTL,INFO)
                KEEP(IPE+N) = IWFR
                METFTN    = 1
                METOPT(1) = 0
                IF (N.LT.50) GO TO 92
                DO 91 I = 1,N
                  IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
                    METOPT(1) = 1
                    METOPT(2) = 3
                    METOPT(3) = 1
                    METOPT(4) = 2
                    METOPT(5) = 0
                    METOPT(6) = 1
                    METOPT(7) = 200
                    METOPT(8) = 1
                    GO TO 92
                  ENDIF
   91           CONTINUE
   92     CALL METIS_NODEND(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                      KEEP(NSTK),KEEP(PERM))
        CALL MA57JD(N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
     +              KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)
                LW = LKEEP - IFCT + 1
        CALL MA57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
     +              KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
                INFO(13) = NCMPA
                NEMIN = ICNTL(12)
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
                NST = KEEP(NSTEPS)
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))
                KEEP(EXPNE) = IWORK(IW5)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)
                IF (FLOAT(INFO(5))/FLOAT(NE).LT.MC47FI) THEN
                  ICNTL6 = 4
                  GO TO 93
                ELSE
C Double groan  ... we will run with MC47 after all
                  ICNTL6=2
C KEEP(IPE) has been corrupted must reset it.
        CALL MA57GD(N,NE,IRN,JCN,KEEP(IFCT),KEEP(IPE),KEEP(COUNT),
     +              KEEP(IW),IWFR,ICNTL,INFO)
                  GO TO 97
                ENDIF
C End of METIS check
              ENDIF
C End of indef case calculation
            ENDIF
C End of logic for ICNTL6 = 5
          ENDIF

   97     IF (ICNTL6.EQ.4) THEN
C Set last pointer in IPE
            KEEP(IPE+N) = IWFR
C Use MeTiS ordering
C Set flag for Fortran-style numbering of arrays
            METFTN    = 1
C This would use only defaults
            METOPT(1) = 0
C Set options for METIS, particularly one for dense columns
C First determine if there are any dense columns
            IF (N.LT.50) GO TO 103
            DO 102 I = 1,N
              IF ((KEEP(IPE+I)-KEEP(IPE+I-1)) .GT. N/10) THEN
C The rest are set to default values
                METOPT(1) = 1
                METOPT(2) = 3
                METOPT(3) = 1
                METOPT(4) = 2
                METOPT(5) = 0
                METOPT(6) = 1
                METOPT(7) = 200
                METOPT(8) = 1
                GO TO 103
              ENDIF
  102       CONTINUE
  103     CALL METIS_NODEND(N,KEEP(IPE),KEEP(IFCT),METFTN,METOPT,
     *                      KEEP(NSTK),KEEP(PERM))
            GO TO 111
          ENDIF


C Obtain ordering using approximate minimum degree ordering.
C Input IPE,IWFR,COUNT,IFCT as from MA57GD.
C Output
C     IPE (- father pointer .. if NV(I) > 0, - subordinate variable
C     pointer if NV(I) = 0)
C     NV(I) for subordinate variables of supervariable, otherwise
C     degree when eliminated.
C     IWFR is set to length required by MC47BD if no compresses.
C     COUNT, IFCT undefined on exit.
C     INVP is inverse permutation and PERM is permutation
C Length of LW set to maximum to avoid compresses in MC47B/BD
          LW = LKEEP-IFCT+1
C ICNTL6 =  0.  MC47 uses code for dealing with dense rows disabled
C In HSL 2002 it was only used in the F90 version.
C ICNTL6 =  2.  MC47 implements code for dealing with dense rows
          IF (ICNTL6 .EQ. 0) ICNT47(4) = -1
          CALL MC47BD(N,LW,KEEP(IPE),IWFR,KEEP(COUNT),
     +                KEEP(IFCT),IWORK(NV),
     +                KEEP(INVP),KEEP(PERM),IWORK(IW1),
     +                IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +                ICNT47,INF47,RINF47)
          INFO(13) = INF47(2)

        ELSE
C End of ICNTL6 .NE. 3

C MA27 ordering being used.  Must insert row lengths in KEEP(IFCT)
C Length of LW set to maximum to avoid compresses.
          LW = LKEEP-IF27H+1
        CALL MA57VD(N,NE,IRN,JCN,KEEP(IF27H),LW,KEEP(IPE),IWORK(IW1),
     *              IWORK(IW2),IWFR,ICNTL,INFO)
C Analyse using minimum degree ordering
          THRESH = FLOAT(ICNTL(14))/100.0
        CALL MA57HD(N,KEEP(IPE),KEEP(IF27H),LW,IWFR,IWORK(NV),
     *              IWORK(IW1),IWORK(IW2),IWORK(IW3),IWORK(IW4),
     +              2139062143,INFO(13),THRESH)
C Set IPE correctly
          DO 110 I = 1,N
            IF (IWORK(NV+I-1).NE.0) GO TO 110
            IN = I
  105       IL = IN
            IN = - KEEP(IPE+IL-1)
            IF (IWORK(NV+IN-1).EQ.0) GO TO 105
C Make subordinate node point to principal node
            KEEP(IPE+I-1) = -IN
  110     CONTINUE
        ENDIF

C End of block for generating ordering
      ENDIF

  111 IF (ICNTL6.EQ.1 .OR. ICNTL6.EQ.4) THEN
C If we have generated ordering using MeTiS then we need to feed
C     permutation as if it were coming from the user as we do not
C     have a tight coupling to MeTiS as for other orderings.

C Sort using given order
        CALL MA57JD(N,NE,IRN,JCN,KEEP(PERM),KEEP(IFCT),KEEP(IPE),
     +              KEEP(COUNT),IWORK(IW1),IWFR,ICNTL,INFO)

C Generating tree using given ordering
C Input:  N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM)
C Output:  KEEP(IPE),IWORK(NV)
        LW = LKEEP - IFCT + 1
        CALL MA57KD(N,KEEP(IPE),KEEP(IFCT),LW,IWFR,KEEP(PERM),
     +              KEEP(INVP),IWORK(NV),IWORK(IW1),NCMPA)
        INFO(13) = NCMPA

      END IF


C Perform depth-first search of assembly tree
C Set NEMIN
      NEMIN = ICNTL(12)
C Input  IPE,NV,NEMIN
C Output
C     IPE .. father and younger brother pointer
C     NV  .. unchanged
C     NE/NSTK/ND defined for nodes of tree
C     PERM
C     IPS(I) position of node I in order
C     LROW(I) is size of frontal matrix at node I
      CALL MA57LD(N,KEEP(IPE),IWORK(NV),KEEP(IPS),IWORK(NELIM),
     +            KEEP(NSTK),KEEP(NODE),KEEP(PERM),
     +            KEEP(NSTEPS),IWORK(FILS),IWORK(FRERE),IWORK(ND),
     +            NEMIN,KEEP(IRNPRM))
      NST = KEEP(NSTEPS)

C Construct map for storing the permuted upper triangle by rows.
C Input N,NE,IRN,JCN,PERM
C Output MAP,LROW,IRNPRM
      CALL MA57MD(N,NE,IRN,JCN,KEEP(MAP),KEEP(IRNPRM),
     +            KEEP(LROW),KEEP(PERM),
     +            IWORK(IW2),IWORK(IW5))

C Set number of entries in expanded input matrix
      KEEP(EXPNE) = IWORK(IW5)

C Evaluate storage and operation counts.
C Input  LROW,NSTK,NELIM,ND
C Output LROW,NSTK (unchanged)
      CALL MA57ND(N,KEEP(LROW),KEEP(NSTK),IWORK(NELIM),
     +            IWORK(ND),NST,IWORK(IW1),IWORK(IW2),
     +            INFO,RINFO)

C Set INFO entry to record ordering used
   93 INFO(36) = ICNTL6
C Add for BIGA
      ALENB    = 1
C Add for Schnabel-Eskow
      IF (ICNTL(7).EQ.4) ALENB = ALENB + N + 5
C Add for scaling
      IF (ICNTL(15).EQ.1) ALENB = ALENB + N

C Allow enough to get started
      INFO(9)  = MAX(INFO(9)+ALENB,ALENB+KEEP(EXPNE)+1)
      INFO(11) = MAX(INFO(11)+ALENB,ALENB+KEEP(EXPNE)+1)
C This is N+5 for starting the factorization, N for first row (max).
      INFO(10) = MAX(INFO(10),KEEP(EXPNE)+N+5)
      INFO(12) = MAX(INFO(12),KEEP(EXPNE)+N+5)

C Needed by MA57B/BD
      IF (ICNTL(15).EQ.1) THEN
        INFO(9) = MAX(INFO(9),ALENB+3*KEEP(EXPNE)+3*N)
        INFO(11) = MAX(INFO(11),ALENB+3*KEEP(EXPNE)+3*N)
C Allow space for integers in computing scaling factors
        INFO(10) = MAX(INFO(10),3*KEEP(EXPNE)+5*N+1)
        INFO(12) = MAX(INFO(12),3*KEEP(EXPNE)+5*N+1)
      ENDIF

C If requested, print parameter values on exit.
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        NZE = KEEP(EXPNE)
        WRITE (MP,99999) INFO(1),NZE,
     *                  (INFO(I),I=3,13),INFO(36),(RINFO(I),I=1,2)
99999 FORMAT (/'Leaving analysis phase (MA57AD) with ...'/
     1    'INFO(1)  Error indicator                      =',I12/
     2    'Number of entries in matrix with diagonal     =',I12/
     2    'INFO(3)  Number of out-of-range indices       =',I12/
     2    'INFO(4)  Number of off-diagonal duplicates    =',I12/
     2    'INFO(5)  Forecast real storage for factors    =',I12/
     3    '----(6)  Forecast integer storage for factors =',I12/
     3    '----(7)  Forecast maximum front size          =',I12/
     4    '----(8)  Number of nodes in assembly tree     =',I12/
     5    '----(9)  Size of FACT without compress        =',I12/
     6    '----(10) Size of IFACT without compress       =',I12/
     5    '----(11) Size of FACT with compress           =',I12/
     5    '----(12) Size of IFACT with compress          =',I12/
     5    '----(13) Number of compresses                 =',I12/
     5    '----(36) Ordering strategy used by code       =',I12/
     9    'RINFO(1) Forecast additions for assembly      =',1P,D12.5/
     9    'RINFO(2) Forecast ops for elimination         =',1P,D12.5)

        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')  'Permutation array:',
     +                  (KEEP(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +        'Number of entries in rows of permuted matrix:',
     +        (KEEP(LROW+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NZE)
        IF (LDIAG.GE.4) K = NZE
        WRITE (MP,'(/A/(5I12))')
     *        'Column indices of permuted matrix:',
     *                           (KEEP(IRNPRM+I-1),I=1,K)
        IF (K.LT.NZE) WRITE (MP,'(16X,A)') '     . . .'
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')
     +    'Tree nodes at which variables eliminated:',
     +    (KEEP(NODE+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(5I12))') 'Map array:',
     *                               (KEEP(I),I=MAP,MAP+K-1)
        IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
      END IF

      RETURN

C Error conditions.
   20 INFO(1) = -1
      INFO(2) = N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'N has value ',INFO(2)
      RETURN

   30 INFO(1) = -2
      INFO(2) = NE
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'NE has value',INFO(2)
       RETURN

   40 INFO(1) = -15
      INFO(2) = LKEEP
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A,I10/A,I10)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'LKEEP has value    ',INFO(2),
     +    'Should be at least ',5*N+NE+MAX(N,NE)+42
       RETURN

   80 INFO(1) = -9
      INFO(2) = I
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A/A,I10,A)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'Invalid permutation supplied in KEEP',
     +    'Component',INFO(2),' is faulty'
      RETURN

   90 INFO(1) = -18
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(/A,I3/A)')
     +    '**** Error return from MA57AD ****  INFO(1) =',INFO(1),
     +    'MeTiS ordering requested but MeTiS not linked'

      END


C--------------------------------------------------------------------
C            HSL 2000 (2000)
C        --
C-             Copyright Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
      SUBROUTINE MA57BD(N, NE, A, FACT, LFACT, IFACT, LIFACT,
     * LKEEP, KEEP, PPOS, ICNTL, CNTL, INFO, RINFO)
C
C Purpose
C =======
C
C
C This subroutine computes the factorization of the matrix input in
C     A using information (in KEEP and IFACT) from MA57AD.
C
      INTEGER N,NE,LFACT,LIFACT,LKEEP
      DOUBLE PRECISION A(NE),FACT(LFACT)
      DOUBLE PRECISION RINFO(20)
C
C     Control parameters: see description in MA57ID
      DOUBLE PRECISION CNTL(5)
      INTEGER ICNTL(20), IFACT(LIFACT)
      INTEGER   INFO(40), KEEP(LKEEP), PPOS(N)
C
C Parameters
C ==========
C
C N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A. It must be unchanged since the
C    last call to MA57AD and is not altered by the
C    subroutine.  Restriction: N > 0
C
C NE is an INTEGER variable which must be set by the user to the
C    number of    entries in the matrix A.  It is not altered by
C    the subroutine.  Restriction: NE >= 0.
C
C A is a REAL (DOUBLE_PRECISION in the D version) array of length NE.
C   It is not altered by the subroutine.
C
C FACT is a REAL (DOUBLE_PRECISION in the D version)
C      array of length LFACT. It need not
C      be set by the user. On exit, entries 1 to INFO(15) of FACT hold
C      the real part of the factors and should be passed unchanged to
C      MA57CD.
C
C LFACT is an INTEGER variable that must be set by the user to
C      the size of array FACT.
C      It should be passed unchanged to MA57CD.
C
C IFACT is an INTEGER array of length LIFACT. It need not
C      be set by the user. On exit, entries 1 to INFO(16) of IFACT hold
C      the integer part of the factors and should be passed unchanged to
C      MA57CD.
C
C LIFACT is an INTEGER variable that must be set by the user to
C      the size of array IFACT.
C      It should be passed unchanged to MA57CD.
C
C LKEEP is an INTEGER that must be set to the length of array KEEP.
C      Restriction: LKEEP >= 5*N+NE+MAX(N,NE)+42.
C
C KEEP is an INTEGER  array of length LKEEP which must be
C       passed unchanged since the last call to MA57AD.  It is not
C       altered by MA57BD.
C
C PPOS is an INTEGER array of length N that is used as workspace.
C
C ICNTL is an INTEGER array of length 20
C            that contains control parameters and must be set
C            by the user.
C            Default values for the components may be set by a call
C            to MA57ID. Details
C            of the control parameters are given in MA57ID
C
C CNTL is a REAL (DOUBLE_PRECISION in the D version) array of length 5
C            that contains control parameters and must be set
C            by the user.
C            Default values for the components may be set by a call
C            to MA57ID. Details
C            of the control parameters are given in MA57ID
C
C RINFO is a REAL (DOUBLE_PRECISION in the D version) array of
C         length 20 which need not be set by the user.
C         We describe below the components of this array modified
C         in the subroutine.
C
C    ______(3)  Number of floating point operations involved
C         during the assembly process
C
C    ______(4)  Number of floating point operations involved
C         during the elimination process
C
C    ______(5)  Number of extra floating point operations caused by
C         use of Level 3 BLAS
C
C
C      .. Error Return ..
C         ============
C
C  A successful return from MA57BD
C  is indicated by a value of INFO(1) positive.
C  Negative values of INFO(1) correspond to
C  error message whereas positive values correspond to
C  warning messages. A negative value of INFO(1) is associated with
C  an error message which will be output on unit ICNTL(1).
C
C
C       .. Local variables ..
C          ===============
C
      INTEGER EXPNE,HOLD,I,IRNPRM,K,LDIAG,LLFACT,LP,LROW,MAP,MM1,MM2,MP
      INTEGER J,JJ,KK,ISCALE,NUM,NE64,IDUP,IMAT,IPT,JLOOP,JNEW,NN,ISING
      INTEGER NSTEPS,NODE,NSTK,PERM,INEW,ALENB,BIGA

      DOUBLE PRECISION ONE,ZERO,RINF,FD15AD,FCT,SMAX,SMIN,REPS
      PARAMETER (ONE = 1.0D0, ZERO=0.0D0)

      INTRINSIC MIN

C
C EXPNE is number of entries of original matrix plus any missing
C    diagonals.
C HOLD points to position in IFACT for first entry of array that holds
C    values for restart.
C MM1,MM2 are used to define start point of arrays for matrix
C    modification.  Needed if LFACT < N + 5.
C NSTEPS holds the number of nodes in the tree.

C
C External Subroutines
C ====================
C

      EXTERNAL MA57OD,MA57UD,FD15AD,MC34AD,MC64WD

C Set RINF to largest positive real number (infinity)
      RINF = FD15AD('H')
C Set REPS to smallest number st 1.0+REPS > 1.0
      REPS = FD15AD('E')

C Set INFO for -3 and -4 error return. Ok that is done on every call
      INFO(17) = 0
      INFO(18) = 0

C
C
C Initialisation of printing controls.
C
      LP     = ICNTL(1)
      MP     = ICNTL(3)
      LDIAG  = ICNTL(5)
C
C??
C     Check if analysis has been effectively performed
C
      IF (N.LE.0)  GO TO 25
      IF (NE.LT.0) GO TO 30
      IF (LKEEP.LT.5*N+NE+MAX(N,NE)+42) GO TO 40
      IF (ICNTL(7).LT.1 .OR. ICNTL(7).GT.4) GO TO 35

C     Partition of array KEEP
      NSTEPS = KEEP(N+1)
      EXPNE  = KEEP(N+2)
      PERM = 1
      HOLD = PERM + N + 2
      LROW = HOLD + 40
      NODE = LROW + N
      NSTK = NODE + N
      MAP  = NSTK + N
      IRNPRM = MAP + MAX(NE,N)

      BIGA = LFACT
      LLFACT = LFACT - 1

      IF (ICNTL(15).EQ.1) THEN
C Matrix is being scaled using MC64SYM
        ISCALE = LLFACT - N + 1
        LLFACT = ISCALE - 1
      ENDIF

      IF (ICNTL(7).EQ.4) THEN
C Schnabel-Eskow modification being used
C Reserve space in FACT for diagonal entries and controls.
        LLFACT = LLFACT - N - 5
C Set MM1 and MM2 to point to first entries in arrays.
        MM1 = LLFACT+6
        MM2 = LLFACT+1
      ELSE
        MM1 = 1
        MM2 = 1
      ENDIF

C One for BIGA
      ALENB = 1
      IF (ICNTL(7).EQ.4)  ALENB = ALENB + N + 5
      IF (ICNTL(15).EQ.1) ALENB = ALENB + N
C +1 because MAP of o-o-r maps to entry 0
      IF (LLFACT.LT.EXPNE+1)   GO TO 85
C The first five positions and room for a whole row are needed
C at the beginning of the MA57O/OD factorization
      IF (LIFACT.LT.EXPNE+N+5)  GO TO 95
C Check that there is enough space for scaling within MA57B/BD.
      IF (ICNTL(15).EQ.1)  THEN
        IF (LFACT .LT. ALENB + 3*EXPNE  + 3*N) GO TO 85
        IF (LIFACT .LT. 3*EXPNE + 5*N + 1) GO TO 95
      ENDIF

C
C PRINTING OF INPUT PARAMETERS
C*****************************
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999)
99999 FORMAT (//'Entering factorization phase (MA57BD) with ...')
        IF (KEEP(HOLD).GT.0) WRITE (MP,99998)
99998 FORMAT ('Re-entry call after call to MA57ED')
        WRITE (MP,99997) N,NE,EXPNE,(ICNTL(I),I=1,5),ICNTL(7),ICNTL(8),
     +         ICNTL(11),ICNTL(15),ICNTL(16), LFACT, LIFACT, NSTEPS,
     +         CNTL(1), CNTL(2), CNTL(4), CNTL(5)
99997 FORMAT ('N       Order of input matrix               =',I12/
     2        'NE      Entries in input matrix             =',I12/
     2        '        Entries in input matrix (inc diags) =',I12/
     6        'ICNTL(1)  Stream for errors                 =',I12/
     7        ' --- (2)  Stream for warnings               =',I12/
     8        ' --- (3)  Stream for monitoring             =',I12/
     9        ' --- (4)  Stream for statistics             =',I12/
     1        ' --- (5)  Level of diagnostic printing      =',I12/
     1        ' --- (7)  Numerical pivoting control        =',I12/
     1        ' --- (8)  Restart or discard factors        =',I12/
     1        ' --- (11) Block size for Level 3 BLAS       =',I12/
     1        ' --- (15) Scaling control (1 on)            =',I12/
     1        ' --- (16) Dropping control (1 on)           =',I12/
     4        'LFACT   Size of real working space          =',I12/
     5        'LIFACT  Size of integer working space       =',I12/
     7        '        Number nodes in assembly tree       =',I12/
     9        'CNTL(1) Value of threshold parameter        =',D12.5/
     9        'CNTL(2) Threshold for zero pivot            =',D12.5/
     9        'CNTL(4) Control for value of static pivots  =',D12.5/
     9        'CNTL(5) Control for number delayed pivots   =',D12.5)
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        IF (NE.GT.0) THEN
          WRITE (MP,'(/A/(3(I6,A,1P,D16.8,A)))') 'Matrix entries:',
     +     (I,': (',A(I),')',I=1,K)
          IF (K.LT.NE) WRITE (MP,'(A)') '     . . .'
        END IF
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE (MP,'(/A/(5I12))')  'Permutation array:',
     +                    (KEEP(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +          'Number of entries in rows of permuted matrix:',
     +          (KEEP(LROW+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        WRITE (MP,'(/A/(5I12))')
     +    'Tree nodes at which variables eliminated:',
     +    (KEEP(NODE+I-1),I=1,K)
        IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NSTEPS)
        IF (LDIAG.GE.4) K = NSTEPS
        IF (K.GT.0) WRITE (MP,'(/A/(5I12))')
     +     'Number of assemblies at each tree node:',
     +     (KEEP(NSTK+I-1),I=1,K)
        IF (K.LT.NSTEPS) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,NE)
        IF (LDIAG.GE.4) K = NE
        WRITE (MP,'(/A/(5I12))') 'Map array:',
     *                               (KEEP(I),I=MAP,MAP+K-1)
        IF (K.LT.NE) WRITE (MP,'(16X,A)') ' . . .'
        K = MIN(10,EXPNE)
        IF (LDIAG.GE.4) K = EXPNE
        WRITE (MP,'(/A/(5I12))')
     *          'Column indices of permuted matrix:',
     *                             (KEEP(IRNPRM+I-1),I=1,K)
        IF (K.LT.EXPNE) WRITE (MP,'(16X,A)') '     . . .'
      ENDIF


C Can't do because all of A isn't there now .. also it is now scaled
C     BIGA = ZERO
C     DO 291 K = 1,NE
C       BIGA = MAX(BIGA,ABS(A(K)))
C 291 CONTINUE

C Jump if it is reentry
      IF (KEEP(HOLD) .GT. 0) GO TO 22

C
C***************************************************
C MAP input nonzeros to appropriate position in FACT
C***************************************************
C
C?? For the moment to handle missing diagonals
      DO 19 K = 1,EXPNE
        FACT(LLFACT-EXPNE+K) = ZERO
   19 CONTINUE

      FACT(BIGA) = ZERO
      DO 20 K = 1,NE
        FACT(BIGA) = MAX(FACT(BIGA),ABS(A(K)))
        FACT(KEEP(MAP+K-1)+LLFACT-EXPNE) = A(K)
   20 CONTINUE
      RINFO(18) = FACT(BIGA)
      DO 21 K = 1,EXPNE
        IFACT(LIFACT-EXPNE+K) = KEEP(IRNPRM+K-1)
   21 CONTINUE
C Invert array PERM
      DO 23 I = 1,N
        PPOS(KEEP(PERM+I-1)) = I
   23 CONTINUE

      IF (ICNTL(15).EQ.1) THEN
C Scaling using MC64.  Matrix must be generated in correct format.

        IPT = 1
        IDUP = IPT+N+1
        IMAT = IDUP+N
        ISING = IMAT + MAX(NE,EXPNE)

C Copy matrix, remove duplicates, and initialize IP array.
        DO 4444 I = 1,N
          IFACT(IDUP+I-1) = 0
 4444   CONTINUE
C       IFACT(IPT)=1
C       DO 9999 I = 1, N
C         IFACT(IPT+I) = IFACT(IPT+I-1)+KEEP(LROW+I-1)
C9999   CONTINUE
C Must use new coordinates to keep matrix (half) symmetric
        IFACT(IPT) = 1
        KK = 1
        K = 1
        DO 3333 J = 1,N
          DO 2222 JJ = 1,KEEP(LROW+J-1)
            I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
            IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
C Duplicate
              FACT(IFACT(IDUP+I-1)) =
     &          FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
            ELSE
C Remove explicit zeros
              IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
                IFACT(IDUP+I-1) = KK
                FACT(KK) = FACT(LLFACT-EXPNE+K)
                IFACT(IMAT-1+KK) = I
                KK = KK+1
              ENDIF
            ENDIF
            K = K + 1
 2222     CONTINUE
          IFACT(IPT+J) = KK
 3333   CONTINUE
C Expand matrix
        CALL MC34AD(N,IFACT(IMAT),IFACT(IPT),.TRUE.,FACT,KEEP(PERM))
        NE64 = IFACT(IPT+N)-1
        DO 75 J = 1,N
          FCT = ZERO
          DO 60 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            FACT(K) = ABS(FACT(K))
            IF (FACT(K).GT.FCT) FCT = FACT(K)
   60     CONTINUE
          FACT(NE64+2*N+J) = FCT
          IF (FCT.NE.ZERO) THEN
            FCT = LOG(FCT)
          ELSE
C This can happen if only if column is null so matrix is singular.
            FCT = RINF/N
          ENDIF
          DO 70 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
CCC
C Note that zeros have been screened out so that FACT(K) always > 0.
C           IF (FACT(K).NE.ZERO) THEN
              FACT(K) = FCT - LOG(FACT(K))
C           ELSE
C             FACT(K) = RINF/N
C           ENDIF
   70     CONTINUE
   75   CONTINUE
C Scale matrix
C B = DW(3N+1:3N+NE); IW(1:5N) and DW(1:2N) are workspaces
        CALL MC64WD(N,NE64,IFACT(IPT),IFACT(IMAT),FACT,KEEP(PERM),NUM,
     &      IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
     &      IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
     &      FACT(NE64+1),FACT(NE64+N+1))
        IF (NUM.EQ.N) THEN
          DO 80 J = 1,N
C           IF (FACT(NE64+2*N+J).NE.ZERO) THEN
              FACT(NE64+N+J) = FACT(NE64+N+J) - LOG(FACT(NE64+2*N+J))
CCC
C           ELSE
C Can't happen because only possible if matrix was singular NUM NE N).
C             FACT(NE64+N+J) = ZERO
C           ENDIF
   80     CONTINUE
C Check size of scaling factors
C     FCT = 0.5*LOG(RINF)
C     DO 86 J = 1,N
C       IF (FACT(NE64+J).LT.FCT .AND. FACT(NE64+N+J).LT.FCT) GO TO 86
C       INF64(1) = 2
C  86 CONTINUE

C Scaling is permuted to scaling on original matrix
          DO 5555 I=1,N
            FACT(ISCALE+PPOS(I)-1) =
     &        SQRT(EXP(FACT(NE64+I)+FACT(NE64+N+I)))
 5555     CONTINUE
        ELSE
C Matrix is singular
C Regenerate PERM and set PPOS to indicate nonsingular block
        K = 0
        DO 3501 I = 1,N
          IF (KEEP(PERM+I-1).LT.0) THEN
            PPOS(I) = -PPOS(I)
            IFACT(ISING+I-1) = 0
          ELSE
            K = K + 1
            IFACT(ISING+I-1) = K
          ENDIF
 3501   CONTINUE
        DO 3502 I = 1,N
          KEEP(PERM+ABS(PPOS(I))-1) = I
 3502   CONTINUE
C Copy matrix, remove duplicates, and initialize IP array.
        DO 3503 I = 1,N
          IFACT(IDUP+I-1) = 0
 3503   CONTINUE
C Must use new coordinates to keep matrix (half) symmetric
        IFACT(IPT) = 1
        KK = 1
        K = 1
        JNEW = 0
        NN = N
        DO 3505 J = 1,N
          IF (PPOS(J).LT.0) THEN
            NN = NN - 1
            K = K + KEEP(LROW+J-1)
            GO TO 3505
          ENDIF
          JNEW = JNEW + 1
          DO 3504 JJ = 1,KEEP(LROW+J-1)
            I = KEEP(PERM+IFACT(LIFACT-EXPNE+K)-1)
            IF (PPOS(I).GT.0) THEN
              IF (IFACT(IDUP+I-1).GE.IFACT(IPT+J-1)) THEN
C Duplicate
                FACT(IFACT(IDUP+I-1)) =
     &            FACT(IFACT(IDUP+I-1)) + FACT(LLFACT-EXPNE+K)
              ELSE
C Remove explicit zeros
                IF (FACT(LLFACT-EXPNE+K).NE.ZERO) THEN
                  IFACT(IDUP+I-1) = KK
                  FACT(KK) = FACT(LLFACT-EXPNE+K)
                  IFACT(IMAT-1+KK) = IFACT(ISING+I-1)
                  KK = KK+1
                ENDIF
              ENDIF
            ENDIF
            K = K + 1
 3504     CONTINUE
          IFACT(IPT+JNEW) = KK
 3505   CONTINUE
      NE64 = IFACT(IPT+NN)-1
        CALL MC34AD(NN,IFACT(IMAT),IFACT(IPT),.TRUE.,FACT,KEEP(PERM))
        NE64 = IFACT(IPT+NN)-1
        DO 3508 J = 1,NN
          FCT = ZERO
          DO 3506 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
            FACT(K) = ABS(FACT(K))
            IF (FACT(K).GT.FCT) FCT = FACT(K)
 3506     CONTINUE
          FACT(NE64+2*N+J) = FCT
CCC
C         IF (FCT.NE.ZERO) THEN
            FCT = LOG(FCT)
C         ELSE
C Can't happen because no null columns and zeros screened out.
C           FCT = RINF/NN
C         ENDIF
          DO 3507 K = IFACT(IPT+J-1),IFACT(IPT+J)-1
C           IF (FACT(K).NE.ZERO) THEN
              FACT(K) = FCT - LOG(FACT(K))
C           ELSE
C Can't happen because zeros screened out.
C             FACT(K) = RINF/NN
C           ENDIF
 3507     CONTINUE
 3508   CONTINUE
        CALL MC64WD(NN,NE64,IFACT(IPT),IFACT(IMAT),FACT,KEEP(PERM),NUM,
     &      IFACT(IDUP),IFACT(IMAT+NE64),IFACT(IMAT+NE64+N),
     &      IFACT(IMAT+NE64+2*N),IFACT(IMAT+NE64+3*N),
     &      FACT(NE64+1),FACT(NE64+N+1))
        DO 3509 J = 1,NN
CCC
C           IF (FACT(NE64+2*N+J).NE.ZERO) THEN
              FACT(NE64+N+J) = FACT(NE64+N+J) - LOG(FACT(NE64+2*N+J))
C           ELSE
C As in original matrix ..this can't happen.
              FACT(NE64+N+J) = ZERO
C           ENDIF
 3509     CONTINUE
C Check size of scaling factors
C     FCT = 0.5*LOG(RINF)
C     DO 86 J = 1,N
C       IF (FACT(NE64+J).LT.FCT .AND. FACT(NE64+N+J).LT.FCT) GO TO 86
C       INF64(1) = 2
C  86 CONTINUE

C Scaling is permuted to scaling on original matrix for scale factors
C for nonsingular block
          K=0
C Loop is on new indices
          DO 3510 I=1,N
            IF (PPOS(I).LT.0) THEN
              K = K + 1
              FACT(ISCALE-PPOS(I)-1) = ZERO
            ELSE
              FACT(ISCALE+PPOS(I)-1) =
     &          SQRT(EXP(FACT(NE64+I-K)+FACT(NE64+N+I-K)))
            ENDIF
 3510     CONTINUE
C Compute scaling on nonsingular part
C Remember that PPOS maps from new to original but is flag on new
          DO 3516 I = 1,N
            KEEP(PERM+ABS(PPOS(I))-1) = I
 3516     CONTINUE
C Looping on new indices
          K = 1
          DO 3514 JJ = 1,N
            J = PPOS(JJ)
            IF (J.GT.0) THEN
              DO 3511 JLOOP = 1,KEEP(LROW+JJ-1)
                I = IFACT(LIFACT-EXPNE+K)
C I is original index so have to map to new to do PPOS test
                INEW = KEEP(PERM+I-1)
                IF (PPOS(INEW).LT.0)
     &            FACT(ISCALE+I-1) = MAX(FACT(ISCALE+I-1),
     &                 ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+J-1))
                K = K + 1
 3511         CONTINUE
            ELSE
              DO 3512 JLOOP = 1,KEEP(LROW+JJ-1)
                I = IFACT(LIFACT-EXPNE+K)
C I is original index so have to map to new to do PPOS test
                INEW = KEEP(PERM+I-1)
C Shouldn't happen otherwise nonsingular block not maximum
C Sorry can happen but entry is implicit zero on diagonal
C Note that J is negative
                IF (I .NE. -J)  THEN
                FACT(ISCALE-J-1) =
     &              MAX(FACT(ISCALE-J-1),
     &              ABS(FACT(LLFACT-EXPNE+K))*FACT(ISCALE+I-1))
                ENDIF
                K = K + 1
 3512         CONTINUE
            ENDIF
 3514     CONTINUE
C Set scaling factors for singular block and reset PPOS
          DO 3513 I = 1,N
            INEW = KEEP(PERM+I-1)
            IF (PPOS(INEW) .LT. 0) THEN
              PPOS(INEW) = - PPOS(INEW)
              IF (FACT(ISCALE+I-1) .EQ. ZERO) THEN
                FACT(ISCALE+I-1) = ONE
              ELSE
                FACT(ISCALE+I-1) = ONE/FACT(ISCALE+I-1)
              ENDIF
            ENDIF
 3513     CONTINUE
        ENDIF
C End of logic for singular matrix

C         DO 8888 I = 1, N
C           FACT(ISCALE+I-1) = ONE
C8888     CONTINUE
          SMAX = FACT(ISCALE)
          SMIN = FACT(ISCALE)
          DO 5566 I = 1,N
            SMAX = MAX(SMAX,FACT(ISCALE+I-1))
            SMIN = MIN(SMIN,FACT(ISCALE+I-1))
 5566     CONTINUE
          RINFO(16) = SMIN
          RINFO(17) = SMAX
C Scale matrix
          K = 1
          FACT(BIGA) = ZERO
          DO 6666 JJ = 1,N
            J = PPOS(JJ)
            DO 7777 JLOOP = 1,KEEP(LROW+JJ-1)
              I = IFACT(LIFACT-EXPNE+K)
              FACT(LLFACT-EXPNE+K) =
     &          FACT(ISCALE+I-1)*FACT(LLFACT-EXPNE+K)*FACT(ISCALE+J-1)
              FACT(BIGA) = MAX(FACT(BIGA), ABS(FACT(LLFACT-EXPNE+K)))
              K = K + 1
 7777       CONTINUE
 6666     CONTINUE
C End of scaling
      ELSE
        RINFO(16) = ONE
        RINFO(17) = ONE
      ENDIF
C
C**********************************
C Numerical Factorization
C**********************************
C Work arrays FACT(MM1/MM2), KEEP(PERM), IFACT(1)
   22 CALL MA57OD(N, EXPNE, FACT, LLFACT, IFACT, LIFACT, KEEP(LROW),
     *            PPOS,
     *            NSTEPS, KEEP(NSTK), KEEP(NODE), FACT(MM1),
     *            FACT(MM2),
     *            KEEP(PERM),
     *            CNTL, ICNTL,
     *            INFO, RINFO, KEEP(HOLD), FACT(BIGA))
      IF (INFO(1).EQ.10 .OR. INFO(1).EQ.11) THEN
        IF (LDIAG.GT.2 .AND. MP.GE.0)  THEN
          IF (INFO(1).EQ.10) WRITE (MP,99982) INFO(1)
99982 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'Factorization suspended because of lack of real space'/
     1  'INFO (1) = ',I3)
          IF (INFO(1).EQ.11) WRITE (MP,99983) INFO(1)
99983 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'Factorization suspended because of lack of integer space'/
     1  'INFO (1) = ',I3)
        ENDIF
        RETURN
      ENDIF
C Regenerate array PERM
      DO 24 I = 1,N
        KEEP(PERM+PPOS(I)-1) = I
   24 CONTINUE
C Compute INFO(17-20)
        INFO(17) = ALENB + INFO(17)
        INFO(19) = ALENB + INFO(19)
C Allow space for scaling in MA57B/BD
      IF (ICNTL(15).EQ.1) THEN
        INFO(17) = MAX(INFO(17),ALENB + 3*EXPNE+3*N)
        INFO(19) = MAX(INFO(19),ALENB + 3*EXPNE+3*N)
        INFO(18) = MAX(INFO(18),3*EXPNE+5*N+1)
        INFO(20) = MAX(INFO(20),3*EXPNE+5*N+1)
      ENDIF
      IF (INFO(1).EQ.-3) GO TO 85
      IF (INFO(1).EQ.-4) GO TO 95
      IF (INFO(1).LT.0) RETURN
      GO TO 100
C************************
C **** Error returns ****
C************************
   25 INFO(1) = -1
      INFO(2) =  N
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'N has value ',INFO(2)
      RETURN
   30 INFO(1) = -2
      INFO(2) = NE
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'NE has value',INFO(2)
      RETURN
   40 INFO(1) = -15
      INFO(2) = LKEEP
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'LKEEP has value    ',INFO(2),
     +    'Should be at least ',5*N+NE+MAX(N,NE)+42
      RETURN
   35 INFO(1) = -10
      INFO(2) = ICNTL(7)
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'ICNTL(7) has value',ICNTL(7)
      RETURN

   85 INFO(1) = -3
      INFO(2) = LFACT
      INFO(17) = MAX(INFO(17), ALENB + EXPNE + 1)
      IF (ICNTL(15).EQ.1) 
     *    INFO(17) = MAX(INFO(17), ALENB + 3*EXPNE + 3*N)
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'Insufficient real space in FACT, LFACT = ',INFO(2)
      RETURN

   95 INFO(1) = -4
      INFO(2) = LIFACT
      INFO(18) = MAX(INFO(18), EXPNE+N+5)
      IF (ICNTL(15).EQ.1) 
     *    INFO(18) = MAX(INFO(18), 3*EXPNE + 5*N + 1)
      IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57BD ****  INFO(1) =',INFO(1),
     +    'Insufficient integer space in IFACT, LIFACT = ',INFO(2)
      RETURN

C****************
C Printing section
C****************
 100  IF (LDIAG.LE.2 .OR. MP.LT.0) RETURN
      WRITE (MP,99980) INFO(1), INFO(2),
     *    (INFO(I),I=14,25),INFO(28),INFO(29)
      WRITE (MP,99984) (INFO(I),I=31,35),RINFO(3), RINFO(4),
     *                 RINFO(5), RINFO(18)
99980 FORMAT (/'Leaving factorization phase (MA57BD) with ...'/
     1  'INFO (1)                                      =',I12/
     2  ' --- (2)                                      =',I12/
     3  ' --- (14) Number of entries in factors        =',I12/
     4  ' --- (15) Real storage for factors            =',I12/
     5  ' --- (16) Integer storage for factors         =',I12/
     6  ' --- (17) Min LFACT with compresses           =',I12/
     7  ' --- (18) Min LIFACT with compresses          =',I12/
     8  ' --- (19) Min LFACT without compresses        =',I12/
     9  ' --- (20) Min LIFACT without compresses       =',I12/
     *  ' --- (21) Order of largest frontal matrix     =',I12/
     1  ' --- (22) Number of 2x2 pivots                =',I12/
     2  ' --- (23) Number of delayed pivots            =',I12/
     3  ' --- (24) Number of negative eigenvalues      =',I12/
     4  ' --- (25) Rank of factorization               =',I12/
     5  ' --- (28) Number compresses on real data      =',I12/
     6  ' --- (29) Number compresses on integer data   =',I12)
      IF (ICNTL(15).EQ.1) WRITE (MP,99985) RINFO(16),RINFO(17)
99985 FORMAT (
     1  'RINFO(16) Minimum value of scaling factor     =  ',1PD10.3/
     2  '-----(17) Maximum value of scaling factor     =  ',1PD10.3)
99984 FORMAT (
     7  ' --- (31) Number of block pivots in factors   =',I12/
     7  ' --- (32) Number of zeros factors triangle    =',I12/
     7  ' --- (33) Number of zeros factors rectangle   =',I12/
     7  ' --- (34) Number of zero cols factors rect    =',I12/
     7  ' --- (35) Number of static pivots             =',I12/
     1  'RINFO(3)  Operations during node assembly     =  ',1PD10.3/
     2  '-----(4)  Operations during node elimination  =  ',1PD10.3/
     3  '-----(5)  Extra operations because of BLAS    =  ',1PD10.3/
     3  '-----(18) Largest modulus of entry in matrix  =  ',1PD10.3)
C     IF (ICNTL(16).EQ.1) WRITE (MP,99986) INFO(37)
C9986 FORMAT (
C    1  'INFO(37)  Number entries dropped (ICNTL(16)=1) = ',I12)
      IF (INFO(27).GT.0) WRITE (MP,99981) INFO(27),RINFO(14),RINFO(15)
99981 FORMAT (/'Matrix modification performed'/
     1  'INFO (27) Step at which matrix first modified =',I12/
     2  'RINFO(14) Maximum value added to diagonal     =  ',1PD10.3/
     2  'RINFO(15) Smallest pivot in modified matrix   =  ',1PD10.3)
C Print out matrix factors from MA57BD.
      CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
C Print scaling factors
      IF (ICNTL(15).NE.1) RETURN
      K = MIN(10,N)
      IF (LDIAG.GE.4) K = N
      WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
     +                    (FACT(ISCALE+I-1),I=1,K)
      IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'

      END
      SUBROUTINE MA57CD(JOB,N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,W,
     *                  LW,IW1,ICNTL,INFO)
C This subroutine uses the factorisation of the matrix in FACT,IFACT to
C     solve a system of equations.
      INTEGER JOB,N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20),INFO(40)
C JOB must be set by the user to determine the coefficient matrix
C     of the equations being solved.  If the factorization is
C                            T  T
C             A =  P  L  D  L  P
C     then coefficient matrix is:
C
C     JOB <= 1   A
C                       T
C     JOB  = 2   P  L  P
C                       T
C     JOB  = 3   P  D  P
C                    T  T
C     JOB >= 4   P  L  P
C
C N must be set to the order of the matrix and must be unchanged since
C     the call to MA57BD. It is not altered.
C FACT holds information on the factors and must be unchanged since
C     the call to MA57BD. It is not altered by MA57CD.
C LFACT must be set to the length of FACT. It is not altered.
C IFACT holds information on the factors and must be unchanged since
C     the call to MA57BD. It is not altered by MA57CD.
C LIFACT must be set to the length of IFACT. It is not altered.
C NRHS is the number of right-hand sides being solved for.
C RHS must be set to the right hand sides for the equations being
C     solved. On exit, this array will hold the solutions.
C LHS must be set to the leading dimension of array RHS.
C W is used as a work array.
C LW  must be set to the length of array W.  It must be at least
C     as large as N*NRHS.  (Actually only INFO(21)*NRHS but no way to
C     check this).
C IW1 is used as a work array.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 0 suppresses output.
C     ICNTL(2) must be set to the stream number for warning output.
C       A value less than 0 suppresses output.
C     ICNTL(3) must be set to the stream number for monitor output.
C       A value less than 0 suppresses output.
C     ICNTL(4) must be set to the stream number for statistics output.
C       A value less than 0 suppresses output.
C     ICNTL(5) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C     > 3  As 2, plus all parameters on entry and exit.
C     ICNTL(6:12) Not referenced.
C     ICNTL(13) Threshold on number of columns in a block for direct
C           addressing using Level 2 and Level 3 BLAS.

C Procedures
      INTRINSIC MIN
      EXTERNAL MA57QD,MA57RD,MA57SD,MA57TD,MA57UD,MA57XD,MA57YD

C
C Local variables
      DOUBLE PRECISION SCALE,ONE
      PARAMETER (ONE = 1.0D0)
      INTEGER I,J,K,LDIAG,LLW,LP,MP,ISCALE
C I  Temporary variable.
C J  Temporary variable.
C K  Temporary variable.
C LDIAG Control for amount of information output.
C LP Stream number for error printing.

C      C MP Stream number for monitor printing.
C

C Set local print control variables
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)

      INFO(1) = 0

C Check input data
      IF (N.LE.0) THEN
        INFO(1) = -1
        INFO(2) = N
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'N has value',N
        GOTO 500
      ENDIF

      IF (NRHS.LT.1) THEN
        INFO(1) = -16
        INFO(2) = NRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of NRHS =',NRHS,' is less than 1'
        GOTO 500
      ENDIF

      IF (LRHS.LT.N) THEN
        INFO(1) = -11
        INFO(2) = LRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of LRHS =',LRHS,' is less than N=',N
        GOTO 500
      ENDIF

      IF (LW.LT.N*NRHS) THEN
        INFO(1) = -17
        INFO(2) = N*NRHS
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I10,A,I10)')
     +    '**** Error return from MA57CD ****  INFO(1) =',INFO(1),
     +    'value of LW =',LW,' is less than', N*NRHS
        GOTO 500
      ENDIF

C If requested, print input parameters
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999) JOB,N,(ICNTL(I),I=1,5),LFACT,LIFACT,NRHS,
     +         LRHS,LW,ICNTL(13)
99999 FORMAT(/'Entering solution phase (MA57CD) with ...'/
     +    'JOB       Control on coefficient matrix       =',I12/
     +    'N         Order of matrix                     =',I12/
     6    'ICNTL(1)  Stream for errors                   =',I12/
     7    ' --- (2)  Stream for warnings                 =',I12/
     8    ' --- (3)  Stream for monitoring               =',I12/
     9    ' --- (4)  Stream for statistics               =',I12/
     1    ' --- (5)  Level of diagnostic printing        =',I12/
     +    'LFACT     Length of array FACT                =',I12/
     +    'LIFACT    Length of array IFACT               =',I12/
     +    'NRHS      Number of right-hand sides          =',I12/
     +    'LRHS      Leading dimension of RHS array      =',I12/
     +    'LW        Leading dimension of work array     =',I12/
     +    'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12)

C Print out matrix factors.
        CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
C Print scaling factors
        IF (ICNTL(15).EQ.1) THEN
          ISCALE = LFACT-N
          K = MIN(10,N)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(/A/(5D12.5))')  'Scaling factors:',
     +                        (FACT(ISCALE+I-1),I=1,K)
          IF (K.LT.N) WRITE (MP,'(16X,A)') ' . . .'
        ENDIF
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        DO 10 J = 1,NRHS
          WRITE(MP,'(/A,I10)') 'Right-hand side',J
          WRITE (MP,'((1P,5D13.3))') (RHS(I,J),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
   10   CONTINUE
      END IF

      LLW = LW/NRHS


C Scale right-hand side
      IF (ICNTL(15).EQ.1) THEN
        ISCALE = LFACT-N
        DO 5555 I = 1, N
          SCALE = FACT(ISCALE+I-1)
          IF (JOB.GE.4) SCALE = ONE/FACT(ISCALE+I-1)
          DO 4444 J = 1, NRHS
            RHS(I,J) = SCALE*RHS(I,J)
 4444     CONTINUE
 5555   CONTINUE
      ENDIF

C Forward substitution
      IF (JOB.LE.2) THEN
        IF (NRHS.EQ.1) THEN
          CALL MA57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ELSE
          CALL MA57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ENDIF
        IF (JOB.EQ.2) GO TO 15
C Back substitution.
        IF (NRHS.EQ.1) THEN
          CALL MA57YD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ELSE
          CALL MA57RD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                W,LLW,IW1,ICNTL)
        ENDIF
      ENDIF
      IF (JOB.EQ.3)
     *  CALL MA57SD(FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *              W,LLW,ICNTL)
      IF (JOB.GE.4)
     *  CALL MA57TD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *              W,LLW,IW1,ICNTL)

C Scale solution
   15 IF (ICNTL(15).EQ.1) THEN
        ISCALE = LFACT-N
        DO 6666 I = 1, N
          SCALE = FACT(ISCALE+I-1)
          IF (JOB.EQ.2) SCALE = ONE/FACT(ISCALE+I-1)
          DO 7777 J = 1, NRHS
            RHS(I,J) = SCALE*RHS(I,J)
 7777     CONTINUE
 6666   CONTINUE
      ENDIF

C
C If requested, print output parameters.
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,'(//A)')
     *       'Leaving solution phase (MA57CD) with ...'
        DO 20 J = 1,NRHS
          WRITE(MP,'(/A,I10)') 'Solution       ',J
          WRITE (MP,'(1P,5D13.3)') (RHS(I,J),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
   20   CONTINUE
      ENDIF

  500 RETURN

      END


      SUBROUTINE MA57QD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
C This subroutine performs forward elimination using the factors
C     stored in FACT/IFACT by MA57BD.
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C     factors. This must be unchanged since the preceding call to
C     MA57BD. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C     in FACT. This must be unchanged since the preceding call to
C     MA57BD. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS on input, must be set to hold the right hand side vector.  On
C     return, it will hold the modified vector following forward
C     elimination.
C LHS must be set to the leading dimension of array RHS.
C W   used as workspace to hold the components of the right hand
C     sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C     larger than INFO(21) as returned from MA57BD.
C IW1 need not be set on entry. On exit IW1(I) (I = 1,NBLK), where
C     NBLK = IFACT(3) is the number of block pivots, will
C     hold pointers to the beginning of each block pivot in array IFACT.
C ICNTL Not referenced except:
C     ICNTL(13) Threshold on number of columns in a block for using
C           addressing using Level 2 and Level 3 BLAS.
C

C Procedures
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV

C Constant
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C Local variables
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,
     +        NCOLS,NROWS
      DOUBLE PRECISION W1
C
C APOS  Current position in array FACT.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C K     Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.

      APOS = 1
      IWPOS = 4
      DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
        IW1(IBLK) = IWPOS
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Load appropriate components of right-hand sides into W.
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 11 J = 1,NRHS
              W(I,J) = RHS(II,J)
   11       CONTINUE
   10     CONTINUE


C Treat diagonal block (direct addressing)
          DO 12 J = 1,NRHS
            CALL DTPSV('L','N','U',NROWS,FACT(APOS),W(1,J),1)
   12     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (direct addressing)
          IF (NCOLS.GT.NROWS) CALL DGEMM('N','N',NCOLS-NROWS,NRHS,NROWS,
     +                                  ONE,FACT(APOS),NCOLS-NROWS,
     +                                  W,LW,ONE,W(NROWS+1,1),LW)
          APOS = APOS + NROWS* (NCOLS-NROWS)

C Reload W back into RHS.
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 36 J = 1,NRHS
              RHS(II,J) = W(I,J)
   36       CONTINUE
   35     CONTINUE

        ELSE

C Perform operations using indirect addressing.

        J1 = IWPOS
        J2 = IWPOS + NROWS - 1


C Treat diagonal block (indirect addressing)
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          DO 101 II = 1,NRHS
            W1 = RHS(ABS(IFACT(J1)),II)
            K = APOS
            DO 100 J = J1+1,J2
              IRHS = ABS(IFACT(J))
              RHS(IRHS,II) = RHS(IRHS,II) - FACT(K)*W1
              K = K + 1
  100       CONTINUE
  101     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE

C Treat off-diagonal block (indirect addressing)
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS
          DO 135 II = 1,NRHS
            K = APOS
            W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)),II)
            DO 133 J = J1,J2
              IRHS = ABS(IFACT(J))
              RHS(IRHS,II) = RHS(IRHS,II) + W1*FACT(K)
              K = K + 1
  133       CONTINUE
  135     CONTINUE
          APOS = K
  136   CONTINUE

      END IF

      IWPOS = IWPOS + NCOLS
  270 CONTINUE

      END


      SUBROUTINE MA57RD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C     using the factors stored in FACT/IFACT by MA57BD.
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C      factors. This must be unchanged since the
C      preceding call to MA57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C      in FACT. This must be unchanged since the preceding call to
C      MA57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C      the forward substitution operations. On exit, holds the
C      solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C      sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C     larger than INFO(21) as returned from MA57BD.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C      number of block pivots, must hold pointers to the beginning of
C      each block pivot in array IFACT, as set by MA57Q/QD.
C      It is not altered.
C ICNTL Not referenced except:
C     ICNTL(13) Threshold on number of columns in a block for using
C           addressing using Level 2 and Level 3 BLAS.

C Procedures
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV

C Constants
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C Local variables.
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,J1,J2,K,KK,LROW,NCOLS,NROWS
      DOUBLE PRECISION W1
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C       pivot.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C       the second.
C K     Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K     Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
      APOS = IFACT(1)
      APOS2 = IFACT(2)
C Run through block pivot rows in the reverse order.
      DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Load latter part of right-hand side into W.
          DO 5 I = NROWS + 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 3 J = 1,NRHS
              W(I,J) = RHS(II,J)
    3       CONTINUE
    5     CONTINUE


C Multiply by the diagonal matrix (direct addressing)
          DO 10 IPIV = NROWS,1,-1
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            APOS = APOS - (NROWS+1-IPIV)
            DO 9 J = 1,NRHS
              W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9       CONTINUE
   10     CONTINUE
          JPIV = -1
          DO 20 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              DO 19 J = 1,NRHS
                W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19         CONTINUE
              IF (JPIV.EQ.1) APOS2 = APOS2 - 1
              JPIV = -JPIV
            END IF

   20     CONTINUE

C Treat off-diagonal block (direct addressing)
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMM('T','N',NROWS,NRHS,K,ONE,
     +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
     +                           W(NROWS+1,1),LW,ONE,W,LW)

C Treat diagonal block (direct addressing)
          DO 22 J = 1,NRHS
            CALL DTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE
C Reload W back into RHS.
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE

        ELSE
C
C Perform operations using indirect addressing.
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1


C Multiply by the diagonal matrix (indirect addressing)
          JPIV = -1
          DO 210 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            LROW = NROWS + 1 - IPIV

            IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
              APOS = APOS - LROW
              DO 65 J = 1,NRHS
                RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65         CONTINUE
            ELSE
C 2 by 2 pivot
              IF (JPIV.EQ.-1) THEN
                IRHS1 = -IFACT(IWPOS+IPIV-2)
                IRHS2 = -IRHS
                APOS = APOS - LROW - LROW - 1
                DO 68 J = 1,NRHS
                  W1 = RHS(IRHS1,J)*FACT(APOS) +
     +                 RHS(IRHS2,J)*FACT(APOS2)
                  RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
     +                           RHS(IRHS2,J)*FACT(APOS+LROW+1)
                  RHS(IRHS1,J) = W1
   68           CONTINUE
                APOS2 = APOS2 - 1
              END IF
              JPIV = -JPIV
            END IF

  210     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (indirect addressing)
          KK = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 218 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = KK
              DO 215 J = J1,J2
                W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  215         CONTINUE
              RHS(IRHS,II) = W1
  218       CONTINUE
            KK = K
  220     CONTINUE

C Treat diagonal block (indirect addressing)
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            DO 240 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = APOS + 1
              DO 230 J = J1,J2
                W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  230         CONTINUE
              RHS(IRHS,II) = W1
  240       CONTINUE
            J1 = J1 - 1
  260     CONTINUE

        END IF

  380 CONTINUE

      END

      SUBROUTINE MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
C     Print out matrix factors from MA57BD or a symbolic representation
C     of them.
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),ICNTL(20)
C FACT   array holding the reals of the factorization.
C        It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  array holding the integers of the factorization. It is not
C     altered.
C LIFACT length of array IFACT. It is not altered.
C ICNTL is not referenced except:
C     ICNTL(3) must be set to the stream number for diagnostic output.
C       A value less than 1 suppresses output.
C     ICNTL(5) must be set to control the amount of output:
C      <3 None.
C       3 First block only.
C       4 All blocks.
C       5 All blocks, but each entry represented by a single character:
C            + for a positive integer
C            - for a negative integer
C            * for a nonzero entry
C            . for a zero entry

C Procedures
      INTRINSIC MIN,SIGN

C Local variables
      CHARACTER*72 LINE
      INTEGER APOS,APOS2,IBLK,ILINE,IROW,IWPOS,J,JPIV,J1,J2,K,
     +        LDIAG,LEN,MP,NBLK,NCOLS,NROWS
C APOS Current position in FACT.
C APOS2 Position in FACT of next off-diagonal entry of 2x2 pivot.
C ILINE Current position in the line.
C IBLK  Current block.
C IROW  Current row.
C IWPOS Current position in IFACT.
C JPIV  has value 1 only for the first row of a 2x2 pivot.
C J     Column index.
C K     Temporary pointer to position in real array.
C J1    Position of last zero in leading part of row.
C J2    Position of last nonzero in leading part of row.
C K     Temporary DO index.
C LDIAG Control for diagnostic printing.
C LEN   1 if pattern only to be printed and 12 if values to be printed.
C LINE  Character variable in which an output line is built.
C MP    Stream number for warning or diagnostic messages
C NBLK  Number of blocks to be printed.
C NCOLS Number of columns in the block.
C NROWS Number of rows in the block.

      CHARACTER*1 PM(-2:2)
      DATA PM/'*','-','.','+','.'/
      DOUBLE PRECISION ZERO,TINY,FD15AD
      PARAMETER (ZERO=0.0D0)

      EXTERNAL FD15AD

C Initialize MP and LDIAG
      MP = ICNTL(3)
      LDIAG = ICNTL(5)

      TINY = FD15AD('T')

      APOS2 = IFACT(1)
      NBLK = IFACT(3)
      IF (LDIAG.EQ.3) NBLK = MIN(1,NBLK)
      LEN = 12
      IF (LDIAG.EQ.5) LEN = 1
      IF (LEN.EQ.12) THEN
        IF (NBLK.EQ.IFACT(3)) THEN
          WRITE (MP,'(/A)')
     +      'For each block, the following information is provided:'

        ELSE
          WRITE (MP,'(/A,A)') 'For the first block only,',
     +      ' the following information is provided:'
        END IF

      END IF

      IF (LEN.EQ.12) WRITE (MP,'(A)')
     +    '   1. Block number, number of rows, number of columns',
     +    '   2. List of indices for the pivot, each negated if part of'
     +    ,'      a 2x2 pivot',
     +    '   3. The factorized block pivot',
     +    '      It has the form',
     +    '            -1  T',
     +    '        L  D   L ',
     +    '                         -1    T',
     +    '      and is printed as D and L  packed together.',
     +    '   4. List of indices for the non-pivot columns',
     +    '   5. The non-pivot part as rectangular block by rows'

      IWPOS = 4
      APOS = 1

      DO 300 IBLK = 1,NBLK
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

        WRITE (MP,'(/4(A,I6))') 'Block pivot',IBLK,' with',NROWS,
     +        ' rows and', NCOLS,' columns'
        IF (LEN.EQ.12) WRITE (MP,'(6I12)')
     +                       (IFACT(K),K=IWPOS,IWPOS+NROWS-1)
        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (PM(SIGN(1,IFACT(K))),
     +      K=IWPOS,IWPOS+NROWS-1)

        JPIV = 0
        DO 30 IROW = 1,NROWS
          IF (JPIV.EQ.1) THEN
            JPIV = 0
          ELSE
            IF (IFACT(IWPOS+IROW-1).LT.0) JPIV = 1
          END IF

          ILINE = 1
          DO 10 J = 1,IROW - 1
            WRITE (LINE(ILINE:ILINE+LEN-1),'(A)') ' '
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
   10     CONTINUE

          DO 20 J = IROW,NROWS
            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +          '(1P,D12.4)') FACT(APOS)
            IF (LEN.EQ.1) THEN
               IF (FACT(APOS).EQ.ZERO) THEN
                  WRITE (LINE(ILINE:ILINE),'(A)') '.'
               ELSE
                  WRITE (LINE(ILINE:ILINE),'(A)') '*'
               END IF
            END IF
            APOS = APOS + 1
            IF (J.EQ.IROW+1) THEN
              IF (JPIV.EQ.1) THEN
                IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +              '(1P,D12.4)') FACT(APOS2)
                IF (LEN.EQ.1) THEN
                    IF (FACT(APOS2).EQ.ZERO) THEN
                       WRITE (LINE(ILINE:ILINE),'(A)') '.'
                    ELSE
                       WRITE (LINE(ILINE:ILINE),'(A)') '*'
                    END IF
                END IF
                APOS2 = APOS2 + 1
              END IF
            END IF
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
   20     CONTINUE

          IF (ILINE.GT.1) THEN
            LINE(ILINE:) = ' '
            WRITE (MP,'(A)') LINE
          END IF

   30   CONTINUE

        IWPOS = IWPOS + NROWS
        IF (LEN.EQ.12) WRITE (MP,'(6I12)') (IFACT(K),K=IWPOS,
     +      IWPOS+NCOLS-NROWS-1)
        IF (LEN.EQ.1) WRITE (MP,'(72A1)') (PM(SIGN(1,IFACT(K))),
     +      K=IWPOS,IWPOS+NCOLS-NROWS-1)

        IWPOS = IWPOS + NCOLS - NROWS
        DO 280 IROW = 1,NROWS
          J1 = NROWS
          J2 = NCOLS
          ILINE = 1
          DO 110 J = J1 + 1,J2
            IF (LEN.EQ.12) WRITE (LINE(ILINE:ILINE+11),
     +          '(1P,D12.4)') FACT(APOS)
            IF (LEN.EQ.1) THEN
               IF (FACT(APOS).EQ.ZERO) THEN
                  WRITE (LINE(ILINE:ILINE),'(A)') '.'
               ELSE
                  WRITE (LINE(ILINE:ILINE),'(A)') '*'
               END IF
            END IF
            APOS = APOS + 1
            ILINE = ILINE + LEN
            IF (ILINE.GT.72) THEN
              WRITE (MP,'(A)') LINE
              ILINE = 1
            END IF
  110     CONTINUE

          IF (ILINE.GT.1) THEN
            LINE(ILINE:) = ' '
            WRITE (MP,'(A)') LINE
          END IF

  280   CONTINUE
  300 CONTINUE
      END

      SUBROUTINE MA57SD(FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,ICNTL)
C This subroutine divides a vector by the block diagonal matrix of
C     the matrix factors using factor entries stored in FACT/IFACT
C      by MA57BD.
      INTEGER LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER ICNTL(20)
C FACT    must be set to hold the real values corresponding to the
C      factors. This must be unchanged since the
C      preceding call to MA57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C      in FACT. This must be unchanged since the preceding call to
C      MA57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C      the forward substitution operations. On exit, holds the
C      solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C      sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C     larger than INFO(21) as returned from MA57BD.
C ICNTL Not referenced except:
C     ICNTL(13) Threshold on number of columns in a block for direct
C           addressing using Level 2 and Level 3 BLAS.

C Procedures
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV

C
C Local variables.
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,NCOLS,NROWS
      DOUBLE PRECISION W1
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C       pivot.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C JPIV  Has the value 1 for the first row of a 2 by 2 pivot and -1 for
C       the second.
C K     Temporary pointer to position in real array.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
      APOS = 1
      APOS2 = IFACT(1)
      IWPOS = 4
C Run through block pivot rows in the reverse order.
      DO 380 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Multiply by the diagonal matrix (direct addressing)
          DO 10 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 9 J = 1,NRHS
              W(IPIV,J) = RHS(IRHS,J)*FACT(APOS)
    9       CONTINUE
            APOS = APOS + (NROWS+1-IPIV)
   10     CONTINUE
          JPIV = 1
          DO 20 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              DO 19 J = 1,NRHS
                W(IPIV,J) = RHS(IRHS1,J)*FACT(APOS2) + W(IPIV,J)
   19         CONTINUE
              IF (JPIV.EQ.-1) APOS2 = APOS2 + 1
              JPIV = -JPIV
            END IF

   20     CONTINUE

C Reload W back into RHS.
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE

        ELSE
C
C Perform operations using indirect addressing.

C Multiply by the diagonal matrix (indirect addressing)
          JPIV = 1
          DO 210 IPIV = 1,NROWS
            IRHS = IFACT(IWPOS+IPIV-1)

            IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
              DO 65 J = 1,NRHS
                RHS(IRHS,J) = RHS(IRHS,J)*FACT(APOS)
   65         CONTINUE
              APOS = APOS + NROWS - IPIV + 1
            ELSE
C 2 by 2 pivot
              IF (JPIV.EQ.1) THEN
                IRHS1 = -IRHS
                IRHS2 = -IFACT(IWPOS+IPIV)
                DO 68 J = 1,NRHS
                  W1 = RHS(IRHS1,J)*FACT(APOS) +
     +                 RHS(IRHS2,J)*FACT(APOS2)
                  RHS(IRHS2,J) = RHS(IRHS1,J)*FACT(APOS2) +
     +                           RHS(IRHS2,J)*FACT(APOS+NROWS-IPIV+1)
                  RHS(IRHS1,J) = W1
   68           CONTINUE
                APOS2 = APOS2 + 1
              END IF
              JPIV = -JPIV
              APOS = APOS + NROWS - IPIV + 1
            END IF

  210     CONTINUE

        END IF

        IWPOS = IWPOS + NCOLS
        APOS = APOS + NROWS*(NCOLS-NROWS)

  380 CONTINUE

      END

      SUBROUTINE MA57TD(N,FACT,LFACT,IFACT,LIFACT,NRHS,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C     using the factors stored in FACT/IFACT by MA57BD.
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),NRHS,LRHS,LW
      DOUBLE PRECISION W(LW,NRHS),RHS(LRHS,NRHS)
      INTEGER IW1(N),ICNTL(20)
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C      factors. This must be unchanged since the
C      preceding call to MA57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C      in FACT. This must be unchanged since the preceding call to
C      MA57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C NRHS must be set to number of right-hand sides.
C RHS  on entry, must be set to hold the right hand side modified by
C      the forward substitution operations. On exit, holds the
C      solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C      sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C     larger than INFO(21) as returned from MA57BD.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C      number of block pivots, must hold pointers to the beginning of
C      each block pivot in array IFACT, as set by MA57Q\QD.
C      It is not altered.
C ICNTL Not referenced except:
C     ICNTL(13) Threshold on number of columns in a block for direct
C           addressing using Level 2 and Level 3 BLAS.

C Procedures
      INTRINSIC ABS
      EXTERNAL DGEMM,DTPSV

C Constants
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C Local variables.
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,
     +        IWPOS,J,J1,J2,K,KK,NCOLS,NROWS
      DOUBLE PRECISION W1
C APOS  Current position in array FACT.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C       the second.
C K     Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K     Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
      APOS = IFACT(1)

C Set IW1 array
      IWPOS = 4
      DO 10 I = 1,IFACT(3)-1
        IW1(I) = IWPOS
        IWPOS = IWPOS + ABS(IFACT(IWPOS))+2
   10 CONTINUE
      IW1(IFACT(3)) = IWPOS

C Run through block pivot rows in the reverse order.
      DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Load right-hand side into W.
          DO 5 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            DO 3 J = 1,NRHS
              W(I,J) = RHS(II,J)
    3       CONTINUE
    5     CONTINUE


C Treat off-diagonal block (direct addressing)
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMM('T','N',NROWS,NRHS,K,ONE,
     +                           FACT(APOS),K,
     +                           W(NROWS+1,1),LW,ONE,W,LW)

          APOS = APOS-(NROWS*(NROWS+1))/2

C Treat diagonal block (direct addressing)
          DO 22 J = 1,NRHS
            CALL DTPSV('L','T','U',NROWS,FACT(APOS),W(1,J),1)
   22     CONTINUE

C Reload W back into RHS.
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            DO 59 J = 1,NRHS
              RHS(II,J) = W(I,J)
   59       CONTINUE
   60     CONTINUE

        ELSE
C
C Perform operations using indirect addressing.
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1

C Treat off-diagonal block (indirect addressing)
          KK = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            DO 218 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = KK
              DO 215 J = J1,J2
                W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  215         CONTINUE
              RHS(IRHS,II) = W1
  218       CONTINUE
            KK = K
  220     CONTINUE

C Treat diagonal block (indirect addressing)
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            DO 240 II = 1,NRHS
              W1 = RHS(IRHS,II)
              K = APOS + 1
              DO 230 J = J1,J2
                W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)),II)
                K = K + 1
  230         CONTINUE
              RHS(IRHS,II) = W1
  240       CONTINUE
            J1 = J1 - 1
  260     CONTINUE

        END IF

  380 CONTINUE


      END
      SUBROUTINE MA57DD(JOB,N,NE,A,IRN,JCN,FACT,LFACT,IFACT,LIFACT,
     *                  RHS,X,RESID,W,IW,ICNTL,CNTL,INFO,RINFO)

C This subroutine solves a single system using one or more steps of
C     iterative refinement.
C If ICNTL(9) = 10 (the default), this subroutine performs iterative
C     refinement using the strategy of Arioli, Demmel, and Duff.
C     IF (ICNTL(9) = 1, then one step of iterative refinement is
C     performed.

      INTEGER JOB,N,NE
      DOUBLE PRECISION A(NE)
      INTEGER IRN(NE),JCN(NE),LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT)
      DOUBLE PRECISION RHS(N),X(N),RESID(N),W(N,*)
      INTEGER IW(N),ICNTL(20)
      DOUBLE PRECISION CNTL(5)
      INTEGER INFO(40)
      DOUBLE PRECISION RINFO(20)

C JOB must be set by the user to determine what action is desired by
C     the user.
C     Values of JOB and their effect are:
C IF ICNTL(9)>1, JOB=0 if no estimate of solution in X; JOB=2 if
C        estimate of solution in X.
C IF ICNTL(9)=1, then:
C     0: Solve Ax=b, calculate residual r=b-Ax and exit.
C        (Note that MA57CD should be used if solution without residual
C         is required)
C     1: Solve Ax=b, calculate residual r=b-Ax, solve A(dx)=r,
C        update solution and exit.
C If JOB > 1, an estimate of the solution must be input in X.
C     2: Calculate residual r=b-Ax, solve A(dx)=r,
C        update solution and exit.
C If JOB > 2, the residual for this estimate must also be input.
C     3: Solve A(dx)=r, update solution and exit.
C     4: Solve A(dx)=r, update solution, calculate residual for new
C        solution and exit.
C N must be set to the order of the matrix and must be unchanged since
C     the call to MA57BD. It is not altered by MA57DD.
C NE must be set to the number of entries in the matrix and must be
C     unchanged since the call to MA57AD. It is not altered by MA57DD.
C A must be set by the user to the values of the matrix as input to
C     MA57AD. It is not altered by MA57DD.
C IRN,JCN must be set to the row and column indices of the entries in A
C     and must be unchanged since the call to MA57AD. They are
C     not altered by MA57DD.
C FACT holds information on the factors and must be unchanged since
C     the call to MA57BD. It is not altered by MA57DD.
C LFACT must be set to the length of FACT. It is not altered by MA57DD.
C IFACT holds information on the factors and must be unchanged since
C     the call to MA57B/BD. It is not altered by MA57D.
C LIFACT must be set to the length of IFACT. It is not altered by
C     A57D/DD.
C RHS is a real array of length N that must be set to the right-hand
C     side for the equation being solved. It is not altered by MA57DD.
C X is a real array of length N. IF JOB >=2, it must be set on entry to
C     an estimated solution. Otherwise, it need not be set by the user.
C     On exit, the improved solution vector is returned in X.
C RESID is a real array of length N. If JOB>2, it must be set on entry
C     to the value of the residual for the current solution estimate
C     held in X.  Otherwise, it need not be set by the user on entry.
C     If JOB=0 or 4 or if ICNTL(9)>1, on exit it will hold the residual
C     vector for the equations being solved. If 1<= JOB <= 3, then
C     RESID will hold on exit the last correction vector added to the
C     solution X.
C W is used as a work array.  If ICNTL(9) = 1, it must be of length at
C     least N.  If ICNTL(9) > 1, it must be of length at least 3*N if
C     ICNTL(10)=0. If ICNTL(9)>1 and ICNTL(10)>0, then W must be of
C     length at least 4*N.
C IW is an integer array of length N that is used as a work array if
C     ICNTL(1) > 9.  It is not accessed if ICNTL(9) = 1.
C ICNTL must be set by the user as follows and is not altered.
C     ICNTL(1)  must be set to the stream number for error messages.
C       A value less than 0 suppresses output.
C     ICNTL(2) must be set to the stream number for warning output.
C       A value less than 0 suppresses output.
C     ICNTL(3) must be set to the stream number for monitor output.
C       A value less than 0 suppresses output.
C     ICNTL(4) must be set to the stream number for statistics output.
C       A value less than 0 suppresses output.
C     ICNTL(5) must be set to control the amount of output:
C       0 None.
C       1 Error messages only.
C       2 Error and warning messages.
C       3 As 2, plus scalar parameters and a few entries of array
C         parameters on entry and exit.
C     > 3  As 2, plus all parameters on entry and exit.
C     ICNTL(6:12) Not referenced.
C     ICNTL(9)  Maximum permitted number of steps of iterative
C               refinement.
C     ICNTL(10) Flag to request calculation of error estimate and
C           condition numbers.
C     ICNTL(13) Threshold on number of columns in a block for direct
C           addressing using Level 2 and Level 3 BLAS.
C CNTL must be set by the user as follows and is not altered.
C     CNTL(3) is the required decrease in the scaled residuals required
C           by the Arioli, Demmel, and Duff iteration.
C INFO is an integer array that need not be set by the user.  On exit,
C     a value of INFO(1) equal to zero indicates success. A failure is
C     indicated by a negative value for INFO.
C RINFO is a real array that need not be set by the user. On exit,
C     If ICNTL(9)>1, RINFO is set to information on the matrix and
C     solution including backward errors.

C     .. Local constants ..
      DOUBLE PRECISION ZERO,ONE
      PARAMETER (ZERO=0.D0,ONE=1.0D0)

C     .. Local variables ..
      DOUBLE PRECISION COND(2),CTAU,DXMAX,ERROR,OLDOMG(2),OLDOM2,
     *                 OMEGA(2),OM2,TAU
      INTEGER I,ICNTLC(20),ITER,J,K,KASE,KK,LDIAG,LP,MP,KEEP71(5)
      LOGICAL LCOND(2)

C
C COND is condition number of system.
C CTAU is set to 1000*machine precision.
C DXMAX used to calculate max norm of solution.
C ERROR used to accumulate error.
C OLDOMG holds value of previous backward errors.
C OLDOM2 holds previous sum of OMEGAs.
C OMEGA used to accumulate backward error.
C OM2 holds sum of OMEGAs.
C TAU is threshold for denominator in scaled residual calculation.
C I  Temporary variable.
C ICNTLC is control array for MA57C/CD.
C ITER is maximum permitted number of steps of iterative refinement.
C J  Temporary variable.
C K  Temporary variable.
C KASE used when calling MC71AD.
C KK Temporary variable.
C LDIAG Control for amount of information output.
C LCOND used as switch when calculating condition number.
C LP Stream number for error printing.
C MP Stream number for monitor printing.
C KEEP71 Work array required by MC71.
C

C Procedures
      INTRINSIC MIN
      EXTERNAL MA57CD,MA57UD,FD15AD,MC71AD
C EPS is the largest real such that 1+EPS is equal to 1.
      DOUBLE PRECISION EPS,FD15AD


C Set local print control variables
      LP = ICNTL(1)
      MP = ICNTL(3)
      LDIAG = ICNTL(5)

      INFO(1) = 0

C Check input data
      IF (N.LE.0) THEN
        INFO(1) = -1
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'N has value',N
        INFO(2) = N
        GOTO 500
      ENDIF

      IF (NE.LT.0) THEN
        INFO(1) = -2
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I3/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'NE has value',NE
        INFO(2) = NE
        GOTO 500
      ENDIF

      IF (ICNTL(9).LT.1) THEN
        INFO(1) = -13
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'ICNTL(9) has value',ICNTL(9)
        INFO(2) = ICNTL(9)
        GOTO 500
      ENDIF

      IF (JOB.LT.0 .OR. JOB.GT.4 .OR. (ICNTL(9).GT.1 .AND.
     *    (JOB.NE.0 .AND. JOB.NE.2)))  THEN
        INFO(1) = -12
        INFO(2) = JOB
        IF (LDIAG.GT.0 .AND. LP.GE.0) WRITE (LP,'(A,I4/A,I12)')
     +    '**** Error return from MA57DD ****  INFO(1) =',INFO(1),
     +    'JOB has value',JOB
        IF (ICNTL(9).GT.1 .AND. LDIAG.GT.0 .AND. LP.GE.0)
     +    WRITE (LP,'(A,I3)') 'and ICNTL(9) =',ICNTL(9)
        GOTO 500
      ENDIF

C If NE = 0, set variables and return
      IF (NE.EQ.0) THEN
        IF (JOB.NE.3) THEN
          DO 8 I = 1,N
            RESID(I) = ZERO
  8       CONTINUE
        ENDIF
        DO 9 I = 1,N
          X(I) = ZERO
  9     CONTINUE
        INFO(30)=0
        DO 10 I = 6,13
          RINFO(I) = ZERO
 10     CONTINUE
        GO TO 500
      ENDIF

C If requested, print input parameters
      IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE (MP,99999) JOB,N,NE,(ICNTL(I),I=1,5),LFACT,LIFACT,
     +   ICNTL(9),ICNTL(10),ICNTL(13),CNTL(3)
99999 FORMAT(/'Entering iterative refinement solution phase ',
     +  '(MA57DD) with ...'/
     +  'JOB       Control for coefficient matrix      =',I12/
     +  'N         Order of matrix                     =',I12/
     +  'NE        Number of entries in matrix         =',I12/
     6  'ICNTL(1)  Stream for errors                   =',I12/
     7  ' --- (2)  Stream for warnings                 =',I12/
     8  ' --- (3)  Stream for monitoring               =',I12/
     9  ' --- (4)  Stream for statistics               =',I12/
     1  ' --- (5)  Level of diagnostic printing        =',I12/
     +  'LFACT     Length of array FACT                =',I12/
     +  'LIFACT    Length of array IFACT               =',I12/
     +  'ICNTL(9)  Number steps iterative refinement   =',I12/
     +  'ICNTL(10) Control for error analysis          =',I12/
     +  'ICNTL(13) Threshold for Level 2 and 3 BLAS    =',I12/
     +  'CNTL(3)   Convergence test for IR             =',1P,D12.4)

C Print out matrix factors.
        CALL MA57UD(FACT,LFACT,IFACT,LIFACT,ICNTL)
        K = MIN(10,N)
        IF (LDIAG.GE.4) K = N
        WRITE(MP,'(/A)') 'Right-hand side'
        WRITE (MP,'((4X, 1P,5D13.3))') (RHS(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
      END IF


      DO 15 I=1,5
        ICNTLC(I) = ICNTL(I)
   15 CONTINUE
      ICNTLC(13) = ICNTL(13)
      ICNTLC(15) = ICNTL(15)
C Switch off monitor printing in MA57C/CD
      ICNTLC(3) = -1

      IF (JOB.LE.2) THEN
        IF (JOB .LE. 1) THEN
C No estimate given in X.
          DO 14 I = 1,N
            X(I) = RHS(I)
            RESID(I) = RHS(I)
   14     CONTINUE
C Solve system Ax=b using MA57C/CD
          CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,X,N,W,N,IW,
     +                ICNTLC,INFO)
        ELSE
C Estimate of solution was input in X.
          DO 13 I = 1,N
            RESID(I) = RHS(I)
   13     CONTINUE
        ENDIF

        IF (ICNTL(9).EQ.1) THEN
C Compute residual
          DO 16 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 16
            RESID(J) = RESID(J) - A(KK)*X(I)
C Matrix is symmetric
            IF (I.NE.J) RESID(I) = RESID(I) - A(KK)*X(J)
   16     CONTINUE
          IF (JOB.EQ.0) GO TO 340
        ELSE
C Calculate values for iterative refinement strategy of Arioli,
C       Demmel and Duff.
          DO 18 I = 1,N
            W(I,1) = ZERO
C Calculate infinity norms of rows of A in W(I,1)
C Sum |a  |, j=1,N (= ||A  ||        )
C       ij               i.  infinity
            W(I,3) = ZERO
   18     CONTINUE
          DO 17 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 17
            RESID(J) = RESID(J) - A(KK)*X(I)
            W(J,1) = W(J,1) + ABS(A(KK)*X(I))
            W(J,3) = W(J,3) + ABS(A(KK))
C Matrix is symmetric
            IF (I.NE.J) THEN
              RESID(I) = RESID(I) - A(KK)*X(J)
              W(I,1) = W(I,1) + ABS(A(KK)*X(J))
              W(I,3) = W(I,3) + ABS(A(KK))
            ENDIF
   17     CONTINUE

C Calculate max-norm of solution
        DXMAX = ZERO
        DO 221 I = 1,N
          DXMAX = MAX(DXMAX,ABS(X(I)))
  221   CONTINUE

C Initialize EPS
      EPS = FD15AD('E')
C Calculate backward errors OMEGA(1) and OMEGA(2)
C CTAU ... 1000 eps (approx)
        CTAU = 1000.*EPS
C tau is (||A  ||         ||x||   + |b| )*n*1000*epsilon
C            i.  infinity      max     i
          OMEGA(1) = ZERO
          OMEGA(2) = ZERO
          DO 231 I = 1,N
            TAU = (W(I,3)*DXMAX+ABS(RHS(I)))*N*CTAU
            IF ((W(I,1)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C       i               i
              OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
     +                   (W(I,1)+ABS(RHS(I))))
              IW(I) = 1
            ELSE
C TAU will be zero if all zero row in A, for example
              IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||        ||x||   )
C       i        i     i.  infinity     max
                OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
     +                     (W(I,1)+W(I,3)*DXMAX))
              END IF
              IW(I) = 2
            END IF
  231     CONTINUE
C
C  Stop the calculations if the backward error is small
C
          OM2 = OMEGA(1) + OMEGA(2)
          ITER = 0
C Statement changed because IBM SP held quantities in registers
C         IF ((OM2+ONE).LE.ONE) THEN
          IF (OM2.LE.EPS) THEN
            GO TO 270
          ENDIF

C Hold current estimate in case needed later
          DO 251 I = 1,N
            W(I,2) = X(I)
  251     CONTINUE
          OLDOMG(1) = OMEGA(1)
          OLDOMG(2) = OMEGA(2)
          OLDOM2 = OM2

        ENDIF

      ENDIF

C At this point JOB >= 1 or ICNTL(9) > 1
C
C Iterative refinement loop
      DO 260 ITER = 1,ICNTL(9)

C Solve system A(dx) = r using MA57C/CD
        CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,RESID,N,W,N,IW,
     +              ICNTLC,INFO)

C Update solution
        DO 141 I = 1,N
          X(I) = X(I) + RESID(I)
  141   CONTINUE

C Exit without computing residual
        IF (JOB.LT.4 .AND. ICNTL(9).EQ.1) GO TO 340

C
C Calculate residual using information in A,IRN,JCN
C If ICNTL(9).GT.1 also calculate |A| |x|
C
        IF (ICNTL(9).EQ.1) THEN
C Now JOB = 4
          DO 151 I = 1,N
            RESID(I) = RHS(I)
  151     CONTINUE
          DO 181 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 181
            RESID(J) = RESID(J) - A(KK)*X(I)
C Matrix is symmetric
            IF (I.NE.J) RESID(I) = RESID(I) - A(KK)*X(J)
  181     CONTINUE
C Jump because only one step of iterative refinement requested.
          GO TO 340
        ELSE
          DO 153 I = 1,N
C b - Ax
            RESID(I) = RHS(I)
C |A||x|
            W(I,1) = ZERO
  153     CONTINUE
          DO 183 KK = 1,NE
            I = IRN(KK)
            J = JCN(KK)
            IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 183
            RESID(J) = RESID(J) - A(KK)*X(I)
            W(J,1) = W(J,1) + ABS(A(KK)*X(I))
C Matrix is symmetric
            IF (I.NE.J) THEN
              RESID(I) = RESID(I) - A(KK)*X(J)
              W(I,1) = W(I,1) + ABS(A(KK)*X(J))
            ENDIF
  183     CONTINUE
        ENDIF

C Calculate max-norm of solution
        DXMAX = ZERO
        DO 220 I = 1,N
          DXMAX = MAX(DXMAX,ABS(X(I)))
  220   CONTINUE

C Calculate omega(1) and omega(2)
C tau is (||A  ||         ||x||   + |b| )*n*1000*epsilon
C            i.  infinity      max     i
        OMEGA(1) = ZERO
        OMEGA(2) = ZERO
        DO 230 I = 1,N
          TAU = (W(I,3)*DXMAX+ABS(RHS(I)))*N*CTAU
          IF ((W(I,1)+ABS(RHS(I))).GT.TAU) THEN
C |Ax-b| /(|A||x| + |b|)
C       i               i
            OMEGA(1) = MAX(OMEGA(1),ABS(RESID(I))/
     +                 (W(I,1)+ABS(RHS(I))))
            IW(I) = 1
          ELSE
C TAU will be zero if all zero row in A, for example
            IF (TAU.GT.ZERO) THEN
C |Ax-b| /(|A||x| + ||A  ||        ||x||   )
C       i        i     i.  infinity     max
              OMEGA(2) = MAX(OMEGA(2),ABS(RESID(I))/
     +                   (W(I,1)+W(I,3)*DXMAX))
            END IF
            IW(I) = 2
          END IF
  230   CONTINUE
C
C  Stop the calculations if the backward error is small
        OM2 = OMEGA(1) + OMEGA(2)
        IF ((OM2+ONE).LE.ONE) THEN
          GO TO 270
        ENDIF
C
C  Check the convergence.
C
        IF (OM2.GT.OLDOM2*CNTL(3)) THEN
C  Stop if insufficient decrease in omega.
          IF (OM2.GT.OLDOM2) THEN
C Previous estimate was better ... reinstate it.
            OMEGA(1) = OLDOMG(1)
            OMEGA(2) = OLDOMG(2)
            DO 240 I = 1,N
              X(I) = W(I,2)
  240       CONTINUE
          END IF
          GO TO 270
        ELSE
C Hold current estimate in case needed later
          DO 250 I = 1,N
            W(I,2) = X(I)
  250     CONTINUE
          OLDOMG(1) = OMEGA(1)
          OLDOMG(2) = OMEGA(2)
          OLDOM2 = OM2
        END IF

  260 CONTINUE
C End of iterative refinement loop.

      INFO(1) = -8
      IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9170) INFO(1),ICNTL(9)
 9170 FORMAT ('Error return from MA57D/DD because of ','nonconvergence',
     +       ' of iterative refinement'/'Error INFO(1) = ',I2,'  with',
     +       ' ICNTL','(9) = ',I10)

C Set the RINFO parameters
  270 RINFO(6)  = OMEGA(1)
      RINFO(7)  = OMEGA(2)
      RINFO(8) = ZERO
      DO 271 I=1,N
        RINFO(8) = MAX(RINFO(8),W(I,3))
  271 CONTINUE
      RINFO(9) = DXMAX
      RINFO(10) = ZERO
      DO 272 I=1,N
        RINFO(10) = MAX(RINFO(10),ABS(RESID(I)))
  272 CONTINUE
      IF (RINFO(8)*RINFO(9).NE.ZERO)
     *RINFO(10) = RINFO(10)/(RINFO(8)*RINFO(9))
      INFO(30) = ITER

      IF (INFO(1).LT.0) GO TO 340

C Jump if estimate of error not requested.
      IF (ICNTL(10).LE.0) GO TO 340

C
C Calculate condition numbers and estimate of the error.
C
C Condition numbers obtained through use of norm estimation
C     routine MC71A/AD.
C
C Initializations
C
      LCOND(1) = .FALSE.
      LCOND(2) = .FALSE.
      ERROR    = ZERO
      DO 280 I = 1,N
        IF (IW(I).EQ.1) THEN
          W(I,1) = W(I,1) + ABS(RHS(I))
C |A||x| + |b|
          W(I,2) = ZERO
          LCOND(1) = .TRUE.
        ELSE
C |A||x| + ||A  ||        ||x||
C             i.  infinity     max

          W(I,2) = W(I,1) + W(I,3)*DXMAX
          W(I,1) = ZERO
          LCOND(2) = .TRUE.
        END IF
  280 CONTINUE
C
C  Compute the estimate of COND
C
      DO 330 K = 1,2
        IF (LCOND(K)) THEN
C MC71A/AD has its own built in limit of 5 to the number of iterations
C    allowed. It is this limit that will be used to terminate the
C    following loop.
          KASE = 0
          DO 310 KK = 1,40

C MC71A/AD calculates norm of matrix
C We are calculating the infinity norm of INV(A).W

C Initialize W(1,3).  W(1,4) is a work array.
            CALL MC71AD(N,KASE,W(1,3),COND(K),W(1,4),IW,KEEP71)
C
C  KASE = 0........ Computation completed
C  KASE = 1........ W * INV(TRANSPOSE(A)) * Y
C  KASE = 2........ INV(A) * W * Y
C                   W is W/W(*,2) .. Y is W(*,3)
C
            IF (KASE.EQ.0) GO TO 320

            IF (KASE.EQ.1) THEN
C Solve system using MA57C/CD.
C W(1,4) is used as workspace
              CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,W(1,3),
     *                    N,W(1,4),N,IW,ICNTLC,INFO)
              DO 290 I = 1,N
                W(I,3) = W(I,K)*W(I,3)
  290         CONTINUE
            END IF

            IF (KASE.EQ.2) THEN
              DO 300 I = 1,N
                W(I,3) = W(I,K)*W(I,3)
  300         CONTINUE
C Solve system using MA57C/CD.
              CALL MA57CD(1,N,FACT,LFACT,IFACT,LIFACT,1,W(1,3),N,
     *                    W(1,4),N,IW,ICNTLC,INFO)
            END IF

  310     CONTINUE

C Error return if MC71AD does not converge
          INFO(1) = -14
          IF (LP.GE.0 .AND. LDIAG.GE.1) WRITE (LP,9160)
 9160 FORMAT ('Error return from MA57D/DD because of ','error in MC71',
     +       'A/AD'/'Error not calculated')

  320     IF (DXMAX.GT.ZERO) COND(K) = COND(K)/DXMAX
          ERROR = ERROR + OMEGA(K)*COND(K)
        ELSE
          COND(K) = ZERO
        ENDIF

  330 CONTINUE

      RINFO(11)  = COND(1)
      RINFO(12)  = COND(2)
      RINFO(13)  = ERROR

C
C If requested, print output parameters.
 340  IF (LDIAG.GE.3 .AND. MP.GE.0) THEN
        WRITE(MP,99980) INFO(1)
99980 FORMAT (/'Leaving iterative refinement solution phase ',
     +  '(MA57DD) with ...'/
     1      'INFO (1)                                      =',I12/)
        IF (INFO(1).LT.0) GO TO 500
        IF (ICNTL(9).GT.1) THEN
          WRITE(MP,99981) INFO(30),(RINFO(I),I=6,10)
99981     FORMAT(
     1     'INFO(30)  Number steps iterative ref   =',I10/
     1     'RINFO(6)  Backward errors  (OMEGA(1))  =',1PD10.3/
     2     '-----(7)  Backward errors  (OMEGA(2))  =',1PD10.3/
     3     '-----(8)  Infinity norm of matrix      =',1PD10.3/
     4     '-----(9)  Infinity norm of solution    =',1PD10.3/
     5     '-----(10) Norm of scaled residuals     =',1PD10.3)
          IF (ICNTL(10).GT.0) WRITE(MP,99979) (RINFO(I),I=11,13)
99979       FORMAT (
     1       'RINFO(11) Condition number (COND(1))   =',1PD10.3/
     1       'RINFO(12) Condition number (COND(2))   =',1PD10.3/
     1       'RINFO(13) Error in solution            =',1PD10.3)
          WRITE(MP,'(/A,I10)') 'Residual'
          K=MIN(N,10)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(1P,5D13.3)') (RESID(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
        ELSE
C ICNTL(9) = 1
          IF (JOB.GE.1 .AND. JOB.LE.3) THEN
            WRITE(MP,'(/A,I10)') 'Correction to solution'
          ELSE
            WRITE(MP,'(/A,I10)') 'Residual'
          ENDIF
          K=MIN(N,10)
          IF (LDIAG.GE.4) K = N
          WRITE (MP,'(1P,5D13.3)') (RESID(I),I=1,K)
          IF (K.LT.N) WRITE (MP,'(A)') '     . . .'
        END IF

C Print solution
        K=MIN(N,10)
        IF (LDIAG.GE.4) K = N
        WRITE(MP,'(/A,I10)') 'Solution'
        WRITE (MP,'(1P,5D13.3)') (X(I),I=1,K)
        IF (K.LT.N) WRITE (MP,'(A)') '     . . .'

      END IF

 500  RETURN

      END
      SUBROUTINE MA57ED(N,IC,KEEP,FACT,LFACT,NEWFAC,LNEW,
     *                  IFACT,LIFACT,NEWIFC,LINEW,INFO)
      INTEGER N,IC,KEEP(*),LFACT,LNEW,LIFACT,LINEW,INFO(40)
      DOUBLE PRECISION FACT(LFACT),NEWFAC(LNEW)
      INTEGER IFACT(LIFACT),NEWIFC(LINEW)
C Local variables
      INTEGER APOSBB,ASTK,HOLD,I,ISTK,IWPOS,MOVE,NFRONT

C HOLD determines part of keep holding saved variables from MA57OD.
      HOLD = N + 3

C Initialize INFO(1) and INFO(2)
      INFO(1) = 0
      INFO(2) = 0

C Test to see whether to map real or integer space
      IF (IC.GE.1) THEN
C Remap of integer space
        IF (LINEW.LE.LIFACT) THEN
          INFO(1) = -7
          INFO(2) = LINEW
          RETURN
        ENDIF
        IWPOS = KEEP(HOLD+7)
        ISTK  = KEEP(HOLD+14)
        NFRONT = KEEP(HOLD+23)
        DO 10 I = 1,IWPOS+NFRONT-1
          NEWIFC(I) = IFACT(I)
   10   CONTINUE
C Move distance
        MOVE = LINEW - LIFACT
        DO 20 I = ISTK+1,LIFACT
          NEWIFC(I+MOVE) = IFACT(I)
   20   CONTINUE
C Reset INPUT, ISTK, PTIRN
          KEEP(HOLD+13) = KEEP(HOLD+13) + MOVE
          KEEP(HOLD+14) = ISTK + MOVE
          KEEP(HOLD+18) = KEEP(HOLD+18) + MOVE
      ENDIF
      IF (IC.NE.1) THEN
C Remap of real space
        IF (LNEW.LE.LFACT) THEN
          INFO(1) = -7
          INFO(2) = LNEW
          RETURN
        ENDIF
C Was .. APOS = KEEP(HOLD+8)
        APOSBB = KEEP(HOLD+9)
        ASTK   = KEEP(HOLD+15)
        DO 60 I = 1, APOSBB-1
          NEWFAC(I) = FACT(I)
   60   CONTINUE
C Move distance
        MOVE = LNEW - LFACT
        DO 70 I = ASTK+1,LFACT
          NEWFAC(I+MOVE) = FACT(I)
   70   CONTINUE
C Reset AINPUT, ASTK, PTRA
        KEEP(HOLD+12) = KEEP(HOLD+12) + MOVE
        KEEP(HOLD+15) = ASTK + MOVE
        KEEP(HOLD+19) = KEEP(HOLD+19) + MOVE
      ENDIF
      RETURN
      END
      SUBROUTINE MA57GD(N,NE,IRN,JCN,IW,IPE,COUNT,FLAG,IWFR,ICNTL,INFO)

C Given the positions of the entries of a symmetric matrix, construct
C     the sparsity pattern of the whole matrix. Either one of a pair
C     (I,J),(J,I) may be used to represent the pair. Duplicates are
C     ignored.

      INTEGER N,NE,IRN(NE),JCN(NE),IW(NE*2+N),IPE(N),COUNT(N),
     +        FLAG(N),IWFR,ICNTL(20),INFO(40)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K),K=1,2,...,NE must be set to the row indices of the entries on
C     input.  IRN is not changed.
C JCN(K),K=1,2,...,NE must be set to the column indices of the entries
C     on input.  JCN is not changed.
C IW need not be set on input. On output it contains lists of
C     column indices.
C IPE need not be set on input. On output IPE(I) points to the start of
C     the entry in IW for row I, I=1,2,...,N.
C COUNT need not be set. On output COUNT(I), I=1,2,..,N, contains the
C     number of off-diagonal entries in row I excluding duplicates.
C FLAG is used for workspace to hold flags to permit duplicate entries
C     to be identified quickly.
C IWFR need not be set on input. On output it points to the first
C     unused location in IW.
C ICNTL Warning messages are printed on stream number ICNTL(2) if
C      ICNTL(2).GT.0 and ICNTL(5).GT.1.
C INFO need not be set on input. On output,
C  INFO(1) has one of the values:
C     0 No out-of-range index or duplicate entry found.
C     1 Out-of-range index found.
C     2 Duplicate entry found.
C     3 Out-of-range index found and duplicate entry found.
C  INFO(3) is set to the number of out-of-range entries.
C  INFO(4) is set to the number of off-diagonal duplicate entries.
C
      INTRINSIC MAX,MIN
C
C Local variables
      INTEGER I,J,K,L,LDIAG,WP
C I Row index
C J Column index
C K Position in IRN, JCN, or IW.
C L Position in IW.
C LDIAG Level of diagnostic printing.
C WP Stream for printing warning messages.

C Set LDIAG and WP
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      IF (WP.LT.0) LDIAG = 0

C Initialize INFO(1)
      INFO(1) = 0

C Count in INFO(3) the number of out-of-range entries, initialize FLAG,
C    and count in COUNT the numbers of entries in the rows including
C    duplicates.
      INFO(3) = 0
      DO 10 I = 1,N
        FLAG(I) = 0
        COUNT(I) = 0
   10 CONTINUE
      DO 20 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          INFO(3) = INFO(3) + 1
          INFO(1) = 1
          IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
     +        '*** Warning message from subroutine MA57AD ***',
     +        ' INFO(1) =',INFO(1)
          IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
     +         K,'th entry (in row',I,' and column',J,') ignored'

        ELSE IF (I.NE.J) THEN
          COUNT(I) = COUNT(I) + 1
          COUNT(J) = COUNT(J) + 1
        END IF

   20 CONTINUE
C
C Accumulate row counts in IPE which is set so that IPE(I) points to
C     position after end of row I.
      IPE(1) = COUNT(1)+1
      DO 30 I = 2,N
        IPE(I) = IPE(I-1) + COUNT(I)
   30 CONTINUE
C
C Run through putting the matrix entries in the right places. IPE is
C     used for holding running pointers and is left holding pointers to
C     row starts.
      DO 40 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N .OR. I.EQ.J) GO TO 40
        IPE(I) = IPE(I) - 1
        IW(IPE(I)) = J
        IPE(J) = IPE(J) - 1
        IW(IPE(J)) = I
   40 CONTINUE
C
C Remove duplicates.
      INFO(4) = 0
C IWFR points to the current position in the compressed set of rows.
      IWFR = 1
C At the start of cycle I of this loop FLAG(J).LT.I for J=1,2,...,N.
C     During the loop FLAG(J) is set to I if an entry in column J is
C     found. This permits duplicates to be recognized quickly.
      DO 60 I = 1,N
        L = IPE(I)
        IPE(I) = IWFR
        DO 50 K = L,L+COUNT(I)-1
          J = IW(K)

          IF (FLAG(J).NE.I) THEN
            FLAG(J) = I
            IW(IWFR) = J
            IWFR = IWFR + 1
          ELSE
C Count duplicates only once.
            IF (I.LT.J) INFO(4) = INFO(4) + 1
          END IF

   50   CONTINUE
C Set COUNT to number without duplicates.
        COUNT(I) = IWFR - IPE(I)
   60 CONTINUE
C
C Test whether duplicates found
      IF (INFO(4).GT.0) THEN
        INFO(1) = INFO(1) + 2
        IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
     +      '*** Warning message from subroutine MA57AD ***',INFO(4),
     +      ' off-diagonal duplicate entries found'
      END IF

      END
      SUBROUTINE MA57JD(N,NE,IRN,JCN,PERM,IW,IPE,COUNT,FLAG,IWFR,
     +                  ICNTL,INFO)

C Given the positions of the entries of a symmetric matrix and a
C     permutation, construct the sparsity pattern of the upper
C     triangular part of the matrix. Either one of a pair
C     (I,J),(J,I) may be used to represent the pair. Duplicates are
C     ignored.

      INTEGER N,NE,IRN(NE),JCN(NE),IW(NE+N),IPE(N),COUNT(N),
     +        PERM(N),FLAG(N),IWFR,ICNTL(20),INFO(40)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K),K=1,2,...,NE must be set to the row indices of the entries on
C     input. If IRN(K) or JCN(K) is out of range, IRN(K) is replaced by
C     zero. Otherwise, IRN(K) is not changed.
C JCN(K),K=1,2,...,NE must be set to the column indices of the entries
C     on input. If IRN(K) or JCN(K) is out of range, JCN(K) is replaced
C     by zero. Otherwise, JCN(K) is not changed.
C PERM must be set so that PERM(I) holds the position of variable I
C     in the permuted order.  Its validity as a permutation will have
C     been checked in MA57A/AD.
C IW need not be set on input. On output it contains lists of
C     column indices, each list being headed by its length.
C IPE need not be set on input. On output IPE(I) points to the start of
C     the entry in IW for row I, I=1,2,...,N.
C COUNT need not be set. On output COUNT(I), I=1,2,..,N, contains the
C     number of off-diagonal entries in row I including duplicates.
C     COUNT(0) contains the number of entries with one or both indices
C     out of range.
C FLAG is used for workspace to hold flags to permit duplicate entries
C     to be identified quickly.
C IWFR need not be set on input. On output it points to the first
C     unused location in IW.
C ICNTL Warning messages are printed on stream number ICNTL(2) if
C      ICNTL(2).GT.0 and ICNTL(3).GT.1.
C INFO need not be set on input. On output, INFO(1) has one of the
C     values:
C     0 No out-of-range index or duplicate entry found.
C     1 Out-of-range index found.
C     2 Duplicate entry found.
C     3 Out-of-range index found and duplicate entry found.
C  INFO(3) is set to the number of faulty entries.
C  INFO(4) is set to the number of off-diagonal duplicate entries.

      INTRINSIC MAX,MIN
C
C Local variables
      INTEGER I,J,K,L,LDIAG,WP
C I Row index
C J Column index
C K Position in IRN, JCN, or IW.
C L Position in IW.
C LDIAG Level of monitor printing.
C WP Stream for printing.

C Set LDIAG and WP
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      IF (WP.LT.0) LDIAG = 0

C Initialize INFO(1), FLAG, and COUNT.
      INFO(1) = 0
      DO 10 I = 1,N
        FLAG(I) = 0
        COUNT(I) = 0
   10 CONTINUE

C Count in INFO(3) the number of out-of-range entries, initialize FLAG,
C    and count in COUNT the numbers of entries in the rows.
      INFO(3) = 0
      DO 30 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          IRN(K) = 0
          JCN(K) = 0
          INFO(3) = INFO(3) + 1
          INFO(1) = 1
          IF (INFO(3).EQ.1 .AND. LDIAG.GT.1) WRITE (WP,'(2A,I2)')
     +        '*** Warning message from subroutine MA57AD ***',
     +        ' INFO(1) =',INFO(1)
          IF (INFO(3).LE.10 .AND. LDIAG.GT.1) WRITE (WP,'(3(I10,A))')
     +        K,'th entry (in row',I,' and column',J,') ignored'

        ELSE IF (PERM(I).LE.PERM(J)) THEN
          COUNT(I) = COUNT(I) + 1
        ELSE
          COUNT(J) = COUNT(J) + 1
        END IF

   30 CONTINUE
C
C Accumulate row counts in IPE ... one added for row length location.
      IPE(1) = COUNT(1) + 1
      DO 40 I = 2,N
        IPE(I) = IPE(I-1) + COUNT(I) + 1
   40 CONTINUE

C Run through putting the matrix entries in the right places. IPE is
C     used for holding running pointers and is left holding pointers to
C     row starts.
      DO 50 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) GO TO 50
        IF (PERM(I).LE.PERM(J)) THEN
          IW(IPE(I)) = J
          IPE(I) = IPE(I) - 1
        ELSE
          IW(IPE(J)) = I
          IPE(J) = IPE(J) - 1
        END IF
   50 CONTINUE

C Remove duplicates
C IWFR points to the current position in the compressed set of rows.
      IWFR = 1
      INFO(4) = 0

C At the start of cycle I of this loop FLAG(J).LT.I for J=1,2,...,N.
C     During the loop FLAG(J) is set to I if an entry in column J is
C     found. This permits duplicates to be recognized quickly.
      DO 70 I = 1,N
        L = IPE(I)
        IPE(I) = IWFR

        DO 60 K = L + 1,L + COUNT(I)
          J = IW(K)
          IF (FLAG(J).NE.I) THEN
            FLAG(J) = I
            IWFR = IWFR + 1
            IW(IWFR) = J
          ELSE
C Count duplicates only once.
            IF (I.LT.J) INFO(4) = INFO(4) + 1
          END IF
   60   CONTINUE

        IF (IWFR.GT.IPE(I)) THEN
          IW(IPE(I)) = IWFR - IPE(I)
          IWFR = IWFR + 1
        ELSE
          IPE(I) = 0
        END IF

   70 CONTINUE

C Test whether duplicates found
      IF (INFO(4).GT.0) THEN
        INFO(1) = INFO(1) + 2
        IF (LDIAG.GT.1 .AND. WP.GE.0) WRITE (WP,'(A/I10,A)')
     +      '*** Warning message from subroutine MA57AD ***',
     +      INFO(4),' off-diagonal duplicate entries found'
      END IF

      END

      SUBROUTINE MA57KD(N, IPE, IW, LW, IWFR, PERM, IPS, NV, FLAG,
     *                  NCMPA)
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
      INTEGER IW(LW), PERM(N), IPS(N), NV(N), FLAG(N)
C
C Using a given pivotal sequence and a representation of the matrix that
C     includes only non-zeros of the strictly upper-triangular part
C     of the permuted matrix, construct tree pointers.
C
C N must be set to the matrix order. It is not altered.
C IPE(I) must be set to point to the position in IW of the
C     start of row I or have the value zero if row I has no off-
C     diagonal non-zeros. during execution it is used as follows.
C     If variable I is eliminated then IPE(I) points to the list
C     of variables for created element I. If element I is
C     absorbed into newly created element J then IPE(I)=-J.
C IW must be set on entry to hold lists of variables by
C     rows, each list being headed by its length. when a variable
C     is eliminated its list is replaced by a list of variables
C     in the new element.
C LW must be set to the length of IW. It is not altered.
C IWFR must be set to the position in IW of the first free variable.
C     It is revised during execution, continuing to have this meaning.
C PERM(K) must be set to hold the position of variable K in the
C     pivot order. It is not altered.
C IPS(I) need not be set by the user and will be used to hold the
C     inverse permutation to PERM.
C NV need not be set. If variable J has not been eliminated then
C     the last element whose leading variable (variable earliest
C     in the pivot sequence) is J is element NV(J). If element J
C     exists then the last element having the same leading
C     variable is NV(J). In both cases NV(J)=0 if there is no such
C     element. If element J has been merged into a later element
C     then NV(J) is the degree at the time of elimination.
C FLAG is used as workspace for variable flags.
C     FLAG(JS)=ME if JS has been included in the list for ME.
C
      INTEGER I,J,ML,MS,ME,IP,MINJS,IE,KDUMMY,JP
      INTEGER LN,JP1,JS,LWFR,JP2,JE

      EXTERNAL MA57FD

C
C Initializations
      DO 10 I=1,N
        FLAG(I) = 0
        NV(I) = 0
        J = PERM(I)
        IPS(J) = I
   10 CONTINUE
      NCMPA = 0
C
C Start of main loop
C
      DO 100 ML=1,N
C ME=MS is the name of the variable eliminated and
C     of the element created in the main loop.
        MS = IPS(ML)
        ME = MS
        FLAG(MS) = ME
C
C Merge row MS with all the elements having MS as leading variable.
C IP points to the start of the new list.
        IP = IWFR
C MINJS is set to the position in the order of the leading variable
C     in the new list.
        MINJS = N
        IE = ME
        DO 70 KDUMMY=1,N
C Search variable list of element IE.
C JP points to the current position in the list being searched.
          JP = IPE(IE)
C LN is the length of the list being searched.
          LN = 0
          IF (JP.LE.0) GO TO 60
          LN = IW(JP)
C
C Search for different variables and add them to list,
C     compressing when necessary
          DO 50 JP1=1,LN
            JP = JP + 1
C Place next variable in JS.
            JS = IW(JP)
C Jump if variable has already been included.
            IF (FLAG(JS).EQ.ME) GO TO 50
            FLAG(JS) = ME
            IF (IWFR.LT.LW) GO TO 40
C Prepare for compressing IW by adjusting pointer to and length of
C     the list for IE to refer to the remaining entries.
            IPE(IE) = JP
            IW(JP) = LN - JP1
C Compress IW.
            CALL MA57FD(N, IPE, IW, IP-1, LWFR, NCMPA)
C Copy new list forward
            JP2 = IWFR - 1
            IWFR = LWFR
            IF (IP.GT.JP2) GO TO 30
            DO 20 JP=IP,JP2
              IW(IWFR) = IW(JP)
              IWFR = IWFR + 1
   20       CONTINUE
   30       IP = LWFR
            JP = IPE(IE)
C Add variable JS to new list.
   40       IW(IWFR) = JS
            MINJS = MIN0(MINJS,PERM(JS)+0)
            IWFR = IWFR + 1
   50     CONTINUE
C Record absorption of element IE into new element.
   60     IPE(IE) = -ME
C Pick up next element with leading variable MS.
          JE = NV(IE)
C Store degree of IE.
          NV(IE) = LN + 1
          IE = JE
C Leave loop if there are no more elements.
          IF (IE.EQ.0) GO TO 80
   70   CONTINUE
   80   IF (IWFR.GT.IP) GO TO 90
C Deal with null new element.
        IPE(ME) = 0
        NV(ME) = 1
        GO TO 100
C Link new element with others having same leading variable.
   90   MINJS = IPS(MINJS)
        NV(ME) = NV(MINJS)
        NV(MINJS) = ME
C Move first entry in new list to end to allow room for length at
C     front. Set pointer to front.
        IW(IWFR) = IW(IP)
        IW(IP) = IWFR - IP
        IPE(ME) = IP
        IWFR = IWFR + 1
  100 CONTINUE

      RETURN
      END
C** end of MA57KD**

      SUBROUTINE MA57FD(N, IPE, IW, LW, IWFR, NCMPA)
C Compress lists held by MA57KD in IW and adjust pointers
C     in IPE to correspond.
      INTEGER N,LW,IWFR,NCMPA
      INTEGER IPE(N)
C
      INTEGER   IW(LW)
C N is the matrix order. It is not altered.
C IPE(I) points to the position in IW of the start of list I or is
C     zero if there is no list I. On exit it points to the new position.
C IW holds the lists, each headed by its length. On output the same
C     lists are held, but they are now compressed together.
C LW holds the length of IW. It is not altered.
C IWFR need not be set on entry. On exit it points to the first free
C     location in IW.
C
      INTEGER I,K1,LWFR,IR,K,K2
      NCMPA = NCMPA + 1
C Prepare for compressing by storing the lengths of the
C     lists in IPE and setting the first entry of each list to
C     -(list number).
      DO 10 I=1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
C
C Compress
C IWFR points just beyond the end of the compressed file.
C LWFR points just beyond the end of the uncompressed file.
      IWFR = 1
      LWFR = IWFR
      DO 60 IR=1,N
        IF (LWFR.GT.LW) GO TO 70
C Search for the next negative entry.
        DO 20 K=LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
C Pick up entry number, store length in new position, set new pointer
C     and prepare to copy list.
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
C Copy list to new position.
        DO 40 K=K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN
      END
C--------------------------------------------------------------------
C            HSL 2000
C        --
C-             Copyright CCLRC Rutherford Appleton Laboratory
C        --
C--------------------------------------------------------------------
      SUBROUTINE MA57LD(N, IPE, NV, IPS, NE, NA, NODE, PERM, NSTEPS,
     *                  FILS, FRERE, ND, NEMIN, SUBORD)
C
C Tree search
C
C Given son to father tree pointers, reorder so that eldest son has
C     smallest degree and perform depth-first
C     search to find pivot order and number of eliminations
C     and assemblies at each stage.
      INTEGER N, NSTEPS
      INTEGER ND(N)
      INTEGER IPE(N), FILS(N), FRERE(N), SUBORD(N)
C
      INTEGER NV(N), IPS(N), NE(N), NA(N), NODE(N), PERM(N)
      INTEGER NEMIN
C N must be set to the matrix order. It is not altered.
C IPE(I) must be set equal to -(father of node I) or zero if
C      node is a root, if NV(I) > 0. If NV(I) = 0, then I is
C      subordinate variable of a supervariable and -IPE(I) points to
C      principal variable.  It is altered to point to its next
C      younger brother if it has one, but otherwise is not changed.
C NV(I) must be set to zero if variable is a subordinate variable
C      of a supervariable and to the degree otherwise.
C      NV is not altered.
C IPS(I) need not be set. It is used temporarily to hold
C      -(eldest son of node I) if it has one and 0 otherwise. It is
C      finally set to hold the position of node I in the order.
C NE(IS) need not be set. It is set to the number of variables
C      eliminated at stage IS of the elimination.
C NA(IS) need not be set. It is set to the number of elements
C      assembled at stage IS of the elimination.
C NODE (I) need not be set before entry. It is used during the code
C    to hold the number of subordinate variables for variable I and
C    on output it holds
C    the node (in dfs ordering) at which variable I is eliminated.
C    It is also defined for subordinate variables.
C PERM is set to the new permutation after dfs of tree.  PERM(I) is
C    the position of variable I in the pivot order.
C ND(IS) need not be set. It is set to the degree at stage IS of
C     the elimination.
C NSTEPS need not be set. It is set to the number of elimination steps.
C NEMIN is used to control the amalgamation process between
C       a son and its father (if the number of fully summed
C       variables of both nodes is smaller than NEMIN).
C SUBORD(I) need not be set. It holds the first subordinate variable
C       for variable I if I
C       is a principal variable and holds the next subordinate variable
C       if otherwise.  It is zero at the end of the chain.
C
      INTEGER I,IF,IS,NR,NR1,INS,INL,INB,INF,INFS,INSW
      INTEGER K,L,ISON,IN,IFSON,INO
      INTEGER INOS,IB,IL,INT
      INTEGER IPERM

C
C Initialize IPS and NE.
      DO 10 I=1,N
        IPS(I) = 0
        NE(I) = 0
        NODE(I) = 0
        SUBORD(I) = 0
   10 CONTINUE
C
C Set IPS(I) to -(eldest son of node I) and IPE(I) to next younger
C     brother of node I if it has one.
      NR = N + 1
      DO 50 I=1,N
        IF = -IPE(I)
        IF (NV(I).EQ.0) THEN
C I is a subordinate node, principal variable is IF
          IF (SUBORD(IF).NE.0) SUBORD(I) = SUBORD(IF)
          SUBORD(IF) = I
          NODE(IF) = NODE(IF)+1
        ELSE
C Node IF is the father of node I.
          IF (IF.NE.0) THEN
C IS is younger brother of node I.
C IPS(IF) will eventually point to - eldest son of IF.
            IS = -IPS(IF)
            IF (IS.GT.0) IPE(I) = IS
            IPS(IF) = -I
          ELSE
C I is a root node
            NR = NR - 1
            NE(NR) = I
          ENDIF
        ENDIF
   50 CONTINUE
C
C We reorganize the tree so that the eldest son has maximum number of
C variables.  We combine nodes when the number of variables in a son
C is greater than or equal to the number of variables in the father.
C If the eldest son has the maximum number of variables,
C and if a combination is possible, it has to be possible with
C the eldest son.
C
C FILS is just used as workspace during this reorganization and is reset
C afterwards.

      DO 999 I=1,N
       FILS(I) = IPS(I)
 999  CONTINUE

      NR1 = NR
      INS = 0
C Jump if all roots processed.
 1000 IF (NR1.GT.N) GO TO 1151
C Get next root
      INS = NE(NR1)
      NR1 = NR1 + 1
C Depth first search through eldest sons.
 1070 INL = FILS(INS)
      IF (INL.LT.0) THEN
       INS = -INL
       GO TO 1070
      ENDIF
C INS is leaf node.

 1080 IF (IPE(INS).LT.0) THEN
C INS is youngest son otherwise IPE value would be positive.
       INS       = -IPE(INS)
C INS is now the father of the reorganized son so we can
C     clear the pointer to the sons.
       FILS(INS) = 0
C Continue backtracking until we encounter node with younger brother.
       GO TO 1080
      ENDIF

      IF (IPE(INS).EQ.0) THEN
C INS is a root, check for next one.
       INS = 0
       GO TO 1000
      ENDIF
C INB is younger brother of INS.
      INB = IPE(INS)

C?? I think this test is the wrong way round
      IF (NV(INB).GE.NV(INS)) THEN
C?? So reversed
C     IF (NV(INS).GE.NV(INB)) THEN
       INS = INB
C Do depth first search from younger brother
       GO TO 1070
      ENDIF
C
C Exchange INB and INS
C Find previous brother of INS (could be the father)
C then we do depth first search with INS = INB
C
      INF = INB
 1090 INF = IPE(INF)
      IF (INF.GT.0) GO TO 1090
C -INF IS THE FATHER
      INF  = -INF
      INFS = -FILS(INF)
C INFS is eldest son of INF
      IF (INFS.EQ.INS) THEN
C  INS is eldest brother .. a role which INB now assumes
        FILS(INF) = -INB
        IPS(INF)  = -INB
        IPE(INS)  = IPE(INB)
        IPE(INB)  = INS
      ELSE
        INSW = INFS
 1100   INFS = IPE(INSW)
        IF (INFS.NE.INS) THEN
          INSW = INFS
          GO TO 1100
        ENDIF
        IPE(INS) = IPE(INB)
        IPE(INB) = INS
        IPE(INSW)= INB
      ENDIF
        INS      = INB
C Depth first search from moved younger brother
        GO TO 1070
C Set FRERE and FILS
 1151 DO 51 I=1,N
       FRERE(I) = IPE(I)
       FILS(I) = IPS(I)
 51   CONTINUE
C
C Depth-first search.
C IL holds the current tree level. Roots are at level N, their sons
C     are at level N-1, etc.
C IS holds the current elimination stage. We accumulate the number
C     of eliminations at stage is directly in NE(IS). The number of
C     assemblies is accumulated temporarily in NA(IL), for tree
C     level IL, and is transferred to NA(IS) when we reach the
C     appropriate stage IS.
      IS = 1
C I is the current node.
      I = 0
C IPERM is used as pointer to setting permutation vector
      IPERM = 1
      DO 160 K=1,N
        IF (I.GT.0) GO TO 60
C Pick up next root.
C Stop if all roots used (needed because of subordinate variables)
        IF (NR.GT.N) GO TO 161
        I = NE(NR)
        NE(NR) = 0
        NR = NR + 1
        IL = N
        NA(N) = 0
C Go to son for as long as possible, clearing father-son pointers
C     in IPS as each is used and setting NA(IL)=0 for all levels
C     reached.
   60   CONTINUE
        DO 70 L=1,N
          IF (IPS(I).GE.0) GO TO 80
          ISON = -IPS(I)
          IPS(I) = 0
          I = ISON
          IL = IL - 1
          NA(IL) = 0
   70   CONTINUE
   80   CONTINUE
C?? Do we want to expand for subordinate variables
C Record position of node I in the order.
        IPS(I) = K
C Add number of subordinate variables to variable I
        NE(IS) = NE(IS) + NODE(I) + 1
        IF (IL.LT.N) NA(IL+1) = NA(IL+1) + 1
        NA(IS) = NA(IL)
        ND(IS) = NV(I)
        NODE(I) = IS
        PERM(I) = IPERM
        IPERM = IPERM + 1
C Order subordinate variables to node I
        IN = I
  777   IF (SUBORD(IN).EQ.0) GO TO 778
          IN = SUBORD(IN)
          NODE(IN) = IS
          PERM(IN) = IPERM
          IPERM = IPERM + 1
          GO TO 777
C Check for static condensation
  778   IF (NA(IS).NE.1) GO TO 90
        IF (ND(IS-1)-NE(IS-1).EQ.ND(IS)) GO TO 100
C Check for small numbers of eliminations in both last two steps.
   90   IF (NE(IS).GE.NEMIN) GO TO 110
        IF (NA(IS).EQ.0) GO TO 110
        IF (NE(IS-1).GE.NEMIN) GO TO 110

C Combine the last two steps
  100   NA(IS-1) = NA(IS-1) + NA(IS) - 1
        ND(IS-1) = ND(IS) + NE(IS-1)
        NE(IS-1) = NE(IS) + NE(IS-1)
        NE(IS) = 0
        NODE(I) = IS-1
C Find eldest son (IFSON) of node I (IS)
C Note that node I must have a son (node IS-1 is youngest)
        IFSON = -FILS(I)
C Now find youngest son INO (he is node IS-1)
        IN = IFSON
 102    INO = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 102
C Cannot be root node .. so points to father
C Merge node IS-1 (INO) into node IS (I)
        NV(INO) = 0
C IPE already set .. was father pointer now principal variable pointer
C Now make subsidiary nodes of INO into subsidiary nodes of I.
C Subordinate nodes of INO become subordinate nodes of I
        IN = I
  888   IF (SUBORD(IN).EQ.0) GO TO 889
        IN = SUBORD(IN)
        NODE(IN) = IS-1
        GO TO 888
  889   SUBORD(IN) = INO
        IN = INO
        IF (SUBORD(IN).EQ.0) GO TO 887
        IN = SUBORD(IN)
        IPE(IN) = -I
  887   CONTINUE

C INOS is eldest son of INO
      INOS = -FILS(INO)

C Find elder brother of node INO
C First check to see if he is only son
      IF (IFSON.EQ.INO) GO TO 107
      IN = IFSON
 105  INS = IN
      IN =  FRERE(IN)
      IF (IN.NE.INO) GO TO 105
C INS is older brother .. make him brother of first son of INO (ie INOS)
C and  make INOS point to I now as father.
C Jump if there is no son of INO
        IF (INOS.EQ.0) THEN
C Elder brother of INO just points to (new) father.
          FRERE(INS) = -I
          GO TO 120
        ELSE
          FRERE(INS) =  INOS
        ENDIF
C INT is youngest brother of INOS.  Make him point to (new) father.
 107    IN = INOS
        IF (IN.EQ.0) GO TO 120
 108    INT = IN
        IN =  FRERE(IN)
        IF (IN.GT.0) GO TO 108
        FRERE(INT) = -I
        GO TO 120
  110   IS = IS + 1
  120   IB = IPE(I)
        IF (IB.GE.0) THEN
C Node I has a younger brother or is a root
          IF (IB.GT.0) NA(IL) = 0
          I = IB
          GO TO 160
        ELSE
C I has no brothers. Go to father of node I
          I = -IB
          IL = IL + 1
        ENDIF
  160 CONTINUE
  161 NSTEPS = IS - 1
      RETURN
      END
      SUBROUTINE MA57MD(N,NE,IRN,JCN,MAP,IRNPRM,
     +                  LROW,PERM,COUNT,IDIAG)
C
C This subroutine is called by MA57A/AD and generates the map that
C     reorders the user's input so that the upper triangle of the
C     permuted matrix is held by rows. No check is made for duplicates.
C
      INTEGER N,NE
C IRNPRM(N+NE) has this dimension to include possibility of no
C    diagonals in input.
      INTEGER IRN(NE),JCN(NE),MAP(NE),IRNPRM(N+NE),LROW(N),PERM(N),
     +        COUNT(N),
     +        IDIAG(N)
C N must be set to the matrix order. It is not altered.
C NE must be set to the number of entries input. It is not altered.
C IRN(K) and JCN(K), K=1,2,...,NE must be set to the row and column
C     numbers of the entries. If entry (IRN(K),JCN(K)) lies in the
C     lower triangular part of the permuted matrix, the values of
C     IRN(K) and JCN(K) are interchanged. Otherwise, these arrays are
C     not changed.
C MAP need not be set on input and on return holds the positions of
C     entries when the permuted upper triangle is ordered by rows.
C LROW need not be set. On return, LROW(I),I=1,N holds the number of
C     entries in row I of the permuted matrix.
C PERM(I) must be set to the position of variable I in the
C     pivot order, I=1,2,...,N.
C COUNT is used for workspace. It is set to row counts and then
C     accumulated row counts.
C IDIAG is used for workspace. It is used as pointer to diagonal entry
C     that is first in the row.
C
C Local variables
      INTEGER EXPNE,I,J,K
C I Row index
C J Column index
C K Position in IRN or JCN.

C Accumulate row counts in COUNT, and interchange row and column
C     numbers where necessary.
      DO 10 I = 1,N
C Set to 1 since diagonal will always be present as first entry.
        COUNT(I) = 1
C IDIAG used first as flag to identify duplicate diagonals
        IDIAG(I) = 0
   10 CONTINUE

C EXPNE counts number of entries plus missing diagonals
      EXPNE = NE + N
      DO 20 K = 1,NE

        I = IRN(K)
        J = JCN(K)

        IF (MAX(I,J).GT.N .OR. MIN(I,J).LT.1) THEN
          EXPNE = EXPNE - 1
          GO TO 20
        ENDIF

C Check for duplicate diagonals
        IF (I.EQ.J) THEN
          I = PERM(I)
          IF (IDIAG(I).GE.1) THEN
            COUNT(I) = COUNT(I) + 1
            IDIAG(I) = IDIAG(I) + 1
          ELSE
            IDIAG(I) = 1
            EXPNE = EXPNE - 1
          ENDIF
          GO TO 20
        ENDIF

        IF (PERM(I).LT.PERM(J)) THEN
          I = PERM(I)
          COUNT(I) = COUNT(I) + 1
        ELSE
          J = PERM(J)
          COUNT(J) = COUNT(J) + 1
        END IF

   20 CONTINUE
C
C Store row counts in LROW and accumulate row counts in COUNT.
      LROW(1) = COUNT(1)
      IDIAG(1) = MAX(IDIAG(1),1)
      DO 30 I = 2,N
        LROW(I) = COUNT(I)
C COUNT(I) set to point to position of last entry in row I of permuted
C       upper triangle.
        COUNT(I) = COUNT(I-1) + LROW(I)
        IDIAG(I) = COUNT(I-1) + MAX(IDIAG(I),1)
   30 CONTINUE

C Set diagonal entries in IRNPRM.  This is done separately because some
C     diagonals may not be present in the users input.
      DO 35 I = 1,N
        K = PERM(I)
        IRNPRM(IDIAG(K)) = I
   35 CONTINUE
C
C Run through putting the entries in the right places. COUNT is used for
C     holding running pointers and is left holding pointers to row
C     starts.
C Count number of entries in expanded matrix (allowing for non-input
C     diagonals)
      DO 40 K = 1,NE
        I = IRN(K)
        J = JCN(K)
        IF (MIN(I,J).LT.1 .OR. MAX(I,J).GT.N) THEN
          MAP(K) = 0
          GO TO 40
        ENDIF
        I = PERM(IRN(K))
        J = PERM(JCN(K))
        IF (I.EQ.J) THEN
          MAP(K) = IDIAG(I)
          IRNPRM(IDIAG(I)) = IRN(K)
          IDIAG(I) = IDIAG(I) - 1
        ELSE
          IF (I.GT.J) THEN
            MAP(K) = COUNT(J)
            IRNPRM(COUNT(J)) = IRN(K)
            COUNT(J) = COUNT(J) - 1
          ELSE
            MAP(K) = COUNT(I)
            IRNPRM(COUNT(I)) = JCN(K)
            COUNT(I) = COUNT(I) - 1
          ENDIF
C         II = MIN(I,J)
C         MAP(K) = COUNT(II)
C         COUNT(I) = COUNT(II) - 1
        ENDIF
   40 CONTINUE

C Set number of entries in expanded matrix
      IDIAG(1) = EXPNE
      RETURN

      END
      SUBROUTINE MA57ND(N,LENR,NA,NE,ND,NSTEPS,LSTKI,LSTKR,
     *                  INFO,RINFO)
C
C Storage and operation count evaluation.
C
C Evaluate number of operations and space required by factorization
C     using MA57BD.  The values given are exact only if no numerical
C     pivoting is performed.
C
C N must be set to the matrix order. It is not altered.
C LENR is number of entries in upper triangle of each row of permuted
C     matrix. It includes diagonal (no duplicates) and duplicates
C     off-diagonal but excludes any out-of-range entries.
C NA,NE,ND must be set to hold, for each tree node, the number of stack
C     elements assembled, the number of eliminations and the size of
C     the assembled front matrix respectively.  They are not altered.
C NSTEPS must be set to hold the number of tree nodes. It is not
C     altered.
C LSTKI is used as a work array by MA57ND.
C LSTKR is used as a work array by MA57ND.
C
C Counts for operations and storage are accumulated in variables
C     OPS,OPSASS,NRLTOT,NIRTOT,NRLNEC,NIRNEC,NRLADU,NIRADU.
C OPS number of multiplications and additions during factorization.
C OPSASS number of multiplications and additions during assembly.
C NRLADU,NIRADU real and integer storage respectively for the
C     matrix factors.
C NRLTOT,NIRTOT real and integer storage respectively required
C     for the factorization if no compresses are allowed.
C NRLNEC,NIRNEC real and integer storage respectively required for
C     the factorization if compresses are allowed.
C MAXFRT is maximum front size

C     .. Scalar Arguments ..
      INTEGER N,NSTEPS
C     ..
C     .. Array Arguments ..
      INTEGER LENR(N),LSTKI(N),LSTKR(N),NA(NSTEPS),
     +        ND(NSTEPS),NE(NSTEPS),INFO(40)
      DOUBLE PRECISION RINFO(20)
C     ..
C     .. Local Scalars ..
      INTEGER I,IORG,ISTKI,ISTKR,ITOP,ITREE,JORG,K,
     +        LSTK,NASSR,NELIM,NFR,NSTK,NTOTPV,NZ1,NZ2
      DOUBLE PRECISION DELIM
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC MAX
C     ..
      DOUBLE PRECISION OPS,OPSASS
      INTEGER NIRADU,NIRNEC,NIRTOT,NRLADU,NRLNEC,NRLTOT,MAXFRT

C     ..
C     .. Executable Statements ..
C
C Accumulate number of nonzeros with indices in range in NZ1.
C     Duplicates on the diagonal are ignored but NZ1 includes any
C     diagonals not present on input and duplicates off the diagonal.
      NZ1 = 0
      DO 40 I = 1,N
        NZ1 = NZ1 + LENR(I)
   40 CONTINUE
      NZ2 = NZ1
C ISTKR,ISTKI Current number of stack entries in
C     real and integer storage respectively.
C OPS,OPSASS,NRLADU,NIRADU,NIRTOT,NRLTOT,NIRNEC,NRLNEC,NZ2 are defined
C     above.
C NZ2 Current number of original matrix entries not yet processed.
C NTOTPV Current total number of rows eliminated.
C ITOP Current number of elements on the stack.
      ISTKI = 0
      ISTKR = 0
      OPS = 0.0D0
      OPSASS = 0.0D0
      NRLADU = 0
C One location is needed to record the number of blocks actually used.
      NIRADU = 3
C Initialize to what is required in MA57BD (as opposed to MA57OD).
      NIRTOT = NZ1+N+5
      NRLTOT = NZ1
      NIRNEC = NZ2+N+5
      NRLNEC = NZ2
      NTOTPV = 0
      ITOP = 0
      MAXFRT = 0
C
C Each pass through this loop processes a node of the tree.
      DO 100 ITREE = 1,NSTEPS
        NELIM = NE(ITREE)
        DELIM = NELIM
        NFR = ND(ITREE)
        MAXFRT = MAX(MAXFRT,NFR)
        NSTK = NA(ITREE)
C Adjust storage counts on assembly of current frontal matrix.
        NASSR = NELIM*(NELIM+1)/2 + NFR*NFR
C Data for no compresses so use original number, NZ1.
        NRLTOT = MAX(NRLTOT,NRLADU+NASSR+ISTKR+NZ1)
C Data for compresses so use current number, NZ2.
        NRLNEC = MAX(NRLNEC,NRLADU+NASSR+ISTKR+NZ2)

C Decrease NZ2 by the number of entries in rows being eliminated at
C     this stage.
        DO 70 IORG = 1,NELIM
          JORG = NTOTPV + IORG
          OPSASS = OPSASS + LENR(JORG)
          NZ2 = NZ2 - LENR(JORG)
   70   CONTINUE

        NTOTPV = NTOTPV + NELIM

C Remove elements from the stack.  There are ITOP elements on the
C     stack with the appropriate entries in LSTKR and LSTKI giving
C     the real and integer storage respectively for each stack
C     element.
        DO 80 K = 1,NSTK
          LSTK = LSTKR(ITOP)
          ISTKR = ISTKR - LSTK
          OPSASS = OPSASS + LSTK
          LSTK = LSTKI(ITOP)
          ISTKI = ISTKI - LSTK
          ITOP = ITOP - 1
   80   CONTINUE

C Accumulate nonzeros in factors and number of operations.
        NRLADU = NRLADU + (NELIM*(NELIM+1))/2 + (NFR-NELIM)*NELIM
        NIRADU = NIRADU + 2 + NFR
        OPS = OPS + (DELIM* (12*NFR+6*NFR*NFR - (DELIM+1)*
     +        (6*NFR+6-(2*DELIM+1))))/6 + DELIM

        IF (NFR.GT.NELIM) THEN
C Stack remainder of element.
          ITOP = ITOP + 1
          LSTKR(ITOP) = ((NFR-NELIM)*(NFR-NELIM+1))/2
          LSTKI(ITOP) = NFR - NELIM + 1
          ISTKI = ISTKI + LSTKI(ITOP)
          ISTKR = ISTKR + LSTKR(ITOP)
        ENDIF

C Adjust integer counts to stack elements and allow for next front.
        IF (ITREE.EQ.NSTEPS) THEN
          NIRTOT = MAX(NIRTOT,NIRADU+ISTKI+NZ1)
          NIRNEC = MAX(NIRNEC,NIRADU+ISTKI+NZ2)
        ELSE
          NIRTOT = MAX(NIRTOT,NIRADU+(N-NTOTPV+2)+ISTKI+NZ1)
          NIRNEC = MAX(NIRNEC,NIRADU+(N-NTOTPV+2)+ISTKI+NZ2)
        ENDIF

  100 CONTINUE
C

C Set INFO and RINFO
      INFO(5)   = NRLADU
      INFO(6)   = NIRADU
      INFO(7)   = MAXFRT
      INFO(8)   = NSTEPS
      INFO(9)   = NRLTOT
      INFO(10)  = NIRTOT
      INFO(11)  = NRLNEC
      INFO(12)  = NIRNEC
      RINFO(1)  = OPSASS
      RINFO(2)  = OPS

      RETURN
      END
      SUBROUTINE MA57OD(N,NE,A,LA,IW,LIW,LROW,PERM,NSTEPS,NSTK,NODE,
     +                  DIAG,SCHNAB,PPOS,CNTL,ICNTL,INFO,RINFO,HOLD,
     +                  BIGA)
C
C Factorization subroutine
C
C This subroutine operates on the input matrix ordered into a tentative
C     pivot order by MA57BD and produces the matrices U and inv(D)
C     of the factorization A = (U trans) D U, where D is a block
C     diagonal matrix with blocks of order 1 and 2. Gaussian elimination
C     is used with pivots of order 1 and 2, chosen according to the
C     tentative pivot order unless stability considerations
C     require otherwise.
      INTEGER N,NE,LA
C SCHNAB is in fact of dimension 5, but must set to * for Fujitsu (and
C     maybe other) compilers on very small matrices.
      DOUBLE PRECISION A(LA),DIAG(N),SCHNAB(*),CNTL(5),RINFO(20),BIGA
      INTEGER LIW,IW(LIW),LROW(N),PERM(N),NSTEPS,NSTK(NSTEPS),
     +        NODE(N),PPOS(N),ICNTL(20),INFO(40),HOLD(40)
C N   must be set to the order of the matrix. It is not altered.
C NE  must be set to the number of entries in the upper triangle of
C     the permuted matrix. It is not altered.
C A   must be set so that the upper triangle of the permuted matrix
C     is held by rows in positions LA-NE+1 to LA. Explicit zero
C     diagonals are stored. Duplicate entries are permitted
C     and are summed. During the computation, active frontal matrices
C     are held in A by rows. The working front is held in full form.
C     Stacked elements are held by rows in packed form.
C     On exit, entries 1 to INFO(10) of A hold real information on the
C     factors and should be passed unchanged to MA57CD. For each block,
C     the factorized block pivot precedes the rows of the out-of-pivot
C     part held as a rectangular matrix by rows.
C          The factorized pivot has the form:
C                -1  T
C            L  D   L
C
C     where L is unit lower triangular, D is block diagonal with blocks
C     of size 1 or 2. L and the diagonal part of D is held packed by
C     columns and the off-diagonal entries of the 2*2 blocks of D are
C     held from position IW(1).
C LA  length of array A. A value for LA sufficient for the
C     tentative pivot sequence will have been provided by MA57AD
C     in variable INFO(11). LA is not altered.
C IW  must be set on input so that IW(I+LIW-NE) holds the column
C     index of the entry in A(I+LA-NE) for I = 1,..., NE.
C     On exit, entries 1 to INFO(16) hold integer information on the
C     factors and should be passed unchanged to MA57CD.
C     IW(1) will be set to one greater than the number of entries in the
C     factors.
C     IW(2) points to the end of the factorization.
C     IW(3) will be set to the number of block pivots actually used;
C     this may be different from NSTEPS since numerical considerations
C     may prevent us choosing pivots at some stages.
C     Integer information on each block pivot row follows. For each
C     block pivot row, we have:
C       * no. of columns,
C       * no. of rows,
C       * list of column indices. The column indices for a
C         2x2 pivot are flagged negative.
C     During the computation, the array is used to hold indexing
C     information on stacked elements.  IW stores the number of
C     variables then a list of variables in the element.
C LIW must be set to the length of array IW. A value sufficient for the
C     tentative pivot sequence will have been provided by MA57AD in
C     variable INFO(12). LIW is not altered.
C LROW  must be set so that LROW(I) holds the number of entries in row
C     I (in permuted order) of the incoming matrix. LROW is not altered.
C PERM  must be set so that PERM(I) holds the variable that is Ith in
C     the tentative pivot order generated by MA57AD. PERM is not
C     altered.
C NSTEPS must be set to the number of nodes in the tree from the
C     analysis. It is the length of array NSTK. Its
C     value will never exceed N. It is not altered.
C NSTK must be set so that NSTK(I) holds the number of
C     stacked elements to be assembled at tree node I.
C NODE must be unchanged since return from MA57AD. NODE(I) gives
C     the tree node at which variable I was eliminated in the
C     analysis phase. It is of length N and is not altered.
C DIAG is only accessed if ICNTL(7) is equal to 4.
C     In that case it must be set to the values of the diagonals of
C     matrix.
C SCHNAB is only accessed if ICNTL(7) is equal to 4. It is used to hold
C     parameters for the Schnabel-Eskow modification and max/min values
C     of current diagonal.  Specifically (using notation of S-E):
C     SCHNAB(1) == GAMMA
C     SCHNAB(2) == TAUBAR (in fact root of TAUBAR)
C     SCHNAB(3) == MU
C     SCHNAB(4) == Max entry on diag
C     SCHNAB(5) == Min entry on diag
C PPOS  is integer work array of dimension N. If I is a variable in
C     the current front, PPOS(I) is used to indicate its position in the
C     front. For any other uneliminated variable, PPOS(I) is set to N+1.
C CNTL must be set (perhaps by MA57ID) so that CNTL(1) holds the
C     pivot threshold and CNTL(2) holds the pivot tolerance.
C ICNTL must be set (perhaps by MA57ID).  Entries of ICNTL accessed
C     by MA57OD are:
C ICNTL(2) is output unit for warning messages.
C ICNTL(7) is used to control pivoting.  With the default value of 1,
C     1 x 1 and 2 x 2 pivots are used subject to passing a threshold
C     tolerance.  If ICNTL(7) is greater than 1 only 1 x 1 pivots will
C     be used. If ICNTL(7) equal to 2,
C     the subroutine will exit immediately a sign change or zero pivot
C     is detected.  If ICNTL(7) is equal to 3, the subroutine will
C     continue the factorization unless a zero
C     pivot is detected.  If ICNTL(7) is equal to 4, the diagonal of
C     the matrix will be modified so that all pivots are of the same.
C ICNTL(8) is used to control whether, on running out of space, the
C     subroutine exits with an error return (ICNTL(8) = 0), or
C     whether it saves some internal variables so that
C     larger arrays can be allocated and the computation restarted
C     from the point at which it failed.
C ICNTL(11) is the block size used by the Level 3 BLAS (default 32).
C RINFO(3) will be set to the number of floating-point operations
C     required for the assembly.
C RINFO(4) will be set to the number of floating-point operations
C     required for the factorization.  RINFO(5) will be set to the
C     number of extra flops needed for the use of GEMM.
C INFO(1)  holds a diagnostic flag. It need not be set on entry. A zero
C     value on exit indicates success. Possible nonzero values are
C        -3  insufficient storage for A.
C        -4  insufficient storage for IW.
C        -5  zero pivot found when ICNTL(7) = 2 or 3.
C        -6  change in sign of pivots when ICNTL(7) = 2.
C        +4  matrix is singular
C       +10  factorizations pauses because insufficient real space
C       +11  factorizations pauses because insufficient integer space
C INFO(40) is used to accumulate number of reals discarded if
C     elimination continues without restart when space is exhausted.
C

C INFO(32) is set to number of zeros in the triangle of the factors
C INFO(33) is set to number of zeros in the rectangle of the factors
C INFO(34) is set to number of zero columns in rectangle of the factors
C Needed to compute these
      INTEGER ZCOL,RPOS

C
C Constants
      DOUBLE PRECISION ZERO,HALF,ONE
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0)
C
C Local variables
      INTEGER AINPUT
      DOUBLE PRECISION AMAX,AMULT1,AMULT2
      INTEGER APOS,APOSA,APOSB,APOSBB,APOSBK,APOSC,APOSI,APOSJ,APOSM,
     +        APOS1,APOS2,APOS3,APOS4,ASTK,ATRASH,BLK
      DOUBLE PRECISION DELTA,DETPIV
      INTEGER ELT
      DOUBLE PRECISION FLOPSA,FLOPSB,FLOPSX
      INTEGER I,I1,IASS,IBEG,IELL,IEND,IEXCH,IINPUT,INTSPA,
     +        IORG,IPIV,IPOS,IROW,ISNPIV,ISTK,ISWOP,IWNFS,IWPOS,
     +        J,JAY,JA1,JCOL,JJ,JJJ,JMAX,J1,J2,K,
     +        KB,KBLK,KCT,KR,KROW,K1,K2,L,LASPIV,LDIAG,LIELL,
     +        LP,LPIV, NBSTATIC
      LOGICAL LASTBK,LTWO
      INTEGER MAXFRT
      DOUBLE PRECISION MAXPIV
      INTEGER NASS,NBLK,NBLOC,NCMPBI,NCMPBR,NEIG,NELL,NFRONT,NIRBDU
      DOUBLE PRECISION NORMJ
      INTEGER NTWO
C LSTAT is .TRUE. is we can use static pivoting
      LOGICAL SCHUR,LSTAT
      INTEGER MPIV,NPIV,NPOTPV,NRLBDU,NSC1,NST,
     +        NSTACK(2),NSTKAC(2),NTOTPV,
     +        NUMORG,OFFDAG,PHASE,PIVBLK
      DOUBLE PRECISION PIVOT
      INTEGER PIVSIZ,POSELT,POSPV1,POSPV2,PTRA,PTRIRN,RLSPA,
     +        SIZBLK,SIZC,SIZF,TRLSPA,TINSPA,TOTSTA(2),WP,ZCOUNT
      DOUBLE PRECISION RMAX,SWOP,TMAX,TOL,UU,ULOC,UTARG,STCTOL

C AINPUT is the position of the first entry of the original matrix
C     reals since the last compress.  It is reset to the current
C     part being processed when MA57PD is called.
C AMAX is used to record the largest entry in a row.
C AMULT1, AMULT2 are used to hold multipliers.
C     Also used as temporary variables.
C APOS is a pointer to the start of the current front in A.
C APOSA holds the index of the current entry of the matrix A used to
C     form the Schur complement as C = A*B
C APOSB holds the index of the current entry of the matrix B used to
C     form the Schur complement as C = A*B
C APOSBB is a pointer to the start of the buffer before current front
C     in A. It is different from APOS because of the way we store the
C     factors and prevents us overwriting when generating them.
C APOSC holds the index of the current entry of the matrix C used to
C     form the Schur complement as C = A*B
C APOSI is set to the beginning of a row in the working front or just
C    ahead of this position.
C APOSJ is set to the beginning of row JMAX in working front.
C APOSM holds the index of the current entry of the multiplier matrix
C     for the matrix B used to form the Schur complement as C = A*B for
C     the part of the pivot rows outside the pivot block.
C APOS1 is a position in the array A.
C APOS2 is a position in the array A.
C APOS3 is the position in A of the first entry in the copy of U.
C APOS4 is the position in A of the first entry in the Schur complement.
C ASTK indicates position immediately before first stack entry.
C ATRASH is used as limit on space in A being set to zero.
C BLK is DO loop variable for blocks in pivot block.
C DELTA is the amount added to diagonal when in Phase 2 of matrix
C     modification.
C DETPIV is the value of the determinant of the 2x2 pivot or
C     candidate pivot.
C ELT is a DO loop variable indicating the element being processed.
C FLOPSA  counts floating-point operations for assembly
C FLOPSB  counts floating-point operations for elimination
C FLOPSX  counts extra flops required by use of GEMM
C I is a DO loop variable.
C I1 is used to hold limit of DO loop.
C IASS is the index of the current tree node.
C IBEG is the position of the beginning of a row in A, used when
C     performing elimination operations on the front.
C IELL Current element being assembled starts in position IELL of IW.
C IEND is the position of the end of a row in A, used when
C     performing elimination operations on the front.
C IEXCH is used to hold the contents of an array entry when doing a
C     swop.
C IINPUT is the position of the first entry of the original matrix
C     integers since the last compress.  It is reset to the current
C     part being processed when MA57PD is called. with REAL = .FALSE.
C INTSPA is amount of integer workspace needed to progress to current
C     point in elimination.
C IORG is a DO loop variable indicating which current row from the
C     incoming matrix is being processed.
C IPIV is the current relative position of the pivot search.
C IPOS is used to hold (temporarily) a position in IW.
C IROW is a DO loop variable used when scanning a row of A.
C ISNPIV is +1 if first pivot is positive and -1 if it is negative.
C ISTK points to the bottom of the stack in IW (needed by compress).
C    It indicates position immediately before first stack entry.
C ISWOP is used when swopping two integers.
C IWNFS points to the first free location for a variable that is not
C     fully summed.
C IWPOS points to the first free position for factors in IW.
C J is a temporary variable.
C JA1 is a temporary index.
C JCOL is used as an index into array A.
C JJ and JJJ are Do loop indices.
C JMAX is the relative column index in the front of the largest
C     off-diagonal in the fully summed part of the prospective pivot
C     row.
C J1 and J2 are pointers to the beginning and end of a row segment
C     in the array IW.  They are also used as running pointers in IW.
C K is a temporary variable.
C KBLK Blocked GEMM is performed on a block KBLK by KBLK matrix.
C KCT counts the number of unchecked candidates in a pivot sweep.
C KR is pointer to the current row in the assembled block being
C     tested for a potential pivot.
C KROW is a DO loop variable.
C K1 and K2 are used as running indices for array A.
C KB is block row index for blocked GEMM.
C L is a temporary variable.
C LASPIV is set to value of NPIV at end of previous block of pivots.
C LIELL is the order of the size of the reduced matrix from the front.
C     This is the order of the stacked matrix.
C LPIV is the number of pivots selected in a block pivot.
C LASTBK is flag to indicate when we are processing last block of
C       pivot block.
C LTWO  is logical variable used to indicate if current pivot is a 2 x 2
C       pivot.
C MAXFRT is the maximum front size encountered so far.
C MAXPIV is the largest of two diagonal entries of a 2 x 2 pivot.
C NASS holds the number of fully assembled variables in
C     the newly created element.
C NBLK is the number of block pivots used.
C NBLOC Number of rows of the Schur complement calculated by each GEMM
C     reference. Set to ICNTL(11).
C NCMPBR, NCMPBI are the number of compresses on real and integer space
C     respectively.
C NEIG is number of negative eigenvalues detected.
C NELL is used to hold the current number of son elements.
C NFRONT is the number of variables in the front.
C NIRBDU is number of integer entries in factors.
C NORMJ is used in matrix modification for 1-norm of off-diagonals.
C NTWO is the number of two by two full pivots used.
C SCHUR if set to .TRUE. then the Schur complement will be
C      generated using Level 3 BLAS.
C MPIV is the number of pivots so far chosen in the current block.
C NPIV is the number of pivots so far chosen at the current node.
C NPOTPV is the total number of potential pivots so far.  Variables
C     PERM(1), PERM(2), .... PERM(NPOTPV) are fully assembled.
C NRLBDU is number of real entries in factors.
C NST temporary to hold NSC1+1 if NSC1 > 0, 0 otherwise.
C NSTACK(I), I = 1,2 hold the number of active entries on the
C     real/integer stack.
C NSTKAC(I), I =1,2 hold the number of entries on the real/integer
C     stack and original matrix after a compress.
C NTOTPV is the total number of pivots selected. This is used
C     to determine whether the matrix is singular.
C NUMORG is the number of variables in the tentative pivot from the
C     incoming original rows.
C OFFDAG is the position in A of the off-diagonal entry of a 2 x 2
C        pivot.
C PIVBLK Number of rows of each block when using GEMM in pivot block.
C     Set to minimum of NBLOC and NASS.
C PIVOT is a temporary variable used to hold the value of the current
C        pivot.
C PIVSIZ is order of current pivot (has value 1 or 2).
C POSELT is a pointer into the current element being assembled.
C POSPV1 is the position in A of a 1 x 1 pivot or the first diagonal
C     of a 2 x 2 pivot.
C POSPV2 is the position in A of the second diagonal of a 2 x 2 pivot.
C PTRA  points to the next original row in A.
C PTRIRN points to the next original row in IW.
C RLSPA is amount of real workspace needed to progress to current
C     point in elimination.
C SIZBLK Number of rows in current block when using GEMM in pivot block.
C     Set to minimum of PIVBLK and NASS - NPIV
C SIZC is set to number of rows in remainder of pivot block after
C     blocking is done.
C SIZF is set to number of rows in remainder of pivot row (to NFRONT)
C     after blocking is done.
C TOTSTA(I), I =1,2 hold the number of entries on the stack and
C     original matrix.
C TRLSPA is amount of real workspace needed to progress to current
C     point in elimination if no compresses on the data are performed.
C RMAX is used to record the largest entry in a row.
C SWOP is used when swopping two reals.
C TMAX is used to record the largest entry in a row.
C TOL is the tolerance against which singularity is judged.
C     If static pivoting is used, then TOL is the value for this.
C UU is a local variable used to hold threshold parameter.  Its value is
C     between 0 and 0.5.
C ZCOUNT is number of "zero" rows in current front (used when
C     ICNTL(16) is equal to 1).
C
C Procedures
C MA57PD compresses arrays.
C MA57WD  adjusts signs in factors, moves the off-diagonal entries of
C      full 2x2 pivots and updates counts.
      DOUBLE PRECISION FD15AD

C?? To identify bug
C     LOGICAL LCASE
C     COMMON /CCASE/LCASE

      INTRINSIC MIN,MAX,ABS
      EXTERNAL DGEMM,FD15AD,MA57PD,MA57WD

C
C Initialization.
      NBLOC = ICNTL(11)
      TOL = CNTL(2)
      LP = ICNTL(1)
      WP = ICNTL(2)
      LDIAG = ICNTL(5)
      INFO(40) = 0
C A local variable UU is used for the threshold parameter, so that
C     CNTL(1) will remain unaltered.
      UU = MIN(CNTL(1),HALF)
      UU = MAX(UU,ZERO)

      LSTAT = .FALSE.
C Check if static pivoting option is on
      IF (CNTL(4).GT.ZERO) THEN
C LSTAT is not now set until number of delayed pivots is CNTL(5)*N
        IF (CNTL(5).EQ.ZERO) LSTAT = .TRUE.
        UTARG = SQRT(UU/CNTL(4))*CNTL(4)
        STCTOL = BIGA*CNTL(4)
C       TOL = STCTOL
      ENDIF

C Action if we are returning in the middle of the factorization
      IF (HOLD(1).GT.0) THEN
        INFO(1) = 0
        NBLK = HOLD(2)
        NTWO = HOLD(3)
        INFO(23) = HOLD(4)
        NCMPBR = 0
        NCMPBI = 0
        NEIG   = HOLD(6)
        MAXFRT = HOLD(7)
C Test compiler by commenting this out
        IWPOS  = HOLD(8)
        APOS   = HOLD(9)
        APOSBB = HOLD(10)
        NSTKAC(1) = HOLD(11)
        NSTKAC(2) = HOLD(12)
        AINPUT  = HOLD(13)
        IINPUT  = HOLD(14)
        ISTK    = HOLD(15)
        ASTK    = HOLD(16)
        INTSPA  = HOLD(17)
        RLSPA   = HOLD(18)
        PTRIRN  = HOLD(19)
        PTRA    = HOLD(20)
        NTOTPV  = HOLD(21)
        NPOTPV  = HOLD(22)
        NUMORG  = HOLD(23)
        NFRONT  = HOLD(24)
        NASS    = HOLD(25)
C       NCOL    = HOLD(26)
        IF (HOLD(1).EQ.1) NELL    = HOLD(27)
        IF (HOLD(1).EQ.2) NPIV    = HOLD(27)
        IASS    = HOLD(28)
        TINSPA  = HOLD(29)
        TRLSPA  = HOLD(30)
        TOTSTA(1) = HOLD(31)
        TOTSTA(2) = HOLD(32)
        NSTACK(1) = HOLD(33)
        NSTACK(2) = HOLD(34)
        INFO(32)  = HOLD(37)
        INFO(33)  = HOLD(38)
        INFO(34)  = HOLD(39)
        NBSTATIC  = HOLD(40)
        IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) ISNPIV = HOLD(35)
        IF (ICNTL(7).EQ.4) PHASE = HOLD(36)
        IF (HOLD(1).EQ.2) NSC1    = NFRONT-NPIV
        FLOPSA = RINFO(3)
        FLOPSB = RINFO(4)
        FLOPSX = RINFO(5)
        IF (HOLD(1).EQ.1) THEN
C Real arrays expanded
          HOLD(1) = 0
          GO TO 333
        ELSE
          IF (HOLD(1).EQ.3) THEN
C We ran out of space when allocating values for zero pivots at the end
C         of the factorization.  This is most likely to happen when we
C         are running with ICNTL(16) equal to 1 (dropping small entries
C         from front).
            HOLD(1) = 0
            GO TO 555
          ELSE
C Integer arrays expanded
            HOLD(1) = 0
            GO TO 444
          ENDIF
        ENDIF
      ENDIF

C NBSTATIC is the number of modified diagonal entries
      NBSTATIC = 0
C NBLK is the number of block pivots used.
      NBLK = 0
C NTWO is the number of 2 x 2 pivots used.
      NTWO = 0
C NCMPBR, NCMPBI are the number of compresses on real and integer space
C     respectively.
      NCMPBR = 0
      NCMPBI = 0
C FLOPSA is the number of floating-point operations for assembly.
      FLOPSA = ZERO
C FLOPSB is the number of floating-point operations for elimination.
      FLOPSB = ZERO
C FLOPSX  counts extra flops required by use of GEMM
      FLOPSX = ZERO
C NEIG is number of negative eigenvalues detected.
      NEIG = 0
C MAXFRT is the maximum front size encountered so far.
      MAXFRT  = 0
C All relevant INFO and RINFO parameters initialized to zero
C     so they have a valid entry on any error return.
      INFO(1) = 0
      INFO(2) = 0
      INFO(14:29) = 0
      INFO(31:35) = 0
      RINFO(3:5) = ZERO
      RINFO(14:15) = ZERO

C Initialization of array indicating positions of variables in front
      DO 10 I = 1,N
        PPOS(I) = N + 1
   10 CONTINUE
C IWPOS is set to position for first index of first block pivot
      IWPOS = 6
C Set first five entries to dummies to avoid unassigned var in MA57ED
      IW(1) = 0
      IW(2) = 0
      IW(3) = 0
      IW(4) = 0
      IW(5) = 0
C APOSBB is a pointer to the next position for storing factors in A
      APOSBB = 1
C Initialize NSTKAC and INTSPA and RLSPA
      NSTACK(1) = 0
      NSTACK(2) = 0
      NSTKAC(1) = NE
      NSTKAC(2) = NE
      TOTSTA(1) = NE
      TOTSTA(2) = NE
      INTSPA = NE+5+N
      RLSPA = NE
      TINSPA = NE+5+N
      TRLSPA = NE
C PTRIRN points to the next original row in IW.
      PTRIRN = LIW - NE + 1
C PTRA  points to the next original row in A.
      PTRA = LA - NE + 1
C ISTK points to the position in IW immediately before the stack.
      ISTK = PTRIRN - 1
C ASTK points to the position in A immediately before the stack.
      ASTK = PTRA - 1
C AINPUT is the position of the first entry of the original matrix
C     reals since the last compress.  It is reset to the current
C     part being processed when MA57PD is called with REAL .TRUE.
      AINPUT = PTRA
C IINPUT is the position of the first entry of the original matrix
C     integers since the last compress.  It is reset to the current
C     part being processed when MA57PD is called. with REAL = .FALSE.
      IINPUT = PTRIRN
C NTOTPV is the total number of pivots selected.
      NTOTPV = 0
C NPOTPV is the total number of potential pivots so far.
      NPOTPV = 0
C In case we run out of space before the first pivot is chosen, we
C     must initialize ISNPIV.
      IF (ICNTL(7).EQ.2 .OR. ICNTL(7).EQ.3) ISNPIV = 0


C Calculate diagonal of matrix and store in DIAG
      IF (ICNTL(7).EQ.4) THEN
        PHASE = 1
        DO 19 I = 1,N
          DIAG(I) = ZERO
   19   CONTINUE
        APOS1 = PTRA-1
        J1 = PTRIRN
        DO 20 I = 1,N
          J2 = J1 + LROW(I) - 1
          DO 25 JJ = J1,J2
            J = IW(JJ)
            APOS1 = APOS1 + 1
            IF (J.EQ.PERM(I)) DIAG(J) = DIAG(J) + A(APOS1)
   25     CONTINUE
          J1 = J2 + 1
   20   CONTINUE
        SCHNAB(1) = ONE
        SCHNAB(5) = ZERO
        DO 21 I = 1,N
          SCHNAB(1) = MAX(SCHNAB(1),ABS(DIAG(I)))
          SCHNAB(5) = MIN(SCHNAB(5),DIAG(I))
   21   CONTINUE
C Set max entry on diag
        SCHNAB(4) = SCHNAB(1)
C S+E has **2/3 .. Nick wants **1/3
        SCHNAB(2) = FD15AD('E')**(1.0/3.0)
        SCHNAB(3) = 0.1
C Initialize RINFO(15) to compute smallest pivot in modified matrix
        RINFO(15) = FD15AD('H')
        DELTA     = ZERO
      ENDIF

C   *****************************************************************
C   * Each pass through this main loop performs all the operations  *
C   * associated with one node of the assembly tree.                *
C   *****************************************************************

      IASS = 1
C     DO 2160 IASS = 1,NSTEPS
 2160 CONTINUE
C Find the frontal variables, ordered with the fully summed variables
C     of the incoming rows first, the fully summed rows from previous
C     steps next, followed by the rest in any order.

C NUMORG is the number of variables in the tentative pivot from the
C     incoming original rows.
C Calculate NUMORG and put indices of these fully summed rows in IW.
        NUMORG = 0
        DO 30 I = NPOTPV + 1,N
C J is Ith variable in tentative pivotal sequence.
          J = PERM(I)
C Jump if we have finished with variables in current node.
          IF (ABS(NODE(J)).GT.IASS) GO TO 40
          IW(IWPOS+NUMORG) = J
          NUMORG = NUMORG + 1
          PPOS(J) = NUMORG
   30   CONTINUE

C NASS will be set to the total number of fully assembled variables in
C     the newly created element. First set it to NUMORG.
   40   NASS = NUMORG
C Add indices of fully summed variables of stacked sons to IW.
        NELL = NSTK(IASS)
        IELL = ISTK + 1
        DO 70 ELT = 1,NELL
          DO 50 JJ = IELL + 1,IELL + IW(IELL)
            J = IW(JJ)
            IF (NODE(J).GT.IASS) GO TO 50
C Jump if variable already included.
            IF (PPOS(J).LE.N) GO TO 50
            IW(IWPOS+NASS) = J
            NASS = NASS + 1
            PPOS(J) = NASS
   50     CONTINUE
          IELL = IELL + IW(IELL) + 1
   70   CONTINUE
C IWNFS points to the first free location for a variable that is not
C     fully summed.
        IWNFS = IWPOS + NASS

C Incorporate original rows.
C J1 is the position of the start of the first original row associated
C     with this node of the assembly tree.
        J1 = PTRIRN
        DO 90 IORG = 1,NUMORG
          J2 = J1 + LROW(NPOTPV+IORG) - 1
C Run through index list of original row.
          DO 80 JJ = J1,J2
            J = IW(JJ)
C Jump if variable already included.
            IF (PPOS(J).LE.N) GO TO 80
            IW(IWNFS) = J
            IWNFS = IWNFS + 1
            PPOS(J) = IWNFS - IWPOS
   80     CONTINUE
          J1 = J2 + 1
   90   CONTINUE

C Now incorporate stacked elements.
C J1 is set to beginning
C J2 is set to end
        IELL = ISTK + 1
        DO 170 ELT = 1,NELL
          J1 = IELL+1
          J2 = IELL+IW(IELL)
          DO 150 JJ = J1,J2
            J = IW(JJ)
C Jump if already assembled
            IF (PPOS(J).LE.N) GO TO 150
            IW(IWNFS) = J
            IWNFS = IWNFS + 1
            PPOS(J) = IWNFS - IWPOS
  150     CONTINUE
          IELL = J2 + 1
  170   CONTINUE

C NFRONT is the number of variables in the front.
        NFRONT = IWNFS - IWPOS

C MAXFRT is the largest front size so far encountered.
        MAXFRT = MAX(MAXFRT,NFRONT)

C Set APOS to the position of first entry in frontal matrix.
        IF (INFO(1).NE.-3) THEN
C Buffer space allocated so that triangular part of pivot can be stored
C   without danger of overwrite.
          APOS = APOSBB + (NASS*(NASS+1))/2
        ELSE
          APOS = 1
        END IF
C
C Assemble reals into frontal matrix.
C
C Accumulate real space needed
        RLSPA  = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+NSTKAC(1))
        TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+TOTSTA(1))

C If necessary, compress A.

  333   IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN


          CALL MA57PD(A,IW,ASTK,AINPUT,PTRA,.TRUE.)

          NCMPBR = NCMPBR + 1
          IF (APOS+NFRONT*NFRONT-1.GT.ASTK) THEN
            IF (ICNTL(8).NE.0) THEN
C Zero part of A to avoid failure in HSL_MA57
              DO 334 I = APOSBB,ASTK
                A(I) = ZERO
  334         CONTINUE
              HOLD(1) = 1
              HOLD(2) = NBLK
              HOLD(3) = NTWO
              HOLD(4) = INFO(23)
              HOLD(5) = NCMPBI
              HOLD(6) = NEIG
              HOLD(7) = MAXFRT
              HOLD(8) = IWPOS
              HOLD(9) = APOS
              HOLD(10) = APOSBB
              HOLD(11) = NSTKAC(1)
              HOLD(12) = NSTKAC(2)
              HOLD(13) = AINPUT
              HOLD(14) = IINPUT
              HOLD(15) = ISTK
              HOLD(16) = ASTK
              HOLD(17) = INTSPA
              HOLD(18) = RLSPA
              HOLD(19) = PTRIRN
              HOLD(20) = PTRA
              HOLD(21) = NTOTPV
              HOLD(22) = NPOTPV
              HOLD(23) = NUMORG
              HOLD(24) = NFRONT
              HOLD(25) = NASS
C             HOLD(26) = NCOL
              HOLD(27) = NELL
              HOLD(28) = IASS
              HOLD(29) = TINSPA
              HOLD(30) = TRLSPA
              HOLD(31) = TOTSTA(1)
              HOLD(32) = TOTSTA(2)
              HOLD(33) = NSTACK(1)
              HOLD(34) = NSTACK(2)
              IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
              IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
              HOLD(37) = INFO(32)
              HOLD(38) = INFO(33)
              HOLD(39) = INFO(34)
              RINFO(3) =FLOPSA
              RINFO(4) =FLOPSB
              RINFO(5) =FLOPSX
              HOLD(40) = NBSTATIC
              INFO(35) = HOLD(40)
              INFO(1) = 10
              RETURN
            ELSE
C INFO(40) accumulates number of discards from factors.
              INFO(40) = INFO(40) + APOS - 1
              APOS = 1
              APOSBB = 1
              INFO(1) = -3
              IF (NFRONT*NFRONT.GT.ASTK) THEN
                INFO(17) = MAX(INFO(17),RLSPA)
                IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
                INFO(2) = LA
                RETURN
              ENDIF
            ENDIF
          ENDIF
        END IF

        ATRASH = APOS + NFRONT*NFRONT - 1
C Zero out appropriate part of A for incoming potential pivot rows.
        DO 210 JJ = APOS,ATRASH
          A(JJ) = ZERO
  210   CONTINUE

C Incorporate reals from original rows.
        J1 = PTRIRN
        DO 230 IORG = 1,NUMORG
C APOSI indicates the position in A just before the beginning of the row
C       being assembled.
          J = PERM(NPOTPV+IORG)
          APOSI = APOS + (PPOS(J)-1)*NFRONT - 1
          J2 = J1 + LROW(NPOTPV+IORG) - 1
          FLOPSA = FLOPSA + J2 - J1 + 1
          DO 220 JJ = J1,J2
            JAY = IW(JJ)
CCC
C Entries always in upper triangle because of ordering.
C Pivot permutations can only affect fully summed variables
C           IF (PPOS(JAY).GE.PPOS(J)) THEN
              APOS2 = APOSI + PPOS(JAY)
C           ELSE
C             APOS2 = APOS + (PPOS(JAY)-1)*NFRONT + PPOS(J) - 1
C           ENDIF
            A(APOS2) = A(APOS2) + A(PTRA)
            PTRA = PTRA + 1
  220     CONTINUE
C         AINPUT = AINPUT + J2 - J1 + 1
          NSTKAC(1) = NSTKAC(1) - J2 + J1 - 1
          J1 = J2 + 1
  230   CONTINUE
C       IINPUT = IINPUT + J2 - PTRIRN
        NSTKAC(2) = NSTKAC(2) - J1 + PTRIRN
        PTRIRN = J1
C Update NPOTPV
        NPOTPV = NPOTPV + NUMORG

C???
C?? Depends if we need lower triangle and whether all entries are
C       already in upper triangle.
C Correct the leading NUMORG*NASS block to allow for having only
C  assembled one of each pair of off-diagonal entries.
C       DO 410 I = 1,NUMORG
C APOS2 points to the diagonal entry of row I.
C         APOS2 = APOS + (NFRONT+1)* (I-1)
C         DO 390 J = 1,NUMORG - I
C           A(APOS2+J*NFRONT) = A(APOS2+J*NFRONT) + A(APOS2+J)
C 390     CONTINUE
C         DO 400 J = 1,NASS - I
C           A(APOS2+J) = A(APOS2+J*NFRONT)
C 400     CONTINUE
C 410   CONTINUE

C Now assemble reals from stacked elements
C POSELT is a running pointer into that element.
        DO 380 ELT = 1,NELL
          POSELT = ASTK + 1
          LIELL = IW(ISTK+1)
          J1 = ISTK + 2
          J2 = ISTK+1 + LIELL
          FLOPSA = FLOPSA + (LIELL*(LIELL+1))/2
          DO 250 JJ = J1,J2
            J = IW(JJ)
            APOS2 = APOS + (PPOS(J)-1)*NFRONT
            APOS1 = POSELT
            DO 240 JJJ=JJ,J2
              JAY = IW(JJJ)
C This part has been modified  (2/11/04)
C???          APOS3 = APOS2 + PPOS(JAY) - 1
C???          A(APOS3) = A(APOS3) + A(APOS1)
C To ensure there is valid entry in upper triangle
C???          APOS5 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
C???          IF (APOS3.NE.APOS5) A(APOS5) = A(APOS5) + A(APOS1)
              IF (PPOS(JAY) .GE. PPOS(J)) THEN
                APOS3 = APOS2 + PPOS(JAY) - 1
              ELSE
                APOS3 = APOS+(PPOS(JAY)-1)*NFRONT+PPOS(J)-1
              ENDIF
              A(APOS3) = A(APOS3) + A(APOS1)
              APOS1 = APOS1 + 1
  240       CONTINUE
            POSELT = POSELT + LIELL - (JJ-J1)
  250     CONTINUE
C ISTK and ASTK updated to point to posn before next element on stack.
          NSTKAC(2) = NSTKAC(2) - (J2-ISTK)
          NSTACK(2) = NSTACK(2) - (J2-ISTK)
          TOTSTA(2) = TOTSTA(2) - (J2-ISTK)
          ISTK = J2
          ASTK = ASTK + (LIELL*(LIELL+1))/2
          NSTKAC(1) = NSTKAC(1) - (LIELL*(LIELL+1))/2
          NSTACK(1) = NSTACK(1) - (LIELL*(LIELL+1))/2
          TOTSTA(1) = TOTSTA(1) - (LIELL*(LIELL+1))/2
  380   CONTINUE

C       IF (LCASE) THEN
C         write(7,'(/A,I8/A,5I8)') '*** Frontal matrix before step',
C    *    IASS,'NFRONT,NASS,NUMORG,APOS,IWPOS',
C    *          NFRONT,NASS,NUMORG,APOS,IWPOS
C         write(7,'(/A/(10I8))') 'IW array',
C    *      (IW(IWPOS+I-1),I=1,NFRONT)
C         DO 1122 J = 1, NFRONT
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      (A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT)
C1122     CONTINUE
C       ENDIF

C ******************
C Each time round this loop, we sweep through the remainder
C      of the assembled part of the front looking for pivots.
C ******************

C Set PIVBLK
        PIVBLK = MIN(NBLOC,NASS)
C Set pointer to first entry in block
        APOSBK = APOS
C NPIV is the number of pivots so far selected at the current node.
        NPIV = 0
C Set local value for U
        ULOC = UU

C Each pass through loop processes one block.
        DO 918 BLK = 1,NASS
C Set last block flag
        IF (NPIV+PIVBLK .GE. NASS) THEN
          LASTBK = .TRUE.
          SIZBLK = NASS - NPIV
        ELSE
          LASTBK = .FALSE.
          SIZBLK = PIVBLK
        ENDIF
C Record number of rows processed to date
        LASPIV = NPIV
C MPIV is the number of pivots so far selected in the current block.
        MPIV = 0
C KR is relative position of current pivot candidate
        KR = 0
CCC Set to following to force 2 by 2 pivots in Nocedal examples
C       KR = NUMORG
C KCT is set to one more than the number of unsearched pivot candidates
        KCT = SIZBLK + 1

C Loop for pivot searches
  920   CONTINUE
C Increment pointer for circular sweep.
          KR = KR + 1
          KCT = KCT - 1
          IF (KCT.EQ.0) GO TO 930
          IF (KR.GT.SIZBLK) KR = MPIV + 1
C Check to see if no pivot was chosen is complete sweep of pivot block.
C We either take the diagonal entry or the 2 by 2 pivot with the
C     largest fully summed off-diagonal at each stage.
C Note that IPIV is the position within the complete current front.
          IPIV = LASPIV + KR
C APOSI is set to the beginning of row IPIV in working front.
            APOSI = APOS + (IPIV-1)*NFRONT
C Set position and value of potential pivot.
            POSPV1 = APOSI + IPIV - 1
            PIVOT = A(POSPV1)

   29       IF (ICNTL(7).EQ.4) THEN
             IF (PHASE.EQ.2) THEN
C In phase 2 of matrix modification (ICNTL(7) = 4)
C Compute quantity to add to diagonal
C Calculate norm of pivot row
              IF (INFO(27).EQ.0) INFO(27) = NTOTPV + 1
              NORMJ = ZERO
              DO 28 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
                NORMJ = NORMJ + ABS(A(I))
   28         CONTINUE
              DELTA = MAX(ZERO,
     *                    - A(POSPV1) + MAX(NORMJ,SCHNAB(2)*SCHNAB(1)))
              A(POSPV1) = A(POSPV1) + DELTA
              IF (A(POSPV1).EQ.ZERO) GO TO 970
              RINFO(15) = MIN(RINFO(15),A(POSPV1))
              DIAG(PERM(NTOTPV+1)) = DELTA
              PIVSIZ = 1
              GO TO 811
             ENDIF
            ENDIF
            IF (ICNTL(7).GT.1) THEN
C Action if no pivoting requested
              IF (ABS(PIVOT).LE.CNTL(2)) THEN
                IF (ICNTL(7).LT.4) GO TO 970
C We are now in phase 2 of matrix modification (ICNTL(7) = 4)
                PHASE = 2
                GO TO 29
              ENDIF
              IF (NTOTPV.EQ.0) THEN
                IF (PIVOT.GT.ZERO) ISNPIV = 1
                IF (PIVOT.LT.ZERO) ISNPIV = -1
              ELSE
                IF (ICNTL(7).EQ.2 .AND. ISNPIV*PIVOT.LT.ZERO) GO TO 980
                IF (ICNTL(7).EQ.3 .AND. ISNPIV*PIVOT.LT.ZERO) THEN
                    INFO(26) = INFO(26) + 1
                    ISNPIV = -ISNPIV
                ENDIF
              ENDIF
              IF (ICNTL(7).EQ.4) THEN
                IF (PIVOT.GE.SCHNAB(1)*SCHNAB(2) .AND.
     *              SCHNAB(5).GE.-SCHNAB(3)*SCHNAB(4)) THEN
C Update and check values of future diagonals
                  SCHNAB(5) = ZERO
                  SCHNAB(4) = ZERO
                  DO 22 I = POSPV1+1,POSPV1+NFRONT-NPIV-1
                    J = IW(IWPOS+NPIV+I-POSPV1)
                    DIAG(J) = DIAG(J) - A(I)*A(I)/PIVOT
                    SCHNAB(5) = MIN(DIAG(J),SCHNAB(5))
                    SCHNAB(4) = MAX(DIAG(J),SCHNAB(4))
                    IF (DIAG(J).LT.-SCHNAB(3)*SCHNAB(1)) THEN
                      PHASE = 2
                      GO TO 29
                    ENDIF
   22             CONTINUE
                  DIAG(PERM(NTOTPV+1)) = ZERO
                  RINFO(15) = MIN(RINFO(15),PIVOT)
                ELSE
                  PHASE = 2
                  GO TO 29
                ENDIF
              ENDIF
              PIVSIZ = 1
              GO TO 811
            ENDIF
C Find largest off-diagonal entry in the part of the row in which we
C     seek a pivot.
            AMAX = ZERO
            JMAX = 0
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
            DO 110 K = 1, IPIV - NPIV - 1
              IF (ABS(A(POSPV1-K*NFRONT)).GT.AMAX) THEN
                AMAX = ABS(A(POSPV1-K*NFRONT))
                JMAX = IPIV - K
              ENDIF
  110       CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
            DO 111 K =  1, MIN(NASS,LASPIV+PIVBLK) - IPIV
              IF (ABS(A(POSPV1+K)).GT.AMAX) THEN
                AMAX = ABS(A(POSPV1+K))
                JMAX = IPIV + K
              ENDIF
  111       CONTINUE
C Do same for the other part.
            RMAX = ZERO

C     restrict partial pivoting check to the fully summed block
C            RMAX = STCTOL
C            DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NASS
C NFRONT-IPIV
C              RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
C  112       CONTINUE

            DO 112 K = MIN(NASS,LASPIV+PIVBLK)-IPIV+1,NFRONT-IPIV
               RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
 112        CONTINUE

C Action taken if matrix is singular.
            IF (MAX(AMAX,RMAX,ABS(PIVOT)).LE.TOL) THEN
C Skip if all of row is zero.
              GO TO 920
            END IF
C Jump if no nonzero entry in row of pivot block
            IF (MAX(AMAX,ABS(PIVOT)).LE.TOL) GO TO 920
            PIVSIZ = 0
            IF (ABS(PIVOT).GT.ULOC*MAX(RMAX,AMAX)) THEN
              PIVSIZ = 1
              A(POSPV1) = PIVOT
C 1 x 1 pivot is chosen
              GO TO 810

            END IF
C If there is only one remaining fully summed row and column exit.
            IF (NPIV+1.EQ.NASS) THEN
              A(POSPV1) = PIVOT
              GO TO 920
            END IF

C Jump if 2 x 2 candidate is diagonal
            IF (AMAX.LE.TOL) GO TO 920

C      Check block pivot of order 2 for stability.
C      Find largest entry in row IPIV outwith the pivot.
            IF (RMAX.LT.AMAX) THEN
              RMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
              DO 113 K = 1, IPIV - NPIV - 1
                IF (IPIV-K.EQ.JMAX) GO TO 113
                RMAX=MAX(RMAX,ABS(A(POSPV1-K*NFRONT)))
  113         CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
              DO 114 K =  1, NFRONT - IPIV
                IF (IPIV+K.EQ.JMAX) GO TO 114
                RMAX = MAX(RMAX,ABS(A(POSPV1+K)))
  114         CONTINUE
            ENDIF

C APOSJ is set to the beginning of row JMAX in working front.
            APOSJ = APOS + (JMAX-1)*NFRONT
C POSPV2 is the position in A of the second diagonal of a 2 x 2 pivot.
C OFFDAG is the position in A of the off-diagonal of a 2 x 2 pivot.
            POSPV2 = APOSJ + JMAX - 1
            IF (IPIV.GT.JMAX) THEN
              OFFDAG = APOSJ + IPIV - 1
            ELSE
              OFFDAG = APOSI + JMAX - 1
            END IF

C       Find largest entry in row JMAX outwith the pivot.
            TMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
            DO 115 K = 1, JMAX - NPIV - 1
              IF (JMAX-K.EQ.IPIV) GO TO 115
              TMAX=MAX(TMAX,ABS(A(POSPV2-K*NFRONT)))
  115       CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
            DO 116 K =  1, NFRONT - JMAX
              IF (JMAX+K.EQ.IPIV) GO TO 116
              TMAX = MAX(TMAX,ABS(A(POSPV2+K)))
  116       CONTINUE


C DETPIV is the value of the determinant of the 2x2 pivot.
            DETPIV = A(POSPV1)*A(POSPV2) - AMAX*AMAX
            MAXPIV = MAX(ABS(A(POSPV1)),ABS(A(POSPV2)))
            IF (MAXPIV.EQ.ZERO) MAXPIV = ONE
            IF (ABS(DETPIV)/MAXPIV.LE.TOL) GO TO 920
            PIVSIZ = 2
C Check pivot for stability
C Jump if pivot fails test
C This is componentwise test
            IF ((ABS(A(POSPV2))*RMAX+AMAX*TMAX)*ULOC.GT.
     +          ABS(DETPIV)) GO TO 920
            IF ((ABS(A(POSPV1))*TMAX+AMAX*RMAX)*ULOC.GT.
     +          ABS(DETPIV)) GO TO 920
C 2 x 2 pivot is chosen

C
C           Pivot has been chosen. It has order PIVSIZ.
  810       LPIV = IPIV
            IF (PIVSIZ.EQ.2) LPIV = MIN(IPIV,JMAX)
C Change made at Stephane's suggestion
CCC         KR = MAX(KR,NPIV+PIVSIZ)
            KR = MAX(KR,MPIV+PIVSIZ)
            KCT = SIZBLK - MPIV - PIVSIZ + 1

C The following loop moves the pivot block to the top left
C           hand corner of the uneliminated frontal matrix.
            DO 860 KROW = NPIV,NPIV + PIVSIZ - 1
C We jump if swop is not necessary.
              IF (LPIV.EQ.KROW+1) GO TO 850

C Swop first part of rows (going down columns)
C JA1 is used as running index for row LPIV
              JA1 = APOS + (LPIV-1)
C J1 is used as running index for row KROW+1
              J1 = APOS + KROW
              DO 820 JJ = 1,KROW
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
                JA1 = JA1 + NFRONT
                J1 = J1 + NFRONT
  820         CONTINUE
C Swop middle part of rows (KROW+1 by rows, LPIV by columns)
              JA1 = JA1 + NFRONT
              J1 = J1 + 1
              DO 830 JJ = 1,LPIV - KROW - 2
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
                JA1 = JA1 + NFRONT
                J1 = J1 + 1
  830         CONTINUE
C Swop diagonals
              SWOP = A(APOS+KROW* (NFRONT+1))
              A(APOS+KROW* (NFRONT+1)) = A(JA1)
              A(JA1) = SWOP
C Swop last part of rows
              DO 840 JJ = 1,NFRONT - LPIV
                JA1 = JA1 + 1
                J1 = J1 + 1
                SWOP = A(JA1)
                A(JA1) = A(J1)
                A(J1) = SWOP
  840         CONTINUE
C Swop integer indexing information
              IPOS = IWPOS + KROW
              IEXCH = IWPOS + LPIV - 1
              ISWOP = IW(IPOS)
              IW(IPOS) = IW(IEXCH)
              IW(IEXCH) = ISWOP
C Set LPIV for the swop of the second row of block pivot.
  850         LPIV = MAX(IPIV,JMAX)
  860       CONTINUE

C
C Set POSPV1 and POSPV2 to new position of pivots.
  811       POSPV1 = APOS + NPIV* (NFRONT+1)
            POSPV2 = POSPV1 + NFRONT + 1
            IF (PIVSIZ.EQ.1) THEN
C Perform the elimination using entry A(POSPV1) as pivot.
C We store U and D inverse.
C Later we store D inverse U which is passed to the solution entry.
              FLOPSB = FLOPSB + ONE
              A(POSPV1) = ONE/A(POSPV1)
              IF (A(POSPV1).LT.ZERO) NEIG = NEIG + 1
              J1 = POSPV1 + 1
              J2 = POSPV1 + NASS - (NPIV+1)
              IBEG = POSPV1 + NFRONT + 1
              IEND = APOS + (NPIV+1)*NFRONT + NFRONT - 1
              DO 880 JJ = J1,J2
C AMULT1 is used to hold the multiplier
                AMULT1 = -A(JJ)*A(POSPV1)
C Hold original entry for GEMM multiply
                IF (.NOT.LASTBK) A(POSPV1+(JJ-J1+1)*NFRONT) = A(JJ)
                JCOL = JJ
                FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
                IF (MPIV+JJ-J1+2.GT.PIVBLK) GO TO 871
C                The following special comment forces vectorization on
C                   Crays.
CDIR$            IVDEP
                DO 870 IROW = IBEG,IEND
                  A(IROW) = A(IROW) + AMULT1*A(JCOL)
                  JCOL = JCOL + 1
  870           CONTINUE
  871           A(JJ) = AMULT1
                IBEG = IBEG + NFRONT + 1
                IEND = IEND + NFRONT
  880         CONTINUE
              NPIV = NPIV + 1
              MPIV = MPIV + 1
              NTOTPV = NTOTPV + 1

C       IF (LCASE) THEN
C         write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ
C         write(7,'(/A/(10I8))') 'IW array',
C    *      (IW(IWPOS+I-1),I=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      ((A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT),J=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      (A(APOS+I-1),I=1,NFRONT*NFRONT)
C       ENDIF

              IF (MPIV.EQ.SIZBLK) GO TO 930
            ELSE
C Perform elimination using block pivot of order two.
C Replace block pivot by its inverse.
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6.0

              SWOP = A(POSPV2)
              IF (DETPIV.LT.ZERO) THEN
                NEIG = NEIG + 1
              ELSE
                IF (SWOP.LT.ZERO) NEIG = NEIG + 2
              END IF

              A(POSPV2) = A(POSPV1)/DETPIV
              A(POSPV1) = SWOP/DETPIV
              A(OFFDAG) = -A(OFFDAG)/DETPIV

              J1 = POSPV1 + 2
              J2 = POSPV1 + NASS - (NPIV+1)
C             J2 = POSPV1 + NFRONT - (NPIV+1)
              IBEG = POSPV2 + NFRONT + 1
              IEND = APOS + (NPIV+2)*NFRONT + NFRONT - 1
              DO 900 JJ = J1,J2
                K1 = JJ
                K2 = JJ + NFRONT
                AMULT1 = - (A(POSPV1)*A(K1)+A(POSPV1+1)*A(K2))
                AMULT2 = - (A(POSPV1+1)*A(K1)+A(POSPV2)*A(K2))
                IF (.NOT.LASTBK) THEN
                  A(POSPV1 + (JJ-J1+2)*NFRONT) = A(K1)
                  A(POSPV1 + (JJ-J1+2)*NFRONT + 1) = A(K2)
                ENDIF
                FLOPSB = FLOPSB + (IEND-IBEG+1)*4 + 6
                IF (MPIV+JJ-J1+3.GT.PIVBLK) GO TO 891
C                The following special comment forces vectorization on
C                  Crays.
CDIR$            IVDEP
                DO 890 IROW = IBEG,IEND
                  A(IROW) = A(IROW) + AMULT1*A(K1) + AMULT2*A(K2)
                  K1 = K1 + 1
                  K2 = K2 + 1
  890           CONTINUE
  891           A(JJ) = AMULT1
                A(JJ+NFRONT) = AMULT2
                IBEG = IBEG + NFRONT + 1
                IEND = IEND + NFRONT
  900         CONTINUE
C Flag column indices of 2 x 2 pivot.
              IPOS = IWPOS + NPIV
              IW(IPOS) = -IW(IPOS)
              IW(IPOS+1) = -IW(IPOS+1)
              NPIV = NPIV + 2
              MPIV = MPIV + 2
              NTOTPV = NTOTPV + 2
              NTWO = NTWO + 1

C       IF (LCASE) THEN
C         write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ',
C    *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ
C         write(7,'(/A/(10I8))') 'IW array',
C    *      (IW(IWPOS+I-1),I=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      ((A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT),J=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      (A(APOS+I-1),I=1,NFRONT*NFRONT)
C       ENDIF

              IF (MPIV.EQ.SIZBLK) GO TO 930
            END IF

        GO TO 920
C 920   CONTINUE

 930    IF (LASTBK) THEN
          IF (NPIV.EQ.NASS) GO TO 935
C Jump if we do not have static pivoting switched on
          IF (.NOT. LSTAT)  GO TO 935
C First reduce local value of threshold to see if we can find pivots
          ULOC = ULOC/10.0D0
          IF (ULOC.LT.UTARG) THEN
C This stops us reducing it beyond UTARG
C At this point we stop looking for pivots in normal way and go to
C         static pivoting option.
            ULOC = ULOC * 10.0D0
            GO TO 9919
          ENDIF
          KCT = SIZBLK + 1 - MPIV
C Search for pivots using new value of ULOC
          GO TO 920
C Some old experiments
C         IF (ULOC.EQ.CNTL(4)) GO TO 9919
C         ULOC = MAX(ULOC*1.0D-2,CNTL(4))
C         IF (ABS(ULOC-CNTL(4))/CNTL(4) .LT. 10) ULOC = CNTL(4)
C         KCT = SIZBLK + 1 - MPIV
C         GO TO 920
        ENDIF

C Check if any pivots chosen from this block. If not, increase PIVBLK
C       and try again.
        IF (MPIV.EQ.0) THEN
          PIVBLK = 2*PIVBLK
          GO TO 918
        ENDIF
C Finished pivoting on block BLK ... now update rest of pivot block
C       using GEMM
        KBLK = (NASS-(LASPIV+PIVBLK))/PIVBLK
        L = NASS - (LASPIV+PIVBLK)
        APOS4 = APOS+(LASPIV+PIVBLK)*(NFRONT+1)
        DO 931 KB = 1,KBLK
          FLOPSX = FLOPSX + PIVBLK*(PIVBLK-1)*MPIV
          CALL DGEMM('N','N',L-(KB-1)*PIVBLK,PIVBLK,MPIV,ONE,
     +               A(APOSBK+PIVBLK*KB),NFRONT,
     +               A(APOSBK+PIVBLK*KB*NFRONT),NFRONT,ONE,
     +               A(APOS4+PIVBLK*(KB-1)*(NFRONT+1)),NFRONT)
C And now process the part of the pivot row outside the fs block
          IF (NFRONT.GT.NASS)
     +    CALL DGEMM('N','T',NFRONT-NASS,PIVBLK,MPIV,ONE,
     +               A(APOSBK+NASS-LASPIV),NFRONT,
     +               A(APOSBK+PIVBLK*KB),NFRONT,ONE,
     +               A(APOSBK+KB*NFRONT*PIVBLK+NASS-LASPIV),NFRONT)

  931   CONTINUE

       SIZC = NASS - (KBLK+1)*PIVBLK - LASPIV
       SIZF = NFRONT - (KBLK+1)*PIVBLK - LASPIV
       APOSA = APOSBK + (KBLK+1)*PIVBLK
       DO 934 K = 1,MPIV
         APOSB = APOSBK + NFRONT*PIVBLK*(KBLK+1) + K - 1
         APOSM = APOSBK + PIVBLK*(KBLK+1) + (K-1)*NFRONT
         APOSC = APOSBK + PIVBLK*(KBLK+1)*(NFRONT+1)
         DO 933 JJ = 1,SIZC
            DO 932 J = JJ,SIZC
              A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSB)
  932       CONTINUE
C And now process the part of the pivot row outside the fs block
            DO 936 J = SIZC+1,SIZF
              A(APOSC+J-1) = A(APOSC+J-1) + A(APOSA+J-1)*A(APOSM)
  936       CONTINUE
            APOSC = APOSC + NFRONT
            APOSB = APOSB + NFRONT
            APOSM = APOSM + 1
  933     CONTINUE
          APOSA = APOSA + NFRONT
  934   CONTINUE

        APOSBK = APOSBK + MPIV*(NFRONT+1)
        LASPIV = NPIV

C       IF (LCASE) THEN
C         write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C    *   'NFRONT,NASS,APOS,IWPOS,NPIV',
C    *    NFRONT,NASS,APOS,IWPOS,NPIV
C         write(7,'(/A,2I8)') 'After blocking .. APOSBK,LASPIV',
C    *    APOSBK,LASPIV
C         write(7,'(/A/(10I8))') 'IW array',
C    *      (IW(IWPOS+I-1),I=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      ((A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT),J=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      (A(APOS+I-1),I=1,NFRONT*NFRONT)
C       ENDIF

  918   CONTINUE
C End of main elimination loop.
CCC
C       IF (LP.GE.0) WRITE(LP,'(A)') '****** BE WORRIED LOOP 918'

C SCHUR if set to .TRUE. then the Schur complement will be generated
C      using level 3 BLAS at this step so an extra copy of pivot row is
C      made.


C *************************
C     Do static pivoting  *
C *************************
 9919      IPIV = LASPIV+MPIV
 9920      IPIV = IPIV + 1
CADD Probably not needed .. use only IPIV
C          IF (KR.GT.SIZBLK) KR = MPIV + 1
C     Note that IPIV is the position within the complete current front.
C          IPIV = LASPIV + KR
C     APOSI is set to the beginning of row IPIV in working front.
           APOSI = APOS + (IPIV-1)*NFRONT
C     Set position and value of potential pivot.
           POSPV1 = APOSI + IPIV - 1
           PIVOT = A(POSPV1)
CADD
C Although theses are not needed just now they are kept for when
C we use 2 x 2 static pivots.
CCC        PIVSIZ = 1
CCC        LPIV = IPIV

C Logic has changed so no need to calculate AMAX
C This is code from earlier experiments
CCC        AMAX = ZERO
C Split loops in two because only upper triangle is held
C Scan lower triangle by scanning up column of upper triangle.
CCC        DO 9876 K = 1, IPIV - NPIV - 1
CCC          AMAX = MAX(AMAX,ABS(A(POSPV1-K*NFRONT)))
CCC76      CONTINUE
C Scan upper triangle by scanning along row from first off-diagonal
CCC        DO 9878 K =  1, NFRONT - IPIV
CCC          AMAX = MAX(AMAX,ABS(A(POSPV1+K)))
CCC78      CONTINUE
C     Check size of 1 x 1 pivot and adjust if necessary
C          IF (ABS(A(POSPV1)).LT.CNTL(4)) THEN
C              PIVOT = CNTL(4)
C          IF (ABS(A(POSPV1)).LT.MAX(ULOC*AMAX,STCTOL)) THEN
C              PIVOT = MAX(ULOC*AMAX,STCTOL)

           IF (ABS(A(POSPV1)).LT.STCTOL) THEN
               PIVOT = STCTOL
              IF (A(POSPV1) .LT. ZERO) THEN
                 A(POSPV1) = -PIVOT
                 PIVOT     = -PIVOT
              ELSE
                 A(POSPV1) = PIVOT
              ENDIF
              NBSTATIC = NBSTATIC + 1
           ENDIF

C Perform the elimination using entry A(POSPV1) as pivot.
C We store U and D inverse.
C Later we store D inverse U which is passed to the solution entry.
           FLOPSB = FLOPSB + ONE
           A(POSPV1) = ONE/A(POSPV1)
           IF (A(POSPV1).LT.ZERO) NEIG = NEIG + 1

           J1 = POSPV1 + 1
           J2 = POSPV1 + NASS - (NPIV+1)
           IBEG = POSPV1 + NFRONT + 1
           IEND = APOSI + 2*NFRONT - 1
           DO 9880 JJ = J1,J2
C AMULT1 is used to hold the multiplier
              AMULT1 = -A(JJ)*A(POSPV1)
              JCOL = JJ
              FLOPSB = FLOPSB + (IEND-IBEG+1)*2 + 1
C     The following special comment forces vectorization on
C     Crays.
CDIR$            IVDEP
              DO 9870 IROW = IBEG,IEND
                 A(IROW) = A(IROW) + AMULT1*A(JCOL)
                 JCOL = JCOL + 1
 9870         CONTINUE
              A(JJ) = AMULT1
              IBEG = IBEG + NFRONT + 1
              IEND = IEND + NFRONT
 9880      CONTINUE
           NPIV = NPIV + 1
           MPIV = MPIV + 1
           NTOTPV = NTOTPV + 1

C     IF (LCASE) THEN
C     write(7,'(/A,I8/A,7I8)') '*** Frontal matrix at step',IASS,
C     *   'NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ',
C     *    NFRONT,NASS,NUMORG,APOS,IWPOS,NPIV,PIVSIZ
C     write(7,'(/A/(10I8))') 'IW array',
C     *      (IW(IWPOS+I-1),I=1,NFRONT)
C     write(7,'(/A/(5D16.8))') 'Frontal matrix',
C     *      ((A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT),J=1,NFRONT)
C     write(7,'(/A/(5D16.8))') 'Frontal matrix',
C     *      (A(APOS+I-1),I=1,NFRONT*NFRONT)
C     ENDIF

C Get next static pivot
           IF (MPIV.LT.SIZBLK) GO TO 9920
C********************************
C End of static pivoting loop   *
C********************************

  935   SCHUR = (NBLOC.LT.(NFRONT-NASS) .AND. NPIV.GE.NBLOC)
        IF (ICNTL(16).EQ.1) THEN
C Remove "zero" rows from fully summed rows within Schur complement
C ZCOUNT is count of "zero" rows
          ZCOUNT = 0
C APOS4 is beginning of block of fully summed uneliminated variables
          APOS4 = APOS + NPIV*NFRONT + NPIV

C Expand block so that lower triangle is included
C APOSB scans lower triangle by rows
C APOSC sweeps upper triangle by columns
          APOSB = APOS4 + NFRONT
          APOSC = APOS4 + 1
          DO 4444 I = 2,NASS-NPIV
            DO 4443 J = 1,I-1
              A(APOSB) = A(APOSC)
              APOSB = APOSB + 1
              APOSC = APOSC + NFRONT
 4443       CONTINUE
            APOSB = APOS4 + NFRONT*I
            APOSC = APOS4 + I
 4444     CONTINUE

C Remove any zero rows by swopping with "first" row so that all zero
C     rows will be swept to beginning of block
C Row I is the row currently being checked
          I = NASS - NPIV
C Also exchange integer information accordingly
 4445     CONTINUE
          IF (ZCOUNT.EQ.I) GO TO 4450
C Check row I for zero
C APOSB is beginning of row I
          APOSB = APOS4 + (I-1)*NFRONT
          DO 4446 J = 1,NFRONT-NPIV
            IF (ABS(A(APOSB+J-1)).GT.TOL) GO TO 4449
 4446     CONTINUE
C Row is all zero
          ZCOUNT = ZCOUNT + 1
C Swop row ZCOUNT with row I
          DO 4447 J = 1,NFRONT-NPIV
            A(APOSB+J-1) = A(APOS4+NFRONT*(ZCOUNT-1)+J-1)
 4447     CONTINUE
C Zero row ZCOUNT
          DO 4448 J = 1,NFRONT-NPIV
            A(APOS4+NFRONT*(ZCOUNT-1)+J-1) = ZERO
 4448     CONTINUE
C Swop integers
          ISWOP = IW(IWPOS+NPIV+ZCOUNT-1)
          IW(IWPOS+NPIV+ZCOUNT-1) = IW(IWPOS+NPIV+I-1)
          IW(IWPOS+NPIV+I-1) = ISWOP
          GO TO 4445
 4449     I = I - 1
          GO TO 4445
 4450     CONTINUE
        ELSE
          ZCOUNT = 0
        ENDIF
C Set order of Schur complement (including rows of delayed pivots)
C But not including "zero" rows
        NSC1 = NFRONT - NPIV - ZCOUNT

C       IF (LCASE) THEN
C         write(7,'(/A,I8/A,5I8)') '*** Frontal matrix fter step',IASS,
C    *   'NFRONT',
C    *    NFRONT
C         write(7,'(/A/(10I8))') 'IW array',
C    *      (IW(IWPOS+I-1),I=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      ((A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT),J=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      (A(APOS+I-1),I=1,NFRONT*NFRONT)
C       ENDIF

C Accumulate total number of delayed pivots in INFO(23)
C Do not accumulate at last step where matrix is singular
        IF (IASS.NE.NSTEPS) INFO(23) = INFO(23) + NASS - NPIV
C SET LSTAT
        IF (CNTL(4).GT.ZERO .AND. INFO(23).GT.CNTL(5)*N) LSTAT = .TRUE.

C Jump if no Schur complement to form ... just store factors.
        IF (NSC1.EQ.0) GO TO 1830

C Save space for factors and Schur complement (if appropriate).

        IF (.NOT.SCHUR) THEN
C We now compute triangular Schur complement not using BLAS.
C This Schur complement is placed directly on the stack.

C Remove these
          RLSPA = MAX(RLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
     +                      NSTKAC(1))
          TRLSPA = MAX(TRLSPA,INFO(40)+APOS+NFRONT*NFRONT-1+
     +                      TOTSTA(1))
          NSTKAC(1) = NSTKAC(1) + ((NSC1+1)*NSC1)/2
          NSTACK(1) = NSTACK(1) + ((NSC1+1)*NSC1)/2
          TOTSTA(1) = TOTSTA(1) + ((NSC1+1)*NSC1)/2

C Initialize Schur complement
C Copying from the back to avoid overwriting important data
          APOSI = APOS + NFRONT*NFRONT - 1
          DO 1370 JJ = 1,NFRONT-NPIV
            J = APOSI
            DO 1360 JJJ = 1,JJ
                A(ASTK) = A(J)
                ASTK = ASTK - 1
                J = J - 1
 1360       CONTINUE
            APOSI = APOSI - NFRONT
 1370     CONTINUE
C APOS4 is the position in A of the first entry in the Schur complement.
          APOS4 = ASTK + 1


C Perform pivoting operations.
C Initialize variables
          J1 = IWPOS
          LTWO = .FALSE.
          POSPV1 = APOS
          DO 1450 I1 = 1,NPIV
            IF (LTWO) GO TO 1440
            APOSI = APOS + (I1-1)*NFRONT + NASS
            J2 = APOS + NFRONT* (I1-1) + NFRONT - 1
CCC What happens here ??
            APOSC = APOS4 +
     *       ((NASS-NPIV-ZCOUNT)*(2*NFRONT-NPIV-ZCOUNT-NASS+1))/2
C Check to see if current pivot is 1 x 1  or 2 x 2.
            IF (IW(J1).GT.0) THEN
              FLOPSB = FLOPSB + (NFRONT-NASS) +
     *                          (NFRONT-NASS)* (NFRONT-NASS+1)
              DO 1410 JJ = APOSI,J2
                AMULT1 = -A(JJ)*A(POSPV1)
                DO 1400 JJJ = JJ,J2
                  A(APOSC) = A(APOSC) + AMULT1*A(JJJ)
                  APOSC = APOSC + 1
 1400           CONTINUE
                A(JJ) = AMULT1
 1410         CONTINUE
              J1 = J1 + 1
            ELSE
              POSPV2 = POSPV1 + NFRONT + 1
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6* (NFRONT-NASS) +
     +                 2* (NFRONT-NASS)* (NFRONT-NASS+1)
              DO 1430 JJ = APOSI,J2
                AMULT1 = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
                AMULT2 = -A(POSPV2)*A(JJ+NFRONT) - A(OFFDAG)*A(JJ)
                DO 1420 JJJ = JJ,J2
                  A(APOSC) = A(APOSC) + AMULT1*A(JJJ) +
     +                       AMULT2*A(JJJ+NFRONT)
                  APOSC = APOSC + 1
 1420           CONTINUE
                A(JJ) = AMULT1
                A(JJ+NFRONT) = AMULT2
 1430         CONTINUE
              J1 = J1 + 2
              POSPV1 = POSPV2
              LTWO = .TRUE.
              GO TO 1450
            END IF

 1440       LTWO = .FALSE.
C Move to next pivot position
            POSPV1 = POSPV1 + NFRONT + 1
 1450     CONTINUE

C         IF (LCASE) THEN
C           write(7,'(A,I8)') 'GEMM not used at stage',IASS
C           write(7,'(A/(5D16.8))') 'Stacking Schur',
C    *           (A(I),I=APOS4,APOS4+(NSC1*(NSC1+1))/2-1)
C         ENDIF

        ELSE

C Action if SCHUR is true, We now use GEMM.

C Since SCHUR is true, copy U,
C     divide factors by D, and generate Schur complement using GEMM.
C     Then compress factors (to upper trapezoidal), and stack
C     half of the Schur complement.

C APOS4 is position in A of first entry of Schur complement to be
C     updated using GEMM.
          APOS4 = APOS+NASS*(NFRONT+1)

C APOS3 is the position in A of the first entry in the copy of U.
        APOS3 = APOS+NASS*NFRONT

C Copy U and divide factors by D
C Initialize variables
        J1 = IWPOS
        LTWO = .FALSE.
        POSPV1 = APOS
          DO 1490 I = 1,NPIV
            IF (LTWO) GO TO 1480
            APOSI = APOS + (I-1)*NFRONT + NASS
C           POSELT = APOS3 + (I-1)* (NFRONT-NASS)
            POSELT = APOS3 + I - 1
C Check to see if current pivot is 1 x 1  or 2 x 2.
            IF (IW(J1).GT.0) THEN
              FLOPSB = FLOPSB + (NFRONT-NASS)
              DO 1460 JJ = APOSI,APOS + NFRONT*I - 1
                A(POSELT) = A(JJ)
                A(JJ) = -A(JJ)*A(POSPV1)
                POSELT = POSELT + NFRONT
 1460         CONTINUE
              J1 = J1 + 1
            ELSE
              POSPV2 = POSPV1 + NFRONT + 1
              OFFDAG = POSPV1 + 1
              FLOPSB = FLOPSB + 6* (NFRONT-NASS)
              DO 1470 JJ = APOSI,APOS + NFRONT*I - 1
                A(POSELT) = A(JJ)
C               A(POSELT+NFRONT-NASS) = A(JJ+NFRONT)
                A(POSELT+1) = A(JJ+NFRONT)
                A(JJ) = - (A(POSPV1)*A(JJ)+A(OFFDAG)*A(JJ+NFRONT))
                A(JJ+NFRONT) = -A(POSPV2)*A(JJ+NFRONT) -
     +                         A(OFFDAG)*A(POSELT)
                POSELT = POSELT + NFRONT
 1470         CONTINUE
              J1 = J1 + 2
              POSPV1 = POSPV2
              LTWO = .TRUE.
              GO TO 1490
            END IF

 1480       LTWO = .FALSE.
C Move to next pivot position
            POSPV1 = POSPV1 + NFRONT + 1
 1490     CONTINUE

C Now create Schur complement by using Level 3 BLAS GEMM.
C Jump if Schur complement is null.
C Increment FLOPSB
          FLOPSB = FLOPSB + NPIV* (NFRONT-NASS)**2 +
     *                      NPIV* (NFRONT-NASS)
C We divide the multiply into blocks to avoid too many extra
C    computations when using GEMM with a symmetric result.
C Block formed by GEMM has NBLOC rows.
          KBLK = ( NFRONT-NASS)/NBLOC
          L =  NFRONT - NASS
          DO 1500 KB = 1,KBLK
C Accumulate extra flops caused by using GEMM
            FLOPSX = FLOPSX + NBLOC* (NBLOC-1)* (NPIV)
            CALL DGEMM('N','N',L-(KB-1)*NBLOC,NBLOC,NPIV,ONE,
     +                 A(APOS+NASS+NBLOC*(KB-1)),NFRONT,
     +                 A(APOS3+NBLOC*(KB-1)*NFRONT),NFRONT,ONE,
     +                 A(APOS4+NBLOC*(NFRONT+1)*(KB-1)),NFRONT)
 1500     CONTINUE

C Calculate the block upper triangular part of the Schur complement.
          DO 1550 I = 1 + KBLK*NBLOC,L
C APOSA holds the index of the current entry of the matrix A used to
C     form the Schur complement as C = A*B
C APOSB holds the index of the current entry of the matrix B used to
C     form the Schur complement as C = A*B
C APOSC holds the index of the current entry of the matrix C used to
C     form the Schur complement as C = A*B
            APOSA = APOS + NASS
            APOSB = APOS3 +(I-1)*NFRONT
            APOSC = APOS4 + (I-1)*NFRONT - 1
            DO 1540 K = 1,NPIV
              DO 1530 J = I,L
                A(APOSC+J) = A(APOSC+J) + A(APOSA+J-1)*A(APOSB)
 1530         CONTINUE
              APOSA = APOSA + NFRONT
              APOSB = APOSB + 1
 1540       CONTINUE
 1550     CONTINUE

C Stack half of Schur complement.

C Stack reals
C Stack in reverse order to avoid need for compresses.
          JA1 = APOS+NFRONT*NFRONT-1
          NSTKAC(1) = NSTKAC(1) + ((NSC1+1)* (NSC1))/2
          NSTACK(1) = NSTACK(1) + ((NSC1+1)* (NSC1))/2
          TOTSTA(1) = TOTSTA(1) + ((NSC1+1)* (NSC1))/2
C Stack by rows
          DO 1710 I = NSC1,1,-1
            DO 1700 JJ = JA1,JA1-(NSC1-I),-1
              A(ASTK) = A(JJ)
              ASTK = ASTK - 1
 1700       CONTINUE
            JA1 = JA1 - NFRONT
 1710     CONTINUE

C         IF (LCASE) THEN
C           write(7,'(A,I8)') 'GEMM used at stage',IASS
C           write(7,'(A/(5D16.8))') 'Stacking Schur',
C    *           (A(I),I=ASTK+1,ASTK+(NSC1*(NSC1+1))/2)
C         ENDIF

C END of SCHUR being true action (started after label 1450)
        END IF

C Stack integers
        NSTKAC(2) = NSTKAC(2) + NSC1 + 1
        NSTACK(2) = NSTACK(2) + NSC1 + 1
        TOTSTA(2) = TOTSTA(2) + NSC1 + 1
C Record space needed to this point
 1830   IF (IASS.EQ.NSTEPS) THEN
          INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+NSTKAC(2))
          TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+TOTSTA(2))
          GO TO 2158
        ELSE
          INTSPA = MAX(INTSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+NSTKAC(2))
          TINSPA = MAX(TINSPA,IWPOS+NFRONT-1+(N-NTOTPV+2)+TOTSTA(2))
        ENDIF

C Check space and compress if necessary
  444   NST = 0
C +1 for length of stacked entry
        IF (NSC1.GT.0) NST = NSC1 + 1
        IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
C Compress integer storage
          CALL MA57PD(A,IW,ISTK,IINPUT,PTRIRN,.FALSE.)
          NCMPBI = NCMPBI + 1

          IF (IWPOS+NFRONT-1+(N-NTOTPV+2)+NST.GT.ISTK) THEN
C Still insufficient space after compress
            IF (ICNTL(8).NE.0) THEN
              HOLD(1) = 2
              HOLD(2) = NBLK
              HOLD(3) = NTWO
              HOLD(4) = INFO(23)
              HOLD(5) = NCMPBI
              HOLD(6) = NEIG
              HOLD(7) = MAXFRT
              HOLD(8) = IWPOS
              HOLD(9) = APOS
              HOLD(10) = APOSBB
              HOLD(11) = NSTKAC(1)
              HOLD(12) = NSTKAC(2)
              HOLD(13) = AINPUT
              HOLD(14) = IINPUT
              HOLD(15) = ISTK
              HOLD(16) = ASTK
              HOLD(17) = INTSPA
              HOLD(18) = RLSPA
              HOLD(19) = PTRIRN
              HOLD(20) = PTRA
              HOLD(21) = NTOTPV
              HOLD(22) = NPOTPV
              HOLD(23) = NUMORG
              HOLD(24) = NFRONT
              HOLD(25) = NASS
C             HOLD(26) = NCOL
              HOLD(27) = NPIV
              HOLD(28) = IASS
              HOLD(29) = TINSPA
              HOLD(30) = TRLSPA
              HOLD(31) = TOTSTA(1)
              HOLD(32) = TOTSTA(2)
              HOLD(33) = NSTACK(1)
              HOLD(34) = NSTACK(2)
              IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
              IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
              HOLD(37) = INFO(32)
              HOLD(38) = INFO(33)
              HOLD(39) = INFO(34)
              NSC1    = NFRONT-NPIV
              RINFO(3) =FLOPSA
              RINFO(4) =FLOPSB
              RINFO(5) =FLOPSX
              INFO(1) = 11
              HOLD(40) = NBSTATIC
              INFO(35) = HOLD(40)
            ELSE
              INFO(1)  = -4
              INFO(2)  = LIW
              INFO(18) = INTSPA
            ENDIF
            RETURN
          END IF
        END IF

        IF (NSC1.GT.0) THEN
          DO 1720 I = 1,NSC1
            IW(ISTK) = IW(IWPOS+NFRONT-I)
            ISTK = ISTK - 1
 1720     CONTINUE
C         write(11,'(A)') 'Stack integers'
          IW(ISTK) = NSC1
C         write(11,'(A/(10I8))') 'IW ....',(IW(I),I=ISTK,ISTK+NSC1)
          ISTK = ISTK - 1
        ENDIF

C       IF (LCASE) THEN
C         write(7,'(/A,I8/A,5I8)') '*** Frontal matrix end step',IASS,
C    *   'NFRONT',
C    *    NFRONT
C         write(7,'(/A/(10I8))') 'IW array',
C    *      (IW(IWPOS+I-1),I=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      ((A(APOS+I-1+(J-1)*NFRONT),I=1,NFRONT),J=1,NFRONT)
C         write(7,'(/A/(5D16.8))') 'Frontal matrix',
C    *      (A(APOS+I-1),I=1,NFRONT*NFRONT)
C       ENDIF

C Reset PPOS.
        DO 1840 JJ = IWPOS + NPIV,IWPOS + NFRONT - 1
          J = ABS(IW(JJ))
          PPOS(J) = N + 1
 1840   CONTINUE


C********************************
C    STORE FACTORS
C********************************
C Complete the integer information in the factors
 2158   IF (NPIV.EQ.0) GO TO 2159
        NBLK = NBLK + 1

        IW(IWPOS-2) = NFRONT
        IW(IWPOS-1) = NPIV
        IWPOS = IWPOS + NFRONT + 2

C Store information on the reals for the factors.
C We copy from A(JA1) to A(APOS2) ... the use of buffer space from
C   APOSBB to APOS ensures no overwrites.
        IF (INFO(1).EQ.-3) THEN
          INFO(40) = INFO(40) + (NPIV * (2*NFRONT-NPIV+1))/2
          GO TO 2159
        END IF

        APOS2 = APOSBB
C Store reals from full pivot.
        DO 2130 I = 1,NPIV
C JA1 points to the diagonal
          JA1 = APOS + (I-1)* (NFRONT+1)
          DO 2120 J = I,NPIV
            A(APOS2) = A(JA1)
            IF (A(APOS2).EQ.ZERO) INFO(32) = INFO(32) + 1
            APOS2 = APOS2 + 1
            JA1 = JA1 + 1
 2120     CONTINUE
 2130   CONTINUE
        RPOS = APOS2
C Store rectangle
        DO 2150 I = 1,NPIV
          JA1 = APOS + (I-1)*NFRONT + NPIV
          DO 2140 J = 1,NFRONT - NPIV
            A(APOS2) = A(JA1)
            APOS2 = APOS2 + 1
            JA1 = JA1 + 1
 2140     CONTINUE
 2150   CONTINUE
C Set APOSBB for next block of factors
        APOSBB = APOS2

C Check rectangle for zeros
        DO 2152 J = 1,NFRONT-NPIV
        APOS2 = RPOS+J-1
        ZCOL = 1
          DO 2151 I = 1,NPIV
            IF (A(APOS2).EQ.ZERO) INFO(33) = INFO(33)+1
            IF (A(APOS2).NE.ZERO) ZCOL = 0
            APOS2 = APOS2 + NFRONT - NPIV
 2151     CONTINUE
        IF (ZCOL.EQ.1) INFO(34) = INFO(34)+1
 2152   CONTINUE

 2159   IASS = IASS + 1
      IF (IASS.LE.NSTEPS) THEN
C Initialize to zero to avoid problem when calling MA57ED
        IW(IWPOS-2) = 0
        IW(IWPOS-1) = 0
        GO TO 2160
      ENDIF
C2160 CONTINUE
C
C End of loop on tree nodes.
C

      INFO(35) = NBSTATIC
      IF (INFO(1).EQ.-3) THEN
        INFO(2)  = LA
        INFO(17) = MAX(INFO(17),RLSPA)
        IF (ICNTL(7).EQ.4) INFO(17) = MAX(INFO(17),RLSPA + N)
        RETURN
      END IF
      GO TO 1000
 970  INFO(1) = -5
      INFO(2) = NTOTPV + 1
      IF (LDIAG.GT.0 .AND. LP.GE.0)
     *    WRITE(LP,99992) INFO(1),PIVOT,CNTL(2),INFO(2),ICNTL(7)
99992 FORMAT (/'*** Error message from routine MA57BD **',
     *       '   INFO(1) = ',I3/'Pivot has value ',D16.8,' when ',
     *       'CNTL(2) has value ',D16.8/
     *       'at stage',I11,2X,'when ICNTL(7) =',I3)
      RETURN
 980  INFO(1) = -6
      INFO(2) = NTOTPV + 1
      IF (LDIAG.GT.0 .AND. LP.GE.0)
     *    WRITE(LP,99993) INFO(1),INFO(2),ICNTL(7)
99993 FORMAT (/'*** Error message from routine MA57BD **',
     *       '   INFO(1) = ',I3/'Change in sign of pivot at stage',
     *       I10,2X,'when ICNTL(7) = ',I3)
      RETURN
 1000 NRLBDU = APOSBB - 1
      NIRBDU = IWPOS - 3
      IF (NTOTPV.NE.N) THEN
        INFO(1) = 4
        IF (LDIAG.GT.0 .AND. WP.GE.0)
     *      WRITE(WP,99994) INFO(1),NTOTPV
99994 FORMAT (/'*** Warning message from routine MA57BD **',
     *         '   INFO(1) =',I2/5X, 'Matrix is singular, rank =', I5)
      ENDIF

C Recent change was to remove condition that ICNTL(16) was equal to 1
C More recent change to remove deficiency test.  This change means that
C we now test there is sufficicent room to move the off-diagonal entries
C of the two by two pivots.
C 555 IF (NTOTPV.NE.N) THEN
C Check space (by this time there is nothing to compress)
  555 NRLBDU = APOSBB - 1
      NIRBDU = IWPOS - 3
      IF (NIRBDU+3*(N-NTOTPV) .GT. LIW
     +    .OR. NRLBDU+(N-NTOTPV)+NTWO .GT. LA) THEN
C I don't think this can happen ... at least I can't make it happen
C It is left in for "safety" :-)
C Still insufficient space after compress
          IF (ICNTL(8).NE.0) THEN
            HOLD(1) = 3
            HOLD(2) = NBLK
            HOLD(3) = NTWO
            HOLD(4) = INFO(23)
            HOLD(5) = NCMPBI
            HOLD(6) = NEIG
            HOLD(7) = MAXFRT
            HOLD(8) = IWPOS
            HOLD(9) = APOS
            HOLD(10) = APOSBB
            HOLD(11) = NSTKAC(1)
            HOLD(12) = NSTKAC(2)
            HOLD(13) = AINPUT
            HOLD(14) = IINPUT
            HOLD(15) = ISTK
            HOLD(16) = ASTK
            HOLD(17) = INTSPA
            HOLD(18) = RLSPA
            HOLD(19) = PTRIRN
            HOLD(20) = PTRA
            HOLD(21) = NTOTPV
            HOLD(22) = NPOTPV
            HOLD(23) = NUMORG
            HOLD(24) = NFRONT
            HOLD(25) = NASS
C             HOLD(26) = NCOL
            HOLD(27) = NPIV
            HOLD(28) = IASS
            HOLD(29) = TINSPA
            HOLD(30) = TRLSPA
            HOLD(31) = TOTSTA(1)
            HOLD(32) = TOTSTA(2)
            HOLD(33) = NSTACK(1)
            HOLD(34) = NSTACK(2)
            IF (ICNTL(7).EQ.2 .OR.ICNTL(7).EQ.3) HOLD(35) = ISNPIV
            IF (ICNTL(7).EQ.4) HOLD(36) = PHASE
            HOLD(37) = INFO(32)
            HOLD(38) = INFO(33)
            HOLD(39) = INFO(34)
            NSC1    = NFRONT-NPIV
            RINFO(3) =FLOPSA
            RINFO(4) =FLOPSB
            RINFO(5) =FLOPSX
            IF (NRLBDU+(N-NTOTPV)+NTWO .GT. LA) INFO(1) = 10
            IF (NIRBDU+3*(N-NTOTPV) .GT. LIW)  INFO(1) = 11
            HOLD(40) = NBSTATIC
            INFO(35) = HOLD(40)
          ELSE
            IF (NIRBDU+3*(N-NTOTPV) .GT. LIW) THEN
              INFO(1)  = -4
              INFO(2)  = LIW
              INFO(18) = MAX(INTSPA,NIRBDU+3*(N-NTOTPV))
            ELSE
              INFO(1)  = -3
              INFO(2) = LA
              INFO(17) = MAX(INFO(17),RLSPA,NRLBDU+(N-NTOTPV)+NTWO)
              IF (ICNTL(7).EQ.4) INFO(17) =
     +          MAX(INFO(17),RLSPA + N,NRLBDU+(N-NTOTPV)+NTWO)
            ENDIF
          ENDIF
          RETURN
        ENDIF
C     ENDIF

C Add explicit entries in factors for zero pivots (now set to 1.0)
C Initialize flag array to identify indices of zero pivots
      IF (N.NE.NTOTPV) THEN
        DO 3331 I = 1,N
          PPOS(I) = 0
 3331   CONTINUE
        IWPOS = 4
        DO 3332 I = 1,NBLK
          NFRONT = IW(IWPOS)
          NPIV = IW(IWPOS+1)
          DO 3330 J = IWPOS+2,IWPOS+NPIV+1
            PPOS(ABS(IW(J))) = 1
 3330     CONTINUE
          IWPOS = IWPOS + NFRONT + 2
 3332   CONTINUE
        K= 0
        DO 3333 I=1,N
          IF (PPOS(I).EQ.0) THEN
            K=K+1
            NBLK = NBLK + 1
            NRLBDU = NRLBDU+1
            A(NRLBDU) = ONE
            IW(NIRBDU+1) = 1
            IW(NIRBDU+2) = 1
            IW(NIRBDU+3) = I
            NIRBDU = NIRBDU+3
          ENDIF
 3333   CONTINUE
      ENDIF

C
      INFO(14) = NRLBDU
C     Move the off-diagonal entries of the 2x2 pivots within full blocks
C     to the end of A.
      IW(1) = NRLBDU + 1
      IW(2) = NRLBDU + NTWO
      INFO(15) = IW(2)
      IW(3) = NBLK
      INFO(31) = NBLK
C     Negate the entries of L, move the off-diagonal entries of the 2x2
C     pivots within full blocks to the end of A and update NRLBDU to
C     correspond
      CALL MA57WD(A,LA,IW,LIW,NRLBDU)
      INFO(16) = NIRBDU
      INFO(18) = INTSPA
      INFO(20) = TINSPA
      INFO(17) = RLSPA
      INFO(19) = TRLSPA
      INFO(21) = MAXFRT
      INFO(22) = NTWO
C     INFO(23)  .. computed as sum of NASS-NPIV
      INFO(24) = NEIG
      INFO(25) = NTOTPV
      INFO(28) = NCMPBR
      INFO(29) = NCMPBI
      RINFO(3) = FLOPSA
      RINFO(4) = FLOPSB
      RINFO(5) = FLOPSX
      IF (INFO(27).GT.0) THEN
        RINFO(14) = ZERO
        DO 332 I = 1,N
          RINFO(14) = MAX(RINFO(14),DIAG(I))
 332    CONTINUE
      ENDIF

      RETURN

      END


      SUBROUTINE MA57PD(A,IW,J1,J2,ITOP,REAL)
C This subroutine performs a very simple compress (block move).
C     Entries J1+1 to J2-1 (incl.) in A or IW as appropriate are moved
C     to occupy the positions immediately prior to position ITOP.
C A/IW hold the array being compressed.
C J1/J2 define the entries being moved.
C ITOP defines the position immediately after the positions to which
C     J1 to J2 are moved.
C REAL must be set by the user to .TRUE. if the move is on array A,
C     a value of .FALSE. will perform the move on A.
C     .. Scalar Arguments ..
      INTEGER ITOP,J1,J2
      LOGICAL REAL
C     ..
C     .. Array Arguments ..
      DOUBLE PRECISION A(*)
      INTEGER IW(*)
C     ..
C     .. Local Scalars ..
      INTEGER IPOS,JJ
C     ..
C     .. Executable Statements ..
      IF (J2.EQ.ITOP) GO TO 50
      IPOS = ITOP - 1
      IF (REAL) THEN
        DO 10 JJ = J2-1,J1+1,-1
          A(IPOS) = A(JJ)
          IPOS = IPOS - 1
   10   CONTINUE
      ELSE
        DO 20 JJ = J2-1,J1+1,-1
          IW(IPOS) = IW(JJ)
          IPOS = IPOS - 1
   20   CONTINUE
      ENDIF
      J2 = ITOP
      J1 = IPOS
   50 RETURN
      END
      SUBROUTINE MA57WD(A,LA,IW,LIW,NRLBDU)
C     Negate the entries of L, move the off-diagonal entries of the 2x2
C     pivots within full blocks to the end of A and update NRLBDU to
C     correspond.
      INTEGER LA,LIW
      DOUBLE PRECISION A(LA)
      INTEGER IW(LIW)

      INTEGER NRLBDU

C Constants
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D0)
C Local variables
      INTEGER APOS,IBLK,IROW,IWPOS,J,JPIV,NCOLS,NROWS
C APOS  Position in A of current diagonal entry.
C IBLK  Current block.
C IROW  Current row.
C IWPOS Current position in IW.
C J     Do loop variable
C JPIV  Used as a flag so that IPIV is incremented correctly after the
C       use of a 2 by 2 pivot.
C NCOLS Number of columns in the block.
C NROWS Number of rows in the block.

      APOS = 1
      IWPOS = 6
      DO 40 IBLK = 1,IW(3)
        NCOLS = IW(IWPOS-2)
        NROWS = IW(IWPOS-1)
        JPIV = 1
        DO 30 IROW = 1,NROWS
          JPIV = JPIV - 1
          IF (JPIV.EQ.1) GO TO 10
          IF (IW(IWPOS+IROW-1).LT.0) THEN
            JPIV = 2
            NRLBDU = NRLBDU + 1
            A(NRLBDU) = A(APOS+1)
            A(APOS+1) = ZERO
          END IF

   10     DO 20 J = APOS + 1,APOS + NROWS - IROW
            A(J) = -A(J)
   20     CONTINUE
          APOS = APOS + NROWS - IROW + 1
   30   CONTINUE
C Negate entries in rectangular block (was done earlier by MA47OD)
C       DO 35 J = APOS,APOS+NROWS*(NCOLS-NROWS)-1
C         A(J) = -A(J)
C  35   CONTINUE
        APOS = APOS + NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + NCOLS + 2
   40 CONTINUE
      END
      SUBROUTINE MA57XD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
C This subroutine performs forward elimination using the factors
C     stored in FACT/IFACT by MA57BD.
C It is designed for efficiency on one right-hand side.
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS)
      INTEGER IW1(N),ICNTL(20)
C N   must be set to the order of the matrix. It is not altered.
C FACT   must be set to hold the real values corresponding to the
C     factors. This must be unchanged since the preceding call to
C     MA57BD. It is not altered.
C LFACT  length of array FACT. It is not altered.
C IFACT  holds the integer indexing information for the matrix factors
C     in FACT. This must be unchanged since the preceding call to
C     MA57BD. It is not altered.
C LIFACT length of array IFACT. It is not altered.
C RHS on input, must be set to hold the right hand side vector.  On
C     return, it will hold the modified vector following forward
C     elimination.
C LHS must be set to the leading dimension of array RHS.
C W   used as workspace to hold the components of the right hand
C     sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C     larger than INFO(21) as returned from MA57BD.
C IW1 need not be set on entry. On exit IW1(I) (I = 1,NBLK), where
C     NBLK = IFACT(3) is the number of block pivots, will
C     hold pointers to the beginning of each block pivot in array IFACT.
C ICNTL Not referenced except:
C     ICNTL(13) Threshold on number of columns in a block for using
C           addressing using Level 2 and Level 3 BLAS.
C

C Procedures
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV

C Constant
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C Local variables
      INTEGER APOS,I,IBLK,II,IPIV,IRHS,IWPOS,J,J1,J2,K,K1,K2,
     +        NCOLS,NROWS
      DOUBLE PRECISION W1,W2
C
C APOS  Current position in array FACT.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C K     Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.

      APOS = 1
      IWPOS = 4
      DO 270 IBLK = 1,IFACT(3)

C Find the number of rows and columns in the block.
        IW1(IBLK) = IWPOS
        NCOLS = IFACT(IWPOS)
        NROWS = IFACT(IWPOS+1)
        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Load appropriate components of right-hand sides into W.
          DO 10 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = RHS(II)
   10     CONTINUE


C Treat diagonal block (direct addressing)
          CALL DTPSV('L','N','U',NROWS,FACT(APOS),W,1)
          APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (direct addressing)
C         IF (NCOLS.GT.NROWS) CALL DGEMM('N','N',NCOLS-NROWS,1,NROWS,
C    +                                  ONE,FACT(APOS),NCOLS-NROWS,
C    +                                  W,LW,ONE,W(NROWS+1),LW)
          IF (NCOLS.GT.NROWS) CALL DGEMV('N',NCOLS-NROWS,NROWS,
     +                                  ONE,FACT(APOS),NCOLS-NROWS,
     +                                  W,1,ONE,W(NROWS+1),1)
          APOS = APOS + NROWS* (NCOLS-NROWS)

C Reload W back into RHS.
          DO 35 I = 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   35     CONTINUE

        ELSE

C Perform operations using indirect addressing.

        J1 = IWPOS
        J2 = IWPOS + NROWS - 1


C Treat diagonal block (indirect addressing)
        DO 130 IPIV = 1,NROWS
          APOS = APOS + 1
          W1 = RHS(ABS(IFACT(J1)))
          K = APOS
          DO 100 J = J1+1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) - FACT(K)*W1
            K = K + 1
  100     CONTINUE
          APOS = K
          J1 = J1 + 1
  130   CONTINUE

C Treat off-diagonal block (indirect addressing)
C       J2 = IWPOS + NCOLS - 1
C       DO 136 IPIV = 1,NROWS
C         K = APOS
C         W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
C         DO 133 J = J1,J2
C           IRHS = ABS(IFACT(J))
C           RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
C           K = K + 1
C 133     CONTINUE
C         APOS = K
C 136   CONTINUE

C Loop unrolling
        J2 = IWPOS + NCOLS - 1
        DO 136 IPIV = 1,NROWS-1,2
          K1 = APOS
          K2 = APOS+NCOLS-NROWS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          W2 = RHS(ABS(IFACT(IWPOS+IPIV)))
          DO 133 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K1) + W2*FACT(K2)
            K1 = K1 + 1
            K2 = K2 + 1
  133     CONTINUE
          APOS = K2
  136   CONTINUE

        IF (MOD(NROWS,2).EQ.1) THEN
          K = APOS
          W1 = RHS(ABS(IFACT(IWPOS+IPIV-1)))
          DO 137 J = J1,J2
            IRHS = ABS(IFACT(J))
            RHS(IRHS) = RHS(IRHS) + W1*FACT(K)
            K = K + 1
  137     CONTINUE
          APOS = K
        ENDIF
      END IF

      IWPOS = IWPOS + NCOLS
  270 CONTINUE

      END


      SUBROUTINE MA57YD(N,FACT,LFACT,IFACT,LIFACT,RHS,LRHS,
     *                  W,LW,IW1,ICNTL)
C This subroutine performs backward elimination operations
C     using the factors stored in FACT/IFACT by MA57BD.
C It is designed for efficiency on one right-hand side.
      INTEGER N,LFACT
      DOUBLE PRECISION FACT(LFACT)
      INTEGER LIFACT,IFACT(LIFACT),LRHS,LW
      DOUBLE PRECISION W(LW),RHS(LRHS)
      INTEGER IW1(N),ICNTL(20)
C N    must be set to the order of the matrix. It is not altered.
C FACT    must be set to hold the real values corresponding to the
C      factors. This must be unchanged since the
C      preceding call to MA57BD. It is not altered.
C LFACT   length of array FACT. It is not altered.
C IFACT   holds the integer indexing information for the matrix factors
C      in FACT. This must be unchanged since the preceding call to
C      MA57BD. It is not altered.
C LIFACT  length of array IFACT. It is not altered.
C RHS  on entry, must be set to hold the right hand side modified by
C      the forward substitution operations. On exit, holds the
C      solution vector.
C LHS must be set to the leading dimension of array RHS.
C W    used as workspace to hold the components of the right hand
C      sides corresponding to current block pivotal rows.
C LW  must be set as the leading dimension of array W.  It need not be
C     larger than INFO(21) as returned from MA57BD.
C IW1  on entry IW1(I) (I = 1,NBLK), where  NBLK = IFACT(3) is the
C     number of block pivots, must hold pointers to the beginning of
C     each block pivot in array IFACT, as set by MA57X/XD. It is not
C     altered.
C ICNTL Not referenced except:
C     ICNTL(13) Threshold on number of columns in a block for using
C           addressing using Level 2 and Level 3 BLAS.

C Procedures
      INTRINSIC ABS
      EXTERNAL DGEMV,DTPSV

C Constants
      DOUBLE PRECISION ONE
      PARAMETER (ONE=1.0D0)
C
C Local variables.
      INTEGER APOS,APOS2,I,IBLK,II,IPIV,IRHS,IRHS1,
     +        IRHS2,IWPOS,J,JPIV,J1,J2,K,K2,LROW,NCOLS,NROWS
      DOUBLE PRECISION W1,W2
C APOS  Current position in array FACT.
C APOS2 Current position in array FACT for off-diagonal entry of 2x2
C       pivot.
C I     Temporary DO index
C IBLK  Index of block pivot.
C II    Temporary index.
C IPIV  Pivot index.
C IRHS  RHS index.
C IRHS1 RHS index.
C IRHS2 RHS index.
C IWPOS Position in IFACT of start of current index list.
C J     Temporary DO index
C JPIV  Has the value -1 for the first row of a 2 by 2 pivot and 1 for
C       the second.
C K     Temporary pointer to position in real array.
C J1    Position in IFACT of index of leading entry of row.
C J2    Position in IFACT of index of trailing entry of row.
C K     Temporary variable.
C LROW  Length of current row.
C NCOLS Number of columns in the block pivot.
C NROWS Number of rows in the block pivot.
C W1    RHS value.
C
      APOS = IFACT(1)
      APOS2 = IFACT(2)
C Run through block pivot rows in the reverse order.
      DO 380 IBLK = IFACT(3),1,-1

C Find the number of rows and columns in the block.
        IWPOS = IW1(IBLK)
        NCOLS = ABS(IFACT(IWPOS))
        NROWS = ABS(IFACT(IWPOS+1))
        APOS = APOS - NROWS* (NCOLS-NROWS)
        IWPOS = IWPOS + 2

        IF (NROWS.GT.4 .AND. NCOLS.GT.ICNTL(13)) THEN

C Perform operations using direct addressing.

C Load latter part of right-hand side into W.
          DO 5 I = NROWS + 1,NCOLS
            II = ABS(IFACT(IWPOS+I-1))
            W(I) = RHS(II)
    5     CONTINUE


C Multiply by the diagonal matrix (direct addressing)
          DO 10 IPIV = NROWS,1,-1
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            APOS = APOS - (NROWS+1-IPIV)
            W(IPIV) = RHS(IRHS)*FACT(APOS)
   10     CONTINUE
          JPIV = -1
          DO 20 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            IF (IRHS.LT.0) THEN
              IRHS1 = -IFACT(IWPOS+IPIV-1+JPIV)
              W(IPIV) = RHS(IRHS1)*FACT(APOS2) + W(IPIV)
              IF (JPIV.EQ.1) APOS2 = APOS2 - 1
              JPIV = -JPIV
            END IF

   20     CONTINUE

C Treat off-diagonal block (direct addressing)
          K = NCOLS - NROWS
          IF (K.GT.0) CALL DGEMV('T',K,NROWS,ONE,
     +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
     +                           W(NROWS+1),1,ONE,W,1)
C         IF (K.GT.0) CALL DGEMM('T','N',NROWS,1,K,ONE,
C    +                           FACT(APOS+(NROWS*(NROWS+1))/2),K,
C    +                           W(NROWS+1),LW,ONE,W,LW)

C Treat diagonal block (direct addressing)
          CALL DTPSV('L','T','U',NROWS,FACT(APOS),W,1)
C Reload W back into RHS.
          DO 60 I = 1,NROWS
            II = ABS(IFACT(IWPOS+I-1))
            RHS(II) = W(I)
   60     CONTINUE

        ELSE
C
C Perform operations using indirect addressing.
          J1 = IWPOS
          J2 = IWPOS + NCOLS - 1


C Multiply by the diagonal matrix (indirect addressing)
          JPIV = -1
          DO 210 IPIV = NROWS,1,-1
            IRHS = IFACT(IWPOS+IPIV-1)
            LROW = NROWS + 1 - IPIV

            IF (IRHS.GT.0) THEN
C 1 by 1 pivot.
              APOS = APOS - LROW
              RHS(IRHS) = RHS(IRHS)*FACT(APOS)
            ELSE
C 2 by 2 pivot
              IF (JPIV.EQ.-1) THEN
                IRHS1 = -IFACT(IWPOS+IPIV-2)
                IRHS2 = -IRHS
                APOS = APOS - LROW - LROW - 1
                W1 = RHS(IRHS1)*FACT(APOS) +
     +               RHS(IRHS2)*FACT(APOS2)
                RHS(IRHS2) = RHS(IRHS1)*FACT(APOS2) +
     +                         RHS(IRHS2)*FACT(APOS+LROW+1)
                RHS(IRHS1) = W1
                APOS2 = APOS2 - 1
              END IF
              JPIV = -JPIV
            END IF

  210     CONTINUE
          APOS = APOS + (NROWS* (NROWS+1))/2

C Treat off-diagonal block (indirect addressing)
C         KK = APOS
C         J1 = IWPOS + NROWS
C         DO 220 IPIV = 1,NROWS
C           IRHS = ABS(IFACT(IWPOS+IPIV-1))
C           W1 = RHS(IRHS)
C           K = KK
C           DO 215 J = J1,J2
C             W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
C             K = K + 1
C 215       CONTINUE
C           RHS(IRHS) = W1
C           KK = K
C 220     CONTINUE

C Loop unrolling
          K = APOS
          J1 = IWPOS + NROWS
          DO 220 IPIV = 1,NROWS-1,2
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W1 = RHS(IRHS)
            IRHS1 = ABS(IFACT(IWPOS+IPIV))
            W2 = RHS(IRHS1)
            K2 = K+(NCOLS-NROWS)
            DO 215 J = J1,J2
              II = ABS(IFACT(J))
              W1 = W1 + FACT(K)*RHS(II)
              W2 = W2 + FACT(K2)*RHS(II)
              K = K + 1
              K2 = K2 + 1
  215       CONTINUE
            RHS(IRHS) = W1
            RHS(IRHS1) = W2
            K = K2
  220     CONTINUE

          IF (MOD(NROWS,2).EQ.1) THEN
            IRHS = ABS(IFACT(IWPOS+IPIV-1))
            W1 = RHS(IRHS)
            DO 216 J = J1,J2
              W1 = W1 + FACT(K)*RHS(ABS(IFACT(J)))
              K = K + 1
  216       CONTINUE
            RHS(IRHS) = W1
          ENDIF

C Treat diagonal block (indirect addressing)
          J2 = IWPOS + NROWS - 1
          DO 260 IPIV = 1,NROWS
            IRHS = ABS(IFACT(J1-1))
            APOS = APOS - IPIV
            W1 = RHS(IRHS)
            K = APOS + 1
            DO 230 J = J1,J2
              W1 = W1 - FACT(K)*RHS(ABS(IFACT(J)))
              K = K + 1
  230       CONTINUE
            RHS(IRHS) = W1
            J1 = J1 - 1
  260     CONTINUE

        END IF

  380 CONTINUE

      END

      SUBROUTINE MA57VD(N,NZ,IRN,ICN,IW,LW,IPE,IQ,FLAG,IWFR,
     +                 ICNTL,INFO)
C Is identical to subroutine MA27GD.  Internal version for MA57.
C
C SORT PRIOR TO CALLING ANALYSIS ROUTINE MA27H/HD (internal MA57
C subroutine MA57H/HD).
C
C GIVEN THE POSITIONS OF THE OFF-DIAGONAL NON-ZEROS OF A SYMMETRIC
C     MATRIX, CONSTRUCT THE SPARSITY PATTERN OF THE OFF-DIAGONAL
C     PART OF THE WHOLE MATRIX (UPPER AND LOWER TRIANGULAR PARTS).
C EITHER ONE OF A PAIR (I,J),(J,I) MAY BE USED TO REPRESENT
C     THE PAIR. DIAGONAL ELEMENTS AND DUPLICATES ARE IGNORED.
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C NZ MUST BE SET TO THE NUMBER OF NON-ZEROS INPUT. IT IS NOT
C     ALTERED.
C IRN(I),I=1,2,...,NZ MUST BE SET TO THE ROW NUMBERS OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C ICN(I),I=1,2,...,NZ MUST BE SET TO THE COLUMN NUMBERS OF THE
C     NON-ZEROS ON INPUT. IT IS NOT ALTERED UNLESS IT IS EQUIVALENCED
C     TO IW (SEE DESCRIPTION OF IW).
C IW NEED NOT BE SET ON INPUT. ON OUTPUT IT CONTAINS LISTS OF
C     COLUMN INDICES, EACH LIST BEING HEADED BY ITS LENGTH.
C     IRN(1) MAY BE EQUIVALENCED TO IW(1) AND ICN(1) MAY BE
C     EQUIVALENCED TO IW(K), WHERE K.GT.NZ.
C LW MUST BE SET TO THE LENGTH OF IW. IT MUST BE AT LEAST 2*NZ+N.
C     IT IS NOT ALTERED.
C IPE NEED NOT BE SET ON INPUT. ON OUTPUT IPE(I) POINTS TO THE START OF
C     THE ENTRY IN IW FOR ROW I, OR IS ZERO IF THERE IS NO ENTRY.
C IQ NEED NOT BE SET.  ON OUTPUT IQ(I),I=1,N CONTAINS THE NUMBER OF
C     OFF-DIAGONAL N0N-ZEROS IN ROW I INCLUDING DUPLICATES.
C FLAG IS USED FOR WORKSPACE TO HOLD FLAGS TO PERMIT DUPLICATE ENTRIES
C     TO BE IDENTIFIED QUICKLY.
C IWFR NEED NOT BE SET ON INPUT. ON OUTPUT IT POINTS TO THE FIRST
C     UNUSED LOCATION IN IW.
C ICNTL is an INTEGER array of assumed size.
C INFO is an INTEGER array of assumed size.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NZ
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),ICN(*),IPE(N),IQ(N),IRN(*),IW(LW)
      INTEGER ICNTL(*),INFO(*)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,J,JN,K,K1,K2,L,LAST,LR,N1,NDUP
C     ..
C     .. Executable Statements ..
C
C INITIALIZE INFO(2) AND COUNT IN IPE THE
C     NUMBERS OF NON-ZEROS IN THE ROWS AND MOVE ROW AND COLUMN
C     NUMBERS INTO IW.
      INFO(2) = 0
      DO 10 I = 1,N
        IPE(I) = 0
   10 CONTINUE
      LR = NZ
      IF (NZ.EQ.0) GO TO 120
      DO 110 K = 1,NZ
        I = IRN(K)
        J = ICN(K)
        IF (I.LT.J) THEN
          IF (I.GE.1 .AND. J.LE.N) GO TO 90
        ELSE IF (I.GT.J) THEN
          IF (J.GE.1 .AND. I.LE.N) GO TO 90
        ELSE
          IF (I.GE.1 .AND. I.LE.N) GO TO 80
        END IF
        INFO(2) = INFO(2) + 1
        INFO(1) = 1
        IF (INFO(2).LE.1 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=60) INFO(1)
        END IF

   60   FORMAT (' *** WARNING MESSAGE FROM SUBROUTINE MA57AD',
     +          '  *** INFO(1) =',I2)

        IF (INFO(2).LE.10 .AND. ICNTL(2).GT.0) THEN
          WRITE (ICNTL(2),FMT=70) K,I,J
        END IF

   70   FORMAT (I6,'TH NON-ZERO (IN ROW',I6,' AND COLUMN',I6,
     +         ') IGNORED')

   80   I = 0
        J = 0
        GO TO 100

   90   IPE(I) = IPE(I) + 1
        IPE(J) = IPE(J) + 1
  100   IW(K) = J
        LR = LR + 1
        IW(LR) = I
  110 CONTINUE
C
C ACCUMULATE ROW COUNTS TO GET POINTERS TO ROW STARTS IN BOTH IPE AND IQ
C     AND INITIALIZE FLAG
  120 IQ(1) = 1
      N1 = N - 1
      IF (N1.LE.0) GO TO 140
      DO 130 I = 1,N1
        FLAG(I) = 0
        IF (IPE(I).EQ.0) IPE(I) = -1
        IQ(I+1) = IPE(I) + IQ(I) + 1
        IPE(I) = IQ(I)
  130 CONTINUE
  140 LAST = IPE(N) + IQ(N)
      FLAG(N) = 0
      IF (LR.GE.LAST) GO TO 160
      K1 = LR + 1
      DO 150 K = K1,LAST
        IW(K) = 0
  150 CONTINUE
  160 IPE(N) = IQ(N)
      IWFR = LAST + 1
      IF (NZ.EQ.0) GO TO 230
C
C RUN THROUGH PUTTING THE MATRIX ELEMENTS IN THE RIGHT PLACE
C     BUT WITH SIGNS INVERTED. IQ IS USED FOR HOLDING RUNNING POINTERS
C     AND IS LEFT HOLDING POINTERS TO ROW ENDS.
      DO 220 K = 1,NZ
        J = IW(K)
        IF (J.LE.0) GO TO 220
        L = K
        IW(K) = 0
        DO 210 ID = 1,NZ
          IF (L.GT.NZ) GO TO 170
          L = L + NZ
          GO TO 180

  170     L = L - NZ
  180     I = IW(L)
          IW(L) = 0
          IF (I.LT.J) GO TO 190
          L = IQ(J) + 1
          IQ(J) = L
          JN = IW(L)
          IW(L) = -I
          GO TO 200

  190     L = IQ(I) + 1
          IQ(I) = L
          JN = IW(L)
          IW(L) = -J
  200     J = JN
          IF (J.LE.0) GO TO 220
  210   CONTINUE
  220 CONTINUE
C
C RUN THROUGH RESTORING SIGNS, REMOVING DUPLICATES AND SETTING THE
C     MATE OF EACH NON-ZERO.
C NDUP COUNTS THE NUMBER OF DUPLICATE ELEMENTS.
  230 NDUP = 0
      DO 280 I = 1,N
        K1 = IPE(I) + 1
        K2 = IQ(I)
        IF (K1.LE.K2) GO TO 240
C ROW IS EMPTY. SET POINTER TO ZERO.
        IPE(I) = 0
        IQ(I) = 0
        GO TO 280
C ON ENTRY TO THIS LOOP FLAG(J).LT.I FOR J=1,2,...,N. DURING THE LOOP
C     FLAG(J) IS SET TO I IF A NON-ZERO IN COLUMN J IS FOUND. THIS
C     PERMITS DUPLICATES TO BE RECOGNIZED QUICKLY.
  240   DO 260 K = K1,K2
          J = -IW(K)
          IF (J.LE.0) GO TO 270
          L = IQ(J) + 1
          IQ(J) = L
          IW(L) = I
          IW(K) = J
          IF (FLAG(J).NE.I) GO TO 250
          NDUP = NDUP + 1
          IW(L) = 0
          IW(K) = 0
  250     FLAG(J) = I
  260   CONTINUE
  270   IQ(I) = IQ(I) - IPE(I)
        IF (NDUP.EQ.0) IW(K1-1) = IQ(I)
  280 CONTINUE
      IF (NDUP.EQ.0) GO TO 310
C
C COMPRESS IW TO REMOVE DUMMY ENTRIES CAUSED BY DUPLICATES.
      IWFR = 1
      DO 300 I = 1,N
        K1 = IPE(I) + 1
        IF (K1.EQ.1) GO TO 300
        K2 = IQ(I) + IPE(I)
        L = IWFR
        IPE(I) = IWFR
        IWFR = IWFR + 1
        DO 290 K = K1,K2
          IF (IW(K).EQ.0) GO TO 290
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
  290   CONTINUE
        IW(L) = IWFR - L - 1
  300 CONTINUE
  310 RETURN

      END
      SUBROUTINE MA57HD(N,IPE,IW,LW,IWFR,NV,NXT,LST,IPD,FLAG,IOVFLO,
     +                 NCMPA,FRATIO)
C Was identical to subroutine MA27HD.  Internal version for MA57.
C     Changes made in September 2009 because of bug in compress control
C     found by Nick.
C
C ANALYSIS SUBROUTINE
C
C GIVEN REPRESENTATION OF THE WHOLE MATRIX (EXCLUDING DIAGONAL)
C     PERFORM MINIMUM DEGREE ORDERING, CONSTRUCTING TREE POINTERS.
C     IT WORKS WITH SUPERVARIABLES WHICH ARE COLLECTIONS OF ONE OR MORE
C     VARIABLES, STARTING WITH SUPERVARIABLE I CONTAINING VARIABLE I FOR
C     I = 1,2,...,N. ALL VARIABLES IN A SUPERVARIABLE ARE ELIMINATED
C     TOGETHER. EACH SUPERVARIABLE HAS AS NUMERICAL NAME THAT OF ONE
C     OF ITS VARIABLES (ITS PRINCIPAL VARIABLE).
C
C N MUST BE SET TO THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) MUST BE SET TO POINT TO THE POSITION IN IW OF THE
C     START OF ROW I OR HAVE THE VALUE ZERO IF ROW I HAS NO OFF-
C     DIAGONAL NON-ZEROS. DURING EXECUTION IT IS USED AS FOLLOWS. IF
C     SUPERVARIABLE I IS ABSORBED INTO SUPERVARIABLE J THEN IPE(I)=-J.
C     IF SUPERVARIABLE I IS ELIMINATED THEN IPE(I) EITHER POINTS TO THE
C     LIST OF SUPERVARIABLES FOR CREATED ELEMENT I OR IS ZERO IF
C     THE CREATED ELEMENT IS NULL. IF ELEMENT I
C     IS ABSORBED INTO ELEMENT J THEN IPE(I)=-J.
C IW MUST BE SET ON ENTRY TO HOLD LISTS OF VARIABLES BY
C     ROWS, EACH LIST BEING HEADED BY ITS LENGTH.
C     DURING EXECUTION THESE LISTS ARE REVISED AND HOLD
C     LISTS OF ELEMENTS AND SUPERVARIABLES. THE ELEMENTS
C     ALWAYS HEAD THE LISTS. WHEN A SUPERVARIABLE
C     IS ELIMINATED ITS LIST IS REPLACED BY A LIST OF SUPERVARIABLES
C     IN THE NEW ELEMENT.
C LW MUST BE SET TO THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR MUST BE SET TO THE POSITION IN IW OF THE FIRST FREE VARIABLE.
C     IT IS REVISED DURING EXECUTION AND CONTINUES TO HAVE THIS MEANING.
C NV(JS) NEED NOT BE SET. DURING EXECUTION IT IS ZERO IF
C     JS IS NOT A PRINCIPAL VARIABLE AND IF IT IS IT HOLDS
C     THE NUMBER OF VARIABLES IN SUPERVARIABLE JS. FOR ELIMINATED
C     VARIABLES IT IS SET TO THE DEGREE AT THE TIME OF ELIMINATION.
C NXT(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE NEXT
C     SUPERVARIABLE HAVING THE SAME DEGREE AS JS, OR ZERO
C     IF IT IS LAST IN ITS LIST.
C LST(JS) NEED NOT BE SET. DURING EXECUTION IT IS THE
C     LAST SUPERVARIABLE HAVING THE SAME DEGREE AS JS OR
C     -(ITS DEGREE) IF IT IS FIRST IN ITS LIST.
C IPD(ID) NEED NOT BE SET. DURING EXECUTION IT
C     IS THE FIRST SUPERVARIABLE WITH DEGREE ID OR ZERO
C     IF THERE ARE NONE.
C FLAG IS USED AS WORKSPACE FOR ELEMENT AND SUPERVARIABLE FLAGS.
C     WHILE THE CODE IS FINDING THE DEGREE OF SUPERVARIABLE IS
C     FLAG HAS THE FOLLOWING VALUES.
C     A) FOR THE CURRENT PIVOT/NEW ELEMENT ME
C           FLAG(ME)=-1
C     B) FOR VARIABLES JS
C           FLAG(JS)=-1 IF JS IS NOT A PRINCIPAL VARIABLE
C           FLAG(JS)=0 IF JS IS A SUPERVARIABLE IN THE NEW ELEMENT
C           FLAG(JS)=NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C                 ELEMENT THAT HAS BEEN COUNTED IN THE DEGREE
C                 CALCULATION
C           FLAG(JS).GT.NFLG IF JS IS A SUPERVARIABLE NOT IN THE NEW
C                 ELEMENT THAT HAS NOT BEEN COUNTED IN THE DEGREE
C                 CALCULATION
C     C) FOR ELEMENTS IE
C           FLAG(IE)=-1 IF ELEMENT IE HAS BEEN MERGED INTO ANOTHER
C           FLAG(IE)=-NFLG IF ELEMENT IE HAS BEEN USED IN THE DEGREE
C                 CALCULATION FOR IS.
C           FLAG(IE).LT.-NFLG IF ELEMENT IE HAS NOT BEEN USED IN THE
C                 DEGREE CALCULATION FOR IS
C IOVFLO should be set to a high legitimate integer.  It is used as a
C        flag.
C NCMPA number of compresses.
C FRATIO is set to ICNTL(14)/100 and is the density of rows regarded as
C     dense.
C
C     .. Scalar Arguments ..
      DOUBLE PRECISION FRATIO
      INTEGER IWFR,LW,N,IOVFLO,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER FLAG(N),IPD(N),IPE(N),IW(LW),LST(N),NV(N),NXT(N)
C     ..
C     .. Local Scalars ..
      INTEGER I,ID,IDL,IDN,IE,IP,IS,JP,JP1,JP2,JS,K,K1,K2,KE,KP,KP0,KP1,
     +        KP2,KS,L,LEN,LIMIT,LN,LS,LWFR,MD,ME,ML,MS,NEL,NFLG,NP,
     +        NP0,NS,NVPIV,NVROOT,ROOT
C LIMIT  Limit on number of variables for putting node in root.
C NVROOT Number of variables in the root node
C ROOT   Index of the root node (N+1 if none chosen yet).
C     ..
C     .. External Subroutines ..
      EXTERNAL MA57ZD
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,MIN
C     ..
C If a column of the reduced matrix has relative density greater than
C CNTL(2), it is forced into the root. All such columns are taken to
C have sparsity pattern equal to their merged patterns, so the fill
C and operation counts may be overestimated.
C
C IS,JS,KS,LS,MS,NS are used to refer to supervariables.
C IE,JE,KE are used to refer to elements.
C IP,JP,KP,K,NP are used to point to lists of elements
C     or supervariables.
C ID is used for the degree of a supervariable.
C MD is used for the current minimum degree.
C IDN is used for the no. of variables in a newly created element
C NEL is used to hold the no. of variables that have been
C     eliminated.
C ME=MS is the name of the supervariable eliminated and
C     of the element created in the main loop.
C NFLG is used for the current flag value in array FLAG. It starts
C     with the value IOVFLO and is reduced by 1 each time it is used
C     until it has the value 2 when it is reset to the value IOVFLO.
C
C     .. Executable Statements ..
C Initializations
      DO 10 I = 1,N
        IPD(I) = 0
        NV(I) = 1
        FLAG(I) = IOVFLO
   10 CONTINUE
      MD = 1
      NCMPA = 0
      NFLG = IOVFLO
      NEL = 0
      ROOT = N+1
      NVROOT = 0
C
C Link together variables having same degree
      DO 30 IS = 1,N
        K = IPE(IS)
        IF (K.GT.0) THEN
          ID = IW(K) + 1
          NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
        ELSE
C We have a variable that can be eliminated at once because there is
C     no off-diagonal nonzero in its row.
          NEL = NEL + 1
          FLAG(IS) = -1
          NXT(IS) = 0
          LST(IS) = 0
        ENDIF
   30 CONTINUE

C
C Start of main loop
C
      DO 340 ML = 1,N

C Leave loop if all variables have been eliminated.
        IF (NEL+NVROOT+1.GE.N) GO TO 350
C
C Find next supervariable for elimination.
        DO 40 ID = MD,N
          MS = IPD(ID)
          IF (MS.GT.0) GO TO 50
   40   CONTINUE
   50   MD = ID
C Nvpiv holds the number of variables in the pivot.
        NVPIV = NV(MS)
C
C Remove chosen variable from linked list
        NS = NXT(MS)
        NXT(MS) = 0
        LST(MS) = 0
        IF (NS.GT.0) LST(NS) = -ID
        IPD(ID) = NS
        ME = MS
        NEL = NEL + NVPIV
C IDN holds the degree of the new element.
        IDN = 0
C
C Run through the list of the pivotal supervariable, setting tree
C     pointers and constructing new list of supervariables.
C KP is a pointer to the current position in the old list.
        KP = IPE(ME)
        FLAG(MS) = -1
C IP points to the start of the new list.
        IP = IWFR
C LEN holds the length of the list associated with the pivot.
        LEN = IW(KP)
        DO 140 KP1 = 1,LEN
          KP = KP + 1
          KE = IW(KP)
C Jump if KE is an element that has not been merged into another.
          IF (FLAG(KE).LE.-2) GO TO 60
C Jump if KE is an element that has been merged into another or is
C     a supervariable that has been eliminated.
          IF (FLAG(KE).LE.0) THEN
             IF (IPE(KE).NE.-ROOT) GO TO 140
C KE has been merged into the root
             KE = ROOT
             IF (FLAG(KE).LE.0) GO TO 140
          END IF
C We have a supervariable. Prepare to search rest of list.
          JP = KP - 1
          LN = LEN - KP1 + 1
          IE = MS
          GO TO 70
C Search variable list of element KE, using JP as a pointer to it.
   60     IE = KE
          JP = IPE(IE)
          LN = IW(JP)
C
C Search for different supervariables and add them to the new list,
C     compressing when necessary. This loop is executed once for
C     each element in the list and once for all the supervariables
C     in the list.
   70     DO 130 JP1 = 1,LN
            JP = JP + 1
            IS = IW(JP)
C Jump if IS is not a principal variable or has already been counted.
            IF (FLAG(IS).LE.0) THEN
               IF (IPE(IS).EQ.-ROOT) THEN
C IS has been merged into the root
                  IS = ROOT
                  IW(JP) = ROOT
                  IF (FLAG(IS).LE.0) GO TO 130
               ELSE
                  GO TO 130
               END IF
            END IF
            FLAG(IS) = 0

C To fix Nick bug need to add one here to store (eventually) length
C     of new row
            IF (IWFR .GE. LW-1) THEN
C Logic was previously as below
CCC         IF (IWFR.LT.LW) GO TO 100
C Prepare for compressing IW by adjusting pointers and
C     lengths so that the lists being searched in the inner and outer
C     loops contain only the remaining entries.
              IPE(MS) = KP
              IW(KP) = LEN - KP1
              IPE(IE) = JP
              IW(JP) = LN - JP1
C Compress IW
              CALL MA57ZD(N,IPE,IW,IP-1,LWFR,NCMPA)
C Copy new list forward
              JP2 = IWFR - 1
              IWFR = LWFR
              IF (IP.GT.JP2) GO TO 90
              DO 80 JP = IP,JP2
                IW(IWFR) = IW(JP)
                IWFR = IWFR + 1
   80         CONTINUE
C Adjust pointers for the new list and the lists being searched.
   90         IP = LWFR
              JP = IPE(IE)
              KP = IPE(ME)
            ENDIF

C Store IS in new list.
            IW(IWFR) = IS
            IDN = IDN + NV(IS)
            IWFR = IWFR + 1
C Remove IS from degree linked list
            LS = LST(IS)
            LST(IS) = 0
            NS = NXT(IS)
            NXT(IS) = 0
            IF (NS.GT.0) LST(NS) = LS
            IF (LS.LT.0) THEN
              LS = -LS
              IPD(LS) = NS
            ELSE IF (LS.GT.0) THEN
              NXT(LS) = NS
            END IF
  130     CONTINUE
C Jump if we have just been searching the variables at the end of
C     the list of the pivot.
          IF (IE.EQ.MS) GO TO 150
C Set tree pointer and flag to indicate element IE is absorbed into
C     new element ME.
          IPE(IE) = -ME
          FLAG(IE) = -1
  140   CONTINUE

C Store the degree of the pivot.
  150   NV(MS) = IDN + NVPIV

C Jump if new element is null.
        IF (IWFR.EQ.IP) THEN
          IPE(ME) = 0
          GO TO 340
        ENDIF

        K1 = IP
        K2 = IWFR - 1
C
C Run through new list of supervariables revising each associated list,
C     recalculating degrees and removing duplicates.
        LIMIT = NINT(FRATIO*(N-NEL))
        DO 310 K = K1,K2
          IS = IW(K)
          IF (IS.EQ.ROOT) GO TO 310
          IF (NFLG.GT.2) GO TO 170
C Reset FLAG values to +/-IOVFLO.
          DO 160 I = 1,N
            IF (FLAG(I).GT.0) FLAG(I) = IOVFLO
            IF (FLAG(I).LE.-2) FLAG(I) = -IOVFLO
  160     CONTINUE
          NFLG = IOVFLO
C Reduce NFLG by one to cater for this supervariable.
  170     NFLG = NFLG - 1
C Begin with the degree of the new element. Its variables must always
C     be counted during the degree calculation and they are already
C     flagged with the value 0.
          ID = IDN
C Run through the list associated with supervariable IS
          KP1 = IPE(IS) + 1
C NP points to the next entry in the revised list.
          NP = KP1
          KP2 = IW(KP1-1) + KP1 - 1
          DO 220 KP = KP1,KP2
            KE = IW(KP)
C Test whether KE is an element, a redundant entry or a supervariable.
            IF (FLAG(KE).EQ.-1) THEN
              IF (IPE(KE).NE.-ROOT) GO TO 220
C KE has been merged into the root
              KE = ROOT
              IW(KP) = ROOT
              IF (FLAG(KE).EQ.-1) GO TO 220
            END IF
            IF (FLAG(KE).GE.0) GO TO 230
C Search list of element KE, revising the degree when new variables
C     found.
            JP1 = IPE(KE) + 1
            JP2 = IW(JP1-1) + JP1 - 1
            IDL = ID
            DO 190 JP = JP1,JP2
              JS = IW(JP)
C Jump if JS has already been counted.
              IF (FLAG(JS).LE.NFLG) GO TO 190
              ID = ID + NV(JS)
              FLAG(JS) = NFLG
  190       CONTINUE
C Jump if one or more new supervariables were found.
            IF (ID.GT.IDL) GO TO 210
C Check whether every variable of element KE is in new element ME.
            DO 200 JP = JP1,JP2
              JS = IW(JP)
              IF (FLAG(JS).NE.0) GO TO 210
  200       CONTINUE
C Set tree pointer and FLAG to indicate that element KE is absorbed
C     into new element ME.
            IPE(KE) = -ME
            FLAG(KE) = -1
            GO TO 220
C Store element KE in the revised list for supervariable IS and flag it.
  210       IW(NP) = KE
            FLAG(KE) = -NFLG
            NP = NP + 1
  220     CONTINUE
          NP0 = NP
          GO TO 250
C Treat the rest of the list associated with supervariable IS. It
C     consists entirely of supervariables.
  230     KP0 = KP
          NP0 = NP
          DO 240 KP = KP0,KP2
            KS = IW(KP)
            IF (FLAG(KS).LE.NFLG) THEN
               IF (IPE(KS).EQ.-ROOT) THEN
                  KS = ROOT
                  IW(KP) = ROOT
                  IF (FLAG(KS).LE.NFLG) GO TO 240
               ELSE
                  GO TO 240
               END IF
            END IF
C Add to degree, flag supervariable KS and add it to new list.
            ID = ID + NV(KS)
            FLAG(KS) = NFLG
            IW(NP) = KS
            NP = NP + 1
  240     CONTINUE
C Move first supervariable to end of list, move first element to end
C     of element part of list and add new element to front of list.
  250     IF (ID.GE.LIMIT) GO TO 295
          IW(NP) = IW(NP0)
          IW(NP0) = IW(KP1)
          IW(KP1) = ME
C Store the new length of the list.
          IW(KP1-1) = NP - KP1 + 1
C
C Check whether row is is identical to another by looking in linked
C     list of supervariables with degree ID at those whose lists have
C     first entry ME. Note that those containing ME come first so the
C     search can be terminated when a list not starting with ME is
C     found.
          JS = IPD(ID)
          DO 280 L = 1,N
            IF (JS.LE.0) GO TO 300
            KP1 = IPE(JS) + 1
            IF (IW(KP1).NE.ME) GO TO 300
C JS has same degree and is active. Check if identical to IS.
            KP2 = KP1 - 1 + IW(KP1-1)
            DO 260 KP = KP1,KP2
              IE = IW(KP)
C Jump if IE is a supervariable or an element not in the list of IS.
              IF (ABS(FLAG(IE)+0).GT.NFLG) GO TO 270
  260       CONTINUE
            GO TO 290

  270       JS = NXT(JS)
  280     CONTINUE
C Supervariable amalgamation. Row IS is identical to row JS.
C Regard all variables in the two supervariables as being in IS. Set
C     tree pointer, FLAG and NV entries.
  290     IPE(JS) = -IS
          NV(IS) = NV(IS) + NV(JS)
          NV(JS) = 0
          FLAG(JS) = -1
C Replace JS by IS in linked list.
          NS = NXT(JS)
          LS = LST(JS)
          IF (NS.GT.0) LST(NS) = IS
          IF (LS.GT.0) NXT(LS) = IS
          LST(IS) = LS
          NXT(IS) = NS
          LST(JS) = 0
          NXT(JS) = 0
          IF (IPD(ID).EQ.JS) IPD(ID) = IS
          GO TO 310
C Treat IS as full. Merge it into the root node.
  295     IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IW(K) = ROOT
            IPE(IS) = -ROOT
            NV(ROOT) = NV(ROOT) + NV(IS)
            NV(IS) = 0
            FLAG(IS) = -1
          END IF
          NVROOT = NV(ROOT)
          GO TO 310
C Insert IS into linked list of supervariables of same degree.
  300     NS = IPD(ID)
          IF (NS.GT.0) LST(NS) = IS
          NXT(IS) = NS
          IPD(ID) = IS
          LST(IS) = -ID
          MD = MIN(MD,ID)
  310   CONTINUE

C
C Reset flags for supervariables in newly created element and
C     remove those absorbed into others.
        DO 320 K = K1,K2
          IS = IW(K)
          IF (NV(IS).EQ.0) GO TO 320
          FLAG(IS) = NFLG
          IW(IP) = IS
          IP = IP + 1
  320   CONTINUE

        FLAG(ME) = -NFLG
C Move first entry to end to make room for length.
        IW(IP) = IW(K1)
        IW(K1) = IP - K1
C Set pointer for new element and reset IWFR.
        IPE(ME) = K1
        IWFR = IP + 1

C  End of main loop
  340 CONTINUE
C

C Absorb any remaining variables into the root
  350 DO 360 IS = 1,N
        IF(NXT(IS).NE.0 .OR. LST(IS).NE.0) THEN
          IF (NVROOT.EQ.0) THEN
            ROOT = IS
            IPE(IS) = 0
          ELSE
            IPE(IS) = -ROOT
          END IF
          NVROOT = NVROOT + NV(IS)
          NV(IS) = 0
         END IF
  360 CONTINUE
C Link any remaining elements to the root
      DO 370 IE = 1,N
        IF (IPE(IE).GT.0) IPE(IE) = -ROOT
  370 CONTINUE
      IF(NVROOT.GT.0)NV(ROOT)=NVROOT
      END
      SUBROUTINE MA57ZD(N,IPE,IW,LW,IWFR,NCMPA)
C Is identical to subroutine MA27UD.  Internal version for MA57.
C COMPRESS LISTS HELD BY MA27H/HD (MA57H/HD) IN IW AND ADJUST POINTERS
C     IN IPE TO CORRESPOND.
C N IS THE MATRIX ORDER. IT IS NOT ALTERED.
C IPE(I) POINTS TO THE POSITION IN IW OF THE START OF LIST I OR IS
C     ZERO IF THERE IS NO LIST I. ON EXIT IT POINTS TO THE NEW POSITION.
C IW HOLDS THE LISTS, EACH HEADED BY ITS LENGTH. ON OUTPUT THE SAME
C     LISTS ARE HELD, BUT THEY ARE NOW COMPRESSED TOGETHER.
C LW HOLDS THE LENGTH OF IW. IT IS NOT ALTERED.
C IWFR NEED NOT BE SET ON ENTRY. ON EXIT IT POINTS TO THE FIRST FREE
C     LOCATION IN IW.
C     ON RETURN IT IS SET TO THE FIRST FREE LOCATION IN IW.
C NCMPA is number of compresses.
C
C     .. Scalar Arguments ..
      INTEGER IWFR,LW,N,NCMPA
C     ..
C     .. Array Arguments ..
      INTEGER IPE(N),IW(LW)
C     ..
C     .. Local Scalars ..
      INTEGER I,IR,K,K1,K2,LWFR
C     ..
C     .. Executable Statements ..
      NCMPA = NCMPA + 1
C PREPARE FOR COMPRESSING BY STORING THE LENGTHS OF THE
C     LISTS IN IPE AND SETTING THE FIRST ENTRY OF EACH LIST TO
C     -(LIST NUMBER).
      DO 10 I = 1,N
        K1 = IPE(I)
        IF (K1.LE.0) GO TO 10
        IPE(I) = IW(K1)
        IW(K1) = -I
   10 CONTINUE
C
C COMPRESS
C IWFR POINTS JUST BEYOND THE END OF THE COMPRESSED FILE.
C LWFR POINTS JUST BEYOND THE END OF THE UNCOMPRESSED FILE.
      IWFR = 1
      LWFR = IWFR
      DO 60 IR = 1,N
        IF (LWFR.GT.LW) GO TO 70
C SEARCH FOR THE NEXT NEGATIVE ENTRY.
        DO 20 K = LWFR,LW
          IF (IW(K).LT.0) GO TO 30
   20   CONTINUE
        GO TO 70
C PICK UP ENTRY NUMBER, STORE LENGTH IN NEW POSITION, SET NEW POINTER
C     AND PREPARE TO COPY LIST.
   30   I = -IW(K)
        IW(IWFR) = IPE(I)
        IPE(I) = IWFR
        K1 = K + 1
        K2 = K + IW(IWFR)
        IWFR = IWFR + 1
        IF (K1.GT.K2) GO TO 50
C COPY LIST TO NEW POSITION.
        DO 40 K = K1,K2
          IW(IWFR) = IW(K)
          IWFR = IWFR + 1
   40   CONTINUE
   50   LWFR = K2 + 1
   60 CONTINUE
   70 RETURN

      END
