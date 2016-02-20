* STRSY2 solves the real Sylvester matrix equation (unblocked algorithm)
*
* This routine is an exact copy of LAPACK's strsyl.f.
* It serves as an unblocked kernel in the recursive algorithms.
      SUBROUTINE STRSY2( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $                   LDC, SCALE, INFO )
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      REAL               SCALE
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
      REAL               ZERO, ONE
      PARAMETER          ( ZERO = 0.0E+0, ONE = 1.0E+0 )
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            IERR, J, K, K1, K2, KNEXT, L, L1, L2, LNEXT
      REAL               A11, BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,
     $                   SMLNUM, SUML, SUMR, XNORM
      REAL               DUM( 1 ), VEC( 2, 2 ), X( 2, 2 )
      LOGICAL            LSAME
      REAL               SDOT, SLAMCH, SLANGE
      EXTERNAL           LSAME, SDOT, SLAMCH, SLANGE
      EXTERNAL           SLABAD, SLALN2, SLASY2, SSCAL, XERBLA
      INTRINSIC          ABS, MAX, MIN, REAL
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'T' ) .AND. .NOT.
     $    LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'T' ) .AND. .NOT.
     $         LSAME( TRANB, 'C' ) ) THEN
         INFO = -2
      ELSE IF( ISGN.NE.1 .AND. ISGN.NE.-1 ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'STRSY2', -INFO )
         RETURN
      END IF
      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
      EPS = SLAMCH( 'P' )
      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL SLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*REAL( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*SLANGE( 'M', M, M, A, LDA, DUM ),
     $       EPS*SLANGE( 'M', N, N, B, LDB, DUM ) )
      SGN = ISGN
      IF( NOTRNA .AND. NOTRNB ) THEN
         LNEXT = 1
         DO 70 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 70
            IF( L.EQ.N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L ).NE.ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
            KNEXT = M
            DO 60 K = M, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 60
               IF( K.EQ.1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 ).NE.ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 10 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   10                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 20 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   20                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                         C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 40 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   40                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
                  CALL SLASY2( .FALSE., .FALSE., ISGN, 2, 2,
     $                         A( K1, K1 ), LDA, B( L1, L1 ), LDB, VEC,
     $                         2, SCALOC, X, 2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 50 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   50                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
   60       CONTINUE
   70    CONTINUE
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
         LNEXT = 1
         DO 130 L = 1, N
            IF( L.LT.LNEXT )
     $         GO TO 130
            IF( L.EQ.N ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L+1, L ).NE.ZERO ) THEN
                  L1 = L
                  L2 = L + 1
                  LNEXT = L + 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L + 1
               END IF
            END IF
            KNEXT = 1
            DO 120 K = 1, M
               IF( K.LT.KNEXT )
     $            GO TO 120
               IF( K.EQ.M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K ).NE.ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 80 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   80                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 90 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
   90                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 100 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  100                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K1, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L1 ), 1 )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( L1-1, C( K2, 1 ), LDC, B( 1, L2 ), 1 )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
                  CALL SLASY2( .TRUE., .FALSE., ISGN, 2, 2, A( K1, K1 ),
     $                         LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                         2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 110 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  110                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
  120       CONTINUE
  130    CONTINUE
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
         LNEXT = N
         DO 190 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 190
            IF( L.EQ.1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 ).NE.ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
            KNEXT = 1
            DO 180 K = 1, M
               IF( K.LT.KNEXT )
     $            GO TO 180
               IF( K.EQ.M ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K+1, K ).NE.ZERO ) THEN
                     K1 = K
                     K2 = K + 1
                     KNEXT = K + 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K + 1
                  END IF
               END IF
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC,
     $                         B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 140 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  140                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  CALL SLALN2( .TRUE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 150 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  150                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 160 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  160                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K1 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( K1-1, A( 1, K2 ), 1, C( 1, L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                         B( L2, MIN(L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
                  CALL SLASY2( .TRUE., .TRUE., ISGN, 2, 2, A( K1, K1 ),
     $                         LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                         2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 170 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  170                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
  180       CONTINUE
  190    CONTINUE
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
         LNEXT = N
         DO 250 L = N, 1, -1
            IF( L.GT.LNEXT )
     $         GO TO 250
            IF( L.EQ.1 ) THEN
               L1 = L
               L2 = L
            ELSE
               IF( B( L, L-1 ).NE.ZERO ) THEN
                  L1 = L - 1
                  L2 = L
                  LNEXT = L - 2
               ELSE
                  L1 = L
                  L2 = L
                  LNEXT = L - 1
               END IF
            END IF
            KNEXT = M
            DO 240 K = M, 1, -1
               IF( K.GT.KNEXT )
     $            GO TO 240
               IF( K.EQ.1 ) THEN
                  K1 = K
                  K2 = K
               ELSE
                  IF( A( K, K-1 ).NE.ZERO ) THEN
                     K1 = K - 1
                     K2 = K
                     KNEXT = K - 2
                  ELSE
                     K1 = K
                     K2 = K
                     KNEXT = K - 1
                  END IF
               END IF
               IF( L1.EQ.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( M-K1, A( K1, MIN(K1+1, M ) ), LDA,
     $                   C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L1, C( K1, MIN( L1+1, N ) ), LDC,
     $                         B( L1, MIN( L1+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SCALOC = ONE
                  A11 = A( K1, K1 ) + SGN*B( L1, L1 )
                  DA11 = ABS( A11 )
                  IF( DA11.LE.SMIN ) THEN
                     A11 = SMIN
                     DA11 = SMIN
                     INFO = 1
                  END IF
                  DB = ABS( VEC( 1, 1 ) )
                  IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                     IF( DB.GT.BIGNUM*DA11 )
     $                  SCALOC = ONE / DB
                  END IF
                  X( 1, 1 ) = ( VEC( 1, 1 )*SCALOC ) / A11
                  IF( SCALOC.NE.ONE ) THEN
                     DO 200 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  200                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
               ELSE IF( L1.EQ.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, A( K1, K1 ),
     $                         LDA, ONE, ONE, VEC, 2, -SGN*B( L1, L1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 210 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  210                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K2, L1 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.EQ.K2 ) THEN
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                         C( MIN( K1+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = SGN*( C( K1, L1 )-( SUML+SGN*SUMR ) )
                  SUML = SDOT( M-K1, A( K1, MIN( K1+1, M ) ), LDA,
     $                         C( MIN( K1+1, M ), L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = SGN*( C( K1, L2 )-( SUML+SGN*SUMR ) )
                  CALL SLALN2( .FALSE., 2, 1, SMIN, ONE, B( L1, L1 ),
     $                         LDB, ONE, ONE, VEC, 2, -SGN*A( K1, K1 ),
     $                         ZERO, X, 2, SCALOC, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 220 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  220                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 2, 1 )
               ELSE IF( L1.NE.L2 .AND. K1.NE.K2 ) THEN
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 1 ) = C( K1, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K1, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K1, MIN( L2+1, N ) ), LDC,
     $                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 1, 2 ) = C( K1, L2 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L1 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                         B( L1, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 1 ) = C( K2, L1 ) - ( SUML+SGN*SUMR )
                  SUML = SDOT( M-K2, A( K2, MIN( K2+1, M ) ), LDA,
     $                         C( MIN( K2+1, M ), L2 ), 1 )
                  SUMR = SDOT( N-L2, C( K2, MIN( L2+1, N ) ), LDC,
     $                         B( L2, MIN( L2+1, N ) ), LDB )
                  VEC( 2, 2 ) = C( K2, L2 ) - ( SUML+SGN*SUMR )
                  CALL SLASY2( .FALSE., .TRUE., ISGN, 2, 2, A( K1, K1 ),
     $                         LDA, B( L1, L1 ), LDB, VEC, 2, SCALOC, X,
     $                         2, XNORM, IERR )
                  IF( IERR.NE.0 )
     $               INFO = 1
                  IF( SCALOC.NE.ONE ) THEN
                     DO 230 J = 1, N
                        CALL SSCAL( M, SCALOC, C( 1, J ), 1 )
  230                CONTINUE
                     SCALE = SCALE*SCALOC
                  END IF
                  C( K1, L1 ) = X( 1, 1 )
                  C( K1, L2 ) = X( 1, 2 )
                  C( K2, L1 ) = X( 2, 1 )
                  C( K2, L2 ) = X( 2, 2 )
               END IF
  240       CONTINUE
  250    CONTINUE
      END IF
      RETURN
      END
