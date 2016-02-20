* CTRSY2 solves the complex Sylvester matrix equation (unblocked algorithm)
*
* This routine is an exact copy of LAPACK's ctrsyl.f.
* It serves as an unblocked kernel in the recursive algorithms.
      SUBROUTINE CTRSY2( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $                   LDC, SCALE, INFO )
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      REAL               SCALE
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            J, K, L
      REAL               BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,
     $                   SMLNUM
      COMPLEX            A11, SUML, SUMR, VEC, X11
      REAL               DUM( 1 )
      LOGICAL            LSAME
      REAL               CLANGE, SLAMCH
      COMPLEX            CDOTC, CDOTU, CLADIV
      EXTERNAL           LSAME, CLANGE, SLAMCH, CDOTC, CDOTU, CLADIV
      EXTERNAL           CSSCAL, SLABAD, XERBLA
      INTRINSIC          ABS, AIMAG, CMPLX, CONJG, MAX, MIN, REAL
      NOTRNA = LSAME( TRANA, 'N' )
      NOTRNB = LSAME( TRANB, 'N' )
      INFO = 0
      IF( .NOT.NOTRNA .AND. .NOT.LSAME( TRANA, 'C' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRNB .AND. .NOT.LSAME( TRANB, 'C' ) ) THEN
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
         CALL XERBLA( 'CTRSY2', -INFO )
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
      SMIN = MAX( SMLNUM, EPS*CLANGE( 'M', M, M, A, LDA, DUM ),
     $       EPS*CLANGE( 'M', N, N, B, LDB, DUM ) )
      SGN = ISGN
      IF( NOTRNA .AND. NOTRNB ) THEN
         DO 30 L = 1, N
            DO 20 K = M, 1, -1
               SUML = CDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,
     $                C( MIN( K+1, M ), L ), 1 )
               SUMR = CDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
               SCALOC = ONE
               A11 = A( K, K ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 10 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
   10             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
   20       CONTINUE
   30    CONTINUE
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
         DO 60 L = 1, N
            DO 50 K = 1, M
               SUML = CDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = CDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
               SCALOC = ONE
               A11 = CONJG( A( K, K ) ) + SGN*B( L, L )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
   50       CONTINUE
   60    CONTINUE
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
         DO 90 L = N, 1, -1
            DO 80 K = 1, M
               SUML = CDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = CDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,
     $                B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
               SCALOC = ONE
               A11 = CONJG( A( K, K )+SGN*B( L, L ) )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 70 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
   70             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
   80       CONTINUE
   90    CONTINUE
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
         DO 120 L = N, 1, -1
            DO 110 K = M, 1, -1
               SUML = CDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,
     $                C( MIN( K+1, M ), L ), 1 )
               SUMR = CDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,
     $                B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*CONJG( SUMR ) )
               SCALOC = ONE
               A11 = A( K, K ) + SGN*CONJG( B( L, L ) )
               DA11 = ABS( REAL( A11 ) ) + ABS( AIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( REAL( VEC ) ) + ABS( AIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = CLADIV( VEC*CMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 100 J = 1, N
                     CALL CSSCAL( M, SCALOC, C( 1, J ), 1 )
  100             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
  110       CONTINUE
  120    CONTINUE
      END IF
      RETURN
      END
