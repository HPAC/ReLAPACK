* ZTRSY2 solves the complex Sylvester matrix equation (unblocked algorithm)
*
* This routine is an exact copy of LAPACK's ctrsyl.f.
* It serves as an unblocked kernel in the recursive algorithms.
      SUBROUTINE ZTRSY2( TRANA, TRANB, ISGN, M, N, A, LDA, B, LDB, C,
     $                   LDC, SCALE, INFO )
      CHARACTER          TRANA, TRANB
      INTEGER            INFO, ISGN, LDA, LDB, LDC, M, N
      DOUBLE PRECISION   SCALE
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      LOGICAL            NOTRNA, NOTRNB
      INTEGER            J, K, L
      DOUBLE PRECISION   BIGNUM, DA11, DB, EPS, SCALOC, SGN, SMIN,
     $                   SMLNUM
      COMPLEX*16         A11, SUML, SUMR, VEC, X11
      DOUBLE PRECISION   DUM( 1 )
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANGE
      COMPLEX*16         ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           LSAME, DLAMCH, ZLANGE, ZDOTC, ZDOTU, ZLADIV
      EXTERNAL           DLABAD, XERBLA, ZDSCAL
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
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
         CALL XERBLA( 'ZTRSY2', -INFO )
         RETURN
      END IF
      SCALE = ONE
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
      EPS = DLAMCH( 'P' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SMLNUM*DBLE( M*N ) / EPS
      BIGNUM = ONE / SMLNUM
      SMIN = MAX( SMLNUM, EPS*ZLANGE( 'M', M, M, A, LDA, DUM ),
     $       EPS*ZLANGE( 'M', N, N, B, LDB, DUM ) )
      SGN = ISGN
      IF( NOTRNA .AND. NOTRNB ) THEN
         DO 30 L = 1, N
            DO 20 K = M, 1, -1
               SUML = ZDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,
     $                C( MIN( K+1, M ), L ), 1 )
               SUMR = ZDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
               SCALOC = ONE
               A11 = A( K, K ) + SGN*B( L, L )
               DA11 = ABS( DBLE( A11 ) ) + ABS( DIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( DBLE( VEC ) ) + ABS( DIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 10 J = 1, N
                     CALL ZDSCAL( M, SCALOC, C( 1, J ), 1 )
   10             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
   20       CONTINUE
   30    CONTINUE
      ELSE IF( .NOT.NOTRNA .AND. NOTRNB ) THEN
         DO 60 L = 1, N
            DO 50 K = 1, M
               SUML = ZDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = ZDOTU( L-1, C( K, 1 ), LDC, B( 1, L ), 1 )
               VEC = C( K, L ) - ( SUML+SGN*SUMR )
               SCALOC = ONE
               A11 = DCONJG( A( K, K ) ) + SGN*B( L, L )
               DA11 = ABS( DBLE( A11 ) ) + ABS( DIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( DBLE( VEC ) ) + ABS( DIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 40 J = 1, N
                     CALL ZDSCAL( M, SCALOC, C( 1, J ), 1 )
   40             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
   50       CONTINUE
   60    CONTINUE
      ELSE IF( .NOT.NOTRNA .AND. .NOT.NOTRNB ) THEN
         DO 90 L = N, 1, -1
            DO 80 K = 1, M
               SUML = ZDOTC( K-1, A( 1, K ), 1, C( 1, L ), 1 )
               SUMR = ZDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,
     $                B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*DCONJG( SUMR ) )
               SCALOC = ONE
               A11 = DCONJG( A( K, K )+SGN*B( L, L ) )
               DA11 = ABS( DBLE( A11 ) ) + ABS( DIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( DBLE( VEC ) ) + ABS( DIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 70 J = 1, N
                     CALL ZDSCAL( M, SCALOC, C( 1, J ), 1 )
   70             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
   80       CONTINUE
   90    CONTINUE
      ELSE IF( NOTRNA .AND. .NOT.NOTRNB ) THEN
         DO 120 L = N, 1, -1
            DO 110 K = M, 1, -1
               SUML = ZDOTU( M-K, A( K, MIN( K+1, M ) ), LDA,
     $                C( MIN( K+1, M ), L ), 1 )
               SUMR = ZDOTC( N-L, C( K, MIN( L+1, N ) ), LDC,
     $                B( L, MIN( L+1, N ) ), LDB )
               VEC = C( K, L ) - ( SUML+SGN*DCONJG( SUMR ) )
               SCALOC = ONE
               A11 = A( K, K ) + SGN*DCONJG( B( L, L ) )
               DA11 = ABS( DBLE( A11 ) ) + ABS( DIMAG( A11 ) )
               IF( DA11.LE.SMIN ) THEN
                  A11 = SMIN
                  DA11 = SMIN
                  INFO = 1
               END IF
               DB = ABS( DBLE( VEC ) ) + ABS( DIMAG( VEC ) )
               IF( DA11.LT.ONE .AND. DB.GT.ONE ) THEN
                  IF( DB.GT.BIGNUM*DA11 )
     $               SCALOC = ONE / DB
               END IF
               X11 = ZLADIV( VEC*DCMPLX( SCALOC ), A11 )
               IF( SCALOC.NE.ONE ) THEN
                  DO 100 J = 1, N
                     CALL ZDSCAL( M, SCALOC, C( 1, J ), 1 )
  100             CONTINUE
                  SCALE = SCALE*SCALOC
               END IF
               C( K, L ) = X11
  110       CONTINUE
  120    CONTINUE
      END IF
      RETURN
      END
