** ZHETRF_REC2 computes a partial factorization of a complex Hermitian indefinite matrix using the Bunch-Kaufman diagonal pivoting method
*
* This routine is a minor modification of LAPACK's zlahef.
* It serves as an unblocked kernel in the recursive algorithms.
* The blocked BLAS Level 3 updates were removed and moved to the
* recursive algorithm.
      SUBROUTINE ZHETRF_REC2( UPLO, N, NB, KB, A, LDA, IPIV, W, LDW,
     $                        INFO )
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), W( LDW, * )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      INTEGER            IMAX, J, JJ, JMAX, JP, K, KK, KKW, KP,
     $                   KSTEP, KW
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, R1, ROWMAX, T
      COMPLEX*16         D11, D21, D22, Z
      LOGICAL            LSAME
      INTEGER            IZAMAX
      EXTERNAL           LSAME, IZAMAX
      EXTERNAL           ZCOPY, ZDSCAL, ZGEMM, ZGEMV, ZLACGV, ZSWAP
      INTRINSIC          ABS, DBLE, DCONJG, DIMAG, MAX, MIN, SQRT
      DOUBLE PRECISION   CABS1
      CABS1( Z ) = ABS( DBLE( Z ) ) + ABS( DIMAG( Z ) )
      INFO = 0
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
      IF( LSAME( UPLO, 'U' ) ) THEN
         K = N
   10    CONTINUE
         KW = NB + K - N
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
         KSTEP = 1
         CALL ZCOPY( K-1, A( 1, K ), 1, W( 1, KW ), 1 )
         W( K, KW ) = DBLE( A( K, K ) )
         IF( K.LT.N ) THEN
            CALL ZGEMV( 'No transpose', K, N-K, -CONE, A( 1, K+1 ), LDA,
     $                  W( K, KW+1 ), LDW, CONE, W( 1, KW ), 1 )
            W( K, KW ) = DBLE( W( K, KW ) )
         END IF
         ABSAKK = ABS( DBLE( W( K, KW ) ) )
         IF( K.GT.1 ) THEN
            IMAX = IZAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = CABS1( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
            A( K, K ) = DBLE( A( K, K ) )
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
               KP = K
            ELSE
               CALL ZCOPY( IMAX-1, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
               W( IMAX, KW-1 ) = DBLE( A( IMAX, IMAX ) )
               CALL ZCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                     W( IMAX+1, KW-1 ), 1 )
               CALL ZLACGV( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               IF( K.LT.N ) THEN
                  CALL ZGEMV( 'No transpose', K, N-K, -CONE,
     $                        A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW,
     $                        CONE, W( 1, KW-1 ), 1 )
                  W( IMAX, KW-1 ) = DBLE( W( IMAX, KW-1 ) )
               END IF
               JMAX = IMAX + IZAMAX( K-IMAX, W( IMAX+1, KW-1 ), 1 )
               ROWMAX = CABS1( W( JMAX, KW-1 ) )
               IF( IMAX.GT.1 ) THEN
                  JMAX = IZAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, KW-1 ) ) )
               END IF
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                  KP = K
               ELSE IF( ABS( DBLE( W( IMAX, KW-1 ) ) ).GE.ALPHA*ROWMAX )
     $                   THEN
                  KP = IMAX
                  CALL ZCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
               ELSE
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
            KK = K - KSTEP + 1
            KKW = NB + KK - N
            IF( KP.NE.KK ) THEN
               A( KP, KP ) = DBLE( A( KK, KK ) )
               CALL ZCOPY( KK-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               CALL ZLACGV( KK-1-KP, A( KP, KP+1 ), LDA )
               IF( KP.GT.1 )
     $            CALL ZCOPY( KP-1, A( 1, KK ), 1, A( 1, KP ), 1 )
               IF( K.LT.N )
     $            CALL ZSWAP( N-K, A( KK, K+1 ), LDA, A( KP, K+1 ),
     $                        LDA )
               CALL ZSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
            IF( KSTEP.EQ.1 ) THEN
               CALL ZCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               IF( K.GT.1 ) THEN
                  R1 = ONE / DBLE( A( K, K ) )
                  CALL ZDSCAL( K-1, R1, A( 1, K ), 1 )
                  CALL ZLACGV( K-1, W( 1, KW ), 1 )
               END IF
            ELSE
               IF( K.GT.2 ) THEN
                  D21 = W( K-1, KW )
                  D11 = W( K, KW ) / DCONJG( D21 )
                  D22 = W( K-1, KW-1 ) / D21
                  T = ONE / ( DBLE( D11*D22 )-ONE )
                  D21 = T / D21
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = D21*( D11*W( J, KW-1 )-W( J, KW ) )
                     A( J, K ) = DCONJG( D21 )*
     $                           ( D22*W( J, KW )-W( J, KW-1 ) )
   20             CONTINUE
               END IF
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
               CALL ZLACGV( K-1, W( 1, KW ), 1 )
               CALL ZLACGV( K-2, W( 1, KW-1 ), 1 )
            END IF
         END IF
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K-1 ) = -KP
         END IF
         K = K - KSTEP
         GO TO 10
   30    CONTINUE
         J = K + 1
   60    CONTINUE
            JJ = J
            JP = IPIV( J )
            IF( JP.LT.0 ) THEN
               JP = -JP
               J = J + 1
            END IF
            J = J + 1
            IF( JP.NE.JJ .AND. J.LE.N )
     $         CALL ZSWAP( N-J+1, A( JP, J ), LDA, A( JJ, J ), LDA )
         IF( J.LT.N )
     $      GO TO 60
         KB = N - K
      ELSE
         K = 1
   70    CONTINUE
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
         KSTEP = 1
         W( K, K ) = DBLE( A( K, K ) )
         IF( K.LT.N )
     $      CALL ZCOPY( N-K, A( K+1, K ), 1, W( K+1, K ), 1 )
         CALL ZGEMV( 'No transpose', N-K+1, K-1, -CONE, A( K, 1 ), LDA,
     $               W( K, 1 ), LDW, CONE, W( K, K ), 1 )
         W( K, K ) = DBLE( W( K, K ) )
         ABSAKK = ABS( DBLE( W( K, K ) ) )
         IF( K.LT.N ) THEN
            IMAX = K + IZAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = CABS1( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
            A( K, K ) = DBLE( A( K, K ) )
         ELSE
            IF( ABSAKK.GE.ALPHA*COLMAX ) THEN
               KP = K
            ELSE
               CALL ZCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1 )
               CALL ZLACGV( IMAX-K, W( K, K+1 ), 1 )
               W( IMAX, K+1 ) = DBLE( A( IMAX, IMAX ) )
               IF( IMAX.LT.N )
     $            CALL ZCOPY( N-IMAX, A( IMAX+1, IMAX ), 1,
     $                        W( IMAX+1, K+1 ), 1 )
               CALL ZGEMV( 'No transpose', N-K+1, K-1, -CONE, A( K, 1 ),
     $                     LDA, W( IMAX, 1 ), LDW, CONE, W( K, K+1 ),
     $                     1 )
               W( IMAX, K+1 ) = DBLE( W( IMAX, K+1 ) )
               JMAX = K - 1 + IZAMAX( IMAX-K, W( K, K+1 ), 1 )
               ROWMAX = CABS1( W( JMAX, K+1 ) )
               IF( IMAX.LT.N ) THEN
                  JMAX = IMAX + IZAMAX( N-IMAX, W( IMAX+1, K+1 ), 1 )
                  ROWMAX = MAX( ROWMAX, CABS1( W( JMAX, K+1 ) ) )
               END IF
               IF( ABSAKK.GE.ALPHA*COLMAX*( COLMAX / ROWMAX ) ) THEN
                  KP = K
               ELSE IF( ABS( DBLE( W( IMAX, K+1 ) ) ).GE.ALPHA*ROWMAX )
     $                   THEN
                  KP = IMAX
                  CALL ZCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
               ELSE
                  KP = IMAX
                  KSTEP = 2
               END IF
            END IF
            KK = K + KSTEP - 1
            IF( KP.NE.KK ) THEN
               A( KP, KP ) = DBLE( A( KK, KK ) )
               CALL ZCOPY( KP-KK-1, A( KK+1, KK ), 1, A( KP, KK+1 ),
     $                     LDA )
               CALL ZLACGV( KP-KK-1, A( KP, KK+1 ), LDA )
               IF( KP.LT.N )
     $            CALL ZCOPY( N-KP, A( KP+1, KK ), 1, A( KP+1, KP ), 1 )
               IF( K.GT.1 )
     $            CALL ZSWAP( K-1, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL ZSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
            IF( KSTEP.EQ.1 ) THEN
               CALL ZCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  R1 = ONE / DBLE( A( K, K ) )
                  CALL ZDSCAL( N-K, R1, A( K+1, K ), 1 )
                  CALL ZLACGV( N-K, W( K+1, K ), 1 )
               END IF
            ELSE
               IF( K.LT.N-1 ) THEN
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / DCONJG( D21 )
                  T = ONE / ( DBLE( D11*D22 )-ONE )
                  D21 = T / D21
                  DO 80 J = K + 2, N
                     A( J, K ) = DCONJG( D21 )*
     $                           ( D11*W( J, K )-W( J, K+1 ) )
                     A( J, K+1 ) = D21*( D22*W( J, K+1 )-W( J, K ) )
   80             CONTINUE
               END IF
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
               CALL ZLACGV( N-K, W( K+1, K ), 1 )
               CALL ZLACGV( N-K-1, W( K+2, K+1 ), 1 )
            END IF
         END IF
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -KP
            IPIV( K+1 ) = -KP
         END IF
         K = K + KSTEP
         GO TO 70
   90    CONTINUE
         J = K - 1
  120    CONTINUE
            JJ = J
            JP = IPIV( J )
            IF( JP.LT.0 ) THEN
               JP = -JP
               J = J - 1
            END IF
            J = J - 1
            IF( JP.NE.JJ .AND. J.GE.1 )
     $         CALL ZSWAP( J, A( JP, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GT.1 )
     $      GO TO 120
         KB = K - 1
      END IF
      RETURN
      END
