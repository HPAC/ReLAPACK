* DSYTRF_ROOK_REC2 computes a partial factorization of a real symmetric matrix using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
*
* This routine is a minor modification of LAPACK's dlasyf_rook.
* It serves as an unblocked kernel in the recursive algorithms.
* The blocked BLAS Level 3 updates were removed and moved to the
* recursive algorithm.
      SUBROUTINE DSYTRF_ROOK_REC2( UPLO, N, NB, KB, A, LDA, IPIV, W,
     $                             LDW, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, KB, LDA, LDW, N, NB
      INTEGER            IPIV( * )
      DOUBLE PRECISION   A( LDA, * ), W( LDW, * )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   EIGHT, SEVTEN
      PARAMETER          ( EIGHT = 8.0D+0, SEVTEN = 17.0D+0 )
      LOGICAL            DONE
      INTEGER            IMAX, ITEMP, J, JJ, JMAX, JP1, JP2, K, KK,
     $                   KW, KKW, KP, KSTEP, P, II
      DOUBLE PRECISION   ABSAKK, ALPHA, COLMAX, D11, D12, D21, D22,
     $                   DTEMP, R1, ROWMAX, T, SFMIN
      LOGICAL            LSAME
      INTEGER            IDAMAX
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, IDAMAX, DLAMCH
      EXTERNAL           DCOPY, DGEMM, DGEMV, DSCAL, DSWAP
      INTRINSIC          ABS, MAX, MIN, SQRT
      INFO = 0
      ALPHA = ( ONE+SQRT( SEVTEN ) ) / EIGHT
      SFMIN = DLAMCH( 'S' )
      IF( LSAME( UPLO, 'U' ) ) THEN
         K = N
   10    CONTINUE
         KW = NB + K - N
         IF( ( K.LE.N-NB+1 .AND. NB.LT.N ) .OR. K.LT.1 )
     $      GO TO 30
         KSTEP = 1
         P = K
         CALL DCOPY( K, A( 1, K ), 1, W( 1, KW ), 1 )
         IF( K.LT.N )
     $      CALL DGEMV( 'No transpose', K, N-K, -ONE, A( 1, K+1 ),
     $                  LDA, W( K, KW+1 ), LDW, ONE, W( 1, KW ), 1 )
         ABSAKK = ABS( W( K, KW ) )
         IF( K.GT.1 ) THEN
            IMAX = IDAMAX( K-1, W( 1, KW ), 1 )
            COLMAX = ABS( W( IMAX, KW ) )
         ELSE
            COLMAX = ZERO
         END IF
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
            CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
         ELSE
            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN
               KP = K
            ELSE
               DONE = .FALSE.
   12          CONTINUE
                  CALL DCOPY( IMAX, A( 1, IMAX ), 1, W( 1, KW-1 ), 1 )
                  CALL DCOPY( K-IMAX, A( IMAX, IMAX+1 ), LDA,
     $                        W( IMAX+1, KW-1 ), 1 )
                  IF( K.LT.N )
     $               CALL DGEMV( 'No transpose', K, N-K, -ONE,
     $                           A( 1, K+1 ), LDA, W( IMAX, KW+1 ), LDW,
     $                           ONE, W( 1, KW-1 ), 1 )
                  IF( IMAX.NE.K ) THEN
                     JMAX = IMAX + IDAMAX( K-IMAX, W( IMAX+1, KW-1 ),
     $                                     1 )
                     ROWMAX = ABS( W( JMAX, KW-1 ) )
                  ELSE
                     ROWMAX = ZERO
                  END IF
                  IF( IMAX.GT.1 ) THEN
                     ITEMP = IDAMAX( IMAX-1, W( 1, KW-1 ), 1 )
                     DTEMP = ABS( W( ITEMP, KW-1 ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF
                  IF( .NOT.(ABS( W( IMAX, KW-1 ) ).LT.ALPHA*ROWMAX ) )
     $            THEN
                     KP = IMAX
                     CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
                     DONE = .TRUE.
                  ELSE IF( ( P.EQ.JMAX ) .OR. ( ROWMAX.LE.COLMAX ) )
     $            THEN
                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  ELSE
                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                     CALL DCOPY( K, W( 1, KW-1 ), 1, W( 1, KW ), 1 )
                  END IF
               IF( .NOT. DONE ) GOTO 12
            END IF
            KK = K - KSTEP + 1
            KKW = NB + KK - N
            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN
               CALL DCOPY( K-P, A( P+1, K ), 1, A( P, P+1 ), LDA )
               CALL DCOPY( P, A( 1, K ), 1, A( 1, P ), 1 )
               CALL DSWAP( N-K+1, A( K, K ), LDA, A( P, K ), LDA )
               CALL DSWAP( N-KK+1, W( K, KKW ), LDW, W( P, KKW ), LDW )
            END IF
            IF( KP.NE.KK ) THEN
               A( KP, K ) = A( KK, K )
               CALL DCOPY( K-1-KP, A( KP+1, KK ), 1, A( KP, KP+1 ),
     $                     LDA )
               CALL DCOPY( KP, A( 1, KK ), 1, A( 1, KP ), 1 )
               CALL DSWAP( N-KK+1, A( KK, KK ), LDA, A( KP, KK ), LDA )
               CALL DSWAP( N-KK+1, W( KK, KKW ), LDW, W( KP, KKW ),
     $                     LDW )
            END IF
            IF( KSTEP.EQ.1 ) THEN
               CALL DCOPY( K, W( 1, KW ), 1, A( 1, K ), 1 )
               IF( K.GT.1 ) THEN
                  IF( ABS( A( K, K ) ).GE.SFMIN ) THEN
                     R1 = ONE / A( K, K )
                     CALL DSCAL( K-1, R1, A( 1, K ), 1 )
                  ELSE IF( A( K, K ).NE.ZERO ) THEN
                     DO 14 II = 1, K - 1
                        A( II, K ) = A( II, K ) / A( K, K )
   14                CONTINUE
                  END IF
               END IF
            ELSE
               IF( K.GT.2 ) THEN
                  D12 = W( K-1, KW )
                  D11 = W( K, KW ) / D12
                  D22 = W( K-1, KW-1 ) / D12
                  T = ONE / ( D11*D22-ONE )
                  DO 20 J = 1, K - 2
                     A( J, K-1 ) = T*( (D11*W( J, KW-1 )-W( J, KW ) ) /
     $                             D12 )
                     A( J, K ) = T*( ( D22*W( J, KW )-W( J, KW-1 ) ) /
     $                           D12 )
   20             CONTINUE
               END IF
               A( K-1, K-1 ) = W( K-1, KW-1 )
               A( K-1, K ) = W( K-1, KW )
               A( K, K ) = W( K, KW )
            END IF
         END IF
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -P
            IPIV( K-1 ) = -KP
         END IF
         K = K - KSTEP
         GO TO 10
   30    CONTINUE
         J = K + 1
   60    CONTINUE
            KSTEP = 1
            JP1 = 1
            JJ = J
            JP2 = IPIV( J )
            IF( JP2.LT.0 ) THEN
               JP2 = -JP2
               J = J + 1
               JP1 = -IPIV( J )
               KSTEP = 2
            END IF
            J = J + 1
            IF( JP2.NE.JJ .AND. J.LE.N )
     $         CALL DSWAP( N-J+1, A( JP2, J ), LDA, A( JJ, J ), LDA )
            JJ = J - 1
            IF( JP1.NE.JJ .AND. KSTEP.EQ.2 )
     $         CALL DSWAP( N-J+1, A( JP1, J ), LDA, A( JJ, J ), LDA )
         IF( J.LE.N )
     $      GO TO 60
         KB = N - K
      ELSE
         K = 1
   70   CONTINUE
         IF( ( K.GE.NB .AND. NB.LT.N ) .OR. K.GT.N )
     $      GO TO 90
         KSTEP = 1
         P = K
         CALL DCOPY( N-K+1, A( K, K ), 1, W( K, K ), 1 )
         IF( K.GT.1 )
     $      CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE, A( K, 1 ),
     $                  LDA, W( K, 1 ), LDW, ONE, W( K, K ), 1 )
         ABSAKK = ABS( W( K, K ) )
         IF( K.LT.N ) THEN
            IMAX = K + IDAMAX( N-K, W( K+1, K ), 1 )
            COLMAX = ABS( W( IMAX, K ) )
         ELSE
            COLMAX = ZERO
         END IF
         IF( MAX( ABSAKK, COLMAX ).EQ.ZERO ) THEN
            IF( INFO.EQ.0 )
     $         INFO = K
            KP = K
            CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
         ELSE
            IF( .NOT.( ABSAKK.LT.ALPHA*COLMAX ) ) THEN
               KP = K
            ELSE
               DONE = .FALSE.
   72          CONTINUE
                  CALL DCOPY( IMAX-K, A( IMAX, K ), LDA, W( K, K+1 ), 1)
                  CALL DCOPY( N-IMAX+1, A( IMAX, IMAX ), 1,
     $                        W( IMAX, K+1 ), 1 )
                  IF( K.GT.1 )
     $               CALL DGEMV( 'No transpose', N-K+1, K-1, -ONE,
     $                           A( K, 1 ), LDA, W( IMAX, 1 ), LDW,
     $                           ONE, W( K, K+1 ), 1 )
                  IF( IMAX.NE.K ) THEN
                     JMAX = K - 1 + IDAMAX( IMAX-K, W( K, K+1 ), 1 )
                     ROWMAX = ABS( W( JMAX, K+1 ) )
                  ELSE
                     ROWMAX = ZERO
                  END IF
                  IF( IMAX.LT.N ) THEN
                     ITEMP = IMAX + IDAMAX( N-IMAX, W( IMAX+1, K+1 ), 1)
                     DTEMP = ABS( W( ITEMP, K+1 ) )
                     IF( DTEMP.GT.ROWMAX ) THEN
                        ROWMAX = DTEMP
                        JMAX = ITEMP
                     END IF
                  END IF
                  IF( .NOT.( ABS( W( IMAX, K+1 ) ).LT.ALPHA*ROWMAX ) )
     $            THEN
                     KP = IMAX
                     CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
                     DONE = .TRUE.
                  ELSE IF( ( P.EQ.JMAX ) .OR. ( ROWMAX.LE.COLMAX ) )
     $            THEN
                     KP = IMAX
                     KSTEP = 2
                     DONE = .TRUE.
                  ELSE
                     P = IMAX
                     COLMAX = ROWMAX
                     IMAX = JMAX
                     CALL DCOPY( N-K+1, W( K, K+1 ), 1, W( K, K ), 1 )
                  END IF
               IF( .NOT. DONE ) GOTO 72
            END IF
            KK = K + KSTEP - 1
            IF( ( KSTEP.EQ.2 ) .AND. ( P.NE.K ) ) THEN
               CALL DCOPY( P-K, A( K, K ), 1, A( P, K ), LDA )
               CALL DCOPY( N-P+1, A( P, K ), 1, A( P, P ), 1 )
               CALL DSWAP( K, A( K, 1 ), LDA, A( P, 1 ), LDA )
               CALL DSWAP( KK, W( K, 1 ), LDW, W( P, 1 ), LDW )
            END IF
            IF( KP.NE.KK ) THEN
               A( KP, K ) = A( KK, K )
               CALL DCOPY( KP-K-1, A( K+1, KK ), 1, A( KP, K+1 ), LDA )
               CALL DCOPY( N-KP+1, A( KP, KK ), 1, A( KP, KP ), 1 )
               CALL DSWAP( KK, A( KK, 1 ), LDA, A( KP, 1 ), LDA )
               CALL DSWAP( KK, W( KK, 1 ), LDW, W( KP, 1 ), LDW )
            END IF
            IF( KSTEP.EQ.1 ) THEN
               CALL DCOPY( N-K+1, W( K, K ), 1, A( K, K ), 1 )
               IF( K.LT.N ) THEN
                  IF( ABS( A( K, K ) ).GE.SFMIN ) THEN
                     R1 = ONE / A( K, K )
                     CALL DSCAL( N-K, R1, A( K+1, K ), 1 )
                  ELSE IF( A( K, K ).NE.ZERO ) THEN
                     DO 74 II = K + 1, N
                        A( II, K ) = A( II, K ) / A( K, K )
   74                CONTINUE
                  END IF
               END IF
            ELSE
               IF( K.LT.N-1 ) THEN
                  D21 = W( K+1, K )
                  D11 = W( K+1, K+1 ) / D21
                  D22 = W( K, K ) / D21
                  T = ONE / ( D11*D22-ONE )
                  DO 80 J = K + 2, N
                     A( J, K ) = T*( ( D11*W( J, K )-W( J, K+1 ) ) /
     $                           D21 )
                     A( J, K+1 ) = T*( ( D22*W( J, K+1 )-W( J, K ) ) /
     $                             D21 )
   80             CONTINUE
               END IF
               A( K, K ) = W( K, K )
               A( K+1, K ) = W( K+1, K )
               A( K+1, K+1 ) = W( K+1, K+1 )
            END IF
         END IF
         IF( KSTEP.EQ.1 ) THEN
            IPIV( K ) = KP
         ELSE
            IPIV( K ) = -P
            IPIV( K+1 ) = -KP
         END IF
         K = K + KSTEP
         GO TO 70
   90    CONTINUE
         J = K - 1
  120    CONTINUE
            KSTEP = 1
            JP1 = 1
            JJ = J
            JP2 = IPIV( J )
            IF( JP2.LT.0 ) THEN
               JP2 = -JP2
               J = J - 1
               JP1 = -IPIV( J )
               KSTEP = 2
            END IF
            J = J - 1
            IF( JP2.NE.JJ .AND. J.GE.1 )
     $         CALL DSWAP( J, A( JP2, 1 ), LDA, A( JJ, 1 ), LDA )
            JJ = J + 1
            IF( JP1.NE.JJ .AND. KSTEP.EQ.2 )
     $         CALL DSWAP( J, A( JP1, 1 ), LDA, A( JJ, 1 ), LDA )
         IF( J.GE.1 )
     $      GO TO 120
         KB = K - 1
      END IF
      RETURN
      END
