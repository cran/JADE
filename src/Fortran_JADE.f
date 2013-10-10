      SUBROUTINE FGALG(A, RN, B, IP, K, IFAULT, NF, EPSF)
C
      DOUBLE PRECISION A(IP, IP, K), H(IP, IP, K), B(IP, IP), RN(K)
      DOUBLE PRECISION ZERO, ONE, P8, P6, G1, G2, EPSF
      DATA ZERO /0.0/, ONE /1.0/, P8 /0.8/, P6 /0.6/, P1 /0.1/ 
C
C        CHECK INPUT PARAMETERS
C
      IFAULT = 1
      IF (K .LT. 1 .OR. K .GT. K) RETURN
      IFAULT = 2
      IF (IP .LT. 2 .OR. IP .GT. IP) RETURN
      IFAULT = 3
      DO 3 I = 1, K
      IF (RN(I) .LE. ZERO) RETURN
    3 CONTINUE
      DO 4 I = 1, K
      IFAULT = -I
      DO 4 L = 2, IP
      JEND = L - 1
      DO 4 J = 1, JEND
      IF (ABS(A(L, J, I) - A(J, L, I)) .GE. P1) RETURN
    4 CONTINUE
      IFAULT = 0
C
      DO 5 L = 1, IP
      DO 5 J = 1, IP
      B(L, J) = ZERO
      IF (L .EQ. J) B(L, J) = ONE
    5 CONTINUE
C
      DO 7 I = 1, K
      DO 7 L = 1, IP
      DO 7 J = 1, IP
    7 H(L, J, I) = A(L, J, I)
      CALL FALG(A, RN, B, IP, K, IFAULT, NF, EPSF)
      IF (IFAULT .NE. 0 .OR. NF .GT. 1) RETURN
C
C        IF F-ALGORITHM STOPPED AFTER ONLY 1 ITERATION, TRY
C        ANOTHER INITIAL APPROXIMATION TO THE ORTHOGONAL MATRIX B.
C
      DO 8 I = 1, K
      DO 8 L = 1, IP
      DO 8 J = 1, IP
    8 A(L, J, I) = H(L, J, I)
      IP1 = IP - 1
      DO 6 L = 1, IP1
      DO 6 J = 1, IP
      G1 = P8 * B(L, J) + P6 * B(L + 1, J)
      G2 = P8 * B(L + 1, J) - P6 * B(L, J)
      B(L, J) = G1
      B(L + 1, J) = G2
    6 CONTINUE
      CALL FALG(A, RN, B, IP, K, IFAULT, NF, EPSF)
      RETURN
      END
C
      SUBROUTINE FALG(A, RN, B, IP, K, IFAULT, NF, EPSF)
      DOUBLE PRECISION A(IP, IP, K), Q(2, 2), B(IP, IP), BOLD(IP, IP)
      DOUBLE PRECISION F(IP, IP, K), AUX(IP, IP), RN(K), T(2, 2, K)
      DOUBLE PRECISION EPS, EPSF, DIFF, ZERO, ONE, TWO
      DOUBLE PRECISION B1, B2, C, C2, R, S, T1, T2
      DATA ZERO /0.0/, ONE/1.0/, TWO/2.0/
C
      MAXF = NF
      NF = 0
      IFAULT = 0
      IP1 = IP - 1
    4 CONTINUE
C
C        COMPUTE MATRICES F(I), USING THE CURRENT B
C
      DO 2 I = 1, K
      DO 1 M = 1, IP
      DO 1 J = 1, IP
    1 AUX(M, J) = A(M, J, I)
      CALL MULT(AUX, B, IP, IP)
      DO 2 M = 1, IP
      DO 2 J = 1, IP
    2 F(M, J, I) = AUX(M, J)
C
C         STEP F1
C
      NF = NF + 1
      DO 5 M = 1, IP
      DO 5 J = 1, IP
    5 BOLD(M, J) = B(M, J)
C
C         STEP F2 (START SWEEP of F-ALGORITHM) 
C
      DO 8 M = 1, IP1
      JSTART = M + 1
C
C         STEP F21
C
      DO 8 J = JSTART, IP
      DO 6 I = 1, K
      T(1, 1, I) = F(M, M, I)   
      T(1, 2, I) = F(M, J, I)
      T(2, 1, I) = F(J, M, I)   
      T(2, 2, I) = F(J, J, I)
    6 CONTINUE
C
C         STEP F22
C
      CALL GALG(T, Q, RN, K, IFAULT)
      IF (IFAULT .NE. 0) RETURN
      C = Q(1, 1)
      S = Q(2, 1)
      R = S / (ONE + C)
C
C         STEP F23
C  
      T1 = S / C
      T2 = T1 * T1
      C2 = C * C
      DO 9 I = 1, K
      F(M, M, I) = C2 * (T(1, 1, I) + TWO * T1 * T(1, 2, I)
     *           + T2 * T(2, 2, I))  
      F(J, J, I) = C2 * (T2 * T(1, 1, I) - TWO * T1 * T(1, 2, I)
     *           + T(2, 2, I))
      F(M, J, I) = C2 * (T1 * (T(2, 2, I) - T(1, 1, I))
     *           + (ONE - T2) * T(1, 2, I))
      F(J, M, I) = F(M, J, I) 
      DO 9 L = 1, IP
      IF (L .EQ. J .OR. L. EQ. M) GOTO 9
      B1 = F(L, M, I)
      B2 = F(L, J, I)
      F(L, M, I) = B1 + S * (B2 - R * B1)
      F(L, J, I) = B2 - S * (B1 + R * B2)
      F(M, L, I) = F(L, M, I)
      F(J, L, I) = F(L, J, I)
    9 CONTINUE
C
C        STEP F24
C 
      DO 7 L = 1, IP
      B1 = B(L, M)
      B2 = B(L, J) 
      B(L, M) = B1 + S * (B2 - R * B1)
      B(L, J) = B2 - S * (B1 + R * B2)
    7 CONTINUE
    8 CONTINUE
C
C        FOR LARGE PROBLEMS, RE-ORTHOGONALIZE THE MATRIX B HERE, USING
C        THE MODIFIED GRAM-SCHMIDT ALGORITHM (SEE GOLUB & VAN LOAN,  
C        1983, PP. 151-152)
C
C        STEP F3
C
      EPS = ZERO 
      DO 10 M = 1, IP 
      DO 10 J = 1, IP
      DIFF = ABS(B(M, J) - BOLD(M, J))
      IF (DIFF .GT. EPS) EPS = DIFF
   10 CONTINUE
      IF(EPS .GE. EPSF .AND. NF .LT. MAXF) GOTO 4
      IF(EPS .GE. EPSF) IFAULT = 4
C
C         COPY MATRICES F(I) INTO A(I) BEFORE RETURNING TO FGALG
C  
      DO 11 I = 1, K
      DO 11 M = 1, IP
      DO 11 J = 1, IP
   11 A(M, J, I) = F(M, J, I)
      RETURN
      END
C
C
C
      SUBROUTINE GALG(T, Q, RN, K, IFAULT)
      DOUBLE PRECISION T(2, 2, K),  RN(K),  U(2, 2), DELTA(K, 2)
      DOUBLE PRECISION Q(2, 2), COEF(K)
      DOUBLE PRECISION EPSG, ZERO, ONE, TWO, C, S, C2, S2, CS, SOLD
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, EPSG /0.0001/
      MAXG = 5
C
C        STEP G0
C
      IF (K .EQ. 1) MAXG = 1
      C = ONE
      S = ZERO
C
C        STEP G1
C
      DO 1 NG = 1, MAXG
      SOLD = S
      C2 = C * C
      S2 = S * S
      CS = TWO * C * S
C
C        STEP G2 
C
      DO 2 I = 1, K
      DELTA(I, 1) = C2 * T(1, 1, I) + S2 * T(2, 2, I) + CS * T(1, 2, I)
      DELTA(I, 2) = C2 * T(2, 2, I) + S2 * T(1, 1, I) - CS * T(1, 2, I)
C      
C        CHECK FOR NON-POSITIVE DEFINITENESS  
C
      IF (DELTA(I, 1) .GT. ZERO .AND. DELTA(I, 2) .GT. ZERO) GOTO 4
      IFAULT = 5
      RETURN
    4 COEF(I) = RN(I) * (DELTA(I, 1) - DELTA(I, 2)) / (DELTA(I, 1) 
     *        * DELTA(I, 2))
    2 CONTINUE
      DO 3 M = 1, 2
      DO 3 J = 1, 2
      U(M, J) = ZERO
      DO 3 I = 1, K
      U(M, J) = U(M, J) + COEF(I) * T(M, J, I)
    3 CONTINUE
C
C        STEP G3
C
      CALL EIGVEC(U, Q)
      C = Q(1, 1)
      S = Q(2, 1)
C
C        STEP G4       
C
      IF (ABS(S - SOLD) .LT. EPSG) RETURN
    1 CONTINUE
      RETURN
      END
C
      SUBROUTINE EIGVEC(U, Q)
      DOUBLE PRECISION U(2, 2), Q(2, 2)
      DOUBLE PRECISION ZERO, ONE, TWO, FOUR, RATIO, DISCR, C, S
      DOUBLE PRECISION T1, T2, T, EP
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/, FOUR /4.0/, EP /1.E-10/
C
      C = ONE
      S = ZERO
      IF (ABS(U(1, 2)) .LT. EP) RETURN
      RATIO = (U(2, 2) - U(1, 1)) / U(1, 2)
      DISCR = SQRT(RATIO * RATIO + FOUR)
      T1 = (RATIO + DISCR) / TWO
      T2 = (RATIO - DISCR) / TWO
      T = T1
      IF(ABS(T1) .GT. ABS(T2)) T = T2
      C = ONE / SQRT(ONE + T * T)
      S = T * C
      Q(1, 1) = C
      Q(2, 1) = S
      Q(1, 2) = -S
      Q(2, 2) = C  
      RETURN
      END
C
      SUBROUTINE MULT(A, B, IP, IQ)
      DOUBLE PRECISION A(IP, IP), B(IP, IP), H(IP, IP)
      DOUBLE PRECISION ZERO
      DATA ZERO /0.0/
C
      DO 1 I = 1, IQ
      DO 1 J = 1, IP
      H(I, J) = ZERO
      DO 1 L = 1, IP
      H(I, J) = H(I, J) + B(L, I) * A(L, J)
    1 CONTINUE
      DO 2 I = 1, IQ
      DO 2 J = 1, I
      A(I, J) = ZERO
      DO 3 L = 1, IP
    3 A(I, J) = A(I, J) + H(I, L) * B(L, J)
      A(J, I) = A(I, J)
    2 CONTINUE
      RETURN
      END

C
C--------------------------------------------------------------------------------
C
      SUBROUTINE FRJD(A, RN, B, IP, K, IFAULT, NF, EPSF)
      INTEGER I, IFAULT, IP, IP1, J, JSTART, K, L, M, NF, MAXF, CC
      DOUBLE PRECISION A(IP, IP, K), AUX(IP, IP), A1, A2, A3, A4 
      DOUBLE PRECISION B(IP, IP), B1, B2, BOLD(IP, IP)
      DOUBLE PRECISION C, Q(2, 2), RN(K), S, T(K, 2, 2)
      DOUBLE PRECISION ZERO, ONE, TWO, EPSF
      DATA ZERO /0.0/, ONE /1.0/, TWO /2.0/
C
      MAXF = NF
      NF = 0
      IFAULT = 0
      IP1 = IP - 1
C
C        INITIAL MULTIPLICATION
C
      DO 2 I = 1, K
      DO 1 L = 1, IP
      DO 1 J = 1, IP
    1 AUX(L, J) = A(L, J, I)
      CALL MULTR(AUX, B, IP, IP)
      DO 2 L = 1, IP
      DO 2 J = 1, IP
    2 A(L, J, I) = AUX(L, J)
C
    4 NF = NF + 1
      CC = 0
      DO 5 L = 1, IP
      DO 5 J = 1, IP
    5 BOLD(L, J) = B(L, J)
      DO 8 L = 1, IP1
      JSTART = L + 1
      DO 8 J = JSTART, IP
      DO 6 I = 1, K
      T(I, 1, 1) = A(L, L, I)
      T(I, 1, 2) = A(L, J, I)
      T(I, 2, 1) = A(J, L, I)
      T(I, 2, 2) = A(J, J, I)
    6 CONTINUE
C
C        GET ROTATION SINE AND COSINE
C
      CALL GALGR(T, RN, K, Q)
C
C        ROTATE B
C
      C = Q(1, 1)
      S = Q(1, 2)
      IF (abs(S) .GT. EPSF) CC = 1 
      DO 7 M = 1, IP
      B1 = B(L, M)
      B2 = B(J, M)
      B(L, M) = C*B1 + S*B2
      B(J, M) = C*B2 - S*B1
    7 CONTINUE
C
C        UPDATE THE A MATRICES
C
      DO 12 I = 1, K
      DO 10 M = 1, IP
      A1 = A(L, M, I)
      A2 = A(J, M, I)
      A(L, M, I) = C*A1 + S*A2
      A(J, M, I) = C*A2 - S*A1
   10 CONTINUE
      DO 11 M = 1, IP
      A3 = A(M, L, I)
      A4 = A(M, J, I)
      A(M, L, I) = C*A3 + S*A4
      A(M, J, I) = C*A4 - S*A3
   11 CONTINUE
   12 CONTINUE 
    8 CONTINUE
C
C        CHECK FOR CONVERGENCE
C
      IF (CC .NE. 0 .AND. NF .LT. MAXF) GOTO 4
      IF (CC .NE. 0) IFAULT = 4
      RETURN
      END
c
c---------------------------------------------------------------------
c
      SUBROUTINE GALGR (T, RN, K, Q)
      DOUBLE PRECISION C, CO, SI, U, WT, ALPHA, BETA, THETA
      DOUBLE PRECISION Q(2,2), RN(K), T(K,2,2), ZERO, TWO, FOUR 
C
      DATA ZERO/0.0/, TWO/2.0/, FOUR/4.0/
C
      ALPHA = ZERO
      BETA = ZERO
C 
      DO 1  I=1, K
      U  = T(I,1,1)-T(I,2,2)
      C  = T(I,1,2)
      WT = RN(I)
      ALPHA = ALPHA + WT*(U*U-FOUR*C*C)
      BETA = BETA + WT*FOUR*U*C
    1 CONTINUE

      THETA = ATAN(BETA/(ALPHA+SQRT(ALPHA*ALPHA+BETA*BETA)))/TWO

      SI     = SIN(THETA)
      CO     = COS(THETA)
      Q(1,1) = CO
      Q(1,2) = SI
      Q(2,1) = -SI
      Q(2,2) = CO
      RETURN
      END
C   
      SUBROUTINE AT (A, B, RES)
      DOUBLE PRECISION A, B, RES
C
      RES = ATAN(B/(A+SQRT(A*A+B*B)))      
      RETURN
      END



      SUBROUTINE MULTR(A, B, IP, IQ)
      DOUBLE PRECISION A(IP, IP), B(IP, IP), H(IP, IP)
      DOUBLE PRECISION ZERO
      DATA ZERO /0.0/
C
      DO 1 I = 1, IQ
      DO 1 J = 1, IP
      H(I, J) = ZERO
      DO 1 L = 1, IP
      H(I, J) = H(I, J) + B(I, L) * A(L, J)
    1 CONTINUE
      DO 2 I = 1, IQ
      DO 2 J = 1, I
      A(I, J) = ZERO
      DO 3 L = 1, IP
    3 A(I, J) = A(I, J) + H(I, L) * B(J, L)
      A(J, I) = A(I, J)
    2 CONTINUE
      RETURN
      END



