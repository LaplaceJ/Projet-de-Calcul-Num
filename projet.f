      PROGRAM PROJET
* La matrice A et le vecteur b sont surdimensionnés
      INTEGER MMAX, NMAX
      PARAMETER (MMAX = 100)
      PARAMETER (NMAX = 10)
      DOUBLE PRECISION A(MMAX,NMAX), B(MMAX)
* Les dimensions mathématiques de A (m x n) et b (m) sont lues sur l'entrée standard
      INTEGER M, N
*
* Lecture de A et b.
      CALL DREAD_MPL (M, N, A, MMAX)
      CALL DREAD_MPL (M, 1, B, MMAX)
* Affichage de A et b
      CALL DPRINT_MPL ('A', M, N, A, MMAX)
      CALL DPRINT_MPL ('b', M, 1, B, MMAX)
* Résolution par l'algorithme de Cholesky
      CALL CHOLESKY_SOLVE (M, N, A, MMAX, B)
* Résolution par la factorisation QR
      CALL FACT_QR_SOLVE (M, N, A, MMAX, B)
      END PROGRAM

      SUBROUTINE CHOLESKY_SOLVE (M, N, A, LDA, B)
* Les paramètres formels (LDA = Leading Dimension of A)
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,*), B(*)
* Variables locales
      DOUBLE PRECISION ATA(N,N), ATB(N)
      
*
      WRITE (*,*) 'Résolution par la méthode du Cdt Cholesky'
* ATA = Transpose (A) . A   (DGEMM est une BLAS pour le pdt de matrices)
      CALL DGEMM ('Transpose', 'No Transpose', 
     $             N, N, M, 1D0, A, LDA, A, LDA, 0D0, ATA, N)
      CALL DPRINT_MPL ('A**T . A', N, N, ATA, N)
* ATB = Transpose (A) . b   (DGEMV est une BLAS pour le pdt matrice . vecteur)
      CALL DGEMV ('Transpose', M, N, 1D0, A, LDA, B, 1, 0D0, ATB, 1)
      CALL DPRINT_MPL ('A**T . b', N, 1, ATB, N)
* Factorisation de Cholesky : ATA = L . Transpose (L)
* DPOTRF est une fonction LAPACK pour la méthode de Cholesky
      !CALL DPOTRF2 (0, N, ATA, N)
      CALL DPOTRF2 (0, N, ATA, N)
      !IF (INFO .NE. 0) STOP 'Erreur DPOTRF'
      CALL DPRINT_MPL ('Après DPOTRF, A', N, N, ATA, N)
* Résolution L . y = ATB avec résultat dans ATB.
* DTRSV est une BLAS pour la substitution avant/arrière
      CALL DTRSV2 (0,0,N, ATA, N, ATB)
* Résolution Transpose (L) . x = y (résultat dans ATB)
      CALL DTRSV2 (1,1,N, ATA, N, ATB)
      CALL DPRINT_MPL ('solution', N, 1, ATB, N)
      END SUBROUTINE

      SUBROUTINE FACT_QR_SOLVE (M, N, A, LDA, B)
* Les paramètres formels (LDA = Leading Dimension of A)
      INTEGER M, N, LDA
      DOUBLE PRECISION A(LDA,*), B(*) ,  TMP(LDA ,N)
* Variables locales
      DOUBLE PRECISION  QTB(M)
     
*
      WRITE (*,*) 'Résolution par la factorisation QR'
* Factorisation A = Q . R
* DGEQRF2 est une fonction pour la factorisation QR (Householder) sur A 
      CALL DGEQRF2 (M, N, A, LDA , TMP,LDA ,N )
      
      CALL DPRINT_MPL ('Après DGEQRF, A', M, N, A, LDA)
* QTB = b   (DCOPY est une BLAS qui copie un vecteur dans un autre)
      
      CALL DCOPY (M, B, 1, QTB, 1)
* QTB = Transpose (Q) . b
* DORMQR est une fonction LAPACK pour multiplier par Q après DGEQRF sur b 
      
      CALL DORMQR2 ( LDA , M , N , TMP , QTB, M)
      
      CALL DPRINT_MPL ('QTB', N, 1, QTB, N)
* Résolution R . x = QTB (résultat dans QTB)
      CALL DTRSV2 (1,0,N, A, LDA, QTB)
      CALL DPRINT_MPL ('solution', N, 1, QTB, N)
      END SUBROUTINE

