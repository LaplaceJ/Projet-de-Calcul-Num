* Le répertoire cnum est supposé placé dans le répertoire principal
* $HOME/cnum/gfor factorisation_LU_real.f
* cat factorisation_LU.dat
* ./a.out < factorisation_LU.dat

      PROGRAM FAC_LU
      INTEGER DMAX
      PARAMETER (DMAX=10)
      REAL A(DMAX,DMAX), B(DMAX)
      INTEGER IPIV(DMAX), INFO, M, N
*
      WRITE (*,*) 'Entrer A'
* Lit une matrice simple precision A de dimension M x N
* M, N et A sont modifiés.
      CALL SREAD_MPL (M, N, A, DMAX)
      IF (M .NE. N) STOP 'Matrice carree attendue'
      WRITE (*,*) 'Entrer b'
      CALL SREAD_MPL (M, 1, B, DMAX)
* Factorisation A = P . L . U
      CALL SGETRF (M, M, A, DMAX, IPIV, INFO)
      IF (INFO .NE. 0) STOP 'Erreur SGETRF'
* Affiche A, matrice M x M simple precision
      CALL SPRINT_MPL ('A', M, M, A, DMAX)
* Affiche le vecteur colonne d'entiers IPIV de dimension M
      CALL IPRINT_MPL ('IPIV', M, 1, IPIV, DMAX)
* Résolution
      CALL SGETRS('No Transpose', M, 1, A, DMAX, IPIV, B, DMAX, INFO)
      IF (INFO .NE. 0) STOP 'Erreur SGETRS'
* Affiche le vecteur colonne simple precision B de dimension M.
      CALL SPRINT_MPL ('Solution', M, 1, B, DMAX)
      END PROGRAM
