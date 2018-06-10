* Le répertoire cnum est supposé placé dans le répertoire principal
* $HOME/cnum/gfor factorisation_LU_double.f
* cat factorisation_LU.dat
* ./a.out < factorisation_LU.dat

      PROGRAM FAC_LU
      INTEGER DMAX
      PARAMETER (DMAX=10)
      DOUBLE PRECISION A(DMAX,DMAX), B(DMAX)
      INTEGER IPIV(DMAX), INFO, M, N
*
      WRITE (*,*) 'Entrer A'
* Lit une matrice double precision A de dimension M x N
* M, N et A sont modifiés.
      CALL DREAD_MPL (M, N, A, DMAX)
      IF (M .NE. N) STOP 'Matrice carree attendue'
      WRITE (*,*) 'Entrer b'
* Lit un vecteur colonne double precision B de dimension M.
* M et B sont modifiés.
      CALL DREAD_MPL (M, 1, B, DMAX)
* Factorisation A = P . L . U
      CALL DGETRF (M, M, A, DMAX, IPIV, INFO)
      IF (INFO .NE. 0) STOP 'Erreur DGETRF'
* Affiche A, matrice M x M double precision
      CALL DPRINT_MPL ('A', M, M, A, DMAX)
* Affiche le vecteur colonne d'entiers IPIV de dimension M
      CALL IPRINT_MPL ('IPIV', M, 1, IPIV, DMAX)
* Résolution
      CALL DGETRS('No Transpose', M, 1, A, DMAX, IPIV, B, DMAX, INFO)
      IF (INFO .NE. 0) STOP 'Erreur DGETRS'
* Affiche le vecteur colonne double precision B de dimension M.
      CALL DPRINT_MPL ('Solution', M, 1, B, DMAX)
      END PROGRAM
