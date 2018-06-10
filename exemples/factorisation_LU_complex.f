* Le répertoire cnum est supposé placé dans le répertoire principal
* $HOME/cnum/gfor factorisation_LU_double.f
* cat factorisation_LU.dat
* ./a.out < factorisation_LU.dat

      PROGRAM FAC_LU
      INTEGER DMAX
      PARAMETER (DMAX=10)
      COMPLEX A(DMAX,DMAX), B(DMAX)
      INTEGER IPIV(DMAX), INFO, M, N
*
      WRITE (*,*) 'Entrer A'
* Lit une matrice complexe A de dimension M x N
* M, N et A sont modifiés.
      CALL CREAD_MPL (M, N, A, DMAX)
      IF (M .NE. N) STOP 'Matrice carree attendue'
      WRITE (*,*) 'Entrer b'
* Lit un vecteur colonne complexe B de dimension M.
* M et B sont modifiés.
      CALL CREAD_MPL (M, 1, B, DMAX)
* Factorisation A = P . L . U
      CALL CGETRF (M, M, A, DMAX, IPIV, INFO)
      IF (INFO .NE. 0) STOP 'Erreur CGETRF'
* Affiche A, matrice M x M complexe
      CALL CPRINT_MPL ('A', M, M, A, DMAX)
* Affiche le vecteur colonne d'entiers IPIV de dimension M
      CALL IPRINT_MPL ('IPIV', M, 1, IPIV, DMAX)
* Résolution
      CALL CGETRS('No Transpose', M, 1, A, DMAX, IPIV, B, DMAX, INFO)
      IF (INFO .NE. 0) STOP 'Erreur CGETRS'
* Affiche le vecteur colonne complexe B de dimension M.
      CALL CPRINT_MPL ('Solution', M, 1, B, DMAX)
      END PROGRAM
