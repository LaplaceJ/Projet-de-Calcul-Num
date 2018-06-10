* Le répertoire cnum est supposé placé dans le répertoire principal
* $HOME/cnum/gfor matrice_complexe.f
* cat matrice_complexe.dat
* ./a.out < matrice_complexe.dat
* ./a.out < matrice_complexe.dat | ./a.out

      PROGRAM COMP_MAT
      INTEGER DMAX
      PARAMETER (DMAX = 10)
      INTEGER M, N
      COMPLEX A(DMAX,DMAX)
* Lit une matrice complexe A de dimension M x N.
* M, N, A sont modifiés
      CALL CREAD_MPL (M, N, A, DMAX)
* Affiche la matrice lue.
      CALL CPRINT_MPL ('A', M, N, A, DMAX)
      END PROGRAM
