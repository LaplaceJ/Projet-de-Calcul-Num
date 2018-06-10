* Fonction faisant une remontée ou une descente 
* sur une matrice à transposer ou non.
* 
* Entier : 
* UPLO matrice supérieur == 1 (inférieur == 0 )  
* TRANS transposé la matrice lors des calculs == 1 
* si non 0 
* N nombre de colonne de la matrice A 
* LDA La dimension principale de A
* i,j indices de boucle 
*
* DOUBLE PRECISION :  
* Matrice A(LDA,N) de double 
* X(N) vecteur solution  
* 
* Obligation : 
* - Le nombre de colonnes de A == nombre de 
* colonnes de X
* - A peut être rectangle 
* - A doit être triangulaire (inférieur ou supérieur) 
* - les paramètres sont considérés comme bon 


     
         SUBROUTINE DTRSV2 (UPLO, TRANS,  N, A, LDA, X )
       IMPLICIT NONE

      INTEGER LDA , N  , UPLO, i, j , TRANS    
      DOUBLE PRECISION A(LDA,N),X(N)
      
      IF (UPLO.EQ.1) THEN 
      ! remonté 
      
      X(n) = X(n) / A(n,n) 
      
      DO i = n-1 , 1 , -1
       DO j = i + 1  , n 
       
          !pas de transposer
          IF (TRANS.EQ.0) THEN 
           X(i) = X(i) -  X(j) * A(i,j)
          ELSE !transposer
            X(i) = X(i) -  X(j) * A(j,i)
          ENDIF
       ENDDO
       ! 1 sur la diagonal
       X(i) = X(i) / A(i,i) 
      ENDDO 
      
      ELSE
      ! decescente 
      ! A(n,n) = 1 
      X(1) = X(1) / A(1,1) 
      
      DO i = 2 , n 
       DO j = 1  , i - 1  
        !pas de transposer
        IF (TRANS.EQ.0) THEN 
           X(i) = X(i) -  X(j) * A(i,j)
          ELSE!transposer
            X(i) = X(i) -  X(j) * A(j,i)
          ENDIF
       ENDDO
       ! 1 sur la diagonal 
       X(i) = X(i) / A(i,i) 
      ENDDO 

      ENDIF
      END SUBROUTINE
      
