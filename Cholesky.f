* fonction factorisant une matrice sous la 
* forme de cholesky.  A = btb
*
* Cette fonction calcul la matrice b diagonal supérieur 
* ou inférieur.
*
*
*
*
*
* Entier : 
* UPLO matrice supérieur(1) ou inférieur(0)
* N nombre de colonnes de la matrice 
* LDA La dimension principale de A
* i, j, k indices de boucle
*
* DOUBLE PRECISION
* A matrice à factoriser
* tmp variable tampon


      SUBROUTINE DPOTRF2 ( UPLO, N, A, LDA)
      IMPLICIT NONE
      
      INTEGER LDA , N  ,UPLO , i, j , k   
      DOUBLE PRECISION A(LDA,N) , tmp
      
      IF (UPLO.EQ.1) THEN 
      !UPER 
      
      DO i = 1 , N 
       !calcul de A(i,i) boucle 
        tmp = 0 
       DO j = 1 , i - 1 
         tmp = tmp + A(j,i)**2 
       ENDDO
    
       
       A(i,i) = sqrt(A(i,i) - tmp )
       
       
       ! calcul des éléments de la ligne i
       
       DO j = i + 1 , N 
        !calcule de la ligne i , j 
        
        DO k  = 1 , i -1
         A(i,j) = A(i,j) - A( k, i) * A( k,j )
        ENDDO
        
        A(i,j)   = A(i,j)  / A(i,i)  
       ENDDO
      ENDDO
      
      ELSE
* LOWER 
      
      DO i = 1 , N 
       !calcul de A(i,i) boucle 
        tmp = 0 
       DO j = 1 , i - 1 
         tmp = tmp + A(i,j)**2 
       ENDDO
       
       
       A(i,i) = sqrt(A(i,i) - tmp )
       
       
       ! calcul des éléments de la ligne i
       
       DO j = i + 1 , N 
        !calcule de la ligne i , j 
        
        DO k  = 1 , i -1
         A(j,i) = A(j,i) - A( i, k) * A( j,k )
        ENDDO
        
        A(j,i)   = A(j,i)  / A(i,i)  
       ENDDO
      ENDDO
      
      ENDIF
      END SUBROUTINE
