* Ces deux fonctions permettent de faire une factorisation
* QR.
*
*
* La première subroutine calcule la factorisation de A est
* une matrice permettant de calculer plusieurs B pour résoudre
* plusieurs systems différents. 
* 
* La seconde calcule un vecteur B avec la factorisation QR de A. 
* Il ne reste plus qu'à faire une remontée pour obtenir le
* vecteur d'inconnues X
* 
*
* Fonctionnement: 
* 
* Nous utilisons une matrice (TMP) de la forme suivante
* |v1,1|0   |0   |0   |
* |... |v2,1|0   |0   |
* |... |... |... |0   |
* |... |... |... |0   |
* |v1,n|v2,n|vi,n|vn,n|
*
* Ou Vi est un vecteur utiliser lors d'une itération de 
* householder sur la matrice A 
* 
* NB: Une méthode plus optimale pour faire ceci est de 
* stocker les Vi sous la matrice A; où les 0 apparaissent 
* grâce à la méthode QR. Il nous suffirait de stocker v1 
* appart. Ainsi nous aurions économisé de la place en 
* mémoire.
*
* Obligation : 
* - Le nombre de colonnes de A == nombre de colonnes de B
* - A peut être rectangle 
* - les paramètres sont considérés comme bon 


* Entier : 
* LDT La dimension principale de TMP
* MT Nombre de lignes de TMP
* NT Nombre de colonne de TMP
* MB Nombre de lignes de QTB
* j,k indices de boucle 
*
* DOUBLE PRECISION :
* QTB Vecteur B à transformer 
* TMP matrice des vecteurs i
* vtvSurDeux variable tampon du calcul
* VtV / 2 
* scal variable tampon du calcul VB 
*
      SUBROUTINE DORMQR2( LDT, MT, NT, TMP,  QTB,  MB)
       IMPLICIT NONE 

        INTEGER MT, NT , MB  , j , k  , LDT
        DOUBLE PRECISION QTB(MB) , TMP(LDT, NT) , vtvSurDeux, scal
         
         do j =  1 , MB  
          scal = 0.0 
          vtvSurDeux = 0.0 
          
          ! calcul de scal
          do k = j , MT 
           scal = TMP(k,j) * QTB(k) + scal 
          end do 
          
          !calcul de vtvSurDeux
          do k = j , MT
            vtvSurDeux = TMP(k,j)*TMP(k,j) + vtvSurDeux
          enddo
          
          vtvSurDeux = vtvSurDeux / 2 
          
          do k = j , MT
           QTB(k) = QTB(k) - scal * TMP(k,j) / vtvSurDeux  
          end do 
         end do 
        END SUBROUTINE  

 
 
* Entier : 
* LDA La dimension principale de A
* LDA Nombre de lignes de A 
* N Nombre de colonne de A
* M Nombre de lignes de A
* LDtmp La dimension principale de TMP
* Ntmp Nombre de lignes de TMP 
* i, j, k indices de boucle 
*
* DOUBLE PRECISION :
* A matrice avec laquelle ou on applique
* la méthode QR 
* TMP matrice des vecteurs i
* vtvSurDeux variable tampon du calcul
* VtV / 2 
* scal variable tampon du calcul V A i+1 
* norm norme 2 de Ai
      SUBROUTINE DGEQRF2 (M, N, A, LDA, TMP, LDtmp ,Ntmp )
      IMPLICIT NONE
      
      INTEGER LDA , N  , i, j , k  , M  , LDtmp ,Ntmp

      DOUBLE PRECISION A(LDA,N) , TMP(LDtmp ,Ntmp )  
      DOUBLE PRECISION vtvSurDeux , norm , scal
       do i = 1 , N 
        
         ! calcul de la norme²
         norm = 0.0 
         do j = i , M 
          norm = A(j,i)**2  + norm
         end do 
         
         
         ! calcul de vt v  / 2 
         vtvSurDeux = norm 
         norm = sqrt(norm) 
         
         
         if(A(i,i).LT.0) then
          ! si A(i,i) < 0 alors on minimise TMP(i,i): 
          vtvSurDeux = vtvSurDeux - norm * A(i,i)
          norm = - norm 
         else 
         ! si A(i,i) > 0 alors on maximise TMP(i,i): 
          vtvSurDeux = vtvSurDeux + norm * A(i,i)
         endif
    
         !formation du vecteur v
         TMP(i,i) = (A(i,i) + norm )  
         
         !reste du vecteur v 
         do j = i + 1 , M 
          TMP(j,i) = A(j,i) 
         end do 
         
        
        
         
         ! formation de Hn * An
         A(i,i) = -norm
         do j = i + 1 , N 
          scal = 0.0 
          do k = i , M 
           scal = TMP(k,i) * A(k,j) + scal 
          end do 
          do k = i , M 
           A(k,j) = A(k,j) - scal * TMP(k,i) / vtvSurDeux
          end do 
         end do 
        enddo
        
       
        
      END SUBROUTINE    
          
 
   
