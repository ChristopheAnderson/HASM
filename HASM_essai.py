import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import pandas as pd








def insert_matrix(main_matrix, inserted_matrix, row_start, col_start):
    # Vérifier que les dimensions de la matrice à insérer sont valides
    if row_start + inserted_matrix.shape[0] > main_matrix.shape[0] or col_start + inserted_matrix.shape[1] > main_matrix.shape[1]:
        raise ValueError("Les dimensions de la matrice à insérer sont trop grandes pour la matrice principale.")
    
    # Copier la matrice principale pour éviter les modifications inattendues
    result_matrix = np.copy(main_matrix)
    
    # Insérer la matrice à partir de la position spécifiée
    result_matrix[row_start:row_start+inserted_matrix.shape[0], col_start:col_start+inserted_matrix.shape[1]] = inserted_matrix
    
    return result_matrix


def tri_diagonal_matrix(n):
    # Créer une matrice remplie de zéros
    A = np.zeros((n, n))

    # Remplir la diagonale principale avec des 1
    np.fill_diagonal(A, -2)

    # Remplir les diagonales supérieure et inférieure avec des 1
    np.fill_diagonal(A[1:], 1)
    np.fill_diagonal(A[:, 1:], 1)

    return A




def nombre_de_couples_uniques(matrice):
    # Utiliser un ensemble pour stocker les couples uniques
    couples_uniques = set()
    
    # Parcourir chaque ligne de la matrice et ajouter les couples à l'ensemble
    for ligne in matrice:
        couple = (ligne[0], ligne[1])
        couples_uniques.add(couple)
    
    # Retourner le nombre de couples uniques
    return len(couples_uniques)




def HASM(nom_fichier,interval_x,interval_y,fichier_fonction,forma, lambd=100):

    # Charger le fichier Excel sans en-têtes
    df = pd.read_excel(nom_fichier, header=None)
    # Accéder à chaque colonne par son indice
    colonne0 = (df[0]).values
    colonne1 = (df[1]).values
    colonne2 = (df[2]).values
    # Convertir le vecteur en entiers
    colonne0= colonne0.astype(int)
    colonne1= colonne1.astype(int)
    colonne2= colonne2.astype(float)
    
    I,J=max(colonne0),max(colonne1)
    I1,J1=(I-2),(J-2)
    n=I1*J1
    # lambd=1    valeur de lambda par défaut 
    
    colonne0= colonne0-1
    colonne1= colonne1-1
    
    matrice = np.vstack((colonne0, colonne1, colonne2))
    matrice=matrice.T
    
    # Obtenir le nombre de couples uniques
    nombre_uniques = nombre_de_couples_uniques(matrice)
    
    if(nombre_uniques!=len(matrice)):
        print("Attention : Il y a au moins une répétition d'un point de coordonnée dans la base de donnée")
        return
    else:
        
        
        
        
        # Paramètres de discrétisation
        
        x_range = interval_x
        y_range = interval_y
        num_points_x = I
        num_points_y = J
        
        x_vals = np.linspace(x_range[0], x_range[1], num_points_x)
        y_vals = np.linspace(y_range[0], y_range[1], num_points_y)
        h1=x_vals[1]-x_vals[0]
        h2=y_vals[1]-y_vals[0]
        
        if(h1!=h2):
            print("Les pas pour les deux axes ne sont pas égaux")
            return
        
        h=h1
        
        
        
        def coro(i,j):
            retour=i*(J)+j
            if retour<0 or i<0:
                retour=0
            else :
                # Filtrer la matrice pour récupérer les lignes correspondantes
                resultat = matrice[(matrice[:, 0] == i) & (matrice[:, 1] == j)]
                # Récupérer les valeurs de la troisième colonne
                valeurs_col3 = resultat[:, 2]
                retour=valeurs_col3[0]
            return retour
        
        
        #Valeurs initiales
        fo = np.zeros((I*J,1))
        
        # Valeurs de la fonction au J+1 sampling point 
        
        fonc_init = np.zeros((n,1))
        x= colonne0
        y= colonne1
        #mask = (matrice[:, 0] != I-1) & (matrice[:, 1] !=J-1) & (matrice[:, 0] != 0) & (matrice[:, 1] != 0)
        matrice1=matrice[(matrice[:, 0] != I-1) & (matrice[:, 1] !=J-1) & (matrice[:, 0] != 0) & (matrice[:, 1] != 0)]
        nombre=len(matrice1)
        x1 , y1 = matrice1[:,0] , matrice1[:,1]
        x1 = x1.astype(int)
        y1 = y1.astype(int)
        x1=x1-1
        y1=y1-1
        for i in range(nombre):
            fonc_init[x1[i]*(J-2)+y1[i],0]=coro(x1[i]+1, y1[i]+1)
            fo[(x1[i]+1)*J+(y1[i]+1),0] =coro(x1[i]+1, y1[i]+1) 
        
        
        for i in [0,I-1]:
            for j in range(J):
                fo[i*J+j,0]=coro(i,j)
        
        for j in [0,J-1]:
            for i in range(I):
                fo[i*J+j,0]=coro(i,j)
                
        
        def cor(i,j):
            retour=i*(J)+j
            if retour<0 or i<0:
                retour=0
            else :
                retour=fo[retour,0]
            return retour
        
        # Use of all sampling point
        #Sampling = np.zeros((n, n))
        #for i in range(n):
            #Sampling[i,i]=1
            
        
        
        # Valeurs de la fonction au sampling point 
        # k=0
        # fonc_init = np.zeros((n,1))
        # for i in range(1,I-1):
        #     for j in range(1,J-1) :
        #         fonc_init[k,0]=cor(i, j)
        #         k=k+1
        
        
        
        
        # The sampling points are the first J+1 points
        Sampling = np.zeros((n, n))
        for i in range(nombre):
            Sampling[x1[i]*(J-2)+y1[i],x1[i]*(J-2)+y1[i]]=1
            
        
                
        
        # Matrice principale
        main_matrix = np.zeros((n, n))
        matrix_B= np.zeros((n, n))
        # Première matrice à insérer
        inserted_matrix=np.zeros((J1,J1))
        np.fill_diagonal(inserted_matrix , 1)
        # Deuxième matrice à insérer
        matrix_b=tri_diagonal_matrix(J1)
        
        for i in range(I1):
            row_start = i*J1
            col_start = i*J1
            result_matrix = insert_matrix(main_matrix, -2*inserted_matrix, row_start, col_start)
            matrix_B = insert_matrix(matrix_B, matrix_b, row_start, col_start)
            if i>=1:
                result_matrix[J1:][row_start-J1:row_start-J1+inserted_matrix.shape[0], col_start-J1:col_start-J1+inserted_matrix.shape[1]] = inserted_matrix
                result_matrix[:,J1:][row_start-J1:row_start-J1+inserted_matrix.shape[0], col_start-J1:col_start-J1+inserted_matrix.shape[1]] = inserted_matrix
        
            main_matrix = np.copy(result_matrix)  
            
            
        
        
        #Les coefficients de la premiere equation diff et coef de Christofell
        
        
        
        
        # Ajout de array2 à la suite de array1 avec vstack
        #fonc_init= np.vstack([fonc_init, np.zeros((3*J1,1))])
        
        
        E= np.zeros((I,J))
        G= np.zeros((I,J))
        F= np.zeros((I,J))
        L= np.zeros((I,J))
        N= np.zeros((I,J))
        M= np.zeros((I,J))
        #f=fonc_init
        
        for i in range(1,I-1):
            for j in range(1,J-1):
               
            
                E[i,j] = 1 + ((cor(i+1, j) - cor(i-1, j)) / (2 * h))**2
                G[i,j] = 1 + ((cor(i, j+1) - cor(i, j-1)) / (2 * h))**2
                F[i,j] = ((cor(i+1, j) - cor(i-1, j)) / (2 * h))*((cor(i, j+1) - cor(i, j-1)) / (2 * h))
                L[i,j] =((cor(i+1, j) - 2*cor(i, j) +  cor(i-1, j))/h**2)/math.sqrt(1+ ((cor(i+1, j) - cor(i-1, j)) / (2 * h))**2 + ((cor(i, j+1) - cor(i, j-1)) / (2 * h))**2 )
                N[i,j] =((cor(i, j+1)- 2*cor(i, j) +  cor(i, j-1))/h**2)/math.sqrt(1+ ((cor(i+1, j) - cor(i-1, j)) / (2 * h))**2 + ((cor(i, j+1) - cor(i, j-1)) / (2 * h))**2 )
                M[i,j] =(( ((cor(i+1, j+1) - cor(i+1, j-1))) / 4* h**2) -( ((cor(i-1, j+1) - cor(i-1, j-1))) / 4* h**2)  )/math.sqrt(1+ (( cor(i+1, j) - cor(i-1, j)) / (2 * h))**2 + ((cor(i, j+1) - cor(i, j-1)) / (2 * h))**2 )
                
        
        
        
        
        G111= np.zeros((I,J))
        G211= np.zeros((I,J))
        G122= np.zeros((I,J))
        G222= np.zeros((I,J))
        G112= np.zeros((I,J))
        G212= np.zeros((I,J))
        
        
        for i in range(1,I-1):
            for j in range(1,J-1):
                
            
                G111[i,j] = (  G[i,j]*(E[i+1,j]-E[i-1,j]) -2*F[i,j]*(F[i+1,j]-F[i-1,j]) + F[i,j]*(E[i,j+1]-E[i,j-1])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
                G211[i,j] = (  2*E[i,j]*(F[i+1,j]-F[i-1,j]) -E[i,j]*(E[i,j+1]-E[i,j-1]) - F[i,j]*(E[i+1,j]-E[i-1,j])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
                G122[i,j] = (  2*G[i,j]*(F[i,j+1]-F[i,j-1]) -G[i,j]*(G[i+1,j]-G[i-1,j]) - F[i,j]*(G[i,j+1]-G[i,j-1])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
                G222[i,j] = (  E[i,j]*(G[i,j+1]-G[i,j-1]) -2*F[i,j]*(F[i,j+1]-F[i,j-1]) + F[i,j]*(G[i+1,j]-G[i-1,j])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
                G112[i,j] = (  G[i,j]*(E[i,j+1]-E[i,j-1])  - F[i,j]*(G[i+1,j]-G[i-1,j])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
                G212[i,j] = (  E[i,j]*(G[i+1,j]-G[i-1,j])  - F[i,j]*(F[i,j+1]-F[i,j-1])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
        
        
        
        
        # Construire le premier vecteur de second membre
        
        D=np.empty((0,1))
        mon_dictionnaire = {}
        for i in range(1,I-1):
            cle = f'd{i-1}'
            valeur = np.zeros((J1,1))
            a5,a6=0,0
            if (i-1)==0 :
                a5=1
            if (i+1)==(I-1) :
                a6=1
            for j in range(1,J-1):
                valeur[j-1,0]=(1/2)*(cor(i+1, j) - cor(i-1, j))*G111[i,j]*h + (1/2)*(cor(i, j+1)-cor(i, j-1))*G211[i,j]*h + L[i,j]*h**2/math.sqrt(E[i,j]+G[i,j]-1) - a5*cor(i-1, j) - a6*cor(i+1, j)
             
            mon_dictionnaire[cle]=valeur
            D = np.vstack((D, valeur))
            
        
        
        
        # Construire le deuxieme vecteur de second membre
        
        Q=np.empty((0,1))
        mon_dictionnaire2 = {}
        for i in range(1,I-1):
            cle = f'q{i}'
            valeur = np.zeros((J1,1))
           
            for j in range(1,J-1):
                a5,a6=0,0
                if (j-1)==0 :
                    a5=1
                if (j+1)==(J-1) :
                    a6=1
                valeur[j-1,0]=(1/2)*(cor(i+1, j)-cor(i-1, j))*G122[i,j]*h + (1/2)*(cor(i, j+1)-cor(i, j-1))*G222[i,j]*h + N[i,j]*h**2/math.sqrt(E[i,j]+G[i,j]-1) - a5*cor(i, j-1) - a6*cor(i, j+1)
             
            mon_dictionnaire2[cle]=valeur
            Q = np.vstack((Q, valeur))
        
        
        
        
        #fonc_init = f
        #fonc_init = fonc_init[:-3*J]
        
        # Supprimer la première et la derniere ligne
        E = E[1:I-1]
        G = G[1:I-1]
        F = F[1:I-1]
        L = L[1:I-1]
        N = N[1:I-1]
        M = M[1:I-1]
        # Supprimer la première et la dernière colonne
        E = E[:, 1:J-1]
        G = G[:, 1:J-1]
        F = F[:, 1:J-1]
        L = L[:, 1:J-1]
        N = N[:, 1:J-1]
        M = M[:, 1:J-1]
        
        
        
        # Résolution du problème
        
        # Méthode directe
        Z_solu=fonc_init
        nbr_ite=1
        for i in range(nbr_ite):
            inverse_prod = np.linalg.inv( np.dot(main_matrix.T, main_matrix) + np.dot(matrix_B.T, matrix_B) + lambd**2*np.dot(Sampling.T, Sampling) )
            second_prod = ( np.dot(main_matrix.T, D ) + np.dot(matrix_B.T, Q ) + lambd**2*np.dot(Sampling.T, Z_solu ) )
            Z_solu = np.dot( inverse_prod , second_prod )
        
        Solu = np.reshape(Z_solu, ( I1 , J1 )) # redimensionnement de la taille
        
        
        #Méthode itérative
        
        # first = np.dot(main_matrix.T, main_matrix) + np.dot(matrix_B.T, matrix_B) + lambd**2*np.dot(Sampling.T, Sampling)
        # diagonale = np.diag(np.diag(first))
        # supérieure = -(np.triu(first)- np.diag(np.diag(first)))
        # inférieure = -(np.tril(first) - np.diag(np.diag(first)))
        
        # for i in range(8):
        #     Z_final =  np.dot(np.dot( np.linalg.inv(diagonale-inférieure) , supérieure ), fonc_init)  +  np.dot( np.linalg.inv(diagonale-inférieure) , second_prod )
        #     fonc_init=Z_final
        
        # Métrique de performance
            #RMSE
        
        
        # Charger le fichier Excel sans en-têtes
        Z = pd.read_excel(fichier_fonction, header=None)
        Z=np.array(Z)
        
        moyenne=np.mean(np.array(Z[1:I1+1,1:J1+1]))
        dénominateur=0
        RMSE=0
        MAE=0
        for i in range(I1):
            for j in range(J1):
                MAE=MAE+np.abs(Z[i+1,j+1]-Solu[i,j])
                RMSE=RMSE+(Z[i+1,j+1]-Solu[i,j])**2
                dénominateur=dénominateur+(Z[i+1,j+1]-moyenne)**2
        numérateur=RMSE
        R=1-(numérateur/dénominateur)
        RMSE=math.sqrt(RMSE/(I1*J1))  
        MAE=(MAE/(I1*J1))
        print("RMSE=  ", RMSE," MAE=  ", MAE)    
        
        RMSE_rela=0
        for i in range(nombre):
                RMSE_rela=RMSE_rela+(Z[x1[i]+1, y1[i]+1]-Solu[x1[i],y1[i]])**2
        RMSE_rela=math.sqrt(RMSE_rela/(nombre))
        print("RMSE relative au points échantillonné      RMSE=  ", RMSE_rela ,"    R=  ",R)
        
        
        x = np.linspace(x_range[0], x_range[1], I)
        y = np.linspace(y_range[0], y_range[1], J)
        Y, X = np.meshgrid(y, x)
        
        # Supprimer la première et la derniere ligne
        X = X[1:I-1]
        Y = Y[1:I-1]
        
        # Supprimer la première et la dernière colonne
        X = X[:, 1:J-1]
        Y = Y[:, 1:J-1]
        
        
        

        
        # Affichage en 2D avec un nuage de points
        fig=plt.figure(figsize=(8, 6))
        #plt.scatter(X, Y, c=Solu, cmap='viridis')
        #plt.colorbar(label='Valeur de la fonction')
        #contours = plt.contour(X, Y, Solu, cmap='viridis')
        #plt.colorbar(contours, label='Valeur de la fonction')
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Solu, cmap='viridis')
        ax.set_title('Représentation 3D de la fonction')
        plt.title('Surface plot of f(x, y) =')
        #Représentation 2D de la fonction avec un nuage de points
        ax.set_zlim(np.min(Z)+1, np.max(Z)+1)
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.grid(True)
        plt.show()
        
        # Insérer la matrice à partir de la position spécifiée
        
        #print("Matrice après insertion :\n", result_matrix)
        
        if(forma==1):
            # donnée initiale à l'intérieur domaine
            donnee_initiale=Z[1:I1+1,1:J1+1]
            return (X,Y,Solu,donnee_initiale)
        if(forma==2):
            # donnée initiale à l'intérieur domaine et non au bord du domaine
            donnee_initiale=Z[1:I1+1,1:J1+1]
            X = np.reshape(X, ( I1*J1 , 1 )) # redimensionnement de la taille
            Y = np.reshape(Y, ( I1*J1 , 1 )) # redimensionnement de la taille
            Solu = np.reshape(Solu, ( I1*J1 , 1 )) # redimensionnement de la taille
            donnee_initiale = np.reshape(donnee_initiale, ( I1*J1 , 1 )) # redimensionnement de la taille
            
            return (X,Y,Solu,donnee_initiale)
        else:
            return "Attention : Le cinquième paramètre qui est le format ne prend que deux valeurs \n Le format 1 renvoie le résultat de la prédiction sous forme de matrice \n le format 2 renvoie le résultat de la prédiction sous forme de vecteur"


nom_fichier='matrice.xlsx'
fichier_fonction='fonction.xlsx'
x_range = (0, 0.49)
y_range = (0, 0.59)
forma=2
# Le résultat de la prédiction de la méthode HASM renvoie sous le format indiqué dans l'ordre :
# 1- les valeurs des abscisses dicrétisés 
# 2- les valeurs des ordonnées dicrétisés 
# 3- les valeurs de prédiction de l'attributs aux points de coordonées 
#    dicrétisées ci-dessus
# 4- les valeurs des données exactes fournies de l'attibuts à l'intérieur du 
#    du domaine et non au bord du domaine
sortie_prediction= HASM(nom_fichier,x_range,y_range,fichier_fonction,forma,1000)
X,Y,Solu,donnee_initiale=sortie_prediction
# Convertir la matrice en DataFrame

if (forma==1):
    mon_fichier1= pd.DataFrame(X)
    mon_fichier2= pd.DataFrame(Y)
    mon_fichier3= pd.DataFrame(Solu)
    mon_fichier4= pd.DataFrame(donnee_initiale)
    #Enregistrer la DataFrame dans un fichier Excel
    mon_fichier1.to_excel("Sortie_prediction_format1/X1.xlsx", index=False, header=False)
    mon_fichier2.to_excel("Sortie_prediction_format1/Y1.xlsx", index=False, header=False)
    mon_fichier3.to_excel("Sortie_prediction_format1/prédiction1.xlsx", index=False, header=False)
    mon_fichier4.to_excel("Sortie_prediction_format1/Z_donnée1.xlsx", index=False, header=False)
    
    
 
if (forma==2):
    
    # Combiner les matrices en une seule matrice à plusieurs colonnes
    combined_matrix = np.column_stack((X, Y, Solu, donnee_initiale))
    
    # Transformer la matrice en DataFrame pandas
    mon_fichier= pd.DataFrame(combined_matrix, columns=['X1', 'Y1', 'prediction1','Z_donnee1'])
    
    # Spécifier le chemin complet où le fichier Excel sera enregistré
    chemin_fichier = 'Sortie_prediction_format2/Sortie_format1_prediction1.xlsx'
    
    # Enregistrer le DataFrame dans un fichier Excel
    mon_fichier.to_excel(chemin_fichier, index=False)

