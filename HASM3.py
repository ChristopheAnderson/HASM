import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import pandas as pd



I,J=50,50
I1,J1=(I-2),(J-2)
n=I1*J1
lambd=1 # prenom la valeur de lambda=1


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



def discretize_2d_function(func, x_range, y_range, num_points_x, num_points_y):
    # Discrétisation des domaines x et y
    x_vals = np.linspace(x_range[0], x_range[1], num_points_x)
    y_vals = np.linspace(y_range[0], y_range[1], num_points_y)
    
    # Création d'une grille de points 2D
    Y, X = np.meshgrid(y_vals, x_vals)
    
    # Calcul des valeurs de la fonction en chaque point
    Z = func(X, Y)
    
    # Création du dictionnaire pour stocker les valeurs de la fonction avec leur numérotation horizontale
    points_dict = {}
    for i, x in enumerate(x_vals):
        for j, y in enumerate(y_vals):
            point_num = i * num_points_y + j
            points_dict[point_num] = {'x': x, 'y': y, 'value': Z[i, j]}
    
    return X, Y, Z, points_dict

# Fonction exemple
def example_function(x, y):
    return np.cos(10*y) + np.sin(10*(x-y))

# Paramètres de discrétisation

x_range = (0, 1)
y_range = (0, 1)
num_points_x = I
num_points_y = J

x_vals = np.linspace(x_range[0], x_range[1], num_points_x)
h=x_vals[1]-x_vals[0]

# Discrétisation et récupération des valeurs de la fonction
X, Y, Z, points_dict = discretize_2d_function(example_function, x_range, y_range, num_points_x, num_points_y)

# Affichage des résultats
#print("Valeurs de la fonction en chaque point de la grille:")
#print(Z)
#print("\nDictionnaire des points de discrétisation avec leurs valeurs:")
#for point_num, info in points_dict.items():
    #print(f"Point {point_num}: (x={info['x']}, y={info['y']}, value={info['value']})")


def coro(i,j):
    retour=i*(J)+j
    if retour<0 or i<0:
        retour=0
    else :
        retour=points_dict[retour]['value']
    return retour



#Valeurs initiales
fo = np.zeros((I*J,1))

# Valeurs de la fonction au J+1 sampling point 

fonc_init = np.zeros((n,1))
x= random.sample(range(0,I-2),I1)
y= random.sample(range(0,J-2),J1)
for i in x:
    for j in y:
        fonc_init[i*(J-2)+j,0]=coro(i+1, j+1)
        fo[(i+1)*J+(j+1),0] =coro(i+1, j+1) 


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
for i in x:
    for j in y:
        Sampling[i*(J-2)+j,i*(J-2)+j]=1
    

        

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
RMSE=0
MAE=0
for i in range(I1):
    for j in range(J1):
        MAE=MAE+np.abs(Z[i+1,j+1]-Solu[i,j])
        RMSE=RMSE+(Z[i+1,j+1]-Solu[i,j])**2
RMSE=math.sqrt(RMSE/(I1*J1))  
MAE=(MAE/(I1*J1))
print("RMSE= ", RMSE," MAE=", MAE)    

RMSE_rela=0
for i in x:
    for j in y:
        RMSE_rela=RMSE_rela+(Z[i+1, j+1]-Solu[i,j])**2
RMSE_rela=math.sqrt(RMSE_rela/(len(x)*len(y)))
print("RMSE relative au points échantillonné      RMSE= ", RMSE_rela)

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
plt.title('Représentation 2D de la fonction avec un nuage de points')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()

# Insérer la matrice à partir de la position spécifiée

#print("Matrice après insertion :\n", result_matrix)
