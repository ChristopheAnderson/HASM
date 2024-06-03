import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math


I,J=10,10
n=I*J
lambd=2 # prenom la valeur de lambda=1

def cor(i,j):
    retour=i*(J-1)+j
    return retour

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
    X, Y = np.meshgrid(x_vals, y_vals)
    
    # Calcul des valeurs de la fonction en chaque point
    Z = func(X, Y)
    
    # Création du dictionnaire pour stocker les valeurs de la fonction avec leur numérotation horizontale
    points_dict = {}
    for i, x in enumerate(x_vals):
        for j, y in enumerate(y_vals):
            point_num = i * num_points_y + j
            points_dict[point_num] = {'x': x, 'y': y, 'value': Z[j, i]}
    
    return X, Y, Z, points_dict

# Fonction exemple
def example_function(x, y):
    return np.sin(x) + np.cos(y)

# Paramètres de discrétisation

x_range = (0, 2*np.pi)
y_range = (0, 2*np.pi)
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







# The sampling points are the first J points
Sampling = np.zeros((n, n))
for i in range(J):
    Sampling[i,i]=1
# Valeurs de la fonction au sampling point
fonc_init=np.zeros((n,1))
for i in range(n):
    if i < J :
        fonc_init[i,0]=points_dict[i]['value']



# Matrice principale
main_matrix = np.zeros((n, n))
matrix_B= np.zeros((n, n))
# Première matrice à insérer
inserted_matrix=np.zeros((J,J))
np.fill_diagonal(inserted_matrix , 1)
# Deuxième matrice à insérer
matrix_b=tri_diagonal_matrix(J)

for i in range(I):
    row_start = i*J
    col_start = i*J
    result_matrix = insert_matrix(main_matrix, -2*inserted_matrix, row_start, col_start)
    matrix_B = insert_matrix(matrix_B, matrix_b, row_start, col_start)
    if i>=1:
        result_matrix[J:][row_start-J:row_start-J+inserted_matrix.shape[0], col_start-J:col_start-J+inserted_matrix.shape[1]] = inserted_matrix
        result_matrix[:,J:][row_start-J:row_start-J+inserted_matrix.shape[0], col_start-J:col_start-J+inserted_matrix.shape[1]] = inserted_matrix

    main_matrix = np.copy(result_matrix)  
    
    


#Les coefficients de la premiere equation diff et coef de Christofell




# Ajout de array2 à la suite de array1 avec vstack
fonc_init= np.vstack([fonc_init, np.zeros((3*J,1))])



E= np.zeros((I,J))
G= np.zeros((I,J))
F= np.zeros((I,J))
L= np.zeros((I,J))
N= np.zeros((I,J))
M= np.zeros((I,J))
f=fonc_init

for i in range(I):
    for j in range(J):
        a1,a2,a3,a4=1,1,1,1
        if (i+1)>=I :
            a1=0
        if (i-1)<0 :
            a2=0
        if (j+1)>=J :
            a3=0
        if (j-1)<0 :
            a4=0
    
        E[i,j] = 1 + ((a1*f[cor(i+1, j), 0] - a2*f[cor(i-1, j), 0]) / (2 * h))**2
        G[i,j] = 1 + ((a3*f[cor(i, j+1), 0] - a4*f[cor(i, j-1), 0]) / (2 * h))**2
        F[i,j] = ((a1*f[cor(i+1, j), 0] - a2*f[cor(i-1, j), 0]) / (2 * h))*((a3*f[cor(i, j+1), 0] - a4*f[cor(i, j-1), 0]) / (2 * h))
        L[i,j] =((a1*f[cor(i+1, j), 0] - 2*f[cor(i, j), 0] +  a2*f[cor(i-1, j), 0])/h**2)/math.sqrt(1+ ((a1*f[cor(i+1, j), 0] - a2*f[cor(i-1, j), 0]) / (2 * h))**2 + ((a3*f[cor(i, j+1), 0] - a4*f[cor(i, j-1), 0]) / (2 * h))**2 )
        N[i,j] =((a3*f[cor(i, j+1), 0] - 2*f[cor(i, j), 0] +  a4*f[cor(i, j-1), 0])/h**2)/math.sqrt(1+ ((a1*f[cor(i+1, j), 0] - a2*f[cor(i-1, j), 0]) / (2 * h))**2 + ((a3*f[cor(i, j+1), 0] - a4*f[cor(i, j-1), 0]) / (2 * h))**2 )
        M[i,j] =(( ((a1*a3*f[cor(i+1, j+1), 0] - a1*a4*f[cor(i+1, j-1), 0])) / 4* h**2) -( ((a2*a3*f[cor(i-1, j+1), 0] - a2*a4*f[cor(i-1, j-1), 0])) / 4* h**2)  )/math.sqrt(1+ (( a1*f[cor(i+1, j), 0] - a2*f[cor(i-1, j), 0]) / (2 * h))**2 + ((a3*f[cor(i, j+1), 0] - a4*f[cor(i, j-1), 0]) / (2 * h))**2 )
        

        
        
#fonc_init=f
# Ajout de array2 à la suite de array1 avec vstack
fonc_init= np.vstack([fonc_init, np.zeros((3*J,1))])
E= np.vstack([E, np.zeros((2,J))])
G= np.vstack([G, np.zeros((2,J))])
F= np.vstack([F, np.zeros((2,J))])
L= np.vstack([L, np.zeros((2,J))])
N= np.vstack([N, np.zeros((2,J))])
M= np.vstack([M, np.zeros((2,J))])
# Ajout de tableau2 à la suite de tableau1 avec hstack
E= np.hstack((E, np.zeros((I+2,2))))
G= np.hstack((G, np.zeros((I+2,2))))
F= np.hstack((F, np.zeros((I+2,2))))
L= np.hstack((L, np.zeros((I+2,2))))
N= np.hstack((N, np.zeros((I+2,2))))
M= np.hstack((M, np.zeros((I+2,2))))



G111= np.zeros((I,J))
G211= np.zeros((I,J))
G122= np.zeros((I,J))
G222= np.zeros((I,J))
G112= np.zeros((I,J))
G212= np.zeros((I,J))


for i in range(I):
    for j in range(J):
        a1,a2,a3,a4=1,1,1,1
        if (i+1)>=I :
            a1=0
        if (i-1)<0 :
            a2=0
        if (j+1)>=J :
            a3=0
        if (j-1)<0 :
            a4=0
    
        G111[i,j] = (  G[i,j]*(a1*E[i+1,j]-a2*E[i-1,j]) -2*F[i,j]*(a1*F[i+1,j]-a2*F[i-1,j]) + F[i,j]*(a3*E[i,j+1]-a4*E[i,j-1])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
        G211[i,j] = (  2*E[i,j]*(a1*F[i+1,j]-a2*F[i-1,j]) -E[i,j]*(a3*E[i,j+1]-a4*E[i,j-1]) - F[i,j]*(a1*F[i+1,j]-a2*F[i-1,j])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
        G122[i,j] = (  2*G[i,j]*(a3*F[i,j+1]-a4*F[i,j-1]) -G[i,j]*(a1*G[i+1,j]-a2*G[i-1,j]) - F[i,j]*(a3*G[i,j+1]-a4*G[i,j-1])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
        G222[i,j] = (  E[i,j]*(a3*G[i,j+1]-a4*G[i,j-1]) -2*F[i,j]*(a3*F[i,j+1]-a4*F[i,j-1]) + F[i,j]*(a1*G[i+1,j]-a2*G[i-1,j])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
        G112[i,j] = (  G[i,j]*(a3*E[i,j+1]-a4*E[i,j-1])  - F[i,j]*(a1*G[i+1,j]-a2*G[i-1,j])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )
        G212[i,j] = (  E[i,j]*(a1*G[i+1,j]-a2*G[i-1,j])  - F[i,j]*(a3*F[i,j+1]-a4*F[i,j-1])  ) /( 4*h*(E[i,j]*G[i,j]-F[i,j]**2  ) )




# Construire le premier vecteur de second membre

D=np.empty((0,1))
mon_dictionnaire = {}
for i in range(I):
    cle = f'd{i}'
    valeur = np.zeros((J,1))
    a5,a6=0,0
    if (i-1)==-1 :
        a5=1
    if (i+1)==I :
        a6=1
    for j in range(J):
        a1,a2,a3,a4=1,1,1,1
        if (i+1)>=I :
            a1=0
        if (i-1)<0 :
            a2=0
        if (J+1)>=J :
            a3=0
        if (j-1)<0 :
            a4=0
        valeur[j,0]=(1/2)*(a1*f[cor(i+1, j), 0]-a2*f[cor(i-1, j), 0])*G111[i,j]*h + (1/2)*(a3*f[cor(i, j+1), 0]-a4*f[cor(i, j-1), 0])*G211[i,j]*h + L[i,j]*h**2/math.sqrt(E[i,j]+G[i,j]-1) - a5*f[cor(i-1, j), 0] - a6*f[cor(i+1, j), 0]
     
    mon_dictionnaire[cle]=valeur
    D = np.vstack((D, valeur))
    



# Construire le deuxieme vecteur de second membre

Q=np.empty((0,1))
mon_dictionnaire2 = {}
for i in range(I):
    cle = f'q{i}'
    valeur = np.zeros((J,1))
   
    for j in range(J):
        a1,a2,a3,a4,a5,a6=1,1,1,1,0,0
        if (i+1)>=I :
            a1=0
        if (i-1)<0 :
            a2=0
        if (J+1)>=J :
            a3=0
        if (j-1)<0 :
            a4=0
        if (j-1)==-1 :
            a5=1
        if (j+1)==J :
            a6=1
        valeur[j,0]=(1/2)*(a1*f[cor(i+1, j), 0]-a2*f[cor(i-1, j), 0])*G122[i,j]*h + (1/2)*(a3*f[cor(i, j+1), 0]-a4*f[cor(i, j-1), 0])*G222[i,j]*h + N[i,j]*h**2/math.sqrt(E[i,j]+G[i,j]-1) - a4*a5*f[cor(i, j-1), 0] - a3*a6*f[cor(i, j+1), 0]
     
    mon_dictionnaire2[cle]=valeur
    Q = np.vstack((Q, valeur))















# Supprimer la(es) derniere(s) ligne(s) des tableaux
fonc_init = f
fonc_init = fonc_init[:-3*J]
E = E[:-2]
G = G[:-2]
F = F[:-2]
L = L[:-2]
N = N[:-2]
M = M[:-2]
# Supprimer les deux dernieres lignes des tableaux
E = E[:, :-2]
G = G[:, :-2]
F = F[:, :-2]
L = L[:, :-2]
N = N[:, :-2]
M = M[:, :-2]



# Résolution du problème

inverse_prod = np.linalg.inv( np.dot(main_matrix.T, main_matrix) + np.dot(matrix_B.T, matrix_B) + lambd**2*np.dot(Sampling.T, Sampling) )
second_prod = ( np.dot(main_matrix.T, D ) + np.dot(matrix_B.T, Q ) + lambd**2*np.dot(Sampling.T, fonc_init ) )

first = np.dot(main_matrix.T, main_matrix) + np.dot(matrix_B.T, matrix_B) + lambd**2*np.dot(Sampling.T, Sampling)
diagonale = np.diag(np.diag(first))
supérieure = -(np.triu(first)- np.diag(np.diag(first)))
inférieure = -(np.tril(first) - np.diag(np.diag(first)))

for i in range(8):
    Z_final =  np.dot(np.dot( np.linalg.inv(diagonale-inférieure) , supérieure ), fonc_init)  +  np.dot( np.linalg.inv(diagonale-inférieure) , second_prod )
    fonc_init=Z_final

Z_solu = np.dot( inverse_prod , second_prod )

Solu = np.reshape(Z_solu, ( I , J ))



fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
plt.contourf(X, Y, Solu, cmap='viridis')  # Utilisez contourf pour afficher une carte de chaleur
surf = ax.plot_surface(X, Y, Solu, cmap='viridis', edgecolor='none')
fig.colorbar(surf, shrink=0.5, aspect=5)
ax.set_title('Représentation 3D de la fonction')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.show()

# Insérer la matrice à partir de la position spécifiée

#print("Matrice après insertion :\n", result_matrix)
