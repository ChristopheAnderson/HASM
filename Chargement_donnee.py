import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import pandas as pd

# Définir la fonction
def f(x, y):
    return np.cos(10*y) + np.sin(10*(x-y))

#points initiaux
a,b=0,0
# nombre de points sur les deux axes
I,J=50,60
# pas commun des deux axes
h=0.01
# Créer des données pour les axes x et y
x=[ a + i*h for i in range (I) ]
y=[ b + j*h for j in range (J) ]

x_range=(a , x[I-1])
y_range=(b , y[J-1])
# pgcd = math.gcd(I, J)
# pgcd = pgcd/100
# ppcm=(I*J)/pgcd

Y, X = np.meshgrid(y, x)

# Calculer les valeurs de z
Z = f(X, Y)
# X_num=np.linspace(1, I, I)
# y_num=np.linspace(1, J, J)

np.random.seed(42)
point_number=2*I+2*(J-2)+46*56
inner_point=point_number-2*I-2*(J-2)


def generer_couples_aleatoires_numpy(Nx, Ny):
    # Créer une liste de toutes les paires possibles (x, y)
    toutes_paires = np.array([(x, y) for x in range(2,Nx) for y in range(2,Ny)])
    return toutes_paires

def choix(N):
    
    if(N==1):
        
        # Générer les couples aléatoires sans répétition
        couples_aleatoires = generer_couples_aleatoires_numpy(I, J)
        # Mélanger les paires de manière aléatoire
        np.random.shuffle(couples_aleatoires)
        x_val=couples_aleatoires[0:inner_point,0]
        y_val=couples_aleatoires[0:inner_point,1]
            
    
    if(N==2):
        
        np.random.seed(42)
        x_val=np.random.randint(2 , I , size=inner_point)
        y_val=np.random.randint(2 , J , size=inner_point)
    
    return(x_val,y_val)

x_val , y_val = choix(1)

z_val=np.empty(point_number)
incre=0
#For 2*I boundary conditions
x_inter=np.empty(0)
y_inter=np.empty(0)
for i in [0 , I-1]:
    for j in range(J):
        z_val[incre]=Z[i, j]
        # Ajouter la valeur au début du vecteur
        x_inter= np.append(x_inter, i+1)
        y_inter= np.append(y_inter, j+1)
        incre=incre + 1

#For 2*(J-2) boundary conditions
for j in [0 , J-1]:
    for i in range(1,I-1):
        z_val[incre]=Z[i, j]
        x_inter= np.append(x_inter, i+1)
        y_inter= np.append(y_inter, j+1)
        incre=incre + 1


for i in range(inner_point):
    z_val[incre]=Z[x_val[i]-1, y_val[i]-1]
    x_inter= np.append(x_inter, x_val[i])
    y_inter= np.append(y_inter, y_val[i])
    incre=incre + 1


matrice = np.vstack((x_inter, y_inter, z_val))
matrice=matrice.T

# Convertir la matrice en DataFrame
fichier_excel = pd.DataFrame(matrice)
# Enregistrer la DataFrame dans un fichier Excel
fichier_excel.to_excel("matrice.xlsx", index=False, header=False)
# Pareil pour les vecteurs de la fonction
fichier_fonction = pd.DataFrame(Z)
# Enregistrer la DataFrame dans un fichier Excel
fichier_fonction.to_excel("fonction.xlsx", index=False, header=False)


# Créer une figure et un axe 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tracer la surface
ax.plot_surface(X, Y, Z, cmap='viridis')

ax.set_zlim(np.min(Z)+1, np.max(Z)+1)


# Ajouter des étiquettes aux axes
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

# Titre
ax.set_title('Surface plot of f(x, y) =')
# 1/3*np.exp(-(81/4)*((x-0.5)**2+(y-0.5)**2

# Afficher le graphique
plt.show()

