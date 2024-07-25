import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random
import pandas as pd

# Définir la fonction
def f(x, y):
    return x**2 - y**2
I=50
J=50
# Créer des données pour les axes x et y
x = np.linspace(0, 1, I)
y = np.linspace(0, 1, J)
Y, X = np.meshgrid(y, x)

# Calculer les valeurs de z
Z = f(X, Y)
# X_num=np.linspace(1, I, I)
# y_num=np.linspace(1, J, J)

np.random.seed(42)
point_number=2*I+2*(J-2)+46*46
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

x_val , y_val = choix(2)

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


# Créer une figure et un axe 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Tracer la surface
ax.plot_surface(X, Y, Z, cmap='viridis')

# Ajouter des étiquettes aux axes
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')

# Titre
ax.set_title('Surface plot of f(x, y) = 1/3*np.exp(-(81/4)*((x-0.5)**2+(y-0.5)**2')

# Afficher le graphique
plt.show()

