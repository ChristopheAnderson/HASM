import numpy as np
import matplotlib.pyplot as plt
import math

# Fonction exemple
def example_function(x, y):
    return 1/3*np.exp(-(81/4)*((x-0.5)**2+(y-0.5)**2) )

# Paramètres de discrétisation
x_range = (0, 1)
y_range = (0, 1)
num_points_x = 50  # Exemple de nombre de points x
num_points_y = 50 # Exemple de nombre de points y

# Discrétisation des variables x et y
x_vals = np.linspace(x_range[0], x_range[1], num_points_x)
y_vals = np.linspace(y_range[0], y_range[1], num_points_y)
X, Y = np.meshgrid(x_vals, y_vals)

# Calcul des valeurs de la fonction sur la grille
Z = example_function(X, Y)

# Affichage en 2D avec un nuage de points
fig=plt.figure(figsize=(8, 6))
#plt.scatter(X, Y, c=Z, cmap='viridis')
#plt.colorbar(label='Valeur de la fonction')
#contours = plt.contour(X, Y, Z, cmap='viridis')
#plt.colorbar(contours, label='Valeur de la fonction')
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z, cmap='viridis')
plt.title('Représentation 2D de la fonction avec un nuage de points')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(True)
plt.show()
