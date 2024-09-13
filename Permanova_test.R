# Installer le package vegan si ce n'est pas déjà fait
if (!require(vegan)) install.packages("vegan")


library(vegan)

# Exemple de création de données
set.seed(123)
n <- 90 * 1000  # Total de 90 ensembles de données avec 1000 lignes chacun

data <- data.frame(
  Méthode = rep(c("BME", "HASM"), each = n / 2),
  Taille = rep(rep(c("64", "80", "169", "250", "400"), each = n / (2 * 5)), 2),
  Asymétrie = rep(rep(c("négatif", "symétrique", "positif"), each = n / (2 * 5 * 3)), 2),
  Dépendance = rep(rep(c("faible", "modéré", "fort"), each = n / (2 * 5 * 3 * 3)), 2),
  RMSE = rnorm(n, mean = 5, sd = 2),  # Exemple de données générées
  MAE = rnorm(n, mean = 3, sd = 1)   # Exemple de données générées
)

# Créer un tableau des variables explicatives pour PERMANOVA
factors <- data[, c("Méthode", "Taille", "Asymétrie", "Dépendance")]

# Calculer la matrice de distances pour RMSE et MAE
dist_rmse <- dist(data$RMSE)
dist_mae <- dist(data$MAE)



# PERMANOVA pour RMSE
permanova_rmse <- adonis(dist_rmse ~ Méthode * Taille * Asymétrie * Dépendance, data = factors, permutations = 999)
print(permanova_rmse)


# PERMANOVA pour MAE
permanova_mae <- adonis(dist_mae ~ Méthode * Taille * Asymétrie * Dépendance, data = factors, permutations = 999)
print(permanova_mae)



