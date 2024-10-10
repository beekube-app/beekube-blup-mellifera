library(jsonlite)
library(dplyr)
library(pedigreemm)
library(Matrix)
library(lme4)
library(brms)
library(rstan)
library(lm.beta)


# Fonction pour trier le pédigrée
sortPed <- function(ped) {
  ped$generation <- 0
  ped$processed <- FALSE

  while (any(!ped$processed)) {
    for (i in which(!ped$processed)) {
      if (is.na(ped$dam[i]) && is.na(ped$sire[i])) {
        ped$processed[i] <- TRUE
        next
      }

      dam_gen <- if (!is.na(ped$dam[i])) ped$generation[ped$id == ped$dam[i]] else -1
      sire_gen <- if (!is.na(ped$sire[i])) ped$generation[ped$id == ped$sire[i]] else -1

      if (length(dam_gen) == 0) dam_gen <- -1
      if (length(sire_gen) == 0) sire_gen <- -1

      if (all(ped$processed[ped$id %in% c(ped$dam[i], ped$sire[i])] | is.na(c(ped$dam[i], ped$sire[i])))) {
        ped$generation[i] <- max(dam_gen, sire_gen) + 1
        ped$processed[i] <- TRUE
      }
    }
  }

  ped <- ped[order(ped$generation, ped$id),]
  rownames(ped) <- NULL
  return(ped)
}

# Fonction pour charger et préparer les données
load_beekube_data <- function(json_file) {
  data <- fromJSON(json_file)

  criteres <- as.character(data$evaluate)
  criteres_eliminatoires <- as.character(data$evaluate_elimination)

  df <- as.data.frame(data$data)

  # Remplacer "Unknown" et NA par "0" pour le pédigrée
  df$queenbee_parent[df$queenbee_parent == "Unknown" | is.na(df$queenbee_parent)] <- "0"
  df$drone_parent[df$drone_parent == "Unknown" | is.na(df$drone_parent)] <- "0"

  # Convertir les autres colonnes en facteurs si nécessaire
  df <- df %>%
    mutate(across(c(queenbee, drone_parent, apiary, born), as.factor))

  list(df = df, criteres = criteres, criteres_eliminatoires = criteres_eliminatoires)
}

# Fonction pour calculer le BLUP
calculate_blup <- function(data, criteres) {
   blups <- list()
  methods <- list()

  # Préparer les données de pédigrée
  ped_data <- data.frame(
    id = as.character(data$queenbee),
    sire = as.character(data$drone_parent),
    dam = as.character(data$queenbee_parent)
  )

  # Nettoyer le pédigrée
  ped_data$sire[ped_data$sire == "0" | ped_data$sire == "Unknown" | is.na(ped_data$sire)] <- NA
  ped_data$dam[ped_data$dam == "0" | ped_data$dam == "Unknown" | is.na(ped_data$dam)] <- NA
  ped_data <- ped_data[ped_data$id != "0" & ped_data$id != "Unknown",]

  # Ajouter les parents manquants
  all_parents <- unique(c(ped_data$sire, ped_data$dam))
  missing_parents <- setdiff(all_parents, ped_data$id)
  missing_parents <- missing_parents[!is.na(missing_parents)]

  if(length(missing_parents) > 0) {
    ped_data <- rbind(ped_data,
                      data.frame(id = missing_parents,
                                 sire = NA,
                                 dam = NA))
  }

  # Trier le pédigrée
  sorted_ped <- sortPed(ped_data)

  # Créer l'objet pédigrée
  ped <- pedigree(sire = sorted_ped$sire, dam = sorted_ped$dam, label = sorted_ped$id)



  if(inherits(ped, "try-error")) {
    stop("Impossible de créer l'objet pédigrée. Vérifiez la structure du pédigrée.")
  }

  # Calculer la matrice de parenté
  A <- as(getA(ped), "matrix")

  for (critere in criteres) {
    tryCatch({
      # Préparer les données pour le modèle
      model_data <- data.frame(
        id = data$queenbee,
        value = data[[critere]],
        apiary = data$apiary,
        born = data$born
      )

      # Ajuster la matrice A pour correspondre aux données du modèle
      A_subset <- A[model_data$id, model_data$id]

      # Méthode 1: BLUP
      blup_model <- lmer(value ~ (1 | id) + (1 | apiary) + (1 | born),
                             data = model_data,
                             weights = diag(A_subset))

      if (!inherits(blup_model, "try-error")) {
        blups[[critere]] <- ranef(blup_model)$id[, 1]
        methods[[critere]] <- "BLUP"
      } else {
        # Méthode 2: Modèle bayésien
        # Extraire les constantes
        NUM_CHAINS <- 2
        NUM_ITER <- 2000
        WARMUP_ITER <- 1000
        FORMULA <- value ~ (1 | id)
        FAMILY <- gaussian()
        PRIORS <- c(prior(normal(0, 1), class = "Intercept"), prior(cauchy(0, 2), class = "sd"))
        COVARIATES <- list(id = A_subset)

        # Nom plus descriptif pour le modèle
        bayesian_model <-
          brm(
            FORMULA,
            data = model_data,
            family = FAMILY,
            prior = PRIORS,
            cov = COVARIATES,
            chains = NUM_CHAINS,
            iter = NUM_ITER,
            warmup = WARMUP_ITER
          )

        if (!inherits(bayesian_model, "try-error")) {
          blups[[critere]] <- ranef(bayesian_model)$id[, , "Intercept"]
          methods[[critere]] <- "Bayesien"
        } else {
          # Méthode 3: Régression linéaire simple
          lm_model <- lm(value ~ 1, data = model_data)
          blups[[critere]] <- rep(coef(lm_model)[1], nrow(model_data))
          methods[[critere]] <- "Régression linéaire"
        }
      }

      print(paste("BLUP calculé pour le critère", critere, "avec la méthode", methods[[critere]]))

    }, error = function(e) {
      message(paste("Impossible de calculer le BLUP pour le critère", critere, ":", e$message))
      print(str(model_data))
      blups[[critere]] <- rep(NA, nrow(data))
      methods[[critere]] <- "Échec"
    })
  }

  list(blups = blups, methods = methods)
}


# Utilisation
data <- load_beekube_data("/Users/tups/Sites/beekube-blup/lineage-blup.json")

# Afficher la structure des données pour le diagnostic
str(data$df)

results <- calculate_blup(data$df, data$criteres)

# Ajouter les BLUPs et les méthodes aux données originales
for (critere in data$criteres) {
  data$df[[paste0("blup_", critere)]] <- results$blups[[critere]]
  data$df[[paste0("method_", critere)]] <- results$methods[[critere]]
}

# Afficher les résultats
print(head(data$df, 10))

# Exporter les résultats en CSV
write.csv(data$df, "/Users/tups/Sites/beekube-blup/resultats_index.csv", row.names = FALSE)
