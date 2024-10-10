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
    mutate(across(c(queenbee, drone_parent, born), as.factor))

  # Préparer les données d'évaluation
  eval_data <- lapply(criteres, function(critere) {
  eval_list <- df[[critere]]
  if (length(eval_list) > 0) {
    eval_df <- do.call(rbind, lapply(seq_along(eval_list), function(i) {
      if (length(eval_list[[i]]) > 0) {
        data.frame(
          queenbee = df$queenbee[i],
          critere = critere,
          note = eval_list[[i]]$note,
          user = eval_list[[i]]$user,
          apiary = eval_list[[i]]$apiary,
          beehiveType = eval_list[[i]]$beehiveType,
          month = eval_list[[i]]$month,
          year = eval_list[[i]]$year,
          stringsAsFactors = FALSE
        )
      }
    }))
    return(eval_df)
  }
})

  print(eval_data);

  # Combiner tous les data frames en un seul
  eval_data <- do.call(rbind, eval_data)

  # Convertir les colonnes en facteurs si nécessaire
  eval_data <- eval_data %>%
    mutate(across(c(queenbee, critere, user, apiary, beehiveType, month, year), as.factor))

  # Ajouter cette ligne pour vérifier la structure des données
  print(str(eval_data))

  list(df = df, eval_data = eval_data, criteres = criteres, criteres_eliminatoires = criteres_eliminatoires)
}

# Fonction pour calculer le BLUP
# Fonction pour calculer le BLUP
calculate_blup <- function(data, eval_data, criteres) {
  blups <- list()
  methods <- list()

  for (critere in criteres) {
    tryCatch({
      # Préparer les données pour le modèle
      model_data <- eval_data %>%
        filter(critere == !!critere) %>%
        left_join(data, by = "queenbee")

      # Vérifier quelles colonnes sont disponibles
      available_columns <- names(model_data)
      print(paste("Colonnes disponibles pour le critère", critere, ":", paste(available_columns, collapse = ", ")))

      # Construire la formule du modèle en fonction des colonnes disponibles
      random_effects <- c("queenbee", "user", "apiary", "beehiveType", "month", "year")
      formula_parts <- c("note ~ (1|queenbee)")
      for (effect in random_effects[-1]) {
        if (effect %in% available_columns) {
          formula_parts <- c(formula_parts, paste0("(1|", effect, ")"))
        }
      }
      formula_str <- paste(formula_parts, collapse = " + ")
      print(paste("Formule du modèle pour le critère", critere, ":", formula_str))

      # Ajuster la matrice A pour correspondre aux données du modèle
      A_subset <- A[model_data$queenbee, model_data$queenbee]

      # Méthode 1: BLUP avec lmer
      blup_model <- lmer(as.formula(formula_str), data = model_data, weights = diag(A_subset))

      if (!inherits(blup_model, "try-error")) {
        blups[[critere]] <- ranef(blup_model)$queenbee[, 1]
        methods[[critere]] <- "BLUP"
      } else {
        # Méthode 2: Régression linéaire simple
        lm_model <- lm(note ~ queenbee, data = model_data)
        blups[[critere]] <- coef(lm_model)[-1]  # Exclure l'intercept
        methods[[critere]] <- "Régression linéaire"
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
data <- load_beekube_data("/Users/tups/Sites/beekube-blup/lineage.json")

# Afficher la structure des données pour le diagnostic


results <- calculate_blup(data$df, data$eval_data, data$criteres)


# Préparation des résultats pour l'export CSV
export_df <- data$df

# Supprimer les colonnes de type liste
list_columns <- sapply(export_df, is.list)
export_df <- export_df[, !list_columns]

# Ajouter les BLUPs et les méthodes
for (critere in data$criteres) {
  export_df[[paste0("blup_", critere)]] <- results$blups[[critere]]
  export_df[[paste0("method_", critere)]] <- results$methods[[critere]]
}

# Exporter les résultats en CSV
write.csv(export_df, "/Users/tups/Sites/beekube-blup/resultats_index.csv", row.names = FALSE)