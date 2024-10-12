library(jsonlite)
library(dplyr)
library(pedigreemm)
library(Matrix)
library(lme4)
library(brms)
library(rstan)
library(lm.beta)

options(mc.cores = parallel::detectCores())

# Chemins des fichiers
current_dir <- getwd()
input_file <- file.path(current_dir, "data/input.json")
output_file <- file.path(current_dir, "data/resultats_index.json")

# Fonction pour trier correctement le pedigree
sort_pedigree <- function(ped, max_iterations = 1000) {
  ped$sorted <- FALSE
  sorted_ped <- data.frame()
  iteration <- 0

  while(nrow(ped) > 0 && iteration < max_iterations) {
    iteration <- iteration + 1

    # Trouver les individus sans parents ou dont les parents sont déjà triés
    to_add <- ped[is.na(ped$sire) & is.na(ped$dam) |
                  (is.na(ped$sire) | ped$sire %in% sorted_ped$id) &
                  (is.na(ped$dam) | ped$dam %in% sorted_ped$id), ]

    if(nrow(to_add) == 0) {
      # Diagnostic: identifier les individus problématiques
      problem_individuals <- ped[1:min(5, nrow(ped)), ]
      print("Individus problématiques:")
      print(problem_individuals)

      # Vérifier les cycles potentiels
      for(i in seq_len(nrow(problem_individuals))) {
        individual <- problem_individuals$id[i]
        sire <- problem_individuals$sire[i]
        dam <- problem_individuals$dam[i]

        if(!is.na(sire) && sire == individual) {
          stop(paste("Cycle détecté: L'individu", individual, "est son propre père"))
        }
        if(!is.na(dam) && dam == individual) {
          stop(paste("Cycle détecté: L'individu", individual, "est sa propre mère"))
        }
      }

      stop("Impossible de trier le pedigree. Il pourrait y avoir des cycles ou des incohérences.")
    }

    sorted_ped <- rbind(sorted_ped, to_add)
    ped <- ped[!ped$id %in% to_add$id, ]
  }

  if(iteration == max_iterations) {
    stop("Nombre maximum d'itérations atteint. Le pedigree pourrait être trop grand ou contenir des cycles.")
  }

  return(sorted_ped)
}

# Fonction pour vérifier et nettoyer le pedigree
check_and_clean_pedigree <- function(ped) {
  # Convertir toutes les colonnes en caractères
  ped$id <- as.character(ped$id)
  ped$sire <- as.character(ped$sire)
  ped$dam <- as.character(ped$dam)

  # Remplacer les valeurs vides par NA
  ped$sire[ped$sire == ""] <- NA
  ped$dam[ped$dam == ""] <- NA

  # Vérifier les doublons
  if(any(duplicated(ped$id))) {
    warning("Des IDs en double ont été détectés. Ils seront supprimés.")
    ped <- ped[!duplicated(ped$id), ]
  }

  # Vérifier les parents manquants
  missing_sires <- setdiff(ped$sire, c(ped$id, NA))
  missing_dams <- setdiff(ped$dam, c(ped$id, NA))

  if(length(c(missing_sires, missing_dams)) > 0) {
    warning("Certains parents ne sont pas présents dans la liste des IDs. Ils seront remplacés par NA.")
    ped$sire[ped$sire %in% missing_sires] <- NA
    ped$dam[ped$dam %in% missing_dams] <- NA
  }

  return(ped)
}

# Fonction récursive pour gérer les NA
replace_na_with_null <- function(x) {
  if (is.list(x)) {
    return(lapply(x, replace_na_with_null))
  } else if (is.vector(x) && !is.null(x)) {
    # Pour les vecteurs, on garde les NA tels quels
    return(x)
  } else {
    # Pour les valeurs individuelles, on retourne NULL si c'est NA
    return(if (is.na(x)) NULL else x)
  }
}

# fonction pour vérifier si un facteur a suffisamment de niveaux
has_enough_levels <- function(factor, min_levels = 2) {
  length(unique(factor)) >= min_levels
}

# Fonction pour calculer la matrice de parenté adaptée aux abeilles
calculate_bee_kinship <- function(ped, n_drones = 12, has_drone_info = TRUE) {
  # Vérifier si l'information sur les drones est disponible
  if (!has_drone_info) {
    # Si pas d'info sur les drones, utiliser une matrice de parenté standard
    A <- as(getA(ped), "matrix")
    return(A)
  }

  # Calcul basé sur Bienefeld et al. (2007) si l'info sur les drones est disponible
  A <- as(getA(ped), "matrix")

  # Ajustement pour l'accouplement multiple
  A_adjusted <- A * (1 + 1/(4 * n_drones))

  # Ajuster les coefficients de parenté entre reines soeurs
  queen_ids <- unique(ped$id[!is.na(ped$dam)])
  for (i in 1:(length(queen_ids) - 1)) {
    for (j in (i + 1):length(queen_ids)) {
      if (ped$dam[ped$id == queen_ids[i]] == ped$dam[ped$id == queen_ids[j]]) {
        A_adjusted[queen_ids[i], queen_ids[j]] <- A_adjusted[queen_ids[j], queen_ids[i]] <- 0.25 + 0.5 / n_drones
      }
    }
  }

  return(A_adjusted)
}

# Fonction pour charger et préparer les données
# Modification de la fonction load_beekube_data
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
  eval_data <- do.call(rbind, lapply(criteres, function(critere) {
    eval_list <- df[[critere]]
    do.call(rbind, lapply(seq_along(eval_list), function(i) {
      if (length(eval_list[[i]]) > 0) {
        data.frame(
          queenbee = df$queenbee[i],
          critere = critere,
          note = unlist(eval_list[[i]]$note),
          user = unlist(eval_list[[i]]$user),
          apiary = unlist(eval_list[[i]]$apiary),
          beehiveType = unlist(eval_list[[i]]$beehiveType),
          month = unlist(eval_list[[i]]$month),
          year = unlist(eval_list[[i]]$year),
          stringsAsFactors = FALSE
        )
      }
    }))
  }))


  # Convertir les colonnes en facteurs si nécessaire
  eval_data <- eval_data %>%
    mutate(across(c(queenbee, critere, user, apiary, beehiveType, month, year), as.factor))

  list(df = df, eval_data = eval_data, criteres = criteres, criteres_eliminatoires = criteres_eliminatoires)
}

# Fonction principale pour le calcul BLUP
calculate_blup <- function(data, eval_data, criteres, n_drones = 12) {
  blups <- list()
  methods <- list()

  # Préparer les données de pédigrée
  ped_data <- data.frame(
    id = as.character(data$queenbee),
    sire = as.character(data$drone_parent),
    dam = as.character(data$queenbee_parent)
  )

  # Nettoyer le pédigrée
  ped_data$sire[ped_data$sire == "0" |
                  ped_data$sire == "Unknown" |
                  is.na(ped_data$sire)] <- NA
  ped_data$dam[ped_data$dam == "0" |
                 ped_data$dam == "Unknown" |
                 is.na(ped_data$dam)] <- NA
  ped_data <- ped_data[ped_data$id != "0" & ped_data$id != "Unknown",]

  # Ajouter les parents manquants
  all_parents <- unique(c(ped_data$sire, ped_data$dam))
  missing_parents <- setdiff(all_parents, ped_data$id)
  missing_parents <- missing_parents[!is.na(missing_parents)]

  if (length(missing_parents) > 0) {
    ped_data <- rbind(ped_data,
                      data.frame(id = missing_parents,
                                 sire = NA,
                                 dam = NA))
  }

  ped_data <- check_and_clean_pedigree(ped_data)
  # Trier le pedigree
  sorted_ped <- sort_pedigree(ped_data)

  # Création de l'objet pedigree
  ped <- pedigree(sire = sorted_ped$sire, dam = sorted_ped$dam, label = sorted_ped$id)

    # Vérifier si l'information sur les drones est disponible
  has_drone_info <- !all(is.na(data$drone_parent) | data$drone_parent == "0")

  # Calcul de la matrice de parenté adaptée aux abeilles
  A <- calculate_bee_kinship(ped, n_drones, has_drone_info)

  for (critere in criteres) {
    tryCatch({
      model_data <- eval_data %>%
        filter(critere == !!critere) %>%
        left_join(data, by = "queenbee")

      if (nrow(model_data) < 1) {
        message(paste("Pas assez de données pour le critère", critere))
        blups[[critere]] <- rep(NA, nrow(data))
        methods[[critere]] <- "Pas assez de données"
        next
      }

      # Vérifier quelles colonnes sont disponibles
      available_columns <- names(model_data)

      # Construction de la formule du modèle
      random_effects <- c("user", "apiary", "beehiveType", "month", "year")
      formula_parts <- c("note ~ (1|queenbee) + (1|drone_parent)")
       for (effect in random_effects) {
        if (effect %in% available_columns && has_enough_levels(model_data[[effect]])) {
          formula_parts <- c(formula_parts, paste0("(1|", effect, ")"))
        }
      }
      formula_str <- paste(formula_parts, collapse = " + ")

      # Ajuster la matrice A pour correspondre aux données du modèle
      A_subset <- A[as.character(model_data$queenbee), as.character(model_data$queenbee)]

      # Méthode 1: BLUP avec lmer
      blup_model <- lmer(as.formula(formula_str), data = model_data, weights = diag(A_subset))


      if (!inherits(blup_model, "try-error")) {
        tryCatch({
          blup_values <- ranef(blup_model)$queenbee[, 1] + ranef(blup_model)$drone_parent[, 1]
          # Obtenir les noms des reines du modèle
          queen_names <- rownames(ranef(blup_model)$queenbee)
          names(blup_values) <- queen_names

          # Créer un vecteur aligné avec toutes les reines
          aligned_blups <- rep(NA, nrow(data))
          names(aligned_blups) <- as.character(data$queenbee)

          # Utiliser une méthode de correspondance plus robuste
        common_names <- intersect(names(aligned_blups), names(blup_values))
        aligned_blups[common_names] <- blup_values[common_names]

        blups[[critere]] <- aligned_blups
        methods[[critere]] <- "BLUP"
        }, error = function(e) {
          message(paste("Erreur lors de l'extraction des BLUPs pour le critère", critere, ":", e$message))
          blups[[critere]] <- rep(NA, nrow(data))
          methods[[critere]] <- "Échec (extraction)"
        })
      } else {
        # Méthode alternative si BLUP échoue
        lm_model <- lm(as.formula(paste(critere, "~ queenbee")), data = model_data)
        blup_values <- coef(lm_model)[-1]
        aligned_blups <- rep(NA, nrow(data))
        names(aligned_blups) <- as.character(data$queenbee)
        aligned_blups[names(blup_values)] <- blup_values
        blups[[critere]] <- aligned_blups
        methods[[critere]] <- "Régression linéaire"
      }
    }, error = function(e) {
      message(paste("Impossible de calculer le BLUP pour le critère", critere, ":", e$message))
      blups[[critere]] <- rep(NA, nrow(data))
      methods[[critere]] <- "Échec"
    })
  }

  list(blups = blups, methods = methods)
}

# Utilisation
data <- load_beekube_data(input_file)
results <- calculate_blup(data$df, data$eval_data, data$criteres)


# Préparation des résultats pour l'export JSON
export_list <- list()

for (i in seq_len(nrow(data$df))) {
  queen_data <- as.list(data$df[i, ])

  # Supprimer les colonnes de type liste
  queen_data <- queen_data[!sapply(queen_data, is.list)]

  # Ajouter les BLUPs et les méthodes
  blups <- list()
  methods <- list()
  for (critere in data$criteres) {
    blups[[critere]] <- results$blups[[critere]][i]
    methods[[critere]] <- results$methods[[critere]][i]
  }

  queen_data$blups <- blups
  queen_data$methods <- methods

  # Appliquer la fonction pour gérer les NA
  queen_data <- replace_na_with_null(queen_data)

  export_list[[i]] <- queen_data
}

# Convertir la liste en JSON
json_output <- toJSON(export_list, pretty = TRUE, auto_unbox = TRUE, na = "null")

# Exporter les résultats en JSON
write(json_output, output_file)

print("Les résultats ont été exportés en JSON avec les NA remplacés par null.")