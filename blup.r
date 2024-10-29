library(jsonlite)
library(dplyr)
library(pedigreemm)
library(Matrix)
library(lme4)
library(brms)
library(rstan)
library(lm.beta)

options(mc.cores = parallel::detectCores())

if (requireNamespace("rstan", quietly = TRUE)) {
  rstan::rstan_options(auto_write = TRUE)
  rstan::rstan_options(threads_per_chain = 1)
}

# File Paths
current_dir <- getwd()
input_file <- file.path(current_dir, "data/input.json")
output_file <- file.path(current_dir, "data/resultats_index.json")


# Function to calculate bee-adapted heritability
calculate_bee_heritability <- function(data, trait, pedigree, n_drones = 12, has_drone_info = TRUE) {
  # Vérifier la présence des colonnes nécessaires
  required_columns <- c("queenbee", trait)
  if (has_drone_info) {
    required_columns <- c(required_columns, "drone_parent")
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(paste("Colonnes manquantes dans les données:", paste(missing_columns, collapse = ", ")))
  }

  # Vérifier si 'apiary' est présent, sinon utiliser une valeur par défaut
  if (!"apiary" %in% names(data)) {
    warning("Colonne 'apiary' non trouvée. Utilisation d'une valeur par défaut.")
    data$apiary <- factor(1)  # Assigner une valeur par défaut
  }

  # Ajuster la matrice de parenté pour les abeilles
  A <- as(getA(pedigree), "CsparseMatrix") # Utilise la nouvelle syntaxe recommandée
  A_adjusted <- if (has_drone_info) {
    as.matrix(A) * (1 + 1 / (4 * n_drones)) # Convertit explicitement en matrice dense
  } else {
    as.matrix(A)
  }

  # Assure-toi que les poids sont un vecteur numérique
  weights <- diag(A_adjusted)

  # print(A_adjusted);

  # Préparer la formule du modèle
  formula <- if (has_drone_info) {
    as.formula(paste(trait, "~ (1|queenbee) + (1|drone_parent) + (1|apiary)"))
  } else {
    as.formula(paste(trait, "~ (1|queenbee) + (1|apiary)"))
  }

  tryCatch({
    # Ajuster le modèle mixte
    model <- lmer(formula,
                  data = data,
                  control = lmerControl(check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.nRE = "ignore"),
                  weights = weights)

    # Extraire les composantes de la variance
    vc <- VarCorr(model)
    v_a <- vc$queenbee[1]  # Variance génétique additive de la reine
    v_c <- vc$apiary[1]  # Variance due à l'effet de la colonie/rucher
    v_d <- if (has_drone_info) vc$drone_parent[1] else 0  # Variance due aux drones
    v_e <- attr(vc, "sc")^2  # Variance résiduelle

    # Calculer l'héritabilité
    h2 <- (v_a + v_d) / (v_a + v_d + v_c + v_e)

    # Calculer l'erreur standard de l'héritabilité
    se_h2 <- sqrt(var(h2))

    return(list(heritability = h2, se = se_h2,
                v_additive_queen = v_a, v_additive_drone = v_d,
                v_colony = v_c, v_residual = v_e))

  }, error = function(e) {
    message(paste("Erreur dans calculate_bee_heritability:", e$message))
    return(list(heritability = NA, se = NA,
                v_additive_queen = NA, v_additive_drone = NA,
                v_colony = NA, v_residual = NA))
  })
}

# Function to properly sort the pedigree
sort_pedigree <- function(ped, max_iterations = 1000) {
  ped$sorted <- FALSE
  sorted_ped <- data.frame()
  iteration <- 0

  while (nrow(ped) > 0 && iteration < max_iterations) {
    iteration <- iteration + 1

    # Find individuals without parents or whose parents are already sorted
    to_add <- ped[is.na(ped$sire) & is.na(ped$dam) |
                    (is.na(ped$sire) | ped$sire %in% sorted_ped$id) &
                      (is.na(ped$dam) | ped$dam %in% sorted_ped$id),]

    if (nrow(to_add) == 0) {
      # Diagnosis: Identify problematic individuals
      problem_individuals <- ped[1:min(5, nrow(ped)),]
      print("Individus problématiques:")
      print(problem_individuals)

      # Verify potential cycles
      for (i in seq_len(nrow(problem_individuals))) {
        individual <- problem_individuals$id[i]
        sire <- problem_individuals$sire[i]
        dam <- problem_individuals$dam[i]

        if (!is.na(sire) && sire == individual) {
          stop(paste("Cycle detected: The individual", individual, "is his own father"))
        }
        if (!is.na(dam) && dam == individual) {
          stop(paste("Cycle detected: The individual", individual, "is their own mother"))
        }
      }

      stop("Unable to sort the pedigree. There might be cycles or inconsistencies.")
    }

    sorted_ped <- rbind(sorted_ped, to_add)
    ped <- ped[!ped$id %in% to_add$id,]
  }

  if (iteration == max_iterations) {
    stop("Maximum number of iterations reached. The pedigree might be too large or contain cycles.")
  }

  return(sorted_ped)
}

# Function to check and clean the pedigree
check_and_clean_pedigree <- function(ped) {
  # Convert all columns to characters
  ped$id <- as.character(ped$id)
  ped$sire <- as.character(ped$sire)
  ped$dam <- as.character(ped$dam)

  # Replace empty values with NA
  ped$sire[ped$sire == ""] <- NA
  ped$dam[ped$dam == ""] <- NA

  # Check for duplicates
  if (any(duplicated(ped$id))) {
    warning("Duplicate IDs have been detected. They will be removed.")
    ped <- ped[!duplicated(ped$id),]
  }

  # Check for missing parents
  missing_sires <- setdiff(ped$sire, c(ped$id, NA))
  missing_dams <- setdiff(ped$dam, c(ped$id, NA))

  if (length(c(missing_sires, missing_dams)) > 0) {
    warning("Some parents are not present in the list of IDs. They will be replaced by NA.")
    ped$sire[ped$sire %in% missing_sires] <- NA
    ped$dam[ped$dam %in% missing_dams] <- NA
  }

  return(ped)
}

# Recursive function to handle NA
replace_na_with_null <- function(x) {
  if (is.list(x)) {
    return(lapply(x, replace_na_with_null))
  } else if (is.vector(x) && !is.null(x)) {
    # For vectors, keep the NAs as they are
    return(x)
  } else {
    # For individual values, return NULL if it is NA
    return(if (is.na(x)) NULL else x)
  }
}

# function to check if a factor has enough levels
has_enough_levels <- function(factor, min_levels = 2) {
  length(unique(factor)) >= min_levels
}

# Function to calculate the kinship matrix adapted for bees
calculate_bee_kinship <- function(ped, n_drones = 12, has_drone_info = TRUE) {
  # Check if information about drones is available
  if (!has_drone_info) {
    # If no information on drones, use a standard kinship matrix
    A <- as(getA(ped), "matrix")
    return(A)
  }

  # Calculation based on Bienefeld et al. (2007) if drone information is available
  A <- as(getA(ped), "matrix")

  # Adjustment for multiple mating
  A_adjusted <- A * (1 + 1 / (4 * n_drones))

  # Adjust the kinship coefficients between sister queens
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

# Function to load and prepare the data
# Modification of the function load_beekube_data
load_beekube_data <- function(json_file) {
  data <- fromJSON(json_file)

  criteres <- as.character(data$evaluate)
  criteres_eliminatoires <- as.character(data$evaluate_elimination)

  df <- as.data.frame(data$data)

  # Replace "Unknown" and NA with "0" for pedigree
  df$queenbee_parent[df$queenbee_parent == "Unknown" | is.na(df$queenbee_parent)] <- "0"
  df$drone_parent[df$drone_parent == "Unknown" | is.na(df$drone_parent)] <- "0"

  # Convert other columns to factors if necessary
  df <- df %>%
    mutate(across(c(queenbee, drone_parent, born), as.factor))

  # Prepare evaluation data
  eval_data <- do.call(rbind, lapply(criteres, function(critere) {
    eval_list <- df$evaluate[[critere]]
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


  # Convert columns to factors if necessary
  eval_data <- eval_data %>%
    mutate(across(c(queenbee, critere, user, apiary, beehiveType, month, year), as.factor))

  list(df = df, eval_data = eval_data, criteres = criteres, criteres_eliminatoires = criteres_eliminatoires)
}

# Main function for BLUP calculation
calculate_blup <- function(data, eval_data, criteres, n_drones = 12) {
  blups <- list()
  methods <- list()
  heritabilities <- list()

  # Prepare pedigree data
  ped_data <- data.frame(
    id = as.character(data$queenbee),
    sire = as.character(data$drone_parent),
    dam = as.character(data$queenbee_parent)
  )

  # Clean the pedigree
  ped_data$sire[ped_data$sire == "0" |
                  ped_data$sire == "Unknown" |
                  is.na(ped_data$sire)] <- NA
  ped_data$dam[ped_data$dam == "0" |
                 ped_data$dam == "Unknown" |
                 is.na(ped_data$dam)] <- NA
  ped_data <- ped_data[ped_data$id != "0" & ped_data$id != "Unknown",]

  # Add missing parents
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
  # Sort the pedigree
  sorted_ped <- sort_pedigree(ped_data)

  # Creation of the pedigree object
  ped <- pedigree(sire = sorted_ped$sire, dam = sorted_ped$dam, label = sorted_ped$id)


  # Check if information about drones is available
  has_drone_info <- !all(is.na(data$drone_parent) | data$drone_parent == "0")

  # Calculation of the kinship matrix adapted for bees
  A <- calculate_bee_kinship(ped, n_drones, has_drone_info)

  for (critere in criteres) {
    tryCatch({
      model_data <- eval_data %>%
        filter(critere == !!critere) %>%
        left_join(data, by = "queenbee")

      if (nrow(model_data) < 1) {
        message(paste("Not enough data for the criterion", critere))
        blups[[critere]] <- rep(NA, nrow(data))
        methods[[critere]] <- "Not enough data"
        next
      }

      # Check which columns are available
      available_columns <- names(model_data)

      # Constructing the model formula
      random_effects <- c("queenbee", "user", "apiary", "beehiveType", "month", "year", "drone_parent")

      if (has_drone_info) {
        formula_parts <- c("note ~ (1|queenbee) + (1|drone_parent)")
      } else {
        formula_parts <- c("note ~ (1|queenbee)")
      }

      for (effect in random_effects[-1]) {
        if (effect %in% available_columns && has_enough_levels(model_data[[effect]])) {
          formula_parts <- c(formula_parts, paste0("(1|", effect, ")"))
        }
      }
      formula_str <- paste(formula_parts, collapse = " + ")

      # Adjust the matrix A to match the model data
      A_subset <- A[as.character(model_data$queenbee), as.character(model_data$queenbee)]

      # Method 1: BLUP with lmer
      blup_model <- lmer(as.formula(formula_str), data = model_data, weights = diag(A_subset))


      if (!inherits(blup_model, "try-error")) {
        tryCatch({
          if (has_drone_info) {
            blup_values <- ranef(blup_model)$queenbee[, 1] + ranef(blup_model)$drone_parent[, 1]
          } else {
            blup_values <- ranef(blup_model)$queenbee[, 1]
          }

          # Get the names of the queens from the model
          queen_names <- rownames(ranef(blup_model)$queenbee)
          names(blup_values) <- queen_names

          # Create an aligned vector with all the queens
          aligned_blups <- rep(NA, nrow(data))
          names(aligned_blups) <- as.character(data$queenbee)

          # Use a more robust matching method
          common_names <- intersect(names(aligned_blups), names(blup_values))
          aligned_blups[common_names] <- blup_values[common_names]

          blups[[critere]] <- aligned_blups
          methods[[critere]] <- "BLUP"
        }, error = function(e) {
          message(paste("Error extracting BLUPs for the criterion", critere, ":", e$message))
          blups[[critere]] <- rep(NA, nrow(data))
          methods[[critere]] <- "Failure (extraction)"
        })
      } else {
        # Alternative method if BLUP fails
        lm_model <- lm(note ~ queenbee, data = model_data)
        blup_values <- coef(lm_model)[-1]
        aligned_blups <- rep(NA, nrow(data))
        names(aligned_blups) <- as.character(data$queenbee)
        aligned_blups[names(blup_values)] <- blup_values
        blups[[critere]] <- aligned_blups
        methods[[critere]] <- "Linear regression"
      }

      # Calcul de l'héritabilité
      h2_result <- calculate_bee_heritability(model_data, "note", ped, n_drones, has_drone_info)
      heritabilities[[critere]] <- list(
        heritability = h2_result$heritability,
        se = h2_result$se,
        v_additive_queen = h2_result$v_additive_queen,
        v_additive_drone = h2_result$v_additive_drone,
        v_colony = h2_result$v_colony,
        v_residual = h2_result$v_residual
      )

    }, error = function(e) {
      message(paste("Unable to calculate the BLUP for the criterion", critere, ":", e$message))
      blups[[critere]] <- rep(NA, nrow(data))
      methods[[critere]] <- "Failure"
    })
  }

  list(blups = blups, methods = methods, heritabilities = heritabilities)
}


data <- load_beekube_data(input_file)
results <- calculate_blup(data$df, data$eval_data, data$criteres)


# Preparing results for JSON export
export_list <- list()

for (i in seq_len(nrow(data$df))) {
  queen_data <- as.list(data$df[i,])

  # Remove columns of type list
  queen_data <- queen_data[!sapply(queen_data, is.list)]

  # Add the BLUPs and methods
  blups <- list()
  methods <- list()
  for (critere in data$criteres) {
    blups[[critere]] <- results$blups[[critere]][i]
    methods[[critere]] <- results$methods[[critere]][i]
  }

  queen_data$blups <- blups
  queen_data$methods <- methods

  # Apply the function to handle NAs
  queen_data <- replace_na_with_null(queen_data)

  export_list[[i]] <- queen_data
}


# Preparing results for JSON export
export <- list()

# Ajout des héritabilités aux résultats
export$blup <- export_list
export$heritabilities <- results$heritabilities

# Convert the list to JSON
json_output <- toJSON(export, pretty = TRUE, auto_unbox = TRUE, na = "null")

# Export results to JSON
write(json_output, output_file)

print("The results were exported to JSON with NAs replaced by null.")