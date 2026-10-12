suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(lme4)
})



METHOD_BLUP <- "BLUP"
METHOD_FALLBACK <- "lm"
METHOD_NO_DATA <- "no_data"
METHOD_FAILED <- "failed"

empty_heritability <- function() {
  list(heritability = NA_real_, se = NA_real_, v_additive_queen = NA_real_,
       v_additive_drone = NA_real_, v_colony = NA_real_, v_residual = NA_real_)
}

# Function to calculate bee-adapted heritability
calculate_bee_heritability <- function(data, trait, has_drone_info = TRUE) {
  required_columns <- c("queenbee", trait)
  if (has_drone_info) {
    required_columns <- c(required_columns, "drone_parent")
  }

  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0) {
    stop(paste("Missing columns in the data:", paste(missing_columns, collapse = ", ")))
  }

  # Check if 'apiary' is present, otherwise use a default value
  if (!"apiary" %in% names(data)) {
    warning("Column 'apiary' not found. Using a default value.")
    data$apiary <- factor(1)
  }

  data$queenbee <- factor(data$queenbee)
  if (has_drone_info) {
    data$drone_parent <- factor(data$drone_parent)
  }
  data$apiary <- factor(data$apiary)
  data[[trait]] <- as.numeric(data[[trait]])
  data <- data[complete.cases(data[, c(trait, "queenbee", "apiary")]), ]

  if (nrow(data) < 10) {  # arbitrary minimum number
    warning("Pas assez de données pour calculer l'héritabilité")
    return(empty_heritability())
  }

  tryCatch({
    formula <- if (has_drone_info) {
      as.formula(paste(trait, "~ (1|queenbee) + (1|drone_parent) + (1|apiary)"))
    } else {
      as.formula(paste(trait, "~ (1|queenbee) + (1|apiary)"))
    }
    model <- lme4::lmer(formula, data = data,
                        control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                                    check.nobs.vs.nRE = "ignore"))
    vc <- lme4::VarCorr(model)
    v_a <- as.numeric(vc$queenbee[1])
    v_c <- as.numeric(vc$apiary[1])
    v_d <- if (has_drone_info && !is.null(vc$drone_parent)) as.numeric(vc$drone_parent[1]) else 0
    v_e <- attr(vc, "sc")^2
    h2 <- (v_a + v_d) / (v_a + v_d + v_c + v_e)

    # Pas d'erreur standard sans profil de vraisemblance ou bootstrap : null (ROADMAP item 2)
    list(heritability = h2, se = NA_real_, v_additive_queen = v_a,
         v_additive_drone = v_d, v_colony = v_c, v_residual = v_e)
  }, error = function(e) {
    warning(paste("Error in calculating heritability:", e$message))
    empty_heritability()
  })
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

EVAL_COLUMNS <- c("queenbee", "critere", "note", "user", "apiary", "beehiveType", "month", "year")

empty_eval_data <- function() {
  data.frame(queenbee = character(0), critere = character(0), note = numeric(0),
             user = character(0), apiary = character(0), beehiveType = character(0),
             month = character(0), year = character(0), stringsAsFactors = FALSE)
}

# Notes des critères demandés, une ligne par note. Beekube envoie une clé (éventuellement
# vide) pour chaque critère non éliminatoire ; le code tolère aussi bien les clés absentes
# que les listes vides.
collect_evaluations <- function(df, criteres) {
  rows <- lapply(as.character(criteres), function(critere) {
    eval_list <- df$evaluate[[critere]]
    do.call(rbind, lapply(seq_along(eval_list), function(i) {
      ev <- eval_list[[i]]
      if (!is.list(ev) || NROW(ev) == 0) {
        return(NULL)
      }
      data.frame(
        queenbee = as.character(df$queenbee[i]),
        critere = critere,
        note = as.numeric(unlist(ev$note)),
        user = as.character(unlist(ev$user)),
        apiary = as.character(unlist(ev$apiary)),
        beehiveType = as.character(unlist(ev$beehiveType)),
        month = as.character(unlist(ev$month)),
        year = as.character(unlist(ev$year)),
        stringsAsFactors = FALSE
      )
    }))
  })
  out <- do.call(rbind, rows)
  if (is.null(out)) empty_eval_data() else out[, EVAL_COLUMNS]
}

# Reines éliminées : au moins une note sur un critère éliminatoire.
# NB : Beekube envoie aujourd'hui [] pour ces critères (EvaluateBLUP.php:84-88) : liste vide en pratique.
eliminated_queens <- function(df, criteres_eliminatoires) {
  if (length(criteres_eliminatoires) == 0) {
    return(character(0))
  }
  unique(as.character(collect_evaluations(df, criteres_eliminatoires)$queenbee))
}

# Function to load and prepare the data
load_beekube_data <- function(json_file) {
  data <- fromJSON(json_file)

  criteres <- as.character(data$evaluate)
  criteres_eliminatoires <- as.character(data$evaluate_elimination)

  df <- as.data.frame(data$data)

  # Unknown parents are coded "0" (character even when every value is null)
  df$queenbee_parent <- as.character(df$queenbee_parent)
  df$drone_parent <- as.character(df$drone_parent)
  df$queenbee_parent[df$queenbee_parent == "Unknown" | is.na(df$queenbee_parent)] <- "0"
  df$drone_parent[df$drone_parent == "Unknown" | is.na(df$drone_parent)] <- "0"

  eval_data <- collect_evaluations(df, criteres)

  df <- df %>%
    mutate(across(c(queenbee, drone_parent, born, variety), as.factor))

  list(df = df, eval_data = eval_data, criteres = criteres,
       criteres_eliminatoires = criteres_eliminatoires,
       eliminated = eliminated_queens(df, criteres_eliminatoires))
}

# Notes d'un critère jointes au seul champ utile de la reine (drone_parent).
# Sélection explicite : un champ `apiary` côté reine créerait apiary.x / apiary.y.
build_model_data <- function(eval_data, df, critere) {
  target <- as.character(critere)
  queens <- data.frame(
    queenbee = as.character(df$queenbee),
    drone_parent = as.character(df$drone_parent),
    stringsAsFactors = FALSE
  )
  eval_data[as.character(eval_data$critere) == target, , drop = FALSE] %>%
    mutate(queenbee = as.character(queenbee)) %>%
    left_join(queens, by = "queenbee") %>%
    mutate(across(c(queenbee, drone_parent, user, apiary, beehiveType, month, year), as.factor))
}

# L'effet drone_parent n'entre dans le modèle que s'il a au moins 2 niveaux parmi les reines notées
use_drone_effect <- function(model_data, has_drone_info) {
  isTRUE(has_drone_info) && has_enough_levels(model_data$drone_parent)
}

build_formula <- function(model_data, use_drone) {
  parts <- "note ~ (1|queenbee)"
  if (use_drone) {
    parts <- c(parts, "(1|drone_parent)")
  }
  for (effect in c("user", "apiary", "beehiveType", "month", "year")) {
    if (effect %in% names(model_data) && has_enough_levels(model_data[[effect]])) {
      parts <- c(parts, paste0("(1|", effect, ")"))
    }
  }
  paste(parts, collapse = " + ")
}

# BLUP lmer : effet reine + effet de son drone_parent (non pondéré)
fit_mixed_lmer <- function(model_data, use_drone) {
  model <- lme4::lmer(
    stats::as.formula(build_formula(model_data, use_drone)),
    data = model_data,
    control = lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  )
  effects <- lme4::ranef(model)
  queen_effects <- effects$queenbee
  values <- stats::setNames(queen_effects[, 1], rownames(queen_effects))

  # The drone effect is indexed by drone_parent level, not by queen
  drone_effects <- effects$drone_parent
  if (use_drone && !is.null(drone_effects)) {
    queen_drone <- as.character(model_data$drone_parent[match(names(values), as.character(model_data$queenbee))])
    drone_values <- drone_effects[queen_drone, 1]
    drone_values[is.na(drone_values)] <- 0
    values <- values + drone_values
  }
  values
}

# Repli : moyenne par reine centrée = coefficients de lm(note ~ 0 + queenbee) centrés
fit_fallback_lm <- function(model_data) {
  d <- model_data
  d$queenbee <- droplevels(factor(as.character(d$queenbee)))
  if (nlevels(d$queenbee) < 2) {
    stop("fallback lm needs at least 2 rated queens")
  }
  coefs <- stats::coef(stats::lm(note ~ 0 + queenbee, data = d))
  names(coefs) <- sub("^queenbee", "", names(coefs))
  coefs - mean(coefs, na.rm = TRUE)
}

# Calcule un critère et RENVOIE ses résultats (aucune affectation perdue dans un gestionnaire)
fit_criterion <- function(eval_data, df, critere, has_drone_info,
                          fit_mixed = fit_mixed_lmer, fit_fallback = fit_fallback_lm) {
  queen_ids <- as.character(df$queenbee)
  blups <- stats::setNames(rep(NA_real_, length(queen_ids)), queen_ids)
  methods <- stats::setNames(rep(METHOD_NO_DATA, length(queen_ids)), queen_ids)

  in_criterion <- as.character(eval_data$critere) == as.character(critere)
  rated <- intersect(queen_ids, unique(as.character(eval_data$queenbee[in_criterion])))
  if (length(rated) == 0) {
    message(sprintf("No data for criterion %s", critere))
    return(list(blups = blups, methods = methods, heritability = empty_heritability()))
  }

  tryCatch({
    model_data <- build_model_data(eval_data, df, critere)
    use_drone <- use_drone_effect(model_data, has_drone_info)

    method <- METHOD_BLUP
    values <- tryCatch(fit_mixed(model_data, use_drone), error = function(e) {
      message(sprintf("Mixed model failed for criterion %s: %s", critere, conditionMessage(e)))
      NULL
    })
    if (is.null(values)) {
      method <- METHOD_FALLBACK
      values <- tryCatch(fit_fallback(model_data), error = function(e) {
        message(sprintf("Fallback lm failed for criterion %s: %s", critere, conditionMessage(e)))
        NULL
      })
    }

    methods[rated] <- METHOD_FAILED
    if (!is.null(values)) {
      ok <- intersect(rated, names(values)[!is.na(values)])
      blups[ok] <- values[ok]
      methods[ok] <- method
    }

    list(blups = blups, methods = methods,
         heritability = calculate_bee_heritability(model_data, "note", use_drone))
  }, error = function(e) {
    message(sprintf("Unable to calculate the BLUP for criterion %s: %s", critere, conditionMessage(e)))
    methods[rated] <- METHOD_FAILED
    list(blups = stats::setNames(rep(NA_real_, length(queen_ids)), queen_ids),
         methods = methods, heritability = empty_heritability())
  })
}

# Main function for BLUP calculation
calculate_blup <- function(data, eval_data, criteres) {
  has_drone_info <- !all(is.na(data$drone_parent) | as.character(data$drone_parent) == "0")

  # Listes nommées même vides : toJSON les écrit {} et non []
  blups <- stats::setNames(list(), character(0))
  methods <- stats::setNames(list(), character(0))
  heritabilities <- stats::setNames(list(), character(0))

  for (critere in criteres) {
    fit <- fit_criterion(eval_data, data, critere, has_drone_info)
    blups[[critere]] <- fit$blups
    methods[[critere]] <- fit$methods
    heritabilities[[critere]] <- fit$heritability
  }

  list(blups = blups, methods = methods, heritabilities = heritabilities)
}


has_any_index <- function(blup_values, queen) {
  any(vapply(blup_values, function(x) !is.na(x[queen]), logical(1)))
}

# Reines candidates (mère comme partenaire) : vivantes, non éliminées, indexées sur au moins un critère
mating_candidates <- function(df, blup_values, eliminated = character(0)) {
  ids <- as.character(df$queenbee)
  died <- if (is.null(df$died)) rep(FALSE, length(ids)) else as.logical(df$died)
  died[is.na(died)] <- FALSE
  indexed <- vapply(ids, function(q) has_any_index(blup_values, q), logical(1), USE.NAMES = FALSE)
  ids[indexed & !died & !(ids %in% eliminated)]
}

# Function to suggest optimal matings based on genetic values and diversity
suggest_optimal_matings <- function(df, blup_results, eliminated = character(0), n_suggestions = 10) {
  queens <- as.character(df$queenbee)
  blup_values <- blup_results$blups
  candidates <- mating_candidates(df, blup_values, eliminated)

  # Create pedigree data
  ped_data <- data.frame(
    id = queens,
    sire = as.character(df$drone_parent),
    dam = as.character(df$queenbee_parent),
    stringsAsFactors = FALSE
  )

  # Unknown parents are coded "0" upstream: they must not count as a shared ancestor
  unknown_ids <- c("0", "Unknown", "")
  ped_data$sire[ped_data$sire %in% unknown_ids] <- NA
  ped_data$dam[ped_data$dam %in% unknown_ids] <- NA

  # Function to calculate genetic similarity
  calc_similarity <- function(queen1, queen2, ped_data, depth = 3) {
    # Fonction récursive pour trouver tous les ancêtres jusqu'à une certaine profondeur
    get_ancestors <- function(id, current_depth = 0) {
      if (current_depth >= depth || is.na(id)) {
        return(character(0))
      }
      parents <- c(
        ped_data$sire[ped_data$id == id],
        ped_data$dam[ped_data$id == id]
      )
      parents <- parents[!is.na(parents)]
      ancestors <- parents
      for (parent in parents) {
        ancestors <- c(ancestors, get_ancestors(parent, current_depth + 1))
      }
      unique(ancestors)
    }

    # Lignée de chaque reine : elle-même + ses ancêtres, pour qu'une mère
    # et sa fille partagent la mère
    ancestors1 <- unique(c(as.character(queen1), get_ancestors(queen1)))
    ancestors2 <- unique(c(as.character(queen2), get_ancestors(queen2)))

    common_ancestors <- intersect(ancestors1, ancestors2)
    total_ancestors <- union(ancestors1, ancestors2)
    if (length(total_ancestors) == 0) {
      return(1) # Diversité maximale si pas d'ancêtres communs connus
    }
    1 - length(common_ancestors) / length(total_ancestors)
  }

  empty_suggestions <- data.frame(
    mate = character(0), total_score = numeric(0), blup_score = numeric(0),
    diversity_score = numeric(0), queen_blup = numeric(0), mate_blup = numeric(0),
    stringsAsFactors = FALSE
  )

  suggestions <- list()
  for (queen in queens) {
    suggestions[[queen]] <- empty_suggestions
    if (!(queen %in% candidates)) {
      next
    }

    mate_details <- lapply(setdiff(candidates, queen), function(mate) {
      queen_blup <- mean(unlist(lapply(blup_values, function(x) x[queen])), na.rm = TRUE)
      mate_blup <- mean(unlist(lapply(blup_values, function(x) x[mate])), na.rm = TRUE)
      blup_score <- mean(c(queen_blup, mate_blup), na.rm = TRUE)

      # Normalisation du blup_score (plage initiale supposée : -0.5 à 0.5)
      blup_rescaled <- (blup_score + 0.5) / 1
      diversity_score <- calc_similarity(queen, mate, ped_data)
      total_score <- 0.7 * blup_rescaled + 0.3 * diversity_score

      data.frame(mate = mate, total_score = total_score, blup_score = blup_score,
                 diversity_score = diversity_score, queen_blup = queen_blup,
                 mate_blup = mate_blup, stringsAsFactors = FALSE)
    })

    if (length(mate_details) > 0) {
      mate_df <- do.call(rbind, mate_details)
      mate_df <- mate_df[order(mate_df$total_score, decreasing = TRUE), ]
      suggestions[[queen]] <- head(mate_df, n_suggestions)
    }
  }

  suggestions
}


# Convertit le data.frame de suggestions d'une reine en liste JSON
format_suggestions <- function(suggestions) {
  if (is.null(suggestions) || nrow(suggestions) == 0) {
    return(list())
  }
  lapply(seq_len(nrow(suggestions)), function(j) {
    list(
      queen = as.character(suggestions$mate[j]),
      total_score = round(suggestions$total_score[j], 3),
      blup_score = round(suggestions$blup_score[j], 3),
      diversity_score = round(suggestions$diversity_score[j], 3),
      queen_blup = round(suggestions$queen_blup[j], 3),
      mate_blup = round(suggestions$mate_blup[j], 3)
    )
  })
}

# Assemble la sortie JSON du contrat v1 : list(blup = [...], heritabilities = {...})
build_export <- function(data, results, mating_suggestions) {
  export_list <- list()
  for (i in seq_len(nrow(data$df))) {
    queen_data <- as.list(data$df[i, ])
    queen_id <- as.character(queen_data$queenbee)

    # Remove columns of type list
    queen_data <- queen_data[!sapply(queen_data, is.list)]
    queen_data <- replace_na_with_null(queen_data)

    # Objets JSON même sans critère : {} et non []
    blups <- stats::setNames(list(), character(0))
    methods <- stats::setNames(list(), character(0))
    for (critere in data$criteres) {
      blups[[critere]] <- unname(results$blups[[critere]][i])
      methods[[critere]] <- unname(results$methods[[critere]][i])
    }
    queen_data$blups <- blups
    queen_data$methods <- methods
    queen_data$mating_suggestions <- format_suggestions(mating_suggestions[[queen_id]])

    export_list[[i]] <- queen_data
  }
  list(blup = export_list, heritabilities = results$heritabilities)
}

# Calcule le BLUP v1 de input_file et écrit le JSON dans output_file
run_blup <- function(input_file, output_file) {
  data <- load_beekube_data(input_file)
  results <- calculate_blup(data$df, data$eval_data, data$criteres)
  mating_suggestions <- suggest_optimal_matings(data$df, results, data$eliminated)
  export <- build_export(data, results, mating_suggestions)
  json_output <- toJSON(export, pretty = TRUE, auto_unbox = TRUE, na = "null")
  write(json_output, output_file)
  invisible(export)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) != 2) {
    stop("Usage: Rscript blup.r <input.json> <output.json>")
  }
  run_blup(args[[1]], args[[2]])
  message("The results were exported to ", args[[2]])
}

if (sys.nframe() == 0L) {
  main()
}
