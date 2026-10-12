# Tests R de blup.r (API v1).
# Lancement : docker compose build api && docker compose run --rm api Rscript tests/r/test_blup.R
suppressPackageStartupMessages(source("blup.r"))

# ---- Helpers ----
empty_named_list <- function() stats::setNames(list(), character(0))

rec <- function(note, user = 1L, month = "04", apiary = 100L) {
  list(note = note, user = user, apiary = apiary, beehiveType = "dadant-10",
       month = month, year = "24")
}

make_queen <- function(id, parent = NULL, drone = NULL, died = FALSE,
                       evaluate = empty_named_list()) {
  list(queenbee = id, queenbee_parent = parent, drone_parent = drone,
       apiary_default = 100L, born = "2023", evaluate = evaluate,
       variety = "Buckfast", died = died)
}

standard_queens <- function(eliminated_queen = NULL) {
  elim <- function(id) if (identical(id, eliminated_queen)) list(rec(1)) else list()
  list(
    make_queen(1L, evaluate = list(
      "10" = list(rec(0.25, 1L, "03"), rec(0.5, 2L, "04")), "11" = list(rec(0.5)), "20" = elim(1L))),
    make_queen(2L, parent = 1L, evaluate = list(
      "10" = list(rec(0.5, 1L, "03"), rec(0.75, 2L, "04")), "11" = list(rec(0.75)), "20" = elim(2L))),
    make_queen(3L, parent = 1L, evaluate = list(
      "10" = list(rec(0.75, 1L, "03"), rec(1, 2L, "04")), "11" = list(rec(1)), "20" = elim(3L))),
    make_queen(4L, parent = 2L, evaluate = list(
      "10" = list(rec(1, 1L, "03"), rec(0.75, 2L, "04")), "11" = list(), "20" = elim(4L))),
    make_queen(5L, parent = 2L, evaluate = list(
      "10" = list(), "11" = list(), "20" = elim(5L))),
    make_queen(6L, parent = 3L, died = TRUE, evaluate = list(
      "10" = list(rec(0.5, 1L, "03"), rec(0.5, 2L, "04")), "11" = list(), "20" = elim(6L)))
  )
}

write_input <- function(queens, evaluate = c(10L, 11L, 12L), evaluate_elimination = 20L) {
  path <- tempfile(fileext = ".json")
  payload <- list(exploitation = 1L, evaluate = I(evaluate),
                  evaluate_elimination = I(evaluate_elimination), data = queens)
  writeLines(jsonlite::toJSON(payload, auto_unbox = TRUE, null = "null", digits = NA), path)
  path
}

load_standard <- function(...) load_beekube_data(write_input(standard_queens(...)))

run_on <- function(queens, ...) {
  out <- tempfile(fileext = ".json")
  run_blup(write_input(queens, ...), out)
  jsonlite::fromJSON(out, simplifyVector = FALSE)
}

tests <- list()

# ---- Tests ----
tests$blup_r_is_sourceable <- function() {
  stopifnot(is.function(run_blup), is.function(main), is.function(build_export))
}

tests$run_blup_writes_output_file <- function() {
  out <- run_on(standard_queens())
  stopifnot(length(out$blup) == 6, !is.null(out$heritabilities))
}

tests$heavy_packages_not_loaded <- function() {
  stopifnot(!any(c("brms", "rstan", "lm.beta") %in% loadedNamespaces()))
}

tests$ranef_is_called_from_lme4 <- function() {
  code <- paste(readLines("blup.r"), collapse = "\n")
  stopifnot(!grepl("(^|[^:])ranef\\(", code), !grepl("(^|[^:])VarCorr\\(", code),
            !grepl("(^|[^:])lmer\\(", code))
}

# drone_parent renseigné (ids 1 à 5, toujours < i) : sans info mâle, l'ancien
# calculate_bee_kinship sautait la double boucle et le test passerait à tort
big_queens <- function(n = 2000L, rated = 50L) {
  lapply(seq_len(n), function(i) {
    parent <- if (i > 1L) as.integer(i %/% 2L) else NULL
    drone <- if (i > 10L) as.integer(i %% 5L + 1L) else NULL
    ev <- if (i > n - rated) {
      list("10" = list(rec(0.25 + 0.5 * (i %% 2), 1L, "03"), rec(0.5, 2L, "04")), "20" = list())
    } else {
      list("10" = list(), "20" = list())
    }
    make_queen(as.integer(i), parent = parent, drone = drone, evaluate = ev)
  })
}

tests$kinship_code_removed <- function() {
  stopifnot(!exists("calculate_bee_kinship"), !exists("sort_pedigree"),
            !exists("check_and_clean_pedigree"),
            !("pedigreemm" %in% loadedNamespaces()),
            !("pedigree" %in% names(formals(calculate_bee_heritability))))
}

tests$large_pedigree_runs_quickly <- function() {
  setTimeLimit(elapsed = 60, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf), add = TRUE)
  out <- run_on(big_queens(), evaluate = 10L)
  stopifnot(length(out$blup) == 2000)
}

tests$cyclic_pedigree_does_not_abort <- function() {
  qs <- list(
    make_queen(1L, parent = 2L, evaluate = list("10" = list(rec(0.5, 1L, "03"), rec(0.75, 2L, "04")), "20" = list())),
    make_queen(2L, parent = 1L, evaluate = list("10" = list(rec(1, 1L, "03"), rec(0.25, 2L, "04")), "20" = list()))
  )
  out <- run_on(qs, evaluate = 10L)
  stopifnot(length(out$blup) == 2)
}

tests$model_data_keeps_record_apiary <- function() {
  qs <- lapply(standard_queens(), function(q) {
    q$apiary <- 999L  # ancien format : rucher au niveau de la reine
    q
  })
  d <- load_beekube_data(write_input(qs))
  md <- build_model_data(d$eval_data, d$df, "10")
  stopifnot(nrow(md) == 10,
            "apiary" %in% names(md),
            !any(c("apiary.x", "apiary.y") %in% names(md)),
            all(as.character(md$apiary) == "100"),
            all(c("queenbee", "drone_parent", "note", "user", "beehiveType", "month", "year") %in% names(md)))
}

tests$missing_criterion_keys_are_tolerated <- function() {
  qs <- standard_queens()
  for (k in 4:6) qs[[k]]$evaluate[["11"]] <- NULL
  d <- load_beekube_data(write_input(qs))
  stopifnot(setequal(as.character(d$eval_data$queenbee[d$eval_data$critere == "11"]),
                     c("1", "2", "3")))
}

tests$no_records_at_all_runs <- function() {
  qs <- lapply(1:3, function(i) make_queen(i, evaluate = list("10" = list(), "20" = list())))
  out <- run_on(qs)
  stopifnot(length(out$blup) == 3,
            all(vapply(out$blup, function(q) is.null(q$blups[["10"]]), logical(1))))
}

fit_standard <- function(critere, ...) {
  d <- load_standard()
  fit_criterion(d$eval_data, d$df, critere, FALSE, ...)
}

tests$one_note_per_queen_gets_blup <- function() {
  f <- fit_standard("11")
  stopifnot(identical(unname(f$methods[c("1", "2", "3")]), rep("BLUP", 3)),
            !any(is.na(f$blups[c("1", "2", "3")])),
            all(f$methods[c("4", "5", "6")] == "no_data"),
            all(is.na(f$blups[c("4", "5", "6")])))
}

tests$criterion_without_records_is_no_data <- function() {
  f <- fit_standard("12")
  stopifnot(all(f$methods == "no_data"), all(is.na(f$blups)),
            is.na(f$heritability$heritability))
}

tests$methods_cover_every_queen <- function() {
  d <- load_standard()
  r <- calculate_blup(d$df, d$eval_data, d$criteres)
  stopifnot(identical(names(r$methods), c("10", "11", "12")),
            identical(names(r$blups), c("10", "11", "12")),
            identical(names(r$heritabilities), c("10", "11", "12")))
  for (k in names(r$methods)) {
    stopifnot(length(r$methods[[k]]) == 6, length(r$blups[[k]]) == 6,
              identical(names(r$methods[[k]]), as.character(1:6)),
              all(r$methods[[k]] %in% c("BLUP", "lm", "no_data", "failed")))
  }
  stopifnot(r$methods[["10"]][["5"]] == "no_data")
}

tests$fallback_lm_when_mixed_model_fails <- function() {
  f <- fit_standard("10", fit_mixed = function(...) stop("boom"))
  d <- load_standard()
  notes <- d$eval_data[d$eval_data$critere == "10", ]
  m <- tapply(notes$note, notes$queenbee, mean)
  expected <- m - mean(m)
  rated <- c("1", "2", "3", "4", "6")
  stopifnot(all(f$methods[rated] == "lm"), f$methods[["5"]] == "no_data",
            # as.numeric : tapply renvoie un tableau 1-d, all.equal comparerait aussi l'attribut dim
            isTRUE(all.equal(unname(f$blups[rated]), as.numeric(expected[rated]), tolerance = 1e-8)))
}

tests$failed_when_every_method_fails <- function() {
  boom <- function(...) stop("boom")
  f <- fit_standard("10", fit_mixed = boom, fit_fallback = boom)
  stopifnot(all(f$methods[c("1", "2", "3", "4", "6")] == "failed"),
            f$methods[["5"]] == "no_data", all(is.na(f$blups)))
}

tests$single_rated_queen_is_failed_not_dropped <- function() {
  qs <- list(
    make_queen(1L, evaluate = list("10" = list(rec(0.5, 1L, "03"), rec(0.75, 2L, "04")), "20" = list())),
    make_queen(2L, evaluate = list("10" = list(), "20" = list()))
  )
  d <- load_beekube_data(write_input(qs, evaluate = 10L))
  f <- fit_criterion(d$eval_data, d$df, "10", FALSE)
  stopifnot(f$methods[["1"]] == "failed", f$methods[["2"]] == "no_data")
}

tests$mixed_model_is_unweighted <- function() {
  d <- load_standard()
  md <- build_model_data(d$eval_data, d$df, "10")
  values <- fit_mixed_lmer(md, FALSE)
  ref <- lme4::lmer(stats::as.formula(build_formula(md, FALSE)), data = md,
                    control = lme4::lmerControl(check.nobs.vs.nlev = "ignore",
                                                check.nobs.vs.nRE = "ignore"))
  re <- lme4::ranef(ref)$queenbee
  stopifnot(isTRUE(all.equal(unname(values[rownames(re)]), re[, 1], tolerance = 1e-8)))
}

tests$drone_effect_needs_two_levels <- function() {
  one <- data.frame(drone_parent = factor(c("7", "7", "7")))
  two <- data.frame(drone_parent = factor(c("7", "8", "7")))
  stopifnot(!use_drone_effect(one, TRUE), use_drone_effect(two, TRUE),
            !use_drone_effect(two, FALSE),
            !grepl("drone_parent", build_formula(two, FALSE)),
            grepl("(1|drone_parent)", build_formula(two, TRUE), fixed = TRUE))
}

tests$heritability_se_is_null <- function() {
  md <- data.frame(
    queenbee = factor(rep(1:6, each = 2)),
    apiary = factor(rep(c(100, 200), times = 6)),
    note = c(0.2, 0.4, 0.5, 0.6, 0.7, 0.9, 0.3, 0.2, 0.8, 0.9, 0.5, 0.4)
  )
  h <- calculate_bee_heritability(md, "note", FALSE)
  stopifnot(identical(names(h), c("heritability", "se", "v_additive_queen",
                                  "v_additive_drone", "v_colony", "v_residual")),
            is.na(h$se), !is.na(h$heritability))
}

blups_fixture <- function() {
  list(
    "10" = stats::setNames(c(NA, 0.1, 0.2, NA, NA, 0.3), as.character(1:6)),
    "11" = stats::setNames(c(0.4, NA, NA, NA, NA, NA), as.character(1:6))
  )
}

tests$candidates_use_any_criterion <- function() {
  d <- load_standard()
  # 1 n'a d'index que sur le 2e critère ; 6 est morte ; 4 et 5 n'ont aucun index
  stopifnot(identical(mating_candidates(d$df, blups_fixture()), c("1", "2", "3")))
}

tests$candidates_exclude_eliminated <- function() {
  d <- load_standard()
  stopifnot(identical(mating_candidates(d$df, blups_fixture(), eliminated = "2"), c("1", "3")))
}

tests$elimination_records_mark_queen <- function() {
  d <- load_beekube_data(write_input(standard_queens(eliminated_queen = 2L)))
  stopifnot(identical(d$eliminated, "2"))
}

tests$empty_elimination_lists_mark_nobody <- function() {
  stopifnot(length(load_standard()$eliminated) == 0)
}

tests$suggestions_skip_eliminated_queen <- function() {
  d <- load_standard()
  s <- suggest_optimal_matings(d$df, list(blups = blups_fixture()), eliminated = "2")
  stopifnot(nrow(s[["2"]]) == 0, nrow(s[["6"]]) == 0,
            identical(s[["1"]]$mate, "3"), identical(s[["3"]]$mate, "1"))
}

tests$run_blup_never_suggests_eliminated_queen <- function() {
  out <- run_on(standard_queens(eliminated_queen = 2L))
  mates <- unlist(lapply(out$blup, function(q) {
    vapply(q$mating_suggestions, function(m) m$queen, character(1))
  }))
  stopifnot(length(out$blup[[2]]$mating_suggestions) == 0,
            !("2" %in% mates), length(mates) > 0)
}

tests$export_without_criteria_has_objects <- function() {
  out <- tempfile(fileext = ".json")
  run_blup(write_input(standard_queens(), evaluate = integer(0)), out)
  txt <- paste(readLines(out), collapse = "\n")
  stopifnot(grepl('"heritabilities":\\s*\\{\\s*\\}', txt),
            grepl('"blups":\\s*\\{\\s*\\}', txt),
            grepl('"methods":\\s*\\{\\s*\\}', txt))
}

# ---- Runner (laisser en fin de fichier ; insérer les nouveaux tests au-dessus) ----
run_all <- function(tests) {
  failures <- 0L
  for (name in names(tests)) {
    outcome <- tryCatch({
      tests[[name]]()
      "PASS"
    }, error = function(e) paste("FAIL:", conditionMessage(e)))
    cat(sprintf("%-50s %s\n", name, outcome))
    if (!identical(outcome, "PASS")) failures <- failures + 1L
  }
  if (failures > 0L) {
    cat(sprintf("%d R test(s) failed\n", failures))
    quit(status = 1)
  }
  cat("All R tests passed\n")
}

run_all(tests)
