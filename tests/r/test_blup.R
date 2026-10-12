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
