# Classify BayesPrism vs EcoTyper concordance patterns.
# Methods template — consumes already-computed per-cohort coefficients.
# Spec: bayesprism_ecotyper_vs_tacstd2_cldn4.md §B6.2 and §B7.
#
# Input rows (one per cohort x estimand x outcome):
#   cohort, seed (TACSTD2|CLDN4), estimand (E3|E5|E6|E7),
#   outcome (clr_cd8|tls|clr_ce1|clr_ce9),
#   beta, se, inside_empirical_null (logical)
#
# This function does not compute associations. It only assigns a pattern
# label from pre-declared rules so the label cannot be invented after seeing plots.

pattern_from_signs <- function(sign_e3_t, sign_e3_c, sign_e7_ce1, sign_e7_ce9,
                               e4_agree, either_inside_null) {
  if (isTRUE(either_inside_null)) return("P-null")
  seeds_agree <- isTRUE(e4_agree) && sign_e3_t == sign_e3_c && sign_e3_t != 0
  # CE1 expected same direction as "immune-cold" reading; CE9 opposite.
  # Direction itself is not pre-specified; only *mutual* alignment is.
  ce_aligned <- (sign_e7_ce1 != 0 && sign_e7_ce9 != 0 && sign_e7_ce1 == -sign_e7_ce9)
  if (seeds_agree && ce_aligned && sign_e3_t == sign_e7_ce1) return("P-agree")
  if (!seeds_agree && (sign_e3_t != 0 || sign_e3_c != 0)) return("P-gene-split")
  if (seeds_agree && !ce_aligned) return("P-method-split")
  "P-null"
}

sign_of <- function(beta, se, z_floor = 1.96) {
  if (!is.finite(beta) || !is.finite(se) || se <= 0) return(0L)
  z <- beta / se
  if (abs(z) < z_floor) return(0L)
  as.integer(sign(beta))
}
