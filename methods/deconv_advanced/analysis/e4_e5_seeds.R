# E4 seed concordance and E5 unique-information residuals.
# Methods template — no cohort I/O, no immune association.
# Confirm helpers against env/versions.lock before use.
# Spec: bayesprism_ecotyper_vs_tacstd2_cldn4.md §B4.4

lin_ccc <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  if (length(x) < 4L) return(NA_real_)
  mx <- mean(x); my <- mean(y)
  vx <- var(x); vy <- var(y)
  sxy <- cov(x, y)
  2 * sxy / (vx + vy + (mx - my)^2)
}

e4_concordance <- function(e3_t, e3_c) {
  ok <- is.finite(e3_t) & is.finite(e3_c)
  list(
    spearman = if (sum(ok) >= 4L)
      suppressWarnings(cor(e3_t[ok], e3_c[ok], method = "spearman"))
    else NA_real_,
    ccc = lin_ccc(e3_t, e3_c),
    n = sum(ok)
  )
}

e5_residuals <- function(e3_t, e3_c) {
  ok <- is.finite(e3_t) & is.finite(e3_c)
  if (sum(ok) < 4L) {
    return(list(
      tacstd2_given_cldn4 = rep(NA_real_, length(e3_t)),
      cldn4_given_tacstd2 = rep(NA_real_, length(e3_c))
    ))
  }
  r_t <- r_c <- rep(NA_real_, length(e3_t))
  r_t[ok] <- resid(lm(e3_t[ok] ~ e3_c[ok]))
  r_c[ok] <- resid(lm(e3_c[ok] ~ e3_t[ok]))
  list(
    tacstd2_given_cldn4 = r_t,
    cldn4_given_tacstd2 = r_c
  )
}
