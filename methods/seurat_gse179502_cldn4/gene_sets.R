# Pre-specified mouse gene sets for GSE179502 (Cldn4-only).
# A8 Hallmark / MHC / TJ families are human; mapped to mm10 symbols.
# Cldn4 and Stk11 are held out of TJ. Tacstd2 is audit-only, never a gate.

AT2 <- c(
  "Sftpc", "Sftpb", "Sftpa1", "Lamp3", "Slc34a2", "Napsa", "Lyz2",
  "Abca3", "Cxcl15", "Hc", "Lpcat1", "Etv5", "Cebpa", "Nkx2-1"
)

IFN_CORE <- c(
  "Ifit1", "Ifit2", "Ifit3", "Isg15", "Mx1", "Stat1", "Irf7",
  "Cxcl10", "Bst2", "Rsad2"
)

MHC_CORE <- c(
  "H2-K1", "H2-D1", "H2-Q4", "B2m", "Tap1", "Tap2", "Psmb8", "Psmb9", "Nlrc5"
)

TJ_CORE <- c(
  "Cldn3", "Cldn7", "Cldn18", "Ocln", "Tjp1", "Tjp2", "F11r",
  "Marveld2", "Marveld3", "Cdh1", "Cgn", "Pard3"
)

HLA_TO_H2 <- list(
  "HLA-A" = "H2-K1",
  "HLA-B" = "H2-D1",
  "HLA-C" = "H2-Q4",
  "HLA-E" = "H2-T23",
  "HLA-F" = "H2-Q10",
  "HLA-G" = "H2-Q7",
  "HLA-DMA" = "H2-DMa",
  "HLA-DQA1" = "H2-Aa",
  "HLA-DRB1" = "H2-Eb1"
)

ALIASES <- list(
  OAS1 = c("Oas1a", "Oas1g"),
  OASL = c("Oasl1", "Oasl2"),
  OAS2 = "Oas2",
  OAS3 = "Oas3",
  MX1 = "Mx1",
  MX2 = "Mx2",
  IFI27 = c("Ifi27", "Ifi27l2a"),
  IFI44L = "Ifi44",
  WARS1 = "Wars",
  C1R = "C1ra",
  C1S = "C1s1",
  FCGR1A = "Fcgr1",
  SAMD9 = "Samd9l",
  SAMD9L = "Samd9l",
  MARCHF1 = c("Marchf1", "March1"),
  TENT5A = c("Tent5a", "Fam46a"),
  PALS1 = "Mpp5",
  PATJ = "Inadl",
  PNP = "Pnp",
  MT2A = "Mt2"
)

human_to_mouse_candidates <- function(symbol) {
  s <- trimws(symbol)
  if (!nzchar(s)) return(character())
  if (s %in% names(HLA_TO_H2)) return(unname(unlist(HLA_TO_H2[[s]])))
  if (s %in% names(ALIASES)) return(unname(unlist(ALIASES[[s]])))
  if (startsWith(s, "MIR")) return(character())
  if (grepl("-", s, fixed = TRUE) && identical(s, toupper(s))) {
    parts <- strsplit(s, "-", fixed = TRUE)[[1]]
    return(paste0(paste0(toupper(substring(parts[1], 1, 1)), tolower(substring(parts[1], 2))), "-", paste(parts[-1], collapse = "-")))
  }
  if (identical(s, toupper(s))) {
    return(paste0(toupper(substring(s, 1, 1)), tolower(substring(s, 2))))
  }
  paste0(toupper(substring(s, 1, 1)), substring(s, 2))
}

map_human_set <- function(human_genes, present) {
  out <- character()
  seen <- character()
  for (g in human_genes) {
    for (m in human_to_mouse_candidates(g)) {
      if (m %in% present && !(m %in% seen)) {
        seen <- c(seen, m)
        out <- c(out, m)
      }
    }
  }
  out
}

present_subset <- function(genes, present) genes[genes %in% present]

load_a8_mouse <- function(json_path, present) {
  a8 <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  sets <- a8$sets
  ifn_h <- sort(unique(c(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE, sets$HALLMARK_INTERFERON_ALPHA_RESPONSE)))
  mhc_h <- sets$CUSTOM_MHC_I_ANTIGEN_PRESENTATION
  tj_h <- sort(unique(c(sets$KEGG_TIGHT_JUNCTION, sets$GOBP_TIGHT_JUNCTION_ORGANIZATION)))
  foc <- a8$focal_genes
  if (!is.null(foc)) {
    krt <- if (!is.null(sets$KRT_EPITHELIAL)) sets$KRT_EPITHELIAL else character()
    tj_h <- unique(c(tj_h, setdiff(foc, krt)))
  }
  ifn <- map_human_set(ifn_h, present)
  mhc <- map_human_set(mhc_h, present)
  tj <- map_human_set(tj_h, present)
  tj <- setdiff(tj, c("Cldn4", "Stk11"))
  list(IFN_A8 = ifn, MHC_A8 = mhc, TJ_A8 = tj)
}
