#' --------------------
#' Helpers specific to working
#' with sdmTMB and GAMs
#' --------------------


####
## Create basis set and add to data
## See extensive conversation at https://github.com/sdmTMB/sdmTMB/issues/509
####

# Build smoother basis and split into:
# Xs: unpenalized component
# Zs: penalized component

make_basis <- function(smooth_formula, dat, basis_prev = NULL){
  sm <- sdmTMB:::parse_smoothers(smooth_formula, data = dat, basis_prev = basis_prev)
  
  sx <- as.data.frame(sm$Xs)
  colnames(sx) <- paste0("SX", seq_len(ncol(sx)))
  sz <- as.data.frame(do.call(cbind, sm$Zs))
  colnames(sz) <- paste0("SZ", seq_len(ncol(sz)))
  v <- as.data.frame(sm$basis_out[[1]][[1]]$X)
  colnames(v) <- paste0("SV", seq_len(ncol(v)))
  dat <- cbind(dat, sx, sz, v)
  
  return(list(sx = sx, sz = sz, dat = dat, sm = sm, v = v))
}
