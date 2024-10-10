
# install packages
############################################################
# install.packages("bsvars")")
# install.packages("bsvarSIGNs")

# Gross domestic product (GDP); Chain volume
rgdp_dwnld      = readrba::read_rba(series_id = "GGDPCVGDP")
rgdp_tmp        = xts::xts(rgdp_dwnld$value, rgdp_dwnld$date, tclass = 'yearqtr')
drgdp           = na.omit(400 * diff(log(rgdp_tmp)))
drgdp           = xts::to.quarterly(drgdp, OHLC = FALSE)

# Consumer price index; All groups; Quarterly change (in per cent)
picpi_dwnld     = readrba::read_rba(series_id = "GCPIAGSAQP")
pi              = 4 * xts::xts(picpi_dwnld$value, picpi_dwnld$date, tclass = 'yearqtr')
pi              = xts::to.quarterly(pi, OHLC = FALSE)

# Interbank Overnight Cash Rate
cr_dwnld        = readrba::read_rba(series_id = "FIRMMCRID")   # Cash Rate Target
cr_tmp          = xts::xts(cr_dwnld$value, cr_dwnld$date)
cr              = xts::to.quarterly(cr_tmp, OHLC = FALSE)

# Real Trade-Weighted Index
rtwi_dwnld      = readrba::read_rba(series_id = "FRERTWI")
rtwi_tmp        = xts::xts(rtwi_dwnld$value, rtwi_dwnld$date, tclass = 'yearqtr')
rtwi            = 100 * na.omit(diff(log(rtwi_tmp)))
rtwi            = xts::to.quarterly(rtwi, OHLC = FALSE)

y               = na.omit(merge(drgdp, pi, cr, rtwi))
y
plot(y, 
     main = "Australian monetary system",
     legend.loc = "bottomleft", 
     col = c("#FF00FF","#990099","#ff69b4","#330033")
)

# setup
############################################################
library(bsvars)
set.seed(123)

N       = ncol(y)
p       = 9
S_burn  = 1e4
S       = 2e4
thin    = 2

# structural matrix - extended model
############################################################
B_LR = matrix(TRUE, N, N)
B_LR[upper.tri(B_LR)] = FALSE
B_LR[3,4] = TRUE

# estimation - lower-triangular model
############################################################
spec_bsvar0     = specify_bsvar$new(as.matrix(y), p = p, stationary = rep(TRUE, N))
spec_bsvar0 |> 
  estimate(S = S_burn) |> 
  estimate(S = S, thin = thin) -> soe_bsvar0

ir_bsvar0       = compute_impulse_responses(soe_bsvar0, horizon = 20)
ir_bsvar0[,3,,] = 0.25 * ir_bsvar0[,3,,] / mean(ir_bsvar0[3,3,1,])
plot(ir_bsvar0)

# estimation - lower-triangular model - MLE prior for A
############################################################
spec_bsvar      = specify_bsvar$new(as.matrix(y), p = p, stationary = rep(TRUE, N))

A_mle           = t(solve(
  tcrossprod(spec_bsvar$data_matrices$X), 
  tcrossprod(spec_bsvar$data_matrices$X, spec_bsvar$data_matrices$Y)
))
spec_bsvar$prior$A = A_mle

spec_bsvar |> 
  estimate(S = S_burn) |> 
  estimate(S = S, thin = thin) -> soe_bsvar

ir_bsvar        = compute_impulse_responses(soe_bsvar, horizon = 20)
ir_bsvar[,3,,]  = 0.25 * ir_bsvar[,3,,] / mean(ir_bsvar[3,3,1,])
plot(ir_bsvar)

# estimation - extended model - MLE prior for A
############################################################
spec_bsvar_lr   = specify_bsvar$new(as.matrix(y), p = p, B = B_LR, stationary = rep(TRUE, N))
spec_bsvar_lr$prior$A = A_mle

spec_bsvar_lr |> 
  estimate(S = S_burn) |> 
  estimate(S = S, thin = thin) -> soe_bsvar_lr

ir_bsvar_lr        = compute_impulse_responses(soe_bsvar_lr, horizon = 20)
ir_bsvar_lr[,3,,]  = 0.25 * ir_bsvar_lr[,3,,] / mean(ir_bsvar_lr[3,3,1,])
plot(ir_bsvar_lr)

# estimation - lower-triangular SV heteroskedastic model
############################################################
spec_bsvar_sv = specify_bsvar_sv$new(as.matrix(y), p = p, stationary = rep(TRUE, N))
spec_bsvar_sv$prior$A = A_mle

spec_bsvar_sv |>
  estimate(S = S_burn) |> 
  estimate(S = S, thin = thin) -> soe_bsvar_sv

ir_bsvar_sv        = compute_impulse_responses(soe_bsvar_sv, horizon = 20)
ir_bsvar_sv[,3,,]  = 0.25 * ir_bsvar_sv[,3,,] / mean(ir_bsvar_sv[3,3,1,])
plot(ir_bsvar_sv)

# estimation - extended SV heteroskedastic model
############################################################
spec_bsvar_lr_sv = specify_bsvar_sv$new(as.matrix(y), p = p, B = B_LR, stationary = rep(TRUE, N))
spec_bsvar_lr_sv$prior$A = A_mle

spec_bsvar_lr_sv |>
  estimate(S = S_burn) |> 
  estimate(S = S, thin = thin) -> soe_bsvar_lr_sv

ir_bsvar_lr_sv        = compute_impulse_responses(soe_bsvar_lr_sv, horizon = 20)
ir_bsvar_lr_sv[,3,,]  = 0.25 * ir_bsvar_lr_sv[,3,,] / mean(ir_bsvar_lr_sv[3,3,1,])
plot(ir_bsvar_lr_sv)

# save the estimation results
############################################################
save(
  spec_bsvar0,
  spec_bsvar,
  spec_bsvar_lr,
  spec_bsvar_sv,
  spec_bsvar_lr_sv,
  soe_bsvar0,
  soe_bsvar,
  soe_bsvar_lr,
  soe_bsvar_sv,
  soe_bsvar_lr_sv,
  file = "soe_bsvar.rda"
)
