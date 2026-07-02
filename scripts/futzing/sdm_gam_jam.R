library(mgcv)

sdm_time <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                               0, 
                             spatial = "off",
                   time_varying_type = "rw0",
                   time_varying =~1,
                   knots = list(year = 5),
                           #  mesh = basic_mesh,
                             spatiotemporal = "off",
                             time = "year",
                             data = nwa_dat,
                             extra_time = c(1950, 1951, 1952, 1953, 1954, 1955, 1956,
                                            1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965,
                                            1966, 1967, 1969, 1970, 1971, 1972, 1974, 1976))

gam_time <- bam(ln_focal_std_by_ecoregion ~ 
                       s(year) ,
                     data = nwa_dat)

modelbased::estimate_relation(gam_time, 
                              by = "year",
                              length = 200) |> plot()

p <- predict(sdm_time, 
             newdata = data.frame(year = 1970:2025),
            return_tmb_object = TRUE)

p1 <- get_index(p, bias_correct = TRUE)
ggplot(p1, aes(x = year, y = est)) + geom_line()

pred <- predict(sdm_time, 
                newdata = data.frame(year = 1970:2025),
                se.fit = TRUE)

pred$gam <- predict(gam_time, newdata = data.frame(year = 1970:2025))

ggplot(pred, aes(x = year, y = est)) + geom_line() + 
  geom_line(aes(y = gam), color = "red")



sdm_gam <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                     s(year, bs = "tp"), 
                   spatial = "on",
                  mesh = basic_mesh,
                  spatial_varying = ~ 0 + year,
                   data = nwa_dat)

###
sm <- parse_smoothers(ln_focal_std_by_ecoregion ~ s(year, bs = "cr", k = 10), nwa_dat) 

sm$basis_out[[1]][[1]] |> names()
dim(sm$basis_out[[1]][[1]]$X)
matplot(t(sm$basis_out[[1]][[1]]$X), type = "l")

#1 s(year) CRS   NA    1     -0.0358   2008
### 
b <- basis(s(year, bs = "cr", k = 10), 
           data = nwa_dat,
           constraint = TRUE)

ggplot(b, aes(x = year, y = .value, color = .bf)) + 
  geom_line()

smat <- sm$basis_out[[1]][[1]]$X |> as_tibble()
nwa_dat <- bind_cols(nwa_dat, smat)


gam_time <- gam(ln_focal_std_by_ecoregion ~ 
                  s(year, bs = "cr", k = 10) ,
                data = nwa_dat)


basis_time <- gam(ln_focal_std_by_ecoregion ~ 
                  V1+V2+V3+V4+V5+V6+V7+V8+V9,
                data = nwa_dat)

plot(predict(gam_time), predict(basis_time))

pred_dat <- basis(s(year, bs = "cr", k = 10), 
                       data = nwa_dat, 
                       at = data.frame(year = 1950:2025),
                  constraint = TRUE) |>
  tidyr::pivot_wider(names_from = ".bf", 
                     values_from = ".value", 
              names_prefix = "V")

plot(predict(gam_time, newdata = pred_dat), 
     predict(basis_time, newdata = pred_dat))


sdm_st <- sdmTMB(ln_focal_std_by_ecoregion ~ 
                   1+V1+V2+V3+V4+V5+V6+V7+V8+V9, 
                   spatial = "on",
                   spatial_varying = 1+V1+V2+V3+V4+V5+V6+V7+V8+V9,
                   mesh = basic_mesh,
                   spatiotemporal = "off",
                   time = "year",
                   data = nwa_dat,
                   extra_time = c(1950, 1951, 1952, 1953, 1954, 1955, 1956,
                                  1957, 1958, 1959, 1960, 1961, 1962, 1963, 1964, 1965,
                                  1966, 1967, 1969, 1970, 1971, 1972, 1974, 1976))

