##############################################################################
##  This script contains the source code for a function I've modified       ## 
##  from the `cffdrs` package so `terra` is used instead of `raster`.       ##
##                                                                          ##
##  At time of manuscript revision in May 2024, `cffdrs` had been           ##
##  updated to use `terra` instead of `raster`. So the user can utilize     ##
##  `cffdrs`'s  function `fwiRaster` if preferred.                          ##
##                                                                          ##
##  This script is sourced in '00_aggregate_fwi.R'.                         ##
##############################################################################

fwiRaster_terra <- function (input, init = c(ffmc = 85, dmc = 6, dc = 15), mon = 7, 
                             out = "all", lat.adjust = TRUE, uppercase = TRUE) 
  {
    ell01 <- c(6.5, 7.5, 9, 12.8, 13.9, 13.9, 12.4, 10.9, 9.4, 
               8, 7, 6)
    ell02 <- c(7.9, 8.4, 8.9, 9.5, 9.9, 10.2, 10.1, 9.7, 9.1, 
               8.6, 8.1, 7.8)
    ell03 <- c(10.1, 9.6, 9.1, 8.5, 8.1, 7.8, 7.9, 8.3, 8.9, 
               9.4, 9.9, 10.2)
    ell04 <- c(11.5, 10.5, 9.2, 7.9, 6.8, 6.2, 6.5, 7.4, 8.7, 
               10, 11.2, 11.8)
    fl01 <- c(-1.6, -1.6, -1.6, 0.9, 3.8, 5.8, 6.4, 5, 2.4, 0.4, 
              -1.6, -1.6)
    fl02 <- c(6.4, 5, 2.4, 0.4, -1.6, -1.6, -1.6, -1.6, -1.6, 
              0.9, 3.8, 5.8)
    if (!is.na(charmatch("input", search()))) {
      detach(input)
    }
    names(input) <- tolower(names(input))
    temp <- input$temp
    prec <- input$prec
    ws <- input$ws
    rh <- input$rh
    if(units(temp) == "K") {
      temp <- temp - 273.15
      units(temp) <- "C"
    } else {
      temp <- temp
    }
    if ("lat" %in% names(input)) {
      lat <- input$lat
    }
    else {
      lat <- temp
      terra::values(lat) <- 55
    }
    if (!exists("temp") | is.null(temp)) 
      stop("temperature (temp) is missing!")
    if (!exists("prec") | is.null(prec)) 
      stop("precipitation (prec) is missing!")
    if (!is.na(sum(values(prec)[values(prec) < 0]))) 
      stop("precipiation (prec) cannot be negative!")
    if (!exists("ws") | is.null(ws)) 
      stop("wind speed (ws) is missing!")
    if (!is.na(sum(values(ws)[values(ws) < 0])))
      stop("wind speed (ws) cannot be negative!")
    if (!exists("rh") | is.null(rh)) 
      stop("relative humidity (rh) is missing!")
    if (!is.na(sum(values(rh)[values(rh) < 0])))
      stop("relative humidity (rh) cannot be negative!")
    names(init) <- tolower(names(init))
    if (is.numeric(init)) {
      if (is.null(names(init))) {
        names(init) <- c("ffmc", "dmc", "dc")
      }
      ffmc_yda <- dmc_yda <- dc_yda <- temp
      terra::values(ffmc_yda) <- init[["ffmc"]]
      terra::values(dmc_yda) <- init[["dmc"]]
      terra::values(dc_yda) <- init[["dc"]]
    }
    else {
      ffmc_yda <- init$ffmc
      dmc_yda <- init$dmc
      dc_yda <- init$dc
    }
    rh[rh >= 100] <- 99.9999
    wmo <- 147.27723 * (101 - ffmc_yda)/(59.5 + ffmc_yda)
    ra1 <- prec
    ra1[ra1 <= 0.5] <- NA
    ra1 <- ra1 - 0.5
    ra2 <- prec
    ra2[ra2 > 0.5] <- NA
    ra <- terra::cover(ra1, ra2, values = NA)
    wmo1 <- terra::mask(wmo, ra1)
    wmo2 <- terra::mask(wmo, ra2)
    wmo11 <- wmo1
    wmo11[wmo11 <= 150] <- NA
    ra11 <- ra1
    ra11[wmo1 <= 150] <- NA
    wmo11 <- wmo11 + 0.0015 * (wmo11 - 150) * (wmo11 - 150) * 
      sqrt(ra11) + 42.5 * ra11 * exp(-100/(251 - wmo11)) * 
      (1 - exp(-6.93/ra11))
    wmo12 <- wmo1
    wmo12[wmo12 > 150] <- NA
    ra12 <- ra1
    ra12[wmo1 > 150] <- NA
    wmo12 <- wmo12 + 42.5 * ra12 * exp(-100/(251 - wmo12)) * 
      (1 - exp(-6.93/ra12))
    wmo1 <- terra::cover(wmo11, wmo12, values = NA)
    wmo <- terra::cover(wmo1, wmo2, values = NA)
    wmo[wmo > 250] <- 250
    rm(ra1, ra11, ra12, ra2, wmo1, wmo2, wmo11, wmo12)
    ed <- 0.942 * (rh^0.679) + (11 * exp((rh - 100)/10)) + 0.18 * 
      (21.1 - temp) * (1 - 1/exp(rh * 0.115))
    ew <- 0.618 * (rh^0.753) + (10 * exp((rh - 100)/10)) + 0.18 * 
      (21.1 - temp) * (1 - 1/exp(rh * 0.115))
    z0 <- terra::lapp(c(wmo, ed, ew), fun = function(x, y, z) {return(x < y & x < z)})
    z0[z0 == 0] <- NA
    rh0 <- terra::mask(rh, z0)
    ws0 <- terra::mask(ws, z0)
    z <- 0.424 * (1 - (((100 - rh0)/100)^1.7)) + 0.0694 * sqrt(ws0) * 
      (1 - ((100 - rh0)/100)^8)
    z[is.na(z)] <- 0
    z <- terra::mask(z, temp)
    rm(rh0, ws0, z0)
    x <- z * 0.581 * exp(0.0365 * temp)
    z0 <- terra::lapp(c(wmo, ed, ew), fun = function(x, y, z) {return(x < y & x < z)})
    z0[z0 == 0] <- NA
    ew0 <- terra::mask(ew, z0)
    x0 <- terra::mask(x, z0)
    wmo0 <- terra::mask(wmo, z0)
    wmo1 <- ew0 - (ew0 - wmo0)/(10^x0)
    wmo2 <- wmo
    wmo2[!is.na(wmo0)] <- NA
    wm <- terra::cover(wmo1, wmo2, values = NA)
    rm(z0, ew0, x0, wmo0, wmo1, wmo2)
    z0 <- terra::lapp(c(wmo, ed), fun = function(x, y) {return(x > y)})
    z0[z0 == 0] <- NA
    rh0 <- terra::mask(rh, z0)
    ws0 <- terra::mask(ws, z0)
    z0 <- 0.424 * (1 - (rh0/100)^1.7) + 0.0694 * sqrt(ws0) * 
      (1 - (rh0/100)^8)
    z1 <- z
    z1[!is.na(z0)] <- NA
    z <- terra::cover(z0, z1, values = NA)
    rm(rh0, ws0)
    x <- z * 0.581 * exp(0.0365 * temp)
    ed0 <- mask(ed, z0)
    wmo0 <- terra::mask(wmo, z0)
    x0 <- terra::mask(x, z0)
    wm0 <- ed0 + (wmo0 - ed0)/(10^x0)
    wm1 <- terra::mask(wm, z1)
    wm <- terra::cover(wm0, wm1, values = NA)
    rm(ed0, x0, wm0, wm1, wmo0)
    ffmc <- (59.5 * (250 - wm))/(147.2 + wm)
    ffmc[ffmc > 101] <- 101
    ffmc[ffmc < 0] <- 0
    t0 <- temp
    t0[t0 < -1.1] <- -1.1
    rk <- 1.894 * (t0 + 1.1) * (100 - rh) * ell01[mon] * 1e-04
    if (lat.adjust) {
      rk[lat <= 30 & lat > 10] <- 1.894 * (t0[lat <= 30 & lat > 
                                                10] + 1.1) * (100 - rh[lat <= 30 & lat > 10]) * ell02[mon] * 
        1e-04
      rk[lat <= -10 & lat > -30] <- 1.894 * (t0[lat <= -10 & 
                                                  lat > -30] + 1.1) * (100 - rh[lat <= -10 & lat > 
                                                                                  -30]) * ell03[mon] * 1e-04
      rk[lat <= -30 & lat >= -90] <- 1.894 * (t0[lat <= -30 & 
                                                   lat >= -90] + 1.1) * (100 - rh[lat <= -30 & lat >= 
                                                                                    -90]) * ell04[mon] * 1e-04
      rk[lat <= 10 & lat > -10] <- 1.894 * (t0[lat <= 10 & 
                                                 lat > -10] + 1.1) * (100 - rh[lat <= 10 & lat > -10]) * 
        9 * 1e-04
    }
    ra <- prec
    rw <- 0.92 * ra - 1.27
    wmi <- 20 + 280/exp(0.023 * dmc_yda)
    b <- dmc_yda
    b[dmc_yda <= 33] <- 100/(0.5 + 0.3 * dmc_yda[dmc_yda <= 33])
    if (!is.null(dmc_yda[dmc_yda > 33 & dmc_yda <= 65])) {
      b[dmc_yda > 33 & dmc_yda <= 65] <- 14 - 1.3 * log(dmc_yda[dmc_yda > 
                                                                  33 & dmc_yda <= 65])
    }
    if (!is.null(dmc_yda[dmc_yda > 65])) {
      b[dmc_yda > 65] <- 6.2 * log(dmc_yda[dmc_yda > 65]) - 
        17.2
    }
    wmr <- wmi + 1000 * rw/(48.77 + b * rw)
    pr0 <- 43.43 * (5.6348 - log(wmr - 20))
    pr <- pr0
    pr[prec <= 1.5] <- dmc_yda[prec <= 1.5]
    pr[pr < 0] <- 0
    dmc <- pr + rk
    dmc[dmc < 0] <- 0
    t0[temp < (-2.8)] <- -2.8
    pe <- (0.36 * (t0 + 2.8) + fl01[mon])/2
    if (lat.adjust) {
      pe[lat <= -10] <- (0.36 * (t0[lat <= -10] + 2.8) + fl02[mon])/2
      pe[lat > -10 & lat <= 10] <- (0.36 * (t0[lat > -10 & 
                                                 lat <= 10] + 2.8) + 1.4)/2
    }
    ra <- prec
    rw <- 0.83 * ra - 1.27
    smi <- 800 * exp(-1 * dc_yda/400)
    dr0 <- dc_yda - 400 * log(1 + 3.937 * rw/smi)
    dr0[dr0 < 0] <- 0
    dr <- dr0
    dr[prec <= 2.8] <- dc_yda[prec <= 2.8]
    dc <- dr + pe
    dc[dc < 0] <- 0
    fW <- exp(0.05039 * ws)
    fm <- 147.27723 * (101 - ffmc)/(59.5 + ffmc)
    fF <- 91.9 * exp(-0.1386 * fm) * (1 + (fm^5.31)/49300000)
    isi <- 0.208 * fW * fF
    bui <- 0.8 * dc * dmc/(dmc + 0.4 * dc)
    bui[dmc == 0 & dc == 0] <- 0
    p <- (dmc - bui)/dmc
    p[dmc == 0] <- 0
    cc <- 0.92 + ((0.0114 * dmc)^1.7)
    bui0 <- dmc - cc * p
    bui0[bui0 < 0] <- 0
    bui[bui < dmc] <- bui0[bui < dmc]
    bb <- 0.1 * isi * (0.626 * (bui^0.809) + 2)
    bb[bui > 80] <- 0.1 * isi[bui > 80] * (1000/(25 + 108.64/exp(0.023 * 
                                                                   bui[bui > 80])))
    fwi <- exp(2.72 * ((0.434 * log(bb))^0.647))
    fwi[bb <= 1] <- bb[bb <= 1]
    dsr <- 0.0272 * (fwi^1.77)
    if (out == "fwi") {
      new_FWI <- c(ffmc, dmc, dc, isi, bui, fwi, 
                               dsr)
      names(new_FWI) <- c("ffmc", "dmc", "dc", "isi", "bui", 
                          "fwi", "dsr")
      if (uppercase) {
        names(new_FWI) <- toupper(names(new_FWI))
      }
    }
    else {
      if (out == "all") {
        new_FWI <- c(input, ffmc, dmc, dc, isi, 
                                 bui, fwi, dsr)
        names(new_FWI) <- c(names(input), "ffmc", "dmc", 
                            "dc", "isi", "bui", "fwi", "dsr")
        if (uppercase) {
          names(new_FWI) <- toupper(names(new_FWI))
        }
      }
    }
    return(new_FWI)
  }