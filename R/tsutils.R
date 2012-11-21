setClass("seasonality", representation = "list", S3methods = T)

#' Periodicity
#'
#' Determine if a univariate or multivariate time series is periodic
#' @param mts either a univariate or multivariate time series
#' @export
is.periodic = function(mts) {
    if (is.vector(mts)) {
        spec = spectrum(mts)
        return (F)
    }
    else apply(mts, 2, function(j) {
        is.periodic(j)
    })
}

#' Deseason
#'
#' Deseasonalize a time series so that models can target the
#' unexplainable components.
#' @param mts either a vector representing 1 time series or a
#'   a data frame or matrix representing multiple time series
#' @frequency the number of observations per period
#' @export
deseason = function(mts, frequency) {
    periodic = is.periodic(mts)
    if (is.vector(mts)) {
        if (!periodic) {
            warning("time series is not periodic")
            return (structure(list(seasonal = rep(0, length(mts)),
                                   remaining = mts, freq = frequency),
                              class = "seasonality"))
        }
        y = stl(ts(mts, frequency=frequency), na.action=na.omit,
                s.window="periodic")
        seasonal = y$time.series[,1]
        remaining = apply(y$time.series[,-1], 1, sum)
        return (structure(list(seasonal=seasonal, remaining=remaining,
                               freq = frequency), class = "seasonality"))
    }
    else {
        periodic.cols = which(periodic)
        if (length(periodic.cols) == 0) {
            warning("no periodic time series found")
            return (structure(list(seasonal = matrix(0, nrow(mts), ncol(mts)),
                                   remaining = mts, freq = frequency),
                              class = "seasonality"))
        }
        decomp = apply(mts[,periodic.cols], 2, function(j) {
            j.ts = ts(j, frequency=frequency)
            y = stl(j.ts, na.action=na.omit, s.window="periodic")
            seasonal = y$time.series[,1]
            remaining = apply(y$time.series[,-1], 1, sum)
            return (list(seasonal=seasonal, remaining=remaining))
        })
        seasonal = do.call('rbind', lapply(decomp, function(i) i$seasonal))
        remaining = do.call('rbind', lapply(decomp, function(i) i$remaining))
        return (structure(list(seasonal=seasonal, remaining=remaining,
                               freq = frequency), class = "seasonality"))
    }
}

lastPeriod = function(x) {UseMethod("lastPeriod")}

lastPeriod.seasonality = function(seasonality) {
    season = seasonality$seasonal
    if (is.vector(season)) {
        n = length(season)
        season[(n - seasonality$freq + 1) : n]
    } else {
        n = nrow(season)
        season[(n - seasonality$freq + 1) : n,]
    }
}

findPeriod = function(x) {
    return (24)
}
