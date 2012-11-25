setClass("seasonality", representation = "list", S3methods = T)

#' Periodicity
#'
#' Determine if a univariate or multivariate time series is periodic
#' using Fisher's G Test
#' @param mts either a univariate or multivariate time series
#' @export
is.periodic = function(mts) {
    if (is.vector(mts)) {
        spec = spectrum(mts, plot = F)
        g = .fishers.g.test(spec$spec)
        if (g$p <= 0.05) T else F
    }
    else apply(mts, 2, function(j) {
        is.periodic(j)
    })
}

.fishers.g.test = function(spec) {
    g = max(spec) / sum(spec)
    q = length(spec)
    p = 1 - sum(sapply(0:q, function(j) { (-1)^j * choose(q,j) *
                                          max(0, (1-j*g))^(q-1) }))
    list(g = max(spec) / sum(spec),
         p = p)   
}

#' Deseason
#'
#' Deseasonalize a time series so that models can target the
#' unexplainable components.
#' @param mts either a vector representing 1 time series or a
#'   a data frame or matrix representing multiple time series
#' @frequency the number of observations per period - either
#'   a single numeric in the univariate case or a vector in the
#'   multivariate case
#' @export
deseason = function(mts, frequency = NA, auto = F) {
    if (is.na(frequency) & !auto) return (list(seasonal = 0, remaining = mts, freq = NA))
    if (is.na(frequency) & auto) frequency = findPeriod(mts)
    bad.freq = which((2*frequency) > nrow(mts))
    if (length(bad.freq) > 0) {
        warning (paste("in order to deseasonalize a time series, it should have
                        more than 2 periods worth of measurements.  Treating
                        the following columns as not seasonal:", bad.freq))
        frequency[bad.freq] = NA
    }
    if (is.vector(mts)) {
        if (is.na(frequency)) {
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
        periodic.cols = which(!is.na(frequency))
        if (length(periodic.cols) == 0) {
            warning("no periodic time series found")
            return (structure(list(seasonal = matrix(0, nrow(mts), ncol(mts)),
                                   remaining = mts, freq = frequency),
                              class = "seasonality"))
        }
        decomp = apply(rbind(1:ncol(mts),mts)[,periodic.cols], 2, function(j) {
            colIndex = j[1]
            j.orig = j[-1]
            j.ts = ts(j.orig, frequency=frequency[colIndex])
            y = stl(j.ts, na.action=na.omit, s.window="periodic")
            seasonal = y$time.series[,1]
            remaining = apply(y$time.series[,-1], 1, sum)
            return (list(seasonal=seasonal, remaining=remaining))
        })
        remaining = seasonal = matrix(0, nrow(mts), ncol(mts))
        decomp.seasonal = do.call('cbind', lapply(decomp, function(i) i$seasonal))
        decomp.remaining = do.call('cbind', lapply(decomp, function(i) i$remaining))
        seasonal[,periodic.cols] = decomp.seasonal
        remaining[,periodic.cols] = decomp.remaining
        return (structure(list(seasonal=seasonal, remaining=remaining,
                               freq = frequency), class = "seasonality"))
    }
}

#' Last Period of a Time Series
#'
#' Get the last complete period.  If a time series has 36 measurements
#' with a period of 24, then this function should only return the first
#' 24 measurements
#' @export
lastPeriod = function(x) {UseMethod("lastPeriod")}

lastPeriod.seasonality = function(seasonality) {
    if (length(seasonality$freq) > 1) {
        freq.indices = which(!is.na(seasonality$freq))
        if (length(freq.indices) == 0) {
            warning("object was not seasonal")
            return (NA)
        }
        n = nrow(seasonality$seasonal)
        apply(rbind(1:ncol(seasonality$seasonal), seasonality$seasonal), 2, function(j) {
            freq = seasonality$freq[colIndex]
            if (is.na(freq)) return (NA)
            colIndex = j[1]
            j.orig = j[-1]
            endIndex = floor(n / freq) * freq
            startIndex = endIndex - freq + 1
            j.orig[startIndex : endIndex]
        })
    }
    else {
        freq = seasonality$freq
        if (is.na(freq)) {
            warning("object was not seasonal")
            return (NA)
        }
        n = length(seasonality$seasonal)
        endIndex = floor(n / freq)
        startIndex = endIndex - freq + 1
        seasonality$seasonal[startIndex : endIndex]
    }
}

findPeriod = function(x) {
    if (is.vector(x)) {
        spec = spectrum(x, plot = F)
        spec.max = max(spec$spec)
        spec.max.index = which.max(spec$spec)
        p = 1 - (1 - exp(spec.max))^length(spec$spec)
        if (p <= .05) return (floor(spec$freq[spec.max] * length(spec$spec)))
        else {warning("signal is not periodic") return (NA)}
    }
    else {
        apply(x, 2, findPeriod)
    }
}

periodIndex = function(freq, i) {i %% freq}
