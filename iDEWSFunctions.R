################################################
########## iDEWS Data Functions ################
################################################

# Author: Hugo Pedder



# Calculate mean duration in cases for which dadmin and ddisch is reported
# Subtract mean duration (could even do this using SD as well if it is large) from ddisch if dadmin is missing
impute.dadmin <- function(df) {
  # Estimate parameters for negative binomial distribution of days in hospital

  temp <- df[is.na(df$dadmin) & !is.na(df$ddisch),]
  print(paste0("Number of potential entries for which to impute date of admission from date of discharge: ", nrow(temp)))

  temp <- df[!is.na(df$dadmin) & !is.na(df$ddisch),]
  stay <- difftime(df$ddisch, df$dadmin, units = "days")
  stay <- subset(stay, !is.na(stay) & stay>=0)
  stay <- as.numeric(stay)

  # Drop values for stay > 1000
  stay <- subset(stay, stay<1000)

  print(paste0("Number of entries with which to estimate days spent in hospital: ", length(stay)))

  # Fit negative binomial dist USING MASS
  fit <- MASS::fitdistr(stay, "Negative Binomial")

  print("Parameters from negative binomial distribution")
  print(fit$estimate)
  return(fit$estimate)
}



impute.model <- function(df) {

  # temp <- df[is.na(df$dadmin) & !is.na(df$ddisch),]
  # print(paste0("Number of potential entries for which to impute date of admission from date of discharge: ", nrow(temp)))

  df$stay[!is.na(df$dadmin) & !is.na(df$ddisch)] <- 1
  df$stay[!is.na(df$stay)] <- difftime(df$ddisch[!is.na(df$stay)], df$dadmin[!is.na(df$stay)], units = "days")

  # Drop values for stay where >6 months
  df$stay[df$stay>=183 | df$stay<0] <- NA

  print(paste0("Number of entries with which to estimate days spent in hospital: ", nrow(df[!is.na(df$stay),])))
  print(paste0("Number of entries imputed for days spent in hospital: ", nrow(df[is.na(df$stay),])))

  # Using MICE (with predicted mean matching)
  impute.df <- df[, c("stay", "avtemp", "avhumid", "month", "week", "year", "season", "day", "maxtemp", "mintemp", "maxhumid",
                      "minhumid", "yrday")]
  temp <- mice(impute.df, m = 1, method = "pmm")

  for (i in seq_along(temp$imp)) {
    df[[names(temp$imp)[i]]][is.na(df[[names(temp$imp)[i]]])] <- temp$imp[[i]][["1"]]
  }


  # Using handmade imputation
  # stay.df <- subset(df, !is.na(stay))
  #
  # stay.df$stayzero <- ifelse(stay.df$stay==0, 1, 0)
  #
  # zeromod <- randomForest::randomForest(stayzero ~ avtemp + avhumid + month + week + year + season + day + maxtemp + mintemp +
  #                                        maxhumid + minhumid + yrday, data=stay.df, na.action=na.omit)
  #
  # countmod <- randomForest::randomForest(stay ~ avtemp + avhumid + tot + month + week + year + season + day + maxtemp + mintemp +
  #                                         maxhumid + minhumid + yrday, data=subset(stay.df, stayzero==0), na.action=na.omit)
  #
  # df$stayzero[is.na(df$stay)] <- round(predict(zeromod, df[is.na(df$stay),]))
  # df$stay[is.na(df$stay)] <- round(predict(countmod, df[is.na(df$stay),]))
  #
  # df$stay[df$stayzero==1] <- 0


  return(df)

}






prepare.outcome <- function(df, imputeadmin=TRUE, missingzero=TRUE, seed=890421, missingdate=c("2005-12-31", "2007-10-10"),
                            cutdate=c("2002-06-13", "2015-12-31"), agecut=5, imputemodel=FALSE, age="pyoung") {

  # df: data frame of disease/outcome data
  # imputeadmin: Whether or not to impute the date of admission using the date of discharge and drawing from a
  #distribution of days in hospital
  # missingzero: Whether to code missing counts as zero (excluding date in `missingdate` range)
  # missingdate: specifies date range for which to set counts as NA
  # cutdate: exclude dates outside of this range from the dataset
  # age: can be either as a proportion of young ("pyoung") or stratified ("strat")



  if (imputeadmin==TRUE & any(is.na(df$dadmin))) {

    print(paste0("Total number of admissions in dataset: ", nrow(df)))

    if (imputemodel==FALSE) {
      # Estimate paramters for distribution of days in hospital
      params <- impute.dadmin(df)

      # Impute missing dadmin using params
      temp <- df$ddisch[is.na(df$dadmin) & !is.na(df$ddisch)]
      x <- vector()

      if (!is.null(seed)) {
        set.seed(890421)
      }
      for (i in seq_along(temp)) {
        day <- rnbinom(1, size=params[1], mu=params[2])
        x <- append(x, day)
        temp[i] <- temp[i] - day
      }
      #hist(x)

      df$dadmin[is.na(df$dadmin) & !is.na(df$ddisch)] <- temp

    } else if (imputemodel==TRUE) {

      df <- impute.model(df)
      df$dadmin[is.na(df$dadmin) & !is.na(df$stay)] <- df$ddisch[is.na(df$dadmin) & !is.na(df$stay)] -
        df$stay[is.na(df$dadmin) & !is.na(df$stay)]

    }

    print(paste0("Hospital entries that are missing date of admission (even after imputation): ", nrow(df[is.na(df$dadmin),])))
  }

  # Drop missing dadmin
  df <- df[!is.na(df$dadmin),]

  # Add age category
  df$agecat[df$ageyrs>=agecut] <- 0
  df$agecat[df$ageyrs<agecut] <- 1

  # Count by count
  df <- arrange(df, dadmin)

  if (age=="pyoung") {
    df <- df %>% group_by(dadmin) %>% mutate(count=n())
  } else if (age=="strat") {
    df <- df[!is.na(df$ageyrs),]
    df <- df %>% group_by(dadmin, agecat) %>% mutate(count=n())
  }


  if (age=="pyoung") {
    # Calculate proportion of young (age < agecut) - ASSUMES MISISNG AGE IS MAR
    df <- df %>% group_by(dadmin, !is.na(agecat)) %>%
      mutate(n=n()) %>% mutate(pyoung=sum(agecat, na.rm=TRUE)/n)
    df <- df %>% group_by(dadmin) %>% arrange(dadmin, agecat) %>%
      mutate(pyoung=pyoung[1])

    counts.df <- df[, c("dadmin", "count", "pyoung")]
    counts.df <- unique(counts.df)

    counts.df$dadmin <- as.Date(counts.df$dadmin)
    counts.df <- arrange(counts.df, dadmin)

    if (missingzero==TRUE) {
      # Add in zeros for missing days (except missing chunk)

      alldates <- tibble("dadmin"=seq(as.Date(cutdate[1]), as.Date(cutdate[2]), by="days"))
      alldates$count <- NA
      alldates$pyoung <- NA
      alldates$count[!is.na(match(alldates$dadmin, counts.df$dadmin))] <-
        counts.df$count
      alldates$pyoung[!is.na(match(alldates$dadmin, counts.df$dadmin))] <-
        counts.df$pyoung

      alldates$count[is.na(alldates$count)] <- 0
      alldates$count[(alldates$dadmin > missingdate[1] & alldates$dadmin < missingdate[2])] <- NA # Set count=NA for missing years

      counts.df <- alldates
    }


  } else if (age=="strat") {

    counts.df <- df[, c("dadmin", "count", "agecat")]
    counts.df <- unique(counts.df)

    counts.df$dadmin <- as.Date(counts.df$dadmin)
    counts.df <- arrange(counts.df, dadmin, agecat)

    if (missingzero==TRUE) {
      # Add in zeros for missing days (except missing chunk)

      alldates <- tibble("dadmin"=rep(seq(as.Date(cutdate[1]), as.Date(cutdate[2]), by="days"), 2))
      alldates <- arrange(alldates, dadmin)
      alldates$agecat <- rep(c(0,1), nrow(alldates)/2)

      alldates <- full_join(alldates, counts.df)

      alldates$count[is.na(alldates$count)] <- 0
      alldates$count[(alldates$dadmin > missingdate[1] & alldates$dadmin < missingdate[2])] <- NA # Set count=NA for missing years

      counts.df <- alldates
    }
  }


  counts.df <- arrange(counts.df, dadmin)

  return(counts.df)
}





prepare.analysis <- function(df, time="day", var="prevtemp", fulljoin=TRUE) {
  # time: can take "day" or "week"

  # Load disease and weather data
  counts.df <- df

  # Combine into a single data frame#
  if (fulljoin==TRUE) {
    df <- full_join(counts.df, heat, by="dadmin")
    df <- df %>% rename(avtemp=Average, maxtemp=Max, mintemp=Min)

    df <- full_join(df, rain, by="dadmin")

    df <- full_join(df, humid, by="dadmin")
    df <- df %>% rename(avhumid=AVG, maxhumid=MX, minhumid=MN)
  } else {
    # df <- left_join(counts.df, heat, by="dadmin")
    # df <- df %>% rename(avtemp=Average, maxtemp=Max, mintemp=Min)

    df <- counts.df

    df$avtemp <- heat$Average[match(df$dadmin, heat$dadmin, nomatch=match(df$ddisch, heat$dadmin))]
    df$maxtemp <- heat$Max[match(df$dadmin, heat$dadmin, nomatch=match(df$ddisch, heat$dadmin))]
    df$mintemp <- heat$Min[match(df$dadmin, heat$dadmin, nomatch=match(df$ddisch, heat$dadmin))]

    # df <- left_join(df, rain, by=c("dadmin", "ddisch"))

    df$tot <- rain$tot[match(df$dadmin, rain$dadmin, nomatch=match(df$ddisch, rain$dadmin))]

    # df <- left_join(df, humid, by=c("dadmin", "ddisch"))
    # df <- df %>% rename(avhumid=AVG, maxhumid=MX, minhumid=MN)

    df$avhumid <- humid$AVG[match(df$dadmin, humid$dadmin, nomatch=match(df$ddisch, humid$dadmin))]
    df$maxhumid <- humid$MX[match(df$dadmin, humid$dadmin, nomatch=match(df$ddisch, humid$dadmin))]
    df$minhumid <- humid$MN[match(df$dadmin, humid$dadmin, nomatch=match(df$ddisch, humid$dadmin))]
  }


  df <- arrange(df, dadmin)

  # Drop columns with a "."
  if (any(grep("\\.", colnames(df)))) {
    df <- df[, -grep("\\.", colnames(df))]
  }

  # Select counts from 2007 to end 2016
  #df <- subset(df, dadmin>="2007-10-10" & dadmin < "2017-01-01")

  # Add season
  df$month <- month(df$dadmin)
  df$week <- as.numeric(strftime(df$dadmin, format = "%V"))
  df$year <- as.numeric(strftime(df$dadmin, format = "%y"))

  df$season[df$week %in% c(1:8, 48:53)] <- "DJF"
  df$season[df$week %in% c(9:21)] <- "MAM"
  df$season[df$week %in% c(22:34)] <- "JJA"
  df$season[df$week %in% c(35:47)] <- "SON"
  df$season <- factor(df$season)


  # Add day of the week (assumes first row is a Thursday) starting on Sunday
  df$day <- factor(weekdays(df$dadmin), levels=c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday"))
  df$weekend[df$day %in% c("Saturday", "Sunday")] <- 1
  df$weekend[!df$day %in% c("Saturday", "Sunday")] <- 0

  # Add day of the year
  df$yrday <- yday(df$dadmin)

  # Drop week 53 (this will lead to strange values since week will not have full counts)
  df <- df[df$week!=53 | is.na(df$week),]


  if (time=="week") {
    # Cluster by week (probably only for malaria)
    week.df <- df
    week.df <- week.df %>% group_by(week, year) %>% mutate(count=sum(count, na.rm = FALSE))

    week.df <- week.df %>% group_by(week, year) %>% mutate(tot=sum(tot, na.rm = FALSE))

    week.df <- week.df %>% group_by(week, year) %>% mutate(avtemp=mean(avtemp, na.rm = TRUE))
    week.df <- week.df %>% group_by(week, year) %>% mutate(mintemp=min(mintemp, na.rm = TRUE))
    week.df <- week.df %>% group_by(week, year) %>% mutate(maxtemp=max(maxtemp, na.rm = TRUE))

    week.df <- week.df %>% group_by(week, year) %>% mutate(avhumid=mean(avhumid, na.rm = TRUE))
    week.df <- week.df %>% group_by(week, year) %>% mutate(minhumid=min(minhumid, na.rm = TRUE))
    week.df <- week.df %>% group_by(week, year) %>% mutate(maxhumid=max(maxhumid, na.rm = TRUE))

    week.df <- week.df %>% group_by(week, year) %>% mutate(pyoung=sum(count * pyoung, na.rm = TRUE) / sum(count))
    week.df <- unique(week.df[, names(week.df) %in% c("count", "tot", "week", "season", "year", "pyoung",
                                                      "agecat", "avtemp", "avhumid", "mintemp", "minhumid",
                                                      "maxtemp", "maxhumid")])
    week.df$id <- c(1:nrow(week.df))
    df <- week.df
  }

  # Add previous day/week's and temperature range temp/humidity
  df$diftemp <- df$maxtemp - df$mintemp
  df$prevtemp <- NA
  df$prevtemp[2:nrow(df)] <- df$avtemp[2:nrow(df)] - df$avtemp[1:(nrow(df)-1)]

  df$difhumid <- df$maxhumid - df$minhumid
  df$prevhumid <- NA
  df$prevhumid[2:nrow(df)] <- df$avhumid[2:nrow(df)] - df$avhumid[1:(nrow(df)-1)]

  # Add indicator for daily temperature range (+ve or -ve)
  # df$diftemp[!is.na(df$`Min Tm`) & !is.na(df$`Max Tm`) & df$`Min Tm` > df$`Max Tm`] <-
  #   df$diftemp[!is.na(df$`Min Tm`) & !is.na(df$`Max Tm`) & df$`Min Tm` > df$`Max Tm`] * -1

  # Add lags
  if (!is.null(var)) {
    addlag <- function(df, var, lag=1) {
      lim1 <- nrow(df) - (lag-1)
      lim2 <- nrow(df)
      df[[paste0("lag", lag)]] <- c(rep(NA, lag), df[[var]][-c(lim1:lim2)])

      return(df)
    }

    for (i in 1:10) {
      #df <- addlag(df, lag=i, var="tot")
      df <- addlag(df, lag=i, var=var)
    }
  }

  return(df)
}



prepare.weather <- function(impute=TRUE) {
  ########### LOAD WEATHER DATA #############

  ####### Load Thohoyandou data ######
  humid <- read_excel("~/WASH/Climate data 20-09-2017.xls",
                      skip = 24, sheet = 2)
  humid$dadmin <- as.Date(humid$Date)
  humid$Date <- NULL

  rain <- read_excel("~/WASH/Climate data 20-09-2017.xls",
                     skip = 14, sheet = 3)
  rain$tot <- as.numeric(rain$tot)
  rain$dadmin <- as.Date(rain$DATE)
  rain$DATE <- NULL

  heat <- read_excel("~/WASH/Climate data 20-09-2017.xls",
                     skip = 22, sheet = 4)
  heat$dadmin <- as.Date(heat$Date)
  heat$Date <- NULL


  windy <- read_excel("~/WASH/Climate data 20-09-2017.xls",
                      skip = 46, sheet = 5)
  windy$dadmin <- as.Date(windy$Date)
  windy$Date <- NULL


  windsp <- read_excel("~/WASH/Climate data 20-09-2017.xls",
                       skip = 22, sheet = 6)
  windsp$dadmin <- as.Date(windsp$Date)
  windsp$Date <- NULL



  ####### Load Giyani data ######
  giyani <- read_excel("~/WASH/Climate data 20-09-2017.xls",
                       sheet = 7, col_names=FALSE)

  start <- grep("HOURLY DATA", giyani$...1)
  end <- grep("^(avg)|(Max)|(tot)$", giyani$...1)

  # Drop wind direction start points
  start <- start[start<=7605 | start>=10193]

  weather <- list()
  weather.lab <- vector()
  for (i in seq_along(start)) {

    month <- giyani[(start[i]+4):end[i]-1,]
    names(month) <- as.vector(giyani[(start[i]+2),])

    temp <- giyani$...1[start[i]]
    temp <- gsub("(HOURLY DATA( )?\\:)(.+$)", "\\3", temp)

    type <- gsub("(^.+)(-)(.+$)", "\\1", temp)
    type <- trimws(type)
    weather.lab <- append(weather.lab, type)

    date <- gsub("(^.+)(-)(.+$)", "\\3", temp)
    date <- trimws(date)

    dat <- c("type"=type, "date"=date, "month"=list(month))

    weather[[length(weather)+1]] <- dat
  }
  weather.lab <- unique(weather.lab)


  # Extract different weather stats
  g.rain <- data.frame("DD"=NA, "tot"=NA, "date"=NA)
  g.wind <- data.frame("DD"=NA, "avg"=NA, "mx"=NA, "tm"=NA, "date"=NA)
  g.humid <- data.frame("DD"=NA, "avg"=NA, "mx"=NA, "tm"=NA, "mn"=NA, "tm"=NA, "date"=NA)
  g.temp <- data.frame("DD"=NA, "avg"=NA, "mx"=NA, "tm"=NA, "mn"=NA, "tm"=NA, "date"=NA)
  for (i in seq_along(weather)) {
    month <- weather[[i]]$month

    if (weather[[i]]$type=="Humidity(%)") {
      month <- month[,c(1:6)]
      month$date <- weather[[i]]$date
      g.humid <- rbind(g.humid, month)

    } else if (weather[[i]]$type=="Rain(mm)") {
      month <- month[,c(1:2)]
      month$date <- weather[[i]]$date
      g.rain <- rbind(g.rain, month)

    } else if (weather[[i]]$type=="Temperature (C)") {
      month <- month[,c(1:6)]
      month$date <- weather[[i]]$date
      g.temp <- rbind(g.temp, month)

    } else if (weather[[i]]$type=="Average Wind Speed (m/s)") {
      month <- month[,c(1:4)]
      month$date <- weather[[i]]$date
      g.wind <- rbind(g.wind, month)

    } else {
      stop("None of the weather types")
    }

  }
  g.rain <- g.rain[-1,]
  g.wind <- g.wind[-1,]
  g.temp <- g.temp[-1,]
  g.humid <- g.humid[-1,]

  g.rain$tot <- gsub("=", "", g.rain$tot)
  g.rain[,-ncol(g.rain)] <- apply(g.rain[,-ncol(g.rain)], MARGIN = 2, FUN=function(x) {
    as.numeric(x)
  })
  g.rain$dadmin <- as.Date(paste(g.rain$DD, g.rain$date), format="%d %B %Y")

  g.wind[,-ncol(g.wind)] <- apply(g.wind[,-ncol(g.wind)], MARGIN = 2, FUN=function(x) {
    as.numeric(x)
  })
  g.wind$dadmin <- as.Date(paste(g.wind$DD, g.wind$date), format="%d %B %Y")

  g.temp[,-ncol(g.temp)] <- apply(g.temp[,-ncol(g.temp)], MARGIN = 2, FUN=function(x) {
    as.numeric(x)
  })
  g.temp$dadmin <- as.Date(paste(g.temp$DD, g.temp$date), format="%d %B %Y")

  g.humid[,-ncol(g.humid)] <- apply(g.humid[,-ncol(g.humid)], MARGIN = 2, FUN=function(x) {
    as.numeric(x)
  })
  g.humid$dadmin <- as.Date(paste(g.humid$DD, g.humid$date), format="%d %B %Y")


  ########## Merge Giyani and Thoyanadou #############

  g.rain$g.tot <- g.rain$tot
  rain <- full_join(rain, g.rain, by="dadmin")
  rain$tot <- rain$tot.x

  names(g.wind)[2:4] <- paste0("g.", names(g.wind)[2:4])
  windsp <- full_join(windsp, g.wind, by="dadmin")

  names(g.temp)[2:6] <- paste0("g.", names(g.temp)[2:6])
  heat <- full_join(heat, g.temp, by="dadmin")

  names(g.humid)[2:6] <- paste0("g.", names(g.humid)[2:6])
  humid <- full_join(humid, g.humid, by="dadmin")

  if (impute==TRUE) {
    library(mice)

    weather.df <- cbind(rain, heat, windsp, humid)
    weather.df <- data.frame(rain$dadmin, rain$tot, heat$Average, heat$Max, heat$Min, humid$AVG, humid$MX, humid$MN,
                             windsp$avg, windsp$mx)

    imp <- mice(weather.df, maxit = 2, m = 10, seed = 1)

    mistot <- apply(imp$imp$rain.tot, MARGIN=1, FUN=function(x) {mean(x)})

    misheatavg <- apply(imp$imp$heat.Average, MARGIN=1, FUN=function(x) {mean(x)})
    misheatmax <- apply(imp$imp$heat.Max, MARGIN=1, FUN=function(x) {mean(x)})
    misheatmin <- apply(imp$imp$heat.Min, MARGIN=1, FUN=function(x) {mean(x)})

    mishumidavg <- apply(imp$imp$humid.AVG, MARGIN=1, FUN=function(x) {mean(x)})
    mishumidmax <- apply(imp$imp$humid.MX, MARGIN=1, FUN=function(x) {mean(x)})
    mishumidmin <- apply(imp$imp$humid.MN, MARGIN=1, FUN=function(x) {mean(x)})

    rain$tot[is.na(rain$tot)] <- mistot

    heat$Average[is.na(heat$Average)] <- misheatavg
    heat$Max[is.na(heat$Max)] <- misheatmax
    heat$Min[is.na(heat$Min)] <- misheatmin

    humid$AVG[is.na(humid$AVG)] <- mishumidavg
    humid$MX[is.na(humid$MX)] <- mishumidmax
    humid$MN[is.na(humid$MN)] <- mishumidmin

  }

  return(list(rain=rain, windsp=windsp, heat=heat, humid=humid))
}








# Pools multiple zero-inflated negbin models (e.g. from multiple imputation) using Rubin's Rules
pool.rubin <- function(imps=list(), param="count") {

  ncoef <- length(imps[[1]]$coefficients[[param]])
  coefmat <- matrix(ncol=length(imps), nrow=ncoef)
  covarray <- array(dim=c(ncoef, ncoef, length(imps)))
  for (i in seq_along(imps)) {
    coefmat[,i] <- imps[[i]]$coefficients[[param]]

    matnam <- rownames(imps[[i]]$vcov)
    index <- which(grepl(paste0("^", param), matnam))
    covarray[,,i] <- imps[[i]]$vcov[index, index]
  }
  coef <- apply(coefmat, MARGIN=1, FUN=function(x) {mean(x)})
  names(coef) <- names(imps[[1]]$coefficients[[param]])

  within <- apply(covarray, MARGIN=c(1,2), FUN=function(x) {mean(x)})
  between <- apply(coefmat, MARGIN=1, FUN=function(x) {var(x)})

  vcov <- within + between/length(imps)
  #vcov <- within + between/length(imps) + between
  #vcov <- within + (1+1/length(imps)) * (between)
  colnames(vcov) <- names(imps[[1]]$coefficients[[param]])
  rownames(vcov) <- names(imps[[1]]$coefficients[[param]])

  return(list(coef=coef, vcov=vcov))
}






# Run multiple imputation for negative binomial imputation of days in hospital
# IF MODEL CHANGES THEN THESE NEED TO BE CHANGED HERE TOO!
idews.mi <- function(pneum.df, reps=5, nb.size=0.4700373, nb.mu=16.6820675, seed=890421, imputemodel=FALSE,
                     sensitivity=NULL, tempknot=4, templag=3, humidknot=3, humidlag=3, diftempknot=4) {

  set.seed(seed)

  imps <- list()
  imps.df <- list()
  for (i in 1:reps) {

    if (imputemodel==FALSE) {
      df <- prepare.outcome(pneum.df, seed=NULL)

      df <- prepare.analysis(df, time="day", var="count")
      df <- df[df$dadmin>="2007-10-10" & df$dadmin<"2016-01-01",]
    } else {
      df <- prepare.analysis(pneum.df, time="day", fulljoin=FALSE)
      df <- prepare.outcome(df, imputemodel = imputemodel)
      df <- prepare.analysis(df, time="day", var="count")
      df <- df[df$dadmin>="2007-10-10" & df$dadmin<"2016-01-01",]
    }

    if (!is.null(sensitivity)) {
      if (sensitivity=="mintemp") {
        #### Basis matrices ####
        mintempmat <- crossbasis(df$mintemp, lag=21, # Lag of 21 days
                                argvar=list(fun="ns", df=tempknot), # Natural spline with 4 knots for mean temperature
                                arglag=list(fun="ns", df=templag)) # Natural spline with 3 knots for mean temperature lag

        humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                               argvar=list(fun="ns", df=humidknot), # Natural spline with 3 knots for mean temperature
                               arglag=list(fun="ns", df=humidlag)) # Natural spline with 3 knots for mean temperature lag

        diftempmat <- crossbasis(df$diftemp, # No lag
                                 argvar=list(fun="ns", df=diftempknot)) # Natural spline with 3 knots for daily temperature range

        lagmat <- crossbasis(df$lag1, lag=3, # Lag of 21 days
                             argvar=list(fun="lin"))#, # Natural spline with 3 knots for mean temperature

        # BINOMIAL COMPONENT ONLY AFFECTED BY DAY OF THE WEEK AND LAG
        model <- zeroinfl(count ~ mintempmat + humidmat + diftempmat +
                            sin(2*pi*(yrday/(365))) + cos(2*pi*(yrday/(365))) +
                            splines::ns(dadmin, 5) + day + lagmat | day + lagmat,
                          dist="negbin", data=df
        )
      } else if (sensitivity=="maxtemp") {
        #### Basis matrices ####
        maxtempmat <- crossbasis(df$maxtemp, lag=21, # Lag of 21 days
                                 argvar=list(fun="ns", df=tempknot), # Natural spline with 4 knots for mean temperature
                                 arglag=list(fun="ns", df=templag)) # Natural spline with 3 knots for mean temperature lag

        humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                               argvar=list(fun="ns", df=humidknot), # Natural spline with 3 knots for mean temperature
                               arglag=list(fun="ns", df=humidlag)) # Natural spline with 3 knots for mean temperature lag

        diftempmat <- crossbasis(df$diftemp, # No lag
                                 argvar=list(fun="ns", df=diftempknot)) # Natural spline with 3 knots for daily temperature range

        lagmat <- crossbasis(df$lag1, lag=3, # Lag of 21 days
                             argvar=list(fun="lin"))#, # Natural spline with 3 knots for mean temperature

        # BINOMIAL COMPONENT ONLY AFFECTED BY DAY OF THE WEEK AND LAG
        model <- zeroinfl(count ~ maxtempmat + humidmat + diftempmat +
                            sin(2*pi*(yrday/(365))) + cos(2*pi*(yrday/(365))) +
                            splines::ns(dadmin, 5) + day + lagmat | day + lagmat,
                          dist="negbin", data=df
        )
      }
    } else {
      #### Basis matrices ####
      avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                              argvar=list(fun="ns", df=tempknot), # Natural spline with 4 knots for mean temperature
                              arglag=list(fun="ns", df=templag)) # Natural spline with 3 knots for mean temperature lag

      humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                             argvar=list(fun="ns", df=humidknot), # Natural spline with 3 knots for mean temperature
                             arglag=list(fun="ns", df=humidlag)) # Natural spline with 3 knots for mean temperature lag

      diftempmat <- crossbasis(df$diftemp, # No lag
                               argvar=list(fun="ns", df=diftempknot)) # Natural spline with 3 knots for daily temperature range

      lagmat <- crossbasis(df$lag1, lag=3, # Lag of 21 days
                           argvar=list(fun="lin"))#, # Natural spline with 3 knots for mean temperature

      # BINOMIAL COMPONENT ONLY AFFECTED BY DAY OF THE WEEK AND LAG
      model <- zeroinfl(count ~ avtempmat + humidmat + diftempmat +
                          sin(2*pi*(yrday/(365))) + cos(2*pi*(yrday/(365))) +
                          splines::ns(dadmin, 5) + day + lagmat | day + lagmat,
                        dist="negbin", data=df
      )
    }


    imps[[length(imps)+1]] <- model
    imps.df[[length(imps.df)+1]] <- df
  }

  return(list(imps=imps, imps.df=imps.df))
}






plotpred <- function(pred, dens=NULL, adddens=TRUE, xlab="Temp difference", xlim=NULL, logy=FALSE) {

  cols <- brewer.pal(5, "Set1")

  pred.df <- data.frame(x=as.numeric(names(pred$allRRfit)), rr=pred$allRRfit, l95=pred$allRRlow, u95=pred$allRRhigh)

  g <- ggplot(pred.df, aes(x=x, y=rr, ymin=l95, ymax=u95)) +
    geom_ribbon(aes(ymin=l95, ymax=u95), fill=cols[2], alpha=0.5) +
    geom_line(aes(x=x, y=rr), linetype="solid", size=0.75) +
    geom_hline(aes(yintercept=1)) +
    geom_vline(aes(xintercept=pred$cen), linetype="dashed") +
    xlab(xlab) + ylab("RR") +
    theme_classic()

  if (logy==TRUE) {
    g <- g + scale_y_continuous(trans='log2')
  }

  if (!is.null(xlim)) {
    g <- g + xlim(xlim)
  }

  if (adddens==TRUE & !is.null(dens)) {
    g <- g + theme(axis.line.x=element_blank(),axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.x=element_blank())

    g2 <- ggplot(df, aes(df[[dens]])) +
      geom_density(fill=cols[3], size=0.75, alpha=0.5) +
      xlab(xlab) + ylab("Variable density") +
      geom_vline(aes(xintercept=pred$cen), linetype="dashed") +
      theme_classic()

    if (!is.null(xlim)) {
      g2 <- g2 + xlim(xlim)
    }

    grid.newpage()
    grid.draw(rbind(ggplotGrob(g), ggplotGrob(g2), size = "last"))
  } else {
    return(g)
  }

}


plotpred2 <- function(pred, pred.s1=NULL, pred.s2=NULL, labs=c("Final model", "Sensitivity 1", "Sensitivity 2"),
                      dens=NULL, adddens=TRUE, xlab="Temp difference", xlim=NULL, logy=FALSE) {

  cols <- brewer.pal(5, "Set1")

  pred.df <- data.frame(x=as.numeric(names(pred$allRRfit)), rr=pred$allRRfit, l95=pred$allRRlow, u95=pred$allRRhigh, model=1)

  if (!is.null(pred.s1)) {
    pred.df <- rbind(pred.df,
                     data.frame(x=as.numeric(names(pred.s1$allRRfit)), rr=pred.s1$allRRfit,
                                l95=pred.s1$allRRlow, u95=pred.s1$allRRhigh, model=2)
                     )
  }
  if (!is.null(pred.s2)) {
    pred.df <- rbind(pred.df,
                     data.frame(x=as.numeric(names(pred.s2$allRRfit)), rr=pred.s2$allRRfit,
                                l95=pred.s2$allRRlow, u95=pred.s2$allRRhigh, model=3)
    )
  }

  pred.df$Model <- factor(pred.df$model, labels=labs[1:length(unique(pred.df$model))])

  g <- ggplot(pred.df, aes(x=x, y=rr, ymin=l95, ymax=u95)) +
    geom_ribbon(aes(ymin=l95, ymax=u95, fill=Model), alpha=0.5) +
    geom_line(aes(x=x, y=rr, linetype=Model), size=0.75) +
    geom_hline(aes(yintercept=1)) +
    geom_vline(aes(xintercept=pred$cen), linetype="dashed") +
    xlab(xlab) + ylab("RR") +
    theme_classic()

  if (logy==TRUE) {
    g <- g + scale_y_continuous(trans='log2')
  }

  if (!is.null(xlim)) {
    g <- g + xlim(xlim)
  }

  if (is.null(pred.s1)) {
    g <- g + theme(legend.position="none")
  }

  if (adddens==TRUE & !is.null(dens)) {
    g <- g + theme(axis.line.x=element_blank(),axis.text.x=element_blank(),
                   axis.ticks.x=element_blank(),
                   axis.title.x=element_blank())

    g2 <- ggplot(df, aes(df[[dens]])) +
      geom_density(fill=cols[3], size=0.75, alpha=0.5) +
      xlab(xlab) + ylab("Variable density") +
      geom_vline(aes(xintercept=pred$cen), linetype="dashed") +
      theme_classic()

    if (!is.null(xlim)) {
      g2 <- g2 + xlim(xlim)
    }
    if (is.null(pred.s1)) {
      g2 <- g2 + theme(legend.position="none")
    }

    grid.newpage()
    grid.draw(rbind(ggplotGrob(g), ggplotGrob(g2), size = "last"))
  } else {
    return(g)
  }

}



plotlag <- function(pred, vars=c(12,28), facetlab=c("Lagged association at", "°C mean daily temperature")) {

  cols <- brewer.pal(3, "Set1")

  plot.df <- data.frame(x=NA, rr=NA, l95=NA, u95=NA, var=NA)
  for (i in seq_along(vars)) {

    temp.df <- data.frame(x=as.vector(as.numeric(gsub("lag", "", colnames(pred$matRRfit)))),
                          rr=pred$matRRfit[which(rownames(pred$matRRfit)==vars[i]),],
                          l95=pred$matRRlow[which(rownames(pred$matRRfit)==vars[i]),],
                          u95=pred$matRRhigh[which(rownames(pred$matRRfit)==vars[i]),],
                          var=vars[i]
                          )

    plot.df <- rbind(plot.df, temp.df)
  }
  plot.df <- plot.df[-1,]

  plot.df$var <- factor(plot.df$var, labels=paste(facetlab[1], unique(plot.df$var), facetlab[2]))
  #plot.df$var <- factor(plot.df$var)

  g <- ggplot(plot.df, aes(x=x, y=rr, ymin=l95, ymax=u95)) +
    geom_ribbon(aes(ymin=l95, ymax=u95), fill=cols[2], alpha=0.5) +
    geom_line(aes(x=x, y=rr), size=0.75) +
    geom_hline(aes(yintercept=1)) +
    xlab("Lag (days)") + ylab("RR") +
    theme_classic()

  g <- g + facet_wrap(~var, nrow = 2)
  return(g)
}




plotslice <- function(pred, lags=c(0,7,14,21), xlab="Mean daily temperature (°C)", logy=FALSE) {
  cols <- brewer.pal(3, "Set1")

  plot.df <- data.frame(x=NA, rr=NA, l95=NA, u95=NA, lag=NA)
  for (i in seq_along(lags)) {

    temp.df <- data.frame(x=as.vector(as.numeric(rownames(pred$matRRfit))),
                          rr=pred$matRRfit[,lags[i]+1],
                          l95=pred$matRRlow[,lags[i]+1],
                          u95=pred$matRRhigh[,lags[i]+1],
                          lag=lags[i]
    )

    plot.df <- rbind(plot.df, temp.df)
  }
  plot.df <- plot.df[-1,]

  plot.df$lag <- factor(plot.df$lag, labels=paste("Association at", unique(plot.df$lag), "days lag"))

  g <- ggplot(plot.df, aes(x=x, y=rr, ymin=l95, ymax=u95)) +
    geom_ribbon(aes(ymin=l95, ymax=u95), fill=cols[2], alpha=0.5) +
    geom_line(aes(x=x, y=rr), size=0.75) +
    geom_hline(aes(yintercept=1)) +
    geom_vline(aes(xintercept=pred$cen), linetype="dashed") +
    xlab(xlab) + ylab("RR") +
    theme_classic()

  g <- g + facet_wrap(~lag)

  if (logy==TRUE) {
    g <- g + scale_y_continuous(trans='log2')
  }

  return(g)
}






# Attributable risk (From Gasparrini 2014)

attrdl <- function(x,basis,cases,model=NULL,coef=NULL,vcov=NULL,type="af",
                   dir="back",tot=TRUE,cen,range=NULL,sim=FALSE,nsim=5000) {
  ################################################################################
  #
  # CHECK VERSION OF THE DLNM PACKAGE
  if(packageVersion("dlnm")<"2.2.0")
    stop("update dlnm package to version >= 2.2.0")
  #
  # EXTRACT NAME AND CHECK type AND dir
  name <- deparse(substitute(basis))
  type <- match.arg(type,c("an","af"))
  dir <- match.arg(dir,c("back","forw"))
  #
  # DEFINE CENTERING
  if(missing(cen) && is.null(cen <- attr(basis,"argvar")$cen))
    stop("'cen' must be provided")
  if(!is.numeric(cen) && length(cen)>1L) stop("'cen' must be a numeric scalar")
  attributes(basis)$argvar$cen <- NULL
  #
  # SELECT RANGE (FORCE TO CENTERING VALUE OTHERWISE, MEANING NULL RISK)
  if(!is.null(range)) x[x<range[1]|x>range[2]] <- cen
  #
  # COMPUTE THE MATRIX OF
  #   - LAGGED EXPOSURES IF dir="back"
  #   - CONSTANT EXPOSURES ALONG LAGS IF dir="forw"
  lag <- attr(basis,"lag")
  if(NCOL(x)==1L) {
    at <- if(dir=="back") tsModel:::Lag(x,seq(lag[1],lag[2])) else
      matrix(rep(x,diff(lag)+1),length(x))
  } else {
    if(dir=="forw") stop("'x' must be a vector when dir='forw'")
    if(ncol(at <- x)!=diff(lag)+1)
      stop("dimension of 'x' not compatible with 'basis'")
  }
  #
  # NUMBER USED FOR THE CONTRIBUTION AT EACH TIME IN FORWARD TYPE
  #   - IF cases PROVIDED AS A MATRIX, TAKE THE ROW AVERAGE
  #   - IF PROVIDED AS A TIME SERIES, COMPUTE THE FORWARD MOVING AVERAGE
  #   - THIS EXCLUDES MISSING ACCORDINGLY
  # ALSO COMPUTE THE DENOMINATOR TO BE USED BELOW
  if(NROW(cases)!=NROW(at)) stop("'x' and 'cases' not consistent")
  if(NCOL(cases)>1L) {
    if(dir=="back") stop("'cases' must be a vector if dir='back'")
    if(ncol(cases)!=diff(lag)+1) stop("dimension of 'cases' not compatible")
    den <- sum(rowMeans(cases,na.rm=TRUE),na.rm=TRUE)
    cases <- rowMeans(cases)
  } else {
    den <- sum(cases,na.rm=TRUE)
    if(dir=="forw")
      cases <- rowMeans(as.matrix(tsModel:::Lag(cases,-seq(lag[1],lag[2]))))
  }
  #
  ################################################################################
  #
  # EXTRACT COEF AND VCOV IF MODEL IS PROVIDED
  if(!is.null(model)) {
    cond <- paste0(name,"[[:print:]]*v[0-9]{1,2}\\.l[0-9]{1,2}")
    if(ncol(basis)==1L) cond <- name
    model.class <- class(model)
    coef <- dlnm:::getcoef(model,model.class)
    ind <- grep(cond,names(coef))
    coef <- coef[ind]
    vcov <- dlnm:::getvcov(model,model.class)[ind,ind,drop=FALSE]
    model.link <- dlnm:::getlink(model,model.class)
    if(model.link!="log") stop("'model' must have a log link function")
  }
  #
  # IF REDUCED ESTIMATES ARE PROVIDED
  typebasis <- ifelse(length(coef)!=ncol(basis),"one","cb")
  #
  ################################################################################
  #
  # PREPARE THE ARGUMENTS FOR TH BASIS TRANSFORMATION
  predvar <- if(typebasis=="one") x else seq(NROW(at))
  predlag <- if(typebasis=="one") 0 else dlnm:::seqlag(lag)
  #
  # CREATE THE MATRIX OF TRANSFORMED CENTRED VARIABLES (DEPENDENT ON typebasis)
  if(typebasis=="cb") {
    Xpred <- dlnm:::mkXpred(typebasis,basis,at,predvar,predlag,cen)
    Xpredall <- 0
    for (i in seq(length(predlag))) {
      ind <- seq(length(predvar))+length(predvar)*(i-1)
      Xpredall <- Xpredall + Xpred[ind,,drop=FALSE]
    }
  } else {
    basis <- do.call(onebasis,c(list(x=x),attr(basis,"argvar")))
    Xpredall <- dlnm:::mkXpred(typebasis,basis,x,predvar,predlag,cen)
  }
  #
  # CHECK DIMENSIONS
  if(length(coef)!=ncol(Xpredall))
    stop("arguments 'basis' do not match 'model' or 'coef'-'vcov'")
  if(any(dim(vcov)!=c(length(coef),length(coef))))
    stop("arguments 'coef' and 'vcov' do no match")
  if(typebasis=="one" && dir=="back")
    stop("only dir='forw' allowed for reduced estimates")
  #
  ################################################################################
  #
  # COMPUTE AF AND AN
  af <- 1-exp(-drop(as.matrix(Xpredall%*%coef)))
  an <- af*cases
  #
  # TOTAL
  #   - SELECT NON-MISSING OBS CONTRIBUTING TO COMPUTATION
  #   - DERIVE TOTAL AF
  #   - COMPUTE TOTAL AN WITH ADJUSTED DENOMINATOR (OBSERVED TOTAL NUMBER)
  if(tot) {
    isna <- is.na(an)
    af <- sum(an[!isna])/sum(cases[!isna])
    an <- af*den
  }
  #
  ################################################################################
  #
  # EMPIRICAL CONFIDENCE INTERVALS
  if(!tot && sim) {
    sim <- FALSE
    warning("simulation samples only returned for tot=T")
  }
  if(sim) {
    # SAMPLE COEF
    k <- length(coef)
    eigen <- eigen(vcov)
    X <- matrix(rnorm(length(coef)*nsim),nsim)
    coefsim <- coef + eigen$vectors %*% diag(sqrt(eigen$values),k) %*% t(X)
    # RUN THE LOOP
    # pre_afsim <- (1 - exp(- Xpredall %*% coefsim)) * cases # a matrix
    # afsim <- colSums(pre_afsim,na.rm=TRUE) / sum(cases[!isna],na.rm=TRUE)
    afsim <- apply(coefsim,2, function(coefi) {
      ani <- (1-exp(-drop(Xpredall%*%coefi)))*cases
      sum(ani[!is.na(ani)])/sum(cases[!is.na(ani)])
    })
    ansim <- afsim*den
  }
  #
  ################################################################################
  #
  res <- if(sim) {
    if(type=="an") ansim else afsim
  } else {
    if(type=="an") an else af
  }
  #
  return(res)
}
