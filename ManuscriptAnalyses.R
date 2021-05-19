############################################################
####### Run analyses for iDEWS Pneumonia Manuscript ########
############################################################

# Author: HUgo Pedder

library(dlnm)
library(tidyverse)
library(lubridate)
library(knitr)
library(readxl)
library(lmtest)
library(MASS)
library(pscl)
library(RColorBrewer)
library(ggpubr)
library(grid)
library(mice)

####### Functions ########
source("~/WASH/Pneumonia Manuscript/PneumoniaAnalyses/iDEWSFunctions.R")

# Load data
load(file="~/WASH/Pneumonia Manuscript/PneumoniaAnalyses/iDEWSdat.RData")

# Prepare weather
weather <- prepare.weather(impute = FALSE)
for (i in seq_along(weather)) {
  assign(names(weather)[i], weather[[i]])
}



# Prepare data
pneum.df <- hosp.df[hosp.df$pneum %in% c("Yes", "Y") & !is.na(hosp.df$pneum),]
pneum.df <- subset(pneum.df, !(grepl("asp", radmin1) |  # Drop if "aspiration" is specified as reason for admission
                     grepl("asp", radmin2) |
                     grepl("asp", radmin3)))

pneum.df <- subset(pneum.df, !(grepl("pneumothorax", radmin1) |  # Drop if "pneumothorax" is specified as reason for admission
                                 grepl("pneumothorax", radmin2) |
                                 grepl("pneumothorax", radmin3)))

pneum.df <- subset(pneum.df, (!(grepl("chem", radmin1) |  # Drop if "checm" is specified as reason for admission (for chemical pneumonia)
                            grepl("chem", radmin2) |
                            grepl("chem", radmin3))))

# Numbers of missing info
temp <- subset(pneum.df, dadmin>="2007-10-10" & dadmin<"2016-01-01")
table(is.na(temp$ageyrs))[2] / nrow(temp)
table(is.na(temp$sex))[2] / nrow(temp)

df <- prepare.outcome(pneum.df)

df <- prepare.analysis(df, time="day", var="count")
df <- df[df$dadmin>="2007-10-10" & df$dadmin<"2016-01-01",]






# Alternative data prep (using model imputation)
pneum.df <- hosp.df[hosp.df$pneum %in% c("Yes", "Y") & !is.na(hosp.df$pneum),]
pneum.df <- subset(pneum.df, !(grepl("asp", radmin1) |  # Drop if "aspiration" is specified as reason for admission
                                 grepl("asp", radmin2) |
                                 grepl("asp", radmin3)))

pneum.df <- subset(pneum.df, !(grepl("pneumothorax", radmin1) |  # Drop if "pneumothorax" is specified as reason for admission
                                 grepl("pneumothorax", radmin2) |
                                 grepl("pneumothorax", radmin3)))

pneum.df <- subset(pneum.df, (!(grepl("chem", radmin1) |  # Drop if "checm" is specified as reason for admission (for chemical pneumonia)
                                  grepl("chem", radmin2) |
                                  grepl("chem", radmin3))))

pneum.df$dadmin <- as.Date(pneum.df$dadmin)
df <- prepare.analysis(pneum.df, time="day", fulljoin=FALSE, var=NULL)
df <- prepare.outcome(df, imputemodel = TRUE)
df <- prepare.analysis(df, time="day", var="count")
df <- df[df$dadmin>="2007-10-10" & df$dadmin<"2016-01-01",]


# Figure 2
pdf("~/WASH/Pneumonia Manuscript/Graphs/Figure2.pdf", width=15)
par(mfrow=c(2,2))
plot(df$dadmin, df$count, xlab="Date of admission", ylab="Daily admissions", main="Daily hospital admissions for pneumonia over time",
     type="p")
plot(df$dadmin, df$avtemp, xlab="Date of admission", ylab="Mean daily temperature (C)", main="Mean daily temperature over time",
     type="p", ylim=c(0,35))
plot(df$dadmin, df$avhumid, xlab="Date of admission", ylab="Relative humidity (%)", main="Relative humidity over time",
     type="p", ylim=c(0,100))
plot(df$dadmin, df$tot, xlab="Date of admission", ylab="Total daily precipitation (mm)", main="Total daily precipitation over time",
     ylim=c(0,150), type="p")
dev.off()



# Get base model and meteorological matrices
avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=4), # Natural spline with 4 knots for mean temperature
                        arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for mean temperature lag

humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                       argvar=list(fun="ns", df=3), # Natural spline with 3 knots for mean temperature
                       arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for mean temperature lag

diftempmat <- crossbasis(df$diftemp, # No lag
                         argvar=list(fun="ns", df=4)) # Natural spline with 4 knots for daily temperature range

lagmat <- crossbasis(df$lag1, lag=3, # Lag of 3 days
                     argvar=list(fun="lin")) # Linear relationship for previous 3 days' values

model0 <- glm.nb(count ~ avtempmat + humidmat + diftempmat +
                   sin(2*pi*(yrday/(365))) + cos(2*pi*(yrday/(365))) +
                   splines::ns(dadmin, 5) + day + lagmat,
                 data = df)

model <- zeroinfl(count ~ avtempmat + humidmat + diftempmat +
                    sin(2*pi*(yrday/(365))) + cos(2*pi*(yrday/(365))) +
                    splines::ns(dadmin, 5) + day + lagmat | day + lagmat,
                  dist="negbin", data=df
)


# Plot example model fitted values
pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS1.pdf", width=15)
par(mfrow=c(1,1))
plot(df$dadmin[as.numeric(names(model0$y))], df$count[as.numeric(names(model0$y))], ylim=c(0,25),
     xlab="Date of admission", ylab="Daily hospital admissions")
     #main="Model fitted values plotted with daily hospital admission counts")
lines(df$dadmin[as.numeric(names(model0$y))], exp(model$fitted.values), col="red")
legend("topleft", legend=c("Daily counts", "Model fitted values"), col=c("black","red"), pch=c(1,NA), lty=c(NA,1))
dev.off()




# Plot model checks
res <- model$residuals

par(mfrow=c(1,1))
plot(df$dadmin[as.numeric(names(model0$y))], res,
     xlab="Date of admission", ylab="Model residual")

par(mfrow=c(2,2))
plot(df$avtemp[as.numeric(names(model0$y))], res,
     xlab="Mean daily temperature", ylab="Model residual")
plot(df$avhumid[as.numeric(names(model0$y))], res,
     xlab="Relative humidity", ylab="Model residual")
plot(df$diftemp[as.numeric(names(model0$y))], res,
     xlab="Daily temperature range", ylab="Model residual")
plot(df$day[as.numeric(names(model0$y))], res,
     xlab="Day of the week", ylab="Model residual")


par(mfrow=c(1,1))
pacf(model$residuals, main="Partial autocorrelation of model residuals")



# Run multiple imputations
imps <- idews.mi(pneum.df, reps=5, nb.size=0.4700373, nb.mu=16.6820675, seed=890422)
output <- pool.rubin(imps=imps, param="count")

# Run multiple imputations (using model imputations)
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE)
output <- pool.rubin(imps=imps$imps, param="count")
z.output <- pool.rubin(imps=imps$imps, param="zero")

save(imps, output, z.output, file="~/WASH/Pneumonia Manuscript/PneumoniaAnalyses/impmodels.RData")




###################################################################
########### Results ############
###################################################################

load("~/WASH/Pneumonia Manuscript/PneumoniaAnalyses/impmodels.RData")

# Produce table of model estimates
overdisp <- vector()
for (i in seq_along(imps$imps)) {
  overdisp <- append(overdisp, imps$imps[[i]]$theta)
}
mean(overdisp)


coef <- output$coef
names(coef) <- c("Intercept", "Mean Temp v1.l1", "Mean Temp v1.l2", "Mean Temp v1.l3",
                 "Mean Temp v2.l1", "Mean Temp v2.l2", "Mean Temp v2.l3",
                 "Mean Temp v3.l1", "Mean Temp v3.l2", "Mean Temp v3.l3",
                 "Mean Temp v4.l1", "Mean Temp v4.l2", "Mean Temp v4.l3",
                 "Humidity v1.l1", "Humidity v1.l2", "Humidity v1.l3",
                 "Humidity v2.l1", "Humidity v2.l2", "Humidity v2.l3",
                 "Humidity v3.l1", "Humidity v3.l2", "Humidity v3.l3",
                 "Temp Range v1", "Temp Range v2", "Temp Range v3", "Temp Range v4",
                 "sin(Annual)", "cos(Annual)", "LongTerm v1", "LongTerm v2", "LongTerm v3", "LongTerm v4", "LongTerm v5",
                 "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday",
                 "Admission Lag"
)

z.coef <- data.frame("Zero Estimate"=z.output$coef,
                     "Zero SE"=diag(z.output$vcov)^0.5,
                     "Zero p"=(1-pnorm(abs(z.output$coef/(diag(z.output$vcov)^0.5))))*2
                     )
z.coef <- rbind(z.coef[1,],
                as.data.frame(matrix(ncol=3, nrow=length(coef)-nrow(z.coef), dimnames=list(NULL, names(z.coef)))),
                z.coef[2:nrow(z.coef),])

coef <- data.frame("Parameter"=names(coef),
                   "NegBin Estimate"=coef,
                   "NegBin SE"=diag(output$vcov)^0.5,
                   "NegBin p"=(1-pnorm(abs(coef/(diag(output$vcov)^0.5))))*2,
                   row.names = NULL)
coef <- cbind(coef, z.coef, row.names=NULL)

coef[,2:ncol(coef)] <- round(coef[,2:ncol(coef)],3)
coef <- apply(coef, MARGIN=2, FUN=function(x) {trimws(as.character(x))})
coef <- apply(coef, MARGIN=2, FUN=function(x) {gsub("0\\.000", "<0.001", x)})

write.csv(coef, file="~/WASH/Pneumonia Manuscript/TableS1.csv")



# Avetemp
coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)

plotpred(pred, dens="avtemp", xlim=c(10,30), xlab="Mean daily temperature (°C)", logy=TRUE)
plotlag(pred, vars=c(12,28))
plotslice(pred, lags=c(0,7,14,21))
plotslice(pred, lags=c(0,5,10,15,21))

pdf("~/WASH/Pneumonia Manuscript/Graphs/Figure3.pdf")
plotpred(pred, dens="avtemp", xlim=c(10,30), xlab="Mean daily temperature (°C)", logy=TRUE)
dev.off()

pdf("~/WASH/Pneumonia Manuscript/Graphs/Figure4.pdf", width=10)
plotslice(pred, lags=c(0,7,14,21))
dev.off()

pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS2.pdf")
plotlag(pred, vars=c(12,28))
dev.off()


# Humid
coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
pred <- crosspred(humidmat, coef=coef.humid, vcov=vcov.humid, model.link="log", cen=67)

plotpred(pred, dens="avhumid", xlab="Relative humidity (%)", logy=TRUE)
plotlag(pred, vars=c(26,96))
plotslice(pred, lags=c(0,7,14,21), logy=FALSE, xlab = "Relative humidity (%)")


pdf("~/WASH/Pneumonia Manuscript/Graphs/Figure5.pdf")
plotpred(pred, dens="avhumid", xlab="Relative humidity (%)", logy=TRUE)
dev.off()

pdf("~/WASH/Pneumonia Manuscript/Graphs/Figure6.pdf", width=10)
plotslice(pred, lags=c(0,7,14,21), logy=FALSE, xlab = "Relative humidity (%)")
dev.off()

pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS3.pdf")
plotlag(pred, vars=c(24,92))
dev.off()




# Diftemp
coef.diftemp <- output$coef[grepl("diftempmat", names(output$coef))]
vcov.diftemp <- output$vcov[grepl("diftempmat", rownames(output$vcov)), grepl("diftempmat", colnames(output$vcov))]
pred <- crosspred(diftempmat, coef=coef.diftemp, vcov=vcov.diftemp, model.link="log", cen=1.3)

plotpred(pred, dens="diftemp", xlab="Daily temperature range (DTR) (°C)", logy=FALSE)


pdf("~/WASH/Pneumonia Manuscript/Graphs/Figure7.pdf")
plotpred(pred, dens="diftemp", xlab="Daily temperature range (DTR) (°C)", logy=FALSE)
dev.off()





#########################################
######### Attributable risk #########
#########################################


############# For manuscript ###########

af.df <- data.frame("Exposure"=NA, "AF"=NA, "Low"=NA, "High"=NA)

#### AF due to temperature ####
coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]

# For manuscript table
temp <- quantile(attrdl(df$avtemp,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="back", cen=21, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp,
               quantile(attrdl(df$avtemp,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp,
                               range=c(min(df$avtemp), 14), type="af", dir="back", cen=21, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- append(temp,
               quantile(attrdl(df$avtemp,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp,
                               range=c(26, max(df$avtemp)), type="af", dir="back", cen=21, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="Mean daily temperature",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
                           )

temp <- quantile(attrdl(df$avtemp+2,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="forw", cen=21, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp,
               quantile(attrdl(df$avtemp+2,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp,
                               range=c(min(df$avtemp), 14), type="af", dir="forw", cen=21, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- append(temp,
               quantile(attrdl(df$avtemp+2,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp,
                               range=c(26, max(df$avtemp)), type="af", dir="forw", cen=21, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="Mean daily temperature + 2°C",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
                           )

temp <- quantile(attrdl(df$avtemp-2,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="forw", cen=21, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp,
               quantile(attrdl(df$avtemp-2,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp,
                               range=c(min(df$avtemp), 14), type="af", dir="forw", cen=21, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- append(temp,
               quantile(attrdl(df$avtemp-2,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp,
                               range=c(26, max(df$avtemp)), type="af", dir="forw", cen=21, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="Mean daily temperature - 2°C",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
                           )




#### AF due to humidity ####
coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]

temp <- quantile(attrdl(df$avhumid,humidmat,df$count, coef=coef.humid,vcov=vcov.humid ,type="af", dir="back", cen=67, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp,
               quantile(attrdl(df$avhumid,humidmat,df$count, coef=coef.humid,vcov=vcov.humid,
                               range=c(min(df$avhumid), 40), type="af", dir="forw", cen=67, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- append(temp,
               quantile(attrdl(df$avhumid,humidmat,df$count, coef=coef.humid,vcov=vcov.humid,
                               range=c(80, max(df$avhumid)), type="af", dir="forw", cen=67, sim=T,nsim=1000),
                               probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="Relative humidity",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
                           )

temp <- quantile(attrdl(df$avhumid+5,humidmat, df$count, coef=coef.humid,vcov=vcov.humid ,type="af", dir="forw", cen=67, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp,
               quantile(attrdl(df$avhumid+5, humidmat,df$count, coef=coef.humid,vcov=vcov.humid,
                               range=c(min(df$avhumid), 40), type="af", dir="forw", cen=67, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- append(temp,
               quantile(attrdl(df$avhumid+5, humidmat,df$count, coef=coef.humid,vcov=vcov.humid,
                               range=c(80, max(df$avhumid)), type="af", dir="forw", cen=67, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="Relative humidity + 5%",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
)

temp <- quantile(attrdl(df$avhumid-5,humidmat,df$count, coef=coef.humid,vcov=vcov.humid ,type="af", dir="forw", cen=67, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp,
               quantile(attrdl(df$avhumid-5, humidmat,df$count, coef=coef.humid,vcov=vcov.humid,
                               range=c(min(df$avhumid), 40), type="af", dir="forw", cen=67, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- append(temp,
               quantile(attrdl(df$avhumid-5, humidmat,df$count, coef=coef.humid,vcov=vcov.humid,
                               range=c(80, max(df$avhumid)), type="af", dir="forw", cen=67, sim=T,nsim=1000),
                        probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="Relative humidity - 5%",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
)


#### AF due to DTR ####
coef.diftemp <- output$coef[grepl("diftempmat", names(output$coef))]
vcov.diftemp <- output$vcov[grepl("diftempmat", rownames(output$vcov)), grepl("diftempmat", colnames(output$vcov))]
quantile(attrdl(df$diftemp,diftempmat,df$count, coef=coef.diftemp,vcov=vcov.diftemp ,type="af",cen=1.3, sim=T,nsim=1000))

temp <- quantile(attrdl(df$diftemp,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,type="af", dir="back", cen=1.3, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp, quantile(attrdl(df$diftemp,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,
                                     range=c(min(df$diftemp), 1.3), type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                              probs = c(0.025,0.5,0.975)))
temp <- append(temp, quantile(attrdl(df$diftemp,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,
                                     range=c(20, max(df$diftemp)), type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                              probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="DTR",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
)

temp <- quantile(attrdl(df$diftemp+2,diftempmat, df$count, coef=coef.diftemp,vcov=vcov.diftemp ,type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp, quantile(attrdl(df$diftemp+2,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,
                                     range=c(min(df$diftemp), 1.3), type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                              probs = c(0.025,0.5,0.975)))
temp <- append(temp, quantile(attrdl(df$diftemp+2,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,
                                     range=c(20, max(df$diftemp)), type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                              probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="DTR + 2°C",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
)

temp <- quantile(attrdl(df$diftemp-2,diftempmat,df$count, coef=coef.diftemp,vcov=vcov.diftemp ,type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                 probs = c(0.025,0.5,0.975))
temp <- append(temp, quantile(attrdl(df$diftemp-2,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,
                                     range=c(min(df$diftemp), 1.3), type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                              probs = c(0.025,0.5,0.975)))
temp <- append(temp, quantile(attrdl(df$diftemp-2,diftempmat, df$count, coef=coef.diftemp, vcov=vcov.diftemp ,
                                     range=c(20, max(df$diftemp)), type="af", dir="forw", cen=1.3, sim=T,nsim=1000),
                              probs = c(0.025,0.5,0.975)))
temp <- round(temp*100, 1)
af.df <- af.df %>% add_row(Exposure="DTR - 2°C",
                           AF=paste0(temp[2], "% (", temp[1],"%, ", temp[3], "%)"),
                           Low=paste0(temp[5], "% (", temp[4],"%, ", temp[6], "%)"),
                           High=paste0(temp[8], "% (", temp[7],"%, ", temp[9], "%)")
)



# Write to file
af.df <- af.df[-1,]

write.csv(af.df, file="~/WASH/Pneumonia Manuscript/Table1.csv")





############## For exploration ###############

# AF due to temperature
coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]

quantile(attrdl(df$avtemp,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="back", cen=21, sim=T,nsim=1000),
        probs = c(0.025,0.5,0.975))

# due to temperature below 13
quantile(attrdl(df$avtemp,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="back",
                cen=21, range=c(min(df$avtemp), 13), sim=T,nsim=1000),
         probs = c(0.025,0.5,0.975))

# due to temperature >25
quantile(attrdl(df$avtemp,avtempmat,df$count, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="back",
                cen=21, range=c(25, max(df$avtemp)), sim=T,nsim=1000),
         probs = c(0.025,0.5,0.975))

# Due to a change in annual temp of increase 2 degrees (CHECK)
quantile(attrdl(df$avtemp+2,avtempmat,df$count+2, coef=coef.avtemp,vcov=vcov.avtemp ,type="af", dir="forw", cen=21, sim=T,nsim=1000),
         probs = c(0.025,0.5,0.975))



# AF due to humidity
coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
quantile(attrdl(df$avhumid,humidmat,df$count, coef=coef.humid,vcov=vcov.humid ,type="af",cen=21, sim=T,nsim=1000))








##############################################################################
##############################################################################

################ SENSITIVITY ANALYSES ##############

##############################################################################
##############################################################################

########## Using Minimum daily temperature ###########

imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = "mintemp")
output <- pool.rubin(imps=imps$imps, param="count")
z.output <- pool.rubin(imps=imps$imps, param="zero")

mintempmat <- crossbasis(df$mintemp, lag=21, # Lag of 21 days
                         argvar=list(fun="ns", df=4), # Natural spline with 4 knots for min temperature
                         arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for min temperature lag

coef.mintemp <- output$coef[grepl("mintempmat", names(output$coef))]
vcov.mintemp <- output$vcov[grepl("mintempmat", rownames(output$vcov)), grepl("mintempmat", colnames(output$vcov))]
pred <- crosspred(mintempmat, coef=coef.mintemp, vcov=vcov.mintemp, model.link="log", cen = 15)

plotpred(pred, dens="mintemp", xlab="Minimum daily temperature (°C)", logy=TRUE)
plotlag(pred, vars=c(12,28))
plotslice(pred, lags=c(0,7,14,21))

pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS7.pdf")
plotpred(pred, dens="mintemp", xlab="Minimum daily temperature (°C)", logy=TRUE)
dev.off()


########## Using Maximum daily temperature ###########

imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = "maxtemp")

output <- pool.rubin(imps=imps$imps, param="count")
z.output <- pool.rubin(imps=imps$imps, param="zero")

maxtempmat <- crossbasis(df$maxtemp, lag=21, # Lag of 21 days
                         argvar=list(fun="ns", df=4), # Natural spline with 4 knots for max temperature
                         arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for max temperature lag

coef.maxtemp <- output$coef[grepl("maxtempmat", names(output$coef))]
vcov.maxtemp <- output$vcov[grepl("maxtempmat", rownames(output$vcov)), grepl("maxtempmat", colnames(output$vcov))]
pred <- crosspred(maxtempmat, coef=coef.maxtemp, vcov=vcov.maxtemp, model.link="log", cen = 28)

plotpred(pred, dens="maxtemp", xlab="Maximum daily temperature (°C)", logy=TRUE)
plotlag(pred, vars=c(12,28))
plotslice(pred, lags=c(0,7,14,21))

pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS8.pdf")
plotpred(pred, dens="maxtemp", xlab="Maximum daily temperature (°C)", logy=TRUE)
dev.off()



########## Using different numbers of knots for mean daily temperature ###########

# Run with 2 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, tempknot=2)

output <- pool.rubin(imps=imps$imps, param="count")

avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=2), # Natural spline with 2 knots for mean temperature
                        arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for mean temperature lag

coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred2 <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)


# Run with 3 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, tempknot=3)

output <- pool.rubin(imps=imps$imps, param="count")

avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                         argvar=list(fun="ns", df=3), # Natural spline with 3 knots for mean temperature
                         arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for mean temperature lag

coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred3 <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)

# Run with 3 knots and 2 knots lag
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, tempknot=3, templag=2)

output <- pool.rubin(imps=imps$imps, param="count")

avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=3), # Natural spline with 3 knots for mean temperature
                        arglag=list(fun="ns", df=2)) # Natural spline with 2 knots for mean temperature lag

coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred3.2 <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)


# Run with 3 knots and 4 knots lag
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, tempknot=3, templag=4)

output <- pool.rubin(imps=imps$imps, param="count")

avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=3), # Natural spline with 3 knots for mean temperature
                        arglag=list(fun="ns", df=4)) # Natural spline with 4 knots for mean temperature lag

coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred3.4 <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)


# Run with 4 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, tempknot=4)

output <- pool.rubin(imps=imps$imps, param="count")

avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=4), # Natural spline with 4 knots for mean temperature
                        arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for mean temperature lag

coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred4 <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)


# Run with 5 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, tempknot=5)

output <- pool.rubin(imps=imps$imps, param="count")

avtempmat <- crossbasis(df$avtemp, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=5), # Natural spline with 5 knots for mean temperature
                        arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for mean temperature lag

coef.avtemp <- output$coef[grepl("avtempmat", names(output$coef))]
vcov.avtemp <- output$vcov[grepl("avtempmat", rownames(output$vcov)), grepl("avtempmat", colnames(output$vcov))]
pred5 <- crosspred(avtempmat, coef=coef.avtemp, vcov=vcov.avtemp, model.link="log", cen = 21)


g1 <- plotpred2(pred=pred4, pred.s1=pred3, pred.s2=pred5, adddens = FALSE,
          dens="avtemp", xlim=c(10,30), xlab="Mean daily temperature (°C)", logy=TRUE,
          labs=c("4 df", "3 df", "5 df")) +
  ggtitle("Exposure association")

g2 <- plotpred2(pred=pred4, pred.s1=pred3.2, pred.s2=pred3.4, adddens = FALSE,
                dens="avtemp", xlim=c(10,30), xlab="Mean daily temperature (°C)", logy=TRUE,
                labs=c("3 df", "2 df", "4 df")) +
  ggtitle("Lag association")

grid.newpage()
grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))


pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS9.pdf")
grid.newpage()
grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
dev.off()





########## Using different numbers of knots for relative humidity ###########

# Run with 2 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, humidknot=2)

output <- pool.rubin(imps=imps$imps, param="count")

humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                        argvar=list(fun="ns", df=2), # Natural spline with 2 knots for humidity
                        arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for humidity lag

coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
pred2.3 <- crosspred(humidmat, coef=coef.humid, vcov=vcov.humid, model.link="log", cen=67)



# Run with 3 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, humidknot=3)

output <- pool.rubin(imps=imps$imps, param="count")

humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                       argvar=list(fun="ns", df=3), # Natural spline with 3 knots for humidity
                       arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for humidity lag

coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
pred3.3 <- crosspred(humidmat, coef=coef.humid, vcov=vcov.humid, model.link="log", cen=67)



# Run with 4 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, humidknot=4)

output <- pool.rubin(imps=imps$imps, param="count")

humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                       argvar=list(fun="ns", df=4), # Natural spline with 4 knots for humidity
                       arglag=list(fun="ns", df=3)) # Natural spline with 3 knots for humidity lag

coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
pred4.3 <- crosspred(humidmat, coef=coef.humid, vcov=vcov.humid, model.link="log", cen=67)



# Run with 2 lag knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, humidlag=2)

output <- pool.rubin(imps=imps$imps, param="count")

humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                       argvar=list(fun="ns", df=3), # Natural spline with 3 knots for humidity
                       arglag=list(fun="ns", df=2)) # Natural spline with 2 knots for humidity lag

coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
pred3.2 <- crosspred(humidmat, coef=coef.humid, vcov=vcov.humid, model.link="log", cen=67)



# Run with 4 lag knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, humidlag=4)

output <- pool.rubin(imps=imps$imps, param="count")

humidmat <- crossbasis(df$avhumid, lag=21, # Lag of 21 days
                       argvar=list(fun="ns", df=3), # Natural spline with 4 knots for humidity
                       arglag=list(fun="ns", df=4)) # Natural spline with 4 knots for humidity lag

coef.humid <- output$coef[grepl("humidmat", names(output$coef))]
vcov.humid <- output$vcov[grepl("humidmat", rownames(output$vcov)), grepl("humidmat", colnames(output$vcov))]
pred3.4 <- crosspred(humidmat, coef=coef.humid, vcov=vcov.humid, model.link="log", cen=67)




g1 <- plotpred2(pred=pred3.3, pred.s1=pred2.3, pred.s2=pred4.3, adddens = FALSE,
                dens="avhumid", xlab="Relative humidity (%)", logy=TRUE,
                labs=c("3 df", "2 df", "4 df")) +
  ggtitle("Exposure association")

g2 <- plotpred2(pred=pred3.3, pred.s1=pred3.2, pred.s2=pred3.4, adddens = FALSE,
                dens="avhumid", xlab="Relative humidity (%)", logy=TRUE,
                labs=c("3 df", "2 df", "4 df")) +
  ggtitle("Lag association")

grid.newpage()
grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))


pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS10.pdf")
grid.newpage()
grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), size = "last"))
dev.off()





########## Using different numbers of knots for DTR ###########

# Run with 3 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, diftempknot=3)

output <- pool.rubin(imps=imps$imps, param="count")

diftempmat <- crossbasis(df$diftemp, # No lag
                         argvar=list(fun="ns", df=3)) # 3 knots


coef.diftemp <- output$coef[grepl("diftempmat", names(output$coef))]
vcov.diftemp <- output$vcov[grepl("diftempmat", rownames(output$vcov)), grepl("diftempmat", colnames(output$vcov))]
pred3 <- crosspred(diftempmat, coef=coef.diftemp, vcov=vcov.diftemp, model.link="log", cen=1.3)


# Run with 4 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, diftempknot=4)

output <- pool.rubin(imps=imps$imps, param="count")

diftempmat <- crossbasis(df$diftemp, # No lag
                         argvar=list(fun="ns", df=4)) # 4 knots


coef.diftemp <- output$coef[grepl("diftempmat", names(output$coef))]
vcov.diftemp <- output$vcov[grepl("diftempmat", rownames(output$vcov)), grepl("diftempmat", colnames(output$vcov))]
pred4 <- crosspred(diftempmat, coef=coef.diftemp, vcov=vcov.diftemp, model.link="log", cen=1.3)


# Run with 5 knots
imps <- idews.mi(pneum.df, reps=20, nb.size=0.4700373, nb.mu=16.6820675, seed=890422, imputemodel = TRUE,
                 sensitivity = NULL, diftempknot=5)

output <- pool.rubin(imps=imps$imps, param="count")

diftempmat <- crossbasis(df$diftemp, # No lag
                         argvar=list(fun="ns", df=5)) # 5 knots


coef.diftemp <- output$coef[grepl("diftempmat", names(output$coef))]
vcov.diftemp <- output$vcov[grepl("diftempmat", rownames(output$vcov)), grepl("diftempmat", colnames(output$vcov))]
pred5 <- crosspred(diftempmat, coef=coef.diftemp, vcov=vcov.diftemp, model.link="log", cen=1.3)



g1 <- plotpred2(pred=pred4, pred.s1=pred3, pred.s2=pred5, adddens = FALSE,
                dens="diftemp", xlab="Daily temperature range (DTR) (°C)", logy=TRUE,
                labs=c("4 df", "3 df", "5 df")) +
  ggtitle("Exposure association")


pdf("~/WASH/Pneumonia Manuscript/Graphs/FigureS11.pdf", height=5)
plot(g1)
dev.off()
