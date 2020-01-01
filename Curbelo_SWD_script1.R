# Curbelo et al. 2020
# D. suzukii in Hawaii Forest
# 09 December 2019

# getwd
setwd("~/Google Drive/Manuscripts/Curbeloetal_HawaiiDsuzukii/suzukii_analysis/")

# Load libraries
library(car) # diagnostics
library(MuMIn) #model selection
library(date) #for julian date

# list files
list.files()

# set data
df <- read.csv("HavoFinal.csv")
names(df)

# create julian date
df$Date <- mdy.date(df$Mon0, df$Day0, df$Year0)
tmp <- as.POSIXlt(df$Date, format = "%d%b%y")
df$Julianday <- tmp$yday

# create end julian date
df$Date_end <- mdy.date(df$Mon1, df$Day1, df$Year1)
tmp <- as.POSIXlt(df$Date_end, format = "%d%b%y")
df$Julianday_end <- tmp$yday

# export data
write.csv(df, "HavoFinal_v3.csv")

# GLM
test.glm <-glm(df$Suzukii.Total~df$Elevation+df$Site.Type+df$Julianday+df$Day_out,
               family = poisson(), na.action = "na.fail")
test.null <- glm(df$Suzukii.Total~1, family = poisson)
summary(test.glm)

# how many residuals
test.glm$residuals #140

test.glm$deviance/test.glm$df.residual # 16.68396; huge overdispersion
influence.measures(test.glm) # 23, 29, 43, 65, 107, 113, 115, 125, 127, 129, 131, 133; n = `1`
options(max.print = 2000)
plot(test.glm)

# Compare a change in the parameter estimates when observations 50,124 is removed from the model
# use negative sign to indicate you want to update the frog.glm model with observation 47 removed
# not a large change suggesting that inclusion is not a problem
# first round, 12 observations removed, therefore 128 records retained...91% of data used.
compareCoefs(test.glm, update(test.glm, subset=-c(23,29, 43, 65, 107, 113,
                                                  115, 125, 127, 129, 131,
                                                  133)))

# import new data and run model
df <- read.csv("HavoFinal_2.csv")
nrow(df)

# Generate Julian Day
# create julian date
df$Date <- mdy.date(df$Mon0, df$Day0, df$Year0)
tmp <- as.POSIXlt(df$Date, format = "%d%b%y")
df$Julianday <- tmp$yday

# create end julian date
df$Date_end <- mdy.date(df$Mon1, df$Day1, df$Year1)
tmp <- as.POSIXlt(df$Date_end, format = "%d%b%y")
df$Julianday_end <- tmp$yday

# Examines the fit of the model compared to a null model 
## See large difference with summary output but this provides a p-value
anova(test.glm, test.null, test = "Chisq") # significantly different than null
dd <- dredge(test.glm)
dd

### Don't necessarily need to account for overdispersion, but for practice
# Notice that dispersion parameter is set at 2.39
# notice that AIC is set to NA bc log likelihood is not calculated
over.testGLM <- update(test.glm, family=quasipoisson)
summary(over.testGLM)

# function to be used as an argument to extract overdispersion parameter from a model
## (see below example of how it can be used as an argument)
dfun <- function(object) {with(object,sum((weights * residuals^2)[weights > 0])/df.residual)}

# Model Selection using QAICc with function dredge 
## see quasi model selection in R_Boelker.pdf for more details and how to do in other packages
# must change the family name to x.quaasipoisson for MuMIn package
## from ?QAIC - A 'hacked' constructor for quasibinomial family object, 
### that allows for ML estimation
x.quasipoisson <- function(...) {res <- quasipoisson(...)
res$aic <- poisson(...)$aic
res}

# Update model to include new family = x.quasipoisson
ms.over.testGLM <- update(over.testGLM, family = "x.quasipoisson",
                          na.action=na.fail)

# Examine all possible combinations of models based on parameters in global model
## Asking R to use the Quasi likelihood which incorportates variance inflation factor
## also use rank = "QAIC"
## chat = overdispersion parameter -- extracting the value from the original model
dd <- dredge(ms.over.testGLM, rank = "QAICc", chat=dfun(test.glm))
dd

# shown only to illistrate how function is calculated specfic for this example
dfun.test <- function(test.glm) {with(teset.glm,sum((weights * residuals^2)[weights > 0])/df.residual)}

# also just indicate the value for overdispersion based on the output for over.redstartGLM model
# same as dd -- shown as an example
dd2 <- dredge(ms.over.testGLM, rank = "QAICc", chat=11.33288)
dd2

write.csv(dd2, "dd2_quasipoison.csv")

# Model Average Parameter estimates for subset of models
# see MuMln package for details
# use model-avg coefficients list and Std. Error = unconditional se / ignore p-values
top.models <- get.models(dd2, subset = delta <2)
x <-model.avg(top.models)
coefs <- model.avg(dd, subset = delta < 4, revised.var= T)
summary(coefs) # examine full average only, notice the relative variable importance values here!
coefs <-top.models$`15`$coefficients # isolate just parameter estimates
test.glm <-glm(df$Suzukii.Total~df$Elevation+df$Site.Type+df$Julianday,
               family = poisson(), na.action = "na.fail")

# Create plots
# forested site and lava site
nrow(df)
df.for <-df[ which(df$habitat=="Forested"),] 
df.lav <-df[ which(df$habitat=="Non-forested"),] 

# Figure 1 was made in ArcGIS

# Figure 2 generation
tiff(filename = "/Users/jonathankoch/Google Drive/Manuscripts/Curbeloetal_HawaiiDsuzukii/suzukii_analysis/Fig2.tiff",width = 2000, height = 1000, 
     units = "px", pointsize =12, res = 300)
par(mfrow=c(1,2))
plot(df.for$Suzukii.Total~df.for$Elevation, pch = 16, xlab = "Altitude (m)",
     ylab = "D. suzukii abundance")
points(df.lav$Suzukii.Total~df.lav$Elevation, pch = 15, col = "gray")
boxplot(df$Suzukii.Total~df$habitat, xlab = "Habitat type",
        ylab = "D. suzukii abundance")
dev.off()