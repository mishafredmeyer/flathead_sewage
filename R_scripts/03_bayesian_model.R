# rm(list=ls())

library(jagsUI)
library(dplyr)
library(tidyverse)
library(lubridate)
library(PKPDmisc)


#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("flathead_ppcp.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:80){ 
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- beta0 + site[i] + month[i]     
    } 
    
    
     # Priors
    beta0 ~ dnorm(0, 0.0001)
    sigma.e ~ dunif(0, 1)
    sigma.site ~ dunif(0,1)
    sigma.month ~ dunif(0,1)
    
    # Derived quantities
    tau <- pow(sigma.e,-2) 
    tau.site <- pow(sigma.site,-2) 
    tau.month <- pow(sigma.month,-2) 

    
    
    } # end model
    ",fill = TRUE)
sink()


# Read in data for PPCPs
dat <- read.csv("../cleaned_data/ppcp.csv", header = TRUE)

dat_numeric <- dat %>%
  filter(!(SITE %in% c("FI2", "DU2", "HO1"))) %>%
  filter(SITE %in% c("BD", "FLBS", "SA", "SJ", "BO", "LK", "WF",
                     "WS", "DA", "FI")) %>%
  mutate(time_point = ifelse(MONTH == "MAY", "MAY_1", NA),
         time_point = ifelse(MONTH == "JUNE" & between(TIME, 1, 6), "MAY_1", time_point),
         time_point = ifelse(MONTH == "JUNE" & between(TIME, 7, 25), "JUNE_1", time_point),
         time_point = ifelse(MONTH == "JUNE" & between(TIME, 26, 30), "JUNE_2", time_point),
         time_point = ifelse(MONTH == "JULY" & between(TIME, 1, 11), "JUNE_2", time_point),
         time_point = ifelse(MONTH == "JULY" & between(TIME, 12, 23), "JULY_1", time_point),
         time_point = ifelse(MONTH == "JULY" & between(TIME, 24, 31), "JULY_2", time_point),
         time_point = ifelse(MONTH == "AUGUST" & between(TIME, 1, 8), "JULY_2", time_point),
         time_point = ifelse(MONTH == "AUGUST" & between(TIME, 9, 21), "AUGUST_1", time_point),
         time_point = ifelse(MONTH == "AUGUST" & between(TIME, 22, 31), "AUGUST_2", time_point),
         time_point = ifelse(MONTH == "SEPTEMBER" & between(TIME, 1, 5), "AUGUST_2", time_point),
         time_point = ifelse(MONTH == "SEPTEMBER" & between(TIME, 6, 23), "SEPTEMBER_1", time_point),
         time_point = ifelse(MONTH == "SEPTEMBER" & between(TIME, 24, 30), "SEPTEMBER_2", time_point),
         time_point = ifelse(MONTH == "OCTOBER", "SEPTEMBER_2", time_point),
         time_point = as.factor(time_point)) %>%
  group_by(SITE, time_point) %>%
  summarize(y_raw = sum(CONCENTRATION, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(y = log10(y_raw+1)) %>%
  #filter(y <= 2) %>%
  #mutate(y = scale(y_raw, scale = TRUE, center = TRUE)) %>%
  mutate(site = as.numeric(SITE),
         month =as.numeric(time_point)) %>%
  select(-time_point, -SITE) %>%
  #strip_attributes(c("scaled:center", "scaled:scale")) %>%
  mutate(y = as.vector(y)) %>%
  select(-y_raw) %>%
  dplyr::select(y, site, month) %>%
  as.list()

#### Playing with variance ratios

ESS <- dat
  
dat_numeric$time_point <- factor(dat_numeric$time_point, levels = c("MAY_1",
                                                                    "JUNE_1", "JUNE_2",
                                                                    "JULY_1", "JULY_2",
                                                                    "AUGUST_1", "AUGUST_2",
                                                                    "SEPTEMBER_1", "SEPTEMBER_2"))

attr(dat_numeric$y, "attr") <- NULL

dat_numeric <- merTools::stripAttributes(dat_numeric)
# Read in dat for fatty acids
# dat <- read.csv("../cleaned_data/fatty_acids.csv", header = TRUE)
# 
# SAFA <- c("C12.0", "C14.0", "C15.0", "C16.0", "C17.0", "C18.0", "C20.0", "C22.0", "iso.C15.0", "C24.0", "C26.0", "C28.0")
# MUFA <- c( "C15.1", "C14.1n5", "C15.1w7", "C16.1w5", "C16.1w6", "C16.1w7", "C16.1w7c",  "C16.1w8", "C16.1w9", "C17.1n7", "C18.1w7", "C18.1w7c", 
#            "C18.1w9", "C18.1w9c", "C20.1w7", "C20.1w9", "C22.1w7", "C22.1w9", "C22.1w9c")
# LUFA <- c("C16.2", "C16.2w4",  "C16.2w6",  "C16.2w7",  "C16.3w3",  "C16.3w4",  "C16.3w6",  "C18.2w6",  "C18.2w6t", "C18.2w6c",
#           "C18.3w3",  "C18.3w6", "C20.2w6",  "C20.3w3",  "C20.3w6",  "C22.2w6",  "C22.3w3")
# HUFA <- c("C16.4w1", "C16.4w3", "C18.4w3", "C18.4w4", "C18.5w3", "C20.4w2", "C20.4w3", "C20.4w6", "C20.5w3", "C22.4w3", 
#           "C22.4w6", "C22.5w3", "C22.5w6", "C22.6w3")
# SCUFA_LUFA <- c("C16.2", "C16.2w4",  "C16.2w6",  "C16.2w7",  "C16.3w3",  "C16.3w4",  "C16.3w6",  "C18.2w6", "C18.2w6c", 
#                 "C18.2w6t", "C18.3w3",  "C18.3w6")
# LCUFA_LUFA <- c("C20.2w6",  "C20.3w3",  "C20.3w6",  "C22.2w6",  "C22.3w3")
# SCUFA_HUFA <- c("C16.4w1", "C16.4w3", "C18.4w3", "C18.4w4", "C18.5w3")
# LCUFA_HUFA <- c( "C20.4w2", "C20.4w3", "C20.4w6", "C20.5w3", "C22.4w3", 
#                  "C22.4w6", "C22.5w3", "C22.5w6", "C22.6w3")
# SCUFA <- c("C16.2w4",  "C16.2w6",  "C16.2w7",  "C16.3w3",  "C16.3w4",  "C16.3w6",  "C18.2w6",  "C18.2w6t", "C18.2w6t",
#            "C18.3w3",  "C18.3w6", "C16.4w1", "C16.4w3", "C18.4w3", "C18.4w4", "C18.5w3")
# LCUFA <- c( "C20.4w2", "C20.4w3", "C20.4w6", "C20.5w3", "C22.4w3", 
#             "C22.4w6", "C22.5w3", "C22.5w6", "C22.6w3", "C20.2w6",  "C20.3w3",  "C20.3w6",  "C22.2w6",  "C22.3w3")
# 
# C18PUFA <- c("C18.2w6",  "C18.2w6t", "C18.3w3",  "C18.3w6", "C18.4w3", "C18.4w4", "C18.5w3")
# C20PUFA <- c("C20.4w3", "C20.4w6", "C20.5w3", "C20.2w6",  "C20.3w3",  "C20.3w6")
# 
# dat_numeric <- dat %>%
#   filter(!(LOC %in% c("FI2", "DU2", "HO1"))) %>%
#   rename("month" = "MONTH") %>%
#   gather(fatty_acid, conc, C12.0:C28.0) %>%
#   filter(fatty_acid != "C19.0",
#          !(fatty_acid %in% MUFA)) %>%
#   mutate(fatty_acid_group = ifelse(fatty_acid %in% SAFA, "SAFA", NA),
#          fatty_acid_group = ifelse(fatty_acid %in% c(LUFA, HUFA), "PUFA", fatty_acid_group)) %>%
#   group_by(LOC, month, fatty_acid_group) %>%
#   summarize(total_conc = sum(conc)) %>%
#   ungroup() %>%
#   spread(fatty_acid_group, total_conc) %>%
#   mutate(y_raw = PUFA/SAFA) %>% 
#   mutate(site = as.numeric(LOC),
#          month =as.numeric(month),
#          y = as.numeric(y_raw)) %>%
#   select(y, site, month) %>%
#   as.list()


str(dat_numeric)
attach(dat_numeric)



# Initial values
inits <- function (){
  list (beta0 = rnorm(1), sigma.e=runif(1), sigma.site=runif(1),sigma.month=runif(1), 
        sigma.site.month=runif(1))
}


# Parameters monitored
parameters <- c("beta0","sigma.e","sigma.site", "sigma.month", "sigma.site.month",
                "local.b","global.b","a")


# MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 2

start.time = Sys.time()         # Set timer 
# Call JAGS from R 
out <- jags(data = dat_numeric, inits, parameters, "flathead_ppcp.txt", n.chains = nc, 
            n.thin = nt, n.iter = ni, n.burnin = nb,  parallel=T)
end.time = Sys.time()
elapsed.time = round(difftime(end.time, start.time, units='mins'), dig = 2)
cat('Posterior computed in ', elapsed.time, ' minutes\n\n', sep='') 
# Calculate computation time


# Summarize posteriors
print(out, dig = 3)

### Sometimes you have many, many parameters to examine:
# Find which parameters, if any, have Rhat > 1.1
which(out$summary[, c("Rhat")] > 1.1)

# Or see what max Rhat value is
max(out$summary[, c("Rhat")])

# outExp <- out$summary
# write.csv(outExp, "TP_ModelSummary.csv", row.names = T)
mcmcOut <- out$sims.list
saveRDS(mcmcOut, file="../cleaned_data/fatty_acid_mcmc_out.rds")

library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(coda)
out.mcmc <- as.mcmc(out)
S <- ggs(out.mcmc$samples)
# ggs_traceplot(S)
ggmcmc(S, file = "../figures_tables/fatty_acid_mcmc_output.pdf")
# For select parameters

# hist(out$sims.list$tau2)


# sum1 <- out$summary
# write.csv(sum1,"summary_exp.csv", row.names = T)

# traceplot(out,parameters = "a")
# traceplot(out,parameters = "lambda1")
# traceplot(out,parameters = "sigma.year")
# 
# dat$lake <- as.numeric(as.factor(as.numeric(dat$lagoslakeid)))
# dat$year <- as.numeric(as.factor(as.numeric(dat$sampleyear)))
# dat$huc <- as.numeric(as.factor(as.numeric(as.factor(dat$hu4_zoneid))))
# write.csv(dat,"tp_dat_diagnostics.csv",row.names = FALSE)
