# rm(list=ls())

library(jagsUI)
library(dplyr)
library(tidyverse)
library(lubridate)


#################################################################
########## BUGS CODE ############################################
#################################################################

# Define the model in the BUGS language and write a text file
sink("flathead.txt")
cat("
    model {
    
    # Likelihood: 
    # Level-1 of the model
    for (i in 1:74){ 
    y[i] ~ dnorm(mu[i], tau)               
    mu[i] <- beta0 + site[i] + month[i]     
    } 
    
     # Priors
    beta0 ~ dnorm(0, 0.0001)
    sigma.e ~ dunif(0, 1)
    sigma.site ~ dunif(0,2)
    sigma.month ~ dunif(0,1)
    sigma.site.month ~ dunif(0,1)
    
    # Derived quantities
    tau <- pow(sigma.e,-2) 
    tau.site <- pow(sigma.site,-2) 
    tau.month <- pow(sigma.month,-2) 
    tau.site.month <- pow(sigma.site.month,-2)
    
    
    } # end model
    ",fill = TRUE)
sink()


# Read in data
dat <- read.csv("../cleaned_data/ppcp.csv", header = TRUE)

dat_numeric <- dat %>%
  filter(!(SITE %in% c("FI2", "DU2", "HO1"))) %>%
  group_by(SITE, MONTH) %>%
  summarize(y_raw = sum(CONCENTRATION, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(site = as.numeric(SITE),
         month =as.numeric(MONTH),
         y = as.numeric(scale(y_raw, center = TRUE, scale = TRUE))) %>%
  select(y, site, month) %>%
  as.list()

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
ni <- 100000
nt <- 1
nb <- 5000
nc <- 2

start.time = Sys.time()         # Set timer 
# Call JAGS from R 
out <- jags(dat_numeric, inits, parameters, "flathead.txt", n.chains = nc, 
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
saveRDS(mcmcOut, file="../cleaned_data/ppcp_mcmc_out.rds")

library(ggmcmc)
library(gridExtra)
library(ggthemes)
library(coda)
out.mcmc <- as.mcmc(out)
S <- ggs(out.mcmc$samples)
# ggs_traceplot(S)
ggmcmc(S, file = "../figures_tables/mcmc_output.pdf")
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
