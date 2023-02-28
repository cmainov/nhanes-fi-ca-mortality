library( tidyverse )
library( glue )
library( splines )
library( survey )
library( survminer ) # for adjusted survival curves
library( ggsci )

source( "R/utils.R" ) # read in helper functions
source( "R/surv-miner-bug-fix.R" ) # bug fix for generating survival curves with `survminer`

d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds")


# x variables
indices <- c( "fs_enet", "age_enet", "hhs_enet", "fdas_enet", "pc1", "pc2" )




### Analyses on the Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                  "smokstat", "kcal", "weekmetmin", "education_bin",
                  "cci_score", "alc_cat", "ins.status", "binfoodsechh",
                  "foodasstpnowic" ) 


# all-cause mortality
out.res <- list()
fin.res <-data.frame()

for( i in 1:length( indices ) ){
  
  out.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
     subs = "inc == 1",    # subset of data to use
     cuts = 5,             # quantiles to use for categorization
     id.col = "seqn",      # subject id column
     covars = covars.surv, # covariates
     time = "stime",       # survival time column
     mort.ind = "mortstat",
     sample.name = "All Cancer Survivors" )    # mortality indicator column
  
  fin.res <- rbind( fin.res, out.res[[i]]$frame)
}


# cancer mortality
out.res.ca <- list()
fin.res.ca <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.ca[[i]] <- res( df = d, x = indices[i], 
                       subs = "inc == 1", 
                       cuts = 5, 
                       id.col = "seqn", 
                       covars = covars.surv, 
                       time = "stime", 
                       mort.ind = "castat",
                       sample.name = "All Cancer Survivors" ) 
  
  fin.res.ca <- rbind( fin.res.ca, out.res.ca[[i]]$frame)
}


# cvd mortality
out.res.cvd <- list()
fin.res.cvd <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.cvd[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "cvdstat",
                          sample.name = "All Cancer Survivors" ) 
  
  fin.res.cvd <- rbind( fin.res.cvd, out.res.cvd[[i]]$frame )
}




# ---------------------------------------------------------------------------------------------------------------------------------------------------------





### Analyses on the Food Insecure Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv.fi <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                     "smokstat", "kcal", "weekmetmin", "education_bin",
                     "cci_score", "alc_cat", "ins.status",
                     "foodasstpnowic"  ) 

# all-cause mortality
out.res.fi <- list()
fin.res.fi <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.fi[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c("inc == 1", "binfoodsechh == 'Low'"),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.fi, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       sample.name = "Food Insecure Cancer Survivors" )    # mortality indicator column
  
  fin.res.fi <- rbind( fin.res.fi, out.res.fi[[i]]$frame )
  
}



# cancer mortality
out.res.fi.ca <- list()
fin.res.ca.fi <- data.frame()

for( i in 1:length( indices ) ){
  
  out.res.fi.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c("inc == 1", "binfoodsechh == 'Low'"), 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.fi, 
                          time = "stime", 
                          mort.ind = "castat",
                          sample.name = "Food Insecure Cancer Survivors" ) 
  
  fin.res.ca.fi <- rbind( fin.res.ca.fi, out.res.fi.ca[[i]]$frame)
}


# cvd mortality
## NOTE: Too few deaths for this subanalysis, the model does not converge
out.res.fi.cvd <- list()
fin.res.cvd.fi <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.fi.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = c("inc == 1", "binfoodsechh == 'Low'"), 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.fi, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           sample.name = "Food Insecure Cancer Survivors" ) 
  
  fin.res.cvd.fi <- rbind( fin.res.cvd.fi, out.res.fi.cvd[[i]]$frame)
  
}




# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Further Adjust for NHANES ADL Score (Analysis on the Cancer Survivor Population) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv.s <- c( covars.surv, "adl.score" )

# all-cause mortality
out.res.s <- list()
fin.res.s <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.s[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = "inc == 1",    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.s, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       sample.name = "All Cancer Survivors" )    # mortality indicator column
  
  fin.res.s <- rbind( fin.res.s, out.res.s[[i]]$frame)
  
  }


# cancer mortality
out.res.s.ca <- list()
fin.res.s.ca <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.s.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.s, 
                          time = "stime", 
                          mort.ind = "castat",
                          sample.name = "All Cancer Survivors" ) 
  
  fin.res.s.ca <- rbind( fin.res.s.ca, out.res.s.ca[[i]]$frame)
}


# cvd mortality
out.res.s.cvd <- list()
fin.res.s.cvd <- data.frame()
for( i in 1:length( indices ) ){
  
  out.res.s.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = "inc == 1", 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.s, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           sample.name = "All Cancer Survivors" ) 
  
  fin.res.s.cvd <- rbind( fin.res.s.cvd, out.res.s.cvd[[i]]$frame )
  
}



# ---------------------------------------------------------------------------------------------------------------------------------------------------------





### Stratify by Time Since Diagnosis (Only in Cancer Survivor Sample Given Sample Size) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

time.levs <- levels(d$timecafactor) # levels of time since diagnosis factor
out.time.res <- data.frame() # initialie results storage frame

for( j in seq_along( time.levs ) ){
  
  # all-cause mortality
  out.res.time <- list()
  fin.res.time <- data.frame()
  for( i in 1:length( indices ) ){
    
    out.res.time[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                         subs = c( "inc == 1", paste0( "timecafactor == '", time.levs[j], "'" ) ),    # subset of data to use
                         cuts = 5,             # quantiles to use for categorization
                         id.col = "seqn",      # subject id column
                         covars = covars.surv, # covariates
                         time = "stime",       # survival time column
                         mort.ind = "mortstat",
                         sample.name = time.levs[j] )    # mortality indicator column
  
    fin.res.time <- rbind( fin.res.time, out.res.time[[i]]$frame)
    
    }
  
  
  # cancer mortality
  out.res.time.ca <- list()
  fin.res.time.ca <- data.frame()
  for( i in 1:length( indices ) ){
    
    out.res.time.ca[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                                 subs = c( "inc == 1", paste0( "timecafactor == '", time.levs[j], "'" ) ),    # subset of data to use
                                 cuts = 5,             # quantiles to use for categorization
                                 id.col = "seqn",      # subject id column
                                 covars = covars.surv, # covariates
                                 time = "stime",       # survival time column
                                 mort.ind = "castat",
                                 sample.name = time.levs[j] )    # mortality indicator column
    
    fin.res.time.ca <- rbind( fin.res.time.ca, out.res.time.ca[[i]]$frame)
  }
  
  


  # combine results for different types of mortality
  time.res <- bind_rows( data.frame( index = "All-Cause Mortality"),
                         fin.res.time, 
                        data.frame( index = "Cancer-Specific Mortality"),
                        fin.res.time.ca )
  
  # bind results
  out.time.res <- rbind( out.time.res, time.res )
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Sensitivity Analysis: Only Those within 5 Years of a Primary Cancer Diagnosis ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## All Cancer Survivors ##
# all-cause mortality
out.sens.res <- list()
fin.res.sens <- data.frame()
for( i in 1:length( indices ) ){
  
  out.sens.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c( "inc == 1", "timesincecadxmn <= 60" ),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       sample.name = "All Cancer Survivors" )    # mortality indicator column
  
  fin.res.sens <- rbind( fin.res.sens, out.sens.res[[i]]$frame )
}



# cancer mortality
out.sens.res.ca <- list()
fin.res.sens.ca <- data.frame()
for( i in 1:length( indices ) ){
  
  out.sens.res.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c( "inc == 1", "timesincecadxmn <= 60" ), 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "castat",
                          sample.name = "All Cancer Survivors" ) 
  
  fin.res.sens.ca <- rbind( fin.res.sens.ca, out.sens.res.ca[[i]]$frame )
  
}


# cvd mortality
out.sens.res.cvd <- list()
fin.res.sens.cvd <- list()
for( i in 1:length( indices ) ){
  
  out.sens.res.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = c( "inc == 1", "timesincecadxmn <= 60" ), 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           sample.name = "All Cancer Survivors" ) 
  
  fin.res.sens.cvd <- rbind( fin.res.sens.cvd, out.sens.res.cvd[[i]]$frame )
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------





### Assemble Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## All-cause mortality ##

ac.table <- bind_rows( out.res, out.res.fi ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "age_enet", "hhs_enet", 
                                    "fdas_enet", "pc1", "pc2" ) ),
                  sample )


## Cancer-specific mortality ##

ca.table <- bind_rows( out.res.ca, out.res.fi.ca ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "age_enet", "hhs_enet", 
                                    "fdas_enet", "pc1", "pc2" ) ),
           sample )


## CVD-specific mortality ##

cvd.table <- bind_rows( out.res.cvd ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI cvd survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "age_enet", "hhs_enet", 
                                    "fdas_enet", "pc1", "pc2" ) ),
           sample )


## Adjust for ADL Score Models ##
s.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
           fin.res.s, 
           data.frame( index = "Cancer-Specific Mortality"),
           fin.res.s.ca, 
           data.frame( index = "Cardiovascular Disease Mortality"),
           fin.res.s.cvd )

ac.table[ac.table=="fs_enet"] <- "Food Insecurity"
ca.table[ca.table=="fs_enet"] <- "Food Insecurity"
cvd.table[cvd.table=="fs_enet"] <- "Food Insecurity"

ac.table[ac.table=="age_enet"] <- "Age"
ca.table[ca.table=="age_enet"] <- "Age"
cvd.table[cvd.table=="age_enet"] <- "Age"

ac.table[ac.table=="fdas_enet"] <- "Food Assistance (SNAP)"
ca.table[ca.table=="fdas_enet"] <- "Food Assistance (SNAP)"
cvd.table[cvd.table=="fdas_enet"] <- "Food Assistance (SNAP)"

ac.table[ac.table=="hhs_enet"] <- "Household Size"
ca.table[ca.table=="hhs_enet"] <- "Household Size"
cvd.table[cvd.table=="hhs_enet"] <- "Household Size"

ac.table[ac.table=="pc1"] <- "Modified Western"
ca.table[ca.table=="pc1"] <- "Modified Western"
cvd.table[cvd.table=="pc1"] <- "Modified Western"

ac.table[ac.table=="pc2"] <- "Prudent"
ca.table[ca.table=="pc2"] <- "Prudent"
cvd.table[cvd.table=="pc2"] <- "Prudent"


## Generate one table with all causes of death ##

all.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      ca.table, 
                      data.frame( index = "Cardiovascular Disease Mortality"),
                      cvd.table )


## Time since diagnosis supplementary table ##

out.time.res[out.time.res=="fs_enet"] <- "Food Insecurity"

out.time.res[out.time.res=="age_enet"] <- "Age"

out.time.res[out.time.res=="fdas_enet"] <- "Food Assistance (SNAP)"

out.time.res[out.time.res=="hhs_enet"] <- "Household Size"

out.time.res[out.time.res=="pc1"] <- "Modified Western"

out.time.res[out.time.res=="pc2"] <- "Prudent"



## Save tables ##

write.table( ac.table, "04-Tables-Figures/tables/04-table-4-ac.txt", sep = "," )
write.table( ac.table, "04-Tables-Figures/tables/05-table-4-ca.txt", sep = "," )
write.table( ac.table, "04-Tables-Figures/tables/06-table-4-cvd.txt", sep = "," )
write.table( all.table, "04-Tables-Figures/tables/07-table-4-all.txt", sep = "," )
write.table( s.table, "04-Tables-Figures/tables/08-table-s2.txt", sep = "," )
write.table( out.time.res, "04-Tables-Figures/tables/09-table-s3.txt", sep = "," )



### Survival Curves ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## FI pattern survival curves
ggadjustedcurves( fit = out.sens.res[[1]]$q.obj,
                 variable ="fs_enet.q",
                 data = out.sens.res[[1]]$dat,
                 method = "conditional",
                 title = "Food Security Pattern",
                 font.title = c(16, "bold"),
                 legend.title = "Quintile",
                 font.legend = c(10, "bold"),
                 legend = c(0.14,0.25),
                 ylab = "Adjusted Survival Rate",
                 xlab = "Follow-up (Months)",
                 size = 0.6) +
  theme(text=element_text(family="Avenir") ) + 
  scale_color_ordinal()

ggsave( "04-Tables-Figures/figures/02a-fi-surv-curve.png", 
        height = 7.21, 
        width = 6.42 )

# snap pattern survival curves
ggadjustedcurves( fit = out.res[[4]]$q.obj,
                 variable = "fdas_enet.q",
                 data = out.res[[4]]$dat,
                 method = "conditional",
                 title = "SNAP Pattern",
                 font.title = c(16, "bold"),
                 legend.title = "Quintile",
                 font.legend = c(10, "bold"),
                 legend = c(0.14,0.25),
                 ylab = "Adjusted Survival Rate",
                 xlab = "Follow-up (Months)",
                 size = 0.6) +
  theme(text=element_text(family="Avenir") ) + 
  scale_color_ordinal()

ggsave( "04-Tables-Figures/figures/02b-snap-surv-curve.png", 
        height = 7.21, 
        width = 6.42 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

