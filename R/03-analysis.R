library( tidyverse )
library( glue )
library( splines )
library( survey )
library( survminer )

source( "R/utils.R" )

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
for( i in 1:length( indices ) ){
  
  out.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
     subs = "inc == 1",    # subset of data to use
     cuts = 5,             # quantiles to use for categorization
     id.col = "seqn",      # subject id column
     covars = covars.surv, # covariates
     time = "stime",       # survival time column
     mort.ind = "mortstat",
     sample.name = "All Cancer Survivors" )    # mortality indicator column
}

fin.res <- do.call( "rbind", out.res )


# cancer mortality
out.res.ca <- list()
for( i in 1:length( indices ) ){
  
  out.res.ca[[i]] <- res( df = d, x = indices[i], 
                       subs = "inc == 1", 
                       cuts = 5, 
                       id.col = "seqn", 
                       covars = covars.surv, 
                       time = "stime", 
                       mort.ind = "castat",
                       sample.name = "All Cancer Survivors" ) 
}

fin.res.ca <- do.call( "rbind", out.res.ca )


# cvd mortality
out.res.cvd <- list()
for( i in 1:length( indices ) ){
  
  out.res.cvd[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "cvdstat",
                          sample.name = "All Cancer Survivors" ) 
}

fin.res.cvd <- do.call( "rbind", out.res.cvd )



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
for( i in 1:length( indices ) ){
  
  out.res.fi[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c("inc == 1", "binfoodsechh == 'Low'"),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.fi, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       sample.name = "Food Insecure Cancer Survivors" )    # mortality indicator column
}

fin.res.fi <- do.call( "rbind", out.res.fi )



# cancer mortality
out.res.fi.ca <- list()
for( i in 1:length( indices ) ){
  
  out.res.fi.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c("inc == 1", "binfoodsechh == 'Low'"), 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.fi, 
                          time = "stime", 
                          mort.ind = "castat",
                          sample.name = "Food Insecure Cancer Survivors" ) 
}

fin.res.ca.fi <- do.call( "rbind", out.res.fi.ca )


# cvd mortality
## NOTE: Too few deaths for this subanalysis, the model does not converge
out.res.fi.cvd <- list()
for( i in 1:length( indices ) ){
  
  out.res.fi.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = c("inc == 1", "binfoodsechh == 'Low'"), 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.fi, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           sample.name = "Food Insecure Cancer Survivors" ) 
}

fin.res.cvd.fi <- do.call( "rbind", out.res.fi.cvd )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Further Adjust for NHANES ADL Score (Analysis on the Cancer Survivor Population) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv.s <- c( covars.surv, "adl.score" )

# all-cause mortality
out.res.s <- list()
for( i in 1:length( indices ) ){
  
  out.res.s[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = "inc == 1",    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.s, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       sample.name = "All Cancer Survivors" )    # mortality indicator column
}

fin.res.s <- do.call( "rbind", out.res.s )


# cancer mortality
out.res.s.ca <- list()
for( i in 1:length( indices ) ){
  
  out.res.s.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.s, 
                          time = "stime", 
                          mort.ind = "castat",
                          sample.name = "All Cancer Survivors" ) 
}

fin.res.s.ca <- do.call( "rbind", out.res.s.ca )


# cvd mortality
out.res.s.cvd <- list()
for( i in 1:length( indices ) ){
  
  out.res.s.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = "inc == 1", 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.s, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           sample.name = "All Cancer Survivors" ) 
}

fin.res.s.cvd <- do.call( "rbind", out.res.s.cvd )



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


## Generate one table with all causes of death

all.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      ca.table, 
                      data.frame( index = "Cardiovascular Disease Mortality"),
                      cvd.table )


## Save tables ##

write.table( ac.table, "04-Tables-Figures/tables/02-table-4-ac.txt", sep = "," )
write.table( ac.table, "04-Tables-Figures/tables/03-table-4-ca.txt", sep = "," )
write.table( ac.table, "04-Tables-Figures/tables/04-table-4-cvd.txt", sep = "," )
write.table( all.table, "04-Tables-Figures/tables/05-table-4-all.txt", sep = "," )
write.table( s.table, "04-Tables-Figures/tables/06-table-s2.txt", sep = "," )

