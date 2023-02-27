library( tidyverse )
library( glue )
library( splines )
library( survey )
library( survminer )

source( "R/utils.R" )

d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds")

# covariates to include in model
covars.surv <- c( 'race', 'gender', 'age', 'bmxbmi', 'hhsize', "fipr",
                            'smokstat', 'kcal', 'weekmetmin', 'education_bin',
                            'cci_score', "alc_cat", "ins.status", "binfoodsechh",
                  "foodasstpnowic" ) 

# x variables
indices <- c( "fs_enet", "age_enet", "hhs_enet", "fdas_enet", "pc1", "pc2" )




### Analyses on the Cancer population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

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
covars.surv.fi <- c( 'race', 'gender', 'age', 'bmxbmi', 'hhsize', "fipr",
                     'smokstat', 'kcal', 'weekmetmin', 'education_bin',
                     'cci_score', "alc_cat", "ins.status",
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
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.fi, 
                          time = "stime", 
                          mort.ind = "castat",
                          sample.name = "Food Insecure Cancer Survivors" ) 
}

fin.res.ca.fi <- do.call( "rbind", out.res.fi.ca )


# cvd mortality
out.res.fi.cvd <- list()
for( i in 1:length( indices ) ){
  
  out.res.fi.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = "inc == 1", 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.fi, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           sample.name = "Food Insecure Cancer Survivors" ) 
}

fin.res.cvd.fi <- do.call( "rbind", out.res.fi.cvd )
View(fin.res.cvd.fi)

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Assemble Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## All-cause mortality ##
ac.table <- bind_rows( out.res, out.res.fi ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c('fs_enet', 'age_enet', 'hhs_enet', 
                                    'fdas_enet', 'pc1', 'pc2' ) ),
                  sample )
