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
                            'cci_score', "alc_cat", "ins.status", "binfoodsechh" ) 

# x variables
indices <- c( "fs_enet", "age_enet", "hhs_enet", "fdas_enet", "pc1", "pc2" )




### Analyses on the Cancer population
# all-cause mortality
out.res <- list()
for( i in 1:length( indices ) ){
  
  out.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
     subs = "inc == 1",    # subset of data to use
     cuts = 5,             # quantiles to use for categorization
     id.col = "seqn",      # subject id column
     covars = covars.surv, # covariates
     time = "stime",       # survival time column
     mort.ind = "mortstat" )    # mortality indicator column
}

fin.res <- do.call( "rbind", out.res )
View(fin.res)



# cancer mortality
out.res.ca <- list()
for( i in 1:length( indices ) ){
  
  out.res.ca[[i]] <- res( df = d, x = indices[i], 
                       subs = "inc == 1", 
                       cuts = 5, 
                       id.col = "seqn", 
                       covars = covars.surv, 
                       time = "stime", 
                       mort.ind = "castat" ) 
}

# diabetes mortality
out.res.dm <- list()
for( i in 1:length( indices ) ){
  
  out.res.dm[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "dmstat" ) 
}

fin.res.dm <- do.call( "rbind", out.res.dm )
View(fin.res.dm)


# cvd mortality
out.res.cvd <- list()
for( i in 1:length( indices ) ){
  
  out.res.cvd[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "cvdstat" ) 
}

fin.res.cvd <- do.call( "rbind", out.res.cvd )
View(fin.res.cvd)
fin.res.cvd$pcubic
