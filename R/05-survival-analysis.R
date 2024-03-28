###---------------------------------------------------
###   05-VALIDATION ANALYSIS/SURVIVAL ANALYSIS
###---------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we perform the validation survival analysis on the testing
# subsample of the dataset. Note this script makes use of several results-
# generating functions we wrote and that are stored in the utils.R file.
# 
# INPUT DATA FILE: 
# i."03-Data-Rodeo/01-analytic-data.rds"
#
# OUTPUT FILES: 
# Several results files (tables--.txt files--and figures--.png and .tiff files)
#
# Resources: 
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( tidyverse )
library( glue )      # for string gluing
library( splines )   # for spline functions
library( survey )    # survey commands
library( survminer ) # for adjusted survival curves
library( ggsci )     # pallettes
library( latex2exp ) # for latex in plots

source( "R/utils.R" ) # read in helper functions
source( "R/old/surv-miner-bug-fix.R" ) # bug fix for generating survival curves with `survminer`


### (0.0) Read-in Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds")

# specify x variables in the analyses
indices <- c( "fs_enet", "pc1", "pc2", "hei.2015" )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) Analyses on the Entire Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (1.1) Specify covariates to include in model ##
covars.surv <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                  "smokstat", "kcal", "weekmetmin", "education_bin",
                  "cci_score", "alc_cat", "ins.status", "binfoodsechh",
                  "foodasstpnowic" ) 

covars.base <- c( "race", "gender", "age" ) # basic model covariates

covars.null <- c() # null model covariates

covars.list <- list( covars.null, covars.base, covars.surv )

## ---o--- ##


## (1.2) All-cause mortality ##

out.res <- list()
fin.res <- data.frame()

for( i in seq_along( indices ) ){
  
  int.list <- list() # hold intermediary results (i.e., looping over covariate set)
  for( j in 1:length( covars.list ) ){
    
     int.list[[j]] <- res( df = d, x = indices[i],   # data and x variable for model
       subs = "test.set == 1",    # subset of data to use
       cuts = 5,             # quantiles to use for categorization
       id.col = "seqn",      # subject id column
       covars = unlist( covars.list[j] ), # covariates
       time = "stime",       # survival time column
       mort.ind = "mortstat",
       scale.y = 1.5, # for shifting y axis max value
       int.knots = 1, # interior knots for spline models
       sample.name = "All Cancer Survivors",
       model.name = c( "Null Model", "Basic Model", "Full Model" )[j] )    # sequential adjustment
  
     fin.res <- rbind( fin.res, int.list[[j]]$frame )
     
     }
  
  model.index <- which( fin.res$model == "Full Model" )[1]
  out.res[[i]] <- int.list[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
    
}

## ---o--- ##


## (1.3) Cancer mortality ##

out.res.ca <- list()
fin.res.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  int.list.ca <- list() # hold intermediary results (i.e., looping over covariate set)
  for( j in 1:length( covars.list ) ){
    
    int.list.ca[[j]] <- res( df = d, x = indices[i], 
                         subs = "test.set == 1", 
                         cuts = 5, 
                         id.col = "seqn", 
                         covars = unlist( covars.list[j] ), 
                         time = "stime", 
                         mort.ind = "castat",
                         scale.y = 1.5, # for shifting y axis max value
                         int.knots = 1,
                         sample.name = "All Cancer Survivors",
                         model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) # sequential adjustment
    
    fin.res.ca <- rbind( fin.res.ca, int.list.ca[[j]]$frame)
  }
  
  model.index <- which( fin.res.ca$model == "Full Model" )[1]
  out.res.ca[[i]] <- int.list.ca[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) Analyses on the Food Insecure Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Specify covariates to include in model ##

covars.surv.fi <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                     "smokstat", "kcal", "weekmetmin", "education_bin",
                     "cci_score", "ins.status",
                     "foodasstpnowic"  ) 

covars.list.fi <- list( covars.null, covars.base, covars.surv.fi )

## ---o--- ##


## (2.2) All-cause mortality ##
out.res.fi <- list()
fin.res.fi <- data.frame()

for( i in seq_along( indices ) ){
  
  int.list.fi <- list() # hold intermediary results (i.e., looping over covariate set)
  for( j in 1:length( covars.list ) ){
    
    int.list.fi[[j]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c("test.set == 1", "binfoodsechh == 'Low'"),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = unlist( covars.list.fi[j] ), # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 1,
                       sample.name = "Food Insecure Cancer Survivors",
                       model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) # sequential adjustment

    fin.res.fi <- rbind( fin.res.fi, int.list.fi[[j]]$frame)
  }
  
  model.index <- which( fin.res.fi$model == "Full Model" )[1]
  out.res.fi[[i]] <- int.list.fi[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}

## ---o--- ##


## (2.3) Cancer mortality ##

out.res.fi.ca <- list()
fin.res.fi.ca <- data.frame()

for( i in seq_along( indices ) ){
  
  int.list.fi.ca <- list() # hold intermediary results (i.e., looping over covariate set)
  for( j in 1:length( covars.list ) ){
    
    int.list.fi.ca[[j]] <- res( df = d, x = indices[i],   # data and x variable for model
                             subs = c("test.set == 1", "binfoodsechh == 'Low'"),    # subset of data to use
                             cuts = 5,             # quantiles to use for categorization
                             id.col = "seqn",      # subject id column
                             covars = unlist( covars.list.fi[j] ), # covariates
                             time = "stime",       # survival time column
                             mort.ind = "castat",
                             scale.y = 1.3, # for shifting y axis max value
                             int.knots = 1,
                             sample.name = "Food Insecure Cancer Survivors",
                             model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) # sequential adjustment
    
    fin.res.fi.ca <- rbind( fin.res.fi.ca, int.list.fi.ca[[j]]$frame)
  }
  
  model.index <- which( fin.res.fi.ca$model == "Full Model" )[1]
  out.res.fi.ca[[i]] <- int.list.fi.ca[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Further Adjust for NHANES ADL Score (Analysis on the Cancer Survivor Population) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Specify covariates to include in model ##
covars.surv.s <- c( covars.surv, "adl.score" )

## ---o--- ##


## (3.2) All-cause mortality ##

out.res.s <- list()
fin.res.s <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.s[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = "test.set == 1",    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.s, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 1,
                       sample.name = "All Cancer Survivors",
                       model.name = "Full Model" )    # full model
  
  fin.res.s <- rbind( fin.res.s, out.res.s[[i]]$frame)
  
  }

## ---o--- ##


## (3.3) Cancer mortality ##
out.res.s.ca <- list()
fin.res.s.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.s.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = "test.set == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.s, 
                          time = "stime", 
                          mort.ind = "castat",
                          scale.y = 1.3, # for shifting y axis max value
                          int.knots = 1,
                          sample.name = "All Cancer Survivors",
                          model.name = "Full Model") 
  
  fin.res.s.ca <- rbind( fin.res.s.ca, out.res.s.ca[[i]]$frame)
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Sensitivity Analysis: Only Those within 4 Years of a Primary Cancer Diagnosis ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## All Cancer Survivors ##

## (4.1 ) All-cause mortality ##

out.sens.res <- list()
fin.res.sens <- data.frame()
for( i in seq_along( indices ) ){
  
  out.sens.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c( "test.set == 1", "timesincecadxmn <= 48" ),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 1,
                       sample.name = "All Cancer Survivors",
                       model.name = "Full Model" )    # full model
  
  fin.res.sens <- rbind( fin.res.sens, out.sens.res[[i]]$frame )
}

## ---o--- ##


## (4.2) Cancer mortality ##

out.sens.res.ca <- list()
fin.res.sens.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  out.sens.res.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c( "test.set == 1", "timesincecadxmn <= 48" ), 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "castat",
                          scale.y = 1.3, # for shifting y axis max value
                          int.knots = 1,
                          sample.name = "All Cancer Survivors",
                          model.name = "Full Model" ) 
  
  fin.res.sens.ca <- rbind( fin.res.sens.ca, out.sens.res.ca[[i]]$frame )
  
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) Assemble Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) All-cause mortality results ##

# full model results only
ac.table <- bind_rows( fin.res %>% filter( model == "Full Model" ), 
                       fin.res.fi %>% filter( model == "Full Model" ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
                  sample )

## ---o--- ##


## (5.2) Cancer-specific mortality ##

ca.table <- bind_rows( fin.res.ca %>% filter( model == "Full Model" ), 
                       fin.res.fi.ca %>% filter( model == "Full Model" ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
           sample )

## ---o--- ##


## (5.3) Adjust for ADL Score Models ##
s.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
           fin.res.s, 
           data.frame( index = "Cancer-Specific Mortality"),
           fin.res.s.ca )

## ---o--- ##


## (5.4) Sensitivity analysis excluding those >60 mos since dx ##
sens.48.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                            fin.res.sens, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      fin.res.sens.ca )

## ---o--- ##


## (5.5) Rename patterns in table ##
ac.table[ac.table == "fs_enet"] <- "Pattern #1"
ca.table[ca.table == "fs_enet"] <- "Pattern #1"
s.table[s.table == "fs_enet"] <- "Pattern #1"
sens.48.table[sens.48.table == "fs_enet"] <- "Pattern #1"

ac.table[ac.table == "pc1"] <- "Pattern #2"
ca.table[ca.table == "pc1"] <- "Pattern #2"
s.table[s.table == "pc1"] <- "Pattern #2"
sens.48.table[sens.48.table == "pc1"] <- "Pattern #2"

ac.table[ac.table == "pc2"] <- "Pattern #3"
ca.table[ca.table == "pc2"] <- "Pattern #3"
s.table[s.table == "pc2"] <- "Pattern #3"
sens.48.table[sens.48.table == "pc2"] <- "Pattern #3"

ac.table[ac.table == "hei.2015"] <- "HEI-2015"
ca.table[ca.table == "hei.2015"] <- "HEI-2015"
s.table[s.table == "hei.2015"] <- "HEI-2015"
sens.48.table[sens.48.table == "hei.2015"] <- "HEI-2015"

## ---o--- ##


## (5.6) Generate one table (main analysis) with all causes of death ##

all.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      ca.table )

## ---o--- ##


## (5.7) Generate one table with all causes of death for basic and null models only ##

# null and basic model results only
ac.table.nb <- bind_rows( fin.res %>% filter( model %in% c( "Null Model", "Basic Model" ) ), 
                       fin.res.fi %>% filter( model %in% c( "Null Model", "Basic Model" ) ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
           sample )


# cancer-specific mortality 

ca.table.nb <- bind_rows( fin.res.ca %>% filter( model %in% c( "Null Model", "Basic Model" ) ), 
                       fin.res.fi.ca %>% filter( model %in% c( "Null Model", "Basic Model" ) ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
           sample )

# change names inside table
ac.table.nb[ac.table.nb == "fs_enet"] <- "Pattern #1"
ca.table.nb[ca.table.nb == "fs_enet"] <- "Pattern #1"

ac.table.nb[ac.table.nb == "pc1"] <- "Pattern #2"
ca.table.nb[ca.table.nb == "pc1"] <- "Pattern #2"

ac.table.nb[ac.table.nb == "pc2"] <- "Pattern #3"
ca.table.nb[ca.table.nb == "pc2"] <- "Pattern #3"

ac.table.nb[ac.table.nb == "hei.2015"] <- "HEI-2015"
ca.table.nb[ca.table.nb == "hei.2015"] <- "HEI-2015"

# row bind results for this supplementary table
all.table.nb <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table.nb, 
                        data.frame( index = "Cancer-Specific Mortality"),
                        ca.table.nb )

## ---o--- ##


## (5.8) Save tables ##

write.table( ac.table, "04-Tables-Figures/tables/04-table-4-ac.txt", sep = "," )
write.table( ca.table, "04-Tables-Figures/tables/05-table-4-ca.txt", sep = "," )
write.table( all.table, "04-Tables-Figures/tables/06-table-4-all.txt", sep = "," )
write.table( s.table, "04-Tables-Figures/tables/09-table-s4.txt", sep = "," ) # ADL score sensitivity analysis results
write.table( sens.48.table, "04-Tables-Figures/tables/08-table-s3.txt", sep = "," ) # limit to those with cancer < 4 yrs ago results
write.table( all.table.nb, "04-Tables-Figures/tables/07-table-s1.txt", sep = "," ) # null and basic model results

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (6.0) Survival Curves and Spline Plots ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (6.1) Spline plots ##

model.index <- 1:4 
diet.names <- c( "Pattern #1$^{\\dagger a}$", "Pattern #2$^{\\ddagger b}$", 
                 "Pattern #3$^{\\ddagger c}$", "HEI-2015$^d$" ) # x-axis labels

sp.plots.list <- Map( function( x, y ){
  
  prud1.sp <- out.res[[x]]$spline.plot +
    xlab( unname( TeX( y ) ) ) +
    theme( legend.position = "none",
           text = element_text( family = "Avenir" ),
           axis.title.y = element_text( size = 13 ) ,
           axis.text.y = element_text( size = 10, color = "grey30" ),
           axis.title.x = element_text( size = 15 ) ,
           axis.text.x = element_text( size = 10, color = "grey30" ) ) +
    ylab( "" ) +
    coord_cartesian( ylim = c( 0.2, max = 1.8 ) ) 
}, 
x = model.index, y = diet.names )

## ---o--- ##


## (6.2) Survival curves ##
vars.q <- c( "fs_enet.q", "pc1.q", "pc2.q", "hei.2015.q" ) # cut variable names

leg.pos <- list( c(0.15,0.35), "none", "none", "none" ) # legend positions

sc.plots.list <- Map( function( x, y, z, l ){
  ggadjustedcurves( fit =  out.res[[x]]$q.obj,
                    variable = y,
                    data = out.res[[x]]$dat,
                    method = "conditional",
                    title = unname( TeX( z ) ),
                    font.title = c(16, "bold"),
                    legend.title = "Quintile",
                    font.legend = c(10, "bold"),
                    legend = c(0.14,0.31),
                    ylab = "",
                    xlab = "Follow-up (Months)",
                    size = 0.6) +
    theme(text=element_text(family="Avenir") ) + 
    scale_color_ordinal() + 
    theme( legend.position = l,
           text = element_text( family = "Avenir" ),
           axis.title.y = element_text( size = 13 ) ,
           axis.text.y = element_text( size = 10, color = "grey30" ),
           axis.title.x = element_text( size = 13 ) ,
           axis.text.x = element_text( size = 10, color = "grey30" ) )
},
x = model.index,
y = vars.q,
z = diet.names,
l = leg.pos)

## ---o--- ##


## (6.3) Arrange into final figure ##
( sp.sc.comb <- ggarrange( ggarrange( sc.plots.list[[1]] + ylab( "Survival Probability" ) + theme( axis.title.y = element_text(margin = margin( t = 0, r = 9, b = 0, l = 0 ) ) ), # add y axis label for first plot in row with spacing between axis text and title 
                                      sp.plots.list[[1]] + ylab( "Hazard Ratio (HR)" ), # add y axis label for first plot in row
                                      nrow = 2, labels = list( "A", "B" ) ),
           ggarrange( sc.plots.list[[2]], sp.plots.list[[2]], nrow = 2 ),
           ggarrange( sc.plots.list[[3]], sp.plots.list[[3]], nrow = 2 ),
           ggarrange( sc.plots.list[[4]], sp.plots.list[[4]], nrow = 2 ),
           nrow = 1, ncol = 4 ) )

## ---o--- ##


## (6.4) Save Figures ##
ggsave( "04-Tables-Figures/figures/03-surv-spline-comb.png",
        height = 8, 
        width = 13,
        plot = sp.sc.comb,
        dpi = 800 )

ggsave( "04-Tables-Figures/figures/03-surv-spline-comb.tiff",
        height = 8, 
        width = 13,
        plot = sp.sc.comb,
        dpi = 500 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

