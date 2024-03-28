###-----------------------------------------------------
###   06-PROPENSITY SCORE MATCHING SENSITIVITY ANALYSIS
###-----------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we perform the propensity score matching sensitivity analysis
# by generating matched sample based on a series of covariate and then rerunning
# the main analyses. We use the `MatchIt` function to carry out aspects of this
# analysis.
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
library( MatchIt )   # propensity score matching

source( "R/utils.R" ) # read in helper functions
source( "R/old/surv-miner-bug-fix.R" ) # bug fix for generating survival curves with `survminer`


### (0.0) Matching ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (0.1) Read in data and create binary indicators for high low scores for propensity score models ##

d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds") %>%
  mutate( fs.enet.bin = ifelse( .$fs_enet >= median( .$fs_enet, na.rm = TRUE ), 1,
                                ifelse( .$fs_enet < median( .$fs_enet, na.rm = TRUE ), 0, NA ) ),
          
          pc1.bin = ifelse( .$pc1 >= median( .$pc1, na.rm = TRUE ), 1,
                                ifelse( .$pc1 < median( .$pc1, na.rm = TRUE ), 0, NA ) ),
          pc2.bin = ifelse( .$pc2 >= median( .$pc2, na.rm = TRUE ), 1,
                                ifelse( .$pc2 < median( .$pc2, na.rm = TRUE ), 0, NA ) ),
          hei.2015.bin = ifelse( .$hei.2015 >= median( .$hei.2015, na.rm = TRUE ), 1,
                                ifelse( .$hei.2015 < median( .$hei.2015, na.rm = TRUE ), 0, NA ) )) %>%
  select( -contains( "agedx" ) )

# these are the covariates we removed subjects on if they were missing
these.covs <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                 "smokstat", "kcal", "weekmetmin", "education_bin",
                 "cci_score", "alc_cat", "ins.status", "binfoodsechh",
                 "foodasstpnowic" ) 

## ---o--- ##


## (0.2) Sample inclusion indicators ##

# create new indicator variable identifying those that made it to the final sample
# and those that could have been eligible for final sample had they not had missing values

# check which covariates are complete so that we can fit a model on those
covs.missing <- sapply( d[ which(d$test.set == 1) ,these.covs], function( x ) sum(is.na( x ) ) )

covs.not.missing <- covs.missing[ which( covs.missing == 0 ) ] %>% names()

## ---o--- ##


## (0.2) Logit model evaluating Pr[ Selection | covariates ] ##

# fit logistic model evaluating probability of selection given covariates

# "treatment" variables
treats <- c( "fs.enet.bin","pc1.bin", "pc2.bin", "hei.2015.bin" )
trts <- treats %>% str_remove_all(., ".bin" )

d.out <- d # copy dataframe for looping and minimize overwriting

for( i in seq_along( treats ) ){
  
  f.1 <- paste0( treats[i], " ~", paste0( covs.not.missing, collapse = " + " ) )
  
  d.2 <- d %>%
    filter( test.set == 1 ) %>% # subset testing sample
    mutate( new.wts = wtdr18yr / mean( .$wtdr18yr ) ) # normalize weights
  
  # fit propensity score model and conduct nearest neighbor matching
  s.out1 <- matchit( formula( f.1 ), data = d.2,
                    s.weights = "new.wts", # use normalized weights
                    distance = "glm",      # calculates propesnity score with logit model
                    discard = "both" )      # define area of common overlap in propensity scores across treated and untreated groups before matching
  
  
  
  d.3 <- d.2 %>%
    mutate( rowid = rownames(.) ) %>%
    filter( rowid %in% c( s.out1$match.matrix[!is.na( s.out1$match.matrix ) ], 
                          names( s.out1$match.matrix[ !is.na( s.out1$match.matrix) , 1 ] ) ) ) %>%
    mutate( !!paste0( "match.in.", trts[i] ) := 1 ) %>%
    select( -rowid )
  
  d <- left_join( d, d.3 %>% select( seqn, !!paste0( "match.in.", trts[i] ) ),
                  by = "seqn")

}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) Analyses on the Entire Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# x variables
indices <- c( "fs_enet", "pc1", "pc2", "hei.2015" )

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

out.res.ps <- list()
fin.res <- data.frame()

for( i in seq_along( indices ) ){
  
  int.list.ps <- list() # hold intermediary results (i.e., looping over covariate set)
  for( j in 1:length( covars.list ) ){
    
    int.list.ps[[j]] <- res( df = d, x = indices[i],   # data and x variable for model
                          subs = paste0( paste0( "match.in.", trts[i] )," == 1" ),    # subset of data to use which is the column with the matched subset
                          cuts = 5,             # quantiles to use for categorization
                          id.col = "seqn",      # subject id column
                          covars = unlist( covars.list[j] ), # covariates
                          time = "stime",       # survival time column
                          mort.ind = "mortstat",
                          scale.y = 1.5, # for shifting y axis max value
                          int.knots = 1, # interior knots for spline models
                          sample.name = "All Cancer Survivors",
                          model.name = c( "Null Model", "Basic Model", "Full Model" )[j] )    # sequential adjustment
    
    fin.res <- rbind( fin.res, int.list.ps[[j]]$frame )
    
  }
  
  model.index <- which( fin.res$model == "Full Model" )[1]
  out.res.ps[[i]] <- int.list.ps[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}


## ---o--- ##


## (1.3) Cancer mortality ##
out.res.ps.ca <- list()
fin.res.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  int.list.ps.ca <- list() # hold intermediary results (i.e., looping over covariate set)
  for( j in 1:length( covars.list ) ){
    
    int.list.ps.ca[[j]] <- res( df = d, x = indices[i], 
                             subs = paste0( paste0( "match.in.", trts[i] )," == 1" ), 
                             cuts = 5, 
                             id.col = "seqn", 
                             covars = unlist( covars.list[j] ), 
                             time = "stime", 
                             mort.ind = "castat",
                             scale.y = 1.5, # for shifting y axis max value
                             int.knots = 1,
                             sample.name = "All Cancer Survivors",
                             model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) # sequential adjustment
    
    fin.res.ca <- rbind( fin.res.ca, int.list.ps.ca[[j]]$frame)
  }
  
  model.index <- which( fin.res.ca$model == "Full Model" )[1]
  out.res.ps.ca[[i]] <- int.list.ps.ca[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
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
                             subs = c(paste0( paste0( "match.in.", trts[i] )," == 1" ), "binfoodsechh == 'Low'"),    # subset of data to use
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
                                subs = c(paste0( paste0( "match.in.", trts[i] )," == 1" ), "binfoodsechh == 'Low'"),    # subset of data to use
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



### (3.0) Assemble Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) All-cause mortality results ##

# full model results only
ac.table.match <- bind_rows( fin.res %>% filter( model == "Full Model" ), 
                       fin.res.fi %>% filter( model == "Full Model" ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
           sample )

## ---o--- ##


## (3.2) Cancer-specific mortality ##

ca.table.match <- bind_rows( fin.res.ca %>% filter( model == "Full Model" ), 
                       fin.res.fi.ca %>% filter( model == "Full Model" ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
           sample )


# rename patterns in tables
ac.table.match[ac.table.match == "fs_enet"] <- "Food Insecurity"
ca.table.match[ca.table.match == "fs_enet"] <- "Food Insecurity"

ac.table.match[ac.table.match == "pc1"] <- "Prudent #1"
ca.table.match[ca.table.match == "pc1"] <- "Prudent #1"

ac.table.match[ac.table.match == "pc2"] <- "Prudent #2"
ca.table.match[ca.table.match == "pc2"] <- "Prudent #2"

ac.table.match[ac.table.match == "hei.2015"] <- "HEI-2015"
ca.table.match[ca.table.match == "hei.2015"] <- "HEI-2015"

## ---o--- ##


## (3.6) Generate one table (main analysis) with all causes of death ##

all.table.match <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table.match, 
                        data.frame( index = "Cancer-Specific Mortality"),
                        ca.table.match )

write.table( all.table.match, "04-Tables-Figures/tables/10-table-ps-match-all.txt", sep = "," )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Survival Curves and Spline Plots ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.1) Spline plots ##

model.index <- 1:4 
diet.names <- c( "FI Pattern", "Western Pattern", 
                 "Mixed Pattern", "HEI-2015$^a$" ) # x-axis labels

sp.plots.list <- Map( function( x, y ){
  
  prud1.sp <- out.res.ps[[x]]$spline.plot +
    xlab( unname( TeX( y ) ) ) +
    theme( legend.position = "none",
           text = element_text( family = "Avenir" ),
           axis.title.y = element_text( size = 13 ) ,
           axis.text.y = element_text( size = 10, color = "grey30" ),
           axis.title.x = element_text( size = 13 ) ,
           axis.text.x = element_text( size = 10, color = "grey30" ) ) +
    ylab( "" ) +
    coord_cartesian( ylim = c( 0.2, max = 1.8 ) ) 
}, 
x = model.index, y = diet.names )

## ---o--- ##


## (4.2) Survival curves ##

vars.q <- c( "fs_enet.q", "pc1.q", "pc2.q", "hei.2015.q" ) # cut variable names

leg.pos <- list( c(0.15,0.35), "none", "none", "none" ) # legend positions

sc.plots.list.ps <- Map( function( x, y, z, l ){
  ggadjustedcurves( fit =  out.res.ps[[x]]$q.obj,
                    variable = y,
                    data = out.res.ps[[x]]$dat,
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


## (4.3) Arrange into final figure ##

( sp.sc.comb.ps <- ggarrange( ggarrange( sc.plots.list.ps[[1]], sp.plots.list[[1]], nrow = 2, labels = list( "A", "B" ) ),
                           ggarrange( sc.plots.list.ps[[2]], sp.plots.list[[2]], nrow = 2 ),
                           ggarrange( sc.plots.list.ps[[3]], sp.plots.list[[3]], nrow = 2 ),
                           ggarrange( sc.plots.list.ps[[4]], sp.plots.list[[4]], nrow = 2 ),
                           nrow = 1, ncol = 4 ) )

ggsave( "04-Tables-Figures/figures/05-surv-spline-comb-ps-match.png",
        height = 8, 
        width = 13,
        plot = sp.sc.comb,
        dpi = 400 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


