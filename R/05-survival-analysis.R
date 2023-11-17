library( tidyverse )
library( glue )
library( splines )
library( survey )
library( survminer ) # for adjusted survival curves
library( ggsci )
library( latex2exp )

source( "R/utils.R" ) # read in helper functions
source( "R/old/surv-miner-bug-fix.R" ) # bug fix for generating survival curves with `survminer`

d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds")


# x variables
indices <- c( "fs_enet", "pc1", "pc2", "hei.2015" )


### Analyses on the Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                  "smokstat", "kcal", "weekmetmin", "education_bin",
                  "cci_score", "alc_cat", "ins.status", "binfoodsechh",
                  "foodasstpnowic" ) 

covars.base <- c( "race", "gender", "age" ) # basic model covariates

covars.null <- c() # null model covariates

covars.list <- list( covars.null, covars.base, covars.surv )

# all-cause mortality
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
       model.name = c( "Null Model", "Basic Model", "Full Model" )[j] )    # mortality indicator column
  
     fin.res <- rbind( fin.res, int.list[[j]]$frame )
     
     }
  
  model.index <- which( fin.res$model == "Full Model" )[1]
  out.res[[i]] <- int.list[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
    
}


# cancer mortality
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
                         model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) 
    
    fin.res.ca <- rbind( fin.res.ca, int.list.ca[[j]]$frame)
  }
  
  model.index <- which( fin.res.ca$model == "Full Model" )[1]
  out.res.ca[[i]] <- int.list[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### Analyses on the Food Insecure Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv.fi <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                     "smokstat", "kcal", "weekmetmin", "education_bin",
                     "cci_score", "ins.status",
                     "foodasstpnowic"  ) 

covars.list.fi <- list( covars.null, covars.base, covars.surv.fi )

# all-cause mortality
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
                       model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) 

    fin.res.fi <- rbind( fin.res.fi, int.list.fi[[j]]$frame)
  }
  
  model.index <- which( fin.res.fi$model == "Full Model" )[1]
  out.res.fi[[i]] <- int.list[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}



# cancer mortality
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
                             model.name = c( "Null Model", "Basic Model", "Full Model" )[j]) 
    
    fin.res.fi.ca <- rbind( fin.res.fi.ca, int.list.fi.ca[[j]]$frame)
  }
  
  model.index <- which( fin.res.fi.ca$model == "Full Model" )[1]
  out.res.fi.ca[[i]] <- int.list[[model.index]] # save full model results (to obtain spline model and create tables for null/basic and full models separately)
  
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### Further Adjust for NHANES ADL Score (Analysis on the Cancer Survivor Population) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv.s <- c( covars.surv, "adl.score" )

# all-cause mortality
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
                       model.name = "Full Model" )    # mortality indicator column
  
  fin.res.s <- rbind( fin.res.s, out.res.s[[i]]$frame)
  
  }


# cancer mortality
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



### Sensitivity Analysis: Only Those within 5 Years of a Primary Cancer Diagnosis ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## All Cancer Survivors ##

# all-cause mortality
out.sens.res <- list()
fin.res.sens <- data.frame()
for( i in seq_along( indices ) ){
  
  out.sens.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c( "test.set == 1", "timesincecadxmn <= 60" ),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 1,
                       sample.name = "All Cancer Survivors",
                       model.name = "Full Model" )    # mortality indicator column
  
  fin.res.sens <- rbind( fin.res.sens, out.sens.res[[i]]$frame )
}



# cancer mortality
out.sens.res.ca <- list()
fin.res.sens.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  out.sens.res.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c( "test.set == 1", "timesincecadxmn <= 60" ), 
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



### Assemble Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## All-cause mortality ##

# full model results only
ac.table <- bind_rows( fin.res %>% filter( model == "Full Model" ), 
                       fin.res.fi %>% filter( model == "Full Model" ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
                  sample )


## Cancer-specific mortality ##

ca.table <- bind_rows( fin.res.ca %>% filter( model == "Full Model" ), 
                       fin.res.fi.ca %>% filter( model == "Full Model" ) ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "pc1", "pc2", "hei.2015" ) ),
           sample )



## Adjust for ADL Score Models ##
s.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
           fin.res.s, 
           data.frame( index = "Cancer-Specific Mortality"),
           fin.res.s.ca )


## Sensitivity analysis excluding those >60 mos since dx ##
sens.60.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                            fin.res.sens, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      fin.res.sens.ca )


ac.table[ac.table == "fs_enet"] <- "Food Insecurity"
ca.table[ca.table == "fs_enet"] <- "Food Insecurity"
s.table[s.table == "fs_enet"] <- "Food Insecurity"
sens.60.table[sens.60.table == "fs_enet"] <- "Food Insecurity"

ac.table[ac.table == "pc1"] <- "Prudent #1"
ca.table[ca.table == "pc1"] <- "Prudent #1"
s.table[s.table == "pc1"] <- "Prudent #1"
sens.60.table[sens.60.table == "pc1"] <- "Prudent #1"

ac.table[ac.table == "pc2"] <- "Prudent #2"
ca.table[ca.table == "pc2"] <- "Prudent #2"
s.table[s.table == "pc2"] <- "Prudent #2"
sens.60.table[sens.60.table == "pc2"] <- "Prudent #2"

ac.table[ac.table == "hei.2015"] <- "HEI-2015"
ca.table[ca.table == "hei.2015"] <- "HEI-2015"
s.table[s.table == "hei.2015"] <- "HEI-2015"
sens.60.table[sens.60.table == "hei.2015"] <- "HEI-2015"


## Generate one table (main analysis) with all causes of death ##

all.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      ca.table )

## Generate one table with all causes of death for basic and null models only ##

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
ac.table.nb[ac.table.nb == "fs_enet"] <- "Food Insecurity"
ca.table.nb[ca.table.nb == "fs_enet"] <- "Food Insecurity"

ac.table.nb[ac.table.nb == "pc1"] <- "Prudent #1"
ca.table.nb[ca.table.nb == "pc1"] <- "Prudent #1"

ac.table.nb[ac.table.nb == "pc2"] <- "Prudent #2"
ca.table.nb[ca.table.nb == "pc2"] <- "Prudent #2"

ac.table.nb[ac.table.nb == "hei.2015"] <- "HEI-2015"
ca.table.nb[ca.table.nb == "hei.2015"] <- "HEI-2015"

# row bind results for this supplementary table
all.table.nb <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table.nb, 
                        data.frame( index = "Cancer-Specific Mortality"),
                        ca.table.nb )

## Save tables ##

write.table( ac.table, "04-Tables-Figures/tables/04-table-4-ac.txt", sep = "," )
write.table( ca.table, "04-Tables-Figures/tables/05-table-4-ca.txt", sep = "," )
write.table( all.table, "04-Tables-Figures/tables/06-table-4-all.txt", sep = "," )
write.table( s.table, "04-Tables-Figures/tables/09-table-s4.txt", sep = "," ) # ADL score sensitivity analysis results
write.table( sens.60.table, "04-Tables-Figures/tables/08-table-s3.txt", sep = "," ) # limit to those with cancer < 5 yrs ago results
write.table( all.table.nb, "04-Tables-Figures/tables/07-table-s1.txt", sep = "," ) # null and basic model results

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### Survival Curves and Spline Plots ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## spline plots ##
model.index <- 1:4 
diet.names <- c( "FI Pattern", "Western Pattern", 
                 "Mixed Pattern", "HEI-2015$^a$" ) # x-axis labels

sp.plots.list <- Map( function( x, y ){
  prud1.sp <- out.res[[x]]$spline.plot +
    xlab( unname( TeX( y ) ) ) +
    theme( legend.position = "none",
           text = element_text( family = "Avenir" ),
           axis.title.y = element_text( size = 13 ) ,
           axis.text.y = element_text( size = 10, color = "grey30" ),
           axis.title.x = element_text( size = 13 ) ,
           axis.text.x = element_text( size = 10, color = "grey30" ) ) +
    ylab( "" ) +
    coord_cartesian( ylim = c( 0.2, max = 1.8 ) ) 
}, x = model.index, y = diet.names )


## survival curves ##
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

# arrange into final figure
( sp.sc.comb <- ggarrange( ggarrange( sc.plots.list[[1]], sp.plots.list[[1]], nrow = 2, labels = list( "A", "B" ) ),
           ggarrange( sc.plots.list[[2]], sp.plots.list[[2]], nrow = 2 ),
           ggarrange( sc.plots.list[[3]], sp.plots.list[[3]], nrow = 2 ),
           ggarrange( sc.plots.list[[4]], sp.plots.list[[4]], nrow = 2 ),
           nrow = 1, ncol = 4 ) )

ggsave( "04-Tables-Figures/figures/03-surv-spline-comb.png",
        height = 8, 
        width = 13,
        plot = sp.sc.comb,
        dpi = 400 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

