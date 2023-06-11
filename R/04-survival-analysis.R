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
indices <- c( "fs_enet", "age_enet", "hhs_enet", "fdas_enet", "pc1", "pc2", "hei.2015" )




### Analyses on the Cancer Survivor Population ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# covariates to include in model
covars.surv <- c( "race", "gender", "age", "bmxbmi", "hhsize", "fipr",
                  "smokstat", "kcal", "weekmetmin", "education_bin",
                  "cci_score", "alc_cat", "ins.status", "binfoodsechh",
                  "foodasstpnowic" ) 


# all-cause mortality
out.res <- list()
fin.res <- data.frame()

for( i in seq_along( indices ) ){
  
  out.res[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
     subs = "inc == 1",    # subset of data to use
     cuts = 5,             # quantiles to use for categorization
     id.col = "seqn",      # subject id column
     covars = covars.surv, # covariates
     time = "stime",       # survival time column
     mort.ind = "mortstat",
     scale.y = 1.5, # for shifting y axis max value
     int.knots = 2,
     sample.name = "All Cancer Survivors" )    # mortality indicator column
  
  fin.res <- rbind( fin.res, out.res[[i]]$frame )
}


# cancer mortality
out.res.ca <- list()
fin.res.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.ca[[i]] <- res( df = d, x = indices[i], 
                       subs = "inc == 1", 
                       cuts = 5, 
                       id.col = "seqn", 
                       covars = covars.surv, 
                       time = "stime", 
                       mort.ind = "castat",
                       scale.y = 1.5, # for shifting y axis max value
                       int.knots = 2,
                       sample.name = "All Cancer Survivors" ) 
  
  fin.res.ca <- rbind( fin.res.ca, out.res.ca[[i]]$frame)
}


# cvd mortality
out.res.cvd <- list()
fin.res.cvd <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.cvd[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "cvdstat",
                          scale.y = 1.5, # for shifting y axis max value
                          int.knots = 2,
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
for( i in seq_along( indices ) ){
  
  out.res.fi[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = c("inc == 1", "binfoodsechh == 'Low'"),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.fi, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 2,
                       sample.name = "Food Insecure Cancer Survivors" )    # mortality indicator column
  
  fin.res.fi <- rbind( fin.res.fi, out.res.fi[[i]]$frame )
  
}



# cancer mortality
out.res.fi.ca <- list()
fin.res.ca.fi <- data.frame()

for( i in seq_along( indices ) ){
  
  out.res.fi.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c("inc == 1", "binfoodsechh == 'Low'"), 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.fi, 
                          time = "stime", 
                          mort.ind = "castat",
                          scale.y = 1.3, # for shifting y axis max value
                          int.knots = 2,
                          sample.name = "Food Insecure Cancer Survivors" ) 
  
  fin.res.ca.fi <- rbind( fin.res.ca.fi, out.res.fi.ca[[i]]$frame)
}


# cvd mortality
## NOTE: Too few deaths for this subanalysis, the model does not converge
out.res.fi.cvd <- list()
fin.res.cvd.fi <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.fi.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = c("inc == 1", "binfoodsechh == 'Low'"), 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.fi, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           scale.y = 1.3, # for shifting y axis max value
                           int.knots = 2,
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
for( i in seq_along( indices ) ){
  
  out.res.s[[i]] <- res( df = d, x = indices[i],   # data and x variable for model
                       subs = "inc == 1",    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv.s, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 2,
                       sample.name = "All Cancer Survivors" )    # mortality indicator column
  
  fin.res.s <- rbind( fin.res.s, out.res.s[[i]]$frame)
  
  }


# cancer mortality
out.res.s.ca <- list()
fin.res.s.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.s.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = "inc == 1", 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv.s, 
                          time = "stime", 
                          mort.ind = "castat",
                          scale.y = 1.3, # for shifting y axis max value
                          int.knots = 2,
                          sample.name = "All Cancer Survivors" ) 
  
  fin.res.s.ca <- rbind( fin.res.s.ca, out.res.s.ca[[i]]$frame)
}


# cvd mortality
out.res.s.cvd <- list()
fin.res.s.cvd <- data.frame()
for( i in seq_along( indices ) ){
  
  out.res.s.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = "inc == 1", 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv.s, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           scale.y = 1.3, # for shifting y axis max value
                           int.knots = 2,
                           sample.name = "All Cancer Survivors" ) 
  
  fin.res.s.cvd <- rbind( fin.res.s.cvd, out.res.s.cvd[[i]]$frame )
  
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
                       subs = c( "inc == 1", "timesincecadxmn <= 60" ),    # subset of data to use
                       cuts = 5,             # quantiles to use for categorization
                       id.col = "seqn",      # subject id column
                       covars = covars.surv, # covariates
                       time = "stime",       # survival time column
                       mort.ind = "mortstat",
                       scale.y = 1.3, # for shifting y axis max value
                       int.knots = 2,
                       sample.name = "All Cancer Survivors" )    # mortality indicator column
  
  fin.res.sens <- rbind( fin.res.sens, out.sens.res[[i]]$frame )
}



# cancer mortality
out.sens.res.ca <- list()
fin.res.sens.ca <- data.frame()
for( i in seq_along( indices ) ){
  
  out.sens.res.ca[[i]] <- res( df = d, x = indices[i], 
                          subs = c( "inc == 1", "timesincecadxmn <= 60" ), 
                          cuts = 5, 
                          id.col = "seqn", 
                          covars = covars.surv, 
                          time = "stime", 
                          mort.ind = "castat",
                          scale.y = 1.3, # for shifting y axis max value
                          int.knots = 2,
                          sample.name = "All Cancer Survivors" ) 
  
  fin.res.sens.ca <- rbind( fin.res.sens.ca, out.sens.res.ca[[i]]$frame )
  
}


# cvd mortality
out.sens.res.cvd <- list()
fin.res.sens.cvd <- list()
for( i in seq_along( indices ) ){
  
  out.sens.res.cvd[[i]] <- res( df = d, x = indices[i], 
                           subs = c( "inc == 1", "timesincecadxmn <= 60" ), 
                           cuts = 5, 
                           id.col = "seqn", 
                           covars = covars.surv, 
                           time = "stime", 
                           mort.ind = "cvdstat",
                           scale.y = 1.3, # for shifting y axis max value
                           int.knots = 2,
                           sample.name = "All Cancer Survivors" ) 
  
  fin.res.sens.cvd <- rbind( fin.res.sens.cvd, out.sens.res.cvd[[i]]$frame )
}


# ---------------------------------------------------------------------------------------------------------------------------------------------------------





### Assemble Tables ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## All-cause mortality ##

ac.table <- bind_rows( fin.res, fin.res.fi ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "age_enet", "hhs_enet", 
                                    "fdas_enet", "pc1", "pc2" ) ),
                  sample )


## Cancer-specific mortality ##

ca.table <- bind_rows( fin.res.ca, fin.res.ca.fi ) %>%
  data.frame() %>%
  
  # arrange tables first by custom order and second by sample so that ALL Survivors are situated next to estimates for FI CA survivors for a given diet index
  arrange( factor(index, levels = c("fs_enet", "age_enet", "hhs_enet", 
                                    "fdas_enet", "pc1", "pc2" ) ),
           sample )


## CVD-specific mortality ##

cvd.table <- bind_rows( fin.res.cvd ) %>%
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


## Sensitivity analysis excluding those >60 mos since dx ##
sens.60.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                            fin.res.sens, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      fin.res.sens.ca, 
                      data.frame( index = "Cardiovascular Disease Mortality"),
                      fin.res.sens.cvd )

# rename the dietary patterns in these tables
ac.table[ac.table=="fs_enet"] <- "Food Insecurity"
ca.table[ca.table=="fs_enet"] <- "Food Insecurity"
cvd.table[cvd.table=="fs_enet"] <- "Food Insecurity"
s.table[s.table=="fs_enet"] <- "Food Insecurity"
sens.60.tablesens.60.table <- "Food Insecurity"


ac.table[ac.table=="age_enet"] <- "Age"
ca.table[ca.table=="age_enet"] <- "Age"
cvd.table[cvd.table=="age_enet"] <- "Age"
s.table[s.table=="age_enet"] <- "Age"
sens.60.table[sens.60.table=="age_enet"] <- "Age"

ac.table[ac.table=="fdas_enet"] <- "Food Assistance (SNAP)"
ca.table[ca.table=="fdas_enet"] <- "Food Assistance (SNAP)"
cvd.table[cvd.table=="fdas_enet"] <- "Food Assistance (SNAP)"
s.table[s.table=="fdas_enet"] <- "Food Assistance (SNAP)"
sens.60.table[sens.60.table=="fdas_enet"] <- "Food Assistance (SNAP)"

ac.table[ac.table=="hhs_enet"] <- "Household Size"
ca.table[ca.table=="hhs_enet"] <- "Household Size"
cvd.table[cvd.table=="hhs_enet"] <- "Household Size"
s.table[s.table=="hhs_enet"] <- "Household Size"
sens.60.table[sens.60.table=="hhs_enet"] <- "Household Size"

ac.table[ac.table=="pc1"] <- "Prudent #1"
ca.table[ca.table=="pc1"] <- "Prudent #1"
cvd.table[cvd.table=="pc1"] <- "Prudent #1"
s.table[s.table=="pc1"] <- "Prudent #1"
sens.60.table[sens.60.table=="pc1"] <- "Prudent #1"

ac.table[ac.table=="pc2"] <- "Prudent #2"
ca.table[ca.table=="pc2"] <- "Prudent #2"
cvd.table[cvd.table=="pc2"] <- "Prudent #2"
s.table[s.table=="pc2"] <- "Prudent #2"
sens.60.table[sens.60.table=="pc2"] <- "Prudent #2"


## Generate one table (main analysis) with all causes of death ##

all.table <- bind_rows( data.frame( index = "All-Cause Mortality"),
                        ac.table, 
                      data.frame( index = "Cancer-Specific Mortality"),
                      ca.table, 
                      data.frame( index = "Cardiovascular Disease Mortality"),
                      cvd.table )



## Save tables ##

write.table( ac.table, "04-Tables-Figures/tables/04-table-4-ac.txt", sep = "," )
write.table( ac.table, "04-Tables-Figures/tables/05-table-4-ca.txt", sep = "," )
write.table( ac.table, "04-Tables-Figures/tables/06-table-4-cvd.txt", sep = "," )
write.table( all.table, "04-Tables-Figures/tables/07-table-4-all.txt", sep = "," )
write.table( s.table, "04-Tables-Figures/tables/08-table-s3.txt", sep = "," )
write.table( sens.60.table, "04-Tables-Figures/tables/09-table-s4.txt", sep = "," )



### Survival Curves ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## FI pattern survival curves
( fi.sc <- ggadjustedcurves( fit = out.res[[1]]$q.obj,
                 variable ="fs_enet.q",
                 data = out.res[[1]]$dat,
                 method = "conditional",
                 title = "Food Insecurity (FI) Pattern",
                 font.title = c(16, "bold"),
                 legend.title = "Quintile",
                 font.legend = c(10, "bold"),
                 legend = c(0.14,0.34),
                 ylab = "Adjusted Survival Rate",
                 xlab = "Follow-up (Months)",
                 size = 0.6) +
  theme( text = element_text( family = "Avenir" ),
         axis.title.y = element_text( size = 13, margin = unit(c(0,6,6,0), 'pt') ) ,
         axis.text.y = element_text( size = 10, color = "grey30" ),
         axis.title.x = element_text( size = 13 ) ,
         axis.text.x = element_text( size = 10, color = "grey30" ) )+ 
  scale_color_ordinal() )

ggsave( "04-Tables-Figures/figures/02a-fi-surv-curve.png", 
        height = 7.21, 
        width = 6.42 )

# snap pattern survival curves
( snap.sc <- ggadjustedcurves( fit =  out.res[[4]]$q.obj,
                 variable = "fdas_enet.q",
                 data = out.res[[4]]$dat,
                 method = "conditional",
                 title = "Food Assistance (SNAP) Pattern",
                 font.title = c(16, "bold"),
                 legend.title = "Quintile",
                 font.legend = c(10, "bold"),
                 legend = c(0.14,0.31),
                 ylab = "",
                 xlab = "Follow-up (Months)",
                 size = 0.6) +
  theme(text=element_text(family="Avenir") ) + 
  scale_color_ordinal() + 
    theme( legend.position = "none",
           text = element_text( family = "Avenir" ),
           axis.title.y = element_text( size = 13 ) ,
           axis.text.y = element_text( size = 10, color = "grey30" ),
           axis.title.x = element_text( size = 13 ) ,
           axis.text.x = element_text( size = 10, color = "grey30" ) ) )

ggsave( "04-Tables-Figures/figures/02b-snap-surv-curve.png", 
        height = 7.21, 
        width = 6.42 )


## Arrange with Spline Plots (Only Food Insecurity and SNAP Patterns) ##

# FI pattern
fi.sp <- out.res[[1]]$spline.plot +
  xlab( unname( TeX( "FI Pattern Score$^a$" ) ) ) +
  theme( 
         text = element_text( family = "Avenir" ),
         axis.title.y = element_text( size = 13 ) ,
         axis.text.y = element_text( size = 10, color = "grey30" ),
         axis.title.x = element_text( size = 13 ) ,
         axis.text.x = element_text( size = 10, color = "grey30" ),
         legend.text = element_text( size = 10 ) ) 

# SNAP pattern
snap.sp <- out.res[[4]]$spline.plot +
  xlab( unname( TeX( "SNAP Pattern Score$^a$" ) ) ) +
  theme( legend.position = "none",
         text = element_text( family = "Avenir" ),
         axis.title.y = element_text( size = 13 ) ,
         axis.text.y = element_text( size = 10, color = "grey30" ),
         axis.title.x = element_text( size = 13 ) ,
         axis.text.x = element_text( size = 10, color = "grey30" ) ) +
  ylab( "" ) +
  coord_cartesian( ylim = c( 0.7, max = 2.8 ) ) 
  
ggarrange( ggarrange( fi.sc, fi.sp, nrow = 2, labels = list( "A", "B" ) ),
           ggarrange( snap.sc, snap.sp, nrow = 2,labels = list( "C", "D" ) ),
           nrow = 1, ncol = 2 )

## arrange into large ggarrange object

ggsave( "04-Tables-Figures/figures/03-surv-spline-comb.png",
        height = 7.4, 
        width = 8.96 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

