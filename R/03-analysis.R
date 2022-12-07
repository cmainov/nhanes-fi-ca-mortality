library( tidyverse )
library(glue)
library(splines)

d <- readRDS( "03-Data-Rodeo/01-analytic-data.rds")

d.1 <- d %>%
  filter( inc == 1) %>%
  mutate( fs.enet.q = as.factor( quant_cut( var = "fs_enet", x = 5, df = . ) ) ) %>%
  bind_rows( d %>% filter( is.na(inc) ), . )

# covariates to include in model
covars.surv <- tolower( c( 'race', 'Gender', 'Age', 'BMXBMI', 'HHSize',
                            'SmokStat', 'fipr', 'KCAL', 'WeekMetMin', 'Education_bin',
                            'CCI_Score', "alc_cat" ) )

indices <- c( "fs_enet", "age_enet", "hhs_enet", "fdas_enet", "pc1", "pc2" )

# all-cause mortality
out.res <- list()
for( i in 1:length( indices ) ){
  
  out.res[[i]] <- res( df = d, x = indices[i], 
     subs = "inc == 1", 
     cuts = 5, 
     id.col = "seqn", 
     covars = covars.surv, 
     time = "stime", 
     mort.ind = "mortstat" ) 
}

fin.res <- do.call( "rbind", out.res)
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

fin.res.ca <- do.call( "rbind", out.res.ca )
View(fin.res.ca)

