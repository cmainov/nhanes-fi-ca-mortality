###------------------------------------------------------------
###   01-INCLUSIONS AND EXCLUSIONS
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we run through the exclusion and inclusion criteria for the 
# analysis that is then reflected in Figure 1 of the manuscript. We generate
# a new dataset that has an indicator variable that indicates membership in the
# final analytic sample.
#
# INPUT DATA FILES: 
# i."02-Data-Wrangled/01-covariate-mortality-linkage.rds"
#
# OUPUT FILES:
# i."02-Data-Wrangled/02-inclusions-exclusions.rds"
#
# Resources:
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( tidyverse ) 
library( survey )    # for complex survey data model-fitting

source( "R/utils.R" )


### (0.0) Read in Working Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

( d <- readRDS( "02-Data-Wrangled/01-covariate-mortality-linkage.rds" ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ) ) 

# n.total unique
# 1   78964  78964

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) Step 1 ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (1.0) Age >= 20 and response of "1" ( "yes" ) to the MCQ220 probe on cancer history ##
( ex.1 <- as.numeric( ( step1.data <- d %>%
                          filter( age >= 20 & ca == 1 ) ) %>%
                        summarise( subjects.remaining = n( ) ) ) ) # 4079 subjects remain

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) Step 2 ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Exclude those without complete dietary data/missing weights from the 24-hr recall subsample ##

( step2a.data <- step1.data %>%
    filter( is.na( wtdr18yr ) == F & wtdr18yr != 0 ) ) %>%
  summarise( subjects.remaining = n( ) , subjects.excluded = ex.1 - n( ) )

# subjects.remaining subjects.excluded
#               3959              120

ex.2a <- nrow( step2a.data )

## ---o--- ##


## (2.2) Exclude those with a history of non-melanoma skin cancer as their only cancer diagnosis or ##
## outright missing cancer diagnosis data ##

( step2b.data <- step2a.data %>%
    filter( !( primaryca == "Non-Melanoma Skin" ) ) ) %>%
  summarise( subjects.remaining = n( ) , subjects.excluded = ex.2a - n( ) )

# subjects.remaining subjects.excluded
#               3383               576

ex.2b <- nrow( step2b.data )

## ---o--- ##


## (2.3) Missing age at diagnosis of first cancer ##

( step2c.data <- step2b.data %>% 
    rowwise( ) %>%
    mutate( first.age = suppressWarnings( min( agedxa, agedxb, agedxc, na.rm = T ) ),
            first.age = ifelse( is.infinite( first.age ), NA, first.age ) ) %>%
    ungroup( ) %>%
    filter( is.na( first.age ) == F ) ) %>%
  data.frame( ) %>%
  summarise( subjects.remaining = n( ), subjects.excluded = ex.2b - n( ) )

# subjects.remaining subjects.excluded
#               3379                 4

ex.2c <- nrow( step2c.data )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Step 3 ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Missing covariate and survival data


( step3.data <- step2c.data %>%
    filter( !is.na( foodasstpnowic ) & 
              !is.na( binfoodsechh ) &
              !is.na( hhsize ) &
              !is.na( agecat ) &
              !is.na( race ) &
              !is.na( ins.status ) &
              !is.na( age ) &
              !is.na( bmxbmi ) &
              !is.na( smokstat ) &
              !is.na( fipr ) &
              !is.na( kcal ) &
              !is.na( weekmetmin ) &
              !is.na( education_bin ) &
              !is.na( cci_score ) &
              !is.na( alc_cat ) &
              !is.na( stime ) &
              !is.na( mortstat ) ) %>%
    mutate( inc = 1 ) ) %>% # final sample indicator
  summarise( subjects.remaining = n( ) , subjects.excluded = ex.2c - n( ) )



# subjects.remaining subjects.excluded
#               2493                884
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### (4.0) Bind Indicated Data Back and Save ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

( d.1 <- d %>%
  filter( seqn %notin% step3.data$seqn ) %>%
  bind_rows( ., step3.data ) )%>%
  summarise( final.n = n() )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### (5.0) Save ###  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( d.1, "02-Data-Wrangled/02-inclusions-exclusions.rds" )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


