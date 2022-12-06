

# Resources:
#   i. NCHS NHIS/NHANES Mortality Data Linkage Codebook: https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf


library( survey )
library(tidyverse)

# read in helper functions
source( "R/utils.R")


## Read in data from previous analysis ##
url.raw <- "https://github.com/cmainov/NHANES-Diet-Penalized-Regression/blob/main/03-Data-Rodeo/04-Analytic-Data.rds?raw=true"

d <- readRDS( url( url.raw ) ) %>%
  rename_all( tolower ) # lower case all column names



### Read in Mortality Files ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## Download raw files (Accessed 14 November 2022) ##

yrs <- c( "1999_2000", "2001_2002", "2003_2004", "2005_2006", "2007_2008",
          "2009_2010", "2011_2012", "2013_2014", "2015_2016", "2017_2018")

for( i in 1:length( yrs ) ) {
  
  # create URL for each cycle
  url.in <- paste0( "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_", 
                    yrs[i], 
        "_MORT_2019_PUBLIC.dat" )
  
  # download file and save in directory
  download.file( url.in, 
              method = "curl", 
              destfile = paste0( "01-Data-Raw/mort-", 
                                 str_replace( yrs[i],"\\_", "-"),
                                 ".dat" ) )

}



## Use provided code to piece together data ##

dir.files <- paste0( "01-Data-Raw/", dir( "01-Data-Raw" ) )

datout<-data.frame() #initialize dataframe
for( i in 1:length( dir.files ) ){
  
  # read in the fixed-width format ASCII file
  dsn <- read_fwf(file = dir.files[i],
                  col_types = "ciiiiiiiddii",
                  fwf_cols(publicid = c(1,14),
                           eligstat = c(15,15),
                           mortstat = c(16,16),
                           ucod_leading = c(17,19),
                           diabetes = c(20,20),
                           hyperten = c(21,21),
                           dodqtr = c(22,22),
                           dodyear = c(23,26),
                           wgt_new = c(27,34),
                           sa_wgt_new = c(35,42),
                           permth_int = c(43,45),
                           permth_exm = c(46,48)
                  ),
                  na = "."
  )
  
  # create the ID (SEQN) for the NHANES surveys
  dsn$seqn <- substr(dsn$publicid,1,5)
  # NOTE:   SEQN is the unique ID for NHANES.
  
  #Drop NHIS variables
  dsn <- dsn%>%select(-publicid)
  dsn <- select(dsn, -dodqtr)
  dsn <- select(dsn, -dodyear)
  dsn <- select(dsn, -wgt_new)
  dsn <- select(dsn, -sa_wgt_new)
  
  
  # Structure and contents of data
  str(dsn)
  
  
  # Variable frequencies
  
  #ELIGSTAT: Eligibility Status for Mortality Follow-up
  table(dsn$eligstat)
  #1 = "Eligible"
  #2 = "Under age 18, not available for public release"
  #3 = "Ineligible"
  
  #MORTSTAT: Final Mortality Status
  table(dsn$mortstat, useNA="ifany")
  # 0 = Assumed alive
  # 1 = Assumed deceased
  # <NA> = Ineligible or under age 18
  
  #UCOD_LEADING: Underlying Cause of Death: Recode
  table(dsn$ucod_leading, useNA="ifany")
  # 1 = Diseases of heart (I00-I09, I11, I13, I20-I51)
  # 2 = Malignant neoplasms (C00-C97)
  # 3 = Chronic lower respiratory diseases (J40-J47)
  # 4 = Accidents (unintentional injuries) (V01-X59, Y85-Y86)
  # 5 = Cerebrovascular diseases (I60-I69)
  # 6 = Alzheimer's disease (G30)
  # 7 = Diabetes mellitus (E10-E14)
  # 8 = Influenza and pneumonia (J09-J18)
  # 9 = Nephritis, nephrotic syndrome and nephrosis (N00-N07, N17-N19, N25-N27)
  # 10 = All other causes (residual)
  # <NA> = Ineligible, under age 18, assumed alive, or no cause of death data
  
  #DIABETES: Diabetes Flag from Multiple Cause of Death (MCOD)
  table( dsn$diabetes, useNA = "ifany" )
  # 0 = No - Condition not listed as a multiple cause of death
  # 1 = Yes - Condition listed as a multiple cause of death
  # <NA> = Assumed alive, under age 18, ineligible for mortality follow-up, or MCOD not available
  
  #HYPERTEN: Hypertension Flag from Multiple Cause of Death (MCOD)
  table( dsn$hyperten, useNA = "ifany" )
  # 0 = No - Condition not listed as a multiple cause of death
  # 1 = Yes - Condition listed as a multiple cause of death
  # <NA> = Assumed alive, under age 18, ineligible for mortality follow-up, or MCOD not available
  
  # Re-name the dataset, DSN, to the short survey name then remove other R objects
  datout <- rbind( datout, dsn )
}

# ---------------------------------------------------------------------------------------------------------------------------------------------------------





### Merge, Wrangle, and Save ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


dat <- left_join( d, datout %>% mutate( seqn = as.numeric( seqn ) ), by = "seqn" ) %>%
  
  # survival time computed as time since diagnosis
  mutate( stime = permth_int + timesincecadxmn,
          
          ## indicators for cause-specific mortality analyses ##
          # death from malignant neoplasm indicator
          castat = ifelse( ucod_leading == 2, 1, 
                           ifelse( !is.na( mortstat ) & !is.na( ucod_leading ) & ucod_leading != 2, 0, NA ) ),
          
          # death from cardiovascular disease indicator (heart and cerebrovascular dz)
          cvdstat = ifelse( ucod_leading %in% c( 1, 5 ), 1, 
                            ifelse( !is.na( mortstat ) & !is.na( ucod_leading ) & ucod_leading %notin% c( 1, 5 ), 0, NA ) ),
          
          # death from diabetes mellitus indicator
          dmstat = ifelse( ucod_leading == 2, 1, 
                           ifelse( !is.na( mortstat ) & !is.na( ucod_leading ) & ucod_leading != 2, 0, NA ) )) %>%
  
  # clean up those not recognized in above script
  mutate( castat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                             is.na( castat ), 0, castat ),
          
          cvdstat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                              is.na( cvdstat ), 0, cvdstat ),
          
          dmstat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                             is.na( dmstat ), 0, dmstat ) )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### Fix Issues with Duplicates in Mortality-Linked Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# separate those with mortstat == 1 and mortstat == 0 and mortstat == NA

( d.01 <- dat %>%
    filter( mortstat == 1 ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )


( d.00 <- dat %>%
    filter( mortstat == 0 ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

( d.na <- dat %>%
    filter( is.na( mortstat) ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )


## Mortstat == 0 first ##


# look at columns of interest first
d.00 %>%
  filter( duplicated( seqn ) ) %>%
  select( seqn, mortstat, permth_int, hyperten, ucod_leading, age, diabetes ) 

# looks like issue is with follow up time: there are several follow up times provided for some individuals
# we will take the largest follow up time from the set to rectify this


( d.00.1 <- d.00 %>%
    group_by( seqn ) %>%
    filter( permth_int == max( permth_int, na.rm = T ) ) %>%
    ungroup() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )
# 67 duplicates remain

d.00.1 %>%
  filter( duplicated( seqn ) ) %>%
  select( seqn, mortstat, permth_int, hyperten, ucod_leading, age, diabetes, stime, castat, cvdstat, dmstat, eligstat )

# it looks like the remaining duplicates have same value for permth_int in their duplications. Thus, we should be able to use 
# distinc() to solve the issue

( d.00.2 <- d.00.1 %>%
    distinct() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )
# 35 remain duplicated

d.00.2 %>%
  filter( duplicated( seqn ) ) 

# remaining issue appears to be in the permth_exm column which we are not using, so we should be able to remove this column and 
# then use distinct() again

( d.00.3 <- d.00.2 %>%
    select( -permth_exm ) %>%
    distinct() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# resolved


## Mortstat == 1 ## 

# look at columns of interest first
d.01 %>%
  filter( duplicated( seqn ) ) %>%
  select( seqn, mortstat, permth_int, hyperten, ucod_leading, age, diabetes, permth_exm ) 

# looks like it is similar issues with permth_int and permth_exm columns. we will use similar approach as above

( d.01.1 <- d.01 %>%
    group_by( seqn ) %>%
    select( -permth_exm ) %>%
    filter( permth_int == max( permth_int, na.rm = T) ) %>%
    ungroup() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# duplicates resolved


## Mortstat == NA ##

# look at columns of interest first
d.na %>%
  filter( duplicated( seqn ) ) %>%
  select( seqn, mortstat, permth_int, hyperten, ucod_leading, age, diabetes, permth_exm ) 

# doesn't appear to be any differences within duplicates. will try using distinct()

( d.na.1 <- d.na %>%
    distinct() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# 36 remain duplicated. let's look at those columns again
d.na.1.dups <- d.na.1[ duplicated( d.na.1$seqn ), "seqn"]

d.na.1 %>%
  filter( seqn %in% d.na.1.dups ) %>% arrange( seqn ) %>%
  select( seqn, mortstat, permth_int, hyperten, ucod_leading, age, diabetes, permth_exm, eligstat ) 
# issue is with the eligstat column. some individuals have both 2 and 3 (both of which indicate ineligible)
# we will keep only one

( d.na.2 <- d.na.1 %>%
    filter( seqn %in% d.na.1.dups ) %>%
    group_by( seqn ) %>%
    filter( eligstat == min( eligstat ) ) %>%
    ungroup() %>%
    rbind( ., d.na.1 %>% filter( seqn %notin% d.na.1.dups ) ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )
# resolved


## Bind data with non-duplicates back together ##

( d.1 <- bind_rows( d.na.2, d.01.1, d.00.3 ) %>% data.frame() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# we still have some duplicates meaning that there must be individuals with different mortstat indicators

d.1.dups <- d.1[ duplicated( d.1$seqn ), "seqn"]

(d.1.1 <- d.1 %>%
    filter( seqn %in% d.1.dups ) %>%
    group_by( seqn )%>%
    filter( mortstat == max( mortstat, na.rm = T ) ) %>%
    ungroup()) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# bind non-duplicates with rectified duplicates data frame
( d.2 <- ( d.1.2 <- d.1 %>%
            filter( seqn %notin% d.1.dups ) ) %>%
    rbind( ., d.1.1 ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )
# duplicates have been rectified and d.2 is final
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### Save ###  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( d.2, "02-Data-Wrangled/01-covariate-mortality-linkage.rds")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


