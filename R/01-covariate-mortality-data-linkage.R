###------------------------------------------------------------
###   01-DATA MERGE AND WRANGLING
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will import the NHANES survey study data and merge it with
# the linked mortality data from the CDC/NCHS linkage file. We use the provided
# R code file for running the merge. 
#
# INPUT DATA FILES: 
# NHANES data from previous analysis stored in a GitHub repository linked in code below.
#
# OUPUT FILES:
# i."02-Data-Wrangled/01-covariate-mortality-linkage.rds"
#
# Resources:
#   i. NCHS NHIS/NHANES Mortality Data Linkage Codebook: https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf
#   ii. `hei` source code: https://github.com/timfolsom/hei
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( survey )     # survey commands
library( tidyverse )
library( RNHANES )    # for downloading NHANES data
library( haven )      # for reading/saving alternative formats of data
library( hei ) # we modify the `hei` code to support the column names in our data

source( "R/utils.R" ) # helper functions


## (0.0) Read in Data from Previous Analysis ##
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

url.raw <- "https://github.com/cmainov/NHANES-Diet-Penalized-Regression/blob/main/03-Data-Rodeo/04-Analytic-Data.rds?raw=true"

d <- readRDS( url( url.raw ) ) %>%
  rename_all( tolower ) # lower case all column names

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) Read in Mortality Files ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


## (1.1) Download raw files (Accessed 14 November 2022) ##

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

## ---o--- ##


## (1.2) Use NCHS-provided code to piece together data ##

dir.files <- paste0( "01-Data-Raw/", dir( "01-Data-Raw" ) )

datout <- data.frame() #initialize dataframe

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



### (2.0) Merge, Wrangle, and Save ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Clean up mortality variables ##

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
          dmstat = ifelse( ucod_leading == 7, 1, 
                           ifelse( !is.na( mortstat ) & !is.na( ucod_leading ) & ucod_leading != 7, 0, NA ) ),
          
          # death from influenza/pneumonia indicator
          pnstat = ifelse( ucod_leading == 8, 1, 
                           ifelse( !is.na( mortstat ) & !is.na( ucod_leading ) & ucod_leading != 8, 0, NA ) ))%>%
  
  # clean up those not recognized in above script
  mutate( castat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                             is.na( castat ), 0, castat ),
          
          cvdstat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                              is.na( cvdstat ), 0, cvdstat ),
          
          dmstat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                             is.na( dmstat ), 0, dmstat ),
          
          pnstat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                             is.na( pnstat ), 0, pnstat ) )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Fix Issues with Duplicates in Mortality-Linked Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Separate those with mortstat == 1 and mortstat == 0 and mortstat == NA ##

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

## ---o--- ##


## (3.2) Mortstat == 0 first ##


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

## ---o--- ##


## (3.3) Mortstat == 1 ## 

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

## ---o--- ##


## (3.4) Mortstat == NA ##

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
d.na.1.dups <- d.na.1[ duplicated( d.na.1$seqn ), "seqn" ]

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

## ---o--- ##


## (3.5) Bind data with non-duplicates back together ##

( d.1 <- bind_rows( d.na.2, d.01.1, d.00.3 ) %>% data.frame() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# we still have some duplicates meaning that there must be individuals with different mortstat indicators

d.1.dups <- d.1[ duplicated( d.1$seqn ), "seqn"]

( d.1.1 <- d.1 %>%
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



### (4.0) Health Insurance Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.1) Set-up data query details ##

yrs.nh <- c( "1999-2000", "2001-2002", "2003-2004", "2005-2006", "2007-2008",
             "2009-2010", "2011-2012", "2013-2014" )

these.i <- seq_along(yrs.nh)

yrs.short <- c( 99,1,3,5,7,9,11,13 ) 
yrs.short<- ifelse( yrs.short < 10 , paste0( "0", yrs.short ), yrs.short )

hiq.surveys <- c( "HIQ", paste0( "HIQ_", LETTERS[2:8] ) )

## ---o--- ##


## (4.2) Use RHANES package to import data ##

l.hiq <- lapply( these.i,
                  
                  function(x){
                    
                    nhanes_load_data( hiq.surveys[x], yrs.nh[x], demographics = FALSE ) 
                    
                  }
)

l.hiq[[9]] <- read_xpt( file = "https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/HIQ_I.XPT" ) # manual 2015 data
l.hiq[[10]] <- read_xpt( file = "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/HIQ_J.XPT" ) # manual 2017 data

# do `bind_rows` so we can see which variables are not consistent across cycles
d.hiq <- do.call( "bind_rows", l.hiq ) %>%
  rename( seqn = SEQN ) %>%
  mutate( ins.status = ifelse( cycle %in% c( "1999-2000", "2001-2002", "2003-2004" ) & HID010 == 1, "Yes",
                  ifelse( cycle %in% c( "1999-2000", "2001-2002", "2003-2004" ) & HID010 == 2, "No",
                          ifelse( cycle %in% c( "1999-2000", "2001-2002", "2003-2004" ) & HID010 %in% c( 7, 9 ), NA, 
                                  ifelse( cycle %notin% c( "1999-2000", "2001-2002", "2003-2004" ) & HIQ011 == 1, "Yes",
                                          ifelse( cycle %notin% c( "1999-2000", "2001-2002", "2003-2004" ) & HIQ011 == 2, "No",
                                                  ifelse( cycle %notin% c( "1999-2000", "2001-2002", "2003-2004" ) & HIQ011 %in% c( 7, 9 ), NA,
                                                          NA )))))))


## ---o--- ##


## (4.3) Merge with working data ##

( d.3 <- d.2 %>%
    left_join(., d.hiq %>%
                select( seqn, ins.status ) ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### (5.0) Physical Functioning Survey (ADL Scale) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) Set-up data query details ##

pfq.surveys <- c( "PFQ", paste0( "PFQ_", LETTERS[2:8] ) )

## ---o--- ##


## (5.2) Use RHANES package to import data ##

l.pfq <- lapply( these.i,
                 
                 function(x){
                   
                   nhanes_load_data( pfq.surveys[x], yrs.nh[x], demographics = FALSE ) 
                   
                 }
)

## ---o--- ##


## (5.3) 2015-2018 data via manual download ##

# mot currently available in RHANES package functions

l.pfq[[9]] <- read_xpt( file = "https://wwwn.cdc.gov/Nchs/Nhanes/2015-2016/PFQ_I.XPT" ) # manual 2015 data
l.pfq[[10]] <- read_xpt( file = "https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/PFQ_J.XPT" ) # manual 2017 data

## ---o--- ##


## (5.4) Clean up imported data and merge to working data ##

# make 1999 names upper case bc they are not
colnames( l.pfq[[1]] ) <- toupper( colnames( l.pfq[[1]]  ) )

# columns to extract (column names of 19-items in the ADL scale)
nms.adl.99 <- paste0( "PFQ060", LETTERS[1:19] ) #1999-2002
nms.adl.05 <- paste0( "PFQ061", LETTERS[1:19] ) #2003-

adl.out <- data.frame() #initialize empty data.frame for combining the results across cycles
for( i in 1:length( l.pfq ) ){
  
  if ( i %in% 1:2 ) {
    
    
    these.99 <- setNames( l.pfq[[i]][c( "SEQN", nms.adl.99 )],
                          c( "SEQN", nms.adl.05 ) )
    adl.out <- rbind( adl.out, these.99 )
    
  }
  
  if ( i > 3 ) {
    
    adl.out <- rbind( adl.out, l.pfq[[i]][c( "SEQN", nms.adl.05 )])
  }
                         
}

# set to missing values other than 1:4 for these columns
adl.out[, 2:ncol(adl.out) ] <- sapply( 2:ncol(adl.out),
       function(x){
         adl.out[,x]<-ifelse( adl.out[,x] %notin% c(1:4), NA, adl.out[,x] )
       }
       ) %>%
  data.frame()


# merge
d.4 <- adl.out %>%
  rename( seqn = SEQN ) %>%
  mutate( adl.score = PFQ061A + PFQ061B + PFQ061C + PFQ061D + 
            PFQ061E + PFQ061F + PFQ061G + PFQ061H + PFQ061I + 
            PFQ061J + PFQ061K + PFQ061L + PFQ061M + PFQ061N + 
            PFQ061O + PFQ061P + PFQ061Q + PFQ061R + PFQ061S ) %>%
  select( seqn, adl.score ) %>%
  left_join( d.3, . )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (6.0) HEI-2015 Scores ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (6.1) Clean up column names ##

d.5 <- d.4 %>% # rename columns for use in `hei` function
  rename( TKCAL = kcal,
          TSFAT = sfat,
          TMFAT = mfat,
          TPFAT = pfat,
          TSODI = sodium,
          TALCO = alco,
          T_PF_SEAFD_LOW = fish_lo,
          T_F_CITMLB = f_citmelber,
          T_PF_EGGS = eggs,
          T_F_CITMLB = f_citmelber,
          T_F_OTHER = fruitother,
          T_PF_MPS_TOTAL = mpstotal,
          T_PF_SOY = soy,
          T_PF_SEAFD_HI = fish_hi,
          T_PF_NUTSDS = nuts,
          T_SOLID_FATS = solidfats,
          T_ADD_SUGARS = addedsugars,
          T_V_TOTAL = vegtotal,
          T_V_LEGUMES = legumes,
          T_F_TOTAL = fruittotal,
          T_D_TOTAL = dtotal,
          T_V_DRKGR = greenleafy,
          T_G_WHOLE = wholegrain,
          T_G_REFINED = refinedgrain )

# fped data object
fped.hei <- d.5 %>%
  select( SEQN = seqn,
          RIDAGEYR = age,
          T_PF_SEAFD_LOW ,
          T_F_CITMLB,
          T_PF_EGGS,
          T_F_CITMLB,
          T_F_OTHER,
          T_PF_MPS_TOTAL,
          T_PF_SOY,
          T_PF_SEAFD_HI,
          T_PF_NUTSDS,
          T_SOLID_FATS,
          T_ADD_SUGARS,
          T_V_TOTAL,
          T_V_LEGUMES,
          T_F_TOTAL,
          T_D_TOTAL,
          T_V_DRKGR,
          T_G_WHOLE,
          T_G_REFINED )

fped.hei$DRSTZ <- 1

# diet data object
diet.hei <- d.5 %>%
  select( SEQN = seqn,
          TKCAL,
          TSFAT,
          TMFAT,
          TPFAT,
          TSODI,
          TALCO )

# demo data object
dem.hei <- d.5 %>%
  select( SEQN = seqn,
          SDDSRVYR = cycle )

## ---o--- ##


## (6.2) Use `hei` to compute hei-2015 score and merge ##

d.6 <- hei( fped = fped.hei, diet = diet.hei, demograph = dem.hei ) %>%
  select( seqn = SEQN,
          hei.2015 = HEI ) %>%
  left_join( d.4, . ) %>%
  
  # remove elastic net dietary patterns from previous analysis since they will be recomputed in this analysis
  select( -c( fs_enet, fdas_enet, age_enet, hhs_enet,
              pc1, pc2 ) )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (7.0) Save ###  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( d.6, "02-Data-Wrangled/01-covariate-mortality-linkage.rds")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


