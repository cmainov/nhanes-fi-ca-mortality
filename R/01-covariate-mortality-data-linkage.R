

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


left_join( d, datout %>% mutate( seqn = as.numeric( seqn ) ), by = "seqn" ) %>%
  
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
                             is.na( dmstat ), 0, dmstat ) ) %>%
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### Save ###  
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( "02-Data-Wrangled/01-covariate-mortality-linkage.rds")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



d.1 <- left_join( d, datout %>% mutate( seqn = as.numeric( seqn ) ), 
                  by = "seqn" ) %>%
  
  # create survival time variable as time since diagnosis
  mutate( stime = permth_int + timesincecadxmn,
          
          # generate cancer mortality indicator for cancer-mortality analysis
          castat = ifelse( ucod_leading == 2, 1, 
                           ifelse( !is.na( mortstat ) & !is.na( ucod_leading ) & ucod_leading != 2, 0, NA ) ) ) %>%
  
  # address remaining un labeled subjects
  mutate( castat = ifelse( !is.na( mortstat ) & is.na( ucod_leading ) &
                             is.na( castat ), 0, castat ) ) 

###ID model for later use in GGplot construction--see end of program
modelfullaid<-svycoxph(Surv(stime,mortstat)~factor(FS_ENet_q)+Race+Gender+Age+BMXBMI+SmokStat+fipr+
                         KCAL+WeekMetMin+Education2+PrimaryCAGroup+CCI_Score+SEQN,design=mod1)
#bind rownumber and ID #
idmod<-data.frame(id=modelfullaid$model[,14],rownm=rownames(modelfullaid$model))

mortdat$include<-ifelse(mortdat$SEQN %in% idmod$id,1,0)


d.1 <- left_join( d, datout %>% mutate( seqn = as.numeric( seqn ) ), by = "seqn" ) %>%

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


d.try <- d.1 %>%
  filter( !is.na( wtmec18yr))

mod1<-svydesign(id = ~sdmvpsu, weights = ~wtmec18yr, strata = ~sdmvstra, 
                nest = TRUE, survey.lonely.psu = "adjust", data = d.try)
mod1<-subset(mod1,ca==1)#inclusions



covars.logit <- tolower( c( 'race', 'Gender', 'Age', 'BMXBMI', 'HHSize',
                   'SmokStat', 'fipr', 'KCAL', 'WeekMetMin', 'Education_bin',
                   'CCI_Score', "alc_cat" ) )



d.1$fs_enet <- as.vector(d.1$fs_enet)
jj<-as.vector(d.1$fs_enet)

mod.id<-svycoxph(formula( paste0( "Surv(stime,mortstat) ~ binfoodsechh + seqn +", paste0( covars.logit, collapse = " + ") ) ),
         design=mod1)

idmod<-data.frame(id=mod.id$model[,3],rownm=rownames(mod.id$model))

d.2 <- d.1%>%
  mutate( include = ifelse( seqn %in% idmod$id, 1, 0 ) ) %>%
  filter( include == 1 & !is.na( fs_enet) & !is.na( wtdr18yr)) %>%
  mutate( fs.enet.q = as.factor( quant_cut( var = "hhs_enet", x = 5, df = . ) ) ) %>%
  select(fs.enet.q, seqn, include) %>%
  left_join(d.1, ., by = "seqn") %>%
  filter(!is.na(wtdr18yr))



mod2<-svydesign(id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                nest = TRUE, survey.lonely.psu = "adjust", data = d.2)
mod2<-subset(mod2,include == 1 & !is.na(age_enet) & !is.na( wtdr18yr) &
               cycle %in% 1:10)#inclusions


m.2 <- svycoxph(formula( paste0( "Surv(stime,mortstat) ~ fs.enet.q +cycle +", paste0( covars.logit, collapse = " + ") ) ),
         design=mod2)

