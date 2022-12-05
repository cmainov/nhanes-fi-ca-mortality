library( tidyverse ) 
library( survey )    # for complex survy data model-fitting



( d <- readRDS( "02-Data-Wrangled/01-covariate-mortality-linkage.rds" ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ) ) 

# n.total unique
# 1   86538  78964
# NOTE: there are some duplicates associated with the merging of the mortality data


### Exclusions Based on Linked Mortality Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## Remove those ineligible for survival analysis due to : ##

# i. < 18 yrs of age
# ii. insufficient identifying data, per the NCHS ( see: https://www.cdc.gov/nchs/data/datalinkage/public-use-linked-mortality-files-data-dictionary.pdf)
# we will also filter by our age exclusion parameter at this point as well: age < 20

( d.1 <- d %>%
  filter( eligstat == 1 | age >= 20 ) ) %>%
  summarise( n.total = n(), 
             unique = length( unique( .$seqn ) ),
             removed = nrow( d %>% distinct() ) - unique,
             duplicated = n.total - unique ) 

# n.total unique removed
# 1   78250  73180   10732
# NOTE: some duplicates remain

# diagnose the duplicate issue:

# look at columns of interest first
d.1 %>%
  filter( duplicated( seqn ) & !is.na( mortstat ) ) %>%
  select( seqn, mortstat, hyperten, ucod_leading, age, diabetes ) 
# looks like the duplicates all have missing data for the columns we will be using
# except for mortstat

# double check they all have missing day in the pertinent columns
sapply( d.1 %>%
          filter( duplicated( seqn ) ) %>%
          select(seqn, mortstat, hyperten, ucod_leading, age, diabetes ),
        function( x ) sum( is.na( x ) ) )

# decision: we will group by seqn and keep the rows with the "max" value among
# NA, 0, or 1 ( meaning that we should keep the rows that have a value)
( d.2 <- d.1 %>% 
    filter( duplicated( seqn ) ) %>%
    group_by( seqn ) %>%
    filter( mortstat == max( mortstat, na.rm = T ) ) %>%
    ungroup() %>%
  bind_rows( d.1 %>%
               filter( ! duplicated( seqn ) ), . ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             duplicated = n.total - unique )

# n.total unique
# 1   75170  73180  # 1990 still duplicated

# recheck missings
sapply( d.2 %>%
          filter( duplicated( seqn ) ) %>%
          select(seqn, mortstat, hyperten, ucod_leading, age, diabetes ),
        function( x ) sum( is.na( x ) ) )
# still have 51 duplicates without missing data in the mortality columns

# remove those with missing cause of death
( d.3 <- d.2 %>%
  filter( !is.na( ucod_leading ) ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             removed = nrow( d.2 %>% distinct() ) - unique,
             duplicated = n.total - unique )
# n.total unique removed duplicated
# 1    6743   6725   68350         18

# investigate further
d.3 %>%
  filter( duplicated( seqn ) ) %>%
  select( seqn, eligstat, mortstat, ucod_leading, permth_int,
          permth_exm, stime ) %>%
    arrange( seqn )

# it looks like the remaining duplicate rows have a row with shorter follow-up
# time and one with longer follow-up time. We will retain the row with the longer
# follow up time


( d.4 <- d.3 %>%
  group_by( seqn ) %>%
  filter( permth_int == max( permth_int, na.rm = T ) ) %>% 
  ungroup() ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
             removed = nrow( d.3 %>% distinct() ) - unique,
             duplicated = n.total - unique )

# no duplicates remain
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


### Exclusions Based on Covariates ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------




mod.id$model %>% filter( duplicated( seqn ))


( d.4 %>% 
rowwise( ) %>%
  mutate( first.age = suppressWarnings( min( agedxa, agedxb, agedxc, na.rm = T ) ),
          first.age = ifelse( is.infinite( first.age ), NA, first.age ) ) %>%
  ungroup( ) %>%
  filter( is.na( first.age ) == F ) ) %>%
  summarise( n.total = n(), unique = length( unique( .$seqn ) ),
           removed = nrow( d.3 %>% distinct() ) - unique,
           duplicated = n.total - unique )


d.5 <- d.4 %>%
mutate( stime = permth_int + timesincecadxmn,
        
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
  filter( ca == 1) %>%
  mutate( fs.enet.q = as.factor( quant_cut( var = "hhs_enet", x = 5, df = . ) ) ) %>%
  select( seqn,fs.enet.q, castat, cvdstat, dmstat, stime) %>%
  left_join(d.4, .) %>% filter(is.na( wtdr18yr))




d.try <- d.1 %>%
  filter( !is.na( wtmec18yr))

mod1<-svydesign(id = ~sdmvpsu, weights = ~wtmec18yr, strata = ~sdmvstra, 
                nest = TRUE, survey.lonely.psu = "adjust", data = d.5)
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
                nest = TRUE, survey.lonely.psu = "adjust", data = d.5)
mod2<-subset(mod2,!is.na(age_enet) & ca == 1)#inclusions


m.2 <- svycoxph(formula( paste0( "Surv(stime,mortstat) ~ fs.enet.q +cycle +", paste0( covars.logit, collapse = " + ") ) ),
                design=mod2)

