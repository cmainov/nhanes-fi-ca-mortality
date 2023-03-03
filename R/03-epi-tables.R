library( survey )
library( tidyverse )
library( jtools ) # for svycor function

source( "R/utils.R" )



dat <- readRDS( "03-Data-Rodeo/01-analytic-data.rds") %>%
  dplyr::filter( is.na( wtdr18yr ) == F ) %>% # subset to those having non -missing weights 
  
  dplyr::mutate( 
    hhsize.bin = ifelse( hhsize >= 5, 1,
                                        ifelse( hhsize < 5, 0, NA ) ),
    # concatenate cause of death into a single variable for table
    cod = ifelse( castat == 1, "Cancer",
                  ifelse( cvdstat == 1, "Cardiovascular Disease",
                          ifelse( mortstat == 1 & castat != 1 & cvdstat != 1, "Other",
                                  NA ) ) ),
    timesince.cat.5yr = ifelse( timesincecadxmn <= 60, "<= 5 yrs", ">5 yrs"))
  
  # designs
  nhc <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra,
                    nest = TRUE, survey.lonely.psu = "adjust", data = dat )
  
  # create three subsets, food insecure with cancer ( fiw ), food secure with cancer ( fsw ), and combined 
  # food secure and insecure cancer survivors
  fiw <- subset( nhc, binfoodsechh == "Low" & inc == 1 )
  fsw <- subset( nhc, binfoodsechh == "High" & inc == 1 )
  gen <- subset( nhc, inc == 1 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

### Person-Years and -Months Calculations ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
  
( pm <- sum( gen$variables$stime ) ) # person-months

pm/12 # person-years 

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
  
  
### Table 1 ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# write function to return table 1 with specific variables

cafs_table1 <- function( design, df ){
  
  sex <- epitab( var = "gender", data.fr = df, des = design, table.var = "Gender" )
  alcuse <- epitab( var = "alc_cat", data.fr = df, des = design, table.var = "Alcohol Use" )
  race <- epitab( var = "race", data.fr = df, des = design, table.var = "Race/Ethnicity" )
  smokstat <- epitab( var = "smokstat", data.fr = df, des = design, table.var = "Smoking Status" )
  income <- epitab( var = "fipr", data.fr = df, des = design, table.var = "Income:Poverty" )
  hhsize <- epitab( var = "hhsize.bin", data.fr = df, des = design, table.var = "Household Size" )
  education <- epitab( var = "education_bin", data.fr = df, des = design, table.var = "Education Attained" )
  bmi <- epitab.means( cont.var = "bmxbmi", des = design, table.var = "BMI", dig = 2 )
  metmins <- epitab.means( cont.var = "weekmetmin", des = design, table.var = "MET Minutes", dig = 2 )
  calor <- epitab.means( cont.var = "kcal", des = design, table.var = "Calories", dig = 2 )
  age <- epitab.means( cont.var = "age", des = design, table.var = "Age", dig = 2 )
  cci <- epitab.means( cont.var = "cci_score", des = design, table.var = "CCI", dig = 2 )
  site <- epitab( var = "primarycagroup", data.fr = df, des = design, table.var = "Cancer Site" )
  snap <- epitab( var = "foodasstpnowic", data.fr = df, des = design, table.var = "SNAP Assistance" )
  time <- epitab( var = "timesince.cat.5yr", data.fr = df, des = design, table.var = "Years Since Diagnosis" )
  insur <- epitab( var = "ins.status", data.fr = df, des = design, table.var = "Health Insurance Status" )
  disab <- epitab.means( cont.var = "adl.score", des = design, table.var = "NHANES ADL Score", dig = 2 )
  ac.mort <- epitab( var = "cod", data.fr = df, des = design, table.var = "Cause of Death" )

  table1 <- rbind( age, sex, race, education, income, insur, hhsize, disab, bmi, metmins, 
                   calor, cci, snap, time, smokstat, alcuse, ac.mort )
  
  return( table1 )
  
}



# generate table columns for each of the subsets described above
fiw.tab <- cafs_table1( design = fiw, df = dat )
fsw.tab <- cafs_table1( design = fsw, df = dat )
gen.tab <- cafs_table1( design = gen, df = dat )

# merge columns into table
final.tab <- cbind( gen.tab, fiw.tab, fsw.tab )


## add column for p value for t tests and chi square test ##

# vector of strings containing elements that maps to a given row in the table
chi  <- c( "Smoking", "Alcohol", "Gender", "Income", "Size", "Education", "Race",
           "Years", "SNAP", "Insurance", "Cause of" )
these  <- c( "smokstat", "alc_cat", "gender", "fipr", "hhsize.bin", "education_bin", "race", 
             "timesince.cat.5yr", "foodasstpnowic", "ins.status", "cod" )

# chi square
for ( i in 1:length( chi ) ){
  
  final.tab[ which( str_detect( final.tab[ , 1 ], chi[i] ) ), 7 ] <- ifelse( svychisq( as.formula( paste0( "~binfoodsechh + ", these[i] ) ), design = gen )$p.value < 0.01, "< 0.01",
                                                                             round( svychisq( as.formula( paste0( "~binfoodsechh + ", these[i] ) ), design = gen )$p.value, digits = 2 ) )
  
} 

# t test variables
tt  <- c( "Age", "BMI", "HHSize", "MET", "Calories", "CCI", "ADL" )
these.b  <- c( "age", "bmxbmi", "hhsize", "weekmetmin", "kcal", "cci_score", "adl.score" ) 


# t test
for ( i in 1:length( tt ) ){
  
  final.tab[ which( str_detect( final.tab[ , 1 ], tt[i] ) ), 7 ] <- ifelse( svyttest( as.formula( paste0( these.b[i], "~binfoodsechh" ) ), design = gen )$p.value < 0.01, "< 0.01",
                                                                            round( svyttest( as.formula( paste0( these.b[i], "~binfoodsechh" ) ), design = gen )$p.value, digits = 2 ) )
  
}

# remove empty "( )" in table
final.tab[ final.tab == "  ( )" ] <- ""

# remove redundant "characteristics" columns
not.these <- vector()
for( i in 2:ncol( final.tab ) ){
  
  nt <- sum(str_detect( final.tab[,i], "[:alpha:]" ) )
  not.these[i] <- ifelse( is.na( nt ), 0,
                          ifelse( i == 1, 0,
                                  ifelse( nt != 0, 1, 0 )))
}

these.cols <- which( not.these != 1 | is.na( not.these) )

# sample sizes for column headings
gen.n <- nrow( gen$variables)
fiw.n <- nrow( fiw$variables)
fsw.n <- nrow( fsw$variables)

t.1 <- setNames( final.tab[, these.cols ],
                       c( "Characteristic",
                          paste0( "Combined Sample (n = ", gen.n ),
                          paste0( "Food Insecure (n = ", fiw.n ),
                          paste0( "Food Secure (n = ", fsw.n ),
                          "p" ) )

# save
write.table( t.1, "04-Tables-Figures/tables/01-table-1.txt", sep = ", ", row.names = FALSE )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Correlation Matrix (Diet Patterns--(Table 2) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# subset individuals meeting criteria
cordata <- dat[ which( dat$inc == 1 ), ]

# food group column indices
fdgrp.columns <- which( colnames( cordata ) %in% c( "processedmts", "meat", "poultry", 
                                                    "fish_hi", "fish_lo", "eggs", 
                                                    "solidfats", "oils", "milk", 
                                                    "yogurt", "cheese", "alcohol", 
                                                    "fruitother", "f_citmelber", 
                                                    "tomatoes", "greenleafy", 
                                                    "darkylveg", "otherveg", 
                                                    "potatoes", "otherstarchyveg", 
                                                    "legumes", "soy", "refinedgrain", 
                                                    "wholegrain", "nuts", "addedsugars" ) )

fdgrp.columns <- fdgrp.columns[ c( 1, 26, 2:25 ) ] # re-arrange so that Meat column index is second column index


# center and scale food group variables before correlation analysis
for ( j in fdgrp.columns ){# ensure proper variables are indicated by the column index in this line of code before proceeding
  cordata[ , j ] <- ( cordata[ , j ] - mean( cordata[ , j ], na.rm = T ) ) / sd( cordata[ , j ], na.rm = T )
}

# re-join adjusted, centered, and scaled variables to original data
cordat2 <- left_join( dat[ , -fdgrp.columns ], cordata[ , c( 1, fdgrp.columns ) ] )

# subset since svy procedures require all rows in data
# to have weights and no missing weights
mod1 <- cordat2 %>%
  filter( is.na( wtdr18yr ) == F) %>%
  svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
             nest = TRUE, survey.lonely.psu = "adjust", data = . ) %>%
  subset( ., inc == 1 )# inclusions



fdgrp.diet.names <- c( "processedmts", "meat", "poultry", "fish_hi", "fish_lo", 
                       "eggs", "solidfats", "oils", "milk", "yogurt", "cheese", 
                       "alcohol", "fruitother", "f_citmelber", "tomatoes", 
                       "greenleafy", "darkylveg", "otherveg", "potatoes", 
                       "otherstarchyveg", "legumes", "soy", "refinedgrain", 
                       "wholegrain", "nuts", "addedsugars", "fs_enet", "age_enet", 
                       "fdas_enet", "hhs_enet", "pc1", "pc2" )

diet.patt.names <- c( "fs_enet", "age_enet", 
                      "fdas_enet", "hhs_enet", "pc1", "pc2" )

## Loop to generate correlation matrix using svycor function ##

# initialize matrix
corr.matrix <- matrix( NA, ncol = length( diet.patt.names ), nrow = length( fdgrp.diet.names ) )

for ( g in 1:length( diet.patt.names ) ){
  
  loadings.vector <- vector( )
  
  for ( i in 1:length( fdgrp.diet.names ) ){
    
    loadings.vector[ i ] <- round( svycor( as.formula( paste0( "~", fdgrp.diet.names[ i ], "+", diet.patt.names[ g ] ) ), design = mod1 )$cors[ 2 ], digits = 2 )
    
  }
  
  corr.matrix[ , g ] <- loadings.vector
  
}


## text-process correlation matrix ##

# assign column names and rownames to matrix

colnames( corr.matrix ) <- diet.patt.names
rownames( corr.matrix ) <- fdgrp.diet.names

# fix significant digits and trailing zeros:
corr.matrix.b <- print( formatC( corr.matrix, digits = 2, format ="fg", flag ="#" ) )

# replace instances of "0" with "0.00"
corr.matrix.b[ corr.matrix.b =="0" ] <- "0.00"

# eliminate excess trailing zeros
for ( i in 1:ncol( corr.matrix.b ) ){
  
  corr.matrix.b[ , i ] <- str_replace( corr.matrix.b[ , i ], "( ?<=\\.\\d\\d )0", "" )
  
}

# replace NA"s with "--"
corr.matrix.b[ corr.matrix.b ==" NA" ] <- "--"

# save table
write.table( corr.matrix.b, "04-Tables-Figures/tables/02-table-2-corr.txt", sep =",", row.names = FALSE )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Diet Patterns Across Food Security Status ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# diet pattern variable names
diet.patt.names <- c( "fs_enet", "age_enet", 
                      "fdas_enet", "hhs_enet", "pc1", "pc2" )

# how we want them presented in the table
diet.patt.names.table <- c( "Food Security Pattern", "Age Pattern", 
                      "SNAP Pattern", "Household Size Pattern", "Modified Western Pattern",
                      "Prudent Pattern" )

diets.t <- data.frame() # initialize frame
for ( i in seq_along( diet.patt.names ) ){
  
  # we use the three survey design objects already specified above( `fiw` `gen`, `fsw`)
  fi.patt <- epitab.means( cont.var = diet.patt.names[i], des = fiw, table.var = diet.patt.names.table[i], dig = 2 )
  fs.patt <- epitab.means( cont.var = diet.patt.names[i], des = fsw, table.var = diet.patt.names.table[i], dig = 2 )
  gen.patt <- epitab.means( cont.var = diet.patt.names[i], des = gen, table.var = diet.patt.names.table[i], dig = 2 )
  
  all.three.c <- cbind( gen.patt, fi.patt, fs.patt )

  diets.t <- rbind( diets.t, all.three.c )
}


# remove redundant "characteristics" columns
not.these <- vector()
for( i in 2:ncol( diets.t ) ){
  
  nt <- sum(str_detect( diets.t[,i], "[:alpha:]" ) )
  not.these[i] <- ifelse( is.na( nt ), 0,
                          ifelse( i == 1, 0,
                                  ifelse( nt != 0, 1, 0 )))
}

these.cols <- which( not.these != 1 | is.na( not.these) )

diets.t <- diets.t[, these.cols ] # subset and keep only non redundant columns


# t test to compare across food security status
for ( i in seq_along( diet.patt.names ) ){
  
  diets.t[ which( str_detect( diets.t[ , 1 ], diet.patt.names.table[i] ) ), 5 ] <- ifelse( svyttest( as.formula( paste0( diet.patt.names[i], "~binfoodsechh" ) ), design = gen )$p.value < 0.01, "< 0.01",
                                                                            round( svyttest( as.formula( paste0( diet.patt.names[i], "~binfoodsechh" ) ), design = gen )$p.value, digits = 2 ) )
  
}

# save table
write.table( diets.t, "04-Tables-Figures/tables/03-table-3-diet-fi.txt", sep =",", row.names = FALSE )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
