###------------------------------------------------------------
###   04-DESCRIPTIVES
###------------------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will generate table 1 (epidemiological characteristics of the 
# study sample), a correlation matrix and a radar chart with the dietary pattern 
# correlations with the dietary patterns loadings, and a radar chart. Also, we
# will generate a table examining differences in the dietary pattern scores
# across food security status.
#
# INPUT DATA FILES: 
# "03-Data-Rodeo/01-analytic-data.rds"
#
# OUPUT FILES:
# i."04-Tables-Figures/figures/04-ggradar-all.png"
# ii. "04-Tables-Figures/tables/01-table-1.txt"
# iii. "04-Tables-Figures/tables/02-table-2-corr.txt"
# iv. "04-Tables-Figures/tables/03-table-3-diet-fi.txt"
#
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( survey )    # survey commands
library( tidyverse )
library( jtools ) # for svycor function
# library( ggradar ) # for radar chart plotting; see note below with modified and debugged code
library( latex2exp ) # for latex in plots

source( "R/old/ggradar-newoption.R" ) # modified ggradar code so that we don't actually need to load in the library above

# helper functions
source( "R/utils.R" )


### 0.0 Read in Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

dat <- readRDS( "03-Data-Rodeo/01-analytic-data.rds" ) %>%
  dplyr::filter( is.na( wtdr18yr ) == F ) %>% # subset to those having non -missing weights 
  
  dplyr::mutate( 
    hhsize.bin = ifelse( hhsize >= 5, 1,
                                        ifelse( hhsize < 5, 0, NA ) ),
    # concatenate cause of death into a single variable for table
    cod = ifelse( castat == 1, "Cancer",
                  ifelse( cvdstat == 1, "Cardiovascular Disease",
                          ifelse( mortstat == 1 & castat != 1 & cvdstat != 1, "Other",
                                  ifelse( mortstat == 0, "Censored",
                                  NA ) ) ) ),
    timesince.cat.4yr = ifelse( timesincecadxmn <= 48, "<= 4 yrs", ">4 yrs"))
  
  # designs
  nhc <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra,
                    nest = TRUE, survey.lonely.psu = "adjust", data = dat )
  
  # create three subsets, food insecure with cancer ( fiw ), food secure with cancer ( fsw ), and combined 
  # food secure and insecure cancer survivors
  fiw <- subset( nhc, binfoodsechh == "Low" & inc == 1 )
  fsw <- subset( nhc, binfoodsechh == "High" & inc == 1 )
  gen <- subset( nhc, inc == 1 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------

  
  
### (1.0) Person-Years and -Months Calculations ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
  
( pm <- sum( gen$variables$stime ) ) # person-months

pm / 12 # person-years 

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
  
  

### (2.0) Table 1 ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Write function to return table 1 with specific variables ##

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
  time <- epitab( var = "timesince.cat.4yr", data.fr = df, des = design, table.var = "Years Since Diagnosis" )
  insur <- epitab( var = "ins.status", data.fr = df, des = design, table.var = "Health Insurance Status" )
  disab <- epitab.means( cont.var = "adl.score", des = design, table.var = "NHANES ADL Score", dig = 2 )
  ac.mort <- epitab( var = "cod", data.fr = df, des = design, table.var = "Cause of Death" )

  table1 <- rbind( age, sex, race, education, income, insur, hhsize, disab, bmi, metmins, 
                   calor, cci, snap, time, smokstat, alcuse, ac.mort )
  
  return( table1 )
  
}

## ---o--- ##


## (2.2) Generate table columns for each of the subsets described above ##

fiw.tab <- cafs_table1( design = fiw, df = dat )
fsw.tab <- cafs_table1( design = fsw, df = dat )
gen.tab <- cafs_table1( design = gen, df = dat )

# merge columns into table
final.tab <- cbind( gen.tab, fiw.tab, fsw.tab )

## ---o--- ##


## (2.3) Add column for p value for t tests and chi square tests ##

# vector of strings containing elements that maps to a given row in the table
chi  <- c( "Smoking", "Alcohol", "Gender", "Income", "Size", "Education", "Race",
           "Years", "SNAP", "Insurance", "Cause of" )
these  <- c( "smokstat", "alc_cat", "gender", "fipr", "hhsize.bin", "education_bin", "race", 
             "timesince.cat.4yr", "foodasstpnowic", "ins.status", "cod" )

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

## ---o--- ##


## (2.4) Final table clean up ##

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
                          paste0( "Combined Sample (n = ", gen.n, ")" ),
                          paste0( "Food Insecure (n = ", fiw.n, ")" ),
                          paste0( "Food Secure (n = ", fsw.n, ")" ),
                          "p" ) )

# save
write.table( t.1, "04-Tables-Figures/tables/01-table-1.txt", sep = ", ", row.names = FALSE )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Correlation Matrix (Diet Patterns--(Table 2) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Prep data ##

# subset individuals meeting criteria
cordata <- dat[ which( dat$inc == 1 & dat$test.set == 1 ), ]

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
  subset( ., inc == 1 & test.set == 1 ) # inclusions AND performed only on testing set



fdgrp.diet.names <- c( "processedmts", "meat", "poultry", "fish_hi", "fish_lo", 
                       "eggs", "solidfats", "oils", "milk", "yogurt", "cheese", 
                       "alcohol", "fruitother", "f_citmelber", "tomatoes", 
                       "greenleafy", "darkylveg", "otherveg", "potatoes", 
                       "otherstarchyveg", "legumes", "soy", "refinedgrain", 
                       "wholegrain", "nuts", "addedsugars", "fs_enet", 
                       "pc1", "pc2", "hei.2015" )

diet.patt.names <- c( "fs_enet", "pc1", "pc2", "hei.2015" )

## ---o--- ##


## (3.2) Loop to generate correlation matrix using `svycor` function ##

# initialize matrix
corr.matrix <- matrix( NA, ncol = length( diet.patt.names ), nrow = length( fdgrp.diet.names ) )

for ( g in 1:length( diet.patt.names ) ){
  
  loadings.vector <- vector( )
  
  for ( i in 1:length( fdgrp.diet.names ) ){
    
    loadings.vector[ i ] <- round( svycor( as.formula( paste0( "~", fdgrp.diet.names[ i ], "+", diet.patt.names[ g ] ) ), design = mod1 )$cors[ 2 ], digits = 2 )
    
  }
  
  corr.matrix[ , g ] <- loadings.vector
  
}

## ---o--- ##


## (3.3) Text-process correlation matrix ##

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




### (4.0) Radar Chart for Select Dietary Patterns ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## NOTE: code below requires that some of the above code is run prior to running this one.
# specifically, the code chunk in the (3.0) Correlation Matrix (Diet Patterns--(Table 2)
# section needs to be run

## (4.1) Prep names and labels for plot ##

# axis labels
fg.only <- c( "Processed\n Meat", "Meat", "Poultry", "Fish,\n High n-3", "Fish,\n Low n-3", 
              "Eggs", "Solid\n Fats", "Oils", "Milk", "Yogurt", "Cheese", 
              "Alcohol", "Fruit\n Other", "Citrus,\nMelons,\nBerries",  
              "Green\n Leafy\n Veg.", "Dark Ylw.\n Veg.", "Tomatoes", "Other\n Veg.", "Potatoes", 
              "Other\n Starchy\n Veg.", "Legumes ", "Soy", "Refined\n Grain", 
              "Whole\n Grain", "Nuts", "Added\n Sugars" )

fg.only.lo <- c( "processedmts", "meat", "poultry", "fish_hi", "fish_lo", 
                 "eggs", "solidfats", "oils", "milk", "yogurt", "cheese", 
                 "alcohol", "fruitother", "f_citmelber", 
                 "greenleafy", "darkylveg", "tomatoes", "otherveg", "potatoes", 
                 "otherstarchyveg", "legumes", "soy", "refinedgrain", 
                 "wholegrain", "nuts", "addedsugars" )

rp.matrix <- corr.matrix[ fg.only.lo, ]

## ---o--- ##


## (4.2) Radar chart using `ggradar` ##

# prepare frame for plotting (items/variables in columns, dietary patterns in rows)
rp.all.gg <- setNames( data.frame( bind_rows(
  setNames( matrix( rp.matrix[, 1 ], nrow = 1 ), fg.only ),
  setNames( matrix( rp.matrix[, 2 ], nrow = 1 ), fg.only ),
  setNames( matrix( rp.matrix[, 3 ], nrow = 1 ), fg.only ),
  setNames( matrix( rp.matrix[, 4 ], nrow = 1 ), fg.only ), ) ),
  fg.only ) %>%
  
  mutate( across( where( is.numeric ), ~ ( .x + 1 ) / 2 ) ) %>% # map -1-1 scale to 0-1 scale
  data.frame()

# set column names and row names
colnames( rp.all.gg ) <- c ( fg.only.lo )
rownames( rp.all.gg ) <- c( "Pattern #1", "Pattern #2", "Pattern #3",
                            "HEI-2015" )

# need to add the group names as the first column of the dataset you feed to `ggradar` since it looks for those in column #1
# and begins looking at numeric values in column #2 and onward
rp.all.gg <- rp.all.gg %>%
  mutate( group = rownames( . ) ) %>%
  relocate( group, .before = processedmts )

# line colors
colors.line <- c( rgb(0.529,0.808,0.98,0.9), rgb(0,0,0.502,0.9), rgb(0.55,0.10,0.10,0.9) , "goldenrod2" )

# plot
( g.1 <- ggradar( rp.all.gg, values.radar = c("-1.0", "0.0", "1.0"),
                  font.grid.label = "Avenir Next LT Pro Bold", # new option added (modified source code)
                  grid.label.size =  8,
                  font.radar = "Avenir",
                  axis.labels = fg.only,
                  legend.position = c(0.95,0.85),
                  legend.title = "Dietary Pattern",
                  axis.label.size = 7.4,
                  x.centre.range = 1.31 ,
                  legend.labels = c( unname( TeX( "Pattern #1$^{\\dagger a}$" ) ),
                                     unname( TeX("Pattern #2$^{\\ddagger b}$" ) ),
                                     unname( TeX("Pattern #3$^{\\ddagger c}$" ) ),
                                     unname( TeX( "HEI-2015$^d$" ) ) ),
                  group.colours = colors.line ) +
    theme( legend.title = element_text( face = "bold"),
           legend.text = element_text( size = 19, margin = margin(t = 10, b = 10, l = 0.1) ), # margins for legend text
           legend.background = element_rect(fill='transparent'),
           legend.text.align = 0,
           legend.box.margin = unit(c(1,1,1,1.8), 'cm'), # margins for legend box
           plot.margin=unit(c(0.1,1.7,3.9,1), 'cm') ) )  # margins for plot
  

# change axis label text to grey
g.1$layers[[5]]$aes_params <- c( g.1$layers[[5]]$aes_params, colour = "grey40" )


ggsave( "04-Tables-Figures/figures/04-ggradar-all.png",
        width = 14.8, 
        height = 11.9,
        plot = g.1,
        dpi = 800 )

ggsave( "04-Tables-Figures/figures/04-ggradar-all.tiff",
        width = 14.8, 
        height = 11.9,
        plot = g.1,
        dpi = 500 )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) Diet Patterns Across Food Security Status ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (5.1) Variables and labels ##

# diet pattern variable names
diet.patt.names <- c( "fs_enet", "pc1", "pc2", "hei.2015" )

# how we want them presented in the table
diet.patt.names.table <- c( "Pattern #1", "Pattern #2",
                      "Pattern #3", "HEI-2015" )

## ---o--- ##

## (5.2) Generate results with `svy` commands ##

diets.t <- data.frame() # initialize frame
for ( i in seq_along( diet.patt.names ) ){
  
  # we use the three survey design objects already specified above( `fiw` `gen`, `fsw`)
  fi.patt <- epitab.means( cont.var = diet.patt.names[i], des = subset( mod1, binfoodsechh == "Low" ), table.var = diet.patt.names.table[i], dig = 2 )
  fs.patt <- epitab.means( cont.var = diet.patt.names[i], des = subset( mod1, binfoodsechh == "High" ), table.var = diet.patt.names.table[i], dig = 2 )
  gen.patt <- epitab.means( cont.var = diet.patt.names[i], des = mod1, table.var = diet.patt.names.table[i], dig = 2 )
  
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

## ---o--- ##


## (5.3) Cohen's D ##

for ( i in seq_along( diet.patt.names ) ){
  
  cd.col <- round( svycd( x = diet.patt.names[i], design.1 = subset( mod1, binfoodsechh == "Low" ), 
                          design.2 = subset( mod1, binfoodsechh == "High" ) ), 2 )
  
  diets.t[ which( str_detect( diets.t[ , 1 ], diet.patt.names.table[i] ) ), 5 ] <- cd.col
}

## ---o--- ##


## (5.4) T-test to compare across food security status ##
for ( i in seq_along( diet.patt.names ) ){
  
  diets.t[ which( str_detect( diets.t[ , 1 ], diet.patt.names.table[i] ) ), 6 ] <- ifelse( svyttest( as.formula( paste0( diet.patt.names[i], "~binfoodsechh" ) ), design = mod1 )$p.value < 0.01, "< 0.01",
                                                                            round( svyttest( as.formula( paste0( diet.patt.names[i], "~binfoodsechh" ) ), design = mod1 )$p.value, digits = 2 ) )
  
}


# save table
write.table( diets.t, "04-Tables-Figures/tables/03-table-3-diet-fi.txt", 
             sep ="," )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------





