###---------------------------------------------------
###   03-DIETARY PATTERNS EXTRACTION
###---------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will extract dietary patterns using penalized logistic regression and then principal components
# analysis. all dietary patterns scores will be energy adjusted in the final dataset.
# 
# INPUT DATA FILE: "02-Data-Wrangled/02-inclusions-exclusions.rds" 
#
# OUTPUT FILES: "03-Data-Rodeo/01-analytic-data.rds"
#
# Resources: 
# i. reference for penalized logistic regression: https://teazrq.github.io/SMLR/logistic-regression.html
# ii. adjustment for total energy paper: https://doi.org/10.1093/ajcn/65.4.1220S
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( glmnet ) # fit penalized regression models
library( caret ) # control tuning parameters in penalized regression
library( tidyverse )
library( survey ) # complex survey design models
library( jtools ) # svycor function
library( weights )
library( latex2exp ) # to add LaTeX to plots

source( "R/utils.R" ) # read in helper and internal functions

# data read-in
dat <- readRDS( "02-Data-Wrangled/02-inclusions-exclusions.rds" )

# collapse red mt and organ mt to same group give very low intake of organ mts
dat$meat <- dat$redmts + dat$organmts


### (0.0) Prep Data for Penalized Logistic Regression ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# subset data for those included in analysis
caonly <- dat[ which( dat$inc == 1 ), ] # those included after inclusions/exclusions steps

# random sample of included subjects for training the penalized logit models #
set.seed( 0872645 ) # set seed for reproducibility
sz <- ceiling( nrow( caonly ) * 0.5 ) # size of sample (50-50 split into training and test sets)
train <- caonly[ caonly$seqn %in% sample( x = caonly$seqn, size = sz ), ] # training dataset
test <- caonly[ caonly$seqn %notin% train$seqn, ] # test set for later in the analysis

# set categories to numerical dummy variables for glmnet model
train <- train %>%
  mutate( foodasstpnowic = ifelse( foodasstpnowic =='yes', 1,
                                   ifelse( foodasstpnowic =='no', 0, NA ) ) ) %>%
  mutate( agecat = ifelse( agecat =='elderly', 1,
                           ifelse( agecat =='non-elderly', 0, NA ) ) ) %>%
  mutate( binfsh = ifelse( binfoodsechh =='Low', 1,
                           ifelse( binfoodsechh =='High', 0, NA ) ) ) %>%
  mutate( hhsize_bin = ifelse( hhsize>= 5, 1,
                               ifelse( hhsize<5, 0, NA ) ) )

# food groups columns used for the procedure
fdgrp.columns <- which( colnames( train ) %in% c( 'processedmts','meat','poultry',
                                                  'fish_hi','fish_lo','eggs',
                                                  'solidfats','oils','milk',
                                                  'yogurt','cheese','alcohol',
                                                  'fruitother','f_citmelber',
                                                  'tomatoes','greenleafy',
                                                  'darkylveg','otherveg','potatoes',
                                                  'otherstarchyveg','legumes',
                                                  'soy','refinedgrain','wholegrain',
                                                  'nuts','addedsugars') )

fdgrp.columns <- fdgrp.columns[ c( 1, 26, 2:25 ) ] # re-arrange so that meat column index is second column index

# other relevant columns and their indices that we will need for the procedure
fs.outcome.column <- which( colnames( train ) =='binfsh' )
fdas.outcome.column <- which( colnames( train ) =='foodasstpnowic' )
age.outcome.column <- which( colnames( train ) =='agecat' )
hhsize.outcome.column <- which( colnames( train ) =='hhsize_bin' )
weight.column <- which( colnames( train ) =='wtdr18yr' )
kcal.column <- which( colnames( train ) =='kcal' )
seqn.column <- which( colnames( train ) =='seqn' )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) Patterns Extraction with Penalized Logistic Regression ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# normalize variables prior to regression since they are all in different units

for ( j in fdgrp.columns ){ # ensure proper variables are indicated by the column index in this line of code before proceeding
  
  train[ , j ] <- ( train[ , j ] - mean( train[ , j ], na.rm = TRUE ) ) / sd( train[ , j ], na.rm = TRUE )

  }

# **NOTE: we will include the calories column in these procedures so that the coefficients for the food groups are adjusted for calories (standard multivariate method)
# internal function is `enet_pat` that is written in the "utils.R" file

par( mfrow = c( 2, 2 ), mar = c( 3, 3, 2, 1 ) ) # multi panel figure for finding optimizer

## food insecurity binary outcome/dietary pattern ##
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, fs.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights

fsoc <- enet_pat( xmat, yvec, xwts, plot.title = 'Food Insecurity' ) 

# add legend
colorss <- c( "black", "red", "green3", "navyblue",   "cyan",   "magenta", "gold", "gray",
              'pink', 'brown', 'goldenrod' )
legend( "bottomright", 
        legend = c( paste( seq( 0, 1, by = 0.1 ) ), 'Minimizer' ), col = c( colorss, 'grey' ),
        lty = c( rep( 1, 11 ), 2 ), title = TeX( '$\\alpha$' ),
        cex = 0.6, inset = 0, y.intersp = 0.5 )


## age outcome/dietary pattern ##
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, age.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights

ageoc <- enet_pat( xmat, yvec, xwts, plot.title = 'Age' )



## receipt of food assistance outcome/dietary pattern ##
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, fdas.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights

fdasoc <- enet_pat( xmat, yvec, xwts, plot.title = 'Food Assistance ( SNAP )' )


## household size outcome/dietary pattern ##
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, hhsize.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights

hhssoc <- enet_pat( xmat, yvec, xwts, plot.title ='Household Size' )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) Generate Pattern Scores in the Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# generate scores on data
xmatrix <- as.matrix( dat[ which( dat$inc == 1 ), c( fdgrp.columns ) ] ) # subset as a matrix for matrix multiplication in next step
d <- dat[ which( dat$inc == 1 ), ] # keep only those satisfying inclusions/exclusions as dataframe to then merge matrix products back into dataframe format

# center and scale testing data (in the matrix format) before generating scores
for ( i in 1:ncol( xmatrix ) ){
  xmatrix[ , i ] <- ( xmatrix[ , i ] - mean( xmatrix[ , i ], na.rm = T ) ) / sd( xmatrix[ , i ], na.rm = T )
}

# matrix multiplication and add score columns back to dataframe
d$fs_enet <- t( fsoc$coefs[1:length(fdgrp.columns)] %*% t( xmatrix ) )
d$age_enet <- t( ageoc$coefs[1:length(fdgrp.columns)]  %*% t( xmatrix ) )
d$fdas_enet <- t( fdasoc$coefs[1:length(fdgrp.columns)]  %*% t( xmatrix ) )
d$hhs_enet <- t( hhssoc$coefs[1:length(fdgrp.columns)]  %*% t( xmatrix ) )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Dietary Patterns Extraction with PCA ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# to have weights and no missing weights

svy.design <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                         nest = TRUE, survey.lonely.psu = "adjust", data = dat )

varest <- subset( svy.design, diet.ext.ind.pca == 1 &
                    seqn %in% train$seqn ) # inclusions


# pca using svyprcomp
pcaobj <- svyprcomp( ~ processedmts + meat + poultry + 
                     fish_hi + fish_lo + eggs + 
                     solidfats + oils + milk + 
                     yogurt + cheese + alcohol + 
                     fruitother + f_citmelber + 
                     tomatoes + greenleafy + 
                     darkylveg + otherveg + potatoes + 
                     otherstarchyveg + legumes + 
                     soy + refinedgrain + wholegrain + 
                     nuts + addedsugars, design = varest, center = T,
                     scale = T, # center and scale variables using internal arguments
                     scores = TRUE )

# scree plot
plot( pcaobj$sdev, type = 'b',
      main = 'Scree Plot for Diet Patterns PCA', xlab = 'Component', ylab = 'Eigenvalue' )
abline( h = 1, lty = 2 )
# elbow present after the second principal component

# assign factor loading matrix
coefspc <- pcaobj$rotation

# percent of variation accounted for by first two components ( 15.79 % )
sum( pcaobj$sdev[ 1:2 ] / sum( pcaobj$sdev ) ) 

# generate scores
xmatrix <- as.matrix( dat[ which( dat$inc == 1 ), c( fdgrp.columns ) ] )

# center and scale testing data before generating scores
for ( i in 1:ncol( xmatrix ) ){
  xmatrix[ , i ] <- ( xmatrix[ , i ] - mean( xmatrix[ , i ], na.rm = T ) ) / sd( xmatrix[ , i ], na.rm = T )
}

# generate scores using matrix multiplication
d$pc1 <- t( coefspc[ , 1 ] %*% t( xmatrix ) )
d$pc2 <- t( coefspc[ , 2 ] %*% t( xmatrix ) )

# add to original data and save
d.2 <- left_join( dat, d[ , c( "seqn", "fs_enet", "age_enet",
                                       "fdas_enet", "hhs_enet", "pc1", "pc2" ) ] ) %>%
    mutate( test.set = ifelse( seqn %in% test$seqn, 1, 0 ) ) # add indicator variable for membership in test set

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Energy Adjust Principal Component Scores using Residual Method (Willett) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


svy.design.2 <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                         nest = TRUE, survey.lonely.psu = "adjust", data = d.2 )

adj.des <- subset( svy.design.2, inc == 1 ) # survey design object

d.3 <- svy_energy_residual( nutr = c( "pc1", "pc2" ), # columns to be energy adjusted
                            design = adj.des, # design object
                            calories = "kcal", # calories column
                            overwrite = "yes" ) # keep original column names

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) Save Analytic Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( d.3, "03-Data-Rodeo/01-analytic-data.rds" )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



