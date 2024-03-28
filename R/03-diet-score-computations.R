###---------------------------------------------------
###   03-DIETARY PATTERNS EXTRACTION
###---------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will extract dietary patterns using penalized logistic 
# regression and then principal components analysis. We also separate the data 
# into two random subsets, a "training" subset, for extracting the dietary patterns, 
# and a "testing" subset, that is used later for the validation survival analysis. 
# All dietary patterns scores will be energy adjusted in the final dataset.
# 
# INPUT DATA FILE: 
# i."02-Data-Wrangled/02-inclusions-exclusions.rds" 
#
# OUTPUT FILES: 
# i. "03-Data-Rodeo/01-analytic-data.rds"
#
# Resources: 
# i. reference for penalized logistic regression: https://teazrq.github.io/SMLR/logistic-regression.html
# ii. adjustment for total energy paper: https://doi.org/10.1093/ajcn/65.4.1220S
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( glmnet ) # fit penalized regression models
library( caret ) # control tuning parameters in penalized regression
library( tidyverse )
library( survey ) # complex survey design models and commands
library( jtools ) # svycor function
library( weights )
library( latex2exp ) # to add LaTeX to plots

source( "R/utils.R" ) # read in helper and internal functions


### (0.0) Prep Data for Penalized Logistic Regression ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (0.1) Read-in inclusions/exclusions data ##

# data read-in
dat <- readRDS( "02-Data-Wrangled/02-inclusions-exclusions.rds" )

# collapse red mt and organ mt to same group give very low intake of organ mts
dat$meat <- dat$redmts + dat$organmts

# subset data for those included in analysis
caonly <- dat[ which( dat$inc == 1 ), ] # those included after inclusions/exclusions steps

## ---o--- ##


## (0.2) Generate random sample of included subjects for training the penalized logit models ##

set.seed( 0872645 ) # set seed for reproducibility

sz <- ceiling( nrow( caonly ) * 0.3 ) # size of sample (30-70 split into training and test sets)
train <- caonly[ caonly$seqn %in% sample( x = caonly$seqn, size = sz ), ] # training dataset
test <- caonly[ caonly$seqn %notin% train$seqn, ] # test set for later in the analysis

# set categories to numerical dummy variables for glmnet model
train <- train %>%
  mutate( binfsh = ifelse( binfoodsechh =='Low', 1,
                           ifelse( binfoodsechh =='High', 0, NA ) ) )

## ---o--- ##


## (0.3) Set up names and labels of food groups for the extraction analyses ##

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
weight.column <- which( colnames( train ) =='wtdr18yr' )
kcal.column <- which( colnames( train ) =='kcal' )
seqn.column <- which( colnames( train ) =='seqn' )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (1.0) Patterns Extraction with Penalized Logistic Regression ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (1.1) Normalize variables prior to regression since they are all in different units ##

for ( j in fdgrp.columns ){ # ensure proper variables are indicated by the column index in this line of code before proceeding
  
  train[ , j ] <- ( train[ , j ] - mean( train[ , j ], na.rm = TRUE ) ) / sd( train[ , j ], na.rm = TRUE )

  }

# **NOTE: we will include the calories column in these procedures so that the coefficients for the food groups are adjusted for calories (standard multivariate method)
# internal function is `enet_pat` that is written in the "utils.R" file

## ---o--- ##


## (1.2) Setup matrices and vectors to feed to `glmnet` functions ## 

# food insecurity binary outcome/dietary pattern (i.e., Pattern #1)
design.matrix <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, fs.outcome.column, weight.column ) ] ) )
matrix.x <- design.matrix[ , 1:27 ] # food grps and kcal data matrix
vector.y <- design.matrix[ , 28 ] # binary response vector
vector.wts <- design.matrix[ , 29 ] / mean( design.matrix[ , 29 ] ) # generate vector of normalized weights

## ---o--- ##

## (1.3) Use internal function `enet_pat` to extract dietary pattern using penalized logit ## 

fsoc <- enet_pat( xmat = matrix.x, 
                  yvec = vector.y, 
                  wts = vector.wts,
                  plot.title = "Food Insecurity" ) # Pattern # 1

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (2.0) Generate Pattern #1 Scores in the Dataset ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (2.1) Prep Data ##

xmatrix <- as.matrix( dat[ which( dat$inc == 1 ), c( fdgrp.columns ) ] ) # subset as a matrix for matrix multiplication in next step
d <- dat[ which( dat$inc == 1 ), ] # keep only those satisfying inclusions/exclusions as dataframe to then merge matrix products back into dataframe format

## ---o--- ##


## (2.2) Center and scale data (in the matrix format) before generating scores ##
for ( i in 1:ncol( xmatrix ) ){
  xmatrix[ , i ] <- ( xmatrix[ , i ] - mean( xmatrix[ , i ], na.rm = T ) ) / sd( xmatrix[ , i ], na.rm = T )
}

## ---o--- ##


## (2.3) Matrix multiplication and add score columns back to dataframe ##
d$fs_enet <- t( fsoc$coefs[ 1:length(fdgrp.columns) ] %*% t( xmatrix ) )

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (3.0) Dietary Patterns (Patterns #2 and #3) Extraction with PCA ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (3.1) Survey design object ##

svy.design <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                         nest = TRUE, survey.lonely.psu = "adjust", data = dat )

varest <- subset( svy.design, diet.ext.ind.pca == 1 &
                    seqn %in% train$seqn ) # inclusions (training data)

## ---o--- ##


## (3.2) PCA using `svyprcomp` ##

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

## ---o--- ##


## (3.3) Generate index scores in the dataset for patterns #2 and #3 ##

xmatrix <- as.matrix( dat[ which( dat$inc == 1 ), c( fdgrp.columns ) ] )

# center and scale testing data before generating scores
for ( i in 1:ncol( xmatrix ) ){
  xmatrix[ , i ] <- ( xmatrix[ , i ] - mean( xmatrix[ , i ], na.rm = T ) ) / sd( xmatrix[ , i ], na.rm = T )
}

# generate scores using matrix multiplication
d$pc1 <- t( coefspc[ , 1 ] %*% t( xmatrix ) )
d$pc2 <- t( coefspc[ , 2 ] %*% t( xmatrix ) )

# add to original data and save
d.2 <- left_join( dat, d[ , c( "seqn", "fs_enet", "pc1", "pc2" ) ] ) %>%
    mutate( test.set = ifelse( seqn %in% test$seqn, 1, 0 ) ) # add indicator variable for membership in test set

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (4.0) Energy Adjust Principal Component Scores using Residual Method (Willett) ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

## (4.1) Survey design object ##

svy.design.2 <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                         nest = TRUE, survey.lonely.psu = "adjust", data = d.2 )

adj.des <- subset( svy.design.2, inc == 1 ) # survey design object

## ---o--- ##


## (4.2) Use internal helper function (`svy_energy_residual`) for energy adjusting the pattern scores ##

d.3 <- svy_energy_residual( nutr = c( "pc1", "pc2", "hei.2015" ), # columns to be energy adjusted
                            design = adj.des, # design object
                            calories = "kcal", # calories column
                            overwrite = "yes" ) %>% # keep original column names
  select( seqn, fs_enet, pc1, pc2 ) %>%
  left_join( dat, ., by = "seqn" ) %>%
  mutate( test.set = ifelse( seqn %in% test$seqn, 1, 0 ) ) # add indicator variable for membership in test set

# ---------------------------------------------------------------------------------------------------------------------------------------------------------



### (5.0) Save Analytic Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( d.3, "03-Data-Rodeo/01-analytic-data.rds" )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------



