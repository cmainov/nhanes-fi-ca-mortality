###---------------------------------------------------
###   03-DIETARY PATTERNS EXTRACTION
###---------------------------------------------------

# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# 
# In this script, we will extract dietary patterns using penalized logistic regression and then principal components
# analysis.
# 
# INPUT DATA FILE: "02-Data-Wrangled/02-inclusions-exclusions.rds" 
#
# OUTPUT FILES: "03-Data-Rodeo/01-analytic-data.rds"
#
# Resources: 
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

library( glmnet ) # fit penalized regression models
library( caret ) # control tuning parameters in penalized regression
library( tidyverse )
library( survey ) # complex survey design models
library( jtools ) # svycor function
library( weights )
library( latex2exp ) # to add LaTeX to plots


xdata <- readRDS( "02-Data-Wrangled/02-inclusions-exclusions.rds" )


# collapse red mt and organ mt to same group give very low intake of organ mts
xdata$meat <- xdata$redmts + xdata$organmts

# copy
x.data <- xdata


### Dietary Patterns Extraction: Penalized Logistic Regression ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


# function for extracting patterns with penalized logit

## START FUNCTION ##
enet_pat <- function( xmat, yvec, wts, plot.title ){
  
  colorss <- c( "black", "red", "green3", "navyblue",   "cyan",   "magenta", "gold", "gray",
                'pink', 'brown', 'goldenrod' ) # colors for plot
  
  # initialize lists to store outputs
  store <- list( )
  coefsdt <- list( )
  
  alpha.grid <- seq( 0, 1, 0.1 ) # range of alpha values for tuning grid
  
  for ( i in 1:length( alpha.grid ) ){ # set the grid of alpha values
    
    set.seed( 28 ) # seed to reproduce results
    
    # call glmnet with 10-fold cv
    enetr <- cv.glmnet( x = xmat, y = yvec, family ='binomial', weights = wts,
                        nfold = 10, alpha = alpha.grid[ i ] )
    
    # bind values of lambda to cross validation error
    evm <- data.frame( cbind( enetr$lambda, enetr$cvm ) )
    colnames( evm ) <- c( 'lambda', 'av.error' )
    
    # now create a list that stores each coefficients matrix for each value of alpha
    # at the lambda minimizer
    coefsdt[[ i ]] <- list( alpha = paste0( alpha.grid )[ i ], 
                            coefs = coef( enetr, s = "lambda.min" ) )
    
    # create a dataframe that houses the alpha, min lambda, and average error
    resdf <- data.frame( alpha = alpha.grid[ i ], 
                         evm[ which( evm$av.error == min( evm$av.error ) ), 'lambda' ],
                         av.error = evm[ which( evm$av.error == min( evm$av.error ) ), 'av.error' ] )
    colnames( resdf ) <- c( 'alpha', 'lambda', 'av.error' )
    
    store[[ i ]] <- resdf
    
    ## generate plot ##
    
    if ( i == 1 ){ # for the first value of 'i'
      plot( x = enetr$lambda, y = enetr$cvm, type ='l', 
            ylim = c( min( enetr$cvm )-0.02, max( enetr$cvm )-0.02 ),
            xlim = c( min( evm$lambda ), ( resdf$lambda*1.05 ) ), 
            las = 0, 
            cex.axis = 0.7 )
    }
    else if ( i != 1 ){ # each additional line will be superimposed on the plot with a different color
      lines( x = enetr$lambda, 
             y = enetr$cvm, 
             col = colorss[ i ] )
    }
  }
  
  ## superimpose intersecting lines at the minimizer ## 
  cverr <- do.call( 'rbind', store ) # this gives the table of errors for each combination of alpha and lambda
  abline( h = cverr[ which( cverr$av.error == min( cverr$av.error ) ), 'av.error' ],
          lty = 2 )
  abline( v = cverr[ which( cverr$av.error == min( cverr$av.error ) ), 'lambda' ],
          lty = 2 )
  
  
  ## add optimal lambda and alpha values to plot title ## 
  optimall <- cverr[ which( cverr$av.error == min( cverr$av.error ) ), ] # here I extract the optimal combination of
  
  # lambda and alpha
  optlam <- signif( optimall[ 2 ], 2 )
  opta <- optimall[ 1 ]
  title( main = TeX( paste0( plot.title, ' ( $\\lambda_{optimal} =$', optlam, ' and $\\alpha_{optimal} =$', opta, ' )' ) ),
         cex.main = 0.8,
         cex.lab = 0.8,
         xlab = TeX( '$\\lambda$' ), 
         mgp = c( 2, 1, 0 ),
         ylab ='Deviance', mgp = c( 2, 1, 0 ) )
  
  
  
  # the function returns the optimal lambda alpha combo and the set of coefficients that 
  # correspond to that combination of parameters
  return( list( optimall, coefs = as.matrix( coefsdt[[ which( alpha.grid == optimall$alpha ) ]]$coefs )[ -1, ] ) )
}
## END FUNCTION ##



## Subset Data for Penalized Logit Procedure ##
caonly <- x.data[ which( x.data$inc == 1 ), ]

# random sample of included
set.seed( 0872645 ) # set seed for reproducibility
sz <- ceiling( nrow( caonly ) * 0.5 ) # size of sample
train <- caonly[ caonly$seqn %in% sample( x = caonly$seqn, size = sz ), ] # training dataset
test <- caonly[ caonly$seqn %notin% train$seqn, ]

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

fdgrp.columns <- fdgrp.columns[ c( 1, 26, 2:25 ) ] # re-arrange so that Meat column index is second column index

fs.outcome.column <- which( colnames( train ) =='binfsh' )
fdas.outcome.column <- which( colnames( train ) =='foodasstpnowic' )
age.outcome.column <- which( colnames( train ) =='agecat' )
hhsize.outcome.column <- which( colnames( train ) =='hhsize_bin' )
weight.column <- which( colnames( train ) =='wtdr18yr' )
kcal.column <- which( colnames( train ) =='kcal' )
seqn.column <- which( colnames( train ) =='seqn' )



##############################################################################
#################### Run Function for Patterns Extraction #################### 
##############################################################################

for ( j in fdgrp.columns ){# ensure proper variables are indicated by the column index in this line of code before proceeding
  train[ , j ] <- ( train[ , j ] - mean( train[ , j ], na.rm = TRUE ) ) / sd( train[ , j ], na.rm = TRUE )
}

# Save plot figure

# Food insecurity binary outcome
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, fs.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights

fsoc <- enet_pat( xmat, yvec, xwts, plot.title ='Food Insecurity' ) 

# add legend
colorss <- c( "black", "red", "green3", "navyblue",   "cyan",   "magenta", "gold", "gray",
              'pink', 'brown', 'goldenrod' )
legend( "bottomright", 
        legend = c( paste( seq( 0, 1, by = 0.1 ) ), 'Minimizer' ), col = c( colorss, 'grey' ),
        lty = c( rep( 1, 11 ), 2 ), title = TeX( '$\\alpha$' ),
        cex = 0.6, inset = 0, y.intersp = 0.5 )


# Age outcome
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, age.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights
ageoc <- enet_pat( xmat, yvec, xwts, plot.title = 'Age' )



# Receipt of food assistance outcome
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, fdas.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights
fdasoc <- enet_pat( xmat, yvec, xwts, plot.title = 'Food Assistance ( SNAP )' )


# HHSize outcome
mat <- na.omit( as.matrix( train[ , c( fdgrp.columns, kcal.column, hhsize.outcome.column, weight.column ) ] ) )
xmat <- mat[ , 1:27 ] # food grps
yvec <- mat[ , 28 ] # binary response
xwts <- mat[ , 29 ] / mean( mat[ , 29 ] ) # normalize weights
hhssoc <- enet_pat( xmat, yvec, xwts, plot.title ='Household Size' )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Generate Pattern Scores in the Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# generate scores on data
xmatrix <- as.matrix( x.data[ which( x.data$inc == 1 ), c( fdgrp.columns, kcal.column ) ] ) # subset as a matrix
d <- x.data[ which( x.data$inc == 1 ), ] # keep only those satisfying inclusions/exclusions

# remove kcal column
xmatrix <- xmatrix[ , -ncol( xmatrix ) ]

# center and scale testing data before generating scores
for ( i in 1:ncol( xmatrix ) ){
  xmatrix[ , i ] <- ( xmatrix[ , i ] - mean( xmatrix[ , i ], na.rm = T ) ) / sd( xmatrix[ , i ], na.rm = T )
}


# matrix multiplication
d$fs_enet <- t( fsoc$coefs[1:26] %*% t( xmatrix ) )
d$age_enet <- t( ageoc$coefs[1:26]  %*% t( xmatrix ) )
d$fdas_enet <- t( fdasoc$coefs[1:26]  %*% t( xmatrix ) )
d$hhs_enet <- t( hhssoc$coefs[1:26]  %*% t( xmatrix ) )


# ---------------------------------------------------------------------------------------------------------------------------------------------------------




### Patterns Extraction with PCA ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------

# to have weights and no missing weights

pca.design <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                         nest = TRUE, survey.lonely.psu = "adjust", data = x.data )

varest <- subset( pca.design, diet.ext.ind.pca == 1 &
                    seqn %in% train$seqn ) # inclusions


# PCA using svyprcomp
pcaobj <- svyprcomp( ~ processedmts+meat+poultry+
                     fish_hi+fish_lo+eggs+
                     solidfats+oils+milk+
                     yogurt+cheese+alcohol+
                     fruitother+f_citmelber+
                     tomatoes+greenleafy+
                     darkylveg+otherveg+potatoes+
                     otherstarchyveg+legumes+
                     soy+refinedgrain+wholegrain+
                     nuts+addedsugars, design = varest, center = T,
                     scale = T, 
                     scores = TRUE )

# scree plot
# eigenvalues/Scree plot
plot( pcaobj$sdev, type = 'b',
      main = 'Scree Plot for Diet Patterns PCA', xlab = 'Component', ylab = 'Eigenvalue' )
abline( h = 1, lty = 2 )
# elbow present after the second principal component

# assign factor loading matrix
coefspc <- pcaobj$rotation

# save raw loading matrix
# write.table( round( coefspc[ , 1:2 ], digits = 2 ), "04-Manuscript/Tables/PCA-Factorloadings.txt", sep =", ", row.names = FALSE )


sum( pcaobj$sdev[ 1:2 ] / sum( pcaobj$sdev ) ) # percent of variation accounted for by first two components ( 0.1412 )

# generate scores
xmatrix <- as.matrix( x.data[ which( x.data$inc == 1 ), c( fdgrp.columns, kcal.column ) ] )

# remove kcal column
xmatrix <- xmatrix[ , -ncol( xmatrix ) ]

# center and scale testing data before generating scores
for ( i in 1:ncol( xmatrix ) ){
  xmatrix[ , i ] <- ( xmatrix[ , i ] - mean( xmatrix[ , i ], na.rm = T ) ) / sd( xmatrix[ , i ], na.rm = T )
}

# generate scores using matrix multiplication
d$pc1 <- t( coefspc[ , 1 ] %*% t( xmatrix ) )
d$pc2 <- t( coefspc[ , 2 ] %*% t( xmatrix ) )

# add to original data and save

( x.data3 <- left_join( xdata, d[ , c( "seqn", "fs_enet", "age_enet",
                                       "fdas_enet", "hhs_enet", "pc1", "pc2" ) ] ) %>%
    mutate( test.set = ifelse( seqn %in% test$seqn, 1, 0 ) ) ) %>% # add indicator variable for membership in test set
   
  # ---------------------------------------------------------------------------------------------------------------------------------------------------------

### Save Analytic Data ###
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
saveRDS( ., "03-Data-Rodeo/01-analytic-data.rds" )
# ---------------------------------------------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------------------------------------------

