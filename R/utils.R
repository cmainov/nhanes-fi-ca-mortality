###---------------------------------------------------
###   RESULT-GENERATING AND OTHER HELPER FUNCTIONS
###---------------------------------------------------


########## %notin% operator #########
`%notin%` <- Negate( `%in%` )

####################################################################################################
#################################### Quantile Cutting Function #####################################
####################################################################################################

quant_cut<-function(var,x,df){
  
  xvec<-vector() # initialize null vector to store
  
  for (i in 1:x){
    xvec[i]<-i/x
  }
  
  qs<-c(min(df[[var]],na.rm=T), quantile(df[[var]],xvec,na.rm=T))
  
  df[['new']]=x+1 # initialize variable
  
  for (i in 1:(x)){
    df[['new']]<-ifelse(df[[var]]<qs[i+1] & df[[var]]>=qs[i],
                        c(1:length(qs))[i],
                        ifelse(df[[var]]==qs[qs==max(qs)],x,df[['new']]))
  }
  
  return(df[['new']])
}


####################################################################################################
#################################### Trend Variable Function #######################################
####################################################################################################

trend_func<-function(rank.var,cont.var,df,trend.var,x){
  
  df[[trend.var]] = 1
  
  medians<-vector()
  
  for (i in 1:x){
    
    newdf<-df[df[[rank.var]]==i,]
    
    medians[i]<-median(newdf[[cont.var]],na.rm=T)
    
    df[[trend.var]]<-ifelse(df[[rank.var]]==i,medians[i],df[[trend.var]])
    
  }
  
  return(df)
}




####################################################################################################
####################################### Results Function ###########################################
####################################################################################################


res <- function( df, x, subs, cuts, id.col, covars, time, mort.ind ){
  require( tidyverse )
  require( glue )
  require( splines )
  
  these <- which( eval( parse( text = ( paste0( "df$", subs ) ) ) ) )
  
  # compute quantile rank and trend variable
  d.1 <- df[ these, ] %>%
    mutate( !!paste0( x, ".q" ) := as.factor( quant_cut( var = x, x = cuts, df = . ) ) ) %>%
    bind_rows( df[ -these, ], . ) %>%
    trend_func( rank.var = paste0( x, ".q" ), cont.var = x, df = ., trend.var = paste0( x, ".trend" ), x = 5 )
  
  des <- svydesign(id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                   nest = TRUE, survey.lonely.psu = "adjust", data = d.1)
  
  des <- subset( des, eval( parse( text = subs ) ) ) #inclusions
  
  
  ## Important function arguments ##
  
  # standard deviation to scale continuous predictor
  x.scale <- sd( df[ these, x ] )
  
  # knots at the quantiles
  kts <- paste0( levels( as.factor( d.1[[ paste0( x, ".trend" )]] ) ), collapse = ", " )
  
  # levels of cat variable
  cat.l <- levels( as.factor( d.1[[ paste0( x, ".q" ) ]] ) )
  
  ## Fit Models ##
  
  # quantile specification
  m.q <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( x, ".q" ), " + ", paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  sum.m.q <- summary( m.q )$coefficients %>% data.frame()
  ci.m.q <- confint( m.q )
  
  # trend test
  m.t <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( x, ".trend" ), " + ", paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  sum.m.t <- summary( m.t )$coefficients %>% data.frame()
  
  # linear specification
  m.l <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "I( ", x, "/", x.scale, ") + " ), paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  sum.m.l <- summary( m.l )$coefficients %>% data.frame()
  ci.m.l <- confint( m.l )
  
  # cubic polynomial
  m.c <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "I( ", x, "/", x.scale, ") + " ), 
                                    paste0( "I( (", x, "/", x.scale, ")^2) + " ), paste0( "I( (", x, "/", x.scale, ")^3) + " ),
                                    paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  sum.m.c <- summary( m.c )$coefficients %>% data.frame()
  
  # cubic spline 
  m.cs <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "bs(", x,", knots = c(", kts, "), degree = 3 ) + " ), 
                                     paste0( covars, collapse = " + ") ) ),
                    design = des )
  
  
  
  ## Generate Table ##
  
  
  # first table for odds ratios across quantiles
  res.mat.q <- matrix( ncol = 2 )
  res.mat.q[, 1 ] <- x
  res.mat.q[, 2 ] <- "1.00" # referent
  
  for( i in 2:cuts ){
    
    # model object row of interest
    cut.obj <- sum.m.q[ which( str_detect( rownames( sum.m.q ) , paste0( x,".q", cat.l[i]) ) ), ]
    
    # confidence interval object row of interest
    cut.ci <- exp( ci.m.q[ which( str_detect( rownames( sum.m.q ) , paste0( x,".q", cat.l[i]) ) ), ] )
    
    # put together row in table
    res.mat.q <- cbind( res.mat.q, paste0( round( exp( cut.obj[, "coef" ] ), 2 ),
                                           " (",
                                           paste0( round( cut.ci[1], 2 ), "-", round( cut.ci[2], 2 ) ),
                                           ")" ) )
    # asterisk on significant results
    if( cut.obj[, 6] < 0.05 & cut.obj[, 6] >= 0.01 ){
      res.mat.q[, i] <- paste0( res.mat.q[, i], "*" )
    }
    
    if( cut.obj[, 6] < 0.01 ){
      res.mat.q[, i] <- paste0( res.mat.q[, i], "**" )
    }
    
  }
  
  # add trend test p-value
  t.obj <- sum.m.t[ which( str_detect( rownames( sum.m.t ) , paste0( x,".trend" ) ) ), ]
  
  res.mat <- cbind( res.mat.q, round( t.obj[, 6], 2 ) )
  
  # asteriks
  if( t.obj[, 6] < 0.05 & t.obj[, 6] >= 0.01 ){
    res.mat[, ncol( res.mat ) ] <- paste0( res.mat[, ncol( res.mat ) ], "*" )
  }
  
  if( t.obj[, 6] < 0.01 ){
    res.mat[, ncol( res.mat ) ] <-  paste0( "< 0.01**" )
  }
  
  # add linear specification OR and 95% CI
  l.obj <- sum.m.l[ which( str_detect( rownames( sum.m.l ) , paste0( x ) ) ), ]
  l.ci <- exp( ci.m.l[ which( str_detect( rownames( sum.m.l ) , paste0( x ) ) ), ])
  
  res.mat <- cbind( res.mat, paste0( round( exp( l.obj[, "coef" ] ), 2 ),
                            " (",
                            paste0( round( l.ci[1], 2 ), "-", round( l.ci[2], 2 ) ),
                            ")" ) )
  
  # asteriks
  if( l.obj[, 6] < 0.05 & l.obj[, 6] >= 0.01 ){
    res.mat[, ncol( res.mat ) ] <- paste0( res.mat[, ncol( res.mat ) ], "*" )
  }
  
  if( l.obj[, 6] < 0.01 ){
    res.mat[, ncol( res.mat ) ] <- paste0( res.mat[,ncol( res.mat ) ], "**" )
  }
  
  # add cubic specification p-value
  c.obj <- sum.m.c[ which( str_detect( rownames( sum.m.c ) , "\\^3" ) ), ]
  
  res.mat <- cbind( res.mat, round( c.obj[, 6], 2 ) )
  
  # asteriks
  if( c.obj[, 6] < 0.05 & c.obj[, 6] >= 0.01 ){
    res.mat[, ncol( res.mat ) ] <- paste0( res.mat[, ncol( res.mat ) ], "*" )
  }
  
  if( c.obj[, 6] < 0.01 ){
    res.mat[, ncol( res.mat ) ] <-  paste0( "< 0.01**" )
  }
  
  # column names 
  res.frame <- data.frame( res.mat )
  colnames(res.frame) <- c( "index", paste0( "Q", 1:cuts ), "ptrend", "linear", "pcubic" )
  return( res.frame )
}

#res( df = d, x = "fs_enet", subs = "inc == 1", cuts = 5, id.col = "seqn", covars = covars.logit, time = "stime", mort.ind = "mortstat")
