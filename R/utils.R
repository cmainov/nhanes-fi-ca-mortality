##---------------------------------------------------
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


res <- function( df, x, subs, cuts, id.col, covars, time, mort.ind, sample.name ){
  
  require( tidyverse )
  require( glue )
  require( splines )
  
  these <- which( eval( parse( text = ( paste0( "df$", subs, collapse = " & " ) ) ) ) )
  
  # compute quantile rank  and trend variable on subsample of interest and recombine
  d.1 <- df[ these, ] %>%
    mutate( !!paste0( x, ".q" ) := as.factor( quant_cut( var = x, x = cuts, df = . ) ) ) %>%
    bind_rows( df[ -these, ], . ) %>%
    trend_func( rank.var = paste0( x, ".q" ), cont.var = x, df = ., trend.var = paste0( x, ".trend" ), x = 5 )
  
  des <- svydesign(id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                   nest = TRUE, survey.lonely.psu = "adjust", data = d.1)
  
  des <- subset( des, eval( parse( text = paste0( subs, collapse = " & " ) ) ) ) #inclusions
  
  
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
  
  n.table <- nrow( m.l$model )
  
  # cubic polynomial
  m.c <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "I( ", x, "/", x.scale, ") + " ), 
                                    paste0( "I( (", x, "/", x.scale, ")^2) + " ), paste0( "I( (", x, "/", x.scale, ")^3) + " ),
                                    paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  sum.m.c <- summary( m.c )$coefficients %>% data.frame()
  
  # natural cubic spline 
  m.cs <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "ns(", x,", knots = c(", kts, ") ) + " ), 
                                     paste0( covars, collapse = " + ") ) ),
                    design = des )
  
  
  # Non-linearity Likelihood-Ratio test
  p.nl <- pchisq( abs( m.l$ll[2] -  m.cs$ll[2]),
                  df = m.l$degf.resid-m.cs$degf.resid, lower.tail = FALSE )
  
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
      res.mat.q[, (i+1)] <- paste0( res.mat.q[, (i+1)], "*" )
    }
    
    if( cut.obj[, 6] < 0.01 ){
      res.mat.q[, (i+1)] <- paste0( res.mat.q[, (i+1)], "**" )
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
  
  res.mat <- cbind( res.mat, round( p.nl, 2 ) )
  
  # asteriks
  if( p.nl < 0.05 & p.nl >= 0.01 ){
    res.mat[, ncol( res.mat ) ] <- paste0( res.mat[, ncol( res.mat ) ], "*" )
  }
  
  if( p.nl < 0.01 ){
    res.mat[, ncol( res.mat ) ] <-  paste0( "< 0.01**" )
  }
  
  # column names 
  res.frame <- data.frame( res.mat )
  colnames(res.frame) <- c( "index", paste0( "Q", 1:cuts ), "p.trend", "linear", "p.nonlinear" )
  
  
  # add n to table
  res.frame$n <- n.table
  res.frame <- res.frame %>%
    relocate( n, .before = Q1 )
  
  # add subsample name to table
  res.frame$sample <- sample.name
  res.frame <- res.frame %>%
    relocate( sample, .before = n )
  
  ## Significant digits ##
  
  # column indices for odds ratio/ci
  col.ind.q <- c( which( str_detect( colnames( res.frame ), "Q\\d" ) ), 
                  which( str_detect( colnames( res.frame ), "^linear$" ) ) )
  
  # odds ratios and confidence intervals
  for( i in col.ind.q ) {
    res.frame[,i] <- str_replace( res.frame[,i], '(\\(\\d)\\,', "\\1\\.00," ) # match open parenthesis followed by a digit and then a comma. Retain everything except the comma and add ".00,"
    res.frame[,i] <- str_replace( res.frame[,i], "(\\(\\d\\.\\d)\\,","\\10\\,") # match open parenthesis followed by a digit, period, digit, and then a comma. Retain everything except the comma and add ".0,"
    res.frame[,i] <- str_replace( res.frame[,i], "(\\(\\d\\.\\d)\\,", "\\10\\," ) # match open parenthesis followed by digit, period, digit and comma.Retain everything except the comma and add "0,"
    res.frame[,i] <- str_replace( res.frame[,i], "(\\d\\.\\d)\\s", "\\10 " ) # match digit, period, digit and space. Retain everything except space and add "0 "
    res.frame[,i] <- str_replace( res.frame[,i], "^(\\d)\\s\\(", "\\1\\.00 \\(" ) # match beginning of string followed by single digit, followed by a space and parenthesis. Reatin the single digit and and ".00 ()
    res.frame[,i] <- str_replace( res.frame[,i], "(\\d\\.\\d)\\)", "\\10\\)") # match digit, period, digit, close parenthesis. Retain everything except parenthesis and add "0)" to end
    res.frame[,i] <- str_replace( res.frame[,i], "(\\,\\s\\d)\\)", "\\1.00\\)") # match comma, space, digit, close parenthesis. Retain everything except parenthesis and add ".00)" to end
    res.frame[,i] <- str_replace( res.frame[,i], "(\\(\\d\\.\\d)\\-", "\\10\\-") # open parenthesis, digit, period, digit, hypen. Retain everything except hyphen and add "0)" to end
    res.frame[,i] <- str_replace( res.frame[,i], "(\\-\\d)\\)", "\\1.00\\)") # match hyphen, digit, close parenthesis. Retain everything except parenthesis and add ".00)" to end
    res.frame[,i] <- str_replace( res.frame[,i], "(\\(\\d)\\-", "\\1.00\\)") # match hyphen, digit, close parenthesis. Retain everything except parenthesis and add ".00)" to end
  }
  
  # column indices for columns containing p values 
  col.ind.p <- which( str_detect( colnames( res.frame ), "p\\." ) )
  
  for( i in col.ind.p ){
    res.frame[,i] <- str_replace( res.frame[,i], "(\\d\\.\\d)$", "\\10" ) # match digit, period, digit,end and add a 0 before the end
    res.frame[,i] <- str_replace( res.frame[,i], "^1$", "0.99" ) # round down probabilities = 1
  }
  
  return( list( frame = res.frame, q.obj = m.q,
                dat = des$variables ) )
}

# res( df = d, x = "fs_enet", subs = "inc == 1", cuts = 5, id.col = "seqn", covars = covars.logit, time = "stime", mort.ind = "mortstat")



####################################################################################################
################################# Table 1 (Categorical Variables) ##################################
####################################################################################################


epitab <- function(var,data.fr,des,table.var){
  
  attach(data.fr)
  
  typ<-paste0('~',var)
  sumcat=0
  for (i in 1:length(levels(factor(data.fr[[var]])))){
    sumcat<-sumcat+((svytable(formula(typ),design=des)[i]))
  }
  
  wtpct<-vector()
  
  for (i in 1:length(levels(factor(data.fr[[var]])))){
    wtpct[i]<-round((svytable(formula(typ),design=des)[i]/sumcat*100),digits=1)
  }
  wtpct
  wtpct<-c(' ',wtpct,' ')
  total<-vector()
  
  for (i in 1:length(levels(factor(data.fr[[var]])))){
    total[i]<-table(des[["variables"]][var])[i]
  }
  total<-c(' ',total,' ')
  
  levelnames<-c(table.var,levels(as.factor(data.fr[[var]])),' ')
  levelnames<-levelnames[!is.na(levelnames)==T]
  levelnames<-levelnames[!levelnames=='Missing/unknown']
  merged<-data.frame(cbind(levelnames,paste0(total,' (',wtpct,')')))
  
  colnames(merged)<-c('levelnames','mn')
  merged
  detach(data.fr)
  return(merged)
  
}



####################################################################################################
################################## Table 1 (Continuous Variables) ##################################
####################################################################################################

epitab.means <- function(cont.var, des, table.var, dig){ 
  # dig is the number of digits to round to
  mn<-paste0(round(svymean(as.formula(paste0('~',cont.var)),design = des,na.rm=T)[1],digits=dig),
             ' (',round(sqrt(svyvar(as.formula(paste0('~',cont.var)),design = des,na.rm=T))[1],digits=dig),')')
  
  ms2<-data.frame(c('',table.var,''),c('',mn,''))
  
  colnames(ms2)<-c('levelnames','mn')
  
  return(ms2)
}
