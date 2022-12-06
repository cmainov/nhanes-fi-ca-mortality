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
  
  these <- which( eval( parse( text = ( paste0( "df$", subs ) ) ) ) )
  
  # compute quantile rank and trend variable
  d.1 <- df[ these, ] %>%
    mutate( !!paste0( x, ".q" ) := as.factor( quant_cut( var = x, x = cuts, df = . ) ) ) %>%
    bind_rows( df[ -these, ], . ) %>%
    trend_func( rank.var = paste0( x, ".q" ), cont.var = x, df = ., trend.var = paste0( x, ".trend" ), x = 5 )
  
  des <- svydesign(id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra, 
                   nest = TRUE, survey.lonely.psu = "adjust", data = d.1)
  
  des <- subset( des, eval( parse( text = subs ) ) ) #inclusions
  
  # standard deviation to scale continuous predictor
  x.scale <- sd( df[ these, x ] )
  
  # knots at the quintiles
  kts <- paste0( levels( as.factor( d.1[[ paste0( x, ".trend" )]] ) ), collapse = ", " )
  
  
  # quantile specification
  m.q <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( x, ".q" ), " + ", paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  # trend test
  m.t <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( x, ".trend" ), " + ", paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  # linear specification
  m.l <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "I( ", x, "/", x.scale, ") + " ), paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  # cubic polynomial
  m.c <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "I( ", x, "/", x.scale, ") + " ), 
                                    paste0( "I( (", x, "/", x.scale, ")^2) + " ), paste0( "I( (", x, "/", x.scale, ")^3) + " ),
                                    paste0( covars, collapse = " + ") ) ),
                   design = des )
  
  # cubic spline 
  m.cs <- svycoxph( formula( paste0( "Surv(", time, ",",mort.ind," ) ~ ", paste0( "bs(", x,", knots = c(", kts, "), degree = 3 ) + " ), 
                                     paste0( covars, collapse = " + ") ) ),
                    design = des )
  
  
  return( m.q )
}

# res( df = d, x = "fs_enet", subs = "inc == 1", cuts = 5, id.col = "seqn", covars = covars.logit, time = "stime", mort.ind = "mortstat")
