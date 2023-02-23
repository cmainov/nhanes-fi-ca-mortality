library( survey )
library( tidyverse )


source( "R/utils.R" )



dat <- readRDS( "03-Data-Rodeo/01-analytic-data.rds") %>%
  dplyr::filter( is.na( wtdr18yr ) == F ) %>% # subset to those having non -missing weights 
    dplyr::mutate( hhsize.bin = ifelse( hhsize >= 5, 1,
                                        ifelse( hhsize < 5, 0, NA ) ) ) # dichotomize household size column before generating table
  
  # designs
  nhc <- svydesign( id = ~sdmvpsu, weights = ~wtdr18yr, strata = ~sdmvstra,
                    nest = TRUE, survey.lonely.psu = "adjust", data = dat )
  
  # create three subsets, food insecure with cancer ( fiw ), food secure with cancer ( fsw ), and combined 
  # food secure and insecure cancer survivors
  fiw <- subset( nhc, binfoodsechh == "Low" & inc == 1 )
  fsw <- subset( nhc, binfoodsechh == "High" & inc == 1 )
  gen <- subset( nhc, inc == 1 )

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
  bmi <- epitab.means( cont.var = "bmxbmi", des = design, table.var = "BMI" )
  metmins <- epitab.means( cont.var = "weekmetmin", des = design, table.var = "MET Minutes" )
  calor <- epitab.means( cont.var = "kcal", des = design, table.var = "Calories" )
  age <- epitab.means( cont.var = "age", des = design, table.var = "Age" )
  cci <- epitab.means( cont.var = "cci_score", des = design, table.var = "CCI" )
  site <- epitab( var = "primarycagroup", data.fr = df, des = design, table.var = "Cancer Site" )
  snap <- epitab( var = "foodasstpnowic", data.fr = df, des = design, table.var = "SNAP Assistance" )
  time <- epitab( var = "timecafactor", data.fr = df, des = design, table.var = "Years Since Diagnosis" )
  insur <- epitab( var = "ins.status", data.fr = df, des = design, table.var = "Health Insurance Status" )
  
  table1 <- rbind( age, sex, race, education, income, insur, hhsize, bmi, metmins, calor, cci, snap, site, time, smokstat, alcuse )
  
  return( table1 )
  
}



# generate table columns for each of the subsets described above
fiw.tab <- cafs_table1( design = fiw, df = dat )
fsw.tab <- cafs_table1( design = fsw, df = dat )
gen.tab <- cafs_table1( design = gen, df = dat )

# merge columns into table
final.tab <- cbind( gen.tab, fiw.tab, fsw.tab )


## add column for p value for t tests and chi square test ##

# vector of strings containing elements that singal to a given row in the table
chi  <- c( "Smoking", "Alcohol", "Gender", "Income", "Size", "Education", "Race",
           "Years", "Site", "SNAP", "Insurance" )
these  <- c( "smokstat", "alc_cat", "gender", "fipr", "hhsize.bin", "education_bin", "race", 
             "timecafactor", "primarycagroup", "foodasstpnowic", "ins.status" )

# chi square
for ( i in 1:length( chi ) ){
  
  final.tab[ which( str_detect( final.tab[ , 1 ], chi[i] ) ), 7 ] <- ifelse( svychisq( as.formula( paste0( "~binfoodsechh + ", these[i] ) ), design = gen )$p.value < 0.01, "< 0.01",
                                                                             round( svychisq( as.formula( paste0( "~binfoodsechh + ", these[i] ) ), design = gen )$p.value, digits = 2 ) )
  
} # error w/ cancer site

tt  <- c( "Age", "BMI", "HHSize", "MET", "Calories", "CCI" )
these.b  <- c( "age", "bmxbmi", "hhsize", "weekmetmin", "kcal", "cci_score" ) 


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
write.table( final.tab, "04-Tables-Figures/tables/table-1.txt", sep = ", ", row.names = FALSE )
