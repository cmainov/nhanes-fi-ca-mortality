library(fmsb)

corr.matrix[,1]

fg.only <- c( "Processed\n Meat", "Meat", "Poultry", "Fish,\n Hi n-3", "Fish,\n Lo n-3", 
              "Eggs", "Solid\n Fats", "Oils", "Milk", "Yogurt", "Cheese", 
              "Alcohol", "Fruit\n Other", "Citrus,\nMelons,\nBerries",  
              "Green\n Leafy\n Veg.", "Dark Ylw.\n Veg.", "Tomatoes", "Other\n Veg.", "Potatoes", 
              "Other\n Starchy\n Veg.", "Legumes", "Soy", "Refined\n Grain", 
              "Whole\n Grain", "Nuts", "Added\n Sugars" )


fg.only.lo <- c( "processedmts", "meat", "poultry", "fish_hi", "fish_lo", 
                         "eggs", "solidfats", "oils", "milk", "yogurt", "cheese", 
                         "alcohol", "fruitother", "f_citmelber", 
                         "greenleafy", "darkylveg", "tomatoes", "otherveg", "potatoes", 
                         "otherstarchyveg", "legumes", "soy", "refinedgrain", 
                         "wholegrain", "nuts", "addedsugars" )

rp.matrix <- corr.matrix[ fg.only.lo, ]


## radar plot for food insecurity pattern ##

# prepare frame for plotting (first row is max value for each variable and second row is minimum value; third row is value you want to plot)
rp.all <- setNames( data.frame( bind_rows( setNames( matrix( rep( 1, length( fg.only ) ), nrow = 1 ), fg.only ),
                                          setNames( matrix( rep( -1, length( fg.only ) ), nrow = 1 ), fg.only ),
                                          setNames( matrix( rp.matrix[, 1 ], nrow = 1 ), fg.only ),
                                          setNames( matrix( rp.matrix[, 3 ], nrow = 1 ), fg.only ),
                                          setNames( matrix( rp.matrix[, 5 ], nrow = 1 ), fg.only ),
                                          setNames( matrix( rp.matrix[, 6 ], nrow = 1 ), fg.only ), ) ),
                   fg.only )




rownames( rp.all ) <- c( "Max", "Min", "Food Insecurity", "Food Assistance\n (SNAP)", "Prudent #1", "Prudent #2" )
colors.line <- c(rgb(0.529,0.808,0.98,0.9), rgb(0.55,0.10,0.10,0.9) , rgb(0.7,0.5,0.1,0.9), rgb(0,0,0.502,0.9) )
colors.fill <- c(rgb(0.529,0.808,0.98,0.2), rgb(0.55,0.10,0.10,0.2) , rgb(0.7,0.5,0.1,0.2), rgb(0,0,0.502,0.2) )

# plotting parameters
op <- par(family = "Avenir", font = 2, col = "gray38",
          mar = c(2, 2, 2, 2), cex = 1.26) # font, font type, font color, and reduce margins

## plot with above parameters ##

# four select indices overlapping

radar.all <- radarchart( rp.all, axistype=1, axislabcol = 'black',
            caxislabels = c( "-1.0", "-0.5", "0", "0.5", "1.0"),
            plwd=2.5, plty=1, pcol = colors.line,
            pfcol=colors.fill,
            cglty = 3, title = "All Patterns" )

# fi pattern alone
radar.fi <- radarchart( rp.all[ c( 1,2,3), ], axistype=1, axislabcol = 'black',
            caxislabels = c( "-1.0", "-0.5", "0", "0.5", "1.0"),
            plwd=2.5, plty=1, pcol = colors.line[1],
            pfcol=colors.fill[1],
            cglty = 3, title = "Food Insecurity Pattern" )

# SNAP pattern alone
radar.snap <- radarchart( rp.all[ c( 1,2,4), ], axistype=1, axislabcol = 'black',
            caxislabels = c( "-1.0", "-0.5", "0", "0.5", "1.0"),
            plwd=2.5, plty=1, pcol = colors.line[2],
            pfcol=colors.fill[2],
            cglty = 3, title = "Food Assistance (SNAP) Pattern" )

# prudent #1 pattern alone
radar.p1 <- radarchart( rp.all[ c( 1,2,5), ], axistype=1, axislabcol = 'black',
            caxislabels = c( "-1.0", "-0.5", "0", "0.5", "1.0"),
            plwd=2.5, plty=1, pcol = colors.line[3],
            pfcol=colors.fill[3],
            cglty = 3, title = "Prudent Pattern #1")

# prudent #2 pattern alone
radar.p2 <- radarchart( rp.all[ c( 1,2,6), ], axistype=1, axislabcol = 'black',
            caxislabels = c( "-1.0", "-0.5", "0", "0.5", "1.0"),
            plwd=2.5, plty=1, pcol = colors.line[4],
            pfcol=colors.fill[4],
            cglty = 3, title = "Prudent Pattern #2" )


plot_grid(radar.p1, radar.p2)

## ggradar
library(ggradar)
rp.all.gg<-rp.all[3:nrow(rp.all),] %>%
  mutate( across( where( is.numeric ), ~ (.x + 1 ) / 2 ) ) %>% # map -1-1 scale to 0-1 scale
  data.frame()

colnames(rp.all.gg) <- c(fg.only.lo)

# need to add the group names as the first column of the dataset you feed to `ggradar` since it looks for those in column #1
# and begins looking at numeric values in column #2 and onward
rp.all.gg <- rp.all.gg %>%
  mutate( group = rownames( . ) ) %>%
  relocate( group, .before = processedmts )
  
# line colors
colors.line <- c(rgb(0.529,0.808,0.98,0.9), rgb(0,0,0.502,0.9), rgb(0.55,0.10,0.10,0.9) , "goldenrod2" )

# plot
(g.1 <- ggradar( rp.all.gg, values.radar = c("-1.0", "0.0", "1.0"),
         font.radar = "Avenir",
         axis.labels = fg.only,
         legend.position = c(0.93,0.85),
         legend.title = "Dietary Pattern",
         axis.label.size = 6,
         x.centre.range = 1.31 ,
         group.colours = colors.line ))

# change axis label text to grey
g.1$layers[[5]]$aes_params <- c(g.1$layers[[5]]$aes_params, colour = "grey40")
sapply(rp.all.gg, function(x) class(x))


mtcars %>%
  add_rownames( var = "group" ) %>%
  mutate_each(funs(rescale), -group) %>%
  tail(4) %>% select(1:10) -> mtcars_radar

ggradar(mtcars_radar) 
View(mtcars_radar)
