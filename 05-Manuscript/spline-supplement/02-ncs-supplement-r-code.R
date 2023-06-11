# generate some toy data
set.seed( 23 ) # set seed for reproducibility 
x <- rnorm( 20, mean = 35, sd = 5 )
y <- rnorm( 20, mean = 78, sd = 10 )

# specify a natural cubic spline with 2 interior knots (K=4)
#(generate the basis matrix)
k2 <- ns( x, df = 3 )

head( k2, 5 ) # print first 5 rows of the basis matrix

# 1         2           3
# [1,]  0.27013521 0.4949742 -0.31550699
# [2,] -0.07878135 0.4987137 -0.33410946
# [3,]  0.53390462 0.3255878  0.03408727
# [4,] -0.14194707 0.4300676  0.71187952
# [5,]  0.50852945 0.3222180  0.09029069


# linear regression model with `ns`
lm( y ~ ns( x, df = 3 ) )

# Call:
#   lm(formula = y ~ ns(x, df = 3))
# 
# Coefficients:
#   (Intercept)  ns(x, df = 3)1  ns(x, df = 3)2  ns(x, df = 3)3  
# 80.999           6.804         -13.666          -2.166  

# alternative code for specifying the same matrix
k2.v2 <- ns( x, 
             knots = quantile( x, c( 0.33, 0.66 ) ), # interior knots (1st and 2nd tertiles)
             Boundary.knots = quantile( x, c( 0, 1 ) ) ) # boundary knots (min and max)
head( k2.v2, 5 ) # print first 5 rows of the basis matrix

# 1         2           3
# [1,]  0.27454544 0.4970633 -0.31053853
# [2,] -0.08486486 0.5069896 -0.33380781
# [3,]  0.52668014 0.3284054  0.04392163
# [4,] -0.14508287 0.4247300  0.72035287
# [5,]  0.49944729 0.3253491  0.10026887
