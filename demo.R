#Input Data
x <- c(-0.531, -1.266, -0.994, 0.39, 0.204, -1.357, 0.1422, -1.392, 1.003, 0.538,
       -1.270, 0.106, -1.064, -1.022, -0.487, 1.721, -1.265, 1.893, 1.713, -1.172,
       -0.906, 1.470, 1.769, -0.220, 1.044, 0.914, 1.548, -0.084, -0.546, 0.589,
       2.205, 0.901, 3.224, 0.309, 1.275, 1.315, 3.097, 2.265, 1.663, 3.562,
       2.813, 1.882, 1.47, 1.546, 2.684, 2.627, 3.696, 0.415, 3.692, 2.888)
data <- matrix(x, nrow = 10, ncol = 5)
#Calculate sample means
means <- c()
for (i in 1:5){
  means <- c(means, mean(data[,i]))
}

#Shapiro test reveals that this data is not normal, the assumption for Tukey are not met
shapiro.test(data[,2])

#Calculate mean square error
MSE <- mean(apply(data, 1, var)) # Simplified; use pooled variance

#Run outlier-dependent comparison
odc(means,10,MSE)
