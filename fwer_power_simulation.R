simulate_fwer <- function(I, J, sigma = 1, n_sim = 10000, alpha = 0.05) {
  #This function simulates data under the null hypothesis that all group means are the same
  #It then runs the outlier-dependent comparison method to determine how often it returns
  #a false positive
  false_rejections <- 0 #count the number of times the null hypothesis is rejected (false positives)
  for (sim in 1:n_sim) {
    # Generate data under H0 (all means = 0)
    data <- matrix(rnorm(I * J, mean = 0, sd = sigma), nrow = I, ncol = J)
    # Calculate summary statistics
    means <- rowMeans(data) #calculate the group means
    MSE <- mean(apply(data, 1, var)) #calculate the mean square error
    # Simplified; use pooled variance
    
    # Apply the method
    results <- odc(means, J, MSE, alpha)
    # Check if any rejection occurred
    if (any(results$significant)) { #check to see if the method returned any false positives
      false_rejections <- false_rejections + 1
    }
  }
  return(false_rejections / n_sim) #return the percentage of false rejections
}

simulate_power <- function(I, J, delta, sigma = 1, n_sim = 10000, alpha = 0.05) {
  #This function simulates data for groups under the alternative hypothesis that the mean of group 1 is significantly different
  #than the mean of group 4 (this difference is delta). It then runs both the outlier-dependent comparison and Tukey's HSD
  #to test whether these methods identify the difference.
  
  true_rejections_odc <- 0 #counts the number of true rejections of null hypothesis under my method
  true_rejections_tuk <- 0 #counts the number of true rejections of null hypothesis under tukey's method
  
  for (sim in 1:n_sim) {
    # Generate data under Ha:one group differs from the other
    x <- c()
    x <- c(x,rnorm(J*(I-2), mean = 0, sd = sigma))
    x <- c(x,rnorm(J*2, mean = delta, sd = sigma))
    
    data <- matrix(x, nrow = J, ncol = I)
    
    # Calculate summary statistics
    means <- c() #calculate group means
    for (i in 1:5){
      means <- c(means, mean(data[,i]))
    }
    
    #calculate mean squared error
    MSE <- mean(apply(data, 1, var)) # Simplified; use pooled variance
    
    # Apply modified methods to compare groups 1 and 4 which should be different
    true_rejections_odc <- true_rejections_odc + odc2(means,J,MSE)
    true_rejections_tuk <- true_rejections_tuk + tukey(means, J, MSE)
  }
  #return a vector with the percentage of correctly identified differences using outlier-dependent comparison
  #and the percentage of correctly identified differences using Tukey's HSD
  return(c(true_rejections_odc/n_sim, true_rejections_tuk/n_sim))
  
}

odc2 <- function(means, J, MSE, alpha = 0.05) {
#This is a modified version of the ODC method specifically to compare groups 1 and 4 and return whether they are significantly different.
  ID <- 0 # Returns 1 if a difference is identified, 0 if it is not
  I <- length(means) #number of groups
  
  iqr <- IQR(means) #calculate interquartile range
  med <- median(means) #calculate median
  
  #Take groups 1 and 4, the groups of interest
  m1 <- means[1]
  m4 <- means[4]
  n <- choose(length(means),2) #n is the total number of comparisons being made
  
  #If either m1 or m4 is an extreme outlier, compare using bonferroni
  if (m1 < med - 3 * iqr || m1 > med + 3 * iqr || m4 < med - 3 * iqr || m4 > med + 3 * iqr){
    t <- (m1 - m4) / sqrt (2 * MSE / J) #calculte test statistic
    p <- 2*(1-pt(t, I * (J-1))) #calculate p value
    #determine if it is significant with bonferroni correction
    if (p < alpha / n){
      ID <- 1
    }
  }
  
  #If either m1 or m4 is an outlier, compare using modified bonferroni
  else if (m1 < med - 1.5 * iqr || m1 > med + 1.5 * iqr || m4 < med - 1.5 * iqr || m4 > med + 1.5 * iqr){
    t <- (m1 - m4) / sqrt (2*MSE / J) #calculate test statistic
    p <- 2*(1-pt(t, I*(J-1)) )#calculate p value
    #determine if it is significant with modified bonferroni correction
    if (p < alpha / sqrt(n)){
      ID <- 1
    }
  }
  
  #If neither m1 nor m4 is an outlier, compare using Tukey
  else{
    Q <- qtukey(1-alpha, length(means), length(means/2) * (J-1)) #calculate Q
    HSD <- Q * sqrt(MSE / J) #calculate HSD
    if (abs(m1-m4) > HSD){ #determine whether the differences is significant under Tukey
      ID <- 1
    }
  }
  return(ID) #return the result
}

tukey <- function(means, J, MSE, alpha = 0.05) {
  #This is a modified version of the Tukey method specifically to compare groups 1 and 4 and return whether they are significantly different.
  ID <- 0 # Returns 1 if a difference is identified, 0 if it is not
  
  #Take groups 1 and 4, the groups of interest
  m1 <- means[1]
  m4 <- means[4]
  
  Q <- qtukey(1-alpha, length(means), length(means/2) * (J-1)) #calculate Q
  HSD <- Q * sqrt(MSE / J) #calculate HSD
  if (abs(m1-m4) > HSD){ #determine whether the differences is significant under Tukey
    ID <- 1
  }
  return(ID) #return the result
}


#Simulate Power
#set I and J
I <- 5
J <- 10
#create a vector of deltas to test
deltas <- c(.5,1,1.5,2)
power <- c() #create a vector to store the results
for (d in deltas){ #simulate the power for each delta
  result <- simulate_power(I,J,d)
  power <- c(power, result)
}
#Store the results in a data frame
power_results <- data.frame(delta,power)

#Simulate FWER
J <- c(3,5,7,5,5)
I <- c(10,10,10,5,20)
#create a vector to store the results
fwer <- c()
for (k in 1:length(I)){ #Simulate FWER for each I,J combination
  result <- simulate_fwer(I[k], J[k])
  fwer <- c(fwer, result)
}
#Store the results in a data frame
fwer_results <- data.frame(I,J,fwer)
