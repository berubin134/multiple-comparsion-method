# multiple-comparsion-method
My outlier-dependent comparison method is a post-hoc comparison method which identifies significant differences between group means. It combines Tukey's Honestly Significant Difference (HSD) with a modified Bonferroni correction to maximize statistical power while in a more flexible way than either Tukey's HSD or Bonferroni.
The sample means of each group are used to calculate summary statistics (median and interquartile range). Then each sample mean is categorized as an extreme outlier if it falls more than 3 times the IQR from the median, an outlier if it falls between 1.5 and 3 times the IQR from the median, or a non-outlier if it is within 1.5 times the IQR from the median.
The extreme outliers are compared using a t-test with Bonferroni correction of $\frac{\alpha}{n}$. The outliers are compared using a t-test with modified Bonferroni correction of $\frac{\alpha}{\sqrt{n}}$. The non-outliers are compared using Tukey's HSD.

This method is designed to be used following a significant result from an ANOVA test. To use this method, call odc(means, J, MSE, alpha) where means is a vector of sample means from each group, J is the number of observations per group (the method assumes consistent number of observations across groups), MSE is the mean square error from the ANOVA test (Some of squares of the error divided by degrees of freedom of the error), and alpha is the significance level with a default value of 0.05.

Example:

ocd(means = c(14.5, 13.8, 13.3, 14.3, 13.1), J = 9, MSE = .088)

outputs the data frame:

   X1   X2    diff       significant
     
1  14.5 13.8  0.7        TRUE

2  14.5 13.3  1.2        TRUE

3  14.5 14.3  0.2       FALSE

4  14.5 13.1  1.4        TRUE

5  13.8 13.3  0.5        TRUE

6  13.8 14.3  0.5        TRUE

7  13.8 13.1  0.7        TRUE

8  13.3 14.3  1.0        TRUE

9  13.3 13.1  0.2       FALSE

10 14.3 13.1  1.2        TRUE


This method is less powerful than Tukey’s HSD, so Tukey’s would be preferred for cases where the assumptions of normality and equal variance are met, and there are not any extreme values. This method would be preferred over Tukey's when these assumptions are not met. It can maintain a low FWER particularly when there are outliers among the groups being compared. Additionally, it is less conservative and more powerful than the Bonferroni or modified Bonferroni correction.
