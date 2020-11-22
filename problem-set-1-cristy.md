---
title: "Problem Set 1"
author: "(Cristy Almonte )"
date: "due: 11am, Thursday, 12 September 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
```

We have designed these problem sets as a vehicle for learning. And so we want you to pursue completing them in a way that maximizes how much you learn. We assume that there is an optimum for each individual that is a mixture of independent effort, discussion with classmates, and seeking assistance from the TAs or Dr. Kim and myself. So while it is appropriate to seek input from your colleagues, it is not appropriate to copy solutions from each other. You will learn and be the most satisfied with your effort, if you challenge yourself to reach your own solutions. 

Although the recitations are not solely for help with homework, the TAs may spend part of each recitation helping students who may be stuck. You may also email us with questions, but be advised there is no guarantee we will be able to respond in the hours before it's due. And providing some detail in your email about what your problem is and what you have tried will increase the odds of a helpful reply.

> Since this is the first problem set, I will provide two worked out solutions (to Q1a and Q1b) to illustrate how you can fill in your answers in R Markdown. For questions Q1a and Q1b you should look at the R Markdown in the *.Rmd file and see how it's done; you don't need to answer them yourself.

Please view the rendered HTML version of this R Markdown document first. It is much easier to read. In the HTML version the questions are in bold, which in the R Markdown is marked by flanking double asterisks `**`. Be sure to address all sub-parts of each question. Each question has a dedicated space to put your answer.

#### Q1a (example)

**Given the vector of values below, estimate the mean of the vector. Use two approaches. One may use a single function, but also find another approach.**

```{r}
x <- c(8, 0, 11, 15, 5, 3, 9, 20, 0, 9, 5, 0, 18, 2, 1, 1, 44, 11, 7, 9)
# Provide your two solutions below.
```

**A1a:**

*(If a solution can be done in a single code block, I sometimes provide one where you should put your answer. But this is an example problem so I present a solution for you, below.)*

```{r, results="hold"}
# If this was a real homework problem you could put your solution here.
```

**A1a:**

(*Here I offer the solution)*

```{r, results="hold"}
# Two ways to compute a mean in R:
mean(x)
sum(x)/length(x)
```

*The above is a sufficient solution to this problem. But if you're curious, there are other valid answers. There are often multiple ways to do things in R. I show another one of them below. In the homeworks, don't worry if you're doing it the "right" way. Many ways are often right. The point of the homeworks is to learn.*

You could have also provided this solution:

```{r}
# We can think of the mean of N numbers like an expectation with each number having 
# probability 1/N
expectation <- 0
p <- 1/length(x)
for (xi in x) {
  expectation <- expectation + xi
}
print(expectation*p)
```


----

#### Q1b (example)

**In R, both the combine function, `c()`, and the list constructor, `list()` can be used to make ordered data structures of things. Consider the following four lines of code. What are the similarities and differences among these four statements in terms of the resulting data structure, what class the result is, and the way the data is organized? Demonstrate how you can learn about the results with code.**

```{r, eval=FALSE}
c(1, c(1,2))
c(1, list(1,2))
list(1, c(1,2))
list(1, list(1,2))
## returns "character"

```

*For some questions, the answer will be completely blank. In this case you can write your solution as a combination of regular text and code, as warranted.*

**A1b:**

*(put your solution here...actually, normally you would, but I provide a solution below)*

**A1b:**

We can learn about these objects by looking at how the R console formats them when it prints them out.

```{r}
c(1, c(1,2))
```

This produced a numeric vector of length 3.

We can also confirm this by considering it's class and length:

```{r}
length(c(1, c(1,2)))
class(c(1, c(1,2)))
```

If we execute the second example from the question, we get:

```{r}
c(1, list(1,2))
```

We can see that the result is a list, as evidenced by the double bracket indexes. This result is also length 3.

We can check the class as follows:

```{r}
class(c(1, list(1,2)))
```

The next example is:

```{r}
list(1, c(1,2))
```

This one is different in that the result has length 2:

```{r}
length(list(1, c(1,2)))
```

The second item of the list contains a numeric vector containing 1,2. It didn't get flattened out.

The last one is different still. The second item in the list is a list, not a numeric vector:

```{r}
list(1, list(1,2))
```

----

#### Q2

An important part of learning R will be having facility with the built-in help.

**Use the R help to look up the `signif()` method. Use it to print $\pi$ to 5 significant digits. How many digits total are there in the output? How many to the right of the decimal? How is `signif()` different from `round()`, `floor()`, and `ceiling()`. Show some examples and explain.** 
**A2:**

```{r}
signif(pi, digits=5)
round(pi, digits=5)
floor(pi)
ceiling(pi)
```

*(Discuss what you learn here and show code you used to figure out how these functions work.)*

signif() allows for the rounding of the  values to the specified number of significant digits. In this case we specified digits=5, which printed out the first 5 digits in pie. When we applied the round(), it rounded pi to a specified  number of digits (which again was 5 ) but it started rounding after the decimal point. When we applied the floor(), it rounded the largest integer to or smaller than n, which in this case was 3. The ceiling () did the complete opposite, it rounded to the next interger above 3 giving us a value of 4. 



----

#### Q3

In the theory lectures you saw what's called the maximum likelihood estimate (MLE) of the variance (don't worry about what MLE means yet, we'll get to it later):

$$
\hat{\sigma}^2 = \frac{1}{n}\sum_{i}^{n}(x_i - \bar{x})^2
$$

And in the applied lectures you saw an unbiased estimate of the variance:

$$
\hat{\sigma}^2 = \frac{1}{n-1}\sum_{i}^{n}(x_i - \bar{x})^2
$$

These are both reasonable ways to estimate variance.

Consider a sample from a normal distribution, which we can get using the `rnorm()` function:

```{r}
s <- rnorm(10, mean=20, sd=3)
print(s)

```

**Compute both the MLE and unbiased estimates of the variance for a sample without using the `sd()` or `var()` functions. Do the computations using only the `sum()` function and math operators. Then use the `var()` function in R and identify which estimator R is using when computing the variance. What is the true value of the variance?**

**A3:**

*(Show your computations here)*
MLE
```{r}
s <- rnorm(10, mean=20, sd=3)
a= (1)/(10) * sum((s-(sum(s)/10))^2)
print(var(s))
print(a)
```
Unbiased estimate of the variance

```{r}

s <- rnorm(10, mean=20, sd=3)
a= (1)/(10-1) * sum((s-(sum(s)/10))^2)

print(a)
sprintf("this is variance: %f", var(s))
```

The true value for variance is 10.29476 and the unbiased estimate of the variance is the estimator that R is using. 

#### Q4

Consider the code below that samples from a uniform distribution between 0 and 20:

```{r}
x <- runif(10, min=0, max=20)
print(x)
```

Given a random variable, $X$, we can perform a **transformation** that is called computing the **z-score** of that variable. We will be learning about transformations soon in the theory lectures, and this little exercise should get you ready for working with them in R. The z-score transformation subtracts the mean and divides the result by the standard deviation:

$$
Z = \frac{X - \mu_x}{\sigma_x}
$$

This has the effect of giving the transformed variable, $Z$, a mean of 0 and a standard deviation of 1. It is a standard way of making random variables comparable that differ in mean and variance, but otherwise come from similar distributions. It's particularly important for random variables that are somewhat normally distributed.

**Compute the z-score for the sample, `x`.**

Hint: Look in Applied Lesson 2 for the R functions you need to compute $\hat{\mu}$ and $\hat{\sigma}$.

**AQ4:**

```{r}
x <- runif(10, min=0, max=20)
print(x)
ux<- mean(x)
qx<- sd(x)
z= (x-ux)/(qx)
print(z)
```
----

#### Q5a

Given July was the warmest month on record [according Copernicus Climate Change Service](https://climate.copernicus.eu/another-exceptional-month-global-average-temperatures), I thought it would be interesting to have a look at historical global temperatures.

I have included in the `problem-set-1` folder a dataset of historical temperatures. These are global temperatures including both air and sea data. The file is named: `NOAA-temperatures.tsv`.

In order to load this data, you'll need to make sure your R session's working directory is set to the folder containing this data. You can see the current working directory using the `getwd()` function:

```{r}
getwd()
```


If you need to change this, you can either use the `setwd()` function or do so from the RStudio menu:

`Session > Set Working Directory > Choose Directory`

Once you're in the right folder, you can read this data into R as follows:

```{r}
tempData <- read.table(file="NOAA-temperatures.tsv", header=TRUE)
tail(tempData)
```

This is what it looks like:

```{r, echo=FALSE, fig.width=5, fig.height=3.5}
par(mar=c(4,4,1,1))
plot(tempData$year, tempData$tempIndex, cex=0.6, ylab="temperature index", xlab="year")
```

**Compute the median temperature for before the year of your birth and compute the median for after the year of your birth. Also compute the mean temperature for both periods. Try and use some other summaries to to characterize how the temperatures in these two periods were different.**

Here is a hint:

You can index a vector with a logical vector. For example, if we wanted to examine the temperature in years before 1900, we could do:

```{r}
before1900 <- tempData$year < 1990
tempData$tempIndex[before1900]
```

**A5a**

```{r, results="hold"}
# First create a variable containing the year of your birth, for example:
yob <- 1996
# Then create vectors for the temperatures before and after that year.
before1996 <- tempData$year < yob
#print(before1996)

beforeMean <- mean(tempData$tempIndex[before1996])
beforeMedian <- median(tempData$tempIndex[before1996])
print(beforeMean)
print(beforeMedian)

after1996 <- tempData$year > yob
##print(after1996)
afterMean <- mean(tempData$tempIndex[after1996])
afterMedian <-median(tempData$tempIndex[after1996])
print(afterMean)
print(afterMedian)
# Then use functions you've learned to compute summaries of the data.

sdbefore <- sd(tempData$tempIndex[before1996])
varbefore<- var(tempData$tempIndex[before1996])
print (sdbefore)
print (varbefore)
sdafter <- sd (tempData$tempIndex[after1996])
varafter <- var(tempData$tempIndex[after1996])
print (sdafter)
print(varbefore)

```

The before mean and median are very different from one another, while the after mean and median are fairly close together. 

----

#### Q5b

**Determine which years the temperature is above the 90% quantile. Do you see any patterns?**

This problem also uses the NOAA temperature from the above question. Please refer to **Applications Lesson 2** for how to compute quantiles. Your solution should be a subset of the vector, `tempData$year`.

Hint: A good way to solve this problem is to index the vector with a logical vector. Consider using the `>` operator.

**A5b:**

```{r}
temps <-tempData$tempIndex
perc90 <- quantile(temps, probs=.90)

print(perc90)

lessThanPerc90 <- temps > perc90

print(tempData$year[lessThanPerc90])

```

*(Discuss any patterns you see here.)*
 Since 2001, the temeperature is trending higher than the 90th percentile, with exception of a few years (2004, 2008,2011) . This grouping of data is  can be thought as anomalous to the entire data. 
----

#### Q6

Both the full range of the data and the 90% interquantile range describe the spread of the data. However, these summaries are computed from a sample, and so each time you draw a new sample, the summary will be slightly different.

**How does stability (or variation) of the range and the 90% interquantile range compare? Use the standard deviation as a tool for talking about how much estimates of the range vary vs the 90% interquantile range.**

From the lectures we learned that the **interquartile** range is the portion of a sample that ranges from the 25% quantile to the 75% quantile. Here we're interested in the range that spans the middle 90% of the data.

We will be using samples from a normal distribution, like in problem **Q3**.

We can draw 3 samples from a normal distribution, using default values for the mean and standard deviation argument by simply providing the number of samples we want to draw as the single argument:

```{r}
rnorm(3)
```

To answer this question, you can use the `replicate()` function, which allows us to execute the same code several times:

```{r}
replicate(10, rnorm(3))
```

Above we executed the code `rnorm(3)` 10 times, assembling the results in a matrix.

**A6:**

*We provide a partially worked out solution. We show how you can compute the range of many samples. Do the same for computing the 90% interquantile range the same number of times and show how the standard deviation of the estimates inform us about the robustness of the range vs the 90% interquantile range.*

```{r, results="hold"}
n <- 100
m <- 10


print(ranges <- replicate(m, range(rnorm(n))))
print(perc90 <- replicate(m, quantile(range(rnorm(n)), probs=c(.05,.95))))
# We can then see the standard deviation of the lower end of the range for each of the m samples.
cat("sd lower: ", sd(ranges[1,]))
# Also consider the upper end of the range
cat("sd upper: ", sd(ranges[2,]))
# Do the same for computing the interquantile range m times.
cat("sd quantile: ", sd(perc90[1, ]))
cat("sd quantile: ", sd(perc90[2, ]))

# ...
```
The standard deviation for the lower end  of  90% interquantile range is a lot different from the the sd of the lower end of the range. The sd of the upper end of 90% interquantile range is similar to the sd upper to end of the range. 

----

#### Q7

Consider the points visualized in this plot:

```{r, echo=-1, fig.width=3, fig.height=3}
par(mar=c(4,4,1,1))
x <- rnorm(30, mean=2)
y <- exp(rnorm(30, mean=1, sd=1))
plot(x,y, xlab="x", ylab="y", main="", pch=20)
```

**In the theory lecture a statistic was described that is equivalent to *balancing* your data. This is the point such that if each data point had the same weight, you could balance the two-dimensional plane containing the data perfectly at that point. Compute the (x,y) coordinates of this point and use the provided code to visualize it.**

**A7:**

*(edit the following code as indicated in the comments)*

```{r, echo=-1, fig.width=3, fig.height=3}
par(mar=c(4,4,1,1))
# Replace the following line with code to compute the "center of gravity" of these points. This is a tuple containg an x value and a y value.

centerOfGravity <- c(mean(x) , mean(y))
print (centerOfGravity)
# Once you provide a (x,y) tuple in the variable, centerOfGravity, the next two lines will plot it on the graph. You don't need to edit this part. The (NA,NA) tuple won't be visible, but once you compute your answer, you should see a red dot. 
{
  plot(x,y, xlab="x", ylab="y", main="", pch=20)
  points(centerOfGravity[1], centerOfGravity[2], col="red", pch=20)
}
```


----

#### Q8

Here is the code from Applications Lesson 2 that we used to load the transcript model data:

```{r}
d <- read.table(file="human_transcripts-chr1.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")
head(d)
```

**Compute the mean transcript length for transcripts that have just one exon. Do the same for transcripts that have exactly two exons, then the exons, and so on, up to 5 exons. What do you notice about these numbers? Is there a trend? Are there any outliers (i.e., values that seem different from the others)?**

**A8:**
*(place your solution and discussion here)*

The lowest mean length value observed was for two exons with a value of 993.2899 and the highest observed value was for one exon with a value of 1815.119. The one exon appears to be an outlier having a larger mean length when compared to the other values. 
```{r}
data<- read.table(file="human_transcripts-chr1.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")
head(d)

one <- data$exons == 1 
#print(one)
mean(data$length[one])

two <- data$exons == 2 
#print(two)
mean(data$length[two])


three <- data$exons == 3 
#print(three)
mean(data$length[three])

four <- data$exons == 4
#print(four)
mean(data$length[four])


five <- data$exons == 5
#print(five)
mean(data$length[five])

```

----

#### Q9

In this homework problem, we do an estimation game, kind of like how we estimated $\pi$ in Applied Lesson 1. You probably won't be able to look at the code in the Rmd file and fully understand it. But you don't need to. Instead, I am simply asking you add a line or two of code that can answer the question of how much area is in the smiley face.

```{r, echo=FALSE}
# You don't need to read or understand any of this code. DO NOT EDIT THIS BLOCK.
eyeTopY <- function(x) { sqrt(1-x^2)/8+0.8 }
eyeBottomY <- function(x) { -sqrt(1-x^2)/8+0.8 }
plotFace <- function() {
  par(mar=c(1,1,1,1))
  x <- seq(0, 1, length.out=100)
  plot(x, (sin(-x*pi)+1)*0.7, type="l", xlab="", ylab="", ylim=0:1, xlim=0:1, axes=FALSE)
  lines(x, ((sin(-x*pi)+1)/2+0.5)*0.7)
  x2 <- seq(-1, 1, length.out=100)
  lines((x2+2.6)/8, eyeTopY(x2))
  lines((x2+2.6)/8, eyeBottomY(x2))
  lines((8-(x2+2.6))/8, eyeTopY(x2))
  lines((8-(x2+2.6))/8, eyeBottomY(x2))
  box()
}
insideLeftEye <- function(x,y) {
  x > 0.2 & x < 0.45 & y < eyeTopY(x) & y > eyeBottomY(x)
}
insideRightEye <- function(x,y) {
  x > 0.55 & x < 0.8 & y < eyeTopY(x) & y > eyeBottomY(x)
}
insideSmile <- function(x,y) {
  y > (sin(-x*pi)+1)*0.7 & y < ((sin(-x*pi)+1)/2+0.5)*0.7
}
inside <- function(x,y) {
  insideSmile(x,y) | insideLeftEye(x,y) | insideRightEye(x,y)
}
```

I have created a function `plotFace()` that takes no arguments and does the following:

```{r, fig.width=3, fig.height=3}
plotFace()
```

I made the face a complicated shape so that you'd need to solve it by sampling (and not with, trigonometry, for example.). Do you think it's a bit creepy?

Our goal is to estimate the area inside the eyes and mouth combined.

We can draw a sample of (x,y) coordinates like this:

```{r}
x <- runif(50000, min=0, max=1)
y <- runif(50000, min=0, max=1)
```

The horizontal axis of the plot is in the range (0,1) and the vertical axis is also in the range (0,1).

I also provide a black-box function `inside()` that you can use to test whether a vector of x and y coordinates are inside the eyes and mouth.

We can use it like this:

```{r}
inside(head(x), head(y))
```

If we color the points, we can see which points are inside the eyes and mouth:

```{r, fig.width=3, fig.height=3}
{
  plotFace()
  colorIndex <- inside(x,y)+1
  points(x,y, pch=20, col=c("black", "blue")[colorIndex])
}
#print(colorIndex)
```

In the above code, we're converting a logical vector with `TRUE` and `FALSE` values to corresponding 1 and 0 values. We then add 1 so that what were `TRUE` values index the second item in the character vector of color names and what were `FALSE` values index the first item in the color names.

**Your job is to estimate the area of the region consisting of the eyes and mouth combined. Assume the area of the square is 1 square meter.**

**A9:**

*Edit the following code to arrive at and print out your solution*

```{r}
# First draw what you feel draw what you feel are a sufficient number of x and y
# samples. Change the value of n if necessary to get an accurate estimate. You might 
# need to try sampling a few times to get a sense of if your sample size is big enough.
n <- 50000
x <- runif(n, min=0, max=1)
y <- runif(n, min=0, max=1)

insidePoints = sum(inside(x, y))

#outsidePoints = sum(!inside(x, y))

#print(insidePoints)

meterPerSample = 1/50000

finalValue = meterPerSample * insidePoints

# Add a few lines of code here to use your x, y samples to estimate the area inside 
# the eyes and mouth. Remember to use the function, inside(x,y), as part of your 
# solution. Make sure your answer is scaled appropriately to the units.

# Print out your answer. Edit the following code, replacing NA with a variable 
# containing your estimate:
cat("The area of the face is about", finalValue, "square meters")
```

----

#### Q10 (bonus)

Recall our analysis of DNase hypersensitivity data from **Applications Lesson 2**. Here we load in another similar data set for a different gene, JUN:

```{r}
h <- read.table(file="dnase_hypersensitivity_near_JUN.tsv", header=TRUE, stringsAsFactors=FALSE)
head(h)
```

We can visualize this as follows:

**Use anything you've learned so far to identify where the peaks are.**

*This is a hard statistical problem, and it since this is a bonus question, it doesn't matter how far you get. The point of this question is to get you exploring R and playing with a real problem in an open-ended way. Try and have fun with it and learn some new things. It's easy to look at this plot and make a subjective judgment about where the peaks are. But what if you need to do this for 10,000s of genes? In that case we need a computational approach. How would you design a statistic? If you're not sure how to use the R you know to do this, design your statistic as precisely as you can using words.*

**A10:**

*(See what you can come up with, and put your efforts here)*

```{r}
library(pracma)

h <- read.table(file="dnase_hypersensitivity_near_JUN.tsv", header=TRUE, stringsAsFactors=FALSE)
tail(h)
#print(h)

print(length(h$position))

plot(h$position, h$intensity, type="l", xlab="Intensity", ylab="DNase intensity", 
     main="DNase hypersensitivity near the JUN gene")

peaks <- findpeaks(h$intensity,minpeakdistance = 100)

peaks = peaks[order(peaks[,2]),]

x = c(1:length(h$position))

masker <- function(xi) {
  return (xi %in% peaks[,2])
}

y = lapply(x, masker) 

y = unlist(y, use.names=FALSE)

#y

h$position[y]

plot(h, type='l')
points(h$position[y], peaks[,1])
  
```