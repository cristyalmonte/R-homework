---
title: "Problem Set 6"
author: "Cristy Almonte"
date: "due: 5pm, Friday, 15 November 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
options(width=100)
```

---

*Make sure you have read through the full Theory Lecture Notes 10 before working on this problem set.*

#### Tutorial of regression in R

Here I preview some features of R that will help with linear regression in R that we haven't had time to present in class.

In this tutorial we'll use simulated data for $x$ from a normal distribution:

```{r}
x <- runif(20, min=1, max=5)
print(x)
```

We'll design $y$ to be a linear function of $x$ plus some noise from $N(0,1)$:

```{r}
y <- 3.0 + 2.5 * x + rnorm(length(x))
```

Because this is simulated data, we know the true values of the coefficients, $b_0 = 3.0$ and $b_1 = 2.5$.

We can plot our simulated data as follows:

```{r, echo=-1, fig.width=6, fig.height=4}
par(mar=c(4,4,1,1))
plot(x, y)
```

R has a special type of object that is used to represent model formulas. These objects are created using a special operator, `~`, which is the tilde character:

```{r}
response ~ predictor
```

In the above expression, the terms `response` and `predictor` are not referring (yet) to variables that exist in our R session. They are simply names that R will later use to look for the data in either a provided data frame or in as variables in the R session, once we ask it to do something with the `formula`.

We can see that this is of class `formula`:

```{r}
class(response ~ predictor)
```

We use this type of formula when specifying to R how we want to do linear regression.

For example, we can fit a linear model that predicts our `y` dependent variable using our `x` independent variable using the `lm()` function and our formula `y ~ x` as follows:

```{r} 
m <- lm(y ~ x)
print(m)
```

Our true intercept was `3.0` and our true slope was `2.5`. We can see that our fit model didn't quite recover these values, but it got reasonably close.

The object `m` is a special object class:

```{r}
class(m)
```

We can see the what properties the model has by looking at it's `names()`:

```{r}
names(m)
```

If we want the coefficients, we can get those as follows:

```{r}
m$coefficients
```

There also happens to be a `coef()` helper function that let's us extract this information:

```{r}
coef(m)
```

The same is true for the residuals:

```{r}
m$residuals
resid(m)
```


We can get a lot more details about the fit of the linear model by asking R for a summary of that fit model object:

```{r}
summary(m)
```

The result of `summary()` is also an object we can get data out of:

```{r}
names(summary(m))
summary(m)$r.squared
```

If we want to plot our regression line on the scatter plot of our data, we can do so as follows:

```{r, echo=-1, fig.width=6, fig.height=4}
par(mar=c(4,4,1,1))
plot(x, y)
abline(a=m$coefficients[1], b=m$coefficients[2], col="gray", lty=2)
```


#### Data preparation

For the next few questions we'll examine the relationship between lifespan and telomeres for a few different species. 

We'll use some data from this paper:

> [Telomere shortening rate predicts species life span](https://www.pnas.org/content/116/30/15122)<br>
Kurt Whittemore, Elsa Vera, Eva Martínez-Nevado, Carola Sanpera, and Maria A. Blasco<br>
PNAS July 23, 2019 116 (30) 15122-15127

Telomeres are regions of repetitive DNA at the ends of chromosomes that serve as a buffer to protect the information-containing portion of the chromosomes. In humans the sequence `TTAGGG` is repeated thousands of times at the end of each chromosome. During cellular replication, the very end of the chromosome can't be replicated and so each division the telomeres get shorter. A special enzyme, telomerase can lengthen telomeres, but this generally doesn't happen between ordinary somatic divisions. Telomeres have been implicated as important to many aspects of aging and cancer.

We can load the data into R as follows:

```{r}
telo <- read.table(file="telomeres.tsv", header=TRUE, stringsAsFactors=FALSE)
n <- nrow(telo)
print(telo)
```

Here I plot a `dotchart()` for each of the data frame columns (except the first):

```{r, fig.width=9, fig.height=2.2}
par(mar=c(4,1,1,1), mgp=c(2.2, 0.8, 0), mfcol=c(1,ncol(telo)-1))
for (i in 2:ncol(telo)) {
  # Reverse the row order so the order matches the data frame row order
  dotchart(telo[n:1,i], xlab=names(telo)[i], pch=20, cex.lab=1.5)
}
```

You'll notice a few of the variables have one value that is much higher than the rest. This can pose a problem for linear regression if the the variable is not on the natural scale that is relevant to the thing being predicted. When variables vary over several orders of magnitude we can often mitigate this problem by using the log transform. 

Here we built a new data frame comprising log transformed variables:

```{r}
d <- cbind(species=telo[,1], log10(telo[,-1]))
print(d)
```

We can draw the same dotcharts for the log-transformed data:

```{r, fig.width=9, fig.height=2.2}
par(mar=c(4,1,1,1), mgp=c(2.2, 0.8, 0), mfcol=c(1,ncol(d)-1))
for (i in 2:ncol(d)) {
  dotchart(d[n:1,i], xlab=names(d)[i], pch=20, cex.lab=1.5)
}
```

We can also look at how the variables correlate with each other using the `pairs()` plotting function:

```{r, fig.width=7, fig.height=7}
# We exclude the first column by using a negative index
pairs(d[,-1])
```

#### Q1

Simple linear regression with one predictor follows the following model:

$$
y = b_0 + b_1x + \epsilon
$$

Where each $x$ is the single predictor variable, $\epsilon \sim N(0,1)$ is the noise and $b_0, b_1$ are the intercept and slope parameters that we fit (coefficients). 

We can estimate the parameters using the equation Eq 10.2.3 from the theory notes:

$$
\begin{align}
\hat{b}_1 =& \ \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sum_{i=1}^n (x_i - \bar{x})^2} \\
\hat{b}_0 =& \ \bar{y} - \hat{b}_1\bar{x}
\end{align}
$$

With these data we're most interested whether average lifespan can be predicted by either telomere length or how fast the telomeres shorten over time. In this question we'll first consider predicting average lifespan using the telomere shortening rate, as these appear highly correlated in the above `pairs()` plot.

**Compute the coefficients $b_0$ and $b_1$ using the formulas above and when using telomere shortening rate to predict average lifespan. Don't use the `lm()` function yet.**

**A1:**

```{r}
##shortening rate to predict lifespan

c <- sum(d$shortening-mean(d$shortening))* (d$avgLifespan-mean(d$avgLifespan))
g  <- sum((d$shortening-mean(d$shortening))^2)
b1 <- c/g
b1
b0 <- mean(d$avgLifespan)- b1*mean(d$shortening)
b0

```


---

#### Q2

Recall that the predicted value of a linear model is:

$$
\hat{y} = \hat{b}_0 + \hat{b}_1 x 
$$

And the residual is defined as how much our prediction is off from the observed value:

$$
\delta = y - \hat{y}
$$

**Plot a scatter plot of $\hat{y}$ vs $y$ for the linear model predicting average lifespan using telomere shortening rates. Draw your plot so it is square and both the `xlim` and `ylim` are the same. Add the diagonal line (using the `abline()` function) with slope 1 and intercept 0 for reference so we can get a sense of how good the predictions are (if they were perfect all the points would lie on the diagonal). Also compute the residuals. What is the geometric interpretation of the residuals on the plot of $\hat{y}$ vs $y$?**
```{r}

```

**A2:**

*(add the code to compute the predictions and residuals and draw your plot here)*
The residuals helps us understand our regression model.From this plot, we can tell that the residuals are 


```{r}
yhat = b0 + b1*d$shortening
yhat
y <- d$avgLifespan
s= d$avgLifespan - yhat
s

fit <- lm(shortening ~ avgLifespan, data=d)
fit$coefficients
plot(yhat ~ y, data=d) + abline(0, 1, col="red", lty=2) +points(s, pch=20, col="orange")

```


---

#### Q3

Recall that the total sum of squares is:

$$
SS_{tot} = \sum_{i=1}^n(y_i - \bar{y})^2
$$

And that the residual sum of squares is:

$$
SS_{res} = \sum_{i=1}^n(y_i - f(x_i))^2
$$

Where $f(x)$ is the linear function we fit to predict data using our variable, $x$. That is to say, $\hat{y}_i = f(x_i)$.

The proportion of the variance explained by the model is:

$$
R^2 = \frac{SS_{tot} - SS_{res}}{SS_{tot}}
$$

**Compute the proportion of variance explained by the model.**

**A3:**

```{r}
SStot <- sum((d$shortening - mean(d$shortening))^2)
SSres <- sum(fit$residuals^2)
R2 <- (SStot - SSres) / SStot
R2
```

---

#### Q4

As we saw before with the simulated data, R can do these calculations for us. 

For example, we can fit a linear model that predicts the `avgLifespan` column in the data frame `d` with the `shortening` column as follows:

```{r}
# Here our formula is referring to columns in our data frame and so we need to pass lm() the data
# frame. Without a data frame lm() looks in the current R session environment.
mod1 <- lm(avgLifespan ~ shortening, data=d)
summary(mod1)
```

Here I run `summary()` again, but here I print it in a way that we get line numbers, which will be helpful for answering this question.

```{r,}
# I'm using the capture.output() function so we get line numbers
out <- capture.output(summary(mod1))
cat(paste(1:length(out), "#", out, collapse="\n"), "\n")
```

**Give the line number where you find the pieces of information requested in the answer section of this question. Check that each of these numbers matches what you computed in Q1 and Q3.**

**A4:**

*(For each item, replace the character `N` with the line number containing the requested piece of information from the summary listing above.)*

1. estimate of the intercept: 3.24075
2. slope coefficient: -0.73811
3. proportion of the variance explained: 0.9342

---

#### Q5

In this question we'll be performing **multiple linear regression**. That is, regression with multiple predictors.

Recall that multiple linear regression has this form:

$$
y = b_0 + b_1x_1 + b_2x_2 + \cdots + b_kx_k + \epsilon
$$

Where each $x_i$ are separate (ideally independent) variables. $\epsilon \sim N(0,1)$ is the noise and $b_0, \dots, b_k$ are the fit parameters (coefficients). 

We can also write this equation using vector notation. We can think of the vector: $\mathbf{b} = (b_0, b_1, \dots, b_k)$ multiplied by the row vector of predictors with an extra 1 added at the front for the intercept: $\mathbf{x}^T = (1, x_1, x_2, \dots, x_k)$:

$$
y = \mathbf{x}^T\mathbf{b}
$$

Or if we have a feature matrix, $\mathbf{X}$, of data with one row per object and the columns are features, we have:

$$
\mathbf{y} = \mathbf{X}\mathbf{b}
$$

Like the vector, $\mathbf{x}$, where we added a 1 in front for the intercept term, we can add a column of 1's to $\mathbf{X}$.

To indicate that our model has multiple predictors, we can simply give the formula more terms. Here is an example formula with three predictors:

```{r}
y ~ x1 + x2 + x3
```

**Perform multiple linear regression using telomere length, shortening, heart rate, and body weight as predictors for predicting average lifespan. Does this model explain a larger percent of the variance than the simple linear regression (with one predictor).**

**A5:**

This model does explain a larger percent of the variance, when compared to the simple linear regression.

```{r}
#head(d)
fit2 <- lm( avgLifespan ~ telomereLen + shortening+ heartRate + bodyWeight, data=d)
summary(fit2)
```
 
---

#### Q6

In Q5 you performed multiple linear regression.

You'll notice in the output to `summary()` that there is a quantity, `Adjusted R-squared`.

Adjusted $R^2$ is a way to account for the fact that you've added additional terms to your model. Whenever you add more parameters the model fit is going to improve, even if you're adding something that is totally independent of the variable you are trying to predict. This is simply because the least squares fitting procedure is going to assign a coefficient that allows the model to use the noise present in that predictor to slightly reduce the sum of squares and thereby increase the variance explained.

[Adjusted R-squared](https://en.wikipedia.org/wiki/Coefficient_of_determination#Adjusted_R2) is calculated as follows:

$$
\bar{R}^2 = 1 - (1-R^2)\frac{n-1}{n-p-1}
$$

Where $p$ is the number of predictors in the model. For simple linear regression (with just one predictor), $p=1$. For the case when we use all four predictors we have $p=4$. This correction to $R^2$ accounts for the extra parameters we're using to fit the data.

**Compute the model predictions using the coefficients from the model you fit in Q5. Then compute $R^2$ for the multiple linear regression model with the 4 predictors using the equation given in Q3. Then compute the adjusted $R^2$ for both the simple linear regression $R^2$ in Q3 and the multiple linear regression $R^2$ you computed for this question. Check that both adjust $R^2$ values match what you found using `summary()`.**

*Hint: You can get the coefficients out of the `lm` object returned by `lm()` by using the `coef()` method. Don't worry about calculating the coefficients the long way.*
```{r}


##Equation 3

SStot <- sum((d$shortening - mean(d$shortening))^2)
SStot2 <- sum((d$telomereLen - mean(d$telomereLen))^2)
SStot3 <- sum((d$heartRate - mean(d$heartRate))^2)
SStot4 <- sum((d$bodyWeight - mean(d$bodyWeight))^2)

SStotal <- SStot + SStot2 + SStot3 + SStot4
SSres <- sum(fit$residuals^2)

SSre2 <- sum(fit2$residuals^2)
SSre2

R2b <- (SStotal - SSre2) / SStotal
R2b
R2a <-(SStot - SSres) / SStot
R2a

n <- length(d$avgLifespan)
p <- 4
R2bar.multiple <- 1- (1-R2b)* ((n-1)/(n-p-1))
R2bar.multiple


R2bar.simple  <- 1- (1-R2a)* ((n-1)/(n-p-1))
R2bar.simple 

summary(fit2)
summary(fit)

```




---

#### Q7

In this question we'll explore performing **polynomial regression**.

Recall that polynomial regression of order $k$ fits a model of the form:

$$
y = b_0 + b_1x + b_2x^2 + \cdots + b_kx^k + \epsilon
$$

Where $\epsilon \sim N(0,1)$ is the noise and $b_0, \dots, b_k$ are the fit parameters. In the above formula, there is a single predictor variable, $x$, but it is raised to multiple powers.

There is a data set provided with the `MASS` library that comes with R that is from an old study relating how rat muscle tissue contracts when dipped in various concentrations of $\textrm{CaCl}_2$.

Strips of muscle tissue was dissected from rat hearts. Each strip was dipped into one of 6 concentrations of calcium chloride and how much the strip shortened due to contraction of the tissue was measured.

We can view the data as follows:

```{r}
library(MASS)
head(muscle)
```

The `Length` column gives a measurement of how much the muscle contracted (larger is more contraction).

We can see how many observations we have for each concentration as follows:

```{r}
table(muscle$Conc)
```

And we can plot the data as follows:

```{r, echo=-1, fig.width=6, fig.height=4}
par(mar=c(4,4,1,1))
plot(muscle$Conc, muscle$Length, xlab="concentration", ylab="contraction amount")
```

This question has multiple parts:

1. **Fit two models to these data. In the first model fit the data to a line (first order polynomial).**
2. **In the second model fit a quadratic function (a second order polynomial).**
3. **Use the `summary()` function on models and discuss which one has a better adjusted $R^2$. What happens to the significance of the intercept coefficient as you go from a first order polynomial to a second order polynomial?**
4. **Plot the fit line and quadratic curve with the data on the same plot.**
5. **Do you think this quadratic curve would extrapolate well?**

Hints:

1. The easiest way to fit a second order polynomial is to create a second predictor that is the square of the first predictor and then perform multiple linear regression with these two predictors. You can store this predictor as a new column in the `muscle` data frame, or as a separate variable.
2. One way to plot the fit to the first order polynomial is to use the `abline()` function.
3. To plot the quadratic fit, you can create a sequence using `seq()` that goes from 0 up to the maximum salt concentration and then use the coefficients from your polynomial fit to compute the vector of $y$ values that corresponds to the sequence of $x$ values. Then you can use the `lines()` function to plot the curve. If you use a small enough increment for `seq()` it will look like a smooth curve.

**A7:**
As you go from first order polynomial to a second order polynomial the intercept coefficent  increases from -0.34427 to 0.066538. I do think the quadratic curve could extrapole well. 
```{r}


fit1 <- lm(muscle$Length ~ muscle$Conc)
fit2 <- lm(muscle$Length ~ muscle$Conc + I(muscle$Conc^2))
summary(fit1)
summary(fit2)


seq1 <- seq(1,max(muscle$Conc))

co <- fit2$coefficient[seq1]

plot(muscle$Conc, muscle$Length, main="Scatterplot",  xlab="concentration", ylab="contraction amount") + abline(fit1, col= "blue" ) + lines(seq1, co, type="l", col="red", lwd=2)
  
                                                                                                               
```


---