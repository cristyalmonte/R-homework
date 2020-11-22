---
title: "Problem Set 7"
author: "Cristy Almonte"
date: "due: 5pm, Friday, 22 November 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
options(width=100)
```


*Keep in mind you're responsible for responding to each part of each question. Often the bolded part of the question has multiple sentences or questions, each asking something of your solution.*

---

#### Q1

The formula for the slope of a simple linear regression is:

$$
\hat{b}_1 = \ \frac{\sum_{i=1}^n (x_i - \bar{x})(y_i - \bar{y})}{\sum_{i=1}^n (x_i - \bar{x})^2} \\
$$

It is possible to compute the slope coefficient another way using the something of the following form:

$$
\hat{b}_1 = \textrm{cor}(x,y) \frac{c}{d}
$$

Suppose we have the data:

```{r}
x <- c(1,5,3,9,5,2)
y <- c(2.5, 6.1, 1.6, 11.1, 3.9, 1.2)
```

**What are $c$ and $d$? Express the numerator, $c$, and the denominator, $d$, in conceptual terms, like we did for the left part, $\textrm{cor}(x,y)$, not as a summation of individual $x_i$ or $y_i$ values. Next, confirm that you can use this method involving $\textrm{cor}()$ to compute the same result for $\hat{b}_1$ as you get when using the `lm()` method in R using the toy data in `x` and `y`. Can you offer an interpretation for why this formula makes sense that relates how linear regression is similar to, yet different from computing a correlation?**

**A1:**

c is the standard deviation of y and d is the standard deviation of x.  Both linear regression and correlation quantify the direction of the relationship between two variables.Correlation quantifies the direction and strength of the relationship between two numeric variables between -1.0 and 1.0. Simple linear regression relates X to Y through an equation Y = a + bX.

```{r}

c <- sd(y)
d <- sd(x)
b1 <- cor(x,y) * c/d
b1

b1lm <- lm(y ~ x)
b1lm


```

---

#### Q2

Suppose we have a vector of numerical data, `v`:

```{r}
v <- c(0.1, 0.2, 0.8, 0.4, 0.41, 0.48, 0.52, 3.9)
```

And we are interested in the following statistic:

```{r}
sd(v)/mean(v)
```

This statistic is called the **coefficient of variation**.

Next suppose we want to use the `boot()` function from the **Applied Lesson 9** notes to perform bootstrapping.

Consider the following method of generating the distribution of the estimator of the coefficient of variation:

```{r, eval=FALSE}
distr <- replicate(1000, sd(boot(v))/mean(boot(v)))
hist(distr,breaks=50,col="gray80",main="Bootstrapdistribution",xlab="sd(z)/mean(z)")

sd(distr)
```

**What is wrong with this method of doing bootstrap? What code would we write to fix it? Try computing the standard deviation of the statistic using both methods on the given test data, `v`. What do you observe?**

**A2:**
This method of using bootstrap does not include an indices parameter that the boot() function can use to select for each replication. 

For the first method we have a higher standard deviation meaning that the values are spread out over a wider range. The second method had a lower standard deviation indicates that the values tend to be closer to the mean. 
```{r}
boot <- function(z) {
  # Sample from the vector something the same length, with replacement
  z[sample(1:length(z), length(z), replace=TRUE)]
}

bootDistr <- replicate (1000 , { 
  z <- boot(v)
  sd(z)/mean(z)
  })
sd(bootDistr)
hist(bootDistr,breaks=50,col="gray80",main="Bootstrapdistribution",xlab="sd(z)/mean(z)")

```


---

#### Q3

We have discussed using `sample()` to simulate DNA sequence under a uniform model where all bases are equally likely. In this question you need to use the `sample()` method to simulate AT-rich or GC-rich sequences. Look at the help page, `?sample`, to see how this can be done using the `sample()` function.

**Simulate 4 sequences, each of length 100, but each with a different level of GC content: 90% GC, 70% GC, 50% GC and 25% GC. Prob(A) == Prob(T) and Prob(G) == Prob(C). The resulting object should be a character vector with 4 strings, each string length 100. Then compute the observed GC content percentage for each of these sequences.**

**A3:**
```{r}
seq1 <- sample(c("A","T","G","C"), 100,replace=TRUE,prob = c(10,10,90,90))
seq2 <- sample(c("A","T","G","C"), 100,replace=TRUE,prob = c(0.30,0.30,0.7,0.7))
seq3 <- sample(c("A","T","G","C"), 100,replace=TRUE,prob =  c(0.5,0.5,0.5,0.5))
seq4 <- sample(c("A","T","G","C"), 100,replace=TRUE,prob = c(0.75,0.75,0.25,0.25))


gcContent <- function(letters) {
  mean(letters == "G") + mean(letters == "C")
}

print(gcContent(seq1))
print(gcContent(seq2))
print(gcContent(seq3))
print(gcContent(seq4))
```

---

#### Q4

There is one more function that comes with every probability distribution in R, the distribution's **quantile function**. For the normal distribution this function is `qnorm()`. And for all the other distributions the corresponding functions are also prefixed with `q` (e.g., `qpois()`, `qexp()`, etc.). This function serves a similar purpose to the `quantiles()` function, but instead of giving you the quantiles for an observed data set it gives you the theoretical quantile for a member of the distribution family. This allows you to translate a tail probability into a value from the random variable, which is particularly helpful for thinking about hypothesis testing, computing power, type I and type II error, etc.

For example, if we want to know what threshold corresponds to a (one-sided) type I error rate of 5%, for a particular distribution we simply write:

```{r}
qnorm(0.05, mean=0, sd=1)
```

This means that a there is a 5% change of observing a value less than `r round(qnorm(0.05, mean=0, sd=1), 6)`.

If we are interested in the high tail, we can give `qnorm()` the `lower.tail=FALSE` argument:

```{r}
qnorm(0.05, mean=0, sd=1, lower.tail=FALSE)
```

This time it means there is a 5% chance of observing a value above `r round(qnorm(0.05, mean=0, sd=1, lower.tail=FALSE), 6)`.

Notice that `qnorm()` is the inverse of `pnorm()`.

We can plug the quantile we got from `qnorm()` back into `pnorm()` and get back 5%:

```{r}
pnorm(-1.644854, mean=0, sd=1)
```

Suppose we have a null distribution parameterized as follows:

```{r}
muH0 <- 100
sdH0 <- 10
```

If we plot it we get:

```{r, echo=-1, fig.width=7, fig.height=4}
par(mar=c(4,4,1,1))
# Plot null distributions, H0
x <- seq(20, 150, length.out=200)
plot(x, dnorm(x, mean=muH0, sd=sdH0), type="l", xlab="x", ylab="density")
text(110, 0.035, "H0")

```

Now let's consider testing an alternative hypothesis for a process that has a lower mean and an effect size of 28 units (with the same variance).

That is to say, we parameterize our alternative hypothesis like this:

```{r}
muH1 <- 72
sdH1 <- sdH0
```

So, assuming the effect size is -28 units, we can plot both distributions together:

```{r, echo=-1, fig.width=7, fig.height=4}
par(mar=c(4,4,1,1))
# Plot the two hypothesis distributions
x <- seq(20, 150, length.out=200)
plot(x, dnorm(x, mean=muH0, sd=sdH0), type="l", xlab="x", ylab="density")
lines(x, dnorm(x, mean=muH1, sd=sdH1), col="red")
text(62, 0.035, "H1", col="red")
text(110, 0.035, "H0")
```

In the context of testing the hypothesis that $\mu_1 < \mu_0$, we wish to calculate our power. Recall that power is the area under the H1 distribution that is more extreme than the threshold set for the allowable type I error.

Suppose we want to only allow a type I error rate of 1%.

**Use `qnorm()` and `pnorm()` to compute our power for the test that $\mu_1 < \mu_0$ when the effect size is -28 units. If we grow the effect size by 50%, making it more negative, what power do we get?**

**A4:**
The power is the probability that we would reject under the alternative hypothesis. 
```{r}

muH0 <- 100
sdH0 <- 10
muH1 <- 72
sdH1 <- sdH0

L=qnorm(0.01/2,muH0 ,sdH0) 
U=qnorm(1-(0.01/2), muH0,sdH0 )
muH1 <- 72
sdH1 <- sdH0
pL=pnorm(L,muH1 , sdH1) 
pU=1-pnorm(U,muH1 , sdH1) 
pL + pU 

# 50% more negative 

muH0 <- 100
sdH0 <- 10
L=qnorm(0.01/2,muH0 ,sdH0) 
U=qnorm(1-(0.01/2), muH0,sdH0 )
muH1 <- 36
sdH1 <- sdH0

pL=pnorm(L,muH1 , sdH1) 
pU=1-pnorm(U,muH1 , sdH1) 
pL + pU 
```
---

#### Q5

In the context of hypothesis testing we have four possible outcomes:

1. **True positive (TP)** - the alternative hypothesis is true and we detected (classified/concluded) it as such.
2. **True negative (TN)** - the null hypothesis is true and we detected it as such.
3. **False positive (FP)** - the null hypothesis is true, but we mistakenly concluded the alternative was true. Type 1 error
4. **False negative (FN)** - the alternative hypothesis is actually true, but we mistakenly concluded the null hypothesis was true.Type 2 error

Suppose you are screening small molecules for ones that can inhibit a particular transmembrane sodium channel protein. For simplicity, assume that candidate molecules will either block or will not block the channel functionality.

Further suppose that of the 5000 molecules you plan to screen, 50 of them are actually capable of blocking the sodium channel.

**If you have 80% power to detect molecules that inhibit the protein and you are conducting your test at the 1% type I error rate, and you screen all 5000 small molecules, what is the expected number of molecules that will fall into each class: TP, TN, FP, FN? What fraction of the molecules that you declared inhibitors are actually inhibitors? Of what you declare to be inhibitors, were most of them actually inhibitors?**

*Hint: You may want to start by lumping pairs of classes together and thinking about how many are in each pair, and then thinking of how these split into the individual classes.*

**A5:**
```{r}
n= 5000
postive = 50
negative = 4500
alpha =0.01
p = 0.8
beta = 1- p
TP1 = 50
FN1 = 4950

TP <- TP1/(TP1+FN1)*5000
TP
TN <-  TP1/(TP1+FN1)*5000
TN
FP <-FN1/(TP1+FN1)*5000
FP
FN <- TN-1
FN



```


---

#### Q6

Recall that the t statistic for a one-sided t test is defined as follows:

$$
t = \frac{\bar{x} - \mu_0}{s/\sqrt{n}}
$$

Where $\bar{x}$ is the observed sample mean, $\mu_0$ is the null hypothesis, $s$ is the sample standard deviation and $n$ is the sample size.

Suppose you are testing the ability of a new drug to remove amyloid beta plaque in a mouse model of Alzheimer's. You want to test the hypothesis that your drug removes plaque.

**If you observe reduction in plaque of 1.7 units with sample standard deviation of 3.9 units, what is the smallest sample size that would have been big enough to have rejected the null hypothesis at a type I error rate of 5%?**

*Hint: You may find the function `qt()` useful for this question.*

**A6:**
The smallest sample size would be n=14
```{r}
mu0 = 1.7          
s = 3.9             
alpha = 0.05

L=qnorm(alpha/2,mu0 ,s) 
U=qnorm(1-(alpha/2), mu0,s)

z=  1.645
q = 3.9
E = 1.7

n = (z*q/E)^2
n


```


---
