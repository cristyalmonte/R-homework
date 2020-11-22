---
title: "Problem Set 3"
author: "(Cristy Almonte)"
date: "due: 11am, Friday, 27 September 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
```

#### Q1:

In **Theory Lecture 3** we learned about a number of transforms, including these widely used ones:

1. $X^2$
0. $\sqrt{X}$
0. $log(X)$
0. **shift-log**: $log(X + c)$ 
0. **z-transform**: $(X-\mu)/\sigma$ 
0. **logit transform**: $log(X/(1-X))$ 

There is a file, `xyz.tsv` in the problem set folder. You can read it in as follows:

```{r}
puzzle <- read.table(file="xyz.tsv", header=TRUE, sep="\t")
head(puzzle)
```

Two of these three samples were simulated in such a way that if you transform them, you will get back a uniform distribution. The third has not been transformed. But I'm not telling you which is which.

The transforms you'll need to use are $log(a)$ and $a^2$.

**Experiment with transforming the data in these three columns and visualizing the resulting distribution with the `hist()` function until you are able to figure out which column needs to be transformed with $log(a)$ and which needs to be transformed with $a^2$ and which doesn't need to be transformed, such that the transformed distributions all look roughly uniformly distributed.**

Hint: When you are plotting histograms, I recommend you set the `breaks` argument so that you have sufficient resolution to see the distribution, but not too many breaks such that you don't get accurate estimates of the relative bar heights.

**A1:**

```{r}
hist(x = puzzle$x, breaks = 1000)

hist((puzzle$x)^2)
hist(puzzle$y)

#hist(x= puzzle$y,breaks=1000) 

#hist(puzzle$z)

hist(log(puzzle$z))

#hist(x= puzzle$z,breaks=1000) 
```

*(Show your code that you use to visualize the histograms and discuss what you observe. Insert R Markdown code blocks where necessary to show your code. Be sure to clearly state which two columns need to be transformed and which transform should be used in each case.)*

The two plots that need to be transformed are x and z, with square function and the log function respectively. 
----

#### Q2:

Consider the **z-transform** that you learned about in the theory lecture:

$$
z = \frac{x - \mu}{\sigma}
$$

The z-transform is generally used to convert a normally distributed random variable to a standard normal random variable. The original random variable could have any mean and/or variance, and by using this transform the transformed random variable would have mean 0 and standard deviation 1.


But suppose you wanted to go the other direction. Suppose you can only simulate samples from a standard normal distribution, but you really want a sample from a distribution with mean 10 and standard deviation 4.

**How would you go about doing this? What transformation should you use? See if you can transform a sample from a standard normal to a normal distribution with mean 10 and standard deviation 4. Plot a histogram of your resulting sample and compute both the mean and standard deviation of the sample.**

**A2:**

*(Explain your transformation and fill in the following code block)*

```{r}

# Transform s so that it has mean 10 and standard devaition 4
set.seed(1)

s <- rnorm(sd = 1,n= 300)
sd(s)
mean(s)
x = s * 4 + 10
sd(x)
mean(x)
hist(x)


```
For this transformation I used the z score equation to add the mean and multiply the standard deviation. 


----

#### Q3

As covered in the theory notes, the Normal distribution, with mean $\mu$ and variance $\sigma^2$, has the following distribution:

$$
f(x) = \frac{1}{\sqrt{2 \pi \sigma^2}} \exp(-\frac{(x-\mu)^2}{2\sigma^2})
$$

R has a function, `dnorm()` that computes the normal density:

```{r}
set.seed(1)
x <- (-4):4
dnorm(x, mean=1, sd=2)
```

**Write a function in R to compute the normal density for a vector of numbers for a given mean and standard deviation. Use built-in math operations and functions to perform the computation, not the dnorm() function itself, obviously. Confirm that your function evaluates to the same values as R's dnorm() function.**

**A3:**

```{r}

normalDensity <- function(x, u, std) {
  
return (1/(sqrt(2*pi*std^2))*exp((-(x-u)^2)/(2*std^2)))
}

```

```{r}
normalDensity(x, 1, 2)
```




----

#### Q4

In this question we'll work with the cumulative distribution function, or CDF, of the normal distribution.

R also has a function, `pnorm()`, for computing the cumulative distribution function of the normal distribution:

```{r}
pnorm(-0.3, mean=1, sd=1.5)
```

As discussed in the theory lectures, this function computes the following integral:

$$
P\{ X \le a \} = F(a) = \int_{-\infty}^{a}\frac{1}{\sqrt{2 \pi \sigma^2}} \exp(-\frac{(x-\mu)^2}{2\sigma^2})dx
$$

Consider the following three visualizations of the standard normal distribution. Each plot has a region shaded, corresponding to part of the density.

```{r, echo=FALSE, fig.width=8, fig.height=2.2}
par(mar=c(4,4,3,1), mfcol=c(1,3))
normalAreaUnderCurve <- function(low, high, title) {
  xr <- seq(low, high, length.out=100)
  yr <- dnorm(xr)
  x <- seq(-4, 4, length.out=100)
  col <- rgb(1, 0, 0, 0.2)
  plot(x, dnorm(x), type="n", bty="n", ylab="density", xlab="x")
  polygon(c(xr[1], xr, xr[100]), c(0, yr, 0), 
          col=col, border=col)
  lines(x, dnorm(x))
  mtext(title, side=3, adj=0, line=1)
}
normalAreaUnderCurve(-4, -2, "A")
normalAreaUnderCurve(-1, 1, "B")
normalAreaUnderCurve(0, 4, "C")
```

**Compute the probabilities corresponding to the shaded regions in the three plots above, for the standard normal distribution, $X \sim N(0,1)$. Use the pnorm() R function to perform your calculations.**
w
You may want to read the help page for the normal distribution, which can be reached in the console with `?Normal` or by looking up the help for any of these functions, `?pnorm` or `?dnorm`.

In the above plot, only the portion of the distribution in the domain $x \in (-4,4)$ is shown. For plots A and C, assume that the shaded region extends to negative infinity and positive infinity respectively.

*Hint:* In performing these computations, you many need to call `pnorm()` more than once.

**A4:**

```{r}
# Use pnorm() to compute the areas shaded above for plots A, B, and C, separately.


#The area under the curve for A
areaA <- pnorm(-2,mean=0,sd=1)
print(areaA)
#The area under the curve for B
areaB<- pnorm(1,mean=0,sd=1)-pnorm(-1,mean=0,sd=1)
print(areaB)
#The area under the curve for C
areaC <- pnorm(0,mean=0,sd=1)
print(areaC)

```

----

#### Q5

Recall from the theory lectures that the expectations and variances of two independent, random variables, $X_1$ and $X_2$, obey the following algebraic properties (where $a$ is a constant):

$$
\begin{align}
E[X_1 + X_2] &=  E[X_1] + E[X_2] \\
E[a X_1] &=a E[X_1] \\
Var[X_1 + X_2] &= Var[X_1] + Var[X_2] \\
Var[X_1 - X_2] &= Var[X_1] + Var[X_2] \\
Var[a X_1] &= a^2 Var[X_1]
\end{align}
$$

Assume that:

$$
\begin{align}
E[X_1] &= m_1 \\
E[X_2] &= m_2 \\
Var[X_1] &= v_1 \\
Var[X_2] &= v_2
\end{align}
$$

A constant, $b$, has an expectation that is equal to itself, and no variance (because it's constant):

$$
\begin{align}
E[b] &= b \\
Var[b] &= 0
\end{align}
$$

**See if you can use the above properties of expectations and variances to write formula for the following expressions:**

1. $E[2X_1 + 3X_2]$
2. $E[X_1 - X_2]$
3. $E[X_1 + 1]$
4. $Var[\frac{X_1}{3}]$
5. $Var[X_1 + X_1 + X_1]$
6. $Var[X_1 - X_1]$
7. $Var[2X_1 - 3X_2]$
8. $Var[X_1 + 17]$

Note: The random variable $X_1 + X_1$ is a random variable formed by sampling once from the $X_1$ process and then sampling a new sample from the $X_1$ process and adding the two outcomes together. This is not equivalent to sampling once and multiplying the result by 2.

**Bonus:**

Assume we have an infinite number of independent and identically distributed random variables, $Y_i$ such that $i \in {1, 2, ...}$ for all positive integers. This implies that all the variables, $Y_i$ have the same mean and variance, call them $m = E[Y_i]$ and $v = Var[Y_i]$.

What does this equal?

$$
\sum_{i=1}^{\infty}E\bigg[\frac{Y_i}{2^i}\bigg]
$$

How about this?

$$
\sum_{i=1}^{\infty}Var\bigg[\frac{Y_i}{\sqrt{2}^i}\bigg]
$$

*Hint: A problem closely related to these two summations was discussed in class.*

**A5:**

*Fill in expressions for the numbered items above. You don't need to use any R for this question. For the purposes of writing your solutions, you don't need to worry about formatting it pretty in R Markdown. So for $m_1$ you can just write m1, etc.*




1. $E[2X1+3X2] = 2m1+3m2$ 
2. $E[X1−X2]= m1-m2$
3. $E[X1+1] = m1 + 1$
4. $Var[1/3*X] = 1/9v1$
5. $Var[X1+X1+X1] = 3v1$
6. $Var[X1−X1]= 2v1$
7. $Var[2X1−3X2]= 2v1+3v2$
8. $Var[X1+17]=v1$


----

#### Q6a

As mentioned in the theory lectures, a major driver of many advances in biology in the past two decades has been the cost-effectiveness of obtaining sequencing data. 

The first effort to sequence a whole human genome, The Human Genome Project cost \$2.7 billion. By the time the first genome was nearing completion the cost per Mb would put the cost at sequencing another whole genome at about $100 million (6x coverage of Sanger sequencing). Today, sequencing a genome costs about \$1000 (30x coverage of Illumina sequencing).

This enormous reduction in cost has enabled:

1. The sequencing and annotation of whole genomes
0. Profiling the expression of whole transcriptomes
0. Searching large populations for genome-wide associations of disease phenotypes
0. Identification of rare mutations, structural variations, etc.
0. Metagenomic studies of complex microbial ecosystems
0. Epigenetic and other functional profiling of DNA state
0. and much more

[NHGRI](https://www.genome.gov/about-genomics/fact-sheets/DNA-Sequencing-Costs-Data) has kept tabs on how much sequencing has cost over time.

We can read in this data as follows:

```{r}
d <- read.table(file="sequencing-cost.tsv", stringsAsFactors=FALSE, 
                header=TRUE)
head(d)
```

If we plot the raw data we get this:

```{r}
plot(d$year, d$costPerMb, type="l", xlim=c(2000, 2020),
     ylab="cost per Mb ($)", xlab="year",
     main="The cost of sequencing over time")
```

This plot makes it look like the cost has crashed to zero. But is that what's really going on?

If we zoom in by including just data from the flat portion of the above plot, say post 2010, we get a very different picture:

```{r}
sel <- d$year > 2010
plot(d$year[sel], d$costPerMb[sel], type="l", xlim=c(2000, 2020),
     ylab="cost per Mb ($)", xlab="year",
     main="Cost of sequencing since 2010")
```

This suggests that the first visualization is not capturing well what is happening in the data. In the first plot the cost looks flat from 2010 onward. But from the second plot we can see the price is still falling during this time period. It may be we still care greatly about the price. Instead of sequencing megabases per week, we're now sequencing gigabases every day, and so we very much still care about additional price decreases.

As we learned in the theory lectures, it may be transforming the data is warranted.

In R, there are two ways to plot your data on a log scale:

1. Transform your data first, and plot it like usual.
2. Leave your data untransformed and tell the `plot()` function to transform it for you.

We can accomplish the second of these options by adding the option `log="y"` to the `plot()` function to indicate that the y-axis of the plot should be drawn on a log scale. (Note: if we had wanted to transform both axes, we would have said `log="xy"`). An advantage of letting R do the transform is that the original units used on the axis are preserved, only the scaling is changed.

**Plot the cost of sequencing data both ways.**

**A6a:**

Here, transform your data first, and then plot it:

```{r}
x <- d$year
print(x)
plot(d$year, d$costPerMb, type="l", xlim=c(2000, 2020), 
     ylab="cost per Mb ($)", xlab="year",
     main="Cost of sequencing since 2010")
ylab = "cost per Mb (log scale)"

y <- log(d$costPerMb)

plot(x, y,type="l", xlim=c(2000, 2020), ylab=ylab, xlab="year", main="The cost of sequencing over time")


```

Here, leave your data untransformed, and add the `log="y"` option to the plot.

```{r}

y <- d$costPerMb
# Give your y-axis a title that correctly reflects what is shown.
ylab = "cost per Mb (log scale)"

plot(x, y, log="y",type="l", xlim=c(2000, 2020),
     ylab=ylab, xlab="year", main="The cost of sequencing over time")
```

*Hint: Other than the numbers on the y-axis, and the y-axis label, the two plots should look identical.*

----

#### Q6b

Many people have compared the cost of sequencing to Moore's law, which has been variously stated as the number of transistors in integrated circuits doubles every two years (or every year or every 18 months, depending on the formulation).

```{r echo=FALSE, out.width = "75%"}
knitr::include_graphics("moores-law.png")
```

For years, people have been mentioning on Twitter how the cost of sequencing is beating Moore's Law:

```{r echo=FALSE, out.width = "75%"}
knitr::include_graphics("two-tweets.png")
```

But then it was pointed out that in recent years the pace has seemed to stagnate:

```{r echo=FALSE, out.width = "45%"}
knitr::include_graphics("tweet1.png")
```

In this question we are going to consider a transformation of the data that helps us see what's going on.

An exponential process is one where the rate of increase is proportional to the current amount. This can be described by the following differential equation:

$$
\frac{dY}{dt} = cY
$$

Which has the following solution:

$$
Y = y_0 e^{ct}
$$

Here, $t$ is time and $c$ is how much the quantity increases (or decreases if $c < 0$) over time.

Many phenomena in physics, chemistry, biology, economics, and finance follow this sort of process, including population growth, spread of infections, nuclear reactions, compound interest, and many others.

In R we can generate such a dataset as follows:

```{r, echo=-1, fig.height=3, fig.width=5}
par(mar=c(4,4,1,1))
times <- seq(1, 10, length.out=100)
C <- -1
y <- 10000 * exp(C * times)
plot(times, y, type="l", xlab="t", ylab="y", main="")
```

If we transform the data to a log scale we get a linear relationship:

```{r, echo=-1, fig.height=3, fig.width=5}
par(mar=c(4,4,1,1))
plot(times, log(y), type="l", xlab="t", ylab="y", main="")
```

When the time points are evenly spaced, they have the characteristic that the ratio of successive values are all equal to some constant, $k$:

$$
\frac{y_t}{y_{t+1}} = k
$$

For example:

```{r}
c(y[1]/y[2], y[50]/y[51])
```

**Using the cost of sequencing data, `d$costPerMb`, compute the ratio of successive time points. Specifically, for each pair of successive time points, $y_t$ and $y_{t+1}$, compute the value $k_t$. Plot these $k_t$ values as a function of time. Explain what you see in your plot and how this might relate to technical change over time. Do you think Moore's Law has applied consistently over this period of time?**

**A6b:**

```{r}
# Specify what the time axis should be, putting the result in the variable x. 
# There should be one less ratio than we have time points.
x <- d$year[1:length(d$year)-1]
#length(x)
y = d$costPerMb[1:length(d$costPerMb)-1]
yp1 = d$costPerMb[2:length(d$costPerMb)]
#yp1

k = y/yp1

#length(k)
#length(x)
#k

plot(x, k, type="l", xlab="t", ylab="y", main="")

```

----

#### Q7

In this question we will investigate why the normal distribution is so common.

Let's first consider a distribution that is decidedly not normal. Here is roughly the size distribution of trout in Larry's Creek in central Pennsylvania:

```{r, echo=FALSE, fig.height=4, fig.width=5.5}
par(mar=c(4,4,1,1))
s <- seq(0, 0.9, by=0.01)
plot(s, dexp(s, rate=8)*0.032, type="l", xlab="weight (kg)", ylab="fish per square meter")
```

Most of the fish are very small, and there are a few large ones.

This is an example of an exponential distribution.

The exponential distribution is parameterized by a single parameter, $\lambda$, called the rate. It has probability density:

$$
f(x) = \lambda e^{-\lambda x}
$$

Although not as common as the normal distribution, the exponential distribution crops up all sorts of places. For example, the waiting time until a mutation occurs at a particular nucleotide, follows an exponential distribution.

Back to fishing. If we were to come back every day and catch one fish for a whole year, we would get a sample like this:

```{r, echo=-1, fig.height=3.8}
par(mar=c(4,4,1,1))
w <- rexp(365, rate=8)
hist(w, breaks=50, xlab="weight", main="")
```

But what if we caught more than one fish each day? Let the weight of the first fish we catch be $w_1$, the second fish, $w_2$, and so on.

How much does our full catch weigh each day? For example, what is the distribution of $w_1 + w_2$?

We can simulate this process and visualize it as follows:

```{r, echo=-1, fig.height=3.8}
par(mar=c(4,4,1,1))
w1 <- rexp(365, rate=8)
w2 <- rexp(365, rate=8)
hist(w1+w2, breaks=50, xlab="weight", main="")
```

This distribution is starting to look less exponential and have a bit of a hump.

As you will see in answering this question, as you add the weights of more and more fish, the distribution becomes more like a normal distribution.

> As more and more additive effects are combined, a distribution becomes more normally distributed. This is called the **central limit theorem**.

This is a major reason why the normal distribution is so common. Many variables of interest are the result of many separate effects.

For example, quantitative traits like height, weight, HDL, etc. are under the control of many genetic loci and many environmental factors. All of these things contribute largely additively, resulting in a normal distribution observed in the population.

**Try repeating this process, but consider the distribution of the sum of the weights of 10 fish. Plot the distribution for w1+w2+...+w10.**

**A7:**

*Note: You don't need to use a loop for this part. Feel free to just copy and paste to achieve your solution. Though if you want to try using a for loop or sapply(), you are welcome to. We are covering these topics this week.*

```{r}
# Put your solution here, including plotting a histogram of the resulting distribution
par(mar=c(4,4,1,1))
w1 <- rexp(365, rate=8)
w2 <- rexp(365, rate=8)
w3 <- rexp(365, rate=8)
w4 <- rexp(365, rate=8)
w5 <- rexp(365, rate=8)
w6 <- rexp(365, rate=8)
w7 <- rexp(365, rate=8)
w8 <- rexp(365, rate=8)
w9 <- rexp(365, rate=8)
w10 <- rexp(365, rate=8)
hist(w1+w2+w3+w4+w5+w6+w7+w8+w9+w10, breaks=50, xlab="weight", main="")
```

**Bonus: Can you do the same for the weights of 100 fish? For this question I highly recommend a loop or the use of replicate() and sapply(). Plot a histogram of the distribution. Can you contrast this with a sample form a normal distribution with the same mean and variance as your 100 fish weight distribution?**

