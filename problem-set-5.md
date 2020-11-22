---
title: "Problem Set 5"
author: "Cristy Almonte "
date: "due: 5pm, Friday, 8 November 2019"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
```

---

#### Q1

The particular sequence of C followed by a G in a genome is denoted CpG to distinguish this sequence of adjacent nucleotides from the complementary pairing of C with G on the opposite strand. The cytosine nucleotide in a CpG dinucleotide may be methylated, though this can cause it to be easily mutated to thymine. As a result, vertebrate genomes often have a much lower incidence of CpG dinucleotides than would be expected by chance based on base frequencies alone. The exception is in the promoter regions of genes, where the CpG dinucleotides are generally not methylated, and thus not susceptible to this biased mutation process.

Consider a genomic region in a gene that is far from the promoter and which happens to have 60 guanine nucleotides. Suppose the probability of a cytosine base just upstream (5 prime) of each guanine under these circumstances is 0.1. Note that this is substantially lower than the overall frequency of cytosine in the genome because of historical mutations that erased many of the methylated cytosine bases that were part of CpG sites.

**Using the binomial distribution as a model for this data, compute probabilities using R for the following circumstances:**

1. The probability of no CpG sites.
2. The probability that all guanines are CpG sites.
3. The probability that at least 20% of the guanines are CpG sites.
4. The probability that at most 10% of the guanines are CpG sites.
5. The probability that the number of CpG sites is a prime number. 

For these calculations you can ignore the complications that arise from overlapping dinucleotide patterns. Just assume a binomial distribution.

**A1:**

```{r}
# 1.
  # k = 0, n = 60, p=0.1
dbinom(x=0, size=60, prob=0.1)

```



```{r}
# 2.
 # k = 60, n = 60, p=0.1
dbinom(x=60, size=60, prob=0.1)
```



```{r}
# 3.
#k= number of sucess at least 20% of 60 trials is at least 12 success. 
# k = 12:60, n = 60, p=0.1
prob1<- dbinom(12:60, size=60, prob=0.1)
sum(prob1)    

```

```{r}
# 4.
#k= number of sucess at most 10%  of 60 trials is at most 6. 
# k = 0:6, n = 60, p=0.1
prob2 <- dbinom(0:6, size=60, prob=0.1)
sum(prob2)

```

```{r}
# 5.
# Here I provide a list of prime numbers that are <= 60. Use this vector to compute
# your answer below.
primesBelow60 <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59)
prob3 <- dbinom(primesBelow60, size=60, prob=0.1)
sum(prob3)

```
---

#### Q2

The Poisson distribution has the interesting property that the mean is the same as it's variance. The goal of this question is to visualize how this holds for parameterizations of lambda that span several orders of magnitude.

Here are a few values of lambda that we are interested in considering:

```{r}
lambdas <- c(0.1, 0.5, 1, 5, 10, 50, 100)

```

**For each value of lambda in the `lambdas` vector, draw a sample of 1000 draws from the Poisson distribution with that lambda parameter. Then plot the means and variances of these 7 samples against the corresponding values of lambda to show that the variances and means are equal.**

**A2:**


```{r, fig.height=6, fig.width=6}
 
 par(mar=c(4,4,1,1))

# Use sapply() to draw a sample for each value of lambda in the lambdas vector.

k = sapply (1:length(lambdas), function(i) rpois(1000,lambdas[i]))
#k


# Compute the means of the samples. Replace this variable with the means of the samples:

x = sapply(1:ncol(k), function(i) mean(k[,i]))

#x

# Compute the variances of the samples. Replace this vector with the variances:

y =  sapply(1:ncol(k), function(i) var(k[,i]))

#y

plot(x, y, log="xy", ylim=c(0.09, 110), xlim=c(0.09, 110), xlab="mean", ylab="var",
     cex.axis=0.85)
abline(a=0, b=1, lty=2, col="gray")
```

---

#### Q3

Consider the letters:

```{r}
sevenLetters <- letters[1:7]
print(sevenLetters)
```

You'll need to keep in mind these basic facts about combinatorics:

- Sequences are ordered. Thus, `abc` is a different sequence from `bca`.
- Sets are unordered. Thus $\{a,b,c\}$ is the same set as $\{b,c,a\}$. And sets cannot contain the same
item multiple times by definition.
- There are $n!$ ways to arrange $n$ distinct objects in a sequence.
- There are ${n \choose k}$ ways to choose $k$ objects out of $n$.

**Answer the questions below (in the A3 answer section) about permutations and combinations.**

Hint: All of these can be done with computations involving the `factorial()` and `choose()` functions and basic math operations (multiplication, division, addition, etc.). You shouldn't need to enumerate anything in R. But in order to learn how to think about combinatorics you may find it helpful to use a pencil and paper to think about sequences and sets.

**A3:**

**How many different sequences of length 7 are there, using each letter exactly once?**

```{r}

factorial(7)
```

**How many different sequences are there of length 1?**

```{r}

factorial(7) / ( factorial(7-1))
```

**How many different sequences are there of length 2, such that no sequence includes the same letter multiple times?**

```{r}

factorial(7) /factorial(7-2)

```

**How many different sequences are there of length 3, such that no sequence includes the same letter multiple times?**

```{r}
factorial(7) /  factorial(7-3)
```

**How many distinct sets of three letters are there?**

```{r}
factorial(7) /  factorial(7-3)

```

**How many sequences are there whereby the sequence is at least length 1 and at most length 7 and which doesn't contain any letters multiple times?**

```{r}

x <-sapply(1:7, function(i) factorial(7) / ( factorial(7-i)))
sum(x)
```

---

#### Q4

One way the Poisson distribution arises is as an approximation to the binomial distribution when $n$ is large and $p$ is small, but the product, $np$ is intermediate sized. For example, perhaps $p = 0.01$ and $n = 800$.

```{r}
n <- 800
p <- 0.01
```

Since the mean of the binomial distribution with is $np$, we get the following mean:

```{r}
n*p
```

The mean of a distribution is its first moment. If we want to approximate the binomial distribution with the Poisson distribution using just the first moment, then we should set the $\lambda$ (lambda) parameter of the Poisson distribution to: `lambda=n*p`.

We know that the variance of a Poisson distribution is equal to it's mean. This is not exactly true of the binomial distribution, but it's not very var off under these circumstances. The binomial variance is $np(1-p)$. But since $p << 1$, this implies that $(1-p)$ is very nearly $1$, which means:

$$
np(1-p) \approx np
$$

In our case we get the following variance for the binomial distribution:

```{r}
n*p*(1-p)
```

Which is very nearly the value of 8 that is the variance of the Poisson that we're using to approximate this binomial distribution.

**Sample 10,000 values from this binomial distribution, and another 10,000 values from the corresponding Poisson distribution. Compute the sample means and variances for each of the two samples. Then use the barplot() method to visualize each of the two distributions.**

When visualizing the distributions, you'll want to somehow count up how many of each outcome you observe for each distribution. You're welcome to use either `table()` or `sapply()` for this purpose.
**A4:**

*(Place your code to sample from the two distributions and visualize the probability mass functions here. Create code blocks as necessary.)
```{r}
## Sample from binomial distribution and compute mean and variance 

x <- rbinom(10000, size=n, prob=p)
table(x)
barplot(table(x), space=0, cex.axis=0.2, col="gray80")
text((0:n)+0.5, -mean(x==n), labels=0:n, cex=0.2)
mean(x)
var(x)

y <- rpois(10000,n*p)
table(y)
barplot(table(y), space=0, cex.axis=0.2, col="gray80")
text((0:n)+0.5, -mean(y==n), labels=0:n, cex=0.2)
mean(y)
var(y)



```




---

#### Q5

The following diagram shows a simulated phylogeny with 5 species, labeled $\{A, B, C, D, E\}$. The branches are labeled $b_1$ through $b_8$. Mutations occurring on each branch are shown as black dots. *(If you are curious how this tree was simulated, you can check the included file, `simulateTree.R`, though it's not necessary for answering this question.)*

```{r echo=FALSE, out.width = "60%"}
knitr::include_graphics("tree.png")
```

We can read in a table giving branch lengths and counts of mutations along the branches:

```{r}
d <- read.table(file="branch-lengths.tsv", sep="\t", stringsAsFactors=FALSE, header=TRUE)
print(d)
```

Assume that the number of mutations on any given branch follows a Poisson distribution with $\lambda$ equal to the branch length. This means that each branch has a different member of the Poisson family representing the distribution of mutations along that length.

**Compute the probability of the observed number of mutations or less on each branch, using the correct $\lambda$ for each branch. Which branch has the lowest probability under this model of mutation?**

**A5:**

The lowest probabiliy is B2 with a probability of 0.05519909. 

```{r}
ppois(2, 2.72, lower.tail = TRUE, log.p = FALSE)
ppois(2, 6.16, lower.tail = TRUE, log.p = FALSE)
ppois(3, 4.30, lower.tail = TRUE, log.p = FALSE)
ppois(6, 6.52, lower.tail = TRUE, log.p = FALSE)
ppois(8, 5.68, lower.tail = TRUE, log.p = FALSE)
ppois(2, 1.14, lower.tail = TRUE, log.p = FALSE)
ppois(7, 5.96, lower.tail = TRUE, log.p = FALSE)
ppois(5, 3.58, lower.tail = TRUE, log.p = FALSE)


```

----

#### Q6

One of the distributions covered in the theory review of univariate random variables is the log-normal distribution. A random variable that when log-transformed is normally distributed is said to be log-normal.

Here are examples of some phenomena that follow a log-normal distribution (culled from Wikipedia and my own Googling):

1. Number of RNA-seq reads in a genomic interval.
0. Relative species abundance.
0. Distribution of minerals in the Earth's crust.
0. Size of fruits.
0. Latent time from investigation to symptoms.
0. Age of onset for Alzheimer's disease.
0. Rainfall from clouds.
0. The distribution of some pollutants.
0. Crystals in ice cream.
0. Abundance of fungal spores in the air.
0. Neuron firing rates.
0. The length of comments posted in Internet discussion forums.
0. How long people linger on a web page.
0. English word lengths.
0. Lengths of chess games.
0. Income distribution (excluding the highest earners, which follows a power law).

One situation where the log-normal distribution can arise is from multiplicative effects. Similar to how the central limit theorem explains how the sum of many additive effects is often normally distributed, the product of many multiplicative effects is often log-normal.

Transcriptional regulation is a good example of this. Let me outline several reasons why multiplicative effects should make sense in this context.

**Patterning is only possible with multiplicative interactions:** Achieving precise expression patterns requires performing the equivalent of logical computations. Suppose you wanted to only express gene A if regulators B and C are both expressed and D is not expressed. 

If one considers expression levels transformed to be between 0 and 1, then this can be represented by the following equation:

$$
A = BC(1-D)
$$

Gene A will only be maximally expressed if both B and C are near 1 and D is near zero. With additive effects it's not possible to get anything similar to a logical AND operator. 

Consider the additive model:

$$
A = B + C - D
$$

In this model, if you express B but not C you'll just get intermediate expression of A rather than getting A only if both B and C are expressed.

**Probabilities of independent events are multiplicative**. Many physical processes are the result of multiple, independent events coinciding.

As we know from basic probability theory, if event $A$ and event $B$ are independent, that is $A \perp B$, then the probability of them jointly occurring is just the product of their individual probabilities:

$$
P(A \cap B) = P(A)P(B)
$$

Suppose we have one factor that controls when transcription happens, say by recruiting the necessary apparatus, and another factor that controls how many molecules of RNA are transcribed given the apparatus has been recruited. We only get transcription when the apparatus is recruited (event A) and when the apparatus produces a molecule (event B). 

This model is inherently multiplicative because the joint probability of independent events is multiplicative. If the apparatus is recruited $F$ fraction of the time, and when it's there, transcription proceeds at rate $R$. Then the number of molecules produced overall is $FR$. Doubling $F$ necessarily doubles overall transcription, regardless of how many molecules are produced when the apparatus is there. Each variable is independent of the scale of the other.

Here's a picture of how complicated transcription can be with lots of different factors binding. Often the kinetics of this sort of situation involve many multiplicative effects.

```{r echo=FALSE, out.width = "50%"}
knitr::include_graphics("regulators.png")
```

Let the combined transcriptional output, $X$, be some product of multiple effects, $X_i$:

$$
X = X_1 X_2 \cdots X_k
$$

Suppose that we can sample transcriptional effects using this function:

```{r}
transcriptionalEffects = function(n) {
  x <- abs(rnorm(n, mean=1.3, sd=0.25))
}
```

The individual effects have a distribution that looks like this:

```{r, echo=-1, fig.height=3, fig.width=5}
par(mar=c(5,4,1,1))
hist(transcriptionalEffects(10000), breaks=50, main="", probability=TRUE, col="gray")
```

It's basically a normal distribution, but I've taken the absolute value so we can't get anything negative because a negative multiplicative effect doesn't make any sense.

**Use the `transcriptionalEffects()` function to simulate the individual effects, $X_i$, for 10,000 genes for each of $k=5$ separate regulators. Compute the combined output, $X$, for each of the 10k samples and plot the distribution using `hist()`, passing it the `probability=TURE` argument. Also try log transforming the distribution and plotting it. Which one corresponds to a log normal distribution?**

**Bonus: See if you can add the theoretical normal and log-normal distributions to the plot using `seq()` to generate a sequence of x-values, `dlnorm()` and `dnorm()` to compute the theretical distributions, and `lines()` to add a curve to each plot indicating the estimated density. Estimate the parameters for the distributions from the data.**

**A6:**
The original data corresponds to the log normal. 

```{r}

set.seed(0)

hist(transcriptionalEffects(10000), breaks=50, main="", probability=TRUE, col="gray")
hist(log(transcriptionalEffects(10000)), breaks=50, main="", probability=TRUE, col="gray")

```


---

#### Q7

A conditional probability expresses the probability that an event occurs given another event has already occurred, or also occurs.

For example, the probability an individual is left-handed is about 0.1. But the probability an individual is left-handed given both parents are left-handed is more like 0.26. Thus, the additional information that both parents are left-handed, changes our assessment of how likely it is the individual is left-handed.

Conditional probability can be expressed mathematically as follows:

$$
\textrm{Prob}\big(A \big| B\big) = \frac{\textrm{Prob}\big(A \cap B\big)}{\textrm{Prob}\big(B\big)}
$$

This is to say, the probability that event $A$ occurs given $B$ occurs is equal to the joint probability of $A$ and $B$, divided by the probability of $B$ alone.

For the handedness example, we're told that "the probability an individual is left-handed given both parents are left-handed is more like 0.26." This corresponds to the $\textrm{Prob}(A | B)$ term. The $A$ event corresponds to the event "an individual is left-handed" and the $B$ event corresponds to "both parents are left-handed."

We can rearrange the above equation to the following:

$$
\textrm{Prob}\big(A \cap B\big) = \textrm{Prob}\big(A \big| B\big)\textrm{Prob}\big(B\big)
$$

If we assume that people marry each other and have kids independent of handedness, then we can assume that $P(B) = 0.1^2 = 0.01$. And therefore $\textrm{Prob}(A \cap B) = 0.26 * 0.1^2 = 0.0026$. Note that this is somewhat greater than the situation where the handedness of all three individuals (both parents and kid) were independent, which would have a probability of $0.1^3 = 0.001$. Although the precise mechanism is not known, handedness is believed to have a substantial genetic component, explaining the non-independence of these events.

The point of explaining conditional probability is to contextualize the following statement about the exponential distribution:

$$
\textrm{Prob}\big(X > S+T\ \big|\ X > S\big) = \textrm{Prob}\big(X > T\big)
$$

This is called the memoryless property of the exponential distribution. This property is called memoryless because the distribution going forward is independent of the past, that is, it has no memory of it. This makes the most sense when using the exponential distribution to model a waiting time. Here, X is how long before some event occurs, such as the failure of a light bulb.

If we model the distribution of the time it takes the light bulb to fail as an exponential distribution, then it has the peculiar property that conditional on the light bulb lasting at least $S$ months, then the probability that it lasts $T$ more months is equal to the light bulb lasting at least $T$ months in the first place.

Let's first have a look at the exponential distribution:

```{r, fig.height=4, fig.width=6, echo=-1}
par(mar=c(4,4,1,1))
s <- seq(0, 100, length.out=5000)
length(s)
plot(s, dexp(s, rate=1/20), type="l", xlab="t", ylab="f(t)")
abline(v=40, col="red", lty=2)
abline(v=80, col="blue", lty=2)
```

As $t$ grows, the density of the exponential distribution becomes vanishingly small. So it seems a bit counter-intuitive to realize that if the light bulb survives 40 months then it's just as likely to survive another 40 months. It seems like under this logic if the light bulb keeps being lucky, and surviving through each 40-month period, then it might just go on surviving forever.

But we can check that the memoryless property does indeed hold.

We can compute the probability that the light bulb survives at least 40 months using the `pexp()` function, which gives the cumulative distribution function for the exponential distribution:

```{r}
pexp(40, rate=1/20, lower.tail=FALSE)
```

And we can compute the probability that the light bulb survives at least 80 months:

```{r}
pexp(80, rate=1/20, lower.tail=FALSE)
```

If we then check the conditional probability that the light bulb survives an additional 40 months, conditional that it survives the first 40 months, this is:

```{r}
pexp(80, rate=1/20, lower.tail=FALSE) / pexp(40, rate=1/20, lower.tail=FALSE)
```

But notice that this is identical to the probability that the light bulb survived the first 40 months. So, regardless of how old the light bulb gets, it's chances of surviving another 40 months never changes.

In this question we examine the memoryless property with a simulated sample of 10,000 light bulbs.

**First sample 100,000 "times until failure" from an exponential distribution with rate parameter, $\lambda$ = 0.05. Then for each whole number $S$ in `1:80`, estimate the conditional probability, $\textrm{Prob}\big(X > S+40\ \big|\ X > S\big)$, from the sample. Plot 1:80 vs your estimated conditional probabilities using `plot()` with `type="l"` set. What do you notice? Show the theoretical value on your plot using the `abline()` method and using the `h=` option.**

Hint: Estimate the conditional probability by estimating $\textrm{Prob}\big(A \cap B\big)$ and $\textrm{Prob}\big(B\big)$ and then taking their ratio. Keep in mind that if $X > S+40$ then it's also necessarily true that $X > S$.

**A7:**

```{r}
s <- seq(0, 100, length.out=5000)
#length(s)

#plot(s, dexp(s, rate=1/20), type="l", xlab="t", ylab="f(t)")

pexp(40, rate=1/20, lower.tail=FALSE)

p = sapply(1:80, function(i) {pexp(40+i, rate=1/20, lower.tail=FALSE) / pexp(i, rate=1/20, lower.tail=FALSE)})
#p

{ plot(s[1:80], p, type="l", xlab="1:80", ylab="estimate") 
  abline(h=0.1353353, col="red", lty=2)
}



```
All the values are found around 0.1353353. 

---

