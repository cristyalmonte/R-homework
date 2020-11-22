---
title: "Problem Set 2"
author: "(Cristy Almonte)"
date: 'due: 5pm, Friday, 20 September 2019'
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
```

#### Q1

An important part of learning R and programming more generally is to be able to figure things out for yourself from the documentation. The goal of this problem is to use the built-in R help to figure out how to accomplish something with a function we haven't covered in class.

You've seen that you can split a DNA string into individual bases as follows:

```{r}
tataBox <- strsplit("TATATAA", split="")[[1]]
print(tataBox)
```

**Use the built-in help to read about the function `paste()` and then use it to put the string back together into the original sequence, `TATATAA`.**

This sequence is an instance of what's called a TATA box motif, which is commonly found at eukaryotic promoters.

**A1:**

```{r}
paste(tataBox, collapse = "")

```
  
#### Q2

In order to write your own functions in R, it's helpful to have a solid understanding of what the syntax of a function is.

Consider the following function:

```{r}
hypotenuse <- function(a,b) {
  sqrt(a^2 + b^2)
}

```

**Answer the series of questions in the A2 answer section about this function.**

**A2:**

*(Some answers may warrant code, in which case you should insert code blocks into the markdown. Other times a simple explanation in words is fine.)*

1. What is the name of the function? hypotenuse 
2. How many arguments does this function take? 2
3. How would you call this function? Pick any numbers you want to use and demonstrate.
```{r}
hypotenuse(1,2)
```

4. Does this function work with vectors longer than length 1? Demonstrate.
Yes
```{r}
hypotenuse(c(1, 2), c(7, 9))
```

5. The variables used for the arguments of the function can be named anything you like. Rewrite the function so that the arguments have names of your choice. For example, you could call them thelma and louise if you want. Confirm that you can call the function and it gives you back the correct result.
```{r}
hypotenuse <- function(bob,bill) {
  sqrt(bob^2 + bill^2)
}
```


#### Q3

An important part of understanding R and using R is having an accurate mental model for how the code is parsed and how the computations are performed. In this question I hope to help you build your own capabilities for understanding how complicated expressions can be broken down and sequenced like the R interpreter itself operates.

Code naturally has a nesting structure. For example, if the argument for a function is itself an expression that requires computation, the expression comprising the argument must be computed before the result can be passed to the function.

Consider this example:

```{r}
sum(c("A", "T", "C", "C") == "C")
```

Before the function, `sum()`, can be called, it's argument must be evaluated. That argument is this:

```{r}
c("A", "T", "C", "C") == "C"
```

But this itself is an expression that includes the operator `==` which has two operands, `c("A", "T", "C", "C")` and `"C"`.

So before R can evaluate the logical operator, `==`, it must call the combine function, `c()`:

```{r}
c("A", "T", "C", "C")
```

Once R has created the vector of letters, it can evaluate the logical expression. And once it has evaluated the logical expression, it can compute the sum.

The expression, `c("A", "T", "C", "C") == "C"` is said to be **nested** inside the expression, `sum()`.

We could have computed this step-wise and assigned the values to variables:

```{r}
s1 <- c("A", "T", "C", "C")
s2 <- s1 == "C"
sum(s2)
```

This returns the same result as above.

I'll give another example of how to break down a complicated expression into line-by-line expressions assigned to variables.

Here's the full expression:

```{r}
sort(table(strsplit("ACAGACATGAAAGCTAACCGAGGACTTGAGAGACTCAA", "")[[1]]))[1]
```

And here I've broken it down into expressions according to the nesting structure:

```{r}
a1 <- strsplit("ACAGACATGAAAGCTAACCGAGGACTTGAGAGACTCAA", "")
a2 <- a1[[1]]
a3 <- table(a2)
a4 <- sort(a3)
a4[1]
```

**Your goal is to split the following expression into individual lines, assigning each intermediate computation to a variable, like above. Try and explain what it is doing.**

Here's the expression:

```{r}
mean(diff(c(1,3,4,7,10,12,13,19,20,31,34,39,51)) >= 10)
```

I have purposely included a function you haven't seen before because this strategy of breaking down code into pieces is a good way to figure out what's going on. You can look up `diff()` in the help as well as inspect the output of `diff()` in this particular case, compared with it's input.

**A3:**

```{r}
x1 <- c(1,3,4,7,10,12,13,19,20,31,34,39,51)
cat("My vector: ", x1,  "\n") 
x2 <- diff(x1)
cat("Differentiated vector: ", x2, "\n")
x3 <- x2 >= 10
cat("Boolean vector: ", x3, "\n")
x4 <- mean(x3)
# remember the mean function converts TRUE -> 1 and FALSE -> 0
cat("Mean of boolean vector", x4, "\n")



```

First, the program  combines  the elements into a vector. Then, using the diff() it computes the finite difference approximantion($\frac{F2-F1}{\Delta x}$).  Then, the program uses a boolean vector for the values  >= 10, which gave us two TRUE VALUES.  Then the program computes the mean, while implicityl converting  all  TRUE and FALSE values to 1 and 0 respectively. 

#### Q4:

The reverse complement of a DNA sequence is the sequence found on the opposite strand of double-stranded DNA and in the reverse order.

So if we have this sequence, `AGGTCCC`, then the reverse complement is `GGGACCT`. Adenine pairs with thymine and cytosine pairs with guanine.

**Your goal is to write a function that reverse-complements a DNA sequence.**

It should behave in such a way that this statement evaluates to `TRUE`:

```{r, eval=FALSE}
revComp("AGGTCCC") == "GGGACCT" 
```

What follows is some information that should be of use.

The function, `rev()`, will reverse the order of the components in a vector:

```{r}
rev(c(1,2,3))
```

A good way to look up the corresponding nucleotides is to create a `list` that serves as a lookup table:

```{r}
complements <- list(A="T", T="A", G="C", C="G")
```

We can then look up individual DNA base complements using the **double bracket index operator**, `[[]]`:

```{r}
complements[["A"]]
```

Or we can look up multiple ones by using the **single bracket index operator**, `[]`:

```{r}
v <- c("A", "G")
complements[v]
```

Hint: You will likely need to use your insight from Q1 to answer this question.

**A4:**

Fill in your solution in this code block:

```{r}
revComp<-function(x){

  bases <- c("A","C", "G","T") 

  rev_vector <-rev(unlist(strsplit(toupper(x),NULL)))

  #cat("xx is: ", rev_vector, "\n")
  
  paste(unlist(lapply(rev_vector,function(element){
  
    if(element=="A") 
      elementComplement<-"T"
  
    if(element=="C") 
      elementComplement<-"G"
  
    if(element=="G") 
      elementComplement<-"C"
  
    if(element=="T") 
      elementComplement<-"A"
  
    if(!element %in% bases) 
      elementComplement<-"N"
  
    return(elementComplement)

}
)),collapse="")
}
 
  

```

Once you've implemented the above `revComp()` function, these test cases below should all evaluate to `TRUE` instead of `logical(0)`, which is the empty logical vector. That is to say, you should have `TRUE` repeated 6 times below.

```{r, results="hold"}
revComp("AGGTCCC") == "GGGACCT"
revComp("A") == "T"
revComp("") == ""
revComp("GGGGGGGGGT") == "ACCCCCCCCC"
revComp("CAA") != "AAC"
revComp("TAG") != "GAT"
```

#### Q5:

One of the organisms with the lowest GC content is *Plasmodium falciparum*, the protozoan parasite responsible for malaria.

The sequence for one of *P. falciparum's* many stevor genes is included with the problem set in a fasta file, `stevor.fasta`. You can read it into a character vector as follows:

```{r}
stevor <- scan(file="stevor.fasta", what="", skip=1)
class(stevor)
```

In the lecture notes we saw one possible implementation of a function, `gcContent()`, for computing the GC content for a DNA sequence. Another way to implement this function would have been to use the `table()` function.

If we naively try and tabulate DNA bases using `table()` using a character vector, we will miss the zero counts if not all bases are present in the sequence:

```{r, results="hold"}
table(c("A", "T", "G", "C", "C", "T"))
table(c("T", "T", "G", "C"))

```

In **Applied Lesson 3** we saw that factors can be used in conjunction with `table()` for counting instances of each category in nominal data.

This works even when not all the categories are present.

For example, if we want to compute how often each base occurs, we can use a factor and explicitly give a vector for the `levels` parameter, listing all the possible DNA bases:

```{r}
bases <- c("A", "T", "G", "C")
table(factor(c("A", "T", "T", "G"), levels=bases))
```

**See if you can use this strategy to write a new version of the `gcContent()` function that doesn't use the `sum()` function or logical vectors. Use this method to then compute the GC content of the *P. falciparum* stevor gene. If the nucleotides were all equally frequent, what would the GC content have been?**

**A5:**

```{r}
gcContent2 <- function(dna) {
  letters <- strsplit(dna, "")[[1]]
  bases <- c("G", "C", "A", "T")
  content = table(factor(letters, levels=bases))
  l = length(letters)
  
  (content[["C"]] / l) + (content[["G"]] / l) 
}


```

```{r}
#gcContent2("GGATCACTTGAACCCAAGAAGTGGAGGTTAAAAAAAAAAAAAAAAAAGGAATGTTTT")
gcContent2(stevor)

```


#### Q6

I have included a file in the problem set named `smoking.tsv`. 

**See if you can read this file in using `read.table()` and organize the data into a contingency table. Once you have your table, use the `chisq.test()` function to evaluate whether smoking status and strokes are independent or assorted in some non-random pattern.**

This data set details the study participants' status as a smoker and whether the individual had a stroke during the observation period. The variable encoding the participants' smoking status has three categories. These should be treated as nominal, unordered data.

Use the `table()` function to create your contingency table.

**A6:**

```{r}
x <- read.table("smoking.tsv",header=T,sep="\t")
#print(x)

smoker = x$smoker
smokerTable = table(smoker)

haStroke = x$haStroke
haStrokeTable = table(haStroke)

tbl = table(x$smoker, x$haStroke)
tbl

chisq.test(tbl)

```
The Chi-squared test indicated a p-value of 0.1171, which indicates a non significant value. 
#### Q7

Imagine a lonely polar bear wandering in the arctic foraging for food. The bear is tracked with a GPS and so we have (x,y) coordinates giving the location at various times over a 1-week period.

Here I simulate such a data set as a random walk:

```{r}
set.seed(3)
tm <- cumsum(abs(rnorm(100, mean=1, sd=0.5)) + 0.5)
x <- cumsum(rnorm(100))
y <- cumsum(rnorm(100))
path <- data.frame(day=ceiling(tm/24), hour=(tm %% 24), x=x, y=y)
head(path, 10)
tail(path)
```

We can visualize it like so:

```{r, fig.height=4, fig.width=4, echo=-1}
par(mar=c(4,4,1,1))
{
  plot(path$x, path$y, type="l", xlab="<= west  (x)  east =>", 
       ylab="<= south  (y)  north =>")
  points(path$x, path$y, pch=20, col="red", cex=0.5)
  points(path$x[1], path$y[1], pch=20, col="blue", cex=1.5)
}
```

I've tried to indicate cardinal directions (north, south, east, west) on the axes of the plot, with directions indicated with arrows.

This question is designed to practice indexing things and using logical vectors. You may also need basic functions like max(), min(), range(), median(), mean(), head(), tail(), etc.

**Your goal is to write 1-line R expressions for each prompt in the answer section for this question.** 

**A7:**

What is the x coordinate of the bear at the beginning of day 2?
```{r}
path[path$day == 2,]$x[1]
```

How many data points were collected on the 7th day?
```{r}
length(path[path$day == 7,]$day)
```

Between the hours of 10pm and midnight, what is the furthest West the bear travels? (notice the clock is a 24 hour clock)
```{r}
min(path[path$hour >= 22,]$x)
```

There is only one day when the bear is West of -8, which day is this?
```{r}
path[path$x < -8,]$day[1]
```

What is the furthest the bear gets South on the second day?
```{r}
min(path[path$day == 2,]$y)
```

**BONUS:** If we assume the bear always walked exactly in a straight line from point to point, how far did the bear travel on average per day?


```{r}
#go through every day
#find min and max
#take difference
#thats how far you traveled that day
#keep running total of distance traveled
# divide by # of days

totalDistance =0

for (i in 1:7)
{
  x_min = min(path[path$day == i,]$x)
  x_max = max(path[path$day == i,]$x)
  
  y_min = min(path[path$day == i,]$y)
  y_max = max(path[path$day == i,]$y)
  
  totalDistance = totalDistance + sqrt((x_max - x_min)^2 + (y_max - y_min)^2)
  
  #total_distance = total_distance + max + abs(min)
  print(totalDistance)

}
cat("Average distance bear friend traveled per day:", totalDistance/7)

```

