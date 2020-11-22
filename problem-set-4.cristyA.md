---
title: "Problem Set 4"
author: "Cristy Almonte"
date: 'due: 5pm, Friday, 25 October 2019'
output:
  html_document: default
  pdf_document: default
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, fig.align="center")
set.seed(1)
```

#### Q1

In this question we consider how distances work in a non-Euclidean geometry. And although non-Euclidean geometry can get pretty weird, they still should obey the basic rules we laid out in the theory lectures.

I've included the longitude and latitude for four international airports:

```{r}
airports <- read.table(file="airports.tsv", header=TRUE, sep="\t", stringsAsFactors=FALSE)
```

We can use the `kable()` function, with comes with the `knitr` library to render a pretty table in R Markdown:

```{r}
library(knitr)
kable(airports)

```

Here I've drawn the great circle paths between each pair of airports using Google Earth:

```{r echo=FALSE, out.width = "90%"}
knitr::include_graphics("great-circle-paths.png")
```

Here is a function you can use to compute the great circle distance between two sets of coordinates:

```{r}
# NOTE: You don't need to undertand how this function works.
# Source: https://www.movable-type.co.uk/scripts/gis-faq-5.1.html
greatCircleDistance <- function(coords1, coords2) {
  # Make sure we have matrices
  coords1 <- matrix(as.matrix(coords1), ncol=2)
  coords2 <- matrix(as.matrix(coords2), ncol=2)
  # Convert everything to radians and organize as separate variables
  lat1 <- coords1[,1] * pi/180
  lat2 <- coords2[,1] * pi/180
  lon1 <- coords1[,2] * pi/180
  lon2 <- coords2[,2] * pi/180
  # Compute the central angle between the two radii
  a <- sin((lat2 - lat1)/2)^2 + cos(lat1) * cos(lat2) * sin((lon2 - lon1)/2)^2
  # Compute the portion of the circumference between the two points, The min() is there to
  # increase numerical stability.
  c <- 2 * asin(pmin(1,sqrt(a)))
  # Although the radius depends on location because the Earth is an oblate spheroid, it's not a
  # bad approximation to always use a single value.
  earthRadius <- 6367
  earthRadius * c
}
```

The above function computes the great circle distance between two points anywhere on earth. It takes two matrices, each with two columns giving latitude and longitude expressed in degrees (with fractional parts) for the two locations. It can also take two vectors of length 2 or two data frames with two columns each.

For example, we can see the great circle distance between Narita and LAX is about 8,748 km:

```{r}
greatCircleDistance(airports[1, 5:6], airports[3, 5:6])

```

**Use the function, greatCircleDistance(), and the data frame, airports, to compute a 4x4 distance matrix between all pairs of the 4 airports. Then show that these great circle distances obey the triangle inequality for each path from Philadelphia to Tokyo with at most two segments.**

**A1:**

1. You need every combination 
2.Apply the function to it
3.Create a 4x4 matrix 

```{r}
 vector <- c()

for (i in 1:4) {
  for (j in 1:4){
  
    x = greatCircleDistance(airports[i, 5:6], airports[j, 5:6])
  
    vector <- append(vector, x)
  }
}


data = matrix(vector,nrow = 4,ncol = 4)

data

```


---

#### Q2

As you likely know, the genetic code exhibits degeneracy. That is to say, for many amino acids, there are multiple three-base codons that code for the same amino acid. For example, the amino acid Valine is encoded for by the codons GTT, GTC, GTA, and GTG. Notice that it doesn't matter what the third nucleotide is, as long as the bases at the first two positions are G and T, the third position can be any of the four nucleotides and the protein will still have a V (Valine) amino acid at that position. You can read more about the [degenerecy of the genetic code on Wikipedia](https://en.wikipedia.org/wiki/Genetic_code).

In general the codons tend to be most degenerate at the third position, as in the case Valine. Because often changes at the third positions do not affect the protein that is translated from that gene, evolutionary biologists have posited that these positions in the genome are likely under less selective pressure than the positions that do result in protein coding changes.

Given Halloween is around the corner, in this question we'll work with a gene from *Drosophila* called *disembodied*, or *dib* for short. This gene is one in a family of so called *Halloween genes*, which includes *spook*, *spookier*, *phantom*, *shadow*, and *shade*.

Perhaps you're curious, what is the effect of mutating the *disembodied* gene? If we consider null mutations, which are when the gene function is completely disrupted, the fly doesn't even make it through embryogenesis. Compare the left (wild-type) and right (double knockout) phenotypes below:

```{r echo=FALSE, out.width = "40%"}
knitr::include_graphics("disembodied-phenotypes.png")
```

The rows show different ways of staining *Drosophila* embryos during embryogenesis. The wild-type ones on the left look normal, and the *dib* knockout phenotypes are very severely messed up. Perhaps the very strange looking embryo looks disembodied?

I have included sequences for the protein-coding portion of *disembodied* for both *D. melanogaster* and *D. yakuba*. 

We can read them in as follows:



```{r}
mel <- readLines("mel-dib.seq")
yak <- readLines("yak-dib.seq")
nchar(c(mel, yak))

```

The pairwise sequence alignment for these two protein-coding regions contains no gaps, and these sequences are the same length, extending from start codon to stop codon. Thus we can assume that the i'th position of the *D. melanogaster* sequence corresponds to the i'th position of the *D. yakuba* sequence.

**Compute the Hamming distance between these two sequences, separately for the first, second, and third codon positions. Then make a barplot showing these three distances so we can get a sense of their relative sizes. Be sure to give your barplot labels for the bars and for the axes. What do you observe? How do these distances align with the above discussion of how much selective pressure the third positions are under vs the first and second positions?**

Hints: The sequence has 837 nucleotides. There are 279 first positions, 279 second positions, and 279 third positions. In computing the Hamming distance for first positions, you should use all of the 279 first positions from *melanogaster* and compare these to the 279 first positions of *yakuba*. In Applications Lesson 6 you learned about the modulo operator, `%%`, which could come in handy for this problem. You may want to consult the help of the `barplot()` function to learn how to plot bar plots.

**A2:**


```{r}
hammingDist <- function(seq1, seq2) {
  
  bar_data <- c()
  
  v1 <- unlist(strsplit(seq1, ""))
  v2 <- unlist(strsplit(seq2, ""))
  
  for (i in 1:3) {
    if (i == 1) {
      choose = c(TRUE, FALSE, FALSE)
    }
    if (i == 2) {
      choose = c(FALSE, TRUE, FALSE)
    }
    if (i == 3) {
      choose = c(FALSE, FALSE, TRUE)
    }
    
    seq1_codon = v1[choose]
    seq2_codon = v2[choose]
    
    x = sum(sapply(1:length(seq1_codon), function(i) {sum(seq1_codon[[i]] != seq2_codon[[i]])}))
    
    bar_data <- append(bar_data, x)
    
  }
  
  return (bar_data)
}
print(hammingDist(mel,yak))

mel <- readLines("mel-dib.seq")
yak <- readLines("yak-dib.seq")

barplot(hammingDist(mel, yak), beside=TRUE, ylab="Counts",
 main="Hamming Distances", ylim=c(0,279),
 names.arg=c("Codon One","Codon Two","Codon Three"))

```



---

#### Q3

In this question we examine the iris data again:

```{r}
#head(iris)

#tail(iris)
```
```{r}


```

**Write a function that, for any of the three species, identifies which pair of features is most highly correlated. Obviously any feature is perfectly correlated with itself, so these should be excluded from consideration. Use your function once for each species to find which pair of features, is most highly correlated for that species. Thus your answer should be three pairs of features.**

Your function should be passed a single argument, the species, and it should return a character vector of length 2 giving two feature names.

For example, you would call it like this:

```{r, eval=FALSE}
 mostCorrelatedFeatures("setosa")

```

And it would return something like this (this is only an example, not the features most correlated for setosa:

```
[1] "Sepal.Width"  "Petal.Width"
```

Hints:

1. You might consider using a pair of nested `for` loops, or the function `expand.grid()`.
2. I recommend getting things working for one species before turning your code into a function that takes the species as an argument. 
3. Writing this function will be more involved than functions you've written for previous problem sets.

**A3:**

```{r}

# setosa, virginica, versicolor
mostCorrelatedFeatures <- function(flower) {
  
  features = c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width")
  
  current = 0 
  max_val = 0
  corr1 = ""
  corr2 = ""
  
  for (feature_1 in features) {
    for (feature_2 in features) {
      
    vector_1 = iris[which(iris$Species == flower), ][1:4][feature_1]
    vector_2 = iris[which(iris$Species == flower), ][1:4][feature_2]
    
      current = cor(vector_1, vector_2)
     
      
      if (current > max_val & current < .99) {
        max_val = current
        corr1 = feature_1
        corr2 = feature_2
      }
    }
  }

  return (paste(corr1, " ", corr2))
}
print (mostCorrelatedFeatures("setosa"))
print(mostCorrelatedFeatures("versicolor"))
print(mostCorrelatedFeatures("virginica"))


```

---

#### Q4a

Consider the following Cartesian plane with lettered points:

```{r, fig.width=5, fig.height=5, echo=FALSE}
# Don't worry about undertanding this code chunk.
par(mar=c(4,4,1,1))
xPoints <- c(0, -2, 0, 2, -1, -2, 3, 0, 3, 0, 0, -1, -4, 3, -4, 2)
yPoints <- c(-1, 0, -2, 0, 4, 5, 0, 4, 3, 5, 3, 0, -1, -1, 0, -2)
h <- max(abs(c(xPoints, yPoints)))
plot(c(-h, h), c(-h, h), xlab="", ylab="", bty="n", axes=FALSE, type="n")
abline(h=0)
abline(v=0)
points(xPoints, yPoints, cex=4.0, pch=20, bg="white")
text(xPoints, yPoints, labels=LETTERS[1:length(xPoints)], col="white")
```

There are 6 points that are not on either the vertical or horizontal axis.

**For each of the 6 points make a list giving the letters of the points projected onto the x-axis and make another list of the 6 points giving the letters of their projections onto the y-axis.**

Hint: This problem doesn't require any coding. Just visually inspect the plots and write down the answers.

**A4a:**

Use the lists below and fill them by placing a letter after each arrow.

Projections onto the x-axis:

```
E -> L
F -> B
I -> G
M -> O
N -> G
P -> D
```

Projections onto the y-axis:

```
E -> H
F -> J
I -> K
M -> A
N -> A
P -> D
```

---

#### Q4b

This question is similar to the previous one, but now we're showing their projections onto the two principle components.

```{r, fig.width=5, fig.height=5, echo=FALSE}
# Don't worry about undertanding this code chunk.
par(mar=c(4,4,1,1))
rotate <- function(xy, angle) {
  xy %*% matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), ncol=2)
}
sel <- xPoints != 0 & yPoints != 0
xy <- cbind(xPoints, yPoints)
offRes <- prcomp(xy[sel,], center=FALSE)
projections <- t(cbind(t(offRes$x) * c(1,0), t(offRes$x) * c(0,-1)))
#projections <- t(cbind(t(offRes$x) * c(1,0), t(offRes$x) * c(0,1)))
pc1 <- offRes$rotation[,"PC1"]
pc2 <- offRes$rotation[,"PC2"]
a <- acos((t(pc1) %*% c(1,0))[1,1])
rotatedProjections <- rotate(projections, -a)
plot(c(-h, h), c(-h, h), xlab="", ylab="", bty="n", axes=FALSE, type="n")
abline(h=0, col="gray", lty=2)
abline(v=0, col="gray", lty=2)
abline(b=pc1[2]/pc1[1], a=0, col="red")
abline(b=pc2[2]/pc2[1], a=0, col="blue")
points(xPoints[sel], yPoints[sel], cex=4.0, pch=20, col="black")
text(xPoints[sel], yPoints[sel], labels=LETTERS[1:nrow(xy)][sel], col="white")
points(rotatedProjections[1:6,1], rotatedProjections[1:6,2], pch=20, col="red")
text(rotatedProjections[1:6,1]-0.3, rotatedProjections[1:6,2]-0.1,
     labels=rev(LETTERS)[1:6], col="red", cex=0.7)
o <- order(rotatedProjections[7:12,1])
points(rotatedProjections[7:12,1][o], rotatedProjections[7:12,2][o], pch=20, col="blue")
text(rotatedProjections[7:12,1][o], rotatedProjections[7:12,2][o]+c(-0.3,0.3), 
     labels=rev(LETTERS)[7:12], col="blue", cex=0.7)

```

**For each of the 6 black points, make a list giving the letters of the points projected onto the first PC (red) and make another list of the 6 black points giving the letters of their projections onto second PC (blue).**

**A4b:**

Use the lists below and fill them by placing a letter after each arrow.

Projections onto PC1 (red):

```
E -> z
F -> y
I -> w
M -> w
N -> n
P -> u
```

Projections onto PC2 (blue):

```
E -> s
F -> r
I -> o
M -> t
N -> p
P -> r
```

---

#### Q5

In this problem we'll use the projection we found in Applications Lesson 6 of Obama and Trump speeches and apply counts from three George W. Bush speeches to see where they land.

First we'll need to read in the projection itself. I've created a file, `speeches.rds`, in the folder for this problem set. This file contains the `prcomp` object that was returned from the `prcomp()` method in the analysis of Trump and Obama speeches. The file is written in an R-specific data format called RDS.

We can read it in as follows:

```{r}
politics <- readRDS("speeches.rds")
class(politics)
colnames(politics$x)
politics$x[,"PC1"]

fullCounts <- read.table(file="presidential-speech-word-counts.tsv",
 header=TRUE, row.names=1, sep="\t") 

n <- ncol(fullCounts)
q <- nrow(fullCounts) 

totalCounts <- sapply(1:q, function(i) sum(fullCounts[i,])) 

speechesContainingWord <- sapply(1:q, function(i) {
 sum(fullCounts[i,] >= 1)
}) 

wordCounts <- fullCounts[speechesContainingWord >= 4 & totalCounts <= 500,]


label <- factor(c("trump", "obama")[grepl("obama", colnames(wordCounts))+1])

pcPlot <- function(xPC, yPC) {
 plot(politics$x[,xPC], politics$x[,yPC], col=c("blue", "red")[unclass(label)],
 xlab=xPC, ylab=yPC, pch=20)
}


pcPlot("PC2", "PC3") # Right plot
```

I've included word counts from three of President George W. Bush's speeches, which we can read in as follows:

```{r}
bush <- read.table(file="bush-word-counts.tsv", header=TRUE, stringsAsFactors=FALSE,
                   row.names=1)
head(bush)
```

R has a **generic function** called `predict()` that uses an existing trained model of some kind to spit out predictions for a new dataset. There are many types of models that `predict()` can be used with:

```{r}
methods(predict)
```

We can call predict, to get the coordinates of the new speech, as projected, onto the 20 PCs:

```{r}
bushProj <- predict(politics, newdata=t(bush))
print(bushProj)


```

**Plot the 10 Trump and 10 Obama speeches as projected onto PC2 and PC3, as in Applications Lesson 6, with red and blue dots. Then add the three speeches from President Bush to the plot using the `points()` function, as a single point. Make sure we can tell which point are the Bush speeches. Discuss the location of the points relative to the Obama and Trump speeches and what it might mean.**

Hint: You'll likely want to borrow code from the Applications 6 notes for parts of this question.

**A5:**

The locations of the Bush points are closer to the Obama points. Indicating they use similar words. 
```{r}

label <- factor(c("trump", "obama")[grepl("obama", colnames(wordCounts))+1])

pcPlot <- function(xPC, yPC) {
 plot(politics$x[,xPC], politics$x[,yPC], col=c("blue", "red")[unclass(label)],
 xlab=xPC, ylab=yPC, pch=20)
points(bushProj[,xPC], bushProj[,yPC], col=c("green")[unclass(label)])
}

pcPlot("PC2", "PC3") # Right plot
```

---

#### Q6

I've included a data set giving a 3D representation of an object.

It can be read in as follows:

```{r}
obj <- read.table(file="mystery-object.tsv", header=TRUE, sep="\t")
head(obj)
```

This point cloud representation includes more than 40k points:

```{r}
nrow(obj)
```

The object is three dimensional, but we can visualize each possible pair of dimensions.

Here is the x vs y scatter plot:

```{r, fig.width=8, fig.height=4, echo=-1}
par(mar=c(4,4,1,1), mgp=c(2.2,0.8,0))
plot(obj$x, obj$y, pch=20, cex=0.5, col=rgb(0,0,0,0.2), xlab="X", ylab="Y")
```

Here is the x vs z scatter plot:

```{r, fig.width=5, fig.height=4, echo=-1}
par(mar=c(4,4,1,1), mgp=c(2.2,0.8,0))
plot(obj$x, obj$z, pch=20, cex=0.5, col=rgb(0,0,0,0.2), xlab="Y", ylab="Z")
```

And this is the y vs z view:

```{r, fig.width=8, fig.height=4, echo=-1}
par(mar=c(4,4,1,1), mgp=c(2.2,0.8,0))
plot(obj$y, obj$z, pch=20, cex=0.5, col=rgb(0,0,0,0.2), xlab="X", ylab="Z")
```

**Perform PCA on the object matrix and make three plots: (1) PC1 vs PC2, and (2) PC1 vs PC3, and (3) PC2 vs PC3. See if you can explain what is shown in each plot. Based of these projections, can you better tell what object this is? Why do you think we get these particular three views of this object?**

Hint: You'll need to think of whether to use the original matrix or the transpose, t(), of the matrix. Remember, in this problem each object is a point with three features (x, y, and z coordinates). And you want the columns to be features and the rows to be objects.

**A6:**
The object is a car and we get these these different views because of the different dimensions of the object. 

```{r}

head(obj)
pcn <- prcomp(obj)
pcn

#dim(pcn$x)
#colnames(pcn$x)
pcPlot <- function(xPC, yPC) {
 plot(pcn$x[,xPC], pcn$x[,yPC], col=c("blue", "red"),
 xlab=xPC, ylab=yPC, pch=20)
}

pcPlot("PC1", "PC2")
pcPlot("PC1", "PC3")
pcPlot("PC2", "PC3")


```


---



