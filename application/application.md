---
title: Fay-Herriot Model Estimates for an Income and Living Condition Survey
subject: Notes for short course on small area estimation
author: Tobias Schoch (tobias.schoch@fhnw.ch)
version: February 12, 2022
---

#### Application with R

# Fay-Herriot Model Estimates for an Income and Living Condition Survey

Tobias Schoch (February 12, 2022)

[Morales et al.](#References) (2021, Chapter 1.2) have prepared data on the basis of a living conditions survey (LCS).[^(1)^](#Notes) The goal is to estimate average income for 26 small areas.

First, we introduce the LCS survey data. Then, we study the direct estimator of average income by area. Finally, we estimate the Fay-Herriot model.  Along our discussion, you will be asked to solve tasks.

* **Task 1.** Compute the direct estimator of average income;
* **Task 2.** Compute the generalized variance function of the direct estimator's variance (this will be discussed below);
* **Task 3.** Compute estimates of average income by area with Fay-Herriot model (incl. model diagnostics, and estimation of mean square error).

We will come back to the tasks as we go along.

## 1 Survey Data

#### 1.1 LCS Survey

The survey data `datLCS` contains 6 variables, which are measured for individuals living in private households. The households are identified by variable `house`. In total, the dataset contains data on 2512 individuals living in 962 households. The sampling weights `w` refers to the households. There are 26 small areas (identified by the domain indicator `dom` ). The variables of `datLCS` are summarized as follows.

| Variable | Description                                                  |
| -------- | ------------------------------------------------------------ |
| `sex`    | man=1; woman=2                                               |
| `house`  | household identifier                                         |
| `income` | net equivalent income in euros                               |
| `lab`    | labor status (0=child, i.e., age < 16; 1=employed; 2=unemployed; 3=inactive) |
| `dom`    | domain indicator (small area)                                |
| `w`      | sampling weight (level: household; calibrated)               |

The data are stored as `datLCS.txt` file on the https://github.com/tobiasschoch/SAEcourse GitHub repo (and are available in the in zip-archive). We can load the data by

```R
datLCS <- read.table("datLCS.txt", header = TRUE, sep = "\t", dec = ",")
```

The first three lines of `datLCS` are printed below.

```
dom  sex  house         w    income  lab
 27    1     68  3022.840   6262.40    3
 27    2     68  3022.840   6262.40    3
 27    1     68  3022.840   6262.40    3
```

#### 1.2 Direct Estimator

Average income in area $i=1,\ldots,26$ is computed with the Hajek estimator, which is defined as
$$
\begin{equation*}
\bar{y}_i = \frac{1}{\widehat{N}_i}\sum_{j \in s_i} w_j \cdot \mathrm{income}_j, \qquad \text{with} \quad \widehat{N}_i = \sum_{j \in s_i}w_j,
\end{equation*}
$$
where $s_i$ is the part of the sample $s$ that falls into the $i$-th area, and $w_i$ denotes the sampling weight.

There are several equivalent ways to compute the Hajek estimator of income by small area (`dom`) with the R  software ([R Core Team](#References), 2022).[^(2)^](#Notes) We want to stick with the functions of the R `base` package.  First, we split the `datLCS` data into a list by `dom` . 

```R
datLCS_dom <- split(datLCS, datLCS$dom)
```

The object `datLCS_dom` is a list with 26 list entries (one for each area). The list entries consist of the area-specific part of the data `datLCS`.

In the next step, we use `sapply()` to compute the Hajek estimator (i.e., weighted mean) by area.

```R
sapply(datLCS_dom, function(u) weighted.mean(u$income, u$w))
```

A more complete function of the Hajek estimator (which is also capable of computing an approximate variance of the estimator) is given by

```R
hajek <- function(x, w)
{
    avg <- weighted.mean(x, w)
    Nhat <- sum(w)
    var <- sum(w * (w - 1) * (x - avg)^2) / Nhat^2
    c(avg = avg, var = var, Nhat = Nhat)
}
```

Note that function `hajek()` also computes `Nhat`, i.e., the estimated population size $\widehat{N}_i$.

**Task 1. ** Use the `hajek()` function to compute the Hajek estimator and its variance for all areas. Also, compute the coefficient of variation (in %), defined as
$$
\begin{equation*}
cv_i = 100 \cdot \frac{v_i}{\bar{y}_i},
\end{equation*}
$$
for all $i=1,\ldots,n$ areas, where $v_i$ denotes the variance of the Hajek estimator.

## 2 Auxiliary Data

The `auxLCS` contains aggregated auxiliary data for all $i=1,\ldots,26$ areas. The first three observations are printed below

```
dom     TOT      Mwork     Mnowork     Minact         ss
  3   82001  0.3632226  0.12764258  0.3574135  0.5695195
  5  251866  0.3564652  0.15503770  0.3192122  0.4323160
  6  190653  0.3405221  0.15860230  0.3158246  0.5998553
```

where 

* `TOT`:  total number of individuals in area,

+ `Mwork`: domain mean of `lab=1`,
+ `Mnowork`: domain mean of `lab=2`,
+ `Minact`: domain mean of `lab=3`.

We will utilize the auxiliary information as explanatory variables in the Fay-Herriot model. The dataset `auxLCS` have been processed by [Morales et al.](#References) (2021) and are ready to use. In practice, we have to take care of this process ourselves.

## 3 Fay-Herriot Model

#### 3.1 Preparations

Before ; *generalized variance functions*. This is method attempts to compute the variances of an estimator by exploiting a simple mathematical relationship connecting the variance to the expectation of the estimator. To be explicit, let's consider estimating a population characteristics (e.g., mean). Let $\widehat{\theta}$ be an unbiased survey estimator of the population parameter $\theta$ with variance $v(\widehat{\theta})$. We define the relative variance 
$$
\Delta^2 = \frac{v(\widehat{\theta})}{\theta^2}\tag{A}
$$
which is equal to the squared coefficient of variation. In place of (A), we use a model
$$
\begin{equation*}
\Delta^2 = \alpha + \frac{\beta}{\theta}
\end{equation*}
$$

see [Wolter](#References) (2007, Chapter 7).

#### 3.2 Estimation of the Fay-Herriot model



## Notes

^(1)^ The LCS data synthetically generated data that imitate the structure of an income an living condition survey; see [Morales et al.](#References) (2021, p. 4).

^(2)^ Instead of the functions in the R `base` package, we may use R packages `data.table` or `tidyr` to compute the area-specific estimates.

## References

Morales, Esteban, Pérez & Hobza (2021). *A Course in Small Area Estimation and Mixed Models: Methods, Theory and Applications in R*, Cham: Springer Nature.

Wolter (2007). *Introduction to Variance Estimation*, New York: Springer.

R Core Team (2022). *R: A language and environment for statistical computing*. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

