## Synopsis: Timing an rwbs simulations. What takes the most time? Does the
## intervals() function take a long time? What actually goes on during the
## intervals that makes it take a long time?
library(microbenchmark)
source("../main/artificial/artif-helpers.R")
source("../main/artificial/artif-sim-rwbs-timing-helpers.R")

## How much time is spent?
## 1. Forming intervals
n=2000
before.time = microbenchmark({print("here"); intervals_old(n=n, numIntervals=n)}, times=5)
after.time = microbenchmark({print("here"); intervals_new(n=n, numIntervals=n)}, times=5)

## 2. Filling them in
set.seed(0)
old.time = microbenchmark({intv = intervals_old(n=n, numIntervals=n)}, times=1)
mn = rep(0,n)
mn[100:200] = mn[1100:1200] = mn[1800:2000] = lev
set.seed(99)
y = mn + rnorm(n, 0, sigma)

## What about adding criteria?
la("../../biomotors/biomotors")
intv = addcrit.intervals(y=y, intervals=intv)

## How much memory is spent?



