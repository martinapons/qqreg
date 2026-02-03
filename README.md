# qqreg

## Install the package

```r
# Option 1: remotes
install.packages("remotes")
remotes::install_github("martinapons/qqreg")

# Option 2: devtools
install.packages("devtools")
devtools::install_github("martinapons/qqreg")
```

## Description

qqreg implements the Quantile on Quantiles (QQ) estimator suggested in [Pons (2026)](https://martinapons.github.io/files/QQmodel.pdf), a two-step estimator that recovers a two-dimensional quantile surface.  
The surface characterizes heterogeneity jointly along two dimensions:

- **Within-group rank** (e.g., individuals within a region/firm/cohort), and
- **Between-group rank** (e.g., the relative standing of regions/firms/cohorts).

The main workflow is:
1. Fit the model with `qq_fit()` on a user-specified quantile grid.
2. Inspect results with `summary()`.
3. Visualize with `plot()` as surfaces or slices, returning `ggplot2` objects that can be further customized. To see available options, type `?plot.qqfit`.


## Examples

Getting Started
```r
# load package
library(qqreg)
# help files
?qqreg
?qq_fit
?plot.qqfit
```

Generate data with heterogeneous effects
```r
set.seed(123)
n_groups <- 100
n_per_group <- 100
n <- n_groups * n_per_group

# Covariates
x1 <- 1 + rnorm(n) + rep(runif(n_groups), each = n_per_group)
x2 <- rep(rnorm(n_groups), each = n_per_group)

# Shocks: eta (group-level), nu (individual-level)
eta <- rep(rnorm(n_groups), each = n_per_group)
nu <- rnorm(n)

# DGP: y = 1 + x1 + x2 + eta*(1 - 0.1*x1 - 0.1*x2) + nu*(1 + 0.1*x1 + 0.1*x2)
y <- 1 + x1 + x2 + eta * (1 - 0.1*x1 - 0.1*x2) + nu * (1 + 0.1*x1 + 0.1*x2)

data <- data.frame(
  group = rep(1:n_groups, each = n_per_group),
  x1 = x1,
  x2 = x2,
  y = y
)
```

### QQ Model Without Covariates
Estimate the model
```r
# Estimate the model for quantiles 0.1, 0.2,..., 0.9. Estimation runs in parallel by default. 
fit1 <- qq_fit(y ~1, data = data, group = "group",
              taus = seq(0.1, 0.9, 0.1), nboot = 200)

# Print the results
summary(fit1)
```

```text
Quantile on Quantiles Regression Summary
========================================

Call:
qq_fit(formula = y ~ 1, data = data, group = "group", taus = seq(0.1, 
    0.9, 0.1), nboot = 200, parallel = T, parallel_first = T)

Observations: 10000 
Groups: 100 

---
Coefficient: (Intercept) 
---

Point Estimates:
        v=0.1   v=0.2  v=0.3  v=0.4  v=0.5  v=0.6  v=0.7  v=0.8  v=0.9
u=0.1 -0.4953 -0.2300 0.1479 0.2891 0.5874 0.8916 1.2888 1.5478 1.8251
u=0.2 -0.0285  0.3871 0.6144 0.8745 1.1927 1.5002 1.9356 2.1452 2.4727
u=0.3  0.4528  0.9142 1.1392 1.3521 1.7357 1.9993 2.3025 2.6641 2.9950
u=0.4  0.8087  1.3241 1.5062 1.7467 2.1021 2.3704 2.6349 3.0474 3.3767
u=0.5  1.1670  1.5910 1.8590 2.1152 2.4068 2.7266 3.0439 3.3531 3.7477
u=0.6  1.5031  1.9695 2.1967 2.4905 2.7913 3.0956 3.4280 3.7258 4.0054
u=0.7  1.9642  2.2230 2.6630 2.9284 3.2386 3.4824 3.8652 4.1880 4.5327
u=0.8  2.4994  2.7086 3.1678 3.4037 3.7044 4.0393 4.4182 4.6660 5.0111
u=0.9  3.0441  3.4294 3.8588 4.0517 4.4024 4.7171 5.0500 5.5409 5.8252

[...]
```
Surface Plot
```r
# Default settings
plot(fit1, which = "(Intercept)", type = "surface")
```
<img width="731" height="608" alt="image" src="https://github.com/user-attachments/assets/bcbc02b8-ce8e-4bef-8267-42ca3cf2030c" />

```r
# Personalize the plot (change axis labels, title, and colors)
myColors <- c("#7e227fff", "#174aa9ff", "#1ea861ff", "#d8e400ff")
plot(fit1, which = "(Intercept)", type = "surface", xlab = "Between Dimension",
ylab = "Within Dimension", ggtitle = "Two-Dimensional Quantile Function of Y") +
  ggplot2::scale_fill_gradientn(name = " ", colors = myColors)
```
<img width="730" height="608" alt="image" src="https://github.com/user-attachments/assets/74af917f-2f45-419f-ba5a-ca8940c221b8" />


### QQ Model With Covariates

Estimate the model (1000 bootstrap repetitions!) 
```r
fit2 <- qq_fit(y ~ x1 + x2, data = data, group = "group",
              taus = seq(0.1, 0.9, 0.1), nboot = 1000)
# Print the results
summary(fit2)
```
```text

Quantile on Quantiles Regression Summary
========================================

Call:
qq_fit(formula = y ~ x1 + x2, data = data, group = "group", taus = seq(0.1, 
    0.9, 0.1), nboot = 1000)

Observations: 10000 
Groups: 100 

---
Coefficient: (Intercept) 
---

Point Estimates:
        v=0.1   v=0.2   v=0.3   v=0.4   v=0.5  v=0.6  v=0.7  v=0.8  v=0.9
u=0.1 -1.3439 -0.8189 -0.6231 -0.3511 -0.0902 0.1454 0.3902 0.5758 0.9500
u=0.2 -0.8014 -0.3701 -0.1155  0.0773  0.3512 0.5918 0.8237 1.0686 1.3973
u=0.3 -0.4634 -0.0285  0.1758  0.4135  0.6513 0.9297 1.2118 1.3656 1.7742
u=0.4 -0.2244  0.2110  0.4620  0.6222  0.8535 1.1475 1.4803 1.6522 2.0376
u=0.5  0.0012  0.3944  0.7053  0.9281  1.1931 1.4004 1.7073 1.8250 2.2483
u=0.6  0.3399  0.6426  0.9646  1.1694  1.4594 1.7258 1.9406 2.1310 2.5551
u=0.7  0.6384  0.8986  1.2654  1.5032  1.7360 1.9387 2.3067 2.4563 2.8509
u=0.8  0.9364  1.2892  1.5735  1.8214  1.9928 2.2590 2.5398 2.8004 3.1708
u=0.9  1.2797  1.7733  2.0830  2.2303  2.4804 2.6616 2.8983 3.1625 3.6391
[...]
```

Surface Plot 
```r
# Add custom color palette
myColors <- c("#160417ff", "#b51fb0ff", "#f0aa5bff", "#fdebb1ff")
plot(fit2, which = "x1", type = "surface", ggtitle = "Coefficient on X1", 
    ylab = "Within Dimension", xlab = "Between Dimension") + 
    ggplot2::scale_fill_gradientn(name = " ", colors = myColors)
```
<img width="730" height="610" alt="image" src="https://github.com/user-attachments/assets/62a31f12-5520-48dd-9aa3-c8841cef79ff" />

Plot Slices of u
```r
# Plot slices with default settings. Get 1 graph for each value of u (within rank). 
plot(fit2, which = "x1", type = "slice_u")
```
<img width="987" height="609" alt="image" src="https://github.com/user-attachments/assets/ce9f50ea-7191-4dca-939d-e6dc4828cbcf" />

```r
# Plot two slices and force them on a common y-axis and common X label. 
plot(fit2, which = "x1", type = "slice_u",
     taus_slice = c(0.2, 0.8),  panel_prefix = "Within quantile: ", xlab ="Between Dimension",
     shared_ylim = TRUE, ggtitle = " ", line_color = "#ff9d15ff", ribbon_fill = "#f6941cff",
     ribbon_alpha = 0.3) 
```
<img width="986" height="604" alt="image" src="https://github.com/user-attachments/assets/8d104431-c43d-4c38-a4ab-19e626ff6bfe" />

```r

```

# References
Pons M. (2026). Quantiles on Quantles [Link to the paper](https://martinapons.github.io/files/QQmodel.pdf)
