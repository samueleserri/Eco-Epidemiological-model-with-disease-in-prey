# Eco-Epidemiological-model-with-disease-in-prey

This repository contains the files for the project of the course "Sistemi complessi" at the University of Trieste.

The repo includes the implementation of the eco-epidemiological model with disease in prey, as described in the paper:

Hu, Z., Teng, Z., Jia, C. et al. Complex dynamical behaviors in a discrete eco-epidemiological model with disease in prey. Adv Differ Equ 2014, 265 (2014).
https://doi.org/10.1186/1687-1847-2014-265

As well as numerical simulations and visualizations of the model's dynamics.

The model is described by the following system of difference equations:

```math
\begin{align}
S(t + 1) &= S(t) \exp\left\{r\left(1-\frac{S(t) + I(t)}{K}\right) - \beta I(t)\right\} \\
I(t + 1) &= I(t) \exp\left\{\beta S(t) - c - \frac{b Y(t)}{m Y(t) + I(t)}\right\} \\
Y(t + 1) &= Y(t) \exp\left\{\frac{k b I(t)}{m Y(t) + I(t)} - d\right\}
\end{align}
```
