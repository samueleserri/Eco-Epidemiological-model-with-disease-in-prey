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

## Repository layout
- src/

  - EcoEpiModel.py — abstract base class for eco epidemiological model with parameters:
    - $r$: intrinsic birth rate of the prey population.
    - $K$: carrying capacity for the prey population.
    - $\beta$: transmission coefficient of the disease.
    - $c$: death rate of the infected prey.
    - $b$: predation coefficient.
    - $m$: ratio-dependent rate.
    - $k$: coefficient in conversing prey into predator offspring.
    - $d$: death rate of the predator.

  - DiscreteTimeModel.py — implementation of the discrete model along with plottings methods.
  - ContinuousTimeModel.py — Continuous version of the model. Not used in numerical simulations.
  Reference : Xiao, Y, Chen, L: A ratio-dependent predator-prey model with disease in the prey. Appl. Math. Comput. 131, 397-414 (2002) 
  - simulations.ipynb — notebook with numerical simulations and visualizations 


- slides/ — Beamer presentation sources


- requirements.txt — Python dependencies

Minimal Python example:
```python
from DiscreteTimeModel import DiscreteTimeModel

params = {"r":0.2,"b":0.1,"c":1.5,"d":0.2,"k":0.2,"m":0.5,"beta":0.2,"K":10}
model = DiscreteTimeModel(params)
init = {"S0":4.0,"I0":0.5,"Y0":0.1}
res = model.run(init, max_iter=1000)
model.plot_phase_space(init, res=res)
```