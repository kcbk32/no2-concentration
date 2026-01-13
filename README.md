# Bayesian Time-Series Analysis of NO₂ Concentrations

This project analyzes long-term trends in atmospheric nitrogen dioxide (NO₂) using
Bayesian time-series models. NO₂ is a short-lived air pollutant primarily emitted by
traffic and industry and is widely used as an indicator of air quality and emission
policy effectiveness.

The goal of the analysis is to quantify long-term change in NO₂ concentrations while
correctly accounting for temporal dependence and uncertainty.

---

## What was done

Daily NO₂ measurements from a European monitoring station (1997–2025) were modeled
using three Bayesian models of increasing temporal complexity:

1. **Trend + seasonality (independent errors)**  
   Captures long-term change and annual seasonal patterns, assuming observations are
   otherwise independent.

2. **Trend + seasonality with AR(1) errors**  
   Accounts for short-term persistence in daily deviations, a common feature of
   environmental time series.

3. **Trend + seasonality with seasonal autoregressive errors**  
   Allows stochastic dependence across annual lags, beyond deterministic seasonality.

All models were fitted using Hamiltonian Monte Carlo and evaluated using
leave-one-out cross-validation (PSIS-LOO).

---

## Key findings

- **Short-term temporal dependence is essential.**  
  The AR(1) model substantially outperforms both the independent and seasonal
  autoregressive alternatives in predictive accuracy.

- **Additional seasonal stochastic complexity is not supported.**  
  The seasonal autoregressive model is more flexible but performs worse out of sample.

- **NO₂ concentrations show a clear long-term decline.**  
  The best-performing model estimates a persistent annual decrease of roughly
  **6–7% per year**, with uncertainty fully quantified via posterior distributions.

---

## Why this matters

This analysis demonstrates how Bayesian modeling can:
- Separate long-term trends from short-term variability
- Provide principled uncertainty estimates
- Avoid overfitting through predictive model comparison

The same workflow applies broadly to time-series problems in environmental science,
economics, operations, and applied data science.

---

## Repository contents

- `report/` – Full written analysis (PDF)
- `R/` – Data preparation, model fitting, evaluation, and plotting scripts
- `stan/` – Stan model definitions
- `output/` – Generated results and figures (not tracked in version control)

---

## Tools used

- **R**
- **Stan (HMC)**
- **PSIS-LOO cross-validation**
- **Bayesian time-series modeling**

---
