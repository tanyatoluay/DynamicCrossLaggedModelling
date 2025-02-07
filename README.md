# Cross-Lagged Dynamic Bayesian Modelling in Panel Data 

## Introduction
This study examines the relationship between alcohol consumption and depression using simulated data. The analysis applies statistical inference and probabilistic modeling.

## Methodology

### 1. Data Simulation
To ensure controlled conditions, synthetic data is generated:
- **Alcohol Consumption (
$X$)**: Continuous variable modeled as a normal distribution.
- **Depression Score (
$Y$)**: Continuous response variable influenced by \(X\) with added noise.


Mathematically:
X ~ N(μ_X, σ_X²),   Y = β₀ + β₁X + ε,   ε ~ N(0, σ_ε²)


### 2. Statistical Model
We apply linear regression to estimate the effect of alcohol on depression:
$\hat{Y} = \hat{\beta}_0 + \hat{\beta}_1 X$ where $\hat{\beta}_1$ indicates 
the change in depression score per unit increase in alcohol consumption.

### 3. Hypothesis Testing
We test:
- **Null Hypothesis (
$H_0$)**: No significant relationship (
$\beta_1 = 0$)
- **Alternative Hypothesis (
$H_A$)**: A significant relationship exists (
$\beta_1 \neq 0$)

A t-test on $\hat{\beta}_1$ determines statistical significance.

### 4. Bayesian Approach
A Bayesian regression is also implemented:

$$ P(\beta | X, Y) \propto P(Y | X, \beta) P(\beta) $$

where priors are assigned to \(\beta\).

### 5. Model Evaluation
- **Goodness of Fit**: R-squared (
$R^2$) and Adjusted R-squared.
- **Residual Analysis**: Normality and homoscedasticity checked via diagnostics.
- **KL Divergence**: Comparing empirical vs. theoretical distributions.
