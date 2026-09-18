# Bayesian-Hierarchical-Claim-Frequency-Modeling-Prior-Sensitivity-B-hlmann-Credibility
Modeled regional auto insurance claim frequencies across 21 French administrative regions using a hierarchical Poisson GLM fit in Stan, varying sample size from N=500 to 676k policies to study when Bayesian partial pooling adds practical value over classical estimation.

Methods & Technical Work:
        • Fit a hierarchical Poisson GLM with region-level random effects and non-centered parameterization; compared weakly informative vs. informative Half-Normal priors on the regional variance component​
        • Demonstrated that regional variance is identified by the number of groups (21 regions), not total N — an informative prior effectively imposes complete pooling at small N regardless of individual policy count
        • Quantified shrinkage breakdown across sample sizes: at N=500 raw regional rates span 0–0.3, while posterior means collapse toward the global mean; at N=676k the model largely trusts each region's own data
        • Derived Bühlmann credibility equivalence analytically and verified numerically — credibility factors range from 0.30 (Franche-Comté, n=1,323) to 0.98 (Centre, n=160,529)
        • Showed that covariate adjustment separates what classical Bühlmann cannot: Île-de-France's raw rate of 0.132 drops to a posterior of 0.101, explained by driver demographics and density rather than genuine regional risk

Tools: R, Stan (cmdstanr), ggplot2, posterior, loo

Key Finding: Using 676k French motor insurance policies, the model found that regional claim rates vary meaningfully across France — but that a small regional insurer with 500 policies cannot reliably estimate this variation without strong prior assumptions, while a large national insurer with hundreds of thousands of policies needs no prior at all. Prior choice is a practical pricing decision, not a statistical formality.
