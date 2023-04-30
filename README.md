# Estimating time-varying treatment effects in longitudinal studies

This manuscript has been accepted for publication in *Psychological Methods*.

Preprint available at: https://psyarxiv.com/8sc74

A tutorial on implementing g-estimation of an SNMM in lavaan using a simple example has been accepted for publication in *Advances in Methods and Practices in Psychological Science*.

Preprint available at: https://psyarxiv.com/f23gj/

**Abstract**

Longitudinal study designs are frequently used to investigate the effects of a naturally observed predictor (treatment) on an outcome over time. Because the treatment at each time point or wave is not randomly assigned, valid inferences of its causal effects require adjusting for covariates that confound each treatment-outcome association. But adjusting for covariates which are inevitably time-varying is fraught with difficulties. On the one hand, standard regression adjustment for variables affected by treatment can lead to severe bias. On the other hand, omitting time-varying covariates from confounding adjustment precipitates spurious associations that can lead to severe bias. Thus, either including or omitting time-varying covariates for confounding adjustment can lead to incorrect inferences. In this paper, we introduce an estimation strategy from the causal inference literature for evaluating the causal effects of time-varying treatments in the presence of time-varying confounding. G-estimation of the treatment effect at a particular wave proceeds by carefully adjusting for only pre-treatment instances of all variables while dispensing with any post-treatment instances. The introduced approach has various appealing features. Effect modification by time-varying covariates can be investigated using covariate-treatment interactions. Treatment may be either continuous or noncontinuous with any mean model permitted. Unbiased estimation requires correctly specifying a mean model for either the treatment or the outcome, but not necessarily both. The treatment and outcome models can be fitted with standard regression functions. In summary, g-estimation is effective, flexible, robust, and relatively straightforward to implement.

