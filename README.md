# Multinomial logit modal choice models for aggregated data with Box-Cox transformations

## Introduction
This project illustrates the implementation of (conditional) multinomial
logit mode choice models for aggregated data with Box-Cox transform
of the explanatory variable(s).

## Input data
The model can be used when observations are not individual but aggregated. 
This is for instance the case when one observe transported quantities, 
such as tons, between origins and destinations. The provided sample dataset 
(sampleData.Rda) contains, for every origin-destination pair, the 
following data:
-	The ID of the origin ORG and the ID of the destination DST
-	The unit transportation cost between ORG and DST for three modes
-	The total transit time (duration) between ORG and DST for three modes
-	The observed transported quantities between ORG and DST for three modes.

The transportation modes are:
-	1 for Trucks (road network)
-	2 for Barges (inland waterways network)
-	3 for Trains (railroad network)

Note that all transportation modes are not necessarily available between 
all ORG-DST pairs. This is particularly the case for barges, as the 
inland waterways network is not available everywhere. 

## Objective of the models
Given the input data, the proposed multinomial logit models try to 
estimate the quantities transported between every ORG-DST pair using one 
or two explanatory variables, which are Box-Cox transformed.

## Used tools
The models are written in [R](https://www.r-project.org) and use the 
[‘mlogit’](https://cran.r-project.org/web/packages/mlogit/index.html) 
and [‘mnlogit’](https://cran.r-project.org/web/packages/mnlogit/index.html) 
packages.

## Why Box-Cox transformations?
The Box-Cox transformation is a well-known 
[power transform](https://en.wikipedia.org/wiki/Power_transform#Box–Cox_transformation). 
In the framework of modal choice models, it can be useful in order to:
-	Improve the log-likelihood of the model.
-	Obtain the expected signs of the estimators when explanatory variables are correlated.

The Box-Cox transform is rather easy to apply once the “lambda” to use is 
known. Obtaining the optimal value(s) of the lambda(s) is, however, not 
immediate. 

Two models are provided in this project: the first one is an univariate 
model for which only one lambda must be determined. The second is a 
bivariate model, for which a combination of two lambdas’ (one for each 
explanatory variable) must be found.

## Model 1: The univariate case
The "Logit1.R" script illustrates a simple univariate model that uses the
unit transportation cost (or transit time) as only explanatory variable.

"BoxCoxLogit1.R" uses the same explanatory variable, but applies a Box-Cox
transform to it. The aim is here  to increase the log-likelihood of 
the model using an “optimal” lambda. Therefore, all the values in 
the [-2,2] range are tested, with a step of 0.1.

This model is described and used with real datasets in 
[**Jourquin B.**, *Estimating Elasticities for Freight Transport Using a 
Network Model: An Applied Methodological Framework*, Journal of 
Transportation Technologies, 2019, 9, 1-13](https://doi.org/10.4236/jtts.2019.91001).

## Model 2: The bivariate case
The "Logit2.R" script illustrates the case of a model with two explanatory 
variables: the unit transportation cost and the transit time (duration).
It shows that, without a propre transform of the variables, the sign of
an estimator can be wrong. In this case, the estimator for the "duration"
variable is positive. This is not expected as the modal share of a 
transportation mode should decrease when its transit time increases. This
unexpected sign is related to the fact that both explanatory variables
are heavily correlated.

"BoxCoxLogit2.R" applies a Box-Cox transform, looking for the "best" 
lambda's, not only to improve the log-likelihood of the model, but also 
to obtain the expected signs for the estimators. 

This model is thoroughly described and discussed when used with real 
World datasets in [**Jourquin B. and Beuthe M.**, *Cost, transit time 
and speed elasticity calculations for the European continental freight 
transport*, Transport Policy, 83, 1-12, 2019](https://doi.org/10.1016/j.tranpol.2019.08.009).

## Output
Both scripts print the estimated model (values of the estimators, t-values,
log-likelihood...), the 100 first rows of the input data along with the
related estimated quantities and the correlations between observed and 
estimated quantities.

## Use of the estimators in the Nodus software
[Nodus](http://nodus.uclouvain.be) is a transportation network modeling 
software especially designed for multimodal and intermodal freight 
transport. It offers the possibility to develop user-defined modal 
choice plugins. It can, for instance, use the estimated parameters by 
the two models presented in this project.

To start with, the user has to gather the input data needed to estimate 
the parameters (a dataset which structure is similar to the one of 
“sampleData.Rda”, but specific to the case to handle). The unit 
transportation costs and transit times per mode can be obtained by means of 
an assignment with “save paths” checked (“detailed paths” is not needed)
in the Nodus "assignment" panel. 
In the “header” table resulting from this assignment, one can find, for 
each OD pair (and each “group of commodities”) the unit transportation 
costs and transit times (note that these are expressed in seconds and 
not hours). The observed quantities per mode are to be gathered from the 
(modal) OD matrixes that are used.  

Once the input data prepared, the R script can be run to obtain the 
estimators. If "printNodusEstimators" is set to TRUE in the R script,
the estimators are printed in a format that can directly be 
cut&pasted in a Nodus “.costs” file. Example:

    lambda.0=0.7
    (Intercept).2.0=-4.14188571102275
    (Intercept).3.0=-4.34237545509826
    cost.1.0=-0.709553978510067
    cost.2.0=-0.709553978510067
    cost.3.0=-0.709553978510067

The "plugin" directory contains two user-defined modal choice 
plugins for Nodus that can make use of these parameters. 

## Further improvements
The strategies used to find the optimal lambda's in the two models needs 
to be improved as the proposed approaches are simple but time consuming.


