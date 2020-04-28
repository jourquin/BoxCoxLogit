# Copyright (c) 2019 Université catholique de Louvain Center for Operations Research and Econometrics (CORE) http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script showing how to implement an bivariate (conditional) multinomial logit when the following input data is available for an
# origin-destination pair:
# - Origin and a destinataion ID's
# - Unit cost (or duration) for each mode on the OD relation
# - Observed demand (tons for instance) for each mode between O and D
#
# The script uses the obtained parameters for the logit model to estimate the demand (tons) for each mode on each OD relation.
#
# The provided sample dataset is related to freight transport, with 3 modes (1: road, 2: inland waterways, 3: rail). Not all modes
# are available between all OD pairs (NA values in the dataframe). Unit costs are expressed in euros/ton and duration in hours,
# including time needed to load/unload...
#
# Using the provided data set, the sign of the "duration" variable is positive. This is not the expected sign, as the modal share
# of a transportation mode should decrease when its transit time increases. This strange result is due to the fact that the "cost"
# and "duration" variables are (heavily) correlated. A Box-Cox transformation of these variables can help to fix this, as
# illustrated by the BoxCoxLogit2.R script.

library(mlogit)
library(mnlogit)


# The mlogit package can't cope with 0 values
smallQty <- .Machine$double.xmin

# Solves the univariate conditional logit with Box-Cox transform
solveLogit <- function(x) {
  # Replace null quantities with a small one
  x$qty.1[x$qty.1 == 0] <- smallQty
  x$qty.2[x$qty.2 == 0] <- smallQty
  x$qty.3[x$qty.3 == 0] <- smallQty

  # Replace missing costs with an high value
  highValue <- max(x$cost.1, x$cost.2, x$cost.3, na.rm = TRUE) * 1000
  x$cost.1[is.na(x$cost.1)] <- highValue
  x$cost.2[is.na(x$cost.2)] <- highValue
  x$cost.3[is.na(x$cost.3)] <- highValue

  # Replace missing durations with an high value
  highValue <- max(x$duration.1, x$duration.2, x$duration.3, na.rm = TRUE) * 1000
  x$duration.1[is.na(x$duration.1)] <- highValue
  x$duration.2[is.na(x$duration.2)] <- highValue
  x$duration.3[is.na(x$duration.3)] <- highValue

  # Formula for a conditional multinomial logit. (see mlogit package)
  f <- mFormula(mode ~ cost + duration | 1 | 1)

  # Total transported quantity
  x$totQty <- x$qty.1 + x$qty.2 + x$qty.3

  # Create wideData data, with one record per mode for each OD pair (see mlogit documentation)
  wideData <- data.frame()
  for (mode in 1:3) {
    wd <- data.frame(mode = integer(nrow(x)))
    wd$mode <- mode

    wd$cost.1 <- x$cost.1
    wd$cost.2 <- x$cost.2
    wd$cost.3 <- x$cost.3

    wd$duration.1 <- x$duration.1
    wd$duration.2 <- x$duration.2
    wd$duration.3 <- x$duration.3

    wd$qty <- x$qty.1
    if (mode == 2) {
      wd$qty <- x$qty.2
    } else if (mode == 3) {
      wd$qty <- x$qty.3
    }

    wideData <- rbind(wideData, wd)
  }

  # Transform into "long" format data (see mlogit documentation)
  longData <-
    mlogit.data(wideData,
      choice = "mode",
      shape = "wide",
      varying = 2:7
    ) # First column is "mode", variables are in columns 2 to 7.


  # Solve the model, using the mnlogit package, faster (parallelized) that mlogit
  nbCores <- parallel:::detectCores()
  model <- mnlogit(
    f,
    longData,
    weights = wideData$qty, # This is a weighted logit
    na.rm = FALSE,
    ncores = nbCores
  )

  return(list("data" = x, "model" = model))
}

######################################################################
# Real entry point  ##################################################
######################################################################

# Change working directory to the location of this script
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Change the output width
options(width = 200)

load("sampleData.Rda")

r <- solveLogit(sampleData)
model <- r$model
df <- r$data

# Print the resulting model
print(summary(model)) # !!! With the provided data, the estimator for "duration" is positive !!!

# Restore null quantities
df$qty.1[df$qty.1 == smallQty] <- 0
df$qty.2[df$qty.2 == smallQty] <- 0
df$qty.3[df$qty.3 == smallQty] <- 0

#####################################################################
# 3) Compute the estimated quantities using the estimated logit model
#####################################################################

# Numerators of the multinomial logit

df$utility.1 <- coef(model)["cost"] * df$cost.1 + coef(model)["duration"] * df$duration.1
df$utility.2 <- coef(model)["cost"] * df$cost.2 + coef(model)["duration"] * df$duration.2 + coef(model)["(Intercept):2"]
df$utility.3 <- coef(model)["cost"] * df$cost.3 + coef(model)["duration"] * df$duration.3 + coef(model)["(Intercept):3"]

# Denominator of the multinomial logit
df$denominator <- exp(df$utility.1) + exp(df$utility.2) + exp(df$utility.3)

# Compute the estimated quantities for each mode applying the multinomial logit
df$est_qty.1 <- round(df$totQty * exp(df$utility.1) / df$denominator)
df$est_qty.2 <- round(df$totQty * exp(df$utility.2) / df$denominator)
df$est_qty.3 <- round(df$totQty * exp(df$utility.3) / df$denominator)

# View the estimated tonnages
cat("\nFirst rows of input data and estimations:\n")
print(head(df, 50))

# Compute a simple correlation between observed and estimated quantities
cat("\nCorrelations between observed and estimated quantities: ")
cat(paste("Road:", round(cor(df$qty.1, df$est_qty.1), 2), ", ", sep = ""))
cat(paste("IWW:", round(cor(df$qty.2, df$est_qty.2), 2), ", ", sep = ""))
cat(paste("Rail:", round(cor(df$qty.3, df$est_qty.3), 2), "\n\n"))
