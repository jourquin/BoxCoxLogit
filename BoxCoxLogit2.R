# Copyright (c) 2019 Université catholique de Louvain
# Center for Operations Research and Econometrics (CORE)
# http://www.uclouvain.be
# Written by Pr Bart Jourquin, bart.jourquin@uclouvain.be
#
# R script showing how to implement a bivariate (conditional) multinomial logit with a Box-Cox transform of the 2 explanatory
# variables when the following input data is available for an origin-destination pair:
# - Origin and a destinataion ID's
# - Unit cost (or duration) for each mode on the OD relation
# - Observed demand (tons for instance) for each mode between O and D
#
# The script looks for the optimal combination of lambda's to apply for the Box-Cox transformation and uses the obtained parameters
# for the logit model to estimate the demand (tons) for each mode on each OD relation.
#
# The provided sample dataset is related to freight transport, with 3 modes (1: road, 2: inland waterways, 3: rail). Not all modes
# are available between all OD pairs (NA values in the dataframe). Unit costs are expressed in euros/ton and duration in hours,
# including time needed to load/unload...
#

library(mlogit)
library(mnlogit)

# Print the estimators in a Nodus compatilbe format
printNodusEstimators <- FALSE

# The mlogit package can't cope with 0 values
smallQty <- .Machine$double.xmin

# Solves the bivariate conditional logit with Box-Cox transform
solveBoxCoxLogit <- function(x, lambdaCost, lambdaTime) {

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

  # Apply a Box-Cox transform to the two explanatory variables
  if (lambdaCost != 0) {
    x$cost.1 <- (x$cost.1^lambdaCost - 1) / lambdaCost
    x$cost.2 <- (x$cost.2^lambdaCost - 1) / lambdaCost
    x$cost.3 <- (x$cost.3^lambdaCost - 1) / lambdaCost
  }
  else {
    x$cost.1 <- log(x$cost.1)
    x$cost.2 <- log(x$cost.2)
    x$cost.3 <- log(x$cost.3)
  }

  if (lambdaTime != 0) {
    x$duration.1 <- (x$duration.1^lambdaTime - 1) / lambdaTime
    x$duration.2 <- (x$duration.2^lambdaTime - 1) / lambdaTime
    x$duration.3 <- (x$duration.3^lambdaTime - 1) / lambdaTime
  }
  else {
    x$duration.1 <- log(x$duration.1)
    x$duration.2 <- log(x$duration.2)
    x$duration.3 <- log(x$duration.3)
  }

  # Total transported quantity
  x$totQty <- x$qty.1 + x$qty.2 + x$qty.3

  # Formula for a conditional multinomial logit. (see mlogit package)
  f <- mFormula(mode ~ duration + cost | 1 | 1)

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

  # Create the "long" format data (see mlogit documentation)
  longData <-
    mlogit.data(wideData,
      choice = "mode",
      shape = "wide",
      varying = 2:7
    ) # First column is "mode", variables are in columns 2 to 7.

  # mlogit version 1.1 returns a dfidx object, but mnlogit only accepts data frames for now
  # Be sure to have a dataframe
  longData = as.data.frame(longData)

  # Solve the model, using the mnlogit package, faster (parallelized) that mlogit
  nbCores <- parallel:::detectCores()
  model <- mnlogit(
    f,
    choiceVar = "alt",
    longData,
    weights = wideData$qty, # This is a weighted logit
    na.rm = FALSE,
    ncores = nbCores
  )

  return(list("data" = x, "model" = model))
}

# Test if all the estimators are of the expected sign (both must be negative)
signsAreExpected <- function(model) {
  c <- coef(model)
  correctSign <- TRUE
  # Browse de coefficients names (see output of "summary(model)")
  for (j in 1:length(c)) {
    name <- names(c[j])
    if (substring(name, 1, 1) != "(") {
      # "(Intercept)" must not be tested
      if (c[name] > 0) {
        correctSign <- FALSE
        break
      }
    }
  }
  return(correctSign)
}

########################################################################################################## @
# Main entry point
########################################################################################################## @

# Change working directory to the location of this script
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

# Change the output width
options(width = 200)

###############################
# 1) Find the best lambda value
###############################

bestLL <- -10000000000.0
bestLambdaCost <- 10.0
bestLambdaTime <- 10.0

# Set range to search into
initialThreshold <- 2.4
minLambdaCost <- -initialThreshold
maxLambdaCost <- initialThreshold
minLambdaTime <- -initialThreshold
maxLambdaTime <- initialThreshold

# Initial step = 0.8
step <- 0.8
for (j in 1:4) {
  # Steps are 0.8, 0.4, 0.2 and 0.1
  lambdaCost <- minLambdaCost
  while (lambdaCost <= maxLambdaCost) {
    lambdaTime <- minLambdaTime
    while (lambdaTime <= maxLambdaTime) {
      cat(
        paste(
          "Testing for lambdaCost",
          lambdaCost,
          ", lambdaTime",
          lambdaTime,
          "and step",
          step,
          "\n"
        )
      )

      load("sampleData.Rda")

      # Some combinations of lambda's can lead to numerical singularity
      # This code intercepts the error and ignore it before running the
      # next loop
      res <- try(
        {
          r <- solveBoxCoxLogit(sampleData, lambdaCost, lambdaTime)

          if (r$model$logLik > bestLL) { # Better solution ?
            if (signsAreExpected(r$model)) { # And are the signs of coefficients expected ?
              bestLL <- r$model$logLik
              bestLambdaCost <- lambdaCost
              bestLambdaTime <- lambdaTime
            }
          }
        },
        silent = TRUE
      )
      if (inherits(res, "try-error")) {
        # Just ignore error
      }

      lambdaTime <- round(lambdaTime + step, 2)
    }
    lambdaCost <- round(lambdaCost + step, 2)
  }

  # Prepare next iteration
  if (bestLambdaCost > minLambdaCost) {
    minLambdaCost <- bestLambdaCost - step
  }
  if (bestLambdaCost < maxLambdaCost) {
    maxLambdaCost <- bestLambdaCost + step
  }
  if (bestLambdaTime > minLambdaTime) {
    minLambdaTime <- bestLambdaTime - step
  }
  if (bestLambdaTime < maxLambdaTime) {
    maxLambdaTime <- bestLambdaTime + step
  }

  # Decrease step size
  step <- step / 2
}


####################################################
# 2) Solve the model with the best or given lambda's
####################################################

# Retain these values to solve the logit
lambdaCost <- bestLambdaCost
lambdaTime <- bestLambdaTime

cat(paste("\n\nSolving with lambdaCost = ", lambdaCost, " and lambdaTime = ", lambdaTime))

# Solve the logit
load("sampleData.Rda")

r <- solveBoxCoxLogit(sampleData, lambdaCost, lambdaTime)
model <- r$model
df <- r$data

# Print the resulting model
print(summary(model))

# Restore null quantities
df$qty.1[df$qty.1 == smallQty] <- 0
df$qty.2[df$qty.2 == smallQty] <- 0
df$qty.3[df$qty.3 == smallQty] <- 0

if (printNodusEstimators) {
  # The estimated parameters can be used in a user defined modal-split method in Nodus (BoxCox1.java).
  # Therefore, copy&paste the following output into a project costs file, assuming that the estimated
  # values are for group 0.
  cat("The following lines can be pasted in a Nodus '.costs' file:\n")
  group <- 0
  cat(paste("lambda.", group, "=", lambda, "\n", sep = ""))
  c <- coef(model)
  for (j in 1:length(c)) {
    name <- names(c[j])
    mode <- substr(name, nchar(name), nchar(name))
    if (mode == "1" || mode == "2" || mode == "3") {
      mode <- paste(".", mode, ".", group, sep = "")
      name <- substr(name, 1, nchar(name) - 2)
      name <- paste(name, mode, "=", unname(c[j]), "\n", sep = "")
      cat(name)
    } else {
      # Conditional variables are replicated for the three modes
      mode <- paste(".1.", group, sep = "")
      n <- paste(name, mode, "=", unname(c[j]), "\n", sep = "")
      cat(n)
      mode <- paste(".2.", group, sep = "")
      n <- paste(name, mode, "=", unname(c[j]), "\n", sep = "")
      cat(n)
      mode <- paste(".3.", group, sep = "")
      n <- paste(name, mode, "=", unname(c[j]), "\n", sep = "")
      cat(n)
    }
  }
}

#####################################
# 3) Compute the estimated quantities
#####################################

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


# View the estimated tonnages for the first 50 rows
cat("\nFirst rows of input data and estimations. Note that costs and durations are Box-Cox transformed here:\n")
print(head(df, 50))

# Compute a simple correlation between observed and estimated quantities
cat("\nCorrelations between observed and estimated quantities: ")
cat(paste("Road:", round(cor(df$qty.1, df$est_qty.1), 2), ", ", sep = ""))
cat(paste("IWW:", round(cor(df$qty.2, df$est_qty.2), 2), ", ", sep = ""))
cat(paste("Rail:", round(cor(df$qty.3, df$est_qty.3), 2), "\n\n"))
