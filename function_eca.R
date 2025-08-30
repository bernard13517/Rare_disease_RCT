lm_model <- lm(Y ~ X1 + X2,data=histo.data)

all_coefs <- coef(lm_model)
  
  # Extract the full variance-covariance matrix of the coefficients
vcov_matrix <- vcov(lm_model)
target_predictors <- c("X1", "X2")
mean_vector <- all_coefs[target_predictors]
sd_estime <- summary(lm_model)$sigma  
  # Subset the variance-covariance matrix for just our target predictors
  cov_matrix <- vcov_matrix[target_predictors, target_predictors]
  sampled_vector <- MASS::mvrnorm(n = 1, mu = mean_vector, Sigma = cov_matrix) 
sim.ECA <- function(npat,delta,mx,r)
{
set.seed(r)
 N <- npat  # Number of patients
rho <- 0.5  # Correlation between X1 and X2

# Mean vector
mean <- mx

# Covariance matrix
cov_matrix <- matrix(c(1, 5 * rho, 5 * rho, 25), nrow = 2)

# Generate the covariates X1 and X2

covariates <- MASS::mvrnorm(N, mean, cov_matrix)
X1 <- covariates[, 1]
X2 <- covariates[, 2]

# Generate the parameters for the outcome model
alpha <- 1  # Intercept
beta1 <- - 0.5  # Coefficient for X1
beta2 <- 0.3  # Coefficient for X2

# Generate the error term
epsilon <- rnorm(N, 0, 1)

# Calculate the outcome Y
Y <- alpha + beta1 * (X1 - 5) + beta2 * (X2 - 30) + epsilon 
A=rep(0,N)
control <-  data.frame(X1, X2, Y,A)
control.actual <- control[1:(N/2),]

# Generate the covariates X1 and X2

covariates <- MASS::mvrnorm(N, mean, cov_matrix)
X1 <- covariates[, 1]
X2 <- covariates[, 2]

# Generate the parameters for the outcome model
alpha <- 1  # Intercept
beta1 <- - 0.5  # Coefficient for X1
beta2 <- 0.3  # Coefficient for X2

# Generate the error term
epsilon <- rnorm(N, 0, 1)

# Calculate the outcome Y
Y <- delta + alpha + beta1 * (X1 - 5) + beta2 * (X2 - 30) + epsilon
A=rep(1,N)
drug <-  data.frame(X1, X2, Y,A)

# Complement last half of placebo 

ECA <- drug[(1+N/2):N,]
# histor <- histo.data |> mutate(A=0)
# histor <- rbind(histor,drug.to_match)

# ps_model <- glm(A ~ X1 + X2, data = histor, family = binomial)
# histor$pscore <- predict(ps_model, type = "response")

# #### Step 2: Perform PS matching using nearest neighbor
# match_result <- matchit(A ~ X1 + X2, data = histor, method = "nearest")
for (k in 1:(N/2))
  {sampled_vector <- MASS::mvrnorm(n = 1, mu = mean_vector, Sigma = cov_matrix)
   ECA$Y[k] = sampled_vector[1]*ECA$X1[k] + sampled_vector[2]*ECA$X2[k] + sd_estime*rnorm(1) 
  }
  
#### Step 3: Extract matched data
# matched_data <- match.data(match_result)
control <- rbind(control.actual,ECA)

data <- rbind(drug,control)

model <- lm(Y ~ A,data=data)
model_summary <- summary(model)

# Extract coefficients
coefficients <- model_summary$coefficients

coef_values <- coefficients[, "Estimate"]

# P-values only
p_values <- coefficients[, "Pr(>|t|)"]

return(c(coef_values[2],p_values[2]))

}

sim.matched(npat=62,delta=0.8,mx=c(5,30),r=7741)


