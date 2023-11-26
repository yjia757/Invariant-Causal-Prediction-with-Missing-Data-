# Step 1 - Create simulation data 
# It is in low dimension context, therefore l0 norm penalty is not involved here
# install.packages("MASS")
library(MASS)

# By the way we generate, the set of important variables is {1, 2, 3}; after 
# some calculations using the definition of linear spurious variables, the set 
# of linear spurious variables is {7, 8, 9}

# Create data from environment 1 
set.seed(123)
n1 = 300
d = 13 
mu <- rep(0, d)
Sigma <- diag(d)
u1 <- mvrnorm(n1, mu, Sigma)
x11 <- u1[, 1]
x14 <- u1[, 4]
x12 <- sin(x14) + u1[, 2]
x13 <- cos(x14) + u1[, 3]
x15 <- sin(x13 + u1[, 5])
x110 <- 2.5 * x11 + 1.5 * x12 + u1[, 10]
y1 <- 3 * x11 + 2 * x12 - 0.5 * x13 + u1[, 13]
x16 <- 0.8 * y1 * u1[, 6]
x17 <- 0.5 * x13 + y1 + u1[, 7]
x18 <- 0.5 * x17 - y1 + x110 + u1[, 8]
x19 <- tanh(x17) + 0.1 * cos(x18) + u1[, 9]
x111 <- 0.4 * (x17 + x18) * u1[, 11]
x112 <- u1[, 12]
env1 <- data.frame(x11, x12, x13, x14, x15, x16, x17, x18, x19, x110, x111, x112, y1)

# Create data from environment 2
set.seed(456)
n2 = 300
d = 13 
mu <- rep(0, d)
Sigma <- diag(d)
u2 <- mvrnorm(n2, mu, Sigma)
x21 <- u2[, 1]
x24 <- (u2[, 4])^2 - 1 
x22 <- sin(x24) + u2[, 2]
x23 <- cos(x24) + u2[, 3]
x25 <- sin(x23 + u2[, 5])
x210 <- 2.5 * x21 + 1.5 * x22 + u2[, 10]
y2 <- 3 * x21 + 2 * x22 - 0.5 * x23 + u2[, 13]
x26 <- 0.8 * y2 * u2[, 6]
x27 <- 4 * x23 + tanh(y2) + u2[, 7]
x28 <- 0.5 * x27 - y2 + x210 + u2[, 8]
x29 <- tanh(x27) + 0.1 * cos(x28) + u2[, 9]
x211 <- 0.4 * (x27 + x28) * u2[, 11]
x212 <- u2[, 12]

env2 <- data.frame(x21, x22, x23, x24, x25, x26, x27, x28, x29, x210, x211, x212, y2)


# Step 2 - do regular least square using traditional tool for later comparison 
model1 <- lm(y1 ~ . - 1, data = env1)
summary(model1)

model2 <- lm(y2 ~ . - 1, data = env2)
summary(model2)

names_env1 <- names(env1)
names(env1) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "y")
names_env2 <- names(env2)
names(env2) <- c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9", "x10", "x11", "x12", "y")
comb_data <- rbind(env1, env2)
model3 <- lm(y ~ ., data = comb_data)
summary(model3)
names(env1) <- names_env1
names(env2) <- names_env2


# Step 3: find the minimizer of the pooled LS plus regularization J using gradient descent. If gamma = 0 , then it is just pooled LS. 

# set bounds for the parameters
lower_bounds <- c(rep(-Inf, 12), rep(0, 2))  # No lower bound on beta, but omega must be >= 0
upper_bounds <- c(rep(Inf, 12), rep(1, 2)) # No upper bound on beta, but omega each entry must be <= 1

pooled_J_obj <- function(params, env1, env2, gamma_pool = 3, gamma = 0){
  
  # Extract variables

  beta <- params[1:12]
  omega <- params[13:14]
  
  y1 <- env1$y1
  X1 <- as.matrix(env1[, 1:12]) 
  
  y2 <- env2$y2
  X2 <- as.matrix(env2[, 1:12]) 
  
  residual1 <- y1 - X1 %*% beta
  residual2 <- y2 - X2 %*% beta
  
  # Calculate weighted sum of squared errors
  error_pooled <- omega[1] * (1/length(y1)) * sum((residual1)^2) + 
    omega[2] * (1/length(y2)) * sum((residual2)^2)
  
  error_J <- 0
  
  for (j in 1:length(beta)){
    if (beta[j] != 0){
      error_J = error_J + omega[1] * (1/length(y1)^2) * (X1[, j] %*% (y1 - (X1 %*% beta)))^2 + 
        omega[2] * (1/length(y2)^2) * (X2[, j]  %*% (y2 - (X2 %*% beta)))^2
    }
  }
  
  penalty <- abs(sum(omega) - 1)
  
  return(error_pooled + gamma_pool * penalty + gamma * error_J)
}
# Initial guess for the parameter 
params_initial <- c(model3$coefficients[2:13], 0.5, 0.5)
# Call optim to find the minimizer
# Notice gamma_pool cannot be 0 because o.w. it will automatically converge to they be 0 by the construction of the objective
result <- optim(params_initial, pooled_J_obj, method = "L-BFGS-B", lower = lower_bounds, 
                upper = upper_bounds, env1 = env1, env2 = env2, gamma_pool = 3, gamma = 0)
print(result)


# Step 3: find the minimizer of the FIXED pooled LS plus regularization J using gradient descent. If gamma = 0 , then it is just pooled LS.

pooled_J_obj <- function(params, env1, env2, omega, gamma){
  
  beta <- params
  y1 <- env1$y1
  X1 <- as.matrix(env1[, 1:12]) 
  
  y2 <- env2$y2
  X2 <- as.matrix(env2[, 1:12])
  
  residual1 <- y1 - X1 %*% beta
  residual2 <- y2 - X2 %*% beta
  
  # Calculate weighted sum of squared errors
  error_pooled <- omega[1] * (1/length(y1)) * sum((residual1)^2) + 
    omega[2] * (1/length(y2)) * sum((residual2)^2)
  
  error_J <- 0
  
  for (j in 1:length(beta)){
    if (beta[j] != 0){
      error_J = error_J + omega[1] * (1/length(y1)^2) * (X1[, j] %*% (y1 - (X1 %*% beta)))^2 + 
        omega[2] * (1/length(y2)^2) * (X2[, j]  %*% (y2 - (X2 %*% beta)))^2
    }
  }
  
  return(error_pooled + gamma * error_J)
}

# Initial guess for the parameter 
params_initial <- model3$coefficients[2:13]
# Call optim to find the minimizer
# Notice gamma_pool cannot be 0 because o.w. it will automatically converge to they be 0 by the construction of the objective
result <- optim(params_initial, pooled_J_obj, method = "BFGS", env1 = env1, env2 = env2, omega = c(0.5, 0.5), gamma = 0)
print(result)


# Step 4: doing step 4 but using the authors' brute force search method 

# install.packages("combinat")
library(combinat)
p <- 1:12
power_set <- list()

# Generate the power set
for (k in 0:length(p)) {
  # Get all combinations of size k
  combos <- combn(p, k, simplify = FALSE)
  # Append the combinations to the power set
  power_set <- c(power_set, combos)
}

y1 <- env1$y1 
X1 <- as.matrix(env1[, 1:12]) 

y2 <- env2$y2
X2 <- as.matrix(env2[, 1:12])

gamma = 500

coll_obj_value  <- c()
coll_temp_beta <- list()

for (po in 1:length(power_set)){
  
  S <- power_set[po][[1]]
  
  # First find the argmin of the objective function 
  
  term2 <- 0
  
  for (j in 1:length(S)){
    term2 <- term2 + gamma * (1/n1)^2 * t(X1[, S]) %*% X1[, S[j]]  %*% t(X1[, S[j]]) %*% X1[, S] + 
      gamma * (1/n2)^2 * t(X2[, S]) %*% X2[, S[j]]  %*% t(X2[, S[j]]) %*% X2[, S]
  }
    
  term1 <- solve((1/n1) * t(X1[, S]) %*% X1[, S] + (1/n2) * t(X2[, S]) %*% X2[, S] + term2)
  
  term4 <- 0 
  
  for (j in 1:length(S)){
    term4 <- suppressWarnings(term4 + gamma * (1/n1)^2 * (t(X1[, S[j]]) %*% y1) * as.vector((t(X1[, S]) %*% X1[, S[j]])) + 
      gamma * (1/n2)^2 * (t(X2[, S[j]]) %*% y2) * as.vector((t(X2[, S]) %*% X2[, S[j]])))
  }
  
  term3 <- (1/n1) * t(X1[, S]) %*% y1 + (1/n2) * t(X2[, S]) %*% y2 + term4 
  
  temp_beta <- term1 %*% term3
  coll_temp_beta[[po]] <- temp_beta
  
 # Second find the corresponding objective value with beta 
  
  term5 <- 0
  for (j in 1:length(S)){ 
    term5 <- term5 + gamma * (1/2) * (1/n1)^2 * (t(matrix(X1[, S[j]])) %*% (y1 - (X1[, S] %*% temp_beta)))^2 + 
      gamma * (1/2) * (1/n2)^2 * (t(matrix(X2[, S[j]])) %*% (y2 - (X2[, S] %*% temp_beta)))^2
  }
  
  value <- (1/2) * (1/n1) * sum((y1 - X1[, S] %*% temp_beta)^2) + (1/2) * (1/n2) * sum((y2 - X2[, S] %*% temp_beta)^2) + term5
  
  coll_obj_value[po] <- value
}

# Find which subset of [p] gives the minimum objective value, and what is the corresponding beta value 
loss <- min(coll_obj_value)
emp_S <- power_set[which.min(coll_obj_value)]
emp_temp_beta <- coll_temp_beta[[which.min(coll_obj_value)]]
emp_beta <- rep(0, 12)
indices <- as.numeric(sub("x1", "", rownames(emp_temp_beta)))
emp_beta[indices] <- emp_temp_beta[, 1]
emp_beta <- as.matrix(emp_beta, nrow=12, ncol=1)
print(list(loss, emp_S, emp_beta))

# If you are interested in knowing the minimized objective value of estimator beta corresponds to 
# a partciular set S in the power set, you can first use below to find the index of that S in the 
# power set, and then find the minimized objective value by coll_obj_value[po]
target_set <- c(1, 2, 3, 5)
for (po in 1:length(power_set)) {
  if (length(power_set[[po]]) == length(target_set)) {
    if (all(sort(power_set[[po]]) == sort(target_set))){
      print(po)
    }
  }
}

# Step 5: implementing multipler bootstrap for the estimation of beta 

num_set_bootstrap = 100
set_all_b_ci_coll <- list()
for (gbs in 1:num_set_bootstrap){
  
  num_bootstrap = 100
  all_b_beta <- list()
  all_b_obj <- c()
  all_b <- list(all_b_beta, all_b_obj)
  
  for (bs in 1:num_bootstrap){
    
    # Focus on one bootstrap 
    b_coll_obj_value  <- c()
    b_coll_temp_beta <- list()
    
    for (po in 1:length(power_set)){
      
      S <- power_set[po][[1]]
      
      # Generate Rademacher random variables 
      gene_radem <- function(n) {
        sample(c(-1, 1), size = n, replace = TRUE, prob = c(1/2, 1/2))
      }
      
      dw1 <- list()
      for (i in 1:(length(S) + 1)){
        dw1[[i]] <- diag(gene_radem(n1) + 1 )
      }
      
      dw2 <- list()
      for(i in 1:(length(S) + 1)){
        dw2[[i]] <- diag(gene_radem(n2) + 1)
      }
      
      # First find the argmin of the objective function 
      term2 <- 0
      
      for (j in 1:length(S)){
        term2 <- term2 + gamma * (1/n1)^2 * t(X1[, S]) %*% dw1[[j+1]] %*% X1[, S[j]]  %*% t(X1[, S[j]]) %*% X1[, S] + 
          gamma * (1/n2)^2 * t(X2[, S]) %*% dw2[[j+1]] %*% X2[, S[j]]  %*% t(X2[, S[j]]) %*% X2[, S]
      }
      
      term1 <- solve((1/n1) * t(X1[, S]) %*% dw1[[1]] %*% X1[, S] + (1/n2) * t(X2[, S]) %*% dw2[[1]] %*% X2[, S] + term2)
      
      term4 <- 0 
      
      for (j in 1:length(S)){
        
        term4 <- suppressWarnings(term4 + gamma * (1/n1)^2 * (t(dw1[[j+1]] %*% X1[, S[j]]) %*% y1) * as.vector((t(X1[, S]) %*%  dw1[[j+1]] %*% X1[, S[j]])) + 
                                    gamma * (1/n2)^2 * (t(dw2[[j+1]] %*% X2[, S[j]]) %*% y2) * as.vector((t(X2[, S]) %*% dw2[[j+1]] %*% X2[, S[j]])))
      }
      
      term3 <- (1/n1) * t(X1[, S]) %*% dw1[[1]] %*% y1 + (1/n2) * t(X2[, S]) %*% dw2[[1]] %*% y2 + term4 
      
      temp_beta <- term1 %*% term3
      b_coll_temp_beta[[po]] <- temp_beta
      
      # Second find the corresponding objective value with beta
      
      term5 <- 0
      for (j in 1:length(S)){ 
        term5 <- term5 + gamma * (1/2) * (1/n1)^2 * (t(dw1[[j+1]] %*% matrix(X1[, S[j]])) %*% (y1 - (X1[, S] %*% temp_beta)))^2 + 
          gamma * (1/2) * (1/n2)^2 * (t(dw2[[j+1]] %*% matrix(X2[, S[j]])) %*% (y2 - (X2[, S] %*% temp_beta)))^2
      }
      
      value <- (1/2) * (1/n1) * sum(dw1[[1]] %*% (y1 - X1[, S] %*% temp_beta)^2) + (1/2) * (1/n2) * sum(dw2[[1]] %*% (y2 - X2[, S] %*% temp_beta)^2) + term5
      
      b_coll_obj_value[po] <- value
    } 
    
    # Find which subset of [p] gives the minimum objective value, and what is the corresponding beta value
    b_loss <- min(b_coll_obj_value)
    b_emp_S <- power_set[which.min(b_coll_obj_value)]
    b_emp_temp_beta <- b_coll_temp_beta[[which.min(b_coll_obj_value)]]
    b_emp_beta <- rep(0, 12) 
    b_indices <- as.numeric(sub("x1", "", rownames(b_emp_temp_beta)))
    b_emp_beta[b_indices] <- b_emp_temp_beta[, 1]
    b_emp_beta <- as.matrix(b_emp_beta, nrow=12, ncol=1)
    # print(list(b_loss, b_emp_S, b_emp_beta))
    
    all_b[[1]][[bs]] <- b_emp_beta
    all_b[[2]][bs] <- b_loss
  }
  
  # Step 6: compute the pivotal interval for each coefficients 
  
  coll <- list()
  for (j in 1:length(all_b[[1]][[1]])){
    temp_coll <- c()
    for (i in 1:length(all_b[[1]])){
      temp_coll[i] <- all_b[[1]][[i]][[j]]
    }
    coll[[j]] <- temp_coll
  }
  
  quan_coll <- list()
  for (z in 1:length(coll)){
    q_low <- quantile(coll[[z]], probs = 0.025)
    q_up <- quantile(coll[[z]], probs = 0.975)
    temp_quan <- c(q_low, q_up)
    quan_coll[[z]] <- temp_quan
  }
  
  ci_coll <- list()
  for (j in 1:length(quan_coll)){
    ci_low <- 2 * emp_beta[j] - quan_coll[[j]][1]
    ci_up <- 2 * emp_beta[j] - quan_coll[[j]][2]
    temp_ci <- unname(c(ci_low, ci_up))
    ci_coll[[j]] <- temp_ci
  }
  
  set_all_b_ci_coll[[gbs]] <- ci_coll
  
}


# Step 7: compute the coverage probabilities over all the coefficients 

true_beta <- c(3, 2, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0)
cov_coll <- c()

for (c in 1:length(true_beta)){
  temp <- c()
  for(q in 1:length(set_all_b_ci_coll)){
    if (true_beta[c] >= set_all_b_ci_coll[[q]][[c]][1] && true_beta[c] <= set_all_b_ci_coll[[q]][[c]][2]){
      temp[q] <- 1
    } else {
      temp[q] <- 0
    }
  }
  cov_coll[c] <- mean(temp)
}

# To compute the average coverage probability over all the coefficients 
avg_cov_prob <- mean(cov_coll)


# Step 8: compute the average confidence interval widths over all the coefficients 

width_coll <- c() 

for (c in 1:length(true_beta)){
  temp <- c()
  for(q in 1:length(set_all_b_ci_coll)){
    temp[q] <- set_all_b_ci_coll[[q]][[c]][2] - set_all_b_ci_coll[[q]][[c]][1]
  }
  width_coll[c] <- mean(temp)
}

# To compute the average coverage probability over all the coefficients 
avg_width_coll <- mean(width_coll)













