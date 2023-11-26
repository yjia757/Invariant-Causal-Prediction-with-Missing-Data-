# Nov 7th, 2023 
# Author: Yiran Jia 

# Demonstrate the toy example on page 17. 

library(MASS)
library(combinat)

basic <- function(alpha1 = 1, alpha2 = 1, gamma_value = 10, sample_size_n1 = 3000, 
                  sample_size_n2 = 3000, bs_round = 3, true_vec = c(1, 0)) { 
  
  coll_results <- vector("list", length = bs_round)
  
  # Generate bs_round many times data for later taking average result of them
  for (i in 1:bs_round) {
    set.seed(122 + i)
    # Create data from environment 1 
    n1 = sample_size_n1
    d = 4
    mu <- rep(0, d)
    Sigma <- diag(d)
    u1 <- mvrnorm(n1, mu, Sigma)
    
    h1 <- u1[,4]
    x11 <- sqrt(5) * u1[, 2]
    y1 <- 1 * x11 + 1 * h1 + sqrt(0.5) * u1[, 1]  
    x12 <- alpha1 * h1 + u1[, 3]
    
    env1 <- data.frame(x11, x12, y1)
    
    set.seed(455 + i)
    # Create data from environment 2 
    n2 = sample_size_n2
    d = 4
    mu <- rep(0, d)
    Sigma <- diag(d)
    u2 <- mvrnorm(n1, mu, Sigma)
    
    h2 <- u2[,4]
    x21 <- sqrt(0.5) * u2[, 2]
    y2 <- 1 * x21 + 1 * h1 + sqrt(0.5) * u2[, 1]
    x22 <- alpha2 * h2 + u2[, 3]
    
    
    env2 <- data.frame(x21, x22, y2)
    
    p <- 1:2
    power_set <- list()
    # Generate the power set
    for (k in 0:length(p)) {
      # Get all combinations of size k
      combos <- combn(p, k, simplify = FALSE)
      # Append the combinations to the power set
      power_set <- c(power_set, combos)
    }
    
    # Start EILLS algorithm
    y1 <- env1$y1 
    X1 <- as.matrix(env1[, 1:2]) 
    
    y2 <- env2$y2
    X2 <- as.matrix(env2[, 1:2])
    
    gamma = gamma_value
    
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
      
      
      rownames(term1)<-NULL
      rownames(term3)<-NULL
      colnames(term1)<-NULL
      colnames(term3)<-NULL
      
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
    
    loss <- min(coll_obj_value)
    subset_member <- power_set[which.min(coll_obj_value)][[1]]
    estimate <- coll_temp_beta[[which.min(coll_obj_value)]]
    
    if (length(subset_member) == 2) {
      final_estimate <- estimate
    } else {
      temp <- rep(0, 2)
      temp[subset_member] <- estimate
      final_estimate <- temp
    }
      
    rec1 <- cov(x11, 1 * h1 + sqrt(0.5) * u1[, 1]) # notice I update the error 
    rec2 <- cov(x12, 1 * h1 + sqrt(0.5) * u1[, 1]) # notice I update the error 
    rec3 <- c(rec1, rec2)
    rec4 <- cov(x21, 1 * h2 + sqrt(0.5) * u2[, 1]) # notice I update the error 
    rec5 <- cov(x22, 1 * h2 + sqrt(0.5) * u2[, 1]) # notice I update the error 
    rec6 <- c(rec4, rec5)
    rec7 <- cov(x11, x12)
    rec8 <- cov(x21, x22)
    # This is beta_S^e, where S = {2}, and e = environment 1 
    rec9 <- true_vec[2] + solve(mean(x12^2)) * (rec2 + rec7 * 1)
    # This is beta_S^e, where S = {2}, and e = environment 2 
    rec10 <- true_vec[2] + solve(mean(x22^2)) * (rec5 + rec8 * 1)
    # This is beta_S^e, where S = {1,2}, and e = environment 1
    rec11 <- true_vec + solve(cov(X1[,c(1,2)])) %*% c(rec1, rec2)
    row.names(rec11) <- NULL
    # This is beta_S^e, where S = {1,2}, and e = environment 2
    rec12 <- true_vec + solve(cov(X2[,c(1,2)])) %*% c(rec4, rec5)
    row.names(rec12) <- NULL
    
    temp_result <- list(final_estimate, rec1, rec2, rec3, rec4, rec5, rec6, rec7, rec8, rec9, rec10, rec11, rec12)
    coll_results[[i]] <- temp_result
    
  }
  
  # Collect and Calculate the Average Estimation on each variable 
  all_est <- vector("list", length = 2) # Because we have 2 variables 
  for (j in 1:2) {
    for (i in 1:bs_round) {
      all_est[[j]][i] <- coll_results[[i]][[1]][j]
    }
  }
    
  avg_est <- sapply(all_est, mean)
    
    
  # Calculate the l2 error 
  truth <- true_vec
  l2_error <- sqrt(sum((avg_est - truth)^2))
    
  # Collect and Calculate h_S1, where S1 = {2}
  h_S1_temp1 <- sapply(coll_results, function(x) x[[3]]) # rec2 
  h_S1_temp2 <- sapply(coll_results, function(x) x[[6]]) # rec5
  avg_h_S1_temp1 <- mean(h_S1_temp1)
  avg_h_S1_temp2 <- mean(h_S1_temp2)
  h_S1 <- (1/2) * (avg_h_S1_temp1 - mean(avg_h_S1_temp1, avg_h_S1_temp2))^2 + 
    (1/2) * (avg_h_S1_temp2 - mean(avg_h_S1_temp1, avg_h_S1_temp2))^2
  
  # Collect and Calculate h_S2, where S2 = {1,2}
  h_S2_temp1 <- sapply(coll_results, function(x) x[[2]]) # rec1
  h_S2_temp2 <- sapply(coll_results, function(x) x[[5]]) # rec4
  avg_h_S2_temp1 <- mean(h_S2_temp1)
  avg_h_S2_temp2 <- mean(h_S2_temp2)
  h_S2_temp3 <- c(avg_h_S2_temp1, avg_h_S1_temp1)
  h_S2_temp4 <- c(avg_h_S2_temp2, avg_h_S1_temp2)
  h_S2_temp5 <- (h_S2_temp3 + h_S2_temp4) / 2 
  h_S2 <- (1/2) * (norm(h_S2_temp3 - h_S2_temp5,type="2"))^2 + (1/2) * (norm(h_S2_temp4 - h_S2_temp5,type="2"))^2
  
  # Collect and Calculate b_S1, where S1 = {2}
  b_S1 <- ((1/2) * avg_h_S1_temp1 + (1/2) * avg_h_S1_temp2)^2

  # Collect and Calculate b_S2, where S2 = {1,2}
  b_S2 <- norm((1/2) * h_S2_temp3 + (1/2) * h_S2_temp4, type="2")
  
  # Collect and Calculate d_S1, where S1 = {2}
  d_S1_temp1 <- sapply(coll_results, function(x) x[[10]]) # rec9
  d_S1_temp2 <- sapply(coll_results, function(x) x[[11]]) # rec10
  avg_d_S1_temp1 <- mean(d_S1_temp1)
  avg_d_S1_temp2 <- mean(d_S1_temp2)
  d_S1 <- (1/2) * (avg_d_S1_temp1 - mean(avg_d_S1_temp1, avg_d_S1_temp2))^2 + 
    (1/2) * (avg_d_S1_temp2 - mean(avg_d_S1_temp1, avg_d_S1_temp2))^2
  
  # Collect and Calculate d_S2, where S2 = {1,2}
  d_S2_temp1 <- sapply(coll_results, function(x) x[[12]]) # rec9
  d_S2_temp2 <- sapply(coll_results, function(x) x[[13]]) # rec10
  
  sum_d_S2_temp1 <- 0
  for(i in 1:bs_round){
    sum_d_S2_temp1 <- d_S2_temp1[,i] + sum_d_S2_temp1
  }
  avg_d_S2_temp1 <- sum_d_S2_temp1 / bs_round

  sum_d_S2_temp2 <- 0
  for(i in 1:bs_round){
    sum_d_S2_temp2 <- d_S2_temp2[,i] + sum_d_S2_temp2
  }
  avg_d_S2_temp2 <- sum_d_S2_temp2 / bs_round
  
  d_S2_temp3 <- (avg_d_S2_temp1 + avg_d_S2_temp2)/2
    
  d_S2 <- (1/2) * (norm(avg_d_S2_temp1 - d_S2_temp3,type="2"))^2 + (1/2) * (norm(avg_d_S2_temp2 - d_S2_temp3,type="2"))^2
    
  # Compute the quantity b_S1 / h_S1
  
  b_S1_div_h_S1 <- b_S1 / h_S1
  
  # Compute the quantity b_S2 / h_S2
  
  b_S2_div_h_S2 <- b_S2 / h_S2
  
  # Compute the quantity b_S1 / d_S1
  
  b_S1_div_d_S1 <- b_S1 / d_S1
  
  # Compute the quantity b_S2 / d_S2
  
  b_S2_div_d_S2 <- b_S2 / d_S2
  
  # Compute the quantity sup b_S / h_S
  
  sup_b_S_div_h_S <- max(b_S1_div_h_S1, b_S2_div_h_S2)
  
  # Compute the quantity sup b_S / d_S
  sup_b_S_div_d_S <- max(b_S1_div_d_S1, b_S2_div_d_S2)
  
  
  return(list(avg_est, # 1
              l2_error, # 2
              h_S1, # 3
              h_S2, # 4
              b_S1, # 5
              b_S2, # 6
              d_S1, # 7
              d_S2, # 8
              b_S1_div_h_S1, # 9
              b_S2_div_h_S2, # 10
              b_S1_div_d_S1, # 11
              b_S2_div_d_S2, # 12
              sup_b_S_div_h_S, # 13
              sup_b_S_div_d_S # 14
              )) 
}
  
  
define_alpha1_vec <- seq(3, 0, length.out = 30)
define_alpha2_vec <- seq(0, 3, length.out = 30)
length_alpha <- length(define_alpha1_vec)


plot_result <- function(alpha1_vec = define_alpha1_vec, 
                        alpha2_vec = define_alpha2_vec, 
                        gamma_value_vec = seq(50, 50, length.out = length_alpha), 
                        sample_size_n1_vec = seq(3000, 3000, length.out = length_alpha), 
                        sample_size_n2_vec = seq(3000, 3000, length.out = length_alpha), 
                        bs_round_value = 3){
  
  
  result_avg_est_holder <- vector("list", length = length_alpha) # 1 
  result_l2_error_holder <- numeric(length_alpha) # 2 
  result_h_S1_holder <- numeric(length_alpha) # 3 
  result_h_S2_holder <- numeric(length_alpha) # 4
  result_b_S1_holder <- numeric(length_alpha) # 5
  result_b_S2_holder <- numeric(length_alpha) # 6
  result_d_S1_holder <- numeric(length_alpha) # 7
  result_d_S2_holder <- numeric(length_alpha) # 8
  result_b_S1_div_h_S1_holder <- numeric(length_alpha) # 9
  result_b_S2_div_h_S2_holder <- numeric(length_alpha) # 10
  result_b_S1_div_d_S1_holder <- numeric(length_alpha) # 11
  result_b_S2_div_d_S2_holder <- numeric(length_alpha) # 12
  result_sup_b_S_div_h_S_holder <- numeric(length_alpha) # 13
  result_sup_b_S_div_d_S_holder <- numeric(length_alpha) # 14
  
  for (i in 1:length_alpha){
    full_result <- basic(alpha1 = alpha1_vec[i], alpha2 = alpha2_vec[i], 
                         gamma_value = gamma_value_vec[i], sample_size_n1 = sample_size_n1_vec[i], 
                         sample_size_n2 = sample_size_n2_vec[i], bs_round = bs_round_value)
    result_avg_est_holder[i] <- full_result[1] 
    result_l2_error_holder[i] <- full_result[2] 
    result_h_S1_holder[i] <- full_result[3] 
    result_h_S2_holder[i] <- full_result[4] 
    result_b_S1_holder[i] <- full_result[5] 
    result_b_S2_holder[i] <- full_result[6] 
    result_d_S1_holder[i] <- full_result[7] 
    result_d_S2_holder[i] <- full_result[8] 
    result_b_S1_div_h_S1_holder[i] <- full_result[9] 
    result_b_S2_div_h_S2_holder[i] <- full_result[10] 
    result_b_S1_div_d_S1_holder[i] <- full_result[11]
    result_b_S2_div_d_S2_holder[i] <- full_result[12]
    result_sup_b_S_div_h_S_holder[i] <- full_result[13]
    result_sup_b_S_div_d_S_holder[i] <- full_result[14]
  }
  
  # Plot the relationship between the sequence and the estimation - 
  # confirmation of previous presentation with now average of runs 
  plot_1_x_axis <- seq(1:length_alpha)
  plot_1_y_axis <- result_avg_est_holder
  
  plot_1_y_matrix <- do.call(cbind, plot_1_y_axis)
  
  colors <- rainbow(nrow(plot_1_y_matrix))
  lty_seq <- 1:nrow(plot_1_y_matrix)
  
  dev.new()
  plot(plot_1_x_axis, plot_1_y_matrix[1,], type = "n", xlim = range(plot_1_x_axis), ylim = c(-1,5),
       xlab = "seq", ylab = "Avg Estimation", main = "Average Estimation of covariates")
  
  for(i in 1:nrow(plot_1_y_matrix)) {
    lines(plot_1_x_axis, plot_1_y_matrix[i,], col = colors[i], lwd = 2, lty = lty_seq[i])
  }
  
  # Add a legend
  legend("topright", legend = paste("Line", 1:nrow(plot_1_y_matrix)), col = colors, lty = lty_seq, lwd = 2, cex = 0.4)
  
  
  # Plot the relationship between the l2 error and the h_S, S = {2}
  plot_2_x_axis <- as.numeric(result_h_S1_holder)
  plot_2_y_axis <- as.numeric(result_l2_error_holder)
  dev.new()
  plot(plot_2_x_axis, plot_2_y_axis, type = "p", main = "l2 error vs h_S where S = {2}", xlab = "h_S, S={2}", ylab = "l2 error",
       xlim = range(plot_2_x_axis), ylim = c(0,2), cex = 1.5, col = "blue", lwd = 2)

  # Plot the relationship between the ls error and the d_S, S = {2}
  plot_3_x_axis <- as.numeric(result_d_S1_holder)
  plot_3_y_axis <- as.numeric(result_l2_error_holder)
  dev.new()
  plot(plot_3_x_axis, plot_3_y_axis, type = "p", main = "l2 error vs d_S where S = {2}", xlab = "d_S, S={2}", ylab = "l2 error",
       xlim = range(plot_3_x_axis), cex = 1.5, col = "blue", lwd = 2)
  
  # Plot the relationship between sequence and the h_S, S = {2}
  
  plot_4_x_axis <- seq(1:length_alpha)
  plot_4_y_axis <- result_h_S1_holder
  dev.new()
  plot(plot_4_x_axis, plot_4_y_axis, type = "p", main = "result_h_S1_holder", xlab = "seq", ylab = "result_h_S1_holder",
       xlim = range(plot_4_x_axis), cex = 1.5, col = "blue", lwd = 2) 
  
  
  # Plot the relationship between sequence and the d_S, S = {2}
  
  plot_5_x_axis <- seq(1:length_alpha)
  plot_5_y_axis <- result_d_S1_holder
  dev.new()
  plot(plot_5_x_axis, plot_5_y_axis, type = "p", main = "result_d_S1_holder", xlab = "seq", ylab = "result_d_S1_holder",
       xlim = range(plot_5_x_axis), cex = 1.5, col = "blue", lwd = 2) 
  
  # Plot the relationship between sequence and the result_sup_b_S_div_h_S_holder
  plot_6_x_axis <- seq(1:length_alpha)
  plot_6_y_axis <- result_sup_b_S_div_h_S_holder
  dev.new()
  plot(plot_6_x_axis, plot_6_y_axis, type = "p", main = "result_sup_b_S_div_h_S_holder", xlab = "seq", ylab = "result_sup_b_S_div_h_S_holder",
       xlim = range(plot_6_x_axis), cex = 1.5, col = "blue", lwd = 2) 
  
  
  # Plot the relationship between sequence and the result_sup_b_S_div_d_S_holder
  plot_7_x_axis <- seq(1:length_alpha)
  plot_7_y_axis <- result_sup_b_S_div_d_S_holder
  dev.new()
  plot(plot_7_x_axis, plot_7_y_axis, type = "p", main = "result_sup_b_S_div_d_S_holder", xlab = "seq", ylab = "result_sup_b_S_div_d_S_holder",
       xlim = range(plot_6_x_axis), cex = 1.5, col = "blue", lwd = 2) 
  

  # Plot the relationship between l2 error and result_sup_b_S_div_h_S_holder 
  plot_8_x_axis <- as.numeric(result_sup_b_S_div_h_S_holder)
  plot_8_y_axis <- as.numeric(result_l2_error_holder)
  dev.new()
  plot(plot_8_x_axis, plot_8_y_axis, type = "p", main = "l2 error vs b_S / h_S", xlab = "b_S / h_S", ylab = "l2 error",
       xlim = range(plot_8_x_axis), cex = 1.5, col = "blue", lwd = 2)  
  
  
  # Plot the relationship between l2 error and result_sup_b_S_div_d_S_holder 
  plot_9_x_axis <- as.numeric(result_sup_b_S_div_d_S_holder)
  plot_9_y_axis <- as.numeric(result_l2_error_holder)
  dev.new()
  plot(plot_9_x_axis, plot_9_y_axis, type = "p", main = "l2 error vs b_S / d_S", xlab = "b_S / d_S", ylab = "l2 error",
       xlim = range(plot_8_x_axis), cex = 1.5, col = "blue", lwd = 2)  
  
}
    
    