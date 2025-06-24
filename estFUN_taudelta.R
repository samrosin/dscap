estFUN_taudelta <- function(data){
  A = data$A
  R = data$R
  X = matrix(c(1, data$X1, data$X2), nrow = 1, ncol = 3)
  
  function(theta, models){
    py0 <- plogis(X %*% theta[1:3])
    py1 <- ifelse(data$A==1, data$Y, 0)
    py2 <- plogis(X %*% theta[4:6])
    py3 <- plogis(X %*% theta[7:9])
    py4 <- plogis(X %*% theta[10:12])
    py5 <- plogis(X %*% theta[13:15])
    
    ps0 <- X %*% theta[16:18]
    ps1 <- ifelse(A==1, data$S, 0)  
    ps2 <- X %*% theta[19:21]
    ps3 <- X %*% theta[22:24]
    ps4 <- X %*% theta[25:27]
    ps5 <- X %*% theta[28:30]
    
    c(
      if(A==0) grab_psiFUN(models[[1]], data)(theta[1:3]) else rep(0, 3),
      #if(A==1) (data$Y - theta[4]) else 0,
      if(A==2) grab_psiFUN(models[[2]], data)(theta[4:6]) else rep(0, 3),
      if(A==3) grab_psiFUN(models[[3]], data)(theta[7:9]) else rep(0, 3),
      if(A==4) grab_psiFUN(models[[4]], data)(theta[10:12]) else rep(0, 3),
      if(A==5) grab_psiFUN(models[[5]], data)(theta[13:15]) else rep(0, 3),
      
      if(A==0) grab_psiFUN(models[[6]], data)(theta[16:18]) else rep(0, 3),
      #if(A==1) (data$S - theta[20]) else 0,
      if(A==2) grab_psiFUN(models[[7]], data)(theta[19:21]) else rep(0, 3),
      if(A==3) grab_psiFUN(models[[8]], data)(theta[22:24]) else rep(0, 3),
      if(A==4) grab_psiFUN(models[[9]], data)(theta[25:27]) else rep(0, 3),
      if(A==5) grab_psiFUN(models[[10]], data)(theta[28:30]) else rep(0, 3),
      
      (1-R)*(py0 - theta[31]),
      if(A==1) (data$Y - theta[32]) else 0,
      (1-R)*(py2 - theta[33]),
      (1-R)*(py3 - theta[34]),
      (1-R)*(py4 - theta[35]),
      (1-R)*(py5 - theta[36]),
      
      (1-R)*(ps0 - theta[37]),
      if(A==1) (data$S - theta[38]) else 0,
      (1-R)*(ps2 - theta[39]),
      (1-R)*(ps3 - theta[40]),
      (1-R)*(ps4 - theta[41]),
      (1-R)*(ps5 - theta[42]),
      
      (1 - theta[32] / theta[31]) - theta[43], #tau 1
      (1 - theta[33] / theta[31]) - theta[44], #tau 2
      (1 - theta[34] / theta[31]) - theta[45], #tau 3
      (1 - theta[35] / theta[31]) - theta[46], #tau 4
      (1 - theta[36] / theta[31]) - theta[47], #tau 5
      
      (theta[38] - theta[37]) - theta[48], #delta1
      (theta[39] - theta[37]) - theta[49], #delta2
      (theta[40] - theta[37]) - theta[50], #delta3
      (theta[41] - theta[37]) - theta[51], #delta4
      (theta[42] - theta[37]) - theta[52], #delta5
      
      
      # causal association parameter beta -- linear regression slope
      cov(c(theta[43],theta[44],theta[45],theta[46],theta[47]),
          c(theta[48],theta[49],theta[50],theta[51],theta[52])) / 
        var(c(theta[48],theta[49],theta[50],theta[51],theta[52]))  -
        theta[53],
    
    # causal association pearson rho
    cov(c(theta[43],theta[44],theta[45],theta[46],theta[47]),
        c(theta[48],theta[49],theta[50],theta[51],theta[52])) / 
     ( sd(c(theta[43],theta[44],theta[45],theta[46],theta[47])) * 
         sd(c(theta[48],theta[49],theta[50],theta[51],theta[52]))) - 
      theta[54]
      
        
    )
  }
}
