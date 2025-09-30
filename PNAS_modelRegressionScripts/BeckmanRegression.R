beckmanModel <- function(s, 
                         var_T, 
                         K_mut_eff=5.988e-7, 
                         Sites=10338, 
                         D_max=10000, 
                         numPoints=100, 
                         b=0.25, 
                         r_n=1, 
                         a=c(0.1, 0.8, 0.1), 
                         d_n = 0.18) {
  # Test a
  if (sum(a) != 1 | length(a) != 3) {
    stop("a must be a three member list of floats, with the values summing to 1.")
  }
  if (var_T <= 0) {
    stop("var_T must be positive")
  }
  if (var_T != round(var_T)) {
    stop("T must be an integer")
  }
  
  # set up derived constants
  r_1 = r_n + s
  r_2 = r_n
  r_3 = r_n - s
  
  f_r1_T = (exp(log(2) * (r_1 - 1) * b * var_T) - 1)/(log(2) * (r_1 - 1) * b * var_T)
  f_r2_T = 1
  f_r3_T = (exp(log(2) * (r_3 - 1) * b * var_T) - 1)/(log(2) * (r_3 - 1) * b * var_T)
  
  beta_1 = K_mut_eff * f_r1_T
  beta_2 = K_mut_eff * f_r2_T
  beta_3 = K_mut_eff * f_r3_T
  
  # set up depth list
  Depth=seq.int(0, D_max, D_max/numPoints)
  F_1_i1 = Depth * -1*beta_1
  F_2_i1 = Depth * -1*beta_2
  F_3_i1 = Depth * -1*beta_3
  
  F_1 = a[1] * exp(F_1_i1)
  F_2 = a[2] * exp(F_2_i1)
  F_3 = a[3] * exp(F_3_i1)
  
  F_Tot = F_1 + F_2 + F_3
  numMuts = Sites - (F_Tot * Sites)
  output = data.frame(D=Depth, F_app_unmut = F_Tot, Muts=numMuts)
}

beckmanModel2 <- function(s, 
                         var_T, D, 
                         K_mut_eff=5.988e-7, 
                         Sites=10338, 
                         numPoints=100, 
                         b=0.25, 
                         r_n=1, 
                         a=c(0.1, 0.8, 0.1), 
                         d_n = 0.18) {
  # Test a
  if (sum(a) != 1 | length(a) != 3) {
    stop("a must be a three member list of floats, with the values summing to 1.")
  }
  if (var_T <= 0) {
    stop("var_T must be positive")
  }
  if (var_T != round(var_T)) {
    stop("T must be an integer")
  }
  
  # set up derived constants
  r_1 = r_n + s
  r_2 = r_n
  r_3 = r_n - s
  
  f_r1_T = (exp(log(2) * (r_1 - 1) * b * var_T) - 1)/(log(2) * (r_1 - 1) * b * var_T)
  f_r2_T = 1
  f_r3_T = (exp(log(2) * (r_3 - 1) * b * var_T) - 1)/(log(2) * (r_3 - 1) * b * var_T)
  
  beta_1 = K_mut_eff * f_r1_T
  beta_2 = K_mut_eff * f_r2_T
  beta_3 = K_mut_eff * f_r3_T
  
  # set up depth list
  Depth=D
  F_1_i1 = Depth * -1*beta_1
  F_2_i1 = Depth * -1*beta_2
  F_3_i1 = Depth * -1*beta_3
  
  F_1 = a[1] * exp(F_1_i1)
  F_2 = a[2] * exp(F_2_i1)
  F_3 = a[3] * exp(F_3_i1)
  
  F_Tot = F_1 + F_2 + F_3
  numMuts = Sites - (F_Tot * Sites)
  output = data.frame(D=Depth, F_app_unmut = F_Tot, Muts=numMuts)
}

check_fit_model = function(inData, 
                           inS, 
                           inVar_T, 
                           inK_mut_eff=5.988e-7, 
                           inSites=10338, 
                           inB=0.25, 
                           inR_n=1, 
                           inA=c(0.1, 0.8, 0.1), 
                           inD_n = 0.18) {
  inD = inData$X5
  outModel=beckmanModel2(s=inS, 
                      var_T = inVar_T, 
                      D = inD, 
                      K_mut_eff = inK_mut_eff, 
                      Sites = inSites, 
                      b = inB, 
                      r_n = inR_n, 
                      a = inA, 
                      d_n = inD_n)
  resids = inData$lnF - log(outModel$F_app_unmut)
  outRSS = sum(resids^2)
  SSTO = sum((inData$lnF - mean(inData$lnF))^2)
  R_squared = 1 - (outRSS / SSTO)
  output = list(RSS=outRSS, model=outModel, residuals=resids, R_sq = R_squared)
}

beckman_regression = function(sMin, sMax, sStep, 
                              k_exp_min, k_exp_max, 
                              inputData, 
                              inputVar_T, 
                              inputSites=10338, 
                              inputB=0.25, 
                              inputR_n=1, 
                              inputA=c(0.1, 0.8, 0.1), 
                              inputD_n = 0.18) {
  if (sMax <= sMin) {
    stop("sMax must be greater than sMin")
  } 
  testS_vals = c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  s_values = c()
  k_values = c()
  RSS_vals = c()
  rSq = c()
  for (sVal in seq.int(from = sMin, to = sMax, by = sStep)) {
  #for (sVal in testS_vals) {
    for (k_exp in seq.int(from = k_exp_min, to=k_exp_max, by=1)) {
      for (k_main in seq(from = 1,to = 9.99, by = 0.01)) {
        k_mut = k_main * 10 ^ k_exp
        test = check_fit_model(inData = inputData, 
                               inS = sVal, 
                               inVar_T = inputVar_T, 
                               inK_mut_eff = k_mut, 
                               inSites = inputSites, 
                               inB = inputB, 
                               inR_n = inputR_n, 
                               inA = inputA, 
                               inD_n = inputD_n)
        s_values = c(s_values, sVal)
        k_values = c(k_values, k_mut)
        RSS_vals = c(RSS_vals, test$RSS)
        rSq = c(rSq, test$R_sq)
      }
    }
  }
  output = data.frame(s=s_values, k_mut_eff = k_values, RSS = RSS_vals, rSquared = rSq)
}

beckman_regression2 = function(sMin, sMax, sStep, 
                              k_exp_min, k_exp_max, 
                              a1_min, a1_max, a1_step, 
                              inputData, 
                              inputVar_T, 
                              inputSites=10338, 
                              inputB=0.25, 
                              inputR_n=1, 
                              inputD_n = 0.18) {
  if (sMax <= sMin) {
    stop("sMax must be greater than sMin")
  } 
  
  a1_values = c()
  s_values = c()
  k_values = c()
  RSS_vals = c()
  rSq = c()
  for (a1 in seq(from=a1_min, to=a1_max, by=a1_step)) {
    inputA = c(a1, 1-(2*a1), a1)
    for (sVal in seq.int(from = sMin, to = sMax, by = sStep)) {
      print(c(a1,sVal))
      #for (sVal in testS_vals) {
      for (k_exp in seq.int(from = k_exp_min, to=k_exp_max, by=1)) {
        for (k_main in seq(from = 1,to = 9.99, by = 0.01)) {
          k_mut = k_main * 10 ^ k_exp
          test = check_fit_model(inData = inputData, 
                                 inS = sVal, 
                                 inVar_T = inputVar_T, 
                                 inK_mut_eff = k_mut, 
                                 inSites = inputSites, 
                                 inB = inputB, 
                                 inR_n = inputR_n, 
                                 inA = inputA, 
                                 inD_n = inputD_n)
          a1_values = c(a1_values, a1)
          s_values = c(s_values, sVal)
          k_values = c(k_values, k_mut)
          RSS_vals = c(RSS_vals, test$RSS)
          rSq = c(rSq, test$R_sq)
        }
      }
    }
  }
  output = data.frame(a1 = a1_values, s=s_values, k_mut_eff = k_values, RSS = RSS_vals, rSquared = rSq)
}


getMinimums = function(sMin,sMax,sStep,
                       a1_min, a1_max, a1_step,
                       inRegression) {
  sVals = c()
  k_vals=c()
  RSS_vals = c()
  R_squar = c()
  a1_vals = c()
  #for (sVal in c(0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)) {
  for (aVal in seq(from = a1_min, to = a1_max, by = a1_step)) {
    for (sVal in seq.int(from = sMin, to = sMax, by = sStep)) {
      test2 = inRegression[inRegression$a1 == aVal,]
      test3 = test2[test2$s == sVal,]
      test4 = test3[test3$RSS == min(test3$RSS),]
      sVals = c(sVals, test4$s)
      k_vals = c(k_vals, test4$k_mut_eff)
      RSS_vals = c(RSS_vals, test4$RSS)
      R_squar = c(R_squar, test4$rSquared)
      a1_vals = c(a1_vals, test4$a1)
    }
  }
  output = data.frame(a1 = a1_vals, s=sVals, k_mut_eff = k_vals, RSS = RSS_vals, r_squared = R_squar)
}
