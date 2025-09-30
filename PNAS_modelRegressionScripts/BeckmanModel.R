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
