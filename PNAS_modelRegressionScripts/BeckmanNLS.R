testStart = list(c=0, k=1e-5)

varA1 = .1
varS = 0.08
varB = 2
varT = 25*365
f_r1_T = (exp(varS*varB*varT) - 1)/(varS*varB*varT)
f_r2_T = 1
f_r3_T = (exp(-1*varS*varB*varT) - 1)/(-1*varS*varB*varT)

test = nls(
  F ~ exp(c)*A1*exp(-1*k*D*f_r1_T)+exp(c)*(1 - (2*A1))*exp(-1*k*D)+ exp(c)*A1*exp(-1*k*D*f_r3_T), 
  data=myData, 
  start=testStart, 
  control = nls.control(
    maxiter = 10000, tol=1e-5, minFactor=1e-20, warnOnly = TRUE
    )
  )
summary(test)
