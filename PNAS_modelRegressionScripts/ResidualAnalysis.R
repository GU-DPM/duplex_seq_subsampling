library("readr")
testS_Vals= read_csv("/Users/loeblabm11/Desktop/CRC/TakeToRetreat/TestPythonRegression/FromHome/testS_Vals.csv")
testS_Vals2 = read_csv("/Users/loeblabm11/Desktop/CRC/TakeToRetreat/TestPythonRegression/FromHome/testS_Vals2.csv")
testS_Vals_1 = read_csv("/Users/loeblabm11/Desktop/CRC/TakeToRetreat/TestPythonRegression/FromHome/testS_Vals_20180405_1545.csv")
testVals = seq(1,length(testS_Vals_1$sVal),1)
sVals = c(0.2, 0.19, 0.18, 0.17, 0.16, 0.15, 0.14, 0.13, 0.12, 0.11, 0.1, 0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.005, 0.001,.0008,0.0006,0.0004,0.0002,0.0001,0.00001)
ave1 = c()
ave2 = c()
plusCI1 = c()
minusCI1 = c()
plusCI2 = c()
minusCI2 = c()
for (testVal in testVals) {
  selectModel = beckmanModel2(s=testS_Vals_1$sVal[testVal],var_T=428, D=data3070030$X5, a=c(.015,.97,.015), numPoints = 4, K_mut_eff = testS_Vals_1$k[testVal])
  neutralModel = beckmanModel2(s=testS_Vals_1$sVal[testVal],var_T=428, D=data3070030$X5, a=c(0,1,0), numPoints = 4, K_mut_eff = 9.653627e-07)
  selectResids = (log(selectModel$F_app_unmut) - data3070030$lnF)/data3070030$lnF
  neutralResids = (log(neutralModel$F_app_unmut) - data3070030$lnF)/data3070030$lnF
  #lmSelect = lm(log(selectModel$F_app_unmut) ~ selectModel$D)
  #lmNeut = lm(log(neutralModel$F_app_unmut) ~ neutralModel$D)
  averageResidSelect = sum(abs(selectResids))/length(selectResids)
  averageResidNeut = sum(abs(neutralResids))/length(neutralResids)
  selectPlusCI  = averageResidSelect + 1.96 * sd(abs(selectResids))/sqrt(length(selectResids))
  selectMinusCI = averageResidSelect - 1.96 * sd(abs(selectResids))/sqrt(length(selectResids))
  neutPlusCI = averageResidNeut + 1.96 * sd(abs(neutralResids))/sqrt(length(neutralResids))
  neutMinusCI = averageResidNeut - 1.96 * sd(abs(neutralResids))/sqrt(length(neutralResids))
  ave1 = c(ave1,averageResidSelect)
  ave2 = c(ave2, averageResidNeut)
  plusCI1 = c(plusCI1, selectPlusCI)
  minusCI1 = c(minusCI1, selectMinusCI)
  plusCI2 = c(plusCI2, neutPlusCI)
  minusCI2 = c(minusCI2, neutMinusCI)
}

testPlot = ggplot() + 
  geom_line(mapping = aes(x=testS_Vals_1$sVal, y=ave1), color="blue") +
  geom_line(mapping = aes(x=testS_Vals_1$sVal, y=plusCI1), color="blue", linetype = 2) +
  geom_line(mapping = aes(x=testS_Vals_1$sVal, y=minusCI1), color="blue", linetype = 2) + 
  geom_line(mapping = aes(x=testS_Vals_1$sVal, y=ave2), color="red") +
  geom_line(mapping = aes(x=testS_Vals_1$sVal, y=plusCI2), color="red", linetype = 2) +
  geom_line(mapping = aes(x=testS_Vals_1$sVal, y=minusCI2), color="red", linetype = 2) + 
  scale_y_log10() + xlab("Strength of selection (s)") +
  ylab("Average absolute value of residual as proportion of data") + 
  theme_bw() + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(), 
        axis.line.x = element_line(color="black"), 
        axis.line.y = element_line(color="black"),
        legend.position = "none", 
        axis.title.x = element_text(size=7), 
        axis.text.x = element_text(size=6), 
        axis.title.y = element_text(size=7), 
        axis.text.y = element_text(size=6)
        )
testPlot
ggsave(file="SensitivityThreshold.svg", plot=testPlot, width=183, height=90, units='mm')
