library(ggplot2)
myX1 = c(seq(1,9,1), seq(10,90,10), 
        seq(100, 900, 100), 
        seq(1e3, 1e4-1e3, 1e3), 
        seq(1e4, 1e5-1e3, 1e4), 
        seq(1e5, 1e6-1e3, 1e5), 
        seq(1e6, 1e7-1e3, 1e6), 
        seq(1e7, 1e8-1e3, 1e7), 
        seq(1e8, 1e9-1e3, 1e8), 
        seq(1e9, 1e10, 1e9))
yBreaks=c(seq(1e-10,9e-10,1e-10), 
          seq(1e-9,9e-9,1e-9),  
          seq(1e-8,9e-8,1e-8),
          seq(1e-7,9e-7,1e-7), 
          seq(1e-6,9e-6,1e-6),
          seq(1e-5,9e-5,1e-5),
          seq(1e-4,9e-4,1e-4),
          seq(1e-3,9e-3,1e-3),
          seq(1e-2,9e-2,1e-2),
          seq(1e-1,1e-0,1e-1))
majYBreaks=c("1e-10","1e-9","1e-8","1e-7","1e-6","1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1e0")
majorBreaks=c("1e0","1e1","1e2","1e3","1e4","1e5", "1e6", "1e7", "1e8", "1e9", "1e10")
outBreaks = c()
for (breakPoint in majorBreaks) {
  if (breakPoint != "1e10") {
    outBreaks = c(outBreaks, breakPoint, rep("", 8))
  }
  else {outBreaks = c(outBreaks, breakPoint)}
}
outYBreaks = c()
for (breakPoint in majYBreaks) {
  if (breakPoint != "1e0") {
    outYBreaks = c(outYBreaks, breakPoint, rep("", 8))
  }
  else {outYBreaks = c(outYBreaks, breakPoint)}
}
outBreaks
outYBreaks

# Bozic
deaths=0.18
births=0.25
var_delta = deaths/births
bozicY = (1-var_delta^(myX1))/((myX1))
plot(x=myX1, y=bozicY)

myYSottoriva = 1/myX1
myYBeckman = 6e-7/(1-exp(-6e-7*myX1))

testError = abs(myYBeckman - myYSottoriva)/myYSottoriva
testEqualityM <- myYSottoriva == myYBeckman
myErrorBozic = c(0*myX1[myX1 <= 1000], (abs(myYBeckman[myX1 > 1000] - bozicY[myX1 > 1000]))/myYBeckman[myX1 > 1000])
myErrorSottoriva = c((abs(myYSottoriva[myX1 <= 1000] - bozicY[myX1 <= 1000]))/bozicY[myX1 <= 1000], (abs(myYBeckman[myX1 > 1000] - myYSottoriva[myX1 > 1000]))/myYBeckman[myX1 > 1000])
myErrorBeckman = c((abs(bozicY[myX1 <= 1000] - myYBeckman[myX1 <= 1000]))/bozicY[myX1 <= 1000], 0*myYBeckman[myX1 > 1000])

myErrorBozic = myErrorBozic * 100
myErrorSottoriva = myErrorSottoriva * 100
myErrorBeckman = myErrorBeckman * 100
testPlot1=ggplot() + 
  geom_line(mapping=aes(x=myX1, y=bozicY, color="blue"), size=1.5) +
  geom_line(mapping=aes(x=myX1, y=myYBeckman, color="goldenrod"), size=1.5) + 
  geom_line(mapping=aes(x=myX1, y=myYSottoriva, color="red"), linetype=2, size=1.5) +
  geom_ribbon() + 
  scale_y_log10(breaks=yBreaks, labels=outYBreaks) +
  scale_x_log10(breaks=myX1, labels=outBreaks) + 
  xlab("N(t)") + 
  ylab("Average Mutant Allele Frequency of Singlet") +
  theme_bw() + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(), 
        axis.line = element_line(color="black"), 
        legend.position = c(0.75,0.9),
        legend.background = element_blank(),
        axis.text.y = element_text(size=6),
        axis.title.y = element_text(size=7), 
        axis.title.x = element_text(size=7), 
        legend.text = element_text(size=6), 
        axis.text.x = element_text(size=6), 
        ) + 
  scale_color_manual(guide="legend", name="", 
                     values=c("blue"="blue", "goldenrod"="goldenrod", "red"="red"), 
                     labels=c("Stochastic model, infinite sites", 
                              "Deterministic model, non-infinite sites", 
                              "Deterministic model, infinite sites"))
  
testPlot1
ggsave(file="AvgMAF_Sing.eps", plot=testPlot1, width=87, height=87, units='mm')

BozicRange = c(1e0,1e5)
Sottoriva1 = c(1e0,1e2)
Sottoriva2 = c(1e2,1e5)
BeckmanRange = c()
errorPlot = ggplot() + 
  geom_line(mapping=aes(myX1, myErrorBozic, color="blue"), size=1.5) + 
  geom_line(mapping=aes(myX1, myErrorBeckman, color="goldenrod"), size=1.5) + 
  geom_line(mapping=aes(myX1, myErrorSottoriva, color="red"), linetype=2, size=1.5) + 
  ylim(0,101) + 
  geom_hline(yintercept=10) + 
  geom_text(aes(1e8,10, label="10%", vjust=-1, hjust = -2), size=7*.35) + 
  geom_hline(yintercept=1) + 
  geom_text(aes(1e8,1, label="1%", vjust=-1, hjust = -2), size=7*.35) + 
  scale_x_log10(breaks=myX1, labels=outBreaks) +  
  theme_bw() + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(), 
        axis.line = element_line(color="black"), 
        legend.position = c(0.34,0.8),
        legend.background = element_rect(aes(alpha = 0)), 
        axis.title.y = element_text(size=7), 
        axis.title.x = element_text(size=7), 
        legend.text = element_text(size=6), 
        axis.text.x = element_text(size=6), 
        axis.text.y = element_text(size=6)
  ) +
  xlab("N(t)") + 
  ylab("Average Error (%)") + 
  scale_color_manual(guide="legend", name="", 
                     values=c("blue"="blue", "goldenrod"="goldenrod", "red"="red"), 
                     labels=c("Stochastic model, infinite sites", 
                              "Deterministic model, non-infinite sites", 
                              "Deterministic model, infinite sites"))
errorPlot        
ggsave(file="ErrorPlot.eps", plot=errorPlot, width=87, height=45, units='mm')

ApplicabilityPlot = ggplot() +
  geom_rect(aes(xmin=3e4, xmax=3e5, ymin=0, ymax=1, color="blue", fill="blue", alpha=0.5)) +
  geom_rect(aes(xmin=8e0, xmax=2e1, ymin=1.5, ymax=2.5, color="goldenrod", fill="goldenrod", alpha=0.5)) +
  geom_rect(aes(xmin=8e0, xmax=3e5, ymin=3, ymax=4, color="red", fill="red", alpha=0.5)) +
  geom_rect(aes(xmin=1e0, xmax=3e4, ymin=0, ymax=1, color="blue", fill="blue", alpha=1)) +
  geom_rect(aes(xmin=2e1, xmax=1e10, ymin=1.5, ymax=2.5, color="goldenrod", fill="goldenrod", alpha=1)) +
  geom_rect(aes(xmin=2e1, xmax=3e4, ymin=3, ymax=4, color="red", fill="red", alpha=1)) + 
  geom_text(aes(x=2e2, y=0.5, label="Stochastic model, \ninfinite sites", color = "white"), size=7*.35) + 
  geom_text(aes(x=1.5e5, y=2, label="Deterministic model, \nnon-infinite sites"), size=7*.35) +
  geom_text(aes(x=8e2, y=3.5, label="Deterministic model, \ninfinite sites"), size=7*.35) +
  scale_x_log10(breaks=myX1, labels=outBreaks) +  
  theme_bw() + 
  theme(panel.grid.major=element_blank(), 
        panel.grid.minor=element_blank(), 
        panel.border=element_blank(), 
        axis.line.x = element_line(color="black"), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = "none", 
        axis.title.x = element_text(size=7), 
        axis.text.x = element_text(size=6)
  ) + 
  xlab("N(t)") + 
  ylab("") +
  scale_color_manual(values=c("blue"="blue", "goldenrod"="goldenrod", "red"="red", "white"="white"), 
                     labels=c("Stochastic model, infinite sites", 
                              "Deterministic model, non-infinite sites", 
                              "Deterministic model, infinite sites", 
                              "textLabel")) +
  scale_fill_manual(values=c("blue"="blue", "goldenrod"="goldenrod", "red"="red"), 
                    labels=c("Stochastic model, infinite sites", 
                             "Deterministic model, non-infinite sites", 
                             "Deterministic model, infinite sites"))
ggsave(file="AplicabilityPlot.eps", width=87, height=45, units='mm')
ApplicabilityPlot
dev.off()
mmToIn <- function(in_mm) {return(in_mm/25.4)}

