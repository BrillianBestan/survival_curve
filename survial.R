#install.packages("survival")

setwd("C:\\Users\\Administrator\\Desktop\\her")   
library(survival)
rt=read.table("survivalInput.txt",header=T,sep="\t")
rt$futime=rt$futime/365       
a=rt[,"expression"]<median(rt[,"expression"])
diff=survdiff(Surv(futime, fustat) ~a,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=round(pValue,5)
fit <- survfit(Surv(futime, fustat) ~ a, data = rt)
summary(fit)   
pdf(file="survival.pdf")
plot(fit, lty = 2:3,col=c("red","blue"),xlab="time (day)",ylab="surival rate",
     main=paste("surival curve (p=", pValue ,")",sep=""))
legend("topright", c("MAP3K21 high expression", "MAP3K21 low expression"), lty = 2:3, col=c("red","blue"))
dev.off()
