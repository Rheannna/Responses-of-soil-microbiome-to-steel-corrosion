#exploring the env data
env <- read.csv("../20170728土壤样品微量测定.csv",stringsAsFactors = F)
env = env[,-c(1,3)]
##by
by(env$Cd.mg.kg., env$date, sd)
by(env$Cd.mg.kg., env$date, mean)

##aggregate
aggregate(C ~ date,env[3:36,],mean)
aggregate(C ~ date,env[3:36,],sd)

##多组分组
sd = aggregate(env[,c(2:16)], by = list(env$date), sd)
mean = aggregate(env[,c(2:16)], by = list(env$date), mean)


rnorm(40, mean=1.82, sd =0.2390)

# 统计数据
library(pastecs)
s3 <- stat.desc(env)

mean <- s3["mean", ]
sd <- s3["std.dev", ]
