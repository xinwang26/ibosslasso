#these already written in /Library/Frameworks/R.framework/Versions/3.3/Resources/etc/Renviron, no need to run again, simply keep as reference for future
#Sys.setenv(HADOOP_CONF_DIR="/Users/xinwang/Library/hadoop/etc/hadoop")
#Sys.setenv(HADOOP_PREFIX='/Users/xinwang/Library/hadoop')
#Sys.setenv(HADOOP_HOME='/Users/xinwang/Library/hadoop/')
#Sys.setenv(HADOOP_BIN="/Users/xinwang/Library/hadoop/bin/")
#Sys.setenv(HADOOP_CMD='/Users/xinwang/Library/hadoop/bin/hadoop')
#Sys.setenv(HADOOP_STREAMING='/Users/xinwang/Library/hadoop/share/hadoop/tools/lib/hadoop-streaming-2.7.2.jar')
#Sys.setenv()

#sudo ln -f -s $(/usr/libexec/java_home)/jre/lib/server/libjvm.dylib /usr/local/lib
#export PKG_CONFIG_PATH='/usr/local/lib'
#export LD_LIBRARY_PATH='/usr/local/lib'

library(rJava)
library(Rhipe)
rhinit()

map <- expression({
  lapply(seq_along(map.values),function(r){
    x <- runif(map.values[[r]])
    rhcollect(map.keys[[r]],c(n=map.values[[r]],mean=mean(x),sd=sd(x)))
  })
})
## Create a job object
z <- rhwatch(map, ofolder="/tmp/test", inout=c('lapply','sequence'),
          N=10,mapred=list(mapred.reduce.tasks=0),jobname='test')
## Submit the job
rhex(z)
## Read the results
res <- rhread('/tmp/test')
colres  <- do.call('rbind', lapply(res,"[[",2))
colres

#worked
rhwrite(list(list("a",1),list("b",2)), "/tmp/x")
rhread("/tmp/x")

