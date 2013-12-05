#Function to create the Error Diagram
#Parameters Needed
#  Data in the format of first column = time, second column = magnitude
#  startDate, testStartDate, finishDate = ISO data format
#  modelFunction = function written in R for Eqarthquake analysis
#  parameters = list of parameters

#  To install 'pracma' R package:
#  Within R Studio open the window in the lower right corner of your screen.
#  Open _Packages_ > *Install Packages* > search for 'pracma' and install.
#  or you could just run following code -> install.packages('pracma')

library(pracma)
createErrorDiagram <- function(data, startDate, testStartDate, finishDate, modelFunction, parameters) {

  #Extracts Data
  times = data[,1]
  mags = data[,2]
  
  #Sets time period for training, testing, and entire period
  training.period = as.numeric(testStartDate - startDate)
  testing.period = as.numeric(finishDate - testStartDate)
  entire.period = as.numeric(finishDate - startDate)
  timelist = 2:training.period
  n.training = sum(times<training.period)
  n.test = sum(times<entire.period) - n.training
  
  #Dummy Confidence interval object
  CI.dist=rep(NA,length(timelist))
  CI.list=rep(NA,n.training)
  
  n.events = sapply(timelist,function(x){sum(times<x)})
  
  #Extracts Parameters
  mu.hat = parameters["mu.hat"]
  K.hat = parameters["K.hat"]
  alpha.hat = parameters["alpha.hat"]
  c.hat = parameters["c.hat"]
  p.hat = parameters["p.hat"]
  
  for(KK in 1:length(timelist)){
    CI.dist[KK]=etas.CI(timelist[KK],times[1:n.events[KK]],
                        mags[1:n.events[KK]],m0=1,mu=mu.hat,K=K.hat,alpha=alpha.hat/log(10),c=c.hat,p=p.hat)}
  
  for(KK in 1:length(CI.list)){
    CI.list[KK]=etas.CI(times[1+KK],times[1:(KK)],mags[1:(KK)],
                        m0=1,mu=mu.hat,K=K.hat,alpha=alpha.hat/log(10),c=c.hat,p=p.hat)}
  
  nu.CI=sapply(CI.dist,function(x){mean(CI.list<x)})
  nu.CI=sort(nu.CI,decreasing=TRUE)
  
  xx=seq(1,length(nu.CI))
  yy=length(xx)
  if (length(xx) == 2) {
    print("Not Enough Significant Data")
  } else {
    x = seq(0,1,length.out=yy)
    y = nu.CI[xx]
    #Actually plots the data
    plot(x,y,type="l",xlab=expression(tau),ylab=expression(nu),main="Error Diagram",ylim=c(0,1),col="red",lty=1)
  
    #Requires "pracma" package to use trapz function.
    #Gets area under the curve using trap
    AUC = trapz(x,y)
    print(AUC)
  }
}

#Sample example for how the above function works
#data set is available at https://github.com/aparker92/error_analyzers/blob/master/socal.txt
socal.dat = read.table("/Users/bonghyunkim/Downloads/socal.txt",header=T)
start = ISOdate(1984,1,1,0,0,0)
test.start = ISOdate(2004,6,18,0,0,0)
finish = ISOdate(2010,1,1,0,0,0)
etas.CI <- function (time,t.events,mag.events,m0,mu,K,alpha,c,p) {
  mu+sum(K*10^(alpha*(mag.events-m0))/(time-t.events+c)^p)
}

#Example 1

parameters = c(mu.hat = 0, 
               K.hat = 0,
               alpha.hat = 0,
               c.hat = 0,
               p.hat = 0)

createErrorDiagram(socal.dat, start, test.start, finish, etas.CI, parameters)


#Example 2
parameters = c(mu.hat = .329505837595229, 
               K.hat = .0224702963795154,
               alpha.hat = 1.5839343640414,
               c.hat = .037651249192514,
               p.hat = 1.38508560377488)
               
createErrorDiagram(socal.dat, start, test.start, finish, etas.CI, parameters)
