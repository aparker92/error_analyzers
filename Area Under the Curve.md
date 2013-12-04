Area under the curve


I just wanted to update you guys:
The Quakers made R code that calculated area under the curve for their group. I'll work with that code and fiddle with it so that it can apply to our code and output. I'll update you guys when I translate it into something we can use for our project:

      area=rep(0,length(quanlist))
      xx=1:nrow(TuningE.nu)
      for (j in 1:length(quanlist)){
        area[j]=integrate.xy(x=seq(0,1,length.out=length(xx)), TuningE.nu[xx,j])
      }
      
      jpeg(file=paste("AreaCom(ns=",ns,").jpeg",sep=""),1200,800)
      plot(quanlist,area,pch=19,col=rainbow(length(quanlist)),main=paste("Area under Error Curve(",ns,"successors)"))
      dev.off()
      
      #It took us 6 min to test for 21 parameters
      
      #which parameter does the best?
      min(area)
      quanlist[which(area==min(area))] 
