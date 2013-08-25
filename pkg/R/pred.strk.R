pred.strk <- function (temp,
                       newdata, # only STFDF
                       pred.id= 'tempPred',
                       zcol=1,
                       zero.tol=0,
                       dynamic.cov=c(1,2),
                       static.cov=c(1,2),
                       reg.coef=list(tmean= c(-0.1265044154,0.4051734447,0.4943247727,0.0001837527,-0.0189207588), # Intercept, temp_geo,modis,dem,twi 
                                     tmin = c(-0.9825601517,0.5672140021,0.3344561638, 0.0003119777,-0.0243629638),
                                     tmax = c(1.7873573081,0.350228076, 0.5569091092, 0.0002571338,-0.0012988123) ) [[ 'tmean']] ,
                       vgm.model=list(tmean=vgmST("sumMetric",
                                                  space=vgm( 14.13, "Sph", 5903, 1.933),
                                                  time =vgm(0, "Sph",  0.1, 0),
                                                  joint=vgm(9.06, "Sph", 2054, 0.474),
                                                  stAni=497.9),
                                      tmin = vgmST("sumMetric",
                                                  space=vgm( 22.682, "Sph", 5725, 3.695),
                                                  time =vgm(0, "Sph",  0.1, 0),
                                                  joint=vgm(9.457, "Sph",1888, 1.67),
                                                  stAni=485),
                                     tmax = vgmST("sumMetric",
                                                  space=vgm( 8.31, "Sph", 4930, 2.872),
                                                  time =vgm(0, "Sph",  0.1, 0),
                                                  joint=vgm(11.175, "Sph", 2117, 1.75),
                                                  stAni=527) ) [['tmean']] ,
                       tiling= TRUE,
                       ntiles=64,
                       parallel.processing=TRUE,
                       cpus=3,
                       sp.nmax=18,
                       time.nmax=2,
                       longlat= TRUE,
                       fast= FALSE,
                       computeVar=FALSE,
                       do.cv= TRUE,
                       only.cv = FALSE,
                       out.remove = TRUE,
                       treshold.res=15){
  
  temp <- rm.dupl(temp, zcol,zero.tol)
  gg <- newdata
  time <- gg@time
     
  if(!is.null(dynamic.cov)) {
   dyn<- as.data.frame ( gg@data[,dynamic.cov] )
  }else{ dyn<-NA }
  
  if(!is.null(static.cov)  ) {
    sta <- lapply(static.cov, function(i) rep(gg@sp@data[,i],length(time)) )
    sta <- as.data.frame( do.call(cbind,sta) )
    names(sta)<- names(gg@sp@data[,static.cov])
  }else{ sta<- NA }
  
  df= cbind(dyn,sta)
  df <- df[,colSums(is.na(df))<nrow(df)]
  
  if(nrow(df)==0){
     gg$tlm<-0 
     temp$tres <- temp@data[,zcol]- temp$tlm
     } else {
       
       gg$tlm= reg.coef[1] + as.matrix(df )  %*%  reg.coef[-1]
       gg$tlm<- as.vector(gg$tlm)
       
       gg@sp=as(gg@sp,'SpatialPixelsDataFrame')
       ov <- sapply(1:length(time),function(i) over(temp@sp,as(gg[,i,'tlm'],'SpatialPixelsDataFrame' ) ) )
       ov <- do.call('cbind',ov)
       ov <- as.vector(ov)
       
       t1 <- which(as.Date(index(time[1])) == as.Date(index(temp@time)) )
       t2 <- which(as.Date(index(time[ length(time) ])) == as.Date(index(temp@time)) )
       temp <- temp[,t1:t2]
       
       temp$tlm <- ov
       temp$tres <- temp@data[,zcol]- temp$tlm
       
       # count NAs per stations
       numNA <- apply(matrix(temp@data[,'tres'],
                             nrow=length(temp@sp),byrow=F), MARGIN=1,
                      FUN=function(x) sum(is.na(x)))
       
       # Remove stations out of covariates
       rem <- numNA != length(time)
       temp <- temp[rem,]
#        row.names(temp@sp) <- 1: nrow(temp@sp)
       
       
     }
  
  i_1 <- (1:length(time)) - ceiling(time.nmax/2)
  i_1[i_1<1]=1        
  ip1 <- i_1 + floor(time.nmax/2)
  ip1[ip1>length(time)] <- length(time)  
  
  sfInit ( parallel = parallel.processing , cpus =cpus)
  if(parallel.processing) {
  sfLibrary(gstat)
  sfLibrary(zoo)
  sfLibrary(spacetime)
  sfLibrary(sp)
  sfExport("vgm.model" )
  sfExport( "longlat" )
  sfExport( "computeVar" )
  sfExport( "i_1" )
  sfExport( "ip1" )
  sfExport( "time" )
  }
  
#   

  gg@sp=as(gg@sp,'SpatialPointsDataFrame')
  row.names(gg@sp) = 1:nrow(gg@sp)
  
  remove <- rep(FALSE,length(temp@sp))
  robs=NA
  idsOUT= NA
  
  if(do.cv){
    
    N_POINTS <- length(temp@sp@coords[,1])
    
    if(parallel.processing) {
    sfExport( "temp" )
    sfExport( "N_POINTS" )
    sfExport( "sp.nmax" )   }
    
    
    cv = as.list ( rep(NA, N_POINTS)) 
    
    cv =   sfLapply(1:N_POINTS, function(i) 
    {
      st= temp@sp
      st$dist=spDists(temp@sp,temp@sp[i,] , longlat = longlat)
      
      tmp_st<-st[ order(st$'dist') , ]
      
      local_t= row.names(tmp_st[2:sp.nmax,] ) # remove target station
      
      xxx = as.list ( rep(NA, length(time) ) )
      for( ii in 1:length(time) ) {
        xxx[[ii]]=krigeST(as.formula("tres~1"),
                          data=as(temp[local_t, i_1[ii]:ip1[ii],'tres'],"STSDF"), 
                          newdata=STF(as(temp@sp[i,],"SpatialPoints"),
                                      temp@time[ii],  
                                      temp@endTime[ii]),     
                          modelList=vgm.model,
                          computeVar=computeVar)@data[,1]
        
      } # end of  for
      
      unlist(xxx)  } ) # end of sfLapply
    
    cv <- do.call(rbind,cv)
    cv <- as.vector(cv)
    cv.temp <- temp
    cv.temp$pred.cv <- cv + cv.temp$tlm
    cv.temp$resid.cv <- cv.temp$pred.cv  - cv.temp@data[,zcol]
    
    temp<- cv.temp
    
    bb <- cv.temp
    bb <- as.data.frame(bb)
    bb$sp.ID <-as.character(bb$sp.ID )
    bb$abs.res <-abs(bb$resid.cv )
    bb <- bb[order(bb$abs.res,  decreasing = TRUE),]
    idsOUT <- unique(bb[which( bb$abs.res > treshold.res),]$sp.ID) 
#     idsOUT <- as.numeric (idsOUT) 
    #     stplot( temp[idsOUT,,c('tempc','pred.cv')] , mode='tp')

    if(out.remove==TRUE & length(idsOUT)!=0){
      
      while(length(idsOUT)>0){ 
        
        rm.station = which(row.names(temp@sp)==idsOUT[1])
        remove [rm.station ] <- TRUE
        
        rm.station = which(row.names(cv.temp@sp)==idsOUT[1])
        cv.temp <- cv.temp[ - rm.station , ]
        
        N_POINTS <- length(cv.temp@sp@coords[,1])
        
        if(parallel.processing) {
        sfExport( "cv.temp" )
        sfExport( "N_POINTS" )
        sfExport( "sp.nmax" )  }
        
        
        cv = as.list ( rep(NA, N_POINTS)) 
        
        cv =   sfLapply(1:N_POINTS, function(i) 
        {
          st= cv.temp@sp
          st$dist=spDists(cv.temp@sp,cv.temp@sp[i,] , longlat = longlat)
          
          tmp_st<-st[ order(st$'dist') , ]
          
          local_t= row.names(tmp_st[2:sp.nmax,] ) # remove target station
          
          xxx = as.list ( rep(NA, length(time) ) )
          for( ii in 1:length(time) ) {
            xxx[[ii]]=krigeST(as.formula("tres~1"),
                              data=as(cv.temp[local_t, i_1[ii]:ip1[ii],'tres'],"STSDF"), 
                              newdata=STF(as(cv.temp@sp[i,],"SpatialPoints"),
                                          cv.temp@time[ii],  
                                          cv.temp@endTime[ii]),     
                              modelList=vgm.model,
                              computeVar=computeVar)@data[,1]
            
          } # end of  for
          
          unlist(xxx)  } ) # end of sfLapply
        
        cv <- do.call(rbind,cv)
        cv <- as.vector(cv)
        cv.temp$pred.cv <- cv + cv.temp$tlm
        cv.temp$resid.cv <- cv.temp$pred.cv  - cv.temp@data[,zcol]
        
        bb <- cv.temp
        bb <- as.data.frame(bb)
        bb$sp.ID <-as.character(bb$sp.ID )
        bb$abs.res <-abs(bb$resid.cv )
        bb <- bb[order(bb$abs.res,  decreasing = TRUE),]
        idsOUT <- unique(bb[which( bb$abs.res > treshold.res),]$sp.ID) 
#         idsOUT <- as.numeric (idsOUT) 
        
      } # while
      remove <- temp@sp[remove,]
      robs <- temp[remove,]
      temp <- cv.temp
    }# remove
    
    
    if(only.cv){
      sfStop()
      return(list(pred=NA,cv =cv.temp, out= idsOUT, remst =remove, remobs=robs) )
    }
    
  }else{ 
    cv.temp =NA
    idsOUT=NA}
  
  if (!tiling) {
    dimnames(gg@sp@coords)[[2]] <- c('x','y')
    Mpoint=data.frame(x=mean(gg@sp@coords[,1]),y=mean(gg@sp@coords[,2]) )
    coordinates(Mpoint)=~x+y
    
    temp.local<-temp[,,'tres']
    
    if(parallel.processing) {
    sfExport( "temp.local" )
    sfExport( "gg" )    } 
    
    
    xxx<- sfLapply(1:length(time), function(i) {
      
      obs=as(temp.local[,i_1[i]:ip1[i],'tres'],"STSDF")
      pred.var=1
      if(computeVar==TRUE){ pred.var=2}
        
      krigeST(as.formula("tres~1"),
              data=obs, # [,1:6] # in this case I must to limit for the fist few days
              newdata=STF(as(gg@sp,"SpatialPoints"),
                          temp.local@time[i],    # [3]    
                          #srb@data,  # [3,1]
                          temp.local@endTime[i]),    # [3]    
              modelList=vgm.model,
              computeVar=computeVar)@data[,pred.var]
      
    } )

    res=do.call(cbind,xxx)  
    res= as.vector(res)
    gg$temp.pred<- res + gg$tlm
  
    stfdf <-gg


  }else{ 
    dimnames(gg@sp@coords)[[2]] <- c('x','y')
    xy=as.data.frame(gg@sp)
    xy= xy[row.names(gg@sp),]
    nxy= floor(sqrt(ntiles) )
    xy$xg=as.character( cut(xy$x,nxy,labels=paste("x",1:nxy,sep="") ) )
    xy$yg=as.character( cut(xy$y,nxy,labels=paste("y",1:nxy,sep="") )  )
    #  xy=xy[order(-xy$yg, xy$xg),]
    xy$g=as.factor(paste(xy$xg,xy$yg,sep="") )
    xy$xg=NULL
    xy$yg=NULL
    coordinates(xy) = ~ x+y
    xy$index=1:nrow(xy)
    g_list <- split(xy, xy$g)
    
    # g_list= lapply(g_list,function(i) g[i,])
    
    # mAKE CHUNKS OF INITIAL DATA
    Mpoint=data.frame(x=mean(gg@sp@coords[,1]),y=mean(gg@sp@coords[,2]) )
    coordinates(Mpoint)=~x+y
    Mpoint@proj4string <- gg@sp@proj4string
    st=temp@sp
    
    
#     # Find nearest 500 stations
#     #
#     st$dist=spDists(temp@sp,Mpoint)
#     tmp_st<-st[ order(st$'dist') ,]
#     local1= row.names(tmp_st[1:(sp.nmax+1),] )
    
    
    # Middle point for each chunk
    Mpts= sfLapply(g_list,function(i) {
      Mpoint=data.frame(x=mean(i@coords[,1]),y=mean(i@coords[,2]) )
      coordinates(Mpoint)=~x+y
      Mpoint   })
    
    # Local obs.
#     temp.local <-temp[local1, ,'tres']
    temp.local <-temp[ , ,'tres']
    st=temp.local@sp
    
    if(parallel.processing) {
    sfExport( "temp.local" )
    sfExport( "st" )
    sfExport( "Mpts" )
    sfExport( "g_list" )
    sfExport( "sp.nmax" ) }
    
    res =   sfLapply(1:length(g_list), function(i) 
    {
      pred.var = 1
      if(computeVar==TRUE){ pred.var=2}
      
      st$dist=spDists(temp.local@sp,Mpts[[i]],  longlat = longlat)
      
      tmp_st<-st[ order(st$'dist') ,]
      
      local_t= row.names(tmp_st[1:sp.nmax,] )
      
      xxx = as.list ( rep(NA, length(time) ) )
      for( ii in 1:length(time) ) {
        xxx[[ii]]=krigeST(as.formula("tres~1"),
                          data=as(temp.local[local_t, i_1[ii]:ip1[ii],'tres'],"STSDF"), 
                          newdata=STF(as(g_list[[i]],"SpatialPoints"),
                                      temp.local@time[ii],  
                                      temp.local@endTime[ii]),     
                          modelList=vgm.model,
                          computeVar=computeVar)@data[,pred.var]
      } # end of  for
      
      do.call(cbind,xxx)  } ) # end of sfLapply
    
    res=do.call(rbind,res)
    cc=do.call(rbind,g_list)
    g=cc
    
      stfdf= gg[as.character(g$index),]    
      res <- as.numeric(res)
      stfdf$temp.pred <- res +stfdf$tlm
    
    if(!fast){
      stfdf1= stfdf
      
      dimnames(gg@sp@coords)[[2]] <- c('x','y')
      xy=as.data.frame(gg@sp)
      xy= xy[row.names(gg@sp),]
      nxy= floor(sqrt(ntiles+ ntiles) )
      xy$xg=as.character( cut(xy$x,nxy,labels=paste("x",1:nxy,sep="") ) )
      xy$yg=as.character( cut(xy$y,nxy,labels=paste("y",1:nxy,sep="") )  )
      #  xy=xy[order(-xy$yg, xy$xg),]
      xy$g=as.factor(paste(xy$xg,xy$yg,sep="") )
      xy$xg=NULL
      xy$yg=NULL
      coordinates(xy) = ~ x+y
      xy$index=1:nrow(xy)
      g_list <- split(xy, xy$g)
      
      # g_list= lapply(g_list,function(i) g[i,])
      
      # mAKE CHUNKS OF INITIAL DATA
      Mpoint=data.frame(x=mean(gg@sp@coords[,1]),y=mean(gg@sp@coords[,2]) )
      coordinates(Mpoint)=~x+y
      Mpoint@proj4string <- gg@sp@proj4string
      st=temp@sp

      
      
      # Middle point for each chunk
      Mpts= sfLapply(g_list,function(i) {
        Mpoint=data.frame(x=mean(i@coords[,1]),y=mean(i@coords[,2]) )
        coordinates(Mpoint)=~x+y
        Mpoint   })
      
      # Local obs.
      #     temp.local <-temp[local1, ,'tres']
      temp.local <-temp[ , ,'tres']
      st=temp.local@sp
      if(parallel.processing) {
      sfExport( "temp.local" )
      sfExport( "st" )
      sfExport( "Mpts" )
      sfExport( "g_list" )
      sfExport( "sp.nmax" )
      sfExport( "computeVar" ) }
      
      res =   sfLapply(1:length(g_list), function(i) 
      {
        pred.var = 1
        if(computeVar==TRUE){ pred.var=2}
        
        st$dist=spDists(temp.local@sp,Mpts[[i]] ,  longlat = longlat)
        
        tmp_st<-st[ order(st$'dist') ,]
        
        local_t= row.names(tmp_st[1:sp.nmax,] )
        
        xxx = as.list ( rep(NA, length(time) ) )
        for( ii in 1:length(time) ) {
          xxx[[ii]]=krigeST(as.formula("tres~1"),
                            data=as(temp.local[local_t, i_1[ii]:ip1[ii],'tres'],"STSDF"), 
                            newdata=STF(as(g_list[[i]],"SpatialPoints"),
                                        temp.local@time[ii],  
                                        temp.local@endTime[ii]),     
                            modelList=vgm.model,
                            computeVar=computeVar)@data[,pred.var]
        
        } # end of  for
        
        do.call(cbind,xxx)  } ) # end of sfLapply
      
      res=do.call(rbind,res)
      cc=do.call(rbind,g_list)
      g=cc
      
      stfdf= gg[as.character(g$index),]    
      res <- as.numeric(res)
      stfdf$temp.pred <- res +stfdf$tlm
      
      stfdf<- stfdf[ row.names(stfdf1@sp)  , ]
      
      stfdf$temp.pred <- 0.5*stfdf$temp.pred + 0.5*stfdf1$temp.pred
      
      
    }
    
  } 
  
  stfdf@sp <- as(stfdf@sp, class(newdata@sp) )
   
  stfdf <- STFDF(stfdf@sp,stfdf@time, data.frame(stfdf$temp.pred))
  
  names(stfdf@data) <- pred.id 
  

  
  sfStop()
  return(list(pred=stfdf,cv =cv.temp, out= idsOUT, remst =remove, remobs=robs) )
}
  
    

  
  
  