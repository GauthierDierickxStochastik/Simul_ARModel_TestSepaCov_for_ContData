

library(expm)
library(foreach)
library(doParallel)


space_array= c(10)
#space_array = c(4,6,8,10,12,14,20)
rejection_array = 1:6*0
rho=0.5


for(Space in space_array){
  
  Time = 50
  boots_l = 3
  N =100
  bootstrap_reps = 100
  dep_vec = c(rep(boots_l:1)/boots_l, 1:(N-boots_l)*0)
  weight_matrix = sqrtm(toeplitz(dep_vec))
  sp =c( 3, 2, 1, 1) #third entry is c
  separability_parameters = sp
  Psi_mat = matrix(data = NA, nrow=Space, ncol=Space)
  for(s1 in 1:Space){
    for(s2 in 1:Space){
      Psi_mat[s1,s2] = exp(-25*(s1-s2)^2/(Space-1))
    }
  }
  
  
  C = matrix(data =NA, nrow =(Space*Time), ncol =(Space*Time))
  help_matrix = matrix(data =NA, nrow=Time, ncol = Space)
  for(t in 1:Time){
    for(s in 1:Space){
      help_matrix[t,s] =
        sp[4]^2/sqrt(sp[1]*abs(t-1)/Time+1)*exp(-sp[2]^2*(s-1)^2
                                                
                                                /((Space-1)^2*(sp[1]*abs(t-1)/Time+1)^(sp[3])))
    }
  }
  for(s in 1:Space){
    for(ss in 1:Space){
      for(t in 1:Time){
        for(tt in 1:Time){
          C[((s-1)*Time+t),((ss-1)*Time+tt)] =
            help_matrix[(abs(t-tt)+1),(abs(s-ss)+1)]
        }
      }
    }
  }
  Sq_Cov = Re(sqrtm(C))
  #Sq_Covold <- sqrtm(
  #  Cov_Matrix ( separability_parameters, Space = Space, Time= Time )
  #)
  

  parallel::detectCores()
  
  n.cores <- parallel::detectCores() - 1
  
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "FORK"
  )
  
  #check cluster definition (optional)
  print(my.cluster)
  
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  #check if it is registered (optional)
  foreach::getDoParRegistered()
  
  #how many workers are available? (optional)
  foreach::getDoParWorkers()
  
  
  # Stop cluster
  parallel::stopCluster(cl = my.cluster)
  
 
  
  separability_parameters <- function( a, c, beta , sigma  ){
    c(a, c, beta , sigma ) }
  
  





Cov_Model <- function( s, u, t, w, separability_parameters, Space,
                       Time )
{
  (separability_parameters[4]^2/(separability_parameters[1] * abs( t -
                                                                     w) / Time  + 1
  ) ^(1/2) ) * exp( - separability_parameters[2]^2 *(abs( s - u ) /
                                                       (Space-1) )^2/(
                                                         separability_parameters[1] * abs( t - w ) / Time+
                                                           1)^(separability_parameters[3] ) )
}




#TraceApprox <- PartSpace %x% PartTrace / TraceMeasure( Data )


Cov_Matrix <- function( separability_parameters, Space, Time ){
  Dim <- (Space * Time)
  HelpMatrix <- matrix( NA,nrow = Dim, ncol = Dim )
  #return(HelpMatrix)
  Outcome <- integer(7)
  
  parallel::detectCores()
  
  n.cores <- parallel::detectCores() - 1
  
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "FORK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  Output <- foreach( x = 0: ( Dim^2 - 1 ), .combine = 'c' ) %dopar% {
    #calculate quotient and rest with modulo
    RowNr <- (x %/%  Dim ) + 1
    Rest <- (x %% Dim) + 1
    s <-  ( ( RowNr - 1 ) %/% Time ) + 1 #/ Space
    t <- ( ( RowNr - 1 ) %% Time ) + 1  #/ Time
    u <- ( ( Rest - 1 ) %/% Time )  + 1 #/ Space
    w <- ( ( Rest - 1 ) %% Time ) + 1 #/ Time
    Var<- c(x+1, RowNr, Rest, s , t , u, w)
    Outcome<-rbind( Outcome, Var )
    #}
    #return(Outcome)
    #}
    #HelpMatrix[RowNr,Rest] <-
    Cov_Model( s, u, t, w, separability_parameters, Space, Time )
  }
  # Stop cluster
  parallel::stopCluster(cl = my.cluster)
  
  return( matrix ( Output, Space * Time, byrow = TRUE ) )
}







DataGenerationList<-function( N , separability_parameters , Space,
                              Time, Sq_Cov){
  std_vec <-  matrix(data = rnorm( Space * Time*(N+100), 0, 1 ),
                     ncol=(N+100))
  G =  Sq_Cov %*% std_vec
  for(i in 2:((N+100))){
    G[,i]=G[,i]+rho*G[,(i-1)]
  }
  G=G[,101:(N+100)]
  G <- G - rowMeans(G)
  return(G)
}




#
#         CENTER THE DATA
#


MeanOfData<-function( Data )
{
  Mean<- ( Reduce('+', Data) / length(Data) )
  return(Mean)
}

CenteringFunction<-function( X, M )
{
  Output <- ( X - M )
  return(Output)
}

CenteringOfData<-function( Data )
{
  Mean <- MeanOfData( Data )
  OutputData <-lapply( Data, CenteringFunction, M = Mean )
  return(OutputData)
}




#
#           CALCULATE TRACE
#
#         AND
#
#           CALCULATE PARTIAL TRACE
#


TraceMeasure <-function( Data ){
  DummyList <- sapply( Data, function(x) sum( x^2) )
  Output<- Reduce('+', DummyList)
  return( Output)
}

PartialTraceFreq <- function( TwoDimArr ){
  PartialTracesMatFreq<- TwoDimArr %*% t( TwoDimArr)
  return(PartialTracesMatFreq)
}

PartialTraceFreqDataMat <- function( Data ){
  FreqData<-lapply(Data, FUN=PartialTraceFreq )
  PartTraceFreq <- Reduce("+", FreqData)
  return(PartTraceFreq)
}



PartialTraceTime <- function( TwoDimArr ){
  PartialTracesMatTime<- t(TwoDimArr ) %*%  TwoDimArr
  return(PartialTracesMatTime)
}

PartialTraceTimeDataMat <- function( Data ){
  TimeData<-lapply(Data, FUN=PartialTraceTime )
  PartTraceTime <- Reduce("+", TimeData)
  return(PartTraceTime)
}


PartSpParallelold<-function( Matrix, Space , Time ){
  parallel::detectCores()
  n.cores <- parallel::detectCores() - 1
  my.cluster <- parallel::makeCluster(
    n.cores,
    type = "FORK")
  doParallel::registerDoParallel(cl = my.cluster)
  Output <-foreach( i = 0: ( Space^2 - 1 ), .combine = 'c' ) %dopar%
    {
      sum( diag (
        Matrix[
          (
            ( i %/% Space ) * Time + 1
          )
          :
            (
              ( i %/% Space + 1 ) * Time
            )
          ,
          (
            ( i %% Space ) * Time + 1
          )
          :
            (
              ( i %% Space + 1 ) * Time
            )
        ]
      )
      )
    }
  # Stop cluster
  parallel::stopCluster(cl = my.cluster)
  return( matrix( Output, Space, byrow = TRUE ) )
}
################
PartSpParallel<-function( Matrix, Space , Time ){
  B = matrix(data = rep(0, Space^2), nrow=Space)
  for(r in 1:Space){
    for(i in 1:Space){
      B[i,r] =
        sum(sum(diag(Matrix[((i-1)*Time+1):(i*Time),((r-1)*Time+1):(r*Time)])))
    }
  }
  return(B)
}







#Statistic

TestStat<-function( Data, Space, Time )
{
  EmpCovOp <- Data%*%t(Data)/N
  
  PartSpace <- PartSpParallel( EmpCovOp, Space= Space, Time= Time)
  #PartTraces <- PartSpParallel(
  #  TurnSpaceTimeToTimeSpace (EmpCovOp, Space = Space, Time = Time ),
  #  Space = Time, Time= Space)
  h = matrix(data=Data, ncol=(N*Space))
  PartTrace=h%*%t(h)/N
  TraceApprox <- PartSpace %x% PartTrace / sum(diag(EmpCovOp))
  
  Output <- max(abs( EmpCovOp - TraceApprox ))
  
  return(Output )
}




BootTestStat<-function( dummy)
{
  EmpCovOp <- Data%*%t(Data)/N
  
  PartSpace <- PartSpParallel( EmpCovOp, Space= Space, Time= Time)
  #PartTrace <- PartSpParallel(
  #  TurnSpaceTimeToTimeSpace (EmpCovOp, Space = Space, Time = Time ),
  #  Space = Time, Time= Space)
  h = matrix(data=Data, ncol=(N*Space))
  PartTrace=h%*%t(h)/N
  TraceApprox <- PartSpace %x% PartTrace / sum(diag(EmpCovOp))
  Z = weight_matrix%*% rnorm(N,0,1)
  Z_root = sqrt(as.complex(Z))
  Data_comp =  t(t(Data) *matrix(data=rep(Z_root,Space*Time), nrow=N))
  C_boot = Re(Data_comp%*%t(Data_comp)/N)
  EmpCovOpboot <- TraceApprox + C_boot
  
  PartSpaceboot <- PartSpParallel( EmpCovOpboot, Space= Space, Time=
                                     Time)
  h = matrix(data=Data_comp, ncol=(N*Space))
  PartTraceboot=Re(h%*%t(h)/N)+PartTrace
  TraceApproxboot <- PartSpaceboot %x% PartTraceboot /
    sum(diag(EmpCovOpboot))
  return(max(abs(C_boot+TraceApprox-TraceApproxboot)))
}


#######


number_rejections=0

for(k in 1:bootstrap_reps){
  Data = DataGenerationList( N , separability_parameters , Space,
                             Time, Sq_Cov)
  true_stat = TestStat( Data, Space, Time )
  quantsize = 400
  brownarray=array(0, dim=quantsize)
  helparray = mclapply(brownarray, BootTestStat, mc.cores =
                         getOption("mc.cores",7))
  for(quant in 1:quantsize){
    brownarray[quant]=helparray[[quant]]
  }
  sortbrown =sort(brownarray, decreasing = FALSE)
  upperquantnf=(1-0.05)*quantsize
  
  if(sortbrown[upperquantnf]<true_stat){number_rejections=number_rejections+1}
  print(c(Space,bootstrap_reps,k,number_rejections))
}

rejection_array[(Space-2)/2] = number_rejections/bootstrap_reps

}

