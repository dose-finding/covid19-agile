## DOSE COMBINATIONS
##dose expansion based on SAFETY FIRST One dose expansion per week
#When efficacy information becomes available FOR efficacy.start PATIENTS, expand based on that
#efficacy expansion has varable cohort size
#option for efficacy expansion only of max safe dose
#option to only use corresponding dose controls

library("rjags")
library("mvtnorm")
library("survival")

AGILE_simulations_combo<-function(z,esc_saf=1,efficacy.start,
                         toxicity,
                         HR,
                         cohort.control, 
                         cohort.trt,cohort.control.eff, 
                         cohort.trt.eff,share.control=1, starting.dose=2, futile,   
                         promise, max.control, HRalt=1.75, 
                         rho0=0.5, prec=0.7, mrec=14 ,
                         max.week=200,max.stages,eff.only.max=1,
                         correlation=0.80,active.arm.decision=20, 
                         target.increase=0.20,delta=0.05,
                         Prior.Pvec=c(var1=0.60,
                                      var2=0.25,
                                      slope1=0.00,
                                      slope2=-0.00,
                                      spacing1=0.075,
                                      p0.control=0.10,
                                      priorEta=10,
                                      overdose=0.25)){
  
  max.eff.ss=efficacy.start+max.stages*cohort.trt.eff
  doses<-(nrow(HR)-1)*(ncol(HR)-1)+1
  
  max.total.ss=max.eff.ss*(doses-1)
  # if esc_saf=0, we wait for efficacy data before expanding at all
  # if esc_saf=1, we use safety data to expand until efficacy data from efficacy.start patients is available
  # efficacy.start = number of patients whose efficacy data is needed before efficacy expansion
  # max.stages = maximum number of stages we can look at efficacy
  # eff.only.max = 1 if when expanding for efficacy we only expand the maximum safe dose
  # eff.only.max = 0 if when expanding for efficacy we can expand if it is safe, even if not maximum
  
  
  # cohort.control.eff            # Cohort Size for the Control Treatment (Efficacy expansion)
  # cohort.trt.eff                # Cohort Size for the Active Doses (Efficacy expansion)
  # doses   Number of Doses

  
  
  # Allocation
  # cohort.control            # Cohort Size for the Control Treatment
  # cohort.trt                # Cohort Size for the Active Doses
  # starting.dose             # Starting dose (Start from Dose 1)
  # 
  # # Efficacy Analysis
  # futile                 # Futility Bound
  # promise                 # Efficacy Bound
  # max.control               # Cap on the Number of Patients on Control Used in Testing
  # HRalt                  # Alternative Hypothesis for Posterior Predictive Probability
  # rho0                   # Prior on the HR
  # prec                     # Assumed Recovery on Control
  # mrec                     # Assumed Median Recovery
  # 
  # # Stopping Criteria
  # max.eff.ss                # Maximum Number of Patients Per Dose
  # max.week                  # Maximum Duration of the Study
  # max.total.ss             # Maximum Number of Patients Per Investigation Arm (All Doses)
  # 
  # 
  # correlation     # Correlation Between Toxicity/Efficacy Profiles
  # 
  set.seed(z)
  
  
  if((nrow(HR)!=nrow(toxicity))|(ncol(HR)!=ncol(toxicity))) stop("Incompatable HR/toxicity specification")
  
  
  
  model1.string <-"
model {
for (i in 1:m1){
logit(p[i]) <- alpha0[1] + alpha1[1] * sdose1[i]
odds1[i]<- p[i]/(1-p[i])
}

for (j in 1:m2){
logit(q[j]) <- alpha0[1] + alpha1[2] * sdose2[j]
odds2[j]<- q[j]/(1-q[j])
}

for (i in 1:m1){
for (j in 1:m2){
odds0[i,j]<-odds1[i] + odds2[j] + odds1[i] * odds2[j]
odds[i,j]<-odds0[i,j]*exp(eta * sdose1[i] * sdose2[j] )
tox[i,j]<-odds[i,j]/(1+odds[i,j])
s[i,j] ~ dbin(tox[i,j], n[i,j])		
}
}
eta~dnorm(0,priorEta)
for(t in 1:2){
theta[1:2, t] ~ dmnorm(priorMean[1:2, t],
priorPrec[1:2,
1:2,
t])
## extract actual coefficients
alpha0[t] <- theta[1, t]
alpha1[t] <- exp(theta[2, t])
}
}
"
model1.spec<-textConnection(model1.string)

find.x <- function(ptox, alpha ) {
  alpha<-matrix(alpha,ncol=2)
  x <- (qlogis(ptox)-(alpha[,1]))/alpha[,2]
  return( x )
}


#Prior parameters to use
# var1<-Prior.Pvec[1]
# var2<-Prior.Pvec[2]       #  0.35 is calibrated
# slope<-Prior.Pvec[3]  # -0.25 calibrated
# spacing<-Prior.Pvec[4]
# p0.control<-Prior.Pvec[5]
# overdose<-Prior.Pvec[6]
#####

#Prior parameters to use
# var1<-1.00
# var2<-var3<-0.25
# slope1<-0.00
# slope2<--0.25
# spacing1<-spacing2<-0.10
# p0.control<-0.10
# priorEta<-30
# overdose<-0.25



#Prior parameters to use
var1<-Prior.Pvec[1]
var2<-var3<-Prior.Pvec[2] 
slope1<-Prior.Pvec[3]
slope2<-Prior.Pvec[4]
spacing1<-spacing2<-Prior.Pvec[5]
p0.control<-Prior.Pvec[6]
priorEta<-Prior.Pvec[7]
overdose<-Prior.Pvec[8]
#####

doses1<-c(1:(nrow(toxicity)))     # including control  (agent B)
doses2<-c(1:(ncol(toxicity)))       # including control  (agent A)



doses<-c(1:doses)

#combo

p.tox0.1<-c(p0.control/2,p0.control/2 + spacing1* seq(1,length(doses1)-1))
p.tox0.2<-c(p0.control/2,p0.control/2 + spacing2* seq(1,length(doses2)-1))

priorMean<-mat.or.vec(2,2)
priorMean[,1]<-c(log((p0.control)/2/(1-(p0.control)/2)),slope1)  
priorMean[,2]<-c(log((p0.control)/2/(1-(p0.control)/2)),slope2)  


priorPrec<-array(0,dim=c(2,2,2))
priorVar<-array(0,dim=c(2,2,2))
priorVar[,,1]<-matrix(c(var1,0.0,0.0,var2),2,2)  
priorVar[,,2]<-matrix(c(var1,0.0,0.0,var3),2,2)  
priorPrec[,,1]<-solve(priorVar[,,1])
priorPrec[,,2]<-solve(priorVar[,,2])
priorEta<-priorEta


alpha.prior.plug.1<-c(priorMean[1,1],exp(priorMean[2,1]+diag(priorVar[,,1])[2]/2))
alpha.prior.plug.2<-c(priorMean[1,2],exp(priorMean[2,2]+diag(priorVar[,,2])[2]/2))
sdose1<-find.x(p.tox0.1,alpha=alpha.prior.plug.1) # finding the skeleton
sdose2<-find.x(p.tox0.2,alpha=alpha.prior.plug.2) # finding the skeleton




blrm.covid.combo<-function(n,s,nextdose,sdose1,sdose2,priorMean,priorPrec,priorEta,overdose=0.25,delta=0.05,iter=10000){
  
  toxicity<-prob.next<-mat.or.vec(length(sdose1),length(sdose2))
  
  model1.spec<-textConnection(model1.string)
  mydata <- list(n=n,s=s,m1=length(sdose1),m2=length(sdose2),sdose1=sdose1,sdose2=sdose2,priorMean=priorMean,priorPrec=priorPrec,priorEta=priorEta)
  jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iter,quiet=TRUE)
  update(jags, iter,progress.bar="none")
  tt<-jags.samples(jags,c('alpha0','alpha1','eta'),iter,progress.bar="none")
  
  a01<-tt$alpha0[1,,]
  a02<-tt$alpha0[1,,]
  a11<-(tt$alpha1[1,,])
  a12<-(tt$alpha1[2,,])
  e<-(tt$eta[1,,])
  
  
  for (j in 1:1){
    for (k in 1:1){
      help1 <- a01 + a11 * sdose1[j]
      p1<-exp(help1)/(1+exp(help1))
      help2 <- a02 + a12 * sdose2[k]
      p2<-exp(help2)/(1+exp(help2))
      odds1<-p1/(1-p1)
      odds2<-p2/(1-p2)
      odds0<-odds1 + odds2 + odds1*odds2
      odds<-odds0 * exp(e * sdose1[j] * sdose2[k])
      control.distr<-odds/(1+odds)
    }
  }
  
  for (j in 2:length(sdose1)){
    for (k in 2:length(sdose2)){
      help1 <- a01 + a11 * sdose1[j]
      p1<-exp(help1)/(1+exp(help1))
      help2 <- a02 + a12 * sdose2[k]
      p2<-exp(help2)/(1+exp(help2))
      odds1<-p1/(1-p1)
      odds2<-p2/(1-p2)
      odds0<-odds1 + odds2 + odds1*odds2
      odds<-odds0 * exp(e * sdose1[j] * sdose2[k])
      distr<-odds/(1+odds)
      toxicity[j,k]<-mean(distr)
      y<-distr-control.distr
      prob.next[j,k]<-mean(y <=(target.increase+delta) & (y>=target.increase-delta))
      if(mean(y>=(target.increase+delta))>overdose){
        prob.next[j,k]<-0
      }
    }
  }
  
  for (j in 2:length(sdose1)){
    for (k in 2:length(sdose2)){
      if((j>(nextdose[1]+1))| (k>(nextdose[2]+1)) | ((j>=(nextdose[1]+1))& (k>=(nextdose[2]+1)))){
        prob.next[j,k]<-0
      }
    }
  }
  
  return( prob.next)
}




cox_decide_postprob <- function(control,trt,HR=1.5,p1=0.5) {
  #control = control data
  #trt = treatment data
  #HR = HR for alternative hypothesis
  #p1 = prior weight to alternative hypothesis
  data <- rbind(control,trt)
  data$trt <- rep(c(0,1),c(dim(control)[1],dim(trt)[1]))
  data$offset <- data$trt * log(HR)
  fit1 <- coxph(Surv(time,event)~1,data=data)
  fit2 <- coxph(Surv(time,event)~offset(offset),data=data)
  PPP <- exp(fit2$loglik)*p1/(exp(fit2$loglik)*p1 + exp(fit1$loglik)*(1-p1))
  if (is.na(PPP) | is.nan(PPP)) error_dump <<- list(control=control,trt=trt)
  return(PPP)
}

simulate_single <- function(n=1,u,HR,prec=0.7,mrec=14) {
  #n: sample size for a given arm
  #HR: hazard ratio for this arm 
  #prec: proportion recovering by 28 days under the null i.e. HR=1
  #mrec: median recovery time under the null 
  if (mrec > 28 & prec > 0.5) stop("Incompatible null")
  if (mrec < 28 & prec < 0.5) stop("Incompatible null")
  
  #Find the Weibull that ensures S(28)=1-prec and S(mrec)=15
  alp <- log(-log(1-prec)/log(2))/log(28/mrec)
  lam <- log(2)/mrec^alp
  lamH <- lam*HR
  #Simulate from this Weibull
  t <- (-1/lamH * log(u))^(1/alp)
  #Assume all censoring occurs at 28 days
  x <- pmin(t,28) #Follow-up time
  d <- 1*(x<28)  #Event indicator
  data.frame(time=x,event=d)
}


matrix2vector<-function(MAT){
  rows<-nrow(MAT)
  cols<-ncol(MAT)
  VEC<-c(MAT[1,1])
  for(i in 2:rows){
    for(j in 2:cols){
      VEC<-c(VEC,MAT[i,j])
    }
  }
  VEC
}

vector2matrix<-function(VEC,rows,cols){
  MAT<-matrix(rep(0,rows*cols),nrow=rows,ncol=cols)
  MAT[1,1]<-VEC[1]
  INDEX<-1
  for(i in 2:rows){
    for(j in 2:cols){
      INDEX<-INDEX+1
      MAT[i,j]<-VEC[INDEX]
    }
  }
  MAT
}

comat<-matrix(nrow=length(doses),ncol=2)
comat[1,]<-c(1,1)
INDEX<-1

for(i in 2:(length(doses1))){
  for(j in 2:(length(doses2))){
    INDEX<-INDEX+1
    comat[INDEX,]<-c(i,j)
  }
}

vector2matrix_index<-function(INDEX){
  as.vector(comat[INDEX,])
}

matrix2vector_index<-function(INDEX_VECTOR){
  which((comat[,1]==INDEX_VECTOR[1])&(comat[,2]==INDEX_VECTOR[2]))
}


toxicity<-matrix2vector(toxicity)
HR<-matrix2vector(HR)



bad.eff<-bad.tox<-bad.both<-bad.doses<-good.doses<-unknown.doses<-c()
selection<-rep(NA,length(doses))

nextdose<-starting.dose
stop<-0


n.tox<-rep(0,length(doses))
s.tox<-rep(0,length(doses))



p<-array(0,dim=c(sdose1,sdose2,iter=10000))



prob.next<-mat.or.vec(length(doses1),length(doses2))


response.eff.all<-vector("list",length(doses))
time<-vector("list",length(doses))
checked<-vector("list",length(doses))
control_vec<-c()

for(j in 1:(length(doses))){
  response.eff.all[[j]]<-data.frame()
}


week<-1
add.toxicity<-1
add.efficacy<-0
stop<-0

###############
## trial sim ##
###############

#study can only go on for a max of max.week weeks
while(week<=max.week){

  if(add.toxicity==1){

    #if we add a new dose (and looking at tox), then we need a new control group
    response.eff <- data.frame()
    #patients on control
    for (i in 1:cohort.control){
      st.sigma<-rbind(c(1,correlation),c(correlation,1))
      profile<-rmvnorm(1, mean = c(0,0), sigma = st.sigma)
      profile<-pnorm(profile)
      response.tox<-as.numeric(profile[2]<toxicity[1])
      response.eff<-rbind(response.eff,simulate_single(n=1,u=profile[1],HR=HR[1],prec=prec,mrec=mrec))
      n.tox[1]<-n.tox[1]+1
      s.tox[1]<-s.tox[1]+ response.tox
      control_vec<-c(control_vec,nextdose)
    }
    
    response.eff.all[[1]]<-rbind(response.eff.all[[1]],response.eff)
    
    time[[1]]<-rbind(time[[1]],as.matrix(rep(week,cohort.control),1,cohort.control))
    
    response.eff <- data.frame()
    
    #patients in treatment group
    for (i in 1:cohort.trt){
      st.sigma<-rbind(c(1,correlation),c(correlation,1))
      profile<-rmvnorm(1, mean = c(0,0), sigma = st.sigma)
      profile<-pnorm(profile)
      response.tox<-as.numeric(profile[2]<toxicity[nextdose])
      response.eff<-rbind(response.eff,simulate_single(n=1,u=profile[1],HR=HR[nextdose],prec=prec,mrec=mrec))
      n.tox[nextdose]<-n.tox[nextdose]+1
      s.tox[nextdose]<-s.tox[nextdose]+ response.tox
    }
    # browser()
    response.eff.all[[nextdose]]<-rbind(response.eff.all[[nextdose]],response.eff)
    time[[nextdose]]<-rbind(time[[nextdose]],as.matrix(rep(week,cohort.trt),1,cohort.trt))
    
    
  }
  
  if(add.efficacy==1){

    #if we are adding a new cohort looking at efficacy for the dose  
    for (w in 1:length(which.add)){
      
      q<-which.add[w]
      response.eff <- data.frame()
      for (i in 1:cohort.control.eff){
        st.sigma<-rbind(c(1,correlation),c(correlation,1))
        profile<-rmvnorm(1, mean = c(0,0), sigma = st.sigma)
        profile<-pnorm(profile)
        response.tox<-as.numeric(profile[2]<toxicity[1])
        response.eff<-rbind(response.eff,simulate_single(n=1,u=profile[1],HR=HR[1],prec=prec,mrec=mrec))
        n.tox[1]<-n.tox[1]+1
        s.tox[1]<-s.tox[1]+ response.tox
        control_vec<-c(control_vec,q)
      }
      
      response.eff.all[[1]]<-rbind(response.eff.all[[1]],response.eff)
      time[[1]]<-rbind(time[[1]],as.matrix(rep(week,cohort.control.eff),1,4))
      
      response.eff <- data.frame()
      for (i in 1:cohort.trt.eff){
        st.sigma<-rbind(c(1,correlation),c(correlation,1))
        profile<-rmvnorm(1, mean = c(0,0), sigma = st.sigma)
        profile<-pnorm(profile)
        response.tox<-as.numeric(profile[2]<toxicity[q])
        response.eff<-rbind(response.eff,simulate_single(n=1,u=profile[1],HR=HR[q],prec=prec,mrec=mrec))
        n.tox[q]<-n.tox[q]+1
        s.tox[q]<-s.tox[q]+ response.tox
      }
      
      response.eff.all[[q]]<-rbind(response.eff.all[[q]],response.eff)
      time[[q]]<-rbind(time[[q]],as.matrix(rep(week,cohort.trt.eff),1,4))
    }
  }

  #which dose to escalate to (only if safe, this is zero if not safe)
  nextdose<-vector2matrix_index(nextdose)
  mat.prob.next<-blrm.covid.combo(vector2matrix(n.tox,cols=length(doses2),rows=length(doses1)),vector2matrix(s.tox,cols=length(doses2),rows=length(doses1)),nextdose=nextdose,sdose1=sdose1,sdose2=sdose2,priorMean=priorMean,priorPrec=priorPrec,priorEta=priorEta)
  prob.next<-matrix2vector(mat.prob.next)


  if(all(prob.next==0)){

    stop<-1
    break() 
  }else{
    max.safe.dose<-max(which(prob.next>0))
    nextdose<-which.max(prob.next) #(no skipping doses)
  }
  if(esc_saf==0){
    #if we already have patients on the next dose, we don't go to that dose next
    if(n.tox[nextdose]>0){
      add.toxicity<-0
    }else{
      add.toxicity<-1 #we do go to that dose next
    }
  }else{
    
    #browser()
    if((length(time[[nextdose]])<max.eff.ss)&(prob.next[nextdose]>0)){ #if less than max ss
      if(length(time[[nextdose]])>0){
        if(sum((time[[nextdose]])<=week-3)>=efficacy.start){
          add.toxicity<-0
        }else{
          add.toxicity<-1
        }
      }else{
        add.toxicity<-1
      }
      
    }else{
      add.toxicity<-0
    }
  }
  
  
  
  add.efficacy<-0 
  which.add<-c()
  for (dose in 2:length(doses)){
    if(!(dose %in% bad.doses | dose %in% good.doses)){
      if(length(which(time[[dose]]<=(week-3)))>=efficacy.start){##if we have efficacy data for at least 12 patients
        if(length(checked[[dose]])==0){
          condition.check<-all(checked[[dose]]!= unique(time[[dose]])[which(unique(time[[dose]])<=(week-3))])  
        }else{
          if(length(checked[[dose]])==length(unique(time[[dose]])[which(unique(time[[dose]])<=(week-3))])){
            condition.check<-  !all(checked[[dose]]==unique(time[[dose]])[which(unique(time[[dose]])<=(week-3))])
          }else{
            condition.check<- TRUE
          }
        }
        #if condition.check is true, this means we have new cohorts that have ended efficacy window
        #that means we have new information so can make decisions
        if(condition.check){ 
          
          end.week<-week

          checked[[dose]]<-    unique(time[[dose]])[which(unique(time[[dose]])<=(week-3))]
          #for patients whose efficacy window has closed, which week did they enter the trial?
          if(share.control==1){
            if(length(which(time[[1]]<=(week-3)))<=max.control){ 
              #if the number in control group with efficacy window closed is less than pre-specified maximum
              PPP <- cox_decide_postprob(control=response.eff.all[[1]][which(time[[1]]<=(week-3)),],trt=response.eff.all[[dose]][which(time[[dose]]<=(week-3)),],HR=HRalt,p1=rho0)
            }else{
              #only look at the last max.control patients
              index<-which(time[[1]]<=(week-3))[(length(which(time[[1]]<=(week-3)))-max.control+1):length(which(time[[1]]<=(week-3)))]
              PPP <- cox_decide_postprob(control=response.eff.all[[1]][index,],trt=response.eff.all[[dose]][which(time[[dose]]<=(week-3)),],HR=HRalt,p1=rho0)
            }
          }else{

            control.index<-which((time[[1]]<=(week-3))&(control_vec==dose))
            PPP <- cox_decide_postprob(control=response.eff.all[[1]][control.index,],trt=response.eff.all[[dose]][which(time[[dose]]<=(week-3)),],HR=HRalt,p1=rho0)
          }
          
          #decision on whether to escalate
          #if posterior prob is in between bounds and the dose is safe
          if(PPP > futile & PPP < promise & prob.next[dose]>0 ){
            #if max number of patients reached
            if(length(time[[dose]])>=max.eff.ss){
              #browser()
              unknown.doses<-unique(c(unknown.doses,dose) )
              #otherwise go to efficacy
            }else{
              if(eff.only.max==1){
                if(dose==max.safe.dose){
                  add.efficacy<-1   
                  which.add<-c(which.add,dose)
                }
              }else{
                add.efficacy<-1   
                which.add<-c(which.add,dose)
              }
            }
          }else{
            if(PPP<futile | (prob.next[dose]==0 & sum(n.tox[2:length(doses)])>=active.arm.decision ) ){   
              bad.doses<-c(bad.doses,dose)
              unknown.doses<-unknown.doses[-which(unknown.doses==dose)]
              if(PPP<futile & !(prob.next[dose]==0 & sum(n.tox[2:length(doses)])>=active.arm.decision ) ){
                bad.eff<-c(bad.eff,dose)
              }
              if(!(PPP<futile) & (prob.next[dose]==0 & sum(n.tox[2:length(doses)])>=active.arm.decision ) ){
                bad.tox<-c(bad.tox,dose)
              }
              if(PPP<futile & (prob.next[dose]==0 & sum(n.tox[2:length(doses)])>=active.arm.decision ) ){
                bad.both<-c(bad.both,dose)
              }
            }else{
              if(PPP>promise & prob.next[dose]>0 ){
                good.doses<-c(good.doses,dose) 
                unknown.doses<-unknown.doses[-which(unknown.doses==dose)]
              }
            }
          }
        }# condition check
      } #if efficacy window has ended for some patients
    } #if dose isn't labelled as good or bad
  } # for dose

  cond.1<-((length(unknown.doses)+length(bad.doses)+length(good.doses))==(length(doses)-1))
  cond.2<-week>=max.week
  cond.3<-sum(n.tox[2:length(doses)])>=max.total.ss
  cond.4<-(add.efficacy==0 & add.toxicity==0 & max(unlist(time))<week-3)
  if(cond.1 | cond.2 |cond.3 |cond.4){
    break()
  }
  

  week<-week+1
}

samples<-  n.tox

if(stop==0){

  selection[bad.eff]<-(-1)
  selection[bad.tox]<-(-2)
  selection[bad.both]<-(-3)
  selection[good.doses]<-1
  selection[unknown.doses]<-0
}

selection[is.na(selection)]<-99

if(!exists("end.week")){
  end.week<-week
}
browser()
return(c(selection,samples,week))
}