#SINGLE AGENT

##dose expansion based on SAFETY FIRST One dose expansion per week
#When efficacy information becomes available FOR efficacy.start PATIENTS, expand based on that
#efficacy expansion has varable cohort size
#option for efficacy expansion only of max safe dose
#option to only use corresponding dose controls
library("rjags")
library("mvtnorm")
library("survival")


AGILE_simulations_singleagent<-function(z,esc_saf=1,efficacy.start,doses,HR,
                         toxicity,cohort.control, 
                         cohort.trt,cohort.control.eff, 
                         cohort.trt.eff,share.control=1, starting.dose=2, futile,   
                         promise, max.control, HRalt=1.75, 
                         rho0=0.5, prec=0.7, mrec=14 ,
                         max.week=200,max.stages,eff.only.max=1,
                         correlation=0.80,active.arm.decision=20, 
                         target.increase=0.20,delta=0.05,
                         Prior.Pvec=c(var1=1.40, var2=0.35 ,slope=(-0.25),
                                      spacing=0.125, p0.control=0.10,
                                      overdose=0.25)){
  
  max.eff.ss=efficacy.start+max.stages*cohort.trt.eff
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
  set.seed(z)

  
  # Allocation
  # cohort.control            # Cohort Size for the Control Treatment
  # cohort.trt                # Cohort Size for the Active Doses
  # starting.dose             # Starting dose (Start from Dose 1) - labelles as 2 becasue ctrl is 1
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
  
  
  
  if(length(HR)!=doses|length(toxicity)!=doses) stop("Doses incompatable with HR/toxicity specification")
  
  model1.string <-"
model {
for (i in 1:m){
logit(p[i]) <- alpha0[1] + alpha1[1] * sdose[i]
s[i] ~ dbin(p[i], n[i])		
}

theta[1:2] ~ dmnorm(priorMean[1:2],
priorPrec[1:2,
1:2])
## extract actual coefficients
alpha0<- theta[1]
alpha1 <- exp(theta[2])
}
"
model1.spec<-textConnection(model1.string)

find.x <- function(ptox, alpha ) {
  alpha<-matrix(alpha,ncol=2)
  x <- (qlogis(ptox)-(alpha[,1]))/alpha[,2]
  return( x )
}

#Prior parameters to use
var1<-Prior.Pvec[1]
var2<-Prior.Pvec[2]       #  0.35 is calibrated
slope<-Prior.Pvec[3]  # -0.25 calibrated
spacing<-Prior.Pvec[4]
p0.control<-Prior.Pvec[5]
overdose<-Prior.Pvec[6]
#####

# Other Calibrated Prior Parameters
# var1<-1.40
# var2<-0.40       
# slope<-(-0.40)   
# spacing<-0.15  
# p0.control<-0.10
# overdose<-0.30
#####

doses<-c(1:doses)
#target.increase<-0.20
#delta<-0.05


p.tox0<-c(p0.control,p0.control + spacing* seq(1,length(doses)-1))
priorMean<-c(log(p0.control/(1-p0.control)),slope)  
priorVar<-matrix(c(var1,0.0,0.0,var2),2,2)  
priorPrec<-solve(priorVar)
alpha.prior.plug<-c(priorMean[1],exp(priorMean[2]+diag(priorVar)[2]/2))
sdose<-find.x(p.tox0,alpha=alpha.prior.plug) # finding the skeleton



blrm.covid<-function(n,s,sdose,priorMean,priorPrec,overdose=0.25,delta=0.05,iter=10000){
  p<-mat.or.vec(iter,length(doses))
  model1.spec<-textConnection(model1.string)
  mydata <- list(n=n,s=s,m=length(sdose),sdose=sdose,priorMean=priorMean,priorPrec=priorPrec)
  jags <- jags.model(model1.spec,data =mydata,n.chains=1,n.adapt=iter,quiet=TRUE)
  update(jags, iter,progress.bar="none")
  tt<-jags.samples(jags,c('alpha0','alpha1'),iter,progress.bar="none")
  
  a0<-tt$alpha0[1,,]
  a1<-tt$alpha1[1,,]
  
  
  for (j in 1:length(doses)){
    logit <- a0 + a1 * sdose[j]
    p[,j]<-exp(logit)/(1+exp(logit))
  }
  
  prob.next<-mat.or.vec(length(doses),1)
  for (j in 2:length(doses)){
    y<-p[,j]-p[,1]
    prob.next[j]<-mean(y <=(target.increase+delta) & (y>=target.increase-delta))
    if(mean(y>=(target.increase+2*delta))>overdose){
      prob.next[j]<-0
    }
  }
  
  return(prob.next)
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

# doses<-c(1:doses)
selection<-rep(NA,length(doses))
nextdose<-starting.dose  
n.tox<-rep(0,length(doses))
s.tox<-rep(0,length(doses))
bad.eff<-bad.tox<-bad.both<-bad.doses<-good.doses<-unknown.doses<-c()


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
  prob.next<-blrm.covid(n=n.tox,s=s.tox,sdose=sdose,priorMean=priorMean,priorPrec=priorPrec,overdose=0.25,delta=0.05,iter=10000)
  prevdose<- max(which(n.tox>0))
  if(all(prob.next==0)){

    stop<-1
    break() 
  }else{
    max.safe.dose<-max(which(prob.next>0))
    nextdose<-min(prevdose+1,which.max(prob.next)) #(no skipping doses)
  }
  if(esc_saf==0){
    #if we already have patients on the next dose, we don't go to that dose next
    if(n.tox[nextdose]>0){
      add.toxicity<-0
    }else{
      add.toxicity<-1 #we do go to that dose next
    }
  }else{

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
            #browser()
            control.index<-which((time[[1]]<=(week-3))&(control_vec==dose))
             PPP <- cox_decide_postprob(control=response.eff.all[[1]][control.index,],trt=response.eff.all[[dose]][which(time[[dose]]<=(week-3)),],HR=HRalt,p1=rho0)
          }
          
          #decision on whether to escalate
          #if posterior prob is in between bounds and the dose is safe
          if(PPP > futile & PPP < promise & prob.next[dose]>0 ){
            #if max number of patients reached
            if(length(time[[dose]])>=max.eff.ss){

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

return(c(selection,samples,week))
}