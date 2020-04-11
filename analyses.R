#tertemiz

source("source.R")

#### DATA IMPORT ####
impdata <- read.csv("data/COVID_DATA_v13.csv") #import dataset
impdata$Solid.malignancy_enc <- ifelse(impdata$Solid.malignancy==0,0,1)
impdata$malignancy_any <- factor(ifelse(impdata$Solid.malignancy_enc==1 | impdata$Hematologic.malignancy == 1, 1, 0))
impdata$ICU_yesno <- factor(impdata$ICU_yesno, labels = c("no", "yes"))
impdata$Otomatik.gruplama <- factor(impdata$Otomatik.gruplama, labels = c("ICU-free", "ICU"))
impdata$ARB.or.ACEi <- factor(impdata$ARB.or.ACEi, labels = c("other", "ARB or ACEi"))
impdata$ACEi <- factor(impdata$ACEi, labels = c("other or ARB", "ACEi"))
impdata$ARB <- factor(impdata$ARB, labels = c("other or ACEi", "ARB"))
impdata$RT_PCR <- factor(ifelse(impdata$RT_PCR==1,"positive","negative"))
impdata$RT_PCR[is.na(impdata$RT_PCR)] <- "negative"
dim(impdata)

#### EXCLUDE NON-COVID AND DROP-OUTS####
alldata <- subset(impdata, impdata$Makalede.olsun==1)
dim(alldata)

table(alldata$RT_PCR)
table(alldata$CT_typical_covid, alldata$RT_PCR)

#### SUBSET THOSE WHO FINISHED (DISCHARGED OR WENT TO ICU) ####
completedata <- subset(alldata, !is.na(alldata$Otomatik.gruplama))
dim(completedata)
disc_or_icu <- subset(completedata, completedata$ICU_yesno == "yes" | completedata$Taburcu == 1)
dim(disc_or_icu) #Should have the same dimensions as completedata, or there is something wrong

table(disc_or_icu$ICU_yesno)
table(disc_or_icu$ICU_yesno, disc_or_icu$RT_PCR)
table(disc_or_icu$CT_typical_covid, disc_or_icu$RT_PCR)
table(disc_or_icu$CT_typical_covid, disc_or_icu$Otomatik.gruplama)
fisher.test(table(disc_or_icu$CT_typical_covid, disc_or_icu$Otomatik.gruplama))
binom.test(table(disc_or_icu$CT_typical_covid, disc_or_icu$Otomatik.gruplama)[c(2,1),1])




######## REAL DEAL ############
#Abstract
table(alldata$Gender_1male_2female)

disc_or_icu$Days.hospitalized
myquart("Days.hospitalized", varq = "ICU_yesno", mydata = disc_or_icu)
range(disc_or_icu$Days.hospitalized)

# Percentage and confidence interval for RT-PCR positivity rates
t <- table(alldata$RT_PCR)
t
binom.test(t[2], t[1]+t[2])

# Percentage and confidence interval for CT scan positivity rates
t <- table(alldata$CT_typical_covid)
t
binom.test(t[2], t[1]+t[2])

# Distribution of positive RT-PCR tests between ICU-free and ICU groups (should be similar)
table(alldata$RT_PCR, alldata$ICU_yesno)
t <- table(alldata$RT_PCR, alldata$Otomatik.gruplama)
t
fisher.test(t)

# Plots of age distributions (make new_device=F if working on a non-Windows machine)
prism_plot("age", "ICU_yesno", title_add = "Age vs Intensive Care Requirement", mydata = disc_or_icu, new_device = F)
densitycurves("age", data = disc_or_icu, groups = "ICU_yesno", new_device = F, delete_title = F,
              start_from_zero = T, bin_count = 18, breaks_custom = seq(0,100,20))
# Test for normality and calculate skewness of the age distribution, separately for ICU-free and ICU groups
library(e1071)
shapiro.test(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="no")])
skewness(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="no")])
shapiro.test(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="yes")])
skewness(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="yes")])

#### Odds Ratio of dyspnea as a symptom
library(epitools)
citation("epitools")
packageVersion("epitools")

table(disc_or_icu$dyspnea)
epitab(table(disc_or_icu$dyspnea, disc_or_icu$ICU_yesno))

#### Odds Ratio of Dyspnea as a finding (sorry for the inconvenient parameter coding)
table(disc_or_icu$Dyspnea)
epitab(table(disc_or_icu$Dyspnea, disc_or_icu$ICU_yesno))

#### Odds Ratio of Respiratory Rate > 22
table(disc_or_icu$Respiratory.Rate>22, disc_or_icu$ICU_yesno)
epitab(table(disc_or_icu$Respiratory.Rate>22, disc_or_icu$ICU_yesno))

#### Plotting the Respiratory Rate
prism_plot("Respiratory.Rate", varx = "ICU_yesno", new_device =F)
densitycurves("Respiratory.Rate", data = disc_or_icu, new_device = F)

#### 95% Confidence Interval calculations for the frequency of
#### travel history, contact history and being a heathcare worker among all admitted
#### patients
binom.test(table(disc_or_icu$Travel.history)[c(2,1)])
binom.test(table(disc_or_icu$Contact.history)[c(2,1)])
binom.test(table(disc_or_icu$Healthcare)[c(2,1)])

#### Odds Ratios of Hypertension, Diabetes, COPD/Asthma, CAH and CHF in the history:
binom.test(table(disc_or_icu$Hypertension.history)[c(2,1)])
epitab(table(disc_or_icu$Hypertension.history, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$DM.History, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$COPD.asthma.history, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Coronary.Artery.Disease, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Congestive.Heart.Failure, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Congestive.Heart.Failure, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Solid.malignancy_enc, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Hematologic.malignancy, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$malignancy_any, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$smoking, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$alcohol, disc_or_icu$ICU_yesno))


#### HYPERTENSION AND RAAS ####
dim(disc_or_icu)
ht <- subset(disc_or_icu, disc_or_icu$Hypertension.history==1)
dim(ht)

ht$ARB.or.ACEi
prism_plot("age", "ARB.or.ACEi", mydata = ht)
prism_plot("age", "ARB", mydata = ht)
prism_plot("age", "ACEi", mydata = ht)

#### ODDS RATIOS ####
table(ht$ARB.or.ACEi, ht$ICU_yesno)
table(ht$ARB.or.ACEi)
fisher.test(table(ht$ARB.or.ACEi, ht$ICU_yesno))
epitab(ht$ARB.or.ACEi, ht$ICU_yesno)

table(ht$ACEi, ht$ICU_yesno)
table(ht$ACEi)
fisher.test(table(ht$ACEi, ht$ICU_yesno))
epitab(ht$ACEi, ht$ICU_yesno)

table(ht$ARB, ht$ICU_yesno)
table(ht$ARB)
fisher.test(table(ht$ARB, ht$ICU_yesno))
epitab(ht$ARB, ht$ICU_yesno)



#### ARE THE AGES DIFFERENT? ####
prism_plot("age", "ARB.or.ACEi", mydata = ht)
wilcox.test(age~ARB.or.ACEi, data = ht)

prism_plot("age", "ACEi", mydata = ht)
wilcox.test(age~ACEi, data = ht)
k<-lm(ICU_yesno~ACEi+age, data = ht)
k<-lm(D.dimer~age, data = ht)
summary(k)

prism_plot("age", "ARB", mydata = ht)
wilcox.test(age~ARB, data = ht)

plot(disc_or_icu$D.dimer, disc_or_icu$Lymp)


table1 <- function(indep = "ICU_yesno", mydata = disc_or_icu){
  mytbl <- data.frame("Variable"="varx", "group1"="group1", "group2"="group2", "pvalue"="pvalue", "marker"="marker",
                      stringsAsFactors = F)
  gt <- table(mydata$Gender_1male_2female, mydata$ICU_yesno) #for gender
  gt1 <- paste0(gt[2,1],"F:",gt[1,1],"M") # for gender
  gt2 <- paste0(gt[2,2],"F:",gt[1,2],"M") # for gender
  fgt <- fisher.test(gt)
  mytbl<- rbind(
    mytbl, 
    # disc_or_icu$Gender_1male_2female
    c("N", as.character(table(mydata[indep])),""),
    probs2("RT_PCR", indep, mydata = mydata), 
    c("CT_typical_covid", length(mydata$CT_typical_covid), length(mydata$CT_typical_covid), length(mydata$CT_typical_covid)),
    myquart("age", indep, mydata = mydata, k=0),
    c("Gender", gt1, gt2, round(fgt$p.value, digits=4), rownames(gt)[2]),
    myquart("BMI", indep, mydata = mydata, k=1),
    # probs2("Gender_1male_2female", indep, mydata = mydata),
    probs2("Travel.history", indep, mydata = mydata),
    probs2("Contact.history", indep, mydata = mydata),
    probs2("Healthcare", indep, mydata = mydata),
    
    probs2("fever.reported", indep, mydata = mydata),
    probs2("coughing", indep, mydata = mydata),
    probs2("sputum", indep, mydata = mydata),
    probs2("dyspnea", indep, mydata = mydata),
    probs2("fatigue.or.myalgia", indep, mydata = mydata),
    probs2("nausea", indep, mydata = mydata),
    probs2("diarrhea", indep, mydata = mydata),
    probs2("anosmia", indep, mydata = mydata),
    myquart("fever.days", indep, mydata = mydata, k=0),
    myquart("coughing.days", indep, mydata = mydata, k=0),
    myquart("sputum.days", indep, mydata = mydata, k=0),
    myquart("dyspnea.days", indep, mydata = mydata, k=0),
    myquart("fatigue.or.myalgia.days", indep, mydata = mydata, k=0),
    myquart("Days.from.first.symptom.to.hospitalization", indep, mydata = mydata, k=0),
    myquart("Days.hospitalized", indep, mydata = mydata, k=0),
    myquart("Days.since.first.symptom", indep, mydata = mydata, k=0),
    
    probs2("Hypertension.history", indep, mydata = mydata),
    # probs2("ARB.or.ACEi", indep, mydata = subset(mydata, mydata["Hypertension.history"]==1)),
    # probs2("ARB", indep, mydata = subset(mydata, mydata["Hypertension.history"]==1)),
    # probs2("ACEi", indep, mydata = subset(mydata, mydata["Hypertension.history"]==1)),
    probs2("ARB.or.ACEi", indep, mydata = mydata),
    probs2("ARB", indep, mydata = mydata),
    probs2("ACEi", indep, mydata = mydata),
    probs2("DM.History", indep, mydata = mydata),
    probs2("COPD.asthma.history", indep, mydata = mydata),
    probs2("Coronary.Artery.Disease", indep, mydata = mydata),
    probs2("Congestive.Heart.Failure", indep, mydata = mydata),
    probs2("Solid.malignancy_enc", indep, mydata = mydata),
    probs2("Hematologic.malignancy", indep, mydata = mydata),
    probs2("smoking", indep, mydata = mydata),
    probs2("alcohol", indep, mydata = mydata),
    
    myquart("Body.temp", indep, mydata = mydata, k=1),
    myquart("spo2", indep, mydata = mydata, k=0),
    myquart("Systolic", indep, mydata = mydata, k=0),
    myquart("Diastolic", indep, mydata = mydata, k=0),
    myquart("Pulse", indep, mydata = mydata, k=0),
    myquart("Respiratory.Rate", indep, mydata = mydata, k=0),
    probs2("Dyspnea", indep, mydata = mydata, k=0),
    
    myquart("pH", indep, mydata = mydata, k=2),
    myquart("pO2", indep, mydata = mydata, k=0),
    myquart("pCo2", indep, mydata = mydata, k=0),
    myquart("HCO3", indep, mydata = mydata, k=0),
    myquart("Lactate", indep, mydata = mydata, k=1),
    myquart("Hgb", indep, mydata = mydata, k=1),
    myquart("Plt", indep, mydata = mydata, k=0),
    myquart("WBC", indep, mydata = mydata, k=0),
    myquart("Neut", indep, mydata = mydata, k=0),
    myquart("Lymp", indep, mydata = mydata, k=0),
    myquart("Mon", indep, mydata = mydata, k=0),
    # myquart("Eos", indep, mydata = mydata, k=0),
    myquart("BUN", indep, mydata = mydata, k=0),
    myquart("Kre", indep, mydata = mydata, k=1),
    myquart("Na", indep, mydata = mydata, k=0),
    myquart("Cl", indep, mydata = mydata, k=0),
    myquart("K", indep, mydata = mydata, k=1),
    myquart("Glu", indep, mydata = mydata, k=0),
    myquart("AST", indep, mydata = mydata, k=0),
    myquart("ALT", indep, mydata = mydata, k=0),
    myquart("GGT", indep, mydata = mydata, k=0),
    myquart("ALP", indep, mydata = mydata, k=0),
    myquart("LDH", indep, mydata = mydata, k=0),
    myquart("T.prot", indep, mydata = mydata, k=1),
    myquart("Alb", indep, mydata = mydata, k=1),
    myquart("CRP", indep, mydata = mydata, k=0),
    myquart("Prokalsitonin", indep, mydata = mydata, k=2),
    myquart("Ferritin", indep, mydata = mydata, k=0),
    # myquart("Trig", indep, mydata = mydata, k=0),
    myquart("D.dimer", indep, mydata = mydata, k=0),
    myquart("Trop", indep, mydata = mydata, k=1),
    myquart("pro.BNP", indep, mydata = mydata, k=0),
    myquart("Fibrinojen", indep, mydata = mydata, k=0),
    myquart("INR", indep, mydata = mydata, k=1),
    myquart("APTT", indep, mydata = mydata, k=0)
  )
  
  return(mytbl)
}
print(table1())
t <- table1()
write.csv(t, file = "t.csv")

# h <- table1(indep = "ARB.or.ACEi", mydata = ht)
# write.csv(h, file = "h.csv")

by(disc_or_icu$age, disc_or_icu$ICU_yesno, range)
logit <- glm(ICU_yesno ~ age,data=disc_or_icu,family="binomial")
summary(logit)

predict(logit,data.frame(age=40))


disc_or_icu$Days.from.first.symptom.to.ICU

disc_or_icu$Days.from.first.symptom.to.ICU

disc_or_icu$Solid.malignancy
table(disc_or_icu$Solid.malignancy)




