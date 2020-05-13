#### Analyses for Istanbul Faculty of Medicine Corona-center deidentified patient data ####

source(file = "functions.R")
library(epitools)
citation("epitools") 
packageVersion("epitools")
library(e1071)
citation("e1071")
packageVersion("e1071")


#### DATA IMPORT ####
impdata <- read.csv("COVID_ISTANBUL_DATA_12_04_2020.csv") #import dataset
dim(impdata) # Should have 290 rows 12.April.2020
impdata$Solid.malignancy_enc <- ifelse(impdata$Solid.malignancy==0,0,1) #encode solid malignancy as binary
impdata$malignancy_any <- factor(ifelse(impdata$Solid.malignancy_enc==1 | impdata$Hematologic.malignancy == 1, 1, 0)) #encode hemato malignancy as binary
impdata$ICU_yesno <- factor(impdata$ICU_yesno, labels = c("no", "yes")) # ICU_yesno is our main outcome
impdata$Otomatik.gruplama <- factor(impdata$Otomatik.gruplama, labels = c("ICU-free", "ICU")) # NA are patients who are still neither discharged nor ICU 12.04.2020
#### Collapse CT scoring ####
# We originally had subcategories in between main categories. These are not yet standardized that we are working on standardizing
impdata$CT_collapse <- as.character(impdata$CT)
impdata$CT_collapse[impdata$CT_collapse %in% c("non-COVID", "normal")] <- "Normal or non-COVID"
impdata$CT_collapse[impdata$CT_collapse %in% c("very mild", "mild", "mild-moderate")] <- "Mild"
impdata$CT_collapse[impdata$CT_collapse %in% c("moderate", "moderate-diffuse")] <- "Moderate"
impdata$CT_collapse[impdata$CT_collapse %in% c("diffuse", "very diffuse")] <- "Diffuse"
impdata$CT_collapse[impdata$CT_collapse == c("")] <- NA
impdata$CT_collapse <- factor(impdata$CT_collapse)
impdata$CT_collapse <- factor(impdata$CT_collapse, levels(impdata$CT_collapse)[c(4,2,3,1)]) # reorder factors logically
table(impdata$CT_collapse, impdata$Otomatik.gruplama)

# Factoring ACEi's and ARB's
impdata$ARB.or.ACEi <- factor(impdata$ARB.or.ACEi, labels = c("other", "ARB or ACEi"))
impdata$ACEi <- factor(impdata$ACEi, labels = c("other or ARB", "ACEi"))
impdata$ARB <- factor(impdata$ARB, labels = c("other or ACEi", "ARB"))
impdata$RT_PCR <- factor(ifelse(impdata$RT_PCR==1,"positive","not positive"))
impdata$RT_PCR[is.na(impdata$RT_PCR)] <- "not positive"

#### EXCLUDE NON-COVID, DECLINED TO PARTICIPATE AND DROP-OUTS####
# (nobody declined to participate - 12.04.2020)
alldata <- subset(impdata, impdata$Include.Exclude==1)
dim(alldata) # Excluding those ineligible or with missing data 

#### Comparison of RT_PCR and CT grading rates #### 
# How do RT_PCR results compare to CT?
table(alldata$RT_PCR)
table(alldata$CT_collapse, alldata$RT_PCR)
#distribution of positive PCR results are distributed evenly between different CT categories
fisher.test(table(alldata$CT_collapse, alldata$RT_PCR)) 
# However, they are slightly higher in more severe patients. 
# Higher viral load or overdiagnosis?
# Percentage and confidence interval for RT-PCR positivity rates
t <- table(alldata$RT_PCR)
t
binom.test(t[2], t[1]+t[2])
# Percentage and confidence interval for CT scan positivity rates
t <- table(alldata$CT_typical_covid)
t
binom.test(t[2], t[1]+t[2])




#### SUBSET THOSE WHO REACHED THE ENDPOINT as of 12.04.2020 (DISCHARGED OR WENT TO ICU) ####
completedata <- subset(alldata, !is.na(alldata$Otomatik.gruplama))
dim(completedata)
disc_or_icu <- subset(completedata, completedata$ICU_yesno == "yes" | completedata$Discharged == 1)
dim(disc_or_icu) 


#### ANALYSES ON ICU-free vs ICU, all participants ####
# Total patients ICU-free vs ICU
table(disc_or_icu$ICU_yesno)

# Distribution of ICU-free vs ICU shown on Figure 1A
table(disc_or_icu$RT_PCR, disc_or_icu$ICU_yesno)
table(disc_or_icu$CT_collapse, disc_or_icu$ICU_yesno) # Infiltration level on CT

# Days spent in the hospital as of 12.04.2020 
myquart("Days.hospitalized", varq = "ICU_yesno", mydata = disc_or_icu) # ICU-free vs ICU, median (IQR)

# Distribution of positive RT-PCR tests between ICU-free and ICU groups (should be similar)
table(alldata$RT_PCR, alldata$ICU_yesno)
fisher.test(table(alldata$RT_PCR, alldata$ICU_yesno))

# Distribution of positive CT results between ICU-free and ICU groups (they are not similar)
table(disc_or_icu$CT_collapse, disc_or_icu$ICU_yesno)
fisher.test(table(disc_or_icu$CT_collapse, disc_or_icu$ICU_yesno))
epitab(disc_or_icu$CT_collapse, disc_or_icu$ICU_yesno)

# Age distribution
# Plots of age distributions
prism_plot("age", "ICU_yesno", title_add = "Age vs Intensive Care Requirement", mydata = disc_or_icu)
densitycurves("age", data = disc_or_icu, groups = "ICU_yesno", delete_title = F,
              start_from_zero = T, bin_count = 18, breaks_custom = seq(0,100,20))
# Test for normal distribution of Age
shapiro.test(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="no")])
skewness(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="no")])
shapiro.test(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="yes")])
skewness(disc_or_icu$age[which(disc_or_icu$ICU_yesno=="yes")])


# Odds ratios of presenting symptoms
epitab(table(disc_or_icu$fever.reported, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$coughing, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$dyspnea, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$fatigue.or.myalgia, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$nausea, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$diarrhea, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$anosmia, disc_or_icu$ICU_yesno))

#### Odds Ratio of Labored breathing as a finding 
# (I apologize dearly for the inconvenient parameter coding 
#("d"yspnea for shortness of breath, "D"yspnea for labored breathing))
table(disc_or_icu$Dyspnea)
epitab(table(disc_or_icu$Dyspnea, disc_or_icu$ICU_yesno))

#### Odds Ratio of Respiratory Rate >= 22
table(disc_or_icu$Respiratory.Rate>=22, disc_or_icu$ICU_yesno)
epitab(table(disc_or_icu$Respiratory.Rate>=22, disc_or_icu$ICU_yesno)) # I tuned this manually to maximize the odds ratio
prism_plot("Respiratory.Rate", varx = "ICU_yesno") #### Plotting the Respiratory Rate
densitycurves("Respiratory.Rate", data = disc_or_icu)

#### Odds Ratios of Comorbidities:
epitab(table(disc_or_icu$Hypertension.history, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$DM.History, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$COPD.asthma.history, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Coronary.Artery.Disease, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Congestive.Heart.Failure, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Solid.malignancy_enc, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$Hematologic.malignancy, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$malignancy_any, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$smoking, disc_or_icu$ICU_yesno))
epitab(table(disc_or_icu$alcohol, disc_or_icu$ICU_yesno))

# 95%CI of the percentage of Hypertension history
binom.test(table(disc_or_icu$Hypertension.history)[c(2,1)])

#### 95% Confidence Interval calculations for the frequency of
#### travel history, contact history and being a heathcare worker among all admitted
#### patients
binom.test(table(disc_or_icu$Travel.history)[c(2,1)])
binom.test(table(disc_or_icu$Contact.history)[c(2,1)])
binom.test(table(disc_or_icu$Healthcare)[c(2,1)])


################################################################################################
################################ PRIMARY ANALYSES ##############################################
################################################################################################
#### HYPERTENSION AND RAAS ####
# Creating subset of participants with hypertension only 
dim(disc_or_icu)
ht <- subset(disc_or_icu, disc_or_icu$Hypertension.history==1)
dim(ht)

#### ODDS RATIOS OF RAAS USE ####
table(ht$ARB.or.ACEi, ht$ICU_yesno)
epitab(ht$ARB.or.ACEi, ht$ICU_yesno)
table(ht$ACEi, ht$ICU_yesno)
epitab(ht$ACEi, ht$ICU_yesno)
table(ht$ARB, ht$ICU_yesno)
epitab(ht$ARB, ht$ICU_yesno)

#### ARE THE AGES DIFFERENT? ####
prism_plot("age", "ARB.or.ACEi", mydata = ht) # Colors show ICU-free or ICU
prism_plot("age", "ARB", mydata = ht) # Colors show ICU-free or ICU
prism_plot("age", "ACEi", mydata = ht) # Colors show ICU-free or ICU
prism_plot("age", "ARB_ACEi_OTHER", mydata = ht) # Colors show ICU-free or ICU

wilcox.test(age~ARB.or.ACEi, data = ht)

wilcox.test(age~ACEi, data = ht)
t.test(age~ACEi, data = ht)
lm(ICU_yesno~ACEi+age, data = ht)

wilcox.test(age~ARB, data = ht)

