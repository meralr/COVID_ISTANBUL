# Analyses

# Ran using R v4.0.0 - - "Arbor Day"

#### Analyses for Istanbul Faculty of Medicine Corona-center deidentified patient data ####
library(psych)
citation("psych") 
packageVersion("psych") #‘1.9.12.31’

library(epitools)
citation("epitools") 
packageVersion("epitools") # ‘0.5.10.1’

library(e1071)
citation("e1071")
packageVersion("e1071") # ‘1.7.3’

source(file = "functions.R") # Handy home brewed functions for plotting data and building tables
source(file = "data_import.R") # Imports the data

suppressWarnings(dir.create("table_output"))

# imported data set contains all cases that were admitted to the istanbul faculty of medicine
# corona center at any point
dim(impdata)


# CREATING A SUBSET OF ONLY THOSE WHO REACHED AN ENDPOINT
# End-points are:
# # severe     = ICU admission, hospitalization >= 14 days or death
# # not severe = Discharge event-free
### DEFINING THE OUTCOME VARIABLE ###
impdata$outcome <- as.factor(ifelse(impdata$Ex=="Excitus", "severe",
                                    ifelse(impdata$iCU_yesno=="yes", "severe",
                                           ifelse(impdata$Days_hospitalized_on_day_x >= 14, "severe",
                                                  ifelse(impdata$Taburcu==1, "not severe", "did not reach endpoint")))))
impdata$outcome_category <- as.factor(ifelse(impdata$Ex=="Excitus", "Death",
                                             ifelse(impdata$iCU_yesno=="yes", "ICU", 
                                                    ifelse(impdata$Days_hospitalized_on_day_x >= 14, "Long hospitalization",
                                                           ifelse(impdata$Taburcu==1, "not severe", "did not reach endpoint")))))
table(impdata$outcome)
table(impdata$outcome_category)


# Drop cases admitted before 2020-05-11
impdata <- subset(impdata, impdata$Date_admitted<="2020-05-11")
print(dim(impdata))
table(impdata$outcome)


# Drop cases with neither a positive PCR nor Radiologic diagnosis:
table(impdata$CT_arbitrary_score)
impdata <- subset(impdata, impdata$RT_PCR=="positive" | impdata$CT_arbitrary_score != "Normal or non-COVID")
dim(impdata)
table(impdata$Gender_1male_2female)
summary(impdata$age)

# PCR positivity rate
table(impdata$RT_PCR)
prop.table(table(impdata$RT_PCR))
# Is there a link between PCR positivity and sick days before hospitalization?
# (we would expect less positivity in late arrivals per http://jamanetwork.com/article.aspx?doi=10.1001/jama.2020.8259)
prism_plot("Days_from_first_symptom_to_hospitalization", "RT_PCR", mydata = impdata)
wilcox.test(Days_from_first_symptom_to_hospitalization~RT_PCR, data = impdata) # If anything, positives presented later!

#CT positivity
table(impdata$CT_category)
table(impdata$CT_category!="Normal" & !is.na(impdata$CT_category));prop.table(table(impdata$CT_category!="Normal" & !is.na(impdata$CT_category)))
table(impdata$CT_arbitrary_score, impdata$RT_PCR); prop.table(table(impdata$CT_arbitrary_score, impdata$RT_PCR))
table(impdata$CT_category, impdata$RT_PCR)


# Patients early in their disease course are dropped from the analyses
# THESE LINES OF CODE ARE NOW REDUNDANT AS ALL PATIENTS REACHED AN OUTCOME
# (ie. all patients admitted on or before 11th of May were either discharged or classified as having severe disease on 22nd of May)
covdata <- subset(impdata, impdata$outcome != c("did not reach endpoint"))
dim(covdata)
covdata$outcome <- factor(covdata$outcome) #reset factor levels
covdata$outcome_category <- factor(covdata$outcome_category)

#Outcome breakdown
table(covdata$outcome); prop.table(table(covdata$outcome))
table(covdata$outcome_category)
table(covdata$Days_hospitalized_on_day_x>=14, covdata$outcome_category)
table(covdata$Days_hospitalized_on_day_x>=14, covdata$Ex)
table(covdata$iCU_yesno, covdata$outcome)
table(covdata$iCU_yesno, covdata$outcome_category)
table(covdata$Ex, covdata$outcome_category)
table(covdata$iCU_yesno, covdata$Ex) # 9 patients died outside of the ICU
summary(subset(covdata, covdata$iCU_yesno == "no" & covdata$Ex == "Excitus")["age"])
#days spent in the hospital on May 22nd
summary(covdata$Days_hospitalized_on_day_x)
describeBy(covdata$Days_hospitalized_on_day_x, group = covdata$outcome)

# Overall case fatality rate among hospitalized patients at Istanbul Faculty of Medicine, Corona-center
table(covdata$Ex)
prop.table(table(covdata$Ex)) # 8.7% is subject to increase somewhat as some patients are still intubated and are in critical condition 
binom.test(table(covdata$Ex)[c(2,1)])

# age difference
describeBy(covdata$age, covdata$outcome)
describeBy(covdata$age, covdata$iCU_yesno)
describeBy(covdata$age, covdata$Ex) # As young as 33 died and as old as 98 were saved
wilcox.test(age~outcome, data = covdata)
prism_plot("age", "outcome", mydata=covdata)
prism_plot("age", "iCU_yesno", mydata=covdata)
prism_plot("age", "Ex", mydata=covdata)

# Data needed for our "population pyramid" style graph on figure 1
# Graphpad prism v8.3 was used
# Method can be found here: https://www.youtube.com/watch?v=YonH4LiLEdU
age_stack <- table(cut(covdata$age, breaks = seq(0, 100, 5)), covdata$outcome_category)
write.csv(age_stack, "table_output/age_groups.csv")
# Powerpoint file in the misc folder contains plots that will open in prism upon double-clicking
# (A copy of Prism must be installed on your computer)


# Contact history
binom.test(table(covdata$Contact_history)[c(2,1)]) # Only 34% of all patients were aware of any sick contacts


##### ODDS RATIOS #####


# Level of infiltration on CT highly predicts severe disease
epitab(covdata$CT_category_for_odds, covdata$outcome)

# Male gender 
epitab(covdata$Gender_1male_2female, covdata$outcome, rev = "rows")

# Presenting symptoms
epitab(table(covdata$fever_reported, covdata$outcome))
epitab(table(covdata$coughing, covdata$outcome))
epitab(table(covdata$fatigue_or_myalgia, covdata$outcome)) 
    # Patients reporting fatigue or myalgia appear less likely to have severe disease.
    # Perhaps they ignore these because their clinical picture is dominated by dyspnea.
epitab(table(covdata$dyspnea, covdata$outcome))
epitab(table(covdata$Dyspnea, covdata$outcome)) # "D"yspnea encoded with capital D is labored breathing
epitab(table(covdata$Respiratory_Rate>=22, covdata$outcome)) # I tuned this manually to maximize the odds ratio
epitab(table(covdata$nausea, covdata$outcome))
epitab(table(covdata$diarrhea, covdata$outcome))
epitab(table(covdata$anosmia, covdata$outcome))

#### Odds Ratios of Comorbidities:
table(covdata$Hypertension_history);epitab(table(covdata$Hypertension_history, covdata$outcome))
prop.table(table(covdata$Hypertension_history)) # 41% of hospitalized COVID-19 patients had HT
epitab(table(covdata$DM_History, covdata$outcome))
epitab(table(covdata$COPD_asthma_history, covdata$outcome))
epitab(table(covdata$Coronary_Artery_Disease, covdata$outcome))
epitab(table(covdata$Congestive_Heart_Failure, covdata$outcome))
epitab(table(covdata$Solid_malignancy_enc, covdata$outcome))
epitab(table(covdata$Hematologic_malignancy, covdata$outcome))
epitab(table(covdata$malignancy_any, covdata$outcome))

# SMOKING
epitab(table(covdata$smoking, covdata$outcome))
prop.table(table(covdata$smoking))
table(covdata$smoking, covdata$Gender_1male_2female)
prop.table(table(covdata$smoking, covdata$Gender_1male_2female)[,1])
prop.table(table(covdata$smoking, covdata$Gender_1male_2female)[,2])
covdata$COPD_asthma_history
smokey <- matchControls(smoking~age+
                          Gender_1male_2female+
                          Hypertension_history+
                          DM_History+
                          COPD_asthma_history+
                          Coronary_Artery_Disease+
                          Congestive_Heart_Failure,
                        data = covdata, caselabel = "yes", contlabel = "no")
s <- covdata[c(smokey$cases,smokey$controls),]
dim(s)
addmargins(table(s$smoking, s$outcome))
epitab(s$smoking, s$outcome)


#Hypertension percentage
binom.test(table(covdata$Hypertension_history)[c(2,1)])
epitab(covdata$Hypertension_history, covdata$outcome)
# Create a subset of patients with Hypertension
ht <- subset(covdata, covdata$Hypertension_history == "yes")
print(dim(ht))
table(ht$outcome); prop.table(table(ht$outcome))
table(ht$Ex); prop.table(table(ht$Ex))
table(covdata$Hypertension_history, covdata$Ex); fisher.test(table(covdata$Hypertension_history, covdata$Ex))
prism_plot("Systolic", "outcome", mydata = ht); describeBy(x = ht$Systolic, group = ht$outcome)
wilcox.test(Systolic~outcome, data = ht); myquart("Systolic", "outcome", mydata = ht)
prism_plot("Diastolic", "outcome", mydata = ht); describeBy(x = ht$Diastolic, group = ht$outcome)
wilcox.test(Diastolic~outcome, data = ht); myquart("Diastolic", "outcome", mydata = ht)
prism_plot("Pulse", "outcome", mydata = ht); describeBy(x = ht$Pulse, group = ht$outcome)
wilcox.test(Pulse~outcome, data = ht); myquart("Pulse", "outcome", mydata = ht)


# Crude odds ratios of exposre to each RAAS inhibitor class with respect to other classes as controls
table(ht$exposure)
table(ht$exposure, ht$outcome)
epitab(ht$exposure, ht$outcome)
# Crude Odds ratio of exposure to any RAAS inhibitor drug vs outcome 
epitab(ht$RAAS, ht$outcome)


# Drop 
# - Patients using more than 3 anti-hypertensives
ht3 <- subset(ht, ht$number_of_anti_ht_yes_tzd<=3)
print(dim(ht3))

##### MATCHING #####
# 1:1 Match an ARB to each ACEi
# then 1:1 match 1 "Other" to each ACEi
ht3$acei_vs_arb <- factor(ifelse(ht3$exposure=="ACEi", "case", 
                                  ifelse(ht3$exposure=="ARB", "control", NA)))
ht3$acei_vs_control <- factor(ifelse(ht3$exposure=="ACEi", "case", 
                                ifelse(ht3$exposure=="Other", "control", NA)))
set.seed(100)
match_an_ARB <- matchControls(acei_vs_arb~age+
                                Gender_1male_2female+
                                Days_from_first_symptom_to_hospitalization+
                                smoking+
                                DM_History+
                                COPD_asthma_history+
                                Coronary_Artery_Disease+
                                Congestive_Heart_Failure+
                                number_of_anti_ht_yes_tzd+
                                Kre, 
                              data = ht3, caselabel = "case", contlabel = "control")
match_a_control <- matchControls(acei_vs_control~age+
                                   Gender_1male_2female+
                                   Days_from_first_symptom_to_hospitalization+
                                   smoking+
                                   DM_History+
                                   COPD_asthma_history+
                                   Coronary_Artery_Disease+
                                   Congestive_Heart_Failure+
                                   number_of_anti_ht_yes_tzd+
                                   Kre, 
                                 data = ht3, caselabel = "case", contlabel = "control")
m <- ht3[c(match_an_ARB$cases,match_an_ARB$controls,match_a_control$controls),]
m$exposure <- factor(m$exposure) #resetting factor to get rid of blanks
dim(m)
print(table(m$exposure))
addmargins(table(m$exposure, m$outcome))

##### HOW DID THE MATCHING GO? #####
# Checking for distributions before and after matching:

densitycurves("age", data=ht, groups = "exposure")
densitycurves("age", data=m, groups = "exposure") # Looks well matched
prism_plot("age", "exposure", mydata = ht)
print(prism_plot("age", "exposure", mydata = m)) # Looks well matched
kruskal.test(m$age, m$exposure)
summary(aov(age~exposure, m)) # Looks well matched

table(ht$exposure, ht$Gender_1male_2female)
table(m$exposure, m$Gender_1male_2female) # Looks evenly distributed
fisher.test(table(m$exposure, m$Gender_1male_2female))

table(m$Loop_diuretic, m$exposure); fisher.test(table(m$Loop_diuretic, m$exposure)) # Does not look evenly distributed
table(m$Beta_Blocker_enc, m$exposure); fisher.test(table(m$Beta_Blocker_enc, m$exposure)) 

densitycurves("Days_from_first_symptom_to_hospitalization", data=ht, groups = "exposure")
densitycurves("Days_from_first_symptom_to_hospitalization", data=m, groups = "exposure") # Looks well matched
prism_plot("Days_from_first_symptom_to_hospitalization", "exposure", mydata = ht)
prism_plot("Days_from_first_symptom_to_hospitalization", "exposure", mydata = m) # Medians are different
kruskal.test(Days_from_first_symptom_to_hospitalization~exposure, data = m) # Ranks are not statistically different
summary(aov(Days_from_first_symptom_to_hospitalization~exposure, data = m)) # Parametric test also doesn't detect difference

table(ht$exposure, ht$smoking)
table(m$exposure, m$smoking) # There appear to be extremely few smokers in the cohort, but smoking status reportedly was not properly questioned and is unreliable
table(ht$COPD_asthma_history)
table(ht$COPD_asthma_history, ht$smoking)

table(ht$exposure, ht$DM_History)
table(m$exposure, m$DM_History) # Looks evenly distributed

table(ht$exposure, ht$COPD_asthma_history)
table(m$exposure, m$COPD_asthma_history) # Looks evenly distributed

table(ht$exposure, ht$Coronary_Artery_Disease)
table(m$exposure, m$Coronary_Artery_Disease); fisher.test(table(m$exposure, m$Coronary_Artery_Disease)) # Looks OK

table(ht$exposure, ht$Congestive_Heart_Failure)
table(m$exposure, m$Congestive_Heart_Failure) # Looks evenly distributed

table(ht$exposure, ht$number_of_anti_ht_yes_tzd)
table(m$exposure, m$number_of_anti_ht_yes_tzd) # Bad match visually
fisher.test(table(m$exposure, m$number_of_anti_ht_yes_tzd)) # Bad match statistically
table(ht$exposure, ht$Thiazide_enc)
table(m$exposure, m$Thiazide_enc) # Bad match appears to be driven by combination preparations that included thiazides
                                  # ARB and ACEi preparations are commonly combined with thiazides.
# Does thiazide exposure predict severe disease?
epitab(m[m$exposure%in%c("ARB", "ACEi"),]$Thiazide_enc, m[m$exposure%in%c("ARB", "ACEi"),]$outcome) # No

densitycurves("Kre", data=ht, groups = "exposure")
densitycurves("Kre", data=m, groups = "exposure")
prism_plot("Kre", "exposure", mydata = ht)
prism_plot("Kre", "exposure", mydata = m) # Unable to tell visually
kruskal.test(Kre~exposure, data = m) # Good match statistically
summary(aov(Kre~exposure, m)) # Parametric test also not bad

##### PRIMARY ANALYSES #####
# Unmatched
print(addmargins(table(ht$exposure, ht$outcome)))
write.csv(addmargins(table(ht$exposure, ht$outcome)), "table_output/unmatched.csv")
print(epitab(ht$exposure, ht$outcome))
# Matched
print(addmargins(table(m$exposure, m$outcome)))
print(fisher.test(table(m$exposure, m$outcome)))
write.csv(addmargins(table(m$exposure, m$outcome)), "table_output/matched.csv")
print(epitab(m$exposure, m$outcome))





