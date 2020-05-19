# Analyses

#### Analyses for Istanbul Faculty of Medicine Corona-center deidentified patient data ####

source(file = "functions.R")
source(file = "data_import.R")

# imported data set contains all cases that were admitted to the istanbul faculty of medicine
# corona center at any point
dim(impdata)

### DEFINING THE OUTCOME VARIABLE ###
impdata$outcome <- as.factor(ifelse(impdata$iCU_yesno=="yes", "severe", 
                                    ifelse(impdata$Ex=="Excitus", "severe",
                                           ifelse(impdata$Days_hospitalized_on_day_x > 14, "severe",
                                                  ifelse(impdata$Taburcu==1, "not severe", "did not reach endpoint")))))

# Drop cases with neither a positive PCR nor Radiologic diagnosis:
impdata <- subset(impdata, impdata$RT_PCR=="positive" | impdata$CT_arbitrary_score != "Normal or non-COVID")
dim(impdata)
table(impdata$outcome)

# Drop cases until 14 days ago
impdata <- subset(impdata, impdata$Date_admitted<="2020-04-28")
dim(impdata)

# CREATING A SUBSET OF THOSE WHO REACHED AN ENDPOINT
# Patients early in their disease course are dropped from the analyses
covdata <- subset(impdata, impdata$outcome != c("did not reach endpoint"))
dim(covdata)
covdata$outcome <- factor(covdata$outcome) #reset factor levels

# Overall case fatality rate among hospitalized patients at Istanbul Faculty of Medicine, Corona-center
table(covdata$Ex)
prop.table(table(covdata$Ex))
binom.test(table(covdata$Ex)[c(2,1)])

# PCR positivity rate
table(covdata$RT_PCR)
table(covdata$RT_PCR, covdata$outcome)
fisher.test(table(covdata$RT_PCR, covdata$outcome)) # distribution of (+) swabs had a tendency to be different among groups
prop.table(table(covdata$RT_PCR, covdata$outcome)[,1]) # 58% of non-severe cases had (+) swabs
prop.table(table(covdata$RT_PCR, covdata$outcome)[,2]) # whereas 69% of severe cases had (+) swabs

# Distribution of CT findings
table(covdata$CT_arbitrary_score, covdata$RT_PCR)
table(covdata$CT_category, covdata$RT_PCR)
table(covdata$CT_arbitrary_score, covdata$outcome)
table(covdata$CT_category, covdata$outcome)

# Contact history
binom.test(table(covdata$Contact_history)[c(2,1)])

#Hypertension percentage
binom.test(table(covdata$Hypertension_history)[c(2,1)])

# Hypertension subset
ht <- subset(covdata, covdata$Hypertension_history == "yes")
dim(ht)

# Odds ratios
library(epitools)
citation("epitools") 
packageVersion("epitools")




# Odds ratios of presenting symptoms
epitab(table(covdata$fever_reported, covdata$outcome))
epitab(table(covdata$coughing, covdata$outcome))
epitab(table(covdata$dyspnea, covdata$outcome))
epitab(table(covdata$fatigue_or_myalgia, covdata$outcome))
epitab(table(covdata$nausea, covdata$outcome))
epitab(table(covdata$diarrhea, covdata$outcome))
epitab(table(covdata$anosmia, covdata$outcome))


epitab(ht$RAAS, ht$outcome)
epitab(ht$ACEi, ht$outcome)
epitab(ht$ARB, ht$outcome)



table(ht$ARB_or_ACEi)
table(ht$ARB_or_ACEi, ht$outcome)


dim(ht)
ht$exposure_unmatched <- ht$ARB_or_ACEi
levels(ht$exposure_unmatched)[levels(ht$exposure_unmatched)=='ARB and ACEi'] <- NA
levels(ht$exposure_unmatched)[levels(ht$exposure_unmatched)=='Spironolactone'] <- NA
table(ht$exposure_unmatched)
print(epitab(ht$exposure_unmatched, ht$outcome, rev = "rows"))


##### MATCHING #####
library(e1071)
citation("e1071")
packageVersion("e1071")

# Match an ARB to each ACEi
ht$acei_vs_arb <- factor(ifelse(ht$ARB_or_ACEi=="ACEi", "case", 
                                  ifelse(ht$ARB_or_ACEi=="ARB", "control", NA)))
ht$acei_vs_control <- factor(ifelse(ht$ARB_or_ACEi=="ACEi", "case", 
                                ifelse(ht$ARB_or_ACEi=="Other", "control", NA)))
set.seed(100)
match_an_ARB <- matchControls(acei_vs_arb~age+
                                Gender_1male_2female+
                                Days_from_first_symptom_to_hospitalization+
                                smoking+
                                DM_History+
                                COPD_asthma_history+
                                Coronary_Artery_Disease+
                                Congestive_Heart_Failure+
                                Kre, 
                              data = ht, caselabel = "case", contlabel = "control")
match_a_control <- matchControls(acei_vs_control~age+
                                   Gender_1male_2female+
                                   Days_from_first_symptom_to_hospitalization+
                                   smoking+
                                   DM_History+
                                   COPD_asthma_history+
                                   Coronary_Artery_Disease+
                                   Congestive_Heart_Failure+
                                   Kre, 
                                 data = ht, caselabel = "case", contlabel = "control")
m <- ht[c(match_an_ARB$cases,match_an_ARB$controls,match_a_control$controls),]
m$exposure <- factor(m$ARB_or_ACEi)
dim(m)
table(m$exposure)
table(m$exposure, m$outcome)

kruskal.test(age~exposure_unmatched, data = ht)
kruskal.test(age~Days_from_first_symptom_to_hospitalization , data = ht)

table(ht$acei_vs_arb, ht$COPD_asthma_history)
prop.table(table(ht$acei_vs_arb, ht$COPD_asthma_history))
fisher.test(table(ht$acei_vs_arb, ht$COPD_asthma_history))

table(ht$exposure_unmatched, ht$COPD_asthma_history)
prop.table(table(ht$exposure_unmatched, ht$COPD_asthma_history))
fisher.test(table(ht$exposure_unmatched, ht$COPD_asthma_history))


ht$residuals <- NA
ht[which(ht$ARB_or_ACEi=="ACEi"),"residuals"] <- "ACEi"
ht[which(ht$ARB_or_ACEi=="ARB" & (row.names(ht) %in% match_an_ARB$controls)),"residuals"] <- "ARB_matched"
ht[which(ht$ARB_or_ACEi=="ARB" & !(row.names(ht) %in% match_an_ARB$controls)),"residuals"] <- "ARB_residual"
ht[which(ht$ARB_or_ACEi=="Other" & (row.names(ht) %in% match_a_control$controls)),"residuals"] <- "Other_matched"
ht[which(ht$ARB_or_ACEi=="Other" & !(row.names(ht) %in% match_a_control$controls)),"residuals"] <- "Other_residual"
ht$residuals <- factor(ht$residuals)
table(ht$residuals)


# Check match quality:
densitycurves("age", data=m, groups = "exposure")
prism_plot("age", "exposure", mydata = m)
table(m$exposure, m$Gender_1male_2female);fisher.test(table(m$exposure, m$Gender_1male_2female))
densitycurves("Days_from_first_symptom_to_hospitalization", data=m, groups = "exposure")
prism_plot("Days_from_first_symptom_to_hospitalization", "exposure", mydata = m)
kruskal.test(Days_from_first_symptom_to_hospitalization~exposure, data = m)
table(m$exposure, m$smoking)
table(m$exposure, m$DM_History)
table(m$exposure, m$COPD_asthma_history)
table(m$exposure, m$Coronary_Artery_Disease)
table(m$exposure, m$Congestive_Heart_Failure)
densitycurves("Kre", data=m, groups = "exposure")
prism_plot("Kre", "exposure", mydata = m)
kruskal.test(Kre~exposure, data = m)

### PRIMARY ANALYSIS ###
table(m$exposure, m$outcome)
fisher.test(table(m$exposure, m$outcome))
print(epitab(m$exposure, m$outcome, rev="rows"))
table(m$exposure, m$CT_category)
epitab(table(m$exposure, m$CT_category)[,c(1,2)], rev = "rows")

prism_plot("CT_ordinal", "exposure", mydata = m)
kruskal.test(CT_ordinal~exposure, data = m)

prism_plot("Kre", "exposure", mydata=m)


#### Odds Ratio of Labored breathing as a finding 
# (I apologize dearly for the inconvenient parameter coding 
#("d"yspnea for shortness of breath, "D"yspnea for labored breathing))
table(covdata$Dyspnea)
epitab(table(covdata$Dyspnea, covdata$outcome))

#### Odds Ratio of Respiratory Rate >= 22
table(covdata$Respiratory_Rate>=22, covdata$outcome)
epitab(table(covdata$Respiratory_Rate>=22, covdata$outcome)) # I tuned this manually to maximize the odds ratio
prism_plot("Respiratory_Rate", varx = "outcome", new_device =F) #### Plotting the Respiratory Rate
densitycurves("Respiratory_Rate", data = covdata)

#### Odds Ratios of Comorbidities:
epitab(table(covdata$Hypertension_history, covdata$outcome))
epitab(table(covdata$DM_History, covdata$outcome))
epitab(table(covdata$COPD_asthma_history, covdata$outcome))
epitab(table(covdata$Coronary_Artery_Disease, covdata$outcome))
epitab(table(covdata$Congestive_Heart_Failure, covdata$outcome))
epitab(table(covdata$Solid_malignancy_enc, covdata$outcome))
epitab(table(covdata$Hematologic_malignancy, covdata$outcome))
epitab(table(covdata$malignancy_any, covdata$outcome))
epitab(table(covdata$smoking, covdata$outcome))


prism_plot("age", "residuals", mydata=ht)
prism_plot("Days_from_first_symptom_to_hospitalization", "residuals", mydata=ht)
table(ht$DM_History, ht$residuals)
table(ht$COPD_asthma_history, ht$residuals)
table(ht$Coronary_Artery_Disease, ht$residuals)
table(ht$Congestive_Heart_Failure, ht$residuals)
table(ht$Solid_malignancy_enc, ht$residuals)
table(ht$smoking, ht$residuals)
