# Table 1

#### GENERATE TABLE 1 ####

table(covdata$CT_category, covdata$outcome);fisher.test(table(covdata$CT_category, covdata$outcome))
table(covdata$number_of_anti_ht_yes_tzd, covdata$outcome);fisher.test(table(covdata$number_of_anti_ht_yes_tzd, covdata$outcome))

table1 <- function(indep = "outcome", mydata = covdata){
  # This function will automatically stitch together the numbers on table 1 
  mytbl <- data.frame("Variable"="varx", "group1"="group1", "group2"="group2", "pvalue"="pvalue", "marker"="marker",
                      stringsAsFactors = F)
  gt <- table(mydata$Gender_1male_2female, mydata$outcome) #for gender
  gt1 <- paste0(gt[2,1],"F:",gt[1,1],"M") # for gender
  gt2 <- paste0(gt[2,2],"F:",gt[1,2],"M") # for gender
  fgt <- fisher.test(gt)
  mytbl<- rbind(
    mytbl, 
    c("N", as.character(table(mydata[indep])),""),
    probs2("RT_PCR", indep, mydata = mydata), 
    myquart("age", indep, mydata = mydata, k=0),
    c("Gender", gt1, gt2, round(fgt$p.value, digits=4), rownames(gt)[2]),
    myquart("BMi", indep, mydata = mydata, k=1),
    myquart("Bmi_2", indep, mydata = mydata, k=1),
    probs2("Travel_history", indep, mydata = mydata),
    probs2("Contact_history", indep, mydata = mydata),
    probs2("fever_reported", indep, mydata = mydata),
    probs2("coughing", indep, mydata = mydata),
    probs2("sputum", indep, mydata = mydata),
    probs2("dyspnea", indep, mydata = mydata),
    probs2("fatigue_or_myalgia", indep, mydata = mydata),
    probs2("nausea", indep, mydata = mydata),
    probs2("diarrhea", indep, mydata = mydata),
    probs2("anosmia", indep, mydata = mydata),
    myquart("fever_days", indep, mydata = mydata, k=0),
    myquart("coughing_days", indep, mydata = mydata, k=0),
    myquart("sputum_days", indep, mydata = mydata, k=0),
    myquart("dyspnea_days", indep, mydata = mydata, k=0),
    myquart("fatigue_or_myalgia_days", indep, mydata = mydata, k=0),
    myquart("Days_from_first_symptom_to_hospitalization", indep, mydata = mydata, k=0),
    myquart("Days_hospitalized_on_day_x", indep, mydata = mydata, k=0),
    myquart("Days_since_first_symptom", indep, mydata = mydata, k=0),
    probs2("Hypertension_history", indep, mydata = mydata),
    myquart("number_of_anti_ht_yes_tzd", indep, mydata = mydata, k=0),
    probs2("DM_History", indep, mydata = mydata),
    probs2("COPD_asthma_history", indep, mydata = mydata),
    probs2("Coronary_Artery_Disease", indep, mydata = mydata),
    probs2("Congestive_Heart_Failure", indep, mydata = mydata),
    probs2("KBY", indep, mydata = mydata),
    probs2("Solid_malignancy_enc", indep, mydata = mydata),
    probs2("Hematologic_malignancy", indep, mydata = mydata),
    probs2("smoking", indep, mydata = mydata),
    myquart("Body_temp", indep, mydata = mydata, k=1),
    myquart("spo2", indep, mydata = mydata, k=0),
    myquart("Systolic", indep, mydata = mydata, k=0),
    myquart("Diastolic", indep, mydata = mydata, k=0),
    myquart("Pulse", indep, mydata = mydata, k=0),
    myquart("Respiratory_Rate", indep, mydata = mydata, k=0),
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
    myquart("T_prot", indep, mydata = mydata, k=1),
    myquart("Alb", indep, mydata = mydata, k=1),
    myquart("CRP", indep, mydata = mydata, k=0),
    myquart("Prokalsitonin", indep, mydata = mydata, k=2),
    myquart("Ferritin", indep, mydata = mydata, k=0),
    myquart("D_dimer", indep, mydata = mydata, k=0),
    myquart("Trop", indep, mydata = mydata, k=1),
    myquart("pro_BNP", indep, mydata = mydata, k=0),
    myquart("Fibrinojen", indep, mydata = mydata, k=0),
    myquart("iNR", indep, mydata = mydata, k=1),
    myquart("APTT", indep, mydata = mydata, k=0)
  )
  return(mytbl)
}
print(table1())
t1 <- table1()
write.csv(t1, file = "table_output/table_1.csv")