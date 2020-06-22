# Table 2

#### GENERATE TABLE 2 #########

table(m$outcome, m$exposure); fisher.test(table(m$outcome, m$exposure))
table(m$RT_PCR, m$exposure)


prism_plot("CT_category_numeric", "exposure", m)
densitycurves("CT_category_numeric", groups="exposure", data = m)
table(m$CT_category, m$exposure); fisher.test(table(m$CT_category, m$exposure))
fisher.test(table(m$CT_category_numeric>4, m$exposure))
kruskal.test(CT_category_numeric~exposure, data = m)
wilcox.test(CT_category_numeric~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(CT_category_numeric~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(CT_category_numeric~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))

wilcox.test(Pulse~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(Pulse~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(Pulse~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("Pulse", "exposure", m)


wilcox.test(spo2~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(spo2~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(spo2~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("spo2", "exposure", m)


wilcox.test(Lymp~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(Lymp~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(Lymp~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("Lymp", "exposure", m)

wilcox.test(Mon~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(Mon~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(Mon~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("Mon", "exposure", m)

wilcox.test(CRP~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(CRP~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(CRP~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("CRP", "exposure", m)

wilcox.test(Ferritin~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(Ferritin~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(Ferritin~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("Ferritin", "exposure", m)

wilcox.test(APTT~exposure, data = subset(m, m$exposure%in%c("Other", "ACEi")))
wilcox.test(APTT~exposure, data = subset(m, m$exposure%in%c("Other", "ARB")))
wilcox.test(APTT~exposure, data = subset(m, m$exposure%in%c("ACEi", "ARB")))
prism_plot("APTT", "exposure", m)

m$Mero_or_Pip <- ifelse(m$Meropenem_alma == "yes" | m$Pip_tazo_alma == 1, 1, 0)
m$any_drug <- ifelse(m$Meropenem_alma == "yes" | 
                       m$Pip_tazo_alma == 1 |
                       m$Favipiravir_alma == "yes" |
                       m$Toci_alma == "yes"|
                       m$Anakinra_alma == "yes", 1, 0)
m$immuno_modulator <- ifelse(m$Toci_alma == "yes"| m$Anakinra_alma == "yes", 1, 0)
table(m$Mero_or_Pip, m$exposure); fisher.test(table(m$Mero_or_Pip, m$exposure))
table(m$any_drug, m$exposure); fisher.test(table(m$any_drug, m$exposure))
table(m$immuno_modulator, m$exposure); fisher.test(table(m$immuno_modulator, m$exposure))


table(m$CT_category, m$exposure);fisher.test(table(m$CT_category, m$exposure))
table(m$CT_category, m$exposure)[,c(2,3)];fisher.test(table(m$CT_category, m$exposure)[,c(2,3)])
table(m$CT_category, m$exposure)[,c(1,3)];fisher.test(table(m$CT_category, m$exposure)[,c(1,3)])
table(m$CT_category, m$exposure)[,c(1,2)];fisher.test(table(m$CT_category, m$exposure)[,c(1,2)])
table(m$number_of_anti_ht_yes_tzd, m$exposure);fisher.test(table(m$number_of_anti_ht_yes_tzd, m$exposure))
describeBy(m$Systolic, m$exposure)


table(m$Ca_channel_blocker_type, m$exposure)
table(m$Ca_channel_blocker_type_enc, m$exposure)
k <- subset(m, m$Ca_channel_blocker_type_enc=="Yes");table(k$outcome, k$exposure)

table(m$Alpha1_blocker, m$exposure)
table(m$Alpha1_blocker_enc, m$exposure)

table(m$Beta_Blocker, m$exposure)
table(m$Beta_Blocker_enc, m$exposure)
k <- subset(m, m$Beta_Blocker_enc=="Yes");table(k$outcome, k$exposure)

table(m$Thiazide, m$exposure)
table(m$Thiazide_enc, m$exposure)
k <- subset(m, m$Thiazide_enc=="Yes");table(k$outcome, k$exposure)

table(m$ARB_type, m$exposure)
table(m$ACEi_type, m$exposure)

table(m$ACEi_type, m$outcome); fisher.test(table(m$ACEi_type, m$outcome)[c(3,4),]) # Perindopril may have been the better drug, but we did not have the power to detect this
table(m$ARB_type, m$outcome)

table1 <- function(indep = "exposure", mydata = m){
  # This function will automatically stitch together the numbers on table 1 
  mytbl <- data.frame("Variable"="varx", "group1"="group1", "group2"="group2", "group3"="group3", "pvalue"="pvalue", "marker"="marker",
                      stringsAsFactors = F)
  gt <- table(mydata[,"Gender_1male_2female"], mydata[,indep]) #for gender
  gt1 <- paste0(gt[2,1],"F:",gt[1,1],"M") # for gender
  gt2 <- paste0(gt[2,2],"F:",gt[1,2],"M") # for gender
  gt3 <- paste0(gt[2,3],"F:",gt[1,3],"M") # for gender
  fgt <- fisher.test(gt)
  test_type_m="Kruskal"
  mytbl<- rbind(
    mytbl, 
    c("N", as.character(table(mydata[indep])),""),
    probs2("RT_PCR", indep, mydata = mydata, n_columns = 3), 
    myquart("CT_category_numeric", indep, mydata = mydata, k=0, test_type = test_type_m),
    probs2("Loop_diuretic", indep, mydata = mydata, n_columns = 3),
    probs2("Ca_channel_blocker_type_enc", indep, mydata = mydata, n_columns = 3),
    probs2("Alpha1_blocker_enc", indep, mydata = mydata, n_columns = 3), 
    probs2("Beta_Blocker_enc", indep, mydata = mydata, n_columns = 3), 
    probs2("Thiazide_enc", indep, mydata = mydata, n_columns = 3),  
    probs2("Favipiravir_alma", indep, mydata = mydata, n_columns = 3),
    probs2("Toci_alma", indep, mydata = mydata, n_columns = 3),  
    probs2("Anakinra_alma", indep, mydata = mydata, n_columns = 3),  
    probs2("Meropenem_alma", indep, mydata = mydata, n_columns = 3),  
    probs2("Pip_tazo_alma", indep, mydata = mydata, n_columns = 3),  
    myquart("age", indep, mydata = mydata, k=0, test_type = test_type_m),
    c("Gender", gt1, gt2, gt3, round(fgt$p.value, digits=4), rownames(gt)[2]),
    myquart("BMi", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("Bmi_2", indep, mydata = mydata, k=1, test_type = test_type_m),
    probs2("Travel_history", indep, mydata = mydata, n_columns = 3),
    probs2("Contact_history", indep, mydata = mydata, n_columns = 3),
    probs2("fever_reported", indep, mydata = mydata, n_columns = 3),
    probs2("coughing", indep, mydata = mydata, n_columns = 3),
    probs2("sputum", indep, mydata = mydata, n_columns = 3),
    probs2("dyspnea", indep, mydata = mydata, n_columns = 3),
    probs2("fatigue_or_myalgia", indep, mydata = mydata, n_columns = 3),
    probs2("nausea", indep, mydata = mydata, n_columns = 3),
    probs2("diarrhea", indep, mydata = mydata, n_columns = 3),
    probs2("anosmia", indep, mydata = mydata, n_columns = 3),
    myquart("fever_days", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("coughing_days", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("sputum_days", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("dyspnea_days", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("fatigue_or_myalgia_days", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Days_from_first_symptom_to_hospitalization", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Days_hospitalized_on_day_x", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Days_since_first_symptom", indep, mydata = mydata, k=0, test_type = test_type_m),
    probs2("Hypertension_history", indep, mydata = mydata, n_columns = 3),
    probs2("DM_History", indep, mydata = mydata, n_columns = 3),
    probs2("COPD_asthma_history", indep, mydata = mydata, n_columns = 3),
    probs2("Coronary_Artery_Disease", indep, mydata = mydata, n_columns = 3),
    probs2("Congestive_Heart_Failure", indep, mydata = mydata, n_columns = 3),
    probs2("KBY", indep, mydata = mydata, n_columns = 3),
    probs2("Solid_malignancy_enc", indep, mydata = mydata, n_columns = 3),
    probs2("Hematologic_malignancy", indep, mydata = mydata, n_columns = 3),
    probs2("smoking", indep, mydata = mydata, n_columns = 3),
    myquart("Body_temp", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("spo2", indep, mydata = mydata, k=0, test_type = test_type_m),
    probs2("spo2_93", indep, mydata = mydata, n_columns = 3),
    myquart("Systolic", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Diastolic", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Pulse", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Respiratory_Rate", indep, mydata = mydata, k=0, test_type = test_type_m),
    probs2("Dyspnea", indep, mydata = mydata, k=0, n_columns = 3),
    myquart("pH", indep, mydata = mydata, k=2, test_type = test_type_m),
    myquart("pO2", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("pCo2", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("HCO3", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Lactate", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("Hgb", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("Plt", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("WBC", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Neut", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Lymp", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Mon", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("BUN", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Kre", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("Na", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Cl", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("K", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("Glu", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("AST", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("ALT", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("GGT", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("ALP", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("LDH", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("T_prot", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("Alb", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("CRP", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Prokalsitonin", indep, mydata = mydata, k=2, test_type = test_type_m),
    myquart("Ferritin", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("D_dimer", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Trop", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("pro_BNP", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("Fibrinojen", indep, mydata = mydata, k=0, test_type = test_type_m),
    myquart("iNR", indep, mydata = mydata, k=1, test_type = test_type_m),
    myquart("APTT", indep, mydata = mydata, k=0, test_type = test_type_m),
    probs2("outcome", indep, mydata = mydata, n_columns = 3),
    probs2("Ex", indep, mydata = mydata, n_columns = 3),
    probs2("iCU_yesno", indep, mydata = mydata, n_columns = 3),
    probs2("long_hospitalization", indep, mydata = mydata, n_columns = 3)
  )
  return(mytbl)
}
print(table1())
t2 <- table1()
write.csv(t2, file = "table_output/table_2.csv")


