source("functions.R")

#### GENERATE TABLE 1 ####
table1 <- function(indep = "ICU_yesno", mydata = disc_or_icu){
  # This function will automatically stitch together the numbers on table 1 
  mytbl <- data.frame("Variable"="varx", "group1"="group1", "group2"="group2", "pvalue"="pvalue", "marker"="marker",
                      stringsAsFactors = F)
  gt <- table(mydata$Gender_1male_2female, mydata$ICU_yesno) #for gender
  gt1 <- paste0(gt[2,1],"F:",gt[1,1],"M") # for gender
  gt2 <- paste0(gt[2,2],"F:",gt[1,2],"M") # for gender
  fgt <- fisher.test(gt)
  mytbl<- rbind(
    mytbl, 
    c("N", as.character(table(mydata[indep])),""),
    probs2("RT_PCR", indep, mydata = mydata), 
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


######## GENERATE AND SAVE MANY MANY GRAPHS FOR THE EXPLORATION OF DATA (POST-HOC ANALYSES) #########

### Define the parametric columns:
parametric_columns <- c("age", "BMI", "Days.hospitalized", "Days.from.first.symptom.to.hospitalization",
                        "Days.since.first.symptom", "Systolic", "Diastolic", "Pulse", "Respiratory.Rate",
                        "pH", "spo2", "pO2", "pCo2", "HCO3", "Lactate", "Hgb", "Plt", "WBC", "Neut", "Lymp", 
                        "Mon", "Eos", 
                        "BUN", "Kre", "Na", "Cl", "K", "Glu", "AST", "ALT", "GGT", "ALP", "LDH", 
                        "T.prot", "Alb", "CRP", "Prokalsitonin", "Ferritin", "Trig", "ESR", "D.dimer", 
                        "Trop", "pro.BNP", "Fibrinojen", "INR", "APTT", "CT_arbitrary_score")
dir.create("output") # creates an output folder for the graphs to be generated
sapply(parametric_columns, function(x)prism_plot(x, "ICU_yesno", save_graph=T))

#### The graphs are now in the output folder. Good luck and stay healthy!

