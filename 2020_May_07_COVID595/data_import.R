# Data import script

#### DATA IMPORT ####
impdata <- read.csv("COVID_DATA 07.05.20.csv", stringsAsFactors = F, row.names = 1) #import dataset
dim(impdata) # Should have 499 rows 26.April.2020

impdata[impdata==""] <- NA
# names(impdata) <- gsub(x = names(impdata), pattern = "ü", replacement = "u")
names(impdata) <- gsub(x = names(impdata), pattern = "\\.", replacement = "_")
names(impdata) <- gsub(x = names(impdata), pattern = "\\__", replacement = "_")

impdata$Date_admitted <- as.Date(impdata$Date_admitted, "%d.%m.%y")
impdata$Taburcu_tarihi <- as.Date(impdata$Taburcu_tarihi, "%d.%m.%y")
impdata$Ex <- factor(impdata$Ex, levels = c(0,1), labels = c("Alive", "Excitus"))
impdata$Ex_tarihi <- as.Date(impdata$Ex_tarihi, "%d.%m.%y")
impdata$Days_hospitalized_on_day_x <- ifelse(impdata$Taburcu==1,
                                             difftime(impdata$Taburcu_tarihi, impdata$Date_admitted, tz = 'UTC', units = "days"),
                                             ifelse(impdata$Ex == 1,
                                                    difftime(impdata$Ex_tarihi, impdata$Date_admitted, tz = 'UTC', units = "days"),
                                                    difftime("2020-05-07", impdata$Date_admitted, tz = 'UTC', units = "days")))
impdata$iCU_yesno <- factor(impdata$iCU_yesno, levels = c(0,1), labels = c("no", "yes")) # Ever admitted to ICU during this hospital stay?
impdata$Gender_1male_2female <- factor(impdata$Gender_1male_2female, levels = c(1,2), labels = c("male", "female"))
impdata$smoking <- factor(impdata$smoking, levels=c(0,1) ,labels = c("no", "yes"))
impdata$alcohol <- factor(impdata$alcohol, levels = c(0,1), labels = c("no", "yes")) # alcohol consumption is probably poorly collected in this data set
impdata$Travel_history <- factor(impdata$Travel_history, levels = c(0,1), labels = c("no", "yes")) # Interesting fact: travel history is non-existent in more recent data due to country-wide travel bans
impdata$Contact_history <- factor(impdata$Contact_history, levels = c(0,1), labels = c("no", "yes"))


# Encoding RT_PCR test results
impdata$RT_PCR <- factor(ifelse(impdata$RT_PCR==1,"positive","not positive"))
impdata$RT_PCR[is.na(impdata$RT_PCR)] <- "not positive"

# Encoding CT scores:
impdata$CT_arbitrary_score[impdata$CT_arbitrary_score %in% c("", "CHECK TYPING")] <- NA
impdata$CT_category <- as.character(impdata$CT_arbitrary_score)
impdata$CT_ordinal <- as.numeric(impdata$CT_arbitrary_score)
impdata$CT_arbitrary_score <- factor(impdata$CT_arbitrary_score, 
                                     levels = c(0,1,2,3,4,5,6,7),
                                     labels = c("Normal or non-COVID",
                                                "very mild", "mild", "mild-moderate",
                                                "moderate", "moderate-diffuse", "diffuse",
                                                "very diffuse"))
impdata$CT_category[impdata$CT_category %in% c("0", "1", "2", "3")] <- "Normal or Mild"
impdata$CT_category[impdata$CT_category %in% c("4", "5", "6", "7")] <- "Moderate"
impdata$CT_category[impdata$CT_category %in% c("6", "7")] <- "Diffuse"
impdata$CT_category[impdata$CT_category %in% c("", "CHECK TYPING")] <- NA
impdata$CT_category <- factor(impdata$CT_category, levels = c("Normal or Mild", "Moderate", "Diffuse"),
                              labels = c("Normal or Mild", "Moderate", "Diffuse"))


# Factoring symptoms
impdata$fever_reported <- factor(impdata$fever_reported, levels = c(0,1), labels = c("no", "yes"))
impdata$coughing <- factor(impdata$coughing, levels = c(0,1), labels = c("no", "yes"))
impdata$sputum <- factor(impdata$sputum, levels = c(0,1), labels = c("no", "yes"))
impdata$dyspnea <- factor(impdata$dyspnea, levels = c(0,1), labels = c("no", "yes"))
impdata$fatigue_or_myalgia <- factor(impdata$fatigue_or_myalgia, levels = c(0,1), labels = c("no", "yes"))
impdata$nausea <- factor(impdata$nausea, levels = c(0,1), labels = c("no", "yes"))
impdata$diarrhea <- factor(impdata$diarrhea, levels = c(0,1), labels = c("no", "yes"))
impdata$anosmia <- factor(impdata$anosmia, levels = c(0,1), labels = c("no", "yes"))

# Factoring Hypertension and anti-hypertensive medication
impdata$Hypertension_history <- factor(impdata$Hypertension_history, levels = c(0,1), labels = c("no", "yes"))
impdata$ACEi <- factor(impdata$ACEi, levels = c(0,1), labels = c("other or ARB", "ACEi"))
impdata$ARB <- factor(impdata$ARB, levels = c(0,1), labels = c("other or ACEi", "ARB"))
impdata$Spironolactone <- factor(impdata$Spironolactone, levels = c(0,1), labels = c("other", "Spironolactone"))
impdata$ARB_or_ACEi <- ifelse(impdata$ARB=="ARB" & impdata$ACEi=="ACEi","ARB and ACEi",
                              ifelse(impdata$ARB=="ARB", "ARB",
                                     ifelse(impdata$ACEi=="ACEi", "ACEi",
                                            ifelse(impdata$Spironolactone=="Spironolactone","Spironolactone",
                                                   ifelse(impdata$Hypertension_history=="yes", "Other", "")))))
impdata$ARB_or_ACEi <- factor(impdata$ARB_or_ACEi)
impdata$RAAS <- ifelse(impdata$ARB =="ARB" | impdata$ACEi=="ACEi" | impdata$Spironolactone == "Spironolactone",
                       "RAAS",
                       ifelse(impdata$Hypertension_history=="yes", "Other", ""))
impdata$RAAS <- factor(impdata$RAAS, levels = c("Other", "RAAS"), labels = c("Other", "RAAS"))


#Factoring other comorbidities:
impdata$DM_History <- factor(impdata$DM_History, levels = c(0,1), labels = c("no", "yes"))
impdata$COPD_asthma_history <- factor(impdata$COPD_asthma_history, levels = c(0,1), labels = c("no", "yes"))
impdata$Coronary_Artery_Disease <- factor(impdata$Coronary_Artery_Disease, levels = c(0,1), labels = c("no", "yes"))
impdata$Congestive_Heart_Failure <- factor(impdata$Congestive_Heart_Failure, levels = c(0,1), labels = c("no", "yes"))
impdata$Solid_malignancy_enc <- ifelse(impdata$Solid_malignancy==0,0,1) #encode solid malignancy as binary
impdata$Hematologic_malignancy <- ifelse(impdata$Hematologic_malignancy==0,0,1) #encode solid malignancy as binary
impdata$malignancy_any <- factor(ifelse(impdata$Solid_malignancy_enc==1 | impdata$Hematologic_malignancy == 1, 1, 0)) #encode any malignancy as factor

# COVID medication
impdata$Favipiravir_alma <- factor(impdata$Favipiravir_alma, levels = c(0,1), labels = c("no", "yes"))
impdata$Favi_start_date <- impdata$Favi_baslama_gun
impdata$Favi_start_date <- as.Date.character(impdata$Favi_start_date, "%d.%m.%y")
impdata$Favi_delay <- difftime(impdata$Favi_start_date, impdata$Date_admitted, units = "days")+
  impdata$Days_from_first_symptom_to_hospitalization


impdata$Toci_alma <- factor(impdata$Toci_alma, levels = c(0,1), labels = c("no", "yes"))
impdata$Toci_start_date <- as.Date(impdata$Toci_tedavi_gun, "%d.%m.%y")
impdata$Toci_delay <- difftime(impdata$Toci_start_date, impdata$Date_admitted, units = "days")+
  impdata$Days_from_first_symptom_to_hospitalization

impdata$Anakinra_alma <- factor(impdata$Anakinra_alma, levels = c(0,1), labels = c("no", "yes"))
impdata$Anakin_start_date <- as.Date(impdata$Anakinra_baslama_gun, "%d.%m.%y")
impdata$Anakin_delay <- difftime(impdata$Anakin_start_date, impdata$Date_admitted, units = "days")+
  impdata$Days_from_first_symptom_to_hospitalization

impdata$Meropenem_alma <- factor(impdata$Meropenem_alma, levels = c(0,1), labels = c("no", "yes"))

