#### Analyses for Istanbul Faculty of Medicine Corona-center deidentified patient data ####

source(file = "functions.R")

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
dim(alldata)

table(alldata$RT_PCR)
table(alldata$CT_collapse, alldata$RT_PCR)
#distribution of positive PCR results are distributed evenly between different CT categories
fisher.test(table(alldata$CT_collapse, alldata$RT_PCR)) 
# however, they are slightly higher in more severe patients. Higher viral load?

