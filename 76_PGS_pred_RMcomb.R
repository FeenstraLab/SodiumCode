# among thiazide users, compare the risk of hyponatremia (e.g., <=125, <=130, <135) between those at different PGS. E.g. just as continous or comparing by centiles.
#it would be looking at thiazide users with an on-drug measurement within 120 days (and taking the first if multiple records)

library(optparse)
library(ggplot2)
library(dplyr)
# library(gridExtra)
# library(MASS)
# library(nnet)
library(pROC)
parsed_opts <- list(
  make_option("--PGS",
              action = "store",
              default = NA,
              type = 'character',
              help = 'PRS'),
    make_option("--THI",
              action = "store",
              default = NA,
              type = 'character',
              help = 'THI sodium  data'),
    make_option("--CCB",
              action = "store",
              default = NA,
              type = 'character',
              help = 'CCB sodium data'),
    make_option("--RAS",
              action = "store",
              default = NA,
              type = 'character',
              help = 'RAS sodium data'),
    make_option("--out_plot",
              action = "store",
              default = NA,
              type = 'character',
              help = ' AUC plot'),
    make_option("--out_plot2",
              action = "store",
              default = NA,
              type = 'character',
              help = ' odd ratio plot'),
    make_option("--auc",
              action = "store",
              default = NA,
              type = 'character',
              help = ' AUC results'),
    make_option("--effect",
              action = "store",
              default = NA,
              type = 'character',
              help = ' odd ratio results')
  
)








if (interactive()) {

} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
}

# input sodium data
  THI <- read.table(opt$THI,header=T) %>%
  dplyr::select(cpr_enc,EKSD,cohort,sex,on_drug,pre_drug,SAMPLINGDATE_on,SAMPLINGDATE_pre,ln_diff, diff,thiazide_age,PC1,PC2,PC3,PC4,PC5) %>%
   mutate(drug="THI")    

  CCB <- read.table(opt$CCB,header=T) %>%
    dplyr::select(cpr_enc,EKSD,cohort,sex,on_drug,pre_drug,SAMPLINGDATE_on,SAMPLINGDATE_pre,ln_diff, diff,thiazide_age,PC1,PC2,PC3,PC4,PC5) %>%
    mutate(drug="CCB")

  RAS <- read.table(opt$RAS,header=T) %>%
     dplyr::select(cpr_enc,EKSD,cohort,sex,on_drug,pre_drug,SAMPLINGDATE_on,SAMPLINGDATE_pre,ln_diff, diff,thiazide_age,PC1,PC2,PC3,PC4,PC5) %>%
   mutate(drug="RAS")


(nrow(THI)) 
(nrow(CCB)) 
(nrow(RAS)) 

# input PGS data
  PGS <- read.table(opt$PGS,header=T) %>%
  mutate(PGS=-PGS) %>%
    mutate(zscore=(PGS-mean(PGS))/sd(PGS)) %>%
    filter(zscore>=-5 & zscore<=5)
  ## quantiles
  # get_quantile <- function(x, num.quant, quant.ref){
  #   quant <- as.numeric(cut(x,
  #                           breaks = unique(quantile(
  #                             x, probs = seq(0, 1, 1 / num.quant)
  #                           )),
  #                           include.lowest = T))
  #   if(is.null(quant.ref)){
  #     quant.ref <- ceiling(num.quant / 2)       #define reference factor
  #   }
  #   quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
  #   return(quant)
  # }


  # PGS$quantile <- get_quantile(PGS$PGS,num.quant=10, quant.ref=1)


 # q975 <- qnorm(0.975)
 
  ##select lowest 0.1%/1%/10% of PGS
 PGS <- PGS[order(PGS$zscore),]
 PGS$lowest_per25 <- NA
 PGS[1:(nrow(PGS)*0.25),"lowest_per25"] <- 1
 PGS[((nrow(PGS)*0.25)+1):nrow(PGS),"lowest_per25"] <- 2
  
 PGS$lowest_per1 <- NA
 PGS[1:(nrow(PGS)*0.01),"lowest_per1"] <- 1
 PGS[((nrow(PGS)*0.01)+1):nrow(PGS),"lowest_per1"] <- 2
  
 PGS$lowest_per10 <- NA
 PGS[1:(nrow(PGS)*0.1),"lowest_per10"] <- 1
 PGS[((nrow(PGS)*0.1)+1):nrow(PGS),"lowest_per10"] <- 2
  
  
 PGS$lowVShigh_per25 <- NA
 PGS[1:(nrow(PGS)*0.25),"lowVShigh_per25"] <- 1
 PGS[(nrow(PGS)-(nrow(PGS)*0.25)+1):nrow(PGS),"lowVShigh_per25"] <- 2
 
 # PGS$lowVShigh_per1 <- NA
 # PGS[1:(nrow(PGS)*0.01),"lowVShigh_per1"] <- 1
 # PGS[(nrow(PGS)-(nrow(PGS)*0.01)+1):nrow(PGS),"lowVShigh_per1"] <- 2
  
 PGS$lowVShigh_per10 <- NA
 PGS[1:(nrow(PGS)*0.1),"lowVShigh_per10"] <- 1
 PGS[(nrow(PGS)-(nrow(PGS)*0.1)+1):nrow(PGS),"lowVShigh_per10"] <- 2


## merge

## compare EKSD data for 3 drugs and keep only the earliest 
df1 <- rbind(THI[,c("cpr_enc","EKSD","SAMPLINGDATE_on","SAMPLINGDATE_pre","drug")],CCB[,c("cpr_enc","EKSD","SAMPLINGDATE_on","SAMPLINGDATE_pre","drug")])
df2 <- rbind(df1,RAS[,c("cpr_enc","EKSD","SAMPLINGDATE_on","SAMPLINGDATE_pre","drug")])

##remove comorbidites of s_heart_failure_in2y_120d, s_other_CNS_in2y_120d, s_liver_peritonitis_in2y_120d, s_pancreatitis_in2y_120d, s_renal_in2y_120d, s_polydipsia_in2y_120d, alcohol_drug_in2y_120d,


df_mindate <- df2 %>%
  mutate(EKSD=as.Date(EKSD,format="%Y-%m-%d")) %>%
  group_by(cpr_enc) %>%
  arrange(drug) %>%
  summarise(min_date=min(EKSD),drug_no=paste0(unique(drug),collapse="-"))
df3 <- merge(df2,df_mindate,by="cpr_enc") %>%
  mutate(SAMPLINGDATE_on=as.Date(SAMPLINGDATE_on,format="%Y-%m-%d")) %>%
  mutate(SAMPLINGDATE_pre=as.Date(SAMPLINGDATE_pre,format="%Y-%m-%d")) 


##========================== THI
print("how many thiazide user has used THI only")
THIonly <- unique(df3[df3$drug=="THI" & df3$drug_no=="THI","cpr_enc"])
print(length(THIonly))

print("how many thiazide user has used THI and then used CCB/RAS")

THI_first <- df3 %>%
  filter(drug=="THI" & EKSD==min_date & drug_no!="THI" ) %>%
  dplyr::select(cpr_enc) %>%
  unique() 

#If individuals start with THI, and then change to RAS/CCB before sodium test, remove them;
#if individuals start with THI, and then change to RAS/CCB after sodium test, keep it as THI user.
THI_bf_CCBorRAS <- df3 %>%
  filter(cpr_enc %in% THI_first$cpr_enc) %>%
  group_by(cpr_enc) %>%
  mutate(min_samplingdate=min(SAMPLINGDATE_on)) %>%
  filter(drug!="THI" & SAMPLINGDATE_on>min_samplingdate) %>%
  dplyr::select(cpr_enc) %>%
  unique()

## ========================== CCB
print("how many thiazide user has used THI only")
CCBonly <- unique(df3[df3$drug=="CCB" & df3$drug_no=="CCB","cpr_enc"])
print(length(CCBonly))

print("how many thiazide user has used THI and then used CCB/RAS")

CCB_first <- df3 %>%
  filter(drug=="CCB" & EKSD==min_date & drug_no!="CCB" ) %>%
  dplyr::select(cpr_enc) %>%
  unique() 

#If individuals start with THI, and then change to RAS/CCB before sodium test, remove them;
#if individuals start with THI, and then change to RAS/CCB after sodium test, keep it as THI user.
CCB_bf_THIorRAS <- df3 %>%
  filter(cpr_enc %in% CCB_first$cpr_enc) %>%
  group_by(cpr_enc) %>%
  mutate(min_samplingdate=min(SAMPLINGDATE_on)) %>%
  filter(drug!="CCB" & SAMPLINGDATE_on>min_samplingdate) %>%
  dplyr::select(cpr_enc) %>%
  unique()


## ========================== RAS
print("how many thiazide user has used THI only")
RASonly <- unique(df3[df3$drug=="RAS" & df3$drug_no=="RAS","cpr_enc"])
print(length(RASonly))

print("how many thiazide user has used THI and then used THI/CCB")

RAS_first <- df3 %>%
  filter(drug=="RAS" & EKSD==min_date & drug_no!="RAS" ) %>%
  dplyr::select(cpr_enc) %>%
  unique() 

#If individuals start with RAS, and then change to THI/CCB before sodium test, remove them;
#if individuals start with RAS, and then change to THI/CCB after sodium test, keep it as THI user.
RAS_bf_THIorCCB <- df3 %>%
  filter(cpr_enc %in% RAS_first$cpr_enc) %>%
  group_by(cpr_enc) %>%
  mutate(min_samplingdate=min(SAMPLINGDATE_on)) %>%
  filter(drug!="RAS" & SAMPLINGDATE_on>min_samplingdate) %>%
  dplyr::select(cpr_enc) %>%
  unique()
print(nrow(RAS_bf_THIorCCB))



## summarise
THI_ID <- c(THIonly,THI_bf_CCBorRAS$cpr_enc)
CCB_ID <- c(CCBonly,CCB_bf_THIorRAS$cpr_enc)
RAS_ID <- c(RASonly,RAS_bf_THIorCCB$cpr_enc)
THI_ID_2 <- THI_ID[!(THI_ID %in% c(CCB_ID,RAS_ID))]
CCB_ID_2 <- CCB_ID[!(CCB_ID %in% c(THI_ID,RAS_ID))]
RAS_ID_2 <- RAS_ID[!(RAS_ID %in% c(THI_ID,CCB_ID))]
sum(THI_ID_2 %in% CCB_ID_2)
sum(THI_ID_2 %in% RAS_ID_2)
sum(RAS_ID_2 %in% CCB_ID_2)

print("Final number of THI, CCB and RAS")
(length(unique(THI_ID_2)))
(length(unique(CCB_ID_2)))
(length(unique(RAS_ID_2)))

df10 <- rbind(THI[THI$cpr_enc %in% THI_ID_2,],CCB[CCB$cpr_enc %in% CCB_ID_2,])
df20 <- rbind(df10,RAS[RAS$cpr_enc %in% RAS_ID_2,])


dat <- merge(df20,PGS,by="cpr_enc") %>%
    mutate(drug_le130=case_when(on_drug<=130 ~ 1,
    on_drug>130 ~ 0,
    TRUE ~ NA
    ))


##============================binormial logsitic regression for dat_THI
out <- as.data.frame(matrix(NA,4,5))
dat_THI <- dat[dat$drug=="THI",]
dat_THI$drug_le130 <- as.factor(dat_THI$drug_le130)
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(dat_THI), replace=TRUE, prob=c(0.7,0.3))
train <- dat_THI[sample, ]
test <- dat_THI[!sample, ] 


## model


m1 <- glm(drug_le130 ~ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m1))
predicted <- predict(m1, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[1,1] <- "Sex + age + PCs"
out[1,2] <- CI[[1]]
out[1,3] <- CI[[2]]
out[1,4] <- CI[[3]]
out[1,5] <- "THI data"


m2 <- glm(drug_le130 ~ pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m2))
predicted <- predict(m2, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[2,1] <- "+ pre_drug"
out[2,2] <- CI[[1]]
out[2,3] <- CI[[2]]
out[2,4] <- CI[[3]]
out[2,5] <- "THI data"

m3 <- glm(drug_le130 ~ zscore + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m3))
predicted <- predict(m3, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[3,1] <- "+ PGS"
out[3,2] <- CI[[1]]
out[3,3] <- CI[[2]]
out[3,4] <- CI[[3]]
out[3,5] <- "THI data"

m5 <- glm(drug_le130 ~ lowest_per25 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m5))
predicted <- predict(m5, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[5,1] <- "+ PGS25%"
out[5,2] <- CI[[1]]
out[5,3] <- CI[[2]]
out[5,4] <- CI[[3]]
out[5,5] <- "THI data"

m6 <- glm(drug_le130 ~ lowest_per1 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m6))
predicted <- predict(m6, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[6,1] <- "+ PGS1%"
out[6,2] <- CI[[1]]
out[6,3] <- CI[[2]]
out[6,4] <- CI[[3]]
out[6,5] <- "THI data"

m7 <- glm(drug_le130 ~ lowest_per10 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m7))
predicted <- predict(m7, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[7,1] <- "+ PGS10%"
out[7,2] <- CI[[1]]
out[7,3] <- CI[[2]]
out[7,4] <- CI[[3]]
out[7,5] <- "THI data"


m4 <- glm(drug_le130 ~ zscore + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m4))
predicted <- predict(m4, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[4,1] <- "+ PGS + pre_drug"
out[4,2] <- CI[[1]]
out[4,3] <- CI[[2]]
out[4,4] <- CI[[3]]
out[4,5] <- "THI data"

m8 <- glm(drug_le130 ~ lowest_per25 + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m8))
predicted <- predict(m8, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[8,1] <- "+ PGS25% + pre_drug"
out[8,2] <- CI[[1]]
out[8,3] <- CI[[2]]
out[8,4] <- CI[[3]]
out[8,5] <- "THI data"

m9 <- glm(drug_le130 ~ lowest_per1 + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m9))
predicted <- predict(m9, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[9,1] <- "+ PGS1% + pre_drug"
out[9,2] <- CI[[1]]
out[9,3] <- CI[[2]]
out[9,4] <- CI[[3]]
out[9,5] <- "THI data"

m10 <- glm(drug_le130 ~ lowest_per10 + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m10))
predicted <- predict(m10, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[10,1] <- "+ PGS10% + pre_drug"
out[10,2] <- CI[[1]]
out[10,3] <- CI[[2]]
out[10,4] <- CI[[3]]
out[10,5] <- "THI data"



# m11 <- glm(drug_le130 ~ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m11))
# predicted <- predict(m11, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[11,1] <- "Sex + age + PCs"
# out[11,2] <- CI[[1]]
# out[11,3] <- CI[[2]]
# out[11,4] <- CI[[3]]
# out[11,5] <- "THI data"



m12 <- glm(drug_le130 ~ lowVShigh_per25 + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m12))
predicted <- predict(m12, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[12,1] <- "+ PGS25%_LvsH + pre_drug"
out[12,2] <- CI[[1]]
out[12,3] <- CI[[2]]
out[12,4] <- CI[[3]]
out[12,5] <- "THI data"

# m13 <- glm(drug_le130 ~ lowVShigh_per1 + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m13))
# predicted <- predict(m13, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[13,1] <- "+ PGS1%_LvsH + pre_drug"
# out[13,2] <- CI[[1]]
# out[13,3] <- CI[[2]]
# out[13,4] <- CI[[3]]
# out[13,5] <- "THI data"

m14 <- glm(drug_le130 ~ lowVShigh_per10 + pre_drug+ sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m14))
predicted <- predict(m14, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[14,1] <- "+ PGS10%_LvsH + pre_drug"
out[14,2] <- CI[[1]]
out[14,3] <- CI[[2]]
out[14,4] <- CI[[3]]
out[14,5] <- "THI data"

m15 <- glm(drug_le130 ~ lowest_per25 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m15))
predicted <- predict(m15, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[15,1] <- "+ PGS25%_LvsH"
out[15,2] <- CI[[1]]
out[15,3] <- CI[[2]]
out[15,4] <- CI[[3]]
out[15,5] <- "THI data"

m16 <- glm(drug_le130 ~ lowest_per1 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m16))
predicted <- predict(m16, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[16,1] <- "+ PGS1%_LvsH"
out[16,2] <- CI[[1]]
out[16,3] <- CI[[2]]
out[16,4] <- CI[[3]]
out[16,5] <- "THI data"

m17 <- glm(drug_le130 ~ lowest_per10 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m17))
predicted <- predict(m17, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[17,1] <- "+ PGS10%_LvsH"
out[17,2] <- CI[[1]]
out[17,3] <- CI[[2]]
out[17,4] <- CI[[3]]
out[17,5] <- "THI data"


names(out) <- c("Covariates","Lower","AUC","Upper","dataset")
print(out)

OUT <- out

##============================binormial logsitic regression for all dat(THI,RAS,CCB)
out <- as.data.frame(matrix(NA,4,5))
dat$drug_le130 <- as.factor(dat$drug_le130)
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(dat), replace=TRUE, prob=c(0.7,0.3))
train <- dat[sample, ]
test <- dat[!sample, ] 


## model
m1 <- glm(drug_le130 ~ drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m1))
predicted <- predict(m1, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[1,1] <- "Sex + age + PCs + drug"
out[1,2] <- CI[[1]]
out[1,3] <- CI[[2]]
out[1,4] <- CI[[3]]
out[1,5] <- "THI+CCB+RAS"


m2 <- glm(drug_le130 ~ pre_drug+ drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m2))
predicted <- predict(m2, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[2,1] <- "+ pre_drug"
out[2,2] <- CI[[1]]
out[2,3] <- CI[[2]]
out[2,4] <- CI[[3]]
out[2,5] <- "THI+CCB+RAS"

m3 <- glm(drug_le130 ~ zscore + drug + zscore:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m3))
predicted <- predict(m3, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[3,1] <- "+ PGS"
out[3,2] <- CI[[1]]
out[3,3] <- CI[[2]]
out[3,4] <- CI[[3]]
out[3,5] <- "THI+CCB+RAS"


m5 <- glm(drug_le130 ~ lowest_per25 +drug +lowest_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m5))
predicted <- predict(m5, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[5,1] <- "+ PGS25%"
out[5,2] <- CI[[1]]
out[5,3] <- CI[[2]]
out[5,4] <- CI[[3]]
out[5,5] <- "THI+CCB+RAS"

m6 <- glm(drug_le130 ~ lowest_per1 + drug+lowest_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m6))
predicted <- predict(m6, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[6,1] <- "+ PGS1%"
out[6,2] <- CI[[1]]
out[6,3] <- CI[[2]]
out[6,4] <- CI[[3]]
out[6,5] <- "THI+CCB+RAS"

m7 <- glm(drug_le130 ~ lowest_per10 + drug + lowest_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m7))
predicted <- predict(m7, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[7,1] <- "+ PGS10%"
out[7,2] <- CI[[1]]
out[7,3] <- CI[[2]]
out[7,4] <- CI[[3]]
out[7,5] <- "THI+CCB+RAS"

m4 <- glm(drug_le130 ~ zscore + pre_drug+ drug + zscore:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m4))
predicted <- predict(m4, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[4,1] <- "+ PGS + pre_drug"
out[4,2] <- CI[[1]]
out[4,3] <- CI[[2]]
out[4,4] <- CI[[3]]
out[4,5] <- "THI+CCB+RAS"



m8 <- glm(drug_le130 ~ lowest_per25 + pre_drug+ drug +lowest_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m8))
predicted <- predict(m8, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[8,1] <- "+ PGS25% + pre_drug"
out[8,2] <- CI[[1]]
out[8,3] <- CI[[2]]
out[8,4] <- CI[[3]]
out[8,5] <- "THI+CCB+RAS"

m9 <- glm(drug_le130 ~ lowest_per1 + pre_drug+ drug+lowest_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m9))
predicted <- predict(m9, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[9,1] <- "+ PGS1% + pre_drug"
out[9,2] <- CI[[1]]
out[9,3] <- CI[[2]]
out[9,4] <- CI[[3]]
out[9,5] <- "THI+CCB+RAS"

m10 <- glm(drug_le130 ~ lowest_per10 + pre_drug+drug+ lowest_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m10))
predicted <- predict(m10, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[10,1] <- "+ PGS10% + pre_drug"
out[10,2] <- CI[[1]]
out[10,3] <- CI[[2]]
out[10,4] <- CI[[3]]
out[10,5] <- "THI+CCB+RAS"

# m11 <- glm(drug_le130 ~ drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m11))
# predicted <- predict(m11, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[11,1] <- "Sex + age + PCs + drug"
# out[11,2] <- CI[[1]]
# out[11,3] <- CI[[2]]
# out[11,4] <- CI[[3]]
# out[11,5] <- "THI+CCB+RAS"



m12 <- glm(drug_le130 ~ lowVShigh_per25 +drug +lowVShigh_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m12))
predicted <- predict(m12, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[12,1] <- "+ PGS25%_LvsH"
out[12,2] <- CI[[1]]
out[12,3] <- CI[[2]]
out[12,4] <- CI[[3]]
out[12,5] <- "THI+CCB+RAS"

# m13 <- glm(drug_le130 ~ lowVShigh_per1 + drug+lowVShigh_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m13))
# predicted <- predict(m13, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[13,1] <- "+ PGS1%_LvsH"
# out[13,2] <- CI[[1]]
# out[13,3] <- CI[[2]]
# out[13,4] <- CI[[3]]
# out[13,5] <- "THI+CCB+RAS"

m14 <- glm(drug_le130 ~ lowVShigh_per10 + drug + lowVShigh_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m14))
predicted <- predict(m14, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[14,1] <- "+ PGS10%_LvsH"
out[14,2] <- CI[[1]]
out[14,3] <- CI[[2]]
out[14,4] <- CI[[3]]
out[14,5] <- "THI+CCB+RAS"

m15 <- glm(drug_le130 ~ lowVShigh_per25 + pre_drug+ drug +lowVShigh_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m15))
predicted <- predict(m15, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[15,1] <- "+ PGS25%_LvsH + pre_drug"
out[15,2] <- CI[[1]]
out[15,3] <- CI[[2]]
out[15,4] <- CI[[3]]
out[15,5] <- "THI+CCB+RAS"

# m16 <- glm(drug_le130 ~ lowVShigh_per1 + pre_drug+ drug+lowVShigh_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m16))
# predicted <- predict(m16, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[16,1] <- "+ PGS1%_LvsH + pre_drug"
# out[16,2] <- CI[[1]]
# out[16,3] <- CI[[2]]
# out[16,4] <- CI[[3]]
# out[16,5] <- "THI+CCB+RAS"

m17 <- glm(drug_le130 ~ lowVShigh_per10 + pre_drug+drug+ lowVShigh_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m17))
predicted <- predict(m17, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[17,1] <- "+ PGS10%_LvsH + pre_drug"
out[17,2] <- CI[[1]]
out[17,3] <- CI[[2]]
out[17,4] <- CI[[3]]
out[17,5] <- "THI+CCB+RAS"

names(out) <- c("Covariates","Lower","AUC","Upper","dataset")
print(out)


(OUT <- rbind(OUT,out))


##============================binormial logsitic regression for THI+RAS)
out <- as.data.frame(matrix(NA,4,5))
dat_THI_RAS <- dat[dat$drug!="CCB",]
dat_THI_RAS$drug_le130 <- as.factor(dat_THI_RAS$drug_le130)
set.seed(1)

#Use 70% of dataset as training set and remaining 30% as testing set
sample <- sample(c(TRUE, FALSE), nrow(dat_THI_RAS), replace=TRUE, prob=c(0.7,0.3))
train <- dat_THI_RAS[sample, ]
test <- dat_THI_RAS[!sample, ] 


## model
m1 <- glm(drug_le130 ~ drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m1))
predicted <- predict(m1, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[1,1] <- "Sex + age + PCs + drug"
out[1,2] <- CI[[1]]
out[1,3] <- CI[[2]]
out[1,4] <- CI[[3]]
out[1,5] <- "THI+RAS"


m2 <- glm(drug_le130 ~ pre_drug+ drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m2))
predicted <- predict(m2, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[2,1] <- "+ pre_drug"
out[2,2] <- CI[[1]]
out[2,3] <- CI[[2]]
out[2,4] <- CI[[3]]
out[2,5] <- "THI+RAS"

m3 <- glm(drug_le130 ~ zscore + drug + zscore:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m3))
predicted <- predict(m3, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[3,1] <- "+ PGS"
out[3,2] <- CI[[1]]
out[3,3] <- CI[[2]]
out[3,4] <- CI[[3]]
out[3,5] <- "THI+RAS"


m5 <- glm(drug_le130 ~ lowest_per25 +drug +lowest_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m5))
predicted <- predict(m5, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[5,1] <- "+ PGS25%"
out[5,2] <- CI[[1]]
out[5,3] <- CI[[2]]
out[5,4] <- CI[[3]]
out[5,5] <- "THI+RAS"

m6 <- glm(drug_le130 ~ lowest_per1 + drug+lowest_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m6))
predicted <- predict(m6, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[6,1] <- "+ PGS1%"
out[6,2] <- CI[[1]]
out[6,3] <- CI[[2]]
out[6,4] <- CI[[3]]
out[6,5] <- "THI+RAS"

m7 <- glm(drug_le130 ~ lowest_per10 + drug + lowest_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m7))
predicted <- predict(m7, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[7,1] <- "+ PGS10%"
out[7,2] <- CI[[1]]
out[7,3] <- CI[[2]]
out[7,4] <- CI[[3]]
out[7,5] <- "THI+RAS"

m4 <- glm(drug_le130 ~ zscore + pre_drug+ drug + zscore:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m4))
predicted <- predict(m4, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[4,1] <- "+ PGS + pre_drug"
out[4,2] <- CI[[1]]
out[4,3] <- CI[[2]]
out[4,4] <- CI[[3]]
out[4,5] <- "THI+RAS"



m8 <- glm(drug_le130 ~ lowest_per25 + pre_drug+ drug +lowest_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m8))
predicted <- predict(m8, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[8,1] <- "+ PGS25% + pre_drug"
out[8,2] <- CI[[1]]
out[8,3] <- CI[[2]]
out[8,4] <- CI[[3]]
out[8,5] <- "THI+RAS"

m9 <- glm(drug_le130 ~ lowest_per1 + pre_drug+ drug+lowest_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m9))
predicted <- predict(m9, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[9,1] <- "+ PGS1% + pre_drug"
out[9,2] <- CI[[1]]
out[9,3] <- CI[[2]]
out[9,4] <- CI[[3]]
out[9,5] <- "THI+RAS"

m10 <- glm(drug_le130 ~ lowest_per10 + pre_drug+drug+ lowest_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m10))
predicted <- predict(m10, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[10,1] <- "+ PGS10% + pre_drug"
out[10,2] <- CI[[1]]
out[10,3] <- CI[[2]]
out[10,4] <- CI[[3]]
out[10,5] <- "THI+RAS"

# m11 <- glm(drug_le130 ~ drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m11))
# predicted <- predict(m11, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[11,1] <- "Sex + age + PCs + drug"
# out[11,2] <- CI[[1]]
# out[11,3] <- CI[[2]]
# out[11,4] <- CI[[3]]
# out[11,5] <- "THI+RAS"



m12 <- glm(drug_le130 ~ lowVShigh_per25 +drug +lowVShigh_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m12))
predicted <- predict(m12, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[12,1] <- "+ PGS25%_LvsH"
out[12,2] <- CI[[1]]
out[12,3] <- CI[[2]]
out[12,4] <- CI[[3]]
out[12,5] <- "THI+RAS"

# m13 <- glm(drug_le130 ~ lowVShigh_per1 + drug+lowVShigh_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m13))
# predicted <- predict(m13, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[13,1] <- "+ PGS1%_LvsH"
# out[13,2] <- CI[[1]]
# out[13,3] <- CI[[2]]
# out[13,4] <- CI[[3]]
# out[13,5] <- "THI+RAS"

m14 <- glm(drug_le130 ~ lowVShigh_per10 + drug + lowVShigh_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m14))
predicted <- predict(m14, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[14,1] <- "+ PGS10%_LvsH"
out[14,2] <- CI[[1]]
out[14,3] <- CI[[2]]
out[14,4] <- CI[[3]]
out[14,5] <- "THI+RAS"

m15 <- glm(drug_le130 ~ lowVShigh_per25 + pre_drug+ drug +lowVShigh_per25:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m15))
predicted <- predict(m15, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[15,1] <- "+ PGS25%_LvsH + pre_drug"
out[15,2] <- CI[[1]]
out[15,3] <- CI[[2]]
out[15,4] <- CI[[3]]
out[15,5] <- "THI+RAS"

# m16 <- glm(drug_le130 ~ lowVShigh_per1 + pre_drug+ drug+lowVShigh_per1:drug +sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
# (summary(m16))
# predicted <- predict(m16, test, type="response")
# (auc(test$drug_le130, predicted))
# CI <- ci.auc(test$drug_le130, predicted)
# out[16,1] <- "+ PGS1%_LvsH + pre_drug"
# out[16,2] <- CI[[1]]
# out[16,3] <- CI[[2]]
# out[16,4] <- CI[[3]]
# out[16,5] <- "THI+RAS"

m17 <- glm(drug_le130 ~ lowVShigh_per10 + pre_drug+drug+ lowVShigh_per10:drug+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=train)
(summary(m17))
predicted <- predict(m17, test, type="response")
(auc(test$drug_le130, predicted))
CI <- ci.auc(test$drug_le130, predicted)
out[17,1] <- "+ PGS10%_LvsH + pre_drug"
out[17,2] <- CI[[1]]
out[17,3] <- CI[[2]]
out[17,4] <- CI[[3]]
out[17,5] <- "THI+RAS"

names(out) <- c("Covariates","Lower","AUC","Upper","dataset")
print(out)

OUT <- rbind(OUT,out)

OUT <- OUT[!is.na(OUT$AUC),]


#OUT$Upper <- as.numeric(OUT$Upper)



## plot
OUT$groups <- NA
#OUT[OUT$Covariates %in% c("Sex + age + PCs","Sex + age + PCs + drug","Sex + age + PCs + comorbidities","Sex + age + PCs + drug + comorbidities"),"groups"] <- 1
OUT[OUT$Covariates %in% c("Sex + age + PCs","Sex + age + PCs + drug"),"groups"] <- 1

OUT[OUT$Covariates=="+ pre_drug","groups"] <- 2
OUT[OUT$Covariates %in% c("+ PGS","+ PGS25%","+ PGS1%","+ PGS10%"),"groups"] <- 3
OUT[OUT$Covariates %in% c("+ PGS25%_LvsH","+ PGS1%_LvsH","+ PGS10%_LvsH"),"groups"] <- 4
OUT[OUT$Covariates %in% c("+ PGS + pre_drug","+ PGS25% + pre_drug","+ PGS1% + pre_drug","+ PGS10% + pre_drug"),"groups"] <- 5
OUT[OUT$Covariates %in% c("+ PGS25%_LvsH + pre_drug","+ PGS1%_LvsH + pre_drug","+ PGS10%_LvsH + pre_drug"),"groups"] <- 6

OUT$groups <- as.factor(OUT$groups)

#OUT$Covariates <- factor(OUT$Covariates,levels=c("Sex + age + PCs","Sex + age + PCs + comorbidities","Sex + age + PCs + drug","Sex + age + PCs + drug + comorbidities","+ PGS","+ PGS1%","+ PGS10%","+ PGS25%","+ PGS1%_LvsH","+ PGS10%_LvsH","+ PGS25%_LvsH","+ pre_drug","+ PGS + pre_drug","+ PGS1% + pre_drug","+ PGS10% + pre_drug","+ PGS25% + pre_drug","+ PGS1%_LvsH + pre_drug","+ PGS10%_LvsH + pre_drug","+ PGS25%_LvsH + pre_drug"))
OUT$Covariates <- factor(OUT$Covariates,levels=c("Sex + age + PCs","Sex + age + PCs + drug","+ PGS","+ PGS1%","+ PGS10%","+ PGS25%","+ PGS1%_LvsH","+ PGS10%_LvsH","+ PGS25%_LvsH","+ pre_drug","+ PGS + pre_drug","+ PGS1% + pre_drug","+ PGS10% + pre_drug","+ PGS25% + pre_drug","+ PGS1%_LvsH + pre_drug","+ PGS10%_LvsH + pre_drug","+ PGS25%_LvsH + pre_drug"))

#OUT$Covariates <- as.character(OUT$Covariates)


write.table(OUT,opt$auc,col.names=T,row.names=F,quote=F,sep="\t")


myplot <- ggplot(OUT,aes(x=AUC,y=Covariates,xmin=Lower,xmax=Upper,color=groups))+
  geom_point()+
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73","#F0E442", "#0072B2"))+
  geom_errorbarh(height=0.1) +
#  geom_vline(xintercept=0.675,color="grey",linetype="dashed")+
  facet_grid(~dataset,scales="free",space="free")+
  xlab("AUC(95% CI)")+
  ylab("Covariates")
  
ggsave(opt$out_plot,myplot,height=4,width=10)



##===============================================use logistic regression for risk of hyponatremia (<=130) comparing individuals with lowest 0.1%/1%/10% of PGS wih remaining. 
dat_THI <- dat[dat$drug=="THI",]

effect <- as.data.frame(matrix(NA,4,8))
m1 <- glm(drug_le130 ~ zscore + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)
m2 <- glm(drug_le130 ~ lowest_per25 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)
m3 <- glm(drug_le130 ~ lowest_per1 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)
m4 <- glm(drug_le130 ~ lowest_per10 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)
m5 <- glm(drug_le130 ~ lowVShigh_per25 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)
#m6 <- glm(drug_le130 ~ lowVShigh_per1 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)
m7 <- glm(drug_le130 ~ lowVShigh_per10 + sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI)


CI <- as.data.frame(confint.default(m1))
effect[1,1] <- "PGS"
effect[1,2] <- exp(coef(summary(m1))["zscore",1])
effect[1,3:6] <- coef(summary(m1))["zscore",]
effect[1,7:8] <- exp(CI["zscore",])



CI <- as.data.frame(confint.default(m2))
effect[2,1] <- "PGS25%"
effect[2,2] <- exp(coef(summary(m2))["lowest_per25",1])
effect[2,3:6] <- coef(summary(m2))["lowest_per25",]
effect[2,7:8] <- exp(CI["lowest_per25",])

CI <- as.data.frame(confint.default(m3))
effect[3,1] <- "PGS1%"
effect[3,2] <- exp(coef(summary(m3))["lowest_per1",1])
effect[3,3:6] <- coef(summary(m3))["lowest_per1",]
effect[3,7:8] <- exp(CI["lowest_per1",])

CI <- as.data.frame(confint.default(m4))
effect[4,1] <- "PGS10%"
effect[4,2] <- exp(coef(summary(m4))["lowest_per10",1])
effect[4,3:6] <- coef(summary(m4))["lowest_per10",]
effect[4,7:8] <- exp(CI["lowest_per10",])


CI <- as.data.frame(confint.default(m5))
effect[5,1] <- "PGS25%_LvsH"
effect[5,2] <- exp(coef(summary(m5))["lowVShigh_per25",1])
effect[5,3:6] <- coef(summary(m5))["lowVShigh_per25",]
effect[5,7:8] <- exp(CI["lowVShigh_per25",])

# CI <- as.data.frame(confint.default(m6))
# effect[6,1] <- "PGS1%_LvsH"
# effect[6,2] <- exp(coef(summary(m6))["lowVShigh_per1",1])
# effect[6,3:6] <- coef(summary(m6))["lowVShigh_per1",]
# effect[6,7:8] <- exp(CI["lowVShigh_per1",])

CI <- as.data.frame(confint.default(m7))
effect[7,1] <- "PGS10%_LvsH"
effect[7,2] <- exp(coef(summary(m7))["lowVShigh_per10",1])
effect[7,3:6] <- coef(summary(m7))["lowVShigh_per10",]
effect[7,7:8] <- exp(CI["lowVShigh_per10",])


names(effect) <- c("Covariates","OR","beta","se","t","P","exp_Lower","exp_Upper")
effect$drug <- "THI"

EFFECT <- effect


##=====THI+RAS+CCB
effect <- as.data.frame(matrix(NA,4,8))
m1 <- glm(drug_le130 ~ zscore + drug+zscore:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)
m2 <- glm(drug_le130 ~ lowest_per25 + drug+lowest_per25:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)
m3 <- glm(drug_le130 ~ lowest_per1 + drug+lowest_per1:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)
m4 <- glm(drug_le130 ~ lowest_per10 + drug+lowest_per10:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)
m5 <- glm(drug_le130 ~ lowVShigh_per25 + drug+lowVShigh_per25:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)
# m6 <- glm(drug_le130 ~ lowVShigh_per1 + drug+lowVShigh_per1:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)
m7 <- glm(drug_le130 ~ lowVShigh_per10 + drug+lowVShigh_per10:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat)


CI <- as.data.frame(confint.default(m1))
effect[1,1] <- "PGS"
effect[1,2] <- exp(coef(summary(m1))["zscore",1])
effect[1,3:6] <- coef(summary(m1))["zscore",]
effect[1,7:8] <- exp(CI["zscore",])



CI <- as.data.frame(confint.default(m2))
effect[2,1] <- "PGS25%"
effect[2,2] <- exp(coef(summary(m2))["lowest_per25",1])
effect[2,3:6] <- coef(summary(m2))["lowest_per25",]
effect[2,7:8] <- exp(CI["lowest_per25",])

CI <- as.data.frame(confint.default(m3))
effect[3,1] <- "PGS1%"
effect[3,2] <- exp(coef(summary(m3))["lowest_per1",1])
effect[3,3:6] <- coef(summary(m3))["lowest_per1",]
effect[3,7:8] <- exp(CI["lowest_per1",])

CI <- as.data.frame(confint.default(m4))
effect[4,1] <- "PGS10%"
effect[4,2] <- exp(coef(summary(m4))["lowest_per10",1])
effect[4,3:6] <- coef(summary(m4))["lowest_per10",]
effect[4,7:8] <- exp(CI["lowest_per10",])


CI <- as.data.frame(confint.default(m5))
effect[5,1] <- "PGS25%_LvsH"
effect[5,2] <- exp(coef(summary(m5))["lowVShigh_per25",1])
effect[5,3:6] <- coef(summary(m5))["lowVShigh_per25",]
effect[5,7:8] <- exp(CI["lowVShigh_per25",])

# CI <- as.data.frame(confint.default(m6))
# effect[6,1] <- "PGS1%_LvsH"
# effect[6,2] <- exp(coef(summary(m6))["lowVShigh_per1",1])
# effect[6,3:6] <- coef(summary(m6))["lowVShigh_per1",]
# effect[6,7:8] <- exp(CI["lowVShigh_per1",])

CI <- as.data.frame(confint.default(m7))
effect[7,1] <- "PGS10%_LvsH"
effect[7,2] <- exp(coef(summary(m7))["lowVShigh_per10",1])
effect[7,3:6] <- coef(summary(m7))["lowVShigh_per10",]
effect[7,7:8] <- exp(CI["lowVShigh_per10",])


names(effect) <- c("Covariates","OR","beta","se","t","P","exp_Lower","exp_Upper")
effect$drug <- "THI+RAS+CCB"

EFFECT <- rbind(EFFECT,effect)

##=====THI+RAS
dat_THI_RAS <- dat[dat$drug!="CCB",]
effect <- as.data.frame(matrix(NA,4,8))
m1 <- glm(drug_le130 ~ zscore + drug+zscore:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)
m2 <- glm(drug_le130 ~ lowest_per25 + drug+lowest_per25:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)
m3 <- glm(drug_le130 ~ lowest_per1 + drug+lowest_per1:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)
m4 <- glm(drug_le130 ~ lowest_per10 + drug+lowest_per10:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)
m5 <- glm(drug_le130 ~ lowVShigh_per25 + drug+lowVShigh_per25:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)
# m6 <- glm(drug_le130 ~ lowVShigh_per1 + drug+lowVShigh_per1:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)
m7 <- glm(drug_le130 ~ lowVShigh_per10 + drug+lowVShigh_per10:PGS+sex +thiazide_age + PC1+PC2+PC3+PC4+PC5,family = "binomial",data=dat_THI_RAS)


CI <- as.data.frame(confint.default(m1))
effect[1,1] <- "PGS"
effect[1,2] <- exp(coef(summary(m1))["zscore",1])
effect[1,3:6] <- coef(summary(m1))["zscore",]
effect[1,7:8] <- exp(CI["zscore",])



CI <- as.data.frame(confint.default(m2))
effect[2,1] <- "PGS25%"
effect[2,2] <- exp(coef(summary(m2))["lowest_per25",1])
effect[2,3:6] <- coef(summary(m2))["lowest_per25",]
effect[2,7:8] <- exp(CI["lowest_per25",])

CI <- as.data.frame(confint.default(m3))
effect[3,1] <- "PGS1%"
effect[3,2] <- exp(coef(summary(m3))["lowest_per1",1])
effect[3,3:6] <- coef(summary(m3))["lowest_per1",]
effect[3,7:8] <- exp(CI["lowest_per1",])

CI <- as.data.frame(confint.default(m4))
effect[4,1] <- "PGS10%"
effect[4,2] <- exp(coef(summary(m4))["lowest_per10",1])
effect[4,3:6] <- coef(summary(m4))["lowest_per10",]
effect[4,7:8] <- exp(CI["lowest_per10",])


CI <- as.data.frame(confint.default(m5))
effect[5,1] <- "PGS25%_LvsH"
effect[5,2] <- exp(coef(summary(m5))["lowVShigh_per25",1])
effect[5,3:6] <- coef(summary(m5))["lowVShigh_per25",]
effect[5,7:8] <- exp(CI["lowVShigh_per25",])

# CI <- as.data.frame(confint.default(m6))
# effect[6,1] <- "PGS1%_LvsH"
# effect[6,2] <- exp(coef(summary(m6))["lowVShigh_per1",1])
# effect[6,3:6] <- coef(summary(m6))["lowVShigh_per1",]
# effect[6,7:8] <- exp(CI["lowVShigh_per1",])

CI <- as.data.frame(confint.default(m7))
effect[7,1] <- "PGS10%_LvsH"
effect[7,2] <- exp(coef(summary(m7))["lowVShigh_per10",1])
effect[7,3:6] <- coef(summary(m7))["lowVShigh_per10",]
effect[7,7:8] <- exp(CI["lowVShigh_per10",])


names(effect) <- c("Covariates","OR","beta","se","t","P","exp_Lower","exp_Upper")
effect$drug <- "THI+RAS"

(EFFECT <- rbind(EFFECT,effect))

EFFECT <- EFFECT[!is.na(EFFECT$OR),]

write.table(EFFECT,opt$effect,col.names=T,row.names=F,quote=F,sep="\t")


myplot2 <- ggplot(EFFECT[EFFECT$se<1.1,],aes(x=OR,y=Covariates,xmin=exp_Lower,xmax=exp_Upper))+
  geom_point()+
#  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","#009E73"))+
  geom_errorbarh(height=0.1) +
  geom_vline(xintercept=1,color="grey",linetype="dashed")+
  facet_grid(~drug,scales="free",space="free")+
  xlab("OR(95% CI)")+
  ylab("Covariates")
  
ggsave(opt$out_plot2,myplot2,height=4,width=10)

