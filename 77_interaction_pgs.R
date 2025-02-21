# target: plot and tables for effect of PGS on sodium change and TIH ----

## prereq ----
library(optparse)
library(ggplot2)
library(dplyr)
library(tidyverse)
library(grid)
library(gridExtra)
library(ggsci)

#######################################################################


# input sodium data -----
dat <- read.table("for_THI_CCB_RAS_plot_wID.txt", header = T)

## Interaction models for ln_diff and on_drug for table output ---------


## PGS alone ----
model <- print("on_drug~PGS+covs")
comp <- print("THI")
(lr <- summary(lm(on_drug ~ zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("THI"),] %>% 
                    mutate(zscore_rev = -zscore))))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2),]

cov <- rep(NA, 2)
est <- rep(NA, 2)
se <- rep(NA,2)
p <- rep(NA,2)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("PGS", "THI", "INT") 

tab <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
tab

## PGS + drug alone ----
model <- print("on_drug~PGS+drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(lm(on_drug ~ drug + zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("THI","RAS"),]%>% 
                    mutate(zscore_rev = -zscore))))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3),]

cov <- rep(NA, 1)
est <- rep(NA, 1)
se <- rep(NA,1)
p <- rep(NA,1)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("THI","PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)


## model 1: THI vs RAS------
model <- print("ln_diff~PGS:drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,data=dat[dat$drug %in% c("THI","RAS"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

## model 1: THI vs CCB------
model <- print("ln_diff~PGS:drug+covs")
comp <- print("THI vs CCB")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,data=dat[dat$drug %in% c("THI","CCB"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

## model 2: THI vs RAS--------
model <- print("on_drug~PGS:drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,data=dat[dat$drug %in% c("THI","RAS"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

## model 2: THI vs CCB----
model <- print("on_drug~PGS:drug+covs")
comp <- print("THI vs CCB")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,data=dat[dat$drug %in% c("THI","CCB"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

tab

################## binary ----

### PGS alone ---
model <- print("<130~PGS+covs")
comp <- print("THI")
(lr <- summary(glm(hypona_130 ~ zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2),] 
cov <- rep(NA, 2)
est <- rep(NA, 2)
se <- rep(NA,2)
p <- rep(NA,2)
lr_df <- lr_df %>% rbind(
data.frame(cov, est, se, p))

lr_df[,1] <- c("PGS", "THI", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<125~PGS+covs")
comp <- print("THI")
(lr <- summary(glm(hypona_125 ~ zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2),] 
cov <- rep(NA, 2)
est <- rep(NA, 2)
se <- rep(NA,2)
p <- rep(NA,2)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("PGS", "THI", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<135~PGS+covs")
comp <- print("THI")
(lr <- summary(glm(hypona_134 ~ zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2),] 
cov <- rep(NA, 2)
est <- rep(NA, 2)
se <- rep(NA,2)
p <- rep(NA,2)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("PGS", "THI", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)


### PGS+THI not interaction ---
model <- print("<130~PGS+drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(glm(hypona_130 ~ drug +zscore+ sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","RAS"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3),] 
cov <- rep(NA, 1)
est <- rep(NA, 1)
se <- rep(NA,1)
p <- rep(NA,1)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<125~PGS+drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(glm(hypona_125 ~ drug +zscore+ sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","RAS"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3),] 
cov <- rep(NA, 1)
est <- rep(NA, 1)
se <- rep(NA,1)
p <- rep(NA,1)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<135~PGS+drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(glm(hypona_134 ~ drug +zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","RAS"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3),] 
cov <- rep(NA, 1)
est <- rep(NA, 1)
se <- rep(NA,1)
p <- rep(NA,1)
lr_df <- lr_df %>% rbind(
  data.frame(cov, est, se, p))

lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

##

model <- print("<130~PGS:drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","RAS"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0),
                       zscore_rev = -zscore
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<130~PGS:drug+covs")
comp <- print("THI vs CCB")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                  data=dat[dat$drug %in% c("THI","CCB"),] %>% 
                    mutate(
                      hypona_130 = ifelse(on_drug <= 130, 1, 0),
                      hypona_134 = ifelse(on_drug <= 134, 1, 0),
                      hypona_125 = ifelse(on_drug <= 125, 1, 0)
                    )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

tab

model <- print("<125~PGS:drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","RAS"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<125~PGS:drug+covs")
comp <- print("THI vs CCB")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","CCB"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

tab

model <- print("<135~PGS:drug+covs")
comp <- print("THI vs RAS")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","RAS"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("<135~PGS:drug+covs")
comp <- print("THI vs CCB")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("THI","CCB"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

tab
#################


## model 1+2 for HCTZ vs BFZ----
dat <- read.table("for_BFZ_HCTZ_CCB_RAS_plot_wID.txt", header = T)
model <- print("ln_diff~PGS:drug+covs")
comp <- print("HCTZ vs BFZ")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("BFZ","HCTZ"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("on_drug~PGS:drug+covs")
comp <- print("HCTZ vs BFZ")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("BFZ","HCTZ"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)


model <- print("<125~PGS:drug+covs")
comp <- print("HCTZ vs BFZ")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("BFZ","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)

model <- print("<135~PGS:drug+covs")
comp <- print("HCTZ vs BFZ")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("BFZ","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)

model <- print("<130~PGS:drug+covs")
comp <- print("HCTZ vs BFZ")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("BFZ","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)

## model 1+2 for HCTZ vs RAS----
dat <- read.table("for_BFZ_HCTZ_CCB_RAS_plot_wID.txt", header = T)
model <- print("ln_diff~PGS:drug+covs")
comp <- print("HCTZ vs RAS")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("HCTZ","RAS"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("on_drug~PGS:drug+covs")
comp <- print("HCTZ vs RAS")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("RAS","HCTZ"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)



model <- print("<135~PGS:drug+covs")
comp <- print("HCTZ vs RAS")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("RAS","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<130~PGS:drug+covs")
comp <- print("HCTZ vs RAS")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("RAS","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<125~PGS:drug+covs")
comp <- print("HCTZ vs RAS")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("RAS","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<135~PGS:drug+covs")
comp <- print("BFZ vs RAS")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("RAS","BFZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<130~PGS:drug+covs")
comp <- print("BFZ vs RAS")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("RAS","BFZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<125~PGS:drug+covs")
comp <- print("BFZ vs RAS")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("RAS","BFZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)

## model 1+2 for HCTZ vs RAS
dat <- read.table("for_BFZ_HCTZ_CCB_RAS_plot_wID.txt", header = T)
model <- print("ln_diff~PGS:drug+covs")
comp <- print("HCTZ vs CCB")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("HCTZ","CCB"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("on_drug~PGS:drug+covs")
comp <- print("HCTZ vs CCB")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("CCB","HCTZ"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)


model <- print("<135~PGS:drug+covs")
comp <- print("HCTZ vs CCB")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("CCB","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<130~PGS:drug+covs")
comp <- print("HCTZ vs CCB")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("CCB","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<125~PGS:drug+covs")
comp <- print("HCTZ vs CCB")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("CCB","HCTZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<135~PGS:drug+covs")
comp <- print("BFZ vs CCB")
(lr <- summary(glm(hypona_134 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("CCB","BFZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<130~PGS:drug+covs")
comp <- print("BFZ vs CCB")
(lr <- summary(glm(hypona_130 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("CCB","BFZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)


model <- print("<125~PGS:drug+covs")
comp <- print("BFZ vs CCB")
(lr <- summary(glm(hypona_125 ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                   family = "binomial",
                   data=dat[dat$drug %in% c("CCB","BFZ"),] %>% 
                     mutate(
                       hypona_130 = ifelse(on_drug <= 130, 1, 0),
                       hypona_134 = ifelse(on_drug <= 134, 1, 0),
                       hypona_125 = ifelse(on_drug <= 125, 1, 0)
                     )
)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)
lr_df
tab <- tab %>% rbind(lr_df)
## model 1+2 for BFZ vs CCB
dat <- read.table("for_BFZ_HCTZ_CCB_RAS_plot_wID.txt", header = T)
model <- print("ln_diff~PGS:drug+covs")
comp <- print("BFZ vs RAS")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("BFZ","RAS"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("on_drug~PGS:drug+covs")
comp <- print("BFZ vs RAS")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("BFZ","RAS"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

## model 1+2 for BFZ vs CCB
dat <- read.table("for_BFZ_HCTZ_CCB_RAS_plot_wID.txt", header = T)
model <- print("ln_diff~PGS:drug+covs")
comp <- print("BFZ vs CCB")
(lr <- summary(lm(ln_diff ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("BFZ","CCB"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("on_drug~PGS:drug+covs")
comp <- print("BFZ vs CCB")
(lr <- summary(lm(on_drug ~ drug+zscore+drug:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
                  data=dat[dat$drug %in% c("BFZ","CCB"),])))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)
## outputting intaraction table later ------
tab


# Density plot of PGS among thiazid users ----
dat <- read.table("for_THI_CCB_RAS_plot_wID.txt", header = T)

p1 <- dat[dat$drug %in% c("THI"),] %>% 
  mutate(drug=ifelse(drug=="THI", "Thiazide users", "RASi users"),
         drug = factor(drug, levels=c("Thiazide users", "RASi users"))) %>%
  ggplot(aes(x=zscore)) + 
  geom_density(fill="#6690b9", 
               color="#034788", 
               alpha=1,
               adjust=1.2) +
  scale_color_lancet()+ 
  scale_fill_lancet()+ 
  #scale_x_continuous(limits = c(-4,4), breaks = seq(-4,4,1))+
  #scale_y_continuous(limits = c(138,141), breaks = seq(138, 141, 1))+
  theme_classic()+
  labs(y="PGS density", x="Standardized PGS",
       tag="A")+
  theme(legend.position="top",
        legend.title = element_blank()#,
        #legend.text  = element_text(size=10, face="bold")
  ) + 
  coord_cartesian(xlim = c(-4.02,4.02), expand=F)
p1

## Quantile plot -----
get_quantile <- function(x, num.quant, quant.ref){
  quant <- as.numeric(cut(x,
                          breaks = unique(quantile(
                            x, probs = seq(0, 1, 1 / num.quant)
                          )),
                          include.lowest = T))
  if(is.null(quant.ref)){
    quant.ref <- ceiling(num.quant / 2)       #define reference factor
  }
  quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
  return(quant)
}

dat_thi <- dat %>% filter(drug=="THI")
dat_thi$quantile <- get_quantile(dat_thi$zscore,num.quant=10, quant.ref=1)

dbds_phe_sodium <- read.table("dbds_phe.tsv", header=T)
chb_phe_sodium <- read.table("chb_phe.tsv", header=T)
PGS <- read.table("PGS.tsv", header = T) %>%
  mutate(zscore=(PGS-mean(PGS))/sd(PGS)) %>%
  filter(zscore>=-5 & zscore<=5)

dbds <- merge(dbds_phe_sodium,PGS,by=1)
chb <- merge(chb_phe_sodium,PGS,by=1)
# effect of PGS on sodium
dbds_chb <- rbind(dbds,chb)
dbds_chb$quantile_all <- get_quantile(dbds_chb$zscore,num.quant=10, quant.ref=1)
PGS$quantile_all <- get_quantile(PGS$zscore,num.quant=10, quant.ref=1)

dat_thi_residuals <- dat_thi %>%
  left_join(dbds_chb %>% select(c("FID","phe", "quantile_all")), by = c("cpr_enc"="FID"))

lr <- summary(lm(phe ~ -1 + quantile_all ,data=dat_thi_residuals))                 ## TODO should it be on_drug instead??
lr_effect <- lr$coefficients[,1]
lr_sd <- lr$coefficients[,2]
q975 <- qnorm(0.975)
s <- cbind(c(1:10),lr_effect,lr_sd,lr_effect+q975*lr_sd,lr_effect-q975*lr_sd)
s <- as.data.frame(s)
names(s) <- c("Quantile","Effect","SD","CI.U","CI.L")
s$cohort <- "DBDS+CHB"
s
sPGS <- s

# 
# lr <- summary(lm(phe ~ -1 + quantile ,data=dbds_chb))
# lr_effect <- lr$coefficients[,1]
# lr_sd <- lr$coefficients[,2]
# q975 <- qnorm(0.975)
# s <- cbind(c(1:10),lr_effect,lr_sd,lr_effect+q975*lr_sd,lr_effect-q975*lr_sd)
# s <- as.data.frame(s)
# names(s) <- c("Quantile","Effect","SD","CI.U","CI.L")
# s$cohort <- "DBDS+CHB"
# s

## Percentile plot of the residuals of sodium among thiazide users by PGS -----
p2 <- ggplot(s,aes(x=Quantile,
                   y=Effect,                                            
                   ymin=CI.L,
                   ymax=CI.U,
                   color = factor(cohort)))+
  geom_point()+
  geom_errorbar(width=0.1,
                
  ) +
  xlab("Standardized polygenic score, centiles")+
  ylab("Effect (95% CI)")+
  theme_classic() + 
  scale_color_manual(values = c("#034788"))+
  scale_fill_manual(values = "#034788")+
  scale_x_continuous(breaks=seq(1,10,1)) +
  theme(legend.position = "none")+
  labs(tag="B")
p2

#### Interaction pre vs on THI sodium pvalues and plot
datx <- dat_thi_residuals %>% filter(drug == "THI") %>% select(cpr_enc, on_drug, pre_drug, sex, thiazide_age, quantile_all, zscore, PC1:PC5) %>% 
  pivot_longer(cols = c("pre_drug", "on_drug"), names_to= "names", values_to="values")

model <- print("lm(value~PGS:time+covs)")
comp <- print("THI pre vs on")

(lr <- summary(lm(values ~ names+zscore+names:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5 ,
                  data=datx)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,5)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

model <- print("lmer(value~PGS:time+covs|id)")
comp <- print("TIH pre vs on")
library(lme4)
library(lmerTest)
(lr <- summary(lmer(values ~ names+zscore+names:zscore + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5 + (1|cpr_enc),
                    data=datx)))

coefficient_names <- rownames(coef(lr))
lr_df <- data.frame(Coefficient = coefficient_names, Estimate = coef(lr)) 
lr_df <- lr_df[,c(1,2,3,6)]
colnames(lr_df)[1:4] <- c("cov", "est", "se", "p")
lr_df <- lr_df[c(2,3,11),]
lr_df[,1] <- c("THI", "PGS", "INT") 

lr_df <- lr_df %>% pivot_wider(names_from = "cov", values_from= c("est", "se", "p")) %>%
  mutate(model = model,
         comp = comp)

tab <- tab %>% rbind(lr_df)

####
txt <- paste0("P-value for interaction: ", sprintf(("%.2f"),tab[16,9])) # P-value from ln_diff model, adjusted for sex, age + PCs
txt <- as.data.frame(txt)
txt$x <- -2
txt$y <- 141.5
txt$drug <- "Pre-thiazide drug sodium level    "

px <- datx %>% 
  mutate(drug=ifelse(names=="pre_drug", "Pre-thiazide drug sodium level    ", "On-thiazide drug sodium level    "),
         drug = factor(drug, levels=c("On-thiazide drug sodium level    ", "Pre-thiazide drug sodium level    "))) %>%
  ggplot(aes(x=zscore,y=values,color=drug, shape=drug)) + 
  # geom_point(size=0.5) +
  geom_smooth(method=lm, aes(fill=drug),
              se=T, fullrange=TRUE,
              linewidth=0.6,
              alpha=0.4)+
  #scale_shape_manual(values=c(3, 16, 17))+
  geom_text(aes(x=x, y=y, label=txt), data=txt, color="black",
            show.legend = F)+
  scale_color_lancet()+ 
  scale_fill_lancet()+ 
  #scale_x_continuous(limits = c(-4,4), breaks = seq(-4,4,1))+
  #scale_y_continuous(limits = c(138,141), breaks = seq(138, 141, 1))+
  theme_classic()+
  labs(y="Plasma sodium level (mmol/L)", x="Standardized polygenic score"#,
       #tag = "D"
  )+
  theme(legend.position="top",
        legend.title = element_blank()#,
        #legend.text  = element_text(size=10, face="bold")
  ) + 
  coord_cartesian(xlim = c(-4.02,4.02), expand=F)
px


## outputting intaraction table ------
tab
write.table(tab,"THI_RAS_CCB_interaction_models.txt",col.names=T,row.names=F,quote=F,sep="\t")

## Percentile plot of on + pre sodium for thiazide users by PGS -------
lr <- summary(lm(on_drug ~ -1 + quantile_all ,data=dat_thi_residuals)) 
lr_effect <- lr$coefficients[,1]
lr_sd <- lr$coefficients[,2]
q975 <- qnorm(0.975)
s <- cbind(c(1:10),lr_effect,lr_sd,lr_effect+q975*lr_sd,lr_effect-q975*lr_sd)
s <- as.data.frame(s)
names(s) <- c("Quantile","Effect","SD","CI.U","CI.L")
s$cohort <- "On-thiazide drug sodium level    "
ss <- s

lr <- summary(lm(pre_drug ~ -1 + quantile_all ,data=dat_thi_residuals)) 
lr_effect <- lr$coefficients[,1]
lr_sd <- lr$coefficients[,2]
q975 <- qnorm(0.975)
s <- cbind(c(1:10),lr_effect,lr_sd,lr_effect+q975*lr_sd,lr_effect-q975*lr_sd)
s <- as.data.frame(s)
names(s) <- c("Quantile","Effect","SD","CI.U","CI.L")
s$cohort <- "Pre-thiazide drug sodium level    "
ss <- rbind(ss,s) 


txt <- paste0("P-value for interaction: ", sprintf(("%.2f"),tab[16,9])) # P-value from ln_diff model, adjusted for sex, age + PCs
txt <- as.data.frame(txt)
txt$x <- 3
txt$y <- 141
txt$drug <- "Pre-thiazide drug sodium level    "

p3 <- ggplot(ss ,aes(x=Quantile,
                     y=Effect,
                     
                     color = factor(cohort)))+
  geom_point()+
  geom_errorbar(width=0.1, aes(ymin=CI.L,
                               ymax=CI.U)) +
  geom_text(aes(x=x, y=y, label=txt), data=txt, color="black",
           show.legend = F)+
  xlab("Standardized PGS Centiles")+
  ylab("Plasma sodium level (mmol/L)")+
  theme_classic() + 
  scale_color_lancet()+
  scale_x_continuous(breaks=seq(1,10,1)) +
  scale_y_continuous(limits=c(138,141.3),breaks=seq(138,141,1)) +
  theme(
    legend.title = element_blank(),
    legend.position = "top") +
  labs(tag="A")
p3

sPGSss <- rbind(sPGS, ss)
write.table(sPGSss,"PGS_effect_quantile.txt",col.names=T,row.names=F,quote=F,sep="\t")

####
## Interaction-plot THI vs RAS for on drug level ----
txt <- paste0("P-value for interaction: ", sprintf(("%.2f"),tab[3,9])) # P-value from ln_diff model, adjusted for sex, age + PCs
txt <- as.data.frame(txt)
txt$x <- -2
txt$y <- 141.5
txt$drug <- "Thiazide users"


p5 <- dat[dat$drug %in% c("THI","RAS"),] %>% 
  mutate(drug=ifelse(drug=="THI", "Thiazide users", "RASi users"),
         drug = factor(drug, levels=c("Thiazide users", "RASi users"))) %>%
  ggplot(aes(x=zscore,y=on_drug,color=drug, shape=drug)) + 
  # geom_point(size=0.5) +
  geom_smooth(method=lm, aes(fill=drug),
              se=T, fullrange=TRUE,
              linewidth=0.6,
              alpha=0.4)+
  #scale_shape_manual(values=c(3, 16, 17))+
  geom_text(aes(x=x, y=y, label=txt), data=txt, color="black",
            show.legend = F)+
  scale_color_lancet()+ 
  scale_fill_lancet()+ 
  #scale_x_continuous(limits = c(-4,4), breaks = seq(-4,4,1))+
  #scale_y_continuous(limits = c(138,141), breaks = seq(138, 141, 1))+
  theme_classic()+
  labs(y="On-drug plasma sodium concentration (mmol/L)", x="Standardized polygenic score"#,
      # tag = "B"
      )+
  theme(legend.position="top",
        legend.title = element_blank()#,
        #legend.text  = element_text(size=10, face="bold")
  ) + 
  coord_cartesian(xlim = c(-4.02,4.02), expand=F)
p5

## Outputting Fig 1 template ------
p <- grid.arrange(grobs=list(p1, p2, p3, p5), nrow=4, ncol=1, height=c(1,1,1,1))

p <- grid.arrange(grobs=list(p1, p2, p3, p5), nrow=2, ncol=2, height=c(1,1,1,1))

p <- grid.arrange(grobs=list(p3, p5), nrow=2, ncol=1, height=c(1,1,))
px
ggsave("figs_pgs1.png", 
       p, device = "png", dpi = 300,
       width = 7, height = 10, units="cm", scale=2.5)

ggsave("p5_ras_thi_int.png", 
       p5, device = "png", dpi = 300,
       width = 14, height = 10, units="cm", scale=1.2)
ggsave("p5_ras_thi_int.pdf", 
       p5, device = "pdf", dpi = 300,
       width = 14, height = 10, units="cm", scale=1.2)
ggsave("p5_ras_thi_int.tiff", 
       p5, device = "tiff", dpi = 300,
       width = 14, height = 10, units="cm", scale=1.2)

ggsave("px_thi_pre_on_int.png", 
       px, device = "png", dpi = 300,
       width = 14, height = 10, units="cm", scale=1.2)
## Comparing risk of hyponatremia for THI vs RAS in PGS (4) quantiles -----------------
### prepping data
tjex <- dat %>% 
  filter(drug != "CCB",
         !is.na(drug)
  ) %>%
  mutate(
    drug = factor(drug),
    drug = relevel(drug, ref = "RAS"),
    hypona_130 = ifelse(on_drug <= 130, 1, 0),
    hypona_134 = ifelse(on_drug <= 134, 1, 0),
    hypona_125 = ifelse(on_drug <= 125, 1, 0),
    zscore_decile = factor(ntile(zscore, 5)),
    zscore_decile = relevel(zscore_decile, ref = "5"),
    diff = on_drug - pre_drug,
    diff_decile = factor(ntile(-diff, 3)),
    diff_decile = relevel(diff_decile, ref = "3"), 
    diff_lowest_decile = ifelse(diff_decile== "1", 1, 0),
    
    diff_5 = ifelse(diff <=-10, 1,0),
    zscore_rev = -zscore )

tjex$quantile_all <- get_quantile(tjex$zscore,num.quant=4, quant.ref=4)

PGS$quantile_tot <- get_quantile(PGS$zscore,num.quant=4, quant.ref=1)

tjex <- tjex %>%
  left_join(PGS %>% select(c("cpr_enc", "quantile_tot")), by = c("cpr_enc"="cpr_enc")) # Using the overall PGS distribution for quantile categorization

### Bar plot for hyponatremia for THI vs RAS by PGS ----
dat_bar <- tjex %>% group_by(drug) %>% 
  count(quantile_tot, hypona_130) %>%
  pivot_wider(names_from = c( "hypona_130"), names_prefix="q", 
              values_from= c("n")) %>% 
  mutate(sum = q0+q1,
         prop = q1/sum)

dat_bar <- dat_bar %>% mutate(quantile_tot = case_when(quantile_tot == "1" ~ 1,
                                                       quantile_tot == "2" ~ 2,
                                                       quantile_tot == "3" ~ 3,
                                                       quantile_tot == "4" ~ 4,
                                                       T ~ as.numeric(quantile_tot)),
                              drug2=drug,
                              drug = ifelse(drug=="THI", "Thiazide drug users", "RASi drug users"),
                              drug = factor(drug),
                              drug = relevel(drug, ref = "Thiazide drug users"),
)
p6 <- dat_bar %>% ggplot()+
  geom_bar(aes(x=as.factor(quantile_tot), y=prop, fill=drug), 
           stat="identity", position = position_dodge(width = .53), width=0.5, alpha=0.6) +
  scale_y_continuous(labels=scales::percent_format(scale=100, drop0trailing=T))+
  scale_fill_lancet()+
  theme_classic()+
  labs(x="Standardized PGS quantiles",
       y="Proportion of hyponatremia \n (<130 mmol/L) events",
       tag = "A"
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top")
p6
## OR plot for hyponatremia for THI vs RAS by PGS --------
m1 <- glm(hypona_130 ~  quantile_tot + drug:quantile_tot + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
          family = "binomial",
          data= tjex %>% 
            mutate(quantile_tot = relevel(quantile_tot, ref = "4")))
# m1 <- glm(hypona_130 ~  drug + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
#           family = "binomial",
#           data= tjex %>% filter(quantile_all=="4")
# )
coef_table <- data.frame(summary(m1)$coefficients[, c("Estimate", "Std. Error")])
coefficient_names <- rownames(coef_table)

m1_df <- data.frame(coefficient_names,  coef_table) 
m1_df$or <- exp(m1_df$Estimate)
m1_df$lci <- m1_df$Estimate - 1.96*m1_df$Std..Error
m1_df$lci<- exp(m1_df$lci)
m1_df$uci <- m1_df$Estimate + 1.96*m1_df$Std..Error
m1_df$uci <- exp(m1_df$uci)

m1_df <- m1_df[c(12:15), c(1,4:6)]
m1_df <- m1_df %>% mutate(quantile = case_when(str_detect(coefficient_names, "1") ~ 1,
                                               str_detect(coefficient_names, "2") ~ 2,
                                               str_detect(coefficient_names, "3") ~ 3,
                                               str_detect(coefficient_names, "4") ~ 4))
m1_df <- m1_df[,c(5,2:4)] %>% arrange(quantile)
m1_df$drug <- "THI"

p7 <- m1_df %>% ggplot(aes(x=quantile, 
                           y=or,
                           ymin=lci,
                           ymax=uci,
                           color = factor(drug)))+
  geom_point(aes(x=quantile, 
                 y=or),size = 2) +
  geom_errorbar(aes(x=quantile, 
                    ymin=lci,
                    ymax=uci), width = .05)+
  coord_cartesian(ylim=c(.9,7.4)) + 
  geom_hline(yintercept = 1, color="grey") + 
  scale_y_continuous(trans = "log", breaks = seq(1,7,1), expand=c(0,0))+
  
  scale_color_manual(values = c("#034788"))+
  scale_fill_manual(values = "#034788")+
  theme_classic()+
  labs(x="Standardized PGS quantiles",
       y="Hyponatremia (<130 mmol/L) \n OR (95% CI)",
       tag = "B") + 
  theme(legend.position = "none")

p7

### If wanting to combine in 1 figure ------
test <- dat_bar %>% rename("quantile" = quantile_tot) %>%
  mutate(prop = prop*100) %>%
  left_join(m1_df, by=c("drug2"="drug", "quantile"))

pcomb <- test %>%
  ggplot()+
  geom_bar(aes(x=as.factor(quantile), y=prop, fill=drug), 
           stat="identity", position = position_dodge(width = .6), width=0.5, alpha=0.6) +
  geom_point(aes(x=quantile, 
                 y=or),size = 2) +
  geom_errorbar(aes(x=quantile, 
                    ymin=lci,
                    ymax=uci), width = .03)+
  geom_hline(yintercept = 1, color="grey") + 
  scale_y_continuous(name = "Proportion of hyponatremia \n (<130 mmol/L) events (%)", 
                     breaks=seq(0,6,1),
                     sec.axis = sec_axis(~.*1, name = "OR (95% CI)",
                                         breaks=seq(1,6,1)))+
  scale_fill_lancet()+
  theme_classic()+
  labs(x="Standardized PGS quantiles"
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top")

pcomb

ggsave("igs_pgs2_pcomb.png", 
       pcomb, device = "png", dpi = 300,
       width = 7, height = 4, units="cm", scale=2.7)

## Comparing risk of hyponatremia between THI in PGS (4) quantiles ----------------
#TODO (consider modeling including RAS and without RAS, the latter is the current)

### Adjusted OR plot for hyponatremia between THI in PGS quantiles ------
thi_dat <- dat %>% 
  filter(drug == "THI") %>%
  mutate(
    drug = factor(drug),
    hypona_130 = ifelse(on_drug <= 130, 1, 0),
    hypona_134 = ifelse(on_drug <= 134, 1, 0),
    hypona_125 = ifelse(on_drug <= 125, 1, 0),
    zscore_decile = factor(ntile(zscore, 5)),
    zscore_decile = relevel(zscore_decile, ref = "5"),
    diff = on_drug - pre_drug,
    diff_decile = factor(ntile(-diff, 3)),
    diff_decile = relevel(diff_decile, ref = "3"), 
    diff_lowest_decile = ifelse(diff_decile== "1", 1, 0),
    diff_5 = ifelse(diff <=-10, 1,0),
    zscore_rev = -zscore )

PGS$quantile_tot <- get_quantile(PGS$zscore,num.quant=4, quant.ref=4)

thi_dat <- thi_dat %>%
  left_join(PGS %>% select(c("cpr_enc", "quantile_tot")), by = c("cpr_enc"="cpr_enc"))

m130 <- glm(hypona_130 ~ quantile_tot + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
            family = "binomial",
            data= thi_dat)

coef_table <- data.frame(summary(m130)$coefficients[, c("Estimate", "Std. Error")])
coefficient_names <- rownames(coef_table)

m1_df <- data.frame(coefficient_names,  coef_table) 
m1_df$or <- exp(m1_df$Estimate)
m1_df$lci <- m1_df$Estimate - 1.96*m1_df$Std..Error
m1_df$lci<- exp(m1_df$lci)
m1_df$uci <- m1_df$Estimate + 1.96*m1_df$Std..Error
m1_df$uci <- exp(m1_df$uci)

m1_df <- m1_df[c(2:4), c(1,4:6)]
m1_df <- m1_df %>% mutate(quantile = case_when(str_detect(coefficient_names, "1") ~ 1,
                                               str_detect(coefficient_names, "2") ~ 2,
                                               str_detect(coefficient_names, "3") ~ 3,
                                               str_detect(coefficient_names, "4") ~ 4))
m1_df <- m1_df[,c(5,2:4)] %>% arrange(quantile)
m1_df$hyponatremia <- "\u2264 130 mmol/L"
s <- m1_df

m125 <- glm(hypona_125 ~ quantile_tot + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
            family = "binomial",
            data= thi_dat)


coef_table <- data.frame(summary(m125)$coefficients[, c("Estimate", "Std. Error")])
coefficient_names <- rownames(coef_table)

m1_df <- data.frame(coefficient_names,  coef_table) 
m1_df$or <- exp(m1_df$Estimate)
m1_df$lci <- m1_df$Estimate - 1.96*m1_df$Std..Error
m1_df$lci<- exp(m1_df$lci)
m1_df$uci <- m1_df$Estimate + 1.96*m1_df$Std..Error
m1_df$uci <- exp(m1_df$uci)

m1_df <- m1_df[c(2:4), c(1,4:6)]
m1_df <- m1_df %>% mutate(quantile = case_when(str_detect(coefficient_names, "1") ~ 1,
                                               str_detect(coefficient_names, "2") ~ 2,
                                               str_detect(coefficient_names, "3") ~ 3,
                                               str_detect(coefficient_names, "4") ~ 4))
m1_df <- m1_df[,c(5,2:4)] %>% arrange(quantile)
m1_df$hyponatremia <- "\u2264 125 mmol/L"

s <- rbind(s, m1_df)

m134 <- glm(hypona_134 ~ quantile_tot + sex + thiazide_age + PC1+PC2+PC3+PC4+PC5,
            family = "binomial",
            data= thi_dat)

coef_table <- data.frame(summary(m134)$coefficients[, c("Estimate", "Std. Error")])
coefficient_names <- rownames(coef_table)

m1_df <- data.frame(coefficient_names,  coef_table) 
m1_df$or <- exp(m1_df$Estimate)
m1_df$lci <- m1_df$Estimate - 1.96*m1_df$Std..Error
m1_df$lci<- exp(m1_df$lci)
m1_df$uci <- m1_df$Estimate + 1.96*m1_df$Std..Error
m1_df$uci <- exp(m1_df$uci)

m1_df <- m1_df[c(2:4), c(1,4:6)]
m1_df <- m1_df %>% mutate(quantile = case_when(str_detect(coefficient_names, "1") ~ 1,
                                               str_detect(coefficient_names, "2") ~ 2,
                                               str_detect(coefficient_names, "3") ~ 3,
                                               str_detect(coefficient_names, "4") ~ 4))
m1_df <- m1_df[,c(5,2:4)] %>% arrange(quantile)
m1_df$hyponatremia <- "< 135 mmol/L"
s <- rbind(s, m1_df)

s

x <- data.frame(quantile = 4, or = 1.00, lci = NA, uci=NA, hyponatremia="Reference")

ss <- rbind(s,x) %>% 
  mutate(hyponatremia = factor(hyponatremia, levels=c("< 135 mmol/L" , "\u2264 130 mmol/L", "\u2264 125 mmol/L", "Reference"))) 

p9 <- ss %>% ggplot(aes(x=quantile, 
                        y=or,
                        ymin=lci,
                        ymax=uci,
                        color = factor(hyponatremia)))+
  geom_point(aes(x=quantile, 
                 y=or),size = 3, 
             position = position_dodge(width = .5)) +
  geom_errorbar(aes(x=quantile, 
                    ymin=lci,
                    ymax=uci), width = .05,
                position = position_dodge(width = .5))+
  coord_cartesian(ylim=c(.3,4)) + 
  geom_hline(yintercept = 1, color="grey") + 
  scale_y_continuous(trans = "log", breaks = seq(.5,4,.5), expand=c(0,0))+
  scale_x_continuous(breaks = seq(1,4,1), expand=c(0,0), labels=c("<25%", "25%-<50%", "50%-<75%", "\u2265 75%"))+
  #scale_x_discrete() +  
  scale_color_lancet(alpha=0.8) +
  scale_fill_lancet() + 
  theme_classic()+
  labs(x="Standardized polygenic score, quantiles",
       y="Risk of hyponatremia, \n OR (95% CI)",
       tag = "B") + 
  theme(legend.position = "top",
        legend.title = element_blank())

p9

## Bar plot for hyponatremia between THI in PGS quantiles  ------
prop_130 <- thi_dat %>% group_by(drug) %>% 
  count(quantile_tot, hypona_130) %>%
  pivot_wider(names_from = c( "hypona_130"), names_prefix="q", 
              values_from= c("n")) %>% 
  mutate(sum = q0+q1,
         prop_130 = q1/sum) %>%
  select(drug, quantile = quantile_tot, prop_130)


prop_125 <- thi_dat %>% group_by(drug) %>% 
  count(quantile_tot, hypona_125) %>%
  pivot_wider(names_from = c( "hypona_125"), names_prefix="q", 
              values_from= c("n")) %>% 
  mutate(sum = q0+q1,
         prop_125 = q1/sum) %>%
  select(drug, quantile = quantile_tot, prop_125)

prop_134 <- thi_dat %>% group_by(drug) %>% 
  count(quantile_tot, hypona_134) %>%
  pivot_wider(names_from = c( "hypona_134"), names_prefix="q", 
              values_from= c("n")) %>% 
  mutate(sum = q0+q1,
         prop_134 = q1/sum) %>%
  select(drug, quantile = quantile_tot, prop_134)

dat_bar <- prop_130 %>%
  left_join(prop_125, by=c("drug","quantile")) %>%
  left_join(prop_134, by=c("drug","quantile")) %>%
  pivot_longer(cols=c(prop_130, prop_125, prop_134), names_to="hyponatremia", values_to="prop") %>%
  mutate(hyponatremia = case_when(hyponatremia == "prop_130" ~ "\u2264 130 mmol/L",
                                  hyponatremia == "prop_125" ~ "\u2264 125 mmol/L",
                                  hyponatremia == "prop_134" ~ "< 135 mmol/L"),
         hyponatremia = factor(hyponatremia, levels=c("< 135 mmol/L", "\u2264 130 mmol/L", "\u2264 125 mmol/L")))

dat_bar <- dat_bar %>% mutate(quantile = case_when(quantile == "1" ~ 1,
                                                   quantile == "2" ~ 2,
                                                   quantile == "3" ~ 3,
                                                   quantile == "4" ~ 4,
                                                   T ~ as.numeric(quantile))
                              
)
p8 <- dat_bar %>% ggplot()+
  geom_bar(aes(x=as.factor(quantile), y=prop, fill=hyponatremia), 
           stat="identity", position = position_dodge(width = .53), width=0.5, alpha=0.6) +
  scale_y_continuous(labels=scales::percent_format(scale=100, drop0trailing=T), breaks=seq(0,.10,.02))+
  scale_x_discrete(labels=c("<25%", "25%-<50%", "50%-<75%", "\u2265 75%")) +
  scale_fill_lancet()+
  theme_classic()+
  labs(x="Standardized polygenic score, quantiles",
       y="Proportion of hyponatremia events (%)",
       tag = "A"
  ) + 
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.key.size = unit(0.5, "cm"))
p8

pxx <- grid.arrange(grobs=list(p6,p7,p8, p9), nrow=2, ncol=2, height=c(1,1,1,1))

pxx <- grid.arrange(grobs=list(p8, p9), nrow=2, ncol=1, height=c(1,1))

ggsave("figs_pgs2.png", 
       pxx, device = "png", dpi = 300,
       width = 7, height = 8, units="cm", scale=2.7)

ortab <- ss %>% left_join(dat_bar %>% select(-drug))

write.table(ortab, "tab_pgs_orTIH.txt", col.names = T, row.names = F, quote = F)

## END -------

