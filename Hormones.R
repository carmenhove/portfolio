
#PROJECT 1: What are the sex-specific effects of age on hormone production?

#Identify packages 
pkgs = c("nhanesA","tidyverse","magrittr","knitr","patchwork","ggpubr","scales","brms","mgcv","schoenberg","tidymv")

#Install packages
install.packages(pkgs)

#Load packages
inst = lapply(pkgs, library, character.only = TRUE) 

#Install and load Wes Anderson color palette from devtools
devtools::install_github("karthik/wesanderson")
library(wesanderson)

#Load data directly from NHANES
##Identify BMI, demographic, reproductive health, and hormone data files
table.names <- sort(c('DEMO', 'RHQ','BMX', 
                      str_c(c('DEMO','RHQ','BMX'),"_",rep(c('H','I'),each=3)), 
                      'TST_H','TST_I'))

#apply nhanes() to extract SAS data frames, transition to R list
nhanes.l1 <- sapply(table.names, nhanes) 

#remove SAS labels
nhanes.l2 <- sapply(nhanes.l1, remove_all_labels) 

#reduce and merge data files by file type
BMI <- Reduce(full_join,nhanes.l2[grepl("BMX", names(nhanes.l2))])
RHQ <- Reduce(full_join,nhanes.l2[grepl("RHQ", names(nhanes.l2))])
DEMO <- Reduce(full_join,nhanes.l2[grepl("DEMO", names(nhanes.l2))])
HORM <- Reduce(full_join,nhanes.l2[grepl("TST", names(nhanes.l2))])

#Merge, clean dataset
df1 <- full_join(BMI, DEMO) %>% 
  full_join(., RHQ) %>% full_join(., HORM) %>% 
  rename(Sex = RIAGENDR, Age = RIDAGEYR, BMI = BMXBMI, 
         Preg.demo = RIDEXPRG, Preg.rhq = RHD143, #Currently pregnant?
         Estradiol = LBXEST, Testosterone = LBXTST,
         Est.cc = LBDESTLC, Test.cc = LBDTSTLC,SHBG = LBXSHBG) %>% 
  mutate(across(c(Estradiol, Testosterone),abs), #Make sure E and T are absolute values
         Pregnant = case_when(Preg.demo == 1 | Preg.rhq == 1 ~ "Yes",TRUE ~ "No"),
         Sex = as.factor(case_when(Sex == 1 ~ "Male", Sex == 2 ~ "Female")),
         Pregnant.Sex = case_when(Sex == "Male" ~ "Male",
                                  Sex == "Female" & Pregnant == "No" ~ "NP female",
                                  Sex == "Female" & Pregnant == "Yes" ~ "P female")) %>% 
  filter(across(c(Testosterone, Estradiol,Age,Sex), ~ !is.na(.)) &
           (RHD280 == 2 | (Sex == "Female" & Age < 20) | Sex == "Male") & #No hysterectomy (RHD280 ==2), Not measured in any females under the age of 20
           #Remove imputed values due to below-threshold detection limit reported
           Est.cc == 0 & Test.cc == 0) %>% 
  select(Testosterone,Estradiol,SHBG,Pregnant,Sex,Age,BMI,Pregnant.Sex)

df2 <- df1 %>% 
  pivot_longer(Testosterone:SHBG, 
               names_to = "Measure",
               values_to = "Value")

#Structure of dataset
str(df1)

#Check distributions
ggplot(df2, aes(x = log(Value), color = Sex, fill = Sex))+
  geom_histogram(bins = 30, alpha = 0.8)+
  facet_grid(~ Measure,scales = "free")

ggplot(df2, aes(x = Age, y= log(Value),color = Pregnant.Sex))+
  geom_point(position = "jitter") + 
  facet_grid(Measure ~., scales = "free")

ggplot(df1, aes(x = log(Estradiol), 
                y= log(Testosterone),color = Pregnant.Sex))+
  geom_point(position = "jitter") #+ 
  #facet_grid(Measure ~., scales = "free")

#Number of pregnant individuals in sample
table(df1[df1$Sex=="Female",]$Pregnant)

ggplot(df1, aes(x = Estradiol, color = Sex))+
  geom_histogram(bins=45)+
  facet_grid(~ Sex,scales = "free")

#Impact of pregnant individuals on sample
pregplot <- ggplot(df1, aes(x = Age, y = log(Estradiol), 
                                color = Pregnant))+
  geom_point(position = "jitter")

#Weird low values for male Testosterone
ggplot(dataset, aes(x = log(Testosterone), color = Sex))+
  geom_histogram()#+
  #facet_grid(Sex ~.,scales = "free")

#Weird low values for both male and female estradiol

ggplot(df1, aes(x = BMI, y = log(SHBG), color = Sex))+
  geom_point(position = "jitter")

ggplot(dataset, aes(x = Age, y = log(Testosterone), color = Sex))+
  geom_point(position = "jitter")

#Fully Bayesian model using brms
est.modbrms <- brm(formula = log(Estradiol) ~ s(Age, by = Sex),
               data =dataset, family = gaussian(),
               warmup=1000, iter=2500,
               chains=2,cores=2)
summary(est.modbrms)

#glm model
est.modglm <- mgcv::gam(log(Estradiol) ~ s(Age, by = Sex), 
                  data = dataset, method = "REML")
summary(est.modglm)

est.modglm2 <- mgcv::gam(log(Estradiol) ~ s(Age, by = Sex)+BMI*Sex, 
                        data = dataset, method = "REML")
summary(est.modglm)

test <- lm(log(Estradiol) ~ log(Testosterone)*Sex,
   data = df1)
summary(test)

shbg.modglm2 <- mgcv::gam(log(SHBG) ~ s(Age, by = Sex),
                         data = dataset, method = "REML")
summary(shbg.modglm)

install.packages("tidymv")
library(tidymv)
plot_smooths(est.modglm, Age, Sex)
plot_smooths(shbg.modglm2, Age, Sex)
#Plot raw data and model output onto single graph, to show I got the goods






