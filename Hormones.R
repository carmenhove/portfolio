
#PROJECT A

#Identify packages 
pkgs = c("nhanesA","tidyverse","magrittr","knitr",
         "patchwork","ggpubr","scales","brms","mgcv",
        "sjlabelled","schoenberg","tidymv","broom",
        "broom.mixed","tidybayes","calecopal")

#Install packages
install.packages(pkgs)

#Load packages
inst = lapply(pkgs, library, character.only = TRUE) 

#Install and load Wes Anderson color palette from devtools
devtools::install_github("karthik/wesanderson")
library(wesanderson)

#Rstan has been acting squirrely on my OS when I tweak adapt_delta and max treedepth. To get around this, I use cmdstranr on the back-end. 
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2,overwrite = T)
cmdstan_path()
cmdstan_version()
options(brms.backend = "cmdstanr")


#Load data from NHANES
##Identify BMI, demographic, reproductive health, and hormone data files
table.names <- sort(c(str_c(c('DEMO','RHQ','BMX'),"_",rep(c('H','I'),each=3)), 
                      'TST_H','TST_I'))

#apply nhanes() to extract SAS data frames, transition to R list
nhanes.l1 <- map(table.names, nhanes) 

#remove SAS labels
nhanes.l2 <- map(nhanes.l1, remove_all_labels)

#Name list elements
names(nhanes.l2)<-table.names

#Reduce and merge data files by file type
BMI <- Reduce(full_join,nhanes.l2[grepl("BMX", names(nhanes.l2))])
RHQ <- Reduce(full_join,nhanes.l2[grepl("RHQ", names(nhanes.l2))])
DEMO <- Reduce(full_join,nhanes.l2[grepl("DEMO", names(nhanes.l2))])
HORM <- Reduce(full_join,nhanes.l2[grepl("TST", names(nhanes.l2))])

#Merge, clean dataset
PAdf1 <- full_join(BMI, DEMO) %>% 
  full_join(., RHQ) %>% full_join(., HORM) %>% 
  rename(Sex = RIAGENDR, Age = RIDAGEYR, BMI = BMXBMI, 
         Preg.demo = RIDEXPRG, Preg.rhq = RHD143, #Currently pregnant?
         Estradiol = LBXEST, Testosterone = LBXTST,
         Est.cc = LBDESTLC, Test.cc = LBDTSTLC,
         SHBG = LBXSHBG, Hysterectomy = RHD280) %>% 
  mutate(across(c(Estradiol, Testosterone),abs), #Make sure E and T are absolute values
         Pregnant = case_when(Preg.demo == 1 | Preg.rhq == 1 ~ "Yes",TRUE ~ "No"),
         Sex = as.factor(case_when(Sex == 1 ~ "Male", Sex == 2 ~ "Female")),
         Pregnant.Sex = case_when(Sex == "Male" ~ "Male",
                                  Sex == "Female" & Pregnant == "No" ~ "NP female",
                                  Sex == "Female" & Pregnant == "Yes" ~ "P female")) %>% 
  filter(across(c(Testosterone, Estradiol,Age,Sex), ~ !is.na(.)) & #Remove NAs
           #No hysterectomy (not measured in any females under the age of 20)
           (Hysterectomy == 2 | (Sex == "Female" & Age < 20) | Sex == "Male")) %>%  
           #Remove imputed values due to below-threshold detection limit reported
           #Est.cc == 0 & Test.cc == 0) %>% 
  select(Testosterone,Test.cc,Estradiol,Est.cc,Pregnant,Pregnant.Sex,Age,Sex,BMI) %>% 
  mutate(Sex = ordered(Sex, levels = c("Male","Female")),
         Pregnant.Sex = ordered(Pregnant.Sex, levels = c("Male","NP Female","P Female")))

#Structure of dataset
str(PAdf1)

#Pivot longer to create "Measure" variable
PAdf2 <- PAdf1 %>% 
  pivot_longer(c(Testosterone,Estradiol), 
               names_to = "Measure",
               values_to = "Value") %>% 
  filter(!is.na(Value))

str(PAdf2)

#Histogram of outcome variables
ggplot(PAdf2, aes(x = Value, color = Sex))+
  geom_histogram(bins = 30, alpha = 0.8)+
  facet_grid(~ Measure, scales = "free")

#No zero values
sum(PAdf2$Value==0,na.rm = T)

#Log-transformation
ggplot(PAdf2, aes(x = log(Value), color = Sex))+
  geom_histogram(bins = 30, alpha = 0.8)+
  facet_grid(~ Measure, scales = "free")

#Remove imputed values
PAdf3 <- PAdf2 %>%  
  filter(Est.cc == 0 & Test.cc == 0)

#Check distributions
ggplot(PAdf3, aes(x = log(Value), color = Sex))+
  geom_histogram(bins = 30, alpha = 0.8)+
  facet_grid( ~ Measure,scales = "free")

#What about pregnancy?
table(PAdf3$Pregnant) #112 pregnant individuals in this dataset

#Check plots to see how pregnancy relates to Estradiol
ggplot(PAdf3, aes(x = Age, y = log(Value), color = Pregnant.Sex))+
  geom_point()+
  facet_grid(~ Measure, scales = "free")+
  scale_color_manual(values = cal_palette("dudleya"))

#Make the df a list, split by measure and sex (given sex-specific distributions)
PAlist <- split(PAdf3, list(PAdf3$Measure, PAdf3$Sex))

#A function for fully Bayesian models using brms, with cmdstanr on the backend
get.brms <- function(x){
  
  if(unique(x$Measure) == "Testosterone"){
    knot.values = c(6,12,18,20,40,70)
  } else {knot.values = c(6,12,18,20,30,40,50,70)
    }
  model.prior <- get_prior(log(Value) ~ 
                             s(Age, k = length(knot.values)), 
                           knots = list(Age = knot.values),
                           data = x, family = gaussian())
  model.prior$prior[1]<-"normal(0,10"
                           
  model <- brms::brm(formula = log(Value) ~ 
                       s(Age, k = length(knot.values)), 
               knots = list(Age = knot.values),
               data = x, family = gaussian(),
               warmup=1000, iter=2500,prior=model.prior,
               control = list(max_treedepth = 12),
               chains=2, cores=2, seed = 1234)
  model$Measure = unique(x$Measure)
  model$Sex = unique(x$Sex)
  model
}

#Apply brms function to list
brms.models <- map(PAlist,get.brms)

#Create function to gather tidy predicted draws from each model
getbrmspreds <- function(x){
  #Age values from original datasets
  new_data <- data.frame(Age = sort(unique(PAdf3$Age)))
  
  preds <- as_tibble(x %>% fitted(newdata = new_data,
                                  summary = T)) %>% 
    rename(lnValue = Estimate) %>% 
    mutate(Age = sort(unique(PAdf3$Age)), Sex = x$Sex, 
           Measure = x$Measure, Value = exp(lnValue))
  preds
}

#Apply function to list of brms models
brms.predvals <- plyr::ldply(map(brms.models, getbrmspreds),
                             data.frame, .id = NULL) %>% 
  mutate(Model.Type = "BRMS")

#Get gam models using mgcv 
get.gams <- function(x){
  model <- mgcv::gam(log(Value) ~ s(Age), 
                     data = x, method = "REML")
  model$Measure = unique(x$Measure)
  model$Sex = unique(x$Sex)
  model
}

gam.models <- map(PAlist, get.gams)

getgampreds <- function(x){
  preds <- predict_gam(x) %>% 
    rename(lnValue = fit) %>% 
    mutate(Measure = x$Measure,
           Sex= x$Sex)
  preds
}

gam.predvalues <- plyr::ldply(map(gam.models, getgampreds),
                              data.frame, .id = NULL) %>% 
  mutate(Model.Type = "GAM")

str(gam.predvalues)
str(brms.predvals)

colorset1 = c('Male'='#A2A098','Female'='#5E6B7B')
colorset2 = c('Male'='#A2A098','NP female'='#5E6B7B','P female'='#233D3F')
  
ggplot(brms.predvals, aes(x = Age, y = lnValue, color =Sex))+
  geom_point(data = PAdf3, aes(x = Age, y = log(Value),
                               color = Sex))+
  geom_line(size = 0.7)+
  geom_ribbon(aes(ymin = Q2.5, ymax =Q97.5, group = Sex),
              alpha=0.6, linetype = 0) +
  facet_grid(~Measure,scales = "free")+
  scale_color_manual(values = colorset1)+
  labs(title = "predicted values using brms")

ggplot(gam.predvalues, aes(x = Age, y = lnValue, 
                           color = Sex)) + 
  geom_line(size = 0.7)+
  geom_point(data = PAdf3, aes(x = Age, y = log(Value),
                               color = Sex), position = "jitter")+
  geom_ribbon(aes(ymin = lnValue - se.fit, 
                  ymax = lnValue + se.fit, 
                  group = Sex),
              alpha=0.6, linetype = 0,,show.legend = FALSE) +
  facet_grid(~ Measure,scales = "free")+
  scale_color_manual(values = cal_palette("dudleya"))+
  scale_fill_manual(values = cal_palette("dudleya"))+
  labs(title = "predicted values using gam")+
  theme(legend.position = "bottom")

