install.packages("visreg")
library(visreg)
#PROJECT A

#Identify packages 
pkgs = c("nhanesA","tidyverse","magrittr","knitr",
         "patchwork","ggpubr","scales","brms","mgcv",
        "sjlabelled","schoenberg","tidymv","broom",
        "broom.mixed","tidybayes","calecopal","smotefamily",
        "corrplot","car","e1071","InformationValue","class")

#Install packages
install.packages(pkgs)

#Load packages
inst = lapply(pkgs, library, character.only = TRUE) 

#Install and load Wes Anderson color palette from devtools
devtools::install_github("karthik/wesanderson")
library(wesanderson)

colorset1 = c('Male'='#A2A098','Female'='#5E6B7B')
colorset2 = c('Male'='#A2A098','NP female'='#5E6B7B','P female'='#233D3F')

#Rstan has been acting squirrely on my OS when I tweak adapt_delta and max treedepth. To get around this, I use cmdstranr on the back-end. 
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", 
                                       getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain()
install_cmdstan(cores = 2,overwrite = T)
cmdstan_path()
cmdstan_version()
options(brms.backend = "cmdstanr")


#Load data from NHANES
##Identify BMI, demographic, reproductive health, and hormone data files
table.names <- sort(c(str_c(c('DEMO','RHQ','BMX'),"_",
                            rep(c('H','I'),each=3)), 
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
  select(Testosterone,Test.cc,Estradiol,Est.cc,Pregnant,
         Pregnant.Sex,Age,Sex,BMI,) %>% 
  mutate(#Sex = ordered(Sex, levels = c("Male","Female")),
         Pregnant.Sex = ordered(Pregnant.Sex, levels = c("Male","NP female","P female")),
         AgeGroup = case_when(Age < 20 ~ "Under 20",
                              Age >= 20 & Age <=45 ~ "20-45",
                              Age > 45 ~ "45+"),
         AgeGroup = ordered(AgeGroup, levels = c("Under 20",
                            "20-45","45+")))

#Pivot longer to create "Measure" variable
PAdf2 <- PAdf1 %>% 
  pivot_longer(c(Testosterone,Estradiol), 
               names_to = "Measure",
               values_to = "Value") %>% 
  filter(!is.na(Value))

#Remove imputed values
PAdf3 <- PAdf2 %>%  
  filter(Est.cc == 0 & Test.cc == 0)

#Remove pregnant individuals
PAdf4 <- PAdf3 %>% filter(Pregnant == "No")

ggplot(PAdf4,#%>% filter(Measure == "Testosterone"), 
       aes(x = Value, color = Sex, fill = Sex))+
       geom_density(alpha = 0.5)+ 
       facet_grid( ~ Measure, scales = "free")+
       scale_color_manual(values = colorset1)+
       scale_fill_manual(values = colorset1) 

ggplot(PAdf4, aes(x = Age, y = log(Value), color = Sex))+
  geom_point(position = "jitter")+
  facet_grid(Measure ~., scales = "free")+
  scale_color_manual(values = colorset1)+
  scale_fill_manual(values = colorset1)



#Make the df a list, split by measure and sex (given sex-specific distributions)
PAlist <- split(PAdf4, list(PAdf4$Measure,PAdf4$Sex))

##COMPARE LINEAR AND NON-PARAMETRIC MODELS IN GAM (to show that non-linear model is best)
#Get glm models
get.glms <- function(x){
  model <- glm(log(Value) ~ Age+BMI,data = x)
  model$Measure = unique(x$Measure)
  model$Sex = unique(x$Sex)
  model
}

glm.models <- map(PAlist, get.glms)
system.time(result <- myfunction(with, arguments))

#Get gam models using mgcv 
get.gams <- function(x){
  if(unique(x$Measure) == "Testosterone"){
    knot.values = c(6,12,18,20,40,70)
  } else {knot.values = c(6,12,18,20,30,40,50,70)
  }
  model <- mgcv::gam(log(Value) ~ s(Age)+BMI,
                     k =length(knot.values),
                     knots = list(x = knot.values),
                     data = x)
  model$Measure = unique(x$Measure)
  model$Sex = unique(x$Sex)
  model
}

gam.models <- map(PAlist, get.gams)
gam.summaries <- map(gam.models, summary)

#Example model comparison showing that smooth model is more appropriate
anova(glm.models[["Estradiol.Male"]],
      gam.models[["Estradiol.Male"]],test="Chisq")

#Plot showing major heteroscedasticity in glm models
plot(glm.models[["Estradiol.Male"]]$fitted.values, 
     glm.models[["Estradiol.Male"]]$residuals)

#Get gam model predicted values
getgampreds <- function(x){
  preds <- predict_gam(x,values = list(Age = unique(PAdf3$Age),
                                       BMI = median(PAdf3$BMI,
                                                    na.rm = T))) %>% 
    rename(lnValue = fit) %>% 
    mutate(Measure = x$Measure,
           Sex= x$Sex)
  preds
}

#Put predicted values in a dataframe
gam.predvals <- plyr::ldply(map(gam.models, getgampreds),
                            data.frame, .id = NULL) %>% 
  mutate(Model.Type = "GAM")

#Check the structure
str(gam.predvals)

#Plot gam predicted values
ggplot(gam.predvals, aes(x = Age, y = lnValue, 
                           color = Sex, fill = Sex)) + 
  geom_line(size = 0.7)+
  #geom_point(data = PAdf3, aes(x = Age, y = log(Value),
  #                             color = Sex), 
  #           position = "jitter", alpha = 0.3)+
  geom_ribbon(aes(ymin = lnValue - se.fit, 
                  ymax = lnValue + se.fit, 
                  group = Sex),
              alpha=0.6, linetype = 0) +
  facet_grid(Measure ~.,scales = "free")+
  scale_color_manual(values = colorset1)+
  scale_fill_manual(values = colorset1)

##FULLY BAYESIAN NON-PARAMETRIC REGRESSION MODELS
#A function for fully Bayesian models using brms, with cmdstanr on the backend

get.brms <- function(x){
  if(unique(x$Measure) == "Testosterone"){
    knot.values = c(6,12,18,20,40,70)
  } else {knot.values = c(6,12,18,20,30,40,50,70)
    }
  model.prior <- get_prior(log(Value) ~ 
                             s(Age, k = length(knot.values)) + BMI, 
                           knots = list(Age = knot.values),
                           data = x, family = gaussian())
  model.prior$prior[1]<-"normal(0,10"
                           
  model <- brms::brm(formula = log(Value) ~ 
                       s(Age, k = length(knot.values)) + BMI, 
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
brms.summaries <- map(brms.models, summary)

#Create function to gather tidy predicted draws from each model
getbrmspreds <- function(x){
  #Age values from original datasets
  new_data <- data.frame(Age = sort(unique(PAdf3$Age)),
                         BMI = median(PAdf3$BMI,
                                      na.rm = T))
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

#Compare point-estimates
gam.predvals %>% filter(Age == 29)
brms.predvals %>% filter(Age ==29)

#Comparing run time
system.time(gam.models <- map(PAlist, get.gams))
system.time(brms.models <- map(PAlist, get.brms))

