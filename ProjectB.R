
#PROJECT B
install.packages("visreg")
library(visreg)
#PROJECT A

#Identify packages 
pkgs = c("nhanesA",#"tidyverse",
         "magrittr","knitr",
         "patchwork","ggpubr","scales",#"brms",
         "mgcv",
         "sjlabelled","schoenberg","tidymv","broom",
         "broom.mixed","tidybayes","calecopal","smotefamily",
         "corrplot","car","e1071","InformationValue","class","caret")

#Install packages
install.packages(pkgs)

#Load packages
inst = lapply(pkgs, library, character.only = TRUE) 

#Upload data as a csv file
PBdf1 <- read.csv("FullData.csv", stringsAsFactors = FALSE)
str(PBdf1)

#Whitespace is due to utf8 characters - and this is true across all variables
utf8::utf8_print(unique(PBdf1$Ave.Temp), utf8 = FALSE)

#Create a function to replace all utf8 values
replace.values <- function(x) x <- iconv(x, "latin1", "ASCII", sub="")

#NOAA online data source, Seattle Metro area 2000-2021
PBdf2 <- PBdf1 %>% 
  #Remove "T" and "M" values for precipitation, as there is no public documentation to indicate what these values means (and there are only a few instances in the dataset)
  filter(!grepl("T",Precipitation, fixed = T),
         !grepl("M",Precipitation, fixed = T)) %>% 
  #Replace all utf8 values
  mutate_all(replace.values) %>% 
  #Fix dates, which were coded differently for 2021 versus all previous years
  mutate(Date = case_when(grepl("/",Date,fixed = T) ~ 
                            as.Date(Date, format= "%m/%d/%y"),
                          grepl("-",Date,fixed = T) ~ 
                            as.Date(Date, format= "%Y-%m-%d"))) %>% 
  #Separate out Date into month, day, and year
  mutate(Year = lubridate::year(Date), 
         Month = lubridate::month(Date), 
         Day = lubridate::day(Date)) %>% 
  #Make sure relevant variables are in numeric format
  mutate(across(c(Max.Temp:Precipitation,Day,Month,Departure,Year), as.numeric)) %>%
  #Create binary presence/absence variables for freezing and precipitation
  mutate(PrecipB = case_when(Precipitation > 0 ~ 1,
                             Precipitation ==0 ~ 0)) %>% 
  #Select out CDD, which only had a value of zero
  select(PrecipB, ends_with("Temp"),HDD,Departure, Day,Year)

#Splitting the data
#Partition data into training and test datasets using a 75%:25% split ratio
set.seed(1234)
sample_set <- sample(nrow(PBdf2),round(nrow(PBdf2)*0.75), replace = F)
training <- PBdf2[sample_set,]
testing <- PBdf2[-sample_set,]

#Compare class distribution for precipitation Y/N in testing and training datasets
round(prop.table(table(select(PBdf2, PrecipB),exclude = NULL)),4)*100

#Testing and training have similar class distributions
round(prop.table(table(select(training, PrecipB),exclude = NULL)),4)*100
round(prop.table(table(select(testing, PrecipB),exclude = NULL)),4)*100

#Dealing with class imbalance using SMOTE on training dataset
smote <- SMOTE(training, training$PrecipB)
training.adj <- smote$data %>% select(-class)
round(prop.table(table(select(training.adj, PrecipB),exclude = NULL)),4)*100

#Visualizing multicollinearity
PBdf2 %>% keep(is.numeric) %>% cor() %>% corrplot(type = "upper")

#LOGISTICAL REGRESSION
#Training model #1
logitm1 <- glm(data = training.adj, family = binomial, 
                 formula = PrecipB ~ Departure + Day + Year)
summary(logitm1)
vif(logitm1)
exp(coef(logitm1)[c("Departure","Year","Day")])

logitpred1 <- predict(logitm1, testing, type = 'response')

ideal.cutoff <- optimalCutoff(
  actuals = testing$PrecipB,
  predictedScores = logitpred1,
  optimiseFor = "Both"
)

logitpred2 <- ifelse(logitpred1 >= ideal.cutoff,1,0)

logit.ptable <- table(testing$PrecipB, logitpred2)

sum(diag(logit.ptable))/nrow(testing)
#predicting accuracy is 77%

#Naive Bayes
bayes.mod1 <- naiveBayes(PrecipB ~ Departure + Day + Year,
                         data = training.adj, laplace = 1)
bayes.mod1

bayes.pred1 <- predict(bayes.mod1, testing, type = "class")
bayes.pred1.table <- table(testing$PrecipB, bayes.pred1)
bayes.pred1.table
sum(diag(bayes.pred1.table))/nrow(testing)
#78.35% predictive accuracy


#k-Nearest Neighbors
normalize <- function(x){
  return((x-mean(x))/sd(x))
}

#Partition data into training and test datasets using a 75%:25% split ratio
set.seed(1234)
PBdf.knn <- PBdf2 %>% select(PrecipB, Departure, Day, Year) %>% 
  mutate(across(c(Departure, Day, Year), normalize))

knn_sample_set <- sample(nrow(PBdf.knn),round(nrow(PBdf.knn)*0.75), replace = F)
knn_training <- PBdf.knn[knn_sample_set,]
knn_testing <- PBdf.knn[-knn_sample_set,]

knn_smote <- SMOTE(knn_training, knn_training$PrecipB)
knn_training.adj <- knn_smote$data %>% select(-class)

train_labels <- as.factor(pull(knn_training, PrecipB))
test_labels <- as.factor(pull(knn_testing, PrecipB))

kvalues <- list(1,5,10,15,20,25,30,35,40)

get.knearest <- function(x){
  knn.preds <- knn(train = knn_training %>% 
                     select(PrecipB, Departure, Day, Year),
                   test = knn_testing,
                   cl = train_labels,
                   k=x
  )
  knn.predtbl <- table(test_labels, knn.preds)
  pred.acc <- sum(diag(knn.predtbl))/nrow(knn_testing)
  pred.acc
}

knn.df <- plyr::ldply(map(kvalues, get.knearest),data.frame, .id = NULL) %>% 
  rename(Value = X..i..) %>% 
  mutate(Kvalue = as.numeric(kvalues),
         Value = as.numeric(Value))

ggplot(knn.df, aes(x = Kvalue, y = Value))+
  geom_point()+
  geom_line()+
  ylab("Predictive accuracy")+
  xlab("K-value")

#K fold cross validation
set.seed(1234)
sample_set <- createDataPartition(y = PBdf2$PrecipB,
                                  p = 0.75,
                                  list = F)
kfold_train <- PBdf2[sample_set,]
kfold_test <- PBdf2[-sample_set,]

kfold_smote <- SMOTE(kfold_train, kfold_train$PrecipB)
kfold_train.adj <- kfold_smote$data %>% select(-class) %>% 
  mutate(PrecipB = as.factor(PrecipB))

#K-fold cross validation

##GLM
kfoldfit.glm <- train(
  form = PrecipB ~Departure + Day + Year,
  data = kfold_train.adj,
  trControl = trainControl(method = "cv",
                           number = 5),
  method = "glm",
  family = "binomial"
)

kfoldfit.glm$resample %>% 
  arrange(Resample) %>% 
  summarise(AveAccuracy = mean(Accuracy))

##NAive bayes
kfoldfit.nb <- train(
  form = PrecipB ~Departure + Day + Year,
  data = kfold_train.adj,
  trControl = trainControl(method = "cv",
                           number = 5),
  method = "naive_bayes")

kfoldfit.nb$resample %>% 
  arrange(Resample) %>% 
  summarise(AveAccuracy = mean(Accuracy))


##KNN
kfoldfit.knn <- train(
  form = PrecipB ~Departure + Day + Year,
  data = kfold_train.adj,
  trControl = trainControl(method = "cv",
                           number = 5),
  method = "knn")

kfoldfit.knn$resample %>% 
  arrange(Resample) %>% 
  summarise(AveAccuracy = mean(Accuracy))

##Leave One Out Cross-validation
##GLM
LOOCV.glm <- train(
  form = PrecipB ~Departure + Day + Year,
  data = kfold_train.adj,
  trControl = trainControl(method = "LOOCV"),
  method = "glm",
  family = "binomial"
)

LOOCV.glm$results

##Random cross-validation
RCV.glm <- train(
  form = PrecipB ~Departure + Day + Year,
  data = kfold_train.adj,
  trControl = trainControl(method = "LGOCV",
                           p = 0.1,
                           number = 10),
  method = "glm",
  family = "binomial"
)

RCV.glm$resample %>% 
  arrange(Resample)

#Bootstrap Sampling
Bootstrap.glm <- train(
  form = PrecipB ~Departure + Day + Year,
  data = kfold_train.adj,
  trControl = trainControl(method = "boot632",
                           number  = 3),
  method = "glm",
  family = "binomial"
)

Bootstrap.glm$resample %>% 
  arrange(Resample)

##BEYOND PREDICTIVE ACCURACY
Bootstrap.glm$
