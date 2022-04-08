
#PROJECT B
install.packages("smotefamily")
install.packages("car")
library(smotefamily)
library(corrplot)
library(car)
install.packages("e1071")
library(e1071)
install.packages("InformationValue")
library(InformationValue)
library(class)

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
  mutate(Freezing = case_when(Min.Temp > 33 ~ 0,
                              Min.Temp <= 33 ~ 1),
         PrecipB = case_when(Precipitation > 0 ~ 1,
                             Precipitation ==0 ~ 0)) %>% 
  #Select out CDD, which only had a value of zero
  select(PrecipB,Freezing,ends_with("Temp"),HDD,Departure, Day,Year)

#Splitting the data
#Partition data into training and test datasets using a 75%:25% split ratio
set.seed(1234)
sample_set <- sample(nrow(PBdf2),round(nrow(PBdf2)*0.75), replace = F)
training <- PBdf2[sample_set,]
testing <- PBdf2[-sample_set,]

#Compare class distribution for precipitation Y/N in testing and training datasets
round(prop.table(table(select(PBdf2, PrecipB),exclude = NULL)),4)*100

#Testing and training have similiar class distributions
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
train.mod <- glm(data = training.adj, family = binomial, 
                 formula = PrecipB ~ .)
summary(train.mod)
exp(coef(train.mod))#["Max.Temp"]) #Exponentiated coefficient 

#Training model #2
training.adj2 <- training.adj %>% select(-Ave.Temp)
train.mod2 <- glm(data = training.adj2, family = binomial, 
                  formula = PrecipB ~ .)

#Error indicates perfect multicolinearity - which we expected! 
preds.1 <- predict(train.mod2, testing, type = 'response')
#Clearly indicates Max.Temp, Min.Temp,and HDD
summary(train.mod2)
vif(train.mod2)

training.adj3 <- training.adj2 %>% select(-Max.Temp,-Min.Temp)
train.mod3 <- glm(data = training.adj3, family = binomial, 
                  formula = PrecipB ~ Day + Year + HDD)
summary(train.mod3)
vif(train.mod3)

preds.2 <- predict(train.mod3, testing, type = 'response')

ideal.cutoff <- optimalCutoff(
  actuals = testing$PrecipB,
  predictedScores = preds.2,
  optimiseFor = "Both"
)

preds.3 <- ifelse(preds.2 >= ideal.cutoff,1,0)

preds.table <- table(testing$PrecipB, preds.3)
sum(diag(preds.table))/nrow(testing)
#predicting accuracy is 74.46%


#Naive Bayes
bayes.mod1 <- naiveBayes(PrecipB ~ Day + Year + HDD,
                         data = training.adj3, laplace = 1)
bayes.mod1

bayes.pred1 <- predict(bayes.mod1, testing, type = "class")
bayes.pred1.table <- table(testing$PrecipB, bayes.pred1)
bayes.pred1.table
sum(diag(bayes.pred1.table))/nrow(testing)
#78.35% predictive accuracy


#k-Nearest Neighbors
testing.knn <- testing %>% select(-Max.Temp,-Min.Temp,-Ave.Temp)

train_labels <- as.factor(pull(training.adj3, PrecipB))
test_labels <- as.factor(pull(testing.knn, PrecipB))

k.spec <- sqrt(nrow(testing.knn))

knn.preds1 <- knn(train = training.adj3,
                  test = testing.knn,
                  cl = train_labels,
                  k=k.spec
)

head(knn.preds1)

knn.predtbl <- table(test_labels, knn.preds1)
knn.predtbl

sum(diag(knn.predtbl))/nrow(testing.knn)
#75.32% predictive accuracy

#Can change k to test predictive accuracy 
