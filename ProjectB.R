
#PROJECT B

#Upload data as a csv file
df1 <- read.csv("FullData.csv", stringsAsFactors = FALSE)

#Whitespace is due to utf8 characters - and this is true across all variables
utf8::utf8_print(unique(df1$Ave.Temp), utf8 = FALSE)

#Create a function to replace all utf8 values
replace.values <- function(x) x <- iconv(x, "latin1", "ASCII", sub="")

#NOAA online data source, Seattle Metro area 2000-2021
df2 <- df1 %>% 
  #Remove "T" values for precipitation, no code to indicate what this value means
  filter(!grepl("T",Precipitation, fixed = T)) %>% 
  #Replace all utf8 values
  mutate_all(replace.values) %>% 
  #Fix dates, which were coded differently for 2021 versus all previous years
  mutate(Date = case_when(grepl("/",Date,fixed = T) ~ as.Date(Date, format= "%m/%d/%y"),
                          grepl("-",Date,fixed = T) ~ as.Date(Date, format= "%Y-%m-%d"))) %>% 
  #Separate out Date into month,day,year
  mutate_at(vars(Date), tibble::lst(year, month, day)) %>% 
  #Make sure variables are in numeric format
  mutate(across(c(Max.Temp:Precipitation,day,year), as.numeric)) %>% 
  #Define "Freezing" and "PrecipB" binary variables for analysis
  mutate(Freezing = case_when(Min.Temp > 33 ~ 0,
                              Min.Temp <= 33 ~ 1),
         PrecipB = case_when(Precipitation > 0 ~ 1,
                             Precipitation ==0 ~ 0)) 

hist(df2$Max.Temp)

#Test model
log.model <- glm(PrecipB ~ day + year + Ave.Temp, 
                 data = Bdf1, family = 'binomial')
summary(log.model)
plot(log.model)
summary(log.model)$coefficients

predict(log.model, 
        newdata = list(Precipitation = c(0, 0.08,0.36)))
