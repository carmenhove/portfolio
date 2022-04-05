#NOAA online data source, Seattle Metro area 2000-2021

NOAAdf1 <- read.csv("FullData.csv", stringsAsFactors = FALSE) 
colnames(NOAAdf1)
str(NOAAdf1)
#table(NOAAdf1$X.Daily.Data) 

#Where do the mislabeled data come from?
NOAAdf1 %>% filter(X.Daily.Data == "Temperature")
NOAAdf1 %>% filter(Year %in% c("2011","2016"))
#2011 and 2016 = Temperature needs to be relabeled as "Max Temperature"

NOAAdf2 <- NOAAdf1 %>% 
  mutate(President = case_when(Year %in% c(2000:2007) ~ "Bush Jr",
                               Year %in% c(2008:2015) ~ "Obama",
                               Year %in% c(2016:2019) ~ "Trump",
                               Year %in% c(2020:2021) ~ "Biden"),
         Parameter = case_when(X.Daily.Data == "Temperature" ~ 
                                 "Max Temperature",
                               TRUE ~ X.Daily.Data)) %>% 
  mutate(Parameter = gsub(" ", "_", Parameter)) %>% 
  select(Parameter,Observed,Year,President) %>% 
  pivot_wider(names_from = Parameter, values_from = Observed) %>% 
  mutate(across(c(3:6), as.numeric),
         PrecipB = case_when(Precipitation == 0 ~ "No",
                             TRUE ~ "Yes"))
         
read.delim(file, header = TRUE, sep = "\t", dec = ".", ...)
