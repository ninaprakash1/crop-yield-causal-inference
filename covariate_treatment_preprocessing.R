
library(dplyr)       # data wrangling
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees
library(countrycode)
library(data.table)
library(ggplot2)


data<-read.csv(file.path("data", "Data for R.csv"), na.strings = "")

#TODO: decide on handling of outliers
data <- subset(data, Outlier_NT ==0 | Outlier_CT == 0)

data_clean <- data %>% 
  rename(Rotation.Conservation.Crop = R.CT, Irrigation = Irrigation.CT)

### covariates of interest
study_covars <- c('Seq') #id for study 
study_location_covars <- c( "Site.country", 
                            "Latitude",
                            "Longitude")

study_year_covars <- c("Start.year.of.NT...or.first.year.of.experiment.if.missing.", 
                       "Sowing.year",
                       "Harvest.year")

seasonal_covars <- c("Sowing.month", 
                     "Harvesting.month"
)

mgmt_practice_treatments <- c("SC.CT",  #categorical: Yes, No or Mixed
                              "SC.NT",
                              "Rotation.CT",# indicator of any crop rotation
                              "Rotation.NT",
                              "Residue.Management.in.CT.in.last.cropping.season..available.in.current.cropping.season.",# categorical: retained, removed or average
                              "Residue.Management.in.NT.in.last.cropping.season..available.in.current.cropping.season.",
                              "Weed.and.pest.control.CT", #indicator
                              "Weed.and.pest.control.NT",
                              "Fertilization.CT", #categorical: Yes, No or Mix
                              "Fertilization.NT",
                              "Durationof.NT.period.at.sowing..yrs." #duration of treatment (in years)
) 


mgmt_practice_covars <- c("Rotation.Conservation.Crop", #indicator of Crop rotation with >= 3 crops in accordance with FAO def of conservation agriculture 
                          "N.input" #Yes, No or Mix
) 

climatic_soil_vars <- c("pH..surface.layer.", #at experiment site at time of experiment start
                        "P", #Precipitation over the growing season	(mm)
                        "E", #Potential evapotranspiration over the growing season (mm))
                        "PB", #Precipitation balance (mm) := the amount of available water in the growing season for rainfed field
                        "Tave", #Average air temperature during the growing season (C)
                        "Tmax", #Max air temperature during the growing season (C) 
                        "Tmin", #Min air temperature during the growing season (C) 
                        "ST") #Soil texture at experiment  location. Categorical: Sandy Loam; Loam; Silt Loam; Sandy Clay Loam; Clay Loam; Sandy Clay; Clay

crop_covar<- c("Crop")
outcome_vars <- c("Ln_ratio", "Yield.of.CT", "Yield.of.NT")
keep_cols <- c(study_covars,study_location_covars, study_year_covars, seasonal_covars,crop_covar, mgmt_practice_covars, 
               mgmt_practice_treatments, climatic_soil_vars, outcome_vars)
data_clean<- data_clean[,c(keep_cols, 'Yield.of.CT', 'Yield.of.NT')]

#In 2 cases, experimement year is before sowing year, which is likely a data entry error, so we remove these observation pairs:
data_clean <- data_clean[data_clean$Durationof.NT.period.at.sowing..yrs. >= 0,]


#The ph covariate is currently coded as a string and we need to recode this as numeric, handling cases where string contains a range or other edge cases.
#set ph vars for 22 cases where different for treatment and control to Null (also the one case with a comma, "5,2" since unclear what this means-- likely a typo)
msk<- grepl("CT|NT|,", data_clean$pH..surface.layer.)
data_clean[msk,]$pH..surface.layer. <-NA
#remove "with lime application" from end of string in 15 cases
msk<- grepl(" with lime application", data_clean$pH..surface.layer.)
data_clean[msk,]$pH..surface.layer. <- substr(data_clean[msk,]$pH..surface.layer., 1,10)
#remove "or above" from 10 cases where ph var is "\d or above" (treat as \d)
msk<- grepl("above", data_clean$pH..surface.layer.)
data_clean[msk,]$pH..surface.layer. <- gsub("([0-9]+).*", "\\1", data_clean[msk,]$pH..surface.layer.)
# for 164 cases that give a range with "to" or "-": replace "{x1} to/- {x2}" strings with mean of x1 and x2
msk<- grepl("to", data_clean$pH..surface.layer.)
range_strings <- data_clean[msk,]$pH..surface.layer. 
data_clean[msk,]$pH..surface.layer. <- sapply(strsplit(range_strings, split = " to ", fixed = TRUE), function(k) mean(as.numeric(k)))
msk<- grepl("-", data_clean$pH..surface.layer.)
range_strings <- data_clean[msk,]$pH..surface.layer. 
data_clean[msk,]$pH..surface.layer. <- sapply(strsplit(range_strings, split = "-", fixed = TRUE), function(k) mean(as.numeric(k)))


to_factor_cols <- c(c("Seq", "Site.country", "Crop", 'ST'), mgmt_practice_covars)
data_clean[to_factor_cols] <- lapply(data_clean[to_factor_cols], factor)
data_clean$pH..surface.layer. <- as.numeric(data_clean$pH..surface.layer.)


#make factor for for each binary management practice pair
make_treatment_pairs <- function(df, c_treatment_col, t_treatment_col, new_col,  binary_treatments = c("Yes", "No")) {
  new_col_name <- quo_name(new_col)
  df_new <- df %>%
    mutate(!!new_col_name := ifelse(.data[[c_treatment_col]] %in% binary_treatments & .data[[t_treatment_col]] %in% 
                                      binary_treatments & .data[[c_treatment_col]] == .data[[t_treatment_col]],
                                    paste(.data[[c_treatment_col]], .data[[t_treatment_col]]),
                                    NA))
  
  df_new <- df_new[!names(df_new) %in% c(c_treatment_col, t_treatment_col)]
  return (df_new)
}

data_clean <- make_treatment_pairs(df = data_clean, 
                                   c_treatment = "SC.CT",
                                   t_treatment = "SC.NT",
                                   new_col = 'Soil.Cover.Pair')

data_clean <- make_treatment_pairs(df = data_clean, 
                                   c_treatment = "Rotation.CT",
                                   t_treatment = "Rotation.NT",
                                   new_col = 'Crop.Rotation.Pair')

data_clean <- make_treatment_pairs(df = data_clean, 
                                   c_treatment = "Residue.Management.in.CT.in.last.cropping.season..available.in.current.cropping.season.",
                                   t_treatment = "Residue.Management.in.NT.in.last.cropping.season..available.in.current.cropping.season.",
                                   new_col = 'Residue.Management.Pair',
                                   binary_treatments = c("retained", "removed"))

data_clean <- make_treatment_pairs(df = data_clean, 
                                   c_treatment = "Weed.and.pest.control.CT",
                                   t_treatment = "Weed.and.pest.control.NT",
                                   new_col = 'Weed.Pest.Control.Pair')

data_clean <- make_treatment_pairs(df = data_clean, 
                                   c_treatment = "Fertilization.CT",
                                   t_treatment = "Fertilization.NT",
                                   new_col = 'Fertilization.Pair')

dur_max <- max(data_clean$Durationof.NT.period.at.sowing..yrs.)
data_clean<-data_clean %>%
  mutate(treatment_duration_4_grouping = cut(Durationof.NT.period.at.sowing..yrs., 
                                             breaks = c(0,5, 10,20,dur_max), 
                                             include.lowest = TRUE),
         treatment_duration_eq_3 = cut_number(Durationof.NT.period.at.sowing..yrs.,3), 
         treatment_duration_eq_4 = cut_number(Durationof.NT.period.at.sowing..yrs.,4), 
         treatment_duration_eq_5 = cut_number(Durationof.NT.period.at.sowing..yrs.,5))

data_clean <- data_clean[!names(data_clean) %in% c("Durationof.NT.period.at.sowing..yrs.")]



write.csv(data_clean, file.path("data", "data_clean.csv"))
