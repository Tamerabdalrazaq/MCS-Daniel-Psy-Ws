##Function for removing outliers using MAD

Sub_outlier_Detector<- function(data,subject_col_name,
                               Measure_col_name,
                               lower_limit_abs=FALSE,
                               upper_limit_abs=FALSE,
                               lower_limit=3,
                               upper_limit=3,
                               Mad=TRUE){
  
  require(dplyr)
  ## in this kind of function I need: 
  subject_col_name <- enquo(subject_col_name)
  Measure_col_name <- enquo(Measure_col_name)
  
#  stopifnot()
  
  if (Mad==TRUE){
    center_function=median
    spread_function=mad
  } else{
    center_function=mean
    spread_function=sd
  }
  
  
  #  stopifnot
 ### hard-coded thresholds for both upper and lower tresholds
  if ((lower_limit_abs!=FALSE) & (upper_limit_abs!=FALSE)){
    MeasureSubjectwise <- data %>% group_by(!!subject_col_name) %>% 
      summarise(MeanMes = mean(!!Measure_col_name, na.rm = TRUE))
    MeasureSubjectwise_limits <- summarise(group_by(MeasureSubjectwise),
                                            lower = lower_limit_abs,
                                            upper = upper_limit_abs) #calculate MAD
    MeasureSubjectwise <- merge(MeasureSubjectwise,MeasureSubjectwise_limits)
    list_of_outliers_Mes<-MeasureSubjectwise %>% filter(MeanMes<lower | MeanMes>upper) #see outlier subjects
    #RawData <- RawData %>% filter(Subject!= "5f50985cc5f9412102ab71f9" & Subject!= "5f3801b18c88962be7831304" & Subject!= "5982”) #exclude them
    return_list<-list(Subject_Mes_VS_treshold=MeasureSubjectwise,Potential_outliers_list=list_of_outliers_Mes)
    
    if (nrow(return_list$Potential_outliers_list)==0){
      return_list$Potential_outliers_list<-"No Outliers"
    }
    return(return_list)
  }
  
### hard-coded thresholds for only upper threshold
  
  if ((lower_limit_abs==FALSE) & (upper_limit_abs!=FALSE)){
    MeasureSubjectwise <- data %>% group_by(!!subject_col_name) %>% 
      summarise(MeanMes = mean(!!Measure_col_name, na.rm = TRUE))
    MeasureSubjectwise_limits <- summarise(group_by(MeasureSubjectwise),
                                            lower = center_function(MeanMes) - lower_limit*spread_function(MeanMes),
                                            upper = upper_limit_abs) #calculate MAD
    MeasureSubjectwise <- merge(MeasureSubjectwise,MeasureSubjectwise_limits)
    list_of_outliers_Mes<-MeasureSubjectwise %>% filter(MeanMes<lower | MeanMes>upper) #see outlier subjects
    #RawData <- RawData %>% filter(Subject!= "5f50985cc5f9412102ab71f9" & Subject!= "5f3801b18c88962be7831304" & Subject!= "5982”) #exclude them
    return_list<-list(Subject_Mes_VS_treshold=MeasureSubjectwise,Potential_outliers_list=list_of_outliers_Mes)
    
    if (nrow(return_list$Potential_outliers_list)==0){
      return_list$Potential_outliers_list<-"No Outliers"
    }
    return(return_list)
  }
  
### hard-coded thresholds for only lower threshold
  
  if ((lower_limit_abs!=FALSE) & (upper_limit_abs==FALSE)){
    MeasureSubjectwise <- data %>% group_by(!!subject_col_name) %>% 
      summarise(MeanMes = mean(!!Measure_col_name, na.rm = TRUE))
    MeasureSubjectwise_limits <- summarise(group_by(MeasureSubjectwise),
                                            lower = lower_limit_abs,
                                            upper = center_function(MeanMes) + upper_limit*spread_function(MeanMes)) #calculate MAD
    MeasureSubjectwise <- merge(MeasureSubjectwise,MeasureSubjectwise_limits)
    list_of_outliers_Mes<-MeasureSubjectwise %>% filter(MeanMes<lower | MeanMes>upper) #see outlier subjects
    #RawData <- RawData %>% filter(Subject!= "5f50985cc5f9412102ab71f9" & Subject!= "5f3801b18c88962be7831304" & Subject!= "5982”) #exclude them
    return_list<-list(Subject_Mes_VS_treshold=MeasureSubjectwise,Potential_outliers_list=list_of_outliers_Mes)
    
    if (nrow(return_list$Potential_outliers_list)==0){
      return_list$Potential_outliers_list<-"No Outliers"
    }
    return(return_list)
  }
  

  #Code to remove subjects with mean error-rate beyond 3 MADs:
  ##Measure
  MeasureSubjectwise <- data %>% group_by(!!subject_col_name) %>% 
    summarise(MeanMes = mean(!!Measure_col_name, na.rm = TRUE))
  MeasureSubjectwise_limits <- summarise(group_by(MeasureSubjectwise),
                                          lower = center_function(MeanMes) - lower_limit*spread_function(MeanMes),
                                          upper = center_function(MeanMes) + upper_limit*spread_function(MeanMes)) #calculate MAD
  MeasureSubjectwise <- merge(MeasureSubjectwise,MeasureSubjectwise_limits)
  list_of_outliers_Mes<-MeasureSubjectwise %>% filter(MeanMes<lower | MeanMes>upper) #see outlier subjects
  #RawData <- RawData %>% filter(Subject!= "5f50985cc5f9412102ab71f9" & Subject!= "5f3801b18c88962be7831304" & Subject!= "5982”) #exclude them
  return_list<-list(Subject_Mes_VS_treshold=MeasureSubjectwise,Potential_outliers_list=list_of_outliers_Mes)
  
  if (nrow(return_list$Potential_outliers_list)==0){
    return_list$Potential_outliers_list<-"No Outliers"
  }
  return(return_list)
  
}

#Sub_outlier_removal(Data, Subject,RespRT,RespAC)


  