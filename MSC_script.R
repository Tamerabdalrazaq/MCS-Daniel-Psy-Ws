library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
rm(list = ls())

# Todo: Fix coordinates display problem
# Check fixation landing on previous trial's object

DISPLAY_W <- 1280
DISPLAY_H <- 800

MSC_PATH = "C:/Users/averr/OneDrive/Desktop/TAU/Year 4/Psych Workshop/Scripts/1. Object Detection/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"
MSC <- read_csv(MSC_PATH)
MSC_ORIGINAL <- MSC
# Define script settings
settings <- list(
  filter_bad_detections = TRUE,
  filter_inconsecutive_trial = TRUE,
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray", "im_h", "im_w"),
  fliter_absent = TRUE
)

# ****************** Helper Functions  ************************ 
# ================================================================


## Data Preprocessing ----
# ****************** Dataset Preprocessing   ************************ 
# ================================================================

# Conver N/A string to NA number for later type conversion
MSC <<- mutate(MSC,
             OBJECT_X1 = na_if(OBJECT_X1, "N\\A"),
             OBJECT_Y1 = na_if(OBJECT_Y1, "N\\A"),
             OBJECT_X2 = na_if(OBJECT_X2, "N\\A"),
             OBJECT_Y2 = na_if(OBJECT_Y2, "N\\A"))

# Convert featrues to be type int

MSC <<- MSC %>%
  mutate(
    across(c(
      OBJECT_X1, OBJECT_Y1, OBJECT_X2, OBJECT_Y2, 
      CURRENT_FIX_X, CURRENT_FIX_Y, subjectnum, 
      im_h, im_w, CURRENT_FIX_INDEX
    ), as.numeric)
  )


# Convert fixation coordinats from relative to absolute
MSC <<- MSC %>% mutate(OBJECT_X1 = OBJECT_X1 + (DISPLAY_W - im_w)/2,
                      OBJECT_X2 = OBJECT_X2 + (DISPLAY_W - im_w)/2,
                      OBJECT_Y1 = OBJECT_Y1 + (DISPLAY_H - im_h)/2,
                      OBJECT_Y2 = OBJECT_Y2 + (DISPLAY_H - im_h)/2)

# MSC <<- MSC %>% mutate(CURRENT_FIX_X = CURRENT_FIX_X - (DISPLAY_W - im_w)/2, 
#                        CURRENT_FIX_Y = CURRENT_FIX_Y - (DISPLAY_H - im_h)/2)


# merge  (im_h, im_w) -> img_res 
MSC <<- MSC %>%
  mutate(
    img_res = paste(im_w, im_h, sep = ",")
  ) %>%
  relocate(
    img_res, .after = im_w
  )
  

# Remove unnecessary features
MSC <- MSC %>%
   select(-all_of(settings$features_to_remove))


# Sort by subjectnum, TRIAL_INDEX
MSC <- MSC %>%
  arrange(subjectnum, TRIAL_INDEX)


# ================================================================

#MSC%>% filter(subjectnum == 1 & TRIAL_INDEX < 100) %>% ggplot()+geom_point(aes(x=CURRENT_FIX_X, y=CURRENT_FIX_Y))+ylim(0,800) +xlim(0,1280)


# ****************** Adding New Features  ************************ 

# Add hit - is fixation inside object
MSC <<- mutate(MSC, hit = ifelse(OBJECT_X1 <= CURRENT_FIX_X & 
                                       OBJECT_Y1 <= CURRENT_FIX_Y &
                                       OBJECT_X2 >= CURRENT_FIX_X &
                                       OBJECT_Y2 >= CURRENT_FIX_Y
                                     , 1, 0))

# ================================================================ 




# ****************** Clustering  ************************ 


# Cluster similar rows by subjectnum and trial_index. Merge unique data into a single array.
 if (settings$cluster_fixations){
   MSC <<- MSC %>%
    group_by(subjectnum, TRIAL_INDEX) %>%
    summarize(
      catcue = first(catcue),
      condition = first(condition),
      img_res = first(img_res),
      img_name = first(img_name),
      TRIAL_FIXATION_TOTAL = first(TRIAL_FIXATION_TOTAL),
      CURRENT_FIX_X = list(CURRENT_FIX_X),
      CURRENT_FIX_Y = list(CURRENT_FIX_Y),
      CURRENT_FIX_DURATION = list(CURRENT_FIX_DURATION),
      img_name = first(img_name),
      OBJECT_X1 = first(OBJECT_X1),
      OBJECT_Y1 = first(OBJECT_Y1),
      OBJECT_X2 = first(OBJECT_X2),
      OBJECT_Y2 = first(OBJECT_Y2),
      confidence = first(confidence),
      hit = list(hit),
      LAST_BUTTON_PRESSED = first(LAST_BUTTON_PRESSED),
    )
 }


# Another way to do the previous
#pivoted_MC = MSC_ORIGINAL %>%pivot_wider(names_from = CURRENT_FIX_INDEX,
#                                    id_cols = c(subjectnum:TRIAL_FIXATION_TOTAL,LAST_BUTTON_PRESSED),
#                                    values_from = c(CURRENT_FIX_X:CURRENT_FIX_DURATION)
#)
# ================================================================


# ****************** Filter Data  ************************ 

# Filter incorrect image detections
if (settings$filter_bad_detections) {
  MSC <<- MSC %>%
    filter(
      (condition == "present" & !is.na(OBJECT_X1)) |
        (condition == "absent" & is.na(OBJECT_X1))
    )
}


# Filter out inconsecutive trials
if (settings$filter_inconsecutive_trial){
  MSC <- MSC %>%
    group_by(subjectnum) %>%
    filter(
      is.na(lag(TRIAL_INDEX)) |  # Keep the first row in each group
        TRIAL_INDEX == lag(TRIAL_INDEX) + 1  # Keep if consecutive
    ) %>%
    ungroup()  # Remove grouping if needed
}

if (settings$fliter_absent){
  MSC <- MSC %>% filter(condition == "present")
}
# ================================================================

