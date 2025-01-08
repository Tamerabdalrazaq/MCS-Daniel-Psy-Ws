library(dplyr)
library(ggplot2)
library(readr)

MSC_PATH = "C:/Users/averr/OneDrive/Desktop/TAU/Year 4/Psych Workshop/Scripts/1. Object Detection/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"

MSC <- read_csv(MSC_PATH)
# Define script settings
settings <- list(
  filter_bad_detections = TRUE,
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray", "LAST_BUTTON_PRESSED", "im_h", "im_w")
)

# ****************** Helper Functions  ************************ 
# ================================================================



# ****************** Dataset Preprocessing   ************************ 
# ================================================================

# Conver N/A string to NA number for later type conversion
MSC <<- mutate(MSC,
             OBJECT_X1 = na_if(OBJECT_X1, "N\\A"),
             OBJECT_Y1 = na_if(OBJECT_Y1, "N\\A"),
             OBJECT_X2 = na_if(OBJECT_X2, "N\\A"),
             OBJECT_Y2 = na_if(OBJECT_Y2, "N\\A"))
MSC <<- mutate(MSC,
             OBJECT_X1 = as.numeric(OBJECT_X1),
             OBJECT_Y1 = as.numeric(OBJECT_Y1),
             OBJECT_X2 = as.numeric(OBJECT_X2),
            OBJECT_Y2 = as.numeric(OBJECT_Y2),
            CURRENT_FIX_X = as.numeric(CURRENT_FIX_X),
            CURRENT_FIX_Y = as.numeric(CURRENT_FIX_Y))

# Convert featrues to be type int

MSC <<- MSC %>%
  mutate(
    across(c(
      OBJECT_X1, OBJECT_Y1, OBJECT_X2, OBJECT_Y2, 
      CURRENT_FIX_X, CURRENT_FIX_Y, subjectnum, 
      im_h, im_w, CURRENT_FIX_INDEX
    ), as.numeric)
  )


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
MSC <- MSC %>%
  group_by(subjectnum) %>%
  filter(
    is.na(lag(TRIAL_INDEX)) |  # Keep the first row in each group
      TRIAL_INDEX == lag(TRIAL_INDEX) + 1  # Keep if consecutive
  ) %>%
  ungroup()  # Remove grouping if needed
# ================================================================


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
      hit = list(hit)
    )
 }
# ================================================================
