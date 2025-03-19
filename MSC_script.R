library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
rm(list = ls())

# TODO:
# Calculate the baseline probability of hitting a random rectangle with the average area of the detected boxes
# == average(area of detected boxes)/average(area of images) 
# area of detected box = (x2-x1) * (y2-y1)
# area of image = im_w * im_h

DISPLAY_W <- 1280
DISPLAY_H <- 800

MSC_PATH = "C:/Users/danie/OneDrive/Documents/אוניברסיטה/דוקטורט/PhD/Projects/MCS_dataset - Tamer/MCS_dataset/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"
MSC <- read_csv(MSC_PATH)
MSC_ORIGINAL <- MSC
# Define script settings
settings <- list(
  filter_bad_detections = TRUE,
  filter_inconsecutive_trial = TRUE,
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray"),
  fliter_prev_absent = TRUE,
  test_subject = TRUE
)

# ******************   ************************ 
# ===============================Helper Functions=================================


## Data Preprocessing ----
# ******************    ************************ 
# ================================Dataset Preprocessing================================

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



# MSC_ORIGINAL%>% filter( CURRENT_FIX_INDEX == 2) %>% ggplot()+geom_point(aes(x=CURRENT_FIX_X, y=CURRENT_FIX_Y))+ylim(0,800) +xlim(0,1280)


# ******************   ************************ 
# ===========================Adding New Features=====================================

# Add hit - is fixation inside object
MSC <<- mutate(MSC, hit = ifelse(OBJECT_X1 <= CURRENT_FIX_X &
                                   OBJECT_Y1 <= CURRENT_FIX_Y &
                                   OBJECT_X2 >= CURRENT_FIX_X &
                                   OBJECT_Y2 >= CURRENT_FIX_Y
                                 , 1, 0))

MSC <<- mutate(MSC, p_hit = ((OBJECT_X2-OBJECT_X1) * (OBJECT_Y2-OBJECT_Y1))/(im_w * im_h)
               )

# 
# 
# MSC <<- mutate(MSC, prev_hit = ifelse(lag(OBJECT_X1) <= CURRENT_FIX_X &
#                                         lag(OBJECT_Y1) <= CURRENT_FIX_Y &
#                                         lag(OBJECT_X2) >= CURRENT_FIX_X &
#                                         lag(OBJECT_Y2) >= CURRENT_FIX_Y
#                                       , 1, 0))



# ============================Clustering==================================== 


# Cluster similar rows by subjectnum and trial_index. Merge unique data into a single array.
if (settings$cluster_fixations){
  MSC <<- MSC %>%
    group_by(subjectnum, TRIAL_INDEX) %>%
    summarize(
      catcue = first(catcue),
      condition = first(condition),
      img_res = first(img_res),
      im_w = first(im_w),
      im_h = first(im_h),
      img_name = first(img_name),
      TRIAL_FIXATION_TOTAL = first(TRIAL_FIXATION_TOTAL),
      FIRST_FIX_X = first(CURRENT_FIX_X),
      FIRST_FIX_Y = first(CURRENT_FIX_Y),
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
      p_hit = first(p_hit),
      # prev_hit = first(prev_hit),
      LAST_BUTTON_PRESSED = first(LAST_BUTTON_PRESSED),
    )
}



# ============================More Features==================================== 

MSC <- MSC %>%
  mutate(
    PREV_OBJ_X1 = lag(OBJECT_X1, default = NA),
    PREV_OBJ_Y1 = lag(OBJECT_Y1, default = NA),
    PREV_OBJ_X2 = lag(OBJECT_X2, default = NA),
    PREV_OBJ_Y2 = lag(OBJECT_Y2, default = NA),
    )
MSC <- MSC %>%
  mutate(
    #p_prev_hit1 = lag(p_hit, default = NA),
    prev_p_hit = ((PREV_OBJ_X2-PREV_OBJ_X1) * (PREV_OBJ_Y2-PREV_OBJ_Y1))/(im_w * im_h),
    prev_condition = lag(condition, default = NA)
    )


# Custom function
calc_hit <- function(FIX_X, FIX_Y, PREV_X1, PREV_Y1, PREV_X2, PREV_Y2) {
  return(ifelse(PREV_X1 <= FIX_X & PREV_Y1 <= FIX_Y & PREV_X2 >= FIX_X & PREV_Y2 >= FIX_Y, 1, 0))
}

# prev_hit
MSC$prev_hit <- mapply(function(x1, x2, x3, x4, x5, x6) calc_hit(x1, x2, x3, x4, x5, x6),
                    MSC$CURRENT_FIX_X, MSC$CURRENT_FIX_Y, MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
                    MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2, SIMPLIFY = FALSE)



# Another way to do the previous
#pivoted_MC = MSC_ORIGINAL %>%pivot_wider(names_from = CURRENT_FIX_INDEX,
#                                    id_cols = c(subjectnum:TRIAL_FIXATION_TOTAL,LAST_BUTTON_PRESSED),
#                                    values_from = c(CURRENT_FIX_X:CURRENT_FIX_DURATION)
#)

# =========================Filter Data=======================================
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
    filter(
        TRIAL_INDEX == lag(TRIAL_INDEX) + 1
    )
}



if (settings$fliter_prev_absent){
  MSC <- MSC %>% filter(prev_condition == "present")
}

if (settings$test_subject){
  MSC <- MSC %>% filter(subjectnum == 1)
}


MSC = MSC%>%mutate(
  prev_hit_1 = sapply(prev_hit, function(x) x[1]),
  prev_hit_2 = sapply(prev_hit, function(x) x[2]),
  prev_hit_3 = sapply(prev_hit, function(x) x[3]),
  prev_hit_4 = sapply(prev_hit, function(x) x[4]),
  prev_hit_5 = sapply(prev_hit, function(x) x[5]),
  prev_hit_6 = sapply(prev_hit, function(x) x[6]),
  prev_hit_7 = sapply(prev_hit, function(x) x[7]),
  prev_hit_8 = sapply(prev_hit, function(x) x[8])
)

##
table(MSC$condition,MSC$prev_condition)


prop.table(table(MSC$prev_hit_1))
prop.table(table(MSC$prev_hit_2))
prop.table(table(MSC$prev_hit_3))
prop.table(table(MSC$prev_hit_4))
prop.table(table(MSC$prev_hit_5))
prop.table(table(MSC$prev_hit_6))
prop.table(table(MSC$prev_hit_7))
prop.table(table(MSC$prev_hit_8))
#prop.table(table(MSC$prev_hit_2))

## Some indiciation of an effect
MSC%>%summarize(mean(prev_p_hit, na.rm = TRUE))


# ============================Orthogonal Projection==================================== 

get_projection <- function(x, y, x1, y1, x2, y2) {
  if(is.na(x) | is.na(y) |is.na(x1) |is.na(y1) |is.na(x2) |is.na(y2)){
    return (NA)
  }
  if (x1 <= x & x <= x2 & y1 <= y & y <= y2 ){
    return (c(x,y))
  }
  # Case 1: x1 <= x <= x2
  if (x1 <= x & x <= x2) {
    if (y <= y1) {
      return(c(x, y1)) # Project to (x, y1)
    } else if (y >= y2) {
      return(c(x, y2)) # Project to (x, y2)
    }
  }

  # Case 2: y1 <= y <= y2
  if (y1 <= y & y <= y2) {
    if (x <= x1) {
      return(c(x1, y)) # Project to (x1, y)
    } else if (x >= x2) {
      return(c(x2, y)) # Project to (x2, y)
    }
  }

  # Case 3: x' and y' are the closest points
  # x' = argmin_i(|x - x_i|), y' = argmin_i(|y - y_i|)
  x_prime <- ifelse(abs(x - x1) < abs(x - x2), x1, x2)
  y_prime <- ifelse(abs(y - y1) < abs(y - y2), y1, y2)
  return(c(x_prime, y_prime))
}

euc_l2 <- function(x1, y1, x2, y2) {
  # Calculate the Euclidean distance
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}

calc_orthogonal_distance <- function(x, y, obj_x1, obj_y1, obj_x2, obj_y2) {
  projection <- get_projection(x, y, obj_x1, obj_y1, obj_x2, obj_y2)
  distance <- euc_l2(x, y, projection[1], projection[2])
  return(distance)
}

# Accepts lists X,Y and scalars x1,y1,x2,y2 and returns a new list by applying calc_orthogonal_distance(x, y, x1, y1, x2, y2)
calc_orthogonal_distances <- function(X, Y, x1, y1, x2, y2) {
  # Check if lengths of X and Y are the same
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length")
  }
  
  # Initialize an empty vector to store results
  distances <- numeric(length(X))
  
  # Iterate over the elements of X and Y
  for (i in seq_along(X)) {
    distances[i] <- calc_orthogonal_distance(X[i], Y[i], x1, y1, x2, y2)
  }
  
  # Return the list of distances
  return(distances)
}

# Apply the function to the MSC dataset
MSC <- MSC %>%
  mutate(
    PREV_ORTHOGONAL_DISTANCE = ifelse(
        !is.na(PREV_OBJ_X1) & !is.na(PREV_OBJ_Y1) & !is.na(PREV_OBJ_X2) & !is.na(PREV_OBJ_Y2),
      mapply(
        calc_orthogonal_distances,
        CURRENT_FIX_X, CURRENT_FIX_Y,
        PREV_OBJ_X1, PREV_OBJ_Y1, PREV_OBJ_X2, PREV_OBJ_Y2
      ),
      NA  # Assign NA if any required column has NA
    )
  )

##needs fixing
MSC = MSC%>%mutate(
  PREV_ORTHOGONAL_DISTANCE_1 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) unlist(x[1])),
  PREV_ORTHOGONAL_DISTANCE_2 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[2]),
  PREV_ORTHOGONAL_DISTANCE_3 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[3]),
  PREV_ORTHOGONAL_DISTANCE_4 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[4]),
  PREV_ORTHOGONAL_DISTANCE_5 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[5]),
  PREV_ORTHOGONAL_DISTANCE_6 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[6]),
  PREV_ORTHOGONAL_DISTANCE_7 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[7]),
  PREV_ORTHOGONAL_DISTANCE_8 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[8])
)

MSC%>%summarize(mean(PREV_ORTHOGONAL_DISTANCE_1, na.rm = TRUE))

## what is the baseline?
## randomly generate dots (based on the previous target window (fixed) and the size of the current image (fixed))
## calculate the distance between the dots and the (previous target) rectangle using your function 
## Average the distance for each picture

#monte carlo simulations/ agent based simulations

