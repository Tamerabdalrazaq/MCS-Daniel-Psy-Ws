library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
library(purrr)

rm(list = ls())

DISPLAY_W <- 1280
DISPLAY_H <- 800
CENTER <- c(DISPLAY_W, DISPLAY_H)/2

MSC_PATH = "C:/Users/averr/OneDrive/Desktop/TAU/Year 4/Semestar A/Psych Workshop/Scripts/1. Object Detection/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"
MSC <- read_csv(MSC_PATH)
MSC_ORIGINAL <- MSC


# Define script settings
settings <- list(
  filter_bad_detections = TRUE,
  prev_condition = "present", 
  curr_condition = "absent", # # TODO: Run when absent and present and both
  filter_inconsecutive_trial = TRUE,
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray"),
  fliter_prev_absent = TRUE,
  test_subject = 1,
  calc_mean_OG = FALSE,
  N_BASELINE_FIXATIONS = 10 #relevant fixations for the baseline
)

# ================================Dataset Preprocessing================================

# Conver N/A string to NA number for later type conversion
MSC <- mutate(MSC,
               OBJECT_X1 = na_if(OBJECT_X1, "N\\A"),
               OBJECT_Y1 = na_if(OBJECT_Y1, "N\\A"),
               OBJECT_X2 = na_if(OBJECT_X2, "N\\A"),
               OBJECT_Y2 = na_if(OBJECT_Y2, "N\\A"))

# Convert featrues to be type int

MSC <- MSC %>%
  mutate(
    across(c(
      OBJECT_X1, OBJECT_Y1, OBJECT_X2, OBJECT_Y2, 
      CURRENT_FIX_X, CURRENT_FIX_Y, subjectnum, 
      im_h, im_w, CURRENT_FIX_INDEX
    ), as.numeric)
  )


# Convert fixation coordinates from relative to absolute
MSC <- MSC %>% mutate(OBJECT_X1 = OBJECT_X1 + (DISPLAY_W - im_w)/2,
                       OBJECT_X2 = OBJECT_X2 + (DISPLAY_W - im_w)/2,
                       OBJECT_Y1 = OBJECT_Y1 + (DISPLAY_H - im_h)/2,
                       OBJECT_Y2 = OBJECT_Y2 + (DISPLAY_H - im_h)/2)



# Remove unnecessary features
MSC <- MSC %>%
  select(-all_of(settings$features_to_remove))


# Sort by subjectnum, TRIAL_INDEX
MSC <- MSC %>%
  arrange(subjectnum, TRIAL_INDEX)

## ============================Data Grouping==================================== 


# Group similar rows by subjectnum and trial_index. Merge unique data into a single array.
if (settings$cluster_fixations){
  MSC <- MSC %>%
    group_by(subjectnum, TRIAL_INDEX) %>%
    summarize(
      condition = first(condition),
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
    )
}



## ============================Add New Features==================================== 

# Relevant New Features
MSC <- MSC %>%
  mutate(
    PREV_OBJ_X1 = lag(OBJECT_X1, default = NA), # Previous Object Coordinates
    PREV_OBJ_Y1 = lag(OBJECT_Y1, default = NA),
    PREV_OBJ_X2 = lag(OBJECT_X2, default = NA),
    PREV_OBJ_Y2 = lag(OBJECT_Y2, default = NA),
    prev_condition = lag(condition, default = NA), # Previous Condition
    lag_trial_index = lag(TRIAL_INDEX) # Previous Trial Index
  )



# Images Baseline
# Contains distributions of fixations accross subject for each image
# Features: Fixations_X, Fixations_Y - Nxk matrix:
#                      N number of subjects that saw the same image, k is the max number of fixations accross the N subjects)
#           subjectnum: 
#                      Nx1 vector containing the subject number of the relevant subjects
IMAGES_DF <- MSC %>%
  group_by(img_name) %>%
  summarise(
    FIXATIONS_X = {
      # Get all vectors in current group
      vectors <- CURRENT_FIX_X
      max_len <- max(lengths(vectors))  # k
      
      # Create matrix Nxk with proper dimensions
      mat <- matrix(nrow = length(vectors), ncol = max_len)
      
      # Fill matrix row by row
      for (i in seq_along(vectors)) {
        vec <- vectors[[i]]
        mat[i, 1:length(vec)] <- vec
      }
      
      # Return as list to preserve matrix structure
      list(mat)
    },
    FIXATIONS_Y = {
      # Get all vectors in current group
      vectors <- CURRENT_FIX_Y
      max_len <- max(lengths(vectors))
      
      # Create matrix with proper dimensions
      mat <- matrix(nrow = length(vectors), ncol = max_len)
      
      # Fill matrix row by row
      for (i in seq_along(vectors)) {
        vec <- vectors[[i]]
        mat[i, 1:length(vec)] <- vec
      }
      
      # Return as list to preserve matrix structure
      list(mat)
    },
    subjectnum = list(subjectnum),
    .groups = "drop"
  )

## =========================Filter Data=======================================
# Filter incorrect image detections
if (settings$filter_bad_detections) {
  MSC <- MSC %>%
    filter(
      (condition == "present" & !is.na(OBJECT_X1)) |
        (condition == "absent" & is.na(OBJECT_X1))
    )
}

# Filter out nonconsecutive trials
if (settings$filter_inconsecutive_trial){
  MSC <- MSC %>%
    filter(
      TRIAL_INDEX == lag_trial_index + 1
    )
}


### Dataframe for calculating subject fixations baseline 

if (settings$fliter_prev_absent){
  MSC <- MSC %>% filter(prev_condition == settings$prev_condition &
                          condition == settings$curr_condition)
}

FIXATIONS_BASELINE_DF <- MSC

if (settings$test_subject){
  MSC <- MSC %>% filter(subjectnum == settings$test_subject)
}



## =========================Subject Fixations Vectors For Baselines=======================================

# Extract fixation coordinates dynamically
for (i in seq_len(settings$N_BASELINE_FIXATIONS)) {
  FIXATIONS_BASELINE_DF[[paste0("FIX_", i, "_X")]] <- sapply(FIXATIONS_BASELINE_DF$CURRENT_FIX_X, function(x) x[i])
  FIXATIONS_BASELINE_DF[[paste0("FIX_", i, "_Y")]] <- sapply(FIXATIONS_BASELINE_DF$CURRENT_FIX_Y, function(x) x[i])
}


# Transfoms the fixations from absolute to relative to an origin (CENTER)
get_diff_vector <- function(coords, origin) {
  # coords: a vector of length 2, c(c1, c2)
  # origin: a vector of length 2, c(DISPLAY_W, DISPLAY_H)
  
  if (length(coords) != 2 || length(origin) != 2) {
    print(coords)
    print(origin)
    stop("Both 'coords' and 'origin' must be numeric vectors of length 2.")
  }
  
  return(coords - origin)
}

# Merges FIX_i_X and FIX_i_Y to FIX_i_VECTOR after applying the above transformation
for (i in 1:settings$N_BASELINE_FIXATIONS) {
  x_col <- paste0("FIX_", i, "_X")
  y_col <- paste0("FIX_", i, "_Y")
  vec_col <- paste0("FIX_", i, "_VECTOR")
  
  FIXATIONS_BASELINE_DF[[vec_col]] <- mapply(
    function(c1, c2) get_diff_vector(c(c1, c2), CENTER),
    FIXATIONS_BASELINE_DF[[x_col]], FIXATIONS_BASELINE_DF[[y_col]], SIMPLIFY = FALSE)
}

PREV_CURR_ABSENT_FIX_VECTORS <- FIXATIONS_BASELINE_DF %>%
  group_by(subjectnum) %>%
  summarize(
    FIX_1_VECTOR = list(FIX_1_VECTOR),
    FIX_2_VECTOR = list(FIX_2_VECTOR),
    FIX_3_VECTOR = list(FIX_3_VECTOR),
    FIX_4_VECTOR = list(FIX_4_VECTOR),
    FIX_5_VECTOR = list(FIX_5_VECTOR),
    FIX_6_VECTOR = list(FIX_6_VECTOR),
    FIX_7_VECTOR = list(FIX_7_VECTOR),
    FIX_8_VECTOR = list(FIX_8_VECTOR),
  )


# ============================Calculation of measures==================================== 
## ============================Distance==================================== 

# Calculate the orgthogonal projection of (x,y) onto the rectangle (x1, y1, x2, y2)
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


## ============================Angles==================================== 

calc_norm <- function(v) {
  dist <- euc_l2(v[1], v[2], 0, 0)
  return(dist)
}

calc_dot_product <- function(v1, v2, normalize = TRUE, max_norm = NULL) {
  if (length(v1) != length(v2)) {
    stop("Vectors must be the same length.")
  }
  
  if (any(is.na(v1)) || any(is.na(v2))) {
    return(NA)
  }
  
  # Handle vector length normalization based on normalize (i.e len = 1) or max_norm  
  v1_norm = calc_norm(v1)
  v2_norm = calc_norm(v2)
  if(normalize){
    v1 <- v1*(1/v1_norm)
    v2 <- v2*(1/v2_norm)
  } else if (!is.null(max_norm)){
    if(v1_norm > max_norm){
      v1 <- (v1/v1_norm)*max_norm
    }
    if(v2_norm > max_norm){
      v2 <- (v2/v2_norm)*max_norm
    }
    v1 <- v1 * (1/max_norm)
    v2 <- v2 * (1/max_norm)
  }
  return(sum(v1 * v2))
}

measure_dot_product <- function(v1, v2) {
    return (calc_dot_product(v1, v2,
                             normalize = TRUE, max_norm = NULL))
}

get_curr_fix_prev_obj_dot_product <- function(FIX_X, FIX_Y, obj_x1, obj_y1, obj_x2, obj_y2){
  products <- numeric(length(FIX_X))
  obj_center <- c(obj_x1+obj_x2, obj_y1+obj_y2)/2
  obj_vector <- get_diff_vector(obj_center, CENTER)
  
  # Iterate over the elements of X and Y
  for (i in seq_along(FIX_X)) {
    fix_vector <- get_diff_vector(c(FIX_X[i], FIX_Y[i]), CENTER)
    products[i] <- measure_dot_product(fix_vector, obj_vector)
  }
  
  # Return the list of distances
  return(products)
}

MSC <- MSC %>%
  mutate(
    PREV_OBJ_FIX_DOT_PRODUCT = ifelse(
      !is.na(PREV_OBJ_X1) & !is.na(PREV_OBJ_Y1) & !is.na(PREV_OBJ_X2) & !is.na(PREV_OBJ_Y2),
      mapply(
        get_curr_fix_prev_obj_dot_product,
        CURRENT_FIX_X, CURRENT_FIX_Y,
        PREV_OBJ_X1, PREV_OBJ_Y1, PREV_OBJ_X2, PREV_OBJ_Y2
      ),
      NA  # Assign NA if any required column has NA
    )
    
  )

# ============================Calculating Different Baselines==================================== 



# Calculate the average OG distance of a random point
# from the * previous * detection



  # Calculate the OG distance between center and  the * previous * detection
  MSC$CENTER_OG_DISTANCE <- mapply(calc_orthogonal_distance,
                                 DISPLAY_W/2, DISPLAY_H/2,
                                 MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
                                 MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2)
  
  
  
  # Calculate the OG distance between a specific fixation mean and  the * previous * detection
  
  # For previous absent
  n_fixations <- 8
  #Fix 1
  
  calc_dist_i_og_distance <- function(fix_num){
    Map(
      function(subj, x1, y1, x2, y2) {
        idx <- which(PREV_CURR_ABSENT_FIX_VECTORS$subjectnum == subj)
        if (length(idx) == 0) {
          warning("No match for subjectnum")
          return(NA)
        }
        if (length(idx) > 1) {
          warning("Multiple matches for subjectnum")
          return(NA)
        }
        vector_name <- paste0("FIX_", fix_num, "_VECTOR")
        FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS[[vector_name]][[idx]]
        distances <- numeric(length(FIX_VECTORS))
        for (i in seq_along(FIX_VECTORS)) {
          # IMPORTANT! curr_vector is relative to the center and not absolute to the screen 
          curr_vector = FIX_VECTORS[[i]][1] + CENTER 
          x =  curr_vector[1]
          y = curr_vector[2]
          dist =calc_orthogonal_distance(x, y, x1, y1, x2, y2)
          distances[i] = dist
        }
        return(distances)
      },
      MSC$subjectnum,
      MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
      MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
    )
    }
  
  for (i in 1:n_fixations) {
    col_name <- paste0("DIST_FIX_", i, "_PREV_OBJ_OG_DISTANCE")
    avg_col_name <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_OG_DISTANCE")
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calc_dist_i_og_distance(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  
  ### Old Baseline Calculation
  ### It is problematic because we average the fixations before the non-linear calculations

#   MSC$AVG_FIX_1_PREV_OG_DISTANCE <- mapply(
#     function(subj, x1, y1, x2, y2) {
#       idx <- which(PREV_ABSENT_FIX_AVG$subjectnum == subj)
#       if (length(idx) == 0) return(NA)
#       
#       fix_x <- PREV_ABSENT_FIX_AVG$FIX_1_X[idx]
#       fix_y <- PREV_ABSENT_FIX_AVG$FIX_1_Y[idx]
#       
#       calc_orthogonal_distance(fix_x, fix_y, x1, y1, x2, y2)
#     },
#     MSC$subjectnum,
#     MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
#     MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
#   )
#   
# 

#####=== Subject-Based Baseline for dot product vector difference===#####
MSC$DIST_FIX_1_PREV_OBJ_DOT_PRODUCT <- Map(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_CURR_ABSENT_FIX_VECTORS$subjectnum == subj)
    if (length(idx) == 0) {
      warning("No match for subjectnum")
      return(NA)
    }
    if (length(idx) > 1) {
      warning("Multiple matches for subjectnum")
      return(NA)
    }
    obj_center = get_diff_vector(c(x1+x2, y1+y2)/2, CENTER)
    FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_1_VECTOR[[idx]]
    products <- numeric(length(FIX_VECTORS))
    for (i in seq_along(FIX_VECTORS)) {
      product = measure_dot_product(FIX_VECTORS[[i]], obj_center)
      products[i] = product
    }
    return(products)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)

MSC <- MSC %>%
  mutate(AVG_DIST_FIX_1_PREV_OBJ_DOT_PRODUCT = sapply(DIST_FIX_1_PREV_OBJ_DOT_PRODUCT, mean))


MSC$DIST_FIX_2_PREV_OBJ_DOT_PRODUCT <- Map(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_CURR_ABSENT_FIX_VECTORS$subjectnum == subj)
    if (length(idx) == 0) {
      warning("No match for subjectnum")
      return(NA)
    }
    if (length(idx) > 1) {
      warning("Multiple matches for subjectnum")
      return(NA)
    }
    obj_center = get_diff_vector(c(x1+x2, y1+y2)/2, CENTER)
    FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_2_VECTOR[[idx]]
    products <- numeric(length(FIX_VECTORS))
    for (i in seq_along(FIX_VECTORS)) {
      product = measure_dot_product(FIX_VECTORS[[i]], obj_center)
      products[i] = product
    }
    return(products)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)

MSC <- MSC %>%
  mutate(AVG_DIST_FIX_2_PREV_OBJ_DOT_PRODUCT = sapply(DIST_FIX_2_PREV_OBJ_DOT_PRODUCT, mean))


MSC$DIST_FIX_3_PREV_OBJ_DOT_PRODUCT <- Map(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_CURR_ABSENT_FIX_VECTORS$subjectnum == subj)
    if (length(idx) == 0) {
      warning("No match for subjectnum")
      return(NA)
    }
    if (length(idx) > 1) {
      warning("Multiple matches for subjectnum")
      return(NA)
    }
    obj_center = get_diff_vector(c(x1+x2, y1+y2)/2, CENTER)
    FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_3_VECTOR[[idx]]
    products <- numeric(length(FIX_VECTORS))
    for (i in seq_along(FIX_VECTORS)) {
      product = measure_dot_product(FIX_VECTORS[[i]], obj_center)
      products[i] = product
    }
    return(products)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)

MSC <- MSC %>%
  mutate(AVG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT = sapply(DIST_FIX_3_PREV_OBJ_DOT_PRODUCT, mean, na.rm = TRUE))

#####=== Image-Based Baseline for dot product vector difference===#####

# First determine number of fixations (columns in your matrices)
n_fixations <- 8

# Create a function that calculates dot products for a specific fixation number
calculate_dot_products <- function(fix_num) {
  Map(
    function(subjectnum, img_name, x1, y1, x2, y2) {
      idx <- which(IMAGES_DF$img_name == img_name)
      if (length(idx) == 0) {
        warning("No match for img_name: ", img_name)
        return(NA)
      }
      if (length(idx) > 1) {
        warning("Multiple matches for img_name: ", img_name)
        return(NA)
      }
      
      obj_center <- get_diff_vector(c(x1 + x2, y1 + y2)/2, CENTER)
      FIXATIONS_X <- IMAGES_DF$FIXATIONS_X[[idx]]
      FIXATIONS_Y <- IMAGES_DF$FIXATIONS_Y[[idx]]
      
      
      if (!(length(FIXATIONS_X) == length(FIXATIONS_Y))) {
        print(length(FIXATIONS_X))
        print(length(FIXATIONS_Y))
        stop("invalid lengths @ image-based baseline")
      }
      
      n_trials <- nrow(FIXATIONS_X)
      n_fixations_curr <- ncol(FIXATIONS_X)
      products <- numeric(n_trials)
      
      if(n_fixations_curr < fix_num){
        return(NA)
      }
      
      
      for (i in 1:n_trials) {
        if(IMAGES_DF$subjectnum[[idx]][i] == subjectnum){
          products[i] <- NA
          next;
        }
        fix_i_x <- FIXATIONS_X[i, fix_num]  # Get specific fixation for this trial
        fix_i_y <- FIXATIONS_Y[i, fix_num]
        fix_vector <- get_diff_vector(c(fix_i_x,fix_i_y), CENTER)
        products[i] <- measure_dot_product(fix_vector, obj_center)
      }
      
      return(products)
    },
    MSC$subjectnum,
    MSC$img_name,
    MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
    MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
  )
}

for (i in 1:n_fixations) {
  col_name <- paste0("CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
  avg_col_name <- paste0("AVG_CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
  # Calculate dot products for current fixation
  MSC[[col_name]] <- calculate_dot_products(i)
  
  # Calculate mean for each list element in the column
  MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
}



# ============================Analysis==================================== 

# Vector Difference Analysis

MSC <- MSC%>%mutate(
  PREV_OBJ_FIX_DOT_PRODUCT_1 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[1]),
  PREV_OBJ_FIX_DOT_PRODUCT_2 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[2]),
  PREV_OBJ_FIX_DOT_PRODUCT_3 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[3]),
  PREV_OBJ_FIX_DOT_PRODUCT_4 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[4]),
  PREV_OBJ_FIX_DOT_PRODUCT_5 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[5]),
  PREV_OBJ_FIX_DOT_PRODUCT_6 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[6]),
  PREV_OBJ_FIX_DOT_PRODUCT_7 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[7]),
  PREV_OBJ_FIX_DOT_PRODUCT_8 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[8])
)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_1, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_1_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)
mean(MSC$AVG_CURR_IMG_DIST_FIX_1_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_2, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_2_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)
mean(MSC$AVG_CURR_IMG_DIST_FIX_2_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_3, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)
mean(MSC$AVG_CURR_IMG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_4, na.rm = TRUE)
mean(MSC$AVG_CURR_IMG_DIST_FIX_4_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_5, na.rm = TRUE)
mean(MSC$AVG_CURR_IMG_DIST_FIX_5_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)



# OG Distance Comparison
MSC <- MSC%>%mutate(
  PREV_ORTHOGONAL_DISTANCE_1 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[1]),
  PREV_ORTHOGONAL_DISTANCE_2 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[2]),
  PREV_ORTHOGONAL_DISTANCE_3 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[3]),
  PREV_ORTHOGONAL_DISTANCE_4 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[4]),
  PREV_ORTHOGONAL_DISTANCE_5 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[5]),
  PREV_ORTHOGONAL_DISTANCE_6 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[6]),
  PREV_ORTHOGONAL_DISTANCE_7 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[7]),
  PREV_ORTHOGONAL_DISTANCE_8 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[8])
)


mean(MSC$PREV_ORTHOGONAL_DISTANCE_1, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_1_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_2, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_2_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_3, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_3_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_4, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_4_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_5, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_5_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_6, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_6_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_7, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_7_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)

mean(MSC$PREV_ORTHOGONAL_DISTANCE_8, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_8_PREV_OBJ_OG_DISTANCE, na.rm = TRUE)



MSC_DOT_PRODUCTS <- MSC[, c("PREV_OBJ_FIX_DOT_PRODUCT_1", "AVG_DIST_FIX_1_PREV_OBJ_DOT_PRODUCT","PREV_OBJ_FIX_DOT_PRODUCT_2", "AVG_DIST_FIX_2_PREV_OBJ_DOT_PRODUCT", "PREV_OBJ_FIX_DOT_PRODUCT_3", "AVG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT")]


# Draw subject fixation i distribution
vec_list <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_1_VECTOR[[settings$test_subject]]
coords <- do.call(rbind, vec_list)

# Base plot on XY plane
plot(
  coords[,1], coords[,2],
  xlab = "X", ylab = "Y",
  main = "Subject 1 Fix 1 Distribution",
  xlim = c(-180, 180),  # add some padding
  ylim = c(-180, 180),
  pch = 19, col = "blue", asp = 1
)
arrows(0, 0, coords[,1], coords[,2], length = 0.1, col = "red")


MSC <- MSC %>%
  mutate(PO_COORDS = pmap(
    list(PREV_OBJ_X1, PREV_OBJ_X2, PREV_OBJ_Y1, PREV_OBJ_Y2),
    ~ c((..1 + ..2)/2, (..3 + ..4)/2)
  ))

MSC <- MSC %>%
  mutate(
    RESULT = map(PO_COORDS, ~ get_diff_vector(.x, CENTER))
  )


# Draw subject fixation i distribution
vec_list <- MSC$RESULT
coords <- do.call(rbind, vec_list)

# Base plot on XY plane
plot(
  coords[,1], coords[,2],
  xlab = "X", ylab = "Y",
  main = "Subject 1 Images Coords Distribution",
  xlim = c(-380, 380),  # add some padding
  ylim = c(-250, 250),
  pch = 19, col = "blue", asp = 1
)
arrows(0, 0, coords[,1], coords[,2], length = 0.1, col = "red")
# MSC$diff <- -MSC$PREV_ORTHOGONAL_DISTANCE_2 + MSC$MEAN_OG_DISTANCE  

#write.csv(MSC[, c("AVG_CURR_IMG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT", "PREV_OBJ_FIX_DOT_PRODUCT_3", "AVG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT")], file = "dot_products.csv")
# write.csv(MSC[, c("PREV_ORTHOGONAL_DISTANCE_2", "CENTER_OG_DISTANCE", "MEAN_OG_DISTANCE", "AVG_FIX_OG_DISTANCE")], file = "MUTATED_MSC.csv")
# write_json(MSC, "__data.json", pretty = TRUE)



# TODO
# (N-2) effect
# Image-based baseline for distance
# Remove > 25% microwaves trials
# Manipulate filtration with possible combinations of present, absent
# Code Review
# Run all above on clocks
