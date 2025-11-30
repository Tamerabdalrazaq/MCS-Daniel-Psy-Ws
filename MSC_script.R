rm(list = ls()) 

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
library(purrr)
library(crayon)

## fix the priming of last fixation analysis and also the baseline analysis  

# Define script settings
settings <<- list(
  object = 'm',
  filter_bad_detections = TRUE,
  filter_inconsecutive_trial = TRUE,
  filter_curr_big_objects = 0.25,
  filter_prev_big_objects = 0.25,
  prev_condition = "present", 
  curr_condition = "absent",
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray"),
  fliter_prev_absent = TRUE,
  prev_tar_hit = TRUE,
  test_subject = 1,
  N_BASELINE_FIXATIONS = 10, 
  POL_EFFECT_ORDER = 1 # Will be set dynamically by the script
)

DISPLAY_W <<- 1280
DISPLAY_H <<- 800
CENTER <<- c(DISPLAY_W, DISPLAY_H)/2
PATH_TO_MICROWAVES <<- "C:/Users/danie/OneDrive/Documents/אוניברסיטה/דוקטורט/PhD/Projects/MCS_dataset - Tamer/MCS_dataset/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"
PATH_TO_CLOCKS <<- "C:/Users/danie/OneDrive/Documents/אוניברסיטה/דוקטורט/PhD/Projects/MCS_dataset - Tamer/MCS_dataset/MSC_PROCESSED_DETECTIONS_CLOCKS.csv"


run_script <- function() {
  
  if(settings$object == 'm'){
    MSC_PATH <<- PATH_TO_MICROWAVES
  } else {
    MSC_PATH <<- PATH_TO_CLOCKS
  }
    
  MSC <- read_csv(MSC_PATH)
  MSC_ORIGINAL <- MSC
  GROUPED_MSC <- MSC
  
  
  # ================================Dataset Preprocessing================================
  ##calculate accuracy 
  MSC = mutate(MSC, 
               Accuracy= ifelse(expected == LAST_BUTTON_PRESSED,1,0)
  )
  
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
  
  
  # Convert detection coordinates from relative to absolute
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
        accuracy = first(Accuracy),
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
  
  GROUPED_MSC <- MSC
  
  ## ============================Add New Features==================================== 
  prev_features <- c("OBJECT_X1", "OBJECT_X2", "OBJECT_Y1", "OBJECT_Y2", "condition", "TRIAL_INDEX","accuracy")
  get_prev_feature <- function (subjectnum, trial_index){
    idx <- which(MSC$subjectnum == subjectnum & MSC$TRIAL_INDEX == trial_index - settings$POL_EFFECT_ORDER)
    if (length(idx) == 0) {
      return(NA)
    }
    if (length(idx) > 1) {
      stop("2 matching trials @ get_prev_feature")
    }
    prev_trial = MSC[idx, ]
    prev_features_list = lapply(prev_features, function(feature) prev_trial[[feature]])
    return((prev_features_list))
  }
  MSC$PREV_FEATURES <- Map(get_prev_feature, MSC$subjectnum, MSC$TRIAL_INDEX)
  
  for (i in seq_along(prev_features)) {
    feature <- prev_features[[i]]  # get actual string
    prev_feature <- paste0("PREV_", feature)
    
    MSC[[prev_feature]] <- sapply(MSC$PREV_FEATURES, function(x) {
      if (is.list(x)) x[i] else NA
    })
  }
  
  
  MSC$PREV_OBJECT_X1 <- as.numeric(MSC$PREV_OBJECT_X1)
  MSC$PREV_OBJECT_Y1 <- as.numeric(MSC$PREV_OBJECT_Y1)
  MSC$PREV_OBJECT_X2 <- as.numeric(MSC$PREV_OBJECT_X2)
  MSC$PREV_OBJECT_Y2 <- as.numeric(MSC$PREV_OBJECT_Y2)
  MSC$PREV_FEATURES <- NULL
  
  # Images Baseline - explortatory 
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
      ) %>% 
      filter(
        (PREV_condition == "present" & !is.na(PREV_OBJECT_X1)) |
          (PREV_condition == "absent" & is.na(PREV_OBJECT_X1))
      )
  }
  
  
  ##!!##
  # Filter out big objects
  get_obj_ratio <- function(obj_x1, obj_y1, obj_x2, obj_y2, img_w, img_h){
    prev_obj_w = (obj_x2 - obj_x1) 
    prev_obj_h = (obj_y2 - obj_y1)
    prev_obj_area = prev_obj_w*prev_obj_h
    img_area = img_w * img_h
    (prev_obj_area / img_area)
  }
  
  
  
  # Filter out nonconsecutive trials
  if (settings$filter_inconsecutive_trial){
    MSC <- MSC %>%
      filter(
        !is.na(PREV_TRIAL_INDEX)
      )
  }
  
  
  if (settings$fliter_prev_absent){
    MSC <- MSC %>% filter(PREV_condition == settings$prev_condition)
    if(settings$curr_condition != FALSE){
      MSC <- MSC %>% filter(condition == settings$curr_condition)
    }
  }
  
  
  if (settings$filter_prev_big_objects) {
    MSC <- MSC %>% filter(  get_obj_ratio(PREV_OBJECT_X1, PREV_OBJECT_Y1,
                                          PREV_OBJECT_X2, PREV_OBJECT_Y2,
                                          im_w, im_h) < settings$filter_prev_big_objects)
  }
  
  if (settings$filter_curr_big_objects) {
    MSC <- MSC %>% filter(condition == 'absent' | get_obj_ratio(OBJECT_X1, OBJECT_Y1,
                                                                OBJECT_X2, OBJECT_Y2,
                                                                im_w, im_h) < settings$filter_curr_big_objects)
  }
  
  
  if (settings$prev_tar_hit) {
    MSC <- MSC %>% filter(PREV_accuracy == 1)
  }
  
  
  FIXATIONS_BASELINE_DF <- MSC
  
  
  if (settings$test_subject){
    MSC <- MSC %>% filter(subjectnum == settings$test_subject)
  }
  
  print("Done Filtering")
  
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
  
  vector_cols <- paste0("FIX_", 1:settings$N_BASELINE_FIXATIONS, "_VECTOR")
  PREV_CURR_ABSENT_FIX_VECTORS <- FIXATIONS_BASELINE_DF %>%
    group_by(subjectnum) %>%
    summarize(across(all_of(vector_cols), list), .groups = "drop")
  
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
    return(list(distances))
  }
  
  # Apply the function to the MSC dataset
  MSC <- MSC %>%
    mutate(
      PREV_ORTHOGONAL_DISTANCE = ifelse(
        !is.na(PREV_OBJECT_X1) & !is.na(PREV_OBJECT_Y1) & !is.na(PREV_OBJECT_X2) & !is.na(PREV_OBJECT_Y2),
        mapply(
          calc_orthogonal_distances,
          CURRENT_FIX_X, CURRENT_FIX_Y,
          PREV_OBJECT_X1, PREV_OBJECT_Y1, PREV_OBJECT_X2, PREV_OBJECT_Y2
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
    cos = sum(v1 * v2)
    return(cos)
  }
  
  measure_dot_product <- function(v1, v2) {
    return (calc_dot_product(v1, v2,
                             normalize = TRUE, max_norm = NULL))
  }
  

  ##dot product
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
    return(list(products))
  }
  
  MSC <- MSC %>%
    mutate(
      PREV_OBJ_FIX_DOT_PRODUCT = ifelse(
        !is.na(PREV_OBJECT_X1) & !is.na(PREV_OBJECT_Y1) & !is.na(PREV_OBJECT_X2) & !is.na(PREV_OBJECT_Y2),
        mapply(
          get_curr_fix_prev_obj_dot_product,
          CURRENT_FIX_X, CURRENT_FIX_Y,
          PREV_OBJECT_X1, PREV_OBJECT_Y1, PREV_OBJECT_X2, PREV_OBJECT_Y2
        ),
        NA  # Assign NA if any required column has NA
      )
      
    )
  # ============================Calculating Different Baselines==================================== 
  ### ==================================Subject-Based Distance==========================================
  
  # Calculate the OG distance between center and  the * previous * detection
  MSC$CENTER_OG_DISTANCE <- mapply(calc_orthogonal_distance,
                                   DISPLAY_W/2, DISPLAY_H/2,
                                   MSC$PREV_OBJECT_X1, MSC$PREV_OBJECT_Y1,
                                   MSC$PREV_OBJECT_X2, MSC$PREV_OBJECT_Y2)
  
  calc_dist_i_og_distance <- function(fix_num){
    Map(
      function(subj, x1, y1, x2, y2, current_x,current_y) {
        idx <- which(PREV_CURR_ABSENT_FIX_VECTORS$subjectnum == subj)
        if (length(idx) == 0) {
          warning("No match for subjectnum")
          return(NA)
        }
        if (length(idx) > 1) {
          warning("Multiple matches for subjectnum")
          return(NA)
        }
        ## remove 
        vector_name <- paste0("FIX_", fix_num, "_VECTOR")
        FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS[[vector_name]][[idx]]
        current_x = current_x[[idx]][fix_num]
        current_y = current_y[[idx]][fix_num]
        if (!is.na(current_x) & !is.na(current_y)){
        remove_vector = get_diff_vector(c(current_x, current_y), CENTER)
        removal_idx <- which.min(sapply(FIX_VECTORS, function(v) sqrt(sum((v - remove_vector)^2))))
        FIX_VECTORS = FIX_VECTORS[-removal_idx]
        }
        distances <- numeric(length(FIX_VECTORS))
        for (i in seq_along(FIX_VECTORS)) {
        # IMPORTANT! FIX_VECTORS is relative to the center and not absolute to the screen 
          curr_vector = FIX_VECTORS[[i]] + CENTER 
          x = curr_vector[1]
          y = curr_vector[2]
          dist =calc_orthogonal_distance(x, y, x1, y1, x2, y2)
          distances[i] = dist
        }
        return(distances)
      },
      MSC$subjectnum,
      MSC$PREV_OBJECT_X1, MSC$PREV_OBJECT_Y1,
      MSC$PREV_OBJECT_X2, MSC$PREV_OBJECT_Y2,
      MSC$CURRENT_FIX_X,MSC$CURRENT_FIX_Y
    )
  }
  
  for (i in 1:settings$N_BASELINE_FIXATIONS) {
    col_name <- paste0("DIST_FIX_", i, "_PREV_OBJ_OG_DISTANCE")
    avg_col_name <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_OG_DISTANCE")
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calc_dist_i_og_distance(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  
  
  ### ==================================Image-Based Distances==========================================
  # calculates dot products for a specific fixation number
  calculate_fix_i_angles <- function(fix_num) {
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
          products[i] <- calc_orthogonal_distance(fix_i_x, fix_i_y, x1, y1, x2, y2)
        }
        return(products)
      },
      MSC$subjectnum,
      MSC$img_name,
      MSC$PREV_OBJECT_X1, MSC$PREV_OBJECT_Y1,
      MSC$PREV_OBJECT_X2, MSC$PREV_OBJECT_Y2
    )
  }
  
  for (i in 1:settings$N_BASELINE_FIXATIONS) {
    col_name <- paste0("CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DISTANCE")
    avg_col_name <- paste0("AVG_CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DISTANCE")
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calculate_fix_i_angles(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  
  ### ==================================Subject-Based Angles==========================================
  calc_dist_i_dot_product <- function(fix_num){
    Map(
      function(subj, x1, y1, x2, y2,current_x,current_y) {
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
        vector_name <- paste0("FIX_", fix_num, "_VECTOR")
        FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS[[vector_name]][[idx]]
        current_x = current_x[[idx]][fix_num]
        current_y = current_y[[idx]][fix_num]
        if (!is.na(current_x) & !is.na(current_y)){
          remove_vector = get_diff_vector(c(current_x, current_y), CENTER)
          removal_idx <- which.min(sapply(FIX_VECTORS, function(v) sqrt(sum((v - remove_vector)^2))))
          FIX_VECTORS = FIX_VECTORS[-removal_idx]
        }
        products <- numeric(length(FIX_VECTORS))
        for (i in seq_along(FIX_VECTORS)) {
          product = measure_dot_product(FIX_VECTORS[[i]], obj_center)
          products[i] = product
        }
        return(products)
      },
      MSC$subjectnum,
      MSC$PREV_OBJECT_X1, MSC$PREV_OBJECT_Y1,
      MSC$PREV_OBJECT_X2, MSC$PREV_OBJECT_Y2,
      MSC$CURRENT_FIX_X,MSC$CURRENT_FIX_Y
    )
  }

  
  for (i in 1:settings$N_BASELINE_FIXATIONS) {
    col_name <- paste0("DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
    avg_col_name <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calc_dist_i_dot_product(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  ### ==================================Image-Based Angles==========================================
  # calculates dot products for a specific fixation number
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
      MSC$PREV_OBJECT_X1, MSC$PREV_OBJECT_Y1,
      MSC$PREV_OBJECT_X2, MSC$PREV_OBJECT_Y2
    )
  }
  
  for (i in 1:settings$N_BASELINE_FIXATIONS) {
    col_name <- paste0("CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
    avg_col_name <- paste0("AVG_CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calculate_dot_products(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  
  
  # ============================Analysis==================================== 
  MSC <- MSC%>%mutate(
    PREV_OBJ_FIX_DOT_PRODUCT_1 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[1]),
    PREV_OBJ_FIX_DOT_PRODUCT_2 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[2]),
    PREV_OBJ_FIX_DOT_PRODUCT_3 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[3]),
    PREV_OBJ_FIX_DOT_PRODUCT_4 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[4]),
    PREV_OBJ_FIX_DOT_PRODUCT_5 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[5]),
    PREV_OBJ_FIX_DOT_PRODUCT_6 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[6]),
    PREV_OBJ_FIX_DOT_PRODUCT_7 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[7]),
    PREV_OBJ_FIX_DOT_PRODUCT_8 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[8]),
    PREV_OBJ_FIX_DOT_PRODUCT_9 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[9]),
    PREV_OBJ_FIX_DOT_PRODUCT_10 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[10])
  )
  
  df <- data.frame(Measure = character(),
                   Subject_BL = character(),
                   Image_BL = character(),
                   stringsAsFactors = FALSE)
  for (i in 1:settings$N_BASELINE_FIXATIONS){
    measure <- paste0("PREV_OBJ_FIX_DOT_PRODUCT_", i)
    baseline_subjet <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
    baseline_image <- paste0("AVG_CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DOT_PRODUCT")
    mean_1 = (mean(MSC[[measure]], na.rm = TRUE))
    mean_2 = (mean(MSC[[baseline_subjet]], na.rm = TRUE))
    mean_3 = (mean(MSC[[baseline_image]], na.rm = TRUE))
    row = c(mean_1, mean_2, mean_3)
    df <- rbind(df, data.frame(Measure = row[1], Subject_BL = row[2], Image_BL = row[3], 
                               stringsAsFactors = FALSE))
  }
  df <- format(df, digits = 2, nsmall = 2)
  df[] <- lapply(df, as.numeric)
  df_angles <- df
  
  # OG Distance Comparison
  MSC <- MSC%>%mutate(
    PREV_ORTHOGONAL_DISTANCE_1 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[1]),
    PREV_ORTHOGONAL_DISTANCE_2 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[2]),
    PREV_ORTHOGONAL_DISTANCE_3 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[3]),
    PREV_ORTHOGONAL_DISTANCE_4 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[4]),
    PREV_ORTHOGONAL_DISTANCE_5 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[5]),
    PREV_ORTHOGONAL_DISTANCE_6 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[6]),
    PREV_ORTHOGONAL_DISTANCE_7 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[7]),
    PREV_ORTHOGONAL_DISTANCE_8 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[8]),
    PREV_ORTHOGONAL_DISTANCE_9 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[9]),
    PREV_ORTHOGONAL_DISTANCE_10 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[10])
  )
  
  df <- data.frame(Measure = character(),
                   Subject_BL = character(),
                   Image_BL = character(),
                   stringsAsFactors = FALSE)
  for (i in 1:settings$N_BASELINE_FIXATIONS){
    measure <- paste0("PREV_ORTHOGONAL_DISTANCE_", i)
    baseline_subjet <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_OG_DISTANCE")
    baseline_image <- paste0("AVG_CURR_IMG_DIST_FIX_", i, "_PREV_OBJ_DISTANCE")
    mean_1 = (mean(MSC[[measure]], na.rm = TRUE))
    mean_2 = (mean(MSC[[baseline_subjet]], na.rm = TRUE))
    mean_3 = (mean(MSC[[baseline_image]], na.rm = TRUE))
    row = c(mean_1, mean_2, mean_3)
    df <- rbind(df, data.frame(Measure = row[1], Subject_BL = row[2], Image_BL = row[3], stringsAsFactors = FALSE))
  }
  df <- format(df, digits = 2, nsmall = 2)
  df[] <- lapply(df, as.numeric)
  df_distances <- df
  return(list(df_angles=df_angles, df_distances=df_distances))
}

object = list('m', 'c')
results_dfs = list()
for(i in 1:2){
  settings$POL_EFFECT_ORDER = i
  for(obj in object){
    settings$object = obj
    results <- run_script()
    df_angles = results$df_angles
    df_distances = results$df_distances
    cat("\n ========Angles=======\n\n")
    print(df_angles)
    cat("\n ========Distance=======\n\n")
    print(df_distances)
    results_name = paste0(obj, "_N_", i, "_", cond, '_', 'angles')
    results_dfs[[results_name]] = df_angles
    results_name = paste0(obj, "_N_", i, "_", cond, '_', 'distances')
    results_dfs[[results_name]] = df_distances
    }
  }


plot_results <- function(df, name) {
  df$Index <- 1:nrow(df)
  df_long <- pivot_longer(df, cols = c(Measure, Subject_BL, Image_BL),
                          names_to = "Variable", values_to = "Value")
  title = paste0()
  ggplot(df_long, aes(x = Index, y = Value, color = Variable)) +
    geom_line() +
    geom_point() +
    labs(title = name, x = "Fixation", y = "Value") +
    theme_minimal()
}

for(name in names(results_dfs)){
  plt<- plot_results(results_dfs[[name]], name)
  
  ggsave(paste0("results/",name,".jpg"), plot = plt, width = 6, height = 4, dpi = 300)
}

## statistical analysis ----

#example: 
MSC = mutate(MSC, PoL_Effect_angles_1 = PREV_OBJ_FIX_DOT_PRODUCT_1 -AVG_DIST_FIX_1_PREV_OBJ_DOT_PRODUCT,
             PoL_Effect_distance_1 = PREV_ORTHOGONAL_DISTANCE_1 -AVG_DIST_FIX_1_PREV_OBJ_OG_DISTANCE
             )

library(brms)
model = PoL_Effect_saccade_1~ 1 + (1|participant)
get_prior(model, data = MSC)

library(brms)
myprior = prior(normal(0,0.3), class = intercept)
model_1 = brms::brm(model, 
                      prior = myprior,
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 2000, 
                      Iter = 2000,
                      #sample_prior = "only"
)
prior_summary(model_1)
summary(model_1)

## same for distance

## Exploratory: ----
## 3.	Dwell time as a function of PoL ----
MSC <- MSC%>%mutate(
  CURRENT_FIX_DURATION_1 = sapply(CURRENT_FIX_DURATION, function(x) x[1]),
  CURRENT_FIX_DURATION_2 = sapply(CURRENT_FIX_DURATION, function(x) x[2]),
  CURRENT_FIX_DURATION_3 = sapply(CURRENT_FIX_DURATION, function(x) x[3]),
  CURRENT_FIX_DURATION_4 = sapply(CURRENT_FIX_DURATION, function(x) x[4]),
  CURRENT_FIX_DURATION_5 = sapply(CURRENT_FIX_DURATION, function(x) x[5]),
  CURRENT_FIX_DURATION_6 = sapply(CURRENT_FIX_DURATION, function(x) x[6]),
  CURRENT_FIX_DURATION_7 = sapply(CURRENT_FIX_DURATION, function(x) x[7]),
  CURRENT_FIX_DURATION_8 = sapply(CURRENT_FIX_DURATION, function(x) x[8]),
  CURRENT_FIX_DURATION_9 = sapply(CURRENT_FIX_DURATION, function(x) x[9]),
  CURRENT_FIX_DURATION_10 = sapply(CURRENT_FIX_DURATION, function(x) x[10])
)

model = CURRENT_FIX_DURATION_1~ PoL_Effect_distance_1 + (PoL_Effect_distance_1|participant) ## if works, if not remove random slopes
get_prior(model, data = MSC)

library(brms)
model_2 = brms::brm(model, 
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 2000, 
                      Iter = 2000,
                      #sample_prior = "only"
)
prior_summary(model_2)
summary(model_2)

## 4.	Priming of last fixation vs. priming previous target----
GROUPED_MSC$last_fixation_x = sapply(GROUPED_MSC$CURRENT_FIX_X, function(x) {tail(x,n = 1)})
GROUPED_MSC$last_fixation_y = sapply(GROUPED_MSC$CURRENT_FIX_Y, function(x) {tail(x,n = 1)})
GROUPED_MSC$last_fixation_duration = sapply(GROUPED_MSC$CURRENT_FIX_DURATION, function(x) {tail(x,n = 1)})

prev_features <- c("TRIAL_INDEX","accuracy","condition","last_fixation_x","last_fixation_y","last_fixation_duration")
get_prev_feature <- function (subjectnum, trial_index){
  idx <- which(GROUPED_MSC$subjectnum == subjectnum & GROUPED_MSC$TRIAL_INDEX == trial_index - settings$POL_EFFECT_ORDER)
  if (length(idx) == 0) {
    return(NA)
  }
  if (length(idx) > 1) {
    stop("2 matching trials @ get_prev_feature")
  }
  prev_trial = GROUPED_MSC[idx, ]
  prev_features_list = lapply(prev_features, function(feature) prev_trial[[feature]])
  return((prev_features_list))
}
GROUPED_MSC$PREV_FEATURES <- Map(get_prev_feature, GROUPED_MSC$subjectnum, GROUPED_MSC$TRIAL_INDEX)

for (i in seq_along(prev_features)) {
  feature <- prev_features[[i]]  # get actual string
  prev_feature <- paste0("PREV_", feature)
  
  GROUPED_MSC[[prev_feature]] <- sapply(GROUPED_MSC$PREV_FEATURES, function(x) {
    if (is.list(x)) x[i] else NA
  })
}


GROUPED_MSC$PREV_last_fixation_x <- as.numeric(GROUPED_MSC$PREV_last_fixation_x)
GROUPED_MSC$PREV_last_fixation_y <- as.numeric(GROUPED_MSC$PREV_last_fixation_y)
GROUPED_MSC$PREV_last_fixation_duration <- as.numeric(GROUPED_MSC$PREV_last_fixation_duration)

# Filter out nonconsecutive trials
GROUPED_MSC <- GROUPED_MSC %>%
    filter(
      !is.na(PREV_TRIAL_INDEX)
    )

GROUPED_MSC <- GROUPED_MSC %>% filter(PREV_condition == "absent")
GROUPED_MSC <- GROUPED_MSC %>% filter(condition == "absent")

if (settings$prev_tar_hit) {
  GROUPED_MSC <- GROUPED_MSC %>% filter(PREV_accuracy == 1)
}

FIXATIONS_BASELINE_DF <- GROUPED_MSC

if (settings$test_subject){
  GROUPED_MSC <- GROUPED_MSC %>% filter(subjectnum == settings$test_subject)
}

print("Done Filtering")

#Subject Fixations Vectors For Baselines

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

vector_cols <- paste0("FIX_", 1:settings$N_BASELINE_FIXATIONS, "_VECTOR")
PREV_ABSENT_ABSENT_FIX_VECTORS <- FIXATIONS_BASELINE_DF %>%
  group_by(subjectnum) %>%
  summarize(across(all_of(vector_cols), list), .groups = "drop")

# measures
#distances from last fixation

euc_l2 <- function(x1, y1, x2, y2) {
  # Calculate the Euclidean distance
  distance <- sqrt((x2 - x1)^2 + (y2 - y1)^2)
  return(distance)
}

# Accepts lists X,Y and scalars x1,y1,x2,y2 and returns a new list by applying calc_orthogonal_distance(x, y, x1, y1, x2, y2)
calc_distances <- function(X, Y, x1, y1) {
  # Check if lengths of X and Y are the same
  if (length(X) != length(Y)) {
    stop("X and Y must have the same length")
  }
  
  # Initialize an empty vector to store results
  distances <- numeric(length(X))
  
  # Iterate over the elements of X and Y
  for (i in seq_along(X)) {
    distances[i] <- euc_l2(X[i], Y[i], x1, y1)
  }
  
  # Return the list of distances
  return(list(distances))
}

# Apply the function to the MSC dataset
GROUPED_MSC <- GROUPED_MSC %>%
  mutate(
    PREV_DISTANCE_last_fixation = ifelse(
      !is.na(PREV_last_fixation_x) & !is.na(PREV_last_fixation_y),
      mapply(
        calc_distances,
        CURRENT_FIX_X, CURRENT_FIX_Y,
        PREV_last_fixation_x, PREV_last_fixation_y
      ),
      NA  # Assign NA if any required column has NA
    )
    
  )


#### ============================Angles==================================== 

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
  cos = sum(v1 * v2)
  return(cos)
}

measure_dot_product <- function(v1, v2) {
  return (calc_dot_product(v1, v2,
                           normalize = TRUE, max_norm = NULL))
}


##dot product
get_curr_fix_prev_fix_dot_product <- function(FIX_X, FIX_Y, x1, y1){
  products <- numeric(length(FIX_X))
  obj_center <- c(x1, y1)
  obj_vector <- get_diff_vector(obj_center, CENTER)
  
  # Iterate over the elements of X and Y
  for (i in seq_along(FIX_X)) {
    fix_vector <- get_diff_vector(c(FIX_X[i], FIX_Y[i]), CENTER)
    products[i] <- measure_dot_product(fix_vector, obj_vector)
  }
  
  # Return the list of distances
  return(list(products))
}

GROUPED_MSC <- GROUPED_MSC %>%
  mutate(
    PREV_OBJ_FIX_DOT_PRODUCT = ifelse(
      !is.na(PREV_last_fixation_x) & !is.na(PREV_last_fixation_y),
      mapply(
        get_curr_fix_prev_fix_dot_product,
        CURRENT_FIX_X, CURRENT_FIX_Y,
        PREV_last_fixation_x, PREV_last_fixation_y
      ),
      NA  # Assign NA if any required column has NA
    )
    
  )

#### ============================Calculating Different Baselines==================================== 
#### ==================================Subject-Based Distance==========================================


calc_dist_i_og_distance <- function(fix_num){
  Map(
    function(subj, x1, y1, current_x,current_y) {
      idx <- which(PREV_ABSENT_ABSENT_FIX_VECTORS$subjectnum == subj)
      if (length(idx) == 0) {
        warning("No match for subjectnum")
        return(NA)
      }
      if (length(idx) > 1) {
        warning("Multiple matches for subjectnum")
        return(NA)
      }
      ## remove 
      vector_name <- paste0("FIX_", fix_num, "_VECTOR")
      FIX_VECTORS <- PREV_ABSENT_ABSENT_FIX_VECTORS[[vector_name]][[idx]]
      current_x = current_x[[idx]][fix_num]
      current_y = current_y[[idx]][fix_num]
      if (!is.na(current_x) & !is.na(current_y)){
        remove_vector = get_diff_vector(c(current_x, current_y), CENTER)
        removal_idx <- which.min(sapply(FIX_VECTORS, function(v) sqrt(sum((v - remove_vector)^2))))
        FIX_VECTORS = FIX_VECTORS[-removal_idx]
      }
      distances <- numeric(length(FIX_VECTORS))
      for (i in seq_along(FIX_VECTORS)) {
        # IMPORTANT! FIX_VECTORS is relative to the center and not absolute to the screen 
        curr_vector = FIX_VECTORS[[i]] + CENTER 
        x = curr_vector[1]
        y = curr_vector[2]
        dist =calc_distances(x, y, x1, y1)
        distances[i] = dist
      }
      return(distances)
    },
    GROUPED_MSC$subjectnum,
    GROUPED_MSC$PREV_last_fixation_x, GROUPED_MSC$PREV_last_fixation_y,
    GROUPED_MSC$CURRENT_FIX_X,GROUPED_MSC$CURRENT_FIX_Y
  )
}

for (i in 1:settings$N_BASELINE_FIXATIONS) {
  col_name <- paste0("DIST_FIX_", i, "_PREV_fix_OG_DISTANCE")
  avg_col_name <- paste0("AVG_DIST_FIX_", i, "_PREV_fix_OG_DISTANCE")
  # Calculate dot products for current fixation
  GROUPED_MSC[[col_name]] <- calc_dist_i_og_distance(i)
  
  # Calculate mean for each list element in the column
  GROUPED_MSC[[avg_col_name]] <- sapply(GROUPED_MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
}

#### ==================================Subject-Based Angles==========================================

calc_dist_i_dot_product <- function(fix_num){
  Map(
    function(subj, x1, y1,current_x,current_y) {
      idx <- which(PREV_ABSENT_ABSENT_FIX_VECTORS$subjectnum == subj)
      if (length(idx) == 0) {
        warning("No match for subjectnum")
        return(NA)
      }
      if (length(idx) > 1) {
        warning("Multiple matches for subjectnum")
        return(NA)
      }
      obj_center = get_diff_vector(c(x1, y1), CENTER)
      vector_name <- paste0("FIX_", fix_num, "_VECTOR")
      FIX_VECTORS <- PREV_ABSENT_ABSENT_FIX_VECTORS[[vector_name]][[idx]]
      current_x = current_x[[idx]][fix_num]
      current_y = current_y[[idx]][fix_num]
      if (!is.na(current_x) & !is.na(current_y)){
        remove_vector = get_diff_vector(c(current_x, current_y), CENTER)
        removal_idx <- which.min(sapply(FIX_VECTORS, function(v) sqrt(sum((v - remove_vector)^2))))
        FIX_VECTORS = FIX_VECTORS[-removal_idx]
      }
      products <- numeric(length(FIX_VECTORS))
      for (i in seq_along(FIX_VECTORS)) {
        product = measure_dot_product(FIX_VECTORS[[i]], obj_center)
        products[i] = product
      }
      return(products)
    },
    GROUPED_MSC$subjectnum,
    GROUPED_MSC$PREV_last_fixation_x, GROUPED_MSC$PREV_last_fixation_y,
    GROUPED_MSC$CURRENT_FIX_X,GROUPED_MSC$CURRENT_FIX_Y
  )
}

for (i in 1:settings$N_BASELINE_FIXATIONS) {
  col_name <- paste0("DIST_FIX_", i, "_PREV_OBJ_fix_PRODUCT")
  avg_col_name <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_fix_PRODUCT")
  # Calculate dot products for current fixation
  GROUPED_MSC[[col_name]] <- calc_dist_i_dot_product(i)
  
  # Calculate mean for each list element in the column
  GROUPED_MSC[[avg_col_name]] <- sapply(GROUPED_MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
}

#### ============================Analysis==================================== 
GROUPED_MSC <- GROUPED_MSC%>%mutate(
  PREV_OBJ_FIX_DOT_PRODUCT_1 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[1]),
  PREV_OBJ_FIX_DOT_PRODUCT_2 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[2]),
  PREV_OBJ_FIX_DOT_PRODUCT_3 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[3]),
  PREV_OBJ_FIX_DOT_PRODUCT_4 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[4]),
  PREV_OBJ_FIX_DOT_PRODUCT_5 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[5]),
  PREV_OBJ_FIX_DOT_PRODUCT_6 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[6]),
  PREV_OBJ_FIX_DOT_PRODUCT_7 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[7]),
  PREV_OBJ_FIX_DOT_PRODUCT_8 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[8]),
  PREV_OBJ_FIX_DOT_PRODUCT_9 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[9]),
  PREV_OBJ_FIX_DOT_PRODUCT_10 = sapply(PREV_OBJ_FIX_DOT_PRODUCT, function(x) x[10])
)

df <- data.frame(Measure = character(),
                 Subject_BL = character(),
                 stringsAsFactors = FALSE)
for (i in 1:settings$N_BASELINE_FIXATIONS){
  measure <- paste0("PREV_OBJ_FIX_DOT_PRODUCT_", i)
  baseline_subjet <- paste0("AVG_DIST_FIX_", i, "_PREV_OBJ_fix_PRODUCT")
  mean_1 = (mean(GROUPED_MSC[[measure]], na.rm = TRUE))
  mean_2 = (mean(GROUPED_MSC[[baseline_subjet]], na.rm = TRUE))
  row = c(mean_1, mean_2)
  df <- rbind(df, data.frame(Measure = row[1], Subject_BL = row[2], 
                             stringsAsFactors = FALSE))
}
df <- format(df, digits = 2, nsmall = 2)
df[] <- lapply(df, as.numeric)
df_angles <- df

##distance
GROUPED_MSC <- GROUPED_MSC%>%mutate(
  PREV_DISTANCE_last_fixation_1 = sapply(PREV_DISTANCE_last_fixation, function(x) x[1]),
  PREV_DISTANCE_last_fixation_2 = sapply(PREV_DISTANCE_last_fixation, function(x) x[2]),
  PREV_DISTANCE_last_fixation_3 = sapply(PREV_DISTANCE_last_fixation, function(x) x[3]),
  PREV_DISTANCE_last_fixation_4 = sapply(PREV_DISTANCE_last_fixation, function(x) x[4]),
  PREV_DISTANCE_last_fixation_5 = sapply(PREV_DISTANCE_last_fixation, function(x) x[5]),
  PREV_DISTANCE_last_fixation_6 = sapply(PREV_DISTANCE_last_fixation, function(x) x[6]),
  PREV_DISTANCE_last_fixation_7 = sapply(PREV_DISTANCE_last_fixation, function(x) x[7]),
  PREV_DISTANCE_last_fixation_8 = sapply(PREV_DISTANCE_last_fixation, function(x) x[8]),
  PREV_DISTANCE_last_fixation_9 = sapply(PREV_DISTANCE_last_fixation, function(x) x[9]),
  PREV_DISTANCE_last_fixation_10 = sapply(PREV_DISTANCE_last_fixation, function(x) x[10])
)

df <- data.frame(Measure = character(),
                 Subject_BL = character(),
                 stringsAsFactors = FALSE)
for (i in 1:settings$N_BASELINE_FIXATIONS){
  measure <- paste0("PREV_DISTANCE_last_fixation_", i)
  baseline_subjet <- paste0("AVG_DIST_FIX_", i, "_PREV_fix_OG_DISTANCE")
  mean_1 = (mean(GROUPED_MSC[[measure]], na.rm = TRUE))
  mean_2 = (mean(GROUPED_MSC[[baseline_subjet]], na.rm = TRUE))
  row = c(mean_1, mean_2)
  df <- rbind(df, data.frame(Measure = row[1], Subject_BL = row[2], stringsAsFactors = FALSE))
}
df <- format(df, digits = 2, nsmall = 2)
df[] <- lapply(df, as.numeric)
df_distances <- df


#remove example: 
GROUPED_MSC = mutate(GROUPED_MSC, last_fixation_Effect_distance_1=   PREV_DISTANCE_last_fixation_1 - AVG_DIST_FIX_1_PREV_fix_OG_DISTANCE, 
            last_fixation_Effect_angles_1=  PREV_OBJ_FIX_DOT_PRODUCT_1 - AVG_DIST_FIX_1_PREV_OBJ_fix_PRODUCT
)

sum_MSC=MSC%>%group_by(subjectnum)%>%summarize (PoL_Effect_angles_1 = mean(PoL_Effect_angles_1),PoL_Effect_distance_1 = mean(PoL_Effect_distance_1))
sum_groupd_MSC=GROUPED_MSC%>%group_by(subjectnum)%>%summarize (last_fixation_Effect_distance_1 = mean(last_fixation_Effect_distance_1),last_fixation_Effect_angles_1 = mean(last_fixation_Effect_angles_1))

joint_MSC = left_join(sum_MSC,sum_groupd_MSC,
                      by = c("subjectnum"))
joint_MSC = mutate(joint_MSC, unique_PoL_distance = PoL_Effect_distance_1- last_fixation_Effect_distance_1,
                   PoL_Effect_angles_1 - last_fixation_Effect_angles_1)

library(brms)
model = unique_PoL_distance~ 1 + (1|participant)
get_prior(model, data = MSC)

model_3 = brms::brm(model, 
                    prior = myprior,
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 2000, 
                    Iter = 2000,
                    #sample_prior = "only"
)
prior_summary(model_3)
summary(model_3)

#5.	Time spent near the previous target  ----

proxmity_fun = function(x_vec, y_vec, xmin, xmax, ymin, ymax,im_w,im_h) {
  
  # center of the object
  cx <- (xmin + xmax) / 2
  cy <- (ymin + ymax) / 2
  
  roi_w <- 0.5 * im_w
  roi_h <- 0.5 * im_h
  
  xmin_p <- cx - roi_w / 2
  xmax_p <- cx + roi_w / 2
  ymin_p <- cy - roi_h / 2
  ymax_p <- cy + roi_h / 2
  
  # classify each fixation
  inside <- (x_vec >= xmin_p & x_vec <= xmax_p &
               y_vec >= ymin_p & y_vec <= ymax_p)
  
  # convert TRUE/FALSE → 1/0
  return(as.integer(inside))
}

MSC <- MSC %>%
  mutate(
    proximity_vec = mapply(
      {proxmity_fun},
      CURRENT_FIX_X,
      CURRENT_FIX_Y,
      PREV_OBJECT_X1,
      PREV_OBJECT_X2,
      PREV_OBJECT_Y1,
      PREV_OBJECT_Y2,im_w,im_h)
  )


MSC <- MSC %>%
  mutate(
    total_time_inside = mapply(function(dur, prox){
      sum(dur[prox == 1], na.rm = TRUE)
    }, CURRENT_FIX_DURATION, proximity_vec),
    total_duration_time = unlist(lapply(CURRENT_FIX_DURATION,sum)),
    prop_inside = total_time_inside/total_duration_time,
    obj_ratio = get_obj_ratio(PREV_OBJECT_X1, PREV_OBJECT_Y1,
                              PREV_OBJECT_X2, PREV_OBJECT_Y2,
                              im_w, im_h),
    roi_ratio = mapply(function(xmin, xmax, ymin, ymax) {
      # center of the object
      cx <- (xmin + xmax) / 2
      cy <- (ymin + ymax) / 2
      
      # ROI dimensions: 25% of screen
      roi_w <- 0.5 * im_w
      roi_h <- 0.5 * im_h
      
      # ROI box
      xmin_p <- cx - roi_w / 2
      xmax_p <- cx + roi_w / 2
      ymin_p <- cy - roi_h / 2
      ymax_p <- cy + roi_h / 2  
      
      return(get_obj_ratio(xmin_p, ymin_p,
                    xmax_p, ymax_p,
                    im_w, im_h))
      
    }, PREV_OBJECT_X1, PREV_OBJECT_X2,
    PREV_OBJECT_Y1, PREV_OBJECT_Y2)
  )
table(MSC$obj_ratio)
table(MSC$roi_ratio)

mean(MSC$prop_inside)

##baseline calculation 
random_prev_object = MSC

set.seed(15)
nrep = 10000
vec_df = as.data.frame(matrix(rep(NA,nrep*nrow(MSC)), ncol =nrep , nrow = nrow(MSC))) 
for (i in 1:nrep){
  ## match for eccentricity 

  random_prev_object = mutate(random_prev_object,
                              obj_width = PREV_OBJECT_X2-PREV_OBJECT_X1,
                              obj_length = PREV_OBJECT_Y2-PREV_OBJECT_Y1,
                              angle = runif(nrow(MSC), 0, 2 * pi),
                              center_x = CENTER[1] + CENTER_OG_DISTANCE * cos(angle),
                              center_y = CENTER[2] + CENTER_OG_DISTANCE * sin(angle),
                              x1_raw = center_x - obj_width  / 2,
                              y1_raw = center_y - obj_length / 2,
                              PREV_OBJECT_X1_R = x1_raw,
                              PREV_OBJECT_X2_R = PREV_OBJECT_X1_R+obj_width,
                              PREV_OBJECT_Y1_R = y1_raw,
                              PREV_OBJECT_Y2_R = PREV_OBJECT_Y1_R+obj_length
  )
  random_prev_object <- random_prev_object %>%
    mutate(
      proximity_vec = mapply(
        {proxmity_fun},
        CURRENT_FIX_X,
        CURRENT_FIX_Y,
        PREV_OBJECT_X1_R,
        PREV_OBJECT_X2_R,
        PREV_OBJECT_Y1_R,
        PREV_OBJECT_Y2_R,im_w,im_h)
    )
  
  random_prev_object <- random_prev_object %>%
    mutate(
      total_time_inside = mapply(function(dur, prox){
        sum(dur[prox == 1], na.rm = TRUE)
      }, CURRENT_FIX_DURATION, proximity_vec),
      total_duration_time = unlist(lapply(CURRENT_FIX_DURATION,sum)),
      prop_inside = total_time_inside/total_duration_time
    )
  vec_df[,i] = random_prev_object$prop_inside
}

MSC = mutate(MSC,
             random_prop_inside = rowMeans(vec_df,na.rm = TRUE),
             time_spent_pol = prop_inside - random_prop_inside
             )

model = time_spent_pol~ 1 + (1|participant) ## if works, if not remove random slopes
get_prior(model, data = MSC)

library(brms)
model_4 = brms::brm(model, 
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 2000, 
                    Iter = 2000,
                    #sample_prior = "only"
)
prior_summary(model_4)
summary(model_4)


