rm(list = ls()) 

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tidyverse)
library(purrr)
library(crayon)

# Define script settings
settings <<- list(
  object = 'm', # 'm' for microwaves, 'c' for clocks
  filter_bad_detections = TRUE, 
  filter_inconsecutive_trial = TRUE,
  filter_curr_big_objects = 0.25, #if the current target object is bigger than % of the screen (between 0-1)
  filter_prev_big_objects = 0.25, #if the previous target object is bigger than % of the screen (between 0-1)
  prev_condition = "present", 
  curr_condition = "absent",
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray"),
  fliter_prev_absent = TRUE,
  prev_tar_hit = TRUE, ## previous response was correct
  test_subject = 1, ## number of subject to remove as the analysis was pre-registered based on him
  N_BASELINE_FIXATIONS = 10, ## number of relevant fixations  
  POL_EFFECT_ORDER = 3 # n-i effect
)

##other parameters
DISPLAY_W <<- 1280
DISPLAY_H <<- 800
CENTER <<- c(DISPLAY_W, DISPLAY_H)/2
PATH_TO_MICROWAVES <<- "~/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"
PATH_TO_CLOCKS <<- "~/MSC_PROCESSED_DETECTIONS_CLOCKS.csv"


#run_script <- function() {
  
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
  
  GROUPED_MSC <- MSC ## save for later use
  
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
  
  # Images Baseline - exploratory 
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
  
  ## target detection accuracy
  MSC%>%group_by(condition)%>%summarize(sum(is.na(OBJECT_X1)))
  #MSC%>%group_by(condition)%>%summarize(mean(confidence, na.rm = TRUE))
  MSC%>%group_by(condition)%>%summarize(n())
  
  ## participants exclusion
  source("~/Subject_outliers_Detector_function.R")
  Sub_outlier_Detector(MSC,subjectnum,accuracy,Mad=FALSE)
  MSC%>%group_by(subjectnum)%>%summarize(accuracy= mean(accuracy))%>%summarize(accuracy= mean(accuracy))
  MSC%>%group_by(subjectnum)%>%summarize(accuracy= mean(accuracy))%>%summarize(sd=sd(accuracy))
  
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
  
  ## filter out previous trial absent
  if (settings$fliter_prev_absent){
    MSC <- MSC %>% filter(PREV_condition == settings$prev_condition)
    if(settings$curr_condition != FALSE){
      MSC <- MSC %>% filter(condition == settings$curr_condition)
    }
  }
  
  ## filter out previous big targets
  if (settings$filter_prev_big_objects) {
    MSC <- MSC %>% filter(get_obj_ratio(PREV_OBJECT_X1, PREV_OBJECT_Y1,
                                          PREV_OBJECT_X2, PREV_OBJECT_Y2,
                                          im_w, im_h) < settings$filter_prev_big_objects)
  }
  ## filter out previous big current targets (only relevant for present-> present anlayses)  
  if (settings$filter_curr_big_objects) {
    MSC <- MSC %>% filter(condition == 'absent' | get_obj_ratio(OBJECT_X1, OBJECT_Y1,
                                                                OBJECT_X2, OBJECT_Y2,
                                                                im_w, im_h) < settings$filter_curr_big_objects)
  }
  
  ## filter previous incorrect trial (i.e. n-1, n-2 according to the relevant comparison )
  if (settings$prev_tar_hit) {
    MSC <- MSC %>% filter(PREV_accuracy == 1)
  }
  
  ##save for later
  FIXATIONS_BASELINE_DF <- MSC
  
  ## filter the participants from the pre-registeration
  if (settings$test_subject){
    MSC <- MSC %>% filter(subjectnum != settings$test_subject)
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
  ## Calculates the saccadic direction - cosine similarity
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
  
  # Calculate the distance between center and  the * previous * detection (subject-based)
  MSC$CENTER_OG_DISTANCE <- mapply(calc_orthogonal_distance,
                                   DISPLAY_W/2, DISPLAY_H/2,
                                   MSC$PREV_OBJECT_X1, MSC$PREV_OBJECT_Y1,
                                   MSC$PREV_OBJECT_X2, MSC$PREV_OBJECT_Y2)
  
  calc_dist_i_og_distance <- function(fix_num){
    Map(
      function(subj,x1, y1, x2, y2, current_x,current_y) {
        idx <- which(PREV_CURR_ABSENT_FIX_VECTORS$subjectnum == subj)
      #  idx_current = which(MSC$subjectnum== subj)
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
        current_x = current_x[fix_num]
        current_y = current_y[fix_num]
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
    col_name <- paste0("DIST_PREV_OBJ_OG_DISTANCE_FIX_", i)
    avg_col_name <- paste0("AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_", i)
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calc_dist_i_og_distance(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  
  
  ### ==================================Image-Based Distances==========================================
  # calculates distance for a specific fixation number
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
    col_name <- paste0("CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_", i)
    avg_col_name <- paste0("AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_", i)
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
        current_x = current_x[fix_num]
        current_y = current_y[fix_num]
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
    col_name <- paste0("DIST_PREV_OBJ_DOT_PRODUCT_FIX_", i)
    avg_col_name <- paste0("AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_", i)
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
    col_name <- paste0("CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_", i)
    avg_col_name <- paste0("AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_", i)
    # Calculate dot products for current fixation
    MSC[[col_name]] <- calculate_dot_products(i)
    
    # Calculate mean for each list element in the column
    MSC[[avg_col_name]] <- sapply(MSC[[col_name]], function(x) mean(unlist(x), na.rm = TRUE))
  }
  
  
  
  # ============================Analysis==================================== 
  MSC <- MSC%>%mutate( ## for each fixation, calculate the direction to previous target location (cosine similarity)
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
  
  df_angles = MSC %>%
    group_by(subjectnum) %>%
    summarize(
      across(
        c(starts_with("PREV_OBJ_FIX_DOT_PRODUCT_"),contains("_PREV_OBJ_DOT_PRODUCT")),
        mean,
        na.rm = TRUE
      )
    )%>%select(subjectnum,starts_with("PREV"), starts_with("AVG"))
  
  res = df_angles%>%pivot_longer(PREV_OBJ_FIX_DOT_PRODUCT_1:AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_10,names_to = "variable", values_to = "value")
  res$fixation_num = as.numeric(str_extract(res$variable, "\\d+$"))
  res$variable = sub("_[0-9]+$", "", res$variable)
  
  final_angles = res%>%group_by(fixation_num,variable)%>%summarize(value = mean(value,na.rm = TRUE))%>%pivot_wider(names_from = variable, values_from = value)%>%
    rename("Img_baseline" = "AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX", "subject_baseline"= "AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX",
           "measure"  = "PREV_OBJ_FIX_DOT_PRODUCT") 
                                                                                                                                                         
    
  # Calculate the distance for each fixation
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
  
  
  df_distance = MSC %>%
    group_by(subjectnum) %>%
    summarize(
      across(
        c(starts_with("PREV_ORTHOGONAL"),contains("DISTANCE")),
        mean,
        na.rm = TRUE
      )
    )%>%select(subjectnum,starts_with("PREV"), starts_with("AVG"), - "CENTER_OG_DISTANCE")
  
  res_dis = df_distance%>%pivot_longer(PREV_ORTHOGONAL_DISTANCE:AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_10,names_to = "variable", values_to = "value")
  res_dis$fixation_num = as.numeric(str_extract(res_dis$variable, "\\d+$"))
  res_dis$variable = sub("_[0-9]+$", "", res_dis$variable)
  
  final_distance = res_dis%>%filter(!is.na(fixation_num))%>%group_by(fixation_num,variable)%>%summarize(value = mean(value,na.rm = TRUE))%>%pivot_wider(names_from = variable, values_from = value)%>%
    rename("Img_baseline" = "AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX", "subject_baseline"= "AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX",
           "measure"  = "PREV_ORTHOGONAL_DISTANCE") 
  
  
  ## statistical analysis ----
  ## PoL effect direction and distances
  MSC = mutate(MSC, PoL_Effect_angles_1 = PREV_OBJ_FIX_DOT_PRODUCT_1 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_1, 
               PoL_Effect_angles_2 = PREV_OBJ_FIX_DOT_PRODUCT_2 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_2,## first saccade
               PoL_Effect_angles_3 = PREV_OBJ_FIX_DOT_PRODUCT_3 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_3,
               PoL_Effect_angles_4 = PREV_OBJ_FIX_DOT_PRODUCT_4 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_4,
               PoL_Effect_angles_5 = PREV_OBJ_FIX_DOT_PRODUCT_5 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_5,
               PoL_Effect_angles_6 = PREV_OBJ_FIX_DOT_PRODUCT_6 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_6,
               PoL_Effect_angles_7 = PREV_OBJ_FIX_DOT_PRODUCT_7 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_7,
               PoL_Effect_angles_8 = PREV_OBJ_FIX_DOT_PRODUCT_8 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_8,
               PoL_Effect_angles_9 = PREV_OBJ_FIX_DOT_PRODUCT_9 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_9,
               PoL_Effect_angles_10 = PREV_OBJ_FIX_DOT_PRODUCT_10 -AVG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_10,
               PoL_Effect_distance_1 = PREV_ORTHOGONAL_DISTANCE_1 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_1,## first saccade
               PoL_Effect_distance_2 = PREV_ORTHOGONAL_DISTANCE_2 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_2,
               PoL_Effect_distance_3 = PREV_ORTHOGONAL_DISTANCE_3 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_3,             
               PoL_Effect_distance_4 = PREV_ORTHOGONAL_DISTANCE_4 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_4,
               PoL_Effect_distance_5 = PREV_ORTHOGONAL_DISTANCE_5 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_5,             
               PoL_Effect_distance_6 = PREV_ORTHOGONAL_DISTANCE_6 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_6,
               PoL_Effect_distance_7 = PREV_ORTHOGONAL_DISTANCE_7 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_7,
               PoL_Effect_distance_8 = PREV_ORTHOGONAL_DISTANCE_8 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_8,             
               PoL_Effect_distance_9 = PREV_ORTHOGONAL_DISTANCE_9 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_9,
               PoL_Effect_distance_10 = PREV_ORTHOGONAL_DISTANCE_10 -AVG_DIST_PREV_OBJ_OG_DISTANCE_FIX_10             
               
  )
  
  write_path = "~/models_181026/"
  
  library(brms)
  myprior = prior(normal(0,0.3), class = "Intercept")
  model = PoL_Effect_angles_1~ 1 + (1|subjectnum)
  ## the first model is not interesting, as it reflects fixation to the fixation point
  # model_1 = brms::brm(model, 
  #                     prior = myprior,
  #                     data = MSC,
  #                     chain = 4,
  #                     cores = 4,
  #                     backend = "cmdstan",
  #                     warmup = 1000, 
  #                     iter = 4000,
  #                     silent =0,
  #                     thin = 2,
  #                     file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_1","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                     #sample_prior = "only"
  # )
  # prior_summary(model_1)
  # summary(model_1)
  library(bayestestR)
 # describe_posterior(model_1)
#  plot(model_1)

  ## First fixation
  set.seed(15)
  model_2 = brms::brm(PoL_Effect_angles_2~ 1 + (1|subjectnum), 
                      prior = myprior,
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 1000, 
                      iter = 4000,
                      silent =0,
                      thin = 2,
                      file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_2","_",format(Sys.Date(), "%d%m%Y"),".rds")
  )
  describe_posterior(model_2)
  
  set.seed(15)
  model_3 = brms::brm(PoL_Effect_angles_3~ 1 + (1|subjectnum), 
                      prior = myprior,
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 1000, 
                      iter = 4000,
                      silent =0,
                      thin = 2,
  file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                      #sample_prior = "only"
  )
  describe_posterior(model_3)
 # plot(model_3)
  
  set.seed(15)
  model_4 = brms::brm(PoL_Effect_angles_4~ 1 + (1|subjectnum), 
                      prior = myprior,
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 1000, 
                      iter = 4000,
                      silent =0,
                      thin = 2,
                     file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_4","_",format(Sys.Date(), "%d%m%Y"),".rds")
                      #sample_prior = "only"
  )
  describe_posterior(model_4)
#  plot(model_4)
  
  set.seed(15)
  model_5 = brms::brm(PoL_Effect_angles_5~ 1 + (1|subjectnum), 
                      prior = myprior,
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 1000, 
                      iter = 4000,
                      silent =0,
                      thin = 2,
                      file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_5","_",format(Sys.Date(), "%d%m%Y"),".rds")
                      #sample_prior = "only"
  )
  describe_posterior(model_5)
 # plot(model_5)
  
  set.seed(15)
  model_6 = brms::brm(PoL_Effect_angles_6~ 1 + (1|subjectnum), 
                      prior = myprior,
                      data = MSC,
                      chain = 4,
                      cores = 4,
                      backend = "cmdstan",
                      warmup = 1000, 
                      iter = 4000,
                      silent =0,
                      thin = 2,
                      file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_6","_",format(Sys.Date(), "%d%m%Y"),".rds")
                      #sample_prior = "only"
  )
  describe_posterior(model_6)
  
  describe_posterior(model_2)
  describe_posterior(model_3)
  describe_posterior(model_4)
  describe_posterior(model_5)
  describe_posterior(model_6)

  ## Since it wasn't significant for either 4th and 5th saccades, we didn't examine further
  # model_7 = brms::brm(PoL_Effect_angles_7~ 1 + (1|subjectnum), 
  #                     prior = myprior,
  #                     data = MSC,
  #                     chain = 4,
  #                     cores = 4,
  #                     backend = "cmdstan",
  #                     warmup = 1000, 
  #                     iter = 4000,
  #                     silent =0,
  #                     thin = 2,
  #                     file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_7","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                     #sample_prior = "only"
  # )
  # describe_posterior(model_7)
  # plot(model_7)
  # 
  # model_8 = brms::brm(PoL_Effect_angles_8~ 1 + (1|subjectnum), 
  #                     prior = myprior,
  #                     data = MSC,
  #                     chain = 4,
  #                     cores = 4,
  #                     backend = "cmdstan",
  #                     warmup = 1000, 
  #                     iter = 4000,
  #                     silent =0,
  #                     thin = 2,
  #                     file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_8","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                     #sample_prior = "only"
  # )
  # describe_posterior(model_8)
  # plot(model_8)
  # 
  # 
  # model_9 = brms::brm(PoL_Effect_angles_9~ 1 + (1|subjectnum), 
  #                     prior = myprior,
  #                     data = MSC,
  #                     chain = 4,
  #                     cores = 4,
  #                     backend = "cmdstan",
  #                     warmup = 1000, 
  #                     iter = 4000,
  #                     silent =0,
  #                     thin = 2,
  #                     file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_9","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                     #sample_prior = "only"
  # )
  # describe_posterior(model_9)
  # plot(model_9)
  # 
  # 
  # model_10 = brms::brm(PoL_Effect_angles_10~ 1 + (1|subjectnum), 
  #                      prior = myprior,
  #                      data = MSC,
  #                      chain = 4,
  #                      cores = 4,
  #                      backend = "cmdstan",
  #                      warmup = 1000, 
  #                      iter = 4000,
  #                      silent =0,
  #                      thin = 2,
  #                      file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_10","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                      #sample_prior = "only"
  # )
  # describe_posterior(model_10)
  # plot(model_10)
  
  #### Distance!!
  myprior_dis = prior(normal(0,100), class = "Intercept")
  ##distance 
  model = PoL_Effect_distance_1~ 1 + (1|subjectnum)
  library(brms)
  # model_1_dis = brms::brm(model, 
  #                         prior = myprior_dis,
  #                         data = MSC,
  #                         chain = 4,
  #                         cores = 4,
  #                         backend = "cmdstan",
  #                         warmup = 1000, 
  #                         iter = 4000,
  #                         silent =0,
  #                         thin = 2,
  #                         file =paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_1","_",format(Sys.Date(), "%d%m%Y"),".rds") 
  #                         #sample_prior = "only"
  # )
  # describe_posterior(model_1_dis)
  # plot(model_1_dis)
  

  set.seed(15)
  model_2_dis = brms::brm(PoL_Effect_distance_2~ 1 + (1|subjectnum), 
                          prior = myprior_dis,
                          data = MSC,
                          chain = 4,
                          cores = 4,
                          backend = "cmdstan",
                          warmup = 1000, 
                          iter = 4000,
                          silent =0,
                          thin = 2,
                          file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_2","_",format(Sys.Date(), "%d%m%Y"),".rds")
                          #sample_prior = "only"
  )
  describe_posterior(model_2_dis)
  plot(model_2_dis)
  
  set.seed(15)
  model_3_dis = brms::brm(PoL_Effect_distance_3~ 1 + (1|subjectnum), 
                          prior = myprior_dis,
                          data = MSC,
                          chain = 4,
                          cores = 4,
                          backend = "cmdstan",
                          warmup = 1000, 
                          iter = 4000,
                          silent =0,
                          thin = 2,
                          file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                          #sample_prior = "only"
  )
  describe_posterior(model_3_dis)
  plot(model_3_dis)
  
  set.seed(15)
  model_4_dis = brms::brm(PoL_Effect_distance_4~ 1 + (1|subjectnum), 
                          prior = myprior_dis,
                          data = MSC,
                          chain = 4,
                          cores = 4,
                          backend = "cmdstan",
                          warmup = 1000, 
                          iter = 4000,
                          silent =0,
                          thin = 2,
                          file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_4","_",format(Sys.Date(), "%d%m%Y"),".rds")
                          
                          #sample_prior = "only"
  )
  describe_posterior(model_4_dis)
  plot(model_4_dis)
  
  set.seed(15)
  model_5_dis = brms::brm(PoL_Effect_distance_5~ 1 + (1|subjectnum), 
                          prior = myprior_dis,
                          data = MSC,
                          chain = 4,
                          cores = 4,
                          backend = "cmdstan",
                          warmup = 1000, 
                          iter = 4000,
                          silent =0,
                          thin = 2,
                          file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_5","_",format(Sys.Date(), "%d%m%Y"),".rds")
                          
                          #sample_prior = "only"
  )
  describe_posterior(model_5_dis)
  plot(model_5_dis)

  
  set.seed(15)
  model_6_dis = brms::brm(PoL_Effect_distance_6~ 1 + (1|subjectnum), 
                          prior = myprior_dis,
                          data = MSC,
                          chain = 4,
                          cores = 4,
                          backend = "cmdstan",
                          warmup = 1000, 
                          iter = 4000,
                          silent =0,
                          thin = 2,
                          file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_6","_",format(Sys.Date(), "%d%m%Y"),".rds")
                          #sample_prior = "only"
  )
  describe_posterior(model_6_dis)
  plot(model_6_dis)
  
  describe_posterior(model_2_dis)
  describe_posterior(model_3_dis)
  describe_posterior(model_4_dis)
  describe_posterior(model_5_dis)
  describe_posterior(model_6_dis)
  # 
  # model_7_dis = brms::brm(PoL_Effect_distance_7~ 1 + (1|subjectnum), 
  #                         prior = myprior_dis,
  #                         data = MSC,
  #                         chain = 4,
  #                         cores = 4,
  #                         backend = "cmdstan",
  #                         warmup = 1000, 
  #                         iter = 4000,
  #                         silent =0,
  #                         thin = 2,
  #                         file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_7","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                         #sample_prior = "only"
  # )
  # describe_posterior(model_7_dis)
  # plot(model_7_dis)
  # 
  # model_8_dis = brms::brm(PoL_Effect_distance_8~ 1 + (1|subjectnum), 
  #                         prior = myprior_dis,
  #                         data = MSC,
  #                         chain = 4,
  #                         cores = 4,
  #                         backend = "cmdstan",
  #                         warmup = 1000, 
  #                         iter = 4000,
  #                         silent =0,
  #                         thin = 2,
  #                         file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_8","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                         
  #                         #sample_prior = "only"
  # )
  # describe_posterior(model_8_dis)
  # plot(model_8_dis)
  # 
  # 
  # model_9_dis = brms::brm(PoL_Effect_distance_9~ 1 + (1|subjectnum), 
  #                         prior = myprior_dis,
  #                         data = MSC,
  #                         chain = 4,
  #                         cores = 4,
  #                         backend = "cmdstan",
  #                         warmup = 1000, 
  #                         iter = 4000,
  #                         silent =0,
  #                         thin = 2,
  #                         file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_9","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                         
  #                         #sample_prior = "only"
  # )
  # describe_posterior(model_9_dis)
  # plot(model_9_dis)
  # 
  # 
  # model_10_dis = brms::brm(PoL_Effect_distance_10~ 1 + (1|subjectnum), 
  #                          prior = myprior_dis,
  #                          data = MSC,
  #                          chain = 4,
  #                          cores = 4,
  #                          backend = "cmdstan",
  #                          warmup = 1000, 
  #                          iter = 4000,
  #                          silent =0,
  #                          thin = 2,
  #                          file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_distance_10","_",format(Sys.Date(), "%d%m%Y"),".rds")
  #                          #sample_prior = "only"
  # )
  # describe_posterior(model_10_dis)
  # plot(model_10_dis)
  # 
  
  
#}



# Exploratory: ----


## 3.image baselines ----
MSC = mutate(MSC, PoL_Effect_angles_1_img = PREV_OBJ_FIX_DOT_PRODUCT_1 -AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_1,
             PoL_Effect_angles_2_img = PREV_OBJ_FIX_DOT_PRODUCT_2 -AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_2,
             PoL_Effect_angles_3_img = PREV_OBJ_FIX_DOT_PRODUCT_3 -AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_3,
             PoL_Effect_angles_4_img = PREV_OBJ_FIX_DOT_PRODUCT_4 -AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_4,
             PoL_Effect_angles_5_img = PREV_OBJ_FIX_DOT_PRODUCT_5 -AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_5,
             PoL_Effect_angles_6_img = PREV_OBJ_FIX_DOT_PRODUCT_6 -AVG_CURR_IMG_DIST_PREV_OBJ_DOT_PRODUCT_FIX_6,
             PoL_Effect_distance_1_img = PREV_ORTHOGONAL_DISTANCE_1 -AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_1,
             PoL_Effect_distance_2_img = PREV_ORTHOGONAL_DISTANCE_2 -AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_2,
             PoL_Effect_distance_3_img = PREV_ORTHOGONAL_DISTANCE_3 -AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_3,             
             PoL_Effect_distance_4_img = PREV_ORTHOGONAL_DISTANCE_4 -AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_4,
             PoL_Effect_distance_5_img = PREV_ORTHOGONAL_DISTANCE_5 -AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_5,
             PoL_Effect_distance_6_img = PREV_ORTHOGONAL_DISTANCE_5 -AVG_CURR_IMG_DIST_PREV_OBJ_DISTANCE_FIX_6,             
             
)


library(brms)
myprior = prior(normal(0,0.3), class = "Intercept")
#get_prior(model, data = MSC)

myprior_dis = prior(normal(0,100), class = "Intercept")
model = PoL_Effect_angles_1_img~ 1 + (1|subjectnum)
# model_1 = brms::brm(model, 
#                     prior = myprior,
#                     data = MSC,
#                     chain = 4,
#                     cores = 4,
#                     backend = "cmdstan",
#                     warmup = 1000, 
#                     iter = 4000,
#                     silent =0,
#                     thin = 2,
#                     #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_1","_",format(Sys.Date(), "%d%m%Y"),".rds")
#                     #sample_prior = "only"
# )
# prior_summary(model_1)
# summary(model_1)
library(bayestestR)
# describe_posterior(model_1)
# plot(model_1)
## same for distance
set.seed(15)
model_2 = brms::brm(PoL_Effect_angles_2_img~ 1 + (1|subjectnum), 
                    prior = myprior,
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 1000, 
                    iter = 4000,
                    silent =0,
                    thin = 2
                    #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_2","_",format(Sys.Date(), "%d%m%Y"),".rds")
                    #sample_prior = "only"
)
describe_posterior(model_2)
#plot(model_2)

set.seed(15)
model_3 = brms::brm(PoL_Effect_angles_3_img~ 1 + (1|subjectnum), 
                    prior = myprior,
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 1000, 
                    iter = 4000,
                    silent =0,
                    thin = 2,
                    #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                    #sample_prior = "only"
)
describe_posterior(model_3)
#plot(model_3)
set.seed(15)

model_4 = brms::brm(PoL_Effect_angles_4_img~ 1 + (1|subjectnum), 
                    prior = myprior,
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 1000, 
                    iter = 4000,
                    silent =0,
                    thin = 2,
                    #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                    #sample_prior = "only"
)
describe_posterior(model_4)
#plot(model_4)
set.seed(15)

model_5 = brms::brm(PoL_Effect_angles_5_img~ 1 + (1|subjectnum), 
                    prior = myprior,
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 1000, 
                    iter = 4000,
                    silent =0,
                    thin = 2,
                    #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                    #sample_prior = "only"
)
describe_posterior(model_5)
#plot(model_5)
set.seed(15)

model_6 = brms::brm(PoL_Effect_angles_6_img~ 1 + (1|subjectnum), 
                    prior = myprior,
                    data = MSC,
                    chain = 4,
                    cores = 4,
                    backend = "cmdstan",
                    warmup = 1000, 
                    iter = 4000,
                    silent =0,
                    thin = 2,
                    #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                    #sample_prior = "only"
)
describe_posterior(model_6)
#plot(model_6)

describe_posterior(model_2)
describe_posterior(model_3)
describe_posterior(model_4)        
describe_posterior(model_5)
describe_posterior(model_6)

myprior_dis = prior(normal(0,100), class = "Intercept")
# model = PoL_Effect_distance_1_img~ 1 + (1|subjectnum)
# model_1_dis = brms::brm(model, 
#                     prior = myprior_dis,
#                     data = MSC,
#                     chain = 4,
#                     cores = 4,
#                     backend = "cmdstan",
#                     warmup = 1000, 
#                     iter = 4000,
#                     silent =0,
#                     thin = 2,
# #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_1","_",format(Sys.Date(), "%d%m%Y"),".rds")
#                     #sample_prior = "only"
# )
# library(bayestestR)
# describe_posterior(model_1_dis)

set.seed(15)
model_2_dis = brms::brm(PoL_Effect_distance_2_img~ 1 + (1|subjectnum), 
                        prior = myprior_dis,
                        data = MSC,
                        chain = 4,
                        cores = 4,
                        backend = "cmdstan",
                        warmup = 1000, 
                        iter = 4000,
                        silent =0,
                        thin = 2
                        #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_2","_",format(Sys.Date(), "%d%m%Y"),".rds")
                        #sample_prior = "only"
)
describe_posterior(model_2_dis)
#plot(model_2)

set.seed(15)
model_3_dis = brms::brm(PoL_Effect_distance_3_img~ 1 + (1|subjectnum), 
                        prior = myprior_dis,
                        data = MSC,
                        chain = 4,
                        cores = 4,
                        backend = "cmdstan",
                        warmup = 1000, 
                        iter = 4000,
                        silent =0,
                        thin = 2,
                        #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                        #sample_prior = "only"
)
describe_posterior(model_3_dis)
#plot(model_3)
set.seed(15)

model_4_dis = brms::brm(PoL_Effect_distance_4_img~ 1 + (1|subjectnum), 
                        prior = myprior_dis,
                        data = MSC,
                        chain = 4,
                        cores = 4,
                        backend = "cmdstan",
                        warmup = 1000, 
                        iter = 4000,
                        silent =0,
                        thin = 2,
                        #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                        #sample_prior = "only"
)
describe_posterior(model_4_dis)
set.seed(15)

model_5_dis = brms::brm(PoL_Effect_distance_5_img~ 1 + (1|subjectnum), 
                        prior = myprior_dis,
                        data = MSC,
                        chain = 4,
                        cores = 4,
                        backend = "cmdstan",
                        warmup = 1000, 
                        iter = 4000,
                        silent =0,
                        thin = 2,
                        #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                        #sample_prior = "only"
)
describe_posterior(model_5_dis)
set.seed(15)

model_6_dis = brms::brm(PoL_Effect_distance_6_img~ 1 + (1|subjectnum), 
                        prior = myprior_dis,
                        data = MSC,
                        chain = 4,
                        cores = 4,
                        backend = "cmdstan",
                        warmup = 1000, 
                        iter = 4000,
                        silent =0,
                        thin = 2,
                        #                    file = paste0(write_path,"n",settings$POL_EFFECT_ORDER,"_",settings$object,"_PoL_Effect_angles_3","_",format(Sys.Date(), "%d%m%Y"),".rds")
                        #sample_prior = "only"
)
describe_posterior(model_6_dis)

describe_posterior(model_2_dis);describe_posterior(model_3_dis);describe_posterior(model_4_dis);
describe_posterior(model_5_dis);describe_posterior(model_6_dis)


## 4.	Priming of last fixation vs. priming previous target----

## re-extracting all features in relation to the last fixation and not the previous target
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

## This anlaysis is done only on absent -> absent trials!
GROUPED_MSC <- GROUPED_MSC %>% filter(PREV_condition == "absent")
GROUPED_MSC <- GROUPED_MSC %>% filter(condition == "absent")

if (settings$prev_tar_hit) {
  GROUPED_MSC <- GROUPED_MSC %>% filter(PREV_accuracy == 1)
}

FIXATIONS_BASELINE_DF <- GROUPED_MSC

if (settings$test_subject){
  GROUPED_MSC <- GROUPED_MSC %>% filter(subjectnum != settings$test_subject)
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
      current_x = current_x[fix_num]
      current_y = current_y[fix_num]
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
  col_name <- paste0("DIST_PREV_fix_OG_DISTANCE_FIX_", i)
  avg_col_name <- paste0("AVG_DIST_PREV_fix_OG_DISTANCE_FIX_", i)
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
      current_x = current_x[fix_num]
      current_y = current_y[fix_num]
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
  col_name <- paste0("DIST_PREV_OBJ_fix_PRODUCT_FIX_", i)
  avg_col_name <- paste0("AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_", i)
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

df_angles_priming = GROUPED_MSC %>%
  group_by(subjectnum) %>%
  summarize(
    across(
      c(contains("PREV_OBJ_fix_PRODUCT"), starts_with("PREV_OBJ_FIX_DOT_PRODUCT"),-"PREV_OBJ_FIX_DOT_PRODUCT"),
      mean,
      na.rm = TRUE
    )
  )%>%select(subjectnum,starts_with("PREV"), starts_with("AVG"))

res_priming = df_angles_priming%>%pivot_longer(PREV_OBJ_FIX_DOT_PRODUCT_1:AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_10,names_to = "variable", values_to = "value")
res_priming$fixation_num = as.numeric(str_extract(res_priming$variable, "\\d+$"))
res_priming$variable = sub("_[0-9]+$", "", res_priming$variable)

final_angles_priming = res_priming%>%group_by(fixation_num,variable)%>%summarize(value = mean(value,na.rm = TRUE))%>%pivot_wider(names_from = variable, values_from = value)%>%
  rename("subject_baseline"= "AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX",
         "measure"  = "PREV_OBJ_FIX_DOT_PRODUCT") 


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

df_distance_priming = GROUPED_MSC %>%
  group_by(subjectnum) %>%
  summarize(
    across(
      c(contains("PREV_DISTANCE"), contains("OG_DISTANCE_FIX_"),-"PREV_DISTANCE_last_fixation"),
      mean,
      na.rm = TRUE
    )
  )%>%select(subjectnum,starts_with("PREV"), starts_with("AVG"))

res_priming_dis = df_distance_priming%>%pivot_longer(PREV_DISTANCE_last_fixation_1:AVG_DIST_PREV_fix_OG_DISTANCE_FIX_10,names_to = "variable", values_to = "value")
res_priming_dis$fixation_num = as.numeric(str_extract(res_priming_dis$variable, "\\d+$"))
res_priming_dis$variable = sub("_[0-9]+$", "", res_priming_dis$variable)

final_distance_priming = res_priming_dis%>%group_by(fixation_num,variable)%>%summarize(value = mean(value,na.rm = TRUE))%>%pivot_wider(names_from = variable, values_from = value)%>%
  rename("subject_baseline"= "AVG_DIST_PREV_fix_OG_DISTANCE_FIX",
         "measure"  = "PREV_DISTANCE_last_fixation") 


#remove example: 
GROUPED_MSC = mutate(GROUPED_MSC, 
            last_fixation_Effect_distance_1=   PREV_DISTANCE_last_fixation_1 - AVG_DIST_PREV_fix_OG_DISTANCE_FIX_1, 
            last_fixation_Effect_distance_2=   PREV_DISTANCE_last_fixation_2 - AVG_DIST_PREV_fix_OG_DISTANCE_FIX_2,
            last_fixation_Effect_distance_3=   PREV_DISTANCE_last_fixation_3 - AVG_DIST_PREV_fix_OG_DISTANCE_FIX_3,
            last_fixation_Effect_distance_4=   PREV_DISTANCE_last_fixation_4 - AVG_DIST_PREV_fix_OG_DISTANCE_FIX_4,
            last_fixation_Effect_distance_5=   PREV_DISTANCE_last_fixation_5 - AVG_DIST_PREV_fix_OG_DISTANCE_FIX_5,
            last_fixation_Effect_angles_1=  PREV_OBJ_FIX_DOT_PRODUCT_1 - AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_1,
            last_fixation_Effect_angles_2=  PREV_OBJ_FIX_DOT_PRODUCT_2 - AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_2,
            last_fixation_Effect_angles_3=  PREV_OBJ_FIX_DOT_PRODUCT_3 - AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_3,
            last_fixation_Effect_angles_4=  PREV_OBJ_FIX_DOT_PRODUCT_4 - AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_4,
            last_fixation_Effect_angles_5=  PREV_OBJ_FIX_DOT_PRODUCT_5 - AVG_DIST_PREV_OBJ_fix_PRODUCT_FIX_5
)


## test fixation priming
library(brms)
# model_fix_primig1 = brms::brm(last_fixation_Effect_distance_1~ 1 + (1|subjectnum), 
#                               # prior = myprior,
#                               data = GROUPED_MSC,
#                               chain = 4,
#                               cores = 4,
#                               backend = "cmdstan",
#                               warmup = 1000, 
#                               iter = 4000,
#                               silent =0,
#                               thin = 2
#                               #sample_prior = "only"
# )
# #prior_summary(model_3)
# summary(model_fix_primig1)
library(bayestestR)
# describe_posterior(model_fix_primig1)
#nothing here

set.seed(15)
model_fix_primig2 = brms::brm(last_fixation_Effect_distance_2~ 1 + (1|subjectnum) , 
                              # prior = myprior,
                              data = GROUPED_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
summary(model_fix_primig2)
describe_posterior(model_fix_primig2)
#plot(model_fix_primig2)


set.seed(15)
model_fix_primig3 = brms::brm(last_fixation_Effect_distance_3~ 1 + (1|subjectnum) , 
                              # prior = myprior,
                              data = GROUPED_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
describe_posterior(model_fix_primig3)

set.seed(15)
model_fix_primig4 = brms::brm(last_fixation_Effect_distance_4~ 1 + (1|subjectnum) , 
                              # prior = myprior,
                              data = GROUPED_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
describe_posterior(model_fix_primig4)

set.seed(15)
model_fix_primig5 = brms::brm(last_fixation_Effect_distance_5~ 1 + (1|subjectnum) , 
                              # prior = myprior,
                              data = GROUPED_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
#summary(unique_PoL_distance_5)
describe_posterior(model_fix_primig5)
#plot(unique_PoL_distance_5)

####

library(brms)
# set.seed(15)
# model_fix_primig1_ang = brms::brm(last_fixation_Effect_angles_1~  1 + (1|subjectnum), 
#                                   # prior = myprior,
#                                   data = GROUPED_MSC,
#                                   chain = 4,
#                                   cores = 4,
#                                   backend = "cmdstan",
#                                   warmup = 1000, 
#                                   iter = 4000,
#                                   silent =0,
#                                   thin = 2
#                                   #sample_prior = "only"
# )
#prior_summary(model_3)
#summary(model_fix_primig1_ang)
#describe_posterior(model_fix_primig1_ang)
#plot(model_fix_primig1_ang)

set.seed(15)
model_fix_primig2_ang = brms::brm(last_fixation_Effect_angles_2~  1 + (1|subjectnum) , 
                                  # prior = myprior,
                                  data = GROUPED_MSC,
                                  chain = 4,
                                  cores = 4,
                                  backend = "cmdstan",
                                  warmup = 1000, 
                                  iter = 4000,
                                  silent =0,
                                  thin = 2
                                  #sample_prior = "only"
)
#prior_summary(model_3)
summary(model_fix_primig2_ang)
describe_posterior(model_fix_primig2_ang)
#plot(model_fix_primig2)

## the effect is significant and negative
set.seed(15)
model_fix_primig3_ang = brms::brm(last_fixation_Effect_angles_3~  1 + (1|subjectnum) , 
                                  # prior = myprior,
                                  data = GROUPED_MSC,
                                  chain = 4,
                                  cores = 4,
                                  backend = "cmdstan",
                                  warmup = 1000, 
                                  iter = 4000,
                                  silent =0,
                                  thin = 2
                                  #sample_prior = "only"
)
#prior_summary(model_3)
summary(model_fix_primig3_ang)
describe_posterior(model_fix_primig3_ang)

#also here it is significant 
set.seed(15)
model_fix_primig4_ang = brms::brm(last_fixation_Effect_angles_4~  1 + (1|subjectnum) , 
                                  # prior = myprior,
                                  data = GROUPED_MSC,
                                  chain = 4,
                                  cores = 4,
                                  backend = "cmdstan",
                                  warmup = 1000, 
                                  iter = 4000,
                                  silent =0,
                                  thin = 2
                                  #sample_prior = "only"
)

### location priming above and beyond fixaiton priming

sum_MSC=MSC%>%group_by(subjectnum)%>%
  summarize (PoL_Effect_angles_1 = mean(PoL_Effect_angles_1),
             PoL_Effect_angles_2 = mean(PoL_Effect_angles_2),
             PoL_Effect_angles_3 = mean(PoL_Effect_angles_3),
             PoL_Effect_angles_4 = mean(PoL_Effect_angles_4),
             PoL_Effect_angles_5 = mean(PoL_Effect_angles_5),
            PoL_Effect_distance_1 = mean(PoL_Effect_distance_1),
            PoL_Effect_distance_2 = mean(PoL_Effect_distance_2),
            PoL_Effect_distance_3 = mean(PoL_Effect_distance_3),
            PoL_Effect_distance_4 = mean(PoL_Effect_distance_4),
            PoL_Effect_distance_5 = mean(PoL_Effect_distance_5)
             )

sum_groupd_MSC=GROUPED_MSC%>%group_by(subjectnum)%>%
  summarize (last_fixation_Effect_distance_1 = mean(last_fixation_Effect_distance_1),
             last_fixation_Effect_distance_2 = mean(last_fixation_Effect_distance_2),
             last_fixation_Effect_distance_3 = mean(last_fixation_Effect_distance_3),
             last_fixation_Effect_distance_4 = mean(last_fixation_Effect_distance_4),
             last_fixation_Effect_distance_5 = mean(last_fixation_Effect_distance_5),
             last_fixation_Effect_angles_1 = mean(last_fixation_Effect_angles_1),
             last_fixation_Effect_angles_2 = mean(last_fixation_Effect_angles_2),
             last_fixation_Effect_angles_3 = mean(last_fixation_Effect_angles_3),
             last_fixation_Effect_angles_4 = mean(last_fixation_Effect_angles_4),
             last_fixation_Effect_angles_5 = mean(last_fixation_Effect_angles_5)
             )

## 
joint_MSC = left_join(sum_MSC,sum_groupd_MSC,
                      by = c("subjectnum"))
joint_MSC = mutate(joint_MSC, 
                   unique_PoL_distance_1 = PoL_Effect_distance_1- last_fixation_Effect_distance_1,
                   unique_PoL_distance_2 = PoL_Effect_distance_2- last_fixation_Effect_distance_2,
                   unique_PoL_distance_3 = PoL_Effect_distance_3- last_fixation_Effect_distance_3,
                   unique_PoL_distance_4 = PoL_Effect_distance_4- last_fixation_Effect_distance_4,
                   unique_PoL_distance_5 = PoL_Effect_distance_5- last_fixation_Effect_distance_5,
                   unique_PoL_angles_1 = PoL_Effect_angles_1 - last_fixation_Effect_angles_1,
                   unique_PoL_angles_2 = PoL_Effect_angles_2 - last_fixation_Effect_angles_2,
                   unique_PoL_angles_3 = PoL_Effect_angles_3 - last_fixation_Effect_angles_3,
                   unique_PoL_angles_4 = PoL_Effect_angles_4 - last_fixation_Effect_angles_4,
                   unique_PoL_angles_5 = PoL_Effect_angles_5 - last_fixation_Effect_angles_5)

library(brms)

set.seed(15)
model_fix_primig2 = brms::brm(unique_PoL_distance_2~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
summary(model_fix_primig2)
describe_posterior(model_fix_primig2)

set.seed(15)
model_fix_primig3 = brms::brm(unique_PoL_distance_3~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
#summary(unique_PoL_distance_3)
describe_posterior(model_fix_primig3)
#plot(unique_PoL_distance_3)

set.seed(15)
model_fix_primig4 = brms::brm(unique_PoL_distance_4~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
describe_posterior(model_fix_primig4)

set.seed(15)
model_fix_primig5 = brms::brm(unique_PoL_distance_5~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
#summary(unique_PoL_distance_5)
describe_posterior(model_fix_primig5)
#plot(unique_PoL_distance_5)

####

library(brms)
set.seed(15)
model_fix_primig2_ang = brms::brm(unique_PoL_angles_2~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
summary(model_fix_primig2_ang)
describe_posterior(model_fix_primig2_ang)

set.seed(15)
model_fix_primig3_ang = brms::brm(unique_PoL_angles_3~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
summary(model_fix_primig3_ang)
describe_posterior(model_fix_primig3_ang)

set.seed(15)
model_fix_primig4_ang = brms::brm(unique_PoL_angles_4~ 1 , 
                              # prior = myprior,
                              data = joint_MSC,
                              chain = 4,
                              cores = 4,
                              backend = "cmdstan",
                              warmup = 1000, 
                              iter = 4000,
                              silent =0,
                              thin = 2
                              #sample_prior = "only"
)
#prior_summary(model_3)
summary(model_fix_primig4_ang)
describe_posterior(model_fix_primig4_ang)
#plot(unique_PoL_distance_4)



## 5.	Dwell time as a function of PoL ----
## we did not include this exploratory analysis eventually 
# MSC <- MSC%>%mutate(
#   CURRENT_FIX_DURATION_1 = sapply(CURRENT_FIX_DURATION, function(x) x[1]),
#   CURRENT_FIX_DURATION_2 = sapply(CURRENT_FIX_DURATION, function(x) x[2]),
#   CURRENT_FIX_DURATION_3 = sapply(CURRENT_FIX_DURATION, function(x) x[3]),
#   CURRENT_FIX_DURATION_4 = sapply(CURRENT_FIX_DURATION, function(x) x[4]),
#   CURRENT_FIX_DURATION_5 = sapply(CURRENT_FIX_DURATION, function(x) x[5]),
#   CURRENT_FIX_DURATION_6 = sapply(CURRENT_FIX_DURATION, function(x) x[6]),
#   CURRENT_FIX_DURATION_7 = sapply(CURRENT_FIX_DURATION, function(x) x[7]),
#   CURRENT_FIX_DURATION_8 = sapply(CURRENT_FIX_DURATION, function(x) x[8]),
#   CURRENT_FIX_DURATION_9 = sapply(CURRENT_FIX_DURATION, function(x) x[9]),
#   CURRENT_FIX_DURATION_10 = sapply(CURRENT_FIX_DURATION, function(x) x[10])
# )
# 
# model = CURRENT_FIX_DURATION_1~ PoL_Effect_distance_1 + (PoL_Effect_distance_1|subjectnum) ## if works, if not remove random slopes
# get_prior(model, data = MSC)
# 
# library(brms)
# model_dur_1 = brms::brm(model, 
#                    #     prior = myprior_dis,
#                         data = MSC,
#                         chain = 4,
#                         cores = 4,
#                         backend = "cmdstan",
#                         warmup = 1000, 
#                         iter = 4000,
#                         silent =0,
#                         thin = 2
#                       #sample_prior = "only"
# )
# prior_summary(model_dur_1)
# describe_posterior(model_dur_1)
# plot(model_dur_1)
# 
# 
# 
# model_dur_2 = brms::brm(CURRENT_FIX_DURATION_2~ PoL_Effect_distance_2 + (PoL_Effect_distance_2|subjectnum), 
#                         #     prior = myprior_dis,
#                         data = MSC,
#                         chain = 4,
#                         cores = 4,
#                         backend = "cmdstan",
#                         warmup = 1000, 
#                         iter = 4000,
#                         silent =0,
#                         thin = 2
#                         #sample_prior = "only"
# )
# prior_summary(model_dur_2)
# describe_posterior(model_dur_2)
# plot(model_dur_2)
# library(emmeans)
# 
# model_dur_3 = brms::brm(CURRENT_FIX_DURATION_3~ PoL_Effect_distance_3 + (PoL_Effect_distance_3|subjectnum), 
#                         #     prior = myprior_dis,
#                         data = MSC,
#                         chain = 4,
#                         cores = 4,
#                         backend = "cmdstan",
#                         warmup = 1000, 
#                         iter = 4000,
#                         silent =0,
#                         thin = 2
#                         #sample_prior = "only"
# )
# describe_posterior(model_dur_3)
# plot(model_dur_3)
# 
# # not significan here
# 
# model_dur_4 = brms::brm(CURRENT_FIX_DURATION_4~ PoL_Effect_distance_4 + (PoL_Effect_distance_4|subjectnum), 
#                         #     prior = myprior_dis,
#                         data = MSC,
#                         chain = 4,
#                         cores = 4,
#                         backend = "cmdstan",
#                         warmup = 1000, 
#                         iter = 4000,
#                         silent =0,
#                         thin = 2
#                         #sample_prior = "only"
# )
# describe_posterior(model_dur_4)
# plot(model_dur_4)
# ## not significant 
# 
# model_dur_5 = brms::brm(CURRENT_FIX_DURATION_5~ PoL_Effect_distance_5 + (PoL_Effect_distance_5|subjectnum), 
#                         #     prior = myprior_dis,
#                         data = MSC,
#                         chain = 4,
#                         cores = 4,
#                         backend = "cmdstan",
#                         warmup = 1000, 
#                         iter = 4000,
#                         silent =0,
#                         thin = 2
#                         #sample_prior = "only"
# )
# describe_posterior(model_dur_5)
# plot(model_dur_5)

## not significant 



##6.	Time spent near the previous target  ----
# ## we did not include this exploratory analysis eventually 
# proxmity_fun = function(x_vec, y_vec, xmin, xmax, ymin, ymax,im_w,im_h) {
#   
#   # center of the object
#   cx <- (xmin + xmax) / 2
#   cy <- (ymin + ymax) / 2
#   
#   roi_w <- 0.5 * im_w
#   roi_h <- 0.5 * im_h
#   
#   xmin_p <- cx - roi_w / 2
#   xmax_p <- cx + roi_w / 2
#   ymin_p <- cy - roi_h / 2
#   ymax_p <- cy + roi_h / 2
#   
#   # classify each fixation
#   inside <- (x_vec >= xmin_p & x_vec <= xmax_p &
#                y_vec >= ymin_p & y_vec <= ymax_p)
#   
#   # convert TRUE/FALSE → 1/0
#   return(as.integer(inside))
# }
# 
# MSC <- MSC %>%
#   mutate(
#     proximity_vec = mapply(
#       {proxmity_fun},
#       CURRENT_FIX_X,
#       CURRENT_FIX_Y,
#       PREV_OBJECT_X1,
#       PREV_OBJECT_X2,
#       PREV_OBJECT_Y1,
#       PREV_OBJECT_Y2,im_w,im_h)
#   )
# 
# 
# MSC <- MSC %>%
#   mutate(
#     total_time_inside = mapply(function(dur, prox){
#       sum(dur[prox == 1], na.rm = TRUE)
#     }, CURRENT_FIX_DURATION, proximity_vec),
#     total_duration_time = unlist(lapply(CURRENT_FIX_DURATION,sum)),
#     prop_inside = total_time_inside/total_duration_time,
#     obj_ratio = get_obj_ratio(PREV_OBJECT_X1, PREV_OBJECT_Y1,
#                               PREV_OBJECT_X2, PREV_OBJECT_Y2,
#                               im_w, im_h)#,
#     # roi_ratio = mapply(function(xmin, xmax, ymin, ymax) {
#     #   # center of the object
#     #   cx <- (xmin + xmax) / 2
#     #   cy <- (ymin + ymax) / 2
#     #   
#     #   # ROI dimensions: 25% of screen
#     #   roi_w <- 0.5 * im_w
#     #   roi_h <- 0.5 * im_h
#     #   
#     #   # ROI box
#     #   xmin_p <- cx - roi_w / 2
#     #   xmax_p <- cx + roi_w / 2
#     #   ymin_p <- cy - roi_h / 2
#     #   ymax_p <- cy + roi_h / 2  
#     #   
#     #   return(get_obj_ratio(xmin_p, ymin_p,
#     #                 xmax_p, ymax_p,
#     #                 im_w, im_h))
#     #   
#     # }, PREV_OBJECT_X1, PREV_OBJECT_X2,
#     # PREV_OBJECT_Y1, PREV_OBJECT_Y2)
#   )
# table(MSC$obj_ratio)
# #table(MSC$roi_ratio)
# 
# MSC%>%group_by(subjectnum)%>%summarize(prop_inside =mean(prop_inside))%>%summarize(mean(prop_inside))
# #MSC%>%group_by(obj_ratio)%>%summarize(prop_inside =mean(prop_inside))%>%summarize(mean(prop_inside))
# 
# MSC%>%ggplot(aes(x = obj_ratio, y = prop_inside))+geom_point()
# 
# #mean(MSC$prop_inside)
# 
# ##baseline calculation 
# random_prev_object = MSC
# 
# set.seed(15)
# nrep = 100
# vec_df = as.data.frame(matrix(rep(NA,nrep*nrow(random_prev_object)), ncol =nrep , nrow = nrow(random_prev_object))) 
# 
# random_prev_object <- random_prev_object %>%
#   mutate(
#     obj_center_x = (PREV_OBJECT_X1 + PREV_OBJECT_X2) / 2,
#     obj_center_y = (PREV_OBJECT_Y1 + PREV_OBJECT_Y2) / 2,
#     CENTER_OG_DISTANCE_TRUE = sqrt((obj_center_x - CENTER[1])^2 + (obj_center_y - CENTER[2])^2)
#   )
# for (i in 1:nrep){
#   ## match for eccentricity 
# 
#   random_prev_object = mutate(random_prev_object,
#                               obj_width = PREV_OBJECT_X2-PREV_OBJECT_X1,
#                               obj_length = PREV_OBJECT_Y2-PREV_OBJECT_Y1,
#                               angle = runif(n(), 0, 2 * pi),
#                               center_x = CENTER[1] + CENTER_OG_DISTANCE_TRUE * cos(angle),
#                               center_y = CENTER[2] + CENTER_OG_DISTANCE_TRUE * sin(angle),
#                               x1_raw = center_x - obj_width  / 2,
#                               y1_raw = center_y - obj_length / 2,
#                               PREV_OBJECT_X1_R = x1_raw,
#                               PREV_OBJECT_X2_R = PREV_OBJECT_X1_R+obj_width,
#                               PREV_OBJECT_Y1_R = y1_raw,
#                               PREV_OBJECT_Y2_R = PREV_OBJECT_Y1_R+obj_length
#   )
#   random_prev_object <- random_prev_object %>%
#     mutate(
#       proximity_vec = mapply(
#         {proxmity_fun},
#         CURRENT_FIX_X,
#         CURRENT_FIX_Y,
#         PREV_OBJECT_X1_R,
#         PREV_OBJECT_X2_R,
#         PREV_OBJECT_Y1_R,
#         PREV_OBJECT_Y2_R,im_w,im_h)
#     )
#   
#   random_prev_object <- random_prev_object %>%
#     mutate(
#       total_time_inside = mapply(function(dur, prox){
#         sum(dur[prox == 1], na.rm = TRUE)
#       }, CURRENT_FIX_DURATION, proximity_vec),
#       total_duration_time = unlist(lapply(CURRENT_FIX_DURATION,sum)),
#       prop_inside = total_time_inside/total_duration_time
#     )
#   vec_df[,i] = random_prev_object$prop_inside
# }
# 
# MSC = MSC%>%ungroup()%>%mutate(
#              random_prop_inside = rowMeans(vec_df,na.rm = TRUE),
#              time_spent_pol = prop_inside - random_prop_inside
#              )
# 
# 
# model = time_spent_pol~ 1 + (1|subjectnum) ## if works, if not remove random slopes
# 
# 
# library(brms)
# model_time_spent = brms::brm(model, 
#                     data = MSC,
#                     chain = 4,
#                     cores = 4,
#                     backend = "cmdstan",
#                     warmup = 1000, 
#                     iter = 4000,
#                     silent =0,
#                     thin = 2
#                     #sample_prior = "only"
# )
# describe_posterior(model_time_spent)
# beep(5)


