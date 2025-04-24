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
RANDOM_POINTS_COUNT <- 1000

MSC_PATH = "C:/Users/averr/OneDrive/Desktop/TAU/Year 4/Semestar A/Psych Workshop/Scripts/1. Object Detection/MSC_DS_PROCESSED_DETECTINOS_25_12_2024.csv"
MSC <- read_csv(MSC_PATH)
MSC_ORIGINAL <- MSC


# Define script settings
settings <- list(
  filter_bad_detections = TRUE,
  filter_inconsecutive_trial = TRUE,
  cluster_fixations = TRUE,
  features_to_remove = c("expected", "searcharray"),
  fliter_prev_absent = TRUE,
  test_subject = TRUE,
  calc_mean_OG = FALSE
)

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


# ===========================Adding New Features=====================================

# Add hit - is fixation inside object
MSC <<- mutate(MSC, hit = ifelse(OBJECT_X1 <= CURRENT_FIX_X &
                                   OBJECT_Y1 <= CURRENT_FIX_Y &
                                   OBJECT_X2 >= CURRENT_FIX_X &
                                   OBJECT_Y2 >= CURRENT_FIX_Y
                                 , 1, 0))

MSC <<- mutate(MSC, p_hit = ((OBJECT_X2-OBJECT_X1) * (OBJECT_Y2-OBJECT_Y1))/(im_w * im_h)
)

MSC_ORIGINAL %>% filter(CURRENT_FIX_INDEX == 6) %>% ggplot() + geom_point(aes(x=CURRENT_FIX_X, y=CURRENT_FIX_Y)) + ylim(0,800) + xlim(0,1280) 
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
    p_prev_obj_hit = ((PREV_OBJ_X2-PREV_OBJ_X1) * (PREV_OBJ_Y2-PREV_OBJ_Y1))/(im_w * im_h),
    prev_condition = lag(condition, default = NA)
  )


# Fix 
MSC <- MSC %>% mutate(
lag_trial_index = lag(TRIAL_INDEX)
)

# Custom function
calc_hit <- function(FIX_X, FIX_Y, PREV_X1, PREV_Y1, PREV_X2, PREV_Y2) {
  return(ifelse(PREV_X1 <= FIX_X & PREV_Y1 <= FIX_Y & PREV_X2 >= FIX_X & PREV_Y2 >= FIX_Y, 1, 0))
}

# prev_hit
MSC$prev_hit <- mapply(function(x1, x2, x3, x4, x5, x6) calc_hit(x1, x2, x3, x4, x5, x6),
                       MSC$CURRENT_FIX_X, MSC$CURRENT_FIX_Y, MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
                       MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2, SIMPLIFY = FALSE)


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
      TRIAL_INDEX == lag_trial_index + 1
    )
}


### Prev absent DF for baseline

PREV_CURR_ABSENT_DF <- MSC %>%
  filter(prev_condition == "absent" & condition == "absent")

N <- 6  # Change this to your desired number

# Extract fixation coordinates dynamically
for (i in seq_len(N)) {
  PREV_CURR_ABSENT_DF[[paste0("FIX_", i, "_X")]] <- sapply(PREV_CURR_ABSENT_DF$CURRENT_FIX_X, function(x) x[i])
  PREV_CURR_ABSENT_DF[[paste0("FIX_", i, "_Y")]] <- sapply(PREV_CURR_ABSENT_DF$CURRENT_FIX_Y, function(x) x[i])
}

# PREV_CURR_ABSENT_DF$FIX_1_VECTOR <- mapply(
#   function(c1, c2) c(c1 - DISPLAY_W/2, c2 - DISPLAY_H/2),
#   PREV_CURR_ABSENT_DF$FIX_1_X, PREV_CURR_ABSENT_DF$FIX_1_Y, SIMPLIFY = FALSE)

N <- 5  # <-- change this to your actual N

# Loop through each i and compute the corresponding FIX_i_VECTOR
get_diff_vector <- function(coords, origin) {
  # coords: a vector of length 2, c(c1, c2)
  # origin: a vector of length 2, c(DISPLAY_W, DISPLAY_H)
  
  if (length(coords) != 2 || length(origin) != 2) {
    stop("Both 'coords' and 'origin' must be numeric vectors of length 2.")
  }
  
  return(coords - origin)
}


for (i in 1:N) {
  x_col <- paste0("FIX_", i, "_X")
  y_col <- paste0("FIX_", i, "_Y")
  vec_col <- paste0("FIX_", i, "_VECTOR")
  
  PREV_CURR_ABSENT_DF[[vec_col]] <- mapply(
    function(c1, c2) get_diff_vector(c(c1, c2), c(DISPLAY_W, DISPLAY_H)/2),
    PREV_CURR_ABSENT_DF[[x_col]], PREV_CURR_ABSENT_DF[[y_col]], SIMPLIFY = FALSE)
}



# Compute means dynamically (for OG dist)
PREV_ABSENT_FIX_AVG <- PREV_CURR_ABSENT_DF %>%
  group_by(subjectnum) %>%
  summarise(
    across(starts_with("FIX_"), mean, na.rm = TRUE)
  ) %>%
  ungroup()


PREV_CURR_ABSENT_FIX_VECTORS <- PREV_CURR_ABSENT_DF %>%
  group_by(subjectnum) %>%
  summarize(
    FIX_1_VECTOR = list(FIX_1_VECTOR),
    FIX_2_VECTOR = list(FIX_2_VECTOR),
    FIX_3_VECTOR = list(FIX_3_VECTOR),
    FIX_4_VECTOR = list(FIX_4_VECTOR),
    FIX_5_VECTOR = list(FIX_5_VECTOR),
  )



# CONTINUE FILTER ###############

if (settings$fliter_prev_absent){
  MSC <- MSC %>% filter(prev_condition == "present")
}

if (settings$test_subject){
  MSC <- MSC %>% filter(subjectnum == 1)
  # AVG_FIXATION_DF <- AVG_FIXATION_DF %>% filter(subjectnum == 1)
}


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


# ============================Vector Difference Dot Product==================================== 

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

get_curr_fix_prev_obj_dot_product <- function(FIX_X, FIX_Y, obj_x1, obj_y1, obj_x2, obj_y2){
  products <- numeric(length(FIX_X))
  obj_center <- c(obj_x1+obj_x2, obj_y1+obj_y2)/2
  obj_vector <- get_diff_vector(obj_center, CENTER)
  
  # Iterate over the elements of X and Y
  for (i in seq_along(FIX_X)) {
    fix_vector <- get_diff_vector(c(FIX_X[i], FIX_Y[i]), CENTER)
    products[i] <- calc_dot_product(obj_vector, fix_vector,
                                    normalize = TRUE, max_norm = NULL)
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

# The point is inside the image 
calculate_single_mean_distance  <- function(im_w, im_h, obj_x1, obj_y1, obj_x2, obj_y2) {
  # Calculate the Euclidean distance
  x_vector <- runif(RANDOM_POINTS_COUNT, min = (DISPLAY_W - im_w)/2, (DISPLAY_W - im_w)/2 + im_w)
  y_vector <- runif(RANDOM_POINTS_COUNT, min = (DISPLAY_H - im_h)/2, (DISPLAY_H - im_h)/2 + im_h)
  distances <- mapply(calc_orthogonal_distance,
                     x_vector, y_vector,
                     MoreArgs = list(obj_x1 = obj_x1, obj_y1 = obj_y1,
                                     obj_x2 = obj_x2, obj_y2 = obj_y2))
  
  return(mean(distances))
}


if (settings$calc_mean_OG) {
  MSC$MEAN_OG_DISTANCE <- mapply(calculate_single_mean_distance,
                               MSC$im_w, MSC$im_h,
                               MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
                               MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2)

}


# Calculate the OG distance between center and  the * previous * detection
MSC$CENTER_OG_DISTANCE <- mapply(calc_orthogonal_distance,
                               DISPLAY_W/2, DISPLAY_H/2,
                               MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
                               MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2)



# Calculate the OG distance between a specific fixation mean and  the * previous * detection

# For previous absent
MSC$AVG_FIX_1_PREV_OG_DISTANCE <- mapply(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_ABSENT_FIX_AVG$subjectnum == subj)
    if (length(idx) == 0) return(NA)

    fix_x <- PREV_ABSENT_FIX_AVG$FIX_1_X[idx]
    fix_y <- PREV_ABSENT_FIX_AVG$FIX_1_Y[idx]

    calc_orthogonal_distance(fix_x, fix_y, x1, y1, x2, y2)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)
# For previous absent
MSC$AVG_FIX_2_PREV_OG_DISTANCE <- mapply(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_ABSENT_FIX_AVG$subjectnum == subj)
    if (length(idx) == 0) return(NA)

    fix_x <- PREV_ABSENT_FIX_AVG$FIX_2_X[idx]
    fix_y <- PREV_ABSENT_FIX_AVG$FIX_2_Y[idx]

    calc_orthogonal_distance(fix_x, fix_y, x1, y1, x2, y2)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)

MSC$AVG_FIX_3_PREV_OG_DISTANCE <- mapply(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_ABSENT_FIX_AVG$subjectnum == subj)
    if (length(idx) == 0) return(NA)

    fix_x <- PREV_ABSENT_FIX_AVG$FIX_3_X[idx]
    fix_y <- PREV_ABSENT_FIX_AVG$FIX_3_Y[idx]

    calc_orthogonal_distance(fix_x, fix_y, x1, y1, x2, y2)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)

MSC$AVG_FIX_4_PREV_OG_DISTANCE <- mapply(
  function(subj, x1, y1, x2, y2) {
    idx <- which(PREV_ABSENT_FIX_AVG$subjectnum == subj)
    if (length(idx) == 0) return(NA)

    fix_x <- PREV_ABSENT_FIX_AVG$FIX_4_X[idx]
    fix_y <- PREV_ABSENT_FIX_AVG$FIX_4_Y[idx]

    calc_orthogonal_distance(fix_x, fix_y, x1, y1, x2, y2)
  },
  MSC$subjectnum,
  MSC$PREV_OBJ_X1, MSC$PREV_OBJ_Y1,
  MSC$PREV_OBJ_X2, MSC$PREV_OBJ_Y2
)







#####=== Baseline for dot product vector difference===#####
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
    obj_center = c(x1+x2, y1+y2)/2
    FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_1_VECTOR[[idx]]
    products <- numeric(length(FIX_VECTORS))
    for (i in seq_along(FIX_VECTORS)) {
      product = calc_dot_product(FIX_VECTORS[[i]], obj_center, 
                                 normalize = TRUE, max_norm = NULL)
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
    obj_center = c(x1+x2, y1+y2)/2
    FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_2_VECTOR[[idx]]
    products <- numeric(length(FIX_VECTORS))
    for (i in seq_along(FIX_VECTORS)) {
      product = calc_dot_product(FIX_VECTORS[[i]], obj_center,
                                 normalize = TRUE, max_norm = NULL)
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
    obj_center = c(x1+x2, y1+y2)/2
    FIX_VECTORS <- PREV_CURR_ABSENT_FIX_VECTORS$FIX_3_VECTOR[[idx]]
    products <- numeric(length(FIX_VECTORS))
    for (i in seq_along(FIX_VECTORS)) {
      product = calc_dot_product(FIX_VECTORS[[i]], obj_center, 
                                 normalize = TRUE, max_norm = NULL)
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

# ============================Analysis==================================== 

# Probability analysis

# MSC = MSC%>%mutate(
#   prev_hit_1 = sapply(prev_hit, function(x) x[1]),
#   prev_hit_2 = sapply(prev_hit, function(x) x[2]),
#   prev_hit_3 = sapply(prev_hit, function(x) x[3]),
#   prev_hit_4 = sapply(prev_hit, function(x) x[4]),
#   prev_hit_5 = sapply(prev_hit, function(x) x[5]),
#   prev_hit_6 = sapply(prev_hit, function(x) x[6]),
#   prev_hit_7 = sapply(prev_hit, function(x) x[7]),
#   prev_hit_8 = sapply(prev_hit, function(x) x[8])
# )
# 
# ##
# table(MSC$condition,MSC$prev_condition)
# 
# 
# prop.table(table(MSC$prev_hit_1))
# prop.table(table(MSC$prev_hit_2))
# prop.table(table(MSC$prev_hit_3))
# prop.table(table(MSC$prev_hit_4))
# prop.table(table(MSC$prev_hit_5))
# prop.table(table(MSC$prev_hit_6))
# prop.table(table(MSC$prev_hit_7))
# prop.table(table(MSC$prev_hit_8))
# 
# ## Some indiciation of an effect
# MSC%>%summarize(mean(p_prev_obj_hit, na.rm = TRUE))


# OG Distance Analysis

# MSC <- MSC%>%mutate(
#   PREV_ORTHOGONAL_DISTANCE_1 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[1]),
#   PREV_ORTHOGONAL_DISTANCE_2 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[2]),
#   PREV_ORTHOGONAL_DISTANCE_3 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[3]),
#   PREV_ORTHOGONAL_DISTANCE_4 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[4]),
#   PREV_ORTHOGONAL_DISTANCE_5 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[5]),
#   PREV_ORTHOGONAL_DISTANCE_6 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[6]),
#   PREV_ORTHOGONAL_DISTANCE_7 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[7]),
#   PREV_ORTHOGONAL_DISTANCE_8 = sapply(PREV_ORTHOGONAL_DISTANCE, function(x) x[8])
# )
# 
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_1, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_2, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_3, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_4, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_5, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_6, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_7, na.rm = TRUE)
# mean(MSC$PREV_ORTHOGONAL_DISTANCE_8, na.rm = TRUE)
# 
# mean(MSC$MEAN_OG_DISTANCE, na.rm = TRUE)
# mean(MSC$CENTER_OG_DISTANCE, na.rm = TRUE)
# mean(MSC$AVG_FIX_1_PREV_OG_DISTANCE, na.rm = TRUE)
# mean(MSC$AVG_FIX_2_PREV_OG_DISTANCE, na.rm = TRUE)
# mean(MSC$AVG_FIX_3_PREV_OG_DISTANCE, na.rm = TRUE)
# mean(MSC$AVG_FIX_4_PREV_OG_DISTANCE, na.rm = TRUE)


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

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_2, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_2_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_3, na.rm = TRUE)
mean(MSC$AVG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_4, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_5, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_6, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_7, na.rm = TRUE)

mean(MSC$PREV_OBJ_FIX_DOT_PRODUCT_8, na.rm = TRUE)

mean(MSC$MEAN_OG_DISTANCE, na.rm = TRUE)
mean(MSC$CENTER_OG_DISTANCE, na.rm = TRUE)
mean(MSC$AVG_FIX_1_PREV_OG_DISTANCE, na.rm = TRUE)
mean(MSC$AVG_FIX_2_PREV_OG_DISTANCE, na.rm = TRUE)
mean(MSC$AVG_FIX_3_PREV_OG_DISTANCE, na.rm = TRUE)
mean(MSC$AVG_FIX_4_PREV_OG_DISTANCE, na.rm = TRUE)


MSC_DOT_PRODUCTS <- MSC[, c("PREV_OBJ_FIX_DOT_PRODUCT_1", "AVG_DIST_FIX_1_PREV_OBJ_DOT_PRODUCT","PREV_OBJ_FIX_DOT_PRODUCT_2", "AVG_DIST_FIX_2_PREV_OBJ_DOT_PRODUCT", "PREV_OBJ_FIX_DOT_PRODUCT_3", "AVG_DIST_FIX_3_PREV_OBJ_DOT_PRODUCT")]


# MSC$diff <- -MSC$PREV_ORTHOGONAL_DISTANCE_2 + MSC$MEAN_OG_DISTANCE  
# MSC %>% ggplot() + geom_point(aes(x=CENTER_OG_DISTANCE, y=diff))

# write.csv(MSC[, c("PREV_ORTHOGONAL_DISTANCE_2", "CENTER_OG_DISTANCE", "MEAN_OG_DISTANCE", "AVG_FIX_OG_DISTANCE")], file = "MUTATED_MSC.csv")
# write.csv(MSC[, c("PREV_ORTHOGONAL_DISTANCE_2", "CENTER_OG_DISTANCE", "MEAN_OG_DISTANCE", "AVG_FIX_OG_DISTANCE")], file = "MUTATED_MSC.csv")
# write_json(MSC, "__data.json", pretty = TRUE)




# Meeting 2.4.25
# New baseline: 
#1 Calculate the distribution of each subject for the i-th fixation angle - similar to avg_fix_i_OG_distance but with different measure and replacing the avg with a distribution
#2 Taking image pairs XY, seeing across subjects where they tend to look?? 


#fix dot vector baseline  issue:
# datastructures problem - list of list of vector -> list of vector
