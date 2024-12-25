from ultralytics import YOLO
import cv2
import pandas as pd
import math

DB_NAME = "microwave_fixations.csv"
MCS_FOLDER_BATH = r"C:\Users\averr\OneDrive\Desktop\TAU\Year 4\Psych Workshop\MCS_dataset"
IMAGE_GROUP = "train"

START_LIMIT = 0
END_LIMIT = 2000
DETECTION_SAVE_PATH = ".\detections"
SAVED_CSV_FILE_NAME = "filtered_data_YOLOV11_1_2000"

def get_img_name_from_path(path):
    return path.split("\\")[-1]

model = YOLO("yolo11x.pt")

def detect_microwave(image_path):
    # Get class names from the model
    class_names = model.names
    # Find the index of "microwave" in the class names dictionary
    microwave_index = [key for key, value in class_names.items() if value == "microwave"][0]

    # Run detection with filtering to only detect microwaves
    results = model([image_path], classes=[microwave_index])  # Use 'classes' to limit detection to microwaves only
    # Iterate through results and print bounding box pixel ranges
    if len(results[0].boxes.cls ) < 1:
        return None 
    detection_box = results[0].boxes

    x1, y1, x2, y2 = map(int, detection_box.xyxy[0])  # Bounding box coordinates
    confidence = detection_box.conf[0].item()         # Confidence score as float
    # Print microwave bounding box details
    # results[0].show()
    results[0].save(filename= f"{DETECTION_SAVE_PATH}/det_{get_img_name_from_path(image_path)}")
    return (x1, y1, x2, y2, round(confidence, 2))



def get_pure_img_name(img_name):
    return img_name.split("_")[-1].lstrip('0')

def get_image_full_path(img_name, object, condition): 
    cond = 'TP' if condition == 'present' else 'TA'
    return fr"{MCS_FOLDER_BATH}\images\{IMAGE_GROUP}\{object}\{cond}\{img_name}"

def find_microwaves():
    df = pd.read_csv(DB_NAME)[START_LIMIT:END_LIMIT]
    total_rows = END_LIMIT - START_LIMIT
    img_n_a = ("N\A", "N\A","N\A","N\A","N\A")
    processed_images = dict()
    current_percentage = 0
    def process_image(img_name, condition, object, index):
        nonlocal current_percentage
        if math.floor((index/total_rows)*100) >= current_percentage + 1:
            current_percentage = current_percentage + 1
            print( f"***********************") 
            print( f"{current_percentage}%") 
            print( f"***********************") 
        if condition == "absent":
            return img_n_a
        if img_name in processed_images:
            return processed_images[img_name]
        full_path = get_image_full_path(img_name, object, condition)
        detection_res =  detect_microwave(full_path)
        if detection_res == None:
            processed_images[img_name] = img_n_a
            return img_n_a
        processed_images[img_name] = detection_res
        return detection_res
    df['img_name'] = df.apply(lambda row: get_pure_img_name(row['searcharray']), axis=1, result_type='expand')
    df[['OBJECT_X1', 'OBJECT_Y1', 'OBJECT_X2', 'OBJECT_Y2', "confidence"]] = df.apply(lambda row: process_image(row['img_name'], row['condition'], row['catcue'], row.name), axis=1, result_type='expand')
    df.to_csv(f'{SAVED_CSV_FILE_NAME}.csv', index=False)

    print(df)



def merge_chuncks():
    # Paths to the two CSV files
    file1 = "filtered_data_YOLOV11_1_2000_25_12_24.csv"  # Contains rows 1-2000
    file2 = "filtered_data_YOLOV11_25_12_24.csv"  # Contains rows 2001-end

    # Read both CSV files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Concatenate the DataFrames vertically
    merged_df = pd.concat([df1, df2], ignore_index=True)

    # Save the merged DataFrame to a new CSV file
    merged_df.to_csv("merged_dataset.csv", index=False)

    print("Merged dataset saved as 'merged_dataset.csv'")


def main():
    merge_chuncks()

main()
