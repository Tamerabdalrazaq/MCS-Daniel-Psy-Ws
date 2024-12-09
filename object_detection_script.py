from ultralytics import YOLO
import cv2
import pandas as pd

DB_NAME = "microwave_fixations.csv"
MCS_FOLDER_BATH = r"C:\Users\averr\OneDrive\Desktop\TAU\Year 4\Psych Workshop\MCS_dataset"
IMAGE_GROUP = "train"

LIMIT = 2000
DETECTION_SAVE_PATH = ".\detections"

def get_img_name_from_path(path):
    return path.split("\\")[-1]

# # Load a model
# model = YOLO("yolo11n-obb.yaml")  # build a new model from YAML
# model = YOLO("yolo11n-obb.pt")  # load a pretrained model (recommended for training)
# model = YOLO("yolo11n-obb.yaml").load("yolo11n.pt")  # build from YAML and transfer weights

# # Train the model
# results = model.train(data="dota8.yaml", epochs=100, imgsz=640)

model = YOLO("yolo11x.pt")

def detect_microwave(image_path):
    # Get class names from the model
    class_names = model.names
    # Find the index of "microwave" in the class names dictionary
    microwave_index = [key for key, value in class_names.items() if value == "microwave"][0]

    # Load image
    # image = cv2.imread(image_path)

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



# Example usage
# images = ['16963.jpg', '89.jpg', '908.jpg', '1355.jpg', '1526.jpg', '4904.jpg', '4978.jpg', '5345.jpg', '9113.jpg', '10909.jpg']
# images = ['16963.jpg', '89.jpg', '908.jpg',]
# for idx, image in enumerate(images):
#     print(f"** Image #{idx} {image} **")
#     detect_microwave(image)
#     print("\n")


def get_pure_img_name(img_name):
    return img_name.split("_")[-1].lstrip('0')

def get_image_full_path(img_name, object, condition): 
    cond = 'TP' if condition == 'present' else 'TA'
    return fr"{MCS_FOLDER_BATH}\images\{IMAGE_GROUP}\{object}\{cond}\{img_name}"

def find_microwaves():
    processed_images = dict()
    def process_image(img_name, condition, object):
        if img_name in processed_images:
            return processed_images[img_name]
        full_path = get_image_full_path(img_name, object, condition)
        detection_res =  detect_microwave(full_path)
        if detection_res == None:
            processed_images[img_name] = ("N\A", "N\A","N\A","N\A","N\A")
            return ("N\A", "N\A","N\A","N\A","N\A")
        processed_images[img_name] = detection_res
        return detection_res
    df = pd.read_csv(DB_NAME)[:LIMIT]
    df['img_name'] = df.apply(lambda row: get_pure_img_name(row['searcharray']), axis=1, result_type='expand')
    df[['OBJECT_X1', 'OBJECT_Y1', 'OBJECT_X2', 'OBJECT_Y2', "confidence"]] = df.apply(lambda row: process_image(row['img_name'], row['condition'], row['catcue']), axis=1, result_type='expand')
    df.to_csv('filtered_data_YOLOV11.csv', index=False)

    print(df)


find_microwaves()
# print(get_pure_img_name("COCO_train2014_000000142896.jpg"))
