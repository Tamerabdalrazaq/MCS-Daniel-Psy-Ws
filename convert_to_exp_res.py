import cv2
import numpy as np
import pandas as pd

df = pd.read_json("R scripts/__data.json")

MCS_FOLDER_BATH = r"C:\Users\averr\OneDrive\Desktop\TAU\Year 4\Semestar A\Psych Workshop\MCS_dataset" # Absolute path for the folder containing DB_NAME
IMAGE_GROUP = "train" # Group in ["train", "test"]

def get_image_full_path(img_name, object, condition): 
    cond = 'TP' if condition == 'present' else 'TA'
    return fr"{MCS_FOLDER_BATH}\images\{IMAGE_GROUP}\{object}\{cond}\{img_name}"

def add_frame(image_path, output_path, target_width=1280, target_height=800, frame_color=(0, 0, 0), prev_p1 = None, prev_p2 = None, points=None, point_color=(0, 255, 0), point_radius=8):
    # Load the image
    image = cv2.imread(image_path)
    if image is None:
        print(image_path)
        return
        # raise ValueError("Could not load image. Check the file path.")
    
    h, w, _ = image.shape  # Get current dimensions
    
    # Compute necessary padding
    top_pad = max(0, (target_height - h) // 2)
    bottom_pad = max(0, target_height - h - top_pad)
    left_pad = max(0, (target_width - w) // 2)
    right_pad = max(0, target_width - w - left_pad)
    
    # Add padding to the image
    framed_image = cv2.copyMakeBorder(image, top_pad, bottom_pad, left_pad, right_pad, cv2.BORDER_CONSTANT, value=frame_color)
    
    # Adjust points based on padding and mark them
    if points:
        for idx, (x, y) in enumerate(points):
            red = (0, 0, 255)
            adjusted_x = int(x)
            adjusted_y = int(y)
            # Draw the circle with the interpolated color
            cv2.circle(framed_image, (adjusted_x, adjusted_y), point_radius, red, -1)
            
            # Add index label near the point
            cv2.putText(framed_image, str(idx + 1), (adjusted_x + 10, adjusted_y - 10), 
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (255, 255, 255), 1, cv2.LINE_AA)
    
    if prev_p1 is not None and prev_p2 is not None:
        cv2.rectangle(framed_image, prev_p1, prev_p2, (0, 0, 255), 2) 
    # Save the result
    cv2.imwrite(output_path, framed_image)
    # print(f"Framed image with points saved to {output_path}")


for index, row in df.iterrows():
    img_name = row['img_name']
    image_path = "detections\det_"+img_name  # Assuming column for image filename
    points = ([(row['CURRENT_FIX_X'][i],row['CURRENT_FIX_Y'][i]) for i in range(len(row['CURRENT_FIX_Y']))])  # Assuming column with list of (x, y) tuples
    # Call the function
    img_identifier = f"{row["subjectnum"]}_{row["TRIAL_INDEX"]}_{img_name}"
    output_path = f".\\fixations\\{img_identifier}"  # Create unique output filenames
    if(row['condition'] == 'absent'):
        image_path = get_image_full_path(img_name, "microwave", 'absent')
    prev_p1 = None
    prev_p2 = None
    if not np.isnan(row['PREV_OBJ_X1']) and not np.isnan(row['PREV_OBJ_Y1']):
        prev_p1 = (int(row['PREV_OBJ_X1']), int(row['PREV_OBJ_Y1']))
        prev_p2 = (int(row['PREV_OBJ_X2']), int(row['PREV_OBJ_Y2']))
    add_frame(image_path, output_path, points=points, prev_p1=prev_p1, prev_p2=prev_p2)
    # if index > 500:
    #     break
print('done.')