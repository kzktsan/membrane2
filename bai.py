import cv2
import numpy as np

input_name = "/Users/satokazuki/Desktop/moto.png"
output_name = "/Users/satokazuki/Desktop/moto2.png"

img = cv2.imread(input_name, 0)

height = img.shape[0]
width = img.shape[1]

twice_img = cv2.resize(img,(255,255))

cv2.imwrite(output_name, twice_img, 0)