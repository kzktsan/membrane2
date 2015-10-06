import cv2
import numpy as np
import glob
import os

FOLDER_NAME = "/Users/satokazuki/Desktop/zikken5/G3_database/"
OUTPUT_FOLDER = "/Users/satokazuki/Desktop/zikken5/G3_database/"

folder = [os.path.basename(r) for r in glob.glob(FOLDER_NAME + "*.png")]

for file in folder : 
	img = cv2.imread(FOLDER_NAME + file)
	tmp = file.split(".")
	tmp2 = tmp[0].split("_")
	print tmp2
	if int(tmp2[3]) < 3 : 
		flip1 = cv2.flip(img, 0)
		flip2 = cv2.flip(img, 1)

		cv2.imwrite(OUTPUT_FOLDER + tmp[0] + "_f1.png", flip1)
		cv2.imwrite(OUTPUT_FOLDER + tmp[0] + "_f2.png", flip2)
