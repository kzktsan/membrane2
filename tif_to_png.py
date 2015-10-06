import cv2
import numpy as np
import glob
import os

FOLDER_NAME = "/Users/satokazuki/Desktop/watanabe/G3/"
OUTPUT_FOLDER = "/Users/satokazuki/Desktop/zikken5/G3/"

folder = [os.path.basename(r) for r in glob.glob(FOLDER_NAME + "*")]

for file in folder :
	print file
	img = cv2.imread(FOLDER_NAME + file)
	tmp = file.split(".")
	cv2.imwrite(OUTPUT_FOLDER + tmp[0] + ".png", img)
