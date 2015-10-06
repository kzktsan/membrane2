import cv2
import numpy as np
import glob
import os

FOLDER_NAME = "/Users/satokazuki/Desktop/zikken5/G1/"
OUTPUT_FOLDER = "/Users/satokazuki/Desktop/zikken5/G1_kirinuki/G1_"
LENGTH = 600
HEIGHT_STRIDE = 200
WIDTH_STRIDE = 200

folder = [os.path.basename(r) for r in glob.glob(FOLDER_NAME + "*.png")]
folder_conter = 0

for file in folder :

	print file 
	img = cv2.imread(FOLDER_NAME + file)
	height = img.shape[0]
	width = img.shape[1]

	sp_x = 0
	sp_y = 0
	sp_counter = 0

	while True :
		sp_x = 0
		while True : 
			copy = img[sp_y : sp_y + LENGTH , sp_x : sp_x + LENGTH]
			cv2.imwrite(OUTPUT_FOLDER + str(folder_conter) + "_" + str(sp_counter) + ".png", copy)
			sp_counter += 1
			sp_x += WIDTH_STRIDE
			if (sp_x + LENGTH) > width:
				break

		sp_y += HEIGHT_STRIDE
		if (sp_y + LENGTH) > height:
			break
	folder_conter += 1
	
