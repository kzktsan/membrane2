# coding: utf-8

import cv, commands, sys, os

FOLDER_PATH = "/Users/satokazuki/Desktop/membrane/"

img_list = os.listdir(FOLDER_PATH)
print img_list
count = 0

for img_name in img_list : 
	angle = (90, 180, 270)
	item = img_name.split(".")

	if item[1] == "jpg" :
		im_in = cv.LoadImage(FOLDER_PATH + img_name)
		im_ro = cv.CreateImage(cv.GetSize(im_in), cv.IPL_DEPTH_8U, 3)
		rotate_mat = cv.CreateMat(2, 3, cv.CV_32FC1)
		count += 1

		for i in range(0, 3):
			cv.GetRotationMatrix2D((im_in.height/2, im_in.width/2), angle[i], 1, rotate_mat)
			cv.WarpAffine(im_in, im_ro, rotate_mat)
			cv.SaveImage(FOLDER_PATH + item[0] + "_" + str(i) + ".jpg", im_ro)
			count += 1

print count