import cv2
import numpy as np

#img = cv2.imread("/Users/satokazuki/Desktop/lined2.jpg")
#img2 = cv2.imread("/Users/satokazuki/Desktop/ori.jpg")
img = cv2.imread("image/lined3.jpg")
img2 = cv2.imread("image/ori.jpg", 0)
#length = 64
#length = 254
length = 126
half = length/2
#f = open("/Users/satokazuki/Desktop/redlist.txt", 'w')

height = img.shape[0]
width = img.shape[1]

count = 0
'''
for y in range(0, height - 1):
	for x in range(0, width - 1):
		if (img[y, x, 2] > 160 and img[y, x, 1] < 100 and img[y, x, 0] <100):
			if (y - half >= 0 and y + half < height and x - half >= 0 and x + half < width):
				copy = img2[y - half : y + half +1 , x - half : x + half +1]
				cv2.imwrite("/Users/satokazuki/Desktop/127ori/127cellmem/" + str(x) + "_"+str(y)+ ".png", copy)
				count = count + 1
			#f.write(str(x) + " " + str(y) + " " + str(img[y, x, 2]) + "\n")
'''
for y in range(461, 717):
	for x in range(314, 570):
		if (img[y, x, 2] > 200 and img[y, x, 1] < 100 and img[y, x, 0] <100):
			if (y - half >= 0 and y + half < height and x - half >= 0 and x + half < width):
				copy = img2[y - half : y + half +1 , x - half : x + half +1]
				cv2.imwrite("/Users/satokazuki/Desktop/zikken4/zikken4_cellmem/" + str(x) + "_"+str(y)+ ".png", copy)
				count = count + 1

print(count)
