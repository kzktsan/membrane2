import cv2
import numpy as np
import random

#img = cv2.imread("/Users/satokazuki/Desktop/lined2.jpg")
#img2 = cv2.imread("/Users/satokazuki/Desktop/ori.jpg")
img = cv2.imread("image/lined3.jpg")
img2 = cv2.imread("image/ori.jpg", 0)
length = 126
half = length/2
#f = open("/Users/satokazuki/Desktop/redlist.txt", 'w')

height = img.shape[0]
width = img.shape[1]

count = 0

while count < 140000:
	x = random.randint(0,width-1)
	y = random.randint(0,height-1)
	if not(img[y, x, 2] > 160 and img[y, x, 1] < 100 and img[y, x, 0] <100) and not(img[y, x, 2] < 10 and img[y, x, 1] < 10 and img[y, x, 0] <10):
			if (y - half >= 0 and y + half < height and x - half >= 0 and x + half < width):
				copy = img2[y - half : y + half +1 , x - half : x + half +1]
				cv2.imwrite("/Users/satokazuki/Desktop/127non/" + str(x) + "_"+str(y)+ "_not.png", copy)
				count = count + 1
'''

for y in range(755, 1011):
	for x in range(130, 386):
		if not(img[y, x, 2] > 160 and img[y, x, 1] < 100 and img[y, x, 0] <100):
			if (y - half >= 0 and y + half < height and x - half >= 0 and x + half < width):
				copy = img2[y - half : y + half +1 , x - half : x + half +1]
				cv2.imwrite("/Users/satokazuki/Desktop/127non/" + str(x) + "_"+str(y)+ "_not.png", copy)
				count = count + 1
'''

print(count)