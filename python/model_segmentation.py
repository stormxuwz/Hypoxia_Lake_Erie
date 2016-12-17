import numpy as np
import logging

class timeSeriesSegmentation(object):
	def __init__(self,max_error):
		self.segmentList=None
		self.x=None
		self.max_error=max_error

	def fit_predict(self,x):
		pass

	def calculate_error(self,x,y):
		return np.max(np.abs(x-y))

	def plot(self):
		import matplotlib.pyplot as plt
		plt.figure()
		plt.plot(self.x,"r")
		for seg in self.segmentList:
			plt.plot(seg[1],seg[0],"+-")
		plt.show()

	def createLine(self,x,method="regression"):
		n=len(x)
		if method=="simple":
			line = np.linspace(x[0],x[-1],n)

		elif method == "regression":
			line = np.poly1d(np.polyfit(range(n),x,1))(range(n))
		
		elif method == "poly":
			line=np.poly1d(np.polyfit(range(n),x,2))(range(n))

		return(line)


class slidingWindow(timeSeriesSegmentation):
	def fit_predict(self,x):
		n=len(x)
		leftNode=0
		segmentList=[]
		print n
		
		while leftNode<n-1:
			print leftNode
			newSeg = False
			for rightNode in range(leftNode+3,n):
				testSeg=x[leftNode:rightNode]

				testLine = self.createLine(testSeg,"regression")
				segError = self.calculate_error(testSeg, testLine)

				if segError>self.max_error:
					segmentList.append([testLine,range(leftNode,rightNode)])
					leftNode=rightNode
					newSeg = True
					break

			if newSeg is False:
				segmentList.append([testLine,range(leftNode,rightNode)])
				leftNode = n-1
		
		self.segmentList=segmentList
		self.x=x

class bottomUp_slow(timeSeriesSegmentation):
	def fit_predict(self,x):
		n=len(x)
		segmentIndexList = [[i,i+1] for i in range(0,n,2)]
		errorList = [self.mergeCost(x[segmentIndexList[i]], x[segmentIndexList[i+1]]) for i in range(len(segmentIndexList)-1)]
		
		while True:
			minIndex = np.argmin(errorList)
			if len(segmentIndexList) ==3:
				break
			new_segmentIndexList = self.mergeRight(segmentIndexList,minIndex)

			errorList = [self.mergeCost(x[new_segmentIndexList[i]], x[new_segmentIndexList[i+1]]) for i in range(len(new_segmentIndexList)-1)]
			if min(errorList)<self.max_error:
				segmentIndexList= new_segmentIndexList
			else:
				break
		
		logging.info("errorList")
		logging.info(errorList)

		self.x=x
		self.segmentList=[[self.createLine(x[segIndex]),segIndex] for segIndex in segmentIndexList]

	def mergeCost(self, leftSeg,rightSeg):
		allSeg=np.concatenate((leftSeg,rightSeg))
		line = self.createLine(allSeg)
		return self.calculate_error(line, allSeg)


	def mergeRight(self,segList,index):
		segList[index]=segList[index]+segList[index+1] # merge 
		segList.pop(index+1) # pop the right segment
		return(segList)


class bottomUp(timeSeriesSegmentation):
	def fit_predict(self,x):
		n=len(x)
		segmentIndexList = [[i,i+1] for i in range(0,n,2)]
		errorList = [self.mergeCost(x[segmentIndexList[i]], x[segmentIndexList[i+1]]) for i in range(len(segmentIndexList)-1)]

		while True:
			minIndex = np.argmin(errorList)
			new_segmentIndexList = self.mergeRight(segmentIndexList,minIndex)
			
			if len(segmentIndexList) == 3:
				break
			if minIndex > 0:
				errorList[minIndex-1] = self.mergeCost(x[new_segmentIndexList[minIndex-1]],x[new_segmentIndexList[minIndex]])

			if minIndex < len(errorList)-1:
				errorList[minIndex+1] = self.mergeCost(x[new_segmentIndexList[minIndex]], x[new_segmentIndexList[minIndex+1]])

			errorList.pop(minIndex)
			
			if len(errorList)!=len(new_segmentIndexList)-1:
				raise ValueError("error length not right")

			if min(errorList)<self.max_error:
				segmentIndexList= new_segmentIndexList
			else:
				break

		self.x=x
		self.segmentList=[[self.createLine(x[segIndex]),segIndex] for segIndex in segmentIndexList]

	def mergeCost(self, leftSeg,rightSeg):
		allSeg=np.concatenate((leftSeg,rightSeg))
		line = self.createLine(allSeg)
		return self.calculate_error(line, allSeg)

	def mergeRight(self,segList,index):
		segList[index]=segList[index]+segList[index+1] # merge 
		segList.pop(index+1) # pop the right segment
		return(segList)
