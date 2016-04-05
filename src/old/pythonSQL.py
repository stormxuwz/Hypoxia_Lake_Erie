from __future__ import division
import pandas as pd
from sqlalchemy import create_engine
import os

sheet = pd.read_excel("/Users/WenzhaoXu/Developer/Project_IO/DO_LakeErie/2015Data/Station_Info.xlsx",sheetname="2015")
# print sheet

for feature in ["Latitude","Longitude"]:
	for i, item in enumerate(sheet[feature]):
		if isinstance(item,unicode):
			data = item.encode("utf-8")
			if "N" in data:
				# print data
				data=data.split("N")[1].split("\xc2\xb0")
				data=float(data[0])+float(data[1])/60

			elif "W" in data:
				data=data.split("W")[1].split("\xc2\xb0")
				data=float(data[0])+float(data[1])/60
			else:
				data=data.split(" ")
				data=float(data[0])+float(data[1])/60
			# print data
		else:
			data = item
		
		# print sheet[feature][i],data

		sheet[feature][i]=data

sheet.rename(columns={"Longitude":"longitude","Latitude":"latitude","LoggerID":"loggerId","Position":"loggerPosition"},inplace=True)
sheet.longitude=-1*sheet.longitude



def readDO_data():
	for root,dirs,files in os.walk("/Users/WenzhaoXu/Developer/Project_IO/DO_LakeErie/2015Data/"):
		for f in files:
			if f.endswith(".csv"):
				print f
				if "_" in f:
					logger=f.split("_")[0]
				else:
					logger=f.split(".")[0]
				print logger
				data=pd.read_csv(os.path.join(root,f),header=1)
				data=data.iloc[5:-5,]
				colNum = len(data.columns)
				colName=["ind","Time","DO","Temp"]+["tmp_"+str(i) for i in range(colNum-4)]
				data.columns=colName
				data["Time"]=pd.to_datetime(data["Time"],format="%m/%d/%y %I:%M:%S %p")
				data["logger"]=logger
				# print data.columns
				data[["Time","logger","DO","Temp"]].to_sql("loggerData",SQL_engine1,flavor="mysql",if_exists="append",index=False)



SQL_engine1 = create_engine("mysql+mysqldb://root:XuWenzhaO@localhost/DO2015")
# SQL_engine2 = create_engine("mysql+mysqldb://root:XuWenzhaO@localhost/DO")

# sheet.to_sql("loggerInfo",SQL_engine1,flavor="mysql",if_exists="replace")
readDO_data()


## Plot




