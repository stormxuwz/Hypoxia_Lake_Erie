from sqlalchemy import create_engine
import pandas as pd
import numpy as np
from model_segmentation import bottomUp
import matplotlib.pyplot as plt

# dbConfig <- list(dbname = "DO", username="root", password="XuWenzhaO", host="do.cm1qoaxjisxm.us-west-2.rds.amazonaws.com")
# varUnit <- list(DO="DO(mg/L)",Temp="Temperature(C)")


# engine = create_engine('mysql+mysqldb://root:XuWenzhaO@do.cm1qoaxjisxm.us-west-2.rds.amazonaws.com/DO')

# sql = "Select date(Time) as Time, AVG(DO) as DO, logger from loggerData_2014 where (logger = 10523446 OR logger = 10523447 OR logger = 10523439 OR logger = 10523450 OR logger = 10523436 OR logger = 10528849 OR logger = 10528846 OR logger = 10523443 OR logger = 10523437 OR logger = 10523445 OR logger = 10384436 OR logger = 10384445 OR logger = 10384437 OR logger = 10384438 OR logger = 10384443 OR logger = 10384449 OR logger = 10461951 OR logger = 10384439) Group by date(Time),logger"

# a = pd.read_sql_query(sql,engine)


fname = "/Users/WenzhaoXu/Desktop/2014_DO_data_hourly.csv"
model = bottomUp(max_error = 0.5)

data_2014 = pd.read_csv(fname)


# for i in range(data_2014.shape[1]-1):

test = np.array(data_2014.iloc[:,3].dropna())
# print test

if len(test) % 2 !=0:
	test = test[:-1]
model.fit_predict(test)

model.plot()
print len(model.segmentList)





