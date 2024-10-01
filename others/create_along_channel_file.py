import pandas as pd

data = pd.read_csv("along_channel_data.csv")

data["S"] = 20
data.to_csv("_along_channel_data.csv")
# ;x,Daveraged,D90,D50,sigma,rhos,thickness