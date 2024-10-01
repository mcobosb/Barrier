import pandas as pd

data = pd.read_csv("along_channel_data.csv")

sediments = pd.DataFrame(-1.0, index=data[";x"].values, columns=["Daveraged", "D90", "D50", "sigma", "rhos", "thickness"])

for i in sediments.index:
    sediments.loc[i, :] = [0.0141, 0.02791, 0.00754, 3.7347, 2650.0, 10.0]

sediments.to_csv("sediments")
# ;x,Daveraged,D90,D50,sigma,rhos,thickness
