import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

data = nc.Dataset("../cmake-build-debug/test_003.nc")


water_flow = data.variables["Q"]
eta = data.variables["eta"]
area = data.variables["A"]

x = data.variables["x"][:].data
t = data.variables["time"][:].data

time_to_plot = np.linspace(0, eta.shape[0] - 1, 12)
fig, axs = plt.subplots(1, 3)
for i in time_to_plot:
    i = int(i)
    axs[0].plot(x, eta[i, :], label="time: " + str(t[i]))
    axs[1].plot(x, water_flow[i, :])
    axs[2].plot(x, area[i, :])
axs[0].legend(ncols=2)
plt.show()
