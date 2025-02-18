import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

plt.rc("text", usetex=True)


data = nc.Dataset("../cmake-build-debug/test_003.nc")


water_flow = data.variables["Q"]
eta = data.variables["eta"]
area = data.variables["A"]

x = data.variables["x"][:].data
t = data.variables["time"][:].data

print(np.shape(eta))
fig, axs = plt.subplots(1, 1)
for i in np.arange(1, 5) * 50:
    axs.plot(area[:, i], water_flow[:, i], label=str(i))
plt.xlabel(r"$\mathbf{A (m^2)}$")
plt.ylabel(r"$\mathbf{Q (m^3/s)}$")
plt.grid()
plt.legend(title="Section \#")
plt.show()
