import matplotlib.cm as cm
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np

plt.rc("text", usetex=True)
plt.rc("font", family="serif", size=10)
params = {"text.latex.preamble": [r"\usepackage{amsmath}"]}

data0 = nc.Dataset("../cmake-build-debug/test_003.nc")
data1 = nc.Dataset("../cmake-build-debug/test_004.nc")


water_flow0 = data0.variables["Q"][:].data
eta0 = data0.variables["eta"][:].data
area0 = data0.variables["A"][:].data
sal0 = data0.variables["S"][:].data

water_flow1 = data1.variables["Q"][:].data
eta1 = data1.variables["eta"][:].data
area1 = data1.variables["A"][:].data
sal1 = data1.variables["S"][:].data

x = data0.variables["x"][:].data
t = data0.variables["time"][:].data

time_to_plot = np.linspace(0, eta1.shape[0] - 1, 12)
eta_dif = eta1[:, :] - eta0[:, :]
# mask = abs(eta_dif) < 1e-2
# eta_dif[mask] = np.nan

dif_water_flow = water_flow1[:, :] - water_flow0[:, :]
# mask = abs(dif_water_flow) < 20
# dif_water_flow[mask] = np.nan

dif_sal = sal1[:, :] - sal0[:, :]

x, t = np.meshgrid(x, t)
plt.figure()
plt.contourf(t / 3600, x / 1000, eta_dif, 100, cmap=cm.RdYlBu)
plt.colorbar(label=r"$\mathbf{\eta (m)}$")
plt.contour(t / 3600, x / 1000, eta_dif, [-1e-2, 1e-2])
plt.xlabel(r"\textbf{time (h)}")
plt.ylabel(r"\textbf{x (km)")
plt.grid(True)

plt.figure()
plt.contourf(t / 3600, x / 1000, dif_water_flow, 100, cmap=cm.RdYlBu)
plt.colorbar(label=r"\textbf{Q ($m^3/s$)}")
plt.contour(t / 3600, x / 1000, dif_water_flow, [-10, 10])
plt.xlabel(r"\textbf{time (h)}")
plt.ylabel(r"\textbf{x (km)}")
plt.grid(True)

plt.figure()
plt.contourf(t / 3600, x / 1000, area1[:, :] - area0[:, :], 100, cmap=cm.RdYlBu)
plt.xlabel(r"\textbf{time (h)}")
plt.ylabel(r"\textbf{x (km)}")
plt.grid(True)
plt.colorbar(label=r"\textbf{A (m2)}")

plt.figure()
plt.contourf(t / 3600, x / 1000, sal1[:, :] - sal0[:, :], 100, cmap=cm.RdYlBu)
plt.xlabel(r"\textbf{time (h)}")
plt.ylabel(r"\textbf{x (km)}")
plt.colorbar(label=r"\textbf{s (psu)}")
plt.grid(True)
plt.show()
