import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from marinetools.estuaries import utils

dConfig = {"geometryFilename": "Guadalquivir_f21_dx50_h16.txt", "nx": 2098, "nz": 401}


df, db = utils._read_geometry_oldfiles(dConfig)
# data = pd.read_csv("along_channel_data_GRE.csv")

sections = [1]
sections = np.linspace(0, 2097, 10)
sections = [int(section) for section in sections]

fig, axs = plt.subplots(1, 3)
# Plot along channel location
x = df.xutm
y = df.yutm
axs[0].plot(x, y)


# Plot along channel depth
x = df.x
z = df.z
axs[2].plot(x, z)

for section in sections:
    # Plot cross-sections depth
    x_across = np.hstack([db["xr"][section, ::-1], db["xl"][section, 1:]])
    z_across = np.hstack([db["eta"][section, ::-1], db["eta"][section, 1:]])
    axs[1].plot(x_across, z_across + z[section])

    # Plot cross-sections location (-ve values of right margin)
    x = np.hstack(
        [
            df.xutm[section] - np.cos(df.AngMd[section]) * db["xr"][section, -1],
            df.xutm[section] + np.cos(df.AngMi[section]) * db["xl"][section, -1],
        ]
    )
    y = np.hstack(
        [
            df.yutm[section] - np.sin(df.AngMd[section]) * db["xr"][section, -1],
            df.yutm[section] + np.sin(df.AngMi[section]) * db["xl"][section, -1],
        ]
    )
    axs[0].plot(x, y)

plt.show()
