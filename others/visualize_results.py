from marinetools.estuaries import utils
import matplotlib.pyplot as plt
import numpy as np

dConfig ={"geometryFilename": "Guadalquivir_f21_dx50_h16.txt", "nx": 2098, "nz": 401}


df, db = utils._read_geometry_oldfiles(dConfig)

sections = [1]
sections = np.linspace(0, 2097, 25)
sections = [int(section) for section in sections]

fig, axs = plt.subplots(3, 1)
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
    axs[1].plot(x_across, z_across)

    # Plot cross-sections location (-ve values of right margin)
    x = np.hstack([df.xutm[section] - np.cos(df.AngMd[section])* db["xr"][section, -1], df.xutm[section] + np.cos(df.AngMi[section])* db["xl"][section, -1]])
    y = np.hstack([df.yutm[section] - np.sin(df.AngMd[section])* db["xr"][section, -1], df.yutm[section] + np.sin(df.AngMi[section])* db["xl"][section, -1]])
    axs[0].plot(x, y)
    
plt.show()
