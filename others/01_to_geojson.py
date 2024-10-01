from marinetools.estuaries import utils
import numpy as np
import geopandas as gpd
from shapely.geometry import LineString


dConfig ={"geometryFilename": "Guadalquivir_f21_dx50_h16.txt", "nx": 2098, "nz": 401}

df, db = utils._read_geometry_oldfiles(dConfig)
sections = np.arange(2097)

x = df.xutm
y = df.yutm

#Crear la línea a partir de las coordenadas
linea = LineString(df.loc[:, ["xutm", "yutm"]])

# Crear un GeoDataFrame con la línea
gdf = gpd.GeoDataFrame(geometry=[linea])

# Asignar el sistema de referencia WGS84 (EPSG:4326)
gdf.set_crs(epsg=25829, inplace=True)
gdf.to_file("thalweg_location.geojson")

lineas = []
for section in sections:
    x = np.hstack([df.xutm[section] - np.cos(df.AngMd[section])* db["xr"][section, -1], df.xutm[section] + np.cos(df.AngMi[section])* db["xl"][section, -1]])
    y = np.hstack([df.yutm[section] - np.sin(df.AngMd[section])* db["xr"][section, -1], df.yutm[section] + np.sin(df.AngMi[section])* db["xl"][section, -1]])
    linea = np.vstack([x, y]).T
    linea = LineString(linea)
    lineas.append(linea)
    

# Crear un GeoDataFrame con las líneas
gdf = gpd.GeoDataFrame(geometry=lineas)

# Asignar un sistema de referencia (WGS84, EPSG:4326)
gdf.set_crs(epsg=25829, inplace=True)
gdf.to_file("sections.geojson")