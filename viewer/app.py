import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

from datetime import datetime, timedelta

import dash_bootstrap_components as dbc
import dash_leaflet as dl
import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import xarray as xr
from dash import Dash, Input, Output, callback, dcc, html
from scipy.interpolate import interp1d
import shutil
import tempfile
import uuid
import time
import atexit

# Global variables for dataset and data caching
ds_global = None
nc_data_cache = None
df_global = None

def load_netcdf_data():
    """Load and process NetCDF data - keep dataset open globally"""
    global ds_global
    
    if ds_global is not None:
        print("Using existing NetCDF dataset (cached)")
        return ds_global
    
    # Get NetCDF path from environment variable
    primary_path = os.getenv('NETCDF_PATH')
    
    # Fallback paths
    possible_paths = []
    if primary_path:
        possible_paths.append(primary_path)
    
    possible_paths.extend([
        "../tests/test005/test_005.nc",
        "/home/mcobos/Barrier/tests/test005/test_005.nc",
    ])
    
    # Find the source file
    source_path = None
    for nc_path in possible_paths:
        print(f"Checking for NetCDF file at: {nc_path}")
        if os.path.exists(nc_path):
            source_path = nc_path
            print(f"Found NetCDF source file at: {source_path}")
            break
    
    if source_path is None:
        print("ERROR: NetCDF source file not found")
        return None
    
    # Try direct access first
    try:
        print(f"Attempting to open NetCDF file directly: {source_path}")
        ds_global = xr.open_dataset(source_path)
        print(f"Successfully loaded NetCDF file directly")
        print(f"NetCDF dimensions: {dict(ds_global.sizes)}")
        print(f"NetCDF variables: {list(ds_global.data_vars)}")
        return ds_global
        
    except Exception as direct_error:
        print(f"Direct access failed: {direct_error}")
        print("Copying to temporary location...")
        
        # Copy to Windows temp directory
        temp_dir = tempfile.gettempdir()
        unique_id = str(uuid.uuid4())[:8]
        temp_filename = f"guadalquivir_netcdf_{unique_id}.nc"
        temp_path = os.path.join(temp_dir, temp_filename)
        
        try:
            # Remove existing temp file if exists
            if os.path.exists(temp_path):
                os.remove(temp_path)
            
            print(f"Copying from {source_path} to {temp_path}")
            source_size = os.path.getsize(source_path)
            
            # Copy with retry logic
            for attempt in range(3):
                try:
                    shutil.copy2(source_path, temp_path)
                    break
                except Exception as e:
                    if attempt < 2:
                        print(f"Copy attempt {attempt + 1} failed, retrying...")
                        time.sleep(3)
                    else:
                        raise e
            
            # Verify copy
            if os.path.getsize(temp_path) != source_size:
                raise IOError("File copy incomplete")
            
            print(f"✅ File copied successfully")
            
            # Load from temp location
            ds_global = xr.open_dataset(temp_path)
            ds_global.attrs['_temp_file_path'] = temp_path
            
            print(f"✅ Successfully loaded from temporary location")
            return ds_global
            
        except Exception as copy_error:
            print(f"❌ Failed to copy and open file: {copy_error}")
            return None

def calculate_river_distances(river_points):
    """Calculate cumulative distances from Alcalá dam (0) to river mouth"""
    distances = []
    
    # Start with 0 at the dam (first point)
    distances.append(0)
    
    # Calculate distances from dam to mouth
    for i in range(len(river_points) - 1):
        lat1, lon1 = river_points[i]
        lat2, lon2 = river_points[i + 1]
        
        # Haversine formula
        R = 6371000  # Earth radius in meters
        dlat = np.radians(lat2 - lat1)
        dlon = np.radians(lon2 - lon1)
        a = (np.sin(dlat/2)**2 + 
             np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.sin(dlon/2)**2)
        c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1-a))
        distance = R * c
        
        # Add cumulative distance
        distances.append(distances[-1] + distance)
    
    print(f"River distances: 0 m (dam) to {max(distances):.1f} m (mouth)")
    return distances

def interpolate_netcdf_to_river(ds, river_points, river_distances):
    """Interpolate NetCDF data to river points"""
    global nc_data_cache
    
    if nc_data_cache is not None:
        print("Using cached NetCDF interpolation data")
        return nc_data_cache
    
    if ds is None:
        return None
    
    try:
        print("Processing NetCDF interpolation...")
        
        # Define variables to interpolate
        vars_to_interpolate = ['A', 'Q', 'B', 'eta', 'rho', 'U', 'c', 'S', 'Qb', 'Qs']
        
        # Find coordinate names
        coord_names = list(ds.coords)
        distance_coord = None
        time_coord = None
        
        for coord in coord_names:
            if any(x in coord.lower() for x in ['distance', 'x', 's']):
                distance_coord = coord
            elif any(x in coord.lower() for x in ['time', 't']):
                time_coord = coord
        
        if not distance_coord:
            print("ERROR: No distance coordinate found")
            return None
        
        # Get NetCDF distances (from dam, same as our river distances)
        nc_distances_from_dam = ds[distance_coord].values.copy()
        
        print(f"NetCDF distances from dam: {nc_distances_from_dam.min():.1f} to {nc_distances_from_dam.max():.1f} m")
        print(f"River distances from dam: {min(river_distances):.1f} to {max(river_distances):.1f} m")
        
        # Check if distance ranges are compatible
        river_min, river_max = min(river_distances), max(river_distances)
        nc_min, nc_max = nc_distances_from_dam.min(), nc_distances_from_dam.max()
        
        print(f"Distance range comparison:")
        print(f"  River: {river_min:.1f} to {river_max:.1f} m")
        print(f"  NetCDF: {nc_min:.1f} to {nc_max:.1f} m")
        
        if abs(river_max - nc_max) > 5000:  # More than 5km difference
            print(f"WARNING: Large difference in maximum distances.")
            print(f"  River max: {river_max/1000:.1f} km")
            print(f"  NetCDF max: {nc_max/1000:.1f} km")
        
        # Handle time coordinate
        if time_coord and time_coord in ds.coords:
            time_seconds = ds[time_coord].values.copy()
            base_date = datetime(2025, 1, 10, 0, 0, 0)
            time_values = [base_date + timedelta(seconds=float(t)) for t in time_seconds]
            print(f"Time range: {time_values[0]} to {time_values[-1]} ({len(time_values)} steps)")
        else:
            time_values = [datetime(2025, 1, 10, 0, 0, 0)]
        
        # Process each variable
        all_variables_data = {}
        
        for var_name in vars_to_interpolate:
            # Find variable in NetCDF
            nc_var_name = None
            if var_name in ds.data_vars:
                nc_var_name = var_name
            else:
                for ds_var in ds.data_vars:
                    if ds_var.lower() == var_name.lower():
                        nc_var_name = ds_var
                        break
            
            if not nc_var_name:
                continue
                
            data_var = ds[nc_var_name]
            
            if distance_coord not in data_var.dims:
                continue
            
            # Load all data into memory
            if time_coord and time_coord in data_var.dims and len(time_values) > 1:
                all_data = data_var.values.copy()
            else:
                all_data = data_var.values.copy()
                if len(all_data.shape) == 1:
                    all_data = all_data.reshape(1, -1)
            
            var_time_series = []
            
            # Process each time step
            for i, timestamp in enumerate(time_values):
                try:
                    if len(all_data.shape) > 1 and i < all_data.shape[0]:
                        data_values = all_data[i, :]
                    else:
                        data_values = all_data
                    
                    # Ensure we have the right length
                    if len(data_values) > len(nc_distances_from_dam):
                        data_values = data_values[:len(nc_distances_from_dam)]
                    elif len(data_values) < len(nc_distances_from_dam):
                        print(f"Warning: Data length mismatch for {nc_var_name} at time {i}")
                        continue
                    
                    # Handle NaN values
                    if np.any(np.isnan(data_values)):
                        valid_indices = ~np.isnan(data_values)
                        if np.any(valid_indices):
                            data_values = np.interp(
                                np.arange(len(data_values)),
                                np.arange(len(data_values))[valid_indices],
                                data_values[valid_indices]
                            )
                        else:
                            continue
                    
                    # Now both distances are from dam (0) - no need to reverse!
                    # NetCDF distances and river distances are in the same reference frame
                    
                    # Interpolate directly
                    interp_func = interp1d(
                        nc_distances_from_dam,  # NetCDF distances from dam
                        data_values,            # Data values (no reversal needed)
                        kind='linear', 
                        bounds_error=False, 
                        fill_value='extrapolate'
                    )
                    interpolated_values = interp_func(river_distances[-1] - river_distances)
                    
                    # Verify interpolation
                    if np.any(np.isnan(interpolated_values)):
                        print(f"Warning: NaN values in interpolated data for {nc_var_name} at time {i}")
                    
                    var_time_series.append({
                        'timestamp': timestamp,
                        'values': interpolated_values.copy()
                    })
                    
                except Exception as e:
                    print(f"Error processing {nc_var_name} at time {i}: {e}")
                    continue
            
            if var_time_series:
                all_variables_data[var_name] = {
                    'time_series_data': var_time_series,
                    'units': data_var.attrs.get('units', ''),
                    'long_name': data_var.attrs.get('long_name', var_name),
                    'nc_var_name': nc_var_name
                }
                print(f"Successfully processed {var_name} with {len(var_time_series)} time steps")
                
                # Show sample statistics
                if len(var_time_series) > 0:
                    sample_values = var_time_series[0]['values']
                    print(f"  Sample values: min={np.min(sample_values):.3f}, max={np.max(sample_values):.3f}")
        
        # Cache results
        nc_data_cache = all_variables_data
        print(f"Successfully processed {len(all_variables_data)} variables")
        return nc_data_cache
            
    except Exception as e:
        print(f"Error processing NetCDF data: {e}")
        return None

def generate_guadalquivir_data(n_hours=24):
    """Generate river data using NetCDF interpolation"""
    global ds_global, nc_data_cache
    
    # Load river path and calculate distances
    key_points = load_river_path()
    river_distances = calculate_river_distances(key_points)
    river_distances = np.array(river_distances)
    
    # Load NetCDF data
    ds = load_netcdf_data()
    if ds is None:
        return None
    
    # Get interpolated data
    nc_data = interpolate_netcdf_to_river(ds, key_points, river_distances)
    if not nc_data:
        return None
    
    # Get timestamps
    first_var = list(nc_data.values())[0]
    timestamps = [item['timestamp'] for item in first_var['time_series_data']]
    
    # Define variables
    all_vars = ['A', 'Q', 'B', 'eta', 'rho', 'U', 'c', 'S', 'Qb', 'Qs']
    available_vars = list(nc_data.keys())
    
    all_data = []
    
    # Process each time step and point
    for time_idx, timestamp in enumerate(timestamps):
        for point_idx, (lat, lon) in enumerate(key_points):
            row_data = {
                "lat": lat,
                "lon": lon,
                "timestamp": timestamp,
                "point_id": point_idx,
                "distance_from_dam": river_distances[point_idx],  # Changed name to be clearer
            }
            
            max_distance = max(river_distances)
            distance_factor = river_distances[point_idx] / max_distance  # 0 at dam, 1 at mouth
            
            # Add data for each variable
            for var_name in all_vars:
                if var_name in nc_data:
                    if time_idx < len(nc_data[var_name]['time_series_data']):
                        value = nc_data[var_name]['time_series_data'][time_idx]['values'][point_idx]
                        
                        # Apply physical constraints
                        # if var_name == 'A':
                        #     value = max(1.0, float(value))
                        # elif var_name == 'Q':
                        #     value = max(0.1, float(value))
                        # elif var_name == 'B':
                        #     value = max(1.0, float(value))
                        # elif var_name == 'rho':
                        #     value = max(900, min(1100, float(value)))
                        # elif var_name == 'U':
                        #     value = max(0.0, float(value))
                        # elif var_name in ['c', 'S', 'Qb', 'Qs']:
                        #     value = max(0.0, float(value))
                        # else:
                        value = float(value)
                        
                        row_data[var_name] = value
                #     else:
                #         row_data[var_name] = get_default_value(var_name, distance_factor)
                # else:
                #     row_data[var_name] = get_default_value(var_name, distance_factor)
            
            all_data.append(row_data)
    
    if not all_data:
        return None
    
    df = pd.DataFrame(all_data)
    print(f"Generated DataFrame with {len(df)} rows")
    return df

def get_default_value(var_name, distance_factor):
    """Get default values for missing variables"""
    # distance_factor: 0 at dam, 1 at mouth
    defaults = {
        'A': 200.0 + 50.0 * (1 - distance_factor),      # Area decreases towards mouth
        'Q': 100.0,                                       # Constant discharge
        'B': 80.0 + 20.0 * (1 - distance_factor),       # Width decreases towards mouth  
        'eta': 0.0,                                       # Sea level reference
        'rho': 1000.0 + 25.0 * distance_factor ** 1.5,  # Density increases towards mouth (salinity)
        'U': 0.4 + 0.8 * distance_factor,               # Velocity higher towards mouth
        'c': 0.02 + 0.18 * (1 - distance_factor),       # Suspended solids higher upstream
        'S': 35.0 * distance_factor ** 1.5,             # Salinity: 0 at dam, 35 at mouth
        'Qb': 0.001 * (1 - distance_factor),            # Bedload transport higher upstream
        'Qs': 0.1 * (1 - distance_factor),              # Suspended transport higher upstream
    }
    return defaults.get(var_name, 0.0)

def load_river_path():
    """Load river path from GPKG file or use default coordinates"""
    gpkg_path = r"C:\Users\m_cob\OneDrive\Papers\2024\GRE I - Formulation\MDT - Guadalquivir\Thalweg_a_Puntos_500.gpkg"

    try:
        if os.path.exists(gpkg_path):
            gdf = gpd.read_file(gpkg_path)
            if gdf.crs != "EPSG:4326":
                gdf = gdf.to_crs("EPSG:4326")

            key_points = []
            for idx, row in gdf.iterrows():
                if row.geometry is not None:
                    if row.geometry.geom_type == "Point":
                        key_points.append((row.geometry.y, row.geometry.x))
                    elif row.geometry.geom_type in ["LineString", "MultiLineString"]:
                        if row.geometry.geom_type == "MultiLineString":
                            coords = list(row.geometry.geoms[0].coords)
                        else:
                            coords = list(row.geometry.coords)
                        
                        for coord in coords:
                            key_points.append((coord[1], coord[0]))

            if key_points:
                print(f"Loaded {len(key_points)} points from GPKG")
                return key_points

    except Exception as e:
        print(f"Error loading GPKG: {e}")

    # Default coordinates
    return [
        (37.5190, -5.9750),  # Alcalá del Río Dam
        (37.50855, -5.980669),
        (37.506221, -5.985770),
        (37.489709, -6.005108),
        (37.486285, -6.006035),
        (37.483134, -6.004705),
        (37.480735, -6.001652),
        (37.478783, -5.993207),
        (37.467414, -5.982668),
        (37.463008, -5.987293),
        (37.459478, -5.998231),
        (37.458570, -6.007990),
        (37.454328, -6.009197),
        (37.442254, -6.006071),
        (37.375215, -6.022755),
        (37.341555, -6.023062),
        (37.325292, -6.014520),
        (37.315000, -6.015712),
        (37.2200, -6.1200),
        (37.2000, -6.1320),
        (37.1800, -6.1450),
        (37.1600, -6.1580),
        (37.1400, -6.1700),
        (37.1200, -6.1820),
        (37.1000, -6.1950),
        (37.0800, -6.2080),
        (37.0600, -6.2200),
        (37.0400, -6.2320),
        (37.0200, -6.2450),
        (37.0000, -6.2580),
        (36.9800, -6.2700),
        (36.9600, -6.2820),
        (36.9400, -6.2940),
        (36.9200, -6.3060),
        (36.9000, -6.3180),
        (36.8800, -6.3290),
        (36.8600, -6.3380),
        (36.8400, -6.3450),
        (36.8200, -6.3500),
        (36.8000, -6.3520),
        (36.7785, -6.3535),  # Sanlúcar - Atlantic mouth
    ]

def cleanup_global_dataset():
    """Clean up global dataset and temp files"""
    global ds_global, nc_data_cache
    if ds_global is not None:
        try:
            temp_file_path = ds_global.attrs.get('_temp_file_path')
            ds_global.close()
            
            if temp_file_path and os.path.exists(temp_file_path):
                os.remove(temp_file_path)
                print(f"Cleaned up temporary file")
        except Exception as e:
            print(f"Error during cleanup: {e}")
        finally:
            ds_global = None
            nc_data_cache = None
            
atexit.register(cleanup_global_dataset)

# Parameter display names
PARAMETER_DISPLAY = {
    "A": "Cross-sectional Area (m²)",
    "Q": "Water Discharge (m³/s)",
    "B": "Channel Width (m)",
    "eta": "Water Level (m)",
    "rho": "Water Density (kg/m³)",
    "U": "Water Current (m/s)",
    "c": "Suspended Solid Concentration (kg/m³)",
    "S": "Salinity (PSU)",
    "Qb": "Bedload Sediment Transport (kg/s)",
    "Qs": "Suspended Sediment Transport (kg/s)",
}

# Monitoring points
MONITORING_POINTS = {
    "Alcalá del Río": {"target_lat": 37.5190, "point_id": None},
    "Seville": {"target_lat": 37.383, "point_id": None},
    "Gelves": {"target_lat": 37.342, "point_id": None},
    "Coria del Río": {"target_lat": 37.289, "point_id": None},
    "Puebla del Río": {"target_lat": 37.265, "point_id": None},
    "Isla Mayor": {"target_lat": 37.140, "point_id": None},
    "Bonanza": {"target_lat": 36.820, "point_id": None},
    "Sanlúcar": {"target_lat": 36.7785, "point_id": None},
}

def find_closest_thalweg_points():
    """Find closest thalweg point for each monitoring station"""
    key_points = load_river_path()
    
    for name, station in MONITORING_POINTS.items():
        target_lat = station["target_lat"]
        
        min_distance = float('inf')
        closest_point_id = 0
        
        for i, (lat, lon) in enumerate(key_points):
            lat_distance = abs(lat - target_lat)
            if lat_distance < min_distance:
                min_distance = lat_distance
                closest_point_id = i
        
        MONITORING_POINTS[name]["point_id"] = closest_point_id

find_closest_thalweg_points()

# Initialize Dash app
app = Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME])

# App layout (simplified)
app.layout = dbc.Container([
    # Header
    dbc.Row([
        dbc.Col([
            html.Div([
                html.H1([
                    html.I(className="fas fa-water me-3", style={"color": "#17a2b8"}),
                    "Guadalquivir River Estuary Monitoring System"
                ], className="text-center mb-2", style={"color": "#2c3e50", "fontWeight": "bold"}),
                html.P("Data visualization system", className="text-center text-muted mb-4"),
                html.Hr(style={"borderColor": "#17a2b8", "borderWidth": "2px"}),
            ], className="bg-white rounded p-4 shadow-sm mb-4")
        ])
    ]),
    
    # Main content
    dbc.Row([
        dbc.Col([
            dbc.Card([
                dbc.CardHeader([
                    html.H5([
                        html.I(className="fas fa-map me-2"),
                        "Estuarine Main Channel Map"
                    ], className="mb-0", style={"color": "#495057"})
                ]),
                dbc.CardBody([
                    html.Div([
                        # Control panel
                        html.Div(
                            dbc.Card([
                                dbc.CardHeader([
                                    html.H6([
                                        html.I(className="fas fa-cogs me-2"),
                                        "Control Panel"
                                    ], className="mb-0", style={"fontSize": "14px"})
                                ]),
                                dbc.CardBody([
                                    dbc.Label("Parameter:", style={"fontWeight": "bold", "fontSize": "12px"}),
                                    dcc.Dropdown(
                                        id="parameter-dropdown",
                                        options=[
                                            {"label": "🌊 Q: Water Discharge (m³/s)", "value": "Q"},
                                            {"label": "➡️ U: Water Current (m/s)", "value": "U"},
                                            {"label": "🧂 S: Salinity (PSU)", "value": "S"},
                                            {"label": "⚖️ ρ: Water Density (kg/m³)", "value": "rho"},
                                            {"label": "🔲 A: Cross-sectional Area (m²)", "value": "A"},
                                            {"label": "📏 B: Channel Width (m)", "value": "B"},
                                            {"label": "📊 η: Water Level (m)", "value": "eta"},
                                            {"label": "🟤 c: Suspended Solids (kg/m³)", "value": "c"},
                                            {"label": "🪨 Qb: Bedload Transport (kg/s)", "value": "Qb"},
                                            {"label": "💧 Qs: Suspended Transport (kg/s)", "value": "Qs"},
                                        ],
                                        value="U",
                                        style={"marginBottom": "15px", "fontSize": "12px"}
                                    ),
                                    dbc.Label("Time:", style={"fontWeight": "bold", "fontSize": "12px"}),
                                    dcc.Slider(
                                        id="time-slider",
                                        min=0, max=47, value=47, step=1,
                                        marks={},
                                        tooltip={"placement": "bottom", "always_visible": True},
                                        className="mb-3"
                                    ),
                                    html.Div(id="selected-time-display", style={
                                        "textAlign": "center", "fontSize": "11px",
                                        "fontWeight": "bold", "color": "#17a2b8"
                                    })
                                ])
                            ], className="shadow", style={
                                "backgroundColor": "rgba(255, 255, 255, 0.95)",
                                "width": "280px"
                            }),
                            style={
                                "position": "absolute", "top": "20px", "left": "20px",
                                "zIndex": 1000, "pointerEvents": "auto"
                            }
                        ),
                        # Map
                        dl.Map(
                            id="river-map",
                            style={"width": "100%", "height": "850px", "borderRadius": "8px"},
                            center=[37.15, -6.1], zoom=10, zoomControl=False,
                            children=[
                                dl.TileLayer(url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"),
                                dl.LayerGroup(id="data-layer"),
                                dl.ZoomControl(position="bottomright")
                            ]
                        )
                    ], style={"position": "relative", "height": "850px"})
                ])
            ], className="shadow-sm mb-4")
        ], width=7),
        
        dbc.Col([
            # Histogram
            dbc.Card([
                dbc.CardHeader([
                    html.H5([html.I(className="fas fa-chart-bar me-2"), "Distribution"], className="mb-0")
                ]),
                dbc.CardBody([dcc.Graph(id="parameter-histogram", style={"height": "320px"})])
            ], className="shadow-sm mb-4"),
            
            # Multi-point comparison
            dbc.Card([
                dbc.CardHeader([
                    html.H5([html.I(className="fas fa-layer-group me-2"), "Multi-Point Comparison"], className="mb-0")
                ]),
                dbc.CardBody([dcc.Graph(id="multipoint-comparison", style={"height": "420px"})])
            ], className="shadow-sm")
        ], width=5)
    ])
], fluid=True, style={"backgroundColor": "#f8f9fa", "minHeight": "100vh", "padding": "20px"})

# Callbacks
@callback(
    [Output("time-slider", "marks"), Output("time-slider", "max")],
    Input("parameter-dropdown", "value")
)
def update_time_slider(selected_param):
    global df_global
    
    if df_global is None:
        df_global = generate_guadalquivir_data(n_hours=48)
        
        if df_global is None:
            return {"0": {"label": "No Data"}}, 0
    
    timestamps = sorted(df_global["timestamp"].unique())
    
    # Only show start, middle, and end marks
    marks = {
        0: {"label": "Start", "style": {"fontSize": "10px"}},
        len(timestamps) // 2: {"label": "Middle", "style": {"fontSize": "10px"}},
        len(timestamps) - 1: {"label": "End", "style": {"fontSize": "10px"}}
    }
    
    return marks, len(timestamps) - 1

@callback(
    Output("selected-time-display", "children"),
    Input("time-slider", "value")
)
def update_time_display(time_index):
    global df_global
    
    if df_global is None or time_index is None:
        return "No data available"
    
    timestamps = sorted(df_global["timestamp"].unique())
    if time_index < len(timestamps):
        selected_time = timestamps[time_index]
        total_hours = len(timestamps)
        progress = f"({time_index + 1}/{total_hours})"
        
        return html.Div([
            html.Div(f"Time: {selected_time.strftime('%Y-%m-%d %H:%M:%S')}", 
                    style={"fontWeight": "bold", "color": "#17a2b8"}),
            html.Div(f"Step {progress}", 
                    style={"fontSize": "10px", "color": "#6c757d"})
        ])
    return "Invalid selection"

@callback(
    [Output("data-layer", "children"), Output("parameter-histogram", "figure"), Output("multipoint-comparison", "figure")],
    [Input("parameter-dropdown", "value"), Input("time-slider", "value")]
)
def update_visualizations(selected_param, time_index):
    global df_global
    
    if df_global is None:
        df_global = generate_guadalquivir_data(n_hours=48)
        
        if df_global is None:
            empty_fig = go.Figure()
            empty_fig.add_annotation(text="ERROR: Data could not be loaded", 
                                   xref="paper", yref="paper", x=0.5, y=0.5, 
                                   showarrow=False, font=dict(size=16, color="red"))
            return [], empty_fig, empty_fig
    
    # Get selected timestamp
    timestamps = sorted(df_global["timestamp"].unique())
    if time_index is not None and time_index < len(timestamps):
        selected_timestamp = timestamps[time_index]
    else:
        selected_timestamp = timestamps[-1]
    
    df_selected = df_global[df_global["timestamp"] == selected_timestamp]

    # Create map elements
    river_coords = [[row["lat"], row["lon"]] for _, row in df_selected.sort_values("point_id").iterrows()]
    
    markers = [dl.Polyline(positions=river_coords, color="#17a2b8", weight=3, opacity=0.8)]

    # Data points with color coding
    param_min, param_max = df_selected[selected_param].min(), df_selected[selected_param].max()
    if param_max > param_min:
        for _, row in df_selected.iterrows():
            normalized_value = (row[selected_param] - param_min) / (param_max - param_min)
            color = f"hsl({240 - normalized_value * 120}, 70%, 50%)"

            markers.append(dl.CircleMarker(
                center=[row["lat"], row["lon"]],
                radius=4 + normalized_value * 8,
                color=color, fill=True, fillOpacity=0.8, weight=2,
                children=[dl.Tooltip(f"{PARAMETER_DISPLAY[selected_param]}: {row[selected_param]:.2f}")]
            ))

    # Monitoring points
    colors = ["#6f42c1", "#dc3545", "#28a745", "#fd7e14", "#17a2b8", "#e83e8c", "#20c997", "#ffc107"]
    
    for i, (name, point) in enumerate(MONITORING_POINTS.items()):
        thalweg_point = df_selected[df_selected["point_id"] == point["point_id"]]
        
        if not thalweg_point.empty:
            lat, lon = thalweg_point.iloc[0]["lat"], thalweg_point.iloc[0]["lon"]
            distance = thalweg_point.iloc[0]["distance_from_dam"]  # Updated column name
            param_value = thalweg_point.iloc[0][selected_param]
            
            icon_color = colors[i % len(colors)]
            
            # Monitoring station marker
            markers.append(dl.Marker(
                position=[lat, lon],
                children=[dl.Tooltip(html.Div([
                    html.H6(f"📍 {name}", style={"margin": "0", "color": "#2c3e50"}),
                    html.P(f"Distance from mouth: {distance/1000:.1f} km", style={"margin": "2px 0"}),  # Updated label
                    html.P(f"{PARAMETER_DISPLAY[selected_param]}: {param_value:.2f}", style={"margin": "2px 0"})
                ]))],
                icon={
                    "iconUrl": "data:image/svg+xml;base64," + __import__('base64').b64encode(
                        f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" width="24" height="24">
                            <circle cx="12" cy="12" r="10" fill="rgba(255,255,255,0.95)" 
                                    stroke="{icon_color}" stroke-width="3"/>
                            <circle cx="12" cy="12" r="6" fill="{icon_color}"/>
                        </svg>'''.encode()
                    ).decode(),
                    "iconSize": [24, 24], "iconAnchor": [12, 12]
                }
            ))
            
            # Station label
            markers.append(dl.Marker(
                position=[lat, lon],
                icon={
                    "iconUrl": "data:image/svg+xml;base64," + __import__('base64').b64encode(
                        f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 100 20" width="100" height="20">
                            <rect x="0" y="0" width="100" height="20" fill="rgba(255,255,255,0.9)" 
                                  stroke="{icon_color}" stroke-width="1" rx="10"/>
                            <text x="50" y="14" text-anchor="middle" font-family="Arial" 
                                  font-size="11" fill="{icon_color}" font-weight="bold">{name}</text>
                        </svg>'''.encode()
                    ).decode(),
                    "iconSize": [100, 20], "iconAnchor": [110, 10]
                }
            ))

    # Charts
    param_display = PARAMETER_DISPLAY[selected_param]

    # Histogram
    hist_fig = px.histogram(df_selected, x=selected_param, nbins=20, template="plotly_white")
    hist_fig.update_layout(
        title=f"{selected_timestamp.strftime('%H:%M:%S')}",
        xaxis_title=param_display, yaxis_title="Frequency",
        height=280, showlegend=False, margin=dict(l=40, r=40, t=40, b=40)
    )
    hist_fig.update_traces(marker_color="#17a2b8", opacity=0.7)

    # Multi-point comparison
    multipoint_fig = go.Figure()
    
    for i, (name, point) in enumerate(MONITORING_POINTS.items()):
        df_point = df_global[df_global["point_id"] == point["point_id"]].sort_values("timestamp")
        multipoint_fig.add_trace(go.Scatter(
            x=df_point["timestamp"], y=df_point[selected_param],
            mode="lines", name=name,
            line=dict(color=colors[i % len(colors)], width=2)
        ))

    # Selected time line
    y_range = [df_global[selected_param].min(), df_global[selected_param].max()]
    multipoint_fig.add_shape(
        type="line", x0=selected_timestamp, x1=selected_timestamp,
        y0=y_range[0], y1=y_range[1],
        line=dict(color="red", width=2, dash="dash")
    )

    multipoint_fig.update_layout(
        xaxis_title="Time", yaxis_title=param_display,
        height=330, template="plotly_white",
        margin=dict(l=40, r=40, t=40, b=40)
    )

    return markers, hist_fig, multipoint_fig

if __name__ == "__main__":
    app.run(debug=True)
