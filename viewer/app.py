import os
from datetime import datetime, timedelta

import dash_bootstrap_components as dbc
import dash_leaflet as dl
import geopandas as gpd
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from dash import Dash, Input, Output, callback, dcc, html


def load_river_path():
    """Load river path from GPKG file or use default coordinates"""
    gpkg_path = r"C:\Users\m_cob\OneDrive\Papers\2024\GRE I - Formulation\MDT - Guadalquivir\Thalweg_a_Puntos_500.gpkg"

    try:
        # Try to load the GPKG file
        if os.path.exists(gpkg_path):
            gdf = gpd.read_file(gpkg_path)
            print(f"Loaded GPKG with {len(gdf)} features")
            print(f"Original CRS: {gdf.crs}")
            print(f"Geometry types: {gdf.geometry.geom_type.unique()}")
            print(f"Columns: {list(gdf.columns)}")

            # Transform to WGS84 if not already in WGS84
            if gdf.crs != "EPSG:4326":
                print(f"Transforming from {gdf.crs} to WGS84 (EPSG:4326)")
                gdf = gdf.to_crs("EPSG:4326")
                print("Transformation completed")
            else:
                print("Data already in WGS84")

            # Extract coordinates from geometry
            key_points = []
            for idx, row in gdf.iterrows():
                if row.geometry is not None and not row.geometry.is_empty:
                    if row.geometry.geom_type == "Point":
                        # For Point geometry
                        lat, lon = row.geometry.y, row.geometry.x
                        key_points.append((lat, lon))
                    elif row.geometry.geom_type in ["LineString", "MultiLineString"]:
                        # For LineString geometry
                        if row.geometry.geom_type == "MultiLineString":
                            # Take the first linestring if it's a MultiLineString
                            coords = list(row.geometry.geoms[0].coords)
                        else:
                            coords = list(row.geometry.coords)

                        for coord in coords:
                            lon, lat = (
                                coord[0],
                                coord[1],
                            )  # Now in WGS84 (lon, lat)
                            key_points.append((lat, lon))

            if key_points:
                print(f"Successfully loaded {len(key_points)} points from {gpkg_path}")
                print(f"First few points (WGS84): {key_points[:3]}")
                print(f"Last few points (WGS84): {key_points[-3:]}")

                # Validate coordinates are reasonable for Guadalquivir area
                lats = [point[0] for point in key_points]
                lons = [point[1] for point in key_points]
                lat_range = (min(lats), max(lats))
                lon_range = (min(lons), max(lons))

                print(f"Latitude range: {lat_range}")
                print(f"Longitude range: {lon_range}")

                # Check if coordinates are in expected range for Guadalquivir
                if (
                    36.5 <= lat_range[0] <= 38.0
                    and 36.5 <= lat_range[1] <= 38.0
                    and -7.0 <= lon_range[0] <= -5.5
                    and -7.0 <= lon_range[1] <= -5.5
                ):
                    print("✓ Coordinates appear to be correctly transformed to WGS84")
                    return key_points
                else:
                    print("⚠ Warning: Coordinates may not be correctly transformed")
                    print("Expected ranges: Lat 36.5-38.0, Lon -7.0 to -5.5")
                    return key_points
            else:
                print("No valid coordinates found in GPKG file")

    except Exception as e:
        print(f"Error loading GPKG file: {e}")
        print("Using default coordinates...")

    # Default coordinates if file loading fails (already in WGS84)
    print("Using default river path coordinates")
    return [
        (37.5190, -5.9750),  # Alcalá del Río Dam
        (37.50855, -5.980669),  # Downstream from dam
        (37.506221, -5.985770),  # Approaching Seville
        (37.489709, -6.005108),  # Northern Seville suburbs
        (37.486285, -6.006035),  # Northern Seville area
        (37.483134, -6.004705),  # Seville city center (Guadalquivir bend)
        (37.480735, -6.001652),  # Seville port area
        (37.478783, -5.993207),  # Seville exit towards Gelves
        (37.467414, -5.982668),  # Approaching Gelves
        (37.463008, -5.987293),  # Gelves
        (37.459478, -5.998231),  # Between Gelves and Palomares
        (37.458570, -6.007990),  # Palomares del Río
        (37.454328, -6.009197),  # Downstream Palomares
        (37.442254, -6.006071),  # Coria del Río
        (37.375215, -6.022755),  # Between Coria and Puebla
        (37.341555, -6.023062),  # Puebla del Río
        (37.325292, -6.014520),  # Downstream Puebla
        (37.315000, -6.015712),  # Approaching marshlands
        (37.2200, -6.1200),  # Entry to marshlands
        (37.2000, -6.1320),  # Northern marshlands
        (37.1800, -6.1450),  # Deep into marshlands
        (37.1600, -6.1580),  # Central marshlands
        (37.1400, -6.1700),  # Isla Mayor area
        (37.1200, -6.1820),  # Continuing through marshlands
        (37.1000, -6.1950),  # Central Guadalquivir marshlands
        (37.0800, -6.2080),  # Deeper marshlands
        (37.0600, -6.2200),  # Approaching Trebujena influence
        (37.0400, -6.2320),  # Central marshlands zone
        (37.0200, -6.2450),  # Continuing towards coast
        (37.0000, -6.2580),  # Central marshlands
        (36.9800, -6.2700),  # Approaching Trebujena area
        (36.9600, -6.2820),  # Trebujena vicinity
        (36.9400, -6.2940),  # Between Trebujena and Sanlúcar
        (36.9200, -6.3060),  # Approaching Sanlúcar
        (36.9000, -6.3180),  # Towards Sanlúcar channel
        (36.8800, -6.3290),  # Sanlúcar approach channel
        (36.8600, -6.3380),  # Bonanza area
        (36.8400, -6.3450),  # Port approach
        (36.8200, -6.3500),  # Bonanza/Port area
        (36.8000, -6.3520),  # Final approach
        (36.7785, -6.3535),  # Sanlúcar de Barrameda - Atlantic mouth
    ]


# Generate random Guadalquivir river data with temporal evolution
def generate_guadalquivir_data(n_hours=24):
    # Load river path coordinates
    np.random.seed(42)

    # Load real path coordinates from GPKG file or use defaults
    key_points = load_river_path()
    n_points = len(key_points)
    print(f"Using {len(key_points)} key points for river path generation")

    # Generate points interpolating between key points with smooth curves
    river_path = []
    for i in range(n_points):
        t = i / (n_points - 1)  # Parameter from 0 to 1

        # Find corresponding segment
        segment_index = min(int(t * (len(key_points) - 1)), len(key_points) - 2)
        local_t = (t * (len(key_points) - 1)) - segment_index

        # Linear interpolation between key points
        p1 = key_points[segment_index]
        p2 = key_points[segment_index + 1]

        lat = p1[0] + local_t * (p2[0] - p1[0])
        lon = p1[1] + local_t * (p2[1] - p1[1])

        # Add very subtle natural river meandering (much reduced to stay in channel)
        dam_factor = 1 - np.exp(-t * 12)  # Less variation near dam
        meander_amplitude = 0.0015  # Very small meanders

        # Different meander patterns for different river sections
        if t < 0.3:  # Upper river (dam to Seville) - more controlled
            lat += meander_amplitude * 0.5 * np.sin(t * 20 * np.pi) * dam_factor
            lon += meander_amplitude * 0.3 * np.cos(t * 18 * np.pi) * dam_factor
        elif t < 0.7:  # Middle river (Seville to marshlands) - moderate meandering
            lat += meander_amplitude * np.sin(t * 25 * np.pi) * np.exp(-t * 1.5)
            lon += meander_amplitude * 0.8 * np.cos(t * 22 * np.pi) * (1 - t * 0.5)
        else:  # Lower river (marshlands to sea) - wider, more stable channel
            lat += meander_amplitude * 0.7 * np.sin(t * 15 * np.pi) * (1 - t * 0.8)
            lon += meander_amplitude * 0.5 * np.cos(t * 12 * np.pi) * (1 - t * 0.9)

        river_path.append((lat, lon, i))

    lats = np.array([point[0] for point in river_path])
    lons = np.array([point[1] for point in river_path])
    point_distances = np.array([point[2] for point in river_path])

    # Generate time series for the last n_hours
    timestamps = [datetime.now() - timedelta(hours=i) for i in range(n_hours, 0, -1)]

    all_data = []
    for i, ts in enumerate(timestamps):
        # Add temporal variability with trends
        time_factor = i / len(timestamps)

        # Current data with temporal and spatial variation along the river
        # Velocity pattern: higher near dam (controlled release), varies downstream
        dam_influence = np.exp(
            -point_distances / 15
        )  # Dam influence decreases downstream
        base_velocities = (
            1.2 * dam_influence  # High velocity from dam release
            + 0.5
            + 1.5 * (point_distances / (n_points - 1))  # Increasing towards sea
        )
        velocities = base_velocities * (
            1 + 0.3 * np.sin(time_factor * 2 * np.pi)
        ) + np.random.normal(0, 0.1, n_points)

        # Directions: more controlled near dam, more variable downstream
        controlled_flow = 200 * dam_influence  # Controlled direction from dam
        natural_flow = 225 + 30 * (
            point_distances / (n_points - 1)
        )  # Natural downstream flow
        base_directions = controlled_flow + natural_flow * (1 - dam_influence)
        directions = (
            base_directions
            + 20
            * np.sin(time_factor * 4 * np.pi)
            * (1 - dam_influence * 0.8)  # Less variation near dam
            + np.random.normal(0, 10, n_points) * (1 - dam_influence * 0.5)
        )

        # Temperature: influenced by dam release (cooler water from depth)
        dam_cooling = -2 * dam_influence  # Cooler water from dam
        base_temps = 22 - 3 * (point_distances / (n_points - 1)) + dam_cooling
        temperatures = (
            base_temps
            + 3 * np.sin(time_factor * 2 * np.pi)
            + np.random.normal(0, 0.5, n_points)
        )

        # Salinity: very low near dam, increases towards sea
        base_salinity = (
            0.5 + 32 * (point_distances / (n_points - 1)) ** 1.5
        )  # Non-linear increase
        salinities = base_salinity * (
            1 + 0.2 * np.cos(time_factor * 2 * np.pi)
        ) + np.random.normal(0, 1, n_points)

        for j in range(n_points):
            all_data.append(
                {
                    "lat": lats[j],
                    "lon": lons[j],
                    "velocity": max(0.1, velocities[j]),
                    "direction": directions[j] % 360,
                    "temperature": temperatures[j],
                    "salinity": max(0, salinities[j]),
                    "timestamp": ts,
                    "point_id": j,
                    "distance_from_source": point_distances[j],
                }
            )

    return pd.DataFrame(all_data)


# Parameter units dictionary
PARAMETER_UNITS = {
    "velocity": "m/s",
    "temperature": "°C",
    "salinity": "PSU",
    "direction": "°",
}

# Parameter display names with units
PARAMETER_DISPLAY = {
    "velocity": "Current Velocity (m/s)",
    "temperature": "Temperature (°C)",
    "salinity": "Salinity (PSU)",
    "direction": "Direction (°)",
}


# Updated monitoring points to match the actual river channel better
MONITORING_POINTS = {
    "Alcalá del Río": {"lat": 37.5190, "lon": -5.9750, "point_id": 0},
    "Seville": {"lat": 37.3891, "lon": -5.9845, "point_id": 12},
    "Gelves": {"lat": 37.3420, "lon": -6.0280, "point_id": 18},
    "Coria del Río": {"lat": 37.2890, "lon": -6.0620, "point_id": 25},
    "Puebla del Río": {"lat": 37.2650, "lon": -6.0820, "point_id": 30},
    "Isla Mayor": {"lat": 37.1400, "lon": -6.1700, "point_id": 42},
    "Bonanza": {"lat": 36.8200, "lon": -6.3500, "point_id": 55},
    "Sanlúcar": {"lat": 36.7785, "lon": -6.3535, "point_id": 59},
}

# Initialize app with custom theme
app = Dash(
    __name__, external_stylesheets=[dbc.themes.BOOTSTRAP, dbc.icons.FONT_AWESOME]
)

# Custom CSS styling
custom_style = {
    "background": "linear-gradient(135deg, #667eea 0%, #764ba2 100%)",
    "minHeight": "100vh",
    "fontFamily": "Arial, sans-serif",
}

# Main layout
app.layout = dbc.Container(
    [
        # Header section
        dbc.Row(
            [
                dbc.Col(
                    [
                        html.Div(
                            [
                                html.H1(
                                    [
                                        html.I(
                                            className="fas fa-water me-3",
                                            style={"color": "#17a2b8"},
                                        ),
                                        "Guadalquivir River Estuary Hydrodynamic Monitoring System",
                                    ],
                                    className="text-center mb-2",
                                    style={"color": "#2c3e50", "fontWeight": "bold"},
                                ),
                                html.P(
                                    "Real-time hydrodynamic data visualization system",
                                    className="text-center text-muted mb-4",
                                    style={"fontSize": "1.1em"},
                                ),
                                html.Hr(
                                    style={
                                        "borderColor": "#17a2b8",
                                        "borderWidth": "2px",
                                    }
                                ),
                            ],
                            className="bg-white rounded p-4 shadow-sm mb-4",
                        )
                    ]
                )
            ]
        ),
        # Main content - reorganized layout
        dbc.Row(
            [
                # Left column - River map with integrated control panel
                dbc.Col(
                    [
                        # River Map with embedded control panel
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    [
                                        html.H5(
                                            [
                                                html.I(className="fas fa-map me-2"),
                                                "River Map - Alcalá Dam to Atlantic",
                                            ],
                                            className="mb-0",
                                            style={"color": "#495057"},
                                        )
                                    ]
                                ),
                                dbc.CardBody(
                                    [
                                        html.Div(
                                            [
                                                # Control panel floating over the map
                                                html.Div(
                                                    dbc.Card(
                                                        [
                                                            dbc.CardHeader(
                                                                [
                                                                    html.H6(
                                                                        [
                                                                            html.I(className="fas fa-cogs me-2"),
                                                                            "Control Panel",
                                                                        ],
                                                                        className="mb-0",
                                                                        style={
                                                                            "color": "#495057",
                                                                            "fontSize": "14px"
                                                                        },
                                                                    )
                                                                ],
                                                                className="py-2",
                                                            ),
                                                            dbc.CardBody(
                                                                [
                                                                    html.Div(
                                                                        [
                                                                            dbc.Label(
                                                                                [
                                                                                    html.I(
                                                                                        className="fas fa-chart-line me-1"
                                                                                    ),
                                                                                    "Parameter:",
                                                                                ],
                                                                                style={
                                                                                    "fontWeight": "bold",
                                                                                    "color": "#495057",
                                                                                    "fontSize": "12px",
                                                                                },
                                                                            ),
                                                                            dcc.Dropdown(
                                                                                id="parameter-dropdown",
                                                                                options=[
                                                                                    {
                                                                                        "label": "🌊 Velocity",
                                                                                        "value": "velocity",
                                                                                    },
                                                                                    {
                                                                                        "label": "🌡️ Temperature",
                                                                                        "value": "temperature",
                                                                                    },
                                                                                    {
                                                                                        "label": "🧂 Salinity",
                                                                                        "value": "salinity",
                                                                                    },
                                                                                ],
                                                                                value="velocity",
                                                                                style={
                                                                                    "marginBottom": "10px",
                                                                                    "fontSize": "12px",
                                                                                },
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    html.Div(
                                                                        [
                                                                            dbc.Label(
                                                                                [
                                                                                    html.I(
                                                                                        className="fas fa-map-marker-alt me-1"
                                                                                    ),
                                                                                    "Point:",
                                                                                ],
                                                                                style={
                                                                                    "fontWeight": "bold",
                                                                                    "color": "#495057",
                                                                                    "fontSize": "12px",
                                                                                },
                                                                            ),
                                                                            dcc.Dropdown(
                                                                                id="monitoring-point-dropdown",
                                                                                options=[
                                                                                    {"label": name, "value": name}
                                                                                    for name in MONITORING_POINTS.keys()
                                                                                ],
                                                                                value="Seville",
                                                                                style={
                                                                                    "marginBottom": "10px",
                                                                                    "fontSize": "12px",
                                                                                },
                                                                            ),
                                                                        ]
                                                                    ),
                                                                    dbc.Button(
                                                                        [
                                                                            html.I(
                                                                                className="fas fa-sync-alt me-1"
                                                                            ),
                                                                            "Refresh",
                                                                        ],
                                                                        id="update-btn",
                                                                        color="info",
                                                                        className="w-100",
                                                                        size="sm",
                                                                        style={"fontSize": "12px"},
                                                                    ),
                                                                ],
                                                                className="py-2",
                                                            ),
                                                        ],
                                                        className="shadow",
                                                        style={
                                                            "border": "none",
                                                            "backgroundColor": "rgba(255, 255, 255, 0.95)",
                                                            "backdropFilter": "blur(5px)",
                                                            "width": "280px",
                                                        },
                                                    ),
                                                    style={
                                                        "position": "absolute",
                                                        "top": "20px",
                                                        "left": "20px",
                                                        "zIndex": 1000,
                                                        "pointerEvents": "auto",
                                                    },
                                                ),
                                                # Map component
                                                dl.Map(
                                                    id="river-map",
                                                    style={
                                                        "width": "100%",
                                                        "height": "650px",
                                                        "borderRadius": "8px",
                                                    },
                                                    center=[37.15, -6.1],
                                                    zoom=9,
                                                    zoomControl=False,  # Deshabilita el control de zoom por defecto
                                                    children=[
                                                        dl.TileLayer(
                                                            url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png"
                                                        ),
                                                        dl.LayerGroup(id="data-layer"),
                                                        # Agregar control de zoom personalizado
                                                        dl.ZoomControl(position="bottomright"),  # O la posición que prefieras
                                                    ],
                                                ),
                                            ],
                                            style={
                                                "position": "relative",
                                                "width": "100%",
                                                "height": "650px",
                                            },
                                        )
                                    ],
                                    className="p-2",
                                ),
                            ],
                            className="shadow-sm",
                            style={"border": "none"},
                        )
                    ],
                    width=8,
                ),
                # Right column - Distribution and Multi-Point Comparison (sin cambios)
                dbc.Col(
                    [
                        # Distribution chart
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    [
                                        html.H5(
                                            [
                                                html.I(
                                                    className="fas fa-chart-bar me-2"
                                                ),
                                                "Distribution",
                                            ],
                                            className="mb-0",
                                            style={"color": "#495057"},
                                        )
                                    ]
                                ),
                                dbc.CardBody(
                                    [
                                        dcc.Graph(
                                            id="parameter-histogram",
                                            style={"height": "300px"},
                                        )
                                    ],
                                    className="p-2",
                                ),
                            ],
                            className="shadow-sm mb-4",
                            style={"border": "none"},
                        ),
                        # Multi-Point Comparison chart
                        dbc.Card(
                            [
                                dbc.CardHeader(
                                    [
                                        html.H5(
                                            [
                                                html.I(
                                                    className="fas fa-layer-group me-2"
                                                ),
                                                "Multi-Point Comparison",
                                            ],
                                            className="mb-0",
                                            style={"color": "#495057"},
                                        )
                                    ]
                                ),
                                dbc.CardBody(
                                    [
                                        dcc.Graph(
                                            id="multipoint-comparison",
                                            style={"height": "350px"},
                                        )
                                    ],
                                    className="p-2",
                                ),
                            ],
                            className="shadow-sm",
                            style={"border": "none"},
                        )
                    ],
                    width=4,
                ),
            ]
        ),
    ],
    fluid=True,
    style={"backgroundColor": "#f8f9fa", "minHeight": "100vh", "padding": "20px"},
)


@callback(
    [
        Output("data-layer", "children"),
        Output("parameter-histogram", "figure"),
        Output("multipoint-comparison", "figure"),
    ],
    [
        Input("update-btn", "n_clicks"),
        Input("parameter-dropdown", "value"),
        Input("monitoring-point-dropdown", "value"),
    ],
)
def update_visualizations(n_clicks, selected_param, selected_point):
    # Generate new data with time series
    df_full = generate_guadalquivir_data(n_hours=48)
    df_latest = df_full[df_full["timestamp"] == df_full["timestamp"].max()]

    # Create river path line
    river_coords = [
        [row["lat"], row["lon"]]
        for _, row in df_latest.sort_values("point_id").iterrows()
    ]

    # Create map markers
    markers = [
        # River path
        dl.Polyline(positions=river_coords, color="#17a2b8", weight=3, opacity=0.8)
    ]

    # Data points with color coding
    for _, row in df_latest.iterrows():
        normalized_value = (row[selected_param] - df_latest[selected_param].min()) / (
            df_latest[selected_param].max() - df_latest[selected_param].min()
        )
        color = f"hsl({240 - normalized_value * 120}, 70%, 50%)"  # Blue to red gradient

        markers.append(
            dl.CircleMarker(
                center=[row["lat"], row["lon"]],
                radius=4 + normalized_value * 8,
                color=color,
                fill=True,
                fillOpacity=0.8,
                weight=2,
                children=[
                    dl.Tooltip(
                        f"{PARAMETER_DISPLAY[selected_param]}: {row[selected_param]:.2f}"
                    )
                ],
            )
        )

    # Monitoring points with icons
    for name, point in MONITORING_POINTS.items():
        # Different icon based on location type
        if "Dam" in name or "Alcalá" in name:
            icon_name = "tint"  # Dam icon
            icon_color = "#6f42c1"  # Purple for dam
        elif "Seville" in name:
            icon_name = "city"  # City icon
            icon_color = "#dc3545"  # Red for major city
        elif any(x in name for x in ["Bonanza", "Sanlúcar"]):
            icon_name = "anchor"  # Port/coastal icon
            icon_color = "#17a2b8"  # Blue for coastal areas
        else:
            icon_name = "map-marker"  # Default marker
            icon_color = "#28a745"  # Green for other locations

        # Highlight selected point
        if name == selected_point:
            icon_color = "#ffc107"  # Yellow/gold for selected point
            icon_size = 25
        else:
            icon_size = 20

        markers.append(
            dl.Marker(
                position=[point["lat"], point["lon"]],
                children=[
                    dl.Tooltip(
                        f"📊 {name} (Monitoring Station)",
                        direction="top",
                        permanent=False,
                    )
                ],
                icon={
                    "iconUrl": f"https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/svgs/solid/{icon_name}.svg",
                    "iconSize": [icon_size, icon_size],
                    "iconAnchor": [icon_size // 2, icon_size],
                    "popupAnchor": [0, -icon_size],
                    "className": "monitoring-point-icon",
                },
            )
        )

    # Charts with improved styling and units
    chart_template = "plotly_white"
    param_unit = PARAMETER_UNITS[selected_param]
    param_display = PARAMETER_DISPLAY[selected_param]

    # Histogram
    hist_fig = px.histogram(
        df_latest, x=selected_param, nbins=20, template=chart_template
    )
    hist_fig.update_layout(
        title=f"Distribution of {param_display}",
        xaxis_title=param_display,
        yaxis_title="Frequency",
        height=280,
        showlegend=False,
        margin=dict(l=40, r=40, t=40, b=40),
    )
    hist_fig.update_traces(marker_color="#17a2b8", opacity=0.7)

    # Multi-point comparison
    multipoint_fig = go.Figure()
    colors = [
        "#17a2b8",
        "#dc3545",
        "#28a745",
        "#fd7e14",
        "#6f42c1",
        "#e83e8c",
        "#20c997",
        "#ffc107",
    ]

    for i, (name, point) in enumerate(MONITORING_POINTS.items()):
        df_point = df_full[df_full["point_id"] == point["point_id"]].sort_values(
            "timestamp"
        )
        multipoint_fig.add_trace(
            go.Scatter(
                x=df_point["timestamp"],
                y=df_point[selected_param],
                mode="lines",
                name=name,
                line=dict(color=colors[i % len(colors)], width=2),
            )
        )

    multipoint_fig.update_layout(
        title=f"{param_display} Comparison",
        xaxis_title="Time",
        yaxis_title=param_display,
        height=330,
        template=chart_template,
        margin=dict(l=40, r=40, t=40, b=40),
    )

    return markers, hist_fig, multipoint_fig


if __name__ == "__main__":
    app.run(debug=True)
