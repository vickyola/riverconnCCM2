# River Connectivity Analysis Project

This project focuses on analyzing river connectivity using geospatial data and network analysis techniques. The workflow involves data selection, preprocessing, spatial transformations, and connectivity analysis. Below is a structured guide to understand and execute the provided Python code for this project.

## Overview

The project involves the following key steps:

1. **Data Selection**: Choose the river basin of interest.
2. **CRS Matching**: Ensure spatial data are in the same coordinate reference system.
3. **Data Preprocessing**: Prepare and subset data for analysis.
4. **Geospatial Transformations**: Process and transform spatial data.
5. **Network Construction**: Build and refine the river network graph.
6. **Index Calculations**: Compute various connectivity indices.
7. **Visualization**: Plot the results to interpret connectivity metrics.

## Prerequisites

- **Python Libraries**: Ensure you have the required Python libraries installed, including `pandas`, `geopandas`, `networkx`, and `matplotlib`. Use the following command to install any missing packages:

  ```bash
  pip install pandas geopandas networkx matplotlib
Data Files

Place your geospatial datasets in the project directory. The required files include:

    Basin data (basin_file.csv)
    Dam data (dam_file.shp)
    River data (river_file.shp)

## Usage
1. Data Selection

Choose a river basin using its name or unique identifier:

```
import pandas as pd

# Load basin data
basin_df = pd.read_csv('path_to_basin_file.csv')

# Select river basin
basin_name = "Tajo"
basin_id = basin_df[basin_df['NAME'] == basin_name]['WSO_ID'].unique()

# Subset the basin data
basin_data = basin_df[basin_df['WSO_ID'] == basin_id]

```

2. CRS Matching

Ensure all spatial data are in the same CRS:

```

import geopandas as gpd

# Load and transform CRS
dam_data = gpd.read_file('path_to_dam_file.shp')
dam_data = dam_data.to_crs(epsg=4326)

# Load and transform other spatial data
basin_shape = gpd.read_file('path_to_basin_shape_file.shp').to_crs(epsg=4326)

```

3. Data Preprocessing

Prepare and subset the data:

```

# Buffer basin shape
buffered_basin = basin_shape.buffer(1)

# Intersect with dams
dam_subset = gpd.overlay(dam_data, buffered_basin, how='intersection')

```

4. Geospatial Transformations

Transform and simplify geometries:

```

# Load river data
river_data = gpd.read_file('path_to_river_file.shp').to_crs(epsg=4326)

# Simplify river shapes
simplified_rivers = river_data.copy()
simplified_rivers['geometry'] = simplified_rivers['geometry'].simplify(tolerance=0.01)

```

5. Network Construction

Build and refine the river network graph:

```

import networkx as nx

# Create a graph from river data
G = nx.Graph()

# Add nodes and edges
for idx, row in river_data.iterrows():
    G.add_node(idx, **row.drop('geometry').to_dict())
    # Add edges based on river connectivity (example)
    # G.add_edge(node1, node2)

```
6. Index Calculations

Compute connectivity indices:

```

# Example function to compute an index
def compute_connectivity_index(graph):
    # Placeholder function for index calculation
    return nx.algorithms.connectivity.degree_connected_components(graph)

# Compute indices
connectivity_index = compute_connectivity_index(G)

```

7. Visualization

Plot the results:

```

import matplotlib.pyplot as plt

# Plot river network
fig, ax = plt.subplots()
simplified_rivers.plot(ax=ax, color='blue', linewidth=0.5)
dam_subset.plot(ax=ax, color='red', markersize=5)
plt.title('River Network and Dams')
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.show()
```

## Contributing

Feel free to fork the repository and submit pull requests for improvements or bug fixes.
##  License

This project is licensed under the MIT License. See the LICENSE file for details.
