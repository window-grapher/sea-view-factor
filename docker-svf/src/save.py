import os
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
from calc import convert_to_number

def plot_save_xy_indices_raster(xy_indices, num_cells_x, num_cells_y, raster=None, file_name=None):
    """
    Plots and optionally saves a scatter plot of xy indices on a raster grid.
    Parameters:
    xy_indices (list of tuples): List of (x, y) indices to plot.
    num_cells_x (int): Number of cells in the x direction of the raster grid.
    num_cells_y (int): Number of cells in the y direction of the raster grid.
    raster (2D array, optional): Raster data to be plotted as a background. Default is None.
    file_name (str, optional): Name of the file to save the plot. If None, the plot is not saved. Default is 'xy_indices_and_raster'.
    Returns:
    None
    """
    try:
        # Plot the interpolated points on the raster grid
        plt.figure(figsize=(num_cells_x/200, num_cells_y/200))
        # Plot the raster if available
        if raster is not None:
            raster_processed = np.squeeze(raster) if len(raster.shape) == 3 else raster
            plt.imshow(raster_processed, extent=[0, num_cells_x, 0, num_cells_y], origin='upper', cmap='terrain', alpha=0.5)
        # Plot the xy_indices
        plt.scatter([xy_indices[i][0] for i in range(len(xy_indices))], [num_cells_y-xy_indices[i][1] for i in range(len(xy_indices))], s=0.1)
        # Save / Show the plot
        if file_name:
            save_path = file_name + '.png'
            plt.savefig(save_path)
            plt.show()
        else:
            plt.show()
    except Exception as e:
        print(f"[save.py]An error occurred while plotting and saving the xy indices raster: {e}")

def save_svf_result_to_geojson(svf_result_coordinates_df, file_name):
    """
    Saves the Sky View Factor (SVF) results to a GeoJSON file.
    Parameters:
    svf_result_coordinates_df (pandas.DataFrame): DataFrame containing the SVF results with columns 'lon', 'lat', 'svf', and 'line_id'.
    file_name (str): The name of the output GeoJSON file (without the .geojson extension).
    Returns:
    None
    Raises:
    Exception: If an error occurs during the process, it prints an error message with the exception details.
    """
    try:
        features = []
        for index, row in svf_result_coordinates_df.iterrows(): # use iterrows() instead of values() because line_id is a string
            lon = convert_to_number(row['lon'])
            lat = convert_to_number(row['lat'])
            svf = convert_to_number(row['svf'])
            line_id = row['line_id'] # line_id is a string
            # Create a GeoJSON feature for each point
            feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Point",
                    "coordinates": [lon, lat]  # Use [longitude, latitude] as GeoJSON expects
                },
                "properties": {
                    "svf": svf,
                    "line_id": line_id
                }
            }
            features.append(feature)

        geojson = {
            "type": "FeatureCollection",
            "features": features
        }

        with open(f"{file_name}.geojson", 'w') as f:
            json.dump(geojson, f, indent=4)
    except Exception as e:
        print(f"[save.py]An error occurred while saving the SVF result to GeoJSON: {e}")

def plot_svf_results(svf_result_df, num_cells_x, num_cells_y, raster=None, file_name=None):
    """
    Plots the Sea View Factor (SVF) results on a raster grid.
    Parameters:
    svf_result_df (DataFrame): DataFrame containing the SVF results with columns 'x_index', 'y_index', and 'svf'.
    num_cells_x (int): Number of cells in the x-direction of the raster grid.
    num_cells_y (int): Number of cells in the y-direction of the raster grid.
    raster (ndarray, optional): 2D array representing the raster data to be plotted as the background. Default is None.
    file_name (str, optional): Name of the file to save the plot. If None, the plot will not be saved. Default is None.
    Returns:
    None
    Raises:
    Exception: If an error occurs during plotting, it will be caught and printed.
    """
    try:
        # Extract the SVF values and line IDs
        x_indices_array = svf_result_df['x_index'].values
        y_indices_array = svf_result_df['y_index'].values
        svf_array = svf_result_df['svf'].values
        # Plot the SVF values on the raster grid
        plt.figure(figsize=(num_cells_x/200, num_cells_y/200))
        # Plot the raster if available
        if raster is not None:
            plt.imshow(raster, extent=[0, num_cells_x, 0, num_cells_y], origin='upper', cmap='terrain', alpha=0.5)
        # Plot the SVF values
        plt.scatter([x_index for x_index in x_indices_array], [y_index for y_index in y_indices_array], c=[svf for svf in svf_array], cmap='Reds', s=5)
        # Add colorbar and title
        plt.colorbar(label='Sea View Factor')
        plt.title('Sea View Factor Results')
        # Save / Show the plot
        if file_name:
            save_path = file_name + '.png'
            plt.savefig(save_path)
            plt.show()
    except Exception as e:
        print(f"[save.py]An error occurred while plotting the SVF results: {e}")

def save_list(data_list, file_name):
    """
    Saves a list of data to a text file, with each item on a new line.
    Parameters:
        data_list (list): The list of data to save.
        file_name (str): The name of the file to save the data to.
    Returns:
        None
    """
    try:
        with open(f"{file_name}.txt", 'w') as f:
            for item in data_list:
                f.write(f"{item}\n")
    except Exception as e:
        print(f"[save.py]An error occurred while saving the list: {e}")

def save_dataframe(df, file_name):
    """
    Saves a pandas DataFrame to a CSV file.
    Parameters:
        df (pandas.DataFrame): The DataFrame to save.
        file_name (str): The name of the file to save the DataFrame to.
    Returns:
        None
    """
    try:
        df.to_csv(f"{file_name}.csv", index=False)
        print(f"[save.py]DataFrame saved to {file_name}.csv")
    except Exception as e:
        print(f"[save.py]An error occurred while saving the DataFrame: {e}")

def load_list(file_name):
    """
    Loads a list of data from a text file, with each item on a new line.
    Parameters:
        file_name (str): The name of the file to load the data from.
    Returns:
        list: The list of data loaded from the file.
    """
    try:
        with open(f"{file_name}.txt", 'r') as f:
            data_list_str = f.read().splitlines() # Read the file and split the lines into a list
            data_list = []
            for line in data_list_str:
                try:
                    match = re.findall(r"(?:np\.\w+\()?([\d\.\-e]+)\)?", line) # Extract numbers from np.XXXX(number) or simple number
                    if match:
                        parsed = [float(m) if "." in m or "e" in m else int(m) for m in match] # Convert to float if decimal point or scientific notation is present
                        if len(parsed) == 1:
                            data_list.append(parsed[0])
                        else:
                            data_list.append(parsed)
                except Exception as e:
                    print(f"Skipping invalid line: {line}, Error: {e}")
        print(f"[save.py]Loaded data from {file_name}.txt")
        return data_list
    except Exception as e:
        print(f"[save.py]An error occurred while loading the list: {e}")
        return []

def xy_indices_to_x_y_indices(xy_indices):
    x_indices = [int(pair[0]) for pair in xy_indices]
    y_indices = [int(pair[1]) for pair in xy_indices]
    return x_indices, y_indices