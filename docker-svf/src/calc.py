import numpy as np
import pandas as pd
from geopy.distance import geodesic
import json
import math
import os
from sklearn.cluster import DBSCAN
from tqdm import tqdm
from svf import svf

def find_grid_cell(lat, lon, A, B, cell_width_m, cell_height_m, num_cells_x, num_cells_y):
    """
    Calculate which grid cell a given coordinate (lat, lon) belongs to.
    Parameters:
    - lat (float): Latitude of the point.
    - lon (float): Longitude of the point.
    - A (float): Latitude of the bottom-left corner of the grid.
    - B (float): Longitude of the bottom-left corner of the grid.
    - cell_width_m (float): Width of each cell in meters.
    - cell_height_m (float): Height of each cell in meters.
    - num_cells_x (int): Number of cells in the longitude direction.
    - num_cells_y (int): Number of cells in the latitude direction.
    Returns:
    - x_index (int): Cell index in the longitude direction (0 to num_cells_x - 1).
    - y_index (int): Cell index in the latitude direction (0 to num_cells_y - 1).
    """
    try:
        # Validate latitude and longitude ranges
        if not (-90 <= lat <= 90):
            raise ValueError(f"[calc.py]Latitude {lat} is out of range. It must be between -90 and 90 degrees.")
        if not (-180 <= lon <= 180):
            raise ValueError(f"[calc.py]Longitude {lon} is out of range. It must be between -180 and 180 degrees.")
        # Calculate the distance from the origin (A, B) to the point (Lat, Lon)
        distance_lon_m = geodesic((A, B), (A, lon)).meters
        distance_lat_m = geodesic((A, B), (lat, B)).meters
        # Calculate the cell indices
        x_index = int(distance_lon_m / cell_width_m)
        y_index = int(distance_lat_m / cell_height_m)
        # Boundary check (clamp indices to valid range)
        x_index = min(0, num_cells_x - 1) if B > lon else min(x_index, num_cells_x - 1)
        y_index = min(0, num_cells_y - 1) if A > lat else min(y_index, num_cells_y - 1)
        return x_index, y_index
    except Exception as e:
        print(f"[calc.py]An error occurred while calculating grid cell indices: {e}")
        return None, None

try:
    vectorized_find_grid_cell = np.vectorize(find_grid_cell)
except Exception as e:
    print(f"[calc.py]An error occurred while vectorizing the find_grid_cell function: {e}")

def calc_indices(lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids):
    """
    Calculate the grid cell indices for given latitude and longitude coordinates.
    Parameters:
    lats_raster (numpy.ndarray): Array of latitude values of the raster grid.
    lons_raster (numpy.ndarray): Array of longitude values of the raster grid.
    num_cells_x (int): Number of cells along the x-axis (longitude).
    num_cells_y (int): Number of cells along the y-axis (latitude).
    all_interpolated_lats (numpy.ndarray): Array of interpolated latitude values.
    all_interpolated_lons (numpy.ndarray): Array of interpolated longitude values.
    line_ids (numpy.ndarray): Array of line IDs corresponding to the interpolated points.
    Returns:
    tuple: Two numpy arrays containing the x and y indices of the grid cells for each interpolated point.
    """
    try:
        # Extract the latitude and longitude coordinates of the corners of the grid
        A = np.min(lats_raster)  # Latitude of the bottom-left corner of the grid
        C = np.max(lats_raster)  # Latitude of the top-right corner of the grid
        B = np.min(lons_raster)  # Longitude of the bottom-left corner of the grid
        D = np.max(lons_raster)  # Longitude of the top-right corner of the grid
        # Calculate the total width and height of the grid in meters
        total_width_m = geodesic((A, B), (A, D)).meters  # Distance in meters between the left and right edges
        total_height_m = geodesic((A, B), (C, B)).meters  # Distance in meters between the bottom and top edges
        # Calculate the width and height of each cell in meters
        cell_width_m = total_width_m / num_cells_x
        cell_height_m = total_height_m / num_cells_y
        # Apply the vectorized function to the arrays
        x_indices, y_indices = vectorized_find_grid_cell(all_interpolated_lats, all_interpolated_lons, A, B, cell_width_m, cell_height_m, num_cells_x, num_cells_y)
        print(f"[calc.py]Type of x_indices: {type(x_indices)}")
        print(f"[calc.py]Type of y_indices: {type(y_indices)}")
        # Convert to (x, y) format with origin at the upper-left corner
        xy_indices = [[x, num_cells_y - y] for x, y in zip(x_indices, y_indices)]
        print(f"[calc.py]Type of xy_indices: {type(xy_indices)}")
        return x_indices, y_indices, line_ids, xy_indices
    except Exception as e:
        print(f"[calc.py]An error occurred while calculating indices: {e}")
        return None, None, None, None

def calc_svf(settings, array_band, xy_indices, obstacles=None):
    """
    Calculate the Sea View Factor (SVF) for given observer positions.
    Parameters:
    settings (dict): A dictionary containing the following keys:
        - 'delta_angle' (float): The angle increment for SVF calculation.
        - 'max_distance' (float): The maximum distance to consider for obstacles.
        - 'print_info' (bool): Whether to print detailed information during calculation.
    array_band (numpy.ndarray): A 2D array representing the elevation data.
    xy_indices (list of tuples): A list of (x, y) tuples representing observer positions.
    obstacles (numpy.ndarray, optional): A 2D array representing obstacles. If not provided, an array of zeros will be used.
    Returns:
    list: A list of SVF values for each observer position.
    """
    try:
        # Extract relevant settings
        delta_angle = settings['delta_angle']
        max_distance = settings['max_distance']
        print_info = settings['print_info']
        # Initialize the list of SVF values    
        sea_view_factor_list = []
        # Create obstacles array if not provided
        if obstacles is None:
            obstacles = np.zeros((array_band.shape[0], array_band.shape[1]))
        # Calculate the SVF for each observer position
        for observer_pos in tqdm(xy_indices, desc="Calculating SVF", unit="position"):
            if print_info:
                print(f"[calc.py]Calculating SVF for observer position: {observer_pos}")
            # Calculate the SVF for the observer position
            sea_view_factor, visibility = svf(array_band, obstacles, observer_pos, delta_angle, max_distance, print_info)
            if print_info:
                print(f"[calc.py]Sea View Factor: {sea_view_factor:.2f}")
            sea_view_factor_list.append(sea_view_factor) # Append the SVF value (numpy.float64) to the list
        return sea_view_factor_list
    except Exception as e:
        print(f"[calc.py]An error occurred while calculating SVF: {e}")
        return None

def convert_to_number(s):
    """
    Convert a given input to a number (int or float).
    Parameters:
    s (int, float, str): The input value to be converted. It can be an integer, 
                            a float, or a string representing a number.
    Returns:
    int or float: The converted number if the input can be successfully converted.
    np.nan: If the input cannot be converted to a number, returns numpy's NaN.
    Raises:
    ValueError: If the input string cannot be converted to a number.
    Notes:
    - If the input is a string containing a decimal point or an exponent ('e' or 'E'), it will be converted to a float.
    - If the input is a string without a decimal point or exponent, it will be converted to an int.
    - If the conversion fails, an error message will be printed and np.nan will be returned.
    """
    try:
        if isinstance(s, int):
            return int(s)
        elif isinstance(s, float):
            return float(s)
        elif isinstance(s, str):
            # 小数として変換可能かを判定
            if '.' in s or 'e' in s.lower():
                return float(s)
            else:
                # 整数として変換
                return int(s)
    except ValueError as e:
        print(f"An Error occurred when converting to number: {e}")
        print("return np.nan")
        return np.nan

try:
    vectorized_convert_to_number = np.vectorize(convert_to_number)
except Exception as e:
    print(f"[calc.py]An error occurred while vectorizing the convert_to_number function: {e}")

def save_svf_result_to_dataframe(sea_view_factor_list, x_indices, y_indices, all_interpolated_lats, all_interpolated_lons, line_ids):
    """
    Save Sea View Factor (SVF) results to a pandas DataFrame.
    Parameters:
    sea_view_factor_list (list): List of sea view factor values.
    x_indices (list): List of x indices.
    y_indices (list): List of y indices.
    all_interpolated_lats (list): List of interpolated latitudes.
    all_interpolated_lons (list): List of interpolated longitudes.
    line_ids (list): List of line IDs.
    Returns:
    tuple: A tuple containing two pandas DataFrames:
        - result_df: DataFrame with columns 'x_index', 'y_index', 'svf', and 'line_id'.
        - result_coordinates_df: DataFrame with columns 'lon', 'lat', 'svf', and 'line_id'.
    Raises:
    ValueError: If the lengths of the input lists do not match.
    Exception: If an error occurs during the conversion of input lists to numpy arrays or DataFrames.
    """
    try:
        # Check if the lengths of the input lists match
        print(f"[calc.py]data x_indices: {x_indices}", f"[calc.py]data y_indices: {y_indices}", f"[calc.py]data sea_view_factor_list: {sea_view_factor_list}", f"[calc.py]data line_ids: {line_ids}")
        print(f"[calc.py]Length of x_indices: {len(x_indices)}", f"[calc.py]Length of y_indices: {len(y_indices)}", f"[calc.py]Length of sea_view_factor_list: {len(sea_view_factor_list)}", f"[calc.py]Length of line_ids: {len(line_ids)}")
        if len(x_indices) != len(sea_view_factor_list) or len(sea_view_factor_list) != len(line_ids) or len(line_ids) != len(y_indices):
            raise ValueError("[calc.py]Lengths of input lists do not match.")
        # Ensure that xy_indices and sea_view_factor_list can be combined into an array
        try:
            x_indices = np.array(x_indices)
            converted_x_indices = vectorized_convert_to_number(x_indices)
            y_indices = np.array(y_indices)
            converted_y_indices = vectorized_convert_to_number(y_indices)
            sea_view_factor_list = np.array(sea_view_factor_list)
            converted_sea_view_factor_list = vectorized_convert_to_number(sea_view_factor_list)
            line_ids = np.array(line_ids)
            all_interpolated_lats = np.array(all_interpolated_lats)
            converted_all_interpolated_lats = vectorized_convert_to_number(all_interpolated_lats)
            print(f"[calc.py]Data converted_all_interpolated_lats {converted_all_interpolated_lats}")
            all_interpolated_lons = np.array(all_interpolated_lons)
            converted_all_interpolated_lons = vectorized_convert_to_number(all_interpolated_lons)
            print(f"[calc.py]Data converted_all_interpolated_lons {converted_all_interpolated_lons}")
            # 各配列の dtype を表示
            print(f"[calc.py]converted_x_indices dtype: {converted_x_indices.dtype}")
            print(f"[calc.py]converted_y_indices dtype: {converted_y_indices.dtype}")
            print(f"[calc.py]converted_sea_view_factor_list dtype: {converted_sea_view_factor_list.dtype}")
            print(f"[calc.py]line_ids dtype: {line_ids.dtype}")
            print(f"[calc.py]converted_all_interpolated_lats dtype: {converted_all_interpolated_lats.dtype}")
            print(f"[calc.py]converted_all_interpolated_lons dtype: {converted_all_interpolated_lons.dtype}")
        except Exception as e:
            print(f"[calc.py]An error occurred while converting input lists to pd.DataFrame: {e}")
            return None, None
        # Combine the data into a DataFrame instead of a Numpy Array because line_id is a string
        result_df = pd.DataFrame({'x_index': converted_x_indices, 'y_index': converted_y_indices, 'svf': converted_sea_view_factor_list, 'line_id': line_ids})
        print(f"result dtype: {result_df.dtypes}")
        # Save the results to a DataFrame with latitude and longitude
        result_coordinates_df = pd.DataFrame({'lon': converted_all_interpolated_lons, 'lat': converted_all_interpolated_lats, 'svf': converted_sea_view_factor_list, 'line_id': line_ids})
        print(f"result_coordinates dtype: {result_coordinates_df.dtypes}")
        return result_df, result_coordinates_df
    except Exception as e:
        print(f"[calc.py]An error occurred while saving the results: {e}")
        return None, None

def meters_to_radians(distance_meters):
    """
    Convert a distance from meters to radians.
    Parameters:
    distance_meters (float): The distance in meters to be converted.
    Returns:
    float: The distance in radians.
    """
    # Earth radius in meters
    earth_radius_m = 6371000
    # Convert distance in meters to radians
    return distance_meters / earth_radius_m

def filter_and_cluster_svf_result(svf_result_coordinates_df, svf_threshold=0.07, distance_threshold_m=1000):
    """
    Filters SVF (Sky View Factor) results based on a threshold and clusters the filtered points using DBSCAN clustering.
    Parameters:
    svf_result_coordinates_df (pd.DataFrame): DataFrame containing SVF results with 'lon', 'lat', and 'svf' columns.
    svf_threshold (float, optional): Threshold value for filtering SVF results. Points with SVF values greater than or equal to this threshold are considered. Default is 0.07.
    distance_threshold_m (float, optional): Distance threshold in meters for DBSCAN clustering. Points within this distance are considered in the same cluster. Default is 1000 meters.
    Returns:
    dict: A dictionary where keys are cluster labels and values are lists of points (tuples of longitude and latitude) in each cluster. Returns None if an error occurs.
    """
    try:
        # Filter the SVF results based on the threshold value
        points = []
        for index, row in svf_result_coordinates_df.iterrows(): # use iterrows() instead of values() because line_id is a string
            try:
                if row['svf'] >= svf_threshold:
                    points.append((row['lat'], row['lon']))  # (lat, lon)
            except ValueError as e:
                print(f"[calc.py]An error occurred while filtering the SVF results: {e}")
        # DBSCAN clustering # Threshold distance (e.g., points within 1000 meters are considered in the same cluster)
        distance_threshold_radians = meters_to_radians(distance_threshold_m)
        # Perform DBSCAN clustering
        db = DBSCAN(eps=distance_threshold_radians, min_samples=1, algorithm='ball_tree', metric='haversine').fit(np.radians(points))
        # Get the cluster labels
        labels = db.labels_
        # Group the points by cluster
        clusters = {}
        for label, point in zip(labels, points):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append(point)
        return clusters
    except Exception as e:
        print(f"[calc.py]An error occurred while filtering and clustering the SVF results: {e}")
        return None

def calculate_circle(points):
    """
    Calculate the center and radius of the smallest circle that can enclose all given points.
    Parameters:
        points (list of tuple): A list of tuples where each tuple contains the latitude and longitude of a point.
    Returns:
        tuple: A tuple containing the center (latitude, longitude) of the circle and the radius in meters.
                If an error occurs, returns (None, None).
    Raises:
        Exception: If an error occurs during the calculation, it will be caught and printed.
    """
    try:
        # Extract the latitude and longitude coordinates of the points
        latitudes = [point[0] for point in points] # (lat, lon)
        longitudes = [point[1] for point in points]
        # Calculate the center of the cluster
        center_lat = np.mean(latitudes)
        center_lon = np.mean(longitudes)
        # Calculate the center of the cluster
        center = (center_lat, center_lon)
        # Calculate the maximum distance and use it as the radius
        radius = max(geodesic(center, point).meters for point in points)    
        return center, radius
    except Exception as e:
        print(f"[calc.py]An error occurred while calculating the circle: {e}")
        return None, None

def create_circle_polygon_in_geojson(center, radius, num_points=100):
    """
    Create a polygon representing a circle around a given center point.
    Parameters:
        center (tuple): A tuple containing the latitude and longitude of the center point (lat, lon).
        radius (float): The radius of the circle in meters.
        num_points (int, optional): The number of points to generate for the circle polygon. Default is 100.
    Returns:
        list: A list of [longitude, latitude] pairs representing the circle polygon because GeoJSON expects [lon, lat].
                Returns None if an error occurs.
    Raises:
        Exception: If an error occurs during the creation of the circle polygon, an exception is caught and an error message is printed.
    """
    try:
        # Create a polygon representing a circle around the
        lat, lon = center
        radius_in_radians = meters_to_radians(radius)  # Convert meters to degrees using geopy # ここ要修正
        # Create the circle polygon
        polygon = []
        for i in range(num_points):
            angle = math.radians(float(i) / num_points * 360) # Convert to radians
            dx = radius_in_radians * math.cos(angle)
            dy = radius_in_radians * math.sin(angle)
            polygon.append([lon + dx, lat + dy])
        # Close the polygon by adding the first point again
        polygon.append(polygon[0])
        return polygon
    except Exception as e:
        print(f"[calc.py]An error occurred while creating the circle polygon: {e}")
        return None

def build_target_area(svf_result_coordinates):
    """
    Builds a target area from SVF result coordinates.
    This function filters and clusters the SVF (Sky View Factor) results, calculates the center and radius of each cluster,
    and creates a GeoJSON FeatureCollection representing the target area.
    Parameters:
        svf_result_coordinates (list): A list of SVF result coordinates.
    Returns:
        dict: A dictionary representing the target area in GeoJSON format, or None if an error occurs.
    Raises:
        Exception: If an error occurs during the processing of SVF result coordinates.
    """
    try:
        # Filter and cluster the SVF results
        clusters = filter_and_cluster_svf_result(svf_result_coordinates)
        # Calculate the center and radius of each cluster
        circles = [calculate_circle(cluster) for cluster in clusters.values()] # cluster is a list of points
        # Create a GeoJSON FeatureCollection
        features = []
        for center, radius in circles:
            feature = {
            "type": "Point",
            "coordinates": [center[1], center[0]],  # [lon, lat]
            "properties": {
                "radius": radius
                }
            }
            features.append(feature)
        target_area = {"target_area": features}
        return target_area
    except Exception as e:
        print(f"[calc.py]An error occurred while building the target area: {e}")
        return None

def build_target_area_polygons(target_area):
    """
    Converts a target area defined by circles into a GeoJSON FeatureCollection of polygons.
    Parameters:
        target_area (dict): A dictionary containing the target area information. 
                            It should have a key 'target_area' which is a list of features.
                            Each feature should be a dictionary with 'coordinates' (list of [longitude, latitude])
                            and 'properties' (dictionary containing 'radius' and other properties).
    Returns:
        dict: A GeoJSON FeatureCollection containing the polygons created from the circles.
                Returns None if an error occurs during the conversion process.
    Raises:
        Exception: If an error occurs during the conversion process, it will be caught and printed.
    """
    try:
        # Create a new GeoJSON FeatureCollection
        target_area_polygons = {
            "type": "FeatureCollection",
            "features": []
        }
        # Convert each circle to a polygon
        for feature in target_area['target_area']:
            center = (feature['coordinates'][1], feature['coordinates'][0])  # (lat, lon)
            radius = feature['properties']['radius']
            polygon = create_circle_polygon_in_geojson(center, radius)
            polygon_feature = {
                "type": "Feature",
                "geometry": {
                    "type": "Polygon",
                    "coordinates": [polygon]
                },
                "properties": feature['properties']
            }
            target_area_polygons['features'].append(polygon_feature)
        return target_area_polygons
    except Exception as e:
        print(f"[calc.py]An error occurred while building the target area polygons: {e}")
        return None

def save_area_to_geojson(area, file_name):
    """
    Save the given area data to a GeoJSON file.
    Parameters:
    area (dict): The area data to be saved, typically in GeoJSON format.
    file_name (str): The name of the file where the area data will be saved.
    Raises:
    Exception: If an error occurs while saving the area data to the GeoJSON file.
    """
    try:
        # Save the target area to a GeoJSON file
        with open(file_name, 'w') as f:
            json.dump(area, f)
        print(f"[calc.py]GeoJSON file created: {file_name}")
    except Exception as e:
        print(f"[calc.py]An error occurred while saving the target area to GeoJSON: {e}")

def convert_to_target_area(svf_result_coordinates, output_folder):
    """
    Converts the SVF result coordinates to a target area and saves it to GeoJSON files.
    Parameters:
        svf_result_coordinates (list): A list of coordinates representing the SVF results.
        output_folder (str): The folder where the GeoJSON files will be saved. If None, the files will not be saved.
    Returns:
        tuple: A tuple containing the target area and target area polygons. If an error occurs, returns (None, None).
    Raises:
        Exception: If an error occurs during the conversion or file saving process, it will be caught and printed.
    """
    try:
        # Save the target area to a GeoJSON file
        target_area = build_target_area(svf_result_coordinates)
        if output_folder:
            save_area_to_geojson(target_area, os.path.join(output_folder, 'target_area.geojson'))
        # Save the target area polygons to a GeoJSON file
        target_area_polygons = build_target_area_polygons(target_area)
        if output_folder:
            save_area_to_geojson(target_area_polygons, os.path.join(output_folder, 'target_area_polygons.geojson'))
        return target_area, target_area_polygons
    except Exception as e:
        print(f"[calc.py]An error occurred while converting to target area: {e}")
        return None, None