import geopandas as gpd
import json
from shapely.geometry import shape
import zipfile
import pandas as pd
import os
import shutil
from interpolate import geojson_to_interpolated_points, shapes_df_to_interpolated_points, stops_df_to_interpolated_points

def load_filter_railroad_geojson(input_file_path, filter_key, filter_value):
    """
    Load a GeoJSON file, filter its contents based on a specified key-value pair, 
    and save the filtered data to a new GeoJSON file.
    Parameters:
    input_file_path (str): The path to the input GeoJSON file.
    filter_key (str): The key in the GeoJSON properties to filter by.
    filter_value (str): The value in the GeoJSON properties to filter by.
    Returns:
    GeoDataFrame: The filtered GeoDataFrame.
    Raises:
    ValueError: If no data is found for the specified key-value pair, 
                if the GeoJSON file has unsupported types, 
                or if there are issues with the GeoJSON structure.
    KeyError: If the specified filter_key is not found in the GeoJSON properties.
    Exception: For any other unexpected errors during processing.
    """
    try:
        # GeoJSONファイルの読み込み
        data = gpd.read_file(input_file_path)
        # 指定されたキーと値でフィルタリング
        filtered_gdf = None
        if filter_key is not None and filter_value is not None:
            filter_value_list = list(filter_value)
            for filter_value in filter_value_list:
                filtered_gdf_ = data[data[filter_key] == filter_value]
                if filtered_gdf is None:
                    filtered_gdf = filtered_gdf_
                else:
                    filtered_gdf = pd.concat([filtered_gdf, filtered_gdf_])
        # フィルタリングされたデータを新しいGeoJSONファイルに保存
        if filtered_gdf.empty:
            raise ValueError(f"[line.py]No data found for {filter_key} = {filter_value}.")
        else:
            output_file_path = input_file_path.replace('.geojson', f'_filtered.geojson')
            filtered_gdf.to_file(output_file_path, driver='GeoJSON')
            print(filtered_gdf.info())
            print(f"[line.py]Filtered data saved to {output_file_path}")
            return filtered_gdf
    except KeyError as e:
        raise ValueError(f"[line.py]KeyError: {e}. Consider changing key names in GeoJSON properties.")
    except Exception as e:
        # 失敗した場合の処理
        print(f"[line.py]Failed to load GeoDataFrame: {e} Try loading the GeoJSON file.")
        # ファイルを手動処理する
        try:
            with open(input_file_path, 'r') as f:
                geojson = json.load(f)
            # GeoJSON の type に基づく処理
            if geojson['type'] == 'LineString':
                geometry = shape(geojson)
                gdf = gpd.GeoDataFrame([{'geometry': geometry}], geometry='geometry')
            elif geojson['type'] == 'MultiLineString':
                geometries = [shape({'type': 'LineString', 'coordinates': coords}) for coords in geojson['coordinates']]
                gdf = gpd.GeoDataFrame([{'geometry': geom} for geom in geometries], geometry='geometry')
            else:
                raise ValueError(f"[line.py]Unsupported GeoJSON type: {geojson['type']}")
            # 指定されたキーと値でフィルタリング
            filtered_gdf = data[data[filter_key] == filter_value]
            # フィルタリングされたデータを新しいGeoJSONファイルに保存
            if filtered_gdf.empty:
                raise ValueError(f"[line.py]No data found for {filter_key} = {filter_value}.")
            else:
                output_file_path = input_file_path.replace('.geojson', f'_filtered.geojson')
                filtered_gdf.to_file(output_file_path, driver='GeoJSON')
                return filtered_gdf
        except json.JSONDecodeError as jde:
            raise ValueError(f"[line.py]Failed to parse GeoJSON: {jde}")
        except KeyError as ke:
            raise ValueError(f"[line.py]GeoJSON missing required key: {ke}")
        except Exception as inner_e:
            raise ValueError(f"[line.py]An unexpected error occurred during manual processing: {inner_e}")

def load_gtfs_files(zip_file_path):
    """
    Extracts GTFS zip file and checks for the presence of shapes.txt and stops.txt.
    Returns dataframes for shapes.txt and stops.txt (or None if not present).
    Parameters:
    - zip_file_path (str): Path to the GTFS zip file.
    Returns:
    - shapes_df (pd.DataFrame or None): Dataframe for shapes.txt or None if not present.
    - stops_df (pd.DataFrame or None): Dataframe for stops.txt or None if not present.
    """
    try:
        shapes_df = None
        stops_df = None
        with zipfile.ZipFile(zip_file_path, 'r') as zip_ref:
            # Extract all files to a temporary directory
            temp_dir = os.path.join(os.path.dirname(zip_file_path), 'temp_gtfs')
            zip_ref.extractall(temp_dir)
            # Check for shapes.txt and stops.txt
            shapes_file_path = os.path.join(temp_dir, 'shapes.txt')
            stops_file_path = os.path.join(temp_dir, 'stops.txt')
            if os.path.exists(shapes_file_path):
                shapes_df = pd.read_csv(shapes_file_path)
            if os.path.exists(stops_file_path):
                stops_df = pd.read_csv(stops_file_path)
        # Clean up the temporary directory
        shutil.rmtree(temp_dir)
        return shapes_df, stops_df
    except zipfile.BadZipFile as e:
        raise Warning(f"[line.py]BadZipFile: {e}. The provided file is not a valid zip file.")
    except Exception as e:
        raise Warning(f"[line.py]An unexpected error occurred while extracting GTFS files: {e}")

def load_filter_gtfs_files(zip_file_path, filter_key, filter_value):
    """
    Filters GTFS shapes and stops dataframes based on a specified key and value.
    Parameters:
    - shapes_df (pd.DataFrame): Dataframe for shapes.txt.
    - stops_df (pd.DataFrame): Dataframe for stops.txt.
    - filter_key (str): The key to filter the data.
    - filter_value (str): The value to filter the data.  
    Returns:
    - filtered_shapes_df (pd.DataFrame): Filtered dataframe for shapes.txt.
    - filtered_stops_df (pd.DataFrame): Filtered dataframe for stops.txt.
    """
    try:
        # Extract GTFS files
        shapes_df, stops_df = load_gtfs_files(zip_file_path)
        # Filter shapes and stops dataframes
        filtered_shapes_df = None
        filtered_stops_df = None
        if shapes_df is not None:
            if filter_key is not None and filter_value is not None:
                filter_value_list = list(filter_value)
                for filter_value in filter_value_list:
                    print(f"[line.py]Filtering shapes_df by {filter_key} = {filter_value}")
                    filtered_shapes_df_ = shapes_df[shapes_df[filter_key] == filter_value]
                    if filtered_shapes_df is None:
                        filtered_shapes_df = filtered_shapes_df_
                    else:
                        filtered_shapes_df = pd.concat([filtered_shapes_df, filtered_shapes_df_])
                return filtered_shapes_df, None
            else:
                filtered_shapes_df = shapes_df
                return filtered_shapes_df, None
        elif stops_df is not None:
            if filter_key is not None and filter_value is not None:
                filter_value_list = list(filter_value)
                for filter_value in filter_value_list:
                    print(f"[line.py]Filtering stops_df by {filter_key} = {filter_value}")
                    filtered_stops_df_ = stops_df[stops_df[filter_key] == filter_value]
                    if filtered_stops_df is None:
                        filtered_stops_df = filtered_stops_df_
                    else:
                        filtered_stops_df = pd.concat([filtered_stops_df, filtered_stops_df_])
                return None, filtered_stops_df
            else:
                filtered_stops_df = stops_df
                return None, filtered_stops_df
    except zipfile.BadZipFile as e:
        raise Warning(f"[line.py]BadZipFile: {e}. The provided file is not a valid zip file.")
    except Exception as e:
        raise Warning(f"[line.py]An unexpected error occurred while extracting GTFS files: {e}")

def load_filter_interpolate_line(input_file_path, input_file_type, filter_key, filter_value, id_key, distance_m=10):
    """
    Load, filter, and interpolate line data from a specified input file.
    Parameters:
    input_file_path (str): The path to the input file.
    input_file_type (str): The type of the input file. Supported values are 'railway_geojson' and 'gtfs'.
    filter_key (str): The key to filter the data.
    filter_value (str): The value to filter the data.
    distance_m (int, optional): The distance in meters for interpolation. Default is 10.
    Returns:
    tuple: A tuple containing three elements:
        - all_interpolated_lats (list): A list of interpolated latitudes.
        - all_interpolated_lons (list): A list of interpolated longitudes.
        - line_ids (list): A list of line IDs.
    Raises:
    ValueError: If the input file type is invalid or if neither stops.txt nor shapes.txt is found in GTFS files.
    """
    if input_file_type == "railway_geojson":
        gdf = load_filter_railroad_geojson(input_file_path, filter_key, filter_value)
        try:
            all_interpolated_lats, all_interpolated_lons, line_ids = geojson_to_interpolated_points(gdf, distance_m, id_key)
            return all_interpolated_lats, all_interpolated_lons, line_ids
        except KeyError as e:
            raise Warning(f"[line.py]KeyError: {e}. Consider changing key names in GeoJSON properties.")
        except Exception as e:
            raise Warning(f"[line.py]An unexpected error occurred: {e}")
    elif input_file_type == "gtfs":
        filtered_shapes_df, filtered_stops_df = load_filter_gtfs_files(input_file_path, filter_key, filter_value)
        if filtered_shapes_df is not None:
            try:
                all_interpolated_lats, all_interpolated_lons, line_ids = shapes_df_to_interpolated_points(filtered_shapes_df, distance_m, shapes_lon_key='shape_pt_lon', shapes_lat_key='shape_pt_lat', shape_id_key='shape_id')
                return all_interpolated_lats, all_interpolated_lons, line_ids
            except KeyError as e:
                raise Warning(f"[line.py]KeyError: {e}. Consider changing key names in shapes.txt.")
            except Exception as e:
                raise Warning(f"[line.py]An unexpected error occurred: {e}")
        elif filtered_stops_df is not None:
            try:
                all_interpolated_lats, all_interpolated_lons, line_ids = stops_df_to_interpolated_points(filtered_stops_df, distance_m, stops_lon_key='stop_lon', stops_lat_key='stop_lat', stops_id_key='stop_id')
                return all_interpolated_lats, all_interpolated_lons, line_ids
            except KeyError as e:
                raise Warning(f"[line.py]KeyError: {e}. Consider changing key names in stops.txt.")
            except Exception as e:
                raise Warning(f"[line.py]An unexpected error occurred: {e}")
        else:
            raise ValueError("[line.py]Neither of stops.txt nor shapes.txt is found in GTFS files.")
    else:
        raise ValueError(f"[line.py]Invalid input file type: {input_file_type}. Supported file types are 'railway_geojson' and 'gtfs'.")
