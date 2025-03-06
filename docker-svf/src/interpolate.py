from geopy.distance import geodesic
import pandas as pd
import json

def interpolate_points(start, end, distance_m):
    """
    Interpolates points along a line from start to end at specified distance intervals.
    Parameters:
        start (list): The starting point as [latitude, longitude].
        end (list): The ending point as [latitude, longitude].
        distance_m (float): The distance in meters between interpolated points.
    Returns:
        list: A list of interpolated points, each as [latitude, longitude]
    """
    lats_line = [start[0]]
    lons_line = [start[1]]
    total_distance = geodesic(start, end).meters # EPSG:4326でもメートルで計算される
    if total_distance <= distance_m:
        return lats_line, lons_line
    num_interpolated_points = int(total_distance // distance_m) # 端点を除いた間隔の数を定める
    for i in range(1, num_interpolated_points + 1): # distance_m以下の距離での補間を行う
        interpolation_fraction = i / num_interpolated_points
        lat = start[0] + interpolation_fraction * (end[0] - start[0])
        lon = start[1] + interpolation_fraction * (end[1] - start[1])
        lats_line.append(lat)
        lons_line.append(lon)
    return lats_line, lons_line

def interpolate_linestring(coordinates, distance_m):
    """
    Interpolates points along a linestring at a specified distance interval.
    Used instead of the shapes.txt of GTFS Dataset.
    Parameters:
        coordinates (list of tuples): A list of (longitude, latitude) tuples representing the linestring.
        distance_m (float): The distance interval in meters at which to interpolate points.
    Returns:
        list of lists: A list of interpolated points, where each point is represented as a [latitude, longitude] list.
    """
    try:
        interpolated_lats = []
        interpolated_lons = []
        for i in range(len(coordinates) - 1):
            start = [coordinates[i][1], coordinates[i][0]]  # 緯度・経度の順に変換
            end = [coordinates[i+1][1], coordinates[i+1][0]]
            segment_lats, segment_lons = interpolate_points(start, end, distance_m)
            interpolated_lats.extend(segment_lats[:-1])  # 最後のポイントは次のセグメントの開始点になるので除外
            interpolated_lons.extend(segment_lons[:-1])  # 最後のポイントは次のセグメントの開始点になるので除外
        interpolated_lats.append(coordinates[-1][1])  # 最後のポイントを追加
        interpolated_lons.append(coordinates[-1][0])  # 最後のポイントを追加
        return interpolated_lats, interpolated_lons
    except Exception as e:
        raise RuntimeError(f"[interpolate.py]An error occurred while interpolating linestring: {e}")

def geojson_to_interpolated_points(gdf, distance_m, id_key=None):
    """
    Interpolates points along LineString and MultiLineString geometries in a GeoDataFrame.
    Parameters:
    gdf (GeoDataFrame): Input GeoDataFrame containing geometries to be interpolated.
    distance_m (float): Distance in meters between interpolated points.
    id_key (str, optional): Key to identify unique lines. If None, no unique identification is used.
    Returns:
    tuple: A tuple containing three lists:
        - all_interpolated_lats (list): List of interpolated latitudes.
        - all_interpolated_lons (list): List of interpolated longitudes.
        - line_ids (list): List of line IDs corresponding to the interpolated points.
    Raises:
    ValueError: If a required key is missing in the GeoJSON data.
    RuntimeError: If an error occurs while processing the GeoJSON data.
    """
    all_interpolated_lats = []
    all_interpolated_lons = []
    line_ids = []
    current_line_lats = []
    current_line_lons = []
    previous_id = None
    try:
        # 入力はGeoDataFrameを想定
        # FeatureCollectionの場合
        for index, feature in gdf.iterrows(): # GeoDataFrameの各行=各featureを取得
            geometry = feature.geometry # geometryを取得
            geom_type = feature.geometry.geom_type # geometryのタイプを取得
            id = feature[id_key] if feature[id_key] else None  # idの取得(shape_idを想定) # id_keyがNoneの場合はNoneになる
            # idが異なっても1つのリストに追加する # idがNoneの場合はFalseになる
            if id != previous_id and current_line_lats and current_line_lons:
                all_interpolated_lats.extend(current_line_lats) # appendではなくextendにすることでリストを結合して追加
                all_interpolated_lons.extend(current_line_lons) # appendではなくextendにすることでリストを結合して追加
                line_ids.extend([previous_id] * len(current_line_lats))  # shape_idを追加
                current_line_lats = []
                current_line_lons = []
            # 線の補間
            if geom_type == 'LineString':
                coordinates = list(geometry.coords)
                current_line_lats.extend(interpolate_linestring(coordinates, distance_m)[0])
                current_line_lons.extend(interpolate_linestring(coordinates, distance_m)[1])
            elif geom_type == 'MultiLineString':
                for line in geometry.geoms:
                    coordinates = list(line.coords)
                    current_line_lats.extend(interpolate_linestring(coordinates, distance_m)[0])
                    current_line_lons.extend(interpolate_linestring(coordinates, distance_m)[1])
            # その他の場合
            else:
                print(f"[interpolate.py]Line GeoDataFrame Row {index}: Unsupported GeoJSON geometry type: {geom_type}")
            # idの更新 # idがNoneの場合はNoneに更新される
            previous_id = id
        # 最後の線の追加
        if current_line_lats and current_line_lons:
            all_interpolated_lats.extend(current_line_lats)  # 最後の線を追加
            all_interpolated_lons.extend(current_line_lons)  # 最後の線を追加
            line_ids.extend([previous_id] * len(current_line_lats)) # 最後の線のidを追加 # idがない場合はNoneが追加される
        print(f'[interpolate.py]length of all_interpolated_lats: {len(all_interpolated_lats)}', f'length of all_interpolated_lons: {len(all_interpolated_lons)}', f'length of line_ids: {len(line_ids)}')
        return all_interpolated_lats, all_interpolated_lons, line_ids
    except KeyError as e:
        raise ValueError(f"[interpolate.py]Missing key in GeoJSON data: {e}")
    except Exception as e:
        raise RuntimeError(f"[interpolate.py]An error occurred while processing the GeoJSON data: {e}")

def shapes_df_to_interpolated_points(shapes_df, distance_m, shapes_lon_key, shapes_lat_key, shape_id_key):
    """
    Interpolates points along the lines defined by the shapes in the given DataFrame.
    Parameters:
    shapes_df (pd.DataFrame): DataFrame containing shape information with columns for longitude, latitude, and shape ID.
    distance_m (float): Distance in meters between interpolated points.
    stops_lon_key (str, optional): Column name for longitude values in shapes_df. Default is 'stop_lon'.
    stops_lat_key (str, optional): Column name for latitude values in shapes_df. Default is 'stop_lat'.
    shape_id_key (str, optional): Column name for shape ID values in shapes_df. Default is 'shape_id'.
    Returns:
    tuple: A tuple containing three lists:
        - all_interpolated_lats (list): List of interpolated latitude values.
        - all_interpolated_lons (list): List of interpolated longitude values.
        - line_ids (list): List of shape IDs corresponding to each interpolated point.
    """
    all_interpolated_lats = []
    all_interpolated_lons = []
    line_ids = []
    try:
        for shape_id in shapes_df[shape_id_key].unique():
            shape = shapes_df[shapes_df[shape_id_key] == shape_id]
            lats = shape[shapes_lat_key].tolist()
            lons = shape[shapes_lon_key].tolist()
            coordinates = list(zip(lons, lats))
            current_line_lats, current_line_lons = interpolate_linestring(coordinates, distance_m)
            all_interpolated_lats.extend(current_line_lats)
            all_interpolated_lons.extend(current_line_lons)
            line_ids.extend([shape_id] * len(current_line_lats))
        print(f'[interpolate.py]length of all_interpolated_lats: {len(all_interpolated_lats)}', f'length of all_interpolated_lons: {len(all_interpolated_lons)}', f'length of line_ids: {len(line_ids)}')
        return all_interpolated_lats, all_interpolated_lons, line_ids
    except KeyError as e:
        raise ValueError(f"[interpolate.py]Missing key in shapes DataFrame: {e}")
    except Exception as e:
        raise RuntimeError(f"[interpolate.py]An error occurred while processing the shapes DataFrame: {e}")

def stops_df_to_interpolated_points(stops_df, distance_m, stops_lon_key, stops_lat_key, stops_id_key):
    """
    Interpolates points along the lines defined by the stop coordinates in the given DataFrame.
    Parameters:
        stops_df (pd.DataFrame): DataFrame containing stop coordinates.
        distance_m (float): Distance in meters between interpolated points.
        stops_lon_key (str, optional): Column name for longitude values in stops_df. Defaults to 'stop_lon'.
        stops_lat_key (str, optional): Column name for latitude values in stops_df. Defaults to 'stop_lat'.
    Returns:
        tuple: A tuple containing three elements:
            - all_interpolated_lats (list): List of interpolated latitude values.
            - all_interpolated_lons (list): List of interpolated longitude values.
            - line_ids (list): List of None values with the same length as all_interpolated_lats.
    """
    try:
        coordinates = stops_df[[stops_lon_key, stops_lat_key]].values.tolist()
        all_interpolated_lats, all_interpolated_lons = interpolate_linestring(coordinates, distance_m)
        line_ids = stops_df[[stops_id_key]].values.tolist()
        print(f'[interpolate.py]length of all_interpolated_lats: {len(all_interpolated_lats)}', f'length of all_interpolated_lons: {len(all_interpolated_lons)}', f'length of line_ids: {len(line_ids)}')
        return all_interpolated_lats, all_interpolated_lons, line_ids
    except KeyError as e:
        raise ValueError(f"[interpolate.py]Missing key in stops DataFrame: {e}")
    except Exception as e:
        raise RuntimeError(f"[interpolate.py]An error occurred while processing the stops DataFrame: {e}")