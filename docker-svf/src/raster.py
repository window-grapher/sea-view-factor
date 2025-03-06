import rasterio
import numpy as np
import zipfile
import os

def extract_geotiff_from_zipfile(zip_path, extract_to='.'):
    """
    Extracts a GeoTIFF file from a ZIP archive.
    Parameters:
    zip_path (str): The file path to the ZIP archive.
    extract_to (str): The directory to extract the files to. Defaults to the current directory.
    Returns:
    str: The file path to the extracted GeoTIFF file.
    """
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            # Extract all the contents of zip file in current directory
            zip_ref.extractall(extract_to)
            # Find the first GeoTIFF file in the extracted contents
            for file_name in zip_ref.namelist():
                if file_name.lower().endswith('.tif') or file_name.lower().endswith('.tiff'):
                    return os.path.join(extract_to, file_name)
        raise FileNotFoundError("[raster.py]No GeoTIFF file found in the ZIP archive.")
    except Exception as e:
        print(f"[raster.py]An error occurred while extracting the GeoTIFF file from the ZIP archive: {e}")
        return None

def geotiff_to_numpy_array(geotiff_path):
    """
    Convert a GeoTIFF file to a NumPy array.
    Parameters:
    geotiff_path (str): The file path to the GeoTIFF file.
    Returns:
    numpy.ndarray: A NumPy array containing the raster data from the GeoTIFF file.
    If the GeoTIFF has multiple bands, the result will be a 3D array.
    """
    try:
        # Open the GeoTIFF file
        with rasterio.open(geotiff_path) as dataset:
            # Read the raster data into a NumPy array
            array = dataset.read() # This reads all bands # The shape of GeoTIFF data is typically (bands, height, width)
            # Get the metadata of the GeoTIFF file
            metadata = dataset.meta
            print(f'[raster.py]geotiff_to_numpy_array: metadata for all bands: {metadata}')
        return array, metadata
    except Exception as e:
        print(f"[raster.py]An error occurred while processing the GeoTIFF file: {e}")
        return None, None

def extract_lat_lon(array_band, metadata):
    """
    Extracts the latitude and longitude coordinates for each pixel in a raster array.
    Parameters:
    array_band (numpy.ndarray): The raster data array. Can be 2D or 3D.
    metadata (dict): Metadata containing spatial information, including the 'transform' key.
    Returns:
    tuple: Two numpy arrays containing the latitude and longitude coordinates for each pixel.
            Returns (lats_raster, lons_raster) if successful, otherwise (None, None) in case of an error.
    Raises:
    ValueError: If the array_band has an unexpected shape.
    Exception: If any other error occurs during the extraction process.
    """
    try:
        # Extract the spatial information from the metadata
        transform = metadata['transform']
        if array_band.ndim == 3:
            height, width = array_band.shape[1], array_band.shape[2]
        elif array_band.ndim == 2:
            height, width = array_band.shape[0], array_band.shape[1]
        else:
            raise ValueError("[raster.py]Unexpected array shape: {}".format(array_band.shape))
        # Get the coordinates of the center of each pixel
        cols, rows = np.meshgrid(np.arange(width), np.arange(height))
        # Transform the pixel coordinates to geographic coordinates
        xs, ys = rasterio.transform.xy(transform, rows, cols)
        # Convert to NumPy arrays
        lons_raster = np.array(xs)
        lats_raster = np.array(ys)
        return lats_raster, lons_raster
    except Exception as e:
        print(f"[raster.py]An error occurred while extracting latitude and longitude: {e}")
        return None, None

def load_convert_raster(geotiff_path, band=1, print_info=False):
    """
    Processes a GeoTIFF file and extracts raster data along with latitude and longitude information.
    Parameters:
        geotiff_path (str): The file path to the GeoTIFF file.
        band (int, optional): The specific band to process if the GeoTIFF contains multiple bands. Defaults to 1.
        print_info (bool, optional): If True, prints information about the array and extracted latitude/longitude. Defaults to True.
    Returns:
        tuple: A tuple containing:
            - array_band (numpy.ndarray): The processed raster data for the specified band.
            - lats (numpy.ndarray): The latitude values corresponding to the raster data.
            - lons (numpy.ndarray): The longitude values corresponding to the raster data.
    """
    try:
        # Extract the GeoTIFF file if it is inside a ZIP archive
        if geotiff_path.lower().endswith('.zip'):
            geotiff_path = extract_geotiff_from_zipfile(geotiff_path)
        # Convert the GeoTIFF file to a 3d NumPy array
        array, metadata = geotiff_to_numpy_array(geotiff_path)
        if array is None or metadata is None:
            raise ValueError("[raster.py]Failed to load GeoTIFF file.")
        # Access specific bands if the array is multi-band:
        array_band = array[band - 1] # Access the specified band (1-indexed)
        # Determine the number of cells in x and y directions
        # The shape of the array is typically (bands, height, width)
        if array_band.ndim == 2:
            num_cells_x, num_cells_y = array_band.shape[1], array_band.shape[0]
        elif array_band.ndim == 3:
            num_cells_x, num_cells_y = array_band.shape[2], array_band.shape[1]
        else:
            raise ValueError("[raster.py]Unexpected array shape: {}".format(array_band.shape))
        if print_info:
            print("[raster.py]Array shape:", array_band.shape)  # Prints the shape of the array (height, width)
            print("[raster.py]num_cells_x:", num_cells_x)
            print("[raster.py]num_cells_y:", num_cells_y)
            print("[raster.py]Array data type:", array_band.dtype)  # Prints the data type of the array
        # Extract latitude and longitude information
        lats_raster, lons_raster = extract_lat_lon(array_band, metadata)
        # Check if the latitude and longitude arrays were successfully extracted
        if lats_raster is None or lons_raster is None:
            raise ValueError("[raster.py]Failed to extract latitude and longitude.")
        if print_info:
            print("[raster.py]Latitude array shape:", lats_raster.shape)  # Shape will match the raster dimensions (height, width)
            print("[raster.py]Longitude array shape:", lons_raster.shape)
            print("[raster.py]Latitude range:", [np.min(lats_raster), np.max(lats_raster)])
            print("[raster.py]Longitude range:", [np.min(lons_raster), np.max(lons_raster)])
        return array_band, lats_raster, lons_raster, num_cells_x, num_cells_y
    except Exception as e:
        print(f"[raster.py]An error occurred while loading and converting the raster: {e}")
        return None, None, None, None, None