import json
import os
from raster import load_convert_raster
from line import load_filter_interpolate_line
from calc import calc_indices, calc_svf, save_svf_result_to_dataframe, convert_to_target_area
from save import plot_save_xy_indices_raster, save_svf_result_to_geojson, plot_svf_results, save_list, save_dataframe, load_list, xy_indices_to_x_y_indices

"""
Main script to run the processing and calculation steps for the project.
"""

def load_settings(filename):
    """
    Load settings from a JSON file.
    Parameters:
        filename (str): The name of the JSON file to load.
    Returns:
        dict: The settings loaded from the JSON file, or None if an error occurs.
    Raises:
        Exception: If there is an error loading the JSON file.
    """
    try:
        file_path = os.path.join(os.path.dirname(__file__), filename)
        with open(file_path, 'r') as file:
            settings = json.load(file)
        return settings
    except Exception as e:
        print(f"[main.py]Error loading settings: {e}")
        return None

def run_processing(settings):
    """
    Processes raster and line data based on the provided settings.
    Parameters:
        settings (dict): A dictionary containing the following keys:
            - 'raster_geotiff_path' (str): Path to the raster GeoTIFF file.
            - 'raster_geotiff_band' (int): Band number to process from the raster.
            - 'raster_print_info' (bool): Flag to print raster information.
            - 'line_input_file_path' (str): Path to the line input file.
            - 'line_input_file_type' (str): Type of the line input file.
            - 'line_filter_key' (str): Key to filter the line data.
            - 'line_filter_value' (str): Value to filter the line data.
            - 'line_distance_m' (float): Distance in meters for line interpolation.
    Returns:
        tuple: A tuple containing:
            - raster: Processed raster data, or None if an error occurred.
            - line: Processed line data, or None if an error occurred.
    """
    # Initialize variables
    array_band, lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids = None, None, None, None, None, None, None, None
    # Process raster and line data
    try:
        array_band, lats_raster, lons_raster, num_cells_x, num_cells_y = load_convert_raster(settings['raster_geotiff_path'], settings['raster_geotiff_band'], settings['raster_print_info'])
        print("[main.py]Raster processing completed.")
    except Exception as e:
        print(f"[main.py]Error processing raster: {e}")
    try:
        all_interpolated_lats, all_interpolated_lons, line_ids = load_filter_interpolate_line(settings['line_input_file_path'], settings['line_input_file_type'], settings['line_filter_key'], settings['line_filter_value'], settings['line_id_key'], settings['line_distance_m'])
        print("[main.py]Line processing completed.")
    except Exception as e:
        print(f"[main.py]Error processing line: {e}")
    return array_band, lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids # Return the processed raster and line data, if available otherwise None

def run_calculations(settings, array_band, lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids, folder_path, obstacles=None):
    """
    Perform a series of calculations including indices calculation, SVF calculation, and saving results.
    Parameters:
        settings (dict): Configuration settings for the calculations.
        array_band (numpy.ndarray): The raster data array.
        lats_raster (numpy.ndarray): Latitude values of the raster.
        lons_raster (numpy.ndarray): Longitude values of the raster.
        num_cells_x (int): Number of cells in the x-direction.
        num_cells_y (int): Number of cells in the y-direction.
        all_interpolated_lats (list): List of interpolated latitude values.
        all_interpolated_lons (list): List of interpolated longitude values.
        line_ids (list): List of line identifiers.
        folder_path (str): Path to the folder where results will be saved.
    Returns:
        tuple: A tuple containing:
            - svf_result (list): List of SVF results.
            - svf_result_coordinates (list): List of SVF result coordinates.
    """
    try:
        # Calculate indices
        x_indices, y_indices, line_ids, xy_indices = calc_indices(lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids)
        # Plot and save xy indices on raster
        if not folder_path:
            folder_path = './output'
            print("[main.py]Output folder not specified. Saving results to default folder 'output'.")
        file_name = os.path.join(folder_path, 'xy_indices_and_raster')
        plot_save_xy_indices_raster(xy_indices, num_cells_x, num_cells_y, array_band, file_name)
        # Calculate SVF
        if obstacles: # Check if obstacles are provided in the settings # Obstacleのloadを書いていない
            print(f"[main.py]Starting SVF calculation...Type of settings: {type(settings)}", f"Type of array_band: {type(array_band)}, Type of xy_indices: {type(xy_indices)}", f"Type of obstacles: {type(settings['obstacles'])}")
            svf_list = calc_svf(settings, array_band, xy_indices, obstacles=settings['obstacles'])
        else:
            print(f"[main.py]Starting SVF calculation...Type of settings: {type(settings)}", f"Type of array_band: {type(array_band)}, Type of xy_indices: {type(xy_indices)}")
            svf_list = calc_svf(settings, array_band, xy_indices, obstacles=None)
        # Save SVF results to text files
        if not folder_path:
            folder_path = './output'
            print("[main.py]Output folder not specified. Saving results to default folder 'output'.")
        save_list(svf_list, os.path.join(folder_path, 'svf_list'))
        save_list(xy_indices, os.path.join(folder_path, 'xy_indices'))
        print("[main.py]SVF calculation completed.")
        return svf_list, xy_indices, x_indices, y_indices # Return the SVF results and coordinates
    except Exception as e:
        raise Warning(f"[main.py]Error during calculations: {e}") # Raise a warning if an error occurs during calculations

def load_calculation_results(folder_path):
    """
    Load existing calculation results from the specified folder.
    This function attempts to load previously saved calculation results from the given folder path.
    It expects to find 'svf_list' and 'xy_indices' files in the specified folder.
    Parameters:
        folder_path (str): The path to the folder containing the calculation result files.
    Returns:
        tuple: A tuple containing the following elements:
            - svf_list (list): The loaded list of SVF (Sky View Factor) values.
            - xy_indices (list): The loaded list of XY indices.
            - x_indices (list): The X indices extracted from the XY indices.
            - y_indices (list): The Y indices extracted from the XY indices.
    Raises:
        Exception: If there is an error loading the calculation results, an exception is raised with an error message.
    """
    try:
        # Load existing calculation results
        svf_list = load_list(os.path.join(folder_path, 'svf_list'))
        xy_indices = load_list(os.path.join(folder_path, 'xy_indices'))
        x_indices, y_indices = xy_indices_to_x_y_indices(xy_indices)
        return svf_list, xy_indices, x_indices, y_indices
    except Exception as e:
        raise Exception(f"[main.py]Error loading existing calculation results: {e}")

def run_saving_output_plotting(svf_list, x_indices, y_indices, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids, folder_path):
    """
    Saves SVF results to various formats and plots the results.
    Parameters:
    svf_list (list): List of SVF values.
    x_indices (list): List of x indices.
    y_indices (list): List of y indices.
    num_cells_x (int): Number of cells in the x direction.
    num_cells_y (int): Number of cells in the y direction.
    all_interpolated_lats (list): List of all interpolated latitudes.
    all_interpolated_lons (list): List of all interpolated longitudes.
    line_ids (list): List of line IDs.
    folder_path (str): Path to the folder where results will be saved.
    Returns:
    tuple: A tuple containing:
        - svf_result (list): List of SVF results.
        - svf_result_coordinates (list): List of SVF result coordinates.
    Raises:
    Exception: If there is an error in saving output or plotting results.
    """
    try:
        # Save SVF results to lists
        svf_result_df, svf_result_coordinates_df = save_svf_result_to_dataframe(svf_list, x_indices, y_indices, all_interpolated_lats, all_interpolated_lons, line_ids)
        # Save SVF results to text files in different formats
        save_dataframe(svf_result_df, os.path.join(folder_path, 'svf_result'))
        save_dataframe(svf_result_coordinates_df, os.path.join(folder_path, 'svf_result_coordinates'))
        print(f"[main.py]SVF results saved. example of svf_result: ", svf_result_df[:5])
        print(f"[main.py]SVF results coordinates saved. example of svf_result_coordinates: ", svf_result_coordinates_df[:5])
        # Save SVF results to GeoJSON
        file_name_geojson = os.path.join(folder_path, 'svf_result_latlon')
        save_svf_result_to_geojson(svf_result_coordinates_df, file_name_geojson)
        # Plot SVF results
        file_name_plot = os.path.join(folder_path, 'svf_result_latlon')
        plot_svf_results(svf_result_df, num_cells_x, num_cells_y, raster=None, file_name=file_name_plot)
        return svf_result_df, svf_result_coordinates_df
    except Exception as e:
        raise Exception(f"[main.py]Error saving output and plotting results: {e}")

def process_data():
    """
    Process data by loading settings, checking for completed calculations, and running necessary processing and calculation steps.
    The function performs the following steps:
    1. Load calculation settings from './input/settings_calc.json'.
    2. Check if an output folder is specified in the settings; if not, use the default './output' folder.
    3. Check if the calculation has already been completed by verifying the existence of 'svf_list.txt' and 'xy_indices.txt' in the output folder.
    4. If the calculation is completed, load the results, save the output, plot the results, convert SVF results to the target area, and remove temporary files.
    5. If the calculation is not completed, load processing settings from './input/settings_input.json' and run the processing steps.
    6. Load calculation settings again and run the calculation steps, save the output, plot the results, and convert SVF results to the target area.
    Raises:
        ValueError: If loading calculation or processing settings fails.
        Exception: If an error occurs during processing or if calculation is needed but not completed.
    """
    try:
        # Load settings for processing and calculations
        settings_processing = load_settings('./input/settings_input.json')
        # Run processing steps
        if settings_processing is not None:
            print("[main.py]Processing settings loaded.")
            print("[main.py]Starting processing...")
            array_band, lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids = run_processing(settings_processing)
            print("[main.py]Processing completed.")
            # Load settings for calculations
            settings_calc = load_settings('./input/settings_calc.json')
            # Run calculation steps
            if settings_calc is not None:
                # Check if an output folder is specified
                if settings_calc['output_folder']:
                    folder_path = settings_calc['output_folder']
                else:
                    folder_path = './output'
                    print("[main.py]Output folder not specified. Saving results to default folder 'output'.")
                # Check if the calculation has already been completed
                if os.path.exists(os.path.join(folder_path, 'svf_list.txt')) and os.path.exists(os.path.join(folder_path, 'xy_indices.txt')):
                    try:
                        print("[main.py]Existing calculation results found.")
                        # Load existing calculation results
                        svf_list, xy_indices, x_indices, y_indices = load_calculation_results(folder_path)
                        print("[main.py]Existing calculation results loaded.")
                    except Exception as e:
                        print(f"[main.py]Error loading existing calculation results: {e}")
                else:
                    try:
                        # Run calculations and save results
                        svf_list, xy_indices, x_indices, y_indices = run_calculations(settings_calc, array_band, lats_raster, lons_raster, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids, folder_path)
                        print("[main.py]Calculations completed.")
                    except Exception as e:
                        print(f"[main.py]Error during calculations: {e}")
                # Saving Output and Plotting
                svf_result_df, svf_result_coordinates_df = run_saving_output_plotting(svf_list, x_indices, y_indices, num_cells_x, num_cells_y, all_interpolated_lats, all_interpolated_lons, line_ids, folder_path)
                # Convert SVF results to target area
                convert_to_target_area(svf_result_coordinates_df, folder_path)
                print("[main.py]Calculation and saving completed.")
            else:
                raise ValueError("[main.py]Failed to load calculation settings.")
        else:
            raise ValueError("[main.py]Failed to load processing settings.")
    except Exception as e:
        print(f"[main.py]An error occurred during processing: {e}")

def main():
    """
    Main function to prompt the user to add files to the 'input' folder and process them.
    Steps:
    1. Prompts the user to add files to the 'input' folder.
    2. Checks if the 'input' folder exists and creates it if it does not.
    3. Waits for the user to press Enter after adding files.
    4. Prompts the user to ensure that the settings JSON files are correctly modified.
    5. Lists the files in the 'input' folder.
    6. If no files are found, prompts the user to add files and try again.
    7. If files are found, lists the files and proceeds with processing them.
    Exceptions:
    - Catches and prints any exceptions that occur during execution.
    """
    print("[main.py]Please add your files to the 'input' folder. Do not erase the 'input' folder.")
    print("[main.py]The 'input' folder should contain raster files, line files, and settings JSON files.")
    print("[main.py]If you start a new calculation, make sure to clear the 'output' folder. Otherwise, the existing results will be loaded.")
    # Check if the 'input' folder exists and create it if it does not exist
    input_folder = '/app/input' # docker build automatically creates the 'input' folder
    # Prompt the user to press Enter after adding files
    input("[main.py]Press Enter after you have added your files (raster files, line files, settings json files) to the 'input' folder...")
    # Prompt the user to check the settings JSON files
    print("[main.py]Make sure that the settings JSON files are correctly modified.")
    try:
        # List the files in the 'input' folder
        files = os.listdir(input_folder)
        if not files:
            print("[main.py]No files found in the 'input' folder. Please add files and try again.")
        else:
            print(f"[main.py]Found {len(files)} file(s) in the 'input' folder:")
            for file in files:
                print(f" - {file}")            
            # Proceed with processing the files
            process_data()
    except OSError as e:
        print(f"[main.py]An error occurred while listing files in the 'input' folder: {e}")
    except Exception as e:
        print(f"[main.py]An error occurred during execution: {e}")

if __name__ == "__main__":
    main()