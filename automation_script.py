from typing import Dict
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress
import numpy as np

from graph_manager import *

# Constants
pellet_file_id_index_start = 12
pellet_file_id_index_end = 13

time_column_title = 'Time(s)'
weight_column_title = 'Weight'

max_time_plot_s = 60*60
experiment_duration = 100*60
plateau_time = 3*60
hematite_oxygen_pct = 0.300564
hematite_iron_pct = 0.699436

def get_file_data():
    """
    - Read the data files from local directory
    - Format data for Iron Reduction percentage
    returns Dict[FileName:str, Data:DataFrame]
    """
    # Get the data from files
    dir_current = os.path.dirname(os.path.abspath(__file__))
    experiment_data_files:Dict[str, pd.DataFrame] = read_data_files(dir_current)
    return experiment_data_files
    
class LinregressRange:
    def __init__(self, min, max) -> None:
        self.min = min
        self.max = max

class Pellet:
    def __init__(self, initial_mass:float, start_time_s:int, color:str, label:str, iron_content_XRD:float,initial_radius:float,final_mass:float, final_iron_content:float, S_w_to_iron:float, S_m_to_w:float):
        iron_content =0.699436
        oxygen_content = 0.300564
        gangue_content = 0.04
        self.initial_mass = initial_mass
        self.start_time_s = start_time_s
        self.iron = initial_mass*(1-gangue_content)*iron_content
        self.oxy = initial_mass*(1-gangue_content)*oxygen_content
        self.color = color
        self.label = label
        self.iron_content_XRD = iron_content_XRD
        self.initial_radius = initial_radius
        self.final_mass = final_mass
        self.final_iron_content = final_iron_content
        self.S_w_to_iron = S_w_to_iron
        self.S_m_to_w = S_m_to_w

    def get_linregress_y_values(self, x, y):

        # Compute linregress
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        slope:float
        intercept:float

        # Predict Y values based on the linear fit
        predicted_Y = slope * x + intercept

        return predicted_Y, slope, intercept
    
    def plot(self, graph_config:GraphConfig, graph_time: Series, linregress_ranges:list[LinregressRange] = []):
        self._plot(
            graph=graph_config.graph,
            time_data_s=graph_time.iloc[:len(graph_config.y_values())],
            data=graph_config.y_values(),
            linregress_ranges=linregress_ranges
        )

    def _plot(self, graph, time_data_s, data, linregress_ranges:list[LinregressRange] = []):
        # copy label for modifications
        label = self.label

        # do linregress for each range
        for linrange in linregress_ranges:
            filtered_time = time_data_s[(data >= linrange.min) & (data <= linrange.max)]
            filtered_data = data[(data >= linrange.min) & (data <= linrange.max)]
            predicted_y, slope, intercept = self.get_linregress_y_values(
                x=filtered_time,
                y=filtered_data
            )
            # plot linregress
            self.lines, = graph.plot(
                filtered_time / 60, 
                predicted_y, 
                color='black',
            )
            # label += f"\nRange:[{linrange.min:.0f} - {linrange.max:.0f}] R = {slope:.2f}*t + {intercept:.2f}"
            label += f"\nRange:[{linrange.min:.0f} - {linrange.max:.0f}] Slope = {slope:.4f}"

        # plot curve
        graph.plot(
            time_data_s / 60, 
            data, 
            color=self.color, 
            label=label,
            alpha=1
        )

# returns a list of x values closest to respective y_values_to_find values
def get_nearest_x_for_y(y_values_to_find, y_list: Series, x_list):
    result = []
    for y_value in y_values_to_find:
        closest_index = (y_list - y_value).abs().idxmin()
        closest_x_value = x_list[closest_index]
        result.append(closest_x_value)
    return result

# returns a list of x values closest to respective y_values_to_find values
def get_nearest_y_index_from_value(y_value:float, y_list: Series):
    return int((y_list - y_value).abs().idxmin())

def convert_to_pellet_config(pellet_number:int) -> Pellet:
    #Information Pellet 1
    pellet_1 = Pellet(
        initial_mass=9.695,
        start_time_s=4670,
        color='blue',
        label='D=17.1mm',
        iron_content_XRD=0.992,
        initial_radius = 0.0085,
        final_mass= 7.0,
        final_iron_content = 0.928 + (1/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 30,
        S_m_to_w = 30
    )

    #Information Pellet 2
    pellet_2 = Pellet(
        initial_mass=9.2991,
        start_time_s=4483,
        color='firebrick',
        label='D=17.2mm',
        iron_content_XRD=0.991,
        initial_radius = 0.0086,
        final_mass= 6.7,
        final_iron_content = 0.925 + (1.1/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 30,
        S_m_to_w = 30
    )

    #Information Pellet 3
    pellet_3 = Pellet(
        initial_mass=8.057,
        start_time_s=4612,
        color='lightseagreen',
        label='D=16.1mm',
        iron_content_XRD=0.994,
        initial_radius = 0.0080,
        final_mass= 5.8,
        final_iron_content = 0.928 + (0.7/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 30,
        S_m_to_w = 30
    )
    #Information Pellet 5
    pellet_5 = Pellet(
        initial_mass=5.423,
        start_time_s=4550,
        color='red',
        label='D=13.9mm',
        iron_content_XRD=1,
        initial_radius = 0.00695,
        final_mass= 3.5,
        final_iron_content = 0.936 + (1.6/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 30,
        S_m_to_w = 30
    )

    #Information Pellet 6
    pellet_6 = Pellet(
        initial_mass=4.8950,
        start_time_s=4552,
        color='magenta',
        label='D=13.4mm',
        iron_content_XRD=0.987,
        initial_radius = 0.0067,
        final_mass= 3.5,
        final_iron_content = 0.936 + (1.6/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 30,
        S_m_to_w = 30
    )

    #Information Pellet 7
    pellet_7 = Pellet(
        initial_mass=4.1002,
        start_time_s=4493,
        color='orange',
        label='D=12.8mm',
        iron_content_XRD=0.998,
        initial_radius = 0.0064,
        final_mass= 2.9,
        final_iron_content = 0.952 + (0.2/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 30,
        S_m_to_w = 30
    )

    #Information Pellet 9
    pellet_9 = Pellet(
        initial_mass=1.6551,
        start_time_s=4580,
        color='green',
        label='D=9.5mm',
        iron_content_XRD=0.996,
        initial_radius = 0.00475,
        final_mass= 1.2,
        final_iron_content = 0.953 + (0.5/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 50,
        S_m_to_w = 31
    )

    #Information Pellet 10
    pellet_10 = Pellet(
        initial_mass=1.9795,
        start_time_s=4660,
        color='purple',
        label='D=9.8mm',
        iron_content_XRD=1,
        initial_radius = 0.0049,
        final_mass= 1.4,
        final_iron_content = 0.961,
        S_w_to_iron = 50,
        S_m_to_w = 31
    )

    #Information Pellet 12
    pellet_12 = Pellet(
        initial_mass=1.3451,
        start_time_s=4573,
        color='steelblue',
        label='D=8.7mm',
        iron_content_XRD=1,
        initial_radius = 0.00437,
        final_mass= 1.0,
        final_iron_content = 0.94,
        S_w_to_iron = 50,
        S_m_to_w = 31
    )

    #Information Pellet 13 whichis the first Keetac Pellet reduced attracted do magnet and broke after reduction
    pellet_13 = Pellet(
        initial_mass=1.3451,
        start_time_s=4581,
        color='navy',
        label='D=12.4mm (II)',
        iron_content_XRD=0.99,
        initial_radius = 0.0062,
        final_mass= 3.072,
        final_iron_content = 0.976 + (0.8/100)*(0.95*55.85/(0.95*55.85+16))+ (0.7/100)*(3*55.85/(3*55.85+4*16)),
        S_w_to_iron = 50,
        S_m_to_w = 30
    )

    #Information Pellet 14 whichis the first Keetac Pellet reduced attracted do magnet and broke after reduction
    pellet_14 = Pellet(
        initial_mass=1.3451,
        start_time_s=4587,
        color='darkslategrey',
        label='D=12.5mm (II)',
        iron_content_XRD=0.998,
        initial_radius = 0.0063,
        final_mass= 3.132,
        final_iron_content = 0.985 + (0.4/100)*(0.95*55.85/(0.95*55.85+16)),
        S_w_to_iron = 50,
        S_m_to_w = 30
    )
    switcher = {
        1: pellet_1,
        2: pellet_2,
        3: pellet_3,
        5: pellet_5,
        6: pellet_6,
        7: pellet_7,
        9: pellet_9,
        10: pellet_10,
        12: pellet_12,
        13: pellet_13,
        14: pellet_14
    }
    
    if pellet_number not in switcher:
        raise ValueError(f"No file found for Pellet number: {pellet_number}")
 
    # get() method of dictionary data type returns value of passed argument if it is present in dictionary otherwise second argument will be assigned as default value of passed argument
    return switcher.get(pellet_number) # type: ignore

def convert_file_name_to_pellet_config(file_name:str) -> Pellet:
    pellet_number = int(file_name[pellet_file_id_index_start:pellet_file_id_index_end+1])
    return convert_to_pellet_config(pellet_number)


def read_data_files(path: str) -> Dict[str, pd.DataFrame]:
    """
    TXT files are expected to be in the following format -> [Time:Weight]

    Returns: Dictionary with {File Name:DataFrame}
    """
    data_mass_dict: Dict[str, pd.DataFrame] = {}  # dict with file name as keys and DataFrame as values

    # Lists all files in the directory
    files_in_folder = os.listdir(path)

    # Iterates over the files in the folder
    for file_name in files_in_folder:
        # Check if the file is a CSV or TXT file
        if file_name.endswith('.csv') or file_name.endswith('.txt'):
            # Assemble the full file path
            full_path = os.path.join(path, file_name)

            if file_name.endswith('.csv'):
                # # Reads the file skipping the first 8 lines and taking ';' as a separator (adjust as necessary)
                # data_temp.append(pd.read_csv(caminho_completo, skiprows=8, on_bad_lines='skip'))
                pass
            else:
                file_name_without_txt = os.path.splitext(file_name)[0]
                data_mass_dict[file_name_without_txt] = pd.read_csv(full_path, sep='\t')

    return data_mass_dict

def format_file_data(file_data:DataFrame, pellet:Pellet) -> DataFrame:
    """
    file_data should be a DataFrame with the following format -> [Time ; Weight]

    Returns -> Data with excluded Initial Values (pellet.start_time_s) and Final Values (experiment_duration)

    Data has the following collumns:

    [time_column_title ; weight_column_title]
    """
    # Use first two columns -> Time : Weight
    formatted_data:DataFrame = file_data.iloc[:, :2]
    formatted_data.columns = [time_column_title,weight_column_title]

    # exclude initial values and set the new initial time as zero
    formatted_data = formatted_data.iloc[pellet.start_time_s:]
    formatted_data[time_column_title] -= pellet.start_time_s
    formatted_data[weight_column_title] -= formatted_data[weight_column_title][pellet.start_time_s]

    # exclude final values
    return formatted_data.iloc[:experiment_duration]
        

# def plot_data_points(experiment_data_files:DataFrame):
    
#     F_values_to_find = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
#     print(' '.join([' '] + [str(cell) for cell in F_values_to_find]))
#     for file_name, file_data in experiment_data_files.items():
#         pellet_config = convert_file_name_to_pellet_config(file_name)
#         format_file_data(file_data, pellet_config)

#         # Iterate over all graphs in Graphs
#         for graph_config in graph_configurations.get_all_graphs():
#             pellet_config.plot(graph_config)
#         reduct_workbook = ReductionWorkbook(pellet_config.initial_radius, pellet_config.initial_mass)
        
#         params_model_t_pore, _ = curve_fit(
#             reduct_workbook.model_t_pore, 
#             np.clip(graph_configurations.reduct_graph_config.y_values().iloc[:max_time_plot_s].to_numpy(),0,1),
#             graph_equations.time.iloc[:max_time_plot_s].to_numpy(), 
#             p0=[1],
#             bounds=(0,5)
#         )

#         params_model_t_external, _ = curve_fit(
#             reduct_workbook.model_t_external, 
#             np.clip(graph_configurations.reduct_graph_config.y_values().iloc[:max_time_plot_s].to_numpy(),0,1),
#             graph_equations.time.iloc[:max_time_plot_s].to_numpy(), 
#             p0=[1],
#             bounds=(0,10)
#         )

#         params_model_t_mixed, _ = curve_fit(
#             reduct_workbook.model_t_mixed, 
#             np.clip(graph_configurations.reduct_graph_config.y_values().iloc[:max_time_plot_s].to_numpy(),0,1),
#             graph_equations.time.iloc[:max_time_plot_s].to_numpy(), 
#             p0=[1],
#             bounds=(0,5)
#         )
        
#         pellet_config.plot_model(reduct_workbook.model_t_pore, params_model_t_pore)
#         pellet_config.plot_model(reduct_workbook.model_t_external, params_model_t_external)
#         pellet_config.plot_model(reduct_workbook.model_t_mixed, params_model_t_mixed)
#         pass