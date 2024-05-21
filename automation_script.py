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

max_time_plot_s = 25*60
experiment_duration = 100*60
plateau_time = 3*60
hematite_oxygen_pct = 0.300564
hematite_iron_pct = 0.699436

graph_equations = GraphEquations()

def __main__():
    """
    - Read the data files from local directory
    - Format data for Iron Reduction percentage
    - Plot graphs
    """
    # Get the data from files
    dir_current = os.path.dirname(os.path.abspath(__file__))
    experiment_data_files:Dict[str, pd.DataFrame] = read_data_files(dir_current)
    create_graphs(experiment_data_files)
    plt.show()

class GraphsConfigurations:
    def __init__(self):
        self.reduct_graph_config = GraphConfig(plt,'F',y_values_callable=lambda:graph_equations.reduction)
        self.iron_layer_limiting_graph_config = GraphConfig(plt,'$\\frac{1}{2}-\\frac{1}{3}F-\\frac{1}{2}(1-F)^{\\frac{2}{3}}$',y_values_callable=lambda:graph_equations.iron_layer_limiting)
        self.mixed_control_limiting_graph_config = GraphConfig(plt,'$1-(1-F)^\\frac{1}{3}$',y_values_callable=lambda:graph_equations.limiting_mixed_control)
        self.complete_internal_burning_graph_config = GraphConfig(plt,'$ln(1-F)$',y_values_callable=lambda:graph_equations.complete_internal_burning)
        self.external_mass_transfer_graph_config = GraphConfig(plt,'$F$',y_values_callable=lambda:graph_equations.external_mass_transfer)

    def get_all_graphs(self):
        all_graphs:list[GraphConfig] = []
        for _, graph_config in vars(self).items():
            if isinstance(graph_config,GraphConfig):
                all_graphs.append(graph_config)
        return all_graphs
    
    def set_legends(self):
        for graph_config in self.get_all_graphs():
            graph_config.set_legend_loc()

graph_configurations = GraphsConfigurations()
    
class LinregressRange:
    def __init__(self, min, max) -> None:
        self.min = min
        self.max = max

class Pellet:
    def __init__(self, initial_mass:float, start_time_s:int, color:str, label:str, iron_content_XRD:float,initial_radius:float):
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

    def get_linregress_y_values(self, x, y):

        # Compute linregress
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        slope:float
        intercept:float

        # Predict Y values based on the linear fit
        predicted_Y = slope * x + intercept

        return predicted_Y, slope, intercept
    
    def plot(self, graph_config:GraphConfig, linregress_ranges:list[LinregressRange] = []):
        self._plot(
            graph=graph_config.graph,
            time_data_s=graph_equations.time.iloc[:max_time_plot_s],
            data=graph_config.y_values().iloc[:max_time_plot_s],
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
            graph.plot(
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

def convert_to_pellet_config(pellet_number:int) -> Pellet:
    #Information Pellet 1
    pellet_1 = Pellet(
        initial_mass=9.695,
        start_time_s=4670,
        color='blue',
        label='D=17.1mm',
        iron_content_XRD=1,
        initial_radius = 0.0085
    )

    #Information Pellet 2
    pellet_2 = Pellet(
        initial_mass=9.2991,
        start_time_s=4483,
        color='firebrick',
        label='D=17.2mm',
        iron_content_XRD=1,
        initial_radius = 0.0086
    )

    #Information Pellet 5
    pellet_5 = Pellet(
        initial_mass=5.423,
        start_time_s=4550,
        color='red',
        label='D=13.9mm',
        iron_content_XRD=1,
        initial_radius = 0.00695
    )

    #Information Pellet 6
    pellet_6 = Pellet(
        initial_mass=4.8950,
        start_time_s=4552,
        color='magenta',
        label='D=13.4mm',
        iron_content_XRD=1,
        initial_radius = 0.0067
    )

    #Information Pellet 7
    pellet_7 = Pellet(
        initial_mass=4.1002,
        start_time_s=4493,
        color='orange',
        label='D=12.8mm',
        iron_content_XRD=1,
        initial_radius = 0.0064
    )

    #Information Pellet 9
    pellet_9 = Pellet(
        initial_mass=1.6551,
        start_time_s=4580,
        color='green',
        label='D=9.5mm',
        iron_content_XRD=1,
        initial_radius = 0.00475
    )

    #Information Pellet 10
    pellet_10 = Pellet(
        initial_mass=1.9795,
        start_time_s=4660,
        color='purple',
        label='D=9.8mm',
        iron_content_XRD=1,
        initial_radius = 0.0049
    )

    switcher = {
        1: pellet_1,
        2: pellet_2,
        5: pellet_5,
        6: pellet_6,
        7: pellet_7,
        9: pellet_9,
        10: pellet_10,
    }
    if pellet_number not in switcher:
        raise ValueError(f"No file found for Pellet number: {pellet_number}")
 
    # get() method of dictionary data type returns value of passed argument if it is present in dictionary otherwise second argument will be assigned as default value of passed argument
    return switcher.get(pellet_number) # type: ignore

def convert_file_name_to_pellet_config(file_name:str) -> Pellet:
    pellet_number = int(file_name[pellet_file_id_index_start:pellet_file_id_index_end+1])
    return convert_to_pellet_config(pellet_number)

# helper methods
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

def format_file_data(file_data:DataFrame, pellet:Pellet) -> GraphEquations:
    """
    file_data should be a DataFrame with the following format -> [Time ; Weight]

    Returns -> Formatted data with the following collumns:

    [time_column_title ; weight_column_title ; oxygen_in_pellet_title ; 
    reduction_title ; reduction_pct_title ; iron_layer_limiting_title ;
    limiting_mixed_control_title ; complete_internal_burning_title]
    """
    # Use first two columns -> Time : Weight
    formatted_data:DataFrame = file_data.iloc[:, :2]
    formatted_data.columns = [time_column_title,weight_column_title]

        # exclude initial values and set the new initial time as zero
    formatted_data = formatted_data.iloc[pellet.start_time_s:]
    formatted_data[time_column_title] -= pellet.start_time_s

        # exclude final values
    formatted_data = formatted_data.iloc[:experiment_duration]
        
    graph_equations.calculate(time=formatted_data[time_column_title], weight=formatted_data[weight_column_title], pellet=pellet, plateau_time=plateau_time)

def create_graphs(experiment_data_files: Dict[str, pd.DataFrame]):

    plot_data_points(experiment_data_files)
    graph_configurations.set_legends()

def plot_data_points(experiment_data_files:DataFrame):
    
    for file_name, file_data in experiment_data_files.items():
        pellet_config = convert_file_name_to_pellet_config(file_name)
        format_file_data(file_data, pellet_config)

        # Iterate over all graphs in Graphs
        for graph_config in graph_configurations.get_all_graphs():
            pellet_config.plot(graph_config)
        pellet_config._plot(
            graph=graph_configurations.reduct_graph_config.graph,
            data=graph_equations.external_mass_transfer.iloc[:max_time_plot_s],
            time_data_s=graph_equations.time.iloc[:max_time_plot_s],
        )

__main__()