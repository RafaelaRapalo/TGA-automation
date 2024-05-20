from typing import Dict, Tuple
import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress
import numpy as np

# Constants
pellet_file_id_index_start = 12
pellet_file_id_index_end = 13

time_column_title = 'Time(s)'
weight_column_title = 'Weight'
oxygen_in_pellet_title = 'Oxygen_Pct'
reduction_title = 'Reduction'
reduction_pct_title = 'Reduction_Pct'
iron_layer_limiting_title = 'Iron_Layer_LRE'
limiting_mixed_control_title = 'Mixed_Control_LRE'
complete_internal_burning_title = 'Complete_Internal_Burning'
external_mass_transfer_title = 'External_Mass_Transfer_LRE'

max_time_plot_s = 25*60
experiment_duration = 100*60
plateau_time = 3*60
hematite_oxygen_pct = 0.300564
hematite_iron_pct = 0.699436

def __main__():
    """
    - Read the data files from local directory
    - Format data for Iron Reduction percentage
    - Plot graphs
    """
    # Get the data from files
    dir_current = os.path.dirname(os.path.abspath(__file__))
    experiment_data_files:Dict[str, pd.DataFrame] = read_data_files(dir_current)
    graphManager = GraphManager()
    graphManager.create_graphs(experiment_data_files)
    graphManager.show_graphs()

class LinregressRange:
    def __init__(self, min, max) -> None:
        self.min = min
        self.max = max

class Pellet:
    def __init__(self, initial_mass, start_time_s, color:str, label:str, iron_content_XRD:float,initial_radius):
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

    def plot(self, graph, time_data_s:Series, data:Series, linregress_ranges:list[LinregressRange] = []):
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
        
class GraphDataConfig:
    def __init__(self, y_label:str, title:str, x_data, y_data, font_size:int = 20, label_color:str = "black", x_label:str = 'Time (min)') -> None:
        _,self.graph = plt.subplots(figsize=(10, 6))
        self.title = title
        self.graph.set_xlabel(x_label, color=label_color, fontsize=font_size)
        self.graph.set_ylabel(y_label, color=label_color, fontsize=font_size)
        self.graph.tick_params(labelsize=font_size)
        self.graph.legend(loc='lower right')
        self.x_data = x_data
        self.y_data = y_data
        
class FormattedValues:
    def __init__(self, time_data:Series, weight_data:Series, plateau_value:float) -> None:
        self.time_data:Series = time_data
        self.weight_data:Series = weight_data
        self.plateau_value:float = plateau_value

class GraphManager:        
    def show_graphs():
        plt.show()
        
    def create_graphs(experiment_data_files: Dict[str, pd.DataFrame]):

        GraphManager.plot_data_points(experiment_data_files)
        
    def format_file_data(file_data:DataFrame, pellet_config:Pellet) -> FormattedValues:
        """
        file_data should be a DataFrame with the following format -> [Time ; Weight]

        Returns -> FormattedValues with time_data, weight_data and plateau_value
        """
        # Use first two columns -> Time : Weight
        formatted_data:DataFrame = file_data.iloc[:, :2]
        formatted_data.columns = [time_column_title,weight_column_title]

        # exclude initial values and set the new initial time as zero
        formatted_data = formatted_data.iloc[pellet_config.start_time_s:]
        formatted_data[time_column_title] -= pellet_config.start_time_s
        plateau_value = formatted_data[weight_column_title].tail(plateau_time).median()
        
        # exclude final values
        formatted_data = formatted_data.iloc[:experiment_duration]
        
        return FormattedValues(
            time_data=formatted_data[time_column_title],
            weight_data=formatted_data[weight_column_title],
            plateau_value=plateau_value
        )


    def plot_data_points(experiment_data_files:DataFrame):
        
        for file_name, file_data in experiment_data_files.items():
            pellet_config = convert_file_name_to_pellet_config(file_name)
            graph_configs:list[GraphDataConfig] = GraphManager.get_graph_configs(file_data, pellet_config)
            
            for graph_config in graph_configs:
                pellet_config.plot(
                    graph=graph_config.graph,
                    time_data_s=graph_config.x_data.iloc[:max_time_plot_s],
                    data=graph_config.y_data.iloc[:max_time_plot_s]
                )

    def get_graph_configs(self,file_data:Series, pellet_config:Pellet) -> list[GraphDataConfig]:
        formatted_values:FormattedValues = self.format_file_data(file_data, pellet_config)
        
        reduction_values = ((formatted_values.weight_data/formatted_values.plateau_value) * pellet_config.iron_content_XRD)
        iron_layer_limiting_values = 1/2-(1/3)*reduction_values-(1/2)*pow(1-reduction_values,2/3)
        limiting_mixed_control_values = 1-pow(1-reduction_values,1/3)
        complete_internal_burning_values = np.log(1-reduction_values)
        external_mass_transfer_values = 3.7594*0.0008155/(2*pellet_config.initial_radius)*4*np.pi*pow(pellet_config.initial_radius,2)*(11.04921-4.40365)*formatted_values.time_data/(pellet_config.initial_mass/(55.85+1.5*16))
        # # Reduction = ((Ox / Fe)[fe2o3] - (Ox / Fe)[DRI]) / (Ox / Fe)[fe2o3]
        # formatted_data[reduction_title] = ((hematite_oxygen_pct / hematite_iron_pct) - (formatted_data[oxygen_in_pellet_title]/hematite_iron_pct)) / (hematite_oxygen_pct / hematite_iron_pct)
        

        # # formatted_data[reduction_pct_title] = formatted_data[weight_column_title] delta * 97%
        # # Reduction Percentage = ((Ox / Fe)[fe2o3] - (Ox / Fe)[DRI]) * 100 / (Ox / Fe)[fe2o3]
        # formatted_data[reduction_pct_title] = formatted_data[reduction_title] * 100

        # Correcting reduction curves from XRD result
        # plateau_value = formatted_data[weight_column_title].tail(plateau_time).min()
        # formatted_data[weight_column_title] = formatted_data[weight_column_title].rolling(window=10).mean()

        return [
            GraphDataConfig(
                y_label='F',
                title=reduction_title,
                x_data=formatted_values.time_data,
                y_data=reduction_values
            ),
            GraphDataConfig(
                y_label='$\\frac{1}{2}-\\frac{1}{3}F-\\frac{1}{2}(1-F)^{\\frac{2}{3}}$',
                title=iron_layer_limiting_title,
                x_data=formatted_values.time_data,
                y_data=iron_layer_limiting_values
            ),
            GraphDataConfig(
                y_label='$1-(1-F)^\\frac{1}{3}$',
                title=limiting_mixed_control_title,
                x_data=formatted_values.time_data,
                y_data=limiting_mixed_control_values
            ),
            GraphDataConfig(
                y_label='$ln(1-F)$',
                title=complete_internal_burning_title,
                x_data=formatted_values.time_data,
                y_data=complete_internal_burning_values
            ),
            GraphDataConfig(
                y_label='$F$',
                title=external_mass_transfer_title,
                x_data=formatted_values.time_data,
                y_data=external_mass_transfer_values
            ),
        ]


    
def convert_to_pellet_config(pellet_number:int) -> Pellet:

    # switcher for each pellet number
    switcher = {
        1: Pellet(
            initial_mass=9.695,
            start_time_s=4670,
            color='blue',
            label='D=17.1mm',
            iron_content_XRD=1,
            initial_radius = 0.0085
        ),
        2: Pellet(
            initial_mass=9.2991,
            start_time_s=4483,
            color='firebrick',
            label='D=17.2mm',
            iron_content_XRD=1,
            initial_radius = 0.0086
        ),
        5: Pellet(
            initial_mass=5.423,
            start_time_s=4550,
            color='red',
            label='D=13.9mm',
            iron_content_XRD=1,
            initial_radius = 0.00695),
        6: Pellet(
            initial_mass=4.8950,
            start_time_s=4552,
            color='magenta',
            label='D=13.4mm',
            iron_content_XRD=1,
            initial_radius = 0.0067),
        7: Pellet(
            initial_mass=4.1002,
            start_time_s=4493,
            color='orange',
            label='D=12.8mm',
            iron_content_XRD=1,
            initial_radius = 0.0064),
        9: Pellet(
            initial_mass=1.6551,
            start_time_s=4580,
            color='green',
            label='D=9.5mm',
            iron_content_XRD=1,
            initial_radius = 0.00475),
        10: Pellet(
            initial_mass=1.9795,
            start_time_s=4660,
            color='purple',
            label='D=9.8mm',
            iron_content_XRD=1,
            initial_radius = 0.0049),
    }
    
    if pellet_number not in switcher:
        raise ValueError(f"No configuration found for Pellet number: {pellet_number}")
 
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

__main__()