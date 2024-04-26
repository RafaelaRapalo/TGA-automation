from typing import Dict
import pandas as pd
from pandas import DataFrame
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress
import numpy as np

class LinregressRange:
    def __init__(self, min, max) -> None:
        self.min = min
        self.max = max

class Pellet:
    def __init__(self, initial_mass, start_time_s, color:str, label:str, iron_content_XRD:float):
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

    def get_linregress_y_values(self, x, y):

        # Compute linregress
        slope, intercept, r_value, p_value, std_err = linregress(x, y)
        slope:float
        intercept:float

        # Predict Y values based on the linear fit
        predicted_Y = slope * x + intercept

        return predicted_Y, slope, intercept

    def plot(self, graph, time_data_s, data, linregress_ranges:list[LinregressRange] = []):
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

def convert_to_pellet_config(pellet_number) -> Pellet:
    #Information Pellet 1
    pellet_1 = Pellet(
        initial_mass=9.695,
        start_time_s=4670,
        color='blue',
        label='D=17.1mm',
        iron_content_XRD=1
    )

    #Information Pellet 5
    pellet_5 = Pellet(
        initial_mass=5.423,
        start_time_s=4550,
        color='red',
        label='D=13.9mm',
        iron_content_XRD=1
    )

    #Information Pellet 9
    pellet_9 = Pellet(
        initial_mass=1.6551,
        start_time_s=4580,
        color='green',
        label='D=9.5mm',
        iron_content_XRD=1
    )

    #Information Pellet 10
    pellet_10 = Pellet(
        initial_mass=1.9795,
        start_time_s=4660,
        color='purple',
        label='D=9.8mm',
        iron_content_XRD=1
    )

    switcher = {
        1: pellet_1,
        5: pellet_5,
        9: pellet_9,
        10: pellet_10,
    }
 
    # get() method of dictionary data type returns value of passed argument if it is present in dictionary otherwise second argument will be assigned as default value of passed argument
    return switcher.get(pellet_number)

# config variables
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

max_time_plot_s = 25*60
experiment_duration = 100*60
plateau_time = 3*60
hematite_oxygen_pct = 0.300564
hematite_iron_pct = 0.699436

# Get the script's current directory
dir_current = os.path.dirname(os.path.abspath(__file__))

# Lists all files in the directory
files_in_folder = os.listdir(dir_current)

# List to store the DataFrames of each file
data_temp = {}
data_mass_dict: Dict[str,pd.DataFrame] = {} # dict with file name as keys and DataFrame as values

# Iterates over the files in the folder
for file_name in files_in_folder:
    # Check if the file is a CSV or TXT file
    if file_name.endswith('.csv') or file_name.endswith('.txt'):
        # Assemble the full file path
        full_path = os.path.join(dir_current, file_name)

        if file_name.endswith('.csv'):
            # # Reads the file skipping the first 8 lines and taking ';' as a separator (adjust as necessary)
            # data_temp.append(pd.read_csv(caminho_completo, skiprows=8, on_bad_lines='skip'))
            pass
        else:
            file_name_without_txt = os.path.splitext(file_name)[0]
            data_mass_dict[file_name_without_txt] = pd.read_csv(full_path, sep='\t')

###Filtrar

# prep figure
fig1, reduct_graph = plt.subplots(figsize=(10, 6))
fig2, iron_layer_limiting_graph = plt.subplots(figsize=(10, 6))
fig3, mixed_control_limiting_graph = plt.subplots(figsize=(10, 6))
fig4, complete_internal_burning_graph = plt.subplots(figsize=(10, 6))

for file_name, file_data in data_mass_dict.items():
    pellet_number = int(file_name[pellet_file_id_index_start:pellet_file_id_index_end+1])
    pellet_config = convert_to_pellet_config(pellet_number)
    formatted_data:DataFrame = file_data.iloc[:, :2]
    formatted_data.columns = [time_column_title,weight_column_title]

    # exclude initial values and set the new initial time as zero
    formatted_data = formatted_data.iloc[pellet_config.start_time_s:]
    formatted_data[time_column_title] -= pellet_config.start_time_s

    # exclude final values
    formatted_data = formatted_data.iloc[:experiment_duration]
 
    # Oxygen percentage in pellet = (weight + oxygen)/(weight + oxygen + iron)
    formatted_data[oxygen_in_pellet_title] = (formatted_data[weight_column_title] + pellet_config.oxy) / (formatted_data[weight_column_title] + pellet_config.oxy + pellet_config.iron)

    # # Reduction = ((Ox / Fe)[fe2o3] - (Ox / Fe)[DRI]) / (Ox / Fe)[fe2o3]
    # formatted_data[reduction_title] = ((hematite_oxygen_pct / hematite_iron_pct) - (formatted_data[oxygen_in_pellet_title]/hematite_iron_pct)) / (hematite_oxygen_pct / hematite_iron_pct)
    

    # # formatted_data[reduction_pct_title] = formatted_data[weight_column_title] delta * 97%
    # # Reduction Percentage = ((Ox / Fe)[fe2o3] - (Ox / Fe)[DRI]) * 100 / (Ox / Fe)[fe2o3]
    # formatted_data[reduction_pct_title] = formatted_data[reduction_title] * 100

    # Correcting reduction curves from XRD result
    # plateau_value = formatted_data[weight_column_title].tail(plateau_time).min()
    # formatted_data[weight_column_title] = formatted_data[weight_column_title].rolling(window=10).mean()
    plateau_value = formatted_data[weight_column_title].tail(plateau_time).median()
    formatted_data[reduction_title] = (formatted_data[weight_column_title]/plateau_value) * pellet_config.iron_content_XRD
    formatted_data[reduction_pct_title] = formatted_data[reduction_title]*100
        
    
    # Plot reduction percentage curve
    pellet_config.plot(
        graph=reduct_graph,
        time_data_s=formatted_data[time_column_title].iloc[:max_time_plot_s],
        data=formatted_data[reduction_title].iloc[:max_time_plot_s],
        # linregress_ranges=[
        #     LinregressRange(min=40, max=50),
        #     LinregressRange(min=50, max=60),
        #     LinregressRange(min=60, max=70),
        #     LinregressRange(min=70, max=80),
        #     LinregressRange(min=80, max=90),
        #     ]
    )

    # Fit the rate modlesd
    formatted_data[iron_layer_limiting_title] = 1/2-(1/3)*formatted_data[reduction_title]-(1/2)*pow(1-formatted_data[reduction_title],2/3)
    formatted_data[limiting_mixed_control_title] = 1-pow(1-formatted_data[reduction_title],1/3)
    formatted_data[complete_internal_burning_title] = np.log(1-formatted_data[reduction_title])

    # Plot both equations (iron layer limiting rate equation curve + reduction)
    pellet_config.plot(
        graph=iron_layer_limiting_graph,
        time_data_s=formatted_data[time_column_title].iloc[:max_time_plot_s],
        data=formatted_data[iron_layer_limiting_title].iloc[:max_time_plot_s]
    )
    #Plot equation of mixed control limiting rate equation curve
    pellet_config.plot(
        graph=mixed_control_limiting_graph,
        time_data_s=formatted_data[time_column_title].iloc[:max_time_plot_s],
        data=formatted_data[limiting_mixed_control_title].iloc[:max_time_plot_s]
    )

    #Plot equation of mixed control limiting rate equation curve
    pellet_config.plot(
        graph=complete_internal_burning_graph,
        time_data_s=formatted_data[time_column_title].iloc[:max_time_plot_s],
        data=formatted_data[complete_internal_burning_title].iloc[:max_time_plot_s]
    )
    # ax.plot(xlingress/60, predicted_Y, color='gray', label='linear fitting' + pellet_config.label)

#plt.grid(True)
reduct_graph.set_xlabel('Time (min)', color = "black", fontsize = 12)
reduct_graph.set_ylabel('F', color = "black", fontsize = 12)

iron_layer_limiting_graph.set_xlabel('Time (min)', color = "black", fontsize = 12)
iron_layer_limiting_graph.set_ylabel(r'$\\frac{1}{2}-\\frac{1}{3}F-\\frac{1}{2}(1-F)^{\\frac{2}{3}}$', color="black", fontsize=12)

mixed_control_limiting_graph.set_xlabel('Time (min)', color = "black", fontsize = 12)
mixed_control_limiting_graph.set_ylabel(r'$1-(1-F)^\\frac{1}{3}$', color = "black", fontsize = 12)

complete_internal_burning_graph.set_xlabel('Time (min)', color = "black", fontsize = 12)
complete_internal_burning_graph.set_ylabel('$ln(1-F)$', color = "black", fontsize = 12)

# Adicione uma legenda
reduct_graph.legend(loc = 'lower right')
iron_layer_limiting_graph.legend(loc = 'lower right')
mixed_control_limiting_graph.legend(loc = 'lower right')
complete_internal_burning_graph.legend(loc = 'lower right')

plt.show()