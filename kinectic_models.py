from typing import Dict
from matplotlib import pyplot as plt
import pandas as pd
import automation_script
from graph_manager import GraphConfig, GraphEquations

graph_equations = GraphEquations()

class GraphConfigKinectics:
    def getGraphConfig(y_label:str, y_values_callable:lambda:pd.Series, label_font_size:int=20,tick_font_size:int=20, label_color:str="black", time_label:str='Time (min)', legend_loc:str='lower right'):
        _, graph = plt.subplots(figsize=(10, 6))
        return GraphConfig(graph, y_label, y_values_callable, label_font_size,tick_font_size, label_color, time_label, legend_loc)

class GraphsConfigurations:
    def __init__(self):
        self.reduct_graph_config = GraphConfigKinectics.getGraphConfig('F',y_values_callable=lambda:graph_equations.reduction)
        self.iron_layer_limiting_graph_config = GraphConfigKinectics.getGraphConfig('$\\frac{1}{2}-\\frac{1}{3}F-\\frac{1}{2}(1-F)^{\\frac{2}{3}}$',y_values_callable=lambda:graph_equations.iron_layer_limiting)
        self.mixed_control_limiting_graph_config = GraphConfigKinectics.getGraphConfig('$1-(1-F)^\\frac{1}{3}$',y_values_callable=lambda:graph_equations.limiting_mixed_control)
        self.complete_internal_burning_graph_config = GraphConfigKinectics.getGraphConfig('$ln(1-F)$',y_values_callable=lambda:graph_equations.complete_internal_burning)
        self.joint_emt_plus_sc_graph_config = GraphConfigKinectics.getGraphConfig('$F$',y_values_callable=lambda:graph_equations.joint_emt_plus_sc)

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

def create_graphs(experiment_data_files: Dict[str, pd.DataFrame]):

    plot_data_points(experiment_data_files)
    graph_configurations.set_legends()

def plot_data_points(experiment_data_files:pd.DataFrame):
    
    F_values_to_find = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    print(' '.join([' '] + [str(cell) for cell in F_values_to_find]))
    for file_name, file_data in experiment_data_files.items():
        pellet_config = automation_script.convert_file_name_to_pellet_config(file_name)
        formatted_data = automation_script.format_file_data(file_data, pellet_config)
        graph_equations.calculate(time=formatted_data[automation_script.time_column_title], weight=formatted_data[automation_script.weight_column_title], pellet=pellet_config, plateau_time=automation_script.plateau_time)

        # Iterate over all graphs in Graphs
        for graph_config in graph_configurations.get_all_graphs():
            pellet_config.plot(graph_config, graph_equations.time)


create_graphs(automation_script.get_file_data())
plt.show(block=True)
