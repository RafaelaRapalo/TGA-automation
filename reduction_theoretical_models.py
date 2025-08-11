from typing import Dict
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import automation_script
from graph_manager import GraphConfig, GraphEquations
from reduction_workbook import ReductionWorkbook
from scipy.optimize import curve_fit

graph_equations = GraphEquations()
F_values_to_find = [0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]

def create_graphs(experiment_data_files: Dict[str, pd.DataFrame]):

    plot_data_points(experiment_data_files)

def calculate_r_squared(y_data, y_pred):
    residuals = y_data- y_pred
    residual_sum_of_squares = np.sum(residuals**2)
    total_sum_of_squares = np.sum((y_data-np.mean(y_data))**2)
    r_squared = 1 - (residual_sum_of_squares/total_sum_of_squares)
    return r_squared

def plot_data_points(experiment_data_files: Dict[str, pd.DataFrame]):
    
    print(','.join([str(cell) for cell in F_values_to_find]))
    
    for file_name, file_data in experiment_data_files.items():
        pellet_config = automation_script.convert_file_name_to_pellet_config(file_name)
        formatted_data = automation_script.format_file_data(file_data, pellet_config)
        graph_equations.calculate(time=formatted_data[automation_script.time_column_title], weight=formatted_data[automation_script.weight_column_title], pellet=pellet_config, plateau_time=automation_script.plateau_time)
        
        _, graph = plt.subplots(figsize=(10, 6))
        graph_config = GraphConfig(graph,'F',y_values_callable=lambda:graph_equations.reduction)
        pellet_config.plot(
            graph_config= graph_config,
            graph_time=graph_equations.time
        )
        reduct_workbook = ReductionWorkbook(pellet_config)


        last_index_for_fitting_F = automation_script.get_nearest_y_index_from_value(
            y_list=graph_config.y_values(),
            y_value=0.95
        )
        limited_F_values = np.clip(graph_config.y_values().iloc[:last_index_for_fitting_F-pellet_config.start_time_s].to_numpy(),0,1)
        curve_fit_time_values = graph_equations.time.iloc[:last_index_for_fitting_F-pellet_config.start_time_s].to_numpy()
        params_model_t_pore, _ = curve_fit(
            reduct_workbook.model_t_pore, 
            limited_F_values,
            curve_fit_time_values, 
            p0=[1],
            bounds=(0,5)
        )

        params_model_t_external, _ = curve_fit(
            reduct_workbook.model_t_external, 
            limited_F_values,
            curve_fit_time_values, 
            p0=[1],
            bounds=(0,10)
        )

        params_model_t_mixed, _ = curve_fit(
            reduct_workbook.model_t_mixed, 
            limited_F_values,
            curve_fit_time_values, 
            p0=[1],
            bounds=(0,5)
        )

        #r squared fitting
        r_squared_model_t_pore = calculate_r_squared(curve_fit_time_values, reduct_workbook.model_t_pore(limited_F_values, *params_model_t_pore))
        r_squared_model_t_external = calculate_r_squared(curve_fit_time_values, reduct_workbook.model_t_external(limited_F_values, *params_model_t_external))
        r_squared_model_t_mixed = calculate_r_squared(curve_fit_time_values, reduct_workbook.model_t_mixed(limited_F_values, *params_model_t_mixed))
        
        graph.plot(
            reduct_workbook.model_t_mixed(limited_F_values, *params_model_t_mixed) / 60, 
            limited_F_values, 
            color='black', 
            label=f'Limiting mixed control\nS·De = {params_model_t_mixed[0]*pow(10,2):.4f} $·10^{{-2}} m^4/mol·s$\n$R^2$ = {r_squared_model_t_mixed:.4f}',
            alpha=1,
            ls='-.'
        )

        graph.plot(
            reduct_workbook.model_t_pore(limited_F_values, *params_model_t_pore) / 60, 
            limited_F_values, 
            color='black', 
            label=f'Diffusion in porous iron layer\nDe\' = {params_model_t_pore[0]*pow(10,4):.4f} $·10^{{-4}} m^2/s$\n$R^2$ = {r_squared_model_t_pore:.4f}',
            alpha=1,
            ls=':'
        )
        graph.plot(
            reduct_workbook.model_t_external(limited_F_values, *params_model_t_external) / 60, 
            limited_F_values, 
            color='black', 
            label=f'External mass transfer\nm = {params_model_t_external[0]:.4f} $m/s$\n$R^2$ = {r_squared_model_t_external:.4f}',
            alpha=1,
            ls='--'
        )

        graph_config.set_legend_loc()
        pass


create_graphs(automation_script.get_file_data())
plt.show(block=True)
