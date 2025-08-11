import numpy as np
import pandas as pd
import os
from scipy.optimize import curve_fit
from automation_script import Pellet
from reduction_workbook import ReductionWorkbook
from matplotlib import pyplot as plt

paper_pellet_params = Pellet(
    initial_mass=6,
    start_time_s=4550,
    color='red',
    label='Heidari 2025',
    iron_content_XRD=1,
    initial_radius = 0.0135/2,
    final_mass= 6*0.7,
    final_iron_content = 0.936 + 1.6*(0.95*55.85/(0.95*55.85*16))
)

dir_current = os.path.dirname(os.path.abspath(__file__))
file_path = dir_current + '\paper_comparison_reduction_data\heidari_paper_pellet_01.xlsx'
# Read the Excel file
df = pd.read_excel(file_path)
f_values = df['F']
time_values = df['Time (s)']

# Filter data where reduction ≤ 0.95
mask = f_values <= 0.95
limited_f_values = f_values[mask]
limited_t_values = time_values[:len(limited_f_values)]

def calculate_r_squared(y_data, y_pred):
    residuals = y_data- y_pred
    residual_sum_of_squares = np.sum(residuals**2)
    total_sum_of_squares = np.sum((y_data-np.mean(y_data))**2)
    r_squared = 1 - (residual_sum_of_squares/total_sum_of_squares)
    return r_squared

reduct_workbook = ReductionWorkbook(paper_pellet_params)

params_model_t_pore, _ = curve_fit(
    reduct_workbook.model_t_pore, 
    limited_f_values,
    limited_t_values, 
    p0=[1],
    bounds=(0,5)
)

params_model_t_external, _ = curve_fit(
    reduct_workbook.model_t_external, 
    limited_f_values,
    limited_t_values, 
    p0=[1],
    bounds=(0,10)
)

params_model_t_mixed, _ = curve_fit(
    reduct_workbook.model_t_mixed, 
    limited_f_values,
    limited_t_values, 
    p0=[1],
    bounds=(0,5)
)


#r squared fitting
r_squared_model_t_pore = calculate_r_squared(limited_t_values, reduct_workbook.model_t_pore(limited_f_values, *params_model_t_pore))
r_squared_model_t_external = calculate_r_squared(limited_t_values, reduct_workbook.model_t_external(limited_f_values, *params_model_t_external))
r_squared_model_t_mixed = calculate_r_squared(limited_t_values, reduct_workbook.model_t_mixed(limited_f_values, *params_model_t_mixed))

plt.plot(
    (time_values) / 60, 
    f_values, 
    color='red', 
    label=f'Heidari 2025',
    alpha=1,
)

plt.plot(
    reduct_workbook.model_t_pore(limited_f_values, *params_model_t_pore) / 60, 
    limited_f_values, 
    color='black', 
    label=f'Diffusion in porous iron layer\nDe = {params_model_t_pore[0]*10000:.4f} $·10^{{-4}}m^2/s$\n$R^2$ = {r_squared_model_t_pore:.2f}',
    alpha=1,
    ls=':'
)
plt.plot(
    reduct_workbook.model_t_external(limited_f_values, *params_model_t_external) / 60, 
    limited_f_values, 
    color='black', 
    label=f'External mass transfer\nm = {params_model_t_external[0]:.4f} $m/s$\n$R^2$ = {r_squared_model_t_external:.4f}',
    alpha=1,
    ls='--'
)
plt.plot(
    reduct_workbook.model_t_mixed(limited_f_values, *params_model_t_mixed) / 60, 
    limited_f_values, 
    color='black', 
    label=f'Limiting mixed control\n$De·S$ = {params_model_t_mixed[0]*100:.4f}$·10^{{-2}}m^4/(s·mol)$\n$R^2$ = {r_squared_model_t_mixed:.4f}',
    alpha=1,
    ls='-.'
)

plt.xlabel('Time (s)')
plt.ylabel('Degree of Reduction')
plt.legend(loc='lower right')
plt.show()