from matplotlib import pyplot as plt
import numpy as np
from pandas import DataFrame, Series
from automation_script import Pellet

class GraphEquations:
    def __init__(self):
        pass
    def __init__(self, time:Series, weight:Series, pellet:Pellet, plateau_time:int):
        self.calculate(time, weight, pellet, plateau_time)
    def calculate(self, time:Series, weight:Series, pellet:Pellet, plateau_time:int):
        plateau_value = weight.tail(plateau_time).median()
        
        # Oxygen percentage in pellet = (weight + oxygen)/(weight + oxygen + iron)
        self.time = time
        self.oxygen_in_pellet = (weight + pellet.oxy) / (weight + pellet.oxy + pellet.iron)
        self.reduction = (weight/plateau_value) * pellet.iron_content_XRD
        self.reduction_pct = self.reduction * 100
        self.iron_layer_limiting = 1/2-(1/3)*self.reduction-(1/2)*pow(1-self.reduction,2/3)
        self.limiting_mixed_control = 1-pow(1-self.reduction,1/3)
        self.complete_internal_burning = np.log(1-self.reduction)
        self.external_mass_transfer = 3.7594*0.0008155/(2*pellet.initial_radius)*4*np.pi*pow(pellet.initial_radius,2)*(11.04921-4.40365)*time/(pellet.initial_mass/(55.85+1.5*16))
        self.external_mass_transfer = np.where(self.external_mass_transfer > 1, 1, self.external_mass_transfer)
        self.new_equation = 1

        # # Reduction = ((Ox / Fe)[fe2o3] - (Ox / Fe)[DRI]) / (Ox / Fe)[fe2o3]
        # self.reduction = ((hematite_oxygen_pct / hematite_iron_pct) - (self.oxygen_in_pellet/hematite_iron_pct)) / (hematite_oxygen_pct / hematite_iron_pct)        

        # # self.reduction_pct = self.weight delta * 97%
        # # Reduction Percentage = ((Ox / Fe)[fe2o3] - (Ox / Fe)[DRI]) * 100 / (Ox / Fe)[fe2o3]
        # self.reduction_pct = self.reduction * 100

        # Correcting reduction curves from XRD result
        # plateau_value = self.weight.tail(plateau_time).min()
        # self.weight = self.weight.rolling(window=10).mean()
        
        
class GraphConfig:
    def __init__(self, y_label:str, y_values, label_font_size:int=20,tick_font_size:int=20, label_color:str="black", time_label:str='Time (min)', legend_loc:str='upper right'):
        _,graph = plt.subplots(figsize=(10, 6))

        graph.set_xlabel(time_label, color=label_color, fontsize=label_font_size)
        graph.set_ylabel(y_label, color=label_color, fontsize=label_font_size)
        graph.tick_params(labelsize=tick_font_size)

        graph.legend(loc=legend_loc)
        
        self.graph = graph
        self.y_values = y_values
    