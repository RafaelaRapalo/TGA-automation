from typing import Dict
import pandas as pd
from pandas import DataFrame, Series
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress
import numpy as np
from automation_script import Pellet

class ReductionWorkbook:
    def __init__(self, pellet_data:Pellet) -> None:
        self.pellet_radius_m = pellet_data.initial_radius
        #H2-H2O spreadsheet
        X_h2 = 0.9 # E3
        temperature = 900+273 # A12
        tube_diameter_m =0.05
        pellet_diameter_m = self.pellet_radius_m*2
        iron_density = (pellet_data.final_mass * pellet_data.final_iron_content/55.85)/((4/3)*np.pi*pow(self.pellet_radius_m,3))
        self.oxygen_density_hematite = 1.5*iron_density
        self.oxygen_density_wustite = 1.05*iron_density

#external mass transfer
        gas_volume_inlet= 0.000025
        pressure_bar = 1
        molar_mass_h2 = 2
        molar_mass_h2o = 18
        molar_mass_mixture = 2/((1/molar_mass_h2)+(1/molar_mass_h2o))
        X_h2o = 1- X_h2
        molar_mass_average = molar_mass_h2*X_h2 + molar_mass_h2o*X_h2o
        characteristic_lenght_h2_anstrong = 2.827
        characteristic_lenght_h2o_anstrong = 2.641
        characteristic_lenght_mixture = np.mean([characteristic_lenght_h2_anstrong,characteristic_lenght_h2o_anstrong])
        e_k_ratio_h2 = 59.7
        e_k_ratio_h2o = 809.1
        e_k_ratio_mixture = pow(e_k_ratio_h2*e_k_ratio_h2o,1/2)
        kT_e_ratio = temperature/e_k_ratio_mixture
        diffusion_collision_integral_dimentionless = 1.06036/pow(kT_e_ratio, 0.15610)+00.19300/np.exp(0.47635*kT_e_ratio)+1.03587/np.exp(1.52996*kT_e_ratio)+1.76474/np.exp(3.89411*kT_e_ratio)
        omega_v = 1.16145/pow(kT_e_ratio, 0.14874)+0.52487/np.exp(0.7732*kT_e_ratio)+2.16178/np.exp(2.43787*kT_e_ratio)
        diffusion_coefficient_ab_SI = 0.00266*pow(temperature, 3/2)/(pressure_bar*np.sqrt(molar_mass_mixture)*diffusion_collision_integral_dimentionless*pow(characteristic_lenght_mixture,2)*10000)
        diffusion_coefficient_ab_cm = 0.00266*pow(temperature, 1.5)/(pressure_bar*np.sqrt((molar_mass_mixture)*diffusion_collision_integral_dimentionless*pow(characteristic_lenght_mixture,2)))
        viscosity_h2 = (0.000026693*np.sqrt(molar_mass_h2*temperature)/(pow(characteristic_lenght_h2_anstrong,2)*omega_v))/10
        viscosity_h2o = 0.000026693*np.sqrt(molar_mass_h2o*temperature)/(pow(characteristic_lenght_h2o_anstrong,2)*omega_v)
        # viscosity_mixture = viscosity_h2/(1+(X_h2o/X_h2)*pow((1+np.sqrt(viscosity_h2/viscosity_h2o)*pow((molar_mass_h2o/molar_mass_h2),0.25)),2)/(4/np.sqrt(2)*(1+np.sqrt(molar_mass_h2/molar_mass_h2o))))+viscosity_h2o/(1+X_h2/X_h2o*(1+np.sqrt(viscosity_h2o/viscosity_h2)*pow(molar_mass_h2/molar_mass_h2o,0.25))^2/(4/np.sqrt(2)*(1+np.sqrt(molar_mass_h2o/molar_mass_h2))))
        kinematic_viscosity = viscosity_h2*8.314*temperature/((molar_mass_average/1000)*pressure_bar*100000)
        Sc = kinematic_viscosity/diffusion_coefficient_ab_SI
        gas_velocity = (temperature/298)*gas_volume_inlet/(np.pi*pow(tube_diameter_m,2)/4)
        reynolds_number_Re = gas_velocity*pellet_diameter_m/kinematic_viscosity
        sherwood_number_Sh = 2+0.552*np.sqrt(reynolds_number_Re)*pow(Sc,1/3)
        mass_transfer_coeff = sherwood_number_Sh*diffusion_coefficient_ab_SI/pellet_diameter_m

#turkdogan's models

        # logK= A+B/T A= 0.3956 B=-775.72
        K=pow(10,(0.3956-(775.72/temperature)))
        p_eq = pressure_bar/(1+K)
        phi = pow(10,((-6516/temperature)-1.31))
        self.mols_number = pellet_data.initial_mass/(55.85*0.95+16) #K24
        self.density_grams_m_cubic = pellet_data.initial_mass/((16+0.95*55.85))/((4/3)*np.pi*pow(self.pellet_radius_m,3)) # J24
        self.C_h2_bulk = X_h2*pressure_bar*101325/(8.314*temperature) # G28
        self.C_h2_eq = p_eq*101325/(8.314*temperature) # H28
        self.extraction_rate = mass_transfer_coeff*np.pi*pow(self.pellet_radius_m,2)*(self.C_h2_bulk-self.C_h2_eq)
        self.rate_constant_k = temperature*8.2*pow(10,-5)*pow(10,4)*phi/60

        print('diameter {}, mass transfer coeff {}, Re {}, Sh {}, Sc {}, Dab {}, kin.vis {}'.format(pellet_diameter_m,mass_transfer_coeff,reynolds_number_Re,sherwood_number_Sh, Sc, diffusion_coefficient_ab_SI, kinematic_viscosity))
        # # Equations
        # self.t_pore:Series = ((0.5-F/3-0.5*pow(1-F,2/3))*self.mols_density*pow(pellet_radius_m,2)/((De*0.0001)*(self.C_h2_bulk-self.C_h2_eq)))
        # self.t_external:Series = (F*self.mols_number/self.extraction_rate)
        # self.t_mixed:Series = ((1-pow(1-F,1/3))*pellet_radius_m*self.mols_density/((self.C_h2_bulk-self.C_h2_eq)*pow(self.phi*S*self.mols_density*De,1/2)))

    # def compose_functions(self, pore_constant, external_constant, mixed_constant):
    #     return self.t_pore*pore_constant + self.t_external*external_constant + self.t_mixed*mixed_constant

    def model_t_pore(self, F, De):
        value =  ((0.5-F/3-0.5*pow(1-F,2/3))*self.oxygen_density_wustite*pow(self.pellet_radius_m,2)/((De)*(self.C_h2_bulk-self.C_h2_eq)))
        return value
    
    def model_t_external(self, F, m):
        return (F*self.oxygen_density_hematite*self.pellet_radius_m/(3*m*(self.C_h2_bulk-self.C_h2_eq)))
    
    def model_t_mixed(self, F, De, S=1):
        return ((1-pow(1-F,1/3))*self.pellet_radius_m*self.oxygen_density_wustite/((self.C_h2_bulk-self.C_h2_eq)*pow(self.rate_constant_k*S*self.oxygen_density_wustite*De,1/2)))
