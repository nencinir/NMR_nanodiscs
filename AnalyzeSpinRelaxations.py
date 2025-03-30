import sys

import os
import users_settings as us
from datetime import date
import yaml
import numpy as  np
import fnmatch
sys.path.append('needed_modules/')
import manage_files as nm0
import correlation_functions as nm1
import timescales as nm2
import relaxations as nm3
#import matplotlib.pyplot as plt

#from matplotlib.backends.backend_pdf import PdfPages

gammaD=41.695*10**6; #r*s^(-1)*T^(-1)
gammaH=267.513*10**6;
gammaC=67.262*10**6;
gammaN=-27.166*10**6;


if us.magnetic_field_units=='MHz':
    us.magnetic_field=us.magnetic_field*2*np.pi/gammaH*10**6




if us.perform_analysis == 1:
    timescale_report=[]
    spin_report=[]
    compositions=[]
    temperatures=[]
    for filen in os.listdir(us.parent_folder_path):
        folder_path = us.parent_folder_path+os.fsdecode(filen)+"/"
        for system in us.systems:
            if fnmatch.fnmatch(os.fsdecode(filen), "*"+system+"*"):
                print(f' \n \n ########################### \n')
                #print(f' 1) Creating and updating README.yaml for \n    {folder_path} \n')
                #composition,temperature=nm0.go_through_simulation(folder_path)
                #compositions.append(composition)
                #temperatures.append(temperature)
                #nm0.remove_water(folder_path,us.selection,us.compress_xtc)
                #with open(f'{folder_path}/README.yaml') as yaml_file:
                #    readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                #xtcfile=f"{folder_path}/{readme['FILES_FOR_RELAXATION']['xtc']['NAME']}"
                #tprfile=f"{folder_path}/{readme['FILES_FOR_RELAXATION']['tpr']['NAME']}"
                title=filen
                output_name=title  

                #print(f'\n 2) Calculating Correlation functions for \n    {folder_path} \n')
                
                #nm1.calculate_correlation_functions(xtcfile,tprfile,us.output_path_correlation,output_name,us.end,us.begin,us.atom1_atom2_bonds, us.split_groups,title)
                correlations_path=f'{us.output_path_correlation}/{filen}/'
                
                print(f'\n 3) Calculating Timescales for \n    {folder_path} \n')               
                for analyze in list(range(2000,20500,500)):
                    print(f'analyzing {analyze/100} ns of corr. func \n')
                    timescale_return = nm2.get_timescales_for_system(correlations_path,us.OP,us.smallest_corr_time,us.biggest_corr_time, us.N_exp_to_fit,analyze,us.output_path_timescales,output_name,save_yaml=True,save_txt=us.save_timescales_txt)
                if len(timescale_return)>1:
                    timescale_report.append(timescale_return)
                #print(f'\n 4) Calculating Spin relaxation times for \n    {folder_path} \n')
                #timescales_file=f'{us.output_path_timescales}/{output_name}_timescales.yaml'
                 
                #for analyze in list(range(2000,20500,500)):
                #T1s, T2s, NOEs,  residues = nm3.get_spin_relaxation_times(us.magnetic_field,us.OP,us.smallest_corr_time, us.biggest_corr_time, us.N_exp_to_fit,analyze,timescales_file,us.nuclei,us.output_path_relaxations,output_name,save_yaml=True,save_txt=us.save_relaxations_txt)
                #spin_report.append((T1s,T2s,NOEs,residues))
    #nm3.print_report3(timescale_report,spin_report,compositions,temperatures,us.magnetic_field,us.nuclei,us.report_name)

                
if us.perform_analysis == 2:
    timescale_report=[]
    spin_report=[]
    composition,temperature=nm0.get_basic_info(us.tprfile)
    print(f'\n 1) Calculating Correlation functions for \n    {us.title} \n')
    nm1.calculate_correlation_functions(us.xtcfile,us.tprfile,us.output_path_correlation,us.output_name,us.end,us.begin,us.atom1_atom2_bonds, us.split_groups,us.title)
    print(f'\n 2) Calculating Timescales for \n    {us.title} \n')  
    correlations_path=f'{us.output_path_correlation}/{us.output_name}/'   
    timescale_return = nm2.get_timescales_for_system(correlations_path,us.OP,us.smallest_corr_time,us.biggest_corr_time, us.N_exp_to_fit,us.analyze,us.output_path_timescales,us.output_name,save_yaml=True,save_txt=us.save_timescales_txt)
    if len(timescale_return)>1:
        timescale_report.append(timescale_return)
    print(f'\n 4) Calculating Spin relaxation times for \n    {us.title} \n')
    timescales_file=f'{us.output_path_timescales}/{us.output_name}_timescales.yaml'
    T1s, T2s, NOEs,  residues = nm3.get_spin_relaxation_times(us.magnetic_field,us.OP,us.smallest_corr_time, us.biggest_corr_time, us.N_exp_to_fit,us.analyze,timescales_file,us.nuclei,us.output_path_relaxations,us.output_name,save_yaml=True,save_txt=us.save_relaxations_txt)
    spin_report.append((T1s,T2s,NOEs,residues))
    nm3.print_report3(timescale_report,spin_report,[composition],[temperature],us.magnetic_field,us.nuclei,us.report_name)

