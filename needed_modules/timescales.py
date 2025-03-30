import numpy as np 
from scipy import optimize
import os
import re
import yaml

gammaD=41.695*10**6; #r*s^(-1)*T^(-1)
gammaH=267.513*10**6;
gammaC=67.262*10**6;
gammaN=-27.166*10**6;

def read_corr_func(correlation_file):
    # for reading the correlation function data
    opf = open(correlation_file, 'r')
    lines = opf.readlines()
    data_times = []
    data_F = []
    for i,line in enumerate(lines):
        if '#' in line:
            continue
        if '&' in line:
            continue
        if '@' in line:
            continue    
        if 'label' in line:
            continue
        if line == "":
            continue
        parts = line.split()
        if np.shape(parts)[0]==2:
            try:
                data_F.append(float(parts[1]))
                data_times.append(float(parts[0]))
            except:
                print(i)
                break
   
    data_Fout = np.array(data_F)
    times_out = np.array(data_times)
    return data_Fout, times_out

def get_timescales(correlation_file,OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze_until):
    org_corrF, times_out=read_corr_func(correlation_file)
    #analyze_until = round(len(org_corrF)*analyze)
    length=int(len(org_corrF)*(times_out[1]-times_out[0]))

    org_corrF=org_corrF[0:analyze_until]
    times_out=times_out[0:analyze_until]
    
    # normalized correlation fuction
    NcorrF = (org_corrF - OP ** 2) / (1 - OP ** 2);

    
    # Create correlation times from the times and number of exponential specified by the user
    step_exp=(biggest_corr_time-smallest_corr_time)/N_exp_to_fit
    Ctimes = 10 ** np.arange(smallest_corr_time, biggest_corr_time, step_exp)  # units: ps
    Ctimes = np.logspace(smallest_corr_time-12, biggest_corr_time-12, N_exp_to_fit)*10**12
    # First, no forcing the plateou
    # create exponential functions and put them into a matrix, individual exponentials in columns
    #the lengthe of correlationd data to be used is specified by the user
    n = len(times_out)
    m = len(Ctimes)
    Cexp_mat = np.zeros((n, m))

    for i in range(0, n):
        for j in range(0, m):
            Cexp_mat[i, j] = np.exp(-times_out[i] / Ctimes[j])

    #least square solution
    Coeffs, res = optimize.nnls(Cexp_mat, NcorrF[0:n])

    # Effective correlation time from components, in units of sec
    Teff = sum(Coeffs * Ctimes * 0.001 * 10 ** (-9))

    # calculate t_eff from area
    dt = times_out[2] - times_out[1]
    pos = np.argmax(NcorrF[0:n] < 0)

    if pos > 0:
        tau_eff_area = sum(NcorrF[0:pos]) * dt * 0.001 * 10 ** (-9)
    else:
        tau_eff_area = sum(NcorrF[0:n]) * dt * 0.001 * 10 ** (-9)

    # changin the unit of time permanently from [ps] to [s]
    Ctimes = Ctimes * 0.001 * 10 ** (-9)
    
    return Ctimes, Coeffs, Teff, tau_eff_area, org_corrF, times_out, length
    
def get_timescales_for_system(folder_path,OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,output_path,output_prefix,save_yaml=False,save_txt=True):
    output_name=f'{output_path}/{output_prefix}_timescales.yaml'
    try:
        with open(output_name) as yaml_file:
            content = yaml.load(yaml_file, Loader=yaml.FullLoader)
        
    except:
        content={}
   
    info={}
    info["07_OP"]=OP
    info["04_smallest_corr_time_[s]"]=10**(int(smallest_corr_time)-12)
    info["05_biggest_corr_time_[s]"]=10**(int(biggest_corr_time)-12)
    info["03_N_exp_to_fit"]=N_exp_to_fit
    info["06_analyze"]=analyze/100
    
    try:
        with open(f'{folder_path}/Correlation_files_INFO.yaml') as yaml_file:
            corr_readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
        info['10_xtc'] = corr_readme['xtc']
        info['11_tpr'] = corr_readme['tpr']
        info['12_xtc_modified'] = corr_readme['FROM_XTC']
    except:
        pass
    
    new_ana=False
    if len(content)==0:
        new_ana=True
    else:
        help1=False 
        for analysis in content:
            help2=False
            for value in info:
                try:
                    if not content[analysis]['info'][value]==info[value]:
                        help2=True
                except:
                    help2=True
            if not help2:
                help1=True
        if not help1:
            new_ana=True
    if save_txt:
        new_ana=True
    if new_ana:
        if 'Correlation_files_INFO.yaml' in os.listdir(folder_path) or 'README_correl.yaml' in os.listdir(folder_path):
            cr=True
            timescales=[0]*(len(os.listdir(folder_path)))
            eff_times=[0]*(len(os.listdir(folder_path))-1)
            eff_times_area=[0]*(len(os.listdir(folder_path))-1)
            org_data=[0]*(len(os.listdir(folder_path)))
        else:
            cr=False
            timescales=[0]*(len(os.listdir(folder_path))+1)
            eff_times=[0]*(len(os.listdir(folder_path)))
            eff_times_area=[0]*(len(os.listdir(folder_path)))
            org_data=[0]*(len(os.listdir(folder_path))+1)
    
        all_ind=[]
        for file in os.listdir(folder_path):
            if "INFO" not in os.fsdecode(file) and 'README' not in os.fsdecode(file):
                all_ind.append(int(file.split("_")[-1].split(".")[0]))
        all_ind.sort()        
        for file in os.listdir(folder_path):
            if "INFO" not in os.fsdecode(file) and 'README' not in os.fsdecode(file):
                pre_AA_index=int(file.split("_")[-1].split(".")[0])
                AA_index=all_ind.index(pre_AA_index)
                input_corr_file = folder_path+os.fsdecode(file)
                Ctimes, Coeffs, Teff, tau_eff_area, org_corrF, times_out, length=get_timescales(input_corr_file,OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze)
                tims=[]
                for c in Coeffs:
                    tims.append(float(c))
                timescales[AA_index+1]=tims
                eff_times[AA_index]=float(Teff)
                eff_times_area[AA_index]=float(tau_eff_area)
                org_data[AA_index+1]=list(org_corrF)
              
        org_data[0]=list(times_out)
    
        tims=[]
        for c in Ctimes:
            tims.append(float(c))
        timescales[0]=tims
    
        try:
            with open(folder_path+'Correlation_files_INFO.yaml') as yaml_file:
                corr_readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
            residues=corr_readme['RELAXATIONS']
        except:
            try: 
                with open(folder_path+'README_correl.yaml') as yaml_file:
                    corr_readme = yaml.load(yaml_file, Loader=yaml.FullLoader)
                residues=corr_readme['RELAXATIONS']
            except:
                residues={}
        
        artificials=[]
        for i in range(1,len(timescales)):
            if timescales[i][-1]>0:
                if len(residues)>0:
                    print(f'Artificially slow timescale detected for id {i-1}, residue {residues[i-1][0]}, bond {" - ".join(residues[i-1][1:])}')
                    artificials.append(f'Artificially slow timescale detected for id {i-1}, residue {residues[i-1][0]}, bond {" - ".join(residues[i-1][1:])}')
                else:
                    print(f'Artificially slow timescale detected for id {i-1}')
                    artificials.append(f'Artificially slow timescale detected for id {i-1}')
        info["08_corr_func_length_[ps]"]=length
        info["09_saving_frequency_[ps]"]=float((org_data[0][1]-org_data[0][0]))
        if save_txt:
            output_name=f'{output_path}/{output_prefix}'
            save_timescales_txt(output_name,timescales, eff_times, eff_times_area, residues,info,content)
        if save_yaml:
            output_name=f'{output_path}/{output_prefix}_timescales.yaml'
            yaml_timescales = save_timescales_yaml(output_name,timescales, eff_times, eff_times_area, residues,info,content)
        return timescales, eff_times, eff_times_area, residues, org_data, output_prefix, artificials, info
    print('   * Timescales from correlation fuction already calculated')
    return []
    
def save_timescales_yaml(output_name,timescales, eff_times, eff_times_area, residues,info,content):
       
    analysis=f'analysis{len(content)}'
        
        
    content[analysis]={} 
    content[analysis]['info']=info   
    content[analysis]["results"]={}
   
    
    content[analysis]["results"]['timescales']=timescales
    content[analysis]["results"]['eff_times'] = list(eff_times)
    content[analysis]["results"]['eff_times_area'] =  eff_times_area 
    content[analysis]['residues'] =  residues

    with open(output_name, 'w') as f:
        yaml.dump(content,f, sort_keys=True)
    
    return content[analysis]
    

def save_timescales_txt(output_name,timescales, eff_times, eff_times_area, residues,info,content):
    
    ts=np.transpose(np.array(timescales))
    with open(f'{output_name}_timescales.dat','w') as f:
        f.write('# Timescale analysis \n \n')
        for key, data in info.items():
            f.write(f'# {key}: {data} \n')
        f.write('\n')
        f.write('# Residues: \n')
        f.write('#')
        for res, data in residues.items():
            f.write(f'{res}: {data[0]}, {data[1]}, {data[2]}; ')
            
        f.write('\n \n')
        for row in ts:
            for i,column in enumerate(row):
                if i>0:
                    f.write(f'{column:.2f}  ')
                else:
                    f.write(f'{column:.3e}  ')
            f.write('\n')
