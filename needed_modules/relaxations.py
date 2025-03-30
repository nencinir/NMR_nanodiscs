import yaml
import numpy as np

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

gammaD=41.695*10**6; #r*s^(-1)*T^(-1)
gammaH=267.513*10**6;
gammaC=67.262*10**6;
gammaN=-27.166*10**6;


def get_relaxation_D(magnetic_field,Coeffs,Ctimes,OP):
    omega = gammaD * magnetic_field
    
    #initiate spectral densities
    J0 = 0
    J1 = 0
    J2 = 0
    
    m = len(Ctimes)
    for i in range(0, m):
        w=0
        J0 = J0 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = omega
        J1 = J1 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = 2* omega
        J2 = J2 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

    xksi=167000 # quadrupolar coupling constant [Hz]
    R1 = 3 * (xksi  * np.pi) ** 2 / 40.0 * (1 - OP ** 2) * (0 * J0 + 2 * J1 + 8 * J2)
    R2 = 3 * (xksi  * np.pi) ** 2 / 40.0 * (1 - OP ** 2) * (3 * J0 + 5 * J1 + 2 * J2)

    return 1/R1, 1/R2, 0


def get_relaxation_C(magnetic_field,Coeffs,Ctimes,OP):
    omega = gammaD * magnetic_field
    
    wc = gammaC * magnetic_field;
    wh = gammaH * magnetic_field;
        
    #initiate spectral densities
    J0 = 0
    J1 = 0
    J2 = 0
    Jw1 = 0

    m = len(Ctimes)
    for i in range(0, m):
        w = wh - wc
        J0 = J0 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wc
        J1 = J1 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wc + wh
        J2 = J2 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

    # note! R1's are additive. Nh from the Ferreira2015 paper correctly omitted here
    R1 = (22000 * 2 * np.pi) ** 2 / 20.0 * (1 - OP ** 2) * (J0 + 3 * J1 + 6 * J2)


    return 1/R1, 0, 0


def get_relaxation_N(magnetic_field,Coeffs,Ctimes,OP):
    
    
    wh = gammaH * magnetic_field 
    wn = gammaN * magnetic_field 
    
    #initiate spectral densities
    J0 = 0
    JhMn = 0
    JhPn = 0
    Jh = 0
    Jn = 0

    m = len(Ctimes)
    for i in range(0, m):
        w = 0
      
        J0 = J0 + 2 * Coeffs[i] * Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])
        
        w = wh-wn;
        JhMn = JhMn + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wn;
        Jn = Jn + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])
        
        w = wh;
        Jh= Jh + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])

        w = wn+wh;
        JhPn = JhPn + 2 * Coeffs[i]* Ctimes[i] / (1.0 + w * w * Ctimes[i] * Ctimes[i])


    mu = 4 * np.pi * 10**(-7) #magnetic constant of vacuum permeability
    h_planck = 1.055 * 10**(-34); #reduced Planck constant
    rN = 0.101 * 10**(-9); # average cubic length of N-H bond
    d = 1 * (mu * gammaN * gammaH * h_planck) / (4 * np.pi * rN**3); # dipolar coupling constant

    #units were corrected by S.Ollila and E.Mantzari, removed 2*pi from R1 and R2
    R1 = (d**2 / 20) * (1 * JhMn + 3 * Jn + 6 * JhPn) + Jn * (wn * 160 * 10**(-6))**2 / 15   ; 
    R2 = 0.5 * (d**2 / 20) * (4 * J0 + 3 * Jn + 1 * JhMn + 6 * Jh + 6 * JhPn) + (wn * 160 * 10**(-6))**2 / 90 * (4 * J0 + 3 * Jn);
    NOE = 1 + (d**2 / 20) * (6 * JhPn - 1 * JhMn) * gammaH / (gammaN * R1);


           
    return 1/R1, 1/R2, NOE
    


        
         
choose_nuclei = {
    "13C": get_relaxation_C,
    "2H": get_relaxation_D,
    "15N": get_relaxation_N
}




def get_spin_relaxation_times(magnetic_field,OP,smallest_corr_time, biggest_corr_time, N_exp_to_fit,analyze,timescales_file,nuclei,path_relax,prefix_relax,save_yaml=False,save_txt=True):
    run_ana=False

    with open(timescales_file) as yaml_file:
        content = yaml.load(yaml_file, Loader=yaml.FullLoader)
    
    
    for ana in content:
        help2=False
        if not content[ana]['info']['07_OP']==OP:
            help2=True
        if not content[ana]['info']['04_smallest_corr_time_[s]']==10**(int(smallest_corr_time)-12):
            help2=True
        if not content[ana]['info']['05_biggest_corr_time_[s]']==10**(int(biggest_corr_time)-12):
            help2=True
        if not content[ana]['info']['03_N_exp_to_fit']==N_exp_to_fit:
            help2=True
        if not content[ana]['info']['06_analyze']==analyze:
            help2=True
        if not help2:
            T1s, T2s, NOEs= [], [], []
            Ctimes=content[ana]["results"]['timescales'][0]
            artificials=[]
            for i in range(1,len(content[ana]["results"]['timescales'])):
                Coeffs=content[ana]["results"]['timescales'][i]
                if Coeffs[-1]>0:
                    print(f'residue {i-1} has non zero biggest timescale, setting to 0 for spin relaxation times calculations')
                    Coeffs[-1]=0
                    artificials.append(i-1)
                T1, T2, NOE = choose_nuclei[nuclei](magnetic_field,Coeffs,Ctimes,OP) 
                T1s.append(T1)
                T2s.append(T2)
                NOEs.append(NOE) 
            if save_yaml:
              
                relax_yaml=f'{path_relax}/{prefix_relax}_relaxations.yaml'
                save_relaxations_yaml(relax_yaml,OP,smallest_corr_time,biggest_corr_time,N_exp_to_fit,analyze,magnetic_field,nuclei,T1s,T2s,NOEs,content[ana]['residues'],artificials)
            if save_txt:
                relax_txt=f'{path_relax}/{prefix_relax}_relaxations.dat'
                save_relaxations_txt(relax_txt,OP,smallest_corr_time,biggest_corr_time,N_exp_to_fit,analyze,magnetic_field,nuclei,T1s,T2s,NOEs,content[ana]['residues'],artificials)
                
    return T1s, T2s, NOEs,  content[ana]['residues']


def save_relaxations_yaml(relax_yaml,OP,smallest_corr_time,biggest_corr_time,N_exp_to_fit,analyze,magnetic_field,nuclei,T1s,T2s,NOEs,residues,artificials):
    try:
        with open(relax_yaml) as yaml_file:
            content = yaml.load(yaml_file, Loader=yaml.FullLoader)
    except:
        content={}
    info={}
    info['07_OP']=OP
    info['04_smallest_corr_time_[s]']=10**(int(smallest_corr_time)-12)
    info['05_biggest_corr_time_[s]']=10**(int(biggest_corr_time)-12)
    info['03_N_exp_to_fit']=N_exp_to_fit
    info['06_analyze']=analyze
    info["01_magnetic_field_[T]"]=magnetic_field
    info["02_magnetic_field_[MHz]"]=float(np.round(magnetic_field /(2*np.pi/gammaH*10**6),2))
    info["00_nuclei"]=nuclei

    new_save=True

    for ana in content:
        if content[ana]['info']==info:
            new_save=False
    if new_save:
        ana=f'analysis{len(content)}'
        content[ana]={}
        content[ana]['info']=info
        content[ana]['residues']=residues
        content[ana]['T1']=T1s
        content[ana]['T2']=T2s
        content[ana]['hetNOE']=NOEs
        content[ana]['artificial_timescales']=artificials
                        
                    
        with open(relax_yaml, 'w') as f:
            yaml.dump(content,f, sort_keys=True)

def save_relaxations_txt(relax_txt,OP,smallest_corr_time,biggest_corr_time,N_exp_to_fit,analyze,magnetic_field,nuclei,T1s,T2s,NOEs,residues,artificials):
    info={}
    info['07_OP']=OP
    info['04_smallest_corr_time_[s]']=10**(int(smallest_corr_time)-12)
    info['05_biggest_corr_time_[s]']=10**(int(biggest_corr_time)-12)
    info['03_N_exp_to_fit']=N_exp_to_fit
    info['06_analyze']=analyze
    info["01_magnetic_field_[T]"]=magnetic_field
    info["02_magnetic_field_[MHz]"]=float(np.round(magnetic_field /(2*np.pi/gammaH*10**6),2))
    info["00_nuclei"]=nuclei
    
    resi=[]
    for rID in residues:
       resi.append(rID)
    
    ts=np.transpose(np.array([resi,T1s,T2s,NOEs]))
    with open(relax_txt,'w') as f:
        f.write('# Spin relaxation times analysis \n \n')
        for key, data in info.items():
            f.write(f'# {key}: {data} \n')
        f.write('\n')
        f.write('# Residues: \n')
        f.write('#')
        for res, data in residues.items():
            f.write(f'{res}: {data[0]}, {data[1]}, {data[2]}; ')
        f.write('\n \n')
        f.write('# Artificial Timescales at Residues: ')
        f.write(", ".join(map(str,artificials)))
            
        f.write('\n \n')
  
        for row in ts:
            for i,column in enumerate(row):
                if i>0:
                    f.write(f'{column:10.3f}  ')
                else:
                    f.write(f'{int(column):3d}  ')
                
            f.write('\n')
    pass

def exp_sum(xvalues,timescales):
    results=[]
    for res in range(1,len(timescales)):
        corr=np.zeros(len(xvalues))
        for i,timescale in enumerate(timescales[0]):
            corr+=timescales[res][i]*np.exp(-xvalues/(timescale*10**12))
        results.append(corr)
    return results    
    
        
        
def print_report3(timescale_report,spin_report,compositions,temperatures,magnetic_field,nuclei,report_name):               
    magnetic_field_MHz=float(np.round(magnetic_field /(2*np.pi/gammaH*10**6),2))
    with PdfPages(report_name) as pdf:    
        rows=4
        cols=2
        per_page=rows*cols

        plt.rcParams["figure.figsize"] = [8.27, 11.7]
        plt.rcParams["figure.autolayout"] = True
        
        fig, axs = plt.subplots(rows, cols)
        i=0
        for i,system in enumerate(timescale_report):
            plt.rcParams["figure.figsize"] = [8.25, 11.75]
            plt.rcParams["figure.autolayout"] = True
            fig, axs = plt.subplots(rows, cols)
            xvalues=np.array(system[4][0])
            fits=exp_sum(xvalues,system[0])

            for j in range(1,len(system[4])):
                axs[0,0].plot(xvalues/1000,system[4][j],"-",c=f"C{j-1}")
                
                axs[0,0].plot(xvalues/1000,fits[j-1],"--",c=f"C{j-1}")
            
            
            axs[0,0].set_xlabel('Time [ns]')
            axs[0,0].set_ylabel('Correlation function')
            
            for residue in range(1,len(system[0])):
                for j,weight in enumerate(system[0][residue]):
                    ms=int(np.round(weight*20,0))
                    axs[0,1].set_yscale('log')
                    axs[0,1].plot(residue-1,system[0][0][j]*10**9, "o", markersize=ms, markeredgecolor='C0', markerfacecolor='C0')
            
            
            axs[0,1].set_xlabel('Residue ID')
            axs[0,1].set_ylabel('Timescale [ns]')
            
            plt.suptitle(system[5],size=15)
            
            res=spin_report[i][3].keys()
            labels=[]
            for key,data in spin_report[i][3].items():
                labels.append(data[0])
            
            axs[1,0].set_xticks(list(res), labels, rotation='vertical',fontsize=4)
            axs[1,0].plot(res,spin_report[i][0])
            axs[1,0].set_xlabel('Residue')
            axs[1,0].set_ylabel('T1 [s]')
            
            axs[1,1].set_xticks(list(res), labels, rotation='vertical',fontsize=4)
            axs[1,1].plot(res,spin_report[i][1])
            axs[1,1].set_xlabel('Residue')
            axs[1,1].set_ylabel('T2 [s]')
            
            axs[2,0].set_xticks(list(res), labels, rotation='vertical',fontsize=4)
            axs[2,0].plot(res,spin_report[i][2])
            axs[2,0].set_xlabel('Residue')
            axs[2,0].set_ylabel('hetNOE')
            
            axs[2,1].axis('off')
            
            axs[3,0].axis('off')    
            axs[3,0].text(0,0,'\n'.join(system[6]),fontsize=6, fontdict=None)
            
            axs[3,1].axis('off')
            text=""
            for key,data in system[7].items():
                text+=f'{str(" ".join(key.split("_")[1:])):<30s}'
                
                data=str(data).strip()
                for l in range(len(str(data))//40+1):
                    if l==0:
                        text+=f"{data[l*40:(l+1)*40].strip():40s}\n"
                        
                    else:
                        text+=f"{' '*30}{data[l*40:(l+1)*40].strip()}\n"
            text+=f'{"Temperature [K]":<30s}{temperatures[i]}\n'
            text+=f'{"Magnetic field [T]":<30s}{magnetic_field:.2f}\n'
            text+=f'{"Magnetic field [MHz]":<30s}{magnetic_field_MHz}\n'
            text+=f'{"Nuclei":<30s}{nuclei}\n'
            text+=f'Composition: \n'
            for co,nu in compositions[i].items():
                text+=f'    {co:<15s}{nu}\n'
            axs[3,1].text(0,0,text,fontsize=6,  family='monospace')
            
            
            pdf.savefig()  # saves the current figure into a pdf page
                  
        plt.close()                 
