import os 
from datetime import date
import time
import yaml

def create_index_file(xtcfile, tprfile, atom1, atom2, split_groups):
    
    grofile='gro.gro'
    os.system(f'echo System | gmx trjconv -f {xtcfile} -s {tprfile} -dump 0 -o {grofile} >& /dev/null')
    ndxfile="index_"+atom1+"_"+atom2+".ndx"
    
    if split_groups:    
        with open(grofile, 'rt') as gro_file:
            residue=""
            residues=0
            with open(ndxfile, 'w') as fo:
                for line in gro_file:
                    if 'Title' in line or len(line.split())==1 or len(line.split())==3:
                        pass
                    else:    
                    
                        if line.split()[1].startswith(atom1):
                            residue=line.split()[0]
                            if len(line.split())==5:
                                N=line.split()[1][len(atom1):]
                            else:
                                N=int(line.split()[2])
                        if line.split()[1].startswith(atom2):
                            if len(line.split())==5:
                                HN=line.split()[1][len(atom2):]
                            else:
                                HN=int(line.split()[2])
            
                            if residue==line.split()[0]:
                                fo.write("[ {} ]\n {} {}\n".format(residue,N,HN))
                                residues+=1
                                
    else:
        with open(grofile, 'rt') as gro_file:
            residue=""
            residues=1
            with open(ndxfile, 'w') as fo:
                fo.write("[ {}_{} ] \n".format(atom1,atom2))
                for line in gro_file:
                    if 'Title' in line or len(line.split())==1 or len(line.split())==3:
                        pass
                    else:    
                    
                        if line.split()[1].startswith(atom1):
                            residue=line.split()[0]
                            if len(line.split())==5:
                                N=line.split()[1][len(atom1):]
                            else:
                                N=int(line.split()[2])
                        if line.split()[1].startswith(atom2):
                            if len(line.split())==5:
                                HN=line.split()[1][len(atom2):]
                            else:
                                HN=int(line.split()[2])
                            if residue==line.split()[0]:
                                fo.write(" {} {}\n".format(N,HN))
    os.system(f'rm {grofile}')
    return ndxfile, residues
    

def individual_correlation_function(xtcfile,tprfile,output_path,output_name,end,begin,atom1, atom2, split_groups,number_from):
    """ Function to calculate Rotational Correlation functions from MD simulations.
    \n
    Calculates RCF for the enteries in the index file.
    \n
    Takes following arguments:
      xtcfile -  
      tprfile -  
      ndxfile
      output_path - parent folder to store correlation function folders
    \n
    Output:
        Creates a folder at working directory with the name of gro file and saves correlation functions there.
        When README.yaml available, it saves the path of the correlations functions and date of analysis there
    """


    

    ndxfile, residues = create_index_file(xtcfile, tprfile, atom1, atom2, split_groups)
    
    new_folder=output_path+output_name
    end_flag=f' -e {end} '
    if end==-1:
        end_flag=''
    if begin==-1:
        begin=0
    last_frame_should=-42
    
    
    correl={}
    correl['PROBLEMS']=[]
    
    ##### GET CORRELATION FUNCTIONS #####
    all_alright=True
    
    print("   * Number of corelation functions to calculate: {} \n".format(residues))
    for i in range(0+number_from,residues+number_from):
        print("     - Calculating correlation function {}".format(i+1),end=", ")
            
        os.system("echo " + str(i-number_from) + ' | gmx rotacf -f ' + xtcfile + ' -s ' + tprfile + '  -n ' + ndxfile + '  -o ' + new_folder + '/NHrotaCF_' + str(i) + ' -P 2 -d  ' + str(end_flag) + ' -b ' +str(begin)+' 2> corr.log')
        groups=[]
        with open("corr.log", 'rt') as corr_log:
            for line in corr_log:
                if "Reading frame" in line:
                    last_frame=int(float(line.split()[4]))
                if "Last frame" in line:
                    last_frame=int(float(line.split()[4]))
                if "Group" in line:
                    groups.append(line.split()[3])
                if "Done with trajectory" in line:
                    true_last_frame=last_frame
                    last_frame=last_frame_should
         
                  
        
        if not last_frame==last_frame_should:
            all_alright=False
            correl['PROBLEMS'].append(f'Problem at {last_slide}, residue{groups[i-number_from][0:len(groups[i-number_from])-1]}, index group {i-number_from}')
            print(" last frame",last_frame)
        
        else:
            print(" last frame",true_last_frame)
                               
        try:
            os.system("rm corr.log")
        except:
            pass
        
      
        
         
    today = str(date.today())        
    correl["name"]='a'
            
            
    timepre=os.path.getmtime(xtcfile)
    file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
    correl["FROM_XTC"]=file_mod
    correl["xtc"]=xtcfile
    correl["tpr"]=tprfile
        
        
    correl["BEGGIN"]=begin
    correl["END"]=end
       
    correl["RELAXATIONS"]={}  
    for i,group in enumerate(groups):
        correl["RELAXATIONS"][i+number_from]=group[0:len(groups[i])-1]
                      
        
    correl["LENGTH"]=true_last_frame
    correl["ANALYZED"]=today
    
    os.system(f'rm {ndxfile}')
    
    return residues, correl
    
    
def calculate_correlation_functions(xtcfile,tprfile,output_path,output_name,end,begin,atom1_atom2_bonds, split_groups,title):  
    
    new_folder=output_path+output_name
    nomore=False
    if begin==-1:
        begin=0
    if os.path.isdir(new_folder):
        correlF=new_folder+"/Correlation_files_INFO.yaml"
        if os.path.isfile(correlF):
            with open(correlF) as yaml_file:
                old_readme= yaml.load(yaml_file, Loader=yaml.FullLoader)
            nomore=True
            
            if not old_readme['xtc']==xtcfile:
                nomore=False
            if not old_readme['tpr']==tprfile:
                nomore=False
            timepre=os.path.getmtime(xtcfile)
            file_mod = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(timepre))
            if not old_readme["FROM_XTC"]==file_mod:
                nomore=False
        
            if not old_readme["BEGGIN"]==begin:
                nomore=False               
            if not old_readme["END"]==end:
                nomore=False
        if not nomore:            
            os.system("rm -r "+new_folder)
        else:
            print('   * Correlation fuction already calculated')
    if not nomore:
        os.system("mkdir " + new_folder)
    
        number_from=0

        readme={}
        readme['PROBLEMS']=[]
        readme['NAME']=title

        readme['RELAXATIONS']={}

        for bond in atom1_atom2_bonds:
            atom1=bond[0]
            atom2=bond[1]

            residues, correl=individual_correlation_function(xtcfile,tprfile,output_path,output_name,end,begin,atom1, atom2, split_groups, number_from)
            number_from+=residues
            readme['PROBLEMS']+=correl['PROBLEMS']
            for i in correl['RELAXATIONS']:
                readme['RELAXATIONS'][i]=[correl['RELAXATIONS'][i], atom1,atom2]

        readme['FROM_XTC']= correl['FROM_XTC']
        readme['xtc']=correl['xtc']
        readme['tpr']=correl['tpr']
        readme['BEGGIN'] = correl['BEGGIN']
        readme['END'] = correl['END'] 
        readme['LENGTH'] = correl['LENGTH'] 
        readme['ANALYZED'] = str(date.today())  


        with open(output_path+output_name+"/Correlation_files_INFO.yaml", 'w') as f:
            yaml.dump(readme,f, sort_keys=False) 
