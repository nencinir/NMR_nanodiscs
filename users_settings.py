perform_analysis = 1   # 1/2/3   #1 - goes through fodlers in parent folder, takes '*.xtc', '*tpr' from there, calculates correlation functions, timescales, spin relaxation times
                                 #2 - user inserts '*.xtc' and '*.tpr' files, code calculates correlation functions, timescales, spin relaxation times
                                 #3 - user inserts folder path with timescale yaml files, code calculates spin relaxation times    - will come

########################################################################
# Analysis 1 settings
########################################################################

if perform_analysis == 1:
    """Path settings, saving settings"""
    parent_folder_path="/home/nenciric/Documents/NMR_paper/corr_func/" # A parent folder that contains subfolders with simulations    
    systems=[""]                            # select only systems which folder name contains some of these
    output_path_correlation="/home/nenciric/Documents/NMR_paper/corr_func/"   
    output_path_timescales="//home/nenciric/Documents/NMR_paper/timescales/" 
    output_path_relaxations="//home/nenciric/Documents/NMR_paper/relaxations/" 
    
    save_timescales_txt=False # Generate also txt files?, Results are already saved to yaml files             
    save_relaxations_txt=False # Generate also txt files?, Results are already saved to yaml files 
    
    report_name='analysis_report.pdf'
                
    """Handeling of trajectory file"""                                                  
    compress_xtc=False # True/False/"Original" 
                       #     True - creates tpr, xtc, gro with selection "selection"
                       #     False - creates tpr, gro with the selection, 
                       #             assumes that reduced xtc already exists
                       #     "Original" - uses original trajectories for correlation function calculation
                       #                slows significantly down the calculations
                    
    selection="non-Water"  # selection for the compression
                         # at the moment only standard selections that exists in a default index file supported
                         # such as non-Water 
                 

    """Settings for correlation times calculation"""                 
    begin=-1     # trajectory to be analyzed from the time 'begin' [ps], -1 means from the start
    end=-1       # trajectory to be analyzed until the tme 'end'  [ps]


    split_groups=True   # True/False        #  True      - creates an index file with all atom1-atom2 pairs separatelly
                                            #  False     - creates an index file, where all atom1-atom2 pairs are in one group, the correlation function 
                                                         # is an average over these then

    atom1_atom2_bonds=[('N','HN')] # An example for 'Protein'

    #atom1_atom2_bonds=[("C1","H11"),("C2","H21"),("C3","H31"), 
    #("C4","H41"),("C5","H51"),("C6","H61"),
    #("C7","H71"),("C8","H81"),("C9","H91"),
    #("C10","H101"),("C11","H111"),("C12","H121")]  # En example for lipid


    """Settings for timescale analysis"""
    OP=0 # order parameter
    smallest_corr_time=0 # enter in log scale -3 fs; 0 ps; 3 ns; 6 us;
    biggest_corr_time=6 # same as above
    N_exp_to_fit=181 # number of exponential functions to be fitted between the samlles and biggest corr time
    analyze=1/12 # the proportin of correlation data to be used for fitting, ex. 1/2 uses first half of the data
                 # keep in mind that length of the corr. function is already 1/2 of the length of the simulation
    magnetic_field=850
    magnetic_field_units='MHz' # 'MHz'/'T'
    nuclei="15N" #nuclei to calculate: 2H-deutherium; 13C - carbon; 15N - nitrogen 

########################################################################
# Analysis 2 settings
########################################################################

if perform_analysis == 2:
    """Path and saving settings"""
    xtcfile="//media/ricky/Ricky2020/2020/akseli/aa_70_uf_md/aa_70_uf_md.xtc"
    tprfile="//media/ricky/Ricky2020/2020/akseli/aa_70_uf_md/no_water.tpr" 
    title='nanodisks ares good'
    output_name='aa_70_uf'  #used for correlation, timescales, spin relaxation times
    
    save_timescales_txt=False # Generate also txt files?, Results are already saved to yaml files             
    save_relaxations_txt=False # Generate also txt files?, Results are already saved to yaml files 
    
    report_name='aa_70_uf_lipid.pdf'
                


    output_path_correlation="/home/ricky/Documents/school/thesis/66_timescales_lipids_nanodisks/corr_func/"   
    output_path_timescales="//home/ricky/Documents/school/thesis/66_timescales_lipids_nanodisks/timescales_yamls/" 
    output_path_relaxations="//home/ricky/Documents/school/thesis/66_timescales_lipids_nanodisks/relaxations/" 
    
    """Settings for correlation times calculation"""
    end=-1       # trajectory to be analyzed until the tme 'end'  [ps]
    begin=-1     # trajectory to be analyzed from the time 'begin' [ps], -1 means from the start
    atom1_atom2_bonds=[('C12','H2R'),('C13','H3R'),('C14','H4R'),('C15','H5R'),('C16','H6R'),('C17','H7R'),
                      ('C18','H8R'),('C19','H9R'),('C110','H10R'),('C111','H11R'),('C112','H12R'),('C113','H13R')] # An example for 'Protein'

    split_groups=False   # True/False        #  True      - creates an index file with all atom1-atom2 pairs separatelly
                                            #  False     - creates an index file, where all atom1-atom2 pairs are in one group, the correlation function 
                                                         # is an average over these then
    
    
    """Settings for timescale analysis"""
    OP=0 # order parameter
    smallest_corr_time=0 # enter in log scale -3 fs; 0 ps; 3 ns; 6 us;
    biggest_corr_time=11 # same as above
    N_exp_to_fit=141 # number of exponential functions to be fitted between the samlles and biggest corr time
    analyze=1/10 # the proportin of correlation data to be used for fitting, ex. 1/2 uses first half of the data
                 # keep in mind that length of the corr. function is already 1/2 of the length of the simulation
    magnetic_field=850
    magnetic_field_units='MHz' # 'MHz'/'T'
    nuclei="2H" #nuclei to calculate: 2H-deutherium; 13C - carbon; 15N - nitrogen 
   

