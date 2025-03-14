import pandas as pd
import sqlite3
import numpy as np
from IPython.display import display, HTML
import json
import os
sqlite3.register_adapter(np.int64, lambda val: int(val))
sqlite3.register_adapter(np.int32, lambda val: int(val))

def efield_to_SHE(d, pzc, efield):
    she = efield*d+pzc
    return she

def SHE_to_efield(d, pzc, she):
    efield = (she-pzc)/d
    efield = round(efield, 3)
    return efield

def display_df(df, title=None):
    html_table = df.to_html()
    if title:       
        html_with_title = f"<h3>{title}</h3>{html_table}"
    else:
        html_with_title = html_table
    display(HTML(html_with_title))

ts_dict = {
    'CO2-H-ele': {
        'activation_energy': 0.06, #+0.54 eV = 0.60 eV
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 32,
        'ae': 'Ea'
    },
    'COOH-H-ele': {
        'activation_energy': 0.65, #-0.25 eV = 0.40 eV
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 33,
        'ae': 'Ea'
    },
    'H-CO2-ele': {
        'activation_energy': 0.52, #+0.48 eV = 1.00 eV
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 34,
        'ae': 'Ea'
    },
    'HCOO-H-ele': {
        'activation_energy': 0.66, #-0.26 eV = 0.40 eV
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 35,
        'ae': 'Ea'
    },
}

def generate_formation_energies(metal='Cu', facets=['100', '211', '111'], 
                                she_range=(-2.00, -0.40, 0.02),
                                db_path=None,
                                if_efield=True,
                                ts_dict=ts_dict,
                                correction_dict=None,
                                display_potential=False,
                                save_dir=None):
    """
    Generate formation energies for different metal surfaces and facets
    """
    print(f"\n{'='*50}")
    print(f"Starting formation energy generation for {metal}")
    print(f"{'='*50}")
    
    # Dictionary for metal-specific parameters
    metal_params = {
        'Cu': {
            'd': 1.2,  # Angstrom
            'pzc': -0.54,  # V_SHE
            'she_range': she_range,
            'gas_dict': {
                'H': {'formation_energy': 0.0, 'reference': 'This work'},
                'O2': {'formation_energy': 0.0, 'reference': 'Fake ads'},
                'ele': {'formation_energy': 0.0, 'reference': 'This work'}
            },
            'fallback_map': {
                'OCCH': 'CHCO',
                'OCC': 'CCO'
            }
        }
    }

    
    if metal not in metal_params:
        raise ValueError(f"Parameters for {metal} not defined")
    
    params = metal_params[metal]
    total_facets = len(facets)
    
    conn = sqlite3.connect(db_path)
    try:
        for facet_idx, facet in enumerate(facets, 1):
            print(f"\n{'-'*50}")
            print(f"Processing {metal}{facet} ({facet_idx}/{total_facets})")
            print(f"{'-'*50}")
            # Update ads_dict for current facet
            ads_dict = {
                'O2': {
                    'formation_energy': 9.0,
                    'reference': 'Fake ads',
                    'surface_name': metal,
                    'site_name': facet,
                    'frequencies': []
                },
                'OCCH': {
                    'formation_energy': None,
                    'reference': 'This work',
                    'surface_name': metal,
                    'site_name': facet,
                    'frequencies': []
                },
                'OCCHO': {
                    'formation_energy': None,
                    'reference': 'This work',
                    'surface_name': metal,
                    'site_name': facet,
                    'frequencies': []
                }
            }
            # Step 1: Get unique species
            print("\nStep 1: Getting unique species...")
            unique_species = pd.read_sql_query("""
                SELECT DISTINCT Species.species_name, Status.status
                FROM ReactionStepSpecies rss
                JOIN Species ON rss.species_id = Species.species_id
                JOIN Status ON rss.status_id = Status.status_id
            """, conn)
            print(f"Found {len(unique_species)} unique species")
            
            # Step 2: Create initial dataframe
            print("\nStep 2: Creating initial dataframe...")
            df_empty = create_initial_dataframe(
                unique_species, 
                metal, 
                facet, 
                conn, 
                params['fallback_map'],
                params['gas_dict'],
                ads_dict
            )
            #display(df_empty)
            print(f"Initial dataframe created with {len(df_empty)} rows")
            

            
            if if_efield:
                # Step 4: Apply electric field effects
                print("\nStep 3: Applying electric field effects...")
                start, end, step = params['she_range']
                total_potentials = len(np.arange(start, end, step))
                
                for pot_idx, she in enumerate(np.arange(start, end, step), 1):
                    efield = SHE_to_efield(params['d'], params['pzc'], she)
                    print(f"\rProcessing potential {pot_idx}/{total_potentials}: SHE = {she:.2f}V", end='')
                    
                    # Create output directory
                    if save_dir:
                        output_dir = f'{save_dir}/{metal}{facet}/SHE_{she:.2f}V'
                    else:
                        output_dir = f'{metal}{facet}/SHE_{she:.2f}V'
                    os.makedirs(output_dir, exist_ok=True)
                    
                    # Apply efield corrections and save results
                    df_with_efield = apply_efield_corrections(
                        df_empty,
                        efield,
                        metal,
                        facet,
                        conn,
                        params['fallback_map']
                    )
                    # Apply correction to formation_energy
                    # if correction_dict:
                    #     print(f"Applying correction to formation energy")
                    #     df_with_efield = apply_formation_energy_corrections(df_with_efield, correction_dict)

                    df_after_ts = calculate_all_activation_energies(df_with_efield, conn, ts_dict)
                    
                    # Apply correction to formation_energy
                    if correction_dict:
                        print(f"Applying correction to formation energy")
                        df_after_ts = apply_formation_energy_corrections(df_after_ts, correction_dict)


                    # Save results
                    output_path = f'{output_dir}/input.txt'
                    df_after_ts.to_csv(output_path, sep='\t', index=False)
                    
                    # Display specific potentials if requested
                    if display_potential:
                        if abs(she - params['pzc']) < 0.001 or abs(she + 0.42) < 0.001:
                            display_df(df_after_ts, f"SHE_{she:.2f}V/input.txt")
            else:
                # Step 3: Calculate activation energies
                print("\nStep 3: Calculating activation energies...")
                df_after_ts = calculate_all_activation_energies(df_empty, conn, ts_dict)
                ts_count = len(df_after_ts[df_after_ts['status'] == 'ts'])
                print(f"Processed {ts_count} transition states")
                df_with_efield_after_ts = df_after_ts

                os.makedirs(f'{metal}{facet}', exist_ok=True)
                df_after_ts.to_csv(f'{metal}{facet}/input.txt', sep='\t', index=False)
            
            print(f"\nCompleted processing for {metal}{facet}")


    except Exception as e:
        print(f"\nError processing {metal}{facet}: {str(e)}")
        raise
    finally:
        conn.close()
    
    print(f"\n{'='*50}")
    print(f"Completed formation energy generation for {metal}")
    print(f"{'='*50}")
    
def apply_formation_energy_corrections(df, correction_dict):
    """
    Apply formation energy corrections to dataframe
    """
    for species, correction in correction_dict.items():
        df.loc[(df['species_name'] == species) & (df['status'] != 'gas'), 'formation_energy'] += correction
        # # Add note
        # if df.loc[df['species_name'] == species, 'status'] == 'ads':
        #     df.loc[df['species_name'] == species, 'note'] = f"additional correction: {correction} eV"
        # elif df.loc[df['species_name'] == species, 'status'] == 'ts':
        #     df.loc[df['species_name'] == species, 'note'] += f" | {correction} eV"

    return df

def create_initial_dataframe(unique_species, metal, facet, conn, fallback_map, gas_dict, ads_dict):
    """
    Create initial dataframe with formation energies
    
    Parameters:
    -----------
    unique_species : DataFrame
        DataFrame containing unique species and their status
    metal : str
        Metal name (e.g., 'Cu')
    facet : str
        Surface facet (e.g., '100')
    conn : sqlite3.Connection
        Database connection
    fallback_map : dict
        Dictionary mapping species names to alternative names
    gas_dict : dict
        Dictionary containing gas phase species parameters
    ads_dict : dict
        Dictionary containing adsorbate species parameters
    
    Returns:
    --------
    DataFrame
        Initial dataframe with formation energies
    """
    columns = ['status', 'surface_name', 'site_name', 'species_name', 
               'formation_energy', 'frequencies', 'reference', 
               'active_site', 'cell_size','solvation_energy']
    df_empty = pd.DataFrame(columns=columns)
    new_rows_without_ts = []
    print("\tProcessing species:", end=' ')
    total_species = len(unique_species)
    try:
        for idx, row in unique_species.iterrows():
            print(f"\r\tProcessing species: {idx+1}/{total_species}", end='')
            species_name = row['species_name']
            status = row['status']
            
            formation_energy = reference = surface_name = site_name = frequencies = None
            active_site = cell_size = None
            
            # Get frequencies
            freq_query = """
            SELECT Reference.reference, Frequencies.frequencies
            FROM Frequencies
            JOIN Species ON Frequencies.species_id = Species.species_id
            JOIN Status ON Frequencies.status_id = Status.status_id
            JOIN Reference ON Frequencies.reference_id = Reference.reference_id
            WHERE Species.species_name = ? AND Status.status = ?
            """
            
            frequencies = pd.read_sql_query(freq_query, conn, params=(species_name, status))
            
            # Try fallback name for frequencies if needed
            if frequencies.empty:
                fallback_name = fallback_map.get(species_name, species_name)
                if fallback_name != species_name:
                    frequencies = pd.read_sql_query(freq_query, conn, params=(fallback_name, status))
            
            if not frequencies.empty:
                # Prioritize specific references
                if 'Ara' in frequencies['reference'].values:
                    frequencies = frequencies[frequencies['reference'] == 'Ara']
                elif 'PengRole2020' in frequencies['reference'].values:
                    frequencies = frequencies[frequencies['reference'] == 'PengRole2020']
                else:
                    frequencies = frequencies.iloc[[0]]
                
                frequencies = frequencies['frequencies'].tolist()[0]
            else:
                frequencies = []
            
            # Get formation energies
            energy_query = """
            SELECT 
                Catalysts.surface_name, 
                Species.species_name, 
                Status.status, 
                Facets.facet as site_name, 
                fm.formation_energy, 
                r.reference, 
                fm.efield, 
                fm.active_site, 
                fm.cell_size, 
                fm.coverage, 
                r.reference_id
            FROM FormationEnergy fm
            JOIN Species ON fm.species_id = Species.species_id
            JOIN Status ON fm.status_id = Status.status_id
            JOIN Reference r ON fm.reference_id = r.reference_id
            JOIN Facets ON fm.facet_id = Facets.facet_id
            JOIN Catalysts ON fm.surface_id = Catalysts.surface_id
            WHERE Species.species_name = ? 
            AND Status.status = ? 
            AND (efield BETWEEN -0.1 AND 0.1 OR efield IS NULL)
            ORDER BY r.reference_id
            """
            
            energy_data = pd.read_sql_query(energy_query, conn, params=(species_name, status))
            
            if energy_data.empty:
                fallback_name = fallback_map.get(species_name, species_name)
                if fallback_name != species_name:
                    energy_data = pd.read_sql_query(energy_query, conn, params=(fallback_name, status))
                    
            # Process based on status
            if status == 'ads':
                energy_data = energy_data[
                    (energy_data['surface_name'] == metal) & 
                    (energy_data['site_name'] == facet)
                ]
                
                if '-' in species_name:
                    if species_name.endswith('-'):
                        status = 'ads'
                    else:
                        status = 'ts'
                    
                if energy_data.empty and species_name in ads_dict:
                    formation_energy = ads_dict[species_name]['formation_energy']
                    reference = ads_dict[species_name]['reference']
                    surface_name = ads_dict[species_name]['surface_name']
                    site_name = ads_dict[species_name]['site_name']
                    frequencies = ads_dict[species_name]['frequencies']
                    
            elif status == 'gas':
                surface_name = 'None'
                site_name = 'gas'
                if energy_data.empty and species_name in gas_dict:
                    formation_energy = gas_dict[species_name]['formation_energy']
                    reference = gas_dict[species_name]['reference']
            
            # Process energy data if available
            if not energy_data.empty:
                # Select best reference data
                energy_data = select_best_reference_data(energy_data, species_name)
                
                formation_energy = energy_data['formation_energy']
                reference = energy_data['reference']
                active_site = energy_data['active_site']
                cell_size = energy_data['cell_size']
                
                if status == 'ads':
                    surface_name = energy_data['surface_name']
                    site_name = energy_data['site_name']
            
            # Add solvation effects
            solvation_query = """
            SELECT sol.solvation_energy
            FROM SolvationEffect sol
            JOIN Species ON sol.species_id = Species.species_id
            JOIN Status ON sol.status_id = Status.status_id
            WHERE Species.species_name = ? 
            AND Status.status = ? 
            AND sol.approach = 'explicit'
            """
            
            solvation_data = pd.read_sql_query(solvation_query, conn, params=(species_name, status))
            
            sol_note = None
            if not solvation_data.empty:
                solvation_energy = solvation_data['solvation_energy'].iloc[0]
                if formation_energy is not None:
                    formation_energy += solvation_energy
                    sol_note = f'{solvation_energy:.2f} eV'
            
            # Create new row
            new_row = {
                'status': status,
                'surface_name': surface_name,
                'site_name': site_name,
                'species_name': species_name,
                'formation_energy': formation_energy,
                'frequencies': frequencies,
                'reference': reference,
                'active_site': active_site,
                'cell_size': cell_size,
                'solvation_energy': sol_note
            }
            new_rows_without_ts.append(new_row)
            
        # Create final DataFrame
        df = pd.concat([df_empty, pd.DataFrame(new_rows_without_ts)], ignore_index=True)
        
        # Sort by status
        return df.sort_values(
            by='status', 
            key=lambda x: pd.Categorical(x, categories=['gas', 'ads', 'ts'], ordered=True)
        )
        
    except Exception as e:
        print(f"Error in create_initial_dataframe: {str(e)}")
        raise

def select_best_reference_data(energy_data, species_name):
    """Select best reference data based on priority"""
    if 'Ara' in energy_data['reference'].values:
        return energy_data[energy_data['reference'] == 'Ara'].iloc[0]
    elif 'PengRole2020' in energy_data['reference'].values:
        return energy_data[energy_data['reference'] == 'PengRole2020'].iloc[0]
    else:
        return energy_data.iloc[0]
    

def calculate_all_activation_energies(df, conn, ts_dict):
    """
    Calculate activation energies for transition states
    
    Parameters:
    -----------
    df : DataFrame
        Input DataFrame containing species data
    conn : sqlite3.Connection
        Database connection
    
    Returns:
    --------
    DataFrame
        DataFrame with updated activation energies
    """
    new_rows = []
    

    
    try:
        for _, row in df.iterrows():
            
            species_name = row['species_name']
            status = row['status']
            new_row = row.copy()
            # if status == 'ts':
            #     # Get activation energy data
            #     activation_data = get_activation_energy_data(conn, species_name, ts_dict, row)
                
            #     if activation_data:
            #         # Calculate formation energy including activation energy
            #         new_row = calculate_formation_energy_with_activation(
            #             conn, df, row, activation_data
            #         )
                    
            new_row = calculate_activation_energy(row, df, conn, ts_dict)
            
            new_rows.append(new_row)
            
        return pd.DataFrame(new_rows)
        
    except Exception as e:
        print(f"Error in calculate_all_activation_energies: {str(e)}")
        raise


def calculate_activation_energy(row, df_before_ts, conn, ts_dict):
    """
    Calculate activation energy for a transition state
    
    Parameters:
    -----------
    row : Series
        Current row being processed
    df_before_ts : DataFrame
        DataFrame containing species data before transition state
    conn : sqlite3.Connection
        Database connection
    ts_dict : dict
        Dictionary containing transition state data
        
    Returns:
    --------
    Series
        Updated row with activation energy data
    """
    new_row = row.copy()
    species_name = row['species_name']
    status = row['status']
    metal = df_before_ts[(df_before_ts['status'] != 'gas') & (df_before_ts['surface_name'].notna())]['surface_name'].unique()[0]
    facet = df_before_ts[(df_before_ts['status'] != 'gas') & (df_before_ts['site_name'].notna())]['site_name'].unique()[0]
    # Activation Energy
    activation_query = """
        SELECT ae.activation_energy_Ea_eV, ae.activation_energy_Ga_eV, ae.direction, rss.reaction_step_id, ae.pH, r.reference
        FROM ActivationEnergy ae
        JOIN ReactionStepSpecies rss ON ae.reaction_step_id = rss.reaction_step_id
        JOIN Species s ON rss.species_id = s.species_id
        JOIN Status st ON rss.status_id = st.status_id
        JOIN ReactionSteps rs ON ae.reaction_step_id = rs.reaction_step_id
        JOIN ReactionType rt ON rs.reaction_type_id = rt.reaction_type_id
        JOIN Reference r ON ae.reference_id = r.reference_id
        WHERE s.species_name = ? AND st.status = 'ads'
        GROUP BY ae.reaction_step_id
        """

    is_fs_query = """
        SELECT rss.role, s.species_name, st.status, rss.stoichiometry
        FROM ReactionStepSpecies rss
        JOIN Species s ON rss.species_id = s.species_id
        JOIN Status st ON rss.status_id = st.status_id
        WHERE rss.reaction_step_id = ?
        """
    note=""
    if status == 'ts':
        #print(f"{species_name} {status} {metal} {facet} {type(row['frequencies'])}")
        activation_data = pd.read_sql_query(activation_query, conn, params=(species_name,))
        #print(activation_data)

        if not activation_data.empty:
            reference = activation_data['reference'].iloc[0]
            #print(f"reference {reference}")
            #print(f"reaction_step_id {type(activation_data['reaction_step_id'].iloc[0])}")
            is_fs_data = pd.read_sql_query(is_fs_query, conn, params=(activation_data['reaction_step_id'].iloc[0],))
            #print(is_fs_data)
            if str(row['frequencies']) == '[]':
                ae='Ga'
                activation_energy = activation_data['activation_energy_Ga_eV'].iloc[0]
                #print(f"Ga {activation_energy}")
            else:
                ae='Ea'
                activation_energy = activation_data['activation_energy_Ea_eV'].iloc[0]
                #print(f"Ea {activation_energy}")
            direction = activation_data['direction'].iloc[0]
        
        if species_name in ts_dict:
            ae = ts_dict[species_name]['ae']
            activation_energy = ts_dict[species_name]['activation_energy']
            direction = ts_dict[species_name]['direction']
            reaction_step_id = ts_dict[species_name]['reaction_step_id']
            reference = ts_dict[species_name]['reference']
            is_fs_data = pd.read_sql_query(is_fs_query, conn, params=(reaction_step_id,))
        #print(f"\t\t activation energy {species_name} {status} : {activation_energy}")
        #print(is_fs_data)
        #print(f"\t\t activation energy {species_name} {status} : {activation_energy}")
        if direction == 'forward':
            # Store species and their status for IS
            is_species_list = is_fs_data[is_fs_data['role'] == 'IS']['species_name'].tolist()
            is_status_list = is_fs_data[is_fs_data['role'] == 'IS']['status'].tolist()
            is_stoichiometry_list = is_fs_data[is_fs_data['role'] == 'IS']['stoichiometry'].tolist()
            # Calculate formation energy 
            is_energy_sum = 0.0
            for j, (is_species, is_status, is_stoichiometry) in enumerate(zip(is_species_list, is_status_list, is_stoichiometry_list), 1):
                mask = (df_before_ts['species_name']==is_species) & (df_before_ts['status'] == is_status)
                if not df_before_ts[mask].empty:
                    is_formation_energy = df_before_ts[mask]['formation_energy'].iloc[0]
                    is_energy_sum += is_formation_energy*is_stoichiometry
                    #print(f"\t\tIS {is_species} {is_status} : {is_formation_energy}")
                    if is_stoichiometry != 1:
                        note += f"{is_stoichiometry}*{is_species}_{is_status}"
                    else:
                        note += f"{is_species}_{is_status}"
                    if j != len(is_species_list):
                        note += '+'
                else:
                    #print(f"\t\tIS {is_species} {is_status} : not found")
                    pass
            formation_energy = is_energy_sum + activation_energy
            note += f"+({ae}, {activation_energy} eV)"
            #print(f"\t\t formation energy {species_name} {status} : {formation_energy}")
            
        elif direction == 'backward':
            fs_species_list = is_fs_data[is_fs_data['role'] == 'FS']['species_name'].tolist()
            fs_status_list = is_fs_data[is_fs_data['role'] == 'FS']['status'].tolist()
            fs_stoichiometry_list = is_fs_data[is_fs_data['role'] == 'FS']['stoichiometry'].tolist()
            fs_energy_sum = 0.0
            for k, (fs_species, fs_status, fs_stoichiometry) in enumerate(zip(fs_species_list, fs_status_list, fs_stoichiometry_list), 1):
                mask = (df_before_ts['species_name']==fs_species) & (df_before_ts['status'] == fs_status)
                if not df_before_ts[mask].empty:
                    fs_formation_energy = df_before_ts[mask]['formation_energy'].iloc[0]
                    fs_energy_sum += fs_formation_energy*fs_stoichiometry
                    #print(f"\t\tFS {fs_species} {fs_status} : {fs_formation_energy}")
                    if fs_stoichiometry != 1:
                        note += f"{fs_stoichiometry}*{fs_species}_{fs_status}"
                    else:
                        note += f"{fs_species}_{fs_status}"
                    if k != len(fs_species_list):
                        note += '+'
                else:
                    #print(f"\t\tFS {fs_species} {fs_status} : not found")
                    pass
            formation_energy = fs_energy_sum + activation_energy
            note += f"+({ae}, {activation_energy} eV)"
            #print(f"\t\t formation energy {species_name} {status} : {formation_energy}")
        else:
            formation_energy = row['formation_energy']
            #print(f"\t\t formation energy {species_name} {status} : {formation_energy}")

        new_row['formation_energy'] = formation_energy
        new_row['reference'] = reference
        new_row['surface_name'] = metal
        new_row['site_name'] = facet
        new_row['note'] = note
        #print(new_row)
    else:
        # if status is gas, surface_name is 'None' and site_name is 'gas'
        if row['status'] == 'gas':
            new_row['surface_name'] = 'None'
        elif row['status'] == 'slab':
            new_row['surface_name'] = metal
            new_row['site_name'] = facet
    return new_row

def get_efield_coefficients(efield_coeff, row, metal, facet, fallback_map):
    """
    Get electric field coefficients for species
    
    Parameters:
    -----------
    efield_coeff : DataFrame
        DataFrame containing all efield coefficients
    row : Series
        Current row being processed
    metal : str
        Metal name
    facet : str
        Surface facet
    fallback_map : dict
        Dictionary mapping species names to alternative names
        
    Returns:
    --------
    dict or None
        Dictionary containing efield coefficients if found
    """
    try:
        species_name = row['species_name']
        status = row['status']
        reference = row['reference']
        active_site = row['active_site']
        cell_size = row['cell_size']
        
        # First try with original species name
        mask = (
            (efield_coeff['facet'] == facet) &
            (efield_coeff['surface_name'] == metal) &
            (efield_coeff['status'] == status) &
            (efield_coeff['species_name'] == species_name)
        )
        
        coeff_data = efield_coeff[mask]
        
        # If no results, try fallback name
        if coeff_data.empty:
            fallback_name = fallback_map.get(species_name, species_name)
            if fallback_name != species_name:
                mask = (
                    (efield_coeff['facet'] == facet) &
                    (efield_coeff['surface_name'] == metal) &
                    (efield_coeff['status'] == status) &
                    (efield_coeff['species_name'] == fallback_name)
                )
                coeff_data = efield_coeff[mask]
        
        if not coeff_data.empty:
            # Try to match reference
            ref_data = coeff_data[coeff_data['reference'] == reference]
            
            if not ref_data.empty:
                if len(ref_data) == 1:
                    coeff_data = ref_data
                else:
                    # Try to match active_site
                    site_data = ref_data[ref_data['active_site'] == active_site]
                    if len(site_data) == 1:
                        coeff_data = site_data
                    else:
                        # Try to match cell_size
                        cell_data = site_data[site_data['cell_size'] == cell_size]
                        if len(cell_data) == 1:
                            coeff_data = cell_data
            
            # If no exact match found, try 'Ara' reference
            if len(coeff_data) > 1:
                ara_data = coeff_data[coeff_data['reference'] == 'Ara']
                if not ara_data.empty:
                    coeff_data = ara_data.iloc[[0]]
            
            # If still multiple results, take first one
            if len(coeff_data) > 1:
                coeff_data = coeff_data.iloc[[0]]
            
            # Convert to dictionary
            if len(coeff_data) == 1:
                return {
                    'a0': coeff_data['a0'].iloc[0],
                    'a1': coeff_data['a1'].iloc[0],
                    'a2': coeff_data['a2'].iloc[0],
                    'r2': coeff_data['r2'].iloc[0],
                    'reference': coeff_data['reference'].iloc[0],
                    'active_site': coeff_data['active_site'].iloc[0]
                }
        
        return None
        
    except Exception as e:
        print(f"Error in get_efield_coefficients for {row['species_name']}: {str(e)}")
        return None
    
def calculate_efield_formation_energy(formation_energy, efield, coeff_data):
    """Calculate formation energy with electric field correction"""
    return formation_energy + coeff_data['a1']*efield + coeff_data['a2']*efield**2


def apply_efield_corrections(df, efield, metal, facet, conn, fallback_map):
    """
    Apply electric field corrections to formation energies
    
    Parameters:
    -----------
    df : DataFrame
        Input DataFrame containing species data
    efield : float
        Electric field value to apply
    metal : str
        Metal name
    facet : str
        Surface facet
    conn : sqlite3.Connection
        Database connection
    fallback_map : dict
        Dictionary mapping species names to alternative names
    
    Returns:
    --------
    DataFrame
        DataFrame with applied electric field corrections
    """
    new_rows = []
    
    try:
        # Get efield coefficients
        efield_query = """
        SELECT 
            r.reference, 
            f.facet, 
            c.surface_name, 
            st.status, 
            s.species_name, 
            ef.coverage, 
            ef.cell_size, 
            ef.active_site,
            ef.a0, 
            ef.a1, 
            ef.a2, 
            ef.mae, 
            ef.r2
        FROM EfieldCoefficient ef
        JOIN Reference r ON ef.reference_id = r.reference_id
        JOIN Status st ON ef.status_id = st.status_id
        JOIN Species s ON ef.species_id = s.species_id
        JOIN Facets f ON ef.facet_id = f.facet_id
        JOIN Catalysts c ON ef.surface_id = c.surface_id
        """
        
        efield_coeff = pd.read_sql_query(efield_query, conn)
        
        for _, row in df.iterrows():
            new_row = row.copy()
            
            # Skip if gas phase or no formation energy
            if row['status'] == 'gas' or pd.isna(row['formation_energy']):
                new_rows.append(new_row)
                continue
            
            # Get efield coefficients for this species
            coeff_data = get_efield_coefficients(
                efield_coeff, 
                row, 
                metal, 
                facet, 
                fallback_map
            )
            
            if coeff_data is not None:
                # Calculate new formation energy with efield
                formation_energy = calculate_efield_formation_energy(
                    row['formation_energy'],
                    efield,
                    coeff_data
                )
                
                # Update row with new data
                new_row.update({
                    'formation_energy': formation_energy,
                    'note': f"e-field: a1={coeff_data['a1']:.2f}, "
                           f"a2={coeff_data['a2']:.2f}, "
                           f"r2={coeff_data['r2']:.2f} by "
                           f"{coeff_data['reference']} at "
                           f"{coeff_data['active_site']}"
                })
            
            new_rows.append(new_row)
            
        return pd.DataFrame(new_rows)
        
    except Exception as e:
        print(f"Error in apply_efield_corrections: {str(e)}")
        raise
