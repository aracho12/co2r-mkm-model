from core.formation_energy_calculator import generate_formation_energies
ts_dict = {
    'H-ele': {
        'activation_energy': 0.80, 
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 1,
        'ae': 'Ea'
    },
    'H2-ele': {
        'activation_energy': 0.76,
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 2,
        'ae': 'Ea'
    },
    'CO2-H-ele': {
        'activation_energy': 0.60, 
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 32,
        'ae': 'Ga'
    },
    'COOH-H-ele': {
        'activation_energy': 0.40, 
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 33,
        'ae': 'Ga'
    },
    'H-CO2-ele': {
        'activation_energy': 0.7, 
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 34,
        'ae': 'Ga'
    },
    'HCOO-H-ele': {
        'activation_energy': 0.40, 
        'direction': 'forward',
        'reference': 'This work',
        'reaction_step_id': 35,
        'ae': 'Ga'
    },
}

correction_dict = {
    'COH': -0.15,
    'CO-H-ele': -0.1,
}

# Example usage:
DB_NAME = 'co2r-suncat.db'
save_dir = '.'

if __name__ == "__main__":
    # Generate for Cu with multiple facets
    generate_formation_energies(
        db_path=DB_NAME,
        metal='Cu',
        facets=['100','211'],
        she_range=(-2.00, -0.40, 0.02),
        if_efield=True,
        ts_dict=ts_dict,
        correction_dict=correction_dict,
        save_dir=save_dir
    )
