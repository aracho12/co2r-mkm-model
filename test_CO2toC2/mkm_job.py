from for_a_happy_life import plot_setting
from catmap import ReactionModel
import sys
from string import Template
import os
import pandas as pd
from matplotlib.table import Table
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import pickle
import numpy as np
from mpmath import mpf
import subprocess
include_rate_control = False
FED = True

for s in [1.0]:
    # Find .mkm file in the current directory
    mkm_path = [f for f in os.listdir('.') if f.endswith('.mkm')][0]
    model = ReactionModel(setup_file = mkm_path)
    model.output_variables+=['production_rate', 'free_energy', 'selectivity', 'interacting_energy', 'interaction_matrix']
    if include_rate_control:
        model.output_variables += ['rate_control']
    model.run(recalculate=True)

    if 'COR' in mkm_path:
        possible_prdt = ('CH4_g', 'C2H4_g', 'H2_g')
    elif 'CO2R' in mkm_path:
        if 'CO2toCO' in mkm_path:
            possible_prdt = ('CO_g', 'HCOOH_g', 'H2_g')
        else:
            possible_prdt = ('CO_g', 'HCOOH_g', 'H2_g', 'CH4_g', 'C2H4_g')

    model.output_labels['my_selectivity'] = possible_prdt
    prdt_idx = {}
    for prdt in possible_prdt:
        prdt_idx[prdt] = model.output_labels['production_rate'].index(prdt)
    my_selectivity = []
    for descri, rate in model.production_rate_map:
        total = 1.e-99
        for prdt in possible_prdt:
            total += rate[prdt_idx[prdt]]
        selectivity = [rate[prdt_idx[prdt]]/total for prdt in possible_prdt]
        my_selectivity.append([descri, selectivity])
    model.my_selectivity_map = my_selectivity

    from catmap import analyze
    vm = analyze.VectorMap(model)

    # Plot reaction rates
    vm.plot_variable = 'rate'
    vm.log_scale = True
    vm.min = 1e-10
    vm.max = 1e6
    fig = vm.plot(save=False)
    fig.savefig('rate.pdf')

    # Plot production rates
    vm.plot_variable = 'production_rate'
    vm.log_scale = True
    vm.min = 1e-10
    vm.max = 1e6
    fig = vm.plot(save=False)
    fig.savefig('production_rate'+str(s)+'.pdf')

    # Plot surface coverage (linear scale)
    vm = analyze.VectorMap(model)
    vm.log_scale = False
    vm.unique_only = False
    vm.plot_variable = 'coverage'
    vm.min = 0
    vm.max = 1
    fig = vm.plot(save=False)
    fig.savefig('coverage'+str(s)+'.pdf')

    # Plot surface coverage (log scale)
    vm = analyze.VectorMap(model)
    vm.log_scale = True
    vm.unique_only = False
    vm.plot_variable = 'coverage'
    vm.min = 1e-15
    vm.max = 1
    fig = vm.plot(save=False)
    fig.savefig('coverageLog'+str(s)+'.pdf')

    # Plot selectivity (linear scale)
    vm = analyze.VectorMap(model)
    vm.log_scale = False
    vm.plot_variable = 'my_selectivity'
    vm.include_labels = possible_prdt
    vm.min = 0
    vm.max = 1
    fig = vm.plot(save=False)
    fig.savefig('my_selectivity'+str(s)+'.pdf')

    # Plot selectivity (log scale)
    vm = analyze.VectorMap(model)
    vm.plot_variable = 'my_selectivity'
    vm.include_labels = possible_prdt
    vm.min = 0
    vm.max = 1
    vm.log_scale = True
    fig = vm.plot(save=False)
    fig.savefig('my_selectivityLog'+str(s)+'.pdf')

    if include_rate_control:
        mm = analyze.MatrixMap(model)
        mm.plot_variable = 'rate_control'
        mm.log_scale = False
        mm.min = -2
        mm.max = 2
        mm.plot(save='rate_control.pdf')

# Mechanism analysis
ma = analyze.MechanismAnalysis(model)
ma.energy_type = 'free_energy'
label_size = 10
ma.kwarg_dict = {
        'CO2-CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2-HCOO':{'label_positions':None,'label_args':{'color':'b','rotation':45,'ha':'center','size':label_size}},
        'C1_via_CO-H-ele':{'label_positions':None,'label_args':{'color':'b','rotation':45,'ha':'center','size':label_size}},
        'C1_via_H-CO':{'label_positions':None,'label_args':{'color':'c','rotation':45,'ha':'center','size':label_size}},
        'C2_via_OCCOH':{'label_positions':None,'label_args':{'color':'g','rotation':45,'ha':'center','size':label_size}},
        'C2_via_C-CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'C2_via_CH-CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'C2_via_CH2-CH2':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'HER_V-H':{'label_positions':None,'label_args':{'color':'g','rotation':45,'ha':'center','size':label_size}},
        'HER_Tafel':{'label_positions':None,'label_args':{'color':'g','rotation':45,'ha':'center','size':label_size}},
        'CO2R_CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_HCOOH':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_C1_via_CO-H-ele':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_C1_via_H-CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_C2_via_OCCOH':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_C2_via_C-CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_C2_via_CH-CO':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        'CO2R_C2_via_CH2-CH2':{'label_positions':None,'label_args':{'color':'k','rotation':45,'ha':'center','size':label_size}},
        }

if FED:
    for voltage in [-0.7, -0.9, -1.1]:
        ma.subplots_adjust_kwargs = {'top': 0.87, 'bottom': 0.22}
        ma.include_labels = True
        ma.pressure_correction = True
        ma.coverage_correction = True

        # Create figure with mechanism analysis plot
        fig = ma.plot(save=False, plot_variants=[voltage])
        fig.set_figwidth(7.5)
        fig.set_figheight(10)  # Increase total height

        # Create GridSpec for subplots
        gs = GridSpec(2, 1, height_ratios=[2, 1])  # 2 rows, 1 column (graph + table)

        # Adjust position of original plot
        ax = fig.axes[0]  # Get first axis
        ax.set_position(gs[0].get_position(fig))
        ax.set_subplotspec(gs[0])

        # Add subplot area for tables
        table_ax = fig.add_subplot(gs[1])
        table_ax.axis('off')  # Hide subplot axis

        tables = []
        for i, (key, value) in enumerate(ma.data_dict.items()):
            # Match lengths of ∆G and Ga lists
            length = max(len(value[0]), len(value[1]))
            delta_g = value[0] + [None] * (length - len(value[0]))
            ga = value[1] + [None] * (length - len(value[1]))

            df = pd.DataFrame({
                '∆G': [round(val, 2) if val is not None else None for val in delta_g],
                'Ga': [round(val, 2) if val is not None else None for val in ga]
            })

            # Create individual table
            table = plt.table(cellText=df.values,
                              colLabels=df.columns,
                              cellLoc='center',
                              loc='center',
                              bbox=[0.33 * i, 0, 0.3, 1])  # Adjust table position and size
            table.auto_set_font_size(False)
            table.set_fontsize(10)
            tables.append(table)

            # Add title above table
            table_ax.text(0.33 * i + 0.15, 1.05, key, ha='center', va='bottom', fontsize=12, weight='bold')

        fig.savefig(f'FED_int_{voltage}_with_table.pdf')

# Run mkm_plotter.py after completing all calculations
plotter_path = "/Users/aracho/bin/for_a_happy_life/mkm/old/mkm_plotter.py"
subprocess.run(['python', plotter_path], check=True)