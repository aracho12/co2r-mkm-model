import pickle
import matplotlib.pyplot as plt
from catmap import ReactionModel
import plot_setting
import os
import numpy as np
from catmap import analyze
import math
import argparse
import sys
import pandas as pd
#%matplotlib inline

# Define the parser
parser = argparse.ArgumentParser(description='Plot the data from the CatMAP simulations')
parser.add_argument('-i', '--input', type=str, help='Input folder containing the data')
parser.add_argument('-p', '--potentialvs', action='store_true', help='Plot the data against potentials')
parser.add_argument('-t', '--temperaturevs', action='store_true', help='Plot the data against temperatures')
parser.add_argument('-y', '--ylim', type=float, nargs=2, help='Set the y-axis limits')
parser.add_argument('-x', '--xlim', type=float, nargs=2, help='Set the x-axis limits')
parser.add_argument('-l', '--log', action='store_true', help='Use logarithmic scale for the y-axis')
parser.add_argument('-it', '--intercept', action='store_true', help='Find and plot the intersections between lines')
parser.add_argument('-ty', '--target_y', type=float, help='Find and plot the x-values for the given y-value')
args = parser.parse_args()

potentialvs=True
COR=True
files = os.listdir()
graph_colors=plot_setting.colors

if args.input:
    if os.path.exists(args.input):
        os.chdir(args.input)
        print(f"directory: '{args.input}'")
    else:
        print(f"The input folder '{args.input}' does not exist.")
        sys.exit(1)

# Find file containing '.mkm'
template_file = None
for file in files:
    if '.mkm' in file:
        template_file = file
        break  # Use first file found

if template_file is not None:
    # Check if filename contains 'CO2R'
    if 'CO2R' in template_file:
        COR = False 
    else:
        COR = True
    print(f"Found file: {template_file}, COR: {COR}")
else:
    print("No file found with 'template.mkm.txt' in the name.")

if not potentialvs:
    temperaturevs=True
    potentialvs=False
    print("temperature vs")
else:
    temperaturevs=False
# Define the current directory and the pH value
current_dir = os.getcwd()
#print(f"Current directory: {current_dir}")
pH = 7
# e(elementary charge) 𝜌(surface density of active sites)
eq = 237*1/3*0.001 # uC/cm^2
#print(eq)
num_ele_dict={
    # total number of electron transferred
    'CO_g': 2,
    'H2_g':2,
    'HCOOH_g': 2,
    'HCOO_g': 2,
    'CH4_g':8,
    'CH3OH_g':6,
    'C2H4_g': 12,
    'C2H5OH_g': 12,
    'CH3CH2OH_g':12,
    
}

# Function to get folders starting with 'u_she=' and extract SHE values
def get_she_folders():
    folders = [f for f in os.listdir(current_dir) if f.startswith('u_she=')]
    she_values = sorted([float(f.split('=')[-1]) for f in folders], reverse=True)
    if folders == 0:
        folders = ['.']
        she_values = None
    return folders, she_values

# Get the folders and SHE values
folders_to_process, she_values = get_she_folders()
if not she_values:
    folders_paths = ['.']
    she_values = [0]
else:
    folders_paths = [os.path.join(current_dir, f'u_she={she}') for she in she_values]
log_file=[f for f in os.listdir(folders_paths[0]) if '.log' in f]
# template_file의 확장자 제거 한 이름
template_file_name = template_file.replace('.mkm', '')
log_file=f'{template_file_name}.log'
#print(f"Folders to process: {folders_to_process}")
#print(f"Folders paths: {folders_paths}")

# Calculate the interval `d`
if len(she_values) > 1:
    d = abs(she_values[1] - she_values[0])
else:
    d = 0  # In case there's only one folder

def get_data_labels(data_type, folder):
    setup_file = os.path.join(folder, log_file)
    model = ReactionModel(setup_file=setup_file)

    if data_type == 'coverage_map':
        data_map = model.coverage_map
        label_list=model.output_labels['coverage']
    elif data_type == 'production_rate_map':
        data_map = model.production_rate_map
        label_list=model.output_labels['production_rate']
        
    elif data_type in ['current_density', 'faradaic_efficiency', 'my_selectivity_map']:
        if COR:
            excluded_labels = {'O2_g', 'ele_g', 'H2O_g', 'H_g','OH_g','CO2_g', 'CO_g'}
        else:
            excluded_labels = {'O2_g', 'ele_g', 'H2O_g', 'H_g','OH_g','CO2_g'}
        all_labels = model.output_labels['production_rate']
        label_list = [label for label in all_labels if label not in excluded_labels]           
        
        possible_prdt = label_list
        prdt_idx = {prdt: model.output_labels['production_rate'].index(prdt) for prdt in possible_prdt}
        current_density_list = []
        partial_current_density_list = []
        faradaic_efficiency_list = []
        my_selectivity = []

        for descri, rate in model.production_rate_map:
            total_current_density = 1.e-99
            total_rate = 1.e-99
            for prdt in possible_prdt:
                current_density = eq * num_ele_dict[prdt]*rate[prdt_idx[prdt]]
                total_current_density += current_density
                current_density_list.append(current_density)
                total_rate += rate[prdt_idx[prdt]]
            
            selectivity = [rate[prdt_idx[prdt]] / total_rate for prdt in possible_prdt]
            my_selectivity.append([descri, selectivity])

            partial_current_density = [eq * num_ele_dict[prdt]*rate[prdt_idx[prdt]] for prdt in possible_prdt]
            partial_current_density_list.append([descri, partial_current_density])

            faradaic_efficiency = [eq * num_ele_dict[prdt]*rate[prdt_idx[prdt]] / total_current_density for prdt in possible_prdt]
            faradaic_efficiency_list.append([descri, faradaic_efficiency])

        if data_type == 'current_density':
            data_map = partial_current_density_list
        elif data_type == 'faradaic_efficiency':
            data_map = faradaic_efficiency_list
        elif data_type == 'my_selectivity_map':
            data_map = my_selectivity

    return data_map, label_list

def format_label(label):
    label = label.replace('_g', '(g)')
    label = label.replace('_t', '*')
    formatted_label = ""
    for char in label:
        if char.isdigit():
            formatted_label += f"$_{char}$"
        else:
            formatted_label += char
    return formatted_label



def line_equation_from_points(x1, y1, x2, y2, log=False):
    """
    Calculates the slope and intercept of the line passing through two points.
    
    :param x1, y1: First point
    :param x2, y2: Second point
    :return: (slope, intercept) if the line is not vertical, otherwise None
    """
    if x1 == x2:
        # Vertical line
        return None
    
    if log:
        m = np.log10(y2/y1) / (x2 - x1)
        b = np.log10(y1) - m * x1
    else:
        m = (y2 - y1) / (x2 - x1)
        b = y1 - m * x1
    
    return (m, b) # slope, intercept

def find_intersection_of_lines(slope1, intercept1, slope2, intercept2, log=False):
    """
    Finds the intersection point of two lines given their slopes and y-intercepts.

    :param slope1: Slope of the first line
    :param intercept1: Y-intercept of the first line
    :param slope2: Slope of the second line
    :param intercept2: Y-intercept of the second line
    :return: (x-coordinate, y-coordinate) of the intersection point, or None if no intersection exists
    """
    # If the slopes are the same, the lines are parallel or identical
    if slope1 == slope2:
        return None  # No intersection
    
    # Calculate the x-coordinate of the intersection point
    x_intersection = (intercept2 - intercept1) / (slope1 - slope2)
    # Calculate the x-coordinate of the intersection point
    if log:
        log_y_intersection = slope1 * x_intersection + intercept1
        y_intersection = 10**(log_y_intersection)
    else:
        # Regular linear case
        y_intersection = slope1 * x_intersection + intercept1
    
    return x_intersection, y_intersection    


def find_intersection_of_segments(x1, y1, x2, y2, x3, y3, x4, y4, log=False):
    """
    Finds the intersection of two line segments.
    
    :param x1, y1, x2, y2: First line segment (x1, y1) to (x2, y2)
    :param x3, y3, x4, y4: Second line segment (x3, y3) to (x4, y4)
    :return: (x, y) coordinates of the intersection, or None if they don't intersect within the segments
    """
    # Get the equation of the first line
    line1 = line_equation_from_points(x1, y1, x2, y2, log)
    if line1 is None:
        return None
    m1, b1 = line1
    
    # Get the equation of the second line
    line2 = line_equation_from_points(x3, y3, x4, y4, log)
    if line2 is None:
        return None
    m2, b2 = line2
    
    # Check if the lines are parallel
    if m1 == m2:
        return None
    
    # Calculate the x coordinate of the intersection
    x_intersection, y_intersection = find_intersection_of_lines(m1, b1 , m2, b2, log)
    if min(x1, x2) <= x_intersection <= max(x1, x2):
        #y_intersection = 10**(m1 * x_intersection + np.log10(b1))
        return (x_intersection, y_intersection)    

    
    return None

def find_intersections_between_labels(x, y_list, log=False):
    """
    Finds intersections between line segments defined by multiple labels.
    
    :param x: List of x values (shared across all labels)
    :param y_list: List of (label, y values) for each label
    :return: List of intersections (label1, label2, x_intersection, y_intersection)
    """
    intersections = []
    
    # Loop through all pairs of labels
    for i in range(len(y_list)):
        label1, y_values1 = y_list[i]
        for j in range(i + 1, len(y_list)):
            label2, y_values2 = y_list[j]
            print(label1,y_values1,label2,y_values2)
            
            # Check segments between each consecutive x point
            for k in range(len(x) - 1):
                x1, x2 = x[k], x[k + 1]
                y1_label1, y2_label1 = y_values1[k], y_values1[k + 1]
                y1_label2, y2_label2 = y_values2[k], y_values2[k + 1]
                
                # Find intersection between the segments (x1, y1_label1) -> (x2, y2_label1)
                # and (x1, y1_label2) -> (x2, y2_label2)
                #print(f"Checking intersection between {label1} and {label2} at x={x1:.2f} to {x2:.2f}")
                intersection = find_intersection_of_segments(
                    x1, y1_label1, x2, y2_label1,
                    x1, y1_label2, x2, y2_label2,
                    log
                )
                print(intersection)
                if intersection:
                    intersections.append((label1, label2, intersection[0], intersection[1]))
    
    return intersections

def find_x_for_target_y(x1, y1, x2, y2, target_y=1e-4, log=False):
    """
    Finds the x-value for a given y-value (target_y) along the line segment from (x1, y1) to (x2, y2).
    
    :param x1, y1: Coordinates of the first point on the segment
    :param x2, y2: Coordinates of the second point on the segment
    :param target_y: The y-value for which we want to find the corresponding x-value
    :param log: If True, perform the calculations in log scale
    :return: x-coordinate where the line segment intersects the target_y, or None if no solution exists
    """
    # If the y-values don't bracket the target_y, return None
    if not (min(y1, y2) <= target_y <= max(y1, y2)):
        return None
    
    # Get the slope and intercept of the line segment
    line = line_equation_from_points(x1, y1, x2, y2, log)
    if line is None:
        return None  # Vertical line case, no valid intersection with y=target_y
    
    m, b = line
    
    # Calculate the x-value for y = target_y
    if log:
        target_y = np.log10(target_y)
    x_target = (target_y - b) / m
    
    # Check if the x-value is within the segment range
    if min(x1, x2) <= x_target <= max(x1, x2):
        return x_target
    else:
        return None

def find_x_for_target_y_in_labels(x, y_list, target_y=1e-4, log=False):
    """
    Finds the x-values where the line segments defined by multiple labels cross the target_y.
    
    :param x: List of x values (shared across all labels)
    :param y_list: List of (label, y values) for each label
    :param target_y: Target y-value for which to find x-values
    :return: List of (label, x-value) where the line segment crosses the target_y
    """
    x_for_target_y = []
    
    # Loop through each label's y-values
    for label, y_values in y_list:
        for k in range(len(x) - 1):
            x1, x2 = x[k], x[k + 1]
            y1, y2 = y_values[k], y_values[k + 1]
            
            # Find the x-value for the given target_y on this segment
            x_target = find_x_for_target_y(x1, y1, x2, y2, target_y, log)
            
            if x_target is not None:
                x_for_target_y.append((label, x_target))
    
    return x_for_target_y  



# Function to load and filter data based on the adjusted VRHE range
def collect_data(folders, she_values, d, pH, data_type):
    
    combined_data = {}
    combined_potentials = []
    ranges = []

    for i, folder in enumerate(folders):
        os.chdir(folder)
        she_value = she_values[i]
        rhe_value = she_value + 0.059 * pH

        if i < len(folders) - 1:
            lower_bound = rhe_value - 0.5 * d
            upper_bound = rhe_value + 0.5 * d
        else:  # Last folder
            lower_bound = rhe_value - 10*d
            upper_bound = rhe_value + 0.5 * d  # Extending to the end

        if len(folders) == 1:
            if temperaturevs:
                # temperature
                lower_bound, upper_bound = 0,1000
            else:
                lower_bound, upper_bound = -2, 2


        ranges.append((lower_bound, upper_bound))  
        pkl_file_path = os.path.join(folder, log_file)

        if os.path.exists(pkl_file_path):
            
            model = ReactionModel(setup_file=log_file)
            vm = analyze.VectorMap(model)
            vm.get_pts_cols()

            if i == 0:
                data_map, label_list = get_data_labels(data_type, folder)
                combined_data = {label: [] for label in label_list}
            else:
                data_map, _ = get_data_labels(data_type, folder)


            # Filter data within the adjusted VRHE range
            filtered_data = [
                item for item in data_map
                if lower_bound <= item[0][0] <= upper_bound
            ]

            if filtered_data:
                # Extract potentials and data values from filtered data
                potentials = [item[0][0] for item in filtered_data]
                #print(potentials)
                data_values = [item[1] for item in filtered_data]

                # print(f'Potentials: {potentials}')
                # print(f'Data values: {data_values}')

                combined_potentials.extend(potentials)

                # Combine the data for plotting
                for j, label in enumerate(label_list):
                    values = [data_value[j] for data_value in data_values]
                    combined_data[label].extend(values)
            else:
                print(f'No data in the ±0.5*d VRHE range for {folder}')
            
        else:
            print(f'log file not found in {folder}')
        os.chdir(current_dir)

    return combined_potentials, combined_data, ranges

def save_data_to_csv(data_type, sorted_potentials, sorted_combined_data):
    """
    Function to save data to a CSV file
    
    :param data_type: Type of data (coverage_map, production_rate_map etc.)
    :param sorted_potentials: List of sorted potential values
    :param sorted_combined_data: Dictionary of sorted data
    """
    # Prepare data for DataFrame
    df_data = {'potential_or_temperature': sorted_potentials}
    
    # Add data for each label
    for label, values in sorted_combined_data.items():
        formatted_label = format_label(label)
        df_data[formatted_label] = values
    
    # Create DataFrame
    df = pd.DataFrame(df_data)
    
    # Set x-axis label
    if potentialvs:
        x_label = 'Potential_vs_RHE(V)'
    else:
        x_label = 'Temperature(K)'
    df = df.rename(columns={'potential_or_temperature': x_label})
    
    # Save to CSV file
    output_filename = f'{data_type}_data.csv'
    df.to_csv(output_filename, index=False)
    print(f"Data saved to {output_filename}")
    
# Function to plot the combined data
def plot_combined_data(data_type,
                       ylim=None,
                       xlim=None,
                       xlabel=None,
                       ylabel=None,
                       log=False,
                       save_title=None):

    try:
         # Define colors and markers for different folders
        if len(folders_paths) == 1:
            colors = ['white']
        else:
            colors = plt.cm.viridis(np.linspace(0, 1, len(folders_paths)))
        markers = ['o', 's', '^', 'D', 'v', '<', '>']
        # Collect the combined data
        combined_potentials, combined_data, ranges = collect_data(folders_paths, she_values, d, pH, data_type)
        
        sorted_indices = sorted(range(len(combined_potentials)), key=lambda k: combined_potentials[k])
        sorted_potentials = [float(combined_potentials[i]) for i in sorted_indices]
        label_list = list(combined_data.keys())
        sorted_combined_data = {label: [float(combined_data[label][i]) for i in sorted_indices] for label in label_list}
        
        # Save data to CSV file
        save_data_to_csv(data_type, sorted_potentials, sorted_combined_data)
        
        # Plot the combined data
        plt.figure()
        for i, (label, values) in enumerate(sorted_combined_data.items()):
            # Filter out values below 0.01
            if data_type in ['coverage_map', 'production_rate_map']:
                if any(value >= 0.001 for value in values):
                    formatted_label = format_label(label)
                    plt.plot(sorted_potentials, values, label=formatted_label)
                    plt.scatter(sorted_potentials, values, edgecolors='none')
            else:
                formatted_label = format_label(label)
                plt.plot(sorted_potentials, values, label=formatted_label)
                plt.scatter(sorted_potentials, values, edgecolors='none')
        

        if args.intercept:
            print("Finding intersections")
            if data_type == 'current_density':
                # Find and plot intersections
                y_list = [(label, sorted_combined_data[label]) for label in sorted_combined_data.keys()]
                intersections = find_intersections_between_labels(sorted_potentials, y_list, log)
                #print(f"Intersections: {intersections}")
                for label1, label2, x, y in intersections:
                    plt.plot(x, y, 'ro')
                    #plt.annotate(f'({x:.2f}, {y:.2e})', (x, y), textcoords="offset points", 
                    #            xytext=(0,10), ha='center', fontsize=8, color='red')
                    print(f"Intersection between {label1} and {label2}: x={x:.2f}, y={y:.2e}")

        if args.target_y:
            if data_type=='current_density':
                target_y = args.target_y
                print(f"Finding x-values for y={target_y}")
                plt.axhline(y=target_y, color='gray', linestyle='--', linewidth=0.5)
                y_list = [(label, sorted_combined_data[label]) for label in sorted_combined_data.keys()]

                x_for_target_y = find_x_for_target_y_in_labels(sorted_potentials, y_list, target_y, log)
                for label, x in x_for_target_y:
                    plt.plot(x, target_y, 'ro')
                    #plt.annotate(f'({x:.2f}, {target_y:.2e})', (x, target_y), textcoords="offset points", 
                    #            xytext=(0,10), ha='center', fontsize=8, color='red')
                    #print(f"Label: {label}, x={x:.2f}")
                # horizontal line indicating the target y-value
                
                
                # data
                df = pd.DataFrame(x_for_target_y, columns=['species', 'potential (V RHE)'])
                df['potential (V RHE)'] = df['potential (V RHE)'].round(2)
                print(df)



        # Add labels and legend
        if xlabel:
            plt.xlabel(xlabel)
        else:
            if potentialvs:
                plt.xlabel('$E$ vs $RHE (V)$')
            if temperaturevs:
                plt.xlabel('Temperature (K)')

        if ylabel:
            plt.ylabel(ylabel)
        else:
            plt.ylabel(data_type.replace('_', ' ').capitalize())
        if ylim:
            plt.ylim(ylim)

        if xlim:
            plt.xlim(xlim)
        else:
            if potentialvs:
                if args.xlim:
                    plt.xlim(args.xlim[0],args.xlim[1])
                else:
                    plt.xlim(-1.4,-0.1)
            if temperaturevs:
                if args.xlim:
                    plt.xlim(args.xlim[0],args.xlim[1])
                else:
                    plt.xlim(270,415)
        if log:
            plt.yscale('log')
        
        if len(folders_paths) == 1:    
            plt.title(f'{data_type.replace("_", " ").capitalize()} for pH {pH}')
        else:
            plt.title(f'{data_type.replace("_", " ").capitalize()} (Combined) for pH {pH}')

        # Add annotations for each range
        for (lower_bound, upper_bound), color in zip(ranges, colors):
            plt.axvspan(lower_bound, upper_bound, color=color, alpha=0.09)
            # plt.text((lower_bound + upper_bound) / 2, plt.ylim()[1] * 0.9, f'{lower_bound:.2f} to {upper_bound:.2f}', 
            #          horizontalalignment='center', color=color, fontsize=8)

        plt.legend()
        #plt.grid(True)

        # Save the plot as an image file
        if save_title:
            plt.savefig(f'{save_title}.png')
        else:
            plt.savefig(f'combined_{data_type}.png')
        plt.show()
    except Exception as e:
        print(f"An error occurred: {e}")
    



plot_combined_data('coverage_map',ylim=(0,1),ylabel='Coverage')
plot_combined_data('production_rate_map', ylabel='Rate',log=True, ylim=(1e-10, 1e15))
plot_combined_data('my_selectivity_map',ylim=(0,1))
plot_combined_data('current_density',ylabel='$j_{partial}$ $(mA/Cm^2)$',ylim=(0,1))
plot_combined_data('current_density',ylabel='$j_{partial}$ $(mA/Cm^2)$',log=True, ylim=(1e-10, 1e10),save_title='combined_current_density_Log')
plot_combined_data('faradaic_efficiency',ylabel='FE',ylim=(0,1))
#plot_combined_data('tafel_plot',ylabel='FE',ylim=(0,1))
