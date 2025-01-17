import numpy as np
import sys
from read_in_routines import read_in, read_file_to_dict

# Function that does everything
def main():
    """
    Main function to do everything
    """
    if len(sys.argv) != 2:  # Ensure exactly one argument is provided
        print("Usage: python3 just_plot_it.py <file_path>")
        sys.exit(1)

    file_path = sys.argv[1]  # Get the file path from the command-line arguments

    try:
        result = read_file_to_dict(file_path)
        print(result)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()

def just_plot_it(fresh_readin=True):


    # Check if folder where the plots are generated exists; if 
    # not, make one. 
    if not os.path.isdir(path):
        os.makedirs(path)

    # First lets read in the file names and locations
    read_file_to_dict(file_path)

    if fresh_readin:
        saved_dict_names_list = read_in()
    elif not fresh_readin:
        np.load(saved_dict_name)
    else:
        print(f'fresh_readin must be a boolean, currently: {fresh_readin}')
    return