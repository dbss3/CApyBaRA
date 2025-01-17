import numpy as np
from read_in_for_2DCARMA import read_in

def just_plot_it(fresh_readin=True):

    if fresh_readin:
        saved_dict_names_list = read_in()
    elif not fresh_readin:
        np.load(saved_dict_name)
    else:
        print(f'fresh_readin must be a boolean, currently: {fresh_readin}')
    return