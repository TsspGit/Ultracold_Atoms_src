__author__ = '@Tssp'
import numpy as np

def load_potential(path):
    with open(path, 'r') as file:
        count = 0
        for line in file:
            if '       ' in line:
                count +=1
            else:
                break
    positives = np.loadtxt(path, max_rows=count, delimiter='       ')
    negatives = np.loadtxt(path, skiprows=count, delimiter='     ')
    return np.concatenate([positives, negatives])