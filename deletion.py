import numpy as np
import vtk
import pickle
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay

def save_to_pickle(array, filename):
    with open(filename, "wb") as f:
        pickle.dump(array, f)
        
def load_from_pickle(filename):
    with open(filename, "rb") as f:
        loaded = pickle.load(f)
    return loaded

new_mesh = load_from_pickle("remeshed_groups.pkl")
mesh = new_mesh[2]
print(mesh)

