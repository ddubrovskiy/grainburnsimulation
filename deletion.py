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
collisions = load_from_pickle("collisions_final.pkl")

print(len(collisions), len(new_mesh))
mesh = new_mesh[2]
print(mesh)
print("-------------------")
print(collisions[2])
point = collisions[2][0]
print("-----------------------")
for coord in collisions[2]:
    z = np.copy(coord[2])
    z_coord = np.array([0, 0, z])
    print(z_coord, coord)
    print(np.linalg.norm(coord - z_coord))
    
    
#getting a triangle by node that is within this triangle can be achieved by simple if a in s: flag = True logic
    