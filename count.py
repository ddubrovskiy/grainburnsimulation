import numpy
from stl import mesh
import math

def parse_stl(filename):
    stl_mesh = mesh.Mesh.from_file(filename)
    num_faces = stl_mesh.vectors.shape[0]
    vertices = stl_mesh.vectors
    return num_faces, vertices

n, coords = parse_stl("complicated_grain_mesh.stl")

normals = []
coords_n = numpy.array(coords)
all_coords = []
for coord in coords:
    for c in coord: 
        all_coords.append(c)
        
unique_coords = numpy.unique(all_coords, axis = 0)

main_array = numpy.array([[1.5, 2.5, 3.5], [4.5, 5.5, 6.5], [7.5, 8.5, 9.5], [1.5, 2.5, 3.6], [4.5, 5.5, 6.5], [5.5, 3.5, 2.5], [1.5, 2.5, 3.5], [4.5, 5.5, 6.5], [1.3, 1.2, 4.1]])
sub_arrays = numpy.array([[1.5, 2.5, 3.5], [4.5, 5.5, 6.5]])

result = {}
bbl = []
for i in range(len(sub_arrays)):
    matches = numpy.all(main_array == sub_arrays[i], axis=1)
    matching_indexes = numpy.where(matches)[0] 
    bbl.append(matching_indexes)
    
print(bbl)
for i in range(len(bbl)):
    for j in range(len(bbl[i])):
        bbl[i][j] = float(bbl[i][j])
        bbl[i][j] = bbl[i][j]/3
        bbl[i][j] = math.ceil(bbl[i][j])

fi = False
if(fi == True):
    print("Indices:")
    print(bbl)
    print("------------------")


nnn = numpy.array([numpy.zeros(3), numpy.zeros(3)])
ns = numpy.array([[0.1, 0.2, 0.3], [0.2, 0.3, 0.4], [0.3, 0.4, 0.5]])
for i in range(len(sub_arrays)):
    for j in range(len(bbl[i])):
        nnn[i]+=ns[bbl[i][j]]
    nnn[i]/=len(bbl[i])

fm = False
if(fm == True):
    print("Mean normals:")
    print(nnn)
    print("------------------")


for i in range(len(sub_arrays)):
    matches = numpy.all(main_array == sub_arrays[i], axis = 1)
    matching_indexes = numpy.where(matches)[0]
    for m in matching_indexes:
        main_array[m]-=nnn[i]
   
fr = False     
if(fr == True):
    print("Result:")
    print(main_array)
    print("------------------")  

 
#print(sum(counts.values())) #43182
#print(len(all_coords))      #43182
    

for coord in coords_n:
    v1 = coord[1] - coord[0]
    v2 = coord[2] - coord[1]
    normals.append(numpy.cross(v1, v2))
    

counter = 0
for n in normals:
    counter+=1
    l = numpy.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
    n[0]/=l
    n[1]/=l
    n[2]/=l



index = []
for i in range(len(unique_coords)):
    matches = numpy.all(all_coords == unique_coords[i], axis=1)
    matching_indexes = numpy.where(matches)[0]
    index.append(matching_indexes)

for i in range(len(index)):
    for j in range(len(index[i])):
        index[i][j] = index[i][j]
        index[i][j] = index[i][j]/3
        index[i][j] = math.ceil(index[i][j])
        
normals_mean = numpy.zeros((len(unique_coords), 3))
for i in range(len(unique_coords)):
    for j in range(len(index[i])):
        normals_mean[i]+=normals[index[i][j]]
    normals_mean[i]/=len(index[i])

print(unique_coords[0])
print(index[0])
print("----------------------------")
for j in index[0]:
    print(normals[j])
    print("=")
    print(coords[j])
    print("--------------------------------")
print("Mean:")
print(normals_mean[0])

burn = all_coords
rate = 10
for i in range(len(unique_coords)):
    matches = numpy.all(burn == unique_coords[i], axis = 1)
    matching_indexes = numpy.where(matches)[0]
    for m in matching_indexes:
        burn[m]-=(rate*normals_mean[i])
        
numpy.save("iteration_0.npy", burn)
