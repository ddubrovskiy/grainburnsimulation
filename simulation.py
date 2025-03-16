import numpy as np
from stl import mesh
import vtk
import pickle

def save_to_pickle(array, filename):
    with open(filename, "wb") as f:
        pickle.dump(array, f)
        
def load_from_pickle(filename):
    with open(filename, "rb") as f:
        loaded = pickle.load(f)
    return loaded

def nonzero(vec, eps = 1e-9):
    for i in range(len(vec)):
        if ((vec[i] > eps or vec[i] < eps) and vec[i]!=0):
            return i
 
def getNorm(Tr):
    v1 = Tr[1] - Tr[0]
    v2 = Tr[2] - Tr[1]
    n = np.cross(v1, v2)
    l = np.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
    n[0]/=l
    n[1]/=l
    n[2]/=l
    return n

def max_abs_index(array):
    max_value = array[0]
    max_index = 0
    
    for i in range(len(array)):
        if np.abs(array[i]) > np.abs(max_value):
            max_value = array[i]
            max_index = i
    
    return max_index

def parse_stl(filename):
    stl_mesh = mesh.Mesh.from_file(filename)
    num_faces = stl_mesh.vectors.shape[0]
    vertices = stl_mesh.vectors
    return num_faces, vertices

def different_sign(D):
    a = 0
    p_c = 0
    n_c = 0
    for i in range(len(D)):
        if D[i]>0: p_c +=1
        else: n_c +=1
    if p_c > n_c:
        for i in range(len(D)):
            if D[i] < 0:
                a = i
    else:
        for i in range(len(D)):
            if D[i] > 0:
                a = i
    return a

def name_ind(a):
    if a == 0:
        b = 1
        c = 2
    elif a == 1:
        b = 2
        c = 0
    elif a == 2:
        b = 0
        c = 1
    return b, c 
    
def check_position(triangle):
    for coord in triangle:
        if(np.abs(np.sqrt(coord[0]**2 + coord[1]**2) - 122) > 0.1):
            return False
    return True
    
def check_collisions(Tr1, Tr2):
    A_min = np.min(Tr1, axis=0)
    A_max = np.max(Tr1, axis=0)
    B_min = np.min(Tr2, axis=0)
    B_max = np.max(Tr2, axis=0)
    if(B_min[0] > A_max[0] or B_min[1] > A_max[1] or B_min[2] > A_max[2] or A_min[0] > B_max[0] or A_min[1] > B_max[1] or A_min[2] > B_max[2]):
        return 0
    else:
        Y1 = np.array([Tr2[0] - Tr1[0], Tr2[0] - Tr1[1], Tr2[0] - Tr1[2]])
        Y2 = np.array([Tr2[1] - Tr1[0], Tr2[1] - Tr1[1], Tr2[1] - Tr1[2]])
        Y3 = np.array([Tr2[2] - Tr1[0], Tr2[2] - Tr1[1], Tr2[2] - Tr1[2]])
        Y4 = np.array([Tr1[0] - Tr2[0], Tr1[0] - Tr2[1], Tr1[0] - Tr2[2]])
        Y5 = np.array([Tr1[1] - Tr2[0], Tr1[1] - Tr2[1], Tr1[1] - Tr2[2]])
        Y6 = np.array([Tr1[2] - Tr2[0], Tr1[2] - Tr2[1], Tr1[2] - Tr2[2]])
        det1 = np.linalg.det(Y1)
        det2 = np.linalg.det(Y2)
        det3 = np.linalg.det(Y3)
        det4 = np.linalg.det(Y4)
        det5 = np.linalg.det(Y5)
        det6 = np.linalg.det(Y6)
        if(det1 > 0 and det2 > 0 and det3 > 0):
            return 0
        elif(det1 < 0 and det2 < 0 and det3 < 0):
            return 0
        elif(det4 > 0 and det5 > 0 and det6 > 0):
            return 0
        elif(det4 < 0 and det5 < 0 and det6 < 0):
            return 0
        elif(det1 == 0 or det2 == 0 or det3 == 0):
            return 0
        elif(det4 == 0 or det5 == 0 or det6 == 0):
            return 0
        else:          
            D1 = [det1, det2, det3]
            D2 = [det4, det5, det6]
            #print("D1 and D2:\n", D1, D2)
            a1 = different_sign(D2)
            a2 = different_sign(D1)
            b1, c1 = name_ind(a1)
            #print(f"a1 = {a1}, b1 = {b1}, c1 = {c1}")
            b2, c2 = name_ind(a2)
            #print(f"a2 = {a2}, b2 = {b2}, c2 = {c2}")
            Y7 = np.array([Tr2[b2] - Tr1[a1], Tr2[b2] - Tr1[b1], Tr2[b2] - Tr2[a2]])
            Y8 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[c1], Tr2[a2] - Tr2[c2]])
            #print(Y7, Y8)
            det7 = np.linalg.det(Y7)
            det8 = np.linalg.det(Y8)
            #print(det7, det8)
            if (det7 <= 0 and det8 <= 0):
                
                Y9 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[c1], Tr2[a2] - Tr2[b2]])
                Y10 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[b1], Tr2[a2] - Tr2[c2]])
                det9 = np.linalg.det(Y9)
                det10 = np.linalg.det(Y10)
                #return 1 (k-i-l-j)
                #return 2 (k-i-j-l)
                #return 3 (i-k-l-j)
                #return 4 (i-k-j-l)
                #i = a1->b1
                #j = a1->c1
                #k = a2->b2
                #l = a2->c2
                if (det9 > 0 and det10 > 0):
                    return 1
                elif (det9 > 0 and det10 <= 0):
                    return 2
                elif (det9 <= 0 and det10 > 0):
                    return 3
                elif (det9 <= 0 and det10 <= 0):
                    return 4
            else:
                flag = additional_check(Tr1, Tr2, D1, D2)
                return flag
            
def additional_check(Tr1, Tr2, D1, D2):
    #print("D1 and D2:\n", D1, D2)
    a1 = different_sign(D2)
    a2 = different_sign(D1)
    c1, b1 = name_ind(a1)
    #print(f"a1 = {a1}, b1 = {b1}, c1 = {c1}")
    c2, b2 = name_ind(a2)
    #print(f"a2 = {a2}, b2 = {b2}, c2 = {c2}")
    Y7 = np.array([Tr2[b2] - Tr1[a1], Tr2[b2] - Tr1[b1], Tr2[b2] - Tr2[a2]])
    Y8 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[c1], Tr2[a2] - Tr2[c2]])
    #print(Y7, Y8)
    det7 = np.linalg.det(Y7)
    det8 = np.linalg.det(Y8)
    #print(det7, det8)
    if (det7 <= 0 and det8 <= 0):        
        Y9 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[c1], Tr2[a2] - Tr2[b2]])
        Y10 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[b1], Tr2[a2] - Tr2[c2]])
        det9 = np.linalg.det(Y9)
        det10 = np.linalg.det(Y10)
        #return 1 (k-i-l-j)
        #return 2 (k-i-j-l)
        #return 3 (i-k-l-j)
        #return 4 (i-k-j-l)
        #i = a1->b1
        #j = a1->c1
        #k = a2->b2
        #l = a2->c2
        if (det9 > 0 and det10 > 0):
            return 1
        elif (det9 > 0 and det10 <= 0):
            return 2
        elif (det9 <= 0 and det10 > 0):
            return 3
        elif (det9 <= 0 and det10 <= 0):
            return 4
    else:
        return 0

def unique_collisions(triangles):
    collisions = []
    seen = set()
    
    for i in range(len(triangles)):
        for j in range(i+1, len(triangles)):
            Tr1 = triangles[i]
            Tr2 = triangles[j]
            
            if check_collisions(Tr1, Tr2) != 0:
                tr1 = tuple(sorted(map(tuple, Tr1)))
                tr2 = tuple(sorted(map(tuple, Tr2)))
                
                if tr1 not in seen:
                    collisions.append(Tr1)
                    seen.add(tr1)
                    
                if tr2 not in seen:
                    collisions.append(Tr2)
                    seen.add(tr2)
                    
    return collisions

def i_cal(Tr1, n3, P0, a1, b1, c1):
    P2 = Tr1[c1]
    v1 = Tr1[a1] - Tr1[c1]
    v2 = n3
    v12 = np.cross(v1, v2)
    p21 = P0 - P2
    v212 = np.cross(p21, v2)
    #print("i")
    #print("v212", v212)
    #print("v12", v12)
    k = nonzero(v12)
    t = np.linalg.norm(v212)/np.linalg.norm(v12)
    #print("t", t)
    p = P2 + t * v1
    return p

def j_cal(Tr1, n3, P0, a1, b1, c1):
    P2 = Tr1[b1]
    v1 = Tr1[a1] - Tr1[b1]
    v2 = n3
    v12 = np.cross(v1, v2)
    p21 = P0 - P2
    v212 = np.cross(p21, v2)
    #print("j")
    #print("v212", v212)
    #print("v12", v12)
    k = nonzero(v12)
    t = v212[k]/v12[k]
    #print("t", t)
    p = P2 + t * v1
    return p

def k_cal(Tr2, n3, P0, a2, b2, c2):
    P2 = Tr2[b2]
    v1 = Tr2[a2] - Tr2[b2]
    v2 = n3
    v12 = np.cross(v1, v2)
    p21 = P0 - P2
    v212 = np.cross(p21, v2)
    #print("k")
    #print("v212", v212)
    #print("v12", v12)
    k = nonzero(v12)
    t = v212[k]/v12[k]
    #print("t", t)
    p = P2 + t * v1
    return p

def l_cal(Tr2, n3, P0, a2, b2, c2):
    P2 = Tr2[c2]
    v1 = Tr2[a2] - Tr2[c2]
    v2 = n3
    v12 = np.cross(v1, v2)
    p21 = P0 - P2
    v212 = np.cross(p21, v2)
    #print("l")
    #print("v212", v212)
    #print("v12", v12)
    k = nonzero(v12)
    t = v212[k]/v12[k]
    #print("t", t)
    p = P2 + t * v1
    return p

def collision_point(Tr1, Tr2, spec):
    Y1 = np.array([Tr2[0] - Tr1[0], Tr2[0] - Tr1[1], Tr2[0] - Tr1[2]])
    Y2 = np.array([Tr2[1] - Tr1[0], Tr2[1] - Tr1[1], Tr2[1] - Tr1[2]])
    Y3 = np.array([Tr2[2] - Tr1[0], Tr2[2] - Tr1[1], Tr2[2] - Tr1[2]])
    Y4 = np.array([Tr1[0] - Tr2[0], Tr1[0] - Tr2[1], Tr1[0] - Tr2[2]])
    Y5 = np.array([Tr1[1] - Tr2[0], Tr1[1] - Tr2[1], Tr1[1] - Tr2[2]])
    Y6 = np.array([Tr1[2] - Tr2[0], Tr1[2] - Tr2[1], Tr1[2] - Tr2[2]])
    det1 = np.linalg.det(Y1)
    det2 = np.linalg.det(Y2)
    det3 = np.linalg.det(Y3)
    det4 = np.linalg.det(Y4)
    det5 = np.linalg.det(Y5)
    det6 = np.linalg.det(Y6)
    D1 = [det1, det2, det3]
    D2 = [det4, det5, det6]
    a1 = different_sign(D2)
    a2 = different_sign(D1)
    b1, c1 = name_ind(a1)
    b2, c2 = name_ind(a2)
    
    n1 = getNorm(Tr1)
    n2 = getNorm(Tr2)
    
    n3 = np.cross(n1, n2)
    left = np.array([n1, n2, n3])
    right = np.array([np.dot(n1, Tr1[0]), np.dot(n2, Tr2[0]), np.dot(n3, Tr1[0])])
    P_0 = np.linalg.solve(left, right)
    
    if(spec == 1): #i-l
        
        #p_1 = i_count(Tr1, n3, P_0, a1, b1, c1)
        p_1 = i_cal(Tr1, n3, P_0, a1, b1, c1)
        #p_2 = l_count(Tr2, n3, P_0, a2, b2, c2)
        p_2 = l_cal(Tr2, n3, P_0, a2, b2, c2)
            
    elif(spec == 2): #i-j
       
       #p_1 = i_count(Tr1, n3, P_0, a1, b1, c1)
       p_1 = i_cal(Tr1, n3, P_0, a1, b1, c1)
       #p_2 = j_count(Tr1, n3, P_0, a1, b1, c1)
       p_2 = j_cal(Tr1, n3, P_0, a1, b1, c1)
            
    elif(spec == 3): #k-l
       
       #p_1 = k_count(Tr2, n3, P_0, a2, b2, c2)
       p_1 = k_cal(Tr2, n3, P_0, a2, b2, c2)
       #p_2 = l_count(Tr2, n3, P_0, a2, b2, c2)
       p_2 = l_cal(Tr2, n3, P_0, a2, b2, c2)
       
    elif(spec == 4): #k-j
        
        #p_1 = k_count(Tr2, n3, P_0, a2, b2, c2)
        p_1 = k_cal(Tr2, n3, P_0, a2, b2, c2)
        #p_2 = j_count(Tr1, n3, P_0, a1, b1, c1)
        p_2 = j_cal(Tr1, n3, P_0, a1, b1, c1)
    """ 
    p = []
    for i in np.arange(1, 40):
        p.append(P_0+i*n3)
    return p
    """
    return p_1, p_2

n, coords = parse_stl("complicated_grain_mesh.stl")

normals = []
all_coords = []
for coord in coords:
    for c in coord: 
        all_coords.append(c)
        
unique_coords = np.unique(all_coords, axis = 0)

R = 122 #Внешний радиус геометрии заряда
for coord in coords:
    if(check_position(coord) == True):
        normals.append([0, 0, 0])
    else:
        v1 = coord[1] - coord[0]
        v2 = coord[2] - coord[1]
        normals.append(np.cross(v1, v2))


counter = 0
for n in normals:
    counter+=1
    l = np.sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
    if(l!=0):
        n[0]/=l
        n[1]/=l
        n[2]/=l



index = []
for i in range(len(unique_coords)):
    matches = np.all(all_coords == unique_coords[i], axis=1)
    matching_indexes = np.where(matches)[0]
    index.append(matching_indexes)

for i in range(len(index)):
    for j in range(len(index[i])):
        index[i][j] = index[i][j] // 3
        
normals_mean = np.zeros((len(unique_coords), 3))
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
rate = 20
for i in range(len(unique_coords)):
    matches = np.all(burn == unique_coords[i], axis = 1)
    matching_indexes = np.where(matches)[0]
    for m in matching_indexes:
        burn[m]-=(rate*normals_mean[i])
        

grouped = np.array(burn)
grouped = grouped.reshape(-1, 3, 3)
collisions = unique_collisions(grouped)
print(len(collisions))
save_to_pickle(collisions, "collisions_final.pkl")

save_to_pickle(grouped, "step20from0.pkl")









cols = np.load("collisions.npy")
vtk_points = vtk.vtkPoints()

vtk_cells = vtk.vtkCellArray()

burn = np.array(burn)
burn = burn.reshape(-1, 3, 3)

for i, polygon in enumerate(collisions):
    for point in polygon:
        vtk_points.InsertNextPoint(point)
    
    vtk_cells.InsertNextCell(3, [i*3, i*3+1, i*3+2])

poly_data = vtk.vtkPolyData()
poly_data.SetPoints(vtk_points)
poly_data.SetPolys(vtk_cells)

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(poly_data)

actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)
style = vtk.vtkInteractorStyleTrackballCamera()
render_window_interactor.SetInteractorStyle(style)

renderer.AddActor(actor)
renderer.SetBackground(0.1, 0.1, 0.1)  

render_window.Render()
render_window_interactor.Start()
