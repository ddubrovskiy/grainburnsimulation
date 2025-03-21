import vtk
import pickle
from stl import mesh
import numpy as np

def same_point(A, B):
    flag = False
    if (len(A) == len(B)):
        for i in range(len(A)):
            if (A[i] == B[i]):
                flag = True
            else: 
                flag = False
                return flag
        return flag
    else:
        return "no"

def remove_triangles(triangles, triangles_to_remove):
    """
    Удаляет все треугольники из массива 'triangles', которые присутствуют в массиве 'triangles_to_remove'.
    Каждый треугольник представлен как массив из 3 точек (каждая точка - это массив координат).
    
    :param triangles: Список всех треугольников.
    :param triangles_to_remove: Список треугольников, которые нужно удалить.
    :return: Модифицированный список треугольников.
    """
    
    # Преобразуем треугольники в кортежи для быстрого сравнения
    triangles_set = {tuple(sorted(map(tuple, triangle))) for triangle in triangles}
    triangles_to_remove_set = {tuple(sorted(map(tuple, triangle))) for triangle in triangles_to_remove}
    
    # Находим разницу между множествами
    remaining_triangles_set = triangles_set - triangles_to_remove_set
    
    # Преобразуем множество обратно в список
    remaining_triangles = [list(map(list, triangle)) for triangle in remaining_triangles_set]
    
    return remaining_triangles


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
    #print("additional")
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

n, mesh_ = parse_stl("complicated_grain_mesh.stl")

coords = []
for coord in mesh_:
    for c in coord:
        coords.append(c)
        
coords = np.array(coords)

R = 122
outer = []
for coord in coords:
    if(np.abs(np.sqrt(coord[0]**2 + coord[1]**2) - R) < 0.1):
        outer.append(coord)

##################################################### FINDING GROUPED COLLISIONS
def grouped_collisions(collisions):
    grouped_collisions = []
    for i in range(len(collisions)):
        colls = []
        colls.append(collisions[i])
        for j in range(len(collisions)):
            if i!=j:
                if check_collisions(collisions[i], collisions[j]) != 0:
                    colls.append(collisions[j])
        grouped_collisions.append(colls)
    return grouped_collisions

##################################################### PREPARING FOR TRIANGULATION OF A GROUP
def pretriangulate(groups):
    pretriangulate_points = []
    for group in groups:
        pretr = []
        for k in group[0]:
            pretr.append(k)
        for i in range(1, len(group)):
            p1, p2 = collision_point(group[0], group[i], check_collisions(group[0], group[i]))
            pretr.append(p1)
            pretr.append(p2)
        pretriangulate_points.append(pretr)
    
    return pretriangulate_points

trs = load_from_pickle("step20from0.pkl")


collisions = load_from_pickle("collisions_final.pkl")
"""
groups = grouped_collisions(collisions)
save_to_pickle(groups, "grouped_collisions.pkl")
trngltpnts = pretriangulate(groups)
save_to_pickle(trngltpnts, "grouped_triangulation_points.pkl")
"""
groups = load_from_pickle("grouped_collisions.pkl")
trngltpnts = load_from_pickle("grouped_triangulation_points.pkl")

"""
colpnts = []
for i in range(len(collisions)):
    for j in range(i, len(collisions)):
        if(i != j):
            if (check_collisions(collisions[i], collisions[j]) != 0):
                p1, p2 = collision_point(collisions[i], collisions[j], check_collisions(collisions[i], collisions[j]))
                colpnts.append(p1)
                colpnts.append(p2)
save_to_pickle(colpnts, "all_collision_points.pkl")
"""
colpnts = load_from_pickle("all_collision_points.pkl")
vtk_points = vtk.vtkPoints()
more_points = vtk.vtkPoints()

to_show = 1

bnd_points = vtk.vtkPoints()
for i in colpnts:
    bnd_points.InsertNextPoint(i[0], i[1], i[2])

vtk_cells = vtk.vtkCellArray()
for i, polygon in enumerate(trs):
    for point in polygon:
        vtk_points.InsertNextPoint(point)
    
    vtk_cells.InsertNextCell(3, [i*3, i*3+1, i*3+2])

poly_data = vtk.vtkPolyData()
poly_point = vtk.vtkPolyData()
poly_point.SetPoints(bnd_points)
poly_data.SetPoints(vtk_points)
poly_data.SetPolys(vtk_cells)

glyph = vtk.vtkVertexGlyphFilter()
glyph.SetInputData(poly_point)

pointMapper = vtk.vtkPolyDataMapper()
pointMapper.SetInputConnection(glyph.GetOutputPort())

pointActor = vtk.vtkActor()
colors = vtk.vtkNamedColors()
pointActor.SetMapper(pointMapper)
pointActor.GetProperty().SetPointSize(5)
pointActor.GetProperty().SetColor(colors.GetColor3d("DeepPink"))
pointActor.GetProperty().RenderPointsAsSpheresOn()

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(poly_data)


actor = vtk.vtkActor()
actor.SetMapper(mapper)


actor.GetProperty().SetSpecular(0.0)
actor.GetProperty().SetAmbient(0.3)
actor.GetProperty().SetDiffuse(1)
actor.GetProperty().SetColor(0.5, 0.5, 1)

renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)
style = vtk.vtkInteractorStyleTrackballCamera()
render_window_interactor.SetInteractorStyle(style)

renderer.AddActor(actor)
renderer.AddActor(pointActor)
renderer.SetBackground(0.1, 0.1, 0.1)  


render_window.Render()
render_window_interactor.Start()