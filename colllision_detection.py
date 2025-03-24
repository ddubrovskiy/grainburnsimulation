import numpy as np
import vtk
import pickle
from stl import mesh

def save_to_pickle(array, filename):
    with open(filename, "wb") as f:
        pickle.dump(array, f)
        
def load_from_pickle(filename):
    with open(filename, "rb") as f:
        loaded = pickle.load(f)
    return loaded

def parse_stl(filename):
    stl_mesh = mesh.Mesh.from_file(filename)
    num_faces = stl_mesh.vectors.shape[0]
    vertices = stl_mesh.vectors
    return num_faces, vertices

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

def max_abs_index(array):
    max_value = array[0]
    max_index = 0
    
    for i in range(len(array)):
        if np.abs(array[i]) > np.abs(max_value):
            max_value = array[i]
            max_index = i
    
    return max_index

#Нерабочий образец, в случае чего - повторить для j, k, l
def i_count(Tr1, n3, P_0, a1, b1, c1):
    denarr = np.array([n3[0]*(Tr1[b1][1] - Tr1[a1][1]) - n3[1]*(Tr1[b1][0] - Tr1[a1][0]),      #xy
                           n3[0]*(Tr1[b1][2] - Tr1[a1][2]) - n3[2]*(Tr1[b1][0] - Tr1[a1][0]),  #xz
                           n3[2]*(Tr1[b1][1] - Tr1[a1][1]) - n3[1]*(Tr1[b1][2] - Tr1[a1][2])]) #zy
        
    ind = max_abs_index(denarr)
    denom = denarr[ind]
    if(ind == 0):
        t = n3[0]*(P_0[1] - Tr1[a1][1]) + n3[1]*(Tr1[a1][0] - P_0[0])/denom
    elif(ind == 1):
        t = n3[0]*(P_0[2] - Tr1[a1][2]) + n3[2]*(Tr1[a1][0] - P_0[0])/denom
    elif(ind == 2):
        t = n3[2]*(P_0[1] - Tr1[a1][1]) + n3[1]*(Tr1[a1][2] - P_0[2])/denom
        
    point = Tr1[a1] + t * (Tr1[b1] - Tr1[a1])
    
    return point

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

def j_count(Tr1, n3, P_0, a1, b1, c1):
    denarr = np.array([n3[0]*(Tr1[c1][1] - Tr1[a1][1]) - n3[1]*(Tr1[c1][0] - Tr1[a1][0]),      #xy
                           n3[0]*(Tr1[c1][2] - Tr1[a1][2]) - n3[2]*(Tr1[c1][0] - Tr1[a1][0]),  #xz
                           n3[2]*(Tr1[c1][1] - Tr1[a1][1]) - n3[1]*(Tr1[c1][2] - Tr1[a1][2])]) #zy
        
    ind = max_abs_index(denarr)
    denom = denarr[ind]
    if(ind == 0):
        t = n3[0]*(P_0[1] - Tr1[a1][1]) + n3[1]*(Tr1[a1][0] - P_0[0])/denom
    elif(ind == 1):
        t = n3[0]*(P_0[2] - Tr1[a1][2]) + n3[2]*(Tr1[a1][0] - P_0[0])/denom
    elif(ind == 2):
        t = n3[2]*(P_0[1] - Tr1[a1][1]) + n3[1]*(Tr1[a1][2] - P_0[2])/denom
        
    point = Tr1[a1] + t * (Tr1[c1] - Tr1[a1])
    
    return point

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

def k_count(Tr2, n3, P_0, a2, b2, c2):
    
    dxy = n3[0]*(Tr2[b2][1] - Tr2[a2][1]) - n3[1]*(Tr2[b2][0] - Tr2[a2][0])
    dxz = n3[0]*(Tr2[b2][2] - Tr2[a2][2]) - n3[2]*(Tr2[b2][0] - Tr2[a2][0])
    dzy = n3[2]*(Tr2[b2][1] - Tr2[a2][1]) - n3[1]*(Tr2[b2][2] - Tr2[a2][2])
    
    if(abs(dxy) > abs(dxz) and abs(dxy) > abs(dzy)):
        t = n3[0]*(P_0[1] - Tr2[a2][1]) + n3[1]*(Tr2[a2][0] - P_0[0])/dxy
    elif(abs(dxz) > abs(dxy) and abs(dxz) > abs(dzy)):
        t = n3[0]*(P_0[2] - Tr2[a2][2]) + n3[2]*(Tr2[a2][0] - P_0[0])/dxz
    else:
        t = n3[2]*(P_0[1] - Tr2[a2][1]) + n3[1]*(Tr2[a2][2] - P_0[2])/dzy
   
    point = Tr2[a2] + t * (Tr2[b2] - Tr2[a2])
    print(t)
    
    return point

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
    
def l_count(Tr2, n3, P_0, a2, b2, c2):
    denarr = np.array([n3[0]*(Tr2[c2][1] - Tr2[a2][1]) - n3[1]*(Tr2[c2][0] - Tr2[a2][0]),      #xy
                           n3[0]*(Tr2[c2][2] - Tr2[a2][2]) - n3[2]*(Tr2[c2][0] - Tr2[a2][0]),  #xz
                           n3[2]*(Tr2[c2][1] - Tr2[a2][1]) - n3[1]*(Tr2[c2][2] - Tr2[a2][2])]) #zy
        
    ind = max_abs_index(denarr)
    denom = denarr[ind]
    if(ind == 0):
        t = n3[0]*(P_0[1] - Tr2[a2][1]) + n3[1]*(Tr2[a2][0] - P_0[0])/denom
    elif(ind == 1):
        t = n3[0]*(P_0[2] - Tr2[a2][2]) + n3[2]*(Tr2[a2][0] - P_0[0])/denom
    elif(ind == 2):
        t = n3[2]*(P_0[1] - Tr2[a2][1]) + n3[1]*(Tr2[a2][2] - P_0[2])/denom
    
    point = Tr2[c2] + t * (Tr2[a2] - Tr2[c2])
    print(t)
    
    return point

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
            a1 = different_sign(D2)
            a2 = different_sign(D1)
            b1, c1 = name_ind(a1)
            b2, c2 = name_ind(a2)
            Y7 = np.array([Tr2[b2] - Tr1[a1], Tr2[b2] - Tr1[b1], Tr2[b2] - Tr2[a2]])
            Y8 = np.array([Tr2[a2] - Tr1[a1], Tr2[a2] - Tr1[c1], Tr2[a2] - Tr2[c2]])
            det7 = np.linalg.det(Y7)
            det8 = np.linalg.det(Y8)
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

num, grain = parse_stl("complicated_grain_mesh.stl") 
col = []
coordinates = np.load("step20from0.npy")
coor = np.copy(coordinates)
coordinates = coordinates.reshape(-1, 3, 3)     

"""
bnd = []
for i in range(len(clsn_)):
    for j in range(i+1, len(clsn_)):
        var = check_collisions(clsn_[i], clsn_[j])
        if(var!=0):
            p_1, p_2 = collision_point(clsn_[i], clsn_[j], var)
            bnd.append(p_1)
            bnd.append(p_2)
print(bnd)
np.save("bounds.npy", bnd)
"""
lbnd = np.load("bounds.npy")
vtk_points = vtk.vtkPoints()

bnd_points = vtk.vtkPoints()
coord_points = vtk.vtkPoints()
for i in lbnd:
    bnd_points.InsertNextPoint(i[0], i[1], i[2])
for i in coor:
    coord_points.InsertNextPoint(i[0], i[1], i[2])

pklcrd = load_from_pickle("collisions_final.pkl")
trs = load_from_pickle("step20from0.pkl")

vtk_cells = vtk.vtkCellArray()
for i, polygon in enumerate(pklcrd):
    for point in polygon:
        vtk_points.InsertNextPoint(point)
    vtk_cells.InsertNextCell(3, [i*3, i*3+1, i*3+2])

poly_data = vtk.vtkPolyData()
poly_point = vtk.vtkPolyData()
poly_coord = vtk.vtkPolyData()
poly_point.SetPoints(bnd_points)
poly_coord.SetPoints(coord_points)
poly_data.SetPoints(vtk_points)
poly_data.SetPolys(vtk_cells)

# Add random colors to cells
cell_colors = vtk.vtkUnsignedCharArray()
cell_colors.SetNumberOfComponents(3)
cell_colors.SetNumberOfTuples(poly_data.GetNumberOfCells())

import random
for i in range(poly_data.GetNumberOfCells()):
    r = random.randint(0, 255)
    g = random.randint(0, 255)
    b = random.randint(0, 255)
    cell_colors.InsertTuple3(i, r, g, b)

poly_data.GetCellData().SetScalars(cell_colors)

glyph = vtk.vtkVertexGlyphFilter()
glyph.SetInputData(poly_point)
coordglyph = vtk.vtkVertexGlyphFilter()
coordglyph.SetInputData(poly_coord)

pointMapper = vtk.vtkPolyDataMapper()
pointMapper.SetInputConnection(glyph.GetOutputPort())
pointActor = vtk.vtkActor()
colors = vtk.vtkNamedColors()
pointActor.SetMapper(pointMapper)
pointActor.GetProperty().SetPointSize(5)
pointActor.GetProperty().SetColor(colors.GetColor3d("DeepPink"))
pointActor.GetProperty().RenderPointsAsSpheresOn()

coordMapper = vtk.vtkPolyDataMapper()
coordMapper.SetInputConnection(coordglyph.GetOutputPort())
coordActor = vtk.vtkActor()
coordActor.SetMapper(coordMapper)
coordActor.GetProperty().SetPointSize(5)
coordActor.GetProperty().SetColor(colors.GetColor3d("Blue"))
coordActor.GetProperty().RenderPointsAsSpheresOn()

mapper = vtk.vtkPolyDataMapper()
mapper.SetInputData(poly_data)
mapper.SetScalarModeToUseCellData()
mapper.SetColorModeToDirectScalars()

actor = vtk.vtkActor()
actor.SetMapper(mapper)
actor.GetProperty().SetRepresentationToSurface()
actor.GetProperty().SetOpacity(1.0)

def get_triangle_center(triangle_points):
    p1 = triangle_points[0]
    p2 = triangle_points[1]
    p3 = triangle_points[2]
    center = [(p1[0] + p2[0] + p3[0]) / 3, 
              (p1[1] + p2[1] + p3[1]) / 3, 
              (p1[2] + p2[2] + p3[2]) / 3]
    return center

# Добавление текста с индексом треугольника в его центр
def add_text_to_renderer(renderer, center, index):
    if index == 2:
        text = vtk.vtkVectorText()
        text.SetText(str(index))  # Индекс треугольника

        text_mapper = vtk.vtkPolyDataMapper()
        text_mapper.SetInputConnection(text.GetOutputPort())

        text_actor = vtk.vtkActor()
        text_actor.SetMapper(text_mapper)
        text_actor.SetPosition(center[0], center[1], center[2])  # Центр треугольника

        #text_actor.GetPositionCoordinate().SetCoordinateSystemToWorld()  # Мировые координаты
        renderer.AddActor(text_actor)

renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

# Добавление индексов для каждого треугольника
for i in range(len(pklcrd)):  # Для каждого треугольника
    triangle_points = [pklcrd[i][0], pklcrd[i][1], pklcrd[i][2]]
    center = get_triangle_center(triangle_points)  # Находим центр треугольника
    add_text_to_renderer(renderer, center, i)  # Добавляем текст с индексом





render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)
style = vtk.vtkInteractorStyleTrackballCamera()
render_window_interactor.SetInteractorStyle(style)

def update_property(obj, event, property_name):
    slider_widget = obj.GetRepresentation()
    value = slider_widget.GetValue()
    if property_name == "ambient":
        actor.GetProperty().SetAmbient(value)
    elif property_name == "diffuse":
        actor.GetProperty().SetDiffuse(value)
    elif property_name == "specular":
        actor.GetProperty().SetSpecular(value)
    render_window.Render()

def create_slider(title, min_val, max_val, init_val, y_position, property_name):
    slider_rep = vtk.vtkSliderRepresentation2D()
    slider_rep.SetMinimumValue(min_val)
    slider_rep.SetMaximumValue(max_val)
    slider_rep.SetValue(init_val)
    slider_rep.SetTitleText(title)
    slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint1Coordinate().SetValue(0.1, y_position)
    slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
    slider_rep.GetPoint2Coordinate().SetValue(0.4, y_position)
    slider_rep.SetSliderLength(0.02)
    slider_rep.SetSliderWidth(0.03)
    slider_rep.SetEndCapLength(0.02)
    slider_rep.SetEndCapWidth(0.03)
    slider_rep.SetTubeWidth(0.005)
    slider_rep.SetLabelFormat("%.2f")
    slider_rep.SetTitleHeight(0.02)
    slider_rep.SetLabelHeight(0.02)
    slider_widget = vtk.vtkSliderWidget()
    slider_widget.SetInteractor(render_window_interactor)
    slider_widget.SetRepresentation(slider_rep)
    slider_widget.SetAnimationModeToAnimate()
    slider_widget.EnabledOn()
    slider_widget.AddObserver("InteractionEvent", lambda obj, event: update_property(obj, event, property_name))
    return slider_widget

# Create sliders for ambient, diffuse, and specular properties
ambient_slider = create_slider("Ambient", 0.0, 1.0, 0.3, 0.8, "ambient")
diffuse_slider = create_slider("Diffuse", 0.0, 1.0, 1.0, 0.6, "diffuse")
specular_slider = create_slider("Specular", 0.0, 1.0, 0.0, 0.4, "specular")

renderer.AddActor(actor)
#renderer.AddActor(pointActor)
#renderer.AddActor(coordActor)
renderer.SetBackground(0.1, 0.1, 0.1)
render_window.Render()
render_window_interactor.Start()