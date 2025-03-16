import numpy as np
import vtk

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

def check_collisions(Tr1, Tr2):
    A_min = np.min(Tr1, axis=0)
    A_max = np.max(Tr1, axis=0)
    B_min = np.min(Tr2, axis=0)
    B_max = np.max(Tr2, axis=0)
    if(B_min[0] > A_max[0] or B_min[1] > A_max[1] or B_min[2] > A_max[2] or A_min[0] > B_max[0] or A_min[1] > B_max[1] or A_min[2] > B_max[2]):
        return False
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
            return False
        elif(det1 < 0 and det2 < 0 and det3 < 0):
            return False
        elif(det4 > 0 and det5 > 0 and det6 > 0):
            return False
        elif(det4 < 0 and det5 < 0 and det6 < 0):
            return False
        elif(det1 == 0 or det2 == 0 or det3 == 0):
            return False
        elif(det4 == 0 or det5 == 0 or det6 == 0):
            return False
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
                return True
            else:
                return False

coordinates = np.load("step20from0.npy")

grouped_coordinates = coordinates.reshape(-1, 3, 3)
print(grouped_coordinates[0])

colTr = []
for i in range(len(grouped_coordinates)):
    ch = 0
    for j in range(i+1, len(grouped_coordinates)):
        if(check_collisions == True):
            ch+=1
            if(ch==1):
                colTr.append(grouped_coordinates[i])
                colTr.append(grouped_coordinates[j])
            else:
                colTr.append(grouped_coordinates[j])
print(colTr)
                
            

vtk_points = vtk.vtkPoints()

vtk_cells = vtk.vtkCellArray()

for i, polygon in enumerate(grouped_coordinates):
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
