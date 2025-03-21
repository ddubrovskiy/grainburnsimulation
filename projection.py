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

def to_local(points):
    local_points = []
    base = np.array([0, 0])
    local_points.append(base)
    h = np.sqrt((points[1][0] - points[0][0])**2 + (points[1][1] - points[0][1])**2 + (points[1][2] - points[0][2])**2)
    #print(f"h = {h}")
    point_h = np.array([h, 0])
    local_points.append(point_h)
    for k in range(2, len(points)):
        i = ((points[k][0] - points[0][0])*(points[1][0] - points[0][0]) + (points[k][1] - points[0][1])*(points[1][1] - points[0][1]) + (points[k][2] - points[0][2])*(points[1][2] - points[0][2]))/(h)
        #print(f"i = {i}")
        j = np.sign((points[k][0] - points[0][0])**2 + (points[k][1] - points[0][1])**2 + (points[k][2] - points[0][2])**2 - i**2) * np.sqrt(np.abs((points[k][0] - points[0][0])**2 + (points[k][1] - points[0][1])**2 + (points[k][2] - points[0][2])**2 - i**2))
        #print(f"In square root = {(points[k][0] - points[0][0])**2 + (points[k][1] - points[0][1])**2 + (points[k][2] - points[0][2])**2 - i**2}")
        #print(f"j = {j}")
        point_ij = np.array([i, j])
        local_points.append(point_ij)
    #print("Точки в трёхмерных координатах:")
    #print(np.round(points, 3))
    #print("Точки в двумерных координатах:")
    #print(np.round(local_points, 3))
    return local_points


def triangulation_2d(points):
    tri = Delaunay(points)
    simpls = tri.simplices
    #print(simpls)
    return simpls

def triangulation_3d(points, simpls):
    newtrs = []
    for i in range(len(simpls)):
        newtr = []
        for j in simpls[i]:
            newtr.append(points[j])
        newtrs.append(newtr)
    newtrs = np.array(newtrs)
    return newtrs


def triangulate_with_3d_points(points):
    lcl_pnts = []
    for i in range(len(points)):
        lcl = to_local(points[i])
        lcl = np.array(lcl)
        lcl_pnts.append(lcl)
    new_mesh = []
    for i in range(len(lcl_pnts)):
        simpls = triangulation_2d(lcl_pnts[i])
        new_triangles = triangulation_3d(points[i], simpls)
        new_triangles = np.array(new_triangles)
        new_mesh.append(new_triangles)
        
    return new_mesh
    

pnts = load_from_pickle("grouped_triangulation_points.pkl")
print(len(pnts))
remesh = triangulate_with_3d_points(pnts)
save_to_pickle(remesh, "remeshed_groups.pkl")
lens = []
for i in range(len(remesh)):
    lens.append(len(remesh[i]))
print(lens[1539])
lens = np.array(lens)
index = np.argmax(lens)
print(index)
new_triangles = remesh[1539]
"""
pnt = pnts[0]
print(pnt)
lcl_pnts = []
errors = []
for i in range(len(pnts)):
    lcl = to_local(pnts[i])
    if (np.any(np.isnan(lcl))): errors.append(i)
    lcl = np.array(lcl)
    #input()
    lcl_pnts.append(lcl)
print(len(errors))
lcl_pnt = lcl_pnts[0]
print(np.round(lcl_pnt, 1))
"""








vtk_points = vtk.vtkPoints()

vtk_cells = vtk.vtkCellArray()

for i, polygon in enumerate(new_triangles):
    for point in polygon:
        vtk_points.InsertNextPoint(point)
    
    vtk_cells.InsertNextCell(3, [i*3, i*3+1, i*3+2])

poly_data = vtk.vtkPolyData()
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
