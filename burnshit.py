import numpy as np
import vtk

def triangle_plane(triangle):
    """Compute the plane equation (normal vector and d) of a triangle."""
    v0, v1, v2 = triangle
    normal = np.cross(v1 - v0, v2 - v0)  # Normal vector
    normal = normal / np.linalg.norm(normal)  # Normalize
    d = -np.dot(normal, v0)
    return normal, d

def intersect_planes(plane1, plane2):
    """Find the line of intersection of two planes."""
    n1, d1 = plane1
    n2, d2 = plane2

    # Direction vector of the line (cross product of normals)
    line_dir = np.cross(n1, n2)
    if np.allclose(line_dir, 0):  # Parallel or coincident planes
        return None

    # Find a point on the line
    A = np.array([n1, n2, line_dir])
    b = np.array([-d1, -d2, 0])
    point_on_line = np.linalg.solve(A.T, b)  # Solve for a point on the line

    return point_on_line, line_dir

def clip_line_to_triangle(triangle, line_point, line_dir):
    """Clip the line segment to the triangle."""
    v0, v1, v2 = triangle
    edges = [(v0, v1), (v1, v2), (v2, v0)]

    # Intersect line with edges
    segment_points = []
    for edge_start, edge_end in edges:
        edge_dir = edge_end - edge_start
        A = np.vstack([-line_dir, edge_dir]).T
        b = edge_start - line_point
        try:
            t, u = np.linalg.solve(A, b)
            if 0 <= u <= 1:  # Edge segment bounds
                segment_points.append(line_point + t * line_dir)
        except np.linalg.LinAlgError:
            continue

    # Return unique points (sorted)
    if len(segment_points) >= 2:
        segment_points = np.unique(segment_points, axis=0)
        if len(segment_points) == 2:
            return segment_points
    return None

def intersect_triangles(triangle1, triangle2):
    """Find the line segment of intersection between two triangles."""
    plane1 = triangle_plane(triangle1)
    plane2 = triangle_plane(triangle2)

    # Intersect planes
    line = intersect_planes(plane1, plane2)
    if line is None:
        return None
    line_point, line_dir = line

    # Clip line to both triangles
    segment1 = clip_line_to_triangle(triangle1, line_point, line_dir)
    segment2 = clip_line_to_triangle(triangle2, line_point, line_dir)

    if segment1 is not None and segment2 is not None:
        # Find the overlapping segment
        intersection = np.array([np.maximum(segment1[0], segment2[0]),
                                  np.minimum(segment1[1], segment2[1])])
        if np.all(intersection[0] <= intersection[1]):
            return intersection
    return None

# Example Usage
triangle1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
triangle2 = np.array([[0, 0, 0], [1, 1, -1], [0, 1, 1]])
clsn = np.array([triangle1, triangle2])

intersection = intersect_triangles(triangle1, triangle2)
print("Intersection Line Segment:", intersection)

vtk_points = vtk.vtkPoints()

vtk_cells = vtk.vtkCellArray()

for i, polygon in enumerate(clsn):
    for point in polygon:
        vtk_points.InsertNextPoint(point)
    
    vtk_cells.InsertNextCell(3, [i*3, i*3+1, i*3+2])

poly_data = vtk.vtkPolyData()
poly_point = vtk.vtkPolyData()
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

renderer = vtk.vtkRenderer()
render_window = vtk.vtkRenderWindow()
render_window.AddRenderer(renderer)

render_window_interactor = vtk.vtkRenderWindowInteractor()
render_window_interactor.SetRenderWindow(render_window)
"""
axes = vtk.vtkAxesActor()
widget = vtk.vtkOrientationMarkerWidget()
rgba = [0] * 4
colors = vtk.vtkNamedColors()
colors.GetColor('Carrot', rgba)
widget.SetOutlineColor(rgba[0], rgba[1], rgba[2])
widget.SetOrientationMarker(axes)
widget.SetInteractor(render_window_interactor)
widget.SetViewport(0.0, 0.0, 0.4, 0.4)
widget.SetEnabled(1)
widget.InteractiveOn()
"""
style = vtk.vtkInteractorStyleTrackballCamera()
render_window_interactor.SetInteractorStyle(style)


renderer.AddActor(actor)
renderer.AddActor(pointActor)
renderer.SetBackground(0.1, 0.1, 0.1)  


render_window.Render()
render_window_interactor.Start()