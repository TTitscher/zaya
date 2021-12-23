
import vtk
from vtk.util import numpy_support
import numpy as np
import time


class SphereVisualizer:
    def __init__(self, N):
        self.N = N

        self.positions = vtk.vtkPoints()
        self.diameters = vtk.vtkDoubleArray()

        for i in range(N):
            self.positions.InsertNextPoint(0,0,0)
            self.diameters.InsertNextValue(0)
        self.diameters.SetName("diameter")

        self.window = self._sphere_render_window()

    def add_box(self, lx, ly, lz):
        cube = vtk.vtkCubeSource()  
        cube.SetXLength(lx)
        cube.SetYLength(ly)
        cube.SetZLength(lz)
        cube_mapper = vtk.vtkPolyDataMapper()
        cube_mapper.SetInputConnection(cube.GetOutputPort())
        cube_actor = vtk.vtkActor()
        cube_actor.SetMapper(cube_mapper)
        cube_actor.SetPosition(lx/2, ly/2, lz/2)
        cube_actor.GetProperty().SetRepresentationToWireframe()
        self.renderer.AddActor(cube_actor)


    def _reference_sphere(self, resolution = 20):
        sphere = vtk.vtkSphereSource()
        sphere.SetThetaResolution(resolution)
        sphere.SetPhiResolution(int(resolution/2))
        sphere.SetRadius(0.5)
        return sphere

    def _sphere_render_window(self):
        grid = vtk.vtkUnstructuredGrid()
        grid.SetPoints(self.positions)
        grid.GetPointData().AddArray(self.diameters)
        grid.GetPointData().SetActiveScalars("diameter")
        # poly_data.GetPointData().AddArray(self.radii)

        sphere_source = self._reference_sphere()
        glyph= vtk.vtkGlyph3D()
        glyph.SetInputData(grid)
        glyph.SetSourceConnection(sphere_source.GetOutputPort())

        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(glyph.GetOutputPort())
        mapper.SetScalarModeToUsePointFieldData()

        actor = vtk.vtkActor()
        actor.SetMapper(mapper)

        # create a rendering window and renderer
        renderer = vtk.vtkRenderer()
        renderer.SetBackground(1,1,1)
        renderer.AddActor(actor)
        
        self.txt = vtk.vtkTextActor()
        self.txt.GetTextProperty().SetFontFamilyToArial()
        self.txt.GetTextProperty().SetFontSize(18)
        self.txt.GetTextProperty().SetColor(0,0,0)
        self.txt.SetDisplayPosition(20,30)
        renderer.AddActor(self.txt)

        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)
        render_window.SetSize(500,500)

        self.renderer = renderer
        return render_window
         
    def update_data(self, positions, di):
        vtk_p = numpy_support.numpy_to_vtk(positions)
        vtk_d = numpy_support.numpy_to_vtk(di)
        self.positions.SetData(vtk_p)
        self.diameters.SetArray(vtk_d, len(di), 1)

        self.positions.Modified()
        self.diameters.Modified()

    def update_txt(self, txt):
        self.txt.SetInput(txt)


    def show(self):
        interactor = vtk.vtkRenderWindowInteractor()
        interactor.SetRenderWindow(self.window)
        interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        interactor.Initialize()
        interactor.Start()
        interactor.GetRenderWindow().Finalize()


class Animation:
    def __init__(self, window):
        self.window = window
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.window)

        # enable user interface interactor
        self.interactor.SetInteractorStyle(vtk.vtkInteractorStyleTrackballCamera())
        self.interactor.Initialize()

        self.interactor.AddObserver('TimerEvent', self.update_scene)
        self.interactor.CreateRepeatingTimer(20)

        self.t = 0;
        self.dt = 0.02;
        self.timings = {'Update': 0., 'Render' :0.}
    
    def update_scene(self, *args):
        t_start = time.time()
        self.t += self.dt
        self.update(self.t)
        t_update = time.time()

        self.window.Render()
        t_render = time.time()

        self.timings['Update'] +=  t_update - t_start
        self.timings['Render'] +=  t_render - t_update

    def update(self,t):
        raise RuntimeError("implement me!")

    def start(self):
        self.interactor.Start()
        self.interactor.GetRenderWindow().Finalize()

class SphereAnimation(Animation):
    def __init__(self, vvv, positions, radii):
        self.vvv = vvv
        self.positions = positions
        self.radii = radii
        super().__init__(self.vvv.window)
        
    def update(self, t):
        self.positions[:,1] += 0.01 *  np.sin(t)
        self.radii[:] += 0.01 *  np.cos(t)
        self.vvv.update_data(self.positions, self.radii)

if __name__ == "__main__":
    v = SphereVisualizer(2)
    v.add_box(1,1,1)
    v.update_data([[0.2, 0.5, 0.5],[0.7, 0.5, 0.5]], [0.4, 0.6])
    v.show()
    exit()

    N = 100
    vvv = SphereVisualizer(100)

    data = np.zeros((N, 4))
    for i, row in enumerate(data):
        row[0] = i
        row[3] = i*0.01
    p = data[:, 0:3]
    r = data[:, 3]

    animatin = SphereAnimation(vvv, p, r);
    animatin.start()
    print(animatin.timings)
    



