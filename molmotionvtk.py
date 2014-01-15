import sys
import numpy as np
import scipy.spatial.distance as dist
import vtk
from vtk.util import numpy_support as vnp


# np.random.seed(1)

###########################################
# This part sets up the data
###########################################

# We emulate movement of molecules starting with two kinds A and B
# which react to create C.  A + B -> C

particleA = {'icount': 100, # initial count
             'color': np.array((1,0,0), order='C'),   # color for display
             'name': 'A'}
particleB  = {'icount': 50,
              'color': (0,0,1),
              'name': 'B'}
particleC  = {'icount': 0,
              'color': (0,1,0),
              'name': 'C'}

particleA['pos'] = np.zeros((particleA['icount'], 3))
particleB['pos'] = np.zeros((particleB['icount'], 3), order='C')
# Since 1 molecule each of A and B create 1 molecule of C, the maximum
# no. of C is minimum of A and B
particleC['pos'] = np.zeros((min(particleB['icount'], particleA['icount']), 3), order='C')

# Cylindrical compartment inside which all this happens
compartment = {'length': 10e-6,
               'radius': 1e-6}

# Initialize position of A and B. C is not present initially.
# Assign position to A molecules - x position is one sided Gaussian
particleA['pos'][:,0] = np.abs(np.random.normal(0.0, compartment['length']/5, particleA['pos'].shape[0]))
# particleA['pos'][:,0] = np.linspace(0, compartment['length'], particleA['icount'])
# y is Gaussian
particleA['pos'][:,1] = np.random.normal(0.0, compartment['radius']/10, particleA['pos'].shape[0])
# z is uniform positive
particleA['pos'][:,2] = np.random.rand(particleA['pos'].shape[0])*compartment['radius']

# B starts from  the other end of the cylinder
particleB['pos'][:,0] = np.abs(compartment['length'] - 
                               np.abs(np.random.normal(0.0,
                                                       compartment['length']/5,
                                                       particleB['pos'].shape[0])))

# particleB['pos'][:,0] = np.linspace(0, compartment['length'], particleB['pos'].shape[0])
particleB['pos'][:,1] = np.random.normal(0.0,
                                         compartment['radius']/10,
                                         particleB['pos'].shape[0])
particleB['pos'][:,2] = np.random.normal(0.0,
                                         compartment['radius']/10,
                                         particleB['pos'].shape[0])


# Scalars for size of sphere glyphs representing molecules. Setting
# this to 0 will cause the molecule to vanish from view
particleA['ss'] = np.ones(particleA['icount'], order='C')
particleB['ss'] = np.ones(particleB['icount'], order='C')
particleC['ss'] = np.zeros(min(particleB['icount'], particleA['icount']), order='C')
# print particleC['ss']

# If A and B come within bindingRadius, they turn into C
bindingRadius = 0.1e-6

particleScale = 1e-7
particleSpeed = 0.05e-6
particleResolution=8

steps = 1000

###########################################
# Actual visualization
###########################################

                    
def initPoints(particleDict):
    ptsA = vtk.vtkPoints()
    posA = vnp.numpy_to_vtk(particleDict['pos'], deep=True)
    ptsA.SetData(posA)
    polyDataA = vtk.vtkPolyData()
    polyDataA.SetPoints(ptsA)
    particleDict['polydata'] = polyDataA
    sizesA = vnp.numpy_to_vtk(particleDict['ss'], deep=True)
    polyDataA.GetPointData().SetScalars(sizesA)
    sourceA = vtk.vtkSphereSource()
    sourceA.SetRadius(particleScale)
    sourceA.SetThetaResolution(20)
    sourceA.SetPhiResolution(20)
    glyphA = vtk.vtkGlyph3D()
    glyphA.SetSource(sourceA.GetOutput())
    glyphA.SetInput(polyDataA)
    glyphA.ScalingOn()
    glyphA.SetScaleModeToScaleByScalar()
    # glyphA.SetScaleFactor(particleScale)
    mapperA = vtk.vtkPolyDataMapper()
    mapperA.SetInput(glyphA.GetOutput())
    mapperA.ImmediateModeRenderingOn()
    colorFnA = vtk.vtkColorTransferFunction()
    colorFnA.AddRGBPoint(0, *particleDict['color'])
    colorFnA.AddRGBPoint(1, *particleDict['color'])
    mapperA.SetLookupTable(colorFnA)
    actorA = vtk.vtkActor()
    actorA.SetMapper(mapperA)
    actorA.GetProperty().SetOpacity(0.5)
    actorA.GetProperty().SetColor(particleDict['color'])
    particleDict['points'] = ptsA
    particleDict['polydata'] = polyDataA
    particleDict['source'] = sourceA
    particleDict['glyph'] = glyphA
    particleDict['colorfn'] = colorFnA
    particleDict['mapper'] = mapperA
    particleDict['actor'] = actorA
    return particleDict

    
###########################################
# Display
###########################################

def visualize():
    A = initPoints(particleA)
    B = initPoints(particleB)
    C = initPoints(particleC)
    displayDict = {}
    # Set up the renderer and redering window
    renderer = vtk.vtkRenderer()
    renwin = vtk.vtkRenderWindow()
    renwin.SetSize(1600, 1400)
    renwin.StereoCapableWindowOn()
    renwin.StereoRenderOn()
    renwin.SetStereoTypeToCrystalEyes()
    renwin.AddRenderer(renderer)
    ptsComp = vtk.vtkPoints()
    posComp = vnp.numpy_to_vtk(np.array([[0, 0, 0]], order='C'), deep=True)
    ptsComp.SetData(posComp)
    polyDataComp = vtk.vtkPolyData()
    polyDataComp.SetPoints(ptsComp)
    sourceComp = vtk.vtkCylinderSource()
    sourceComp.SetResolution(20)
    sourceComp.SetRadius(compartment['radius'] + particleScale)
    sourceComp.SetHeight(compartment['length'] + 2 * particleScale)
    sourceComp.SetCenter(0, -compartment['length']/2, 0)
    glyphComp = vtk.vtkGlyph3D()
    glyphComp.SetSource(sourceComp.GetOutput())
    glyphComp.SetInput(polyDataComp)
    glyphComp.ScalingOff()
    # glyphComp.SetScaleModeToDefault()
    # glyphComp.SetScaleFactor(particleScale)
    mapperComp = vtk.vtkPolyDataMapper()
    mapperComp.SetInput(glyphComp.GetOutput())
    mapperComp.ImmediateModeRenderingOn()
    actorComp = vtk.vtkActor()
    actorComp.SetMapper(mapperComp)
    actorComp.GetProperty().SetOpacity(0.2)
    actorComp.GetProperty().SetColor(1, 1, 1)
    actorComp.SetOrientation(0, 90, 90)
    renderer.AddActor(actorComp)
    renderer.AddActor(A['actor'])
    renderer.AddActor(B['actor'])
    renderer.AddActor(C['actor'])
    print 'Create camera'
    camera = vtk.vtkCamera()
    camera.SetPosition(300e-6, 200.0e-6, -300.0e-6)
    camera.SetFocalPoint(0, 0, 0)
    camera.ComputeViewPlaneNormal()    
    renderer.SetActiveCamera(camera)
    renderer.ResetCamera()
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renwin)
    callback = TimerCallback(A, B, C)
    interactor.Initialize()
    interactor.AddObserver('TimerEvent', callback.execute)
    timerId = interactor.CreateRepeatingTimer(1)
    print 'Here'
    # renwin.FullScreenOn()
    interactor.Start()

class TimerCallback(object):
    def __init__(self,  dictA, dictB, dictC):
        self.A = dictA
        self.B = dictB
        self.C = dictC
        print 'Created callback'

    def execute(self, obj, event):
        self.A['pos'] += np.random.normal(0, 1, self.A['pos'].shape) * particleSpeed
        # Replace positions outside the X boundary by nearest boundary position
        x1 = np.where((self.A['pos'][:,0] < 0), 0, self.A['pos'][:,0])
        x2 = np.where(x1 > compartment['length'], compartment['length'], x1)
        self.A['pos'][:,0] = x2
        # Fix y-z boundary
        yzlen = np.sqrt(self.A['pos'][:,1]**2 + self.A['pos'][:,2]**2)    # Distance from origin in yz plane
        bad = np.nonzero(yzlen > compartment['radius'])
        # if y0,z0 is the position in yz plane, and angle of the vector is
        # `a` with y axis, then y1 = R * cos(a) = R * y/r and z1 = R *
        # sin(a) = R * z/r, where R is radius of cylinder, r is length of
        # the original position vector.
        self.A['pos'][bad, 1] *= compartment['radius'] / yzlen[bad]
        self.A['pos'][bad, 2] *= compartment['radius'] / yzlen[bad]
        # We want to ignore particles which have been excluded - exclusion
        # is indicated by 'ss'=0
        deadA = np.nonzero(self.A['ss'] <= 0)[0]
        self.A['pos'][deadA,:] = -1e9 # put large value to remove chance of interaction with B

        self.B['pos'] += np.random.normal(0, 1, self.B['pos'].shape) * particleSpeed
        x1 = np.where(self.B['pos'][:,0] < 0, 0, self.B['pos'][:,0])
        x2 = np.where(x1 > compartment['length'], compartment['length'], x1)
        self.B['pos'][:,0] = x2
        yzlen = np.sqrt(self.B['pos'][:,1]**2 + self.B['pos'][:,2]**2)
        bad = np.nonzero(yzlen > compartment['radius'])
        self.B['pos'][bad,1] *= compartment['radius'] / yzlen[bad]
        self.B['pos'][bad,2] *= compartment['radius'] / yzlen[bad]
        deadB = np.nonzero(self.B['ss'] <= 0)[0]
        self.B['pos'][deadB,:] = 1e9
        self.C['pos'] += np.random.normal(0, 1, self.C['pos'].shape) * particleSpeed
        x1 = np.where(self.C['pos'][:,0] < 0, 0, self.C['pos'][:,0])
        x2 = np.where(x1 > compartment['length'], compartment['length'], x1)
        self.C['pos'][:,0] = x2
        yzlen = np.sqrt(self.C['pos'][:,1]**2 + self.C['pos'][:,2]**2)
        bad = np.nonzero(yzlen > compartment['radius'])
        self.C['pos'][bad,1] *= compartment['radius'] / yzlen[bad]
        self.C['pos'][bad,2] *= compartment['radius'] / yzlen[bad]
        # Now calculate the pairwise distances between A and B molecules
        abdist = dist.cdist(self.A['pos'], self.B['pos'], 'euclidean')
        closeA, closeB = np.nonzero(abdist < bindingRadius)
        deadC = np.nonzero(self.C['ss'] <= 0.0)[0]
        if len(closeA)  > 0:
            # Multiple A can overlap with one B and vice versa, only the
            # smaller number participates in reaction
            if closeA.shape[0] < closeB.shape[0]:
                fewReacParticle = self.A
                fewReacIndex = np.unique(closeA)
            else:
                fewReacParticle = self.B
                fewReacIndex = np.unique(closeB)
            # print 'Number of close molecules:', closeA.shape[0], closeB.shape[0], np.unique(closeA).shape[0], np.unique(closeB).shape[0]
            # Fill C positions sequentially - number of entries only
            # increases from 0
            self.C['pos'][deadC[0]:deadC[0]+fewReacIndex.shape[0],:] = fewReacParticle['pos'][fewReacIndex,:]
            self.A['ss'][closeA] = 0.0
            self.B['ss'][closeB] = 0.0
            self.C['ss'][deadC[0]:deadC[0]+fewReacIndex.shape[0]] = 1.0
        posA = vnp.numpy_to_vtk(self.A['pos'], deep=True)
        # print posA
        # print posB
        self.A['points'].SetData(posA)
        sizesA = vnp.numpy_to_vtk(self.A['ss'], deep=True)
        self.A['polydata'].GetPointData().SetScalars(sizesA)
        posB = vnp.numpy_to_vtk(self.B['pos'], deep=True)
        self.B['points'].SetData(posB)
        sizesB = vnp.numpy_to_vtk(self.B['ss'], deep=True)
        self.B['polydata'].GetPointData().SetScalars(sizesB)
        posC = vnp.numpy_to_vtk(self.C['pos'], deep=True)
        self.C['points'].SetData(posC)
        sizesC = vnp.numpy_to_vtk(self.C['ss'], deep=True)
        self.C['polydata'].GetPointData().SetScalars(sizesC)
        obj.GetRenderWindow().Render()
    
if __name__ == '__main__':
    visualize()
