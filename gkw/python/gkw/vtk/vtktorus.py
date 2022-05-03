import numpy as np
import vtk
from math import pi,cos,sin

def threedtorusvtk(eps):
  """
  Output a vtk torus mesh given eps.
  Just a demonstration at the moment.
  """

  R = 1.
  Z = 0.
  npoints = 21

  ## create the points in XYZPoints
  xyz=np.array([1.,1.,1.],dtype=float)
  XYZPoints=vtk.vtkPoints()
  XYZPoints.SetNumberOfPoints(npoints*npoints*npoints)
  l=0
  for i in range(0,npoints):
    for j in range(0,npoints):
      for k in range(0,npoints):

	## construct the point
        theta = 2.*pi*j/(1.*(npoints-1));
	phi   = 2.*pi*k/(1.*(npoints-1));
	r=R*(1.0+eps*cos(theta));
	xyz[0] = R*cos(phi)*(1.0+eps*cos(theta));
        xyz[1] = R*sin(phi)*(1.0+eps*cos(theta));
        xyz[2] = Z + eps*sin(theta);

        ## add the point
        XYZPoints.InsertPoint(l,xyz)
        l=l+1

  ## create the point data in PointData
  PointData = vtk.vtkFloatArray()
  PointData.SetName("Potential")
  for k in range(0,npoints):
    for j in range(0,npoints):
      for i in range(0,npoints):
        _ = PointData.InsertNextValue(1.*(i+j+k))

  ## create a structured grid using the PointData as scalars
  ## and the XYZPoints.
  torusmesh=vtk.vtkStructuredGrid()
  torusmesh.SetDimensions(npoints,npoints,npoints)
  torusmesh.SetPoints(XYZPoints)
  ## (qq is a pointer for the place where the PointData should go)
  qq=torusmesh.GetPointData()
  qq.SetScalars(PointData)

  ## write the structured grid to vtkXML
  writer2=vtk.vtkXMLStructuredGridWriter()
  writer2.SetFileName("Torusmesh2.vts")
  writer2.SetInput(torusmesh)
  writer2.Write()
