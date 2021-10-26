#!/usr/bin/env python3    
 
 
from __future__ import print_function
import vtk
from numpy import random,genfromtxt,size
import numpy as np
import PointCloudFunctions_2016_12_31a as PCF
from vtk.util import numpy_support
import FT_20170106a as FT
import lsf_circle_20160719a as lsfc
from multiprocessing import Process
import os

from tkinter import *
from tkinter import filedialog
from tkinter import ttk
import tkinter as tk

import tkinterEntryWidget2 as dt



#
#python -m numpy.f2py -c --compiler=msvc FT_20160901a.f90 -m FT_20160901a

# TO DO:
# Sphere needs to exist on 2d plots as well as 3d plot.
# Ideally, we would keep track of how much X, Y, and Z rotation has occurred.
# check on distribution of grouping tasks.

# Up and down needs some feedback, potentially number boxes.
# Need labels for overal cutouts
# Display angle 
# Step size is not updating
# Need a counter for the number of features
 
class SphericalIndicator():
    def __init__(self, index = 0, increment = 10):
        self.index = 0
        self.increment = 10
        self.Xpos = 0
        self.Ypos =0
        self.Zpos = 0
     
#http://stackoverflow.com/questions/19973941/creating-2d-array-of-classes-in-python     
class FeatureLog() :
    def __init__(self) :
        #self.FeatureName = np.empty(self.state_max, dtype = "S50")
        self.FeatureName = "Empty"
        self.Description = "Empty"
        self.a = 0.00
        self.b = 0.00
        self.c = 0.00
        self.d = 0.00
        self.e = 0.00
        self.f = 0.00
        
        
        
             
class VtkPointCloud:
    def __init__(self, X, Y, Z, Kmin=-10.0, Kmax=10.0, maxNumPoints=1e7):
        ti = PCF.TimeStamp(0,"Start")
        ti1 = ti
        #print("VTKPC_00")
        #self.Xdisp = X
        #self.Ydisp = Y
        #self.Zdisp = Z

        
        self.X = X
        self.Y = Y
        self.Z = Z
        
        ti2 = PCF.TimeStamp(0,"Start")
        self.BigD = len(X)
        

        self.Xin = np.empty(self.BigD, np.float)
        self.Yin = np.empty(self.BigD, np.float)
        self.Zin = np.empty(self.BigD, np.float)
        
        self.Xin[:] = X[:]
        self.Yin[:] = Y[:]
        self.Zin[:] = Z[:]
        
        self.maxNumPoints = maxNumPoints
        self.vtkPolyData = vtk.vtkPolyData()
        self.vtkPolyData2 = vtk.vtkPolyData()
        self.clearPoints()



        mapper = vtk.vtkPolyDataMapper()
        # Next line has been added on.

        mapper.SetInputData(self.vtkPolyData)
        mapper.SetColorModeToDefault()

        mapper.SetScalarVisibility(1)
        self.vtkActor = vtk.vtkActor()
        self.vtkActor.SetMapper(mapper)
        
        mapper2 = vtk.vtkPolyDataMapper()
        # Next line has been added on.

        mapper2.SetInputData(self.vtkPolyData2)
        mapper2.SetColorModeToDefault()

        mapper2.SetScalarVisibility(1)
        self.vtkActor2 = vtk.vtkActor()
        self.vtkActor2.SetMapper(mapper2)       
        
        ti2 = PCF.TimeStamp(ti2,"After5 Initializations")
        self.IndexArray = np.empty(len(self.X), dtype=np.bool_)
        self.IndexArray0 = np.empty(len(self.X), dtype=np.bool_)
        self.IndexArray1 = np.empty(len(self.X), dtype=np.bool_)
        self.IndexArray2 = np.empty(len(self.X), dtype=np.bool_)
        self.IndexArray3 = np.empty(len(self.X), dtype=np.bool_)
        #self.IndexArray9 = np.empty(len(self.X), dtype=np.bool_)
        
        self.IndexSlice = np.empty(len(self.X), dtype=np.bool_)
        
        self.IndexArray[:] = True
        self.IndexArray0[:] = True
        self.IndexArray1[:] = True
        self.IndexArray2[:] = True
        self.IndexArray3[:] = True
        
        self.IndexSlice[:] = True
        
        #print("Length and sum of Index Array ", len(self.IndexArray), self.IndexArray.sum())
        #print("Length and sum of Index Array0 ", len(self.IndexArray0), self.IndexArray0.sum())
        #print("Length and sum of Index Array1 ", len(self.IndexArray1), self.IndexArray1.sum())
        #print("Length and sum of Index Array2 ", len(self.IndexArray2), self.IndexArray2.sum())
        
        self.IndexCircle = np.empty(len(self.X), dtype=np.bool_)
        self.IndexCircle[:] = True
        self.Index_X = np.empty(len(self.X), dtype=np.bool_)
        self.Index_X[:] = True
        self.Index_Y = np.empty(len(self.X), dtype=np.bool_)
        self.Index_Y[:] = True
        self.Index_Z = np.empty(len(self.X), dtype=np.bool_)
        self.Index_Z[:] = True
        
        
  
        
        
        
        self.IndexDispA = np.empty(len(self.X), dtype=np.bool_)
        self.IndexDispA[:] = True
        self.IndexDisp = np.empty(len(self.X), dtype=np.bool_)
        self.IndexDisp[:] = True
        self.IndexComb = np.empty(len(self.X), dtype=np.bool_)
        self.IndexComb[:] = True
        self.XaxisLoQualifier = np.empty(len(self.X), dtype=np.bool_)
        self.XaxisHiQualifier = np.empty(len(self.X), dtype=np.bool_)
        self.ZaxisHiDensityBand = np.empty(len(self.X), dtype=np.bool_)
        self.XaxisQualifier = np.empty(len(self.X), dtype=np.bool_)
        self.Evaluate = np.empty(len(self.X), dtype=np.bool_)
        self.Rsqd = np.empty(len(self.X), dtype=np.float64)
        self.IndexOrder = np.empty(len(self.X), dtype = np.int64)
        self.PtDisplaySize = 1.0
        self.HistFirstRun = True
        self.HistRange = 10
        self.Rmin = 0.94
        self.G0max = 0.4
        self.G1max = 0.4
        self.Kmin = 0.4 * (self.Z.max() - self.Z.min()) + self.Z.min()
        self.Kmax = 0.6 * (self.Z.max() - self.Z.min()) + self.Z.min()
        self.Sph_index = 0
        self.Sph_increment = 10
        self.Sph_Xpos = 0
        self.Sph_Ypos = 0
        self.Sph_Zpos = 0
        self.Sph_Radius = 6

        self.g0 = np.empty(self.BigD, dtype = np.float64)
        self.g1 = np.empty(self.BigD, dtype = np.float64)
        self.g2 = np.empty(self.BigD, dtype = np.float64)
        self.g3 = np.empty(self.BigD, dtype = np.float64)

        self.i1 = np.empty(self.BigD, dtype = np.float64)
        self.i2 = np.empty(self.BigD, dtype = np.float64)
        self.i3 = np.empty(self.BigD, dtype = np.float64)
        self.i4 = np.empty(self.BigD, dtype = np.float64)
        
        self.g0avg = np.empty(self.BigD, dtype = np.float64)
        self.g1avg = np.empty(self.BigD, dtype = np.float64)
        self.g2avg = np.empty(self.BigD, dtype = np.float64)
        
        self.iproj =  np.empty(self.BigD, dtype = np.float64)
        self.jproj =  np.empty(self.BigD, dtype = np.float64)
        self.kproj =  np.empty(self.BigD, dtype = np.float64)
        
        self.density_big = np.empty(self.BigD, dtype = np.float64)
                        
        self.NextIndex = np.empty(self.BigD, dtype = np.int64)
        self.PtsTotal = np.empty(self.BigD, dtype = np.int64)
        #self.IndexBig = np.empty(self.BigD, dtype = np.int64)
        
        
        
        self.NextIndex[:] = -1
        self.PtsTotal[:] = -1
        
        self.factor = np.int(1)
        self.factor = max(int(self.factor / 10),1)
        self.LittleD_B = self.BigD
        self.i_disp = (self.factor)*np.array(range(np.int(self.BigD/self.factor)))
        self.i_orig = np.array(range(self.BigD))
        
        self.state = 0
        self.state_max = 7
        #This needs to be reset to 0 for full operation of this program.
        self.statenames = np.empty(self.state_max, dtype = "S75")
        self.statenames[0] = "Manually orient component"
        self.statenames[1] = "Highlight main plane"
        self.statenames[2] = "Analytically determine main plane and reorient"
        self.statenames[3] = "Move highlighter to the center of the origin dowel"
        self.statenames[4] = "Rotate part in XY plane to align second dowel"
        self.statenames[5] = "Identify points and orientation that might create a cylinder"
        self.statenames[6] = "Identify a plane and align it to be horizontal in bottom left quadrant."
        
        self.helptext1 = np.empty(self.state_max, dtype = "S75")
        self.helptext1[0] = "Once OK with rough position, hit RETURN"
        self.helptext1[1] = "Use u/U and l/L to adjust limits of slice along dominant axis"
        self.helptext1[2] = "Press space to group points"
        self.helptext1[3] = "Use directional arrows to move cursor to estimate cyl axis to plane intersection"
        self.helptext1[4] = "Rotate part in XY plane to align second dowel"
        self.helptext1[5] = "Align part in bottom left quadrant"
        self.helptext1[6] = "Use slice boundaries to identify points that would contribute to the plane"
        
        self.helptext2 = np.empty(self.state_max, dtype = "S75")
        self.helptext2[0] = "Touch R to position cloud to origin"
        self.helptext2[1] = "Hit Return when done"
        self.helptext2[2] = "Analytically determine main plane and reorient"
        self.helptext2[3] = "Move highlighter to the center of the origin dowel"
        self.helptext2[4] = "Rotate part in XY plane to align second dowel"
        self.helptext2[5] = "Move indicator ball to point where cylinder axis crosses plane"
        self.helptext2[6] = "Once the minimum number of points are highlighted in red, hit 'A' to align"
        
        self.helptext3 = np.empty(self.state_max, dtype = "S75")
        self.helptext3[0] = ""
        self.helptext3[1] = ""
        self.helptext3[2] = "Analytically determine main plane and reorient"
        self.helptext3[3] = "Use u/U and l/L to adjust limits of slice along dominant axis"
        self.helptext3[4] = "Rotate part in XY plane to align second dowel"
        self.helptext3[5] = "Press h or H to highlight points at smaller or larger radius"
        self.helptext3[6] = "Press p to perform a regressive fit of rotation"
        
        self.helptext4 = np.empty(self.state_max, dtype = "S75")
        self.helptext4[0] = "Adjust increment for rotate and translate"
        self.helptext4[1] = ""
        self.helptext4[2] = "Analytically determine main plane and reorient"
        self.helptext4[3] = "Move highlighter to the center of the origin dowel"
        self.helptext4[4] = "Rotate part in XY plane to align second dowel"
        self.helptext4[5] = "Press 'a' or 'A' to project points an analyse"
        self.helptext4[6] = "Press G to group results and Return to fing longest string"
        
        self.helptext5 = np.empty(self.state_max, dtype = "S75")
        self.helptext5[0] = "Use keys x, X, y, Y, z, and Z to rotate or translate about axes"
        self.helptext5[1] = "Use keys x, X, y, Y, z, and Z to rotate or translate about axes"
        self.helptext5[2] = "Analytically determine main plane and reorient"
        self.helptext5[3] = "Move highlighter to the center of the origin dowel"
        self.helptext5[4] = "Rotate part in XY plane to align second dowel"
        self.helptext5[5] = "Press space to record."
        self.helptext5[6] = "Press space to record."
        
        
        ti = PCF.TimeStamp(ti, "Initialize")
        self.MaxWorkingGroupSize = 300000
        self.Grouped = False
        self.radius = 10
        self.Dens_min = 32
        self.DensityOld = True
        self.X_rot = 0.0
        self.Y_rot = 0.0 
        self.Z_rot = 0.0
        self.X_off = 0.0
        self.Y_off = 0.0
        self.Z_off = 0.0
       
        self.rot_inc = 2
        self.off_inc = 1
        self.inc = 1.0
        self.inc_min = 0.01
        self.inc_max = 10
        self.rot_or_tran = 0 # 0 = rotate, 1 = translate
        
        self.outputfilepathandname = ""
        self.GraphStep = 20
        
        self.density_test = 0.0
        
        self.slice_min = 0.00
        self.slice_max = 100.00
        
        self.X_limit_lower = -10000
        self.X_limit_upper = 100000
        
        self.Y_limit_lower = -10000
        self.Y_limit_upper = 100000
        
        self.Z_limit_lower = -10000
        self.Z_limit_upper = 100000
        
        self.feature_index = 1
        
        self.Features = np.empty(50, dtype = FeatureLog)
        for i in range(50) :
            self.Features[i] = FeatureLog()
            
        self.dominant_axis = 2
        self.allow_dominant_axis_change = True
        
        self.i_ctr = -1.0
        self.j_ctr = -1.0
        
        self.x0 = -1.00
        self.y0 = -1.00
        self.z0 = -1.00
        self.RefineState = -1
        
#        self.dist_scale = 0.1
        self.Z_err1_lim = 0.1
        self.Z_err2_lim = 1.0
        
        self.collapse_offset = 0.0
        self.collapse_repeat = 0.00
        self.collapse = False
        self.collapse_offset_last = 0.00
        self.collapse_repeat_last = 0.00
        
    def CalculateDensityProjection(self, Display, xc2, yc2, zc2, xn2, yn2, zn2, length_cam) :
    
        xi = 0.00
        yi = 0.00
        zi = 0.00
        xj = 0.00
        yj = 0.00
        zj = 0.00
        
        coeffs = np.empty(6, dtype = np.float64)

        
#        FT.project_pts_on_plane(self.X, self.Y, self.Z, self.iproj, self.jproj, self.kproj, self.IndexComb, xc2, yc2, zc2, \
#            xn2, yn2, zn2, xi, yi, zi, xj, yj, zj, self.BigD)
            
        FT.project_pts_on_plane(self.X, self.Y, self.Z, self.iproj, self.jproj, self.kproj, self.IndexComb, xc2, yc2, zc2, \
            xn2, yn2, zn2, coeffs, self.BigD)
        



        xi = coeffs[0]
        yi = coeffs[1]
        zi = coeffs[2]
        xj = coeffs[3]
        yj = coeffs[4]
        zj = coeffs[5]
        
            
#       print("Calculate Density xi, yi, zi, xj, yj, zj ", xi, yi, zi, xj, yj, zj)

        self.xi = xi
        self.yi = yi
        self.zi = zi
        self.xj = xj
        self.yj = yj
        self.zj = zj   

#        print("Calculate Density xi, yi, zi, xj, yj, zj ", xi, yi, zi, xj, yj, zj)
        # self.xk = xk
        # self.yk = yk
        # self.zk = zk        
 
        littled = self.IndexComb.sum()
#        print("Length of kproj ", len(self.kproj), " little d = ", littled)
        #print('littled', littled)
        self.iproj_l = np.empty(littled, dtype = np.float64)
        self.jproj_l = np.empty(littled, dtype = np.float64)
        self.kproj_l = np.empty(littled, dtype = np.float64)
        #print("littled ", littled)
        FT.crunch(self.iproj, self.iproj_l, self.IndexComb, self.BigD, littled)
        FT.crunch(self.jproj, self.jproj_l, self.IndexComb, self.BigD, littled)
        FT.crunch(self.kproj, self.kproj_l, self.IndexComb, self.BigD, littled)
        #print("littled ", littled, len(self.iproj_l), len(self.jproj_l))
        xcam1 = xn2 * length_cam + xc2
        ycam1 = yn2 * length_cam + yc2 
        zcam1 = zn2 * length_cam + zc2
        Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1) 
        Display.renderer.GetActiveCamera().SetFocalPoint(xc2, yc2, zc2)        

        
        self.density_big[:] = 0.0

        FT.estimate_density(self.density_big, self.iproj_l, self.jproj_l, self.kproj_l, self.IndexComb, littled, self.BigD)
        #self.pcld.density_big = PCF.EstimateDensity(self.pcld.iproj, self.pcld.jproj, self.pcld.IndexComb)
        density_test = self.density_big.sum()
        
        return(density_test)

        
    def CalculateDensityProjection2(self, Display, xc, yc, zc, xn, yn, zn, xu, yu, zu, length_cam) :
    
        # create unit vectors for current coordinate system
        xj, yj, zj = FT.unitize(xu, yu, zu)
        xk, yk, zk = FT.unitize(xn, yn, zn)
        xi, yi, zi = FT.crossproduct(xj, yj, zj, xk, yk, zk)
        xi, yi, zi = FT.unitize(xi, yi, zi)
        
        print("CDP2 a")
#        FT.project_pts_on_plane2(self.X, self.Y, self.Z, self.iproj, self.jproj, self.kproj, self.IndexComb, \
#            xc, yc, zc, xi, yi, zi, xj, yj, zj, xk, yk, zk, self.BigD)
 
 
        FT.project_pts_on_plane2(self.X, self.Y, self.Z, self.iproj, self.jproj, self.kproj, self.IndexComb, \
            0.0,0.0,0.0, xi, yi, zi, xj, yj, zj, xk, yk, zk, self.BigD)

        print("CDP2 b") 
#       print("Calculate Density xi, yi, zi, xj, yj, zj ", xi, yi, zi, xj, yj, zj)

        self.x0 = xc
        self.y0 = yc
        self.z0 = zc
        self.xi = xi
        self.yi = yi
        self.zi = zi
        self.xj = xj
        self.yj = yj
        self.zj = zj   
     
        print("CDP2 c") 
        
        littled = self.IndexComb.sum()
        self.iproj_l = np.empty(littled, dtype = np.float64)
        self.jproj_l = np.empty(littled, dtype = np.float64)
        self.kproj_l = np.empty(littled, dtype = np.float64)
        #print("littled ", littled)
        FT.crunch(self.iproj, self.iproj_l, self.IndexComb, self.BigD, littled)
        FT.crunch(self.jproj, self.jproj_l, self.IndexComb, self.BigD, littled)
        FT.crunch(self.kproj, self.kproj_l, self.IndexComb, self.BigD, littled)
        #print("littled ", littled, len(self.iproj_l), len(self.jproj_l))
        xcam1 = xn * length_cam + xc
        ycam1 = yn * length_cam + yc 
        zcam1 = zn * length_cam + zc
        Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1) 
        Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)        

        print("CDP2 d")        
        self.density_big[:] = 0.0

        print("CDP2 e")
        FT.estimate_density(self.density_big, self.iproj_l, self.jproj_l, self.kproj_l, self.IndexComb, littled, self.BigD)
        print("CDP2 f")
        #self.pcld.density_big = PCF.EstimateDensity(self.pcld.iproj, self.pcld.jproj, self.pcld.IndexComb)
        density_test = self.density_big.sum()
        
        return(density_test)

        
    def addPoint(self, point):
        #print("add point")
        #if self.vtkPoints.GetNumberOfPoints() < self.maxNumPoints:
        pointId = self.vtkPoints.InsertNextPoint(point[:])
        self.vtkDepth.InsertNextValue(point[2])
        self.vtkCells.InsertNextCell(1)
        self.vtkCells.InsertCellPoint(pointId)
        #self.vtkCells.Modified()
        #self.vtkPoints.Modified()
        #self.vtkDepth.Modified()
        
        #if self.vtkPoints2.GetNumberOfPoints() < self.maxNumPoints:
        pointId2 = self.vtkPoints2.InsertNextPoint(point[:])
        self.vtkDepth2.InsertNextValue(point[2])
        self.vtkCells2.InsertNextCell(1)
        self.vtkCells2.InsertCellPoint(pointId2)
        #self.vtkCells2.Modified()
        #self.vtkPoints2.Modified()
        #self.vtkDepth2.Modified()
        #self.vtkActor.GetProperty().SetPointSize(self.PtDisplaySize)
         
    def clearPoints(self):
        tcp = PCF.TimeStamp(0, "Start")
        self.vtkPoints = vtk.vtkPoints()
        self.vtkCells = vtk.vtkCellArray()
        self.vtkDepth = vtk.vtkDoubleArray()
        self.vtkDepth.SetName('DepthArray')
        self.vtkPolyData.SetPoints(self.vtkPoints)
        self.vtkPolyData.SetVerts(self.vtkCells)
        self.vtkPolyData.GetPointData().SetScalars(self.vtkDepth)
        self.vtkPolyData.GetPointData().SetActiveScalars('DepthArray')
        
        
        self.vtkPoints2 = vtk.vtkPoints()
        self.vtkCells2 = vtk.vtkCellArray()
        self.vtkDepth2 = vtk.vtkDoubleArray()
        self.vtkDepth2.SetName('Array2')
        self.vtkPolyData2.SetPoints(self.vtkPoints2)
        self.vtkPolyData2.SetVerts(self.vtkCells2)
        self.vtkPolyData2.GetPointData().SetScalars(self.vtkDepth2)
        self.vtkPolyData2.GetPointData().SetActiveScalars('Array2')
        
        tcp = PCF.TimeStamp(tcp,"Cleared Points")
        
    def Undo(self):
        #print("Length and sum of Index Array ", len(self.IndexArray), self.IndexArray.sum())
        #print("Length and sum of Index Array0 ", len(self.IndexArray0), self.IndexArray0.sum())
        #print("Length and sum of Index Array1 ", len(self.IndexArray1), self.IndexArray1.sum())
        #print("Length and sum of Index Array2 ", len(self.IndexArray2), self.IndexArray2.sum())
        self.IndexArray[:] = self.IndexArray1[:]
        self.IndexArray0[:] = self.IndexArray1[:]
        self.IndexArray1[:] = self.IndexArray2[:]
        self.IndexArray2[:] = self.IndexArray3[:]
        #print("Length and sum of Index Array ", len(self.IndexArray), self.IndexArray.sum())
        
    def addPoints(self):
 
        tap = PCF.TimeStamp(0, "Start")
        tap2 = tap

        if (not(np.array_equiv(self.IndexArray, self.IndexArray0))) :
            self.IndexArray3[:] = self.IndexArray2[:]
            self.IndexArray2[:] = self.IndexArray1[:]
            self.IndexArray1[:] = self.IndexArray0[:]
            self.IndexArray0[:] = self.IndexArray[:]
           
        #tap2 = PCF.TimeStamp(tap2, "Check arrays")
        
        self.IndexDispA[:] = False
        self.IndexDispA[self.i_disp] = True
        
        #print("length of IndexDispA", self.IndexDispA.sum())
        if (self.X_limit_upper > self.X_limit_lower) :
            self.Index_X = np.logical_and(np.where(np.less(self.X_limit_lower ,self.X), True,False), \
                np.where(np.less(self.X ,self.X_limit_upper), True,False))
        else :
            self.Index_X = np.logical_or(np.where(np.less(self.X_limit_lower ,self.X), True,False), \
                np.where(np.less(self.X ,self.X_limit_upper), True,False))

        self.Index_Y = np.logical_and(np.where(np.less(self.Y_limit_lower ,self.Y), True,False), \
            np.where(np.less(self.Y ,self.Y_limit_upper), True,False))
        if (self.Z_limit_upper > self.Z_limit_lower) :
            self.Index_Z = np.logical_and(np.where(np.less(self.Z_limit_lower ,self.Z), True,False), \
                np.where(np.less(self.Z ,self.Z_limit_upper), True,False))
        else :
            self.Index_Z = np.logical_or(np.where(np.less(self.Z_limit_upper ,self.Z), True,False), \
                np.where(np.less(self.Z ,self.Z_limit_lower), True,False))
                
        #print("length of Index_X", self.Index_X.sum())
        #print("length of Index_Y", self.Index_Y.sum())
        #print("length of Index_Y", self.Index_Z.sum())
        #print("self.X_limit_lower, self.X_limit_upper ", self.X_limit_lower, self.X_limit_upper)
        #print("self.Y_limit_lower, self.Y_limit_upper ", self.Y_limit_lower, self.Y_limit_upper)

        
        #print("length of IndexDispA", self.IndexDispA.sum())
        self.IndexDispA = np.logical_and(self.IndexDispA, \
            self.Index_X)

        #print("length of IndexDispA", self.IndexDispA.sum())            
        self.IndexDispA = np.logical_and(self.IndexDispA, \
            self.Index_Y)
                        
        self.IndexDispA = np.logical_and(self.IndexDispA, \
            self.Index_Z)
            
        #print("length of IndexDispA", self.IndexDispA.sum())

        
#        LittleD_A = len(self.i_disp)
        LittleD_A = np.long(self.IndexDispA.sum())
        # self.XdispA = self.X[self.i_disp]
        # self.YdispA = self.Y[self.i_disp]
        # self.ZdispA = self.Z[self.i_disp]
        
        if (LittleD_A > 0) :
            self.XdispA = np.zeros(LittleD_A)
            self.YdispA = np.zeros(LittleD_A)
            self.ZdispA = np.zeros(LittleD_A)
            FT.crunch_all(self.X, self.Y, self.Z, self.XdispA, self.YdispA, self.ZdispA, self.IndexDispA, self.BigD, LittleD_A)
        
            #print("length of XdispA ", len(self.XdispA))
        #self.XdispA = PCF.AR(self.X, self.i_disp)
        #self.YdispA = PCF.AR(self.Y, self.i_disp)
        #self.ZdispA = PCF.AR(self.Z, self.i_disp)
        
        #LittleD_B = max(self.i_disp)
        #self.XdispA = np.empty(LittleD_B, np.float64)
        #self.YdispA = np.empty(LittleD_B, np.float64)
        #self.ZdispA = np.empty(LittleD_B, np.float64)
             
        #XDA = Process(target = PCF.AR, args = (self.X, self.i_disp, self.XdispA))
        #YDA = Process(target = PCF.AR, args = (self.Y, self.i_disp, self.YdispA))
        #ZDA = Process(target = PCF.AR, args = (self.Z, self.i_disp, self.ZdispA))
        
        #XDA.start()
        #YDA.start()
        #ZDA.start()
        
        #XDA.join()
        #YDA.join()
        #ZDA.join()
        
        #tap2 = PCF.TimeStamp(tap2, "Reduce arrays A")
        #tap3 = tap2
                
        self.IndexDisp[:] = False

        #print("IndexArray.sum ", self.IndexArray.sum())
        #print("IndexComb.sum ", self.IndexComb.sum())
        self.IndexDisp[self.i_disp] = True  
        #print("IndexDisp.sum ", self.IndexDisp.sum())
        #print("len(IndexDisp) ", len(self.IndexDisp))
        
#        self.IndexDisp = np.logical_and(self.IndexDisp, self.Index_X)  
#        self.IndexDisp 
#        self.IndexDisp = np.logical_and(self.IndexArray, self.IndexDispA)         
        self.IndexComb = np.logical_and(self.IndexArray, self.IndexDispA)    
        LittleD_B = np.long(self.IndexComb.sum())
        self.LittleD_B = LittleD_B
        #tap3 = PCF.TimeStamp(tap3, "Combine logical arrays")
        
        self.XdispB = np.zeros(LittleD_B)
        self.YdispB = np.zeros(LittleD_B)
        self.ZdispB = np.zeros(LittleD_B)
        #tap3 = PCF.TimeStamp(tap3, "Setting up arrays")
        #FT.crunch(self.X, self.XdispB, self.IndexComb, self.BigD, LittleD_B)
        #FT.crunch(self.Y, self.YdispB, self.IndexComb, self.BigD, LittleD_B)
        #FT.crunch(self.Z, self.ZdispB, self.IndexComb, self.BigD, LittleD_B)
        
        #print("LittledD_B", LittleD_B )
        if (LittleD_B > 0) :
            FT.crunch_all(self.X, self.Y, self.Z, self.XdispB, self.YdispB, self.ZdispB, self.IndexComb, self.BigD, LittleD_B)
        #FT.crunch2(self.X, self.Y, self.Z, self.XdispB, self.YdispB, self.ZdispB, self.IndexComb, self.BigD, LittleD_B)
        #XFT = Process(target = FT.crunch, args = (self.X, self.XdispB, self.IndexComb, self.BigD, self.LittleD_B))
        #YFT = Process(target = FT.crunch, args = (self.Y, self.YdispB, self.IndexComb, self.BigD, self.LittleD_B))
        #ZFT = Process(target = FT.crunch, args = (self.Z, self.XdispB, self.IndexComb, self.BigD, self.LittleD_B))
            

            
        #XFT.start()
        #YFT.start()
        #ZFT.start()
        #XFT.join()
        #YFT.join()
        #ZFT.join()
        
        
        #tap3 = PCF.TimeStamp(tap3, "Crunch")
        #tap2 = PCF.TimeStamp(tap2, "Reduce arrays, B")
        
            
        LittleD_A = len(self.XdispA)
        #print("Number of points in index array", self.IndexArray.sum())
        data = np.empty((LittleD_A, 3), dtype = np.float64)

        
        data[:,0] = self.XdispA
        data[:,1] = self.YdispA
        data[:,2] = self.ZdispA
        #tap2 = PCF.TimeStamp(tap2, "Build data A")
        
        #numPoints = len(self.X) 
        numPoints = LittleD_A
     
        self.vtkPoints.Allocate(numPoints)
        table = numpy_support.numpy_to_vtk(data, deep = 2, array_type=None)
        #tap2 = PCF.TimeStamp(tap2, "Covert data to table A")
        self.vtkPoints.SetData(table)
        self.vtkActor.GetMapper().SetColorModeToDefault()
        self.vtkActor.GetProperty().SetColor(0.0, 0.0, 1.0)
        self.vtkActor.GetProperty().SetPointSize(self.PtDisplaySize)
        self.vtkPoints.Modified()
        self.vtkPoints.Modified()
        self.vtkDepth.Modified()
        
        #tap2 = PCF.TimeStamp(tap2, "Set up points, A")
        
        #print("self.HistFirstRun = ", self.HistFirstRun)
        #if (not(self.HistFirstRun)) :
        if (LittleD_B > 0) :
            

            #LittleD_B = len(self.XdispB)
            data2 = np.empty((LittleD_B, 3), dtype = np.float64)
            tap2 = PCF.TimeStamp(tap2, "Initialize data2")
            
           #print("Hist:  Sum of Index Array", self.IndexArray.sum())
            data2[:,0] = self.XdispB
            data2[:,1] = self.YdispB
            data2[:,2] = self.ZdispB
            #numPoints = len(self.X2) 
         
            #tap2 = PCF.TimeStamp(tap2, "Build data 2")
            self.vtkPoints2.Allocate(LittleD_B)
            table2 = numpy_support.numpy_to_vtk(data2, deep = 2, array_type=None)
            self.vtkPoints2.SetData(table2)
            self.vtkActor2.GetMapper().SetColorModeToDefault()
            self.vtkActor2.GetProperty().SetColor(1.0, 0.0, 0.0)
            self.vtkActor2.GetProperty().SetPointSize(self.PtDisplaySize+1)
            self.vtkPoints2.Modified()
            self.vtkPoints2.Modified()
            self.vtkDepth2.Modified()
            #tap = PCF.TimeStamp(tap,"Added Points")
            #tap2 = PCF.TimeStamp(tap2, "Set up points B")
            
            return()
            
            LittleD = np.long(self.IndexArray.sum())
            Xsel = np.empty(LittleD, dtype = np.float64)
            Ysel = np.empty(LittleD, dtype = np.float64)
            Zsel = np.empty(LittleD, dtype = np.float64)
            BigD = np.long(len(self.X))
            
            FT.crunch(self.X, Xsel, self.IndexArray, BigD, LittleD)
            FT.crunch(self.Y, Ysel, self.IndexArray, BigD, LittleD)
            FT.crunch(self.Z, Zsel, self.IndexArray, BigD, LittleD)
            
            #print("LittleD =", LittleD)
             
            #print("self.Xsel ", len(Xsel), len(Ysel), len(Zsel))
            G2 = np.zeros(4, np.float64)
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran Fit plane')
            #Gauss.fitplane(Xsel, Ysel, Zsel, G2, LittleD)
            FT.fitplane(Xsel, Ysel, Zsel, G2, LittleD)
            ##print("G0 = ", G2[0], " G1 = ", G2[1], " G2 = ", G2[2])
            #print("Regression analysis: ", G2[3])
                                              
    def CreatePointCloud(self):   
        #print("CreatePointCloud_00")
        point = np.zeros(3)
        tm = PCF.TimeStamp(0, 'Prior to calling add point in loop')
        
        data = np.empty((len(self.X), 3), dtype = np.float64)
        data[:,0] = self.X
        data[:,1] = self.Y
        data[:,2] = self.Z
            
        for k in range(len(self.X)):
            #point[0] = self.X[k]
            #point[1] = self.Y[k]
            #point[2] = self.Z[k]
          
            point = data[k,:]
            self.addPoint(point)   
            
        #print("CreatePointCloud_01b")
        tm = PCF.TimeStamp(tm, 'After to calling add point in loop')
        self.clearPoints
        self.addPoints
        return self
        
    def AxisRotate(self, Axis, Angle_deg):
        ax = PCF.TimeStamp(0,"Start")
        Angle_rad = Angle_deg * np.pi/180
   
        if Axis == 0: FT.rotatepoints(self.Y, self.Z, Angle_rad, self.BigD)   
        if Axis == 1: FT.rotatepoints(self.X, self.Z, Angle_rad, self.BigD)   
        if Axis == 2: FT.rotatepoints(self.X, self.Y, Angle_rad, self.BigD)  
        ax = PCF.TimeStamp(ax,"Rotating about axis")
        return self
        
    # def AxisRotate2(self, Axis, Angle_deg):
        # ax = PCF.TimeStamp(0,"Start")
        
   
        # if Axis == 0: 
            # self.X_rot += Angle_deg
            # self.X_rot = np.mod(self.X_rot, 360)
            # Angle_rad = self.X_rot * np.pi/180
            # FT.rotatepoints2(self.Yin, self.Zin, self.Y, self.Z, self.Y_off, self.Z_off, Angle_rad, self.BigD)   
        # if Axis == 1: 
            # self.Y_rot += Angle_deg
            # self.Y_rot = np.mod(self.Y_rot, 360)
            # Angle_rad = self.Y_rot * np.pi/180
            # FT.rotatepoints2(self.Xin, self.Zin, self.X, self.Z, self.X_off, self.Z_off, Angle_rad, self.BigD)   
        # if Axis == 2: 
            # self.Z_rot += Angle_deg
            # self.Z_rot = np.mod(self.Z_rot, 360)
            # Angle_rad = self.Z_rot * np.pi/180
            # FT.rotatepoints2(self.Xin, self.Yin, self.X, self.Y, self.X_off, self.Y_off, Angle_rad, self.BigD)  
        # ax = PCF.TimeStamp(ax,"Rotating about axis")
        # return self
        
        
    def collapse_X(self) :
        print ("collapse_X", self.collapse_repeat, self.collapse_offset)
        self.X[:] = np.mod(self.Xin[:], self.collapse_repeat) + self.collapse_offset
        
    def SetOnOrigin(self):    
        self.X = self.X - self.X.min()
        self.Y = self.Y - self.Y.min()
        self.Z = self.Z - self.Z.min()
        
    def Save(self):   
    
            #    outputfilepathandname
        if (self.outputfilepathandname != "") :
            PCF.OutputData(self.outputfilepathandname,self.XdispB,self.YdispB,self.ZdispB, len(self.XdispB))

    def ExportFeatures(self):
        if (self.outputfilepathandname != "") :
            outfile = open(self.outputfilepathandname, "w")
            i = 0
            
            while i < self.feature_index :
                if i > 0 :
                    a_str = self.Features[i].FeatureName
                    b_str = self.Features[i].Description
                    c_str = '%.4f' % self.Features[i].a
                    d_str = '%.4f' % self.Features[i].b
                    e_str = '%.4f' % self.Features[i].c
                    f_str = '%.4f' % self.Features[i].d
                    g_str = '%.4f' % self.Features[i].e
                    h_str = '%.4f' % self.Features[i].f
                else :
                    a_str = 'Geometry'
                    b_str = 'Description'
                    c_str = 'a'
                    d_str = 'b'
                    e_str = 'c'
                    f_str = 'd'
                    g_str = 'e'
                    h_str = 'f'

                LineToWrite = a_str + '\t' + \
                    b_str + '\t' + \
                    c_str + '\t' + \
                    d_str + '\t' + \
                    e_str + '\t' + \
                    f_str + '\t' + \
                    g_str + '\t' + \
                    h_str + '\n'
                outfile.write(LineToWrite)
                i = i + 1
            outfile.close()

            
    def ExportProjection(self):
        if (self.outputfilepathandname != "") :
            outfile = open(self.outputfilepathandname, "w")
            i = 0
            
            a_str = '%.4f' % self.x0
            b_str = '%.4f' % self.y0
            c_str = '%.4f' % self.z0
          
            LineToWrite = a_str + '\t' + \
                b_str + '\t' + \
                c_str + '\n'
            outfile.write(LineToWrite)
            
            
            a_str = '%.4f' % self.i_ctr
            b_str = '%.4f' % self.j_ctr
            c_str = ''
          
            LineToWrite = a_str + '\t' + \
                b_str + '\t' + \
                c_str + '\n'
            outfile.write(LineToWrite)
            
            a_str = '%.4f' % self.xi
            b_str = '%.4f' % self.yi
            c_str = '%.4f' % self.zi
          
            LineToWrite = a_str + '\t' + \
                b_str + '\t' + \
                c_str + '\n'
            outfile.write(LineToWrite)
            
            a_str = '%.4f' % self.xj
            b_str = '%.4f' % self.yj
            c_str = '%.4f' % self.zj
          
            LineToWrite = a_str + '\t' + \
                b_str + '\t' + \
                c_str + '\n'
            outfile.write(LineToWrite)
                            
            Length = len(self.iproj_l)
            print("Length = ", Length)
            while i < Length :
                a_str = '%.4f' % self.iproj_l[i]
                b_str = '%.4f' % self.jproj_l[i]
                c_str = '%.4f' % self.kproj_l[i]

                LineToWrite = a_str + '\t' + \
                    b_str + '\t' + \
                    c_str + '\n'
                outfile.write(LineToWrite)
                i = i + 1
            outfile.close()
            
    def SaveGroupResults(self) :
    
        if (self.outputfilepathandname != "") :
            outfile = open(self.outputfilepathandname, "w")
            LittleD = self.IndexArray.sum()
            index_sav = np.empty(LittleD, dtype = np.float64)    
            Xsav = np.empty(LittleD, dtype = np.float64)                
            Xsav = np.empty(LittleD, dtype = np.float64)
            Ysav = np.empty(LittleD, dtype = np.float64)
            Zsav = np.empty(LittleD, dtype = np.float64)
            g0sav = np.empty(LittleD, dtype = np.float64)
            g1sav = np.empty(LittleD, dtype = np.float64)
            g2sav = np.empty(LittleD, dtype = np.float64)
            g3sav = np.empty(LittleD, dtype = np.float64)
            i1sav = np.empty(LittleD, dtype = np.float64)
            i2sav = np.empty(LittleD, dtype = np.float64)
            i3sav = np.empty(LittleD, dtype = np.float64)
            i4sav = np.empty(LittleD, dtype = np.float64)
            g0avgsav = np.empty(LittleD, dtype = np.float64)
            g1avgsav = np.empty(LittleD, dtype = np.float64)
            g2avgsav = np.empty(LittleD, dtype = np.float64)
            NextIndexsav = np.empty(LittleD, dtype = np.float64)
            PtsTotalsav = np.empty(LittleD, dtype = np.float64)
                
                
            BigD = self.BigD
            FT.crunch(self.i_orig, index_sav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.X, Xsav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.Y, Ysav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.Z, Zsav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.g0, g0sav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.g1, g1sav, self.IndexArray, BigD, LittleD) 
            FT.crunch(self.g2, g2sav, self.IndexArray, BigD, LittleD) 
            FT.crunch(self.g3, g3sav, self.IndexArray, BigD, LittleD) 
            FT.crunch(self.i1, i1sav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.i2, i2sav, self.IndexArray, BigD, LittleD) 
            FT.crunch(self.i3, i3sav, self.IndexArray, BigD, LittleD) 
            FT.crunch(self.i4, i4sav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.g0avg, g0avgsav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.g1avg, g1avgsav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.g2avg, g2avgsav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.NextIndex, NextIndexsav, self.IndexArray, BigD, LittleD)
            FT.crunch(self.PtsTotal, PtsTotalsav, self.IndexArray, BigD, LittleD)

#            outfile = open("GroupedResults.xyz", "w")
            i = 0
            LineToWrite = "Index X Y Z G0 G1 G2 G3 i1 i2 i3 i4 G0avg G1avg G2avg NextIndex PtsTotal \n" 
            outfile.write(LineToWrite)
            while i < LittleD:
                a_str = '%.0f' % index_sav[i]
                b_str = '%.4f' % Xsav[i]
                c_str = '%.4f' % Ysav[i]
                d_str = '%.4f' % Zsav[i]
                e_str = '%.4f' % g0sav[i]
                f_str = '%.4f' % g1sav[i]
                g_str = '%.4f' % g2sav[i]
                h_str = '%.4f' % g3sav[i]
                i_str = '%.0f' % i1sav[i]
                j_str = '%.0f' % i2sav[i]
                k_str = '%.0f' % i3sav[i]
                l_str = '%.0f' % i4sav[i]
                m_str = '%.4f' % g0avgsav[i]
                n_str = '%.4f' % g1avgsav[i]
                o_str = '%.4f' % g2avgsav[i]
                p_str = '%.0f' % NextIndexsav[i]
                q_str = '%.0f' % PtsTotalsav[i]

                

                LineToWrite = a_str + '\t' + b_str + '\t' + c_str + '\t' + d_str + \
                    '\t' + e_str + '\t' + f_str + '\t' + g_str + '\t' + h_str + '\t' + \
                    i_str + '\t' + j_str + '\t' + k_str + '\t' + l_str + '\t' + m_str + \
                    '\t' + n_str + '\t' + o_str +\
                    '\t' + p_str + '\t' + q_str + '\n'
                    
                outfile.write(LineToWrite)
                i = i + 1
            outfile.close()
                
                # #Now only select index values that are in the longest chain
                # self.IndexArray[:] = False
                # i = Pts_Max_Index
                # while (i > -1) :
                    # self.IndexArray[i] = True
                    # i = self.NextIndex[i]
                              
            
    
    def PtSizePlus(self):
        if self.PtDisplaySize < 10 :
            self.PtDisplaySize += 1
       #print("self.DisplaySize = ", self.PtDisplaySize)
        self.vtkActor.GetProperty().SetPointSize(self.PtDisplaySize)
        self.vtkActor2.GetProperty().SetPointSize(self.PtDisplaySize+1)
        
    def RangePlus(self):
        if (self.HistRange < len(self.X)/4) : self.HistRange = self.HistRange * 2
        
    def RangeMinus(self) :
        self.HistRange = self.HistRange / 2
       
    def PtSizeMinus(self):
        self.PtDisplaySize = self.PtDisplaySize - 1
        self.PtDisplaySize = max(self.PtDisplaySize, 1.0)
       #print("self.DisplaySize = ", self.PtDisplaySize)
        self.vtkActor.GetProperty().SetPointSize(self.PtDisplaySize)
        self.vtkActor2.GetProperty().SetPointSize(self.PtDisplaySize+1)
                
    def ReduceDisplayFactor(self):
        self.factor = max(int(self.factor / 10),1)
        self.i_disp = (self.factor)*np.array(range(np.int(self.BigD/self.factor)))
        #print("Length of i_disp", self.i_disp)
          
    def IncreaseDisplayFactor(self):
        self.factor = int(self.factor * 10)
        self.factor = min(self.factor, np.int(self.BigD/100))
        self.i_disp = (self.factor)*np.array(range(np.int(self.BigD/self.factor)))
        #print("Length of i_disp", self.i_disp)
       
    def Group(self) :
 
            BigD = len(self.X) 
            Dist = np.empty(BigD, dtype = np.float64)            
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran subroutine distance')
            
            FT.distance(self.X, self.Y, self.Z, Dist, self.IndexArray, BigD)
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran subroutine distance')
            FT.distance2(self.X, self.Y, self.Z, Dist, BigD)
            
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran subroutine distance')

            #group_points(p, q, r, groups, stilltosort, gsize, npoints)
            #print("Hello, in Group")
            IndexArray_Real = np.zeros(len(self.IndexArray), np.float64)
            IndexArray_Real[:] = 0.0
            #G_Results = np.empty(len(self.X), np.float64)
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran subroutine group')
            self.clearPoints()
            self.addPoints()
            ###########
            LittleD_B = self.IndexComb.sum()
            if (np.pi > 4) :
                LittleD_B = self.IndexComb.sum()
                if (LittleD_B > 1) :
                    if (LittleD_B < 1E5 ) :
                        #use XdispB, YdispB, ZdispB and create other arrays reduced indentically by self.IndexComb
                        FT.crunch(self.i_orig, self.i_b, self.BigD, LittleD_B)

                        j = 0
                        self.IndexArray[:] = False
                        for i in range(0,LittleD_B - 1) :
                            if (self.iar[i] > 0.5) :
                                j = self.i_b[i]
                                self.g0[j] = self.g0B[i]
                                self.g1[j] = self.g1B[i]
                                self.g2[j] = self.g2B[i]
                                self.g1[j] = self.g1B[i]
                                self.i1
                            
                    else :
                        print("Try to reduce the number of points.  Over 1E5 points will take too long")
                ###########
            
            #print("length of i2", len(self.i2))
            FT.group_points2(self.X, self.Y, self.Z, self.g0, self.g1, self.g2, self.g3, self.i1, self.i2, self.i3, self.i4, self.IndexComb, IndexArray_Real, \
                self.G0max, self.G1max, self.Rmin, 4, LittleD_B, len(self.X))
                                                    
            #print("Maximum of self.g0", self.g0.max())
            tm = PCF.TimeStamp(tm, 'After calling Fortran subroutine group')
            #print("Returned from group_points")
           #print(self.IndexArray)
           
            self.IndexArray[:] = False
            self.IndexArray = np.where(np.less(0.5, IndexArray_Real), True, False)
            
           #print(self.IndexArray)
            #print("sum of IndexArray_Real and sum of self.IndexArray ", sum(IndexArray_Real), sum(self.IndexArray))
           #print("Length of self.IndexArray", len(self.IndexArray))
            self.clearPoints()
            self.addPoints()
            self.Grouped = True
   
       
    def Group2(self) :
 
            BigD = len(self.iproj) 
            Dist = np.empty(BigD, dtype = np.float64)            
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran subroutine distance')
            
#            FT.distance(self.iproj, self.jproj, self.kproj, Dist, self.IndexArray, BigD)
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran subroutine distance')
            FT.distance2(self.iproj, self.jproj, self.kproj, Dist, BigD)
            
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran subroutine distance')

            #group_points(p, q, r, groups, stilltosort, gsize, npoints)
            #print("Hello, in Group")
            IndexArray_Real = np.zeros(len(self.IndexArray), np.float64)
            IndexArray_Real[:] = 0.0
            #G_Results = np.empty(len(self.X), np.float64)
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran subroutine group')
            self.clearPoints()
            self.addPoints()
            ###########
            LittleD_B = self.IndexComb.sum()

            #print("length of i2", len(self.i2))
            FT.group_points2(self.kproj, self.iproj, self.jproj, self.g0, self.g1, self.g2, self.g3, self.i1, self.i2, self.i3, self.i4, self.IndexComb, IndexArray_Real, \
                self.G0max, self.G1max, self.Rmin, 4, LittleD_B, len(self.X))
                                                    
            #print("Maximum of self.g0", self.g0.max())
            tm = PCF.TimeStamp(tm, 'After calling Fortran subroutine group')
            #print("Returned from group_points")
           #print(self.IndexArray)
           
            self.IndexArray[:] = False
            self.IndexArray = np.where(np.less(0.5, IndexArray_Real), True, False)
            
           #print(self.IndexArray)
            #print("sum of IndexArray_Real and sum of self.IndexArray ", sum(IndexArray_Real), sum(self.IndexArray))
           #print("Length of self.IndexArray", len(self.IndexArray))
            self.clearPoints()
            self.addPoints()
            self.Grouped = True

            
    def Group3_trash(self) :
 
            BigD = len(self.iproj_l) 
            Dist = np.empty(BigD, dtype = np.float64)            
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran subroutine distance')
            
            FT.distance2(self.iproj_l, self.jproj_l, self.kproj_l, Dist, BigD)
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran subroutine distance')
            
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran subroutine distance')

            #group_points(p, q, r, groups, stilltosort, gsize, npoints)
            #print("Hello, in Group")
            IndexArray_Real = np.zeros(len(self.IndexArray), np.float64)
            IndexArray_Real[:] = 0.0
            #G_Results = np.empty(len(self.X), np.float64)
            tm = PCF.TimeStamp(0, 'Prior to calling Fortran subroutine group')
            self.clearPoints()
            self.addPoints()
            ###########
            LittleD_B = self.IndexComb.sum()

            self.g0_l = np.empty(BigD, dtype = np.float64)
            self.g1_l = np.empty(BigD, dtype = np.float64)
            self.g2_l = np.empty(BigD, dtype = np.float64)
            self.g3_l = np.empty(BigD, dtype = np.float64)
            self.i1_l = np.empty(BigD, dtype = np.float64) 
            self.i2_l = np.empty(BigD, dtype = np.float64)
            self.i3_l = np.empty(BigD, dtype = np.float64)            
            self.i4_l = np.empty(BigD, dtype = np.float64)            
            
 
            FT.group_points2(self.iproj_l, self.jproj_l, self.kproj_l, self.g0, self.g1, self.g2, self.g3, self.i1, self.i2, self.i3, self.i4, self.IndexComb, IndexArray_Real, \
                self.G0max, self.G1max, self.Rmin, 4, LittleD_B, len(self.X))
                                                    
            #print("Maximum of self.g0", self.g0.max())
            tm = PCF.TimeStamp(tm, 'After calling Fortran subroutine group')
            #print("Returned from group_points")
           #print(self.IndexArray)
           
            self.IndexArray[:] = False
            self.IndexArray = np.where(np.less(0.5, IndexArray_Real), True, False)
            
           #print(self.IndexArray)
            #print("sum of IndexArray_Real and sum of self.IndexArray ", sum(IndexArray_Real), sum(self.IndexArray))
           #print("Length of self.IndexArray", len(self.IndexArray))
            self.clearPoints()
            self.addPoints()
            self.Grouped = True
             
    def FindLongestString(self):
        LittleD = np.long(self.IndexArray.sum())
        next_index_r = np.empty(self.BigD, np.float64)
        indexarray_r = np.empty(self.BigD, np.float64)
        pts_total_r = np.empty(self.BigD, np.float64)
        #subroutine findlongeststring(p, q, r, g0, g1, g2, g3, g0avg, g1avg, g2avg, next_index, pts_total, &
                                #stilltosort, indexarray_r, dist_max, littled, bigd, index_pt)
                              
        tm = PCF.TimeStamp(0, 'Prior to running large loop')
        FT.findlongeststring(self.X, self.Y, self.Z, self.g0, self.g1, self.g2, self.g3, self.g0avg, self.g1avg, \
                            self.g2avg, next_index_r, pts_total_r, self.IndexArray, indexarray_r, \
                            self.Z_err1_lim, self.Z_err2_lim, \
                            LittleD, self.BigD)
                            
        self.IndexArray = np.where(np.less(0.5, indexarray_r), True, False)
        self.NextIndex = np.int64(next_index_r)
        self.PtsTotal = np.int64(pts_total_r)
        tm = PCF.TimeStamp(tm, 'After running FT.findlongeststring')
        #print("sum of IndexArray_Real and sum of self.IndexArray ", sum(indexarray_r), sum(self.IndexArray))

    def FindLongestString2(self):
        LittleD = np.long(self.IndexArray.sum())
        next_index_r = np.empty(self.BigD, np.float64)
        indexarray_r = np.empty(self.BigD, np.float64)
        pts_total_r = np.empty(self.BigD, np.float64)

                              
        tm = PCF.TimeStamp(0, 'Prior to running large loop')
        # Note that in the command below, the order of the coordinates has changed.  This is because that
        # function expects the third parameter to be the 'height' or the explicitly defined parameter.
        FT.findlongeststring(self.kproj, self.iproj, self.jproj, self.g0, self.g1, self.g2, self.g3, \
            self.g0avg, self.g1avg, self.g2avg, next_index_r, pts_total_r, self.IndexArray, indexarray_r, \
            self.Z_err1_lim, self.Z_err2_lim, \
            LittleD, self.BigD)
                            
        self.IndexArray = np.where(np.less(0.5, indexarray_r), True, False)
        self.NextIndex = np.int64(next_index_r)
        self.PtsTotal = np.int64(pts_total_r)
        tm = PCF.TimeStamp(tm, 'After running FT.findlongeststring')
    
                            
    def Hist(self):
        if (self.HistFirstRun) :
            self.HistFirstRun = False
           #print("Starting Histogram")
            #global XcontextView
            data = np.empty((len(self.X), 3), dtype = np.float64)
            data[:,0] = self.X
            data[:,1] = self.Y
            data[:,2] = self.Z
            
            bins = 20
            #identify a slice of XY data from the Z axis which represents a sudden increase in density
            Xmin = self.X.min()
            Xmax = self.X.max()
            Ymin = self.Y.min()
            Ymax = self.Y.max()
            Kmin = self.Z.min()
            Kmax = self.Z.max()
            Zstep = (Kmax - Kmin)/bins
            #print("Zstep = ", Zstep)
            MassArray = np.empty(len(self.X))
            DensityArray = np.empty(bins-1)
            DensityArray[:] = 0
            #Mass = np.empty((len(self.X), 3), dtype = np.float64)
            #print(Ymin, Ymax)
            tm = PCF.TimeStamp(0, 'Prior to finding heaviest band')
            for i in range(0,bins-1):
                Zctr = (i+1)*Zstep + Kmin
                #print(i, Zctr)
                #bool_test = np.logical_and(np.less(self.Z,Zctr+Zstep/2), np.less(Zctr-Zstep/2,self.Z))
                MassArray = np.where(np.logical_and(np.less(self.Z,Zctr+Zstep/2), np.less(Zctr-Zstep/2,self.Z)), 1, MassArray)
                DensityArray[i] = MassArray.sum()/len(self.X)
                #print(i, MassArray)
                MassArray = 0

            tm = PCF.TimeStamp(tm, 'After finding heaviest band')
            #print(DensityArray)
            #print("Location of maximum = ", np.argmax(DensityArray))
            #print("Location in mm = ",(np.argmax(DensityArray)+1)*Zstep + Kmin)
            
            
            Index = np.argmax(DensityArray)
            Zctr = (Index+1)*Zstep + Kmin
            self.ZaxisHiDensityBand = np.where(np.logical_and(np.less(self.Z,Zctr+Zstep), np.less(Zctr-Zstep,self.Z)), True,False)
            
            tm = PCF.TimeStamp(tm, 'Where loop')
            
            #Create masks which allow only points at the lo and hi ends of the data
            self.XaxisLoQualifier = np.where(np.less(self.X,200), True,False)
            self.XaxisHiQualifier = np.where(np.less(300,self.X), True,False)
            
            self.XaxisQualifier = np.logical_or(self.XaxisLoQualifier, self.XaxisHiQualifier)
            self.IndexArray = np.logical_and(self.ZaxisHiDensityBand , self.XaxisQualifier)
               
        #self.IndexArray[:] = True
           #print("FinishingHistogramStage 1, IndexArray.sum = ", self.IndexArray.sum())
        
        else :
            LittleD = np.long(self.IndexArray.sum())
           #print("LittleD =", LittleD)
            self.Xsel = np.empty(LittleD, )
            self.Ysel = np.zeros(LittleD)
            self.Zsel = np.zeros(LittleD)
            BigD = np.long(len(self.X))
            tm = PCF.TimeStamp(0, 'Starting Crunch')
            FT.crunch(self.X, self.Xsel, self.IndexArray, BigD, LittleD)
            FT.crunch(self.Y, self.Ysel, self.IndexArray, BigD, LittleD)
            FT.crunch(self.Z, self.Zsel, self.IndexArray, BigD, LittleD)
             
            tm = PCF.TimeStamp(tm, 'Completed Crunch x 3')
            
            #(G, Rsqd) = PCF.FitPlane(self.Xsel, self.Ysel, self.Zsel, LittleD)
            tm = PCF.TimeStamp(tm, 'Completed finding plane in PCF')
            #print("G = ", G)
            #print("Rsqd = ", Rsqd)
            #print("LittleD =", LittleD)
                
            G2 = np.zeros(4, np.float64)
            tm = PCF.TimeStamp(tm, 'Prior to calling Fortran Fit plane')
            #Gauss.fitplane(self.Xsel, self.Ysel, self.Zsel, G2, LittleD)
            FT.fitplane(self.Xsel, self.Ysel, self.Zsel, G2, LittleD)
            tm = PCF.TimeStamp(tm, 'Completed finding plane, Fortran')
            #print("G2 = ", G2)
            
            
            Rsqd_max = 0
            j = 0
            self.Rsqd[:] = 0
            for i in range(BigD):
    # #        for i in range(1) :
                if (self.IndexArray[i]) :
                    self.Evaluate[:] = True
                    self.Evaluate[i] = False
                    #print("Length of self.Evaluate ", self.Evaluate.sum())
                    self.Evaluate = np.logical_and(self.IndexArray, self.Evaluate)
                    #print("Length of self.Evaluate ", self.Evaluate.sum())

                    EvenSmallerD = np.long(self.Evaluate.sum())

                    self.Rsqd[i] = 0
                    if EvenSmallerD > 0 :
                        #print("EvenSmallerD", EvenSmallerD)
                        #print("length of Xsel", len(self.Xsel))
                        #print("Evaluate = ", self.Evaluate)
                        #print("Evaluate.sum()", self.Evaluate.sum())
                        self.Xeval = np.empty(EvenSmallerD, dtype=np.float64)
                        self.Yeval = np.empty(EvenSmallerD, dtype=np.float64)
                        self.Zeval = np.empty(EvenSmallerD, dtype=np.float64)
                        FT.crunch(self.X, self.Xeval, self.Evaluate, BigD, EvenSmallerD)
                        FT.crunch(self.Y, self.Yeval, self.Evaluate, BigD, EvenSmallerD)
                        FT.crunch(self.Z, self.Zeval, self.Evaluate, BigD, EvenSmallerD)
                        #print("Xeval, Yeval, Zeval")
                        #print(self.Xeval, self.Yeval, self.Zeval)
                        #(G, Rsqd) = PCF.FitPlane(self.Xeval, self.Yeval, self.Zeval, EvenSmallerD)
                        #Gauss.fitplane(self.Xeval, self.Yeval, self.Zeval, G2, EvenSmallerD)
                        FT.fitplane(self.Xeval, self.Yeval, self.Zeval, G2, EvenSmallerD)
                        #print("G2 = ", G2)
                        self.Rsqd[i] = G2[3]
                        j = j + 1
                        #print("self.Rsqd[i]", self.Rsqd[i], "i = ", i)
                        #print("Calculations left", i, j)
                        #print("self.Rsqd ", self.Rsqd)
            tm = PCF.TimeStamp(tm, 'Completed loop to evaluate planes')

            #print("self.IndexArray", self.IndexArray)
            self.Rsqd = 1-self.Rsqd
            self.IndexOrder = np.argsort(self.Rsqd)
           #print("self.IndexOrder", self.IndexOrder[0:4])

            self.Evaluate[:] = True
            i_range = min(EvenSmallerD,self.HistRange) 
            for i in range(i_range):
                self.Evaluate[self.IndexOrder[i]] = False
           #print("self.Evaluate", self.Evaluate)
           #print(self.Evaluate.sum())
            self.IndexArray = np.logical_and(self.IndexArray, self.Evaluate)
           #print("IndexArraySum", self.IndexArray.sum())
           #print("self.IndexArray", self.IndexArray)
                
                

            self.Evaluate = np.empty(len(self.X), dtype=np.bool_)
            self.Rsqd = np.empty(len(self.X), dtype=np.float64)
            self.IndexOrder = np.empty(len(self.X), dtype = np.int_)
            
            self.length = len(self.X)
        self.clearPoints()
        self.addPoints()

# #        return(IndexArray)

    def InitBound(self) :
        if self.dominant_axis == 0 :
            self.Kmin = 0.4 * (self.X.max() - self.X.min()) + self.X.min()
            self.Kmax = 0.6 * (self.X.max() - self.X.min()) + self.X.min()
        if self.dominant_axis == 1 :
            self.Kmin = 0.4 * (self.Y.max() - self.Y.min()) + self.Y.min()
            self.Kmax = 0.6 * (self.Y.max() - self.Y.min()) + self.Y.min()
        if self.dominant_axis == 2 :
            self.Kmin = 0.4 * (self.Z.max() - self.Z.min()) + self.Z.min()
            self.Kmax = 0.6 * (self.Z.max() - self.Z.min()) + self.Z.min()
        
    def Bound(self):
        
        tm = PCF.TimeStamp(0, 'Prior to finding heaviest band')

        if (self.dominant_axis == 0) : 
            self.IndexArray = np.where(np.logical_and(np.less(self.X,self.Kmax), np.less(self.Kmin,self.X)), True,False)
        if (self.dominant_axis == 1) :
            self.IndexArray = np.where(np.logical_and(np.less(self.Y,self.Kmax), np.less(self.Kmin,self.Y)), True,False)
        if (self.dominant_axis ==2) :
            self.IndexArray = np.where(np.logical_and(np.less(self.Z,self.Kmax), np.less(self.Kmin,self.Z)), True,False)
        
        tm = PCF.TimeStamp(tm, 'Where loop')
        
        #Create masks which allow only points at the lo and hi ends of the data
  
        self.clearPoints()
        self.addPoints()

    def FindPlane(self):
        print("FindPlane")
        LittleD = np.long(self.IndexArray.sum())
        Xsel = np.empty(LittleD, dtype = np.float64)
        Ysel = np.empty(LittleD, dtype = np.float64)
        Zsel = np.empty(LittleD, dtype = np.float64)
        BigD = np.long(len(self.X))
        
        FT.crunch(self.X, Xsel, self.IndexArray, BigD, LittleD)
        FT.crunch(self.Y, Ysel, self.IndexArray, BigD, LittleD)
        FT.crunch(self.Z, Zsel, self.IndexArray, BigD, LittleD)
        
        #print("LittleD =", LittleD)
         
        #print("self.Xsel ", len(Xsel), len(Ysel), len(Zsel))
        G2 = np.zeros(4, np.float64)
        tm = PCF.TimeStamp(0, 'Prior to calling Fortran Fit plane')
        #Gauss.fitplane(Xsel, Ysel, Zsel, G2, LittleD)
        FT.fitplane(Xsel, Ysel, Zsel, G2, LittleD)
        
        #self.pcld.AxisRotate(1 , 1)

        self.Z = self.Z - G2[2]
        Angle_X = np.arctan(G2[1])
        Angle_Y = np.arctan(G2[0])
        print ("Angles = ", Angle_X, Angle_Y)
        self.Y,self.Z=PCF.RotatePoints(self.Y,self.Z,Angle_X)
        self.X,self.Z=PCF.RotatePoints(self.X,self.Z,Angle_Y)
        
        self.clearPoints()
        self.addPoints()
        return self
        
    def FindPlane2(self):
        print("FindPlane2")
        LittleD = np.long(self.IndexArray.sum())
        Xsel = np.empty(LittleD, dtype = np.float64)
        Ysel = np.empty(LittleD, dtype = np.float64)
        Zsel = np.empty(LittleD, dtype = np.float64)
        BigD = np.long(len(self.X))
        
        FT.crunch(self.X, Xsel, self.IndexArray, BigD, LittleD)
        FT.crunch(self.Y, Ysel, self.IndexArray, BigD, LittleD)
        FT.crunch(self.Z, Zsel, self.IndexArray, BigD, LittleD)
        
        #print("LittleD =", LittleD)
         
        #print("self.Xsel ", len(Xsel), len(Ysel), len(Zsel))
        G2 = np.zeros(4, np.float64)
        tm = PCF.TimeStamp(0, 'Prior to calling Fortran Fit plane')
        #Gauss.fitplane(Xsel, Ysel, Zsel, G2, LittleD)
        FT.fitplane(Xsel, Ysel, Zsel, G2, LittleD)
        
        print("Planar regression = ", G2)
        
        #self.pcld.AxisRotate(1 , 1)

        # self.Z = self.Z - G2[2]
        # Angle_X = np.arctan(G2[1])
        # Angle_Y = np.arctan(G2[0])
        # print ("Angles = ", Angle_X, Angle_Y)
        # self.Y,self.Z=PCF.RotatePoints(self.Y,self.Z,Angle_X)
        # self.X,self.Z=PCF.RotatePoints(self.X,self.Z,Angle_Y)
        
        # self.clearPoints()
        # self.addPoints()
        return self
        
        
    def FindCircle(self) :
        #self.IndexCircle
        tc = PCF.TimeStamp(0, 'Finding center of circle')
        #diameter = 30.00 #mm
        radius = self.radius
        
        p = np.empty(self.BigD, np.float64)
        q = np.empty(self.BigD, np.float64)
        r = np.empty(self.BigD, np.float64)
        
        p[:] = self.X[:]
        q[:] = self.Y[:]
        r[:] = 0.0
        
        p = p - self.Sph_Xpos
        q = q - self.Sph_Ypos
        #r = r - self.Sph_Zpos
        
        Dist = np.empty(self.BigD, dtype = np.float64)   
        FT.distance2(p, q, r, Dist, self.BigD)
        Dist = Dist**0.5
        print("Max of Dist ", Dist.max(), "Min of Dist ", Dist.min())
        
        #self.IndexCircle = False
        #self.IndexCircle = np.where((Dist < radius) , True, False)
        self.IndexCircle = np.where(np.logical_and(np.less(radius,Dist + 3), np.less(Dist,radius)), True,False)

        
       #print("Number of points within radius", self.IndexCircle.sum())
        self.i_crunched = np.empty(self.IndexCircle.sum(), np.int)
        #FT.crunch(self.i_orig, self.i_crunched, self.BigD, self.IndexCircle.sum())
        #self.IndexArray[:] = self.IndexCircle[:]
        self.IndexArray[:] = np.logical_and(self.IndexCircle[:], self.IndexSlice[:])
        
        #Find an area defined by an OD - ID of 2mm which causes the greatest change in point count
        #We will use the center of this area as the approximate center.
        #Min_DensityChange
        
        self.clearPoints()
        self.addPoints()
        
       #print("self.IndexCircle.sum() ", self.IndexCircle.sum(), self.IndexArray.sum())
       #print("self.XdispB, len ", len(self.XdispB))
        
    def FindCylinder(self, Display) :
    
        xc = self.Sph_Xpos
        yc = self.Sph_Ypos
        zc = self.Sph_Zpos
        
        self.xc = xc
        self.yc = yc
        self.zc = zc
        
        xn, yn, zn = Display.FindViewNormal() 
        
        # xc = self.xc
        # yc = self.yc
        # zc = self.zc
        # xn = self.xn
        # yn = self.yn
        # zn = self.zn
        
        self.IndexSlice = np.where(np.logical_and(np.less(self.slice_dist,self.slice_max), \
            np.less(self.slice_min,self.slice_dist)), True,False)
        #self.IndexCircle
       #print("Find Cylinder with xc, yc, zc, xn, yn, zn", xc, yc, zc, xn, yn, zn)
        tc = PCF.TimeStamp(0, 'Finding center of circle')
        #diameter = 30.00 #mm
        radius = self.radius
        
        p = np.empty(self.BigD, np.float64)
        q = np.empty(self.BigD, np.float64)
        r = np.empty(self.BigD, np.float64)
        
        p[:] = self.X[:]
        q[:] = self.Y[:]
        r[:] = self.Z[:]
        
        # p = p - self.Sph_Xpos
        # q = q - self.Sph_Ypos
        # #r = r - self.Sph_Zpos
        
        Dist = np.empty(self.BigD, dtype = np.float64)   
        self.IndexArray[:]= True
        FT.dist_pt_to_line(p, q, r, Dist, self.IndexArray, xc, yc, zc, xn, yn, zn, self.BigD)
        #Dist = Dist**0.5
       #print("Max of Dist ", Dist.max(), "Min of Dist ", Dist.min())
        
        #self.IndexCircle = False
        #self.IndexCircle = np.where((Dist < radius) , True, False)
        self.IndexCircle = np.where(np.logical_and(np.less(radius,Dist + 3), np.less(Dist,radius)), True,False)

        
       #print("Number of points within radius", self.IndexCircle.sum())
        self.i_crunched = np.empty(self.IndexCircle.sum(), np.int)
        #FT.crunch(self.i_orig, self.i_crunched, self.BigD, self.IndexCircle.sum())
        self.IndexArray[:] = np.logical_and(self.IndexCircle[:], self.IndexSlice[:])
        
        #Find an area defined by an OD - ID of 2mm which causes the greatest change in point count
        #We will use the center of this area as the approximate center.
        #Min_DensityChange
        
        self.clearPoints()
        self.addPoints()
        
       #print("self.IndexCircle.sum() ", self.IndexCircle.sum(), self.IndexArray.sum())
       #print("self.XdispB, len ", len(self.XdispB))

        #self.sphere.SetCenter(0, 0, 0)
        
        tc = PCF.TimeStamp(tc, 'After finding center of circle')
        
    #def FitCircle(self) :
    
    def HighlightPlane(self, Display) :

        self.IndexArray[:] = self.IndexSlice[:]

        self.clearPoints()
        self.addPoints()



    
    
            
    def temp2(self):
       #print("StartingFindPlane")
        #global XcontextView
        data = np.empty((len(self.X), 3), dtype = np.float64)
        data[:,0] = self.X
        data[:,1] = self.Y
        data[:,2] = self.Z
        
        bins = 40
        #identify a slice of XY data from the Z axis which represents a sudden increase in density
        Xmin = self.X.min()
        Xmax = self.X.max()
        Ymin = self.Y.min()
        Ymax = self.Y.max()
        Kmin = self.Z.min()
        Kmax = self.Z.max()
        Zstep = (Kmax - Kmin)/bins
        #print("Zstep = ", Zstep)
        MassArray = np.empty(len(self.X))
        DensityArray = np.empty(bins-1)
        #Mass = np.empty((len(self.X), 3), dtype = np.float64)
        #print(Ymin, Ymax)
        tm = PCF.TimeStamp(0, 'Prior to finding heaviest band')
        for i in range(0,bins-1):
            Zctr = (i+1)*Zstep + Kmin
            #print(i, Zctr)
            #bool_test = np.logical_and(np.less(self.Z,Zctr+Zstep/2), np.less(Zctr-Zstep/2,self.Z))
            MassArray = np.where(np.logical_and(np.less(self.Z,Zctr+Zstep/2), np.less(Zctr-Zstep/2,self.Z)), 1, MassArray)
            DensityArray[i] = MassArray.sum()/len(self.X)
            #print(i, MassArray)
            MassArray = 0

        tm = PCF.TimeStamp(tm, 'After finding heaviest band')
        #print(DensityArray)
        #print("Location of maximum = ", np.argmax(DensityArray))
        #print("Location in mm = ",(np.argmax(DensityArray)+1)*Zstep + Kmin)
        
        
        Index = np.argmax(DensityArray)
        Zctr = (Index+1)*Zstep + Kmin
        self.ZaxisHiDensityBand = np.where(np.logical_and(np.less(self.Z,Zctr+Zstep/2), np.less(Zctr-Zstep/2,self.Z)), True,False)
        
        tm = PCF.TimeStamp(tm, 'Where loop')
        
        #Create masks which allow only points at the lo and hi ends of the data
        self.XaxisLoQualifier = np.where(np.less(self.X,200), True,False)
        self.XaxisHiQualifier = np.where(np.less(300,self.X), True,False)
        
        self.XaxisQualifier = np.logical_or(self.XaxisLoQualifier, self.XaxisHiQualifier)
        self.IndexArray = np.logical_and(self.ZaxisHiDensityBand , self.XaxisQualifier)
        
        LittleD = np.long(self.IndexArray.sum())
       #print("LittleD =", LittleD)
        self.Xsel = np.zeros(LittleD)
        self.Ysel = np.zeros(LittleD)
        self.Zsel = np.zeros(LittleD)
        BigD = np.long(len(self.X))
        FT.crunch(self.X, self.Xsel, self.IndexArray, BigD, LittleD)
        FT.crunch(self.Y, self.Ysel, self.IndexArray, BigD, LittleD)
        FT.crunch(self.Z, self.Zsel, self.IndexArray, BigD, LittleD)
        tm = PCF.TimeStamp(tm, 'Completed Crunch x 3')
        
        (G, Rsqd) = PCF.FitPlane(self.Xsel, self.Ysel, self.Zsel, LittleD)
        tm = PCF.TimeStamp(tm, 'Completed finding plane')
       #print("G = ", G)
       #print("Rsqd = ", Rsqd)
       #print("LittleD =", LittleD)
    
        Rsqd_max = 0
        j = 0
        self.Rsqd[:] = 1
        for i in range(BigD):
#        for i in range(1) :
            if (self.IndexArray[i]) :
                self.Evaluate[:] = True
                self.Evaluate[i] = False
               #print("Length of self.Evaluate ", self.Evaluate.sum())
                self.Evaluate = np.logical_and(self.IndexArray, self.Evaluate)
               #print("Length of self.Evaluate ", self.Evaluate.sum())

                EvenSmallerD = np.long(self.Evaluate.sum())

                #self.Rsqd[i] = 0
    #            if EvenSmallerD > 0 :

                #print("EvenSmallerD", EvenSmallerD)
                #print("length of Xsel", len(self.Xsel))
                #print("Evaluate = ", self.Evaluate)
                #print("Evaluate.sum()", self.Evaluate.sum())
                self.Xeval = np.empty(EvenSmallerD, dtype=np.float64)
                self.Yeval = np.empty(EvenSmallerD, dtype=np.float64)
                self.Zeval = np.empty(EvenSmallerD, dtype=np.float64)
                FT.crunch(self.X, self.Xeval, self.Evaluate, BigD, EvenSmallerD)
                FT.crunch(self.Y, self.Yeval, self.Evaluate, BigD, EvenSmallerD)
                FT.crunch(self.Z, self.Zeval, self.Evaluate, BigD, EvenSmallerD)
                (G, Rsqd) = PCF.FitPlane(self.Xeval, self.Yeval, self.Zeval, EvenSmallerD)
                self.Rsqd[i] = Rsqd
                j = j + 1
               #print("self.Rsqd[i]", self.Rsqd[i])
               #print("Calculations left", i, j)
       #print("self.Rsqd ", self.Rsqd)

       #print("self.IndexArray", self.IndexArray)
        self.IndexOrder = np.argsort(self.Rsqd)
       #print("self.IndexOrder", self.IndexOrder)
        #self.Rsqd = 1-self.Rsqd
        self.Evaluate[:] = False
        for i in range(3):
            self.Evaluate[self.IndexOrder[i]] = True
       #print("self.Evaluate", self.Evaluate)
       #print(self.Evaluate.sum())
        self.IndexArray = np.logical_and(self.IndexArray, self.Evaluate)
       #print("IndexArraySum", self.IndexArray.sum())
       #print("self.IndexArray", self.IndexArray)
            
            

        #self.Evaluate = np.empty(len(self.X), dtype=np.bool_)
        #self.Rsqd = np.empty(len(self.X), dtype=np.float64)
        #self.IndexOrder = np.empty(len(self.X), dtype = int_)
        
        self.length = len(self.X)
       #print("FinishingFindPlane, IndexArray.sum = ", self.IndexArray.sum())

#        return(IndexArray)
             
    def Temp(self):

        #Xsav = np.zeros(LittleD)
        #Ysav = np.zeros(LittleD)
        #Zsav = np.zeros(LittleD)

        index_big = np.empty(BigD, dtype = np.float64)               
        i = 0
        while i < BigD :
            index_big[i] = np.float64(i)
            i = i + 1
        
        tm = PCF.TimeStamp(tm, 'Prior to running large loop')
        
        i = 0
        Pts_Max_Index = -1
        while i < BigD :
        # Go through the entire set of points
            if (self.IndexArray[i]) :
                #but only make calculations on the points selected.
                Found = False
                j = 0
                while  j <= i :
                    #starting at 0, go through all prior points to see if there are examples
                    # in which the current point is within some margin of the plane described by a previous point
                    #print("in grouping loop, j = ", j)
                    
                    if (self.PtsTotal[j] > 0) and (self.IndexArray[j]): 
                        Z_pred = self.X[i]*self.g0avg[j] + self.Y[i]*self.g1avg[j] + self.g2avg[j]
                        if abs(self.Z[i] - Z_pred) < 2.0 :
                            # add this point to the chain:
                            j_start = j
                            k = j
                            PtsCounter = 1
                            G0_tot = self.g0[k]
                            G1_tot = self.g1[k]
                            G2_tot = self.g2[k]
                            while self.NextIndex[k] > -1 :
                                PtsCounter += 1
                                k = self.NextIndex[k]
                                G0_tot += self.g0[k]
                                G1_tot += self.g1[k]
                                G2_tot += self.g2[k]
                               #print("k = ", k)
                            if (k != i) : self.NextIndex[k] = i
                            self.g0avg[j_start] = G0_tot / PtsCounter
                            self.g1avg[j_start] = G1_tot / PtsCounter
                            self.g2avg[j_start] = G2_tot / PtsCounter
                            self.PtsTotal[j_start] = PtsCounter
                            if (Pts_Max_Index > -1) :
                                if (PtsCounter > self.PtsTotal[Pts_Max_Index]) : Pts_Max_Index = j_start
                            else : 
                                Pts_Max_Index = j_start
                            Found = True
                        #print("Found = ", Found, "i = ", i, " j = ", j, " self.Z[i]", self.Z[i], "Z_pred ", Z_pred)
                    j += 1
                if not(Found) :
                    self.g0avg[i] = self.g0[i]
                    self.g1avg[i] = self.g1[i]
                    self.g2avg[i] = self.g2[i]
                    self.PtsTotal[i] = 1
               #print("in grouping loop, i = ", i)
            i = i + 1
        
        tm = PCF.TimeStamp(tm, 'After running large loop')
 
  

    def SphRev(self):
        increment = self.Sph_increment
        self.Sph_index = np.int(self.Sph_index - increment)
        while not(self.IndexArray[self.Sph_index]) :
            self.Sph_index = np.int(self.Sph_index - increment)
        if (self.Sph_index < 0) : self.Sph_index = 0
        print ("self.Sph_index = ", self.Sph_index)

    def SphFwd(self):
        increment = self.Sph_increment
        self.Sph_index = np.int(self.Sph_index + increment)
        while not(self.IndexArray[self.Sph_index]) :
            self.Sph_index = np.int(self.Sph_index + increment)
        if (self.Sph_index > len(self.IndexArray)) : self.Sph_index = len(self.IndexArray)
        print ("self.Sph_index = ", self.Sph_index)

    def AdvanceState(self) :
        self.state += 1
        self.state = min(self.state, self.state_max - 1)
        self.RefineState = 0
        print("State: ", self.state, " ", self.statenames[self.state])
                        
    def IncSphereSize(self) :
        if self.Sph_Radius < 20 : self.Sph_Radius += 1.00
        
    def DecSphereSize(self) :
        self.Sph_Radius -= 1.00
        if self.Sph_Radius < 1 : self.Sph_Radius = 1.0
        
        
    def IncInc(self) :
        self.inc = self.inc * 1.6
        if self.inc > self.inc_max : self.inc = self.inc_max
        
    def DecInc(self) :
        self.inc = self.inc / 1.6
        if self.inc < self.inc_min : self.inc = self.inc_min
        
class TwoDPlot():


    def __init__(self, view, HchartA,VchartA, HorTitle,VerTitle, GraphMinH, GraphMinV, GraphRange, GraphStep):
       #print("2d plot init_00")
        self.view = view
        self.HchartA = HchartA
        self.VchartA = VchartA
        self.HorTitle = HorTitle
        self.VerTitle = VerTitle
        self.chart = vtk.vtkChartXY()
        self.chart.SetShowLegend(False)
        self.tableA = vtk.vtkTable()  
        self.tableB = vtk.vtkTable()        
        #self.chart = vtk.vtkChartXY()
        #self.GraphMaxOverall = GraphMaxOverall
        #self.GraphMinOverall = GraphMinOverall
        self.GraphMinH = GraphMinH
        self.GraphMinV = GraphMinV
        self.GraphRange = GraphRange
        self.GraphStep = GraphStep
        self.PtDisplaySize = 2.0
        self.scale_Hor = 1
        self.scale_Ver = 1
       #print("2DPI_01")
         
        numPoints = len(HchartA)
        #print("Number of points ", numPoints)
        
        #Xtable = np.empty(numPoints, np.float64)
        #Ytable = np.empty(numPoints, np.float64)
       
        #self.Xtable = self.Xchart
        #self.Ytable = self.Ychart
        
        #self.tableA.SetNumberOfRows(LittleD_A)
        
        self.tableA.SetNumberOfRows(numPoints)
        self.tableB.SetNumberOfRows(numPoints)
       #print("2DPI_02")
         
        self.HarrA = numpy_support.numpy_to_vtk(self.HchartA, deep = 0, array_type=None)
        self.VarrA = numpy_support.numpy_to_vtk(self.VchartA, deep = 0, array_type=None)
       #print("2DPI_03:")
        self.HarrA.SetName('Horizontal')
        self.VarrA.SetName('Vertical')
        self.tableA.AddColumn(self.HarrA)
        self.tableA.AddColumn(self.VarrA)
        self.chart.ClearPlots()
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle(HorTitle)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetTitle(VerTitle)
        
       #print("2DPI_04:")
        self.points = self.chart.AddPlot(vtk.vtkChart.POINTS)
        self.points.SetInputData(self.tableA, 0, 1)
        self.points.SetColor(0, 0, 255, 255)
        self.points.SetWidth(0.1)
        self.points.SetMarkerSize(self.PtDisplaySize)
        self.points.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)
        
        #view = vtk.vtkContextView()
        self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        self.view.GetRenderWindow().SetSize(400, 300)
        self.view.GetScene().ClearItems()
        self.view.GetScene().AddItem(self.chart)
        
       #print("2DPI_05:")
        
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle(self.HorTitle)
        #chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMinimum(0)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetTitle(self.VerTitle)
          
        #view = vtk.vtkContextView()
        #view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        #view.GetRenderWindow().SetSize(400, 300)
        #view.GetScene().AddItem(chart)
        
        self.ResetScale(self.GraphMinH, self.GraphMinV, self.GraphRange, self.GraphStep)

        # #GraphMinH, GraphMinV, GraphRange, GraphStep
        # self.chart.GetAxis(vtk.vtkAxis.LEFT).SetMinimum(self.GraphMinV)
        # self.chart.GetAxis(vtk.vtkAxis.LEFT).SetMaximum(self.GraphMinV + self.GraphRange)
        # #self.chart.GetAxis(vtk.vtkAxis.LEFT).AutoScale()
        # self.chart.GetAxis(vtk.vtkAxis.LEFT).SetBehavior(vtk.vtkAxis.FIXED)
        
        # print("2DPI_06:")
        
        # self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMinimum(self.GraphMinH)
        # self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMaximum(self.GraphMinH + 2*self.GraphRange)
        # #self.chart.GetAxis(vtk.vtkAxis.BOTTOM).AutoScale()
        # self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetBehavior(vtk.vtkAxis.FIXED)
        # print("2DPI_07:")
        #return view
        
    def PtSizePlus(self):
        if self.PtDisplaySize < 10 :
            self.PtDisplaySize += 1
        print("self.DisplaySize = ", self.PtDisplaySize)
       
    def PtSizeMinus(self):
        self.PtDisplaySize = self.PtDisplaySize - 1
        self.PtDisplaySize = max(self.PtDisplaySize, 1.0)
        print("self.DisplaySize = ", self.PtDisplaySize)
        
    def ZoomIn(self):
        #self.GraphMaxOverall = self.GraphMaxOverall * 1.2
        #self.GraphMaxOverall = max(self.GraphMaxOverall, 1.0)
        self.GraphRange_New = min(self.GraphRange * 1.2, self.GraphRange + 2 * self.GraphStep)
        self.GraphRange = self.GraphStep * int(self.GraphRange_New / self.GraphStep)
        
        
        
        
        self.GraphMinV = min(self.GraphMinV - self.GraphStep, self.GraphMinV + self.GraphRange / 2 - self.GraphRange_New / 2)
        self.GraphMinV = self.GraphStep * int(self.GraphRange / self.GraphStep)
        self.GraphMinH = min(self.GraphMinH - self.GraphStep, self.GraphMinH + self.GraphRange / 2 - self.GraphRange_New / 2)
        self.GraphMinH = self.GraphStep * int(self.GraphRange / self.GraphStep)
        
        
       #print("self.GraphMaxOverall = ", self.GraphMaxOverall)
       
       
    def ZoomOut(self):
        self.GraphMaxOverall = self.GraphMaxOverall / 2
       #print("self.GraphMaxOverall = ", self.GraphMaxOverall)
    
    def ResetScale(self, GraphMinH, GraphMinV, GraphRange, GraphStep) :
    
        self.GraphMinH = GraphMinH
        self.GraphMinV = GraphMinV
        self.GraphRange = GraphRange
        self.GraphStep = GraphStep
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetMinimum(self.GraphMinV)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetMaximum(self.GraphMinV + self.GraphRange)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetBehavior(vtk.vtkAxis.FIXED)
        
        GraphRangeH = int(self.GraphRange * 1.2 / self.GraphStep) * self.GraphStep
        
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMinimum(self.GraphMinH)
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMaximum(self.GraphMinH + GraphRangeH)
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetBehavior(vtk.vtkAxis.FIXED)
        
    
#    def Update(self, HchartA, VchartA, HchartB,VchartB,i_disp, SelChart):
    def Update(self, HchartA, VchartA, HchartB,VchartB):
        #print("Update 2 d table")
         
        LittleD_A = len(HchartA)

        self.HchartA = HchartA
        self.VchartA = VchartA
        self.HchartB = HchartB
        self.VchartB = VchartB

        
        self.HarrA = vtk.vtkFloatArray()
        self.VarrA = vtk.vtkFloatArray()
        self.HarrA = numpy_support.numpy_to_vtk(self.HchartA, deep = 0, array_type=None)
        self.VarrA = numpy_support.numpy_to_vtk(self.VchartA, deep = 0, array_type=None)

        self.HarrA.SetName('Horizontal')
        self.VarrA.SetName('Vertical')
        self.tableA.SetNumberOfRows(LittleD_A)
        self.tableA.AddColumn(self.HarrA)
        self.tableA.AddColumn(self.VarrA)
        
        self.chart.ClearPlots()
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetTitle(self.HorTitle)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetTitle(self.VerTitle)
        
        points = self.chart.AddPlot(vtk.vtkChart.POINTS)   
        points.SetInputData(self.tableA, 0, 1)
        points.SetColor(0, 0, 255, 255)
        points.SetWidth(0.1)
        points.SetMarkerSize(self.PtDisplaySize)
        points.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)
        LittleD_B = len(self.HchartB)
        #LittleD = np.long(SelChart.sum())
        if (LittleD_B>0) : 
            
            HarrB = numpy_support.numpy_to_vtk(self.HchartB, deep = 0, array_type=None)
            VarrB = numpy_support.numpy_to_vtk(self.VchartB, deep = 0, array_type=None)
            HarrB.SetName('Hor_sel')
            VarrB.SetName('Ver_sel')
            self.tableB.SetNumberOfRows(LittleD_B)
            self.tableB.AddColumn(HarrB)
            self.tableB.AddColumn(VarrB)
            self.points_sel = self.chart.AddPlot(vtk.vtkChart.POINTS)
            self.points_sel.SetInputData(self.tableB, 0, 1)
            self.points_sel.SetColor(1.0, 0.0, 0.0)
            self.points_sel.SetWidth(0.1)
            self.points_sel.SetMarkerSize(self.PtDisplaySize+1)
              
        self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        self.view.GetRenderWindow().SetSize(400, 300)
        self.view.GetScene().ClearItems()
        self.view.GetScene().AddItem(self.chart)

        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetMinimum(self.GraphMinV)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetMaximum(self.GraphMinV + self.GraphRange)
        #self.chart.GetAxis(vtk.vtkAxis.LEFT).AutoScale()
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetAxisVisible(True)
        self.chart.GetAxis(vtk.vtkAxis.LEFT).SetBehavior(vtk.vtkAxis.FIXED)
        
        GraphRangeH = int(self.GraphRange * 1.2 / self.GraphStep) * self.GraphStep
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMinimum(self.GraphMinH)
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetMaximum(self.GraphMinH + GraphRangeH)
        
        #self.chart.GetAxis(vtk.vtkAxis.BOTTOM).AutoScale()
        self.chart.GetAxis(vtk.vtkAxis.BOTTOM).SetBehavior(vtk.vtkAxis.FIXED)
        
    # def ExportFeatures(self) :
        
        # self.feature_index = 0
        
        # for i in xrange(50) :
            # Features[i] = FeatureLog()


class MyInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):

 
    def __init__(self,parent, Display, pointCloud, SphInt, Sphere):
    #(self.renderWindowInteractor, self, pcl, self.SI, self.sphere))
        print ("MIS_01")
        self.AddObserver("LeftButtonPressEvent",self.leftButtonPressEvent)
        self.AddObserver("LeftButtonReleaseEvent",self.leftButtonReleaseEvent)
        self.AddObserver("MiddleButtonPressEvent",self.middleButtonPressEvent)
        self.AddObserver("MiddleButtonReleaseEvent",self.middleButtonReleaseEvent)
        self.AddObserver("MouseMoveEvent",self.mouseMoveEvent)
        self.AddObserver("KeyPressEvent", self.keyPress)
        #print("self= ", self)
        #print("parent = ", parent)
        self.RWI = parent
        self.pcld = pointCloud
        #self.SphInt = SphInt
        self.sphere = Sphere
        self.Display = Display
        print ("MIS_02")
        self.timer = PCF.TimeStamp(0,"Start")
        self.lastkey = ""
        
        xn, yn, zn = self.Display.renderer.GetActiveCamera().GetViewPlaneNormal()
        self.pcld.xn = xn
        self.pcld.yn = yn
        self.pcld.zn = zn
        self.mbpe = False
        
                
                
    def leftButtonPressEvent(self,obj,event):
        #print("Middle Button pressed")
        xn, yn, zn = self.Display.renderer.GetActiveCamera().GetViewPlaneNormal()
        self.pcld.xn = xn
        self.pcld.yn = yn
        self.pcld.zn = zn
#        self.OnLeftButtonDown()
        self.OnMiddleButtonDown()
        self.Display.calculate()
        return
        
    def leftButtonReleaseEvent(self,obj,event):
        #print("Middle Button pressed")
        xn, yn, zn = self.Display.renderer.GetActiveCamera().GetViewPlaneNormal()
        self.pcld.xn = xn
        self.pcld.yn = yn
        self.pcld.zn = zn
#        self.OnLeftButtonUp()
        self.OnMiddleButtonUp()
        self.Display.calculate()
        return
 
    def middleButtonPressEvent(self,obj,event):
        #print("Middle Button pressed")
#        self.OnMiddleButtonDown()
#        self.RWI.SetKeySym('Shift')
        key = self.RWI.GetKeySym()
        ctrl = self.RWI.GetControlKey()
        print(key)
       #print(ctrl)
#        if (key == 'Control_L') : 
        if (ctrl == 8) :
            self.OnMiddleButtonDown()
            self.mbpe = True
        else : self.OnLeftButtonDown()
        return
 
    def middleButtonReleaseEvent(self,obj,event):
        #print("Middle Button released")
#        self.OnMiddleButtonUp()
        key = self.RWI.GetKeySym()
#        print(key)
        if (self.mbpe) : 
            self.OnMiddleButtonUp()
            self.mbpe = False
        else : self.OnLeftButtonUp()
        return
        
    def mouseMoveEvent(self, obj, event) :
        #print("MouseMove")
        # xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
        # #print ("xu, yu, zu ", xu, yu, zu)
        
        # xn, yn, zn = self.Display.FindViewNormal()
        # xp, yp, zp = self.Display.renderer.GetActiveCamera().GetPosition()
        # xf, yf, zf = self.Display.renderer.GetActiveCamera().GetFocalPoint()
        # xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
        # xd, yd, zd = self.Display.renderer.GetActiveCamera().GetDirectionOfProjection()
        
        # xn1 = xp - xf
        # yn1 = yp - yf
        # zn1 = zp - zf
        # length_cam = np.sqrt(xn1**2 + yn1**2 + zn1**2)

        # xn1 = xn1 / length_cam
        # yn1 = yn1 / length_cam
        # zn1 = zn1 / length_cam
        # print("xn1, yn1, zn1 ", xn1, yn1, zn1)

        # print("xn, yn, zn = ", xn, yn, zn)
        # print("xp, yp, zp = ", xp, yp, zp)
        # print("xf, yf, zf = ", xf, yf, zf)
        # print("xu, yu, zu = ", xu, yu, zu)
        # print("xd, yd, zd", xd, yd, zd)
        # print("length of camera to focal point ", length_cam)
           
           
        self.OnMouseMove()
        
        
        
     
    def AlignPointsCylindrical(self) :
        
        self.pcld.clearPoints()
        self.pcld.addPoints()                        
        littled = self.pcld.IndexComb.sum()
        self.pcld.iproj_l = np.empty(littled, dtype = np.float64)
        self.pcld.jproj_l = np.empty(littled, dtype = np.float64)
        self.pcld.kproj_l = np.empty(littled, dtype = np.float64)

        tdens = PCF.TimeStamp(0, "")

        
        self.pcld.clearPoints()
        self.pcld.addPoints()
        self.Display.renderWindow.Render() 

        xn, yn, zn = self.Display.FindViewNormal()
        
        dotn = xn*xn + yn*yn +zn*zn

        xn = xn / dotn
        yn = yn / dotn
        zn = zn / dotn

        #self.Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)
        xc, yc, zc = self.Display.renderer.GetActiveCamera().GetFocalPoint()
        
        length_cam = self.Display.length_cam 
        
        offset_i_max = 0
        offset_j_max = 0
        offset_max_val = 0.00
        
        outfile = open("sweep_data_003.xyz", "w")
        
        range_max = 10
        #scan_inc = 0.001
        scan_inc = 0.005
        xn3 = xn
        yn3 = yn
        zn3 = zn
        toggle_i_j = 1

        for levels in range(0,4) :
            offset_i = 0.00
            offset_j = 0.00
            for i in range(range_max) :
                if toggle_i_j > 0 : 
                    offset_i = (i -range_max / 2)*scan_inc
                else :
                    offset_j = (i -range_max / 2)*scan_inc
                
                print("offset_i, j", offset_i, offset_j)
                print("xn, yn, zn", xn, yn, zn)
                if (self.pcld.dominant_axis == 0) : 
                    yn3 = yn + offset_i
                    zn3 = zn + offset_j
                if (self.pcld.dominant_axis == 1) : 
                    xn3 = xn + offset_i
                    zn3 = zn + offset_j
                if (self.pcld.dominant_axis == 2) : 
                    xn3 = xn + offset_i
                    yn3 = yn + offset_j

                dotn = xn3*xn3 + yn3*yn3 +zn3*zn3

                xn3 = xn3 / dotn
                yn3 = yn3 / dotn
                zn3 = zn3 / dotn                                
                
#                        tden = PCF.TimeStamp(tden,"1 set up")
                density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn3, yn3, zn3, length_cam)
#                        tden = PCF.TimeStamp(tden,"2 after Calculate Density")
            
                print("range_max", range_max)
#                print("i ", i, "offset_i = ", offset_i, "offset_j", offset_j, " xn + offset ", xn + offset, " Density sum = ", int(density_test))
                self.pcld.density_test = density_test
                self.pcld.xn = xn3
                self.pcld.yn = yn3
                self.pcld.zn = zn3
#                        tden = PCF.TimeStamp(tden,"3")
                self.Display.calculate()
                self.Display.root.update()
                self.Display.renderWindow.Render()
#                        tden = PCF.TimeStamp(tden,"4 after update, render")
                if (density_test > offset_max_val) : 
                    offset_i_max = offset_i
                    offset_j_max = offset_j
                    offset_max_val = density_test
#                        tden_TOT = PCF.TimeStamp(tden0,"5 TOTAL")
                a_str = '%.4f' % offset_i
                b_str = '%.4f' % offset_j
                c_str = '%.4f' % xn3
                d_str = '%.4f' % yn3
                e_str = '%.4f' % zn3
                f_str = '%.4f' % self.pcld.density_test
                LineToWrite = a_str + '\t' + b_str + '\t' + c_str + '\t' + d_str + \
                    '\t' + e_str + '\t' + f_str +'\n'
                outfile.write(LineToWrite)
                print("offset_i, offset_j, xn3, yn3, zn3", LineToWrite)
                               
#                    print("Final offset_max ", offset_max)

            if (self.pcld.dominant_axis == 0) : 
                yn = yn + offset_i_max
                zn = zn + offset_j_max
            if (self.pcld.dominant_axis == 1) : 
                xn = xn + offset_i_max
                zn = zn + offset_j_max
            if (self.pcld.dominant_axis == 2) : 
                xn = xn + offset_i_max
                yn = yn + offset_j_max
                
            a_str = '%.4f' % offset_i_max
            b_str = '%.4f' % offset_j_max
            c_str = '%.4f' % offset_max_val
            d_str = '%.4f' % xn
            e_str = '%.4f' % yn
            f_str = '%.4f' % zn
            LineToWrite = a_str + '\t' + b_str + '\t' + c_str + '\t' + d_str + \
                '\t' + e_str + '\t' + f_str +'\n'
            outfile.write(LineToWrite)
            print("Maximum")
            print(LineToWrite)
            
            toggle_i_j = toggle_i_j * -1
            if toggle_i_j > 0 :
                scan_inc = scan_inc / 10   
                        
            offset_max_val_final = offset_max_val
            offset_max_val = 0
            offset_i_max = 0
            offset_j_max = 0
            
                
            dotn = xn*xn + yn*yn +zn*zn
            xn = xn / dotn
            yn = yn / dotn
            zn = zn / dotn
                   
            xn3 = xn
            yn3 = yn 
            zn3 = zn
            self.pcld.xn = xn3 
            self.pcld.yn = yn3
            self.pcld.zn = zn3
            density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn3, yn3, zn3, length_cam)                    
            self.Display.calculate()
            self.Display.root.update()
            self.Display.renderWindow.Render() 
                           
            xcam1 = xn * length_cam + xc
            ycam1 = yn * length_cam + yc 
            zcam1 = zn * length_cam + zc
            self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1)   

#                  density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn3, yn3, zn3, length_cam)
            self.pcld.density_test = offset_max_val_final
            self.Display.calculate()
            self.Display.root.update()
            self.Display.renderWindow.Render()                    


        outfile.close()
        
        # dotn = xn*xn + yn*yn +zn*zn
        # xn = xn / dotn
        # yn = yn / dotn
        # zn = zn / dotn
        xn, yn, zn = FT.unitize(xn, yn, zn)

        density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn, yn, zn, length_cam)
        

        xcam1 = xn * length_cam + xc
        ycam1 = yn * length_cam + yc 
        zcam1 = zn * length_cam + zc
        self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1)   
        self.Display.renderWindow.Render() 
        

        
        xn, yn, zn = self.Display.FindViewNormal()
        self.pcld.xn = xn
        self.pcld.yn = yn
        self.pcld.zn = zn
        print("recalculated:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
            
        self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
        self.pcld.clearPoints()
        self.pcld.addPoints()

            
     
    def AlignPointsPlanar(self) :
        
#        print("1")
        self.pcld.clearPoints()
        self.pcld.addPoints()  
#        print("2")
        littled = self.pcld.IndexComb.sum()
        self.pcld.iproj_l = np.empty(littled, dtype = np.float64)
        self.pcld.jproj_l = np.empty(littled, dtype = np.float64)
        self.pcld.kproj_l = np.empty(littled, dtype = np.float64)

#        print("3")
        self.pcld.clearPoints()
        self.pcld.addPoints()
        self.Display.renderWindow.Render() 

        xn, yn, zn = self.Display.FindViewNormal()
        xf, yf, zf = self.Display.renderer.GetActiveCamera().GetFocalPoint()
        dotn = xn*xn + yn*yn +zn*zn
     
        xn = xn / dotn
        yn = yn / dotn
        zn = zn / dotn
        
        dominant_axis_n = 0
        if (abs(zn) > abs(xn) and abs(zn) > abs(yn)) : dominant_axis_n = 2
        if (abs(yn) > abs(xn) and abs(yn) > abs(zn)) : dominant_axis_n = 1
        print("4")
        xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
        dotu = xu*xu + yu*yu +zu*zu
     
        xu = xu / dotu
        yu = yu / dotu
        zu = zu / dotu
        
        print("xu, yu, zu ", xu, yu, zu)
        self.pcld.dominant_axis_u = 0
        if (abs(zu) > abs(xu) and abs(zu) > abs(yu)) : self.pcld.dominant_axis_u = 2
        if (abs(yu) > abs(xu) and abs(yu) > abs(zu)) : self.pcld.dominant_axis_u = 1        
        
        print("self.pcld.dominant_axis_u ", self.pcld.dominant_axis_u)
        print("5")

        xc, yc, zc = self.Display.renderer.GetActiveCamera().GetFocalPoint()

        
        length_cam = self.Display.length_cam 

#        self.Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)
        
        offset_i_max = 0
        offset_j_max = 0
        offset_max_val = 0.00
          
        range_max = 10
        scan_inc = 0.5
        xn3 = xn
        yn3 = yn
        zn3 = zn
        toggle_i_j = 1
        
        print("6")

        for levels in range(0,2) :
            self.Display.renderer.GetActiveCamera().Elevation(-1 * scan_inc * range_max/2)

            print("7: loop")
            offset = 0.0
            offset_max_val = 0.0
            for i in range(range_max) :
                print("8: loop")
                self.Display.renderer.GetActiveCamera().Elevation(scan_inc)
                print("8a: loop")                
                offset = offset + scan_inc
                xn, yn, zn = self.Display.FindViewNormal()
                print("8b: loop")
                density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)
#                CalculateDensityProjection2(self, Display, xc, yc, zc, xn, yn, zn, xu, yu, zu, length_cam) :
                print("8c: loop")
                self.pcld.density_test = density_test
                self.pcld.xn = xn
                self.pcld.yn = yn
                self.pcld.zn = zn
                print("8d: loop")
                self.Display.calculate()
                self.Display.root.update()
                self.Display.renderWindow.Render()
                print("8e: loop")

                if (density_test > offset_max_val) : 
                    offset_max = offset
                    offset_max_val = density_test

            scan_inc = scan_inc / 10   
                        
            self.Display.renderer.GetActiveCamera().Elevation(offset_max - offset)  
            xn, yn, zn = self.Display.FindViewNormal()            
            offset_max_val_final = offset_max_val
            offset_max_val = 0
            offset_max = 0


            self.pcld.xn = xn
            self.pcld.yn = yn
            self.pcld.zn = zn
            
            xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
            dotu = xu*xu + yu*yu +zu*zu
     
            xu = xu / dotu
            yu = yu / dotu
            zu = zu / dotu
            
        
            #density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, length_cam)     
            density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)            
            self.Display.calculate()
            self.Display.root.update()
            self.Display.renderWindow.Render() 
                           
            xcam1 = xn * length_cam + xf
            ycam1 = yn * length_cam + yf 
            zcam1 = zn * length_cam + zf
            self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1)   

            self.pcld.density_test = offset_max_val_final
            self.Display.calculate()
            self.Display.root.update()
            self.Display.renderWindow.Render()                    


            
        xn, yn, zn = self.Display.FindViewNormal()            

        self.pcld.xn = xn
        self.pcld.yn = yn
        self.pcld.zn = zn
        
        xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
        dotu = xu*xu + yu*yu +zu*zu
 
        xu = xu / dotu
        yu = yu / dotu
        zu = zu / dotu
        
    
        #density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, length_cam)     
        density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)  
            
        #density_test = self.pcld.CalculateDensityProjection(self.Display, xf, yf, zf, xn, yn, zn, length_cam)
        
        dotn = xn*xn + yn*yn +zn*zn
        xn = xn / dotn
        yn = yn / dotn
        zn = zn / dotn
            
        xcam1 = xn * length_cam + xc
        ycam1 = yn * length_cam + yc 
        zcam1 = zn * length_cam + zc
        self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1)   
        
            
        self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
        self.pcld.clearPoints()
        self.pcld.addPoints()

            

            
        
    def keyPress(self, obj,event):
        #global XcontextView
        #global YcontextView
        #global ZcontextView
        global Selection
        #key = renderWindowInteractor.GetKeySym()
        #code = renderWindowInteractor.GetKeyCode()
        key = self.RWI.GetKeySym()
        code = self.RWI.GetKeyCode()
        #print(key)
        if (key == 'Control_L') : 
#            self.RWI.SetKeySym('Shift')
            return()
        
        xcam = 0.00
        ycam = 0.00
        zcam = 0.00
        
        xcam,ycam,zcam = self.Display.renderer.GetActiveCamera().GetPosition()
        self.xcam = xcam
        self.ycam = ycam
        self.zcam = zcam
        
        print("code = " , code)
        print("key = ", key)
        
        ts = 0
        tl = PCF.TimeStamp(ts, 'Start')
#        if (key == 'Return') :
#            if self.pcld.statenames[self.pcld.state + 1] != "" :
#                self.pcld.state += 1
#                print("State: ", self.pcld.state)


       #print("State: ", self.pcld.state, '\t', self.pcld.statenames[self.pcld.state])
        if (key =='Escape') : exit()


        #self.statenames[0] = "Manually orient component"
        #self.statenames[1] = "Highlight main plane"
        #self.statenames[2] = "Analytically determine main plane and reorient"
        #self.statenmaes[3] = "Align center of origin dowel to origin"
        #self.statenames[4] = "Rotate part in XY plane to align second dowel"
        
        # xn, yn, zn = self.Display.renderer.GetActiveCamera().GetViewPlaneNormal()
        # #print("Plane Normal obtained, ", xn, yn, zn)
        # self.pcld.xn = xn
        # self.pcld.yn = yn
        # self.pcld.zn = zn
        
        xn, yn, zn = self.Display.FindViewNormal()
        print("keypress:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))


        if (self.pcld.dominant_axis == 0) : self.pcld.Sph_Xpos = 0
        if (self.pcld.dominant_axis == 1) : self.pcld.Sph_Ypos = 0
        if (self.pcld.dominant_axis == 2) : self.pcld.Sph_Zpos = 0

        if (np.pi < 3) :
            if (abs(xn) >  max(abs(yn),abs(zn))) : 
                if (self.pcld.dominant_axis == 1) :
                    H_pos = self.pcld.Sph_Xpos
                    V_pos = self.pcld.Sph_Zpos
                    self.pcld.Sph_Xpos = 0.00
                    self.pcld.Sph_Ypos = H_pos
                    self.pcld.Sph_Zpos = V_pos

                if (self.pcld.dominant_axis == 2) :                
                    H_pos = self.pcld.Sph_Xpos
                    V_pos = self.pcld.Sph_Ypos
                    self.pcld.Sph_Xpos = 0.00
                    self.pcld.Sph_Ypos = H_pos
                    self.pcld.Sph_Zpos = V_pos
                    
                self.pcld.dominant_axis = 0
               #print("mostly oriented in the x direction")                
            else :
                if (abs(yn) > max(abs(xn),abs(zn))) :
                    if (self.pcld.dominant_axis == 0) :
                        H_pos = self.pcld.Sph_Ypos
                        V_pos = self.pcld.Sph_Zpos
                        self.pcld.Sph_Xpos = H_pos
                        self.pcld.Sph_Ypos = 0.00
                        self.pcld.Sph_Zpos = V_pos

                    if (self.pcld.dominant_axis == 2) :                
                        H_pos = self.pcld.Sph_Xpos
                        V_pos = self.pcld.Sph_Ypos
                        self.pcld.Sph_Xpos = H_pos
                        self.pcld.Sph_Ypos = 0.00
                        self.pcld.Sph_Zpos = V_pos
                        
                    self.pcld.dominant_axis = 1
                   #print("mostly oriented in the y direction)")

                else :
                    if (self.pcld.dominant_axis == 0) :
                        H_pos = self.pcld.Sph_Ypos
                        V_pos = self.pcld.Sph_Zpos
                        self.pcld.Sph_Xpos = H_pos
                        self.pcld.Sph_Ypos = 0.00
                        self.pcld.Sph_Zpos = V_pos

                    if (self.pcld.dominant_axis == 1) :
                        H_pos = self.pcld.Sph_Xpos
                        V_pos = self.pcld.Sph_Zpos
                        self.pcld.Sph_Xpos = 0.00
                        self.pcld.Sph_Ypos = H_pos
                        self.pcld.Sph_Zpos = V_pos
                        
                    self.pcld.dominant_axis = 2
                   #print("mostly oriented in the z direction")
                        
        self.pcld.xn = xn
        self.pcld.yn = yn
        self.pcld.zn = zn
        xc = self.pcld.Sph_Xpos
        yc = self.pcld.Sph_Ypos
        zc = self.pcld.Sph_Zpos
        
        xp, yp, zp = self.Display.renderer.GetActiveCamera().GetPosition()
        xf, yf, zf = self.Display.renderer.GetActiveCamera().GetFocalPoint()
        xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
        xd, yd, zd = self.Display.renderer.GetActiveCamera().GetDirectionOfProjection()
        
        xn1 = xcam - xf
        yn1 = ycam - yf
        zn1 = zcam - zf
        length_cam = np.sqrt(xn1**2 + yn1**2 + zn1**2)
        self.Display.length_cam = length_cam
        
        if (key == 'asterisk') or (key == '8'): 
            #xn, yn, zn = self.Display.renderer.GetActiveCamera().GetViewPlaneNormal()

            

            xn1 = xn1 / length_cam
            yn1 = yn1 / length_cam
            zn1 = zn1 / length_cam

           
            sign_xn = 1.00
            sign_yn = 1.00
            sign_zn = 1.00
            if xn1 < 0.00 : sign_xn = -1.00
            if yn1 < 0.00 : sign_yn = -1.00
            if zn1 < 0.00 : sign_zn = -1.00
            
            sign_xu = 1.00
            sign_yu = 1.00
            sign_zu = 1.00
            if xn1 < 0.00 : sign_xu = -1.00
            if yn1 < 0.00 : sign_yu = -1.00
            if zn1 < 0.00 : sign_zu = -1.00
            
            #print("signs ", sign_x, sign_y, sign_z)
           #print("xn1, yn1, zn1 ", xn1, yn1, zn1)
           #print("xcam, ycam, zcam = ", xcam, ycam, zcam)
           #print("xn, yn, zn = ", xn, yn, zn)
           #print("xp, yp, zp = ", xp, yp, zp)
           #print("xf, yf, zf = ", xf, yf, zf)
           #print("xu, yu, zu = ", xu, yu, zu)
           #print("xc, yc, zc", xc, yc, zc)
           #print("length of camera to focal point ", length_cam)
            xn1 = 0.00
            yn1 = 0.00
            zn1 = 0.00
            self.Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)
           
           #print("xn ", xn, "abs(xn) ", abs(xn), "abs(yn)", abs(yn), " abs(zn)", abs(zn), "max(abs(yn),abs(zn))", max(abs(yn),abs(zn)))
            if (abs(xn) >  max(abs(yn),abs(zn))) : 
                xn1 = 1.00
               #print("mostly oriented in the x direction")
                if (key == 'asterisk') or (self.lastkey == '8') : 
                    self.Display.renderer.GetActiveCamera().SetViewUp(0,0,1.0)  

                #self.Display.renderer.GetActiveCamera().GetViewPlaneNormal(0,1.00,0)                
            else :
                if (abs(yn) > max(abs(xn),abs(zn))) :
                    yn1 = -1.00
                   #print("mostly oriented in the y direction)")
                    if (key == 'asterisk') or (self.lastkey == '8') : 
                        self.Display.renderer.GetActiveCamera().SetViewUp(0,0,1)  
                    #self.Display.renderer.GetActiveCamera().GetViewPlaneNormal(0,0,-1.0) 
                else :
                    #if (zn > 0.9) :
                    zn1 = 1.00
                   #print("mostly oriented in the z direction")
                    if (key == 'asterisk') or (self.lastkey == '8') : 
                        self.Display.renderer.GetActiveCamera().SetViewUp(0,1,0)  
            
            if (key == '8') : #Simply duplicate the general orientation that the part had been in.
                xn1 = sign_xn * abs(xn1)
                yn1 = sign_yn * abs(yn1)
                zn1 = sign_zn * abs(zn1)
            
            xcam1 = xn1 * length_cam + xc
            ycam1 = yn1 * length_cam + yc 
            zcam1 = zn1 * length_cam + zc
           #print("xcam1, ycam1, zcam1", xcam1, ycam1, zcam1)
            self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1) 
            self.Display.renderer.GetActiveCamera().OrthogonalizeViewUp()

        if (self.pcld.state == 6) :   #Extract plane geometry
        
            # xc = self.pcld.X.sum() / self.pcld.BigD
            # yc = self.pcld.Y.sum() / self.pcld.BigD
            # zc = self.pcld.Z.sum() / self.pcld.BigD
            # self.Display.renderer.GetActiveCamera().SetFocalPoint(xc,yc,zc)
            
            self.pcld.slice_dist = (self.pcld.X) * xu + \
                (self.pcld.Y) * yu + \
                (self.pcld.Z) * zu
                
            if (key == 'U' ) :
                self.pcld.slice_max = self.pcld.slice_max + self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.HighlightPlane(self.Display)
            if (key == 'u' ) :
                self.pcld.slice_max = self.pcld.slice_max - self.pcld.inc
                print("self.pcld.slice_max", self.pcld.slice_max, self.pcld.slice_dist.max(), self.pcld.slice_dist.min(), \
                    len(self.pcld.slice_dist), self.pcld.IndexSlice.sum())
#                self.pcld.IndexSlice = np.where(np.less(self.pcld.slice_dist,self.pcld.slice_max), True,False)
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                print("self.pcld.slice_max", self.pcld.slice_max, self.pcld.slice_dist.max(), self.pcld.IndexSlice.sum())
                self.pcld.HighlightPlane(self.Display)                
            if (key == 'L' ) :
                self.pcld.slice_min = self.pcld.slice_min + self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.HighlightPlane(self.Display)                    
            if (key == 'l' ) :
                self.pcld.slice_min = self.pcld.slice_min - self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.HighlightPlane(self.Display)  
                
            if (key == 'a') or (key == 'A') :
                LittleD_B = self.pcld.IndexComb.sum()
                print("LittleD_B", LittleD_B)
                if (LittleD_B > 2e4) : 
                    print("Too many points highlighted")
                else :
                    print("Realign allowed")
                    reorient = True
                    #if (self.pcld.DensityOld) :

                    print("a:  Highlight Plane")
                    self.pcld.HighlightPlane(self.Display)  
                    print("a:  AlignPointsPlanar coming up")
                    self.AlignPointsPlanar()
                    self.pcld.RefinedState = 1
                    
            if (key == 'p') or (key == 'P') :

            
                # Now that we have a set of points projected onto a plane, we can find the line.
                xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
                density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)
                #density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)
                if (len(self.pcld.iproj_l) > 10) : 

                    
      
#                    self.pcld.i_ctr, self.pcld.j_ctr, R = lsfc.find_circle(self.pcld.iproj_l, self.pcld.jproj_l)
                    slope, intercept = PCF.CreateSlope(self.pcld.iproj_l, self.pcld.jproj_l, len(self.pcld.iproj_l))
                    print("Roll by this amount", -1*180/np.pi*np.arctan(slope))
                    self.Display.renderer.GetActiveCamera().Roll(-180/np.pi*np.arctan(slope))

                    xu, yu, zu = self.Display.renderer.GetActiveCamera().GetViewUp()
                    print("xu, yu, zu ", xu, yu, zu)
                    xf, yf, zf = self.Display.renderer.GetActiveCamera().GetFocalPoint()
                    xn, yn, zn = self.Display.FindViewNormal()
                    density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)
                    slope, intercept = PCF.CreateSlope(self.pcld.iproj_l, self.pcld.jproj_l, len(self.pcld.iproj_l))
                    print("Roll by this amount", -1*180/np.pi*np.arctan(slope))
                    self.Display.renderer.GetActiveCamera().Roll(-180/np.pi*np.arctan(slope))
                    #xp = xf + xu * intercept
                    #yp = yf + yu * intercept
                    #zp = zf + zu * intercept
                    print("xf, yf, zf ", xf, yf, zf)
                    print("xp, yp, zp ", xp, yp, zp)
                    print("xu, yu, zu ", xu, yu, zu)
                    #print("Based on the i-intercept, the plane includes this point", xp, yp, zp)
                    #D = (xp * xu + yp * yu + zp * zu)
                    self.pcld.A = xu
                    self.pcld.B = yu
                    self.pcld.C = zu
                    self.pcld.D = intercept
                    print("up vector ", xu, yu, zu)
                    print("Equation of the plane = ", self.pcld.A, "x ", self.pcld.B, "y ", \
                        self.pcld.C, "z + ", self.pcld.D)
                 
            if (key == 'space') :

                i =self.pcld.feature_index
                self.pcld.feature_index += 1
                self.pcld.Features[i].FeatureName = "Plane"
                self.pcld.Features[i].Description = dt.start()
                self.pcld.Features[i].a = self.pcld.A
                self.pcld.Features[i].b = self.pcld.B
                self.pcld.Features[i].c = self.pcld.C
                self.pcld.Features[i].d = self.pcld.D
                self.pcld.Features[i].e = 0.00
                self.pcld.Features[i].f = 0.00
                
            if (key == 'g' or key == 'G') :
                density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)
                self.pcld.Group2()
                #self.pcld.SaveGroupResults()
                
            if (key == 'Up') :
                self.pcld.Rmin = self.pcld.Rmin**0.5
                print("Min Rsquared value = ", self.pcld.Rmin)
                self.pcld.Undo()
                self.pcld.Group2()
            if (key == 'Down' ) : 
                self.pcld.Rmin = self.pcld.Rmin**2
                print("Min Rsquared value = ", self.pcld.Rmin)
                self.pcld.Undo()
                self.pcld.Group2()
            if (key == "Left") : 
                self.pcld.G0max = self.pcld.G0max * 0.9
                self.pcld.G1max = self.pcld.G1max * 0.9
                print("Max slope ", self.pcld.G0max)
                self.pcld.Undo()
                self.pcld.Group2()
            if (key == "Right") :
                self.pcld.G0max = self.pcld.G0max * 1.1
                self.pcld.G1max = self.pcld.G1max * 1.1
                print("Max slope ", self.pcld.G0max)
                self.pcld.Undo()
                self.pcld.Group2()
                
                
            
            if (key == 'Return') and (self.lastkey != 'Return') :
                density_test = self.pcld.CalculateDensityProjection2(self.Display, xf, yf, zf, xn, yn, zn, xu, yu, zu, length_cam)
                self.pcld.Group2()
                print("Hit Return, Find Longest String2")
                self.pcld.FindLongestString2()
                #self.pcld.FindPlane2() 
                #self.pcld.DensityOld = True                    
                #self.pcld.AdvanceState()
                
                
        if (self.pcld.state == 5) :  #test state
            self.allow_dominant_axis_change = True

            # xc = self.pcld.Sph_Xpos
            # yc = self.pcld.Sph_Ypos
            # zc = self.pcld.Sph_Zpos
            #
            xcam1 = xn * length_cam + xc
            ycam1 = yn * length_cam + yc 
            zcam1 = zn * length_cam + zc
           #print("xcam1, ycam1, zcam1", xcam1, ycam1, zcam1)
            #
            
            reorient = False
            # if (key == 'x') or (key == 'X'):  xcam = xcam + 100
            # if (key == 'y') or (key == 'Y'):  ycam = ycam + 100
            # if (key == 'z') or (key == 'Z'):  zcam = zcam + 100
            

                
            #if (key == '7') : 
                
            
            self.pcld.slice_dist = (self.pcld.X - xc) * xn + \
                (self.pcld.Y - yc ) * yn + \
                (self.pcld.Z - zc) * zn
                                 
            if (key == 'U' ) :
                self.pcld.slice_max = self.pcld.slice_max + self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.FindCylinder(self.Display)
            if (key == 'u' ) :
                self.pcld.slice_max = self.pcld.slice_max - self.pcld.inc
                print("self.pcld.slice_max", self.pcld.slice_max, self.pcld.slice_dist.max(), self.pcld.slice_dist.min(), \
                    len(self.pcld.slice_dist), self.pcld.IndexSlice.sum())
#                self.pcld.IndexSlice = np.where(np.less(self.pcld.slice_dist,self.pcld.slice_max), True,False)
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                print("self.pcld.slice_max", self.pcld.slice_max, self.pcld.slice_dist.max(), self.pcld.IndexSlice.sum())
                self.pcld.FindCylinder(self.Display)                
            if (key == 'L' ) :
                self.pcld.slice_min = self.pcld.slice_min + self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.FindCylinder(self.Display)                    
            if (key == 'l' ) :
                self.pcld.slice_min = self.pcld.slice_min - self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.FindCylinder(self.Display)                    
            
            
            if (key == 'h') :
                print("h was pressed")
                self.pcld.radius = self.pcld.radius * 0.9
                self.pcld.DensityOld = True 
                self.pcld.FindCylinder(self.Display)
                reorient = True
            if (key == 'H') :
                print("H was pressed")
                self.pcld.radius = self.pcld.radius * 1.2
                self.pcld.DensityOld = True 
                self.pcld.FindCylinder(self.Display)
                reorient = True
                
            if (key == 'd') or (key == 'D') :
                self.Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)
                density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn, yn, zn, length_cam)
                self.pcld.density_test = density_test
                self.pcld.xn = xn
                self.pcld.yn = yn
                self.pcld.zn = zn
                self.Display.calculate()
                self.Display.root.update()
                        
            
            if (key == 'a') or (key == 'A') :
  
                reorient = True
                #if (self.pcld.DensityOld) :
                self.pcld.FindCylinder(self.Display)
                self.AlignPointsCylindrical()
                
                xn, yn, zn = self.Display.FindViewNormal()
                self.pcld.xn = xn
                self.pcld.yn = yn
                self.pcld.zn = zn
                print("After AlignPointsCylindrical:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
                
            xn, yn, zn = self.Display.FindViewNormal()
            print("After exiting 'a':  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
                
            if (key == 'c') or (key == 'C') :
                        
                # Now that we have a set of points projected onto a plane, we can find the center.
                density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn, yn, zn, length_cam)
                if (len(self.pcld.iproj_l) > 10) : 

                    
      
                    self.pcld.i_ctr, self.pcld.j_ctr, R = lsfc.find_circle(self.pcld.iproj_l, self.pcld.jproj_l)
                    #self.pcld.X_ctr, self.pcld.Y_ctr, R = lsfc.find_circle(self.pcld.iproj_l, self.pcld.jproj_l)
                    
                    xc2 = self.pcld.i_ctr * self.pcld.xi + self.pcld.j_ctr * self.pcld.xj
                    yc2 = self.pcld.i_ctr * self.pcld.yi + self.pcld.j_ctr * self.pcld.yj
                    zc2 = self.pcld.i_ctr * self.pcld.zi + self.pcld.j_ctr * self.pcld.zj
                    
                    print("self.pcld.i_ctr, j_ctr", self.pcld.i_ctr, self.pcld.j_ctr)
                    print("Recenter, xc2, yc2, zc2", xc2, yc2, zc2)
                    
                    print("self.pcld.xi, yi, zi", self.pcld.xi, self.pcld.yi, self.pcld.zi)
                    print("self.pcld.xj, yj, zj", self.pcld.xj, self.pcld.yj, self.pcld.zj)
                    print("xc2, yc2, zc2", xc2, yc2, zc2)
                    
                    print("dominant aixs", self.pcld.dominant_axis)
                    #  This needs to be adjusted to allow resolving other major planes.

                    if self.pcld.dominant_axis == 0 :
                        m = xc2 / xn
                        xt = 0.00
                        yt = yc2 - m * yn
                        zt = zc2 - m * zn  
                        
                    if self.pcld.dominant_axis == 1 :
                        m = yc2 / yn
                        yt = 0.00
                        xt = xc2 - m * xn
                        zt = zc2 - m * zn                       
                        
                    if self.pcld.dominant_axis == 2 :
                        m = zc2 / zn
                        zt = 0.00
                        xt = xc2 - m * xn                        
                        yt = yc2 - m * yn

                    print("xt, yt, zt", xt, yt, zt)
                    print("xc, yc, zc", xc, yc, zc)
                
                    xt = xt + xc
                    yt = yt + yc
                    zt = zt + zc
                    self.pcld.Sph_Xpos = xt
                    self.pcld.Sph_Ypos = yt
                    self.pcld.Sph_Zpos = zt
                    self.sphere.SetCenter(self.pcld.Sph_Xpos, self.pcld.Sph_Ypos, self.pcld.Sph_Zpos)



                    

            if (key == 'space') :

                i =self.pcld.feature_index
                self.pcld.feature_index += 1
                self.pcld.Features[i].FeatureName = "Cylinder"
                self.pcld.Features[i].Description = dt.start()
                self.pcld.Features[i].a = self.pcld.Sph_Xpos
                self.pcld.Features[i].b = self.pcld.Sph_Ypos
                self.pcld.Features[i].c = self.pcld.Sph_Zpos
                self.pcld.Features[i].d = xn
                self.pcld.Features[i].e = yn
                self.pcld.Features[i].f = zn
                    
                if (self.pcld.dominant_axis == 0) :
                    if (yn != 0) : 
                        Angle_XY = 180/np.pi*np.arctan(xn/yn)
                    else :
                        Angle_XY = 90                   
                    if (zn != 0) : 
                        Angle_XZ = 180/np.pi*np.arctan(xn/zn)
                    else :
                        Angle_XZ = 90
                        
                if (self.pcld.dominant_axis == 1) :
                    if (xn != 0) : 
                        Angle_XY = 180/np.pi*np.arctan(yn/xn)
                    else :
                        Angle_XY = 90                   
                    if (zn != 0) : 
                        Angle_YZ = 180/np.pi*np.arctan(yn/zn)
                    else :
                        Angle_YZ = 90
                                      
                    print("Angle of cylinder in XY plane", Angle_XY)
                    print("Angle of cylinder in YZ plane", Angle_YZ)
                    
                if (self.pcld.dominant_axis == 2) :
                    if (yn != 0) : 
                        Angle_YZ = 180/np.pi*np.arctan(zn/yn)
                    else :
                        Angle_YZ = 90                   
                    if (xn != 0) : 
                        Angle_XZ = 180/np.pi*np.arctan(zn/xn)
                    else :
                        Angle_XZ = 90
                                      
                    print("Angle of cylinder in YZ plane", Angle_YZ)
                    print("Angle of cylinder in XZ plane", Angle_XZ)
                #print("Cylinder axis crosses xy plane at x, y = ", self.pcld.X_ctr, self.pcld.Y_ctr)        
        
                    

            if (reorient) :
                print("reorient")
                xcam1 = xn * length_cam + xc
                ycam1 = yn * length_cam + yc 
                zcam1 = zn * length_cam + zc
            
                self.Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)
                self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1) 

        if (self.pcld.state == 5) or (self.pcld.state == 6) :
            if (key == 'greater') :
                self.pcld.Dens_min = self.pcld.Dens_min * 0.4
                print("New Density_min = ", self.pcld.Dens_min)
                self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
            if (key == 'less') :
                self.pcld.Dens_min = self.pcld.Dens_min * 2
                print("New Density_min = ", self.pcld.Dens_min)
                self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
        
        xn, yn, zn = self.Display.FindViewNormal()
        print("After last state 5 calcs:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
        
        if (self.pcld.state == 4) :  #Align rotation of part
            self.allow_dominant_axis_change = False
        #self.sphere.SetCenter(self.pcld.X_ctr,self.pcld.Y_ctr,0)
            # if (key == 'space') :
                # print("Space was pressed")
                # self.pcld.FindCircle()
                # self.sphere.SetCenter(self.pcld.X_ctr,self.pcld.Y_ctr,0)
            if (key == 'Return') and (self.lastkey != 'Return') :
                #self.sphere.SetCenter(0,0,0)
                Angle = 180 / np.pi * np.arctan(self.pcld.Y_ctr / self.pcld.X_ctr)
                print("Rotate part in the XY plane by angle ", Angle)
                self.pcld.AxisRotate(2 , Angle)
                #self.pcld.AxisRotate(2 , 2)
                self.pcld.AdvanceState()
                
        if (self.pcld.state == 3) : #Align center of bolt hole with origin
            self.allow_dominant_axis_change = False
            if (key == 'Return') and (self.lastkey != 'Return') :
                self.sphere.SetCenter(0,0,0)
                self.pcld.Sph_Xpos = 0
                self.pcld.Sph_Ypos = 0
                self.pcld.Sph_Zpos = 0
                self.pcld.X = self.pcld.X - self.pcld.X_ctr
                self.pcld.Y = self.pcld.Y - self.pcld.Y_ctr
                self.pcld.AdvanceState()
                
        if (self.pcld.state == 3) or (self.pcld.state == 4) :  
        
            if (self.pcld.dominant_axis == 0) : self.pcld.slice_dist = self.pcld.X - xc
            if (self.pcld.dominant_axis == 1) : self.pcld.slice_dist = self.pcld.Y - yc
            if (self.pcld.dominant_axis == 2) : self.pcld.slice_dist = self.pcld.Z - zc

            print("xn, yn, zn ", xn, yn, zn)
                
            
            if (key == 'U' ) :
                self.pcld.slice_max = self.pcld.slice_max + self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.FindCircle()
            if (key == 'u' ) :
                self.pcld.slice_max = self.pcld.slice_max - self.pcld.inc
                print("self.pcld.slice_max", self.pcld.slice_max, self.pcld.slice_dist.max(), self.pcld.slice_dist.min(), \
                    len(self.pcld.slice_dist), self.pcld.IndexSlice.sum())
#                self.pcld.IndexSlice = np.where(np.less(self.pcld.slice_dist,self.pcld.slice_max), True,False)
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                print("self.pcld.slice_max", self.pcld.slice_max, self.pcld.slice_dist.max(), self.pcld.IndexSlice.sum())
                self.pcld.FindCircle()              
            if (key == 'L' ) :
                self.pcld.slice_min = self.pcld.slice_min + self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.FindCircle()                    
            if (key == 'l' ) :
                self.pcld.slice_min = self.pcld.slice_min - self.pcld.inc
                self.pcld.IndexSlice = np.where(np.logical_and(np.less(self.pcld.slice_dist,self.pcld.slice_max), \
                    np.less(self.pcld.slice_min,self.pcld.slice_dist)), True,False)
                self.pcld.FindCircle() 
                
                
            if (key == 'space') :
                if (self.pcld.DensityOld) :
                    self.pcld.FindCircle()
                    self.pcld.clearPoints()
                    self.pcld.addPoints()
                    #self.pcld.X_ctr, self.pcld.Y_ctr, R, self.pcld.Index = PCF.EstimateCenterPoint(self.pcld.XdispB, self.pcld.YdispB, self.pcld.IndexComb, self.pcld.Dens_min)
                    self.pcld.density_big = PCF.EstimateDensity(self.pcld.XdispB, self.pcld.YdispB, self.pcld.IndexComb)
                    self.pcld.DensOld = False
                self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
                self.pcld.clearPoints()
                self.pcld.addPoints()
                if (len(self.pcld.XdispB) > 10) : 
                    self.pcld.X_ctr, self.pcld.Y_ctr, R = lsfc.find_circle(self.pcld.XdispB, self.pcld.YdispB)
                    self.pcld.Sph_Xpos = self.pcld.X_ctr
                    self.pcld.Sph_Ypos = self.pcld.Y_ctr
                    self.sphere.SetCenter(self.pcld.Sph_Xpos, self.pcld.Sph_Ypos, self.pcld.Sph_Zpos)
                    print("X_ctr, Y_ctr", self.pcld.X_ctr, self.pcld.Y_ctr) 
                else :
                    print("Insufficient data points, len(self.pcld.XdispB) = ", len(self.pcld.XdispB))
            if (key == 'greater') :
                self.pcld.Dens_min = self.pcld.Dens_min * 0.8
                print("New Density_min = ", self.pcld.Dens_min)
                self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
            if (key == 'less') :
                self.pcld.Dens_min = self.pcld.Dens_min * 2
                print("New Density_min = ", self.pcld.Dens_min)
                self.pcld.IndexArray = np.where(self.pcld.density_big > self.pcld.Dens_min, True, False)
            if (key == 'c') :
                print("c was pressed")
                self.pcld.radius = self.pcld.radius * 0.9
                self.pcld.DensityOld = True 
                self.pcld.FindCircle()
            if (key == 'C') :
                print("C was pressed")
                self.pcld.radius = self.pcld.radius * 1.2
                self.pcld.DensityOld = True 
                self.pcld.FindCircle()
            
        if (self.pcld.state == 3) or (self.pcld.state ==4) or (self.pcld.state == 5) :   
            V_off = 0.00
            H_off = 0.00
            return_com = 0
            
            cam_pos_adjust = False
            
            if (self.pcld.rot_or_tran == 0) :
                if (key == 'x') :  
                    xn = xn - self.pcld.inc/100
                    cam_pos_adjust = True
                if (key == 'X') :  
                    xn = xn + self.pcld.inc/100
                    cam_pos_adjust = True
                if (key == 'y') :  
                    yn = yn - self.pcld.inc/100
                    cam_pos_adjust = True
                if (key == 'Y') :  
                    yn = yn + self.pcld.inc/100
                    cam_pos_adjust = True
                if (key == 'z') :  
                    zn = zn - self.pcld.inc/100
                    cam_pos_adjust = True
                if (key == 'Z') :  
                    zn = zn + self.pcld.inc/100
                    cam_pos_adjust = True
            
            if (cam_pos_adjust) :
                self.Display.renderer.GetActiveCamera().SetFocalPoint(xc, yc, zc)

                
                dotn = xn*xn + yn*yn +zn*zn
                xn = xn / dotn
                yn = yn / dotn
                zn = zn / dotn
                    
                xcam1 = xn * length_cam + xc
                ycam1 = yn * length_cam + yc 
                zcam1 = zn * length_cam + zc
                self.Display.renderer.GetActiveCamera().SetPosition(xcam1, ycam1, zcam1) 
                
                density_test = self.pcld.CalculateDensityProjection(self.Display, xc, yc, zc, xn, yn, zn, length_cam)
                self.pcld.density_test = density_test
                self.pcld.xn = xn
                self.pcld.yn = yn
                self.pcld.zn = zn
                self.Display.calculate()
                self.Display.root.update()
                
                
            if (key == 'Right') : #self.pcld.SphFwd()
                #Amount 
                H_off += self.pcld.inc
                self.sphere.SetCenter(self.pcld.Sph_Xpos, self.pcld.Sph_Ypos, self.pcld.Sph_Zpos)
                return_com = 1
            if (key == 'Left') : #self.pcld.SphRev()
                H_off-= self.pcld.inc
                self.sphere.SetCenter(self.pcld.Sph_Xpos, self.pcld.Sph_Ypos, self.pcld.Sph_Zpos)
                return_com = 1
            if (key == 'Up') : 
                V_off+= self.pcld.inc
                self.sphere.SetCenter(self.pcld.Sph_Xpos, self.pcld.Sph_Ypos, self.pcld.Sph_Zpos)
                return_com = 1
            if (key == 'Down') :                
                V_off -= self.pcld.inc
                return_com = 1
            
            if return_com == 1 :   
                if (self.pcld.dominant_axis == 0) :
                    self.pcld.Sph_Xpos = 0.00
                    self.pcld.Sph_Ypos+= H_off
                    self.pcld.Sph_Zpos+= V_off
                else :
                    if (self.pcld.dominant_axis == 1) :
                        self.pcld.Sph_Xpos += H_off 
                        self.pcld.Sph_Ypos = 0.00
                        self.pcld.Sph_Zpos += V_off
                    else :
                        if (self.pcld.dominant_axis == 2) :
                            print("Adjust sphere position in XY plane ", H_off, V_off, self.pcld.Sph_Xpos, self.pcld.Sph_Ypos)
                            self.pcld.Sph_Xpos += H_off
                            self.pcld.Sph_Ypos += V_off
                            self.pcld.Sph_Zpos = 0.00
            
                self.sphere.SetCenter(self.pcld.Sph_Xpos, self.pcld.Sph_Ypos, self.pcld.Sph_Zpos)
                self.Display.renderWindow.Render()
                self.Display.calculate()
                self.Display.root.update()
                self.pcld.DensityOld = True               
            
                return
            
                                
        if (self.pcld.state == 2) :  #
            self.allow_dominant_axis_change = False
            if (key == 'Up') :
                self.pcld.Rmin = self.pcld.Rmin**0.5
                print("Min Rsquared value = ", self.pcld.Rmin)
                self.pcld.Undo()
                self.pcld.Group()
            if (key == 'Down' ) : 
                self.pcld.Rmin = self.pcld.Rmin**2
                print("Min Rsquared value = ", self.pcld.Rmin)
                self.pcld.Undo()
                self.pcld.Group()
            if (key == "Left") : 
                self.pcld.G0max = self.pcld.G0max * 0.9
                self.pcld.G1max = self.pcld.G1max * 0.9
                print("Max slope ", self.pcld.G0max)
                self.pcld.Undo()
                self.pcld.Group()
            if (key == "Right") :
                self.pcld.G0max = self.pcld.G0max * 1.1
                self.pcld.G1max = self.pcld.G1max * 1.1
                print("Max slope ", self.pcld.G0max)
                self.pcld.Undo()
                self.pcld.Group()
            if (key == 'space') :
                print("Space was pressed")
                self.pcld.Group()
                #self.pcld.SaveGroupResults()
            if (key == 'Return') and (self.lastkey != 'Return') :
                if (self.pcld.Grouped) :
                    self.pcld.FindLongestString()
                    self.pcld.FindPlane() 
                    self.pcld.DensityOld = True                    
                    self.pcld.AdvanceState()
                    
        if (self.pcld.state == 0) or (self.pcld.state == 1) :
            print("self.pcld.rot_or_tran", self.pcld.rot_or_tran)
            if (self.pcld.rot_or_tran == 0) :
                if (key == 'x') :  self.pcld.AxisRotate(0 , -1*self.pcld.inc)
                if (key == 'X') :  self.pcld.AxisRotate(0 , self.pcld.inc)
                if (key == 'y') :  self.pcld.AxisRotate(1 , -1*self.pcld.inc)
                if (key == 'Y') :  self.pcld.AxisRotate(1 , self.pcld.inc)
                if (key == 'z') :  self.pcld.AxisRotate(2 , -1*self.pcld.inc)
                if (key == 'Z') :  self.pcld.AxisRotate(2 , self.pcld.inc)
            else :
                if (key == 'x') :  
                    self.pcld.X_off -= self.pcld.inc
                    self.pcld.X -= self.pcld.inc
                    print("self.pcld.inc ", self.pcld.inc)
                if (key == 'X') :  
                    self.pcld.X_off += self.pcld.inc
                    self.pcld.X += self.pcld.inc
                if (key == 'y') :  
                    self.pcld.Y_off -= self.pcld.inc
                    self.pcld.Y -= self.pcld.inc
                if (key == 'Y') :  
                    self.pcld.Y_off += self.pcld.inc
                    self.pcld.Y += self.pcld.inc
                if (key == 'z') :  
                    self.pcld.Z_off -= self.pcld.inc
                    self.pcld.Z -= self.pcld.inc
                if (key == 'Z') :  
                    self.pcld.Z_off += self.pcld.inc
                    self.pcld.Z += self.pcld.inc
            
                
        if self.pcld.state == 1 :  #highlight main plane
            self.allow_dominant_axis_change = True
            
            if (key == 'U' ) :
                self.pcld.Kmax = self.pcld.Kmax + self.pcld.inc
                self.pcld.Bound()
            if (key == 'u' ) :
                self.pcld.Kmax = self.pcld.Kmax - self.pcld.inc
                self.pcld.Bound()
            if (key == 'L' ) :
                self.pcld.Kmin = self.pcld.Kmin + self.pcld.inc
                self.pcld.Bound()
            if (key == 'l' ) :
                self.pcld.Kmin = self.pcld.Kmin - self.pcld.inc
                self.pcld.Bound()
            if (key == 'Return') and (self.lastkey != 'Return') :
                if (self.pcld.IndexArray.sum() < self.pcld.MaxWorkingGroupSize) :
                    print("Advancing from state ", self.pcld.state, "to state ", self.pcld.state + 1)
                    self.pcld.AdvanceState()
                else:
                    print("Unable to advance until selected points is reduced below ", self.pcld.MaxWorkingGroupSize)  
               
        if self.pcld.state == 0 :  #manually orient
            self.allow_dominant_axis_change = True
            # if (key == 'x') :  self.pcld.AxisRotate(0 , -1*self.pcld.rot_inc)
            # if (key == 'X') :  self.pcld.AxisRotate(0 , self.pcld.rot_inc)
            # if (key == 'y') :  self.pcld.AxisRotate(1 , -1*self.pcld.rot_inc)
            # if (key == 'Y') :  self.pcld.AxisRotate(1 , self.pcld.rot_inc)
            # if (key == 'z') :  self.pcld.AxisRotate(2 , -1*self.pcld.rot_inc)
            # if (key == 'Z') :  self.pcld.AxisRotate(2 , self.pcld.rot_inc)
            if (key == 'Return') and (self.lastkey != 'Return') :
                
                
                self.pcld.AxisRotate(0 , 0)
                self.pcld.AxisRotate(1 , 0)
                self.pcld.AxisRotate(2 , 0)
                
                self.pcld.InitBound()
                
                #self.pcld.rot_inc = self.rot_inc * 0.1
                self.pcld.AdvanceState()
                
                
        self.lastkey = key
        
                # slider_rounded = (int(100*self.root.sc_01.get())/100)
        # self.root.db_06d.set(slider_rounded)
        # self.pcl.inc = float(slider_rounded)
        
        if (key == 'braceright' ) :
            self.pcld.IncInc()        
            self.Display.root.sc_01.set(self.pcld.inc)
            slider_rounded = (int(100*self.Display.root.sc_01.get())/100)
            self.Display.root.db_06d.set(slider_rounded)
        #if (key == 'bracketright' ) : 
        if (key == 'braceleft' ) : 
            self.pcld.DecInc()
            self.Display.root.sc_01.set(self.pcld.inc)
            slider_rounded = (int(100*self.Display.root.sc_01.get())/100)
            self.Display.root.db_06d.set(slider_rounded)
            
        if (key == 'i') :
            print("Roll part along the horizontal axis")
            self.Display.renderer.GetActiveCamera().Elevation(self.pcld.inc)
            
        if (key == 'I') :
            self.Display.renderer.GetActiveCamera().Elevation(-1 * self.pcld.inc)
            
        if (key == 'm') :
            print("Aximuth")
            self.Display.renderer.GetActiveCamera().Azimuth(self.pcld.inc)
            
        if (key == 'M') :
            self.Display.renderer.GetActiveCamera().Azimuth(-1 * self.pcld.inc)
        
        if (key == 'o') :
            self.Display.renderer.GetActiveCamera().Roll(self.pcld.inc)
            
        if (key == 'O') :
            self.Display.renderer.GetActiveCamera().Roll(-1*self.pcld.inc)
            
        if (key == 'BackSpace') : 
            if (self.pcld.state >= 1) : self.pcld.state -= 1
        if (key == 'Tab') : 
            self.pcld.AdvanceState()
            self.Display.calculate()
            self.Display.root.update()
            return
        if (key == 'R') : self.pcld.SetOnOrigin()
        if (key == 'r') : self.pcld.IndexArray[:] = True
#        if (key == 'p') : self.pcld.FindPlane()
        #if (key == 's') or (key == 'S') : self.pcld.Save()
        # if (key == 'h') or (key == 'H') :
            # self.pcld.Hist()
            # length = len(self.pcld.X)
        # if (key == 'd') : 
            # self.pcld.PtSizeMinus()
        # if (key == 'D') : 
            # self.pcld.PtSizePlus()
        # if (key == 't') : 
            # self.Display.TDP_X.PtSizeMinus()
            # self.Display.TDP_Y.PtSizeMinus()
            # self.Display.TDP_Z.PtSizeMinus()
        # if (key == 'T') : 
            # self.Display.TDP_X.PtSizePlus()
            # self.Display.TDP_Y.PtSizePlus()
            # self.Display.TDP_Z.PtSizePlus()
        if (key == 'v') : 
            self.Display.TDP_X.ZoomIn()
            self.Display.TDP_Y.ZoomIn()
            self.Display.TDP_Z.ZoomIn()
        if (key == 'V') : 
            self.Display.TDP_X.ZoomOut()
            self.Display.TDP_Y.ZoomOut()
            self.Display.TDP_Z.ZoomOut()
        if (key == 'minus') :
            #self.pcld.RangePlus()
            self.pcld.IncreaseDisplayFactor()
            print("Factor Increased to ", self.pcld.factor)
            #print("Increased range being removed from histogram")
        if (key == 'plus') :
            #self.pcld.RangeMinus()
            self.pcld.ReduceDisplayFactor()
            print("Factor Reduced to ", self.pcld.factor)
            #print("Reduced range being removed from histogram")



        if (key == 'w') or (key == 'W') :
            self.pcld.Undo()
 #       if (key == 'f') or (key == 'F') :
 #           self.pcld.FindLongestString()

        if (key == 'Shift_L') : return
        # if (key == 'c' ) : 
            # print("c was pressed")
            # tfc = PCF.TimeStamp(0,"")
            # lsfc.find_circle(self.pcld.XdispB, self.pcld.YdispB)
            # tfc = PCF.TimeStamp(tfc, "Find Circle")
            
        xn, yn, zn = self.Display.FindViewNormal()
        print("Prior to clearpoints, addpoints:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
            
        ts = PCF.TimeStamp(ts, 'Check Selection')
        self.pcld.clearPoints()
        ts = PCF.TimeStamp(ts, 'Add Points')
        self.pcld.addPoints()
        
        xn, yn, zn = self.Display.FindViewNormal()
        print("After clearpoints, addpoints:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
        
        ts = PCF.TimeStamp(ts, 'AddPoints')
        #print("key %s" , key)
        #ts = PCF.TimeStamp(ts, 'Print')
       
        #self.pcld.IndexComb = np.logical_and(self.pcld.IndexArray, self.pcld.IndexDisp)

        self.Display.TDP_X.Update(self.pcld.YdispA, self.pcld.ZdispA, self.pcld.YdispB, self.pcld.ZdispB)
        self.Display.TDP_Y.Update(self.pcld.XdispA, self.pcld.ZdispA, self.pcld.XdispB, self.pcld.ZdispB)
        self.Display.TDP_Z.Update(self.pcld.XdispA, self.pcld.YdispA, self.pcld.XdispB, self.pcld.YdispB)
        ts = PCF.TimeStamp(ts, 'UpdateTwoDPlots')
        
        xn, yn, zn = self.Display.FindViewNormal()
        print("Prior to last render:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))
                
        #renderWindow.Render()
        self.Display.renderWindow.Render() 
        print("self.Display ", self.Display)
        self.Display.calculate()
        self.Display.root.update()
        ts = PCF.TimeStamp(ts, 'Render')
        tl = PCF.TimeStamp(tl, 'Total loop')
        
        xn, yn, zn = self.Display.FindViewNormal()
        print("AfterLoop:  xn, yn, zn", xn, yn, zn, np.sqrt(xn*xn + yn*yn+zn*zn))

         
class FourPanelDisplay():
    
    def __init__(self): 
# Renderer

        #print("FPD_00")
        self.root = Tk()
        

        
        self.root.title("Point cloud console")
        

        mainframe = ttk.Frame(self.root, padding="3 3 12 12")
        mainframe.grid(column=4, row=20, sticky=(N, W, E, S))
        mainframe.columnconfigure(0, weight=1)
        mainframe.rowconfigure(0, weight=1)
        print (type(mainframe))
        self.root.update()
        MessageText = StringVar()
        MessageText.set("Importing subroutines ")
        # ttk.Label(self.root, 
# #            text = "Importing subroutines.",
            # textvariable = MessageText,
            # justify = LEFT,
            # padx = 20).grid(column=1, row=1, sticky=E)    
            
        ttk.Label(self.root, 
#            text = "Importing subroutines.",
            textvariable = MessageText
            #              MessageText
            ).grid(column=1, row=1, sticky=E)    
            
        #ttk.Label(mainframe, textvariable = self.root.db_01a).grid(column=1, row=1, sticky=E)

        self.root.update()
        
        inputfile = ""
        Tk().withdraw()
        self.root.withdraw()
        inputfile = filedialog.askopenfilename(defaultextension=".xyz", filetypes =(("xyz", "*.xyz"),("All Files","*.*")),
            title = "Choose a file."
            )
        
        #
        self.root.deiconify()
        MessageText.set("Loading data and initializing arrays.  This can take some time.")      
        self.root.update()
        #Tk().withdraw()
        
        #ttk.Label(Tk(), textvariable = self.root.db_01a).grid(column=1, row=1, sticky=E)
        
        
        # Label(self.root, 
              # text = "Loading data and initializing arrays.  This can take some time.",
              # justify = LEFT,
              # padx = 20).grid(column=1, row=1, sticky=E)
              
        
        
        
        
        # Label(self.root, 
              # text = "Loading data and initializing arrays.  This can take some time.",
              # justify = LEFT,
              # padx = 20).grid(column=1, row=1, sticky=E)
        
        if (inputfile == "") : return
        #Tk().withdraw() 
                
        #inputfile = filedialog.askopenfilename()
        #inputfile = "C:/Python34/PointCloud/3_LT1_Reduced.xyz"
        (Xpcl,Ypcl,Zpcl) = PCF.PandaLoad(inputfile)
        
        
        self.SI = SphericalIndicator()
        pcl = VtkPointCloud(Xpcl, Ypcl, Zpcl)
        #print("FPD_01")
        pcl.CreatePointCloud()
        self.pcl = pcl
        #print("FPD_02")
        # X_st_min = pcl.X.min()
        # X_st_max = pcl.X.max()
        # Y_st_min = pcl.Y.min()
        # Y_st_max = pcl.Y.max()
        # Z_st_min = pcl.Z.min()
        # Z_st_max = pcl.Z.max()
        # GraphMinOverall = 20*int(min(X_st_min, Y_st_min, Z_st_min)/20)-40
        # GraphMinOverall = min(GraphMinOverall, -20)
        # GraphMaxOverall = 2*(20*int(max(X_st_max, Y_st_max, Z_st_max)/20)+40)
        
        # print("X_st_min = ", X_st_min)
        # print("Y_st_min = ", Y_st_min)
        # print("Z_st_min = ", Z_st_min)
        # print("GraphMinOverall =", GraphMinOverall)
        # print("GraphMaxOverall =", GraphMaxOverall)
        #print("FPD_03")
                
        hor = 800
        ver = 600
        
        self.root.geometry('+' + str(hor+20) + '+10')
        self.root.update()

        
        menu = Menu(self.root)
        self.root.config(menu=menu)

        file = Menu(menu)

        file.add_command(label = 'Save', command = self.save_subset)
        file.add_command(label = 'Export Features', command = self.export_features)
        file.add_command(label = 'Export Projection', command = self.export_projection)
        file.add_command(label = 'Save Group Results', command = self.save_group_results)
        # #if (key == 's') or (key == 'S') : self.pcld.Save()
        file.add_command(label = 'Exit', command = lambda:exit())

        menu.add_cascade(label = 'File', menu = file)
    

        xmins=[0,.5,0,.5]
        xmaxs=[0.5,1,0.5,1]
        ymins=[0,0,.5,.5]
        ymaxs=[0.5,0.5,1,1]
    
        self.renderer = vtk.vtkRenderer()
        self.renderer.AddActor(pcl.vtkActor)
        self.renderer.AddActor(pcl.vtkActor2)
        self.renderer.SetBackground(1.0, 1.0, 1.0)
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().ParallelProjectionOn()
        #print("FPD_04")
 
# Render Window
        self.renderWindow = vtk.vtkRenderWindow()
        data_min=min(pcl.X.min(),pcl.Y.min(),pcl.Z.min())
        data_max=max(pcl.X.max(),pcl.Y.max(),pcl.Z.max())
        cubeAxesActor = CreateCubeAxes(data_min,data_max,"All")
        cubeAxesActor.SetCamera(self.renderer.GetActiveCamera())
        self.renderer.AddActor(cubeAxesActor)
        
        self.sphere = vtk.vtkSphereSource()
        #self.sphere.SetCenter(pcl.X[self.SI.index], pcl.Y[self.SI.index], pcl.Z[self.SI.index])
        self.sphere.SetCenter(0, 0, 0)
        self.sphere.SetRadius(pcl.Sph_Radius)
        self.sphere.SetPhiResolution(20)
        self.sphere.SetThetaResolution(20)

        smapper = vtk.vtkPolyDataMapper()
        smapper.SetInputConnection(self.sphere.GetOutputPort())
        #smapper.ImmediateModeRenderingOn()
        #The above is thought to improve rendering speed of very large files.  I could not think of an objective test
        # for this, but subjectively it made no difference.
        #http://www.vtk.org/Wiki/VTK/FAQ#Are_there_bindings_to_languages_other_than_Tcl.3F
        actor = vtk.vtkActor()
        actor.GetProperty().SetColor(1.0, 1.0, 0.0)
        actor.SetMapper(smapper)
        self.renderer.AddActor(actor)
        #ClipRange = 10000
        #renderer.ResetCameraClippingRange(-1*ClipRange,ClipRange, -1*ClipRange, ClipRange, -1*ClipRange, ClipRange)

        #print("FPD_05")
        origin = vtk.vtkSphereSource()
        origin.SetCenter(0,0,0)
        origin.SetRadius(5)
        origin.SetPhiResolution(20)
        origin.SetThetaResolution(20)

        origin_mapper = vtk.vtkPolyDataMapper()
        origin_mapper.SetInputConnection(origin.GetOutputPort())
        origin_actor = vtk.vtkActor()
        origin_actor.GetProperty().SetColor(0.0, 0.0, 1.0)
        
        origin_actor.SetMapper(origin_mapper)
        self.renderer.AddActor(origin_actor)
        #vtkCamera::SetParallelProjectionOn()), 
        
        self.renderer.ResetCamera()
         
        self.renderWindow = vtk.vtkRenderWindow()
        #renderWindow.SetWindowName("Interactive")
        self.renderWindow.AddRenderer(self.renderer)
        self.renderer.SetViewport(xmins[0],ymins[0],xmaxs[0],ymaxs[0])
        
        self.renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self.renderWindowInteractor.SetRenderWindow(self.renderWindow)
        #print("FPD_06")
        #def __init__(self,parent, Display, pointCloud, SphInt, Sphere):
        self.renderWindowInteractor.SetInteractorStyle(MyInteractorStyle(self.renderWindowInteractor, self, pcl, self.SI, self.sphere))
        #self
        #print("FPD_07")

         
        self.renderWindow.SetSize(hor,ver)

        self.XcontextView = vtk.vtkContextView()
        self.YcontextView = vtk.vtkContextView()
        self.ZcontextView = vtk.vtkContextView()
        self.TDP_Init = False
#2        self.pcl.clearPoints()
#2        self.pcl.addPoints()
        self.ResetGraphExtents()
        #print("FPD_08")
        self.TDP_X = TwoDPlot(self.XcontextView, pcl.Y,pcl.Z,"Y","Z", self.pcl.GraphMinY, self.pcl.GraphMinZ, 
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Y = TwoDPlot(self.YcontextView, pcl.X,pcl.Z,"X","Z", self.pcl.GraphMinX, self.pcl.GraphMinZ,
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Z = TwoDPlot(self.ZcontextView, pcl.X,pcl.Y,"X","Y", self.pcl.GraphMinX, self.pcl.GraphMinY, 
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Init = True
        #print("FPD_09")
        #print(self.TDP_X)
        
        self.Xrenderer = self.TDP_X.view.GetRenderer()
        self.Yrenderer = self.TDP_Y.view.GetRenderer()
        self.Zrenderer = self.TDP_Z.view.GetRenderer()
        
        self.renderWindow.AddRenderer(self.Xrenderer)
        self.Xrenderer.SetViewport(xmins[1],ymins[1],xmaxs[1],ymaxs[1])
        self.renderWindow.AddRenderer(self.Yrenderer)
        self.Yrenderer.SetViewport(xmins[2],ymins[2],xmaxs[2],ymaxs[2])
        self.renderWindow.AddRenderer(self.Zrenderer)
        self.Zrenderer.SetViewport(xmins[3],ymins[3],xmaxs[3],ymaxs[3])

        #print("FPD_10")
        
        self.renderWindow.Render()
        self.renderWindow.SetWindowName("Interactive")
        
        ###############tkinter window
        

        
        self.root.db_01a = StringVar()
        self.root.db_01b = StringVar()
        self.root.db_01c = StringVar()
        self.root.db_01d = StringVar()
        self.root.db_02a = StringVar()
        self.root.db_02b = StringVar()
        self.root.db_02c = StringVar()
        self.root.db_02d = StringVar()
        self.root.db_03a = StringVar()
        self.root.db_03b = StringVar()
        self.root.db_03c = StringVar()
        self.root.db_03d = StringVar()
        self.root.db_04a = StringVar()
        self.root.db_04b = StringVar()
        self.root.db_04c = StringVar()
        self.root.db_04d = StringVar()
        self.root.db_05a = StringVar()
        self.root.db_05b = StringVar()
        self.root.db_05c = StringVar()
        self.root.db_05d = StringVar()
        self.root.db_06a = StringVar()
        self.root.db_06b = StringVar()
        self.root.db_06c = StringVar()
        self.root.db_06d = StringVar()
        self.root.db_07a = StringVar()
        self.root.db_07b = StringVar()
        self.root.db_07c = StringVar()
        self.root.db_07d = StringVar()
        self.root.db_08a = StringVar()
        self.root.db_08b = StringVar()
        self.root.db_08c = StringVar()
        self.root.db_08d = StringVar()
        self.root.db_09a = StringVar()
        self.root.db_09b = StringVar()
        self.root.db_09c = StringVar()
        self.root.db_09d = StringVar()
        self.root.db_10a = StringVar()
        self.root.db_10b = StringVar()
        self.root.db_10c = StringVar()
        self.root.db_10d = StringVar()
        self.root.db_11a = StringVar()
        self.root.db_11b = StringVar()
        self.root.db_11c = StringVar()
        self.root.db_11d = StringVar()
        self.root.db_12a = StringVar()
        self.root.db_12b = StringVar()
        self.root.db_12c = StringVar()
        self.root.db_12d = StringVar()
        self.root.db_13a = StringVar()
        self.root.db_13b = StringVar()
        self.root.db_13c = StringVar()
        self.root.db_13d = StringVar()
        
        self.root.db_14a = StringVar()
        self.root.db_14b = StringVar()
        self.root.db_14c = StringVar()
        self.root.db_14d = StringVar()
        
        self.root.db_15a = StringVar()
        self.root.db_15b = StringVar()
        self.root.db_15c = StringVar()
        self.root.db_15d = StringVar()
        self.root.db_16a = StringVar()
        self.root.db_16b = StringVar()
        self.root.db_16c = StringVar()
        self.root.db_16d = StringVar()
        self.root.db_17a = StringVar()
        self.root.db_17b = StringVar()
        self.root.db_17c = StringVar()
        self.root.db_17d = StringVar()        
        self.root.db_18a = StringVar()
        self.root.db_18b = StringVar()
        self.root.db_18c = StringVar()
        self.root.db_18d = StringVar()
        self.root.db_19a = StringVar()
        self.root.db_19b = StringVar()
        self.root.db_19c = StringVar()
        self.root.db_19d = StringVar()        
        self.root.db_20a = StringVar()
        self.root.db_20b = StringVar()
        self.root.db_20c = StringVar()
        self.root.db_20d = StringVar() 
        self.root.db_21a = StringVar()
        self.root.db_21b = StringVar()
        self.root.db_21c = StringVar()
        self.root.db_21d = StringVar() 
        self.root.db_22a = StringVar()
        self.root.db_22b = StringVar()
        self.root.db_22c = StringVar()
        self.root.db_22d = StringVar()
        
        self.root.db_23a = StringVar()
        self.root.db_23b = StringVar()
        self.root.db_23c = StringVar()
        self.root.db_23d = StringVar()

        self.root.db_24a = StringVar()
        self.root.db_24b = StringVar()
        self.root.db_24c = StringVar()
        self.root.db_24d = StringVar()

        self.root.db_25a = StringVar()
        self.root.db_25b = StringVar()
        self.root.db_25c = StringVar()
        self.root.db_25d = StringVar()








        
        
        self.root.db_01a.set("State")
        self.root.db_01b.set(str(self.pcl.state))
        self.root.db_01c.set("")
        self.root.db_01d.set("")
        self.root.db_02a.set(self.pcl.statenames[self.pcl.state])
        self.root.db_02b.set("")
        self.root.db_02c.set("")
        self.root.db_03a.set(self.pcl.helptext1[self.pcl.state])
        self.root.db_03b.set("")
        self.root.db_03c.set("")
        self.root.db_04a.set(self.pcl.helptext2[self.pcl.state])
        self.root.db_04b.set("")
        self.root.db_04c.set("")
        self.root.db_05a.set(self.pcl.helptext3[self.pcl.state])
        self.root.db_05b.set("")
        self.root.db_05c.set("")
        self.root.db_06a.set(self.pcl.helptext4[self.pcl.state])
        self.root.db_06b.set("")
        self.root.db_06c.set("")
        self.root.db_07a.set(self.pcl.helptext5[self.pcl.state])
        self.root.db_07b.set("")
        self.root.db_07c.set("")
        self.root.db_08a.set("")
        self.root.db_08b.set("")
        self.root.db_08c.set("")
        self.root.db_09a.set(inputfile + "   Decimation (1 = no decimation)")
        self.root.db_09b.set(self.pcl.BigD)
        self.root.db_09c.set("")
        self.root.db_09d.set(str(self.pcl.factor))
        self.root.db_10a.set("Collapse values.  Toggle, offset, pitch")
        self.root.db_10b.set("")
        self.root.db_10c.set("")

        self.root.db_11a.set("X boundaries: ")
        self.root.db_11b.set("")
        self.root.db_11c.set("")
        self.root.db_11d.set("")

        self.root.db_12a.set("Y boundaries: ")
        self.root.db_12b.set("")
        self.root.db_12c.set("")
        self.root.db_12d.set("")

        self.root.db_13a.set("Z boundaries: ")
        self.root.db_13b.set("")
        self.root.db_13c.set("")  
        self.root.db_13d.set("")    

        self.root.db_14a.set("Upper and lower boundaries (slice): ")
        self.root.db_14b.set("")
        self.root.db_14c.set("")  
        self.root.db_14d.set("")         
      
        self.root.db_15a.set("Dominant axis")
        self.root.db_15b.set("")
        self.root.db_15c.set("")
        self.root.db_16a.set("Offset in the X, Y, and Z axis")
        self.root.db_16b.set("")
        self.root.db_16c.set("")
        self.root.db_17a.set("Position of spherical cursor")
        self.root.db_17b.set("")
        self.root.db_17c.set("")
        self.root.db_18a.set("Change size of points in primary cloud")
        self.root.db_18b.set("")
        self.root.db_18c.set("")
        self.root.db_19a.set("Change size of points in secondary cloud")
        self.root.db_19b.set("")
        self.root.db_19c.set("")
        self.root.db_20a.set("Change size of spherical pointer")
        self.root.db_20b.set("")
        self.root.db_20c.set("")
        self.root.db_21a.set("Enter values for find longest string thresholds")
        self.root.db_21b.set("")
        self.root.db_21c.set("")
        self.root.db_22a.set("")
        self.root.db_22b.set("")
        self.root.db_22c.set("")
        self.root.db_23a.set("Zoom in or out with V or v")
        self.root.db_23b.set("")
        self.root.db_23c.set("")
        self.root.db_24a.set("Normal vector for adjustable view port")
        self.root.db_24b.set("")
        self.root.db_24c.set("")
        self.root.db_25a.set("                                                                                 ")
        self.root.db_25b.set("         ")
        self.root.db_25c.set("         ")
        self.root.db_25d.set("         ")
        
        #print("FPD_11")
        
        #print("self.root.db_02a.get", self.root.db_02a.get())
        ttk.Label(mainframe, textvariable = self.root.db_01a).grid(column=1, row=1, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_01b).grid(column=2, row=1, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_01c).grid(column=3, row=1, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_01d).grid(column=4, row=1, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_02a).grid(column=1, row=2, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_02b).grid(column=2, row=2, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_02c).grid(column=3, row=2, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_02d).grid(column=4, row=2, sticky=E)

        ttk.Label(mainframe, textvariable = self.root.db_03a).grid(column=1, row=3, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_03b).grid(column=2, row=3, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_03c).grid(column=3, row=3, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_03d).grid(column=4, row=3, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_04a).grid(column=1, row=4, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_04b).grid(column=2, row=4, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_04c).grid(column=3, row=4, sticky=E)
        self.root.rb_01 = IntVar()
        self.root.rb_01.set(0)
        ttk.Radiobutton(mainframe, text="Rotate", variable=self.root.rb_01, 
            command = self.toggle_rot_trans, value=0).grid(column = 2, columnspan = 3, row = 7, sticky = W)
        ttk.Radiobutton(mainframe, text="Translate", variable=self.root.rb_01, 
            command = self.toggle_rot_trans, value=1).grid(column = 2, columnspan = 3, row = 8, sticky = W)
            
        self.root.rb_02 = IntVar()
        self.root.rb_02.set(2)
        ttk.Radiobutton(mainframe, text="X", variable=self.root.rb_02, 
            command = self.toggle_dominant_axis, value=0).grid(column = 2, columnspan = 1, row = 15, sticky = W)
        ttk.Radiobutton(mainframe, text="Y", variable=self.root.rb_02, 
            command = self.toggle_dominant_axis, value=1).grid(column = 3, columnspan = 1, row = 15, sticky = W)
        ttk.Radiobutton(mainframe, text="Z", variable=self.root.rb_02, 
            command = self.toggle_dominant_axis, value=2).grid(column = 4, columnspan = 1, row = 15, sticky = W)
            
            
        self.root.sc_01 = DoubleVar()  
        self.root.sc_01.set(1.0)
        ttk.Scale(mainframe, variable = self.root.sc_01, from_= 0.01, to = 10, length = 160, 
            orient = HORIZONTAL, command = self.slide1_update).grid(column = 2, columnspan = 3, row = 6, sticky = W)  


        # self.root.sc_11b = DoubleVar()  
        # self.root.sc_11b.set(1)
        # ttk.Scale(mainframe, variable = self.root.sc_11b, from_= -1000, to = 1000, length = 160, 
            # orient = HORIZONTAL, command = self.slide11b_update).grid(column = 2, columnspan = 2, row = 11, sticky = W) 

        # self.root.sc_11c = DoubleVar()  
        # self.root.sc_11c.set(1)
        # ttk.Scale(mainframe, variable = self.root.sc_11c, from_= -1000, to = 1000, length = 160, 
            # orient = HORIZONTAL, command = self.slide11c_update).grid(column = 3, columnspan = 2, row = 11, sticky = W)             
            
        # self.root.sc_12b = DoubleVar()  
        # self.root.sc_12b.set(1)
        # ttk.Scale(mainframe, variable = self.root.sc_12b, from_= -1000, to = 1000, length = 80, 
            # orient = HORIZONTAL, command = self.slide12b_update).grid(column = 2, columnspan = 1, row = 12, sticky = W) 

        # self.root.sc_12c = DoubleVar()  
        # self.root.sc_12c.set(1)
        # ttk.Scale(mainframe, variable = self.root.sc_12c, from_= -1000, to = 1000, length = 80, 
            # orient = HORIZONTAL, command = self.slide12c_update).grid(column = 3, columnspan = 1, row = 12, sticky = W)   
            
            
        # self.root.sc_13b = DoubleVar()  
        # self.root.sc_13b.set(1)
        # ttk.Scale(mainframe, variable = self.root.sc_13b, from_= -1000, to = 1000, length = 80, 
            # orient = HORIZONTAL, command = self.slide13b_update).grid(column = 2, columnspan = 1, row = 13, sticky = W) 

        # self.root.sc_13c = DoubleVar()  
        # self.root.sc_13c.set(1)
        # ttk.Scale(mainframe, variable = self.root.sc_13c, from_= -1000, to = 1000, length = 80, 
            # orient = HORIZONTAL, command = self.slide13c_update).grid(column = 3, columnspan = 1, row = 13, sticky = W)   

            
        self.root.db_11b.set(str(10*int(self.pcl.X_limit_lower)/10))
        self.root.db_11c.set(str(10*int(self.pcl.X_limit_upper)/10))        
        self.root.db_12b.set(str(10*int(self.pcl.Y_limit_lower)/10))   
        self.root.db_12c.set(str(10*int(self.pcl.Y_limit_upper)/10))            
        self.root.db_13b.set(str(10*int(self.pcl.Z_limit_lower)/10))            
        self.root.db_13c.set(str(10*int(self.pcl.Z_limit_upper)/10))            
            
        self.root.db_14b.set(str(10*int(self.pcl.slice_min)/10))
        self.root.db_14c.set(str(10*int(self.pcl.slice_max)/10))        
        
        self.tb_10 = -1 
        
        #print (type(self))
        self.toggle_btn_collapse = tk.Button(mainframe, text="Collapse?", bg = 'red', width=12, command = self.toggle_10)
        self.toggle_btn_collapse.grid(column = 2, row = 10, sticky = E)
#                ttk.Button(mainframe, text="Reset Scale", command=self.ResetGraphExtents).grid(column=2, row=22, sticky=W)
        #print (type(self.toggle_btn))
        #print (dir(self.toggle_btn))
        self.toggle_btn_collapse.config(bg = 'red')
        #toggle_btn = ttk.Button(text="Toggle", width=12, relief="raised", command = self.toggle_10).grid(column = 2, row = 10)
        #toggle_btn.pack(pady=5)             
            
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_10c).grid(column = 3, row = 10, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_10d).grid(column = 4, row = 10, sticky = E)
        
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_11b).grid(column = 2, row = 11, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_11c).grid(column = 3, row = 11, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_12b).grid(column = 2, row = 12, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_12c).grid(column = 3, row = 12, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_13b).grid(column = 2, row = 13, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_13c).grid(column = 3, row = 13, sticky = E)
        
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_14b).grid(column = 2, row = 14, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_14c).grid(column = 3, row = 14, sticky = E)
            

        self.root.sc_18 = DoubleVar()  
        self.root.sc_18.set(1)
        ttk.Scale(mainframe, variable = self.root.sc_18, from_= 1, to = 10, length = 160, 
            orient = HORIZONTAL, command = self.slide18_update).grid(column = 2, columnspan = 3, row = 18, sticky = W) 
            

        self.root.sc_19 = DoubleVar()  
        self.root.sc_19.set(1)
        ttk.Scale(mainframe, variable = self.root.sc_19, from_= 2, to = 10, length = 160, 
            orient = HORIZONTAL, command = self.slide19_update).grid(column = 2, columnspan = 3, row = 19, sticky = W)

        self.root.sc_20 = DoubleVar()  
        self.root.sc_20.set(self.pcl.Sph_Radius)
        ttk.Scale(mainframe, variable = self.root.sc_20, from_= 1, to = 50, length = 160, 
            orient = HORIZONTAL, command = self.slide20_update).grid(column = 2, columnspan = 3, row = 20, sticky = W)            

        # ttk.Scale(mainframe, variable = self.root.sc_01, from_=1, to = 200, resolution = 10, 
            # command = self.slider1_update, orient = HORIZONTAL).grid(column = 3, 
            # row = 6, sticky = E)             
        # w = Scale(master, variable = slider_v, from_=1, to=200, resolution = 10, orient=HORIZONTAL)
        
        ttk.Label(mainframe, textvariable = self.root.db_04d).grid(column=4, row=4, sticky=E)
                
        ttk.Label(mainframe, textvariable = self.root.db_05a).grid(column=1, row=5, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_05b).grid(column=2, row=5, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_05c).grid(column=3, row=5, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_05d).grid(column=4, row=5, sticky=E)
                        
        ttk.Label(mainframe, textvariable = self.root.db_06a).grid(column=1, row=6, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_06b).grid(column=2, row=6, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_06c).grid(column=3, row=6, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_06d).grid(column=4, row=6, sticky=E)
                                
        ttk.Label(mainframe, textvariable = self.root.db_07a).grid(column=1, row=7, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_07b).grid(column=2, row=7, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_07c).grid(column=3, row=7, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_07d).grid(column=4, row=7, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_08a).grid(column=1, row=8, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_08b).grid(column=2, row=8, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_08c).grid(column=3, row=8, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_08d).grid(column=4, row=8, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_09a).grid(column=1, row=9, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_09b).grid(column=2, row=9, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_09c).grid(column=3, row=9, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_09d).grid(column=4, row=9, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_10a).grid(column=1, row=10, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_10b).grid(column=2, row=10, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_10c).grid(column=3, row=10, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_10d).grid(column=4, row=10, sticky=E)
        
        
        ttk.Label(mainframe, textvariable = self.root.db_11a).grid(column=1, row=11, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_11b).grid(column=2, row=11, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_11c).grid(column=3, row=11, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_11d).grid(column=4, row=11, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_12a).grid(column=1, row=12, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_12b).grid(column=2, row=12, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_12c).grid(column=3, row=12, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_12d).grid(column=4, row=12, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_13a).grid(column=1, row=13, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_13b).grid(column=2, row=13, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_13c).grid(column=3, row=13, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_13d).grid(column=4, row=13, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_14a).grid(column=1, row=14, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_14b).grid(column=2, row=14, sticky=E)
#        ttk.Label(mainframe, textvariable = self.root.db_14c).grid(column=3, row=14, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_14d).grid(column=4, row=14, sticky=E)
        
        
        ttk.Label(mainframe, textvariable = self.root.db_15a).grid(column=1, row=15, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_15b).grid(column=2, row=15, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_15c).grid(column=3, row=15, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_15d).grid(column=4, row=15, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_16a).grid(column=1, row=16, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_16b).grid(column=2, row=16, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_16c).grid(column=3, row=16, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_16d).grid(column=4, row=16, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_17a).grid(column=1, row=17, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_17b).grid(column=2, row=17, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_17c).grid(column=3, row=17, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_17d).grid(column=4, row=17, sticky=E)

        ttk.Label(mainframe, textvariable = self.root.db_18a).grid(column=1, row=18, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_18b).grid(column=2, row=18, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_18c).grid(column=3, row=18, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_18d).grid(column=4, row=18, sticky=E)

        ttk.Label(mainframe, textvariable = self.root.db_19a).grid(column=1, row=19, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_19b).grid(column=2, row=19, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_19c).grid(column=3, row=19, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_19d).grid(column=4, row=19, sticky=E)

        ttk.Label(mainframe, textvariable = self.root.db_20a).grid(column=1, row=20, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_20b).grid(column=2, row=20, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_20c).grid(column=3, row=20, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_20d).grid(column=4, row=20, sticky=E)

        ttk.Label(mainframe, textvariable = self.root.db_21a).grid(column=1, row=21, sticky=E)
        # ttk.Label(mainframe, textvariable = self.root.db_21b).grid(column=2, row=21, sticky=E)
        # ttk.Label(mainframe, textvariable = self.root.db_21c).grid(column=3, row=21, sticky=E)
        # ttk.Label(mainframe, textvariable = self.root.db_21d).grid(column=4, row=21, sticky=E)

        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_21b).grid(column = 2, row = 21, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_21c).grid(column = 3, row = 21, sticky = E)
        
        
        ttk.Label(mainframe, textvariable = self.root.db_22a).grid(column=1, row=22, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_22b).grid(column=2, row=22, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_22c).grid(column=3, row=22, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_22d).grid(column=4, row=22, sticky=E)

        ttk.Label(mainframe, textvariable = self.root.db_23a).grid(column=1, row=23, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_23b).grid(column=2, row=23, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_23c).grid(column=3, row=23, sticky=E)
        #ttk.Label(mainframe, textvariable = self.root.db_23d).grid(column=4, row=23, sticky=E)        

        ttk.Label(mainframe, textvariable = self.root.db_24a).grid(column=1, row=24, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_24b).grid(column=2, row=24, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_24c).grid(column=3, row=24, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_24d).grid(column=4, row=24, sticky=E)
        
        ttk.Label(mainframe, textvariable = self.root.db_25a).grid(column=1, row=25, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_25b).grid(column=2, row=25, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_25c).grid(column=3, row=25, sticky=E)
        ttk.Label(mainframe, textvariable = self.root.db_25d).grid(column=4, row=25, sticky=E)
        
        ttk.Button(mainframe, text="Reset Scale", command=self.ResetGraphExtents).grid(column=2, row=22, sticky=W)
        
#2        self.root.db_21c.set(self.pcl.XdispA.min())
#2        self.root.db_21d.set(self.pcl.XdispA.max())
        # Y_st_min = self.pcl.YdispA.min()
        # Y_st_max = self.pcl.YdispA.max()
        # Z_st_min = self.pcl.ZdispA.min()
        # Z_st_max = self.pcl.ZdispA.max(
        
        
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_22c).grid(column = 3, row = 22, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_22d).grid(column = 4, row = 22, sticky = E) 
        ttk.Entry(mainframe,width = 8, textvariable=self.root.db_23b).grid(column = 2, row = 23, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_23c).grid(column = 3, row = 23, sticky = E)
        ttk.Entry(mainframe, width = 8, textvariable=self.root.db_23d).grid(column = 4, row = 23, sticky = E)        

        #self.root.columnconfigure(1,minsize = 150)
        self.root.columnconfigure(2,minsize = 100)
        self.root.columnconfigure(3,minsize = 100)
        self.root.columnconfigure(4,minsize = 100)
        
        for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

#        #print("FPD_12")
        self.pcl.clearPoints()
        self.pcl.addPoints()
        #self.root.bind('<Return>', self.calculate)
        self.root.bind('<Return>', self.downloadformdata)
        
#        #print("FPD_13")
        #self.root.db_22d.set("Check")
        self.calculate()
#        print("Prior to running calculate")
        MessageText.set("")
        self.root.mainloop()
#        print("FPD_14")
        self.calculate()
        #print("Did calculate just run")
 #       print("FPD_15")
        self.Rebuild()
        
        ##########tkinter window

    def FindViewNormal(self) :
        xn, yn, zn = self.renderer.GetActiveCamera().GetViewPlaneNormal()
        #print("Plane Normal obtained, ", xn, yn, zn)
        self.pcl.xn = xn
        self.pcl.yn = yn
        self.pcl.zn = zn
        return xn, yn, zn
        
    def Rebuild(self) :
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render() 
        
        
    def ResetGraphExtents(self) :
        #self.renderWindowInteractor.Start()
        X_st_min = self.pcl.X.min()
        X_st_max = self.pcl.X.max()
        Y_st_min = self.pcl.Y.min()
        Y_st_max = self.pcl.Y.max()
        Z_st_min = self.pcl.Z.min()
        Z_st_max = self.pcl.Z.max()

#2        X_st_min = self.pcl.XdispA.min()
#2        X_st_max = self.pcl.XdispA.max()
#2        Y_st_min = self.pcl.YdispA.min()
#2        Y_st_max = self.pcl.YdispA.max()
#2        Z_st_min = self.pcl.ZdispA.min()
#2        Z_st_max = self.pcl.ZdispA.max()
        
        
        self.pcl.GraphMinX = self.pcl.GraphStep*int(X_st_min/self.pcl.GraphStep)- 2 * self.pcl.GraphStep
        self.pcl.GraphMinY = self.pcl.GraphStep*int(Y_st_min/self.pcl.GraphStep)- 2 * self.pcl.GraphStep
        self.pcl.GraphMinZ = self.pcl.GraphStep*int(Z_st_min/self.pcl.GraphStep)- 2 * self.pcl.GraphStep
        
        self.pcl.GraphMaxX = self.pcl.GraphStep*int(X_st_max/self.pcl.GraphStep)+ 2 * self.pcl.GraphStep
        self.pcl.GraphMaxY = self.pcl.GraphStep*int(Y_st_max/self.pcl.GraphStep)+ 2 * self.pcl.GraphStep
        self.pcl.GraphMaxZ = self.pcl.GraphStep*int(Z_st_max/self.pcl.GraphStep)+ 2 * self.pcl.GraphStep
                
        self.pcl.GraphMinOverall = 20*int(min(X_st_min, Y_st_min, Z_st_min)/20)-40
        self.pcl.GraphMinOverall = min(self.pcl.GraphMinOverall, -20)
        self.pcl.GraphMaxOverall = 2*(20*int(max(X_st_max, Y_st_max, Z_st_max)/20)+40) 
        self.pcl.GraphRange = max(self.pcl.GraphMaxX - self.pcl.GraphMinX, self.pcl.GraphMaxY - self.pcl.GraphMinY,
            self.pcl.GraphMaxZ - self.pcl.GraphMinZ)
        self.pcl.GraphRange = self.pcl.GraphStep * int(self.pcl.GraphRange / self.pcl.GraphStep) + self.pcl.GraphStep

        if (self.TDP_Init) :
            self.TDP_X.ResetScale(self.pcl.GraphMinY, self.pcl.GraphMinZ, 
                self.pcl.GraphRange, self.pcl.GraphStep)
            self.TDP_Y.ResetScale(self.pcl.GraphMinX, self.pcl.GraphMinZ,
                self.pcl.GraphRange, self.pcl.GraphStep)
            self.TDP_Z.ResetScale(self.pcl.GraphMinX, self.pcl.GraphMinY, 
                self.pcl.GraphRange, self.pcl.GraphStep)
            self.renderWindow.Render()
            self.calculate()
            self.root.update()

        #print("self.X_st_min, self.X_st_max ", X_st_min, X_st_max)
        #print("self.Y_st_min, self.Y_st_max ", Y_st_min, Y_st_max)
        #print("self.Z_st_min, self.Z_st_max ", Z_st_min, Z_st_max)
        #print("self.pcl.GraphMinX, Y, Z", self.pcl.GraphMinX, self.pcl.GraphMinY, self.pcl.GraphMinZ)
        #print("self.GraphRange ", self.pcl.GraphRange)
        
    def toggle_10(self):
        self.tb_10 = self.tb_10 * -1
        print("self.tb_10", self.tb_10)
        if self.tb_10 > 0 :
            self.toggle_btn_collapse.config(bg = 'green')
            self.pcl.collapse = True
            self.pcl.collapse_X()
        else :
            self.toggle_btn_collapse.config(bg = 'red')
            self.pcl.X[:] = self.pcl.Xin[:]
            self.pcl.collapse = False
        self.Rebuild()
        
    def slide1_update(self, *args) : 
        #print("Value of slider", self.root.sc_01.get())
        slider_rounded = (int(100*self.root.sc_01.get())/100)
        self.root.db_06d.set(slider_rounded)
        self.pcl.inc = float(slider_rounded)
        
                

    def slide11b_update(self, *args) : 
        #print("Value of slider 11b", self.root.sc_11b.get())
        slider_rounded = (int(100*self.root.sc_11b.get())/100)
        self.pcl.X_limit_lower = float(slider_rounded)
        text_comb = str(int(100*self.pcl.X_limit_lower)/100) + " " + str(int(100*self.pcl.X_limit_upper)/100)
        self.root.db_11d.set(text_comb)
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render()         
        
        
    def slide11c_update(self, *args) : 
        #print("Value of slider 11c", self.root.sc_11c.get())
        slider_rounded = (int(100*self.root.sc_11c.get())/100)
        self.pcl.X_limit_upper = float(slider_rounded)
        text_comb = str(int(100*self.pcl.X_limit_lower)/100) + " " + str(int(100*self.pcl.X_limit_upper)/100)
        self.root.db_11d.set(text_comb)
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render() 
        
    def slide12b_update(self, *args) : 
        #print("Value of slider 12b", self.root.sc_12b.get())
        slider_rounded = (int(100*self.root.sc_12b.get())/100)
        self.pcl.Y_limit_lower = float(slider_rounded)
        text_comb = str(int(100*self.pcl.Y_limit_lower)/100) + " " + str(int(100*self.pcl.Y_limit_upper)/100)
        self.root.db_12d.set(text_comb)
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render()         
        
        
    def slide12c_update(self, *args) : 
        #print("Value of slider 12c", self.root.sc_12c.get())
        slider_rounded = (int(100*self.root.sc_12c.get())/100)
        self.pcl.Y_limit_upper = float(slider_rounded)
        text_comb = str(int(100*self.pcl.Y_limit_lower)/100) + " " + str(int(100*self.pcl.Y_limit_upper)/100)
        self.root.db_12d.set(text_comb)
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render() 
        
        
    def slide13b_update(self, *args) : 
        #print("Value of slider 13b", self.root.sc_13b.get())
        slider_rounded = (int(100*self.root.sc_13b.get())/100)
        self.pcl.Z_limit_lower = float(slider_rounded)
        text_comb = str(int(100*self.pcl.Z_limit_lower)/100) + " " + str(int(100*self.pcl.Z_limit_upper)/100)
        self.root.db_13d.set(text_comb)
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render()         
        
        
    def slide13c_update(self, *args) : 
        #print("Value of slider 13c", self.root.sc_13c.get())
        slider_rounded = (int(100*self.root.sc_13c.get())/100)
        self.pcl.Z_limit_upper = float(slider_rounded)
        text_comb = str(int(100*self.pcl.Z_limit_lower)/100) + " " + str(int(100*self.pcl.Z_limit_upper)/100)
        self.root.db_13d.set(text_comb)
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render() 
        
    def slide14b_update(self, *args) : 
        #print("Value of slider 14b", self.root.sc_14b.get())
        print("update slice_max", self.pcl.slice_max)
        slider_rounded = (int(100*self.root.sc_14b.get())/100)
        self.pcl.slice_min = float(slider_rounded)
        self.pcl.IndexSlice = np.where(np.logical_and(np.less(self.pcl.slice_dist,self.pcl.slice_max), \
            np.less(self.pcl.slice_min,self.pcl.slice_dist)), True,False)
        self.pcl.FindCylinder(self)
        self.Rebuild()
        # self.pcl.Z_limit_lower = float(slider_rounded)
        # text_comb = str(int(100*self.pcl.Z_limit_lower)/100) + " " + str(int(100*self.pcl.Z_limit_upper)/100)
        # self.root.db_14d.set(text_comb)
        # self.pcl.clearPoints()
        # self.pcl.addPoints()
        # self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        # self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        # self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        # self.renderWindow.Render()         
        
        
    def slide14c_update(self, *args) : 
        self.Rebuild()
        #print("Value of slider 14c", self.root.sc_14c.get())
        slider_rounded = (int(100*self.root.sc_14c.get())/100)
        self.pcl.slice_max = float(slider_rounded)
        self.pcl.IndexSlice = np.where(np.logical_and(np.less(self.pcl.slice_dist,self.pcl.slice_max), \
            np.less(self.pcl.slice_min,self.pcl.slice_dist)), True,False)
        self.pcl.FindCylinder(self)
        self.Rebuild()

    def slide18_update(self, *args) : 
        #print("Value of slider", self.root.sc_17.get())
        slider_rounded = (int(100*self.root.sc_18.get())/100)
        self.root.db_17d.set(slider_rounded)
        self.pcl.PtDisplaySize = float(slider_rounded)
        self.pcl.vtkActor.GetProperty().SetPointSize(self.pcl.PtDisplaySize)
        self.pcl.vtkActor2.GetProperty().SetPointSize(self.pcl.PtDisplaySize+1)
        #print("self.pcl.PtDisplaySize", self.pcl.PtDisplaySize)
        self.renderWindow.Render() 
        
    def slide19_update(self, *args) : 
        #print("Value of slider", self.root.sc_19.get())
        slider_rounded = (int(100*self.root.sc_19.get())/100)
        self.TDP_X.PtDisplaySize = slider_rounded
        self.TDP_Y.PtDisplaySize = slider_rounded
        self.TDP_Z.PtDisplaySize = slider_rounded
        # self.TDP_X.points.SetMarkerSize(self.TDP_X.PtDisplaySize)
        # self.TDP_Y.points.SetMarkerSize(self.TDP_Y.PtDisplaySize)
        # self.TDP_Z.points.SetMarkerSize(self.TDP_Z.PtDisplaySize)
        #self.Display.TDP_X.Update(self.pcld.YdispA, self.pcld.ZdispA, self.pcld.YdispB, self.pcld.ZdispB)
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        # self.TDP_Y.points.SetMarkerSize(self.TDP_Y.PtDisplaySize)
        # self.TDP_Z.points.SetMarkerSize(self.TDP_Z.PtDisplaySize)
        #print("self.TDP_Z.PtDisplaySize", self.TDP_Z.PtDisplaySize)
        self.root.db_18d.set(slider_rounded)
        #self.PtDisplaySize = float(slider_rounded)

        self.calculate()
        self.renderWindow.Render() 

        self.root.update()

    def slide20_update(self, *args) : 
        #print("Value of slider", self.root.sc_20.get())
        slider_rounded = (int(100*self.root.sc_20.get())/100)
        self.root.db_20d.set(slider_rounded)
        self.pcl.Sph_Radius = float(slider_rounded)
        #print("self.pcl.Sph_Radius",self.pcl.Sph_Radius)
        self.sphere.SetRadius(self.pcl.Sph_Radius)

        self.renderWindow.Render() 
        self.root.update()
        

    
    def toggle_rot_trans(self, *args) :
        # if (self.root.rb_01 == 1) : 
            # self.root.rb_01 = 0
        # else :
            # self.root.rb_01 = 1
        self.pcl.rot_or_tran = self.root.rb_01.get()
        #print("Toggle", self.root.rb_01.get(), self.pcl.rot_or_tran)
        
    def toggle_dominant_axis(self, *args) :
        if (self.pcl.allow_dominant_axis_change) : 
            self.pcl.dominant_axis = self.root.rb_02.get()
        else :
            self.rb_02.set(self.pcl.dominant_axis)

        #print("Toggle", self.root.rb_02.get(), self.pcl.dominant_axis)
        
        
    #def calculate(*args):
    
    def save_subset(self) :
    
        #self.root.withdraw()
        filepathandname = filedialog.asksaveasfilename(defaultextension=".txt", filetypes =(("xyz", "*.xyz"),("All Files","*.*")),
            title = "Save subset data."
            )
        self.pcl.outputfilepathandname = filepathandname
        #        outputfilepathandname
        self.pcl.Save()
        
    def export_features(self) :
    
        filepathandname = filedialog.asksaveasfilename(defaultextension=".txt", filetypes =(("txt", "*.txt"),("All Files","*.*")),
            title = "Export features."
            )
        self.pcl.outputfilepathandname = filepathandname
        #        outputfilepathandname
        self.pcl.ExportFeatures()
        
    def export_projection(self) :
    
        filepathandname = filedialog.asksaveasfilename(defaultextension=".txt", filetypes =(("txt", "*.txt"),("All Files","*.*")),
            title = "Export projection."
            )
        self.pcl.outputfilepathandname = filepathandname
        #        outputfilepathandname
        self.pcl.ExportProjection()

    def save_group_results(self) :
    
        filepathandname = filedialog.asksaveasfilename(defaultextension=".txt", filetypes =(("txt", "*.txt"),("All Files","*.*")),
            title = "Save Group Results."
            )
        self.pcl.outputfilepathandname = filepathandname
        #        outputfilepathandname
        self.pcl.SaveGroupResults()
        
        
    def downloadformdata(self, *args) :
        #print("downloadformdata")
        
        
#                self.collapse_off = 0.0
#        self.collapse_repeat = 0.00
        

#            self.root.db_10c.set(str('%.3f' % float(self.pcl.collapse_offset)))
#            self.root.db_10d.set(str('%.3f' % float(self.pcl.collapse_repeat)))            
            
        self.pcl.collapse_offset = float(self.root.db_10c.get())
        self.pcl.collapse_repeat = float(self.root.db_10d.get())
        
        self.pcl.X_limit_lower = float(self.root.db_11b.get())
        self.pcl.X_limit_upper = float(self.root.db_11c.get())
        self.pcl.Y_limit_lower = float(self.root.db_12b.get())
        self.pcl.Y_limit_upper = float(self.root.db_12c.get())
        self.pcl.Z_limit_lower = float(self.root.db_13b.get())
        self.pcl.Z_limit_upper = float(self.root.db_13c.get())
        
        self.pcl.slice_min = float(self.root.db_14b.get())
        self.pcl.slice_max = float(self.root.db_14c.get())
       
        self.pcl.Z_err1_lim = float(self.root.db_21b.get())
        self.pcl.Z_err2_lim = float(self.root.db_21c.get())
         
 
        self.pcl.GraphRange = float(self.root.db_22c.get())
        self.pcl.GraphStep = float(self.root.db_22d.get())
        self.pcl.GraphMinX = float(self.root.db_23b.get())
        self.pcl.GraphMinY = float(self.root.db_23c.get())
        self.pcl.GraphMinZ = float(self.root.db_23d.get())
        
        self.TDP_X.ResetScale(self.pcl.GraphMinY, self.pcl.GraphMinZ, 
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Y.ResetScale(self.pcl.GraphMinX, self.pcl.GraphMinZ,
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Z.ResetScale(self.pcl.GraphMinX, self.pcl.GraphMinY, 
            self.pcl.GraphRange, self.pcl.GraphStep)
            
        print("Download form data ", self.pcl.slice_min, self.pcl.slice_max)
        if (self.pcl.state >= 5) :
            print(self.pcl.IndexSlice.sum())
            self.pcl.IndexSlice = np.where(np.logical_and(np.less(self.pcl.slice_dist,self.pcl.slice_max), \
                np.less(self.pcl.slice_min,self.pcl.slice_dist)), True,False)
            print(self.pcl.IndexSlice.sum())
        if (self.pcl.state == 6) : self.pcl.HighlightPlane(self)           
        if (self.pcl.state == 5) : self.pcl.FindCylinder(self)
 
#        self.Rebuild()
        
        if (self.pcl.collapse) :
            if (self.pcl.collapse_offset_last != self.pcl.collapse_offset) or \
                (self.pcl.collapse_repeat_last != self.pcl.collapse_repeat) :
                self.pcl.collapse_offset_last = self.pcl.collapse_offset
                self.pcl.colllapse_repeat_last = self.pcl.collapse_repeat
                self.pcl.collapse_X()

        #self.ResetGraphExtents()
        
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        self.renderWindow.Render() 
        self.root.update()

    
    def calculate(self, *args):
        try:
            self.root.db_01b.set(int(self.pcl.state))
            self.root.db_02a.set(self.pcl.statenames[self.pcl.state])
            self.root.db_03a.set(self.pcl.helptext1[self.pcl.state])
            self.root.db_04a.set(self.pcl.helptext2[self.pcl.state])
            self.root.db_05a.set(self.pcl.helptext3[self.pcl.state])
            self.root.db_06a.set(self.pcl.helptext4[self.pcl.state])
            self.root.db_07a.set(self.pcl.helptext5[self.pcl.state])
            self.root.db_09c.set(int(self.pcl.LittleD_B))
            self.root.db_09d.set(str(int(self.pcl.factor)))     

            if self.pcl.collapse :
                self.toggle_btn_collapse.config(bg = 'green')
            else :
                self.toggle_btn_collapse.config(bg = 'red')  

            self.root.db_10c.set(str('%.3f' % float(self.pcl.collapse_offset)))
            self.root.db_10d.set(str('%.3f' % float(self.pcl.collapse_repeat)))
            
            self.root.db_14b.set(str(10*int(self.pcl.slice_min)/10))
            self.root.db_14c.set(str(10*int(self.pcl.slice_max)/10))

            #self.root.db_15b.set(str('%.2f' % float(self.pcl.X_rot)))
            #self.root.db_15c.set(str('%.2f' % float(self.pcl.Y_rot)))
            #self.root.db_15d.set(str('%.2f' % float(self.pcl.Z_rot)))
            self.root.db_16b.set(str('%.2f' % float(self.pcl.X_off)))
            self.root.db_16c.set(str('%.2f' % float(self.pcl.Y_off)))
            self.root.db_16d.set(str('%.2f' % float(self.pcl.Z_off)))
            self.root.db_17b.set(str('%.2f' % float(self.pcl.Sph_Xpos)))
            self.root.db_17c.set(str('%.2f' % float(self.pcl.Sph_Ypos)))
            self.root.db_17d.set(str('%.2f' % float(self.pcl.Sph_Zpos)))
            #self.root.db_17d.set(str(int(self.pcl.factor)))
            self.root.db_18d.set(str(int(self.root.sc_18.get())))
            self.root.db_19d.set(str(int(self.root.sc_19.get())))
            self.root.db_20d.set(str(int(self.root.sc_20.get())))
            #self.root.db_19c.set(str(self.pcl.rot_inc))
            #self.root.db_19d.set(str(self.pcl.off_inc))
            #self.root.db_20c.set(str(self.pcl.Sph_Radius))
            
            self.root.db_21b.set(str('%.3f' % float(self.pcl.Z_err1_lim)))
            self.root.db_21c.set(str('%.3f' % float(self.pcl.Z_err2_lim)))
            
            self.root.db_22c.set(str(int(self.pcl.GraphRange)))
            self.root.db_22d.set(str(int(self.pcl.GraphStep)))
            self.root.db_23b.set(str(int(self.pcl.GraphMinX)))
            self.root.db_23c.set(str(int(self.pcl.GraphMinY)))
            self.root.db_23d.set(str(int(self.pcl.GraphMinZ)))
            self.root.db_24b.set(str('%.3f' % float(self.pcl.xn)))
            self.root.db_24c.set(str('%.3f' % float(self.pcl.yn)))
            self.root.db_24d.set(str('%.3f' % float(self.pcl.zn)))
            self.root.db_25d.set(str(int(self.pcl.density_test)))
#            print("self.pcl.xn ", self.pcl.xn, " self.root.db_23c.get ", self.root.db_23c.get())



            
#            print("calculate")
            #meters.set((0.3048 * value * 10000.0 + 0.5)/10000.0)
        except ValueError:
            pass
            
def CreateCubeAxes(data_min,data_max,Axis):
    cubeAxesActor = vtk.vtkCubeAxesActor()

    Xr = 1
    Yr = 1
    Zr = 1
    if Axis == "X" : Xr = 0
    if Axis == "Y" : Yr = 0
    if Axis == "Z" : Zr = 0
    Xmin = data_min*Xr
    Xmax = data_max*Xr
    Ymin = data_min*Yr
    Ymax = data_max*Yr
    Kmin = data_min*Zr
    Kmax = data_max*Yr
 
      
    cubeAxesActor.SetBounds(Xmin,Xmax,Ymin,Ymax,Kmin,Kmax)
    
    cubeAxesActor.GetTitleTextProperty(0).SetColor(1.0, 0.0, 0.0)
    cubeAxesActor.GetLabelTextProperty(0).SetColor(1.0, 0.0, 0.0)
     
    cubeAxesActor.GetTitleTextProperty(1).SetColor(0.0, 1.0, 0.0)
    cubeAxesActor.GetLabelTextProperty(1).SetColor(0.0, 1.0, 0.0)
     
    cubeAxesActor.GetTitleTextProperty(2).SetColor(0.0, 0.0, 1.0)
    cubeAxesActor.GetLabelTextProperty(2).SetColor(0.0, 0.0, 1.0)
     
    cubeAxesActor.GetXAxesLinesProperty().SetColor(1.0, 0.0, 0.0)
    cubeAxesActor.GetYAxesLinesProperty().SetColor(0.0, 1.0, 0.0)
    cubeAxesActor.GetZAxesLinesProperty().SetColor(0.0, 0.0, 1.0)
    
    cubeAxesActor.GetXAxesGridlinesProperty().SetColor(1.0, 0.0, 0.0)
    cubeAxesActor.GetYAxesGridlinesProperty().SetColor(0.0, 1.0, 0.0)
    cubeAxesActor.GetZAxesGridlinesProperty().SetColor(0.0, 0.0, 1.0)
    
    if Axis != "X" : 
        Title = "YZ"
        cubeAxesActor.DrawXGridlinesOn()
    if Axis != "Y" : 
        Title = "XZ"
        cubeAxesActor.DrawYGridlinesOn()
        
    if Axis != "Z" : 
        Title = "XY"
        cubeAxesActor.DrawZGridlinesOn()
    cubeAxesActor.SetGridLineLocation(cubeAxesActor.VTK_GRID_LINES_FURTHEST)
     
    cubeAxesActor.XAxisMinorTickVisibilityOff()
    cubeAxesActor.YAxisMinorTickVisibilityOff()
    cubeAxesActor.ZAxisMinorTickVisibilityOff()
    return(cubeAxesActor)
    

def DisplayData() :
    import sys 
    #print("DD_01")


    FPD_out = FourPanelDisplay()

    
    
DisplayData()




