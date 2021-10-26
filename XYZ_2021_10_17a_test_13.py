#!/usr/bin/env python3    
 
 
from __future__ import print_function
import vtk
from numpy import random,genfromtxt,size
import numpy as np
import PointCloudFunctions_2016_12_31a as PCF
from vtk.util import numpy_support
import FT_20170106a as FT
from multiprocessing import Process
import os

from tkinter import *
from tkinter import filedialog
from tkinter import ttk
import tkinter as tk

import tkinterEntryWidget2 as dt

        
class VtkPointCloud:
    def __init__(self, X, Y, Z, Kmin=-10.0, Kmax=10.0, maxNumPoints=1e7):
        ti = PCF.TimeStamp(0,"Start")
        ti1 = ti

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
        
        
        self.IndexArray[:] = True
        self.IndexArray0[:] = True
        self.IndexArray1[:] = True
        self.IndexArray2[:] = True
        self.IndexArray3[:] = True

        self.IndexDispA = np.empty(len(self.X), dtype=np.bool_)
        self.IndexDispA[:] = True
        self.IndexDisp = np.empty(len(self.X), dtype=np.bool_)
        self.IndexDisp[:] = True
        self.IndexComb = np.empty(len(self.X), dtype=np.bool_)
        self.IndexComb[:] = True

        self.PtDisplaySize = 1.0

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


        self.factor = np.int(1)
        self.factor = max(int(self.factor / 10),1)
        self.LittleD_B = self.BigD
        self.i_disp = (self.factor)*np.array(range(np.int(self.BigD/self.factor)))

        
        self.state = 0
        self.state_max = 7

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
       
        self.inc = 1.0

        self.GraphStep = 20

        self.X_limit_lower = -10000
        self.X_limit_upper = 100000
        
        self.Y_limit_lower = -10000
        self.Y_limit_upper = 100000
        
        self.Z_limit_lower = -10000
        self.Z_limit_upper = 100000


    def addPoint(self, point):
        #print("add point")
        #if self.vtkPoints.GetNumberOfPoints() < self.maxNumPoints:
        pointId = self.vtkPoints.InsertNextPoint(point[:])
        self.vtkDepth.InsertNextValue(point[2])
        self.vtkCells.InsertNextCell(1)
        self.vtkCells.InsertCellPoint(pointId)
        #print("size of points = ", self.vtkCells.GetSize())

         
    def clearPoints(self):
        print("clearPoints")
        tcp = PCF.TimeStamp(0, "Start")
        self.vtkPoints = vtk.vtkPoints()
        self.vtkCells = vtk.vtkCellArray()
        self.vtkDepth = vtk.vtkDoubleArray()
        self.vtkDepth.SetName('DepthArray')
        self.vtkPolyData.SetPoints(self.vtkPoints)
        self.vtkPolyData.SetVerts(self.vtkCells)
        self.vtkPolyData.GetPointData().SetScalars(self.vtkDepth)
        self.vtkPolyData.GetPointData().SetActiveScalars('DepthArray')
        
        
        tcp = PCF.TimeStamp(tcp,"Cleared Points")
        

        
    def addPoints(self):
        print("addPoints")
 
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
                

        LittleD_A = np.long(self.IndexDispA.sum())
        # self.XdispA = self.X[self.i_disp]
        # self.YdispA = self.Y[self.i_disp]
        # self.ZdispA = self.Z[self.i_disp]
        
        if (LittleD_A > 0) :
            self.XdispA = np.zeros(LittleD_A)
            self.YdispA = np.zeros(LittleD_A)
            self.ZdispA = np.zeros(LittleD_A)
            FT.crunch_all(self.X, self.Y, self.Z, self.XdispA, self.YdispA, self.ZdispA, self.IndexDispA, self.BigD, LittleD_A)
        
                
        self.IndexDisp[:] = False

        self.IndexDisp[self.i_disp] = True  
        
        self.IndexComb = np.logical_and(self.IndexArray, self.IndexDispA)    
        LittleD_B = np.long(self.IndexComb.sum())
        self.LittleD_B = LittleD_B

        
        self.XdispB = np.zeros(LittleD_B)
        self.YdispB = np.zeros(LittleD_B)
        self.ZdispB = np.zeros(LittleD_B)

        # if (LittleD_B > 0) :
            # FT.crunch_all(self.X, self.Y, self.Z, self.XdispB, self.YdispB, self.ZdispB, self.IndexComb, self.BigD, LittleD_B)

        LittleD_A = len(self.XdispA)
        #print("Number of points in index array", self.IndexArray.sum())
        data = np.empty((LittleD_A, 3), dtype = np.float64)

        
        data[:,0] = self.XdispA
        data[:,1] = self.YdispA
        data[:,2] = self.ZdispA

        numPoints = LittleD_A
     
        self.vtkPoints.Allocate(numPoints)
        table = numpy_support.numpy_to_vtk(data, deep = 2, array_type=None)

        self.vtkPoints.SetData(table)
        self.vtkActor.GetMapper().SetColorModeToDefault()
        self.vtkActor.GetProperty().SetColor(0.0, 0.0, 1.0)
        self.vtkActor.GetProperty().SetPointSize(self.PtDisplaySize)
        self.vtkPoints.Modified()
        self.vtkPoints.Modified()
        self.vtkDepth.Modified()
        
        print("size of points = ", self.vtkCells.GetSize())
        

                                              
    def CreatePointCloud(self):   
        print("CreatePointCloud_00")
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

        
class TwoDPlot():

    def __init__(self, view, HchartA,VchartA, HorTitle,VerTitle, GraphMinH, GraphMinV, GraphRange, GraphStep):
       #print("2d plot init_00")
        self.view = view
        self.HchartA = HchartA
        self.VchartA = VchartA
        self.chart = vtk.vtkChartXY()
        self.chart.SetShowLegend(False)
        self.tableA = vtk.vtkTable()         
        self.GraphMinH = GraphMinH
        self.GraphMinV = GraphMinV
        self.GraphRange = GraphRange
        self.GraphStep = GraphStep
        self.PtDisplaySize = 2.0

         
        numPoints = len(HchartA)
        
        self.tableA.SetNumberOfRows(numPoints)

        self.HarrA = numpy_support.numpy_to_vtk(self.HchartA, deep = 0, array_type=None)
        self.VarrA = numpy_support.numpy_to_vtk(self.VchartA, deep = 0, array_type=None)
       #print("2DPI_03:")
        self.HarrA.SetName('Horizontal')
        self.VarrA.SetName('Vertical')
        self.tableA.AddColumn(self.HarrA)
        self.tableA.AddColumn(self.VarrA)
        self.chart.ClearPlots()

        self.points = self.chart.AddPlot(vtk.vtkChart.POINTS)
        self.points.SetInputData(self.tableA, 0, 1)
        #   self.points.SetColor(0, 0, 255, 255)
        self.points.SetWidth(0.1)
        self.points.SetMarkerSize(self.PtDisplaySize)
        self.points.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)
        
        #view = vtk.vtkContextView()
        self.view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        self.view.GetRenderWindow().SetSize(400, 300)
        self.view.GetScene().ClearItems()
        self.view.GetScene().AddItem(self.chart)

          
        self.ResetScale(self.GraphMinH, self.GraphMinV, self.GraphRange, self.GraphStep)
    
    def ResetScale(self, GraphMinH, GraphMinV, GraphRange, GraphStep) :
        print("ResetScale")
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
        

    def Update(self, HchartA, VchartA, HchartB,VchartB):
        print("Update 2 d table")
         
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

        
        points = self.chart.AddPlot(vtk.vtkChart.POINTS)   
        points.SetInputData(self.tableA, 0, 1)
        print("Set Color")
        points.SetColor(0, 0, 255, 255)
        points.SetWidth(0.1)
        points.SetMarkerSize(self.PtDisplaySize)
        points.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)
        LittleD_B = len(self.HchartB)

              
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
        


class MyInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):

 
    def __init__(self,parent, Display, pointCloud):
    #(self.renderWindowInteractor, self, pcl, self.SI, self.sphere))
        print ("MIS_01")
        self.AddObserver("LeftButtonPressEvent",self.leftButtonPressEvent)
        self.AddObserver("LeftButtonReleaseEvent",self.leftButtonReleaseEvent)
        self.AddObserver("MiddleButtonPressEvent",self.middleButtonPressEvent)
        self.AddObserver("MiddleButtonReleaseEvent",self.middleButtonReleaseEvent)
        self.AddObserver("MouseMoveEvent",self.mouseMoveEvent)
        self.AddObserver("KeyPressEvent", self.keyPress)
        self.RWI = parent
        self.pcld = pointCloud


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
        print("Middle Button pressed")
        #print(obj)
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
         self.OnMouseMove()
             
  
    def keyPress(self, obj,event):

        key = self.RWI.GetKeySym()
        print("key = ", key)
 

class Refreshing():
    def __init__(self, renderer, Xpcl, Ypcl, Zpcl):
        self.renderer = renderer
        #(self.Xpcl,self.Ypcl,self.Zpcl) = PCF.PandaLoad(inputfile)
        self.Xpcl = Xpcl
        self.Ypcl = Ypcl
        self.Zpcl = Zpcl
        
        self.counts = 0
        
    def refresh(self, iren, event):

        print("counts = ", self.counts)
        
        if self.counts == 0:
            print("counts = 0")
            
            pointCloud = VtkPointCloud(self.Xpcl, self.Ypcl, self.Zpcl)
            self.renderer.AddActor(pointCloud.vtkActor)
            pointCloud.clearPoints()

            for k in range(len(self.Xpcl)):
                point = ([self.Xpcl[k], self.Ypcl[k], self.Zpcl[k]])
                pointCloud.addPoint(point)

            #print(pointCloud.vtkPoints[100])
        
            print(dir(pointCloud))
            #print("Point 100 = ", pointCloud[100])

        iren.GetRenderWindow().Render()
        # if self.iterations == 30:
        if self.counts == 0:
            self.renderer.ResetCamera()
        self.counts += 1     

         
class FourPanelDisplay():
    
    def __init__(self): 
        inputfile = 'JustFireFace.xyz'

        self.root = Tk()
  
        self.root.title("Point cloud console")
        (Xpcl,Ypcl,Zpcl) = PCF.PandaLoad(inputfile)
       
        pcl = VtkPointCloud(Xpcl, Ypcl, Zpcl)

        pcl.CreatePointCloud()  
        self.pcl = pcl
         
        hor = 800
        ver = 600

        xmins=[0,.5,0,.5]
        xmaxs=[0.5,1,0.5,1]
        ymins=[0,0,.5,.5]
        ymaxs=[0.5,0.5,1,1] 

        self.renderer = vtk.vtkRenderer()
        self.renderer.AddActor(pcl.vtkActor)

        self.renderer.SetBackground(1.0, 1.0, 1.0)
        self.renderer.ResetCamera()
        self.renderer.GetActiveCamera().ParallelProjectionOn() 
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

        origin = vtk.vtkSphereSource()
        origin.SetCenter(0,0,0)
        origin.SetRadius(5)
        origin.SetPhiResolution(20)
        origin.SetThetaResolution(20)

        origin_mapper = vtk.vtkPolyDataMapper()
        origin_mapper.SetInputConnection(origin.GetOutputPort())
        origin_actor = vtk.vtkActor()
        #origin_actor.GetProperty().SetColor(0.0, 0.0, 1.0)
        
        origin_actor.SetMapper(origin_mapper)
        self.renderer.AddActor(origin_actor)       
        self.renderer.ResetCamera()
         
        self.renderWindow = vtk.vtkRenderWindow()
        #renderWindow.SetWindowName("Interactive")
        self.renderWindow.AddRenderer(self.renderer)
        self.renderer.SetViewport(xmins[0],ymins[0],xmaxs[0],ymaxs[0])
        
        self.renderWindowInteractor = vtk.vtkRenderWindowInteractor()
        self.renderWindowInteractor.SetRenderWindow(self.renderWindow)
        self.renderWindowInteractor.SetInteractorStyle(MyInteractorStyle(self.renderWindowInteractor, self, pcl))

        self.renderWindow.SetSize(hor,ver)

        self.XcontextView = vtk.vtkContextView()
        self.YcontextView = vtk.vtkContextView()
        self.ZcontextView = vtk.vtkContextView()
        self.TDP_Init = False

        self.ResetGraphExtents()

        self.TDP_X = TwoDPlot(self.XcontextView, pcl.Y,pcl.Z,"Y","Z", self.pcl.GraphMinY, self.pcl.GraphMinZ, 
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Y = TwoDPlot(self.YcontextView, pcl.X,pcl.Z,"X","Z", self.pcl.GraphMinX, self.pcl.GraphMinZ,
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Z = TwoDPlot(self.ZcontextView, pcl.X,pcl.Y,"X","Y", self.pcl.GraphMinX, self.pcl.GraphMinY, 
            self.pcl.GraphRange, self.pcl.GraphStep)
        self.TDP_Init = True

        
        self.Xrenderer = self.TDP_X.view.GetRenderer()
        self.Yrenderer = self.TDP_Y.view.GetRenderer()
        self.Zrenderer = self.TDP_Z.view.GetRenderer()
        
        self.renderWindow.AddRenderer(self.Xrenderer)
        self.Xrenderer.SetViewport(xmins[1],ymins[1],xmaxs[1],ymaxs[1])
        self.renderWindow.AddRenderer(self.Yrenderer)
        self.Yrenderer.SetViewport(xmins[2],ymins[2],xmaxs[2],ymaxs[2])
        self.renderWindow.AddRenderer(self.Zrenderer)
        self.Zrenderer.SetViewport(xmins[3],ymins[3],xmaxs[3],ymaxs[3])

        self.renderWindow.Render()

        self.renderWindow.SetWindowName("Interactive")

        self.pcl.clearPoints()
        temp = input('1604\n')
        self.pcl.addPoints()

        temp = input('3000\n')
        #time.sleep(1)
        self.Rebuild()
        temp = input('3200\n')


        timerId = self.renderWindowInteractor.CreateRepeatingTimer(10)    

        #(self.Xpcl,self.Ypcl,self.Zpcl) = PCF.PandaLoad(inputfile)
        refreshing = Refreshing(self.renderer, Xpcl,Ypcl,Zpcl)
        self.renderWindowInteractor.AddObserver('TimerEvent', refreshing.refresh)
        refreshing.timerId = timerId
        
        
        self.root.mainloop()
        ##########tkinter window

        
    def Rebuild(self) :
        self.pcl.clearPoints()
        self.pcl.addPoints()
        self.TDP_X.Update(self.pcl.YdispA, self.pcl.ZdispA, self.pcl.YdispB, self.pcl.ZdispB)
        self.TDP_Y.Update(self.pcl.XdispA, self.pcl.ZdispA, self.pcl.XdispB, self.pcl.ZdispB)
        self.TDP_Z.Update(self.pcl.XdispA, self.pcl.YdispA, self.pcl.XdispB, self.pcl.YdispB)
        print("Render")
        self.renderWindow.Render() 
        
        
    def ResetGraphExtents(self) :
        #self.renderWindowInteractor.Start()
        X_st_min = self.pcl.X.min()
        X_st_max = self.pcl.X.max()
        Y_st_min = self.pcl.Y.min()
        Y_st_max = self.pcl.Y.max()
        Z_st_min = self.pcl.Z.min()
        Z_st_max = self.pcl.Z.max()


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
    
    def calculate(self, *args):
        print("calculate")
        print("self.pcl", self.pcl)
        #input("1000")
        #self.renderWindow.Render()
        self.Rebuild()
            
            
        
            
            
            
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
     
    return(cubeAxesActor)
    

def DisplayData() :
    import sys 
    #print("DD_01")

    FPD_out = FourPanelDisplay()

def main():
    
    DisplayData()

if __name__ == '__main__':
    main()




