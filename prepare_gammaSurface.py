#!/nfs/apps/Compilers/Python/Anaconda/2.7/bin/python
##!/usr/bin/python

############################################################################### 
#                                                                          * F# 
#               DIST: A DIslocation-Simulation Toolkit                     2 R# 
# GNU License - Author: Zongrui Pei                       2015-06-10       0 A# 
# Version 1.0                                                              1 N# 
#                                                                          5 K# 
# Syntax:                                                                  0 F# 
# Please find the syntx in the howto.dat of the examples folder            6 U# 
# and the CPC paper: Zongrui Pei, DIST: A DIslocation-Simulation Toolkit,  1 R# 
# Computer Physics Communications 233(2018)44-50.                          0 T# 
#                                                                          * *# 
###############################################################################

import numpy as np
from numpy import pi,arctan,sin,cos
import sys
from scipy.optimize import leastsq
import matplotlib.pyplot as plt

class fit_gammaSurface():
  """fit gamma surface using the five point method"""
  def __init__(self,filename,task):
    self.filename=filename
    self.task=task
    self.data=[]
    self.triGeomFuns6=np.array([0,0,0,0,0,0]) #initial value for 5 point method
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count>1:
          self.data.append(np.array([float(ll[0]),float(ll[1]),float(ll[2])]))
        count +=1
  def cal_triGeomFuns6(self,X,Y):
    R0_1=1.0
    R1_cos=cos(2*pi*(X-Y))+cos(2*pi*(X+Y))+cos(4*pi*Y)
    R2_cos=cos(2*pi*(X-3*Y))+cos(2*pi*(X+3*Y))+cos(4*pi*X)
    R3_cos=cos(4*pi*(X-Y))+cos(4*pi*(X+Y))+cos(8*pi*Y)
    I1_sin=sin(2*pi*(X-Y))-sin(2*pi*(X+Y))+sin(4*pi*Y)
    I2_sin=sin(4*pi*(X-Y))-sin(4*pi*(X+Y))+sin(8*pi*Y)
    return np.array([R0_1,R1_cos,R2_cos,R3_cos,I1_sin,I2_sin])
  def func_leastsq(self,params,xdata,ydata):
    return (ydata-np.dot(xdata,params))
  def cal_coeff(self):
    M,E=[],np.zeros(len(self.data))
    for i in range(0,len(self.data)):
      M.append(self.cal_triGeomFuns6(self.data[i][0],self.data[i][1]))
      E[i]=self.data[i][2]
    if len(self.data)==6:
      self.triGeomFuns6=np.linalg.inv(M)*E
    elif len(self.data)>6:
      xdata,ydata=M,E
      results=leastsq(self.func_leastsq,self.triGeomFuns6,args=(xdata,ydata))
      self.triGeomFuns6=results[0]
    else: 
      return 0
      print("Insufficient data for fitting!")
  def check_fitting(self):
    fitted=[]
    ydata=[]
    for i in self.data:
      the_sin_coss=self.cal_triGeomFuns6(i[0],i[1])
      fitted.append(np.dot(the_sin_coss,self.triGeomFuns6))    
      ydata.append(i[2])
    plt.subplot(121)
    xmax,xmin=np.max(ydata),np.min(ydata)
    plt.scatter(ydata,fitted)
    plt.plot([xmin,xmax],[xmin,xmax],'r-')
    plt.subplot(122)
    numX=np.linspace(1,len(self.data),len(self.data))
    error=(np.asarray(ydata)-np.asarray(fitted)) #/np.asarray(ydata)
    plt.plot(numX,error)
    plt.show()
  def plot_gamma_contour(self):
    sampling_point=20
    x=np.linspace(0,1,sampling_point)
    y=x
    X,Y=np.meshgrid(x,y)
    Z=[]
    for i in range(0,sampling_point):
      Zj=[]
      for j in range(0,sampling_point):
        Zij=np.dot(self.cal_triGeomFuns6(X[i][j],Y[i][j]),self.triGeomFuns6)
        Zj.append(Zij)
      Z.append(Zj)
    #Z=np.dot(self.cal_triGeomFuns6(X,Y),self.triGeomFuns6)
    plt.contourf(X,Y,Z,40)
    plt.xlabel("x-axis")
    plt.ylabel("y-axis")
    plt.colorbar()
    plt.show()
  def cal_shear_modulus(self,direction):
    c1,c2=self.triGeomFuns6[1],self.triGeomFuns6[2]
    a=3.52e-10 # fcc
    # a=a*2**0.5 #hcp
    b=0.5**0.5*a
    d=(2.0/3)**0.5*b
    pseudoG=-8*pi**2*(c1+c2)*d/b**2
    return pseudoG/1e9
  def post_processing(self):
    self.read_data()
    self.cal_coeff()
    if (self.task=='coeff'):
      print(self.triGeomFuns6)
    if (self.task=='check'):
      self.check_fitting()
    if (self.task=='contour'):
      self.plot_gamma_contour()
    if (self.task=='G'):
      print(self.cal_shear_modulus(1))
if __name__=="__main__":
  gS=fit_gammaSurface(sys.argv[1],sys.argv[2])
  gS.post_processing()
