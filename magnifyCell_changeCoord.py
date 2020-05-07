#!/usr/bin/python

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
from numpy import pi,arctan
import sys

class magnifyCell_changeCoord():
  """magnify a cell and change its coordinate types between direct and cardesion"""
  def __init__(self,filename):
    self.filename=filename
    self.latt_para=1.0
    self.w_coord=1
    self.sys_name="" 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.N=[1,1,1] #default, will read from structural file
    self.n_unit=[]
    self.mag_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    self.mag_atoms_pos=[] 
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1: 
          if len(ll)==5:
            self.sys_name,self.w_coord,self.N=ll[0],int(ll[1]),[int(ll[2]),int(ll[3]),int(ll[4])]
          elif len(ll)==4:
            self.sys_name,self.N=ll[0],[int(ll[1]),int(ll[2]),int(ll[3])]
          elif len(ll)==2:
            self.sys_name,self.w_coord=ll[0],int(ll[1])
          else:
            self.sys_name=ll
        if count==2: self.latt_para=float(ll[0])
        if count>2 and count<6:
          self.coord[count-3]=np.array([float(ll[0]),float(ll[1]),float(ll[2])])
        if count==6:
          for i in ll:
            self.n_unit.append(int(i))
        if 'Cartesian' in line:
           self.coord_type='Cartesian'
           break
        if 'Direct' in line:       
           self.coord_type='Direct'
           break
        count +=1
      for line in in_file:
        if line != '\n':
          ll = line.split()
          ll[0],ll[1],ll[2]=float(ll[0]),float(ll[1]),float(ll[2])
          self.atoms_pos.append(ll[0:3])
          self.mag_atoms_pos.append([ll[0:3]])
  def magnify_cell(self):
    self.read_data()
    k1=0
    for i in self.coord:
      self.mag_coord[k1]=np.array([i[0]*self.N[k1],i[1]*self.N[k1],self.N[k1]*i[2]])*self.latt_para
      k1 +=1
      if k1==1: exit
    n0=len(self.atoms_pos)
    if n0 != sum(self.n_unit): print "Some atomic positions are missing!"
    if self.coord_type == 'Direct':
      for ix in range(0,self.N[0]):
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):
            for k in range(0,n0):
               mag_atom = self.atoms_pos[k]+ix*np.array([1,0,0])+iy*np.array([0,1,0])+iz*np.array([0,0,1])  
               mag_atom = mag_atom/self.N #/[2.0,2.0,1.0]
               mag_atom = np.dot(self.mag_coord.transpose(),mag_atom)
               self.mag_atoms_pos[k].append(mag_atom)
    elif self.coord_type == 'Cartesian':
      for ix in range(0,self.N[0]):            
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):          
            for k in range(0,n0):         
               mag_atom =self.atoms_pos[k]+np.dot(self.coord.transpose(),np.array([ix,iy,iz]))  
               self.mag_atoms_pos[k].append(mag_atom*self.latt_para)
  def print_disl(self):
    self.magnify_cell()
    print self.sys_name
    print 1.0 #self.latt_para
    for i in range(0,3):
      print format(self.mag_coord[i,0],"03f"),"	",format(self.mag_coord[i,1],"03f"),"	",\
            format(self.mag_coord[i,2],"03f")
    for i in self.N[0]*self.N[1]*self.N[2]*np.asarray(self.n_unit):
      print i,
    if self.w_coord == 1:
      print "\nCartesian" #self.coord_type
      for k in range(0,len(self.atoms_pos)):
        self.mag_atoms_pos[k].pop(0)
        for i in self.mag_atoms_pos[k]:
          print format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f") 
    elif self.w_coord == 0:
      print "\nDirect" #self.coord_type
      for k in range(0,len(self.atoms_pos)):
        self.mag_atoms_pos[k].pop(0)
        for i in self.mag_atoms_pos[k]:
          i=np.dot(np.linalg.inv(self.mag_coord.transpose()),i)
          print format(i[0],"03f"),"    ",format(i[1], "03f")," ",format(i[2],"03f")
if __name__=="__main__":
  disl1=magnifyCell_changeCoord(sys.argv[1])#"unit_cell")
  disl1.print_disl()
