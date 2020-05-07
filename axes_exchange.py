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
from numpy.linalg import inv
import sys
import random

class axes_exchange():
  """exchange any two of the three primitive vectors"""
  def __init__(self,filename,filename_coord):
    self.filename=filename
    self.filename_coord=filename_coord
    self.latt_para=1.0
    self.w_coord=100 # the default means coord type is not set up yet
    self.sys_name="" 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.n_unit=[]
    self.mag_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    self.write_atom_pos=[] 
    self.new_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1: 
          if len(ll)==1:
            self.sys_name=ll[0]
          elif len(ll)==2:
            self.sys_name,self.w_coord=ll[0],int(ll[1])
          else: "The setting of first line in POSCAR is wrong! Too much information!"
        if count==2: self.latt_para=float(ll[0])
        if count>2 and count<6:
          self.coord[count-3]=np.array([float(ll[0]),float(ll[1]),float(ll[2])])
        if count==6:
          for i in ll:
            self.n_unit.append(int(i))
        if 'Cartesian' in line:
           self.coord_type='Cartesian'
           if self.w_coord==100: self.w_coord=1
           break
        if 'Direct' in line:       
           self.coord_type='Direct'
           if self.w_coord==100: self.w_coord=0
           break
        count +=1
      for line in in_file:
        if line != '\n':
          ll = line.split()
          ll[0],ll[1],ll[2]=float(ll[0]),float(ll[1]),float(ll[2])
          self.atoms_pos.append(ll[0:3])
  def read_new_coord(self):
    with open(self.filename_coord,'r') as in_file:
      count=0
      for line in in_file:
        if line != '\n':
          ll=line.split()
          ll[0],ll[1],ll[2]=float(ll[0]),float(ll[1]),float(ll[2])
          self.new_coord[count]=np.array([ll[0],ll[1],ll[2]])
        count +=1
  def in_box(self,atom_pos):
     for i in range(0,3):
       while(atom_pos[i]<0.0): atom_pos[i]=atom_pos[i]+1.0
       while( atom_pos[i]>=1.0): atom_pos[i]=atom_pos[i]-1.0
     return atom_pos
  def orthogonalize(self):
     self.read_data()
     self.read_new_coord()
     atom_in_direct=[]
     if self.coord_type=='Direct':
       for i in self.atoms_pos:
         atom_in_direct.append(self.in_box(inv(self.new_coord.transpose()).dot(self.coord.transpose().dot(i))))
     if self.coord_type=='Cartesian':
       for i in self.atoms_pos:
         atom_in_direct.append(self.in_box(inv(self.new_coord.transpose()).dot(i)))    
     if self.w_coord==0:
       self.write_atom_pos=atom_in_direct
     elif self.w_coord==1:
       for i in atom_in_direct:
         i=self.new_coord.transpose().dot(i)
         self.write_atom_pos.append(i)
  def print_new_supercell(self):
    self.orthogonalize()
    print self.sys_name
    print self.latt_para
    for i in range(0,3):
      print format(self.new_coord[i,0],"03f"),"	",format(self.new_coord[i,1],"03f"),"	",\
            format(self.new_coord[i,2],"03f")
    for i in self.n_unit:
      print i,
    if self.w_coord == 1:
      print "\nCartesian" #self.coord_type
    elif self.w_coord == 0:
      print "\nDirect" 
    for i in self.write_atom_pos:
      print format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f") 
if __name__=="__main__":
  dist1=axes_exchange(sys.argv[1],sys.argv[2])#"unit_cell")
  dist1.print_new_supercell()
