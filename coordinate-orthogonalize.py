#!/usr/bin/python

import numpy as np
from numpy import pi,arctan
from numpy.linalg import inv
import sys
import random

class gen_disl():
  """generate a dislocation"""
  def __init__(self,filename,filename_coord):
    self.filename=filename
    self.filename_coord=filename_coord
    self.latt_para=1.0
    self.w_coord=1
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
          if len(ll)==2: self.sys_name,self.w_coord=ll[0],int(ll[1])
          else: self.sys_name=ll[0]
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
         atom_in_direct.append(self.in_box(np.dot(inv(self.new_coord.transpose()),np.dot(self.coord.transpose(),i))))
     if self.coord_type=='Cartesian':
       for i in self.atoms_pos:
         atom_in_direct.append(self.in_box(np.dot(inv(self.new_coord.transpose()),i)))    
     if self.w_coord==0:
       self.write_atom_pos=atom_in_direct
     elif self.w_coord==1:
       for i in atom_in_direct:
         i=np.dot(self.new_coord.transpose(),i)
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

disl1=gen_disl(sys.argv[1],sys.argv[2])#"unit_cell")
disl1.print_new_supercell()
#print disl1.coord
#print disl1.new_coord
