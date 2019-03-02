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
from numpy import pi,arctan
import sys

class gen_disl():
  """generate a dislocation"""
  def __init__(self,filename1, filename2):
    self.filename1=filename1
    self.filename2=filename2
    self.latt_para=1.0
    self.w_coord=1
    self.sys_name="" 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.N=[1,1,1] #default, will read from structural file
    self.n_unit=[]
    self.GSFE_requirements=[]
  def read_data(self):
    with open(self.filename1,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1: self.sys_name=ll[0]
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
    with open(self.filename2,'r') as in_file2:
      count=1
      for line in in_file2:
        ll=line.split(":")
        if (count>2 and count <7):
          lll=ll[1].split(",")
          for i in range(0,len(lll)):
            lll[i]=float(lll[i])
          self.GSFE_requirements.append(lll)
        else:
          self.GSFE_requirements.append(int(ll[1]))
        count +=1
  def new_coord(self,i,j):
    new_coord=[]
    new_coord.append(self.coord[0])
    new_coord.append(self.coord[1])
    stepX=np.asarray(self.GSFE_requirements[2])*float(self.GSFE_requirements[4][0])/(int(self.GSFE_requirements[5][0])-1)
    stepY=np.asarray(self.GSFE_requirements[3])*float(self.GSFE_requirements[4][1])/(int(self.GSFE_requirements[5][1])-1)
    X=i*stepX[0]+j*stepY[0]+self.coord[2][0]
    Y=i*stepX[1]+j*stepY[1]+self.coord[2][1]
    Z=i*stepX[2]+j*stepY[2]+self.coord[2][2]
    new_coord.append([X,Y,Z])
    return new_coord
  def new_atom_pos(self,i,j):
    new_atom_pos=[]
    pos_cut=0.499
    for k in range(0,len(self.atoms_pos)):
      if self.atoms_pos[k][2]<self.coord[2][2]*pos_cut:
        new_atom_pos.append(self.atoms_pos[k])
      else:
        new_pos=[]
        stepX=1.0*np.asarray(self.GSFE_requirements[2])*float(self.GSFE_requirements[4][0])/(int(self.GSFE_requirements[5][0])-1)
        stepY=1.0*np.asarray(self.GSFE_requirements[3])*float(self.GSFE_requirements[4][1])/(int(self.GSFE_requirements[5][1])-1)
        X=i*stepX[0]+j*stepY[0]+self.atoms_pos[k][0]  
        Y=i*stepX[1]+j*stepY[1]+self.atoms_pos[k][1]  
        Z=i*stepX[2]+j*stepY[2]+self.atoms_pos[k][2]  
        new_atom_pos.append([X,Y,Z])
    return new_atom_pos
  def write_file(self,new_coord,new_atoms_pos,wfile):
    wfile.write(self.sys_name+'\n')
    wfile.write(str(self.latt_para)+'\n') #self.latt_para
    for i in range(0,3):
      wfile.write(str(format(new_coord[i][0],"03f"))+"	"+str(format(new_coord[i][1],"03f"))+"	"+\
            str(format(new_coord[i][2],"03f"))+'\n')
    for i in np.asarray(self.n_unit):
      wfile.write(str(i)+" "),
    if self.w_coord == 1:
      wfile.write("\nCartesian\n") #self.coord_type
      for i in new_atoms_pos:
          if self.GSFE_requirements[6] == 0:
            wfile.write(str(format(i[0],"03f"))+"	"+str(format(i[1], "03f"))+"	"+str(format(i[2],"03f"))+'\n')
          elif self.GSFE_requirements[6] == 1:
            wfile.write(str(format(i[0],"03f"))+"	"+str(format(i[1], "03f"))+"	"+str(format(i[2],"03f"))+" F F T"+'\n')
          else: print("Wrong label! Please only put number 0 or 1!")
    elif self.w_coord == 0:
      wfile.write("\nDirect\n") #self.coord_type
      for i in new_atoms_pos:
          i=np.dot(np.linalg.inv(self.coord.transpose()),i)
          if self.GSFE_requirements[5] == 0:               
            wfile.write(str(format(i[0],"03f"))+"	"+str(format(i[1], "03f"))+"	"+str(format(i[2],"03f"))+'\n') 
          elif self.GSFE_requirements[5] == 1:             
            wfile.write(str(format(i[0],"03f"))+"	"+str(format(i[1], "03f"))+"	"+str(format(i[2],"03f"))+" F F T"+'\n')
          else: print("Wrong label! Please only put number 0 or 1!")
  def print_disl(self):
    self.read_data()
    for i in range(0,int(self.GSFE_requirements[5][0])):
      if int(self.GSFE_requirements[0]) == 1:
        if int(self.GSFE_requirements[1]) == 1:
          new_coord=self.new_coord(i,0)
          with open("POSCAR"+str(i),'w') as wfile:
            self.write_file(new_coord,self.atoms_pos,wfile)
        elif int(self.GSFE_requirements[1]) == 2:
          for j in range(0,int(self.GSFE_requirements[5][1])):   
            new_coord=self.new_coord(i,j)           
            with open("POSCAR_"+str(i)+"_"+str(j),'w') as wfile:
              self.write_file(new_coord,self.atoms_pos,wfile)
        else: print("Please put a right number for the dimensionality!")
      elif int(self.GSFE_requirements[0]) == 2:
        if int(self.GSFE_requirements[1]) == 1:
          new_atom_pos=self.new_atom_pos(i,i)
          with open("POSCAR"+str(i),'w') as wfile:
            self.write_file(self.coord,new_atom_pos,wfile)
        elif int(self.GSFE_requirements[1]) == 2:
          for j in range(0,int(self.GSFE_requirements[5][1])):   
            new_atom_pos=self.new_atom_pos(i,j)           
            with open("POSCAR_"+str(i)+"_"+str(j),'w') as wfile:
              self.write_file(self.coord,new_atom_pos,wfile)
        else: print("Please put a right number for the dimensionality!")
      else: print("Please pick 1 or 2 for the number of GSFs!")
disl1=gen_disl(sys.argv[1],sys.argv[2])#"unit_cell")
#disl1.read_data()
disl1.print_disl()
