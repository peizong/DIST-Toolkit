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
import random

class gen_disl():
  """generate a dislocation"""
  def __init__(self,filename,Rc_cutoff):
    self.filename=filename
    self.Rc=float(Rc_cutoff)
    self.latt_para=1.0
    self.w_coord=-1
    self.sys_name="" 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.N=[10,5] #default, will read from structural file
    self.n_unit=[]
    self.mag_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    self.randStruct=[] 
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1: 
          #self.sys_name,self.w_coord,self.N=ll[0],int(ll[1]),[int(ll[2]),int(ll[3])]
          self.sys_name=line
          #------remove comment signs for these lines later------------
          #if len(ll)==1: self.sys_name=ll[0]
          #elif len(ll)==2: self.sys_name,self.w_coord=ll[0],int(ll[1])
          #else: 
          #  print("Error from the first line of input file!")
          #  return None
          #  break
          #--------------------end-------------------------------------
        #  for i in range(0,len(self.N)):
        #    self.randStruct.append([[0,0,0]])
        if count==2: self.latt_para=float(ll[0])
        if count>2 and count<6:
          self.coord[count-3]=np.array([float(ll[0]),float(ll[1]),float(ll[2])])
        if count==6:
          for i in ll:
            self.n_unit.append(int(i))
        if 'Cartesian' in line:
           self.coord_type=1 #'Cartesian'
           if self.w_coord==-1: self.w_coord=1
           break
        if 'Direct' in line:       
           self.coord_type=0 #'Direct'
           if self.w_coord==-1: self.w_coord=0
           break
        count +=1
      for line in in_file:
        if line != '\n':
          ll = line.split()
          ll[0],ll[1],ll[2]=float(ll[0]),float(ll[1]),float(ll[2])
          self.atoms_pos.append(ll[0:3])
  def gen_random(self):
    self.read_data()
    #randomArray=random.sample(range(1,(self.N[0]+1)),self.N[1])
    centers=np.array([[1.0/3,1.0/3],[2.0/3,2.0/3]])
    LxLy=np.array([[self.coord[0][0],self.coord[0][1]],[self.coord[1][0],self.coord[1][1]]])
    newCenters=np.dot(centers,LxLy)
    Rc=self.Rc #16.0
    selected_ID,unselected_ID1, unselected_ID2,record_ID,randomized_ID =[],[],[],[],[]
    for i in range(0,len(self.atoms_pos)):
      dist1=np.linalg.norm(self.atoms_pos[i][0:2]-newCenters[0])
      dist2=np.linalg.norm(self.atoms_pos[i][0:2]-newCenters[1])
      if (dist1 <= Rc) or (dist2 <= Rc):
        selected_ID.append(i)
      else:
        if (i<len(self.atoms_pos)/2.0):
          unselected_ID1.append(i)
        else:
          unselected_ID2.append(i)
    randomized_order=random.sample(range(0,len(selected_ID)),len(selected_ID))
    for i in randomized_order:
      randomized_ID.append(selected_ID[i])
    #print randomized_ID
    #randomArray=random.sample(range(0,np.sum(self.n_unit)),np.sum(self.n_unit))
    randomArray=np.concatenate((unselected_ID1,randomized_ID,unselected_ID2))
    for i in range(0,len(self.atoms_pos)):
      self.randStruct.append(self.atoms_pos[randomArray[i]])
    #for i in range(0,len(self.atoms_pos)):
    #  if i in randomArray:
    #    self.randStruct[0].append(self.atoms_pos[i])
    #  else:
    #    self.randStruct[1].append(self.atoms_pos[i]) 
  def print_random(self):
    self.gen_random()
    print(self.sys_name), #remove "," later
    print(self.latt_para)
    for i in range(0,3):
      print(format(self.coord[i,0],"03f"),"	",format(self.coord[i,1],"03f"),"	",\
            format(self.coord[i,2],"03f"))
    #for i in range(1,len(self.N)):
    #  print self.N[i], 
    #print self.N[0]-sum(self.N[1:len(self.N)])
    for i in range(0,len(self.n_unit)-1):
      print(self.n_unit[i]),
    print(self.n_unit[len(self.n_unit)-1])
    if self.w_coord == 1:
      print("Cartesian") #self.coord_type
      for i in self.randStruct:
        if self.w_coord==self.coord_type:
          print(format(i[0],"03f"),"   ",format(i[1], "03f")," ",format(i[2],"03f"))
        else:
          i=np.dot(self.coord.transpose(),i)
          print(format(i[0],"03f"),"   ",format(i[1], "03f")," ",format(i[2],"03f")) 
     # for k in range(0,len(self.N)):
     #   self.randStruct[k].pop(0)
     #   for i in self.randStruct[k]:
     #     print format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f") 
    elif self.w_coord == 0:
      print("Direct") #self.coord_type
      for i in self.randStruct:
        if self.w_coord==self.coord_type:
          print(format(i[0],"03f"),"    ",format(i[1], "03f")," ",format(i[2],"03f"))
        else:
          i=np.dot(np.linalg.inv(self.coord.transpose()),i)
          print(format(i[0],"03f"),"    ",format(i[1], "03f")," ",format(i[2],"03f"))
     # for k in range(0,len(self.N)):
     #   self.randStruct[k].pop(0)
     #   for i in self.randStruct[k]:
     #     i=np.dot(np.linalg.inv(self.coord.transpose()),i)
     #     print format(i[0],"03f"),"    ",format(i[1], "03f")," ",format(i[2],"03f")

disl1=gen_disl(sys.argv[1],sys.argv[2])#"unit_cell")
disl1.print_random()
