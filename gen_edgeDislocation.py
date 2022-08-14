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
from __future__ import print_function
import numpy as np
from numpy import sign, pi,arctan,arctan2,log
import sys

class gen_disl():
  """generate a dislocation"""
  def __init__(self,filename):
    self.filename=filename
    self.latt_para=1.0
    self.b=1.0
    self.w=1.0 #0.1
    self.d=12 #defult 0
    self.sys_name=""
    self.disl_along_axis=1 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.N=[1,1,1] #default, will read from structural file
    self.n_unit=[]
    self.mag_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    self.mag_atoms_pos=[] 
    self.disl_center=np.array([0,0])
    self.disl_atoms_pos=[]
    self.num_disl=1 #[]
    #self.magnify_cell()
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1:
          if len(ll)==2:
            self.sys_name,self.b=ll[0],float(ll[1])
          if len(ll)==4:
            self.sys_name,self.b,self.disl_along_axis,self.num_disl=ll[0],float(ll[1]),int(ll[2]),int(ll[3])
          if len(ll)==6:
            self.sys_name,self.b,self.disl_along_axis,self.N=ll[0],float(ll[1]),int(ll[2]),[int(ll[3]),int(ll[4]),int(ll[5])]
          if len(ll)==7:
            self.sys_name,self.b,self.disl_along_axis,self.num_disl,self.N=ll[0],float(ll[1]),int(ll[2]),int(ll[3]),[int(ll[4]),int(ll[5]),int(ll[6])]
          self.w=self.w*self.b 
        if count==2: self.latt_para,self.b=float(ll[0]),self.b*float(ll[0])
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
      self.mag_coord[k1]=np.array([i[0]*self.N[0],i[1]*self.N[1],self.N[2]*i[2]])*self.latt_para
      k1 +=1
      if k1==1: exit
    n0=len(self.atoms_pos)
    if n0 != sum(self.n_unit): print("Some atomic positions are missing!")
    if self.coord_type == 'Direct':
      for ix in range(0,self.N[0]):
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):
            for k in range(0,n0):
               mag_atom = self.atoms_pos[k]+ix*np.array([1,0,0])+iy*np.array([0,1,0])+iz*np.array([0,0,1])  
               mag_atom = mag_atom/self.N #/[2.0,2.0,1.0]
               mag_atom = self.mag_coord.transpose().dot(mag_atom)
               self.mag_atoms_pos[k].append(self.move_into_box(mag_atom))
    elif self.coord_type == 'Cartesian':
      for ix in range(0,self.N[0]):            
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):          
            for k in range(0,n0):         
               mag_atom =self.atoms_pos[k]+self.coord.transpose().dot(np.array([ix,iy,iz]))#ix*self.coord[0]+iy*self.coord[1]  
               self.mag_atoms_pos[k].append(self.move_into_box(mag_atom*self.latt_para)) #? /self.N
  def cal_disl_center(self):
    center=np.array([0,0])
    for i in self.coord:
      center = center +  np.array([i[0],i[1]]) #self.atoms_pos[i]
    center = center/2.0 #len(self.atoms_pos)
    self.disl_center= center
  def is_in_void_box(self,x,y):
    x,y=x-self.disl_center[0],y-self.disl_center[1]
    if self.num_disl==1:                           
      if self.disl_along_axis==1:
        #if (y>=0 and (x>-self.b/2.0) and (x<=self.b/2.0)): # "=" is added to the 2nd,3rd condition
        if (y>=0 and (x>-self.b/2.0) and (x<=self.b/2.0)):
          return True
        else: return False
      elif self.disl_along_axis==2:
        if (x>=0 and (y>-self.b/2.0) and (y<=self.b/2.0)): 
          return True
        else: return False
      else: return False
    elif self.num_disl==2:
      if self.disl_along_axis==1:
        if ((y>-self.disl_center[1]/2.0) and (y<=self.disl_center[1]/2.0) and (x>-self.b/2.0) and (x<=self.b/2.0)): 
          return True
        else: return False
      elif self.disl_along_axis==2:                              
        if ((x>-self.disl_center[0]/2.0) and (x<=self.disl_center[0]/2.0) and (y>-self.b/2.0) and (y<=self.b/2.0)): 
          return True
        else: return False
      else: return False
    else: print("Not support dislocations more than two!")
      
  def UxUz(self,x,y):
    #nu=0.3
    #e=1e-8 #avoid 1/0 error
    ux,uz=0.0,0.0
    x,y=x-self.disl_center[0],y-self.disl_center[1]
    if self.num_disl==1:
      if self.disl_along_axis==1:
        uz=0.0
        if (y>0):
          if x>0:
            #ux=self.b/(pi)*arctan(-(x-self.b)/self.w) #+self.b/2.0
            ux=self.b/(pi)*arctan(x/self.w)+self.b/2.0
            ux-=self.b
          else: ux=self.b/(pi)*arctan(x/self.w)+self.b/2.0
          #self.b/(2*pi)*(arctan2(y,x)+x*y/(x**2+y**2+e)/(2*(1-nu)))
        else: ux=0.0
      elif self.disl_along_axis==2:
        ux=0.0
        if (x>0):
          uz=self.b/pi*(self.w/(self.w**2+y**2))*sign(-y)
        else: uz=0.0
      else: ux,uz=0.0,0.0
    elif self.num_disl==2:
      if self.disl_along_axis==1:
        uz=0.0
        if ((y>-self.disl_center[1]/2.0) and (y<=self.disl_center[1]/2.0)): #(y<=self.disl_center[1]*0.5 or y>=-self.disl_center[1]*0.5):
          #ux=self.b/pi*arctan(x/self.w)*sign(-x) #+ self.b/2.
          f_ux=self.b/(2*pi)*sign(-x) #do not forget sign(-x) to adjust the displacement
          ux=self.b/2.+ f_ux*arctan((x-self.b*self.d)/self.w)+f_ux*arctan((x+self.b*self.d)/self.w)
          ux *= sign(-x)
          #for negative SFE configuration only
          f_uz=3**0.5*self.b/(6*pi)
          uz=f_uz*arctan((x+self.b*self.d)/self.w)-f_uz*arctan((x-self.b*self.d)/self.w)
          #ux *= sign(x)
          #ux += self.b/2.0
          #print("ux: ",ux)
        else: ux=0.0
      elif self.disl_along_axis==2:
        ux=0.0
        if ((x>-self.disl_center[0]/2.0) and (x<=self.disl_center[0]/2.0)):
          uz=self.b/pi*(self.w/(self.w**2+y**2))*sign(-y)
      else: ux,uz=0.0,0.0  
    else: 
      ux,uz=0.0,0.0
      print("Please supply correct number of dislocation!")
    #uz=0 #-self.b/(2*pi)*((1-2*nu)/4/(1-nu)*log(x**2+y**2+e)+(x**2-y**2)/(4*(1-nu))/(x**2+y**2+e))
    return [ux,uz] #self.b/(2*pi)*self.angle(x-self.disl_center[0],y-self.disl_center[1])

  def make_dislocation(self):
    # remove atoms and make a void box
    list_to_be_deleted=[]
    for i in range(0,len(self.atoms_pos)):
      if self.is_in_void_box(self.atoms_pos[i][0],self.atoms_pos[i][1]):
        list_to_be_deleted.append(i)
    list_to_be_deleted.sort(reverse=True)
    for i in list_to_be_deleted:
      self.atoms_pos.pop(i)
    # displace atoms to remove the void box
    for i in range(0,len(self.atoms_pos)):                                                                   
      uxz=self.UxUz(self.atoms_pos[i][0], self.atoms_pos[i][1])
      self.atoms_pos[i][0] +=uxz[0] # remove these two #s
      self.atoms_pos[i][2] +=uxz[1]

  def displace_atoms(self):
    #self.magnify_cell()
    #self.read_data()
    #self.cal_disl_center()
    #for i_unitCell in range(0,sum(self.n_unit)):
    #  self.mag_atoms_pos[i_unitCell].pop(0)
    #  for i in self.mag_atoms_pos[i_unitCell]: 
    for i in self.atoms_pos:
        i[0],i[2] = i[0]+self.UxUz(i[0],i[1])[0],i[2]+self.UxUz(i[0],i[1])[1]
        self.disl_atoms_pos.append(i)
  def move_into_box(self,pos):
    for j in range(0,2):
      while (pos[j]>self.coord[j][j]):  #self.mag_coord[j][j]):
        pos[j]-=self.coord[j][j] #self.mag_coord[j][j]
      while (pos[j]<0):
        pos[j]+=self.coord[j][j] #self.mag_coord[j][j]
    return pos
  def move_remove_atoms(self):
    list_to_be_deleted=[]
    for i in range(0,len(self.disl_atoms_pos)):
      self.disl_atoms_pos[i]=self.move_into_box(self.disl_atoms_pos[i])
      for k in range(0,i):
        diff=np.asarray(self.disl_atoms_pos[i])-np.array(self.disl_atoms_pos[k])
        if (round(np.linalg.norm(diff),3) ==0) and (i != k):
          list_to_be_deleted.append(i)
    list_to_be_deleted.sort(reverse=True)
    #print list_to_be_deleted,len(list_to_be_deleted)
    for i in list_to_be_deleted:
      self.disl_atoms_pos.pop(i)
  def print_disl(self):
    self.read_data()
    self.cal_disl_center()
    self.make_dislocation()
    #self.displace_atoms()
    #print self.disl_atoms_pos
    #self.move_remove_atoms()
    print(self.sys_name)
    print(1.0) #self.latt_para
    for i in range(0,3):
      print(format(self.coord[i,0],"03f"),"	",format(self.coord[i,1],"03f"),"	",format(self.coord[i,2],"03f"))
    #for i in self.N[0]*self.N[1]*self.N[2]*np.asarray(self.n_unit):
    #  print i,
    print(len(self.atoms_pos))
    #print "\nCartesian" #self.coord_type
    print("Cartesian")
    for i in self.atoms_pos:
      print(format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f"))
disl1=gen_disl(sys.argv[1])#"unit_cell")
disl1.print_disl()
