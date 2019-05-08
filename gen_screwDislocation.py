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
from numpy import pi,arccos,arctan
from numpy.linalg import inv
import sys

class gen_disl():
  """generate a dislocation"""
  def __init__(self,filename):
    self.filename=filename
    self.latt_para=1.0
    self.b=1.0
    self.num_disl=2
    self.sys_name="" 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.N=[1,1,1] #default, will read from structural file
    self.n_unit=[]
    self.mag_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    self.mag_atoms_pos=[] 
    self.disl_center=np.array([0,0])
    self.disl_atoms_pos=[]
    self.disl_atoms_pos_in=[]
    #self.magnify_cell()
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1: 
          if len(ll)==6:
            self.sys_name,self.b,self.num_disl, self.N=ll[0],float(ll[1]),float(ll[2]),[int(ll[3]),int(ll[4]),int(ll[5])]
          elif len(ll)==5:
            self.sys_name,self.b,self.N=ll[0],float(ll[1]),[int(ll[2]),int(ll[3]),int(ll[4])]
          elif len(ll)==3:
            self.sys_name,self.b,self.num_disl=ll[0],float(ll[1]),int(ll[2])
          elif len(ll)==2:
            self.sys_name,self.b=ll[0],float(ll[1])
          else: 
            return None
            print("Error with the first line of input file!")
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
    k1=0
    for i in self.coord:
      self.mag_coord[k1]=np.array([i[0]*self.N[0],i[1]*self.N[1],self.N[2]*i[2]])*self.latt_para
      k1 +=1
      if k1==3: break #exit
    n0=len(self.atoms_pos)
    if n0 != sum(self.n_unit): print "Some atomic positions are missing!"
    if self.coord_type == 'Direct':
      for ix in range(0,self.N[0]):
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):
            for k in range(0,n0):
               mag_atom = self.atoms_pos[k]+ix*np.array([1,0,0])+iy*np.array([0,1,0])+iz*np.array([0,0,1])  
               mag_atom = mag_atom/self.N #/[2.0,2.0,1.0]
               mag_atom = np.dot(self.mag_coord,mag_atom) #self.mag_coord.transpose().dot(mag_atom)
               self.mag_atoms_pos[k].append(mag_atom)
    elif self.coord_type == 'Cartesian':
      for k in range(0,n0):
        for ix in range(0,self.N[0]):            
          for iy in range(0,self.N[1]):
            for iz in range(0,self.N[2]):          
               mag_atom =np.asarray(self.atoms_pos[k])+np.dot(self.coord,np.array([ix,iy,iz]))  #self.coord.transpose().dot(np.array([ix,iy,iz]))#ix*self.coord[0]+iy*self.coord[1]  
               self.mag_atoms_pos[k].append(mag_atom) #*self.latt_para) #? /self.N
               #print mag_atom, "after:","	",self.mag_atoms_pos[k]
    else: return None
  def cal_disl_pattern(self):
    L0=np.array([self.mag_coord[0][0],self.mag_coord[0][1]])
    L1=np.array([self.mag_coord[1][0],self.mag_coord[1][1]])
    #self.x0=0.5*(np.dot(L0,L0)+np.dot(L0,L1))/\
    #        (np.dot(L0,L0)+np.dot(L1,L1)+4*np.dot(L0,L1))
    # for Left bottom right up
    self.x0=0.5*(np.dot(L0,L1)+np.dot(L1,L1))/\
            (np.dot(L0,L0)+np.dot(L1,L1)+2*np.dot(L0,L1))
    L2,L3=L0,(1-2*self.x0)*(L0+L1)
    self.theta0=arccos(np.dot(L2,L3)/np.dot(L2,L2)**0.5/np.dot(L3,L3)**0.5)
    #self.theta0=np.pi-arccos(np.dot(L2,L3)/np.dot(L2,L2)**0.5/np.dot(L3,L3)**0.5)
    # for Left Up right bottom
    #self.x0=0.5*(np.dot(L1,L1)-np.dot(L0,L1))/\
    #        (np.dot(L0,L0)+np.dot(L1,L1)-2*np.dot(L0,L1))
    #L2,L3=L0,-(1-self.x0)*(L0-L1)
    #self.theta0=arccos(np.dot(L2,L3)/np.dot(L2,L2)**0.5/np.dot(L3,L3)**0.5)
  def angle(self,x,y):
     if x==0:
       if y>=0: return pi/2.0
       if y<0: return 3*pi/2.0
     elif x>0:
       if y >=0 : return arctan(y/x)
       if y <0: return 2*pi+arctan(y/x)
     elif x<0 : return pi+arctan(y/x)
  def adjust_angle(self,x,y,start_angle):
    adjusted_angle=self.angle(x,y)+start_angle
    while (adjusted_angle <0): adjusted_angle += 2*pi
    while (adjusted_angle >2*pi): adjusted_angle -= 2*pi
    return adjusted_angle
  def theta(self,x,y,N,theta0):
    shift_z=0
    if self.num_disl==4: 
      sign=np.array([1,-1,-1,1])
      disl_centers=np.array([[0.25,0.25],[0.25,0.75],[0.75,0.25],[0.75,0.75]])
    elif self.num_disl==2: 
      sign=np.array([1,-1])
      #disl_centers=np.array([[self.x0+1.0/self.mag_coord[0][0]*0.21213*5,self.x0+1.0/self.mag_coord[1][1]*0.12247*5],[1-self.x0,1-self.x0]])
      disl_centers=np.array([[self.x0,self.x0],[1-self.x0,1-self.x0]])
      #disl_centers=np.array([[self.x0,1-self.x0],[1-self.x0,self.x0]])
    elif self.num_disl==1: 
      sign=np.array([1])
      disl_centers=np.array([[0.0,0.0]])
    else: return None
    disl_centers=disl_centers*np.array([self.mag_coord[0][0],self.mag_coord[1][1]])
# tobedeleted
#    disl_centers=np.array([[4.24,4.5],[12.37,13.0]])
    for i in range(0,self.num_disl):
      #shift_z +=self.b/(2*pi)*sign[i]*self.adjust_angle(x-disl_centers[i,0],y-disl_centers[i,1],theta0)
    #return self.b/(2*pi)*self.adjust_angle(x-self.disl_center[0],y-self.disl_center[1],0/2.0)
      for jx in range(0,80):
        dc_x=(jx-40)*self.mag_coord[0][0]+disl_centers[i,0]
        for jy in range(0,80):
          dc_y=(jy-40)*self.mag_coord[1][1]+disl_centers[i,1]
          shift_z +=self.b/(2*pi)*sign[i]*self.adjust_angle(x-dc_x,y-dc_y,theta0)
    return shift_z
# to be deleted
  def ur(self,x,y):
    u0=[-0.25,0.25]
    r=[x,y]
    return -np.dot(u0,r)
  def displace_atoms(self):
    self.read_data()
    self.magnify_cell()
    self.cal_disl_pattern()
    for i_unitCell in range(0,sum(self.n_unit)):
      self.mag_atoms_pos[i_unitCell].pop(0)
      for j in range(0,len(self.mag_atoms_pos[i_unitCell])):
        i=self.mag_atoms_pos[i_unitCell][j]
      #for i in self.mag_atoms_pos[i_unitCell]: 
        i[2] += self.theta(i[0],i[1],self.num_disl,self.theta0) #+self.ur(i[0],i[1]) #use 4) for quadruple poles
        self.disl_atoms_pos.append(i)
  def in_box(self,atom_pos):
     for i in range(0,3):
       while(atom_pos[i]<0.0): atom_pos[i]=atom_pos[i]+1.0
       while( atom_pos[i]>=1.0): atom_pos[i]=atom_pos[i]-1.0
     return atom_pos
  def move_in_box(self):
     atom_in_direct=[]
     for i in self.disl_atoms_pos:
       i=self.in_box(np.dot(inv(self.mag_coord.transpose()),i))
       atom_in_direct.append(i)
#       atom_in_direct.append(self.in_box(np.dot(inv(self.mag_coord.transpose()),i)))    
     #if self.coord_type=='Direct':
     #  self.write_atom_pos=atom_in_direct
     if self.coord_type=='Cartesian':
       for i in atom_in_direct:
         i=np.dot(self.mag_coord.transpose(),i)
         self.disl_atoms_pos_in.append(i)
     else:
       self.disl_atoms_pos_in=atom_in_direct
  def print_disl(self):
    self.displace_atoms()
    self.move_in_box()
    print self.sys_name
    print 1.0 #self.latt_para
    for i in range(0,3):
      print format(self.mag_coord[i,0],"03f"),"	",format(self.mag_coord[i,1],"03f"),"	",format(self.mag_coord[i,2],"03f")
    for i in self.N[0]*self.N[1]*self.N[2]*np.asarray(self.n_unit):
      print i,
#    print '\n'
    print "\nCartesian" #self.coord_type
#    print(self.mag_atoms_pos)
    for i in self.disl_atoms_pos_in:
    # for i in j:
      print format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f") 
disl1=gen_disl(sys.argv[1])#"unit_cell")
disl1.print_disl()
