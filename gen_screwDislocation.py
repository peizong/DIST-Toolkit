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
    self.disl_centers=[]
    self.sub_matr2x2=[]
    self.sign=[]
    self.mag_atoms_pos=[] 
    self.disl_center=np.array([0,0])
    self.disl_atoms_pos=[]
    self.disl_atoms_pos_in=[]
    self.period_N=20 # used to calculate u_err to get absolutely convergent strain field
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
    #need to remove this line!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    self.x0=1.0/3 #-0.017
    L2,L3=L0,(1-2*self.x0)*(L0+L1)
    self.theta0=arccos(np.dot(L2,L3)/np.dot(L2,L2)**0.5/np.dot(L3,L3)**0.5)
    # define self.disl_centers, self.sign, sub_matr2x2
    disl_centers=[]
    if self.num_disl==4:
      self.sign=np.array([1,-1,-1,1])
      disl_centers=np.array([[0.25,0.25],[0.25,0.75],[0.75,0.25],[0.75,0.75]])
    elif self.num_disl==2:
      self.sign=np.array([-1,1]) #np.array([1,-1])
      disl_centers=np.array([[self.x0,self.x0],[1-self.x0,1-self.x0]])
      #disl_centers=np.array([[self.x0,1-self.x0],[1-self.x0,self.x0]])
    elif self.num_disl==1:
      self.sign=np.array([1])
      disl_centers=np.array([[0.0,0.0]])
    else: return None
    self.sub_matr2x2=np.array([[self.mag_coord[0][0],self.mag_coord[0][1]],
                               [self.mag_coord[1][0],self.mag_coord[1][1]]])
    self.disl_centers=np.dot(disl_centers,self.sub_matr2x2)
  def angle(self,x,y):
     if x==0:
       if y>=0: return pi/2.0
       if y<0: return 3*pi/2.0
     elif x>0:
       if y >=0 : return arctan(y/x)
       if y <0: return 2*pi+arctan(y/x)
     elif x<0 : return pi+arctan(y/x)
  def adjust_angle(self,x,y):
    start_angle=self.theta0 # adjust the starting angle, currently not useful, but not harmful
    adjusted_angle=self.angle(x,y)+start_angle
    while (adjusted_angle <0): adjusted_angle += 2*pi
    while (adjusted_angle >2*pi): adjusted_angle -= 2*pi
    return adjusted_angle
  #recover the periodicity of the dislocation strain field by adding ghost dislocation
  def get_u_err(self):
    u_arr=[]
    # the three corners, OO,OA,OB
    corner_pts=np.array([[0.,0.],[self.sub_matr2x2[0][0],self.sub_matr2x2[0][1]],\
                                 [self.sub_matr2x2[1][0],self.sub_matr2x2[1][1]]])
    for pts in corner_pts:
      u_pts=0.
      for i in range(0,self.num_disl):
        for jx in range(0,self.period_N): #20):
          for jy in range(0,self.period_N): #20):
            new_center=(jx-self.period_N/2)*self.sub_matr2x2[0]+(jy-self.period_N/2)*self.sub_matr2x2[1]+self.disl_centers[i]
            dc_x,dc_y=new_center[0],new_center[1]
            u_pts +=self.b/(2*pi)*self.sign[i]*self.adjust_angle(pts[0]-dc_x,pts[1]-dc_y)
      u_arr.append(u_pts)
    OAB_mat=np.array([[corner_pts[0][0],corner_pts[0][1],1.],
                      [corner_pts[1][0],corner_pts[1][1],1.],
                      [corner_pts[2][0],corner_pts[2][1],1.]])
    g_mat=np.dot(np.linalg.inv(OAB_mat),np.asarray(u_arr))
    g_vec,g0=g_mat[0:2],g_mat[2]
    return g_vec, g0
  def theta(self,x,y,N): #,theta0):
    shift_z=0
    for i in range(0,self.num_disl):
      for jx in range(0,self.period_N): #20):
        for jy in range(0,self.period_N): #20):
          new_center=(jx-self.period_N/2)*self.sub_matr2x2[0]+ \
                     (jy-self.period_N/2)*self.sub_matr2x2[1]+self.disl_centers[i]
          dc_x,dc_y=new_center[0],new_center[1]
          shift_z +=self.b/(2*pi)*self.sign[i]*self.adjust_angle(x-dc_x,y-dc_y) 
    return shift_z
  def displace_atoms(self):
    self.read_data()
    self.magnify_cell()
    self.cal_disl_pattern()
    g_vec,g0=self.get_u_err() # get the g_vec, g0 to restore the periodicity
    for i_unitCell in range(0,sum(self.n_unit)):
      self.mag_atoms_pos[i_unitCell].pop(0)
      for j in range(0,len(self.mag_atoms_pos[i_unitCell])):
        i=self.mag_atoms_pos[i_unitCell][j]
      #for i in self.mag_atoms_pos[i_unitCell]: 
        u_err=np.dot(g_vec,np.array([i[0],i[1]]))
        i[2] += self.theta(i[0],i[1],self.num_disl)-u_err 
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
    print "\nCartesian" #self.coord_type
    for i in self.disl_atoms_pos_in:
      print format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f") 
disl1=gen_disl(sys.argv[1])#"unit_cell")
disl1.print_disl()
