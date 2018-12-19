#!/usr/bin/python
# todo: find correct Voterra equation for edge;
# solve problem with incorrect atomic positions under slip plane
import numpy as np
from numpy import pi,arctan,arctan2,log
import sys

class gen_disl():
  """generate a dislocation"""
  def __init__(self,filename):
    self.filename=filename
    self.latt_para=1.0
    self.b=1.0
    self.sys_name="" 
    self.coord_type="" #"Direct" #"Cartesian"
    self.coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) 
    self.atoms_pos=[] 
    self.N=[10,10,1] #default, will read from structural file
    self.n_unit=[]
    self.mag_coord=np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])
    self.mag_atoms_pos=[] 
    self.disl_center=np.array([0,0])
    self.disl_atoms_pos=[]
    #self.magnify_cell()
  def read_data(self):
    with open(self.filename,'r') as in_file:
      count=1
      for line in in_file:
        ll=line.split()
        if count==1: self.sys_name,self.b,self.N=ll[0],float(ll[1]),[int(ll[2]),int(ll[3]),int(ll[4])]
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
    if n0 != sum(self.n_unit): print "Some atomic positions are missing!"
    if self.coord_type == 'Direct':
      for ix in range(0,self.N[0]):
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):
            for k in range(0,n0):
               mag_atom = self.atoms_pos[k]+ix*np.array([1,0,0])+iy*np.array([0,1,0])+iz*np.array([0,0,1])  
               mag_atom = mag_atom/self.N #/[2.0,2.0,1.0]
               mag_atom = self.mag_coord.transpose().dot(mag_atom)
               self.mag_atoms_pos[k].append(mag_atom)
    elif self.coord_type == 'Cartesian':
      for ix in range(0,self.N[0]):            
        for iy in range(0,self.N[1]):
          for iz in range(0,self.N[2]):          
            for k in range(0,n0):         
               mag_atom =self.atoms_pos[k]+self.coord.transpose().dot(np.array([ix,iy,iz]))#ix*self.coord[0]+iy*self.coord[1]  
               self.mag_atoms_pos[k].append(mag_atom*self.latt_para) #? /self.N
  def cal_disl_center(self):
    center=np.array([0,0])
    for i in self.mag_coord:
      center = center +  np.array([i[0],i[1]]) #self.atoms_pos[i]
    center = center/2.0 #len(self.atoms_pos)
    self.disl_center= center
  def UxUy(self,x,y):
    nu=0.3
    e=1e-8 #avoid 1/0 error
    x,y=x-self.disl_center[0],y-self.disl_center[1]
    ux=self.b/(2*pi)*(arctan2(y,x)+x*y/(x**2+y**2+e)/(2*(1-nu)))
    uz=-self.b/(2*pi)*((1-2*nu)/4/(1-nu)*log(x**2+y**2+e)+(x**2-y**2)/(4*(1-nu))/(x**2+y**2+e))
    return [ux,uz] #self.b/(2*pi)*self.angle(x-self.disl_center[0],y-self.disl_center[1])
  def displace_atoms(self):
    self.magnify_cell()
    self.cal_disl_center()
    for i_unitCell in range(0,sum(self.n_unit)):
      self.mag_atoms_pos[i_unitCell].pop(0)
      for i in self.mag_atoms_pos[i_unitCell]: 
        i[0],i[1] = i[0]+self.UxUy(i[0],i[1])[0],i[1]+self.UxUy(i[0],i[1])[1]
        self.disl_atoms_pos.append(i)
  def print_disl(self):
    self.displace_atoms()
    print self.sys_name
    print 1.0 #self.latt_para
    for i in range(0,3):
      print format(self.mag_coord[i,0],"03f"),"	",format(self.mag_coord[i,1],"03f"),"	",format(self.mag_coord[i,2],"03f")
    for i in self.N[0]*self.N[1]*self.N[2]*np.asarray(self.n_unit):
      print i,
    print "\nCartesian" #self.coord_type
    for i in self.disl_atoms_pos:
      print format(i[0],"03f"),"	",format(i[1], "03f"),"	",format(i[2],"03f") 
disl1=gen_disl(sys.argv[1])#"unit_cell")
disl1.print_disl()
