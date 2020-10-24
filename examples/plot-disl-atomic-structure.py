#!/usr/bin/python
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from numpy import pi
#data=np.genfromtxt('Mg.dat',delimiter='',skip_header=3)
#print data
def atan(x):
    return np.arctan(x)

class atomic_structure:
  """plot the atomic structure of dislocation"""
  def __init__(self, N,disl_size,mater_properties, file_name, delimiter_file, skip_header_file):
    self.N=N
    self.disl_size=disl_size
    self.mater_properties=mater_properties
    self.filename=file_name
    self.delimiter=delimiter_file
    self.skip_header=skip_header_file
    self.data=[]
    self.p=[] # parameters for describing dislocation core
    self.pressure=[]
    self.header=[]
  def read_data(self):
    self.data=np.genfromtxt(self.filename, delimiter=self.delimiter, skip_header=self.skip_header)
    #self.c=self.data #[self.num_arctan,self.data[4,0],self.data[4,1] ]

  def parse_c(self,c):
    N = self.N
    ind=0
    loc=np.linspace(0,0,2*2*(2*N+1))
    w=np.linspace(0,0,2*2*(2*N+1))
    A=np.linspace(0,0,2*2*(2*N+1))
    for k in range(0,2):
      for j in range(0,2):
        for i in range(0,2*N+1):
          ind=i+j*(2*N+1)+k*(2*(2*N+1))
          loc[ind]=(c[2*k]+(i-N)*c[2*k+1])*(-1)**j
          w[ind]=c[4+k]
          A[ind]=1.0/(2*(2*N+1))*(1+k)
          if (ind>3*(2*N+1)-1): A[ind]=-A[ind]
    return [A,loc,w]

  def pressure_field(self,x,y):
    p=self.parse_c(self.data)
    A,loc,w=p[0],p[1],p[2]
    G,v,b=self.mater_properties[0],self.mater_properties[1],self.mater_properties[2]
    #G=1.96e9
    #v=0.3316
    #b=3.1886e-10
    const=-G*(1+v)*b/3/np.pi/(1-v)
    tot_pressure = 0 
    for i in range(0,len(A)/2):
      tot_pressure = tot_pressure-const*A[i]*(y+np.sign(y)*w[i])/((x-loc[i])**2+(y+np.sign(y)*w[i])**2)
    self.pressure = tot_pressure
    print self.pressure
  def U_msft(self,x):
    b,az=1,3**0.5
    h=az/6
    p=self.parse_c(self.data)
    A,loc,w=p[0],p[1],p[2]
    ux,uz=b/2.0,0.0
    for j in range(0,2*(2*N+1)):
       ux=ux+b/pi*A[j]*atan((x-loc[j])/w[j])
       uz=uz+h/pi*A[j+2*(2*N+1)]*atan((x-loc[j+2*(2*N+1)])/w[j+2*(2*N+1)])
    ux,uz = ux/b, uz/az
    return [ux,uz]

  def plot_disl(self):
    self.read_data()
    #p0=self.data[4,1]
    N=self.disl_size #10                 
    #b=3.1886e-10         
    atom1X=np.linspace(0,0,2*N*2*N)                      
    atom1Y=np.linspace(0,0,2*N*2*N)    
    #c=np.linspace(0,0,2*N*2*N)    
    atom1=np.array([1, 1])   
    Nx=np.array([1, 0])  
    Ny=np.array([0, 1])  
    for i in range(0,2*N):  
        for j in range(0,2*N):  
            new_atom1=atom1+i*Nx+j*Ny   
            if j<N:      
                u=self.U_msft(i-N)
                new_atom1=new_atom1+np.array([u[0]-0,0]) 
            atom1X[i*(2*N)+j]=new_atom1[0]-N-0.5 
            atom1Y[i*(2*N)+j]=new_atom1[1]-N-0.5  
            #c[i*(2*N)+j]=pressure_field(atom1X[i*(2*N)+j]*b,atom1Y[i*(2*N)+j]*b)                                             
    self.pressure_field(atom1X,atom1Y)
    cm = plt.cm.get_cmap('RdYlBu')
    header=plt.scatter(atom1X,atom1Y,c=self.pressure, vmin=-0.25, vmax=0.15,s=300)
    plt.axis([-10, 10, -10, 10],fontsize=18)
    plt.xlabel('x (1/3[11-20])',fontsize=18)
    plt.ylabel('y ([1-100])',fontsize=18)
    plt.legend(str(self.N), loc='best')
    cb=plt.colorbar(header)
    cb.set_label(label='pressure in Pa',fontsize=15)

N=2 # the N in MEFD
disl_size =10
G=1.96e9
v=0.3316
b=3.1886e-10
mater_properties=[G,v,b]
case1=atomic_structure(N,disl_size,mater_properties,'best1',"	",1)
#case2=atomic_structure(5,'mesh-20points/11arctan.dat',",",1)
#plt.subplot(121)
#case2.plot_disl()
#plt.subplot(122)
case1.plot_disl()
#plt.colorbar(case2.header) 
plt.show()

