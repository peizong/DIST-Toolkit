#!/usr/bin/python

#----------------Description of its function--------
# This code can create a supercell with a designated plane {hkl}(e.g. {1,1,1}) for
# any crystal structure.
# This is particularly important and useful in constructing interfacial defects,
# e.g., generalized stacking faults, twin boundaries, and grain boundaries. After
# such a supercell being created, the construction of a interfacial defect is
# straightforward.
#----------------User instruction-------------------
# The user has to supply basic parameters of the crystal structure, plane index,
# supercell size, etc. The output are data of created supercell stored in h5data
# format.
#-----------The structure of this code------------------
# Three parts: (i) definition of fundamental math functions;
# (ii) definition of supercell related functions based on (i) functions;
# (iii) prepare and print the supercell for gamma surface based on (ii) functions.
#---------------General rules-----------------------
# (i) functions all begin with F, e.g. FVecAdd
# (ii) input data all begin with I, e.g., IBravLatt
# (iii) output data all begin with O, e.g., OBravLatt

#__author__ = 'Zongrui Pei'

#from Structure.structure import CrystalStructure, AtomStructure
#import Utilities.h5data as h5
#import Utilities.utilities as u
import numpy as np
import fractions
#from Utilities.h5data import h5Data
from operator import itemgetter
import sys

#define 1rd level functions- mathematical functions
#---------------------------------------------------
# redefine the mathematical functions
pi=np.pi
array=np.array
FMatTranspose=np.transpose
FMatInv=np.linalg.inv
FMatMultiply=np.dot
#FMatDet=np.linalg.det
FVecCrossMul=np.cross
FRound=np.round
FTrunc=np.trunc

def FVecAdd(va,vb):
    vc = array([va[0]+vb[0],va[1]+vb[1],va[2]+vb[2]])
    return vc
def FVecMinus(va,vb):
    vc = array([va[0]-vb[0],va[1]-vb[1],va[2]-vb[2]])
    return vc
def FVecDot(va,vb):
    return va[0]*vb[0]+va[1]*vb[1]+va[2]*vb[2]
def FVecTriProduct(va,vb,vc):
    return FVecDot(FVecCrossMul(va,vb),vc)
def FVecLength(va):
    return np.sqrt(FVecDot(va,va))
def FVecsAngle(va,vb):
    return np.arccos(FVecDot(va,vb)/(FVecLength(va)*FVecLength(vb)))
#def FVecNorm(va):
#    if FVecLength(va)==0.0:
#        print "Attention! Normalized vector is 0 vector!"
#    else:
#        return va/FVecLength(va)
def FVecsAreEqual(va,vb):
    if (va[0]==vb[0]) & (va[1]==vb[1]) & (va[2]==vb[2]):
        return True
    else:
        return False
def FNewCoord(Coord,Trans_Matrix):
    return FMatMultiply(Coord,Trans_Matrix)
def FBravaisCoord(BravLatt):
    LattCons = BravLatt[0]
    angles = BravLatt[1]
    sca0 = np.square(np.cos(angles[0]))
    sca1 = np.square(np.cos(angles[1]))
    ssa2 = np.square(np.sin(angles[2]))
    ca0 = np.cos(angles[0])
    ca1 = np.cos(angles[1])
    ca2 = np.cos(angles[2])
    sa2 = np.sin(angles[2])
    pp = ca0 - ca1*ca2
    mm = np.sqrt(ssa2 - sca0 - sca1 + 2*ca0*ca1*ca2)
    if sa2 != 0.0:
        xyz = array([[LattCons[0],0,0],[LattCons[1]*ca2, LattCons[1]*sa2,0],
            [LattCons[2]*ca1, LattCons[2]*pp/sa2, LattCons[2]*mm/sa2]])
        return xyz
    else:
        print "Error, the crystal structure does not exist!"
        return 0
def FPlaneDistance(hkl,BravLatt):
    LattCons = BravLatt[0]
    angles= BravLatt[1]
    h,k,l = hkl[0],hkl[1],hkl[2]
    ca0 = np.cos(angles[0])
    ca1 = np.cos(angles[1])
    ca2 = np.cos(angles[2])
    sa0 = np.sin(angles[0])
    sa1 = np.sin(angles[1])
    sa2 = np.sin(angles[2])
    sca0 = np.square(ca0)
    sca1 = np.square(ca1)
    ssa2 = np.square(sa2)
    pp = ca0 - ca1*ca2
    nn = ca1 - ca2*ca0
    qq = ca2 - ca0*ca1
    mm = np.sqrt(ssa2 - sca0 - sca1 + 2*ca0*ca1*ca2)
    shkls = np.square(h/LattCons[0]*sa0) + np.square(k/LattCons[1]*sa1) + np.square(l/LattCons[2]*sa2)
    hklL = h*k*qq/(LattCons[0]*LattCons[1]) + h*l*nn/(LattCons[0]*LattCons[2]) + l*k*pp/(LattCons[1]*LattCons[2])
    return mm/np.sqrt(shkls-2*hklL)

def FGreatestComDivisor(na,nb):
    na,nb = np.absolute(na), np.absolute(nb)
    return fractions.gcd(na,nb)

def FIndicesRefine(va):
    h,k,l = va[0], va[1], va[2]
    Divisor = FGreatestComDivisor(FGreatestComDivisor(h,k),FGreatestComDivisor(k,l))
    if Divisor != 0.0:
        va=array([[h/Divisor,k/Divisor,l/Divisor]])
    return va

# definition of 1st level functions ends here
#--------------------------------------------------------------------
# define the 2nd level functions- main algorithms
# 2.1 find the primitive vectors based on [h,k,l]
def FNewPrimitiveVectors(hkl,BravLatt,CryStru):
    h,k,l = hkl[0], hkl[1], hkl[2]
    BravVectors=FBravaisCoord(BravLatt)
    shortAtomDis = 1
    secShortAtomDis = 1
    if FTrunc(CryStru) == 1:
        shortAtomDis = shortAtomDis
    elif CryStru == 2 or CryStru == 4:
        secShortAtomDis=0.5*2**0.5
    elif FTrunc(CryStru) == 3:
        secShortAtomDis=0.5*3**0.5
    else:
        print "Wrong crystal structure! Please check your input data!"
    #V0 = FVecTriProduct(BravVectors[0],BravVectors[1],BravVectors[2])/VFactor
    # first step: get three perpendicular vectors vs. hkl
    if h>0.0 :
        if k != 0.0 :
            a1 = array([-k,h,0.0])
        else:
            a1 = array([0.0,1.0,0.0])
        if l!= 0.0 :
            a2 = array([-l,0.0,h])
        else:
            a2 = array([0.0,0.0,1.0])
    elif h == 0.0 :
        if k == 0.0 :
            a1 = array([0.0,0.0,0.0])
            a2 = array([1.0,0.0,0.0])
        else:
            if l == 0.0 :
                a1 = array([1.0,0.0,0.0])
                a2 = array([0.0,0.0,0.0])
            else:
                a1 = array([1.0,0.0,0.0])
                a2 = array([1.0,0.0,0.0])
    else:
        if k != 0.0 :
            a1 = array([k,-h,0.0])
        else:
            a1 = array([0.0,1.0,0.0])
        if l !=0.0 :
            a2 = array([l,0.0,-h])
        else:
            a2 = array([0.0,0.0,1.0])
    if k > 0.0 :
        if l != 0.0:
            a3 = array([0.0,-l,k])
        else:
            a3 = array([0.0,0.0,1.0])
    elif k == 0.0:
        a3 = array([0.0,l,-k])
    else:
        if l != 0.0:
            a3 = array([0.0,l,-k])
        else:
            a3 = array([0.0,0.0,1.0])
    # second step: remove the greatest common divisor of the three vectors
    FIndicesRefine(a1)
    FIndicesRefine(a2)
    FIndicesRefine(a3)
    # third step: remove one of the three vectors
    if FVecsAreEqual(a1,array([0.0,0.0,0.0])):
        a,b = a2,a3
    elif FVecsAreEqual(a2,array([0.0,0.0,0.0])):
        a,b = a1,a3
    elif FVecsAreEqual(a3,array([0.0,0.0,0.0])):
        a,b = a1,a2
    elif FVecsAreEqual(a1,a2):
        a,b = a1,a3
    elif FVecsAreEqual(a1,a3):
        a,b = a1,a2
    elif FVecsAreEqual(a2,a3):
        a,b = a1,a3
    else:
        a,b = a1,a2
    # fourth step: make sure the selected vectors are shortest ones available
    # for non-cubic structure, it may not true, but it will not affect the results
    if FRound(FVecLength(a)/shortAtomDis,1)==2.0:
        a = a/2.0
    if FRound(FVecLength(a)/secShortAtomDis,1)==2.0:
        a = a/2.0
    if FRound(FVecLength(b)/shortAtomDis,1)==2.0:
        b = b/2.0
    if FRound(FVecLength(b)/secShortAtomDis,1)==2.0:
        b = b/2.0
    if FVecLength(a) > FVecLength(b):
        medium, a, b = a, b, medium
    aAb = FVecAdd(a,b)
    aMb = FVecMinus(a,b)
    if (FRound(FVecLength(aAb)/shortAtomDis,1)==2.0) or \
       (FRound(FVecLength(aMb)/shortAtomDis,1)==2.0):
        if FVecLength(aAb) > FVecLength(aMb):
            b = aMb/2.0
        else:
            b = aAb/2.0
    # translate the vectors into Cartesian coordinate
    c = array([h,k,l])
    a = FNewCoord(a,BravVectors)
    b = FNewCoord(b,BravVectors)
    c = FNewCoord(c,BravVectors)
    # make sure again that a is shorter than b
    if FVecLength(a) > FVecLength(b):
        m = a
        a = b
        b = m
    if FVecTriProduct(a,b,c) <0:
        c = -c
    elif FVecTriProduct(a,b,c) == 0:
        print "Error, found primitive vectors are not correct!"
    abc = array([a,b,c])
    return abc
# 2.2 Judge the positions of atoms relative to the supercell box: inside or outside
def FIsInBox(AtomPosition, BoxSize, CellType):
    ErrorCorrect = 1e-10
    if CellType==0 or CellType==2:
        delta = -1e-5
    elif CellType==1 or CellType==3:
        delta = 1e-5
    IsInBox_X = ((AtomPosition[0]+BoxSize[0,0]+ErrorCorrect)*
                 (AtomPosition[0]-BoxSize[0,1]-delta) < 0)
    IsInBox_Y = ((AtomPosition[1]+BoxSize[1,0]+ErrorCorrect)*
                 (AtomPosition[1]-BoxSize[1,1]-delta) < 0)
    IsInBox_Z = ((AtomPosition[2]+BoxSize[2,0]+ErrorCorrect)*
                 (AtomPosition[2]-BoxSize[2,1]-delta) < 0)
    if IsInBox_X and IsInBox_Y and IsInBox_Z:
        return True
    else:
        return False


# 2.3 define the space size of raw material which fills at least the box
def FSampleSize(SCellSize, Trans_Coord):
    SCellVertex8 = array([[SCellSize[0][0],SCellSize[1][0],SCellSize[2][0]],
                          [SCellSize[0][1],SCellSize[1][0],SCellSize[2][0]],
                          [SCellSize[0][0],SCellSize[1][1],SCellSize[2][0]],
                          [SCellSize[0][0],SCellSize[1][0],SCellSize[2][1]],
                          [SCellSize[0][1],SCellSize[1][1],SCellSize[2][0]],
                          [SCellSize[0][0],SCellSize[1][1],SCellSize[2][1]],
                          [SCellSize[0][1],SCellSize[1][0],SCellSize[2][1]],
                          [SCellSize[0][1],SCellSize[1][1],SCellSize[2][1]]])
    RestQuantity = ([[2,2,2]])
    Vertex8InSampleCoord = FNewCoord(SCellVertex8,Trans_Coord)
    NewSize1 = array(np.mat(array([Vertex8InSampleCoord.min(0)-RestQuantity,
                                    Vertex8InSampleCoord.max(0)+RestQuantity])))
    NewSize2 = FMatTranspose(FRound(NewSize1))
    return NewSize2

def FBravLattUnit(CryStru):
  """
  define the atomic coordinates for unit cell
  # definition of Bravai Lattices
  # 1-primitive, 1.1-hexagonal, 2-side-centered,3-body-centered, 4-face-centered,
  # 1.1-hexagonal,1.2-Monoclinic,1.3-Rhombic,1.4-Tetragonal,1.5-Trigonal(Rhombohedral),
    1.6-Triclinic,1.7-Cubic(simple lattice)
  # 2.1-Rhombic,2.2-Tetragonal,2.3-Cubic(bcc)
  # 3.1-Rhombic,3.2-Cubic(fcc)
  # 4.1-Monoclinic 4.2-Rhombic
  """
  if CryStru == 1.1:
      return array([[1.0/3,2.0/3,1.0/4],[2.0/3,1.0/3,3.0/4]])
  elif round(CryStru) == 1:
      return array([[0,0,0]])
  elif round(CryStru) == 2:
      return array([[0,0,0],[0.5,0.5,0]])
  elif round(CryStru) == 3:
      return array([[0,0,0],[0.5,0.5,0.5]])
  elif round(CryStru) == 4:
      return array([[0,0,0],[0.5,0.5,0],[0.5,0,0.5],[0,0.5,0.5]])
  else:
      print "other cases will be added!"
      return 0

# find the number of basis vectors;
def FNoOfBasisVector(CryStru):
  if CryStru == 1.1:
    return 2
  elif round(CryStru) == 1:
    return 1
  elif round(CryStru) == 2 or round(CryStru) == 3:
    return 2
  elif round(CryStru) == 4:
    return 4
  else:
    print("Incorrect input!!!")
# transform primitive vectors according to Cell Type
def FPrimitiveVectorsByCellType(CellType,OldPrimitiveVector):
    if CellType == 2 or CellType ==3 :
        LengthPrimitive = array([FVecLength(OldPrimitiveVector[0]),
                                 FVecLength(OldPrimitiveVector[1]),
                                 FVecLength(OldPrimitiveVector[2])])
        AnglePrimitive  = array([FVecsAngle(OldPrimitiveVector[1],OldPrimitiveVector[2]),
                                 FVecsAngle(OldPrimitiveVector[0],OldPrimitiveVector[2]),
                                 FVecsAngle(OldPrimitiveVector[0],OldPrimitiveVector[1])])
        NewBravLatt12 = array([LengthPrimitive,AnglePrimitive])
        return FBravaisCoord(NewBravLatt12)
    elif CellType == 0 or CellType == 1:
        return OldPrimitiveVector
    else:
        print "Wrong Cell Type! Please enter the cell type again!"
def FCoordinateType(CoordType):
    if CoordType == 0:
        return "Direct"
    elif CoordType == 1:
        return "Cartesian"
    else:
        return "unknown_coordinate"
def FAtomPosition(primitiveVecs,hklPrimitiveVecs,size,CryStru,CellType,CoordType):
    SCellSize = array([[0,size[0]],[0,size[1]],[0,size[2]]])
    MatTransTwoCoord = FMatMultiply(hklPrimitiveVecs,FMatInv(primitiveVecs))
    OPrimitiveVecs = FPrimitiveVectorsByCellType(CellType,hklPrimitiveVecs)
    SampleSize = FSampleSize(SCellSize, MatTransTwoCoord)
    IMatTransTwoCoord = FMatInv(MatTransTwoCoord)
    NoOfBasisVector=FNoOfBasisVector(CryStru)
    BravLattUnit=FBravLattUnit(CryStru)
    AtomPosition =[]
    for i in range(int(SampleSize[0][0]),int(SampleSize[0][1])):
        for j in range(int(SampleSize[1][0]),int(SampleSize[1][1])):
            for k in range(int(SampleSize[2][0]),int(SampleSize[2][1])):
                for iNoOfBasisVector in range (0,NoOfBasisVector):
                    point_in_sample_rel = FVecAdd(array([i,j,k]),BravLattUnit[iNoOfBasisVector])
                    point_in_NPrimitive = FNewCoord(point_in_sample_rel,IMatTransTwoCoord)
                    if FIsInBox(point_in_NPrimitive,SCellSize,CellType):
                        if CoordType == 1:
                            AtomPosition.append(FNewCoord(point_in_NPrimitive,OPrimitiveVecs)) #testing hklPrimitiveVecs
                        if CoordType == 0:
                            point_in_NPrimitive = point_in_NPrimitive/size
                            AtomPosition.append(point_in_NPrimitive)
    AtomPosition.sort(key=itemgetter(0))
    AtomPosition.sort(key=itemgetter(1))
    AtomPosition.sort(key=itemgetter(2))
    return AtomPosition
# definition of 2nd level functions end here
#----------------------------------------------------------------------------------
# define the 3rd level class
class prepareOutput:
  def __init__(self,input_file):
  #def __init__(self,SysName,LattPara,BravLatt,hkl,size,CryStru,CellType,CoordType):
    self.filename=input_file
    self.sys_name="System Name"
    self.latt_para=3 #Angstrom
    self.CryStru=3 #BCC
    self.BravLatt=array([[1,1,1],[pi/2.0,pi/2.0,pi/2.0]])
    self.hkl=np.array([1,1,1]) #(111)
    self.size=np.array([1,1,1]) # 1x1x1xN
    self.size33=np.array([[0,0,0],[0,0,0],[0,0,0]])
    self.CellType=2 # supercell for glide
    self.CoordType=0 # Direct
    self.OPrimitiveVectors=np.array([[0,0,0],[0,0,0],[0,0,0]])
    self.OSuperCellSize=1
    self.OCoordType="unknown"
    self.OAtomPosition=[]
  def read_data(self):                  
    with open(self.filename,'r') as in_file:
      count=1                           
      for line in in_file:              
        ll=line.split() 
        if count==1:                    
          self.sys_name=ll[0]           
         # for i in range(0,len(self.N)):
         #   self.randStruct.append([[0,0,0]])
        if count==2: self.latt_para=float(ll[0])
        if count==3: self.CryStru=float(ll[0])
        if count==4: self.BravLatt=np.array([[float(ll[0]),float(ll[1]),float(ll[2])], 
                                             [pi*float(ll[3]),pi*float(ll[4]),pi*float(ll[5])]]) 
        if count==5: self.hkl=np.array([float(ll[0]),float(ll[1]),float(ll[2])])
        if count==6:                    
          self.size=np.array([float(ll[0]),float(ll[1]),float(ll[2])])
          self.size33 = np.array([[self.size[0],self.size[0],self.size[0]],
                                  [self.size[1],self.size[1],self.size[1]],
                                  [self.size[2],self.size[2],self.size[2]]])
        if count==7: self.CellType=float(ll[0]) 
        if count==8: self.CoordType=float(ll[0])
        if count>8: break
        count +=1
  def cal_parameters(self):
    primitiveVectors=FBravaisCoord(self.BravLatt) #BravLatt)
    hklPrimitiveVecs = FNewPrimitiveVectors(self.hkl,self.BravLatt,self.CryStru)
    #self.OSysName = SysName
    #self.OLattPara = LattPara
    #self.size33 = array([[size[0],size[0],size[0]],
    #                     [size[1],size[1],size[1]],
    #                     [size[2],size[2],size[2]]])
    self.OPrimitiveVectors = FPrimitiveVectorsByCellType(self.CellType,hklPrimitiveVecs*self.size33)
    self.OCoordType=FCoordinateType(self.CoordType)
    self.OAtomPosition=FAtomPosition(primitiveVectors,hklPrimitiveVecs,self.size,self.CryStru,self.CellType,self.CoordType)
    self.OSuperCellSize = len(self.OAtomPosition)
  def print_supercell(self):
    ''' 3.3 write structural file'''
    self.read_data()
    self.cal_parameters()
    print self.sys_name
    print self.latt_para
    for i in self.OPrimitiveVectors:
      print i[0],"  ",i[1],"        ",i[2]
    print self.OSuperCellSize
    print self.OCoordType
    for i in self.OAtomPosition:
      print i[0],"      ",i[1],"        ",i[2]

Mg=prepareOutput(sys.argv[1]) #("in.gen-hkl-unit-cell")
Mg.print_supercell()

# Here come the input data: CompSysName, SuperCellSize, hkl, BravLatt, CellType
#------------------------------------------------------------------------------------
# definition of Bravai Lattices
# 1-primitive, 1.1-hexagonal, 2-side-centered,3-body-centered, 4-face-centered,
# 1.1-hexagonal,1.2-Monoclinic,1.3-Rhombic,1.4-Tetragonal,1.5-Trigonal(Rhombohedral),
# 1.6-Triclinic,1.7-Cubic(simple lattice)
# 2.1-Rhombic,2.2-Tetragonal,2.3-Cubic(bcc)
# 3.1-Rhombic,3.2-Cubic(fcc)
# 4.1-Monoclinic 4.2-Rhombic
#ISysName = "Mg"
#ILattPara    = 3.1886 # La in Angstom for pure Mg
#----------for simple cubic----------------------------
#ICryStru  = 1 # check the definition of CryStru
#IBravLatt = array([[1,1,1],[pi/2.0,pi/2.0,pi/2.0]])
#----------for hcp----------------------------
#ICryStru  = 1.1 # check the definition of CryStru
#IBravLatt = array([[1,1,1.6261],[pi/2.0,pi/2.0,pi*2/3.0]])
#----------for simple cubic----------------------------
#ICryStru  = 2 # check the definition of CryStru
#IBravLatt = array([[1,1,1],[pi/2.0,pi/2.0,pi/2.0]])
#----------for bcc-----------------------------
#ICryStru=3
#IBravLatt = array([[1,1,1],[pi/2.0,pi/2.0,pi/2.0]]) # [[1,Lb/La,Lc/La],[alpha,beta,gamma]]
#----------for fcc-----------------------------
#ICryStru=4
#IBravLatt = array([[1,1,1],[pi/2.0,pi/2.0,pi/2.0]]) # [[1,Lb/La,Lc/La],[alpha,beta,gamma]]
#Ihkl = array([1,0,2]) # the Bravai indice {H K L}
#Isize = array([1,1,6]) # the number of period along each direction N1,N2,N3
#ICellType = 2  # 0-Poscar cell, 1-Complete cell, 2-Poscar cell for gliding, 3-Complete cell for gliding
#ICoordType = 0 # 0- Direct; 1- Cartesian
#-------------------------------------------------------------------------------------------------
#---------print out as hdf5 data format
#structure = h5Data( file_name= "POSCAR.pzr", path= u.PATH_DUMP) # PATH_DUMP = C:\Users\Pei\workspace\data\dump
#structure.addArray(name='materialSystem',val=x.OSysName)
#structure.addArray(name='latticeParameter',val=x.OLattPara)
#structure.addArray(name='primitiveVector',val=x.OPrimitiveVectors)
#structure.addArray(name='NOofAtoms',val=x.OSuperCellSize)
#structure.addArray(name='typeofCoordinate',val=x.OCoordinate)
#structure.addArray(name='Positions', val=x.OAtomPosition)
#structure.move_up()
#structure.close()
#exit()

