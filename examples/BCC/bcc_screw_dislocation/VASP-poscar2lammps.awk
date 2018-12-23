#!/bin/awk -f
#############################################################################
#                                                                           #
# GNU License - Author: Pavlin D. Mitev                        2012-10-09   #
# Version 0.1                                                               #
#                                                                           #
# Converts VASP POSCAR to LAMMPS imput structure.                           #
# It can handle non-orthogonal simulation boxes.                            #
# Reads both v4.6 and 5.2 (atom labels) POSCAR.                             #
# Accepts "Direct" and "Cartesian" coordinates in the POSCAR.               #
#                                                                           #
# Syntax:                                                                   #
# VASP-poscar2lammps.awk POSCAR > structure.lammps                          #
#                                                                           #
#############################################################################
BEGIN{
  pi=3.14159265358979;rad2deg=180./pi;
  if(ARGC<=1) { print "Syntax: \n      VASP-poscar2res.awk  POSCAR  ..."; ex=1;exit}
  for(i=2;i<=ARGC;i++){typeT[i-1]=ARGV[i];}
    ARGC=2;
    
  # Read header, scale and the basis -----------  
  getline; title=$0
  getline; scale=$1
  getline; h1[1]=$1*scale; h1[2]=$2*scale; h1[3]=$3*scale;
  getline; h2[1]=$1*scale; h2[2]=$2*scale; h2[3]=$3*scale;
  getline; h3[1]=$1*scale; h3[2]=$2*scale; h3[3]=$3*scale;

  a=norm(h1); b=norm(h2); c=norm(h3);                                # Length of the basis vectors
  alpha=  angle(h2,h3);  beta=  angle(h1,h3); gamma=  angle(h1,h2);  # Angles in degree
  alphar= alpha/rad2deg; betar= beta/rad2deg; gammar= gamma/rad2deg; # Angles in radians
  

  # Check for labels -------------------------
  getline; 
  if ($1*1 != $1) {
    for(i=1;i<=NF;i++){typeT[i]=$i;}
    getline;
  }

  for(i=1;i<=NF;i++) {type[i]=$i; natoms=natoms+$i} ntypes=NF;
  # Advance to the coordinate section
  while((tolower($0) !~ "direct")&&(tolower($0) !~ "cart")) getline;

  if (tolower($0) ~ "direct") fractional=1; # Fractional format identified

  # Rotation of the matrix to comply with LAMMPS standards ========================
  p_a= sqrt(h1[1]**2 + h1[2]**2 + h1[3]**2);
  p_b= sqrt(h2[1]**2 + h2[2]**2 + h2[3]**2);
  p_c= sqrt(h3[1]**2 + h3[2]**2 + h3[3]**2);
  lx=   p_a;
  p_xy= p_b * cos(gammar);
  p_xz= p_c * cos(betar);
  ly=   sqrt(p_b**2 - p_xy**2);
  p_yz= (p_b*p_c*cos(alphar)-p_xy*p_xz)/(ly);
  lz=   sqrt(p_c**2 - p_xz**2 - p_yz**2);
  # The new basis H matrix ------------------------
  H1[1]=lx;   H1[2]= 0.000; H1[3]= 0.00;
  H2[1]=p_xy; H2[2]= ly;    H2[3]= 0.00;
  H3[1]=p_xz; H3[2]= p_yz;  H3[3]= lz;
  # Matrix for conversion from cartesian to fractional in the old basis set (if necessary) ----------
  cfv= sqrt(1. -cos(alphar)**2 -cos(betar)**2 -cos(gammar)**2 + 2.*cos(alphar)*cos(betar)*cos(gammar));
  cf1[1]= 1./a;  cf1[2]= -cos(gammar)/(a*sin(gammar));   cf1[3]= (cos(alphar)*cos(gammar)-cos(betar))/(a*cfv*sin(gammar)); 
  cf2[1]= 0.00; cf2[2]= 1./(b*sin(gammar));   cf2[3]= (cos(betar)*cos(gammar)-cos(alphar)i)/(b*cfv*sin(gammar));
  cf3[1]= 0.00; cf3[2]= 0.00; cf3[3]= sin(gammar)/(c*cfv);

  # ===============================================================================

  print "# Converted from POSCAR to lammps format"
  print ""
  print natoms" atoms"
  print ntypes" atom types"
  print ""
  printf "0.000000  %10.6f   xlo xhi\n",lx
  printf "0.000000  %10.6f   ylo yhi\n",ly
  printf "0.000000  %10.6f   zlo zhi\n",lz
  print ""
  printf "%10.6f  %10.6f  %10.6f   xy xz yz\n", p_xy,p_xz,p_yz
  print ""
  print "Atoms"
  print ""

  iatom=0;
  for(k=1;k<=ntypes;k++){
    for(i=1;i<=type[k];i++){
      getline
      iatom++
      x=$1; y=$2; z=$3;     
      if (!fractional){
        xx=x*cf1[1]+y*cf1[2]+z*cf1[3];
        yy=x*cf2[1]+y*cf2[2]+z*cf2[3];
        zz=x*cf3[1]+y*cf3[2]+z*cf3[3];
        x= xx; y= yy; z= zz;
      }
      xx=x*H1[1]+y*H2[1]+z*H3[1];
      yy=x*H1[2]+y*H2[2]+z*H3[2];
      zz=x*H1[3]+y*H2[3]+z*H3[3];
      printf"%4i  %-4s   %7f %7f %7f\n",iatom, k, xx,yy,zz
    }
  }

}

function asin(a)       { return atan2(a,sqrt(1-a*a)) }
function acos(a)       { return pi/2-asin(a) }
function norm(x)       {return (sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]));}
function dotprod (x,y) {return ( x[1]*y[1] + x[2]*y[2] + x[3]*y[3] );}
function angle (v1,v2) {
  myacos = dotprod(v1,v2)/norm(v1)/norm(v2);
  if (myacos>1.0) myacos = 1.0;
  if (myacos<-1.0) myacos = -1.0;
  return(acos(myacos)*180.0/3.14159265358979);
}
