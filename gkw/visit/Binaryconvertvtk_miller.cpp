// This code W A Hornsby 2008
// Modified 2010 F J Casson
// For GKW data visualisation
// Takes GKW 3D Output and makes Visit VTK files.
// Files needed: 3DOutputparm.dat + 3D binary data files.
// Run GKW with at least 64 s points (128 is better) and den3d/ene3d/phi3d=.true.
// ComPile with "[g/i]cc -lm -O2 Binaryconvertvtk.cpp" and run from the data folder.

// To Do: Add ability to create fake s points by linear interpolation.
// To Do: Output not in ASCII format

// Modified to allow plotting on arbitrary Miller parameterised flux surfaces
// Using params read at runtime.
// For a nice looking artifically constructed "equilibirum", of nested flux surfaces
// I used these params:
//
// eps    triang  shift   elong
// 0.09   0.25    0.10     1.30
// 0.13   0.30    0.09     1.32
// 0.17   0.35    0.08     1.34
// 0.21   0.40    0.07     1.36
// 0.25   0.45    0.06     1.38
// 0.29   0.50    0.05     1.40

#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>


void vtkfile(double* temp,int s,int x,int y,int t,double Lx,double Ly);
void threedtorusvtk(double eps, double delta, double kappa, double alpha);
void potentialtorusvtk(double* temp,double* sgrid,int grids, double* psi,int gridx, double* chi,int gridy,int count,double eps,double q);

using namespace std;

ifstream::pos_type size;
double * memblock;
double rho_star, shat; //ascetic parameters
double delta, kappa, alpha;
double Lx,Ly,dx,dy;
const double Pi = 3.1416;
int repeats;

int main () 
{
  int i,j,count,nframes,ninterval,fframe;
  
  char str[10];
  char *str2;
  double f;

  string prefix;
  string postfix;
  //Modify this to plot different data file.
  prefix = "Poten";
  ostringstream strs;

  //!!!! THE FOLLOWING SECTION READS IN VARIOUS PARAMETERS
  //!!!! FOR THE COORDINATE TRANSFORMATIONS


  ifstream finput("3DOutputParam.dat");
  if(!finput)
    {
      cout << "3DOutputParam.dat error - Program aborted" << endl;
      return -1;
    }
  finput >> str >> f;
  cout << str << " " << f  << endl;
  Lx = f;
  finput >> str >> f;
  cout << str << " " << f  << endl;
  Ly = f;
  finput >> str >> f;
  cout << str << " " << f  << endl;
  int gridzeta = (int)f;
  finput >> str >> f;
  cout << str << " " << f  << endl;
  int gridpsi = (int)f;
  finput >> str >> f;
  cout << str << " " << f  << endl;
  int grids = (int)f;
  finput >> str >> f;
  cout << str << " " << f  << endl;
  double q = f;
  finput >> str >> f;
  cout << str << " " << f  << endl;
  double eps = f;

  //Renormalise Ly with kthnorm
  //Ly=Ly/(2.*Pi*eps/q);

  cout << "Perpendicular Lx, Ly, sizes in" << endl;
  cout << Lx << " " << Ly << endl;

  cout << "Number of toroidal repeats? 0=Use Flux Tube, 1=Stretch flux tube round torus, >2=repeat flux tube and stretch" << endl;
  cin >> repeats;
  cout << "Rho_star? 0.1/LX recommended:" << endl;
  cin >> rho_star;

  //rescale Ly so that Copies are an exact multiple of turns
  if(repeats>0) Ly=1./(rho_star*repeats);

  cout << "Perpendicular Lx, Ly sizes out" << endl;
  cout << Lx << " " << Ly << endl;

  double s[grids];
  double psi[gridpsi];
  double zeta[gridzeta];

  dx = Lx/gridpsi;
  dy = Ly/gridzeta;
  cout << gridpsi << " " << gridzeta << " " << grids << endl;
  cout << dx << " " << dy << endl;

  psi[0]=0.0;
  zeta[0]=0.0;
  cout << psi[gridpsi-1] << " " << zeta[gridzeta-1] << endl;
  for(i=0;i<gridpsi-1;i++)
    {
      psi[i+1]=psi[i]+dx;
    }
  for(i=0;i<gridzeta-1;i++)
    {
      zeta[i+1]=zeta[i]+dy;
    }
  
  for(i=0;i<grids;i++)
    s[i] = -0.5+(i+0.5)/grids;

  cout << "eps" << endl;
  cin >> eps;
  cout << "triangularity: (0 for circular)" << endl;
  cin >> delta;
  cout << "shift:  (0 for circular)" << endl;
  cin >> alpha;
  cout << "elongation. (1 for circular)" << endl;
  cin >> kappa;

  //This produces a reference toroidal grid
  threedtorusvtk(eps,delta,kappa,alpha);

  cout << "Data File Prefix?" << endl;
  cin >> prefix;
  cout << "What is the number of the first frame?" << endl;
  cin >> fframe;
  cout << "Interval between frames?" << endl;
  cin >> ninterval;
  cout << "Number of frames?" << endl;
  cin >> nframes;
  cout << "Magnetic shear?" << endl;
  cin >> shat;

  count=fframe;
  for(i=1;i<=nframes;i++)
    {
      stringstream strs;
      strs << setw(8) << setfill('0') << count;
      string aa = strs.str();
      string a = prefix + aa;
      cout << "Filename" << endl;
      str2 = new char[a.size()]; 
      strcpy(str2,a.c_str());
      cout << str2 << endl;

      ifstream file1 (str2, ios::in|ios::binary|ios::ate);
      if (file1.is_open())
	{
	  size = file1.tellg();
	  memblock = new double [size/sizeof(double)];
	  file1.seekg (0, ios::beg);
	  file1.read ((char*)memblock, size);
	  file1.close();
	  
	  cout << "The file content of " << str2 << " is in memory - it is " << size << "bytes long" << endl;
	  
	  //vtkfile(memblock,grids,gridzeta,gridpsi,count,Lx,Ly);
	  potentialtorusvtk(memblock,s,grids,psi,gridpsi,zeta,gridzeta,count,eps,q);
	  count = count + ninterval;
	  delete[] memblock;
	}else 
	{
	  cout << "Unable to open file -  All frames have been made....probably" << endl;
	}
    }
   
   return 0;
 }
  

void vtkfile(double* temp,int ns,int nx,int ny,int t,double Lx,double Ly)//Reconstructs the distribution function cross the magnetic field....This uses the structured points system of data
{
  //Note that the s coordinate is the fastest varying index, followed by radial modes and then poloidal
  int i,j,k,l;
  double dx,dy,ds;
  char filename[35];
  i=0;
  
  ds=Pi/ns;
  dx = Lx/nx;
  dy = Ly/ny;

  sprintf(filename,"3dPotential%04i.vtk",t);
  cout << "Writing file " << filename << endl;
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "Data for timestep" << " " << t << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_POINTS" << endl;
  outfile << "DIMENSIONS" << " " << ns << " " << nx << " " << ny << endl;
  outfile << "ASPECT_RATIO" << " " << ds << " " << dx << " " << dy << endl;
  outfile << "ORIGIN" << " " << 0 << " " << 0 << " " << 0 << endl;
  outfile << "POINT_DATA" << " " << ns*nx*ny << endl;
  outfile << "SCALARS" << " " << "Potential" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  l=0;
  for(k=0;k<nx;k++)
    {
      for(j=0;j<ny;j++)
	{
	  for(i=0;i<ns;i++)
	    {
	      outfile << temp[l] << " ";
	      l++;
	    } 
	}
      outfile << endl;
    }
  outfile.close();

}

void threedtorusvtk(double eps, double delta, double kappa, double alpha)
{
  //Outputs the distribution function as a function of vperp and mu along the s co-ordinate 
  int i,j,k,n,m;
  const double Pi = acos(-1.0);
  double theta,r;//The spherical coordinates
  double x,y,z;
  char filename[35];
  double phi;
  double Z = 0.0;
  int points = 100;
  double R = 1.0;
  
  sprintf(filename,"Torusmesh.vtk");
  cout << "Writing file " << filename << endl;
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "A reconstruction of the distribution function in velocity space" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_GRID" << endl;
  outfile << "DIMENSIONS" << " " << points+1 << " " << points+1 << " " << points+1 << endl;
  outfile << "POINTS" << " " << (points+1)*(points+1)*(points+1) << " " << "float" << endl;
  for(i=0;i<=points;i++)//Loop over theta and thi
    {
      for(j=0;j<=points;j++)
	{
	  for(k=0;k<=points;k++)
	    {
	      //rho = 0.2;//*i/(double)20;
	      theta = 2.0*Pi*j/(double)points;
	      phi = 2.0*Pi*k/(double)points;
	  
	      R=1.0+alpha+eps*cos(theta+asin(delta)*sin(theta));
	      

	      x = R*cos(phi);
	      y = R*sin(phi);
	      z = Z + kappa*eps*sin(theta);
	      
	      outfile << x << " " << y << " " << z << " ";
	    }
	  outfile << endl;
	}
    }
  outfile << "POINT_DATA" << " " << (points+1)*(points+1)*(points+1) << endl;
  outfile << "SCALARS" << " " << "Potential" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  
  for(k=0;k<=points;k++)
    {
      for(j=0;j<=points;j++)
	{
	  for(i=0;i<=points;i++)
	    {
	      outfile << 1.0*(i+j+k) << " ";
	    }
	  outfile << endl;
	}
    }
  outfile.close();

}

void potentialtorusvtk(double* temp,double* sgrid,int grids, double* psi,int gridx, double* chi,int gridy,int count,double eps,double q)
{
  //Outputs the distribution function as a function of vperp and mu along the s co-ordinate 
  int i,j,k,n,m,l,i2,jt;
  const double Pi = acos(-1.0);
  double theta,r;//The spherical coordinates

  double x,y,z;
  char filename[35];
  double R = 1.0;
  double phi,rho, tmp;
  double Z = 0.0;
  cout << count << endl;
  sprintf(filename,"TorusPotential%04i.vtk",count);
  cout << filename << endl;
  cout << "Writing file " << filename << endl;
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << filename << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_GRID" << endl;
  outfile << "DIMENSIONS" << " " << grids << " " << gridx << " " << repeats*gridy << endl;
  //outfile << "POINTS" << " " << (grids+1)*gridx*gridy << " " << "float" << endl;
  outfile << "POINTS" << " " << repeats*grids*gridx*gridy << " " << "float" << endl;

// The chi (zeta) coordinate is the axisymetric one that can be repeated
for (jt=0; jt<repeats;jt++)
{
    for(j=0;j<gridy;j++)
    {
      for(k=0;k<gridx;k++)
      {
        //for(i=0;i<grids+1;i++)
        for(i=0;i<grids;i++)
          {
            //Repeat last theta point for clean joins. 
            //if (i==grids+1) {i2=0;}
            //else {i2=i;}

            theta = 2.0*Pi*sgrid[i];
            // remove the half factor for inside / outside torus
            rho = eps+rho_star*(psi[k]-(Lx/2.));
            phi = q*theta*(1.0+shat*(rho-eps)) - 2.0*Pi*rho_star*(chi[j]-Ly/2+Ly*jt);
            
            x = R*cos(phi)*(1.0+alpha+rho*cos(theta+asin(delta)*sin(theta)));
            y = R*sin(phi)*(1.0+alpha+rho*cos(theta+asin(delta)*sin(theta)));
            z = Z + kappa*rho*sin(theta);

            outfile << x << " " << y << " " << z << " ";
          
          }
        outfile << endl;
      }
    }
}
  //outfile << "POINT_DATA" << " " << (grids+1)*gridx*gridy << endl;
  outfile << "POINT_DATA" << " " << repeats*grids*gridx*gridy << endl;
  outfile << "SCALARS" << " " << "Potential" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  for (jt=0; jt<repeats;jt++)
    {
    l=0;
    for(k=0;k<gridx;k++)
      {
        for(j=0;j<gridy;j++)
          {
            for(i=0;i<grids;i++)
              {
                //if (i==0) {tmp=temp[l];}
                outfile << temp[l] << " ";
                l++;
              } 
            //Repeat last theta point
            //outfile << tmp << " ";
            outfile << endl;
          }
      }
    }
}
