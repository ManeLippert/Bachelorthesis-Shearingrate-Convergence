// This code W A Hornsby 2008
// Modified 2010  F J Casson
// For GKW data visualisation
// Takes GKW 3D Output and makes Visit VTK files.
// Files needed: 3DOutputparm.dat + 3D binary data files.
// Run GKW with at least 64 s points (128 is better) and den3d/ene3d/phi3d=.true.
// Compile with "[g/i]cc -lm -O2 Binaryconvertvtk.cpp" and run from the data folder.

// To Do: Add ability to create fake s points by linear interpolation.
// To Do: Add ability to repeat flux tube data
// To Do: Output not in ASCII format

#include <complex>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <iomanip>


void vtkfile(double* temp,int s,int x,int y,int t,double Lx,double Ly);
void threedtorusvtk(double eps);
void potentialtorusvtk(double* temp,double* sgrid,int grids, double* psi,int gridx, double* chi,int gridy,int count,double eps,double q);

using namespace std;

ifstream::pos_type size;
double * memblock;
double rho_star, shat; //ascetic parameters
double Lx,Ly,dx,dy;
const double Pi = 3.1416;

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

  dx = Lx/gridpsi;
  dy = Ly/gridzeta;

  double s[grids];
  double psi[gridpsi];
  double zeta[gridzeta];
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

  //This produces a reference toroidal grid
  threedtorusvtk(eps);

  cout << "Data File Prefix?" << endl;
  cin >> prefix;
  cout << "What is the number of the first frame?" << endl;
  cin >> fframe;
  cout << "Interval between frames?" << endl;
  cin >> ninterval;
  cout << "Number of frames?" << endl;
  cin >> nframes;
  cout << "Rho_star? 0.1/LX recommended:" << endl;
  cin >> rho_star;
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

void threedtorusvtk(double eps)
{
  //Outputs the distribution function as a function of vperp and mu along the s co-ordinate 
  int i,j,k,n,m;
  const double pi = acos(-1.0);
  double theta,r;//The spherical coordinates
  double x,y,z;
  char filename[35];
  double phi;
  double Z = 0.0;
  int points = 20;
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
	      theta = 2.0*pi*j/(double)points;
	      phi = 2.0*pi*k/(double)points;
	  
	      r=R*(1.0+eps*cos(theta));
	      
	      x = R*cos(phi)*(1.0+eps*cos(theta));
	      y = R*sin(phi)*(1.0+eps*cos(theta));
	      z = Z + eps*sin(theta);
	      
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
  int i,j,k,n,m,l;
  const double pi = acos(-1.0);
  double theta,r;//The spherical coordinates

  double x,y,z;
  char filename[35];
  double R = 1.0;
  double phi,rho;
  double Z = 0.0;
  cout << count << endl;
  sprintf(filename,"TorusPotential%04i.vtk",count);
  cout << filename << endl;
  cout << "Writing file " << filename << endl;
  ofstream outfile;
  outfile.open(filename);
  outfile << "# vtk DataFile Version 2.0" << endl;
  outfile << "A reconstruction of the distribution function in velocity space" << endl;
  outfile << "ASCII" << endl;
  outfile << "DATASET"<< " ";
  outfile << "STRUCTURED_GRID" << endl;
  outfile << "DIMENSIONS" << " " << grids << " " << gridx << " " << gridy << endl;
  outfile << "POINTS" << " " << grids*gridx*gridy << " " << "float" << endl;

    for(j=0;j<gridy;j++)
      {
        for(k=0;k<gridx;k++)
        {
	  for(i=0;i<grids;i++)
	    {
	      theta = 2.0*pi*sgrid[i];
              // remove the half factor for inside / outside torus
	      rho = eps+rho_star*(psi[k]-(Lx/2.));
              phi = q*theta*(1.0+shat*(rho-eps)) - 2.0*pi*rho_star*(chi[j]-Ly/2);
	      
	      x = R*cos(phi)*(1.0+rho*cos(theta));
	      y = R*sin(phi)*(1.0+rho*cos(theta));
	      z = Z + rho*sin(theta);
	      
	      outfile << x << " " << y << " " << z << " ";
	    
	    }
	  outfile << endl;
	}
    }
  outfile << "POINT_DATA" << " " << grids*gridx*gridy << endl;
  outfile << "SCALARS" << " " << "Potential" << " " << "float" << endl;  
  outfile << "LOOKUP_TABLE" << " " << "default" << endl;
  l=0;
  for(k=0;k<gridx;k++)
    {
      for(j=0;j<gridy;j++)
	{
	  for(i=0;i<grids;i++)
	    {
	      outfile << temp[l] << " ";
	      l++;
	    } 
	  outfile << endl;
	}
    }

}
