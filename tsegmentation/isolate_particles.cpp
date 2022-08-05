/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2009 Sergio Galindo                                    *
 * Copyright (C) 2013 William Oquendo                                   *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/

// Std lib
#include <math.h>

// GSL
// #include <gsl/gsl_linalg.h>

// // MechSys Included in plane_segmentation
// #include <mechsys/dem/domain.h>
// #include <mechsys/dem/distance.h>
// #include <mechsys/util/fatal.h>
// #include <mechsys/mesh/structured.h>
// #include <mechsys/mesh/unstructured.h>

//User
#include "plane_segmentation.h"
#include "strain_energy_field.h"

struct UserData
{
  Array<DEM::Particle *>    p;            // the array of particles at which the force is to be applied
  Array<Vec3_t  >    vm0;          // value of the vectors close to the middle section
  Array<Vec3_t *>    vm;           // pointers to the vectors close to the middle section
  String             test;         // Type of test vibraiton or tension
  double             A;            // Area of the plate for stress calculation
  double             Am;           // vibration amplitude
  double             ome;          // vibration frequency
  double             sz;           // Stress State
  double             Tf;           // Final time
  Vec3_t             L0;           // Initial dimensions
  Array<double>      strainEF;     // Array of strain energy field values for all particles in the domain
  std::ofstream      oss_ss;       // file for stress strain data
};

//Setup is called on EACH TIME STEP, this way we can calculate stresses on specific particles or impose boundary conditions
// Specifically Setup is called after building the particle interctions, this means we can't touch the particle interactions here but we can modify the forces, positions and velocities of particles
void Setup (DEM::Domain & Dom, void * UD)
{
  // std::cout<<"SETUP at time "<<Dom.Time<<"\n";
  //Use flexion as the test
  UserData & dat = (*static_cast<UserData *>(UD));
  //Average strain
  double avgFz=0.;
  for(size_t i=0; i<dat.p.Size();i++){
    avgFz+=dat.p[i]->F(2);//Add the foces in z for the moving plane
  }
  avgFz/=dat.p.Size();
  dat.sz = avgFz/dat.A;
}

// Report is called every dtOut steps before particle initialization (i.e. force assignments) so it uses the values of the timestep immediately before to calculate and report results.
void Report_test (DEM::Domain & Dom, void * UD)
{
  std::cout<<"REPORT at time "<<Dom.Time<<"\n";
  // UserData & dat = (*static_cast<UserData *>(UD));
  Array<Mat3_t> stresses = StressTensor(Dom.Interactons, Dom.BInteractons, Dom.Particles.Size());
  std::cout<<"Stress tensor for particle "<<0<<" : "<<stresses[0]<<"\n";
}

void Report (DEM::Domain & Dom, void * UD)
{
  std::cout<<"REPORT at time "<<Dom.Time<<"\n";
  UserData & dat = (*static_cast<UserData *>(UD));
  // Use flexion as the test
  //Calculate the strain energy field for all the particles in the domain
  Array<Mat3_t> stresses = StressTensor(Dom.Interactons, Dom.BInteractons, Dom.Particles.Size());
  std::cout<<"Stresses calculated for all particles! \n";
  Array<double> strainEF(Dom.Particles.Size());
  double max_strainEF = 0.;
  for(size_t p=0; p<Dom.Particles.Size();p++){ //Calculate strain energy field for each particle
    if (stresses[p](0,0)>10000. or stresses[p](0,0)<-10000.){
      std::cout<<"WARNING: Encountered large stress tensor components, printing interacrons.\n";
      std::cout<<"Particle "<<p<<" with stress tensor: \n "<<stresses[p]<<"\n";
      std::cout<<"Particle properties: \n";
      std::cout<<"Tag: "<<Dom.Particles[p]->Tag;
      std::cout<<"Kn: "<<Dom.Particles[p]->Props.Kn;
      std::cout<<"Kt: "<<Dom.Particles[p]->Props.Kt;
      std::cout<<"Gn: "<<Dom.Particles[p]->Props.Gn;
      std::cout<<"Gt: "<<Dom.Particles[p]->Props.Gt;
      std::cout<<"Mu: "<<Dom.Particles[p]->Props.Mu;
      std::cout<<"Bn: "<<Dom.Particles[p]->Props.Bn;
      std::cout<<"Bt: "<<Dom.Particles[p]->Props.Bt;
      std::cout<<"Bm: "<<Dom.Particles[p]->Props.Bm;
      std::cout<<"eps: "<<Dom.Particles[p]->Props.eps;
    }
    strainEF[p] = StrainEnergyField(stresses[p], /*Poisson's ratio*/0.1);
    if(strainEF[p] > max_strainEF) max_strainEF = strainEF[p];
  }
  std::cout<<"Strain energy field calculated! \n";
  std::cout<<"maximum strainEF: "<<max_strainEF<<"\n";
  String fn;
  fn.Printf    ("%s_%04d", "StrainEnergyField", Dom.idx_out);
  Dom.WriteXDMF_User(strainEF,fn.CStr());
  if (Dom.idx_out==0)
    {
      String fs;
      fs.Printf("%s_walls.res",Dom.FileKey.CStr());
      dat.oss_ss.open(fs.CStr());
      // Output of the current time, the stress state sx, and the strains ex,ey and ez
      dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sz_0" << Util::_8s << "sz_avg" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << Util::_8s << "max_strainEF" << std::endl;
      std::cout << Util::_10_6 << "Time" << Util::_8s << "sz_0" << Util::_8s << "sz_avg" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << Util::_8s << "max_strainEF" << std::endl;
    }
  if (!Dom.Finished)
    {
      // Measure the strains and stresses
      Vec3_t Xmin, Xmax;
      Dom.BoundingBox(Xmin, Xmax);
      double ex = (Xmax(0)-Xmin(0)-dat.L0(0))/dat.L0(0);
      double ey = (Xmax(1)-Xmin(1)-dat.L0(1))/dat.L0(1);
      double ez = (Xmax(2)-Xmin(2)-dat.L0(2))/dat.L0(2);
      double sz_0 = dat.p[0]->F(2)/dat.A;
      dat.oss_ss << Util::_10_6 << Dom.Time << Util::_8s << dat.sz  << Util::_8s << sz_0 << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << Util::_8s << max_strainEF << std::endl;
      std::cout << Util::_10_6 << Dom.Time << Util::_8s << dat.sz  << Util::_8s << sz_0 << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << Util::_8s << max_strainEF << std::endl;
      // ReportParticles(Dom);
    }
  else dat.oss_ss.close();
  // Check that no particle goes beyond the strain energy field threshold and write its index to a file if it does
  Array<size_t> particlesToBreakID(0);
  std::cout<<"Created particlesToBreakID array with size"<<particlesToBreakID.Size()<<"\n";
  double strainEFThr = max_strainEF*0.9; //3.0; // Strain Energy Field breaking threshold
  for(size_t i=0; i<strainEF.Size(); i++) if (strainEF[i]>strainEFThr) particlesToBreakID.Push(i);
  if (particlesToBreakID.Size() == 0){
    std::cout<<"No particles to break... yet\n";
  } else {
    std::cout<<"Particles surpassing breaking theshold:"<<particlesToBreakID<<"\n";
    String _fs;
    _fs.Printf("%s_break_index_%04f.res", "particles", Dom.Time);
    std::ofstream ofbreak;
    ofbreak.open(_fs.CStr());
    // Re calculate stresses since the matrices are broken in diagonalization procedure
    stresses = StressTensor(Dom.Interactons, Dom.Particles.Size());
    for(size_t pB=0; pB<particlesToBreakID.Size(); pB++){
      size_t pID = particlesToBreakID[pB];
      std::cout<<"Breaking particle with index "<<pID<<"\n";
      std::cout<<"Located at: "<<Dom.Particles[pID]->x<<"\n";
      // Build plane using the principal stress component and the particle geometric center
      // Calculate the first component of the stress tensor
      Vec3_t eigvalR = Vec3_t(0.,0.,0.), _ = Vec3_t(0.,0.,0.), stress_v0 = Vec3_t(0.,0.,0.);
      Mat3_t _stress = stresses[pID];//Create a new matrix with the stress tensor since EigNonsymm destroys the matrices it uses
      EigNonsymm(_stress, eigvalR, _, stress_v0, _, _, _, _, _);
      Vec3_t planeNormal = Vec3_t(0.,0.,0.), planeCentroid = Vec3_t(0.,0.,0.);
      planeNormal=stress_v0/norm(stress_v0);
      for(size_t v=0; v<Dom.Particles[pID]->Verts.Size(); v++) planeCentroid+=*(Dom.Particles[pID]->Verts[v]);
      planeCentroid/=Dom.Particles[pID]->Verts.Size();
      std::cout<<"Cutting plane with normal" << planeNormal<<" at "<<planeCentroid<<"\n";
      // Print particle number, strainEF, ID in the domain, cutting plane normal, cutting plane centre
      ofbreak << Util::_2 << pB << Util::_8s << strainEF[pID] << Util::_2 << pID << Util::_8s << planeNormal(0) << Util::_8s << planeNormal(1) << Util::_8s << planeNormal(2) << Util::_8s << planeCentroid(0) << Util::_8s << planeCentroid(1) << Util::_8s << planeCentroid(2) << std::endl;
    }
    ofbreak.close();
  }
  std::cout<<"REPORT finish!\n";
}

int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    size_t Nproc = 1;
    bool voronoi = 1;
    if (argc>=3) Nproc=atoi(argv[2]);
    if (argc==4) voronoi=atoi(argv[3]);
    std::cout<<"Loading from voronoi cell centers? "<<voronoi<<std::endl;
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double rho;         // rho
    size_t Restart;     // Restart number
    double max_time;    // Time at which the maximum strainEF was found in previous run, should match name of the restart file to use, if any HERE WE USE IT AS STARTNG TIME FOR THE NEW SIMULATION
    bool cohesion;      // Whether or not to apply cohesion between particles
    {
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> Restart;      infile.ignore(200,'\n');
        infile >> max_time;      infile.ignore(200,'\n');
        infile >> cohesion;      infile.ignore(200,'\n');
    }

    // user data and domain
    DEM::Domain dom;
    String _fs;
    dom.Dilate= true;
    // Tags of bulk, static, moving particles
    int bTag = -1, sTag =-2, mTag =-3;
    std::cout<<"Restart index:"<<Restart<<std::endl;
    std::cout<<"Applying cohesion?:"<<cohesion<<std::endl;
    if(Restart){
      std::cout<<"####### RESTART SCHEME #######"<<std::endl;
      String _fs;
      // _fs.Printf("run_%04d/%s_initial%04d", Restart-1, filekey.CStr(), Restart-1);
      // dom.Load(_fs.CStr());
      if (voronoi)
        {
          _fs.Printf("run_%04d/initial_points_%i.xyz", Restart-1, Restart-1);
          std::cout<<"Loading Voronoi packing from points file "<<_fs.CStr()<<"..."<<std::endl;
          dom.AddVoroPackFromPoints (/*Tag*/bTag, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz, /*FileKey*/_fs.CStr(), /*Density*/rho,  /*Cohesion*/cohesion,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
        }
      else
        {
          std::cout<<"Loading from saved geometry file\n";
          //// Load full geometry
          _fs.Printf("run_%04d/brick_restart_base_initial%04d", Restart-1,Restart-1);
          std::cout<<"Loading geometry from"<<_fs.CStr()<<std::endl;
          dom.LoadCohesion(_fs.CStr(), cohesion, bTag);
        }
      std::cout<<"Loaded initial domain structure!"<<std::endl;
      // Cut particles according to max strain EF
      // Previous run generated a file containing particle ids and the planes necessary to cut them so use that here.
      Array<DEM::Particle*> segParticles(0);
      Array<size_t> particlesToBreakID(0);
      size_t pID = 0; //Particle to cut
      _fs.Printf("run_%04d/particles_warning_index.res", Restart-1);
      // std::ifstream infile("run_0000/particles_break_index_5.600050.res");
      std::ifstream infile(_fs.CStr());
      DEM::Domain vdom;
      while (infile >> pID){
        std::cout<<"Writing visualization for particle: " << pID<<"\n";
        vdom.Particles.Push(dom.Particles[pID]);
        vdom.Particles.Last()-> Tag = pID;
      }
      vdom.WriteXDMF("particle_warning");
    }
}
MECHSYS_CATCH
