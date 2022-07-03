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
    if (argc==3) Nproc=atoi(argv[2]);
    if (!Util::FileExists(filename)) throw new Fatal("File <%s> not found",filename.CStr());
    ifstream infile(filename.CStr());

    double verlet;      // Verlet distance for optimization
    String ptype;       // Particle type 
    String test;       // Particle type 
    size_t RenderVideo; // Decide is video should be render
    double Kn;          // Normal stiffness
    double Kt;          // Tangential stiffness
    double Gn;          // Normal dissipative coefficient
    double Gt;          // Tangential dissipative coefficient
    double Mu;          // Microscopic friction coefficient
    double Bn;          // Cohesion normal stiffness
    double Bt;          // Cohesion tangential stiffness
    double Bm;          // Cohesion torque stiffness
    double Eps;         // Threshold for breking bonds
    double R;           // Spheroradius
    size_t seed;        // Seed of the ramdon generator
    double dt;          // Time step
    double dtOut;       // Time step for output
    double Tf;          // Final time for the test
    double Lx;          // Lx
    double Ly;          // Ly
    double Lz;          // Lz
    double dx;          // dx
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    double Rc;          // Radius for the cylinders
    double rho;         // rho
    double Am;          // vibration force amplitude
    double ome;         // Frequency of vibration
    double ex;          // Final strain for the tensile test (positive extension, negative compression)
    double SEFthr;      // Strain energy field threshold
    size_t Restart;     // Restart number
    double max_time;    // Time at which the maximum strainEF was found in previous run, should match name of the restart file to use, if any HERE WE USE IT AS STARTNG TIME FOR THE NEW SIMULATION
    size_t idx_init;    // Last idx from previous simulation
    bool cohesion;      // Whether or not to apply cohesion between particles
    {
        infile >> verlet;       infile.ignore(200,'\n');
        infile >> ptype;        infile.ignore(200,'\n');
        infile >> test;         infile.ignore(200,'\n');
        infile >> RenderVideo;  infile.ignore(200,'\n');
        infile >> Kn;           infile.ignore(200,'\n');
        infile >> Kt;           infile.ignore(200,'\n');
        infile >> Gn;           infile.ignore(200,'\n');
        infile >> Gt;           infile.ignore(200,'\n');
        infile >> Mu;           infile.ignore(200,'\n');
        infile >> Bn;           infile.ignore(200,'\n');
        infile >> Bt;           infile.ignore(200,'\n');
        infile >> Bm;           infile.ignore(200,'\n');
        infile >> Eps;          infile.ignore(200,'\n');
        infile >> R;            infile.ignore(200,'\n');
        infile >> seed;         infile.ignore(200,'\n');
        infile >> dt;           infile.ignore(200,'\n');
        infile >> dtOut;        infile.ignore(200,'\n');
        infile >> Tf;           infile.ignore(200,'\n');
        infile >> Lx;           infile.ignore(200,'\n');
        infile >> Ly;           infile.ignore(200,'\n');
        infile >> Lz;           infile.ignore(200,'\n');
        infile >> dx;           infile.ignore(200,'\n');
        infile >> nx;           infile.ignore(200,'\n');
        infile >> ny;           infile.ignore(200,'\n');
        infile >> nz;           infile.ignore(200,'\n');
        infile >> Rc;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> Am;           infile.ignore(200,'\n');
        infile >> ome;          infile.ignore(200,'\n');
        infile >> ex;           infile.ignore(200,'\n');
        infile >> SEFthr;       infile.ignore(200,'\n');
        infile >> Restart;      infile.ignore(200,'\n');
        infile >> max_time;      infile.ignore(200,'\n');
        infile >> idx_init;      infile.ignore(200,'\n');
        infile >> cohesion;      infile.ignore(200,'\n');
    }

    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    DEM::Domain sdom(&dat); //Segmentation domain
    String _fs;
    dom.Alpha = verlet;
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
      _fs.Printf("run_%04d/initial_points_%i.xyz", Restart-1, Restart-1);
      std::cout<<"Loading Voronoi packing from points file "<<_fs.CStr()<<"..."<<std::endl;
      dom.AddVoroPackFromPoints (/*Tag*/bTag, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz, /*FileKey*/_fs.CStr(), /*Density*/rho,  /*Cohesion*/cohesion,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
      std::cout<<"Loaded initial domain structure!"<<std::endl;
      // Cut particles according to max strain EF
      // Previous run generated a file containing particle ids and the planes necessary to cut them so use that here.
      Array<DEM::Particle*> segParticles(0);
      Array<size_t> particlesToBreakID(0);
      size_t _ = 0; //assign data we won't use
      size_t pID = 0; //Particle to cut
      double _strainEF = 0.;
      Vec3_t planeNormal = Vec3_t(0.,0.,0.), planeCentroid = Vec3_t(0.,0.,0.);
      _fs.Printf("run_%04d/particles_break_index_%04f.res", Restart-1, max_time);
      // std::ifstream infile("run_0000/particles_break_index_5.600050.res");
      std::ifstream infile(_fs.CStr());
      while (infile >> _ >> _strainEF >> pID >> planeNormal(0) >> planeNormal(1) >> planeNormal(2) >> planeCentroid(0) >> planeCentroid(1) >> planeCentroid(2)){
        particlesToBreakID.Push(pID);
        std::cout<<"Forming a plane with normal" << planeNormal<<" at "<<planeCentroid<<"\n";
        DEM::Domain vdom;
        //vdom.domID = 999999; //Visualization
        //vdom.domType = 1; //Subdomain
        vdom.Particles.Push(dom.Particles[pID]);
        vdom.Dilate = false;
        vdom.WriteXDMF("plane_no_cut_0");
        vdom.Dilate = true;
        vdom.WriteXDMF("plane_no_cut_dilate_0");
        double arg = dot(planeNormal, OrthoSys::e1);
        std::cout<<"n: "<<planeNormal<<", e0:"<<OrthoSys::e1<<", dot: "<<arg<<", acos(dot): "<<acos(arg)<<"\n";
        vdom.AddPlane(sTag, planeCentroid, 0., 5, 5, 1., /*Angle*/acos(arg), /*Axis*/&OrthoSys::e0);
        std::cout<<"Writing plane cut visualization \n";
        vdom.Dilate = false;
        vdom.WriteXDMF("plane_no_cut_1");
        vdom.Dilate = true;
        vdom.WriteXDMF("plane_no_cut_dilate_1");
        Array<Array<int>> intersectionIxs(0);
        size_t prev_segSize = segParticles.Size();
        std::cout<<"DOM0: before bisecting segParticles Size:"<<prev_segSize<<"\n";
        std::cout<<"DOM0: before bisecting intersectionIxs Size:"<<intersectionIxs.Size()<<"\n";
        BisectPolyhedron(dom.Particles[pID], planeNormal, planeCentroid, segParticles, intersectionIxs, /*mechsysErode*/4, 1e-3, true);
        size_t nParts = segParticles.Size() - prev_segSize;
        std::cout<<"Particle "<<pID<<" bisected into "<<nParts<<" parts. Creating new simulation domain...\n";
        // Make a visualization of just the particle and the plane cutting it
        vdom.Particles.Clear(); //same as Resize(0);
        for(size_t i=prev_segSize; i<segParticles.Size();i++){
          vdom.Particles.Push(segParticles[i]);
          vdom.Dilate = false;
          String _f_vis, _f_vis_dil;
          _f_vis.Printf    ("%s_%04d", "particle_cut", i-prev_segSize); _f_vis_dil.Printf    ("%s_%04d", "particle_cut_dilate", i-prev_segSize);
          vdom.WriteXDMF(_f_vis.CStr());
          vdom.Dilate = true;
          vdom.WriteXDMF(_f_vis_dil.CStr());
        }
        vdom.AddPlane(sTag, planeCentroid, 0., 5, 5, 1., /*Angle*/acos(arg), /*Axis*/&OrthoSys::e0);
        std::cout<<"Writing plane cut visualization \n";
        vdom.Dilate = false;
        vdom.WriteXDMF("particle_cut_2");
        vdom.Dilate = true;
        vdom.WriteXDMF("particle_cut_dilate_2");
        vdom.Particles.Clear(); // Clear particles before exiting so they don't get destroyed on domain class destructor
      }
      std::cout<<"Starting actual new simulation domain. \n";
      // sdom.domID = dom.domID + 1;
      // sdom.domType = 1; //Subdomain
      sdom.Alpha = verlet;
      sdom.Dilate = true;
      sdom.Initialized = false; //Re initialize the particles <- NOTE: Sure we wanna do that?
      double cam_x=0.0, cam_y=2*Lx, cam_z=0.0;
      sdom.CamPos = cam_x, cam_y, cam_z;
      //Copy particles from original domain to new segmented domain
      for(size_t p=0; p<dom.Particles.Size();p++){
        if(!particlesToBreakID.Has(p)){//If this particle was not segmented
          sdom.Particles.Push(dom.Particles[p]);//Just add it as it is
          sdom.Particles[sdom.Particles.Size()-1]->Index = sdom.Particles.Size()-1; //But reindex them so that they can be adressed properly
        }
      }
      for(size_t sp=0;sp<segParticles.Size();sp++) {
        sdom.Particles.Push(segParticles[sp]);
        sdom.Particles[sdom.Particles.Size()-1]->Index = sdom.Particles.Size()-1;
        // sdom.eigenParticles.Push(segParticles[sp]); // Make these be deleted on domain destruction
        // DON'T do this, this will create double free or corruption and subsequent segmentation fault
        // sdom.Particles[sdom.Particles.Size()-1]->Initialize(sdom.Particles.Size()-1);//Initialize this particle
        // sdom.Particles[sdom.Particles.Size()-1]->InitializeVelocity(dat.dt);//Initialize this particle
      }
      std::cout<<"DOM0: Number of particles in previous domain: "<<dom.Particles.Size()<<"\nNumber of particles in new domain: "<<sdom.Particles.Size()<<"\n";
      // Save new points for the next domain
      _fs.Printf("initial_points_%i.xyz", Restart);
      sdom.SavePoints(_fs.CStr(), bTag);
    }
    _fs.Printf("%s_%04d", "brick_geometry",Restart);
    sdom.WriteXDMF(_fs.CStr());
    sdom.WritePOV(_fs.CStr());
    std::cout<<"Saving initial domain structure..."<<std::endl;
    _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
    sdom.Save(_fs.CStr());
    std::cout<<"Saved initial domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    sdom.Solve (Tf,dt,dtOut, NULL, &Report, filekey.CStr(), RenderVideo, Nproc);
    _fs.Printf("%s_final%04d", filekey.CStr(),Restart);
    std::cout<<"Saved final domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    sdom.Save(_fs.CStr());
    _fs.Printf("final_points_%i.xyz", Restart);
    sdom.SavePoints(_fs.CStr(), bTag);
}
MECHSYS_CATCH
