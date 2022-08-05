/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2005 Dorival M. Pedroso, Raul Durand                   *
 * Copyright (C) 2009 Sergio Galindo                                    *
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

#include<random>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>

//User
//#include "plane_segmentationv7.h"
//#include "strain_energy_field.h"
//#include "user_report.h"


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
  size_t             idx_init;      // Last index from previous simulation
};

//Setup is called on EACH TIME STEP, this way we can calculate stresses on specific particles or impose boundary conditions
// Specifically Setup is called after building the particle interctions, this means we can't touch the particle interactions here but we can modify the forces, positions and velocities of particles

int main(int argc, char **argv) try
{
    std::cout<<"I do not understand why are you running the old code..THIS IS THE NEW CODE\n";
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
    size_t ZisUp;       // Is Z the upward axis?
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
    double dx;          // dx, size for brick walls
    size_t nx;          // nx
    size_t ny;          // ny
    size_t nz;          // nz
    size_t ndx;         // dnx
    double Rc;          // Radius for the cylinders
    double rho;         // rho
    double Am;          // vibration force amplitude
    double ome;         // Frequency of vibration
    double ex;          // Final strain for the tensile test (positive extension, negative compression)
    double tol1;        // Angle tolerance for cohesive bonds
    double tol2;        // Distance tolerance for cohesive bonds
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
        infile >> ZisUp;        infile.ignore(200,'\n');
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
        infile >> ndx;          infile.ignore(200,'\n');
        infile >> Rc;           infile.ignore(200,'\n');
        infile >> rho;          infile.ignore(200,'\n');
        infile >> Am;           infile.ignore(200,'\n');
        infile >> ome;          infile.ignore(200,'\n');
        infile >> ex;           infile.ignore(200,'\n');
        infile >> tol1;         infile.ignore(200,'\n');
        infile >> tol2;         infile.ignore(200,'\n');
        infile >> SEFthr;       infile.ignore(200,'\n');
        infile >> Restart;      infile.ignore(200,'\n');
        infile >> max_time;      infile.ignore(200,'\n');
        infile >> idx_init;      infile.ignore(200,'\n');
        infile >> cohesion;      infile.ignore(200,'\n');
    }

    // user data and domain
    //UserData dat;
    DEM::Domain dom; //(&dat);
    dom.Alpha = verlet;
    dom.Dilate= true; // Dilate for visualizations

    // Tags of bulk, static, moving particles
    int bTag = -1;
    //int bTag = -1, sTag =-2, mTag =-3;

    if(Restart){
      std::cout<<"####### RESTART SCHEME #######"<<std::endl;
      String _fs;
      _fs.Printf("run_%04d/%s_final%04d", Restart-1, filekey.CStr(), Restart-1);
      std::cout<<"Loading initial domain structure from "<<_fs.CStr()<<"..."<<std::endl;
      if (cohesion) {
        std::cout<<"Loading initial domain structure from points initial_points.xyz..."<<std::endl;
        // dom.AddVoroPackFromPoints (/*Tag*/-1, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz, /*FileKey*/"initial_points.xyz", /*Density*/rho,  /*Cohesion*/true,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
        //dom.LoadCohesion(_fs.CStr(), true, bTag);
        _fs.Printf    ("%s_%04d", "restart_voropoints", Restart-1);
        dom.WriteXDMF(_fs.CStr());
        // dom.SavePoints("restart_voropoints.xyz", -1);
      } else {
        dom.Load(_fs.CStr());
      }
      std::cout<<"Loaded initial domain structure!"<<std::endl;

      //No need to cut particles according to max strain EF
      // Just restart simulation from when it left off
      std::cout<<"Starting actual new simulation domain. \n";
      //dom.Initialized = false; //Re initialize the particles <- NOTE: Sure we wanna do that?
      double cam_x=0.0, cam_y=2*Lx, cam_z=0.0;
      dom.CamPos = cam_x, cam_y, cam_z;
      //std::cout<<"DOM"<<dom.domID<<": Number of particles in previous domain: "<<dom.Particles.Size()<<"\n";
      std::cout<<"DOM"<<0<<": Number of particles in previous domain: "<<dom.Particles.Size()<<"\n";

      // NOTE: ADD Cohesive interactions here!!! XXX
      /*
      if (cohesion){
        std::cout<<"Adding COHESION interactions between particles of the new domain\n";
        //LocalImposeParticleCohesion(bTag, dom);//Only add it for particles in the bulk (those tagged bTag), not to moving rods or planes
        dom.ImposeParticleCohesion(bTag);//Only add it for particles in the bulk (those tagged bTag), not to moving rods or planes
        std::cout<<"Successfully added COHESION interactions between particles of the new domain!\n";
        // Set some particles to not move, or move at a given speed
      }
      */
      Dict B;
      B.Set(bTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);//Breaking material
      dom.SetProps(B);
      //Run new simulation from here, this is recursive simulation
      String f_vis;
      f_vis.Printf    ("%s_%04d", "brick_planes", Restart-1);
      std::cout<<"Writing visualization files for dom.\n";
      dom.WriteXDMF(f_vis.CStr());
      dom.WritePOV(f_vis.CStr());
      std::cout<<"Saving initial domain structure..."<<std::endl;
      _fs.Printf("%s_initial%04d", filekey.CStr(),Restart-1);
      dom.Save(_fs.CStr());
      std::cout<<"Saved initial domain structure into "<<_fs.CStr()<<"!"<<std::endl;
      _fs.Printf("%s_restart%04d", filekey.CStr(),Restart-1);
      // Set initial and final time for new simulation
      dom.Time = max_time;
      // dat.idx_init = idx_init; //NOTE: Maybe we want the first index to be zero to ensure everything is initialized correctly
      //dom.Solve (Tf+max_time,dt,dtOut, &Setup, &Report, _fs.CStr(), RenderVideo, Nproc, dat.idx_init);//NOTE: Maybe we want the first index to be zero to ensure everything is initialized correctly
      dom.Solve (Tf+max_time,dt,dtOut, NULL, NULL, _fs.CStr(), RenderVideo, Nproc, idx_init);//NOTE: Maybe we want the first index to be zero to ensure everything is initialized correctly
      std::cout<<"Saving final domain structure..."<<std::endl;
      _fs.Printf("%s_final%04d", filekey.CStr(),Restart-1);
      dom.Save(_fs.CStr());
      std::cout<<"Saved final domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    }
    else{//Create particles as per usual
      if (ptype=="voronoi")
        {
          dom.AddVoroPack (/*Tag*/-1, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz,  /*Density*/rho,  /*Cohesion*/cohesion,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
          //dom.SavePoints("initial_points.xyz", -1);
        }
      else if (ptype=="cube")
        {
          Mesh::Structured mesh(3);
          mesh.GenBox (false, nx, ny, nz, Lx, Ly, Lz);
          dom.GenFromMesh (mesh, R, rho, true, false);
        }
      else if (ptype=="tetra")
        {
          Mesh::Unstructured mesh(/*NDim*/3);
          mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
          dom.GenFromMesh (mesh, R, rho, true, false);
        }
      else throw new Fatal("Packing for particle type not implemented yet");
      //set the element properties
      Dict B;
      B.Set(-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
      dom.SetProps(B);

      String _fs;
      _fs.Printf("%s%04d", "brick_planes_test",Restart);
      dom.WriteXDMF(_fs.CStr());
      dom.WritePOV(_fs.CStr());
      std::cout<<"Saving initial domain structure..."<<std::endl;
      _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
      dom.Save(_fs.CStr());
      std::cout<<"Saved initial domain structure into "<<_fs.CStr()<<"!"<<std::endl;
      //dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
      dom.Solve (10.,1.0e-4,0.1, NULL, NULL, filekey.CStr(), RenderVideo, Nproc);
    }
}
MECHSYS_CATCH
