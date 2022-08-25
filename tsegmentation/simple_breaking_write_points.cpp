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

int main(int argc, char **argv) try
{
    if (argc<2) throw new Fatal("This program must be called with one argument: the name of the data input file without the '.inp' suffix.\nExample:\t %s filekey\n",argv[0]);
    String filekey  (argv[1]);
    String filename (filekey+".inp");
    size_t Nproc = 1;
    bool voronoi = 1;
    if (argc>=3) Nproc=atoi(argv[2]);
    if (argc>=4) voronoi=atoi(argv[3]);
    String geometry_file (argv[4]); //File containing the geometry to use for restarting
    String planecut_file (argv[5]); //File containing the ids of the particles to be cut and the cutting planes to do so
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
    DEM::Domain sdom; //Segmentation domain
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
        planeCentroid=dom.Particles[pID]->x;
        std::cout<<"Changing centroid to particle centre: " << planeCentroid<<"\n";
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
      // sdom.Alpha = verlet;
      sdom.Dilate = true;
      // sdom.Initialized = false; //Re initialize the particles <- NOTE: Sure we wanna do that?
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
    sdom.Solve (0.2,0.0005,1.0, NULL, NULL, filekey.CStr(), 0, Nproc);
    _fs.Printf("%s_final%04d", filekey.CStr(),Restart);
    std::cout<<"Saved final domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    sdom.Save(_fs.CStr());
    _fs.Printf("final_points_%i.xyz", Restart);
    sdom.SavePoints(_fs.CStr(), bTag);
}
MECHSYS_CATCH
