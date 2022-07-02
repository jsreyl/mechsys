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
  std::cout<<"SETUP at time "<<Dom.Time<<"\n";
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

void Report (DEM::Domain & Dom, void * UD)
{
  std::cout<<"REPORT at time "<<Dom.Time<<"\n";
  UserData & dat = (*static_cast<UserData *>(UD));
  Array<Mat3_t> stresses = StressTensor(Dom.Interactons, Dom.BInteractons, Dom.Particles.Size());
  std::cout<<"Stress tensor for particle "<<0<<" : "<<stresses[0]<<"\n";
  // if (dat.test=="tensile" || dat.test=="flexion")
  //   {
  //     if (Dom.idx_out==0)
  //       {
  //         String fs;
  //         fs.Printf("%s_walls.res",Dom.FileKey.CStr());
  //         dat.oss_ss.open(fs.CStr());
  //         double avgFz=0.;
  //         for(size_t i=0; i<dat.p.Size();i++){
  //           avgFz+=dat.p[i]->F(2);//Add the foces in z for the moving plane
  //         }
  //         avgFz/=dat.p.Size();
  //         dat.sz = avgFz/dat.A;
  //         // Output of the current time, the stress state sx, and the strains ex,ey and ez
  //         dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sz" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << std::endl;
  //       }
  //     if (!Dom.Finished)
  //       {
  //         // Measure the strains and stresses
  //         Vec3_t Xmin, Xmax;
  //         Dom.BoundingBox(Xmin, Xmax);
  //         double ex = (Xmax(0)-Xmin(0)-dat.L0(0))/dat.L0(0);
  //         double ey = (Xmax(1)-Xmin(1)-dat.L0(1))/dat.L0(1);
  //         double ez = (Xmax(2)-Xmin(2)-dat.L0(2))/dat.L0(2);
  //         dat.oss_ss << Util::_10_6 << Dom.Time << Util::_8s << dat.sz << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << std::endl;
  //       }
  //     else dat.oss_ss.close();
  //   }
  // if (dat.test=="bending")
  //   {
  //     if (Dom.idx_out==0)
  //       {
  //         double tol = dat.L0(2)/20.0;
  //         for (size_t i=0;i<Dom.Particles.Size();i++)
  //           {
  //             for (size_t j=0;j<Dom.Particles[i]->Verts.Size();j++)
  //               {
  //                 if (fabs((*Dom.Particles[i]->Verts[j])(2))<tol)
  //                   {
  //                     dat.vm.Push (Dom.Particles[i]->Verts[j]);
  //                     dat.vm0.Push(Vec3_t(*Dom.Particles[i]->Verts[j]));
  //                   }
  //               }
  //           }
  //       }
  //     String fs;
  //     fs.Printf("%s_%08d.res",Dom.FileKey.CStr(),Dom.idx_out);
  //     dat.oss_ss.open(fs.CStr());
  //     dat.oss_ss << Util::_10_6 << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s <<"ux" << Util::_8s << "uy" << Util::_8s << "uz" << std::endl;
  //     for (size_t i=0;i<dat.vm.Size();i++)
  //       {
  //         dat.oss_ss << Util::_8s <<            dat.vm0[i](0) << Util::_8s <<            dat.vm0[i](1) << Util::_8s <<            dat.vm0[i](2);
  //         dat.oss_ss << Util::_8s << (*dat.vm[i])(0)-dat.vm0[i](0) << Util::_8s << (*dat.vm[i])(1)-dat.vm0[i](1) << Util::_8s << (*dat.vm[i])(2)-dat.vm0[i](2) << std::endl;
  //       }
  //     dat.oss_ss.close();

  //   }
  // std::cout<<"REPORT finish!\n";
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
    dom.Alpha = verlet;
    dom.Dilate= true;
    // Tags of bulk, static, moving particles
    int bTag = -1, sTag =-2, mTag =-3;
    std::cout<<"Restart index:"<<Restart<<std::endl;
    std::cout<<"Applying cohesion?:"<<cohesion<<std::endl;

    //XXX: Removing the else ifs for tetra makes this do a free(): invalid pointer error, I have absolutely no idea why. You can try running voro_error to reproduce the error in a minimal way.
    //I've been chasing this error for half a year and atributing it to other things, but it's as misterious as it gets
    if (ptype=="voronoi")
      {
        dom.AddVoroPack (bTag, R, Lx,Ly,Lz, nx,ny,nz, rho,/*Cohhesion*/ cohesion, false, seed, 1.0);
        dom.SavePoints("initial_points.xyz", -1);
      }
    else if (ptype=="voronoi2")
      {
        dom.AddVoroPackFromPoints (/*Tag*/bTag, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz, /*FileKey*/"initial_points.xyz", /*Density*/rho,  /*Cohesion*/cohesion,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
      }
    else if (ptype=="tetra")
      {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
      }
    else throw new Fatal("Packing for particle type not implemented yet");
    // Specific to flexion test
    std::cout<<"Adding cylinder at "<<Vec3_t(0.,0.,Lz/2.+Rc+R)<<" with radious "<<Rc<<"\n";
    std::cout<<"Lz: "<<Lz<<" Rc: "<<Rc<<" R: "<<R<<"\n";
    //Add a cylinder as the connection of two circles of radius R0 and R1 located at X0 and X1
    dom.AddCylinder(/*Tag*/mTag,/*X0*/Vec3_t(0.,-Lz/2.,Lz/2.+R+Rc),/*R0*/Rc,/*X1*/Vec3_t(0.,Lz/2.,Lz/2.+R+Rc),/*R1*/Rc,/*R*/R,/*rho*/rho);
    dom.AddCylinder(/*Tag*/sTag,/*X0*/Vec3_t(-Lx/2.+dx,-Lz/2.,-Lz/2.-R-Rc),/*R0*/Rc,/*X1*/Vec3_t(-Lx/2.+dx,Lz/2.,-Lz/2.-R-Rc),/*R1*/Rc,/*R*/R,/*rho*/rho);
    //Make the third cylinder randomly smaller to give the fracture a predilect direction
    srand(0);
    double Rb=Rc-(rand()%10+1)*Rc/100.;
    std::cout <<"Rb: "<<Rb<<"\n";
    dom.AddCylinder(/*Tag*/sTag,/*X0*/Vec3_t(Lx/2.-dx,-Lz/2.,-Lz/2.-R-Rc),/*R0*/Rb,/*X1*/Vec3_t(Lx/2.-dx,Lz/2.,-Lz/2.-R-Rc),/*R1*/Rb,/*R*/R,/*rho*/rho);

    // Initialize the UserData structure
    dat.test = test;
    dat.A    = Lx*Ly;
    dat.Am   = Am;
    dat.ome  = ome;
    dat.Tf   = Tf;
    Vec3_t Xmin,Xmax;
    dom.BoundingBox(Xmin,Xmax);
    dat.L0   = Xmax - Xmin;

    //Generate planes at the bottom and top
    // dom.GenBoundingPlane(sTag,R,1.0,true);
    dom.AddPlane(/*Tag*/sTag,/*Position*/Vec3_t(0.0,0.0,Xmin(2)-R),/*R*/R,/*Lx*/1.5*dat.L0(0),/*Ly*/1.5*dat.L0(1),/*rho*/1.0,/*Angle*/M_PI,/*Axis*/&OrthoSys::e0);//Bottom plane
    dom.AddPlane(/*Tag*/mTag,/*Position*/Vec3_t(0.0,0.0,Xmax(2)+R),/*R*/R,/*Lx*/1.5*dat.L0(0),/*Ly*/1.5*dat.L0(1),/*rho*/1.0,/*Angle*/M_PI,/*Axis*/&OrthoSys::e0);//Top plane

    // input
    double cam_x=0.0, cam_y=2*Lx, cam_z=0.0;
    dom.CamPos = cam_x, cam_y, cam_z;

    //identify the moving lid
    // dat.p = dom.GetParticle (mTag);
    // if (test=="tensile")
    // {
    //     dat.p->FixVeloc();
    //     dat.p->v = 0.0, ex*dat.L0(1)/Tf, 0.0;
    // }
    dom.GetParticles(mTag,dat.p);//This is the top plane (and upper cylinder if the test is flexion)
    //if (test=="flexion"){
    for(size_t i=0; i<dat.p.Size();i++){//Set the top plane and upper cylinder moving down
      dat.p[i]->FixVeloc();
      // dat.p[i]->v = 0.0, 0.0, ex*dat.L0(2)/Tf;
      dat.p[i]->v = 0.0, 0.0, ex*Ly/Tf;
    }
    std::cout<<"Length of the system: "<<dat.L0(0)<<" "<<dat.L0(1)<<" "<<dat.L0(2)<<std::endl;
    std::cout<<"Length of the brick: "<<dat.L0(0)<<" "<<dat.L0(1)<<" "<<dat.L0(2)-4.*Rc<<std::endl;
    std::cout<<"Given velocity for top plane: "<<ex*Ly/Tf<<std::endl;
    std::cout<<"Velocity for top plane: "<<dat.p[0]->v<<std::endl;
    //}

    //set the element properties
    Dict B;
    B.Set(bTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(sTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(mTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    dom.SetProps(B);

    // fix -2 particles at the left extreme of the beam
    // fix static (tag -2) particles, this is the bottom plane (and the two lower cylinders if the test is flexion)
    Array<DEM::Particle *> p;
    dom.GetParticles(sTag,p);
    for(size_t i=0;i<p.Size();i++) p[i]->FixVeloc();

    std::cout<<"Number of interactons for domain:"<<dom.Interactons.Size()<<std::endl;
    std::cout<<"Number of BInteractons for domain:"<<dom.BInteractons.Size()<<std::endl;
    std::cout<<"Number of CInteractons for domain:"<<dom.CInteractons.Size()<<std::endl;
    // XXX: For some misterious reason calculating stress tensors here generates a free() allocation error on Cluster calculation. However if the Stress Tensor is calculated after Solve is initiated this error does not occurr. Weird, but it works for now.
    // Array<Mat3_t> stresses = StressTensor(dom.Interactons, dom.Particles.Size());
    // Array<Mat3_t> stresses = StressTensor(dom.Interactons, dom.BInteractons, dom.Particles.Size());
    // std::cout<<"Stress tensor for particle "<<0<<" : "<<stresses[0]<<"\n";
    // Array<double> strainEF(dom.Particles.Size());
    // for(size_t p=0;p<dom.Particles.Size();p++)
    //   strainEF[p] = StrainEnergyField(stresses[p], /*Poisson's ratio*/0.1);

    // dat.strainEF = strainEF;
    // std::cout<<"Strain Energy Field for particle "<<0<<" : "<<strainEF[0]<<"\n";

    dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
    // dom.Solve (Tf,dt,dtOut, NULL, NULL, filekey.CStr(), RenderVideo, Nproc);
}
MECHSYS_CATCH
