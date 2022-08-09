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
#include "strain_energy_fieldv2.h"

struct UserData
{
  Array<DEM::Particle *>    p;            // the array of particles at which the force is to be applied
  Array<DEM::Particle *>    ps;           // the array of particles that should be static at every moment
  Array<DEM::Particle *>    pt;           // the array of particles whose forces we keep track of
  Array<Vec3_t  >    vm0;          // value of the vectors close to the middle section
  Array<Vec3_t *>    vm;           // pointers to the vectors close to the middle section
  String             test;         // Type of test vibraiton or tension
  double             A;            // Area of the plate for stress calculation
  double             Am;           // vibration amplitude
  double             ome;          // vibration frequency
  double             sz;           // Stress State
  double             vz;           // Velocity for compression in flexural strength test
  double             Tf;           // Final time
  double             Tcomp;        // Maximum time for compression; after this stop compressing and equilibrate
  Vec3_t             L0;           // Initial dimensions
  int                bTag;         // Bulk particle tag
  int                sTag;         // Static particle tag
  int                mTag;         // Mobile particle tag
  Array<double>      strainEF;     // Array of strain energy field values for all particles in the domain
  std::ofstream      oss_ss;       // file for stress strain data
  bool               verbose;      // Whether to print check statements or not
  bool               compress;// Perform flexural test by compressing with the upper plane or not
};

//Setup is called on EACH TIME STEP, this way we can calculate stresses on specific particles or impose boundary conditions
// Specifically Setup is called after building the particle interctions, this means we can't touch the particle interactions here but we can modify the forces, positions and velocities of particles
void Setup (DEM::Domain & Dom, void * UD)
{
  // std::cout<<"SETUP at time "<<Dom.Time<<"\n";
  //Use flexion as the test
  UserData & dat = (*static_cast<UserData *>(UD));
  //if((Dom.Time > dat.Tcomp) && dat.compress){ //Stop moving the plane after 1.5 seconds (approximately 25micrometers from the start)
  // if(Dom.Time < dat.Tcomp && !dat.compress){ //Stop moving the plane after 1.5 seconds (approximately 25micrometers from the start)
  //   for(size_t i=0; i<dat.p.Size();i++){//Set the top plane and upper cylinder moving down
  //     dat.p[i]->FixVeloc(0.0, 0.0, dat.vz);
  //     dat.compress = true;
  //     // std::cout<<'a'<<dat.p[i]->v<<std::endl;
  //     // dat.p[i]->v = 0.0, 0.0, dat.vz;
  //     // std::cout<<'b'<<dat.p[i]->v<<std::endl;
  //   }
  // } else if (Dom.Time > dat.Tcomp && dat.compress) {
  //   if(dat.compress) std::cout<<"SETUP: Stopping plane with velocity"<<dat.p[0]->v<<std::endl;
  //   for(size_t i=0; i<dat.p.Size();i++){//Set the top plane and upper cylinder moving down
  //     // dat.p[i]->v = 0.0, 0.0, ex*dat.L0(2)/Tf;
  //     dat.p[i]->v = 0.0, 0.0, 0.0;
  //     std::cout<<'a'<<dat.p[i]->v<<std::endl;
  //     dat.p[i]->FixVeloc();
  //     std::cout<<'b'<<dat.p[i]->v<<std::endl;
  //   }
  //   if(dat.compress) std::cout<<"SETUP: Moving plane velocity after stopping"<<dat.p[0]->v<<std::endl;
  //   if(dat.compress) std::cout<<"SETUP: Final position of the top plane 0: "<<dat.p[0]->x<<std::endl;
  //   if(dat.compress) std::cout<<"SETUP: Final position of the top plane 1: "<<dat.p[1]->x<<std::endl;
  //   dat.compress = false;
  // }
  std::cout<<'a'<<dat.p[0]->v<<std::endl;
  std::cout<<'b'<<dat.p[1]->v<<std::endl;
  /*
  //Average strain
  double avgFz=0.;
  for(size_t i=0; i<dat.p.Size();i++){
    avgFz+=dat.p[i]->F(2);//Add the foces in z for the moving plane
  }
  avgFz/=dat.p.Size();
  dat.sz = avgFz/dat.A;
  */
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
  std::cout<<"Velocity for top plane: "<<dat.p[1]->v<<std::endl;
  std::cout<<"Position for top plane: "<<dat.p[1]->x<<std::endl;
  // Use flexion as the test
  //Calculate the strain energy field for all the particles in the domain
  Array<Mat3_t> stresses = StressTensor(Dom.Interactons, Dom.BInteractons, Dom.Particles.Size());
  std::cout<<"Stresses calculated for all particles! \n";
  Array<double> strainEF(Dom.Particles.Size());
  double max_strainEF = 0.;
  // Array<double> max_strains(Dom.Particles.Size()/10);
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
      std::cout<<"Gv: "<<Dom.Particles[p]->Props.Gv;
      std::cout<<"Gm: "<<Dom.Particles[p]->Props.Gm;
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
      dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sz_0" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << Util::_8s << "max_strainEF" << std::endl;
      std::cout << Util::_10_6 << "Time" << Util::_8s << "sz_0" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << Util::_8s << "max_strainEF" << std::endl;
      std::cout <<"Upper cylinder:"<< Util::_2 << "ID" << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "strainEF" << Util::_2 << "Nc" << std::endl;
      std::cout <<"Upper plane:"<< Util::_2 << "ID" << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "strainEF" << Util::_2 << "Nc" << std::endl;
      std::cout <<"Lower cylinder a:"<< Util::_2 << "ID" << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "strainEF"<< Util::_2 << "Nc" << std::endl;
      std::cout <<"Lower cylinder b:"<< Util::_2 << "ID" << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << Util::_8s << "strainEF"<< Util::_2 << "Nc" << std::endl;
      std::cout <<"Lower plane:"<< Util::_2 << "ID" << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "strainEF" << Util::_2 << "Nc" << std::endl;
      std::cout <<"Upper particle:"<< Util::_2 << "ID" << Util::_10_6 << "Time" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "BFx" << Util::_8s << "BFy" << Util::_8s << "BFz" << Util::_8s << "CFx" << Util::_8s << "CFy" << Util::_8s << "CFz" << Util::_8s << "strainEF" << Util::_2 << "Nc" << std::endl;
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
      Vec3_t F0 = dat.p[0]->F;
      dat.oss_ss << Util::_10_6 << Dom.Time << Util::_8s << sz_0  << Util::_8s << F0(0) << Util::_8s << F0(1) << Util::_8s << F0(2) << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << Util::_8s << max_strainEF << std::endl;
      std::cout << Util::_10_6 << Dom.Time << Util::_8s << sz_0  << Util::_8s << F0(0) << Util::_8s << F0(1) << Util::_8s << F0(2) << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << Util::_8s << max_strainEF << std::endl;
      // Upper cylinder and plane
      Vec3_t xp = dat.p[0]->x, vp = dat.p[0]->v, Fp = dat.p[0]->F;
      size_t pID = dat.p[0]->Index;
      double sp = strainEF[pID];
      size_t Nc = CalculateContacts(Dom.Interactons, Dom.BInteractons, pID);
      std::cout <<"Upper cylinder:"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << sp << Util::_2 << Nc << std::endl;
      xp = dat.p[1]->x; vp = dat.p[1]->v; Fp = dat.p[1]->F; pID = dat.p[1]->Index; sp = strainEF[pID]; Nc = CalculateContacts(Dom.Interactons, Dom.BInteractons, pID);
      std::cout <<"Upper plane:"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << sp << Util::_2 << Nc << std::endl;
      //Lower cylinders and plane
      xp = dat.ps[0]->x; vp = dat.ps[0]->v; Fp = dat.ps[0]->F; pID = dat.ps[0]->Index; sp = strainEF[pID]; Nc = CalculateContacts(Dom.Interactons, Dom.BInteractons, pID);
      std::cout <<"Lower cylinder a:"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << sp << Util::_2 << Nc << std::endl;
      xp = dat.ps[1]->x; vp = dat.ps[1]->v; Fp = dat.ps[1]->F; pID = dat.ps[1]->Index; sp = strainEF[pID]; Nc = CalculateContacts(Dom.Interactons, Dom.BInteractons, pID);
      std::cout <<"Lower cylinder b:"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << sp << Util::_2 << Nc << std::endl;
      xp = dat.ps[2]->x; vp = dat.ps[2]->v; Fp = dat.ps[2]->F; pID = dat.ps[2]->Index; sp = strainEF[pID]; Nc = CalculateContacts(Dom.Interactons, Dom.BInteractons, pID);
      std::cout <<"Lower plane:"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << sp << Util::_2 << Nc << std::endl;
      // Particle to track
      xp = dat.pt[0]->x; vp = dat.pt[0]->v; Fp = dat.pt[0]->F; pID = dat.pt[0]->Index; sp = strainEF[pID]; Nc = CalculateContacts(Dom.Interactons, Dom.BInteractons, pID);
      Vec3_t BFp = CalculateForce(Dom.BInteractons, pID);
      Vec3_t CFp = CalculateForce(Dom.CInteractons, pID);
      std::cout <<"Upper particle:"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << BFp(0) << Util::_8s << BFp(1) << Util::_8s << BFp(2) << Util::_8s << CFp(0) << Util::_8s << CFp(1) << Util::_8s << CFp(2) << Util::_8s << sp << Util::_2 << Nc << std::endl;
    }
  else dat.oss_ss.close();
  if (Dom.Time > dat.Tcomp){ // After the compression start looking into how to break the particles
    if (dat.compress) {
      std::cout<<"REPORT: Stopping plane with velocity"<<dat.p[0]->v<<std::endl;
      for(size_t i=0; i<dat.p.Size();i++){//Set the top plane and upper cylinder moving down
        // dat.p[i]->v = 0.0, 0.0, ex*dat.L0(2)/Tf;
        dat.p[i]->v = 0.0, 0.0, 0.0;
        dat.p[i]->xb = dat.p[i]->x; //Force velocities in the Verlet algorithm to be zero by setting particle displacement to zero
        //dat.p[i]->FixVeloc();
        //std::cout<<'b'<<dat.p[i]->v<<std::endl;
      }
      std::cout<<"REPORT: Moving plane velocity after stopping"<<dat.p[1]->v<<std::endl;
      std::cout<<"REPORT: Final position of the top plane: "<<dat.p[1]->x<<std::endl;
      dat.compress = false;
    }
    // Check that no particle goes beyond the strain energy field threshold and write its index to a file if it does
    // NOTE: Select the top 10% of the particles to break
    Array<size_t> particlesToBreakID(0);
    if(dat.verbose) std::cout<<"Created particlesToBreakID array with size"<<particlesToBreakID.Size()<<"\n";
    // Array<double> max_strains = strainEF;
    // max_strains.Sort(); // sort  strainEFs in ascending order
    // if (dat.verbose) std::cout<<"Strain EF: "<<strainEF<<"\n";
    // if (dat.verbose) std::cout<<"Max strains: "<<max_strains<<"\n";
    // double strainEFThr = 0.1*max_strains.Last(); //3.0; // Strain Energy Field breaking threshold
    // double strainEFThr = max_strains[9*Dom.Particles.Size()/10]; //3.0; // Strain Energy Field breaking threshold
    // std::cout<<"Strain energy field threshold for 90th quantile"<<strainEFThr<<"\n";
    double strainEFThr = max_strainEF*0.1; //3.0; // Strain Energy Field breaking threshold
    for(size_t i=0; i<strainEF.Size(); i++) if (strainEF[i]>strainEFThr && Dom.Particles[i]->Tag == dat.bTag) particlesToBreakID.Push(i); //Add only bulk particles
    if (particlesToBreakID.Size() == 0){
      std::cout<<"No particles to break... yet\n";
    } else {
      std::cout<<"Particles surpassing breaking theshold:"<<particlesToBreakID<<"\n";
      String _fs;
      _fs.Printf("%s_break_index_%04f.res", "particles", Dom.Time);
      std::ofstream ofbreak;
      ofbreak.open(_fs.CStr());
      // Re calculate stresses since the matrices are broken in diagonalization procedure
      // stresses = StressTensor(Dom.Interactons, Dom.Particles.Size());
      stresses = StressTensor(Dom.Interactons, Dom.BInteractons, Dom.Particles.Size());
      for(size_t pB=0; pB<particlesToBreakID.Size(); pB++){
        size_t pID = particlesToBreakID[pB];
        std::cout<<"Breaking particle with index "<<pID<<"\n";
        if(dat.verbose) std::cout<<"Located at: "<<Dom.Particles[pID]->x<<"\n";
        // Build plane using the principal stress component and the particle geometric center
        // Calculate the first component of the stress tensor
        Vec3_t eigvalR = Vec3_t(0.,0.,0.), _ = Vec3_t(0.,0.,0.), stress_v0 = Vec3_t(0.,0.,0.), stress_v1 = Vec3_t(0.,0.,0.), stress_v2 = Vec3_t(0.,0.,0.);
        Mat3_t _stress = stresses[pID];//Create a new matrix with the stress tensor since EigNonsymm destroys the matrices it uses
        EigNonsymm(_stress, eigvalR, _, stress_v0, stress_v1, stress_v2, _, _, _);
        Array<Vec3_t> stress_vs(3);
        stress_vs[0] = stress_v0; stress_vs[1] = stress_v1; stress_vs[2] = stress_v2;
        if(dat.verbose) std::cout<<"CHECK: Stress tensor with eigenvalues "<<eigvalR<<" and eigenvectors "<<stress_vs<<std::endl;
        // EigNonsymm returns unordered eigenvalues and eigenvectors, we want sigma_1>sigma_2>sigma_3
        size_t largest_id = 0;
        double largest_stress = eigvalR(0);
        for(size_t i=1; i<3;i++) {
          if (eigvalR(i) > largest_stress){
            largest_stress = eigvalR(i); largest_id = i;
          }
        }
        if(dat.verbose) std::cout<<"CHECK: Largest stress eigenvalue "<<eigvalR(largest_id)<<" and eigenvector "<<stress_vs[largest_id]<<std::endl;
        if(eigvalR(largest_id)<0.) std::cout<<"WARNING: All eigenvalues are negative for particle "<<pID<<std::endl;
        Vec3_t planeNormal = Vec3_t(0.,0.,0.), planeCentroid = Vec3_t(0.,0.,0.);
        planeNormal=stress_vs[largest_id]/norm(stress_vs[largest_id]);
        for(size_t v=0; v<Dom.Particles[pID]->Verts.Size(); v++) planeCentroid+=*(Dom.Particles[pID]->Verts[v]);
        planeCentroid/=Dom.Particles[pID]->Verts.Size();
        if(dat.verbose) std::cout<<"Cutting plane with normal" << planeNormal<<" at "<<planeCentroid<<"\n";
        // Print particle number, strainEF, ID in the domain, cutting plane normal, cutting plane centre
        ofbreak << Util::_2 << pB << Util::_8s << strainEF[pID] << Util::_2 << pID << Util::_8s << planeNormal(0) << Util::_8s << planeNormal(1) << Util::_8s << planeNormal(2) << Util::_8s << planeCentroid(0) << Util::_8s << planeCentroid(1) << Util::_8s << planeCentroid(2) << std::endl;
      }
      ofbreak.close();
    }
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
    // double Gv;          // Linear velocity dissipative coefficient
    // double Gm;          // Linear torque dissipative coefficient
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
        // infile >> Gv;           infile.ignore(200,'\n');
        // infile >> Gm;           infile.ignore(200,'\n');
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
    String _fs;
    dom.Alpha = verlet;
    dom.Dilate= true;
    dat.Tcomp = max_time;
    dat.verbose = false;
    // Tags of bulk, static, moving particles
    int bTag = -1, sTag =-2, mTag =-3;
    std::cout<<"Restart index:"<<Restart<<std::endl;
    std::cout<<"Applying cohesion?:"<<cohesion<<std::endl;

    //XXX: Removing the else ifs for tetra makes this do a free(): invalid pointer error, I have absolutely no idea why. You can try running voro_error to reproduce the error in a minimal way.
    //I've been chasing this error for half a year and atributing it to other things, but it's as misterious as it gets
    if (ptype=="voronoi")
      {
        dom.AddVoroPack (bTag, R, Lx,Ly,Lz, nx,ny,nz, rho,/*Cohhesion*/ cohesion, false, seed, 1.0);
        _fs.Printf("initial_points_%i.xyz", Restart);
        dom.SavePoints(_fs.CStr(), bTag);
        _fs.Printf("%s_%04d", "brick_geometry",Restart);
        dom.WriteXDMF(_fs.CStr());
        dom.WritePOV(_fs.CStr());
        std::cout<<"Saving initial domain structure..."<<std::endl;
        _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
        dom.Save(_fs.CStr());
        std::cout<<"Saved initial domain structure into"<<_fs.CStr()<<std::endl;
      }
    else if (ptype=="voronoi2")
      {
        _fs.Printf("initial_points_%i.xyz", Restart);
        std::cout<<"Reading points for voronoi packing from "<<_fs.CStr()<<"\n";
        dom.AddVoroPackFromPoints (/*Tag*/bTag, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz, /*FileKey*/_fs.CStr(), /*Density*/rho,  /*Cohesion*/cohesion,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
        _fs.Printf("%s_%04d", "brick_geometry",Restart);
        dom.WriteXDMF(_fs.CStr());
        std::cout<<"Saved initial domain structure into"<<_fs.CStr()<<std::endl;
        _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
        dom.Save(_fs.CStr());
        std::cout<<"Saved initial domain structure into"<<_fs.CStr()<<std::endl;
      }
    else if (ptype=="load")
      {
        _fs.Printf("segment_%04d/brick_restart_geometry_initial%04d", Restart-1,Restart);
        // _fs.Printf("brick_geometry_%i", Restart);
        std::cout<<"Loading geometry from "<<_fs.CStr()<<"\n";
        dom.LoadCohesion(_fs.CStr(),cohesion,bTag);
        std::cout<<"Number of particles added: "<<dom.Particles.Size()<<std::endl;
        _fs.Printf("%s_%04d", "brick_geometry",Restart);
        dom.WriteXDMF(_fs.CStr());
        std::cout<<"Saved initial domain structure into"<<_fs.CStr()<<std::endl;
        _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
        dom.Save(_fs.CStr());
        std::cout<<"Saved initial domain structure into"<<_fs.CStr()<<std::endl;
      }
    else if (ptype=="tetra")
      {
        Mesh::Unstructured mesh(/*NDim*/3);
        mesh.GenBox  (/*O2*/false,/*V*/Lx*Ly*Lz/(0.5*nx*ny*nz),Lx,Ly,Lz);
        dom.GenFromMesh (mesh, R, rho, true, false);
      }
    else
      {
        throw new Fatal("Packing for particle type not implemented yet");
      }
    // Specific to flexion test
    std::cout<<"Adding cylinder at "<<Vec3_t(0.,0.,Lz/2.+Rc+R)<<" with radious "<<Rc<<"\n";
    std::cout<<"Lz: "<<Lz<<" Rc: "<<Rc<<" R: "<<R<<"\n";
    //Add a cylinder as the connection of two circles of radius R0 and R1 located at X0 and X1
    dom.AddCylinder(/*Tag*/mTag,/*X0*/Vec3_t(0.,-Lz/2.,Lz/2.+R+Rc),/*R0*/Rc,/*X1*/Vec3_t(0.,Lz/2.,Lz/2.+R+Rc),/*R1*/Rc,/*R*/R,/*rho*/rho);
    std::cout<<"Upper cylinder with index "<<dom.Particles.Size()-1<<" located at "<<dom.Particles.Last()->x<<std::endl;
    dom.AddCylinder(/*Tag*/sTag,/*X0*/Vec3_t(-Lx/2.+dx,-Lz/2.,-Lz/2.-R-Rc),/*R0*/Rc,/*X1*/Vec3_t(-Lx/2.+dx,Lz/2.,-Lz/2.-R-Rc),/*R1*/Rc,/*R*/R,/*rho*/rho);
    std::cout<<"Lower cylinder with index "<<dom.Particles.Size()-1<<" located at "<<dom.Particles.Last()->x<<std::endl;
    //Make the third cylinder randomly smaller to give the fracture a predilect direction
    srand(0);
    double Rb=Rc; //-(rand()%10+1)*Rc/100.;
    std::cout <<"Rb: "<<Rb<<"\n";
    dom.AddCylinder(/*Tag*/sTag,/*X0*/Vec3_t(Lx/2.-dx,-Lz/2.,-Lz/2.-R-Rc),/*R0*/Rb,/*X1*/Vec3_t(Lx/2.-dx,Lz/2.,-Lz/2.-R-Rc),/*R1*/Rb,/*R*/R,/*rho*/rho);
    std::cout<<"Lower cylinder (Radius Rb) with index "<<dom.Particles.Size()-1<<" located at "<<dom.Particles.Last()->x<<std::endl;

    // Initialize the UserData structure
    dat.test = test;
    dat.A    = Lx*Ly;
    dat.Am   = Am;
    dat.ome  = ome;
    dat.Tf   = Tf;
    dat.bTag = bTag; dat.sTag = sTag; dat.mTag = mTag;
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
      // dat.p[i]->v = 0.0, 0.0, ex*Lz/dat.Tcomp;
      dat.p[i]->v = 0.0, 0.0, ex;
      dat.p[i]->Props.Gv = 0; //Cilinders and lower plane should not be affected by the viscous force
      dat.p[i]->Props.Gm = 0;
    }
    dat.vz = ex;
    std::cout<<"Indices of upper cylinder and plane: "<<dat.p[0]->Index<<" "<<dat.p[1]->Index<<std::endl;
    dat.compress = true; // Start the simulation compressing the system as in a flexural strngth test
    std::cout<<"Length of the system: "<<dat.L0(0)<<" "<<dat.L0(1)<<" "<<dat.L0(2)<<std::endl;
    std::cout<<"Length of the brick: "<<dat.L0(0)<<" "<<dat.L0(1)<<" "<<dat.L0(2)-4.*Rc<<std::endl;
    // std::cout<<"Given velocity for top plane: "<<ex*Lz/Tf<<std::endl;
    std::cout<<"Given velocity for top plane: "<<ex<<std::endl;
    std::cout<<"Velocity for top plane: "<<dat.p[0]->v<<std::endl;
    std::cout<<"Starting position of the top plane: "<<dat.p[1]->x<<std::endl;
    std::cout<<"Area of the brick: "<<dat.A<<std::endl;
    //}

    //set the element properties
    // double Gv = 1.0; double Gm = 1.0;
    Dict B;
    // B.Set(bTag,"Kn Kt Bn Bt Bm Gn Gt Gv Gm Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Gv ,Gm ,Mu ,Eps);
    // B.Set(sTag,"Kn Kt Bn Bt Bm Gn Gt Gv Gm Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Gv ,Gm ,Mu ,Eps);
    // B.Set(mTag,"Kn Kt Bn Bt Bm Gn Gt Gv Gm Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Gv ,Gm ,Mu ,Eps);
    B.Set(bTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(sTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(mTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    dom.SetProps(B);

    // fix -2 particles at the left extreme of the beam
    // fix static (tag -2) particles, this is the bottom plane (and the two lower cylinders if the test is flexion)
    //Array<DEM::Particle *> p;
    dom.GetParticles(sTag,dat.ps);
    for(size_t i=0;i<dat.ps.Size();i++) {
      dat.ps[i]->FixVeloc();
      dat.ps[i]->Props.Gv = 0; //Cilinders and lower plane should not be affected by the viscous force
      dat.ps[i]->Props.Gm = 0;
    }
    std::cout<<"Indices of lower cylinders and plane: "<<dat.ps[0]->Index<<" "<<dat.ps[1]->Index<<" "<<dat.ps[2]->Index<<std::endl;

    std::cout<<"Number of interactons for domain:"<<dom.Interactons.Size()<<std::endl;
    std::cout<<"Number of BInteractons for domain:"<<dom.BInteractons.Size()<<std::endl;
    std::cout<<"Number of CInteractons for domain:"<<dom.CInteractons.Size()<<std::endl;
    // XXX : Keeping track of particles forces
    dat.pt.Push(dom.Particles[204]); //Keep track of the particle with index 204 (Touching the upper cylinder
    // XXX: For some misterious reason calculating stress tensors here generates a free() allocation error on Cluster calculation. However if the Stress Tensor is calculated after Solve is initiated this error does not occurr. Weird, but it works for now.
    // Array<Mat3_t> stresses = StressTensor(dom.Interactons, dom.Particles.Size());
    // Array<Mat3_t> stresses = StressTensor(dom.Interactons, dom.BInteractons, dom.Particles.Size());
    // std::cout<<"Stress tensor for particle "<<0<<" : "<<stresses[0]<<"\n";
    // Array<double> strainEF(dom.Particles.Size());
    // for(size_t p=0;p<dom.Particles.Size();p++)
    //   strainEF[p] = StrainEnergyField(stresses[p], /*Poisson's ratio*/0.1);

    // dat.strainEF = strainEF;
    // std::cout<<"Strain Energy Field for particle "<<0<<" : "<<strainEF[0]<<"\n";

    // dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
    // _fs.Printf("%s%04d", "brick_planes",Restart);
    // dom.WriteXDMF(_fs.CStr());
    // dom.WritePOV(_fs.CStr());
    // std::cout<<"Saving initial domain structure..."<<std::endl;
    // _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
    // dom.Save(_fs.CStr());
    // std::cout<<"Saved initial domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    dom.Solve (Tf,dt,dtOut, NULL, &Report, filekey.CStr(), RenderVideo, Nproc);
    // dom.Solve (Tf,dt,dtOut, NULL, NULL, filekey.CStr(), RenderVideo, Nproc);
    std::cout<<"Saving initial domain structure..."<<std::endl;
    _fs.Printf("%s_final%04d", filekey.CStr(),Restart);
    dom.Save(_fs.CStr());
    std::cout<<"Saved final domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    _fs.Printf("final_points_%i.xyz", Restart);
    dom.SavePoints(_fs.CStr(), bTag);
}
MECHSYS_CATCH
