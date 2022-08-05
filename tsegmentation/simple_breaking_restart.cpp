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
// #include <mechsys/dem/domain.h>
// #include <mechsys/dem/distance.h>
// #include <mechsys/util/fatal.h>
// #include <mechsys/mesh/structured.h>
// #include <mechsys/mesh/unstructured.h>

//User
#include "plane_segmentationv7.h"
#include "strain_energy_field.h"
#include "user_report.h"


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
  // force at particle tagged -3
  //Ff is a force fixed on the particles
  UserData & dat = (*static_cast<UserData *>(UD));
  if (dat.test=="vibration")
    for(size_t i=0; i<dat.p.Size();i++)
      dat.p[i]->Ff=0.0,0.0,dat.Am*sin(dat.ome*Dom.Time);
  if (dat.test=="bending" || dat.test=="center")
    {
      // for(size_t i=0; i<dat.p.Size();i++)
      //   dat.p[i]->Ff=0.0,0.0,dat.Am;
      if (Dom.Time < 0.5*dat.Tf){
        for(size_t i=0; i<dat.p.Size();i++)
          dat.p[i]->Ff=0.0,0.0,dat.Am*2*Dom.Time/dat.Tf;
      }else{
        for(size_t i=0; i<dat.p.Size();i++)
          dat.p[i]->Ff=0.0,0.0,dat.Am;
      }
    }
  //Average strain
  double avgFz=0.;
  for(size_t i=0; i<dat.p.Size();i++){
      avgFz+=dat.p[i]->F(2);//Add the foces in z for the moving plane
    }
    avgFz/=dat.p.Size();
    dat.sz = avgFz/dat.A;
    // std::cout<<"stress sz : "<<dat.sz <<std::endl;
  }

  // Report is called every dtOut steps before particle initialization (i.e. force assignments) so it uses the values of the timestep immediately before to calculate and report results.
  void Report (DEM::Domain & Dom, void * UD)
  {
    // std::cout<<" Report!  T="<<Dom.Time<<"\n";
      UserData & dat = (*static_cast<UserData *>(UD));
      if (dat.test=="tensile" || dat.test=="flexion")
      {
        //Calculate the strain energy field for all the particles in the domain
        Array<Mat3_t> stresses = StressTensor(Dom.Interactons, Dom.Particles.Size());
        std::cout<<"Stresses calculated for all particles! \n";
        Array<double> strainEF(Dom.Particles.Size());
        double max_strainEF = 0.;
        for(size_t p=0;p<Dom.Particles.Size();p++){
          // std::cout<<"Particle "<<p<<" with stress tensor: \n "<<stresses[p]<<"\n";
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
            DEM::Domain _vdom;
            _vdom.domType = 1; //Make a subdomain so the particles added are not deleted on domain destruction
            _vdom.domID = 9999999; //Visalization domain doesn't make part of the recurrent simulation domains
            size_t _i = 0;
            _vdom.Particles.Push(Dom.Particles[p]); // This particle belongs to the domain Dom so don't add it to eigenparticles lest it be destroyed
            String _f_vis;
            Array<double> stress_force_max(0);
            stress_force_max.Push(stresses[p](0,0));
            _f_vis.Printf    ("%s_%04d", "stress_max_error", _i);
            std::cout<<"Writing visualization files for error _vdom.\n";
            _vdom.WriteXDMF_User(stress_force_max, _f_vis.CStr());
            for(size_t j=0; j<Dom.Interactons.Size();j++){
              DEM::Interacton * I = Dom.Interactons[j];
              if(I->P1->Index == p){
                std::cout<<"Interacton between p"<<I->P1->Index<<" and p"<<I->P2->Index<<" with forces F1="<<I->F1<<" and F2="<<I->F2<<"\n";
                std::cout<<"At X1="<<I->P1->x<<" and X2="<<I->P2->x<<"\n";
                _vdom.Particles.Push(Dom.Particles[I->P2->Index]);
                stress_force_max.Push(norm(I->F2));
                _i++;
                _f_vis.Printf    ("%s_%04d", "stress_max_error", _i);
                std::cout<<"Writing visualization files for error _vdom.\n";
                _vdom.WriteXDMF_User(stress_force_max, _f_vis.CStr());
              }else if(I->P2->Index == p){
                std::cout<<"Interacton between p"<<I->P1->Index<<" and p"<<I->P2->Index<<" with forces F1="<<I->F1<<" and F2="<<I->F2<<"\n";
                std::cout<<"At X1="<<I->P1->x<<" and X2="<<I->P2->x<<"\n";
                _vdom.Particles.Push(Dom.Particles[I->P1->Index]);
                stress_force_max.Push(norm(I->F1));
                _i++;
                _f_vis.Printf    ("%s_%04d", "stress_max_error", _i);
                std::cout<<"Writing visualization files for error _vdom.\n";
                _vdom.WriteXDMF_User(stress_force_max, _f_vis.CStr());
              }
            }
          }
          strainEF[p] = StrainEnergyField(stresses[p], /*Poisson's ratio*/0.1);
          if(strainEF[p] > max_strainEF) max_strainEF = strainEF[p];
          // std::cout<<"And strain energy field "<<strainEF[p]<<"\n";
        }
        std::cout<<"Strain energy field calculated! \n";
        std::cout<<"maximum strainEF: "<<max_strainEF<<"\n";
        // // dat.strainEF = strainEF;
        // std::cout<<"Strain energy field assigned to dat! \n";
        String fn;
        fn.Printf    ("%s_%04d", "StrainEnergyField", Dom.idx_out);
        Dom.WriteXDMF_User(strainEF,fn.CStr());
        // WriteXDMF_SEF(Dom, strainEF, fn.CStr());
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
        double strainEFThr = max_strainEF/10.; //3.0; // Strain Energy Field breaking threshold
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
            //Build plane using the principal stress component and the particle geometric center
            // Calculate the first component of the stress tensor
            Vec3_t eigvalR = Vec3_t(0.,0.,0.), _ = Vec3_t(0.,0.,0.), stress_v0 = Vec3_t(0.,0.,0.);
            Mat3_t _stress = stresses[pID];//Create a new matrix with the stress tensor since EigNonsymm destroys the matrices it uses
            EigNonsymm(_stress, eigvalR, _, stress_v0, _, _, _, _, _);
            Vec3_t planeNormal = Vec3_t(0.,0.,0.), planeCentroid = Vec3_t(0.,0.,0.);
            planeNormal=stress_v0/norm(stress_v0);
            for(size_t v=0; v<Dom.Particles[pID]->Verts.Size(); v++) planeCentroid+=*(Dom.Particles[pID]->Verts[v]);
            planeCentroid/=Dom.Particles[pID]->Verts.Size();
            std::cout<<"Cutting plane with normal" << planeNormal<<" at "<<planeCentroid<<"\n";
            // ofbreak << Util::_10_6 << Dom.Time << Util::_8s << max_strainEF << Util::_8s << particlesToBreakID << std::endl;
            // Print particle number, strainEF, ID in the domain, cutting plane normal, cutting plane centre
            ofbreak << Util::_2 << pB << Util::_8s << strainEF[pID] << Util::_2 << pID << Util::_8s << planeNormal(0) << Util::_8s << planeNormal(1) << Util::_8s << planeNormal(2) << Util::_8s << planeCentroid(0) << Util::_8s << planeCentroid(1) << Util::_8s << planeCentroid(2) << std::endl;
        }
        ofbreak.close();
      }
    }
    if (dat.test=="bending")
    {
        if (Dom.idx_out==0)
        {
            double tol = dat.L0(2)/20.0;
            for (size_t i=0;i<Dom.Particles.Size();i++)
            {
                for (size_t j=0;j<Dom.Particles[i]->Verts.Size();j++)
                {
                    if (fabs((*Dom.Particles[i]->Verts[j])(2))<tol)
                    {
                        dat.vm.Push (Dom.Particles[i]->Verts[j]);
                        dat.vm0.Push(Vec3_t(*Dom.Particles[i]->Verts[j]));
                    }
                }
            }
        }
        String fs;
        fs.Printf("%s_%08d.res",Dom.FileKey.CStr(),Dom.idx_out);
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_10_6 << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s <<"ux" << Util::_8s << "uy" << Util::_8s << "uz" << std::endl;
        for (size_t i=0;i<dat.vm.Size();i++)
        {
            dat.oss_ss << Util::_8s <<            dat.vm0[i](0) << Util::_8s <<            dat.vm0[i](1) << Util::_8s <<            dat.vm0[i](2);
            dat.oss_ss << Util::_8s << (*dat.vm[i])(0)-dat.vm0[i](0) << Util::_8s << (*dat.vm[i])(1)-dat.vm0[i](1) << Util::_8s << (*dat.vm[i])(2)-dat.vm0[i](2) << std::endl;
        }
        dat.oss_ss.close();
    }
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
    double max_time;    // Time at which the maximum strainEF was found in previous run, should match name of the restart file to use, if any
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
    }

    // user data and domain
    UserData dat;
    DEM::Domain dom(&dat);
    dom.Alpha = verlet;
    dom.Dilate= true; // Dilate for visualizations

    // Tags of bulk, static, moving particles
    int bTag = -1, sTag =-2, mTag =-3;

    if(Restart){
      std::cout<<"####### RESTART SCHEME #######"<<std::endl;
      String _fs;
      _fs.Printf("run_%04d/%s_initial%04d", Restart-1, filekey.CStr(), Restart-1);
      std::cout<<"Loading initial domain structure from "<<_fs.CStr()<<"..."<<std::endl;
      dom.Load(_fs.CStr());
      std::cout<<"Loaded initial domain structure!"<<std::endl;

      //Cut particles according to max strain EF
      // Previous run generated a file containing particle ids and the planes necessary to cut them so use that here.
      Array<DEM::Particle*> segParticles(0);
      Array<size_t> particlesToBreakID(0);
      size_t _ = 0; //assign data we won't use
      size_t pID = 0; //Particle to cut
      double _strainEF = 0.;
      Vec3_t planeNormal = Vec3_t(0.,0.,0.), planeCentroid = Vec3_t(0.,0.,0.);
      // double _nx, _ny, _nz, _x, _y, _z; // storing coordinates and vector components for plane normal and centroid
      _fs.Printf("run_%04d/particles_break_index_%04f.res", Restart-1, max_time);
      // std::ifstream infile("run_0000/particles_break_index_5.600050.res");
      std::ifstream infile(_fs.CStr());
      while (infile >> _ >> _strainEF >> pID >> planeNormal(0) >> planeNormal(1) >> planeNormal(2) >> planeCentroid(0) >> planeCentroid(1) >> planeCentroid(2)){
        particlesToBreakID.Push(pID);
        std::cout<<"Forming a plane with normal" << planeNormal<<" at "<<planeCentroid<<"\n";
        DEM::Domain vdom;
        vdom.domID = 999999; //Visualization
        vdom.domType = 1; //Subdomain
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
        std::cout<<"DOM"<<dom.domID<<": before bisecting segParticles Size:"<<prev_segSize<<"\n";
        std::cout<<"DOM"<<dom.domID<<": before bisecting intersectionIxs Size:"<<intersectionIxs.Size()<<"\n";
        BisectPolyhedron(dom.Particles[pID], planeNormal, planeCentroid, segParticles, intersectionIxs, /*mechsysErode*/4, 1e-3, true);
        size_t nParts = segParticles.Size() - prev_segSize;
        std::cout<<"Particle "<<pID<<" bisected into "<<nParts<<" parts. Creating new simulation domain...\n";
        // Make a visualization of just the particle and the plane cutting it
        vdom.Particles.Clear(); //same as Resize(0);
        for(size_t i=prev_segSize; i<segParticles.Size();i++){
          vdom.Particles.Push(segParticles[i]);
          vdom.Dilate = false;
          String _f_vis, _f_vis_dil;
          _f_vis.Printf    ("%s_%04d", "particle_cut_", i-prev_segSize); _f_vis_dil.Printf    ("%s_%04d", "particle_cut_dilate_", i-prev_segSize);
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
      }
      std::cout<<"Starting actual new simulation domain. \n";
      DEM::Domain sdom(&dat);
      sdom.domID = dom.domID + 1;
      sdom.domType = 1; //Subdomain
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
        sdom.eigenParticles.Push(segParticles[sp]); // Make these be deleted on domain destruction
        // DON'T do this, this will create double free or corruption and subsequent segmentation fault
        // sdom.Particles[sdom.Particles.Size()-1]->Initialize(sdom.Particles.Size()-1);//Initialize this particle
        // sdom.Particles[sdom.Particles.Size()-1]->InitializeVelocity(dat.dt);//Initialize this particle
      }
      std::cout<<"DOM"<<sdom.domID<<": Number of particles in previous domain: "<<dom.Particles.Size()<<"\nNumber of particles in new domain: "<<sdom.Particles.Size()<<"\n";
      std::cout<<"DOM"<<sdom.domID<<": Number of eigenparticles: "<<sdom.eigenParticles.Size()<<"\n";
      std::cout<<"DOM"<<sdom.domID<<": Particles: \n"<<sdom.Particles<<"\nEigenParticles: \n"<<sdom.eigenParticles<<"\n";

      // Initialize the UserData structure
      dat.test = test;
      dat.A  = Lx*Ly;
      dat.Am   = Am;
      dat.ome  = ome;
      dat.Tf   = Tf;
      Vec3_t Xmin,Xmax;
      dom.BoundingBox(Xmin,Xmax);
      dat.L0   = Xmax - Xmin;
      // NOTE: ADD Cohesive interactions here!!! XXX
      // std::cout<<"Adding cohesive interactions between particles of the new domain\n";
      // LocalImposeParticleCohesion(dat.bTag, sdom);//Only add it for particles in the bulk (those tagged bTag), not to moving rods or planes
      // Set some particles to not move, or move at a given speed
      //identify the moving lid
      sdom.GetParticles(mTag,dat.p);//This is the top plane (and upper cylinder if the test is flexion)
      Vec3_t CompV = Vec3_t(0.0, 0.0, ex*Ly/Tf); // NOTE: There is a difference in how the compression velocity is calculated here and in simple_breaking3
      if (dat.test=="tensile")
        {
          dat.p[0]->FixVeloc();
          dat.p[0]->v = CompV;
          std::cout<<"-> New Sim: Velocity for top plane: "<<dat.p[0]->v<<std::endl;
        }
      if (dat.test=="flexion"){
        for(size_t i=0; i<dat.p.Size();i++){//Set the top plane and upper cylinder moving down
          dat.p[i]->FixVeloc();
          dat.p[i]->v = CompV;
        }
        std::cout<<"-> New Sim: Velocity for top plane: "<<dat.p[0]->v<<std::endl;
      }
      // fix -2 particles, this is the bottom plane (and the two lower cylinders if the test is flexion)
      Array<DEM::Particle *> p;
      sdom.GetParticles(sTag,p);
      for(size_t i=0;i<p.Size();i++) p[i]->FixVeloc();
      // Set material properties
      Dict B;
      B.Set(bTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);//Breaking material
      B.Set(sTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);//Static parts
      B.Set(mTag,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);//Moving planes
      sdom.SetProps(B);
      //Run new simulation from here, this is recursive simulation
      String f_vis;
      f_vis.Printf    ("%s_%04d", "brick_planes", Restart);
      std::cout<<"Writing visualization files for sdom.\n";
      sdom.WriteXDMF(f_vis.CStr());
      sdom.WritePOV(f_vis.CStr());
      std::cout<<"Saving initial domain structure..."<<std::endl;
      _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
      sdom.Save(_fs.CStr());
      std::cout<<"Saved initial domain structure into "<<_fs.CStr()<<"!"<<std::endl;
      _fs.Printf("%s_restart%04d", filekey.CStr(),Restart);
      sdom.Solve (Tf,dt,dtOut, &Setup, &Report, _fs.CStr(), RenderVideo, Nproc);
      std::cout<<"Saving final domain structure..."<<std::endl;
      _fs.Printf("%s_final%04d", filekey.CStr(),Restart);
      sdom.Save(_fs.CStr());
      std::cout<<"Saved final domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    }
    else{//Create particles as per usual
      if (ptype=="voronoi") dom.AddVoroPack (/*Tag*/-1, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,Lz,  /*Number of cells*/nx,ny,nz,  /*Density*/rho,  /*Cohesion*/true,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
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

      if(test=="flexion"){//Add cilynders at the bottom and top of the block
        std::cout<<"Adding cylinder at "<<Vec3_t(0.,0.,Lz/2.+Rc+R)<<" with radious "<<Rc<<"\n";
        std::cout<<"Lz: "<<Lz<<" Rc: "<<Rc<<" R: "<<R<<"\n";
        //Add a cylinder as the connection of two circles of radius R0 and R1 located at X0 and X1
        dom.AddCylinder(/*Tag*/-3,/*X0*/Vec3_t(0.,-Lz/2.,Lz/2.+R+Rc),/*R0*/Rc,/*X1*/Vec3_t(0.,Lz/2.,Lz/2.+R+Rc),/*R1*/Rc,/*R*/R,/*rho*/rho);
        dom.AddCylinder(/*Tag*/-2,/*X0*/Vec3_t(-Lx/2.+dx,-Lz/2.,-Lz/2.-R-Rc),/*R0*/Rc,/*X1*/Vec3_t(-Lx/2.+dx,Lz/2.,-Lz/2.-R-Rc),/*R1*/Rc,/*R*/R,/*rho*/rho);
        //Make the third cylinder randomly smaller to give the fracture a predilect direction
        srand(0);
        double Rb=Rc-(rand()%10+1)*Rc/100.;
        std::cout <<"Rb: "<<Rb<<"\n";
        dom.AddCylinder(/*Tag*/-2,/*X0*/Vec3_t(Lx/2.-dx,-Lz/2.,-Lz/2.-R-Rc),/*R0*/Rb,/*X1*/Vec3_t(Lx/2.-dx,Lz/2.,-Lz/2.-R-Rc),/*R1*/Rb,/*R*/R,/*rho*/rho);
      }
      // Initialize the UserData structure
      dat.test = test;
      dat.A  = Lx*Ly;
      dat.Am   = Am;
      dat.ome  = ome;
      dat.Tf   = Tf;
      Vec3_t Xmin,Xmax;
      dom.BoundingBox(Xmin,Xmax);
      dat.L0   = Xmax - Xmin;

      //Generate planes at the bottom and top
      // dom.GenBoundingPlane(-2,R,/*Rescaling factor*/1.0,/*Cohesion*/true);
      dom.AddPlane(/*Tag*/-2,/*Position*/Vec3_t(0.0,0.0,Xmin(2)-R),/*R*/R,/*Lx*/1.5*dat.L0(0),/*Ly*/1.5*dat.L0(1),/*rho*/1.0,/*Angle*/M_PI,/*Axis*/&OrthoSys::e0);//Bottom plane
      dom.AddPlane(/*Tag*/-3,/*Position*/Vec3_t(0.0,0.0,Xmax(2)+R),/*R*/R,/*Lx*/1.5*dat.L0(0),/*Ly*/1.5*dat.L0(1),/*rho*/1.0,/*Angle*/M_PI,/*Axis*/&OrthoSys::e0);//Top plane
    // input
    double cam_x=0.0, cam_y=2*Lx, cam_z=0.0;
    dom.CamPos = cam_x, cam_y, cam_z;

    //identify the moving lid
    dom.GetParticles(-3,dat.p);//This is the top plane (and upper cylinder if the test is flexion)
    if (test=="tensile")
    {
        dat.p[0]->FixVeloc();
        dat.p[0]->v = 0.0, ex*dat.L0(1)/Tf, 0.0;
    }
    if (test=="flexion"){
      for(size_t i=0; i<dat.p.Size();i++){//Set the top plane and upper cylinder moving down
        dat.p[i]->FixVeloc();
        // dat.p[i]->v = 0.0, 0.0, ex*dat.L0(2)/Tf;
        dat.p[i]->v = 0.0, 0.0, ex*Ly/Tf;
      }
      std::cout<<"Length of the system: "<<dat.L0(0)<<" "<<dat.L0(1)<<" "<<dat.L0(2)<<std::endl;
      std::cout<<"Length of the brick: "<<dat.L0(0)<<" "<<dat.L0(1)<<" "<<dat.L0(2)-4.*Rc<<std::endl;
      std::cout<<"Given velocity for top plane: "<<ex*Ly/Tf<<std::endl;
      std::cout<<"Velocity for top plane: "<<dat.p[0]->v<<std::endl;
    }

    //set the element properties
    Dict B;
    B.Set(-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(-2,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    B.Set(-3,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
    dom.SetProps(B);

    // fix -2 particles, this is the bottom plane (and the two lower cylinders if the test is flexion)
    Array<DEM::Particle *> p;
    dom.GetParticles(-2,p);
    for(size_t i=0;i<p.Size();i++) p[i]->FixVeloc();

    Array<Mat3_t> stresses = StressTensor(dom.Interactons, dom.Particles.Size());
    std::cout<<"Stress tensor for particle "<<0<<" : "<<stresses[0]<<"\n";
    Array<double> strainEF(dom.Particles.Size());
    for(size_t p=0;p<dom.Particles.Size();p++)
      strainEF[p] = StrainEnergyField(stresses[p], /*Poisson's ratio*/0.1);

    dat.strainEF = strainEF;
    std::cout<<"Strain Energy Field for particle "<<0<<" : "<<strainEF[0]<<"\n";
    String _fs;
    _fs.Printf("%s%04d", "brick_planes",Restart);
    dom.WriteXDMF(_fs.CStr());
    dom.WritePOV(_fs.CStr());
    std::cout<<"Saving initial domain structure..."<<std::endl;
    _fs.Printf("%s_initial%04d", filekey.CStr(),Restart);
    dom.Save(_fs.CStr());
    std::cout<<"Saved initial domain structure into "<<_fs.CStr()<<"!"<<std::endl;
    //dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
    }
}
MECHSYS_CATCH
