#include "plane_segmentation.h"
#include "strain_energy_fieldv3.h"

/*
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/mesh/unstructured.h>
*/
//User data struct
struct UserData
{
    /* data */
    std::ofstream ofcubes;
};


void Setup(DEM::Domain &Dom, void * UD){
    UserData & dat = (*static_cast<UserData *>(UD));
    if (Dom.idx_out == 0.){
        //String fs;
        //fs.Printf("%s_walls.res",Dom.FileKey.CStr());
        dat.ofcubes.open("parinfo.dat");
    }
    if (!Dom.Finished){
        //Print data
        //dat.ofcubes<<
    }
    else  dat.ofcubes.close();
}

void Report(DEM::Domain & Dom){
    std::cout<<"REPORT! T="<<Dom.Time<<std::endl;
    std::cout<<"Dilating polyhedrons? "<<Dom.Dilate<<std::endl;
    std::cout <<"Interaction log: Interactons->"<<Dom.Interactons.Size()<<" CInteractons->"<<Dom.CInteractons.Size()<<" BInteractons->"<<Dom.BInteractons.Size()<<std::endl;
    Vec3_t xp, vp, Fp, BFp, CFp; int Nc, pID;
    for(size_t i=0; i<Dom.Particles.Size();i++){
        DEM::Particle * p = Dom.Particles[i];
        xp = p->x; vp = p->v; Fp = p->F; pID = p->Index; Nc = CalculateContacts(Dom.Interactons, pID);//sp = strainEF[pID]; 
        BFp = CalculateForce(Dom.BInteractons, pID);
        CFp = CalculateForce(Dom.CInteractons, pID);
        std::cout <<"Particle "<<i<<":"<< Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << Util::_8s << BFp(0) << Util::_8s << BFp(1) << Util::_8s << BFp(2) << Util::_8s << CFp(0) << Util::_8s << CFp(1) << Util::_8s << CFp(2) << Util::_2 << Nc << std::endl; //<< Util::_8s << sp 
    }
}

void Report2(DEM::Domain & Dom, void * UD){
    UserData & dat = (*static_cast<UserData *>(UD));
    if (Dom.idx_out == 0.){
        //String fs;
        //fs.Printf("%s_walls.res",Dom.FileKey.CStr());
        dat.ofcubes.open("parinfo.dat");
    }
    if (!Dom.Finished){
        //Print data
        Vec3_t xp, vp, Fp, BFp, CFp; int Nc, pID;
        DEM::Particle * p = Dom.Particles[0];
        xp = p->x; vp = p->v; Fp = p->F; pID = p->Index; Nc = CalculateContacts(Dom.Interactons, pID);//sp = strainEF[pID]; 
        BFp = CalculateForce(Dom.BInteractons, pID);
        CFp = CalculateForce(Dom.CInteractons, pID);
        //std::cout <<"Interaction log: Interactons->"<<Dom.Interactons.Size()<<" CInteractons->"<<Dom.CInteractons.Size()<<" BInteractons->"<<Dom.BInteractons.Size()<<std::endl;
        dat.ofcubes << Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << \
        Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << \
        Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << \
        Util::_8s << BFp(0) << Util::_8s << BFp(1) << Util::_8s << BFp(2) << \
        Util::_8s << CFp(0) << Util::_8s << CFp(1) << Util::_8s << CFp(2) << \
        Util::_2 << Nc << Util::_2 << Dom.Interactons <<Util::_2 << Dom.BInteractons <<Util::_2 << Dom.CInteractons <<std::endl; //<< Util::_8s << sp 
    }
    else  dat.ofcubes.close();
}


int main(int argc, char **argv) try
{
    double dx = 0., dx2=0., vy=0., Fy=0., Tf = 5.0, dt=1e-3;
    int Render = 0;
    if (argc>1) dx = atof(argv[1]);
    if (argc>2) dx2 = atof(argv[2]);
    if (argc>3) vy = atof(argv[3]);
    if (argc>4) Tf = atof(argv[4]);
    if (argc>5) dt = atof(argv[5]);
    if (argc>6) Fy = atof(argv[6]);
    if (argc>7) Render = atoi(argv[7]);
    std::cout<< "Initial cube intersection:"<<dx<<std::endl;
	//--- Particle properties ---
    int pTag=-1;
    double R = 0.1, Lx = 1., rho=3.;
    //bool cohesion=true;
    //--- Domain definition ---
    double verlet = 0.05;
    UserData dat;
    DEM::Domain dom(&dat);
    //DEM::Domain dom; //(&dat);
    //String _fs;
    dom.Alpha = verlet;
    dom.Dilate = true;
    //---Particle addition---
    //Cube centers
    Vec3_t X0 (0., -Lx/2.-R+dx, Lx/2.);
    Vec3_t X1 (0., Lx/2.+R-dx, Lx/2.);
    // Axis for rotations
    Vec3_t ax (0., 0., 1.);
    dom.AddCube(pTag, X0, R, Lx, rho, /*Angle*/0., &ax);
    dom.AddCube(pTag, X1, R, Lx, rho, /*Angle*/0., &ax);
    for(size_t i=0; i<dom.Particles.Size();i++){
        dom.Particles[i]->v = 0.,0.,0.;
        dom.Particles[i]->w = 0.,0.,0.;
        dom.Particles[i]->Eroded=true;
    }
    if (fabs(vy)>0.){ //Positive for particles separating
        dom.Particles[0]->v = 0., -vy, 0.;
        dom.Particles[1]->v = 0., vy, 0.;
    }
    dom.Particles[0]->Ff = 0., Fy, 0.; //Compression
    dom.Particles[1]->Ff = 0., -Fy, 0.; //Compression
    //Add beam force between particles
    LocalImposeParticleCohesion(pTag, dom, 1e-8, 1e-3, dx2);
    std::cout <<"Interaction log: Interactons->"<<dom.Interactons.Size()<<" CInteractons->"<<dom.CInteractons.Size()<<" BInteractons->"<<dom.BInteractons.Size()<<std::endl;
    
    dom.WriteXDMF("test_cubes_initial");
    dom.Dilate= true;
    dom.Solve(/*Tf*/Tf, /*dt*/dt, /*dtOut*/dt*10., /*Setup*/ NULL, /*Report*/ &Report2, "test_cubes" , /*Render_video*/Render); //2 for XDMF
    
}
MECHSYS_CATCH
