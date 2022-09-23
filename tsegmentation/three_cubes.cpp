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
    double tprint;
    int Render;
};


void Setup(DEM::Domain &Dom){
    std::cout<<"TIME:"<<Dom.Time<<std::endl;
    std::cout <<"Interaction log: Interactons->"<<Dom.Interactons.Size()<<" CInteractons->"<<Dom.CInteractons.Size()<<" BInteractons->"<<Dom.BInteractons.Size()<<std::endl;
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
    std::cout<<"REPORT! T="<<Dom.Time<<std::endl;
    UserData & dat = (*static_cast<UserData *>(UD));
    if (Dom.idx_out == 0.){
        //String fs;
        //fs.Printf("%s_walls.res",Dom.FileKey.CStr());
        dat.ofcubes.open("center_parinfo.dat");
    }
    if (!Dom.Finished){
        //Print data
        Vec3_t xp, vp, Fp, BFp, CFp; Mat3_t Sp; int Nc, pID;
        DEM::Particle * p = Dom.Particles[0];
        xp = p->x; vp = p->v; Fp = p->F; pID = p->Index; Nc = CalculateContacts(Dom.Interactons, pID);//sp = strainEF[pID]; 
        Mat3_t s0 = p->S;
        //std::cout <<"Interaction log: Interactons->"<<Dom.Interactons.Size()<<" CInteractons->"<<Dom.CInteractons.Size()<<" BInteractons->"<<Dom.BInteractons.Size()<<std::endl;
        std::cout<<"Stress tensor for particle 0 s0:"<<s0<<"\n";
        BFp = CalculateForce(Dom.BInteractons, pID);
        CFp = CalculateForce(Dom.CInteractons, pID);
        //Array<Mat3_t> stress = StressTensor(Dom.Interactons, Dom.Particles.Size());
        //Mat3_t s0 = stress[0];
        Array<double> strainEF(Dom.Particles.Size());
        //#pragma omp parallel for
        for(size_t i=0; i<Dom.Particles.Size();i++) strainEF[i] = StrainEnergyField(Dom.Particles[i]->S, 0.2);
        std::cout<<"Stress tensor for particle 0 after strainEF calc:"<<Dom.Particles[pID]->S<<"\n";
        std::cout<<"Stress tensor for particle 0 after strainEF calc:"<<p->S<<"\n";
        std::cout<<"Stress tensor for particle 0 after strainEF calc s0:"<<s0<<"\n";
        if(Dom.Time >= dat.tprint) {
            //std::cout<<"Printing XDMF for strain energy field? "<<dat.Render<<" at time"<<Dom.Time<<"\n";
            if(dat.Render){
                String fn;
                fn.Printf    ("%s_%04d", "three_cubes", Dom.idx_out);
                Dom.WriteXDMF_User(strainEF,fn.CStr());
            }
            String _fs;
            _fs.Printf("%s_break_index_%04f.res", "particles", Dom.Time);
            std::ofstream ofbreak;
            ofbreak.open(_fs.CStr());
            Array<Vec3_t> stress_vs(3);
            Vec3_t eigvalR = Vec3_t(0.,0.,0.), _ = Vec3_t(0.,0.,0.);
            Mat3_t _stress = Dom.Particles[pID]->S;//Create a new matrix with the stress tensor since EigNonsymm destroys the matrices it uses
            EigNonsymm(_stress, eigvalR, _, stress_vs[0], stress_vs[1], stress_vs[2], _, _, _);
            //stress_vs[0] = stress_v0; stress_vs[1] = stress_v1; stress_vs[2] = stress_v2;
            std::cout<<"CHECK: Stress tensor with eigenvalues "<<eigvalR<<" and eigenvectors "<<stress_vs<<std::endl;
            // EigNonsymm returns unordered eigenvalues and eigenvectors, we want sigma_1>sigma_2>sigma_3
            size_t largest_id = 0; double largest_stress = eigvalR(0);
            for(size_t i=1; i<3;i++) {
                //if (fabs(eigvalR(i)) > largest_stress){ //absolute value considers both traction (positive) and compression (negative) stresses
                if (eigvalR(i) > largest_stress){ //consider the highest traction stress
                    largest_stress = eigvalR(i); largest_id = i;
                }
            }
            std::cout<<"CHECK: Largest stress eigenvalue "<<eigvalR(largest_id)<<" and eigenvector "<<stress_vs[largest_id]<<std::endl;
            if(eigvalR(largest_id)<0.) std::cout<<"WARNING: All eigenvalues are negative for particle "<<pID<<std::endl;
            Vec3_t planeNormal = Vec3_t(0.,0.,0.), planeCentroid = Vec3_t(0.,0.,0.);
            planeNormal=stress_vs[largest_id]/norm(stress_vs[largest_id]);
            for(size_t v=0; v<Dom.Particles[pID]->Verts.Size(); v++) planeCentroid+=*(Dom.Particles[pID]->Verts[v]);
            planeCentroid/=Dom.Particles[pID]->Verts.Size();
            std::cout<<"Cutting plane with normal" << planeNormal<<" at "<<planeCentroid<<"\n";
            // Print particle number, strainEF, ID in the domain, cutting plane normal, cutting plane centre
            ofbreak << Util::_2 << 0 << Util::_8s << strainEF[pID] << Util::_2 << pID << Util::_8s << planeNormal(0) << Util::_8s << planeNormal(1) << Util::_8s << planeNormal(2) << Util::_8s << planeCentroid(0) << Util::_8s << planeCentroid(1) << Util::_8s << planeCentroid(2) << std::endl;
            ofbreak.close();

            dat.tprint += 0.1;
        }
        //std::cout <<"Interaction log: Interactons->"<<Dom.Interactons.Size()<<" CInteractons->"<<Dom.CInteractons.Size()<<" BInteractons->"<<Dom.BInteractons.Size()<<std::endl;
        dat.ofcubes << Util::_2 << pID << Util::_10_6 << Dom.Time << Util::_8s << xp(0) << Util::_8s << xp(1) << Util::_8s << xp(2) << \
        Util::_8s << vp(0) << Util::_8s << vp(1) << Util::_8s << vp(2) << \
        Util::_8s << Fp(0) << Util::_8s << Fp(1) << Util::_8s << Fp(2) << \
        Util::_8s << BFp(0) << Util::_8s << BFp(1) << Util::_8s << BFp(2) << \
        Util::_8s << CFp(0) << Util::_8s << CFp(1) << Util::_8s << CFp(2) << \
        Util::_2 << Nc << Util::_2 << Dom.Interactons.Size() <<Util::_2 << Dom.BInteractons.Size() <<Util::_2 << Dom.CInteractons.Size() <<\
        Util::_8s << s0(0,0) << Util::_8s << s0(0,1) << Util::_8s << s0(0,2) <<
        Util::_8s << s0(1,0) << Util::_8s << s0(1,1) << Util::_8s << s0(1,2) <<
        Util::_8s << s0(2,0) << Util::_8s << s0(2,1) << Util::_8s << s0(2,2) <<
        Util::_8s << strainEF[pID] << std::endl; //<< Util::_8s << sp 
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
    dat.tprint = 0.;
    dat.Render = Render;
    //DEM::Domain dom; //(&dat);
    //String _fs;
    dom.Alpha = verlet;
    dom.Dilate = true;
    //---Particle addition---
    //Cube centers
    Vec3_t X0 (0., 0., Lx/2.);
    Vec3_t X1 (0., Lx+2*R-dx, Lx/2.);
    Vec3_t X2 (0., -Lx-2*R+dx, Lx/2.);
    // Axis for rotations
    Vec3_t ax (0., 0., 1.);
    dom.AddCube(pTag, X0, R, Lx, rho, /*Angle*/0., &ax);
    dom.AddCube(pTag, X1, R, Lx, rho, /*Angle*/0., &ax);
    dom.AddCube(pTag, X2, R, Lx, rho, /*Angle*/0., &ax);
    for(size_t i=0; i<dom.Particles.Size();i++){
        dom.Particles[i]->v = 0.,0.,0.;
        dom.Particles[i]->w = 0.,0.,0.;
        dom.Particles[i]->Eroded=true;
    }
    if (fabs(vy)>0.){ //Positive for particles separating
        dom.Particles[1]->v = 0., vy, 0.;
        dom.Particles[2]->v = 0., -vy, 0.;
    }
    dom.Particles[1]->Ff = 0., -Fy, 0.; //Compression
    dom.Particles[2]->Ff = 0., Fy, 0.; //Compression
    //Add beam force between particles
    LocalImposeParticleCohesion(pTag, dom, 1e-8, 1e-3, dx2);
    std::cout <<"Interaction log: Interactons->"<<dom.Interactons.Size()<<" CInteractons->"<<dom.CInteractons.Size()<<" BInteractons->"<<dom.BInteractons.Size()<<std::endl;

    std::cout<<"Saving initial domain structure..."<<std::endl;
    dom.Save("three_cubes_initial0");
    std::cout<<"Saved initial domain structure into "<<"three_cubes_initial0"<<"!"<<std::endl;
    dom.WriteXDMF("three_cubes_initial");
    //dom.Solve(/*Tf*/Tf, /*dt*/dt, /*dtOut*/dt*10., /*Setup*/ &Setup, /*Report*/ &Report2, "three_cubes" , /*Render_video*/0); //2 for XDMF
    dom.Solve(/*Tf*/Tf, /*dt*/dt, /*dtOut*/dt*10., /*Setup*/ NULL, /*Report*/ &Report2, "three_cubes" , /*Render_video*/0); //2 for XDMF
    
}
MECHSYS_CATCH