/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang
 * Copyright (C) 2020 Siqi Sun                                         *
 *                                                                      *
 * This program is free software: you can redistribute it and/or modify *
 * it under the terms of the GNU General Public License as published by *
 * the Free Software Foundation, either version 3 of the License, or    *
 * any later version.                                                   *
 *                                                                      *
 *
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Elastic beam static bending


// MechSys
#include <mechsys/mpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    Array<MPM::Particle *> ForcePar; //Forced particles
    Array<MPM::Particle *> PlaneBeamPar0; //Particles at each plane of the beam
    Array<double >         x0;
    Array<double >         x2; 
    double                       Fb; //Maximun applied body force
    double                       Tf; //Final simulation time
    size_t                     ndiv;
    double                       Lx;
    double                       E;
    double                       Lz;
    std::ofstream                oss_ss; // The output file
};

void Setup (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));

    double F = std::min(10.0*dat.Fb*dom.Time/dat.Tf,dat.Fb)/dat.ForcePar.Size();

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
    for (size_t ip=0; ip<dat.ForcePar.Size(); ip++)
    {
        dat.ForcePar[ip]->b = Vec3_t(0.0,0.0,-F);
    }
       
}

void Report (MPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD)); 

    double Dx   = 1.5*dat.Lx/dat.ndiv;
    

    if (dom.idx_out==0)
    {
        String fs;
        fs.Printf("tmpm01_Comparison.res");
        dat.oss_ss.open(fs.CStr());
        dat.oss_ss << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "y(mpm)" << "\n";
    }


    if (dom.idx_out==99)
    {
        
       
        for (size_t i=0;i<dat.x0.Size();i++)
        {
            double x = dat.x0[i]-Dx;
            double z = dat.x2[i]-dat.PlaneBeamPar0[i]->x(2);
            double yt =  2.0*dat.Fb*pow(x,2)*(3.0*(dat.Lx-Dx)-x)/(dat.E*pow(dat.Lz,4));
            dat.oss_ss << Util::_8s << x << Util::_8s << yt <<Util::_8s << z<< "\n";

        }
                
    }
}

int main(int argc, char **argv) try
{
    MPM::Domain dom;
    UserData dat;
    dom.UserData = &dat;
   
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);
    size_t ndiv = 50; //number of divisions per x lenght
     
    double K	= 10.0e4;
    double Nu	= 0.3;
    double	E	= (3.0*(1.0-2.0*Nu))*K;
    double	G	= E/(2.0*(1.0+Nu));
    double rho  = 3000.0;
    double Cs	= sqrt(E/rho);
    double h    = 1.0/ndiv;
    double dt   = 1.0*h/Cs;
    double Lx   = 1.0;   
    double Ly   = 0.16;
    double Lz   = 0.16;
    double Dx   = 2.0*Lx/ndiv;
    double Bc   = Dx;

    dat.ndiv    = ndiv;
    dat.Lx      = Lx;
    dat.Lz      = Lz;
    dat.E       = E;
    


    dom.AddRectangularBeam(-1, Vec3_t(0.0,0.0,0.0), Vec3_t(Lx,Ly,Lz), rho, ndiv);
    dom.ResizeDomain(Vec3_t(-0.5*Lx,-0.5*Lx,-0.5*Lx),Vec3_t(1.5*Lx,1.5*Lx,1.5*Lx), Dx);

    //Setting properties for the material points
    for (size_t ip=0; ip < dom.Particles.Size(); ip++)
    {
        dom.Particles[ip]->K = K;
        dom.Particles[ip]->G = G;
        if (dom.Particles[ip]->x(0)<Bc)
        {
            dom.Particles[ip]->FixVeloc();
            dom.Particles[ip]->Tag = -2;
        }
        if (dom.Particles[ip]->x(0)>Lx*(ndiv-1.0)/ndiv)
        {
            dom.Particles[ip]->Tag = -3;
            dat.ForcePar.Push(dom.Particles[ip]);
        }

        if (dom.Particles[ip]->x(2)<0.6*Lz && dom.Particles[ip]->x(2)>0.4*Lz)
        {
            dat.PlaneBeamPar0.Push(dom.Particles[ip]);
            dat.x0.Push(dom.Particles[ip]->x(0));
            dat.x2.Push(dom.Particles[ip]->x(2));
        }
    }

    //Setting properties for nodes
    for (size_t in=0; in < dom.Nnodes; in++)
    {
        Vec3_t xn;
        dom.NodePosition(in,xn);
        if (xn(0)<Bc)
        {
            dom.Nodes[in].FixVeloc();
        }
    }

    //dom.ParticleToNode();
    //dom.NodeToParticle();
    dom.Gn = 1.0e0;
    double Tf = 20.0; // Final Time
    dat.Tf = Tf;     
    dat.Fb = 1.0e0;  // Force applied at the end

    dom.Solve(Tf,dt,Tf/100,Setup,Report,"tmpm01",true,Nproc);

}
MECHSYS_CATCH

