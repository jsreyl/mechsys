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
 * This program is distributed in the hope that it will be useful,      *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of       *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the         *
 * GNU General Public License for more details.                         *
 *                                                                      *
 * You should have received a copy of the GNU General Public License    *
 * along with this program. If not, see <http://www.gnu.org/licenses/>  *
 ************************************************************************/
// Elastic beam.


// MechSys
#include <mechsys/lbmmpm/Domain.h>
#include <math.h>
#include <iostream>
#include <fstream>

struct UserData
{
    Vec3_t fluidforce;
    double rhomax;
    double rhomin;
};

void Setup (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
    //#pragma omp parallel for schedule(static) num_threads(dom.Nproc)
	//for (size_t ix=0; ix<dom.LBMDOM.Ndim(0); ++ix)
	//for (size_t iy=0; iy<dom.LBMDOM.Ndim(1); ++iy)
	//for (size_t iz=0; iz<dom.LBMDOM.Ndim(2); ++iz)
    //{
        //dom.LBMDOM.BForce[0][ix][iy][iz] = dom.LBMDOM.Rho[0][ix][iy][iz]*dat.fluidforce;
    //}
    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
	for (size_t i=0; i<dom.LBMDOM.Ndim(1); ++i)
	for (size_t j=0; j<dom.LBMDOM.Ndim(2); ++j)
	{
        if (dom.LBMDOM.IsSolid[0][0][i][j]) continue;
        double * f = dom.LBMDOM.F[0][0][i][j];
        
        f[1] = 1.0/3.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-f[2]-2*f[3]-2*f[4]-2*f[5]-2*f[6]-4*f[8]+2*dat.rhomax);
        f[7] = 1.0/24.0*(-2*f[0]-4*f[10]-4*f[12]-4*f[14]-4*f[2] +f[3]-5*f[4]  +f[5]-5*f[6]+20*f[8]+2*dat.rhomax);
        f[9] = 1.0/24.0*(-2*f[0]+20*f[10]-4*f[12]-4*f[14]-4*f[2]+f[3]-5*f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);
        f[11]= 1.0/24.0*(-2*f[0]-4*f[10]+20*f[12]-4*f[14]-4*f[2]-5*f[3]+f[4]  +f[5]-5*f[6]-4*f[8]+2*dat.rhomax);
        f[13]= 1.0/24.0*(-2*f[0]-4*f[10]-4 *f[12]+20*f[14]-4*f[2]-5*f[3]+  f[4]-5*f[5]+f[6]-4*f[8]+2*dat.rhomax);

        dom.LBMDOM.Vel[0][0][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][0][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k];
            dom.LBMDOM.Vel[0][0][i][j] +=  dom.LBMDOM.F[0][0][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][0][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][0][i][j];
	}

    #pragma omp parallel for schedule(static) num_threads(dom.Nproc)
	for (size_t i=0; i<dom.LBMDOM.Ndim(1); ++i)
	for (size_t j=0; j<dom.LBMDOM.Ndim(2); ++j)
	{
        if (dom.LBMDOM.IsSolid[0][dom.LBMDOM.Ndim(0)-1][i][j]) continue;
        double * f = dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j];

        f[2] = 1/3.0* (-2*f[0]-f[1]-2*(2*f[11]+2*f[13]+f[3]+f[4]+f[5]+f[6]+2*f[7]+2*f[9]-dat.rhomin));
        f[8] = 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] - 5*f[5] + f[6] +20*f[7] - 4*f[9] + 2*dat.rhomin);
        f[10]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] - 4*f[13] - 5*f[3] + f[4] + f[5] - 5*f[6] - 4*f[7] + 20*f[9] + 2*dat.rhomin) ;
        f[12]= 1/24.0*(-2*f[0] - 4*f[1] + 20*f[11] - 4*f[13] + f[3] - 5*f[4] - 5*f[5] + f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        f[14]= 1/24.0*(-2*f[0] - 4*f[1] - 4*f[11] + 20*f[13] + f[3] - 5*f[4] + f[5] - 5*f[6] -  4*f[7] - 4*f[9] + 2*dat.rhomin);
        
        dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] = OrthoSys::O;
        dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j] = 0.0;
        for (size_t k=0;k<dom.LBMDOM.Nneigh;k++)
        {
            dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j] +=  dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j][k];
            dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] +=  dom.LBMDOM.F[0][dom.LBMDOM.Ndim(0)-1][i][j][k]*dom.LBMDOM.C[k];
        }
        dom.LBMDOM.Vel[0][dom.LBMDOM.Ndim(0)-1][i][j] *= dom.LBMDOM.Cs/dom.LBMDOM.Rho[0][dom.LBMDOM.Ndim(0)-1][i][j];
	}
}

void Report (LBMMPM::Domain & dom, void * UD)
{
    UserData & dat = (*static_cast<UserData *>(UD));
}

int main(int argc, char **argv) try
{
    //Number of cores
    size_t Nproc = 0.75*omp_get_max_threads();
    if (argc>1) Nproc = atoi(argv[1]);

    //Properties of LBM
    size_t nx = 200;
    size_t ny = 100;
    size_t nz = 100;
    double nu = 1.0e-6;
    double dx = 0.05;

    //Properties of MPM
    size_t ndiv = 5; //number of divisions per x lenght
    double K	= 10.0e4;
    double Nu	= 0.3;
    double	E	= (3.0*(1.0-2.0*Nu))*K;
    double	G	= E/(2.0*(1.0+Nu));
    double rho  = 3000.0;
    double Lx   = 1.0;   
    double Ly   = 1.0;
    double Lz   = 3.0;
    double Cs	= sqrt(E/rho);
    double h    = Ly/ndiv;
    double dt   = 0.1*h/Cs;
    double Dx   = 2.0*Ly/ndiv;
    double Bc   = Dx;
    //double Bc   = 1.1*Lz;

    LBMMPM::Domain dom(D3Q15,nu,iVec3_t(nx,ny,nz),dx,dt);
    UserData dat;
    dom.UserData = &dat;
    dat.rhomax = 1010.0;
    dat.rhomax =  990.0;

    dom.MPMDOM.AddRectangularBeam(-1, Vec3_t(0.5*nx*dx-0.5*Lx,0.5*ny*dx-0.5*Ly,0.0), Vec3_t(0.5*nx*dx+0.5*Lx,0.5*ny*dx+0.5*Ly,Lz), rho, ndiv);
    dom.MPMDOM.ResizeDomain(Vec3_t(0.5*nx*dx-2.0*Lx,0.5*ny*dx-2.0*Ly,0.5*nz*dx-2.0*Lz), Vec3_t(0.5*nx*dx+2.0*Lx,0.5*ny*dx+2.0*Ly,0.5*nz*dx+2.0*Lz),Dx);
    
    //Setting properties for the material points
    double v0 = 1.0;
    for (size_t ip=0; ip < dom.MPMDOM.Particles.Size(); ip++)
    {
        dom.MPMDOM.Particles[ip]->K = K;
        dom.MPMDOM.Particles[ip]->G = G;
        if (dom.MPMDOM.Particles[ip]->x(2)<Bc)
        {
            dom.MPMDOM.Particles[ip]->FixVeloc();
            dom.MPMDOM.Particles[ip]->Tag = -2;
        }
        //dom.MPMDOM.Particles[ip]->v(0) = v0*(dom.MPMDOM.Particles[ip]->x(2)-Dx)/(Lz-Dx);
    }
    dom.MPMDOM.Gn = 1.0e0;

    //Setting properties for nodes
    for (size_t in=0; in < dom.MPMDOM.Nnodes; in++)
    {
        Vec3_t xn;
        dom.MPMDOM.NodePosition(in,xn);
        if (xn(2)<Bc)
        {
            dom.MPMDOM.Nodes[in].FixVeloc();
        }
    }

    //Setting intial conditions of fluid
    for (size_t ix=0;ix<nx;ix++)
    for (size_t iy=0;iy<ny;iy++)
    for (size_t iz=0;iz<nz;iz++)
    {
        Vec3_t v(0.0,0.0,0.0);
        iVec3_t idx(ix,iy,iz);
        dom.LBMDOM.Initialize(0,idx,1000.0/*rho*/,v);
        if (iz==0)
        {
            dom.LBMDOM.IsSolid[0][ix][iy][iz] = true;
        }
    }  
    
    double Tf = 1.0; // Final Time
    dom.Solve(Tf,Tf/200,Setup,NULL,"tlbmmpm01",true,Nproc);
}
MECHSYS_CATCH

