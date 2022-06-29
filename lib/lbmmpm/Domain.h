/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
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


#ifndef MECHSYS_LBMMPM_DOMAIN_H
#define MECHSYS_LBMMPM_DOMAIN_H

// Mechsys
#include <mechsys/flbm/Domain.h>
#include <mechsys/mpm/Domain.h>

namespace LBMMPM
{

class Domain
{
public:

    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    // Constructor
    Domain (LBMethod      Method, ///< Type of array, for example D2Q9
    double                    nu, ///< Viscosity of the fluid
    iVec3_t               Ndim,   ///< Cell divisions per side
    double                dx,     ///< Space spacing
    double                dt);    ///< Time step

    //Methods
    void Reset();                    ///< Reset LBM grid
    void ImprintLattice();              ///< Imprint the MPM material point into the LBM grid
    void Solve(double Tf, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics
    
#ifdef USE_HDF5    
    void WriteXDMF         (char const * FileKey);                                      ///< Save a xdmf file for visualization
#endif

    //Variables
    double           dt;                ///< Time Step
    size_t        Nproc;                ///< Number of processors
    size_t      idx_out;                ///< The discrete time step for output
    String      FileKey;                ///< File Key for output files
    void *     UserData;                ///< User Data
    double         Time;                ///< Simulation time variable
    FLBM::Domain LBMDOM;                ///< The LBM domain
    MPM::Domain  MPMDOM;                ///< The MPM domain
};

inline Domain::Domain(LBMethod TheMethod, double Thenu, iVec3_t TheNdim, double Thedx, double Thedt)
{
    Nproc  = 1;
    dt = Thedt;
    Time = 0.0;
    LBMDOM = FLBM::Domain(TheMethod, Thenu, TheNdim, Thedx, Thedt);
    MPMDOM = MPM ::Domain();

    LBMDOM.Omeis = new double *** [TheNdim(0)];
    LBMDOM.Gamma = new double **  [TheNdim(0)];
    for (size_t ix=0; ix< TheNdim(0); ix++)
    {
        LBMDOM.Omeis[ix] = new double ** [TheNdim(1)];
        LBMDOM.Gamma[ix] = new double *  [TheNdim(1)];
        for (size_t iy=0; iy< TheNdim(1); iy++)
        {
            LBMDOM.Omeis[ix][iy] = new double * [TheNdim(2)];
            LBMDOM.Gamma[ix][iy] = new double   [TheNdim(2)];
            for (size_t iz=0; iz< TheNdim(2); iz++)
            {
                LBMDOM.Omeis[ix][iy][iz] = new double [LBMDOM.Nneigh];
            }
        }
    }
}

void Domain::ImprintLattice()
{
    double Cs = LBMDOM.Cs;
    double ld = LBMDOM.dx;
    double Tau = LBMDOM.Tau[0];

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < MPMDOM.Particles.Size(); ip++)
    {   
        iVec3_t nlow,nhigh;
        Vec3_t Xp = MPMDOM.Particles[ip]->x;
        Vec3_t Lp = MPMDOM.Particles[ip]->Dx;
        Vec3_t Vp = MPMDOM.Particles[ip]->v;

        nlow (0) = (size_t)std::max((trunc((Xp(0)-Lp(0)-ld)/ld)-1),0.0);
        nlow (1) = (size_t)std::max((trunc((Xp(1)-Lp(1)-ld)/ld)-1),0.0);
        nlow (2) = (size_t)std::max((trunc((Xp(2)-Lp(2)-ld)/ld)-1),0.0);
        nhigh(0) = (size_t)std::min((ceil ((Xp(0)+Lp(0)+ld)/ld)+1),(double)LBMDOM.Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((Xp(1)+Lp(1)+ld)/ld)+1),(double)LBMDOM.Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((Xp(2)+Lp(2)+ld)/ld)+1),(double)LBMDOM.Ndim(2)-1.0);
        
        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            Vec3_t Xn(nx*LBMDOM.dx,ny*LBMDOM.dx,nz*LBMDOM.dx);
            double dx = std::max((std::min(Xn(0)+ld,Xp(0)+0.5*Lp(0)+1.0*ld) - std::max(Xn(0),Xp(0)-0.5*Lp(0)-1.0*ld))/ld,0.0);
            double dy = std::max((std::min(Xn(1)+ld,Xp(1)+0.5*Lp(1)+1.0*ld) - std::max(Xn(1),Xp(1)-0.5*Lp(1)-1.0*ld))/ld,0.0);
            double dz = std::max((std::min(Xn(2)+ld,Xp(2)+0.5*Lp(2)+1.0*ld) - std::max(Xn(2),Xp(2)-0.5*Lp(2)-1.0*ld))/ld,0.0);
            double gamma = dx*dy*dz;
            double Bn  = (gamma*(Tau-0.5))/((1.0-gamma)+(Tau-0.5));
            double rho = LBMDOM.Rho[0][nx][ny][nz];
            //LBMDOM.Gamma[nx][ny][nz] = std::max(gamma,LBMDOM.Gamma[nx][ny][nz]);
            if (gamma >  LBMDOM.Gamma[nx][ny][nz])
            {
                LBMDOM.Gamma[nx][ny][nz] = gamma;
                for (size_t k = 0; k < LBMDOM.Nneigh; k++)
                {
                    double Fvpp     = LBMDOM.Feq(LBMDOM.Op[k],rho,Vp);
                    double Fvp      = LBMDOM.Feq(k           ,rho,Vp);
                    double Omega = LBMDOM.F[0][nx][ny][nz][LBMDOM.Op[k]] - Fvpp - (LBMDOM.F[0][nx][ny][nz][k] - Fvp);
                    MPMDOM.Particles[ip]->h += -Bn*Omega*Cs*Cs*ld*ld*LBMDOM.C[k];
                    LBMDOM.Omeis[nx][ny][nz][k] = Omega;
                }
            }
        }
    }
}

void Domain::Reset()
{
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ix=0; ix < LBMDOM.Ndim(0); ix++)
    for (size_t iy=0; iy < LBMDOM.Ndim(1); iy++)
    for (size_t iz=0; iz < LBMDOM.Ndim(2); iz++)
    {
        LBMDOM.Gamma[ix][iy][iz] = (double) LBMDOM.IsSolid[0][ix][iy][iz];
        for (size_t k=0; k<LBMDOM.Nneigh; k++)
        {
            LBMDOM.Omeis[ix][iy][iz][k] = 0.0;
        }
    }

    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip<MPMDOM.Particles.Size(); ip++)
    {
        MPMDOM.Particles[ip]->h = OrthoSys::O;
    }
}

void Domain::WriteXDMF(char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    //Save LBM domain
    size_t  Nx = LBMDOM.Ndim[0]/LBMDOM.Step;
    size_t  Ny = LBMDOM.Ndim[1]/LBMDOM.Step;
    size_t  Nz = LBMDOM.Ndim[2]/LBMDOM.Step;

    for (size_t j=0;j<LBMDOM.Nl;j++)
    {
        // Creating data sets
        float * Density   = new float[  Nx*Ny*Nz];
        float * Gamma     = new float[  Nx*Ny*Nz];
        float * Vvec      = new float[3*Nx*Ny*Nz];

        size_t i=0;
        for (size_t m=0;m<LBMDOM.Ndim(2);m+=LBMDOM.Step)
        for (size_t l=0;l<LBMDOM.Ndim(1);l+=LBMDOM.Step)
        for (size_t n=0;n<LBMDOM.Ndim(0);n+=LBMDOM.Step)
        {
            double rho    = 0.0;
            double gamma  = 0.0;
            Vec3_t vel    = OrthoSys::O;

            for (size_t ni=0;ni<LBMDOM.Step;ni++)
            for (size_t li=0;li<LBMDOM.Step;li++)
            for (size_t mi=0;mi<LBMDOM.Step;mi++)
            {
                rho    += LBMDOM.Rho    [j][n+ni][l+li][m+mi];
                //vel    += LBMDOM.Vel    [j][n+ni][l+li][m+mi];
                vel    += LBMDOM.Vel    [j][n+ni][l+li][m+mi]*(1.0-LBMDOM.Gamma[n+ni][l+li][m+mi]);
                gamma  += std::max(LBMDOM.Gamma[n+ni][l+li][m+mi], (double) LBMDOM.IsSolid[j][n+ni][l+li][m+mi]);
            }
            rho  /= LBMDOM.Step*LBMDOM.Step*LBMDOM.Step;
            gamma/= LBMDOM.Step*LBMDOM.Step*LBMDOM.Step;
            vel  /= LBMDOM.Step*LBMDOM.Step*LBMDOM.Step;
            Density [i]  = (float) rho;
            Gamma   [i]  = (float) gamma;
            Vvec[3*i  ]  = (float) vel(0);
            Vvec[3*i+1]  = (float) vel(1);
            Vvec[3*i+2]  = (float) vel(2);
            i++;
        }
        
        //Writing data to h5 file
        hsize_t dims[1];
        dims[0] = Nx*Ny*Nz;
        String dsname;
        dsname.Printf("Density_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Density );
        if (j==0)
        {
            dsname.Printf("Gamma");
            H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Gamma   );
        }
        dims[0] = 3*Nx*Ny*Nz;
        dsname.Printf("Velocity_%d",j);
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vvec    );
        dims[0] = 1;
        int N[1];
        N[0] = Nx;
        dsname.Printf("Nx");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Ny;
        dsname.Printf("Ny");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);
        dims[0] = 1;
        N[0] = Nz;
        dsname.Printf("Nz");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,N);

        delete [] Density ;
        delete [] Gamma   ;
        delete [] Vvec    ;
    }

    //Storing MPM data

    // Storing center of mass data
    float * Posvec  = new float[3*MPMDOM.Particles.Size()];
    float * Velvec  = new float[3*MPMDOM.Particles.Size()];
    int   * Tag     = new int  [  MPMDOM.Particles.Size()];
    float * Sigma	= new float[6*MPMDOM.Particles.Size()];
    float * Strain	= new float[6*MPMDOM.Particles.Size()];

    for (size_t i=0;i<MPMDOM.Particles.Size();i++)
    {

        Posvec[3*i  ] = float(MPMDOM.Particles[i]->x(0));
        Posvec[3*i+1] = float(MPMDOM.Particles[i]->x(1));
        Posvec[3*i+2] = float(MPMDOM.Particles[i]->x(2));
        Velvec[3*i  ] = float(MPMDOM.Particles[i]->v(0));
        Velvec[3*i+1] = float(MPMDOM.Particles[i]->v(1));
        Velvec[3*i+2] = float(MPMDOM.Particles[i]->v(2));
        Sigma [6*i  ] = float(MPMDOM.Particles[i]->S(0,0));
        Sigma [6*i+1] = float(MPMDOM.Particles[i]->S(0,1));
        Sigma [6*i+2] = float(MPMDOM.Particles[i]->S(0,2));
        Sigma [6*i+3] = float(MPMDOM.Particles[i]->S(1,1));
        Sigma [6*i+4] = float(MPMDOM.Particles[i]->S(1,2));
        Sigma [6*i+5] = float(MPMDOM.Particles[i]->S(2,2));
        Strain[6*i  ] = float(MPMDOM.Particles[i]->E(0,0));
        Strain[6*i+1] = float(MPMDOM.Particles[i]->E(0,1));
        Strain[6*i+2] = float(MPMDOM.Particles[i]->E(0,2));
        Strain[6*i+3] = float(MPMDOM.Particles[i]->E(1,1));
        Strain[6*i+4] = float(MPMDOM.Particles[i]->E(1,2));
        Strain[6*i+5] = float(MPMDOM.Particles[i]->E(2,2));
        Tag   [i]     = int  (MPMDOM.Particles[i]->Tag);  
    }

    hsize_t dims[1];
    dims[0] = 3*MPMDOM.Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dims[0] = MPMDOM.Particles.Size();
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );
    dims[0] = 6*MPMDOM.Particles.Size();
    dsname.Printf("Sigma");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Sigma);
    dsname.Printf("Strain");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Strain);


    delete [] Posvec;
    delete [] Velvec;
    delete [] Tag;
    delete [] Sigma;
    delete [] Strain;

    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);

	// Writing xmf fil
    std::ostringstream oss;

    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"LBM_Mesh\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"3DCoRectMesh\" Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\"/>\n";
    oss << "     <Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> 0.0 0.0 0.0\n";
    oss << "       </DataItem>\n";
    oss << "       <DataItem Format=\"XML\" NumberType=\"Float\" Dimensions=\"3\"> " << LBMDOM.Step*LBMDOM.dx << " " << LBMDOM.Step*LBMDOM.dx  << " " << LBMDOM.Step*LBMDOM.dx  << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    for (size_t j=0;j<LBMDOM.Nl;j++)
    {
    oss << "     <Attribute Name=\"Density_" << j << "\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Density_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity_" << j << "\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity_" << j << "\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    }
    oss << "     <Attribute Name=\"Gamma\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Nz << " " << Ny << " " << Nx << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Gamma\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << "   <Grid Name=\"MPM_points\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << MPMDOM.Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << MPMDOM.Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << MPMDOM.Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << MPMDOM.Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << MPMDOM.Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << MPMDOM.Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Strain \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    oss << " </Domain>\n";
    oss << "</Xdmf>\n";
    
    fn = FileKey;
    fn.append(".xmf");
    std::ofstream of(fn.CStr(), std::ios::out);
    of << oss.str();
    of.close();

}

inline void Domain::Solve(double Tf, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{
    idx_out     = 0;
    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;
    LBMDOM.Nproc = TheNproc;
    MPMDOM.Nproc = TheNproc;
    MPMDOM.Dt    = dt;
    
    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    printf("%s  C parameter                      =  %g%s\n"       ,TERM_CLR2, LBMDOM.Cs                             , TERM_RST);
    for (size_t i=0;i<LBMDOM.Nl;i++)
    {
    printf("%s  Tau of Lattice %zd                 =  %g%s\n"       ,TERM_CLR2, i, LBMDOM.Tau[i]                    , TERM_RST);
    }


    //Initializing MPM particles
    for (size_t ip=0;ip<MPMDOM.Particles.Size();ip++)
    {
        MPMDOM.Particles[ip]->xb = MPMDOM.Particles[ip]->x - MPMDOM.Particles[ip]->v*dt;
    }


    Reset();
    ImprintLattice();
    double tout = Time;

    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            //std::cout << "3" << std::endl;
            //std::cout << "4" << std::endl;
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    #else
                    //WriteVTK (fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }
        //Coupling between LBM and MPM
        Reset();
        ImprintLattice();

        //The LBM dynamics
        LBMDOM.CollideMPM();
        LBMDOM.StreamSC();

        //The MPM dynamics
        MPMDOM.OneStepUSF();

        //Reset domain
        

        Time += dt;
        //std::cout << Time << std::endl;
    }
}
}
#endif
