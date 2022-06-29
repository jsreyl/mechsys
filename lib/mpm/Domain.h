/************************************************************************
 * MechSys - Open Library for Mechanical Systems                        *
 * Copyright (C) 2020 Sergio Galindo                                    *
 * Copyright (C) 2020 Pei Zhang                                         *
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


#ifndef MECHSYS_MPM_DOMAIN_H
#define MECHSYS_MPM_DOMAIN_H

// Hdf5
#ifdef USE_HDF5
#include <hdf5.h>
#include <hdf5_hl.h>
#endif

// Std lib
#ifdef USE_OMP
#include <omp.h>
#endif

//STD
#include<iostream>

// boost => to read json files
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/foreach.hpp>

// Mechsys
#include <mechsys/util/util.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/util/stopwatch.h>
#include <mechsys/util/numstreams.h>
#include <mechsys/mpm/Node.h>
#include <mechsys/mpm/shape_functions.h>
#include <mechsys/dem/domain.h>


enum MPMethod
{
    Solid,  ///< Simulaiton of flexible solids
    Soil,   ///< Simulation of soil
    Fluid,  ///< Simulation of fluid
};

namespace MPM
{
struct MtData
{
    size_t                        ProcRank; ///< Rank of the thread
    size_t                          N_Proc; ///< Total number of threads
    Array<size_t >                   Valid; ///< Array of valid nodes
};

inline size_t Pt2idx(iVec3_t & iv, iVec3_t & Dim) // Calculates the index of the cell at coordinates iv for a cubic lattice of dimensions Dim
{
    return iv(0) + iv(1)*Dim(0) + iv(2)*Dim(0)*Dim(1);
}

inline void   idx2Pt(size_t n, iVec3_t & iv, iVec3_t & Dim) // Calculates the coordinates from the index
{
    iv(0) = n%Dim(0);
    iv(1) = (n/Dim(0))%(Dim(1));
    iv(2) = n/(Dim(0)*Dim(1));
}

class Domain
{
public:
    //Constructor
    Domain() {Nproc = 1; Gn = 0.001;};
    //Domain(MPMethod Method, ///< Type of material to solve
    //iVec3_t Ndim            ///< Vector of integers with the dimensions
    //);

    //typedefs
    typedef void (*ptDFun_t) (Domain & Dom, void * UserData);

    //Data
    void *               UserData;              ///< User Data
    size_t                  Nproc;              ///< Number of procesors for OMP computation
    size_t                idx_out;              ///< The discrete time step for output
    String                FileKey;              ///< File Key for output files
    double                     Dx;              ///< Grid size in physical units
    double                     Gn;              ///< Dissipative constant 
    double                     Dt;              ///< Time step in physical units
    double                   Time;              ///< Total simulation time
    Vec3_t                   Xlow;              ///< Lower bound of cubic grid
    Vec3_t                  Xhigh;              ///< Upper bound of cubic grid
    Node *                  Nodes;              ///< Array of Nodes
    iVec3_t                  Ndim;              ///< Array with grid dimensions    
    size_t                 Nnodes;              ///< Total number of nodes
    Array<Particle *>   Particles;              ///< Array of particles
    Array<size_t >         VNodes;              ///< Array of valid nodes

    //Methods
    void AddFromJson (int Tag, char const * Filename, double rho, double scale, size_t nx); ///< Loads a Jason mesh into a set of MPM particles
    void AddRectangularBeam(int Tag, Vec3_t const & Xmin, Vec3_t const & Xmax, double rho, size_t nx);   ///< Creates a solid rectangular beam of particles between Xmin and Xmax with a number of divisions nx per x lenght and density rho
    void ResizeDomain(Vec3_t const & Xmin, Vec3_t const & Xmax, double deltax);    ///< Creates a mesh between Xmin and Xmax limits with deltax as the grid size
    void NodePosition(size_t in, Vec3_t & Xn);                                      ///< Gives a posiiton of a node in the coordinate system
    void BoundingBox(Vec3_t & minX, Vec3_t & maxX);                                 ///< Give dimensions of the bounding box enclosing the particles
    void ParticleToNode();                                                          ///< Imprinting particle information into the mesh
    void NodeToParticle();                                                          ///< Transfering back node information into particles
    void OneStepUSF();                                                             ///< One time step in the USF formulation
    void Solve(double Tf, double dt, double dtOut, ptDFun_t ptSetup=NULL, ptDFun_t ptReport=NULL,
    char const * FileKey=NULL, bool RenderVideo=true, size_t Nproc=1);            ///< Solve the Domain dynamics

#ifdef USE_HDF5    
    void WriteXDMF         (char const * FileKey);                                      ///< Save a xdmf file for visualization
#endif
};

inline void Domain::AddFromJson (int Tag, char const * Filename, double rho, double scale, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing Mesh --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    DEM::Domain demdom;
    demdom.AddFromJson(-1, Filename, 0.0, 1.0, 1.0);
    Vec3_t Xmin,Xmax;
    demdom.BoundingBox(Xmin,Xmax);
    double deltax = (Xmax(0)-Xmin(0))/nx;
    size_t Nx = nx;
    size_t Ny = (Xmax(1)-Xmin(1))/deltax;
    size_t Nz = (Xmax(2)-Xmin(2))/deltax;
    size_t ip = 0;
    for (size_t ix=0;ix<Nx;ix++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t iz=0;iz<Nz;iz++)
    {
        Vec3_t X0 = Xmin + deltax*Vec3_t(0.5+ix,0.5+iy,0.5+iz);
        if (demdom.Particles[0]->IsInside(X0)&&demdom.Particles[0]->IsInsideAlt(X0))
        {
            Particles.Push(new Particle(Tag, X0, OrthoSys::O, rho*deltax*deltax*deltax, deltax*deltax*deltax));
            ip++;
        }
    }
    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,ip,TERM_RST);
}

inline void Domain::AddRectangularBeam(int Tag, Vec3_t const & Xmin, Vec3_t const & Xmax, double rho, size_t nx)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Introducing rectangular Beam --------------------------------------------%s\n",TERM_CLR1,TERM_RST);
    double deltax = (Xmax(0)-Xmin(0))/nx;
    size_t Nx = nx;
    size_t Ny = (Xmax(1)-Xmin(1))/deltax;
    size_t Nz = (Xmax(2)-Xmin(2))/deltax;
    size_t Np = Nx*Ny*Nz;
    //double Vol = (Xmax(0)-Xmin(0))*(Xmax(1)-Xmin(1))*(Xmax(2)-Xmin(2));
    double Vol = Np*deltax*deltax*deltax;
    double Mp = rho*Vol/Np;
    size_t ip = 0;
    for (size_t ix=0;ix<Nx;ix++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t iz=0;iz<Nz;iz++)
    {
        Vec3_t X0 = Xmin + deltax*Vec3_t(0.5+ix,0.5+iy,0.5+iz);
        Particles.Push(new Particle(Tag, X0, OrthoSys::O, Mp, Vol/Np));
        ip++;
    }
    printf("%s  Num of material points   = %zd%s\n",TERM_CLR2,Np,TERM_RST);
}

inline void Domain::ResizeDomain(Vec3_t const & Xmin, Vec3_t const & Xmax, double deltax)
{
    Util::Stopwatch stopwatch;
    printf("\n%s--- Initializing MPM Domain --------------------------------------------%s\n",TERM_CLR1,TERM_RST);

    double Lx = (Xmax(0)-Xmin(0));
    double Ly = (Xmax(1)-Xmin(1));
    double Lz = (Xmax(2)-Xmin(2));
    if (Lx<deltax||Ly<deltax||Lz<deltax)
    {
        throw new Fatal("Recheck mesh dimensions and grid size since no cells can be produced");
    }
    size_t Nx = Lx/deltax;
    size_t Ny = Ly/deltax;
    size_t Nz = Lz/deltax;
    Ndim      = iVec3_t(Nx,Ny,Nz);
    Nnodes    = Nx*Ny*Nz;
    Dx        = deltax;
    Nodes     = new Node [Nnodes];
    for (size_t iz=0;iz<Nz;iz++)
    for (size_t iy=0;iy<Ny;iy++)
    for (size_t ix=0;ix<Nx;ix++)
    {
        iVec3_t Pt(ix,iy,iz);
        size_t idx = Pt2idx(Pt,Ndim);
        Nodes[idx].X   = Pt;
        Nodes[idx].Idx = idx;
    }

    //Assigning domain limits
    Xlow  = Xmin;
    Xhigh = Xmax;
   
    printf("%s  Num of nodes   = %zd%s\n",TERM_CLR2,Nx*Ny*Nz,TERM_RST);
}

inline void Domain::NodePosition(size_t in, Vec3_t & Xn)
{
    iVec3_t Pt;
    idx2Pt(in,Pt,Ndim);
    Xn = Vec3_t(Dx*Pt(0),Dx*Pt(1),Dx*Pt(2));
    Xn += Xlow;
}

inline void Domain::BoundingBox(Vec3_t & minX, Vec3_t & maxX)
{
    if (Particles.Size()==0) throw new Fatal("MPM::Domain::BoundingBox: There are no particles to build the bounding box");
    minX = Vec3_t(Particles[0]->x(0)-Particles[0]->Dx(0), Particles[0]->x(1)-Particles[0]->Dx(1), Particles[0]->x(2)-Particles[0]->Dx(2));
    maxX = Vec3_t(Particles[0]->x(0)+Particles[0]->Dx(0), Particles[0]->x(1)+Particles[0]->Dx(1), Particles[0]->x(2)+Particles[0]->Dx(2));
    for (size_t i=1; i<Particles.Size(); i++)
    {
        if (minX(0)>Particles[i]->x(0)-Particles[i]->Dx(0)) minX(0) = Particles[i]->x(0)-Particles[i]->Dx(0);
        if (maxX(0)<Particles[i]->x(0)+Particles[i]->Dx(0)) maxX(0) = Particles[i]->x(0)+Particles[i]->Dx(0);
        if (minX(1)>Particles[i]->x(1)-Particles[i]->Dx(1)) minX(1) = Particles[i]->x(1)-Particles[i]->Dx(1);
        if (maxX(1)<Particles[i]->x(1)+Particles[i]->Dx(1)) maxX(1) = Particles[i]->x(1)+Particles[i]->Dx(1);
        if (minX(2)>Particles[i]->x(2)-Particles[i]->Dx(2)) minX(2) = Particles[i]->x(2)-Particles[i]->Dx(2);
        if (maxX(2)<Particles[i]->x(2)+Particles[i]->Dx(2)) maxX(2) = Particles[i]->x(2)+Particles[i]->Dx(2);
    }
}

#ifdef USE_HDF5    
inline void Domain::WriteXDMF (char const * FileKey)
{
    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Storing center of mass data
    float * Posvec  = new float[3*Particles.Size()];
    float * Velvec  = new float[3*Particles.Size()];
    int   * Tag     = new int  [  Particles.Size()];
    float * Sigma	= new float[6*Particles.Size()];
    float * Strain	= new float[6*Particles.Size()];

    for (size_t i=0;i<Particles.Size();i++)
    {

        Posvec[3*i  ] = float(Particles[i]->x(0));
        Posvec[3*i+1] = float(Particles[i]->x(1));
        Posvec[3*i+2] = float(Particles[i]->x(2));
        Velvec[3*i  ] = float(Particles[i]->v(0));
        Velvec[3*i+1] = float(Particles[i]->v(1));
        Velvec[3*i+2] = float(Particles[i]->v(2));
        Sigma [6*i  ] = float(Particles[i]->S(0,0));
        Sigma [6*i+1] = float(Particles[i]->S(0,1));
        Sigma [6*i+2] = float(Particles[i]->S(0,2));
        Sigma [6*i+3] = float(Particles[i]->S(1,1));
        Sigma [6*i+4] = float(Particles[i]->S(1,2));
        Sigma [6*i+5] = float(Particles[i]->S(2,2));
        Strain[6*i  ] = float(Particles[i]->E(0,0));
        Strain[6*i+1] = float(Particles[i]->E(0,1));
        Strain[6*i+2] = float(Particles[i]->E(0,2));
        Strain[6*i+3] = float(Particles[i]->E(1,1));
        Strain[6*i+4] = float(Particles[i]->E(1,2));
        Strain[6*i+5] = float(Particles[i]->E(2,2));
        Tag   [i]     = int  (Particles[i]->Tag);  
    }

    hsize_t dims[1];
    dims[0] = 3*Particles.Size();
    String dsname;
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dims[0] = Particles.Size();
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );
    dims[0] = 6*Particles.Size();
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
    

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    oss << "   <Grid Name=\"MPM_points\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Sigma\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Sigma \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Strain\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"10\" Format=\"HDF\">\n";
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
#endif

inline void Domain::ParticleToNode()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].Reset();
    }
    VNodes.Resize(0);

    //Transfering info from particle to node
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        iVec3_t nlow,nhigh;
        Vec3_t Xp = Particles[ip]->x;
        Vec3_t Lp = Particles[ip]->Dx;
        Mat3_t Vsp = -Particles[ip]->V*Particles[ip]->S;
        Vec3_t Fex = Particles[ip]->b + Particles[ip]->h;

        //nlow (0) = (size_t(trunc((Xp(0)-Lp(0)-Xlow(0))/Dx)-1+Ndim(0)))%Ndim(0);
        //nlow (1) = (size_t(trunc((Xp(1)-Lp(1)-Xlow(1))/Dx)-1+Ndim(1)))%Ndim(1);
        //nlow (2) = (size_t(trunc((Xp(2)-Lp(2)-Xlow(2))/Dx)-1+Ndim(2)))%Ndim(2);
        //nhigh(0) = (size_t(ceil ((Xp(0)+Lp(0)-Xlow(0))/Dx)+1+Ndim(0)))%Ndim(0);
        //nhigh(1) = (size_t(ceil ((Xp(1)+Lp(1)-Xlow(1))/Dx)+1+Ndim(1)))%Ndim(1);
        //nhigh(2) = (size_t(ceil ((Xp(2)+Lp(2)-Xlow(2))/Dx)+1+Ndim(2)))%Ndim(2);

        nlow (0) = (size_t)std::max((trunc((Xp(0)-0.5*Lp(0)-Xlow(0))/Dx)-1),0.0);
        nlow (1) = (size_t)std::max((trunc((Xp(1)-0.5*Lp(1)-Xlow(1))/Dx)-1),0.0);
        nlow (2) = (size_t)std::max((trunc((Xp(2)-0.5*Lp(2)-Xlow(2))/Dx)-1),0.0);
        nhigh(0) = (size_t)std::min((ceil ((Xp(0)+0.5*Lp(0)-Xlow(0))/Dx)+1),(double)Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((Xp(1)+0.5*Lp(1)-Xlow(1))/Dx)+1),(double)Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((Xp(2)+0.5*Lp(2)-Xlow(2))/Dx)+1),(double)Ndim(2)-1.0);

        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            double Sf=0.0; //Shape function value
            Vec3_t Gs=OrthoSys::O; //Shape function gradient
            //Vec3_t Xn((nx+0.5)*Dx,(ny+0.5)*Dx,(nz+0.5)*Dx);
            Vec3_t Xn(nx*Dx,ny*Dx,nz*Dx);
            Xn += Xlow;
            iVec3_t Pt(nx,ny,nz);
            GIMP3D(Xp,Xn,Dx,Lp,Sf,Gs);
            if (Sf>0.0)
            {
                size_t nc = Pt2idx(Pt,Ndim);
                Particles[ip]->Nodes.Push(nc);
                Particles[ip]->Sfunc.Push(Sf);
                Particles[ip]->GSf  .Push(Gs);
                double ms = Sf*Particles[ip]->m;
                Vec3_t VspGs;
                Mult(Vsp,Gs,VspGs);
                Vec3_t dF = Sf*Fex + VspGs;
                omp_set_lock(&Nodes[nc].lck);
                Nodes[nc].Mass += ms;
                Nodes[nc].Mn   += ms*Particles[ip]->v + Dt*dF;
                Nodes[nc].Fn   += dF;
                //Nodes[nc].Sn    = Nodes[nc].Sn + ms*Particles[ip]->S;
                Nodes[nc].valid = true;
                omp_unset_lock(&Nodes[nc].lck);
            }
        }
        if (Particles[ip]->Nodes.Size()<2)
        {
            std::cout << Time << " " << ip << " " << Particles[ip]->x << " " << Particles[ip]->Nodes.Size() << " " << Particles[ip]->V0 << " " << Particles[ip]->Dx << std::endl;
            std::cout << nlow << " " << nhigh << std::endl;
            throw new Fatal("MPM:Domain:ParticleToNode No mesh cells intersected by particles");
        }
    }

    //Producing an array of only valid nodes
    for (size_t in=0; in<Nnodes; in++)
    {
        if (Nodes[in].valid) VNodes.Push(Nodes[in].Idx);
    }
}

inline void Domain::NodeToParticle()
{
    //Solving the equations at the nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(Gn,Dt);
        Nodes[in].FixNode();
    }

    //TRansfering info from node to particle
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        set_to_zero(Particles[ip]->L);
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            double sf = Particles[ip]->Sfunc[nn];
            Vec3_t gs = Particles[ip]->GSf[nn];
            Vec3_t an = Nodes[nc].Fn/Nodes[nc].Mass;
            Particles[ip]->v += sf*an*Dt;
            Particles[ip]->x += sf*Nodes[nc].Vn*Dt;
            Mat3_t Gv;
            //Dyad(Nodes[nc].Vn,gs,Gv);
            Dyad(gs,Nodes[nc].Vn,Gv);
            Particles[ip]->L = Particles[ip]->L + Gv;
        }
        Particles[ip]->CalcVol(Dt);
        Particles[ip]->CalcStress(Dt);
        Particles[ip]->Reset(Dt);
    }

}

inline void Domain::OneStepUSF()
{
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].Reset();
    }
    VNodes.Resize(0);
    
    //Transfering info from particle to node
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        iVec3_t nlow,nhigh;
        Vec3_t Xp = Particles[ip]->x;
        Vec3_t Lp = Particles[ip]->Dx;

        //nlow (0) = (size_t(trunc((Xp(0)-Lp(0)-Xlow(0))/Dx)-1+Ndim(0)))%Ndim(0);
        //nlow (1) = (size_t(trunc((Xp(1)-Lp(1)-Xlow(1))/Dx)-1+Ndim(1)))%Ndim(1);
        //nlow (2) = (size_t(trunc((Xp(2)-Lp(2)-Xlow(2))/Dx)-1+Ndim(2)))%Ndim(2);
        //nhigh(0) = (size_t(ceil ((Xp(0)+Lp(0)-Xlow(0))/Dx)+1+Ndim(0)))%Ndim(0);
        //nhigh(1) = (size_t(ceil ((Xp(1)+Lp(1)-Xlow(1))/Dx)+1+Ndim(1)))%Ndim(1);
        //nhigh(2) = (size_t(ceil ((Xp(2)+Lp(2)-Xlow(2))/Dx)+1+Ndim(2)))%Ndim(2);

        nlow (0) = (size_t)std::max((trunc((Xp(0)-0.5*Lp(0)-Xlow(0))/Dx)-1),0.0);
        nlow (1) = (size_t)std::max((trunc((Xp(1)-0.5*Lp(1)-Xlow(1))/Dx)-1),0.0);
        nlow (2) = (size_t)std::max((trunc((Xp(2)-0.5*Lp(2)-Xlow(2))/Dx)-1),0.0);
        nhigh(0) = (size_t)std::min((ceil ((Xp(0)+0.5*Lp(0)-Xlow(0))/Dx)+1),(double)Ndim(0)-1.0);
        nhigh(1) = (size_t)std::min((ceil ((Xp(1)+0.5*Lp(1)-Xlow(1))/Dx)+1),(double)Ndim(1)-1.0);
        nhigh(2) = (size_t)std::min((ceil ((Xp(2)+0.5*Lp(2)-Xlow(2))/Dx)+1),(double)Ndim(2)-1.0);

        for (size_t nx=nlow(0); nx<=nhigh(0); nx++)
        for (size_t ny=nlow(1); ny<=nhigh(1); ny++)
        for (size_t nz=nlow(2); nz<=nhigh(2); nz++)
        {
            double Sf=0.0; //Shape function value
            Vec3_t Gs=OrthoSys::O; //Shape function gradient
            //Vec3_t Xn((nx+0.5)*Dx,(ny+0.5)*Dx,(nz+0.5)*Dx);
            Vec3_t Xn(nx*Dx,ny*Dx,nz*Dx);
            Xn += Xlow;
            iVec3_t Pt(nx,ny,nz);
            GIMP3D(Xp,Xn,Dx,Lp,Sf,Gs);
            if (Sf>0.0)
            {
                size_t nc = Pt2idx(Pt,Ndim);
                Particles[ip]->Nodes.Push(nc);
                Particles[ip]->Sfunc.Push(Sf);
                Particles[ip]->GSf  .Push(Gs);
                double ms = Sf*Particles[ip]->m;
                omp_set_lock(&Nodes[nc].lck);
                Nodes[nc].Mass += ms;
                Nodes[nc].Mn   += ms*Particles[ip]->v;
                //Nodes[nc].Sn    = Nodes[nc].Sn + ms*Particles[ip]->S;
                Nodes[nc].valid = true;
                omp_unset_lock(&Nodes[nc].lck);
            }
        }
        if (Particles[ip]->Nodes.Size()<2)
        {
            std::cout << Time << " " << ip << " " << Particles[ip]->x << " " << Particles[ip]->Nodes.Size() << " " << Particles[ip]->V0 << " " << Particles[ip]->Dx << std::endl;
            std::cout << nlow << " " << nhigh << std::endl;
            throw new Fatal("MPM:Domain::OneStepUSF No mesh cells intersected by particles");
        }
    }

    //Producing an array of only valid nodes
    for (size_t in=0; in<Nnodes; in++)
    {
        if (Nodes[in].valid) VNodes.Push(Nodes[in].Idx);
    }

    //Fixing Dirichlet nodes
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(0.0,0.0);
        Nodes[in].FixNode();
    }

    //Transfering info from node to particle
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        set_to_zero(Particles[ip]->L);
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            Vec3_t gs = Particles[ip]->GSf[nn];
            Mat3_t Gv;
            Dyad(gs,Nodes[nc].Vn,Gv);
            Particles[ip]->L = Particles[ip]->L + Gv;
        }
        
        Particles[ip]->CalcVol(Dt);
        Particles[ip]->CalcStress(Dt);
        
        Mat3_t Vsp = -Particles[ip]->V*Particles[ip]->S;
        Vec3_t Fex = Particles[ip]->b + Particles[ip]->h;
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            double Sf = Particles[ip]->Sfunc[nn];
            Vec3_t Gs = Particles[ip]->GSf[nn];
            Vec3_t VspGs;
            Mult(Vsp,Gs,VspGs);
            Vec3_t dF = Sf*Fex + VspGs;
            omp_set_lock(&Nodes[nc].lck);            
            Nodes[nc].Fn   += dF;
            Nodes[nc].Mn   += Dt*dF;
            omp_unset_lock(&Nodes[nc].lck);
        }
    }

    //Solve equations at the nodes    
    #pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t nn=0; nn<VNodes.Size(); nn++)
    {
        size_t in = VNodes[nn];
        Nodes[in].UpdateNode(Gn,Dt);
        Nodes[in].FixNode();
    }

    //Updating the positions of the MPM particles
	#pragma omp parallel for schedule(static) num_threads(Nproc)
    for (size_t ip=0; ip < Particles.Size(); ip++)
    {
        for (size_t nn=0; nn<Particles[ip]->Nodes.Size(); nn++)
        {
            size_t nc = Particles[ip]->Nodes[nn];
            double sf = Particles[ip]->Sfunc[nn];
            Vec3_t gs = Particles[ip]->GSf[nn];
            Vec3_t an = Nodes[nc].Fn/Nodes[nc].Mass;
            Particles[ip]->v += sf*an*Dt;
            Particles[ip]->x += sf*Nodes[nc].Vn*Dt;
        }
        Particles[ip]->Reset(Dt);
    }
}

inline void Domain::Solve(double Tf, double dt, double dtOut, ptDFun_t ptSetup, ptDFun_t ptReport,
                          char const * TheFileKey, bool RenderVideo, size_t TheNproc)
{

    FileKey.Printf("%s",TheFileKey);
    Nproc = TheNproc;

    Util::Stopwatch stopwatch;
    printf("\n%s--- Solving ---------------------------------------------------------------------%s\n",TERM_CLR1    , TERM_RST);
    printf("%s  Time step                        =  %g%s\n"       ,TERM_CLR2, dt                                    , TERM_RST);
    
    idx_out     = 0;
    Time        = 0.0;
    double tout = Time;
    Dt          = dt;
    //Initializing particles
    for (size_t ip=0;ip<Particles.Size();ip++)
    {
        Particles[ip]->xb = Particles[ip]->x - Particles[ip]->v*Dt;
    }

    while (Time < Tf)
    {
        if (ptSetup!=NULL) (*ptSetup) ((*this), UserData);
        if (Time >= tout)
        {
            if (TheFileKey!=NULL)
            {
                String fn;
                fn.Printf    ("%s_%04d", TheFileKey, idx_out);
                if ( RenderVideo) 
                {
                    #ifdef USE_HDF5
                    WriteXDMF(fn.CStr());
                    #endif
                }
                if (ptReport!=NULL) (*ptReport) ((*this), UserData);
            }
            tout += dtOut;
            idx_out++;
        }
        //ParticleToNode();
        //NodeToParticle();
        OneStepUSF();
        

        Time += dt;
        //std::cout << Time << std::endl;
    }
}
}
#endif //MECHSYS_MPM_DOMAIN_H
