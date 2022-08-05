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
#include <gsl/gsl_linalg.h>

// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/dem/distance.h>
#include <mechsys/util/fatal.h>

using std::cout;
using std::endl;

struct UserData
{
    DEM::Particle *    p;            // the array of particles at which the force is to be applied
    Array<Vec3_t  >    vm0;          // value of the vectors close to the middle section
    Array<Vec3_t *>    vm;           // pointers to the vectors close to the middle section
    String             test;         // Type of test vibraiton or tension
    double             A;            // Area of the plate for stress calculation
    double             Am;           // vibration amplitude
    double             ome;          // vibration frequency
    double             sy;           // Stress State
    double             Tf;           // Final time
    Vec3_t             L0;           // Initial dimensions
    std::ofstream      oss_ss;       // file for stress strain data
};

void Setup (DEM::Domain & Dom, void * UD)
{
    // force at particle tagged -3
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.test=="vibration")    dat.p->Ff=0.0,0.0,dat.Am*sin(dat.ome*Dom.Time);
    if (dat.test=="bending")
    {
        if (Dom.Time < 0.5*dat.Tf) dat.p->Ff=0.0,0.0,dat.Am*2*Dom.Time/dat.Tf;
        else                       dat.p->Ff=0.0,0.0,dat.Am;
    }
    dat.sy = dat.p->F(1)/dat.A;
}

void Report (DEM::Domain & Dom, void * UD)
{ 
    UserData & dat = (*static_cast<UserData *>(UD));
    if (dat.test=="tensile")
    {
        if (Dom.idx_out==0)
        {
            String fs;
            fs.Printf("%s_walls.res",Dom.FileKey.CStr());
            dat.oss_ss.open(fs.CStr());
            dat.sy = dat.p->F(1)/dat.A;
            // Output of the current time, the stress state sx, and the strains ex,ey and ez
            dat.oss_ss << Util::_10_6 << "Time" << Util::_8s << "sy" << Util::_8s << "ex" << Util::_8s << "ey" << Util::_8s << "ez" << std::endl;
        }
        if (!Dom.Finished)
        {
            // Measure the strains and stresses
            Vec3_t Xmin, Xmax;
            Dom.BoundingBox(Xmin, Xmax);
            double ex = (Xmax(0)-Xmin(0)-dat.L0(0))/dat.L0(0);
            double ey = (Xmax(1)-Xmin(1)-dat.L0(1))/dat.L0(1);
            double ez = (Xmax(2)-Xmin(2)-dat.L0(2))/dat.L0(2);
            dat.oss_ss << Util::_10_6 << Dom.Time << Util::_8s << dat.sy << Util::_8s << ex << Util::_8s << ey << Util::_8s << ez << std::endl;
        }
        else dat.oss_ss.close();
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
bool wayToSort(double i, double j) { return i > j; }

int main(int argc, char **argv) try
{
  // DEM::Domain dom;
  // dom.Alpha = 0.05;
  // dom.Dilate= true;
  // // dom.AddVoroPack (-1, R, Lx,Ly,Lz, nx,ny,nz, rho,/*Cohhesion*/ true, false, seed, 1.0);
  // dom.AddVoroPack (-1, 0.1, 5.0,50.0,5.0, 5,50,5, 0.3,/*Cohhesion*/ true, false, 1000, 1.0);

  // //set the element properties
  // // Dict B;
  // // B.Set(-1,"Kn Kt Bn Bt Bm Gn Gt Mu Eps",Kn ,Kt ,Bn ,Bt ,Bm ,Gn ,Gt ,Mu ,Eps);
  // // dom.SetProps(B);
  // // dom.Solve (Tf,dt,dtOut, &Setup, &Report, filekey.CStr(), RenderVideo, Nproc);
  // dom.Solve (10.0,5.0e-5,0.1, NULL, NULL, "voro_error", 2, 6);
  Array<double> test(0);
  test.Push(1.2); test.Push(1.4); test.Push(1.6); test.Push(1.4); test.Push(2.0);
  Array<double> test2 = test;
  // test2.Sort();
  // test2.Sort(/*ascending*/std::greater<double>());
  test2.Sort(/*ascending*/wayToSort);
  std::cout<<"test "<<test<<"\ntest2 "<<test2<<std::endl;
}
MECHSYS_CATCH
