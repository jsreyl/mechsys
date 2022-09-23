#ifndef BRICK_H
#define BRICK_H

#include <math.h>
#include <random>

#include <gsl/gsl_linalg.h>
// MechSys
#include <mechsys/dem/domain.h>
#include <mechsys/util/fatal.h>
#include <mechsys/util/util.h>
#include <mechsys/mesh/unstructured.h>
#include <mechsys/mesh/structured.h>
#include <mechsys/linalg/matvec.h>

using std::cout;
using std::endl;

//------------------------------Define functions------------------------------
//Utility
void LocalImposeParticleCohesion(int Tag, DEM::Domain & Dom, double tol1=1.0e-8, double tol2=1.0e-3);
void LocalImposeParticleCohesionBetween(int Tag1, int Tag2, DEM::Domain & Dom, double tol1=1.0e-8, double tol2=1.0e-3);
void GetParticlesByPos(int Tag, double xpos, double err, DEM::Domain & Dom, Array<DEM::Particle*> & P);
//Particle addition
void AddBrickOneHole(DEM::Domain & dom, size_t iReg, size_t iTag, double Lx, double Ly, double Lz, double dx, double Vmax, double R, double rho, Vec3_t delta=OrthoSys::O);
void AddBrickOneHoleVoro(DEM::Domain & dom, size_t iReg, size_t iTag, double Lx, double Ly, double Lz, double dx, size_t nx, size_t ny, size_t nz, size_t ndx, size_t seed, double Vmax, double R, double rho, Vec3_t delta=OrthoSys::O);
void AddBrickSixHolesVoro(DEM::Domain & dom, size_t iReg, size_t iTag, double Lx, double Ly, double Lz, double dx, size_t nx, size_t ny, size_t nz, size_t ndx, size_t seed, double Vmax, double R, double rho, size_t ZisUp=0, Vec3_t delta=OrthoSys::O);
void AddBrickSixHoles(DEM::Domain & dom, size_t iReg, size_t iTag, double iLx, double iLy, double iLz, double dx, double Vmax, double R, double rho, size_t ZisUp=0, Vec3_t delta=OrthoSys::O);
void AddVoroPackAt(DEM::Domain & dom, int Tag, double R, double Lx,double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho, bool Cohesion,  bool Periodic, size_t seed, double fraction, Vec3_t delta=OrthoSys::O);

//-------------------IMPLEMENTATION----------------------------------------
//----------------------------UTILITY---------------------------------
inline  void LocalImposeParticleCohesion(int Tag, DEM::Domain & Dom, double tol1, double tol2){
  /*
   *This function imposes a cohesive interaction between all particles in the given domain
   *Inputs:
   - Tag, of the particles
   - Dom, domain of the particles
   - tolerances for angle and distance separation between faces
   *Outputs: None
   */
  Util::Stopwatch stopwatch; //Count how much time is spent generating the cohesions
  for (size_t i=0;i<Dom.Particles.Size()-1;i++)//Loop over all particles
    {
      DEM::Particle * P1 = Dom.Particles[i];
      for (size_t j=i+1;j<Dom.Particles.Size();j++)
        {
          DEM::Particle * P2 = Dom.Particles[j];
          double R =0.5*(P1->Props.R+P2->Props.R);
          //Debugging
          //std::cout<<"We're particles index "<<P1->Index<<" and "<<P2->Index<<std::endl;
          //std::cout<<"And our radiuses are "<<P1->Props.R<<" and "<<P2->Props.R<<std::endl;
          //std::cout<<"Therefore our mean radius is "<<R<<std::endl;
          //std::cout<<"And our Dmaxes are "<<P1->Dmax<<" and "<<P2->Dmax<<std::endl;
          //std::cout<<"And our Positions are "<<P1->x<<" and "<<P2->x<<std::endl;
          if (DEM::Distance(P1->x,P2->x)<P1->Dmax+P2->Dmax)//If the particles are close together
            {
              //std::cout<<"And we're in contact!"<<std::endl;
              for (size_t k=0;k<P1->Faces.Size();k++)//Look every face
                {
                  DEM::Face * F1 = P1->Faces[k];
                  Vec3_t n1,c1;
                  F1->Normal  (n1);
                  F1->Centroid(c1);
                  bool found = false;
                  //std::cout<<"I'm face "<<k<<" from Particle "<<P1->Index<<std::endl;
                for (size_t l=0;l<P2->Faces.Size();l++)
                    {
                      DEM::Face * F2 = P2->Faces[l];
                      Vec3_t n2,c2;
                      F2->Normal  (n2);
                      F2->Centroid(c2);
                      Vec3_t n = 0.5*(n1-n2);
                      n/=norm(n);
                      // std::cout<<"We're Face "<<k<<" from Particle "<<P1->Index<<" and Face"<<l<<" from Particle "<<P2->Index<<std::endl;
                      // std::cout<<"And we're angled one to the other? "<<(fabs(dot(n1,n2)+1.0)<tol1)<<std::endl;
                      // std::cout<<"And Face 2 is close to the centroid of Face 1? "<<(fabs(Distance(c1,*F2)-2*R)<tol2)<<" in fact we're "<<fabs(Distance(c1,*F2)-2*R)<<" apart"<<std::endl;
                      // std::cout<<"And Face 1 is close to the centroid of Face 2? "<<(fabs(Distance(c2,*F1)-2*R)<tol2)<<" in fact we're "<<fabs(Distance(c2,*F1)-2*R)<<" apart"<<std::endl;
                      if ((fabs(dot(n1,n2)+1.0)<tol1)
                          &&(fabs(Distance(c1,*F2)-2*R)<tol2)
                          &&(fabs(Distance(c2,*F1)-2*R)<tol2))//check for distance and angle treshholds
                        {
                          // std::cout<<"Adding cohesion between P1:"<<P1->Index<<" and P2:"<<P2->Index<<std::endl;
                          Dom.BInteractons.Push(new DEM::BInteracton(P1,P2,k,l));
                          found = true;
                          break;
                        }
                    }
                  if (found) break;
                }
            }
        }
    }
  std::cout<<"Total number of cohesion interactions: "<<Dom.BInteractons.Size()<<std::endl;
}

inline  void LocalImposeParticleCohesionBetween(int Tag1, int Tag2, DEM::Domain & Dom, double tol1, double tol2){
  /*
   *This function imposes a cohesive interaction between particles of Tag1 in contact with particles of Tag2 in the given domain
   *Inputs:
   - Tags, of the particles
   - Dom, domain of the particles
   - tolerances for angle and distance separation between faces
   *Outputs: None
   */
  Util::Stopwatch stopwatch; //Count how much time is spent generating the cohesions
  Array<DEM::Particle *> P1s, P2s;
  Dom.GetParticles(Tag1,P1s);
  Dom.GetParticles(Tag2,P2s);
  for (size_t i=0;i<P1s.Size();i++)//Loop over all Tag1 particles
    {
      DEM::Particle * P1 = P1s[i];
      for (size_t j=0;j<P2s.Size();j++)//Loop over all Tag2 particles and search it they're in contact with a Tag1 particle
        {
          DEM::Particle * P2 = P2s[j];
          double R =0.5*(P1->Props.R+P2->Props.R);
          //Debugging
          //std::cout<<"We're particles index "<<P1->Index<<" and "<<P2->Index<<std::endl;
          //std::cout<<"And our radiuses are "<<P1->Props.R<<" and "<<P2->Props.R<<std::endl;
          //std::cout<<"Therefore our mean radius is "<<R<<std::endl;
          //std::cout<<"And our Dmaxes are "<<P1->Dmax<<" and "<<P2->Dmax<<std::endl;
          //std::cout<<"And our Positions are "<<P1->x<<" and "<<P2->x<<std::endl;
          if (DEM::Distance(P1->x,P2->x)<P1->Dmax+P2->Dmax)//If the particles are close together
            {
              //std::cout<<"And we're in contact!"<<std::endl;
              for (size_t k=0;k<P1->Faces.Size();k++)//Look every face
                {
                  DEM::Face * F1 = P1->Faces[k];
                  Vec3_t n1,c1;
                  F1->Normal  (n1);
                  F1->Centroid(c1);
                  bool found = false;
                  //std::cout<<"I'm face "<<k<<" from Particle "<<P1->Index<<std::endl;
                for (size_t l=0;l<P2->Faces.Size();l++)
                    {
                      DEM::Face * F2 = P2->Faces[l];
                      Vec3_t n2,c2;
                      F2->Normal  (n2);
                      F2->Centroid(c2);
                      Vec3_t n = 0.5*(n1-n2);
                      n/=norm(n);
                      // std::cout<<"We're Face "<<k<<" from Particle "<<P1->Index<<" and Face"<<l<<" from Particle "<<P2->Index<<std::endl;
                      // std::cout<<"And we're angled one to the other? "<<(fabs(dot(n1,n2)+1.0)<tol1)<<std::endl;
                      // std::cout<<"And Face 2 is close to the centroid of Face 1? "<<(fabs(Distance(c1,*F2)-2*R)<tol2)<<" in fact we're "<<fabs(Distance(c1,*F2)-2*R)<<" apart"<<std::endl;
                      // std::cout<<"And Face 1 is close to the centroid of Face 2? "<<(fabs(Distance(c2,*F1)-2*R)<tol2)<<" in fact we're "<<fabs(Distance(c2,*F1)-2*R)<<" apart"<<std::endl;
                      if ((fabs(dot(n1,n2)+1.0)<tol1)
                          &&(fabs(Distance(c1,*F2)-2*R)<tol2)
                          &&(fabs(Distance(c2,*F1)-2*R)<tol2))//check for distance and angle treshholds
                        {
                          // std::cout<<"Adding cohesion between P1:"<<P1->Index<<" and P2:"<<P2->Index<<std::endl;
                          Dom.BInteractons.Push(new DEM::BInteracton(P1,P2,k,l));
                          found = true;
                          break;
                        }
                    }
                  if (found) break;
                }
            }
        }
    }
  std::cout<<"Total number of cohesion interactions: "<<Dom.BInteractons.Size()<<std::endl;
}

inline void GetParticlesByPos(int Tag, double xpos, double err, DEM::Domain & Dom, Array<DEM::Particle*> & P){
  /*
   *This function returns the particles within err distance of a position xpos.
   */
  P.Resize(0);
  int n_center=0;
  for(size_t i=0; i<Dom.Particles.Size();i++){
    if((Dom.Particles[i]->Tag==Tag) && (fabs(Dom.Particles[i]->x(1) - xpos) < err)){
      P.Push(Dom.Particles[i]);
      n_center++;
    }
  }
  std::cout<<"Number of particles in the beam's center: "<<n_center<<std::endl;
}

//---------------------------PARTICLE GENERATION----------------------

void AddVoroPackAt(DEM::Domain & dom, int Tag, double R, double Lx,double Ly, double Lz, size_t nx, size_t ny, size_t nz, double rho, bool Cohesion,  bool Periodic, size_t seed, double fraction, Vec3_t delta){
  dom.AddVoroPack(Tag-1000, R, Lx, Ly, Lz, nx, ny, nz, rho, Cohesion, Periodic, seed, fraction);
  Array<DEM::Particle *> Pvoro;
  dom.GetParticles(/*Tag*/Tag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(delta);
    Pvoro[i]->Tag=Tag;
  }
}
inline void AddBrickOneHole(DEM::Domain & dom, size_t iReg, size_t iTag, double Lx, double Ly, double Lz, double dx, double Vmax, double R, double rho, Vec3_t delta){
  /*
   *This function generates a brick with one hole and divides it in spherothetraedral DEM particles
   *INPUTS:
   - dom: DEM domain for the simulation
   - iReg: integer for specifying a tetgen region
   - iTag: tag for the DEM particles generated
   - Lx, Ly, Lz: dimensions of the brick
   - dx: thickness of the brick walls
   - Vmax: maximum volume for each cell, this is the mesh definition
   - R: spheroradious for the particles
   - rho: density for the particles
   - delta: origin of the block, (0,0,0) by default
   OUTPUTS:
   - None, generates DEM particles in the domain
   */
  Mesh::Unstructured mesh(/*NDim*/3);
  mesh.Set(16/*number of points to define the goemetry*/,16/*number of faces to define it*/,1/*number of regions*/,0/*number of holes*/);

  // double Vmax = 0.001; //maximun volume of cells, controls resolution
  mesh.SetReg(iReg,iTag,Vmax,0.01/*Lx*/,0.01/*Ly*/,0.01/*Lz*/); // The values of L must be apoint inside the region to be meshed

  //Defining the points of the geometry
  Vec3_t x00 = Vec3_t(0.0  , 0.0  , 0.0)+delta;
  Vec3_t x01 = Vec3_t(Lx   , 0.0  , 0.0)+delta;
  Vec3_t x02 = Vec3_t(Lx   , Ly   , 0.0)+delta;
  Vec3_t x03 = Vec3_t(0.0  , Ly   , 0.0)+delta;
  Vec3_t x04 = Vec3_t(dx   , dx   , 0.0)+delta;
  Vec3_t x05 = Vec3_t(Lx-dx, dx   , 0.0)+delta;
  Vec3_t x06 = Vec3_t(Lx-dx, Ly-dx, 0.0)+delta;
  Vec3_t x07 = Vec3_t(dx   , Ly-dx, 0.0)+delta;
  Vec3_t x08 = Vec3_t(0.0  , 0.0  , Lz )+delta;
  Vec3_t x09 = Vec3_t(Lx   , 0.0  , Lz )+delta;
  Vec3_t x10 = Vec3_t(Lx   , Ly   , Lz )+delta;
  Vec3_t x11 = Vec3_t(0.0  , Ly   , Lz )+delta;
  Vec3_t x12 = Vec3_t(dx   , dx   , Lz )+delta;
  Vec3_t x13 = Vec3_t(Lx-dx, dx   , Lz )+delta;
  Vec3_t x14 = Vec3_t(Lx-dx, Ly-dx, Lz )+delta;
  Vec3_t x15 = Vec3_t(dx   , Ly-dx, Lz )+delta;
  mesh.SetPnt( 0/*index of the point*/,-1/*tag of the point*/,x00(0),x00(1),x00(2)); 
  mesh.SetPnt( 1/*index of the point*/,-1/*tag of the point*/,x01(0),x01(1),x01(2)); 
  mesh.SetPnt( 2/*index of the point*/,-1/*tag of the point*/,x02(0),x02(1),x02(2)); 
  mesh.SetPnt( 3/*index of the point*/,-1/*tag of the point*/,x03(0),x03(1),x03(2)); 
  mesh.SetPnt( 4/*index of the point*/,-1/*tag of the point*/,x04(0),x04(1),x04(2)); 
  mesh.SetPnt( 5/*index of the point*/,-1/*tag of the point*/,x05(0),x05(1),x05(2)); 
  mesh.SetPnt( 6/*index of the point*/,-1/*tag of the point*/,x06(0),x06(1),x06(2)); 
  mesh.SetPnt( 7/*index of the point*/,-1/*tag of the point*/,x07(0),x07(1),x07(2)); 
  mesh.SetPnt( 8/*index of the point*/,-1/*tag of the point*/,x08(0),x08(1),x08(2)); 
  mesh.SetPnt( 9/*index of the point*/,-1/*tag of the point*/,x09(0),x09(1),x09(2)); 
  mesh.SetPnt(10/*index of the point*/,-1/*tag of the point*/,x10(0),x10(1),x10(2)); 
  mesh.SetPnt(11/*index of the point*/,-1/*tag of the point*/,x11(0),x11(1),x11(2)); 
  mesh.SetPnt(12/*index of the point*/,-1/*tag of the point*/,x12(0),x12(1),x12(2)); 
  mesh.SetPnt(13/*index of the point*/,-1/*tag of the point*/,x13(0),x13(1),x13(2)); 
  mesh.SetPnt(14/*index of the point*/,-1/*tag of the point*/,x14(0),x14(1),x14(2)); 
  mesh.SetPnt(15/*index of the point*/,-1/*tag of the point*/,x15(0),x15(1),x15(2)); 

  //Setting the faces by indexes of the points
  mesh.SetFac( 0, -2, Array<int>(0,1,9,8)/*array of indexes of the points defined before*/);
  mesh.SetFac( 1, -2, Array<int>(0,3,11,8)/*array of indexes of the points defined before*/);
  mesh.SetFac( 2, -2, Array<int>(1,9,10,2)/*array of indexes of the points defined before*/);
  mesh.SetFac( 3, -2, Array<int>(2,3,11,10)/*array of indexes of the points defined before*/);
  mesh.SetFac( 4, -2, Array<int>(4,5,13,12)/*array of indexes of the points defined before*/);
  mesh.SetFac( 5, -2, Array<int>(4,7,15,12)/*array of indexes of the points defined before*/);
  mesh.SetFac( 6, -2, Array<int>(5,6,14,13)/*array of indexes of the points defined before*/);
  mesh.SetFac( 7, -2, Array<int>(6,7,15,14)/*array of indexes of the points defined before*/);
  mesh.SetFac( 8, -2, Array<int>(0,4,7,3)/*array of indexes of the points defined before*/);
  mesh.SetFac( 9, -2, Array<int>(0,1,5,4)/*array of indexes of the points defined before*/);
  mesh.SetFac(10, -2, Array<int>(1,2,6,5)/*array of indexes of the points defined before*/);
  mesh.SetFac(11, -2, Array<int>(2,3,7,6)/*array of indexes of the points defined before*/);
  mesh.SetFac(12, -2, Array<int>(8,12,15,11)/*array of indexes of the points defined before*/);
  mesh.SetFac(13, -2, Array<int>(8,9,13,12)/*array of indexes of the points defined before*/);
  mesh.SetFac(14, -2, Array<int>(9,10,14,13)/*array of indexes of the points defined before*/);
  mesh.SetFac(15, -2, Array<int>(10,11,15,14)/*array of indexes of the points defined before*/);

  //Generate the mesh
  mesh.Generate();
  //Translate mesh into DEM particles
  //double R = 0.001;          // spheroradius
  //double rho = 1600.0;      // density of material
  dom.GenFromMesh (mesh,/*R*/R,/*rho*/rho,true/*Cohesion*/,/*MC*/false);
}

inline void AddBrickOneHoleVoro(DEM::Domain & dom, size_t iReg, size_t iTag, double Lx, double Ly, double Lz, double dx, size_t nx, size_t ny, size_t nz, size_t ndx, size_t seed, double Vmax, double R, double rho, Vec3_t delta){
  /*
   *This function generates a brick with one hole and divides it in voronoï spheropoliedra DEM particles, more specifically we generte four walls and stick them together
   *INPUTS:
   - dom: DEM domain for the simulation
   - iReg: integer for specifying a tetgen region (not used but left for consistency)
   - iTag: initial tag for the DEM particles generated, up to iTag-3 for each of the for walls
   - Lx, Ly, Lz: dimensions of the brick
   - dx: thickness of the brick walls
   - nx, ny, nz, ndx: number of particles for the diferent dimensions
   - seed: random seed used for voronoi particle generation
   - Vmax: maximum volume for each cell, this is the mesh definition (not used but left for consistency)
   - R: spheroradious for the particles
   - rho: density for the particles
   - delta: origin of the block, (0,0,0) by default
   OUTPUTS:
   - None, generates DEM particles in the domain
  */

  //Add four Voronoi walls and move them so they build a brick
  //Note we use Cohesion=false so we don't add cohesion twice when using LocalImposeParticleCohesion
  //WALL 1
  dom.AddVoroPack (/*Tag*/iTag, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,dx,  /*Number of cells*/nx,ny,ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  Array<DEM::Particle *> Pvoro;
  dom.GetParticles(/*Tag*/iTag,Pvoro);
  Vec3_t trans=delta;//Traslation vector
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans);
  }
  //WALL 2
  dom.AddVoroPack (/*Tag*/iTag-1, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,Lz-2*dx,  /*Number of cells*/nx,ndx,nz-2*ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1,Pvoro);
  trans=Vec3_t(0.,-(Ly-dx)/2.,(Lz-dx)/2.)+delta;//Traslation vector
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans);
  }
  //WALL 3
  dom.AddVoroPack (/*Tag*/iTag-2, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,Lz-2*dx,  /*Number of cells*/nx,ndx,nz-2*ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  dom.GetParticles(/*Tag*/iTag-2,Pvoro);
  trans = Vec3_t(0.,(Ly-dx)/2.,(Lz-dx)/2.)+delta;//Traslation vector
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans);
  }
  //WALL 4
  dom.AddVoroPack (/*Tag*/iTag-3, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,dx,  /*Number of cells*/nx,ny,ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  dom.GetParticles(/*Tag*/iTag-3,Pvoro);
  trans=Vec3_t(0.,0.,Lz-dx)+delta;//Traslation vector
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans);
  }
  //Impose cohesion
  LocalImposeParticleCohesion(-1,dom);
}

inline void AddBrickSixHolesVoro(DEM::Domain & dom, size_t iReg, size_t iTag, double Lx, double Ly, double Lz, double dx, size_t nx, size_t ny, size_t nz, size_t ndx, size_t seed, double Vmax, double R, double rho, size_t ZisUp, Vec3_t delta){
  /*
   *This function generates a brick with one hole and divides it in voronoï spheropoliedra DEM particles, more specifically we generte four walls and stick them together
   *INPUTS:
   - dom: DEM domain for the simulation
   - iReg: integer for specifying a tetgen region (not used but left for consistency)
   - iTag: initial tag for the DEM particles generated, up to iTag-3 for each of the for walls
   - Lx, Ly, Lz: dimensions of the brick
   - dx: thickness of the brick walls
   - nx, ny, nz, ndx: number of particles for the diferent dimensions
   - seed: random seed used for voronoi particle generation
   - Vmax: maximum volume for each cell, this is the mesh definition (not used but left for consistency)
   - R: spheroradious for the particles
   - rho: density for the particles
   - ZisUp: brick orientation, if true the brick will be standing with the z axis as upward direction (holes on the zx plane)
   - delta: origin of the block, (0,0,0) by default
   OUTPUTS:
   - None, generates DEM particles in the domain
  */

  //Add four Voronoi walls and move them so they build a brick
  //Note we use Cohesion=false so we don't add cohesion twice when using LocalImposeParticleCohesion
  double dLy=(Ly-4*dx)/3.;//vertical size of the holes
  // size_t dny=(ny-4*ndx)/3;
  double dLz=(Lz-3*dx)/2.;//horizontal size of the holes
  size_t dnz=(nz-3*ndx)/2;
  // delta=delta+Vec3_t(0.0,0.0,-(Lz-dx)/2.); //This centers the brick
  //Define the lengths and positions of each of the nine walls that make up the brick
  double w1Lx=Lx, w1Ly=Ly, w1Lz=dx; size_t w1nx=nx, w1ny=ny, w1nz=ndx;
  // Vec3_t trans1 = delta;
  Vec3_t trans1 = Vec3_t(0.0,0.0,-(Lz-dx)/2.)+delta;
  double w2Lx=Lx, w2Ly=dx, w2Lz=Lz-2*dx; size_t w2nx=nx, w2ny=ndx, w2nz=nz-2*ndx;
  // Vec3_t trans2 = Vec3_t(0.,-(Ly-dx)/2.,(Lz-dx)/2.)+delta;
  Vec3_t trans2 = Vec3_t(0.,-(Ly-dx)/2.,0.)+delta;
  double w3Lx=Lx, w3Ly=dx, w3Lz=Lz-2*dx; size_t w3nx=nx, w3ny=ndx, w3nz=nz-2*ndx;
  // Vec3_t trans3 = Vec3_t(0.,(Ly-dx)/2.,(Lz-dx)/2.)+delta;
  Vec3_t trans3 = Vec3_t(0.,(Ly-dx)/2.,0.)+delta;
  double w4Lx=Lx, w4Ly=Ly, w4Lz=dx; size_t w4nx=nx, w4ny=ny,w4nz=ndx;
  // Vec3_t trans4 = Vec3_t(0.,0.,Lz-dx)+delta;
  Vec3_t trans4 = Vec3_t(0.,0.,(Lz-dx)/2.)+delta;
  double w5Lx=Lx, w5Ly=Ly-2*dx, w5Lz=dx; size_t w5nx=nx, w5ny=ny-2*ndx,w5nz=ndx;
  // Vec3_t trans5 = Vec3_t(0.,0.,(Lz-dx)/2.)+delta;
  Vec3_t trans5 = delta;
  double w6Lx=Lx, w6Ly=dx, w6Lz=dLz; size_t w6nx=nx, w6ny=ndx,w6nz=dnz;
  // Vec3_t trans6 = Vec3_t(0.,dLy+dx-(Ly-dx)/2.,(dLz+dx)/2.)+delta;
  Vec3_t trans6 = Vec3_t(0.,dLy+dx-(Ly-dx)/2.,(dLz-Lz)/2.+dx)+delta;
  double w7Lx=Lx, w7Ly=dx, w7Lz=dLz; size_t w7nx=nx, w7ny=ndx,w7nz=dnz;
  // Vec3_t trans7 = Vec3_t(0.,dLy+dx-(Ly-dx)/2.,3*(dLz+dx)/2.)+delta;
  Vec3_t trans7 = Vec3_t(0.,dLy+dx-(Ly-dx)/2.,(3*dLz-Lz)/2.+2*dx)+delta;
  double w8Lx=Lx, w8Ly=dx, w8Lz=dLz; size_t w8nx=nx, w8ny=ndx,w8nz=dnz;
  // Vec3_t trans8 = Vec3_t(0.,-dLy-dx+(Ly-dx)/2.,(dLz+dx)/2.)+delta;
  Vec3_t trans8 = Vec3_t(0.,-dLy-dx+(Ly-dx)/2.,(dLz-Lz)/2.+dx)+delta;
  double w9Lx=Lx, w9Ly=dx, w9Lz=dLz; size_t w9nx=nx, w9ny=ndx,w9nz=dnz;
  // Vec3_t trans9 = Vec3_t(0.,-dLy-dx+(Ly-dx)/2.,3*(dLz+dx)/2.)+delta;
  Vec3_t trans9 = Vec3_t(0.,-dLy-dx+(Ly-dx)/2.,(3*dLz-Lz)/2.+2*dx)+delta;
  //If the orientation is on Z, change the values accordingly
  if (ZisUp){
    //Walls L and n: change y for z
    w1Lx=Lx, w1Ly=dx, w1Lz=Ly; w1nx=nx, w1ny=ndx, w1nz=ny;
    trans1=Vec3_t(trans1(0),trans1(2),trans1(1));
    w2Lx=Lx, w2Ly=Lz-2*dx, w2Lz=dx; w2nx=nx, w2ny=nz-2*ndx, w2nz=ndx;
    trans2=Vec3_t(trans2(0),trans2(2),trans2(1));
    w3Lx=Lx, w3Ly=Lz-2*dx, w3Lz=dx; w3nx=nx, w3ny=nz-2*ndx, w3nz=ndx;
    trans3=Vec3_t(trans3(0),trans3(2),trans3(1));
    w4Lx=Lx, w4Ly=dx, w4Lz=Ly; w4nx=nx, w4ny=ndx, w4nz=ny;
    trans4=Vec3_t(trans4(0),trans4(2),trans4(1));
    w5Lx=Lx, w5Ly=dx, w5Lz=Ly-2*dx; w5nx=nx, w5ny=ndx,w5nz=ny-2*ndx;
    trans5=Vec3_t(trans5(0),trans5(2),trans5(1));
    w6Lx=Lx, w6Ly=dLz, w6Lz=dx; w6nx=nx, w6ny=dnz, w6nz=ndx;
    trans6=Vec3_t(trans6(0),trans6(2),trans6(1));
    w7Lx=Lx, w7Ly=dLz, w7Lz=dx; w7nx=nx, w7ny=dnz, w7nz=ndx;
    trans7=Vec3_t(trans7(0),trans7(2),trans7(1));
    w8Lx=Lx, w8Ly=dLz, w8Lz=dx; w8nx=nx, w8ny=dnz, w8nz=ndx;
    trans8=Vec3_t(trans8(0),trans8(2),trans8(1));
    w9Lx=Lx, w9Ly=dLz, w9Lz=dx; w9nx=nx, w9ny=dnz, w9nz=ndx;
    trans9=Vec3_t(trans9(0),trans9(2),trans9(1));
  }

  //WALL 1
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w1Lx,w1Ly,w1Lz,  /*Number of cells*/w1nx,w1ny,w1nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  Array<DEM::Particle *> Pvoro;
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans1);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 2
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w2Lx,w2Ly,w2Lz,  /*Number of cells*/w2nx,w2ny,w2nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans2);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 3
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w3Lx,w3Ly,w3Lz,  /*Number of cells*/w3nx,w3ny,w3nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans3);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 4
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w4Lx,w4Ly,w4Lz,  /*Number of cells*/w4nx,w4ny,w4nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans4);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 5
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w5Lx,w5Ly,w5Lz,  /*Number of cells*/w5nx,w5ny,w5nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans5);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 6
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w6Lx,w6Ly,w6Lz,  /*Number of cells*/w6nx,w6ny,w6nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans6);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 7
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w7Lx,w7Ly,w7Lz,  /*Number of cells*/w7nx,w7ny,w7nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans7);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 8
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w8Lx,w8Ly,w8Lz,  /*Number of cells*/w8nx,w8ny,w8nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans8);
    Pvoro[i]->Tag=iTag;
  }
  //WALL 9
  dom.AddVoroPack (/*Tag*/iTag-1000, /*Spheroradious*/R,  /*Dimentions*/w9Lx,w9Ly,w9Lz,  /*Number of cells*/w9nx,w9ny,w9nz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  //Move the wall so it doesn't intersect with the first one
  dom.GetParticles(/*Tag*/iTag-1000,Pvoro);
  for(size_t i=0;i<Pvoro.Size();i++){
    Pvoro[i]->Translate(trans9);
    Pvoro[i]->Tag=iTag;
  }

  // //WALL 1
  // dom.AddVoroPack (/*Tag*/iTag, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,dx,  /*Number of cells*/nx,ny,ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // Array<DEM::Particle *> Pvoro;
  // dom.GetParticles(/*Tag*/iTag,Pvoro);
  // Vec3_t trans=delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  // }
  // //WALL 2
  // dom.AddVoroPack (/*Tag*/iTag-1, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,Lz-2*dx,  /*Number of cells*/nx,ndx,nz-2*ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // //Move the wall so it doesn't intersect with the first one
  // dom.GetParticles(/*Tag*/iTag-1,Pvoro);
  // trans=Vec3_t(0.,-(Ly-dx)/2.,(Lz-dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 3
  // dom.AddVoroPack (/*Tag*/iTag-2, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,Lz-2*dx,  /*Number of cells*/nx,ndx,nz-2*ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // dom.GetParticles(/*Tag*/iTag-2,Pvoro);
  // trans = Vec3_t(0.,(Ly-dx)/2.,(Lz-dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 4
  // dom.AddVoroPack (/*Tag*/iTag-3, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly,dx,  /*Number of cells*/nx,ny,ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // dom.GetParticles(/*Tag*/iTag-3,Pvoro);
  // trans=Vec3_t(0.,0.,Lz-dx)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 5
  // dom.AddVoroPack (/*Tag*/iTag-4, /*Spheroradious*/R,  /*Dimentions*/Lx,Ly-2*dx,dx,  /*Number of cells*/nx,ny-2*ndx,ndx,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // dom.GetParticles(/*Tag*/iTag-4,Pvoro);
  // trans=Vec3_t(0.,0.,(Lz-dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 6
  // dom.AddVoroPack (/*Tag*/iTag-5, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,dLz, /*Number of cells*/nx,ndx,dnz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // //Move the wall so it doesn't intersect with the first one
  // dom.GetParticles(/*Tag*/iTag-5,Pvoro);
  // trans=Vec3_t(0.,dLy+dx-(Ly-dx)/2.,(dLz+dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 7
  // dom.AddVoroPack (/*Tag*/iTag-6, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,dLz, /*Number of cells*/nx,ndx,dnz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // //Move the wall so it doesn't intersect with the first one
  // dom.GetParticles(/*Tag*/iTag-6,Pvoro);
  // trans=Vec3_t(0.,dLy+dx-(Ly-dx)/2.,3*(dLz+dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 8
  // dom.AddVoroPack (/*Tag*/iTag-7, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,dLz, /*Number of cells*/nx,ndx,dnz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // //Move the wall so it doesn't intersect with the first one
  // dom.GetParticles(/*Tag*/iTag-7,Pvoro);
  // trans=Vec3_t(0.,-dLy-dx+(Ly-dx)/2.,(dLz+dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  // //WALL 9
  // dom.AddVoroPack (/*Tag*/iTag-8, /*Spheroradious*/R,  /*Dimentions*/Lx,dx,dLz, /*Number of cells*/nx,ndx,dnz,  /*Density*/rho,  /*Cohesion*/false,  /*Periodic*/false,  /*Rand seed*/seed,  /*fraction*/1.0);
  // //Move the wall so it doesn't intersect with the first one
  // dom.GetParticles(/*Tag*/iTag-8,Pvoro);
  // trans=Vec3_t(0.,-dLy-dx+(Ly-dx)/2.,3*(dLz+dx)/2.)+delta;//Traslation vector
  // for(size_t i=0;i<Pvoro.Size();i++){
  //   Pvoro[i]->Translate(trans);
  //   Pvoro[i]->Tag=iTag;
  // }
  //Impose cohesion
  LocalImposeParticleCohesion(iTag,dom);
}

inline void AddBrickSixHoles(DEM::Domain & dom, size_t iReg, size_t iTag, double iLx, double iLy, double iLz, double dx, double Vmax, double R, double rho, size_t ZisUp, Vec3_t delta){
  /*
   *This functions generates a brick with 6 holes from a mesh defined by the points and their connection into faces. This mesh is then converted into DEM sphero-thetraedral particles.
   INPUT:
   dom: DEM domain where particles are to be added
   iReg, iTag: id and tag for the region, this tag will be used for the particles generated.
   Lx, Ly, Lz: brick dimentions
   dx: spacing between brick holes
   Vmax: maximum volume for the mesh cells, this controls resolution
   R: spheroradius for particles
   rho: density for particles
   delta: displacement vector to add for the mesh
   ZisUp: bool, whether to use z as vertical axis or y (default)
   */
  //size_t IIndex = dom.Particles.Size();//Number of particles before generation
  std::cout<< "Input dimentions for tetra mesh generation: Lx "<<iLx << " Ly "<< iLy << " Lz " << iLz << " Vmax " << Vmax<<"\n";
  double Lx = iLz; double Ly = iLy; double Lz = iLx;

  Mesh::Unstructured mesh(/*NDim*/3);
  mesh.Set(56/*number of points to define the geometry*/,46/*number of faces to define it*/,1/*number of regions*/,0/*number of holes*/);

  //double Vmax = 0.0001; //maximun volume of cells, controls resolution
  Vec3_t origin = Vec3_t(-Lz/2.,-Ly/2.,-Lx/2.)+delta; // Origin vector to displace all vertices by so the center of the brick is at delta
  if (ZisUp)
    origin = Vec3_t(origin(0),origin(2),origin(1));
  mesh.SetReg(/*Id*/iReg,/*Tag*/iTag,Vmax,0.01+origin[0]/*Lx*/,0.01+origin[1]/*Ly*/,0.01+origin[2]/*Lz*/); // The values of L must be a point inside the region to be meshed

  //Defining the points of the geometry
  //double Lx=1., Ly=3., Lz=1.;
  //double dx=0.1;//Size of the divisions
  double dLy=(Ly-4.*dx)/3.;//vertical size of the holes.
  Vec3_t x00 = Vec3_t(0.0 , 0.0        , 0.0       )+origin;
  Vec3_t x01 = Vec3_t(0.0 , 0.0        , Lx        )+origin;
  Vec3_t x02 = Vec3_t(0.0 , Ly         , Lx        )+origin;
  Vec3_t x03 = Vec3_t(0.0 , Ly         , 0.0       )+origin;
  Vec3_t x04 = Vec3_t(0.0 , dx         , dx        )+origin;
  Vec3_t x05 = Vec3_t(0.0 , dx         , (Lx-dx)/2.)+origin;
  Vec3_t x06 = Vec3_t(0.0 , dx         , (Lx+dx)/2.)+origin;
  Vec3_t x07 = Vec3_t(0.0 , dx         , Lx-dx     )+origin;
  Vec3_t x08 = Vec3_t(0.0 , dx+dLy     , Lx-dx     )+origin;
  Vec3_t x09 = Vec3_t(0.0 , dx+dLy     , (Lx+dx)/2.)+origin;
  Vec3_t x10 = Vec3_t(0.0 , dx+dLy     , (Lx-dx)/2.)+origin;
  Vec3_t x11 = Vec3_t(0.0 , dx+dLy     , dx        )+origin;
  Vec3_t x12 = Vec3_t(0.0 , 2*dx+dLy   , dx        )+origin;
  Vec3_t x13 = Vec3_t(0.0 , 2*dx+dLy   , (Lx-dx)/2.)+origin;
  Vec3_t x14 = Vec3_t(0.0 , 2*dx+dLy   , (Lx+dx)/2.)+origin;
  Vec3_t x15 = Vec3_t(0.0 , 2*dx+dLy   , Lx-dx     )+origin;
  Vec3_t x16 = Vec3_t(0.0 , 2*dx+2*dLy , Lx-dx     )+origin;
  Vec3_t x17 = Vec3_t(0.0 , 2*dx+2*dLy , (Lx+dx)/2.)+origin;
  Vec3_t x18 = Vec3_t(0.0 , 2*dx+2*dLy , (Lx-dx)/2.)+origin;
  Vec3_t x19 = Vec3_t(0.0 , 2*dx+2*dLy , dx        )+origin;
  Vec3_t x20 = Vec3_t(0.0 , 3*dx+2*dLy , dx        )+origin;
  Vec3_t x21 = Vec3_t(0.0 , 3*dx+2*dLy , (Lx-dx)/2.)+origin;
  Vec3_t x22 = Vec3_t(0.0 , 3*dx+2*dLy , (Lx+dx)/2.)+origin;
  Vec3_t x23 = Vec3_t(0.0 , 3*dx+2*dLy , Lx-dx     )+origin;
  Vec3_t x24 = Vec3_t(0.0 , 3*dx+3*dLy , Lx-dx     )+origin;
  Vec3_t x25 = Vec3_t(0.0 , 3*dx+3*dLy , (Lx+dx)/2.)+origin;
  Vec3_t x26 = Vec3_t(0.0 , 3*dx+3*dLy , (Lx-dx)/2.)+origin;
  Vec3_t x27 = Vec3_t(0.0 , 3*dx+3*dLy , dx        )+origin;
  Vec3_t x28 = Vec3_t(Lz  , 0.0        , 0.0       )+origin;
  Vec3_t x29 = Vec3_t(Lz  , 0.0        , Lx        )+origin;
  Vec3_t x30 = Vec3_t(Lz  , Ly         , Lx        )+origin;
  Vec3_t x31 = Vec3_t(Lz  , Ly         , 0.0       )+origin;
  Vec3_t x32 = Vec3_t(Lz  , dx         , dx        )+origin;
  Vec3_t x33 = Vec3_t(Lz  , dx         , (Lx-dx)/2.)+origin;
  Vec3_t x34 = Vec3_t(Lz  , dx         , (Lx+dx)/2.)+origin;
  Vec3_t x35 = Vec3_t(Lz  , dx         , Lx-dx     )+origin;
  Vec3_t x36 = Vec3_t(Lz  , dx+dLy     , Lx-dx     )+origin;
  Vec3_t x37 = Vec3_t(Lz  , dx+dLy     , (Lx+dx)/2.)+origin;
  Vec3_t x38 = Vec3_t(Lz  , dx+dLy     , (Lx-dx)/2.)+origin;
  Vec3_t x39 = Vec3_t(Lz  , dx+dLy     , dx        )+origin;
  Vec3_t x40 = Vec3_t(Lz  , 2*dx+dLy   , dx        )+origin;
  Vec3_t x41 = Vec3_t(Lz  , 2*dx+dLy   , (Lx-dx)/2.)+origin;
  Vec3_t x42 = Vec3_t(Lz  , 2*dx+dLy   , (Lx+dx)/2.)+origin;
  Vec3_t x43 = Vec3_t(Lz  , 2*dx+dLy   , Lx-dx     )+origin;
  Vec3_t x44 = Vec3_t(Lz  , 2*dx+2*dLy , Lx-dx     )+origin;
  Vec3_t x45 = Vec3_t(Lz  , 2*dx+2*dLy , (Lx+dx)/2.)+origin;
  Vec3_t x46 = Vec3_t(Lz  , 2*dx+2*dLy , (Lx-dx)/2.)+origin;
  Vec3_t x47 = Vec3_t(Lz  , 2*dx+2*dLy , dx        )+origin;
  Vec3_t x48 = Vec3_t(Lz  , 3*dx+2*dLy , dx        )+origin;
  Vec3_t x49 = Vec3_t(Lz  , 3*dx+2*dLy , (Lx-dx)/2.)+origin;
  Vec3_t x50 = Vec3_t(Lz  , 3*dx+2*dLy , (Lx+dx)/2.)+origin;
  Vec3_t x51 = Vec3_t(Lz  , 3*dx+2*dLy , Lx-dx     )+origin;
  Vec3_t x52 = Vec3_t(Lz  , 3*dx+3*dLy , Lx-dx     )+origin;
  Vec3_t x53 = Vec3_t(Lz  , 3*dx+3*dLy , (Lx+dx)/2.)+origin;
  Vec3_t x54 = Vec3_t(Lz  , 3*dx+3*dLy , (Lx-dx)/2.)+origin;
  Vec3_t x55 = Vec3_t(Lz  , 3*dx+3*dLy , dx        )+origin;

  if (ZisUp){
    x00 = Vec3_t(0.0 , 0.0        , 0.0       )+origin;
    x01 = Vec3_t(0.0 , Lx         , 0.0       )+origin;
    x02 = Vec3_t(0.0 , Lx         , Ly        )+origin;
    x03 = Vec3_t(0.0 , 0.0        , Ly        )+origin;
    x04 = Vec3_t(0.0 , dx         , dx        )+origin;
    x05 = Vec3_t(0.0 , (Lx-dx)/2. , dx        )+origin;
    x06 = Vec3_t(0.0 , (Lx+dx)/2. , dx        )+origin;
    x07 = Vec3_t(0.0 , Lx-dx      , dx        )+origin;
    x08 = Vec3_t(0.0 , Lx-dx      , dx+dLy    )+origin;
    x09 = Vec3_t(0.0 , (Lx+dx)/2. , dx+dLy    )+origin;
    x10 = Vec3_t(0.0 , (Lx-dx)/2. , dx+dLy    )+origin;
    x11 = Vec3_t(0.0 , dx         , dx+dLy    )+origin;
    x12 = Vec3_t(0.0 , dx         , 2*dx+dLy  )+origin;
    x13 = Vec3_t(0.0 , (Lx-dx)/2. , 2*dx+dLy  )+origin;
    x14 = Vec3_t(0.0 , (Lx+dx)/2. , 2*dx+dLy  )+origin;
    x15 = Vec3_t(0.0 , Lx-dx      , 2*dx+dLy  )+origin;
    x16 = Vec3_t(0.0 , Lx-dx      , 2*dx+2*dLy)+origin;
    x17 = Vec3_t(0.0 , (Lx+dx)/2. , 2*dx+2*dLy)+origin;
    x18 = Vec3_t(0.0 , (Lx-dx)/2. , 2*dx+2*dLy)+origin;
    x19 = Vec3_t(0.0 , dx         , 2*dx+2*dLy)+origin;
    x20 = Vec3_t(0.0 , dx         , 3*dx+2*dLy)+origin;
    x21 = Vec3_t(0.0 , (Lx-dx)/2. , 3*dx+2*dLy)+origin;
    x22 = Vec3_t(0.0 , (Lx+dx)/2. , 3*dx+2*dLy)+origin;
    x23 = Vec3_t(0.0 , Lx-dx      , 3*dx+2*dLy)+origin;
    x24 = Vec3_t(0.0 , Lx-dx      , 3*dx+3*dLy)+origin;
    x25 = Vec3_t(0.0 , (Lx+dx)/2. , 3*dx+3*dLy)+origin;
    x26 = Vec3_t(0.0 , (Lx-dx)/2. , 3*dx+3*dLy)+origin;
    x27 = Vec3_t(0.0 , dx         , 3*dx+3*dLy)+origin;
    x28 = Vec3_t(Lz  , 0.0        , 0.0       )+origin;
    x29 = Vec3_t(Lz  , Lx         , 0.0       )+origin;
    x30 = Vec3_t(Lz  , Lx         , Ly        )+origin;
    x31 = Vec3_t(Lz  , 0.0        , Ly        )+origin;
    x32 = Vec3_t(Lz  , dx         , dx        )+origin;
    x33 = Vec3_t(Lz  , (Lx-dx)/2. , dx        )+origin;
    x34 = Vec3_t(Lz  , (Lx+dx)/2. , dx        )+origin;
    x35 = Vec3_t(Lz  , Lx-dx      , dx        )+origin;
    x36 = Vec3_t(Lz  , Lx-dx      , dx+dLy    )+origin;
    x37 = Vec3_t(Lz  , (Lx+dx)/2. , dx+dLy    )+origin;
    x38 = Vec3_t(Lz  , (Lx-dx)/2. , dx+dLy    )+origin;
    x39 = Vec3_t(Lz  , dx         , dx+dLy    )+origin;
    x40 = Vec3_t(Lz  , dx         , 2*dx+dLy  )+origin;
    x41 = Vec3_t(Lz  , (Lx-dx)/2. , 2*dx+dLy  )+origin;
    x42 = Vec3_t(Lz  , (Lx+dx)/2. , 2*dx+dLy  )+origin;
    x43 = Vec3_t(Lz  , Lx-dx      , 2*dx+dLy  )+origin;
    x44 = Vec3_t(Lz  , Lx-dx      , 2*dx+2*dLy)+origin;
    x45 = Vec3_t(Lz  , (Lx+dx)/2. , 2*dx+2*dLy)+origin;
    x46 = Vec3_t(Lz  , (Lx-dx)/2. , 2*dx+2*dLy)+origin;
    x47 = Vec3_t(Lz  , dx         , 2*dx+2*dLy)+origin;
    x48 = Vec3_t(Lz  , dx         , 3*dx+2*dLy)+origin;
    x49 = Vec3_t(Lz  , (Lx-dx)/2. , 3*dx+2*dLy)+origin;
    x50 = Vec3_t(Lz  , (Lx+dx)/2. , 3*dx+2*dLy)+origin;
    x51 = Vec3_t(Lz  , Lx-dx      , 3*dx+2*dLy)+origin;
    x52 = Vec3_t(Lz  , Lx-dx      , 3*dx+3*dLy)+origin;
    x53 = Vec3_t(Lz  , (Lx+dx)/2. , 3*dx+3*dLy)+origin;
    x54 = Vec3_t(Lz  , (Lx-dx)/2. , 3*dx+3*dLy)+origin;
    x55 = Vec3_t(Lz  , dx         , 3*dx+3*dLy)+origin;
  }

  mesh.SetPnt( 0/*index of the point*/,-1/*tag of the point*/,x00(0),x00(1),x00(2)); 
  mesh.SetPnt( 1/*index of the point*/,-1/*tag of the point*/,x01(0),x01(1),x01(2)); 
  mesh.SetPnt( 2/*index of the point*/,-1/*tag of the point*/,x02(0),x02(1),x02(2)); 
  mesh.SetPnt( 3/*index of the point*/,-1/*tag of the point*/,x03(0),x03(1),x03(2)); 
  mesh.SetPnt( 4/*index of the point*/,-1/*tag of the point*/,x04(0),x04(1),x04(2)); 
  mesh.SetPnt( 5/*index of the point*/,-1/*tag of the point*/,x05(0),x05(1),x05(2)); 
  mesh.SetPnt( 6/*index of the point*/,-1/*tag of the point*/,x06(0),x06(1),x06(2)); 
  mesh.SetPnt( 7/*index of the point*/,-1/*tag of the point*/,x07(0),x07(1),x07(2)); 
  mesh.SetPnt( 8/*index of the point*/,-1/*tag of the point*/,x08(0),x08(1),x08(2)); 
  mesh.SetPnt( 9/*index of the point*/,-1/*tag of the point*/,x09(0),x09(1),x09(2)); 
  mesh.SetPnt(10/*index of the point*/,-1/*tag of the point*/,x10(0),x10(1),x10(2)); 
  mesh.SetPnt(11/*index of the point*/,-1/*tag of the point*/,x11(0),x11(1),x11(2)); 
  mesh.SetPnt(12/*index of the point*/,-1/*tag of the point*/,x12(0),x12(1),x12(2)); 
  mesh.SetPnt(13/*index of the point*/,-1/*tag of the point*/,x13(0),x13(1),x13(2)); 
  mesh.SetPnt(14/*index of the point*/,-1/*tag of the point*/,x14(0),x14(1),x14(2)); 
  mesh.SetPnt(15/*index of the point*/,-1/*tag of the point*/,x15(0),x15(1),x15(2)); 
  mesh.SetPnt(16/*index of the point*/,-1/*tag of the point*/,x16(0),x16(1),x16(2)); 
  mesh.SetPnt(17/*index of the point*/,-1/*tag of the point*/,x17(0),x17(1),x17(2)); 
  mesh.SetPnt(18/*index of the point*/,-1/*tag of the point*/,x18(0),x18(1),x18(2)); 
  mesh.SetPnt(19/*index of the point*/,-1/*tag of the point*/,x19(0),x19(1),x19(2)); 
  mesh.SetPnt(20/*index of the point*/,-1/*tag of the point*/,x20(0),x20(1),x20(2)); 
  mesh.SetPnt(21/*index of the point*/,-1/*tag of the point*/,x21(0),x21(1),x21(2)); 
  mesh.SetPnt(22/*index of the point*/,-1/*tag of the point*/,x22(0),x22(1),x22(2)); 
  mesh.SetPnt(23/*index of the point*/,-1/*tag of the point*/,x23(0),x23(1),x23(2)); 
  mesh.SetPnt(24/*index of the point*/,-1/*tag of the point*/,x24(0),x24(1),x24(2)); 
  mesh.SetPnt(25/*index of the point*/,-1/*tag of the point*/,x25(0),x25(1),x25(2)); 
  mesh.SetPnt(26/*index of the point*/,-1/*tag of the point*/,x26(0),x26(1),x26(2)); 
  mesh.SetPnt(27/*index of the point*/,-1/*tag of the point*/,x27(0),x27(1),x27(2)); 
  mesh.SetPnt(28/*index of the point*/,-1/*tag of the point*/,x28(0),x28(1),x28(2)); 
  mesh.SetPnt(29/*index of the point*/,-1/*tag of the point*/,x29(0),x29(1),x29(2)); 
  mesh.SetPnt(30/*index of the point*/,-1/*tag of the point*/,x30(0),x30(1),x30(2)); 
  mesh.SetPnt(31/*index of the point*/,-1/*tag of the point*/,x31(0),x31(1),x31(2)); 
  mesh.SetPnt(32/*index of the point*/,-1/*tag of the point*/,x32(0),x32(1),x32(2)); 
  mesh.SetPnt(33/*index of the point*/,-1/*tag of the point*/,x33(0),x33(1),x33(2)); 
  mesh.SetPnt(34/*index of the point*/,-1/*tag of the point*/,x34(0),x34(1),x34(2)); 
  mesh.SetPnt(35/*index of the point*/,-1/*tag of the point*/,x35(0),x35(1),x35(2)); 
  mesh.SetPnt(36/*index of the point*/,-1/*tag of the point*/,x36(0),x36(1),x36(2)); 
  mesh.SetPnt(37/*index of the point*/,-1/*tag of the point*/,x37(0),x37(1),x37(2)); 
  mesh.SetPnt(38/*index of the point*/,-1/*tag of the point*/,x38(0),x38(1),x38(2)); 
  mesh.SetPnt(39/*index of the point*/,-1/*tag of the point*/,x39(0),x39(1),x39(2)); 
  mesh.SetPnt(40/*index of the point*/,-1/*tag of the point*/,x40(0),x40(1),x40(2)); 
  mesh.SetPnt(41/*index of the point*/,-1/*tag of the point*/,x41(0),x41(1),x41(2)); 
  mesh.SetPnt(42/*index of the point*/,-1/*tag of the point*/,x42(0),x42(1),x42(2)); 
  mesh.SetPnt(43/*index of the point*/,-1/*tag of the point*/,x43(0),x43(1),x43(2)); 
  mesh.SetPnt(44/*index of the point*/,-1/*tag of the point*/,x44(0),x44(1),x44(2)); 
  mesh.SetPnt(45/*index of the point*/,-1/*tag of the point*/,x45(0),x45(1),x45(2)); 
  mesh.SetPnt(46/*index of the point*/,-1/*tag of the point*/,x46(0),x46(1),x46(2)); 
  mesh.SetPnt(47/*index of the point*/,-1/*tag of the point*/,x47(0),x47(1),x47(2)); 
  mesh.SetPnt(48/*index of the point*/,-1/*tag of the point*/,x48(0),x48(1),x48(2)); 
  mesh.SetPnt(49/*index of the point*/,-1/*tag of the point*/,x49(0),x49(1),x49(2)); 
  mesh.SetPnt(50/*index of the point*/,-1/*tag of the point*/,x50(0),x50(1),x50(2)); 
  mesh.SetPnt(51/*index of the point*/,-1/*tag of the point*/,x51(0),x51(1),x51(2)); 
  mesh.SetPnt(52/*index of the point*/,-1/*tag of the point*/,x52(0),x52(1),x52(2)); 
  mesh.SetPnt(53/*index of the point*/,-1/*tag of the point*/,x53(0),x53(1),x53(2)); 
  mesh.SetPnt(54/*index of the point*/,-1/*tag of the point*/,x54(0),x54(1),x54(2)); 
  mesh.SetPnt(55/*index of the point*/,-1/*tag of the point*/,x55(0),x55(1),x55(2)); 

  //Setting faces by indexes of the points
  //External faces
  mesh.SetFac( 0, -2, Array<int>(0,1,29,28)/*array of indexes of the points defined before*/);
  mesh.SetFac( 1, -2, Array<int>(0,3,31,28)/*array of indexes of the points defined before*/);
  mesh.SetFac( 2, -2, Array<int>(1,29,30,2)/*array of indexes of the points defined before*/);
  mesh.SetFac( 3, -2, Array<int>(2,3,31,30)/*array of indexes of the points defined before*/);

  //Internal faces
  mesh.SetFac( 4, -2, Array<int>(4,5,33,32)/*array of indexes of the points defined before*/);
  mesh.SetFac( 5, -2, Array<int>(4,32,39,11)/*array of indexes of the points defined before*/);
  mesh.SetFac( 6, -2, Array<int>(5,10,38,33)/*array of indexes of the points defined before*/);
  mesh.SetFac( 7, -2, Array<int>(10,11,39,38)/*array of indexes of the points defined before*/);
  mesh.SetFac( 8, -2, Array<int>(6,7,35,34)/*array of indexes of the points defined before*/);
  mesh.SetFac( 9, -2, Array<int>(6,9,37,34)/*array of indexes of the points defined before*/);
  mesh.SetFac(10, -2, Array<int>(7,8,36,35)/*array of indexes of the points defined before*/);
  mesh.SetFac(11, -2, Array<int>(8,9,37,36)/*array of indexes of the points defined before*/);
  mesh.SetFac(12, -2, Array<int>(12,13,41,40)/*array of indexes of the points defined before*/);
  mesh.SetFac(13, -2, Array<int>(12,40,47,19)/*array of indexes of the points defined before*/);
  mesh.SetFac(14, -2, Array<int>(13,18,46,41)/*array of indexes of the points defined before*/);
  mesh.SetFac(15, -2, Array<int>(18,19,47,46)/*array of indexes of the points defined before*/);
  mesh.SetFac(16, -2, Array<int>(14,15,43,42)/*array of indexes of the points defined before*/);
  mesh.SetFac(17, -2, Array<int>(14,42,45,17)/*array of indexes of the points defined before*/);
  mesh.SetFac(18, -2, Array<int>(15,16,44,43)/*array of indexes of the points defined before*/);
  mesh.SetFac(19, -2, Array<int>(16,17,45,44)/*array of indexes of the points defined before*/);
  mesh.SetFac(20, -2, Array<int>(20,21,49,48)/*array of indexes of the points defined before*/);
  mesh.SetFac(21, -2, Array<int>(20,48,55,27)/*array of indexes of the points defined before*/);
  mesh.SetFac(22, -2, Array<int>(21,26,54,49)/*array of indexes of the points defined before*/);
  mesh.SetFac(23, -2, Array<int>(26,27,55,54)/*array of indexes of the points defined before*/);
  mesh.SetFac(24, -2, Array<int>(22,23,51,50)/*array of indexes of the points defined before*/);
  mesh.SetFac(25, -2, Array<int>(22,50,53,25)/*array of indexes of the points defined before*/);
  mesh.SetFac(26, -2, Array<int>(23,24,52,51)/*array of indexes of the points defined before*/);
  mesh.SetFac(27, -2, Array<int>(24,25,53,52)/*array of indexes of the points defined before*/);

  mesh.SetFac(28, -2, Array<int>(0,4,11,12,19,20,27,3)/*array of indexes of the points defined before*/);
  mesh.SetFac(29, -2, Array<int>(0,1,7,6,5,4)/*array of indexes of the points defined before*/);
  mesh.SetFac(30, -2, Array<int>(1,2,24,23,16,15,8,7)/*array of indexes of the points defined before*/);
  mesh.SetFac(31, -2, Array<int>(2,3,27,26,25,24)/*array of indexes of the points defined before*/);
  mesh.SetFac(32, -2, Array<int>(5,6,9,14,17,22,25,26,21,18,13,10)/*array of indexes of the points defined before*/);
  mesh.SetFac(33, -2, Array<int>(19,18,21,20)/*array of indexes of the points defined before*/);
  mesh.SetFac(34, -2, Array<int>(17,16,23,22)/*array of indexes of the points defined before*/);
  mesh.SetFac(35, -2, Array<int>(11,10,13,12)/*array of indexes of the points defined before*/);
  mesh.SetFac(36, -2, Array<int>(9,8,15,14)/*array of indexes of the points defined before*/);

  mesh.SetFac(37, -2, Array<int>(28,32,39,40,47,48,55,31)/*array of indexes of the points defined before*/);
  mesh.SetFac(38, -2, Array<int>(28,29,35,34,33,32)/*array of indexes of the points defined before*/);
  mesh.SetFac(39, -2, Array<int>(29,30,52,51,44,43,36,35)/*array of indexes of the points defined before*/);
  mesh.SetFac(40, -2, Array<int>(30,31,55,54,53,52)/*array of indexes of the points defined before*/);
  mesh.SetFac(41, -2, Array<int>(33,34,37,42,45,50,53,54,49,46,41,38)/*array of indexes of the points defined before*/);
  mesh.SetFac(42, -2, Array<int>(47,46,49,48)/*array of indexes of the points defined before*/);
  mesh.SetFac(43, -2, Array<int>(45,44,51,50)/*array of indexes of the points defined before*/);
  mesh.SetFac(44, -2, Array<int>(39,38,41,40)/*array of indexes of the points defined before*/);
  mesh.SetFac(45, -2, Array<int>(37,36,43,42)/*array of indexes of the points defined before*/);

  std::cout<<"Generating tetrahedral mesh...\n";
  //Generate the mesh
  mesh.Generate();

  //Translate mesh into DEM particles
  //double R = 0.001;          // spheroradius
  //double rho = 1600.0;      // density of material
  //std::cout<<"Tag of the first mesh cell: "<<mesh.Cells[0]->Tag<<std::endl;
  std::cout<<"Generating particles from mesh...\n";
  dom.GenFromMesh (mesh,/*R*/R,/*rho*/rho,true/*Cohesion*/,/*MC*/false);
  //std::cout<<"Tag of the first mesh gen particle: "<<dom.Particles[IIndex]->Tag<<std::endl;
}

#endif
