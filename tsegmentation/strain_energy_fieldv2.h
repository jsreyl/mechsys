#ifndef STRAIN_ENERGY_FIELD
#define STRAIN_ENERGY_FIELD

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

template<class Interactons_T>
Vec3_t CalculateForce(Interactons_T Interactons, size_t pID);
//Vec3_t CalculateForce(Array<DEM::Interacton*> Interactons, size_t pID);
//Vec3_t CalculateForce(Array<DEM::BInteracton*> Interactons, size_t pID);
size_t CalculateContacts(Array<DEM::Interacton*> Interactons, size_t pID);
size_t CalculateContacts(Array<DEM::BInteracton*> Interactons, size_t pID);
size_t CalculateContacts(Array<DEM::Interacton*> Interactons, Array<DEM::BInteracton*> BInteractons, size_t pID);
Array <Mat3_t> StressTensor(Array<DEM::Interacton*> Interactons, size_t n_p);//Calculate all Stress Tensors for the domains interactions
Array <Mat3_t> StressTensor(Array<DEM::BInteracton*> Interactons, size_t n_p);//Calculate all Stress Tensors for the domain's cohesive interactions
Array <Mat3_t> StressTensor(Array<DEM::Interacton*> Interactons, Array<DEM::BInteracton*> BInteractons, size_t n_p);//Calculate all Stress Tensors for the domain's cohesive interactions
double StrainEnergyField(Mat3_t S, double nu, bool verbose = false);//Calculate the strain energy field given a Stress Tensor matrix

template<class Interactons_T>
inline Vec3_t CalculateForce(Interactons_T Interactons, size_t pID) {
  Vec3_t F = Vec3_t(0., 0., 0.);
  for(size_t i=0; i<Interactons.Size();i++){
	  if (Interactons[i]->I1 == pID ) F +=Interactons[i]->F1;
	  if (Interactons[i]->I2 == pID ) F +=Interactons[i]->F2;
  }
  return F;
}

inline size_t CalculateContacts(Array<DEM::Interacton*> Interactons, size_t pID) {
  size_t Nc = 0;
  for(size_t i=0; i<Interactons.Size();i++) if (Interactons[i]->I1 == pID || Interactons[i]->I2 == pID) Nc++;
  return Nc;
}
inline size_t CalculateContacts(Array<DEM::BInteracton*> Interactons, size_t pID) {
  size_t Nc = 0;
  for(size_t i=0; i<Interactons.Size();i++) if (Interactons[i]->I1 == pID || Interactons[i]->I2 == pID) Nc++;
  return Nc;
}
inline size_t CalculateContacts(Array<DEM::Interacton*> Interactons, Array<DEM::BInteracton*> BInteractons, size_t pID) {
  size_t Nc = 0;
  for(size_t i=0; i<Interactons.Size();i++) if (Interactons[i]->I1 == pID || Interactons[i]->I2 == pID) Nc++;
  for(size_t i=0; i<BInteractons.Size();i++) if (BInteractons[i]->I1 == pID || BInteractons[i]->I2 == pID) Nc++;
  return Nc;
}

inline Array<Mat3_t> StressTensor(Array<DEM::Interacton*> Interactons, size_t n_p){
  /*
   *Calculates the stress tensor for all particles in the domain
   INPUTS:
     - I: Interactons of the domain
     - n_p: number of particles in the domain
   */
  Array<Mat3_t> S(n_p);
  for(size_t i=0; i<n_p;i++) set_to_zero(S[i]); //Set everything to zero
  //Now calculate the components
  size_t n_interacton = Interactons.Size();
  for(size_t i=0; i<n_interacton;i++){
    DEM::Interacton * I = Interactons[i];
    Vec3_t r1 = (I->P1->x - I->Xc)/(I->P1->Props.V); //Position between the two particles
    //Particle i1, sum the dyadic product r1 x F
    S[I->I1](0,0) += r1(0)*I->F1(0); S[I->I1](0,1) += r1(0)*I->F1(1); S[I->I1](0,2) += r1(0)*I->F1(2);
    S[I->I1](1,0) += r1(1)*I->F1(0); S[I->I1](1,1) += r1(1)*I->F1(1); S[I->I1](1,2) += r1(1)*I->F1(2);
    S[I->I1](2,0) += r1(2)*I->F1(0); S[I->I1](2,1) += r1(2)*I->F1(1); S[I->I1](2,2) += r1(2)*I->F1(2);
    //Particle i2, sum the dyadic product F x r2
    Vec3_t r2 = (I->P2->x - I->Xc)/(I->P2->Props.V); //Position between the two particles
    //Particle i2, sum the dyadic product r2 x F
    S[I->I2](0,0) += r2(0)*I->F2(0); S[I->I2](0,1) += r2(0)*I->F2(1); S[I->I2](0,2) += r2(0)*I->F2(2);
    S[I->I2](1,0) += r2(1)*I->F2(0); S[I->I2](1,1) += r2(1)*I->F2(1); S[I->I2](1,2) += r2(1)*I->F2(2);
    S[I->I2](2,0) += r2(2)*I->F2(0); S[I->I2](2,1) += r2(2)*I->F2(1); S[I->I2](2,2) += r2(2)*I->F2(2);
  }
  // Make ridiculous things like 3.2e-311 into zeros so we can use a numeric solver on the stress tensors
  for(size_t p=0; p<n_p;p++)
    for(size_t i=0;i<3;i++)
      for(size_t j=0;j<3;j++){
        double s = S[p](i,j);
        S[p](i,j) = (fabs(s)>1e-4) ? s : 0.;
      }
  return S;
}

inline Array<Mat3_t> StressTensor(Array<DEM::BInteracton*> Interactons, size_t n_p){
  /*
   *Calculates the stress tensor for all particles in the domain
   INPUTS:
     - I: Interactons of the domain
     - n_p: number of particles in the domain
   */
  Array<Mat3_t> S(n_p);
  for(size_t i=0; i<n_p;i++) set_to_zero(S[i]); //Set everything to zero
  //Now calculate the components
  size_t n_interacton = Interactons.Size();
  for(size_t i=0; i<n_interacton;i++){
    DEM::BInteracton * I = Interactons[i];
    Vec3_t r1 = (I->P1->x - I->Xc)/(I->P1->Props.V); //Position between the two particles
    Vec3_t r2 = (I->P2->x - I->Xc)/(I->P2->Props.V); //Position between the two particles
    //Particle i1, sum the dyadic product r1 x F
    S[I->I1](0,0) += r1(0)*I->F1(0); S[I->I1](0,1) += r1(0)*I->F1(1); S[I->I1](0,2) += r1(0)*I->F1(2);
    S[I->I1](1,0) += r1(1)*I->F1(0); S[I->I1](1,1) += r1(1)*I->F1(1); S[I->I1](1,2) += r1(1)*I->F1(2);
    S[I->I1](2,0) += r1(2)*I->F1(0); S[I->I1](2,1) += r1(2)*I->F1(1); S[I->I1](2,2) += r1(2)*I->F1(2);
    //Particle i2, sum the dyadic product r2 x F
    S[I->I2](0,0) += r2(0)*I->F2(0); S[I->I2](0,1) += r2(0)*I->F2(1); S[I->I2](0,2) += r2(0)*I->F2(2);
    S[I->I2](1,0) += r2(1)*I->F2(0); S[I->I2](1,1) += r2(1)*I->F2(1); S[I->I2](1,2) += r2(1)*I->F2(2);
    S[I->I2](2,0) += r2(2)*I->F2(0); S[I->I2](2,1) += r2(2)*I->F2(1); S[I->I2](2,2) += r2(2)*I->F2(2);
  }
  // Make ridiculous things like 3.2e-311 into zeros so we can use a numeric solver on the stress tensors
  for(size_t p=0; p<n_p;p++)
    for(size_t i=0;i<3;i++)
      for(size_t j=0;j<3;j++){
        double s = S[p](i,j);
        S[p](i,j) = (fabs(s)>1e-4) ? s : 0.;
      }
  return S;
}

inline Array<Mat3_t> StressTensor(Array<DEM::Interacton*> Interactons,Array<DEM::BInteracton*> BInteractons, size_t n_p){
  /*
   *Calculates the stress tensor for all particles in the domain
   INPUTS:
     - I: Interactons of the domain
     - n_p: number of particles in the domain
   */
  Array<Mat3_t> S(n_p);
  for(size_t i=0; i<n_p;i++) set_to_zero(S[i]); //Set everything to zero
  //Now calculate the components
  size_t n_interacton = Interactons.Size();
  for(size_t i=0; i<n_interacton;i++){
    DEM::Interacton * I = Interactons[i];
    Vec3_t r1 = (I->P1->x - I->Xc)/(I->P1->Props.V); //Position between the two particles
    //Particle i1, sum the dyadic product r1 x F
    S[I->I1](0,0) += r1(0)*I->F1(0); S[I->I1](0,1) += r1(0)*I->F1(1); S[I->I1](0,2) += r1(0)*I->F1(2);
    S[I->I1](1,0) += r1(1)*I->F1(0); S[I->I1](1,1) += r1(1)*I->F1(1); S[I->I1](1,2) += r1(1)*I->F1(2);
    S[I->I1](2,0) += r1(2)*I->F1(0); S[I->I1](2,1) += r1(2)*I->F1(1); S[I->I1](2,2) += r1(2)*I->F1(2);
    //Particle i2, sum the dyadic product r2 x F
    Vec3_t r2 = (I->P2->x - I->Xc)/(I->P2->Props.V); //Position between the two particles
    S[I->I2](0,0) += r2(0)*I->F2(0); S[I->I2](0,1) += r2(0)*I->F2(1); S[I->I2](0,2) += r2(0)*I->F2(2);
    S[I->I2](1,0) += r2(1)*I->F2(0); S[I->I2](1,1) += r2(1)*I->F2(1); S[I->I2](1,2) += r2(1)*I->F2(2);
    S[I->I2](2,0) += r2(2)*I->F2(0); S[I->I2](2,1) += r2(2)*I->F2(1); S[I->I2](2,2) += r2(2)*I->F2(2);
  }
  n_interacton = BInteractons.Size();
  for(size_t i=0; i<n_interacton;i++){
    DEM::BInteracton * I = BInteractons[i];
    Vec3_t r1 = (I->P1->x - I->Xc)/(I->P1->Props.V); //Position between the two particles
    //Particle i1, sum the dyadic product r1 x F
    S[I->I1](0,0) += r1(0)*I->F1(0); S[I->I1](0,1) += r1(0)*I->F1(1); S[I->I1](0,2) += r1(0)*I->F1(2);
    S[I->I1](1,0) += r1(1)*I->F1(0); S[I->I1](1,1) += r1(1)*I->F1(1); S[I->I1](1,2) += r1(1)*I->F1(2);
    S[I->I1](2,0) += r1(2)*I->F1(0); S[I->I1](2,1) += r1(2)*I->F1(1); S[I->I1](2,2) += r1(2)*I->F1(2);
    //Particle i2, sum the dyadic product r2 x F
    Vec3_t r2 = (I->P2->x - I->Xc)/(I->P2->Props.V); //Position between the two particles
    S[I->I2](0,0) += r2(0)*I->F2(0); S[I->I2](0,1) += r2(0)*I->F2(1); S[I->I2](0,2) += r2(0)*I->F2(2);
    S[I->I2](1,0) += r2(1)*I->F2(0); S[I->I2](1,1) += r2(1)*I->F2(1); S[I->I2](1,2) += r2(1)*I->F2(2);
    S[I->I2](2,0) += r2(2)*I->F2(0); S[I->I2](2,1) += r2(2)*I->F2(1); S[I->I2](2,2) += r2(2)*I->F2(2);
  }
  // Make ridiculous things like 3.2e-311 into zeros so we can use a numeric solver on the stress tensors
  for(size_t p=0; p<n_p;p++)
    for(size_t i=0;i<3;i++)
      for(size_t j=0;j<3;j++){
        double s = S[p](i,j);
        S[p](i,j) = (fabs(s)>1e-4) ? s : 0.;
      }
  return S;
}

inline double StrainEnergyField(Mat3_t S, double nu, bool verbose){
  /*
   *Calculates the strain energy field according to Jiang, Y., Alonso-Marroquín, F., Herrmann, H.J. et al. Particle fragmentation based on strain energy field. Granular Matter 22, 69 (2020). https://doi.org/10.1007/s10035-020-01038-6
   INPUTS:
     - S : stress tensor of the particle
     - nu: Poisson's ratio of the material
   OUTPUTS:
     - strainEF: a scalar representing the Strain Energy Field
   */
  //Get the eigenvalues of the stress tensor
  Vec3_t eigvalR = Vec3_t(0.,0.,0.), eigvalI = Vec3_t(0.,0.,0.);//Since the stress tensor is nonsymmetric in general it's eigenvalues are complex numbers in general
  if (verbose) if(fabs(S(1,0)-S(0,1))<1e-2 && fabs(S(2,0)-S(0,2))<1e-2 && fabs(S(2,1)-S(1,2))<1e-2) std::cout<<"Stress tensor is symmetric!\n"<<S<<std::endl;
  try {
    EigNonsymm(S, eigvalR, eigvalI); //NOTE: This function destroys the original matrix S since it uses the Schur form transformation to get the eigenvalues
    // Also NOTE the eigenvalues are unordered, but this does not matter for hydrostatic and von Misses stress calculation
  } catch (...){
    eigvalR = Vec3_t(0.,0.,0.), eigvalI = Vec3_t(0.,0.,0.);
  }
  Vec3_t s = eigvalR;
  // Calculate the hydrostatic and Von-Mises stress scalars
  double stress_hydro2 = std::pow(s(0)+s(1)+s(2),2)/9.; //Hydrostatic stress squared
  double stress_vm2 = (std::pow(s(0)-s(1),2) + \
                       std::pow(s(1)-s(2),2) + \
                       std::pow(s(2)-s(0),2))/2.; // Von-Mises stress squared
  /* std::cout<<"Hydrostatic stress^2: "<<stress_hydro2<<"\n"; */
  /* std::cout<<"Von-Mises stress^2: "<<stress_vm2<<"\n"; */
  // Calculate the Strain Energy Field
  double strainEF = std::sqrt(2.*(1.+nu)*stress_vm2/3. + \
                              3.*(1.-2.*nu)*stress_hydro2);
  return strainEF;
}

#endif
