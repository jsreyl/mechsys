#ifndef USER_REPORT
#define USER_REPORT

#include <mechsys/dem/domain.h>
#include <boost/lexical_cast.hpp>

void ReportParticles(DEM::Domain Dom);
void WriteXDMF_SEF (DEM::Domain d, Array<double> strainEF, char const * FileKey);
void VisualizeParticles(DEM::Domain Dom, char const * FileKey, Array<DEM::Particle *> Parts);

inline void VisualizeParticles(DEM::Domain _vdom, char const * FileKey, Array<DEM::Particle *>Parts, Array<double> stress_vector){

  for(size_t p =0; p<Parts.Size();p++){
    DEM::Particle * P = Parts[p];
    Array<Vec3_t> Verts(0);
    for(size_t v=0; v<P->Verts.Size();v++){
      Verts.Push(*P->Verts[v]);
    }
    vdom.Particles.Push(new DEM::Particle(/*Tag*/P->Tag,
                                          /*vertices*/Verts,
                                          /*edges*/P->EdgeCon,
                                          /*faces*/P->FaceCon,
                                          /*init vel*/P->v,
                                          /*init omega*/P->w,
                                          /*Radius*/P->Props.R,
                                          /*density*/P->Props.rho));
  }
  std::cout<<"Writing visualization files.\n";
  _vdom.WriteXDMF_User(stress_vector, FileKey);
}

inline void ReportParticles(DEM::Domain Dom){
  // Create a file to store particle positions on each Report step
  std::ofstream particleFile;
  String filename;
  filename = "particles_"+std::to_string(Dom.Time)+".txt";
  particleFile.open(filename);
  particleFile << Util::_10_6 << "Time" << Util::_10_6 << "ID" << Util::_10_6 << "Tag" << Util::_8s << "x" << Util::_8s << "y" << Util::_8s << "z" << Util::_8s << "vx" << Util::_8s << "vy" << Util::_8s << "vz" << Util::_8s << "wx" << Util::_8s << "wy" << Util::_8s << "wz" << Util::_8s << "Fx" << Util::_8s << "Fy" << Util::_8s << "Fz" << Util::_8s << "Tx" << Util::_8s << "Ty" << Util::_8s << "Tz" << std::endl;
  for(size_t i=0; i<Dom.Particles.Size(); i++){
    particleFile << Util::_10_6 << Dom.Time << Util::_10_6 << i << Util::_10_6 << Dom.Particles[i]->Tag << Util::_8s << Dom.Particles[i]->x(0) << Util::_8s << Dom.Particles[i]->x(1) << Util::_8s << Dom.Particles[i]->x(2) << Util::_8s << Dom.Particles[i]->v[0] << Util::_8s << Dom.Particles[i]->v[1] << Util::_8s << Dom.Particles[i]->v[2] << Util::_8s << Dom.Particles[i]->w[0] << Util::_8s << Dom.Particles[i]->w[1] << Util::_8s << Dom.Particles[i]->w[2] << Util::_8s << Dom.Particles[i]->F(0) << Util::_8s << Dom.Particles[i]->F(1) << Util::_8s << Dom.Particles[i]->F(2) << Util::_8s << Dom.Particles[i]->T(0) << Util::_8s << Dom.Particles[i]->T(1) << Util::_8s << Dom.Particles[i]->T(2) << std::endl;
  }
  particleFile.close();
  // Also create a file to store interactions on each Report step
  std::ofstream interactionsFile;
  filename = "interactions_"+std::to_string(Dom.Time)+".txt";
  interactionsFile.open(filename);
  interactionsFile << Util::_10_6 << "Time" << Util::_10_6 << "Interacton_ID" << Util::_10_6 << "P1_ID" << Util::_10_6 << "P2_ID" << Util::_10_6 << "P1_Tag" << Util::_10_6 << "P2_Tag" << Util::_8s << "F1x" << Util::_8s << "F1y" << Util::_8s << "F1z" << Util::_8s << "F2x" << Util::_8s << "F2y" << Util::_8s << "F2z" << Util::_8s << "T1x" << Util::_8s << "T1y" << Util::_8s << "T1z" << Util::_8s << "T2x" << Util::_8s << "T2y" << Util::_8s << "T2z" << std::endl;
  for(size_t i=0; i<Dom.Interactons.Size(); i++){
    interactionsFile << Util::_10_6 << Dom.Time << Util::_10_6 << i << Util::_10_6 << Dom.Interactons[i]->I1 << Util::_10_6 << Dom.Interactons[i]->I2 << Util::_10_6 << Dom.Interactons[i]->P1->Tag << Util::_10_6 << Dom.Interactons[i]->P2->Tag << Util::_8s << Dom.Interactons[i]->F1(0) << Util::_8s << Dom.Interactons[i]->F1(1) << Util::_8s << Dom.Interactons[i]->F1(2) << Util::_8s << Dom.Interactons[i]->F2(0) << Util::_8s << Dom.Interactons[i]->F2(1) << Util::_8s << Dom.Interactons[i]->F2(2) << Util::_8s << Dom.Interactons[i]->T1(0) << Util::_8s << Dom.Interactons[i]->T1(1) << Util::_8s << Dom.Interactons[i]->T1(2) << Util::_8s << Dom.Interactons[i]->T2(0) << Util::_8s << Dom.Interactons[i]->T2(1) << Util::_8s << Dom.Interactons[i]->T2(2) << std::endl;
  }
  interactionsFile.close();
}

inline void WriteXDMF_SEF (DEM::Domain d, Array<double> strainEF, char const * FileKey)
{
  /*
   *Writes .h5 files saving the domain propierties in binary format and .xmf files for visualization using VisIt
   INPUTS:
   * d : DEM::Domain containing the particles and all their information
   * strainEF: Array containing the strain energy field scalar calculated for all particles in the domain. Must be of size Particles.Size()
   * FileKey: name for the .h5 and .xmf files
   */
    size_t N_Faces = 0;
    size_t N_Verts = 0;
    for (size_t i=0; i<d.Particles.Size(); i++)
    {
        for (size_t j=0;j<d.Particles[i]->Faces.Size();j++)
        {
            N_Faces += d.Particles[i]->Faces[j]->Edges.Size();
        }
        N_Verts += d.Particles[i]->Verts.Size() + d.Particles[i]->Faces.Size();
    }

    String fn(FileKey);
    fn.append(".h5");
    hid_t     file_id;
    file_id = H5Fcreate(fn.CStr(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    if (N_Faces>0)
    {

        //Geometric information
        float  * Verts   = new float [3*N_Verts];
        int    * FaceCon = new int   [3*N_Faces];
        
        //Atributes
        int    * Tags    = new int   [  N_Faces];
        int    * Clus    = new int   [  N_Faces];
        float  * Vel     = new float [  N_Faces];
        float  * Ome     = new float [  N_Faces];
        float  * sEF     = new float [  N_Faces]; //Paint the sEF on the faces, not the centres
        //float  * Stress  = new float [9*N_Faces];

        size_t n_verts = 0;
        size_t n_faces = 0;
        size_t n_attrs = 0;
        //size_t n_attrv = 0;
        //size_t n_attrt = 0;
        for (size_t i=0;i<d.Particles.Size();i++)
        {
            DEM::Particle * Pa = d.Particles[i];
            size_t n_refv = n_verts/3;
            Array<Vec3_t> Vtemp(Pa->Verts.Size());
            Array<Vec3_t> Vres (Pa->Verts.Size());
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                Vtemp[j] = *Pa->Verts[j];
                Vres [j] = *Pa->Verts[j];
            }
            double multiplier = 0.0;
            if (d.Dilate&&Pa->Eroded&&Pa->Faces.Size()>=4)
            {
                DEM::Dilation(Vtemp,Pa->EdgeCon,Pa->FaceCon,Vres,Pa->Props.R);
                multiplier = 1.0;
            }
            for (size_t j=0;j<Pa->Verts.Size();j++)
            {
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(0);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(1);
                //Verts[n_verts++] = (float) (*Pa->Verts[j])(2);
                Verts[n_verts++] = float(Vres[j](0));
                Verts[n_verts++] = float(Vres[j](1));
                Verts[n_verts++] = float(Vres[j](2));
            }
            size_t n_reff = n_verts/3;
            for (size_t j=0;j<Pa->FaceCon.Size();j++)
            {
                Vec3_t C,N;
                Pa->Faces[j]->Centroid(C);
                Pa->Faces[j]->Normal(N);
                Verts[n_verts++] = float(C(0) + multiplier*Pa->Props.R*N(0));
                Verts[n_verts++] = float(C(1) + multiplier*Pa->Props.R*N(1));
                Verts[n_verts++] = float(C(2) + multiplier*Pa->Props.R*N(2));
                //Verts[n_verts++] = (float) C(0);
                //Verts[n_verts++] = (float) C(1);
                //Verts[n_verts++] = (float) C(2);
                for (size_t k=0;k<Pa->FaceCon[j].Size();k++)
                {
                    size_t nin = Pa->FaceCon[j][k];
                    size_t nen = Pa->FaceCon[j][(k+1)%Pa->FaceCon[j].Size()];
                    FaceCon[n_faces++] = int(n_reff + j);
                    FaceCon[n_faces++] = int(n_refv + nin);
                    FaceCon[n_faces++] = int(n_refv + nen);

                    //Writing the attributes
                    Tags  [n_attrs] = int(Pa->Tag);
                    Clus  [n_attrs] = size_t(Pa->Cluster);
                    Vel   [n_attrs] = float(norm(Pa->v));
                    Ome   [n_attrs] = float(norm(Pa->w));
                    sEF   [n_attrs] = float(strainEF[i]);
                    n_attrs++;

                    //Vel [n_attrv  ] = (float) Pa->v(0);
                    //Vel [n_attrv+1] = (float) Pa->v(1);
                    //Vel [n_attrv+2] = (float) Pa->v(2);
                    //n_attrv += 3;

                    //Stress[n_attrt  ] = (float) Pa->M(0,0);
                    //Stress[n_attrt+1] = (float) Pa->M(1,0);
                    //Stress[n_attrt+2] = (float) Pa->M(2,0);
                    //Stress[n_attrt+3] = (float) Pa->M(0,1);
                    //Stress[n_attrt+4] = (float) Pa->M(1,1);
                    //Stress[n_attrt+5] = (float) Pa->M(2,1);
                    //Stress[n_attrt+6] = (float) Pa->M(0,2);
                    //Stress[n_attrt+7] = (float) Pa->M(1,2);
                    //Stress[n_attrt+8] = (float) Pa->M(2,2);
                    //n_attrt += 9;
                }
            }
            std::cout<<"No problemo saving sEF for particle "<<i<<"\n";
        }

        //Write the data
        hsize_t dims[1];
        String dsname;
        dims[0] = 3*N_Verts;
        dsname.Printf("Verts");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Verts);
        dims[0] = 3*N_Faces;
        dsname.Printf("FaceCon");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,FaceCon);
        dims[0] = N_Faces;
        dsname.Printf("Tag");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Tags   );
        dims[0] = N_Faces;
        dsname.Printf("Cluster");
        H5LTmake_dataset_int(file_id,dsname.CStr(),1,dims,Clus   );
        dims[0] = N_Faces;
        dsname.Printf("Velocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Vel);
        dims[0] = N_Faces;
        dsname.Printf("AngVelocity");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ome);
        dims[0] = N_Faces;
        dsname.Printf("StrainEnergyField");
        H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,sEF);
        
        //Erasing the data
        delete [] Verts;
        delete [] FaceCon;
        delete [] Tags;
        delete [] Clus;
        delete [] Vel;
        delete [] Ome;
        delete [] sEF;
        //delete [] Stress;
    }
    // Storing center of mass data
    
    float * Radius = new float[  d.Particles.Size()];
    float * Posvec = new float[3*d.Particles.Size()];
    float * Velvec = new float[3*d.Particles.Size()];
    float * Omevec = new float[3*d.Particles.Size()];
    float * Aomvec = new float[3*d.Particles.Size()];
    float * Amovec = new float[3*d.Particles.Size()];
    float * Ufovec = new float[3*d.Particles.Size()];
    float * Inertm = new float[6*d.Particles.Size()];
    float * Ekin   = new float[  d.Particles.Size()];
    int   * Tag    = new int  [  d.Particles.Size()];

    for (size_t i=0;i<d.Particles.Size();i++)
    {
        Mat3_t Inertia,Inertiar,R,Rt,t;
        Inertia(0,0) = d.Particles[i]->I(0); Inertia(0,1) = 0.0; Inertia(0,2) = 0.0;
        Inertia(1,1) = d.Particles[i]->I(1); Inertia(1,0) = 0.0; Inertia(1,2) = 0.0;
        Inertia(2,2) = d.Particles[i]->I(2); Inertia(2,0) = 0.0; Inertia(2,1) = 0.0;

        RotationMatrix(d.Particles[i]->Q,R);
        Rt = ~R;

        Mult(R,Inertia,t);
        Mult(t,Rt,Inertiar);

        Vec3_t Ao,Ome,L,t1,t2;
        Rotation(d.Particles[i]->w,d.Particles[i]->Q,Ome);
        Rotation(d.Particles[i]->wa,d.Particles[i]->Q,Ao);
        t1 = d.Particles[i]->I(0)*d.Particles[i]->w(0),d.Particles[i]->I(1)*d.Particles[i]->w(1),d.Particles[i]->I(2)*d.Particles[i]->w(2);
        Rotation (t1,d.Particles[i]->Q,t2);
        L = d.Particles[i]->Props.m*cross(d.Particles[i]->x,d.Particles[i]->v)+t2;


        d.Particles[i]->Verts.Size()==1 ? Radius[i] = float(d.Particles[i]->Dmax) : Radius[i] = 0.0;
        Posvec[3*i  ] = float(d.Particles[i]->x(0));
        Posvec[3*i+1] = float(d.Particles[i]->x(1));
        Posvec[3*i+2] = float(d.Particles[i]->x(2));
        Velvec[3*i  ] = float(d.Particles[i]->v(0));
        Velvec[3*i+1] = float(d.Particles[i]->v(1));
        Velvec[3*i+2] = float(d.Particles[i]->v(2));
        Omevec[3*i  ] = float(Ome(0));
        Omevec[3*i+1] = float(Ome(1)); 
        Omevec[3*i+2] = float(Ome(2)); 
        Aomvec[3*i  ] = float(Ao(0));
        Aomvec[3*i+1] = float(Ao(1)); 
        Aomvec[3*i+2] = float(Ao(2)); 
        Amovec[3*i  ] = float(L(0));
        Amovec[3*i+1] = float(L(1)); 
        Amovec[3*i+2] = float(L(2)); 
        Ufovec[3*i  ] = (float) d.Particles[i]->F(0);
        Ufovec[3*i+1] = (float) d.Particles[i]->F(1); 
        Ufovec[3*i+2] = (float) d.Particles[i]->F(2); 
        Inertm[6*i  ] = (float) Inertiar(0,0);
        Inertm[6*i+1] = (float) Inertiar(0,1);
        Inertm[6*i+2] = (float) Inertiar(0,2);
        Inertm[6*i+3] = (float) Inertiar(1,1);
        Inertm[6*i+4] = (float) Inertiar(1,2);
        Inertm[6*i+5] = (float) Inertiar(2,2);
        Ekin  [i]     = float(d.Particles[i]->Ekin+d.Particles[i]->Erot);
        Tag   [i]     = int  (d.Particles[i]->Tag);  
    }

    hsize_t dims[1];
    dims[0] = 6*d.Particles.Size();
    String dsname;
    dsname.Printf("Inertia");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Inertm);
    dims[0] = 3*d.Particles.Size();
    dsname.Printf("Position");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Posvec);
    dsname.Printf("PVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Velvec);
    dsname.Printf("PAngVelocity");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Omevec);
    dsname.Printf("PAngacceleration");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Aomvec);
    dsname.Printf("PAngMomentum");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Amovec);
    dsname.Printf("PUForce");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ufovec);
    dims[0] = d.Particles.Size();
    dsname.Printf("Radius");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Radius);
    dsname.Printf("PEkin");
    H5LTmake_dataset_float(file_id,dsname.CStr(),1,dims,Ekin);
    dsname.Printf("PTag");
    H5LTmake_dataset_int  (file_id,dsname.CStr(),1,dims,Tag   );


    delete [] Radius;
    delete [] Posvec;
    delete [] Velvec;
    delete [] Omevec;
    delete [] Aomvec;
    delete [] Amovec;
    delete [] Ufovec;
    delete [] Inertm;
    delete [] Ekin;
    delete [] Tag;


    //Closing the file
    H5Fflush(file_id,H5F_SCOPE_GLOBAL);
    H5Fclose(file_id);
    

    //Writing xmf file
    std::ostringstream oss;
    oss << "<?xml version=\"1.0\" ?>\n";
    oss << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
    oss << "<Xdmf Version=\"2.0\">\n";
    oss << " <Domain>\n";
    if(N_Faces>0)
    {
    oss << "   <Grid Name=\"DEM_Faces\">\n";
    oss << "     <Topology TopologyType=\"Triangle\" NumberOfElements=\"" << N_Faces << "\">\n";
    oss << "       <DataItem Format=\"HDF\" DataType=\"Int\" Dimensions=\"" << N_Faces << " 3\">\n";
    oss << "        " << fn.CStr() <<":/FaceCon \n";
    oss << "       </DataItem>\n";
    oss << "     </Topology>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << N_Verts << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Verts \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Tag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Cluster\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Cluster \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Velocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVelocity\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/AngVelocity \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"StrainEnergyField\" AttributeType=\"Scalar\" Center=\"Cell\">\n";
    oss << "       <DataItem Dimensions=\"" << N_Faces << "\" NumberType=\"Float\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/StrainEnergyField \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "   </Grid>\n";
    }
    oss << "   <Grid Name=\"DEM_Center\" GridType=\"Uniform\">\n";
    oss << "     <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << d.Particles.Size() << "\"/>\n";
    oss << "     <Geometry GeometryType=\"XYZ\">\n";
    oss << "       <DataItem Format=\"HDF\" NumberType=\"Float\" Precision=\"4\" Dimensions=\"" << d.Particles.Size() << " 3\" >\n";
    oss << "        " << fn.CStr() <<":/Position \n";
    oss << "       </DataItem>\n";
    oss << "     </Geometry>\n";
    oss << "     <Attribute Name=\"Radius\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Radius \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Ekin\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << "\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PEkin \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Tag\" AttributeType=\"Scalar\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << "\" NumberType=\"Int\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PTag \n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngVel\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngVelocity\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Angacc\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngacceleration\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"AngMom\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PAngMomentum\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Uforce\" AttributeType=\"Vector\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << " 3\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/PUForce\n";
    oss << "       </DataItem>\n";
    oss << "     </Attribute>\n";
    oss << "     <Attribute Name=\"Inertia\" AttributeType=\"Tensor6\" Center=\"Node\">\n";
    oss << "       <DataItem Dimensions=\"" << d.Particles.Size() << " 6\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n";
    oss << "        " << fn.CStr() <<":/Inertia\n";
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
