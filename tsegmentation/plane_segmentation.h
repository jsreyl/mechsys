#ifndef PLANE_SEGMENTATION
#define PLANE_SEGMENTATION

//plane_segmentationv7.h in the old directory

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

// User
//#include <algorithm> //std::find
//#include <iterator> //std::begin std::end std::distance
#include<assert.h>

using std::cout;
using std::endl;

void LocalImposeParticleCohesion(int Tag, DEM::Domain & Dom, double tol1=1.0e-8, double tol2=1.0e-3, double L0=0.);
Array<size_t> GetLargest3(Array<double> Fs);
size_t CalculateIntersection(Vec3_t planeNormal, Vec3_t planePoint, Vec3_t X0, Vec3_t X1, double dist0, double dist1, Vec3_t & intersectionVertex, double tolerance = 1e-5);
bool FindVectorInList(Vec3_t vertex, Array<Vec3_t> vertList, size_t & position, double tolerance = 1e-5);
bool FindIndexInConList(int index, Array<int> conList, size_t & position);
Array< Array <int> > ReorderEdges(Array<Array<int>> edgeConList);
bool FindEdgePairInList(Array<Array<int>> edgeConList, Array<int>edge);
Array<int> ReorderFaceCon(Array<Vec3_t> & Verts, Array<Array<int>> EdgeCon, Array<int> preFaceCon, bool debug=false);
void AddVertexToConLists(Vec3_t vertex, Array<Vec3_t> & vertsList, Array<int> & edgeCon, Array<int> & faceCon);
void IntrudeFace(Array<Vec3_t> & Verts, Array<int> FaceCon, double R);
void SimplifyGeometry(Array<int>& intersectionsIxs, Array<Vec3_t>& Verts, Array<Array<int>>&EdgeCon,Array<Array<int>>& FaceCon, double threshold=0.1, bool verbose=false);
void SimplifyGeometryNew(Array<int>& intersectionsIxs, Array<Vec3_t>& Verts, Array<Array<int>>&EdgeCon,Array<Array<int>>& FaceCon, double threshold=0.1, bool verbose=false);
void SimplifyGeometry_avg(Array<Vec3_t>& GlobalIntersections, Array<Vec3_t>& Verts, Array<Array<int>>&EdgeCon,Array<Array<int>>& FaceCon, double threshold=0.1, bool verbose=false);
inline void PrintGeometry(Array<Vec3_t> Verts,  Array<Array<int>> EdgeCon, Array<Array<int>> FaceCon, const char* optmsg = "", const char* optmsg2="");
void BisectGeometry(Array<Vec3_t> Verts, Array<Array<int>> EdgeCon, Array<Array<int>> FaceCon, Vec3_t planeNormal, Vec3_t planePoint, Array<Vec3_t>& oVertsPos, Array<Array<int>>&oEdgeConPos,Array<Array<int>>& oFaceConPos, Array<int> & intersectionIxsPos, Array<Vec3_t>& oVertsNeg, Array<Array<int>>&oEdgeConNeg,Array<Array<int>>& oFaceConNeg, Array<int> & intersectionIxsNeg, double tolerance = 1e-5, bool debug = false);
void BisectPolyhedron(DEM::Particle * ogParticle, Vec3_t planeNormal, Vec3_t planePoint, Array<DEM::Particle*> & segParticles, Array<Array<int>> & intersectionIxs, int mechsysErode = 1, double tolerance = 1e-5, bool debug = false);

/*
************************IMPLEMENTATION*******************************
 */
inline  void LocalImposeParticleCohesion(int Tag, DEM::Domain & Dom, double tol1, double tol2, double L0){
  /*
   *This function imposes a cohesive interaction between all particles in the given domain
   *Inputs:
   - Tag, of the particles
   - Dom, domain of the particles
   - tolerances for angle and distance separation between faces
   - L0, additional leverage distance for particle interaction
   *Outputs: None
   */
  Util::Stopwatch stopwatch; //Count how much time is spent generating the cohesions
  std::cout<<"Previous number of cohesion interactions: "<<Dom.BInteractons.Size()<<std::endl;
  std::cout<<"Previous number of interactons: "<<Dom.Interactons.Size()<<std::endl;
  //Get particles of the desired Tag
  Array<DEM::Particle *> p;
  Dom.GetParticles(Tag, p);
  if (p.Size()<2) { std::cout<<"Less than two particles of the desired tag... No cohesion can be added\n"; return;}
  for (size_t i=0;i<p.Size()-1;i++)//Loop over all particles
    {
      DEM::Particle * P1 = p[i];
      for (size_t j=i+1;j<p.Size();j++)
        {
          DEM::Particle * P2 = p[j];
          //double R =0.5*(P1->Props.R+P2->Props.R);
          //Debugging
          //std::cout<<"We're particles index "<<P1->Index<<" and "<<P2->Index<<std::endl;
          //std::cout<<"And our radiuses are "<<P1->Props.R<<" and "<<P2->Props.R<<std::endl;
          //std::cout<<"Therefore our mean radius is "<<R<<std::endl;
          //std::cout<<"And our Dmaxes are "<<P1->Dmax<<" and "<<P2->Dmax<<std::endl;
          //std::cout<<"And our Positions are "<<P1->x<<" and "<<P2->x<<std::endl;
          if (DEM::Distance(P1->x,P2->x)<P1->Dmax+P2->Dmax+L0)//If the particles are close together
            {
              double R =0.5*(P1->Props.R+P2->Props.R);
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
                          &&(fabs(Distance(c1,*F2)-2*R)<tol2+L0)
                          &&(fabs(Distance(c2,*F1)-2*R)<tol2+L0))//check for distance and angle treshholds
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
  std::cout<<"Total number of interactions: "<<Dom.Interactons.Size()<<std::endl;
}

inline  size_t CalculateIntersection(Vec3_t planeNormal, Vec3_t planePoint, Vec3_t X0, Vec3_t X1, double dist0, double dist1, Vec3_t & intersectionVertex, double tolerance){
  /*
   *This function calculates the intersection point between an edge X0-X1 and a plane defined by a normal vector, the intersection case is returned and if there is only one intersection point it is passed to the intertectionVertex vector
   *Inputs:
   - planeNormal, normal vector of the intersecting plane
   - X0, X1, vertices of an edge
   - intersectionVertex, reference vector to return intersection coordinates
   *Outputs:
   - intersectionState: 0 if no intersection,  1 if one intersection, 2 if two intersections (both X0 and X1 lie on the plane)
   */
  size_t intersectionState = 0; //No intersections until proven otherwise
  //dist0 = dot(planeNormal, X0-planePoint); //distance from the plane to point X0, notice this is signed so - means to the 'left' and + to the 'right'
  //dist1 = dot(planeNormal, X1-planePoint);//distance from the plane to point X1
  //if any of the distances is 0 it means that point is on the plane
  // for effective purposes there is always a computational error so use a low tolerance value
  if(abs(dist0) < tolerance){
    intersectionVertex = X0;
    //dist0 = 0;
    intersectionState++;
  }
  if(abs(dist1) < tolerance){
    intersectionVertex = X1;
    //dist1 = 0;
    intersectionState++;
  }
  if (dist0*dist1 < 0){//If one of the distances is positive and the other is negative, the vertices of the edge are at opposite sides of the plane and there is an intersection point between them
    Vec3_t edgeDirection = (X1-X0)/norm(X1-X0);
    double cosTheta = dot(edgeDirection, planeNormal);//Cosine of the angle between the plane normal and edge direction
    intersectionVertex = X1 - edgeDirection*dist1/cosTheta;
    intersectionState++;
  }
  // if (dist0*dist1>0) //both points are at the same side of the plane and no intersection is found
  return intersectionState;
}

inline bool FindVectorInList(Vec3_t vertex, Array<Vec3_t> vertList, size_t & position, double tolerance){
  bool foundVertex = false;
  //Find if vertex is in vertList, O(n) algorithm
  for(size_t i=0; i<vertList.Size(); i++){
    if(norm(vertList[i]-vertex) < tolerance) {//If the distance is 0 they are the same vector
      position = i;
      foundVertex = true;
      //std::cout<<"Found vector "<<vertex<<" on list "<<vertList<<" at position "<<position<<"\n";
      //std::cout<<vertList[position]<<"\n";
      break;
    }
  }
  return foundVertex;
}

inline bool FindIndexInConList(int index, Array<int> conList, size_t & position){
  bool foundIndex = false;
  //Find if vertex is in vertList, O(n) algorithm
  for(size_t i=0; i<conList.Size(); i++){
    if(conList[i] == index) {
      position = i;
      foundIndex = true;
      break;
    }
  }
  return foundIndex;
}

inline Array<Array<int>> ReorderEdges(Array<Array<int>> edgeConList){
  //First reorder into a list of edge pairs
  Array<Array<int>>edgePairList;
  for(size_t i=0;i<edgeConList.Size();i++){
    assert(edgeConList[i].Size()%2 == 0);//edges must be PAIRS of vertices, if the list is odd there's a vertex with no pair
    for(size_t j = 0; j<edgeConList[i].Size();j = j + 2){
      //std::cout<<"i: "<<i<<" j: "<<j<<" edge[i][j]: "<<edgeConList[i][j]<<"\n";
      Array<int> _temp(2);
      if(edgeConList[i][j]!=edgeConList[i][j+1]){
        _temp = edgeConList[i][j], edgeConList[i][j+1];
      }else{ //When a plane cuts exactly the vertices the edge connections might shifted by one ciclicaly to turn [0,0,1,1,2,2] into [0,1,1,2,2,0]
        _temp = edgeConList[i][j+1], edgeConList[i][(j+2)%edgeConList[i].Size()];
      }
      //std::cout<<"temp: "<<_temp<<"\n";
      //Push only unique edges
      if(!FindEdgePairInList(edgePairList, _temp))edgePairList.Push(_temp);
    }
  }
  return edgePairList;
}

inline bool FindEdgePairInList(Array<Array<int>> edgeConList, Array<int>edge){
  bool found = false;
  Array<int>edgeSwap(2);
  edgeSwap = edge[1], edge[0];
  for(size_t i=0; i<edgeConList.Size();i++ ){
    if(edgeConList[i] == edge or edgeConList[i] == edgeSwap) found = true;
  }
  return found;
}

inline Array<int> ReorderFaceCon(Array<Vec3_t> & Verts, Array<Array<int>> EdgeCon, Array<int> preFaceCon, bool debug){
  /*Reorders a FaceCon to ensure the order of the vertices is well connected and follows the right hand rule for cycling to vertices
   */
  size_t n = preFaceCon.Size();
  //First vertex is gonna be the same for the new array
  Array<int> newFaceCon;
  if(debug) std::cout<<"ReorderFaceCon from preFaceCon: "<<preFaceCon<<"\n";
  newFaceCon.Push(preFaceCon[0]);
  if(debug) std::cout<<"newFaceCon: "<<newFaceCon<<"\n";
  Array<int> searchFaceCon = preFaceCon;
  searchFaceCon.DelItem(0);//Delete item with index 0 as this one is already on the newFaceCon
  if(debug) std::cout<<"searchFaceCon: "<<searchFaceCon<<"\n";
  // size_t _ = -1;
  while(searchFaceCon.Size()>0){
    if(debug) std::cout<<"newFaceCon: "<<newFaceCon<<"\n";
    if(debug) std::cout<<"searchFaceCon: "<<searchFaceCon<<"\n";
    for(size_t i = 0; i<searchFaceCon.Size();i++){//Look for a connection with each element of the searh list until you find a valid one
      Array<int> temp(2);
      temp = newFaceCon.Last(), searchFaceCon[i];
      if(debug) std::cout<<"Trying pair "<<temp<<"\n";
      if(FindEdgePairInList(EdgeCon, temp)){
        if(debug) std::cout<<"Found connection, adding.\n";
        newFaceCon.Push(searchFaceCon[i]);
        searchFaceCon.DelItem(i);
        break;}
    }
  }
  if(debug) std::cout<<"newFaceCon: "<<newFaceCon<<"\n";
  if(debug) std::cout<<"searchFaceCon: "<<searchFaceCon<<"\n";
  // for(size_t i=0;i<n;i++){//Go through each of the n vertices we have to fill in the new list, we already added the 0th so n-1 to go
  //   for(size_t j=0;j<n;j++){//And search for a connected vertex in the previous list
  //     std::cout<<"i: "<<i<<" j: "<<j<<"\n";
  //     std::cout<<"preFaceCon: "<<preFaceCon<<"\n";
  //     std::cout<<"newFaceCon: "<<newFaceCon<<"\n";
  //     if(!FindIndexInConList(preFaceCon[j], newFaceCon, _)){//If the vertex not already in our new list
  //       std::cout<<"index "<<preFaceCon[j]<<"not found in newFaceCon, searching for a pair...\n";
  //       //see if pair ij is connected and add j if it is
  //       Array<int> temp(2);
  //       temp = preFaceCon[i], preFaceCon[j];
  //       std::cout<<"Trying pair "<<temp<<"\n";
  //       if(FindEdgePairInList(EdgeCon, temp)){
  //         std::cout<<"Found connection, adding.\n";
  //         newFaceCon.Push(preFaceCon[j]);
  //         break;}
  //     }
  //   }
  // }

  //Make sure the normal vector of intersection face that we just reordered points outside the polyhedron
  //Calculate the geometric centre of the polyhedron
  Vec3_t CP = Vec3_t(0.0,0.0,0.0);
  for(size_t v=0; v<Verts.Size();v++) CP += Verts[v];
  CP/=Verts.Size();
  //Calculate the centre of the intersection face
  Vec3_t CF = Vec3_t(0.0, 0.0, 0.0);
  for(size_t v=0; v<newFaceCon.Size();v++) CF += Verts[newFaceCon[v]];
  CF/=newFaceCon.Size();
  Vec3_t out = CF-CP; // And build a vector that points outside the polyhedron towards the face
  //And calculate the normal vector of the intersection face from two edges
  Vec3_t v1 = Verts[newFaceCon[0]]-Verts[newFaceCon[1]];
  Vec3_t v2 = Verts[newFaceCon[2]]-Verts[newFaceCon[1]];
  Vec3_t normal = cross(v2, v1);
  double direction = dot(out, normal);//If the dot product is negative the normal vector points inside the polyhedron and we need to rearrange the vertices on the other direction
  if(direction < 0){//Rearrange
    if(debug) std::cout<<"Rearranging from FaceCon: "<<newFaceCon<<"\n";
    Array<int> reFaceCon(n);//New aray with reordered faces
    for(size_t i=0; i<n;i++){
      if(debug) std::cout<<"newFacecon[n-i-1="<<n-i-1<<"] ="<<newFaceCon[n-i-1]<<"\n";
      reFaceCon[i] = newFaceCon[n-i-1];
    }
    // std::cout<<"Rearranged FaceCon: "<<reFaceCon<<"\n";
    newFaceCon = reFaceCon;
  }
  return newFaceCon;
}

inline void IntrudeFace(Array<Vec3_t> & Verts, Array<int> FaceCon, double R){
  /*Pushes face vertices towards the geometrical centre of a polyhedron, similar to eroding only a face
   * INPUTS:
   * - Verts: List of all the vertex vectors that form a polyhedron
   * - FaceCon: List of the indices of the vertices to be pushed (for instance the vertices that form a face)
   * - R: Distance to intrude the vertices
   * OUTPUTS:
   * - intVerts: List of vertex vectors after intrusion
   */
  //Calculate the geometrical centre
  Vec3_t CP = Vec3_t(0.,0.,0.);
  for(size_t v=0; v<Verts.Size();v++) CP += Verts[v];
  CP/=Verts.Size();
  //Calculate the centre of the intersection face
  Vec3_t CF = Vec3_t(0.0, 0.0, 0.0);
  for(size_t v=0; v<FaceCon.Size();v++) CF += Verts[FaceCon[v]];
  CF/=FaceCon.Size();
  Vec3_t in = CP-CF;
  double x = R/norm(in);//How much to intrude
  //Push verts in FaceCon towards the geometrical centre
  for(size_t i=0; i<FaceCon.Size();i++){
    Vec3_t dir = Verts[FaceCon[i]] - CP;
    Verts[FaceCon[i]]= Verts[FaceCon[i]] - x*dir;
  }
}

inline void AddVertexToConLists(Vec3_t vertex, Array<Vec3_t> & vertsList, Array<int> & edgeCon, Array<int> & faceCon){
  /*
   * Adds a vertex to a list given it's not there and adds its indices to connection lists accordingly
   * INPUTS:
   * - vertex: vector to add to the list if it's not there already
   * - vertsList: list of vectors to search and add vertices to
   * - edgeCon: list of edge connections to store vertex indices in the right order
   * - faceCon: list of face connections to store vertex indices in the right order
   */
  //Find if the vertex is already on the VertsList
  size_t _, v_ix = -1; //position index in case vertex is in the list
  if(!FindVectorInList(vertex, vertsList, v_ix)){//If it isn't, add it at the end of the list
    vertsList.Push(vertex);
    edgeCon.Push(vertsList.Size()-1);
    faceCon.Push(vertsList.Size()-1);
  } else{// If it already exists in the verts list don't add it but consider it in the edges and faces connections
    edgeCon.Push(v_ix);
    if(!FindIndexInConList(v_ix,faceCon, _)) faceCon.Push(v_ix);//If it isn't already in FaceCon add it aswell
  }
}

inline void PrintGeometry(Array<Vec3_t> Verts,  Array<Array<int>> EdgeCon, Array<Array<int>> FaceCon, const char* optmsg, const char* optmsg2){
  std::cout<<optmsg2;
  for(size_t v=0; v<Verts.Size();v++) std::cout<<optmsg<<" Vertex "<<v<<" : "<<Verts[v]<<std::endl;
  for(size_t e=0; e<EdgeCon.Size();e++) std::cout<<optmsg<<" Edge "<<e<<" : "<<EdgeCon[e]<<std::endl;
  for(size_t f=0; f<FaceCon.Size();f++) std::cout<<optmsg<<" Face "<<f<<" : "<<FaceCon[f]<<std::endl;
}

inline void BisectGeometry(Array<Vec3_t> Verts, Array<Array<int>> EdgeCon, Array<Array<int>> FaceCon, Vec3_t planeNormal, Vec3_t planePoint, Array<Vec3_t>& oVertsPos, Array<Array<int>>&oEdgeConPos,Array<Array<int>>& oFaceConPos, Array<int> & intersectionIxsPos, Array<Vec3_t>& oVertsNeg, Array<Array<int>>&oEdgeConNeg,Array<Array<int>>& oFaceConNeg, Array<int> & intersectionIxsNeg, double tolerance, bool debug){
  /*
   * BisectGeometry: Bisect a polyhedron into two using a cutting plane defined using a normal vector and a point in the plane. The direction of the plane's normal vector defines a positive and negative side of the 3D space which can be measured with the non-normalized distance. Therefore we can extract a polyhedron from the positive side and a polyhedron from the negative side of the plane.
   * Assumes the Face connectivity of ogParticle is done cyclically (default for mechsys) and all normals face outwards
   * INPUTS:
   * - Verts, EdgeCon, FaceCon: Starting vertices and their connection
   * - planeNormal: normal vector of the cutting plane
   * - planePoint: point on the cutting plane
   * - tolerance: how far away is it still considered a 0 distance from the plane (default 1e-5)
   * RETURNS:
   * oVertsPos, oEdgeCon, oFaceConPos: Vertices and connections for the geometry on the positive side
   * oVertsNeg, oEdgeCon, oFaceConNeg: Vertices and connections for the geometry on the negative side
   */
  std::cout<<"Bisecting geometry using plane normal: "<<planeNormal<<" with centroid at "<<planePoint<<"\n";
  //Calculate the signed distances to each vertex and count how many positive and negative verts there are
  size_t n_verts = Verts.Size();
  if (debug) PrintGeometry( Verts,  EdgeCon, FaceCon, "og ", "");

  Array<double> dist(n_verts);
  size_t n_pos = 0, n_neg = 0;
  for(size_t v=0; v<n_verts;v++){
    // dist[v] = dot(planeNormal, *ogParticle->Verts[v]-planePoint);
    dist[v] = dot(planeNormal, Verts[v]-planePoint);
    if(abs(dist[v])<tolerance){ dist[v] = 0.; if(debug) std::cout<<"Vrtex "<<Verts[v]<<" is an intersection\n";}
    if(dist[v]>= 0. ) n_pos++;
    if(dist[v]<= 0. ) n_neg++;
  }
  if(debug) std::cout<<"Number of vertices in the positive side: " <<n_pos << "\nNumber of vertices in the negative side: "<< n_neg<<"\n";
  //if all the vertices are positive or all are negative then the plane does not cut the particle and we can return it
  if(n_pos == n_verts){
    std::cout<<"Plane does not cut the geometry.\n";
    oVertsPos = Verts; oEdgeConPos = EdgeCon; oFaceConPos = FaceCon;
    oVertsNeg.Clear(); oEdgeConNeg.Clear(); oFaceConNeg.Clear();
    return;
  } else if (n_neg == n_verts){
    std::cout<<"Plane does not cut the geometry.\n";
    oVertsNeg = Verts; oEdgeConNeg = EdgeCon; oFaceConNeg = FaceCon;
    oVertsPos.Clear(); oEdgeConPos.Clear(); oFaceConPos.Clear();
    return;
  }

  Array<Vec3_t>       VertsPos(0), VertsNeg(0); // Lists of vertex coordinates for the positive and negative polyhedra
  Array<Array <int> > EdgeConPos(0), EdgeConNeg(0);// Edge connectivity arrays for positive and negative polyhedron
  Array<Array <int> > FaceConPos(0), FaceConNeg(0);// Face connectivity arrays for positive and negative polyhedron
  Array<Vec3_t> globalIntersections(0); // List of intersection vertex coordinates across the whole geometry
  size_t _ = -1; // temporary index, used to avoid declarations when it's not necessary for future calculations

  //BASIC LOGIC: loop over every edge of every face and look for intersections, you can store vertices in pos or neg lists depending on their distance to the plane and intersection points belong to both
  for(size_t i=0; i < FaceCon.Size(); i++){
    Array <int> current_face_cons = FaceCon[i];
    if(debug) std::cout<<"Calculating intersections for face "<<i<<"\t"<<current_face_cons<<"\n";
    // Arrays for face and edge connectivity obtained from the current face
    Array <int> new_FaceConPos(0), new_FaceConNeg(0);
    Array<int> tmp_EdgeConPos(0), tmp_EdgeConNeg(0); //Edge arrays to push if two vertices are both in pos/neg
    Array <Vec3_t> localIntersections; // store intersections for each face, we can use these vertices to close faces and edge connections
    //Here we ASSUME the faces are connected cyclically
    for(size_t j=0; j< current_face_cons.Size();j++){
      Array<int> new_EdgeConPos(0), new_EdgeConNeg(0); //Edge arrays to push if two vertices are both in pos/neg
      //Detect intersections and add them to the new arrays
      int vertex_ix1 = current_face_cons[j];//indices
      int vertex_ix2 = current_face_cons[(j+1)%current_face_cons.Size()];
      Vec3_t vertex1 = Verts[vertex_ix1];//Coordinates
      Vec3_t vertex2 = Verts[vertex_ix2];
      //Calculate intersection between a plane and an edge composed of two vertices, the number of intersecitons is returned and, if an intersection point is found, it is calculated to a vector passed by reference
      Vec3_t intersectionVertex(0., 0., 0.); //coordinates of the intersection point (if any)
      double dist1 = dist[vertex_ix1], dist2 = dist[vertex_ix2]; //signed distances from each vertex to the cutting plane
      size_t intersection_number = CalculateIntersection(planeNormal, planePoint, vertex1,vertex2, dist1, dist2, intersectionVertex, tolerance);
      if(debug) std::cout<<"Intersection number: "<<intersection_number<<" for vertices "<<vertex1<<" "<<vertex2<<"\n";
      //ADD VERTEX1 TO THE CORRESPONDING LIST
      if(dist1 > 0){
        AddVertexToConLists(vertex1, VertsPos, new_EdgeConPos, new_FaceConPos);
      }else if(dist1 <0){
        AddVertexToConLists(vertex1, VertsNeg, new_EdgeConNeg, new_FaceConNeg);
      }

      // ADD INTERSECTIONS
      // In the special case of 2 intersections (i.e. edge lies on the cutting plane) add them one by one
      if (intersection_number == 1){// One intersection, add vertex to both lists and to the local and global intersection lists
        if(debug) std::cout<<"Intersection point: "<<intersectionVertex<<"\n";
        if(!FindVectorInList(intersectionVertex, localIntersections, _)) localIntersections.Push(intersectionVertex);
        if(!FindVectorInList(intersectionVertex, globalIntersections, _)) globalIntersections.Push(intersectionVertex);
        AddVertexToConLists(intersectionVertex, VertsPos, new_EdgeConPos, new_FaceConPos);
        AddVertexToConLists(intersectionVertex, VertsNeg, new_EdgeConNeg, new_FaceConNeg);
      }

      // ADD VERTEX 2 TO THE CORRESPONDING LIST
      if(dist2 > 0){
        AddVertexToConLists(vertex2, VertsPos, new_EdgeConPos, new_FaceConPos);
      }else if(dist2 <0){
        AddVertexToConLists(vertex2, VertsNeg, new_EdgeConNeg, new_FaceConNeg);
      }
      //Push edges in case an edge can be formed
      if(new_EdgeConPos.Size()==2) {tmp_EdgeConPos.Push(new_EdgeConPos[0]); tmp_EdgeConPos.Push(new_EdgeConPos[1]);}
      if(new_EdgeConNeg.Size()==2) {tmp_EdgeConNeg.Push(new_EdgeConNeg[0]); tmp_EdgeConNeg.Push(new_EdgeConNeg[1]);}
    }
    // Use interactions to close edge connection lists
    if(localIntersections.Size()>1){
      for(size_t i=0; i<localIntersections.Size();i++){
        size_t ix_pos, ix_neg = -1;
        if(FindVectorInList(localIntersections[i], VertsPos, ix_pos)) tmp_EdgeConPos.Push(ix_pos);
        if(FindVectorInList(localIntersections[i], VertsNeg, ix_neg)) tmp_EdgeConNeg.Push(ix_neg);
      }
    }
    // Push faces and edges
    if(new_FaceConPos.Size()>2){
      FaceConPos.Push(new_FaceConPos);
      EdgeConPos.Push(tmp_EdgeConPos);
    }
    if(new_FaceConNeg.Size()>2){
      FaceConNeg.Push(new_FaceConNeg);
      EdgeConNeg.Push(tmp_EdgeConNeg);
    }
  }
  // Use global intersections to build the face that is cut, edge connections for the intersections should already be in the list
  Array<int> intersectionFaceConPos(0), intersectionFaceConNeg(0);
  size_t n_intersections = globalIntersections.Size();
  if(n_intersections>2){
    for(size_t i=0; i<n_intersections;i++){
      size_t ix_pos = -1, ix_neg = -1;
      if(FindVectorInList(globalIntersections[i], VertsPos, ix_pos)) intersectionFaceConPos.Push(ix_pos);
      if(FindVectorInList(globalIntersections[n_intersections-1-i], VertsNeg, ix_neg)) intersectionFaceConNeg.Push(ix_neg);
    }
    FaceConPos.Push(intersectionFaceConPos);
    FaceConNeg.Push(intersectionFaceConNeg);
  }

  //Reorder intersection face to keep normals outside
  if(FaceConPos.Size()>3){
    if (debug) PrintGeometry( VertsPos,  EdgeConPos, FaceConPos, "Positive ", "BEFORE reordering\n");
    //Reorder edges so they are arrays of size (n_edges,2)
    EdgeConPos = ReorderEdges(EdgeConPos);
    //Reorder the intersected face in case the vertices are not connected cyclically
    FaceConPos[FaceConPos.Size()-1] = ReorderFaceCon(VertsPos, EdgeConPos, FaceConPos[FaceConPos.Size()-1]);
    if (debug) PrintGeometry( VertsPos,  EdgeConPos, FaceConPos, "Positive ", "AFTER reordering\n");
  }
  if(FaceConNeg.Size()>3){
    if (debug) PrintGeometry( VertsNeg,  EdgeConNeg, FaceConNeg, "Negative ", "BEFORE reordering\n");
    //Reorder edges so they are arrays of size (n_edges,2)
    EdgeConNeg = ReorderEdges(EdgeConNeg);
    //Reorder the intersected face in case the vertices are not connected cyclically
    FaceConNeg[FaceConNeg.Size()-1] = ReorderFaceCon(VertsNeg, EdgeConNeg, FaceConNeg[FaceConNeg.Size()-1]);
    if (debug) PrintGeometry( VertsNeg,  EdgeConNeg, FaceConNeg, "Negative ", "AFTER reordering\n");
  }
  // Vertex and Faces should be done!
  std::cout<<"DONE bisecting geometry! \n";

  oVertsPos = VertsPos; oEdgeConPos = EdgeConPos; oFaceConPos = FaceConPos;
  oVertsNeg = VertsNeg; oEdgeConNeg = EdgeConNeg; oFaceConNeg = FaceConNeg;
  intersectionIxsPos = intersectionFaceConPos;  intersectionIxsNeg = intersectionFaceConNeg;
}

inline void BisectPolyhedron(DEM::Particle * ogParticle, Vec3_t planeNormal, Vec3_t planePoint, Array<DEM::Particle*> & segParticles, Array<Array<int>> & intersectionIxs, int mechsysErode, double tolerance, bool debug){
  /*
   * BisectPolyhedron: Bisect a polyhedron into two using a cutting plane defined using a normal vector and a point in the plane. The direction of the plane's normal vector defines a positive and negative side of the 3D space which can be measured with the non-normalized distance. Therefore we can extract a polyhedron from the positive side and a polyhedron from the negative side of the plane.
   * Assumes the Face connectivity of ogParticle is done cyclically (default for mechsys)
   * INPUTS:
   * - ogParticle: original Particle to bisect
   * - planeNormal: normal vector of the cutting plane
   * - planePoint: point on the cutting plane
   * - tolerance: how far away is it still considered a 0 distance from the plane (default 1e-5)
   * - mechsysErode:  0: Don't erode resulting particles, 1: use mechsys Erosion, 2: use internal face intrusion, 3: Simplify and then erosion
   * RETURNS:
   * segParticles: Array containing the two particles resulting from cutting the ogParticle, one on the postive side of the plane and one on the negative side
   */
  std::cout<<"Bisecting particle using plane normal: "<<planeNormal<<" with centroid at "<<planePoint<<"\n";
  //If the particle is eroded, dilate it so we can use its actual vertices to create actual particles and then erode from base mechsys
  Array<Vec3_t> Vtemp(ogParticle->Verts.Size());
  Array<Vec3_t> Vres (ogParticle->Verts.Size());
  Array<Vec3_t> Vog (ogParticle->Verts.Size());
  for (size_t j=0;j<ogParticle->Verts.Size();j++)
    {
      Vtemp[j] = *ogParticle->Verts[j];
      Vres [j] = *ogParticle->Verts[j];
      Vog [j] = *ogParticle->Verts[j];
    }
  if (ogParticle->Eroded&&ogParticle->Faces.Size()>=4&&mechsysErode!=4)
    {
      // std::cout<<"e_dilation\n";
      DEM::Dilation(Vtemp,ogParticle->EdgeCon,ogParticle->FaceCon,Vres,ogParticle->Props.R);
      // std::cout<<"e_dilation_after\n";
      std::cout<<"Eroded particle, dilating, previous vertices: \n"<<Vog<<"\n";
      // std::cout<<"Tmp vertices: \n"<<Vtemp<<"\n";
    }
  if(debug) PrintGeometry( Vres,  ogParticle->EdgeCon, ogParticle->FaceCon, "og", "Starting geometry (possible dilation): \n");

  //Now we can use the resulting vertices coordinates to measure distances
  Array<Vec3_t>       VertsPos(0), VertsNeg(0); // Lists of vertex coordinates for the positive and negative polyhedra
  Array<Array <int> > EdgeConPos(0), EdgeConNeg(0);// Edge connectivity arrays for positive and negative polyhedron
  Array<Array <int> > FaceConPos(0), FaceConNeg(0);// Face connectivity arrays for positive and negative polyhedron
  Array <int> intersectionIxsPos(0), intersectionIxsNeg(0);// Intersection indices for positive and negative polyhedron
  BisectGeometry(Vres, ogParticle->EdgeCon, ogParticle->FaceCon,
                 planeNormal, planePoint,
                 VertsPos, EdgeConPos, FaceConPos, intersectionIxsPos,
                 VertsNeg, EdgeConNeg, FaceConNeg, intersectionIxsNeg,
                 tolerance, debug);

  //----------------- Build the new particles -------------------
  //Push the positive particle given it's a 3D polyhedra
  std::cout<<"Pushing particles\n";
  // DEM::Domain vdom;
  // vdom.domID = 999999; //Visualization
  // vdom.domType = 1; //Subdomain
  // std::cout<<"Pushing Pos to visuzalization domain...\n";
  // vdom.Particles.Clear();
  // vdom.Particles.Push(new DEM::Particle(/*Tag*/-1,
  //                                       /*vertices*/VertsPos,
  //                                       /*edges*/EdgeConPos,
  //                                       /*faces*/FaceConPos,
  //                                       /*init vel*/ogParticle->v,
  //                                       /*init omega*/ogParticle->w,
  //                                       /*Radius*/ogParticle->Props.R,
  //                                       /*density*/ogParticle->Props.rho));
  // vdom.Dilate = false;
  // vdom.WriteXDMF("particle_cut_before_simplify1");
  // std::cout<<"Pushing Neg to visuzalization domain...\n";
  // vdom.Particles.Clear();
  // vdom.Particles.Push(new DEM::Particle(/*Tag*/-2,
  //                                       /*vertices*/VertsNeg,
  //                                       /*edges*/EdgeConNeg,
  //                                       /*faces*/FaceConNeg,
  //                                       /*init vel*/ogParticle->v,
  //                                       /*init omega*/ogParticle->w,
  //                                       /*Radius*/ogParticle->Props.R,
  //                                       /*density*/ogParticle->Props.rho));
  // vdom.WriteXDMF("particle_cut_before_simplify2");

  if(FaceConPos.Size()>3){
    //Simplify geometry by removing intersections that are very close to original vertices
    if (mechsysErode == 3){
      PrintGeometry( VertsPos,  EdgeConPos, FaceConPos, "Previous Pos", "Simplifying vertices... \n");
      std::cout<<"Previous R:\n"<<ogParticle->Props.R<<"\n";
      SimplifyGeometryNew(intersectionIxsPos, VertsPos, EdgeConPos, FaceConPos, ogParticle->Props.R, /*debug*/ true);
    }
    if (debug) PrintGeometry( VertsPos,  EdgeConPos, FaceConPos, "Positive", "AFTER simplify... \n");
    // vdom.Particles.Clear(); //Subdomain
    // std::cout<<"Pushing to visuzalization domain...\n";
    // vdom.Particles.Push(new DEM::Particle(/*Tag*/-1,
    //                                     /*vertices*/VertsPos,
    //                                     /*edges*/EdgeConPos,
    //                                     /*faces*/FaceConPos,
    //                                     /*init vel*/ogParticle->v,
    //                                     /*init omega*/ogParticle->w,
    //                                     /*Radius*/ogParticle->Props.R,
    //                                     /*density*/ogParticle->Props.rho));
    // vdom.Dilate = false;
    // vdom.WriteXDMF("particle_cut_before_erosion1");
  }
  if(FaceConNeg.Size()>3){
    //Simplify geometry by removing intersections that are very close to original vertices
    if (mechsysErode == 3){
      PrintGeometry( VertsNeg,  EdgeConNeg, FaceConNeg, "Previous Neg", "Simplifying vertices... \n");
      std::cout<<"Previous R:\n"<<ogParticle->Props.R<<"\n";
      SimplifyGeometryNew(intersectionIxsNeg, VertsNeg, EdgeConNeg, FaceConNeg, ogParticle->Props.R, /*debug*/ true);
    }
    if (debug) PrintGeometry( VertsNeg,  EdgeConNeg, FaceConNeg, "Negative ", "AFTER simplify \n");
    // vdom.Particles.Clear();
    // std::cout<<"Pushing to visuzalization domain...\n";
    // vdom.Particles.Push(new DEM::Particle(/*Tag*/-2,
    //                                     /*vertices*/VertsNeg,
    //                                     /*edges*/EdgeConNeg,
    //                                     /*faces*/FaceConNeg,
    //                                     /*init vel*/ogParticle->v,
    //                                     /*init omega*/ogParticle->w,
    //                                     /*Radius*/ogParticle->Props.R,
    //                                     /*density*/ogParticle->Props.rho));
    // vdom.Dilate = false;
    // vdom.WriteXDMF("particle_cut_before_erosion2");
  }

  //Erode or push intersection face towards the geometrical centre to avoid overlpping with other segmented particles
  if(FaceConPos.Size()>3){
    if (mechsysErode == 1 || mechsysErode == 3){
      PrintGeometry( VertsPos,  EdgeConPos, FaceConPos, "Previous ", "Eroding PosParticle via mechsys... \n");
      std::cout<<"Previous R:\n"<<ogParticle->Props.R<<"\n";
      DEM::Erosion(VertsPos, EdgeConPos, FaceConPos, ogParticle->Props.R); //, /*debug*/ false);
    }
    else if(mechsysErode == 2) IntrudeFace(VertsPos, FaceConPos.Last(), 1.00*ogParticle->Props.R);
    else if(mechsysErode == 4){
      if(debug) std::cout<<"Eroding PosParticle via internal bisection...\n";
      //Use the plane to bisect to get the particle with no overlaps
      Vec3_t planePoint2 = planePoint + ogParticle->Props.R*planeNormal;
      // We only care about the positive Vertices, 
      Array<Vec3_t> _v(0); Array<Array <int> > _e(0), _f(0); Array <int> _i(0);
      BisectGeometry(VertsPos, EdgeConPos, FaceConPos,
                     planeNormal, planePoint2,
                     VertsPos, EdgeConPos, FaceConPos, intersectionIxsPos,
                     _v, _e, _f, _i,
                     tolerance, debug);
      //Now the new positive vertices should be one speroradius away from the negative particle
    }
    if (debug) PrintGeometry( VertsPos,  EdgeConPos, FaceConPos, "Positive ", "AFTER erosion \n");
    // vdom.Particles.Clear(); //Subdomain
    // std::cout<<"Pushing to visuzalization domain...\n";
    // vdom.Particles.Push(new DEM::Particle(/*Tag*/-1,
    //                                       /*vertices*/VertsPos,
    //                                       /*edges*/EdgeConPos,
    //                                       /*faces*/FaceConPos,
    //                                       /*init vel*/ogParticle->v,
    //                                       /*init omega*/ogParticle->w,
    //                                       /*Radius*/ogParticle->Props.R,
    //                                       /*density*/ogParticle->Props.rho));
    // vdom.Dilate = false;
    // vdom.WriteXDMF("particle_cut_after_erosion1");
    std::cout<<"Pushing to segParts list...\n";
    segParticles.Push(new DEM::Particle(/*Tag*/ogParticle->Tag,
                                        /*vertices*/VertsPos,
                                        /*edges*/EdgeConPos,
                                        /*faces*/FaceConPos,
                                        /*init vel*/ogParticle->v,
                                        /*init omega*/ogParticle->w,
                                        /*Radius*/ogParticle->Props.R,
                                        /*density*/ogParticle->Props.rho));
    std::cout<<"Calculating props...\n";
    segParticles[segParticles.Size()-1]->poly_calc_props(VertsPos);// Initialize the physical values of the particle
    std::cout<<"Props calc'd...\n";
    if (ogParticle->Eroded) {
      segParticles[segParticles.Size()-1]->Eroded      = true;
      std::cout<<"Particle assigned as eroded\n";
    }
  }

  if(FaceConNeg.Size()>3){
    //Push intersection face towards the geometrical centre to avoid overlpping with other segmented particles
    if (mechsysErode == 1 || mechsysErode == 3){
      PrintGeometry( VertsNeg,  EdgeConNeg, FaceConNeg, "Previous ", "Eroding NegParticle via mechsys... \n");
      std::cout<<"Previous R:\n"<<ogParticle->Props.R<<"\n";
      DEM::Erosion(VertsNeg, EdgeConNeg, FaceConNeg, ogParticle->Props.R);//, /*debug*/ false);
    }
    else if(mechsysErode == 2) IntrudeFace(VertsNeg, FaceConNeg[FaceConNeg.Size()-1], 1.00*ogParticle->Props.R);
    else if(mechsysErode == 4){
      if(debug) std::cout<<"Eroding NegParticle via internal bisection...\n";
      //Use the plane to bisect to get the particle with no overlaps
      Vec3_t planePoint2 = planePoint - ogParticle->Props.R*planeNormal;
      // We only care about the negative Vertices,
      Array<Vec3_t> _v(0); Array<Array <int> > _e(0), _f(0); Array <int> _i(0);
      BisectGeometry(VertsNeg, EdgeConNeg, FaceConNeg,
                     planeNormal, planePoint2,
                     _v, _e, _f, _i,
                     VertsNeg, EdgeConNeg, FaceConNeg, intersectionIxsNeg,
                     tolerance, debug);
      //Now the new negative vertices should be one speroradius away from the positive particle
    }
    if (debug) PrintGeometry( VertsNeg,  EdgeConNeg, FaceConNeg, "Negitive ", "AFTER erosion \n");
    //vdom.Particles.Clear(); //Subdomain
    // std::cout<<"Pushing to visuzalization domain...\n";
    // vdom.Particles.Push(new DEM::Particle(/*Tag*/-2,
    //                                       /*vertices*/VertsNeg,
    //                                       /*edges*/EdgeConNeg,
    //                                       /*faces*/FaceConNeg,
    //                                       /*init vel*/ogParticle->v,
    //                                       /*init omega*/ogParticle->w,
    //                                       /*Radius*/ogParticle->Props.R,
    //                                       /*density*/ogParticle->Props.rho));
    // vdom.Dilate = false;
    // vdom.WriteXDMF("particle_cut_after_erosion3");
    // vdom.Particles.Clear(); //Subdomain
    // vdom.Particles.Push(new DEM::Particle(/*Tag*/-2,
    //                                       /*vertices*/VertsNeg,
    //                                       /*edges*/EdgeConNeg,
    //                                       /*faces*/FaceConNeg,
    //                                       /*init vel*/ogParticle->v,
    //                                       /*init omega*/ogParticle->w,
    //                                       /*Radius*/ogParticle->Props.R,
    //                                       /*density*/ogParticle->Props.rho));
    // vdom.Dilate = false;
    // vdom.WriteXDMF("particle_cut_after_erosion2");
    std::cout<<"Pushing...\n";
    segParticles.Push(new DEM::Particle(/*Tag*/ogParticle->Tag,
                                        /*vertices*/VertsNeg,
                                        /*edges*/EdgeConNeg,
                                        /*faces*/FaceConNeg,
                                        /*init vel*/ogParticle->v,
                                        /*init omega*/ogParticle->w,
                                        /*Radius*/ogParticle->Props.R,
                                        /*density*/ogParticle->Props.rho));
    std::cout<<"Calculating props...\n";
    segParticles[segParticles.Size()-1]->poly_calc_props(VertsNeg);// Initialize the physical values of the particle
    std::cout<<"Props calc'd...\n";
    if(ogParticle->Eroded) segParticles[segParticles.Size()-1]->Eroded      = true;
  }
  intersectionIxs.Clear();
  intersectionIxs.Push(intersectionIxsPos);
  intersectionIxs.Push(intersectionIxsNeg);
}

inline void SimplifyGeometry(Array<int>& intersectionIxs, Array<Vec3_t>& Verts, Array<Array<int>>&EdgeCon,Array<Array<int>>& FaceCon, double threshold, bool debug){
  /*
   *SimplifyGeometry: Given a list of vertices it searches in the list for the vertices closest to the intersections and replaces the intersections by those vertices, rewriting the new geometry with (possibly) fewer vertices
   */

  Array<int> newIntersectionIxs(0);
  size_t _ix;
  size_t _ix2;
  //1. Closest vector to the intersection within threshold that is not an intersection will be the one replacing it
  for(size_t i=0; i<intersectionIxs.Size(); i++){ //For every intersection
    int interIx = intersectionIxs[i];
    Vec3_t intersection = Verts[interIx];
    double min_dist = DBL_MAX;
    int min_vIx = INT_MAX;
    for(size_t v=0; v<Verts.Size(); v++){ //Look for the closest vertex for replacement
      if(FindIndexInConList(v, intersectionIxs, _ix)){ continue;} //Ignoring other intersections
      double v_dist = DEM::Distance(intersection, Verts[v]);
      if (v_dist< min_dist){
        min_dist = v_dist;
        min_vIx = v;
      }
    }
    if (min_dist < threshold) { //If the closest vertex is under the simplification threshold
      //save it under new intersections
      newIntersectionIxs.Push(min_vIx);
    } else {
      newIntersectionIxs.Push(interIx);
    }
  }
  if (debug) std::cout<<"Previous intersection ixs: "<<intersectionIxs<<"\nNew intersection ixs: "<<newIntersectionIxs<<"\n";
  //2. Use new intersections to generate new vertices
  Array<Vec3_t> newVerts;
  Array<size_t> vertsMap;
  Array<int> usedIxs(0);
  for(size_t v=0; v<Verts.Size(); v++){
    if(FindIndexInConList(v, intersectionIxs, _ix)){ //If vertex is an old intersection
      //Find it's new corresponding vertex
      int newInterIx = newIntersectionIxs[_ix];
      if (FindIndexInConList(newInterIx, usedIxs, _ix)){
        vertsMap.Push(usedIxs[_ix]);
      } else {//and it's not in the list yet
        usedIxs.Push(newInterIx);
        newVerts.Push(Verts[newInterIx]);
        vertsMap.Push(newVerts.Size()-1);
      }
    } else {
      if (FindIndexInConList(v, usedIxs, _ix2)){
        vertsMap.Push(usedIxs[_ix2]);
      } else {
        usedIxs.Push(v);
        newVerts.Push(Verts[v]);
        vertsMap.Push(newVerts.Size()-1);
      }
    }
  }

  if ( debug ) std::cout<<"SimplifyGeometry:: oldVerts "<<Verts<<"of size "<<Verts.Size()<<" \nnewVerts"<<newVerts<<"of size "<<newVerts.Size()<<"\n";
  if ( debug ) std::cout<<"SimplifyGeometry:: vertsMap "<<vertsMap<<"of size "<<vertsMap.Size()<<"\nUsed Indices "<<usedIxs<<"of size "<<usedIxs.Size()<<"\n";
  //3. Rewrite faces and edges, removing duplicates
  Array<Array<int>>newFaceCon;Array<Array<int>>newEdgeCon;
  for(size_t f=0; f<FaceCon.Size();f++){
    Array<int> _newFace(0);
    for(size_t v=0;v<FaceCon[f].Size();v++){
      int mappedIx = vertsMap[FaceCon[f][v]]; //Get the new index for the old vertex
      if(!FindIndexInConList(mappedIx,_newFace,_ix)){ //Add it if not there already
        _newFace.Push(mappedIx);
      }
    }
    if(_newFace.Size()<3){
      if (debug) std::cout<<"SimplifyGeometry::Less than 3 vertices, deleting face. oldFace: "<<FaceCon[f]<<"newFace"<<_newFace<<"\n";
    } else {
      newFaceCon.Push(_newFace);
    }
  }
  if ( debug ) std::cout<<"SimplifyGeometry:: oldFaceCon "<<FaceCon<<" of size "<<FaceCon.Size()<<" \nnewFaceCon"<<newFaceCon<<" of size "<<newFaceCon.Size()<<"\n";
  //Edges
  for(size_t e=0; e<EdgeCon.Size();e++){
    Array<int> _newEdge(0);
    for(size_t v=0;v<EdgeCon[e].Size();v++){
      int mappedIx = vertsMap[EdgeCon[e][v]]; //Get the new index for the old vertex
      if(!FindIndexInConList(mappedIx,_newEdge,_ix)){ //Add it if not there already
        _newEdge.Push(mappedIx);
      }
    }
    if(_newEdge.Size()<2){
      if (debug) std::cout<<"SimplifyGeometry::Less than 2 vertices, deleting edge. oldEdge: "<<EdgeCon[e]<<"newEdge"<<_newEdge<<"\n";
    } else {
      newEdgeCon.Push(_newEdge);
    }
  }
  if ( debug ) std::cout<<"SimplifyGeometry:: oldEdgeCon "<<EdgeCon<<" of size "<<EdgeCon.Size()<<" \nnewEdgeCon"<<newEdgeCon<<" of size "<<newEdgeCon.Size()<<"\n";

  Array<int> _newIntersectionIxs(newIntersectionIxs.Size());
  for(size_t i=0; i<newIntersectionIxs.Size();i++){
    _newIntersectionIxs[i] = vertsMap[newIntersectionIxs[i]]; //Get the new index for the old intersection
    }
  //Re assign to reference variables for passing
  Verts = newVerts;
  FaceCon = newFaceCon;
  EdgeCon = newEdgeCon;
  intersectionIxs = _newIntersectionIxs;
  //return 0;
}

inline void SimplifyGeometryNew(Array<int>& intersectionIxs, Array<Vec3_t>& Verts, Array<Array<int>>&EdgeCon,Array<Array<int>>& FaceCon, double threshold, bool debug){
  /*
   *SimplifyGeometry: Given a list of vertices it searches in the list for the vertices closest to the intersections and replaces the intersections by those vertices, rewriting the new geometry with (possibly) fewer vertices
   */

  Array<int> newVertsIxs(0);
  size_t _ix;
  // size_t _ix2;
  //1. Closest vector to the intersection within threshold that is not an intersection will be the one replacing it
  for(size_t w=0; w<Verts.Size(); w++){ //For every old vertex
    if(FindIndexInConList(w, intersectionIxs, _ix)){newVertsIxs.Push(w); continue;} //Ignoring intersections
    // int interIx = intersectionIxs[i];
    // Vec3_t intersection = Verts[interIx];
    double min_dist = DBL_MAX;
    int min_vIx = INT_MAX;
    for(size_t v=0; v<Verts.Size(); v++){ //Look for the closest vertex for replacement
      double v_dist = DEM::Distance(Verts[w], Verts[v]);
      if (v_dist< min_dist){
        min_dist = v_dist;
        min_vIx = v;
      }
    }
    if (min_dist < threshold) { //If the closest vertex is under the simplification threshold
      //save it under new intersections
      newVertsIxs.Push(min_vIx);
    } else {
      newVertsIxs.Push(w);
    }
  }
  if (debug) std::cout<<"New verts ixs: "<<newVertsIxs<<" of size "<<newVertsIxs.Size()<<"\n";
  //2. Use new intersections to generate new vertices
  Array<Vec3_t> newVerts;
  Array<size_t> vertsMap;
  Array<int> usedIxs(0);
  for(size_t v=0; v<Verts.Size(); v++){
    //Find it's new corresponding vertex
    int newVertIx = newVertsIxs[v];
    if (FindIndexInConList(newVertIx, usedIxs, _ix)){
      vertsMap.Push(usedIxs[_ix]);
    } else {//if it's not in the list yet
      usedIxs.Push(newVertIx);
      newVerts.Push(Verts[newVertIx]);
      vertsMap.Push(newVerts.Size()-1);
    }
  }

  if ( debug ) std::cout<<"SimplifyGeometry:: oldVerts "<<Verts<<"of size "<<Verts.Size()<<" \nnewVerts"<<newVerts<<"of size "<<newVerts.Size()<<"\n";
  if ( debug ) std::cout<<"SimplifyGeometry:: vertsMap "<<vertsMap<<"of size "<<vertsMap.Size()<<"\nUsed Indices "<<usedIxs<<"of size "<<usedIxs.Size()<<"\n";
  //3. Rewrite faces and edges, removing duplicates
  Array<Array<int>>newFaceCon;Array<Array<int>>newEdgeCon;
  for(size_t f=0; f<FaceCon.Size();f++){
    Array<int> _newFace(0);
    for(size_t v=0;v<FaceCon[f].Size();v++){
      int mappedIx = vertsMap[FaceCon[f][v]]; //Get the new index for the old vertex
      if(!FindIndexInConList(mappedIx,_newFace,_ix)){ //Add it if not there already
        _newFace.Push(mappedIx);
      }
    }
    if(_newFace.Size()<3){
      if (debug) std::cout<<"SimplifyGeometry::Less than 3 vertices, deleting face. oldFace: "<<FaceCon[f]<<"newFace"<<_newFace<<"\n";
    } else {
      newFaceCon.Push(_newFace);
    }
  }
  if ( debug ) std::cout<<"SimplifyGeometry:: oldFaceCon "<<FaceCon<<" of size "<<FaceCon.Size()<<" \nnewFaceCon"<<newFaceCon<<" of size "<<newFaceCon.Size()<<"\n";
  //Edges
  for(size_t e=0; e<EdgeCon.Size();e++){
    Array<int> _newEdge(0);
    for(size_t v=0;v<EdgeCon[e].Size();v++){
      int mappedIx = vertsMap[EdgeCon[e][v]]; //Get the new index for the old vertex
      if(!FindIndexInConList(mappedIx,_newEdge,_ix)){ //Add it if not there already
        _newEdge.Push(mappedIx);
      }
    }
    if(_newEdge.Size()<2){
      if (debug) std::cout<<"SimplifyGeometry::Less than 2 vertices, deleting edge. oldEdge: "<<EdgeCon[e]<<"newEdge"<<_newEdge<<"\n";
    } else {
      newEdgeCon.Push(_newEdge);
    }
  }
  if ( debug ) std::cout<<"SimplifyGeometry:: oldEdgeCon "<<EdgeCon<<" of size "<<EdgeCon.Size()<<" \nnewEdgeCon"<<newEdgeCon<<" of size "<<newEdgeCon.Size()<<"\n";

  Array<int> _newIntersectionIxs(intersectionIxs.Size());
  for(size_t i=0; i<intersectionIxs.Size();i++){
    _newIntersectionIxs[i] = vertsMap[intersectionIxs[i]]; //Get the new index for the old intersection
    }
  //Re assign to reference variables for passing
  Verts = newVerts;
  FaceCon = newFaceCon;
  EdgeCon = newEdgeCon;
  intersectionIxs = _newIntersectionIxs;
  //return 0;
}
inline void SimplifyGeometry_avg(Array<Vec3_t>& GlobalIntersections, Array<Vec3_t>& Verts, Array<Array<int>>&EdgeCon,Array<Array<int>>& FaceCon, double threshold, bool verbose){
  /*
   *SimplifyGeometry: Given a list of vertices it averages the vertices that are closer than a threshold and rewrites the new geometry with (possibly) fewer vertices
   */
  size_t _ix;
  //1. Build simplified GlobalIntersections
  Array<Vec3_t> newGlobalIntersections;
  for(size_t v=0; v<GlobalIntersections.Size(); v++){
    Array<Vec3_t> matches(0);
    Array<size_t> match_ixs(0);
    matches.Push(GlobalIntersections[v]);
    match_ixs.Push(v);
    for(size_t w=0; w<GlobalIntersections.Size(); w++) {
      if(DEM::Distance(GlobalIntersections[v], GlobalIntersections[w])< threshold) {
        matches.Push(GlobalIntersections[w]);
        match_ixs.Push(w);
      }
    }
    if (matches.Size()>1){//Vertex v is close to those in match_ixs
      Vec3_t avg = Vec3_t(0.,0.,0.);
      for(size_t i=0; i<matches.Size();i++) avg += matches[i];
      avg/=matches.Size();
      newGlobalIntersections.Push(avg);
    } else {
      if(!FindVectorInList(GlobalIntersections[v], newGlobalIntersections, _ix)) newGlobalIntersections.Push(GlobalIntersections[v]);
    }
  }
  if (verbose) std::cout<<"SimplifyGeometry:: oldGlobalIntersections "<<GlobalIntersections<<"of size "<<GlobalIntersections.Size()<<" \nnewGlobalIntersections"<<newGlobalIntersections<<"of size "<<newGlobalIntersections.Size()<<"\n";
  //Loop through the original vector list and create a new one with a mapping list
  //2. Build new vector list and mapping to old vertices
  Array<Vec3_t> newVerts;
  Array<size_t> vertsMap;
  size_t _ix2;
  for(size_t v=0; v<Verts.Size(); v++){
    if (FindVectorInList(Verts[v], newGlobalIntersections, _ix, threshold)){
      if (FindVectorInList(Verts[v], newVerts, _ix2, threshold)){
        vertsMap.Push(_ix2);
      }else{
        newVerts.Push(newGlobalIntersections[_ix]);
        vertsMap.Push(newVerts.Size()-1);
      }
    }else{
      newVerts.Push(Verts[v]);
      vertsMap.Push(newVerts.Size()-1);
    }
  }
  if (verbose) std::cout<<"SimplifyGeometry:: oldVerts "<<Verts<<"of size "<<Verts.Size()<<" \nnewVerts"<<newVerts<<"of size "<<newVerts.Size()<<"\n";
  if (verbose) std::cout<<"SimplifyGeometry:: vertsMap "<<vertsMap<<"of size "<<vertsMap.Size()<<"\n";
  //3. Rewrite faces and edges, removing duplicates
  Array<Array<int>>newFaceCon;Array<Array<int>>newEdgeCon;
  for(size_t f=0; f<FaceCon.Size();f++){
    Array<int> _newFace(0);
    for(size_t v=0;v<FaceCon[f].Size();v++){
      int mappedIx = vertsMap[FaceCon[f][v]]; //Get the new index for the old vertex
      if(!FindIndexInConList(mappedIx,_newFace,_ix)){ //Add it if not there already
        _newFace.Push(mappedIx);
      }
    }
    if(_newFace.Size()<3){
      std::cout<<"SimplifyGeometry::Less than 3 vertices, deleting face. oldFace: "<<FaceCon[f]<<"newFace"<<_newFace<<"\n";
    } else {
      newFaceCon.Push(_newFace);
    }
  }
  if (verbose) std::cout<<"SimplifyGeometry:: oldFaceCon "<<FaceCon<<" of size "<<FaceCon.Size()<<" \nnewFaceCon"<<newFaceCon<<" of size "<<newFaceCon.Size()<<"\n";
  //Edges
  for(size_t e=0; e<EdgeCon.Size();e++){
    Array<int> _newEdge(0);
    for(size_t v=0;v<EdgeCon[e].Size();v++){
      int mappedIx = vertsMap[EdgeCon[e][v]]; //Get the new index for the old vertex
      if(!FindIndexInConList(mappedIx,_newEdge,_ix)){ //Add it if not there already
        _newEdge.Push(mappedIx);
      }
    }
    if(_newEdge.Size()<2){
      std::cout<<"SimplifyGeometry::Less than 3 vertices, deleting edge. oldEdge: "<<EdgeCon[e]<<"newEdge"<<_newEdge<<"\n";
    } else {
      newEdgeCon.Push(_newEdge);
    }
  }
  if (verbose) std::cout<<"SimplifyGeometry:: oldEdgeCon "<<EdgeCon<<" of size "<<EdgeCon.Size()<<" \nnewEdgeCon"<<newEdgeCon<<" of size "<<newEdgeCon.Size()<<"\n";

  //Re assign to reference variables for passing
  Verts = newVerts;
  FaceCon = newFaceCon;
  EdgeCon = newEdgeCon;
  //return 0;
}
inline Array<size_t> GetLargest3(Array<double> Fs){
  if (Fs.Size()<3){
    std::cout<<"GetLargest3: Array size is smaller than three, can't construct plane with fewer than three poinst\n";
    throw new Fatal("Number of array is less than 3, cannot get the three highest");
  }
  //1. Set spots to -inf
  double third, second, first;
  third = second = first = DBL_MIN;
  Array<size_t> LargestIDs(3);
  for(size_t i=0; i<Fs.Size(); i++){
    //If current element is greater than the first, reorder
    if(Fs[i]>first){
      third = second;
      second = first;
      first = Fs[i];
      LargestIDs[0]=i;
    }
    else if(Fs[i]>second){//If current element is greater than the second reorder only the third
      third = second;
      second = Fs[i];
      LargestIDs[1]=i;
    }
    else if (Fs[i]>third){
      third = Fs[i];
      LargestIDs[2]=i;
    }
  }
  return LargestIDs;
}
#endif
