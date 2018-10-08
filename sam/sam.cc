#include "sam.h"
#include "samSz.h"
#include <apfGeometry.h>
#include <apfField.h>
#include <apfShape.h>
#include <pcu_util.h>
#include <algorithm>
#include <stdio.h>
#include <math.h>

namespace sam {

apf::Field* errorThreshold(apf::Mesh* m, const char* fieldName,
    const unsigned idx, const double limit, const double factor)
{
  apf::Field* f = m->findField(fieldName);
  PCU_ALWAYS_ASSERT(f);
  apf::synchronize(f);
  apf::Field* newSz = apf::createFieldOn(m,"specifiedIso",apf::SCALAR);
  apf::Field* curSz = samSz::isoSize(m);
  double* vals = new double[f->countComponents()];
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    apf::getComponents(f, vtx, 0, vals);
    double h = apf::getScalar(curSz,vtx,0);
    if ( vals[idx] > limit )
      h *= factor;
    apf::setScalar(newSz,vtx,0,h);
  }
  m->end(itr);
  apf::destroyField(curSz);
  delete [] vals;
  return newSz;
}

apf::Field* computeHessian(apf::Mesh* m, const char* fieldName, const unsigned idx)
{
  apf::Field* f = m->findField(fieldName);
  PCU_ALWAYS_ASSERT(f);
  apf::synchronize(f);
  double* vals = new double[f->countComponents()];
  apf::Field *speedF = apf::createLagrangeField(m, "hess_speed", apf::SCALAR, 1);
  apf::MeshEntity* v;
  apf::MeshIterator* it = m->begin(0);
  apf::Vector3 vel_vect;
  while( (v = m->iterate(it)) ) { 
    apf::getComponents(f, v, 0, vals);
    vel_vect = apf::Vector3(vals[idx],vals[idx+1],vals[idx+2]);
    double speed = vel_vect.getLength(); //speed is velocity magnitude scalar
    apf::setScalar(speedF, v, 0, speed); 
  }
  m->end(it);
  delete [] vals;
  apf::Field *gradSpeed = apf::recoverGradientByVolume(speedF);
  apf::Field *grad2Speed = apf::recoverGradientByVolume(gradSpeed); //second derivative of speed field
  apf::Field *metricf = createLagrangeField(m, "hess_metric", apf::MATRIX, 1);
  it = m->begin(0);
  while ((v = m->iterate(it)))
  {
    apf::Matrix3x3 g2S;
    apf::getMatrix(grad2Speed, v, 0, g2S);
    apf::Matrix3x3 g2St = apf::transpose(g2S);
    apf::Matrix3x3 metric = (g2S + g2St) / 2; //Matrix is symmetrized here.
    apf::setMatrix(metricf, v, 0, metric);
  }
  m->end(it);
  apf::destroyField(speedF);
  apf::destroyField(gradSpeed);
  apf::destroyField(grad2Speed); 
  return metricf;  
}

struct SortingStruct
{
  apf::Vector3 v;
  double wm;
  bool operator<(const SortingStruct &other) const
  {
    return wm < other.wm;
  }
};

void averageToEntity(apf::Field *ef, apf::Field *vf, apf::MeshEntity *ent)
  /*
    Function used to convert a region/element field into a vertex field via averaging.
    For each vertex, the average value of the adjacent elements is taken
    Input:
      ef is the element field
      ent is the target vertex
    Output:
      vf is the vertex field
  */
{
  apf::Mesh *m = apf::getMesh(ef);
  apf::Adjacent elements;
  m->getAdjacent(ent, m->getDimension(), elements);
  double s = 0; 
  for (std::size_t i = 0; i < elements.getSize(); ++i) 
    s += apf::getScalar(ef, elements[i], 0);
  s /= elements.getSize();
  apf::setScalar(vf, ent, 0, s);
  return;
}

void computeAnisoSzFromHessian(apf::Mesh* m, apf::Field* hessian, apf::Field* sz_scale, apf::Field* sz_frame)
{
  /*
    Both scale and frames fields calcualted in this function are attached
    back to the mesh. The existing mesh size field should be evaluated first in
    order to create a new size field baesd on Hessian metric information.
  */
  //estimate the existing size field for each element.
  int nsd = m->getDimension();
  apf::Field *size_reg = apf::createField(m, "reg_size", apf::SCALAR, apf::getConstant(nsd));//place a node on every region
  apf::Field *size_vtx = apf::createLagrangeField(m, "vtx_size", apf::SCALAR, 1);
  apf::MeshIterator *it = m->begin(nsd); //region is at the dimension 3 (given vertex is 0)
  apf::MeshEntity *r;
  while (r = m->iterate(it))
  {
    double h = apf::computeShortestHeightInTet(m,r);
    apf::setScalar(size_reg, r, 0, h);
  }
  m->end(it);
  //lift the size field from elements to vertices 
  it = m->begin(0);
  apf::MeshEntity* v;
  while ((v = m->iterate(it)))
  {
    averageToEntity(size_reg, size_vtx, v);
  }
  m->end(it);

  apf::Vector3 scale;
  it = m->begin(0);
  while(v = m->iterate(it)) {
  apf::Matrix3x3 metric;
  apf::getMatrix(hessian, v, 0, metric);
  apf::Vector3 eigenVectors[3];
  double eigenValues[3];
  apf::eigen(metric, eigenVectors, eigenValues);
  /*
    Sort the eigenvalues and corresponding vectors.
    Larger eigenvalues means a need for a finer mesh.
  */
  SortingStruct ssa[3];
  for (int i = 0; i < 3; ++i) 
  {    
    ssa[i].v = eigenVectors[i];
    ssa[i].wm = std::fabs(eigenValues[i]);
  }    
  std::sort(ssa, ssa + 3);
  PCU_ALWAYS_ASSERT(ssa[2].wm >= ssa[1].wm);
  PCU_ALWAYS_ASSERT(ssa[1].wm >= ssa[0].wm);
  double firstEigenvalue = ssa[2].wm;
  PCU_ALWAYS_ASSERT(firstEigenvalue > 1e-12);
  //produce scale field based on eigenvalus. 
  double h_ref = apf::getScalar(size_vtx, v, 0);
  scale[0] = h_ref; //samllest size is associated with the direction with largest eigenvalue
  scale[1] = sqrt(ssa[2].wm / ssa[1].wm) * scale[0];
  scale[2] = sqrt(ssa[2].wm / ssa[0].wm) * scale[0]; 
  apf::setVector(sz_scale, v, 0, scale);
  //Some protection may be needed here to limit the aspect ratio.

  //produce frames field based on normalized eigenvectors.
  apf::Matrix3x3 frame(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0);
  frame[0] = ssa[2].v;
  frame[1] = ssa[1].v;
  frame[2] = ssa[0].v;
  for (int i = 0; i < 3; ++i)
    frame[i] = frame[i].normalize();
  frame = apf::transpose(frame);
  apf::setMatrix(sz_frame, v, 0, frame); 
  }
  m->end(it);
  apf::destroyField(size_reg);
  apf::destroyField(size_vtx);
  apf::synchronize(sz_scale);
  apf::synchronize(sz_frame);
}

apf::Field* compareIsoSF(apf::Mesh* m, const char* desiredSzFld, int method)
{
  apf::Field* f = m->findField(desiredSzFld);
  PCU_ALWAYS_ASSERT(f);
  PCU_ALWAYS_ASSERT(f->countComponents() == 1);
  apf::synchronize(f);
  apf::Field* newSz = apf::createFieldOn(m,"compareIsoSF",apf::SCALAR);
  apf::Field* curSz = samSz::isoSize(m);
  double dsr = 0.0;
  double cur = 0.0;
  double rsl = 0.0;
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    dsr = apf::getScalar(f, vtx, 0);
    cur = apf::getScalar(curSz,vtx,0);
    switch (method) {
      case 1:  rsl = (cur-dsr) / dsr; break; // ratio of difference
      default: rsl = cur / dsr;
    }
    apf::setScalar(newSz,vtx,0,rsl);
  }
  m->end(itr);
  apf::destroyField(curSz);
  return newSz;
}

apf::Field* specifiedIso(apf::Mesh* m, const char* fieldName, const unsigned idx)
{
  apf::Field* f = m->findField(fieldName);
  PCU_ALWAYS_ASSERT(f);
  apf::synchronize(f);
  apf::Field* newSz = apf::createFieldOn(m,"specifiedIso",apf::SCALAR);
  double* vals = new double[f->countComponents()];
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    apf::getComponents(f, vtx, 0, vals);
    double h = vals[idx];
    apf::setScalar(newSz,vtx,0,h);
  }
  m->end(itr);
  delete [] vals;
  return newSz;
}

void multiplySF(apf::Mesh* m, apf::Field* sf, double factor) {
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    double h = apf::getScalar(sf,vtx,0);
    apf::setScalar(sf,vtx,0,h*factor);
  }
  m->end(itr);
}

void multiplySFBox(apf::Mesh* m, apf::Field* sf, double factor, double* box) {
  PCU_ALWAYS_ASSERT(box[3] > 0 && box[4] > 0 && box[5] > 0);
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  apf::Vector3 points;
  apf::Vector3 center = apf::Vector3(box[0],box[1],box[2]);
  apf::Vector3 size   = apf::Vector3(box[3],box[4],box[5]);
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    m->getPoint(vtx, 0, points);
    if (apf::withinBox(points, center, size)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

void multiplySFCyl(apf::Mesh* m, apf::Field* sf, double factor, double* cyl) {
  /* cylinder is defined as {center_x, center_y, center_z,
     normal_x, normal_y, normal_z, half_height, radius}   */
  PCU_ALWAYS_ASSERT(cyl[6] > 0 && cyl[7] > 0);
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  apf::Vector3 points;
  apf::Vector3 center = apf::Vector3(cyl[0],cyl[1],cyl[2]);
  apf::Vector3 normal = apf::Vector3(cyl[3],cyl[4],cyl[5]);
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    m->getPoint(vtx, 0, points);
    if (apf::withinCyl(points, center, cyl[6], cyl[7], normal)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

void multiplySFRegion(apf::Mesh* m, apf::Field* sf, double factor, int tag) {
  /* tag should be a tag of a region */
  PCU_ALWAYS_ASSERT(m->findField(apf::getName(sf)));
  int type = 3; // 3D region by default
  apf::MeshEntity* vtx;
  apf::MeshIterator* itr = m->begin(0);
  while( (vtx = m->iterate(itr)) ) {
    apf::ModelEntity* m_vtx = m->toModel(vtx);
    apf::ModelEntity* region = m->findModelEntity(type, tag);
    if (m_vtx == region ||
	    m->isInClosureOf(m_vtx, region)) {
      double h = apf::getScalar(sf,vtx,0);
      apf::setScalar(sf,vtx,0,h*factor);
	}
  }
  m->end(itr);
}

}
