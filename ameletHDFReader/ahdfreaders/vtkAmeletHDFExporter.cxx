/**
 * @file   vtkAmeletHDFExporter.cxx
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Fri Nov 13 14:40:34 2015
 * 
 * @brief  
 * 
 * 
 */

#include "vtkAmeletHDFExporter.h"

#include "vtkSmartPointer.h"
// #include "vtkErrorCode.h"
// #include "vtkInformation.h"

#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkTriangleStrip.h"

#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <utility>


char AH5Umesh::AddGroup(
    const std::string &name,
    vtkPoints *pts, vtkCellArray *polys, vtkCellArray *strips) {
  vtkIdType npts = 0;
  vtkIdType *indx = 0;
  double point[3];

  std::vector<float>::size_type nodes_offset = nodes.size();
  nodes.reserve(nodes_offset +  pts->GetNumberOfPoints());
  for (vtkIdType iPoint = 0; iPoint < pts->GetNumberOfPoints(); ++iPoint) {
    pts->GetPoint(iPoint, point);
    for (unsigned idx = 0; idx < 3; ++idx)
      nodes.push_back(static_cast<float>(point[idx]));
  }
  nodes_offset /= 3;

  // Decompose any triangle strips into triangles
  vtkSmartPointer<vtkCellArray> polyStrips =
      vtkSmartPointer<vtkCellArray>::New();
  if (strips->GetNumberOfCells() > 0) {
    vtkIdType *ptIds = 0;
    for (strips->InitTraversal(); strips->GetNextCell(npts, ptIds);) {
      vtkTriangleStrip::DecomposeStrip(npts, ptIds, polyStrips);
    }
  }

  //
  std::vector<char>::size_type eles_offset = elementtypes.size();

  // Convert out triangle strips
  for (polyStrips->InitTraversal(); polyStrips->GetNextCell(npts, indx); ) {
    elementtypes.push_back(UELE_TRI3);
    for (unsigned idx = 0; idx < 3; ++idx)
      elementnodes.push_back(nodes_offset + indx[idx]);
  }

  for (polys->InitTraversal(); polys->GetNextCell(npts, indx); ) {
    if (npts > 3) {
      return EXIT_FAILURE;
    }

    elementtypes.push_back(UELE_TRI3);
    for (unsigned idx = 0; idx < 3; ++idx)
      elementnodes.push_back(nodes_offset + indx[idx]);
  }

  std::vector<int> group(elementtypes.size() - eles_offset);
  std::vector<int>::iterator x;
  unsigned i;
  for (x = group.begin(), i = eles_offset; x != group.end(); ++i, ++x)
    *x = i;
  groups.insert(std::pair< std::string, std::vector<int> >(name, group));

  return EXIT_SUCCESS;
}


char AH5Umesh::AddInput(const std::string &name, vtkDataObject *data) {
  vtkPolyData* poly = vtkPolyData::SafeDownCast(data);

  vtkPoints *pts;
  vtkCellArray *polys;
  vtkCellArray *strips;

  polys = poly->GetPolys();
  strips = poly->GetStrips();
  pts = poly->GetPoints();

  return AddGroup(name, pts, polys, strips);
}


char AH5Umesh::AH5Dump(AH5_umesh_t *mesh) {
  // Define umesh sizes.
  hsize_t nb_groupgroups = 0;
  hsize_t nb_som_tables = 0;

  AH5_init_umesh(
      mesh,
      elementnodes.size(), elementtypes.size(),
      nodes.size() / 3, groups.size(),
      nb_groupgroups, nb_som_tables);

  std::copy(nodes.begin(), nodes.end(), mesh->nodes);
  std::copy(elementnodes.begin(), elementnodes.end(), mesh->elementnodes);
  std::copy(elementtypes.begin(), elementtypes.end(), mesh->elementtypes);

  unsigned grp_idx = 0;
  std::map< std::string, std::vector<int> >::const_iterator iGroup;
  for (iGroup = groups.begin(), grp_idx = 0; iGroup != groups.end(); ++iGroup, ++grp_idx) {
    AH5_init_umsh_group(
      mesh->groups + grp_idx,
      iGroup->first.c_str(), iGroup->second.size(),
      GROUP_ELEMENT, GROUP_FACE);
    std::copy(iGroup->second.begin(), iGroup->second.end(),
              mesh->groups[grp_idx].groupelts);
  }

  return EXIT_SUCCESS;
}

