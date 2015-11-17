/**
 * @file   vtkAmeletHDFExporter.h
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Fri Nov 13 14:31:42 2015
 * 
 * @brief  At end intended to implement vtkExporter.
 */

#ifndef VTKAMELETHDFEXPORTER_H
#define VTKAMELETHDFEXPORTER_H

#include <ah5.h>

#include <vector>
#include <map>
#include <string>


class vtkCellArray;
class vtkPoints;
class vtkPolyData;
class vtkDataObject;

struct AH5Umesh {
  std::vector<float> nodes;
  std::vector<char> elementtypes;
  std::vector<int> elementnodes;
  std::map< std::string, std::vector<int> > groups;

  char AddGroup(
      const std::string &name,
      vtkPoints *pts, vtkCellArray *polys, vtkCellArray *strips);

  char AddInput(const std::string &name, vtkDataObject *data);

  char AH5Dump(AH5_umesh_t *mesh);
};


#endif /* VTKAMELETHDFEXPORTER_H */
