/**
 * @file   vtkAmeletHDFWriter.cxx
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Mon Nov  9 17:46:32 2015
 * 
 * @brief  
 * 
 * 
 */

#include "vtkAmeletHDFWriter.h"

#include "vtkSmartPointer.h"
#include "vtkObjectFactory.h"
#include "vtkErrorCode.h"
#include "vtkInformation.h"

#include "vtkPolyData.h"
#include "vtkCellArray.h"
#include "vtkTriangle.h"
#include "vtkTriangleStrip.h"

#include <vector>
#include <algorithm>
#include <string>

vtkStandardNewMacro(vtkAmeletHDFWriter);


vtkAmeletHDFWriter::vtkAmeletHDFWriter() : FileName(NULL) {}


void vtkAmeletHDFWriter::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);

  os << indent << "FileName: "
     << ((this->GetFileName() == NULL) ?
         "(none)" : this->GetFileName()) << std::endl;
  os << indent << "Input: " << this->GetInput() << std::endl;
}


vtkPolyData* vtkAmeletHDFWriter::GetInput() {
  return vtkPolyData::SafeDownCast(this->GetInput(0));
}


vtkPolyData* vtkAmeletHDFWriter::GetInput(int port) {
  return vtkPolyData::SafeDownCast(this->Superclass::GetInput(port));
}


std::string vtkAmeletHDFWriter::GetMeshName() {
  if (FileName) {
    std::string fname(FileName);
    std::size_t found = fname.find_last_of(".");
    return fname.substr(0, found);
  }
  return std::string();
}


void vtkAmeletHDFWriter::WriteData() {
  // Check data and call converter.
  vtkPoints *pts;
  vtkCellArray *polys;
  vtkCellArray *strips;
  vtkPolyData *input = this->GetInput();

  polys = input->GetPolys();
  strips = input->GetStrips();
  pts = input->GetPoints();
  if (pts == NULL || polys == NULL) {
    vtkErrorMacro(<<"No data to write!");
    this->SetErrorCode(vtkErrorCode::UnknownError);
    return;
  }

  if (this->FileName == NULL) {
    vtkErrorMacro(<< "Please specify FileName to write");
    this->SetErrorCode(vtkErrorCode::NoFileNameError);
    return;
  }

  this->WriteAH5(pts, polys, strips);
  if (this->ErrorCode == vtkErrorCode::OutOfDiskSpaceError) {
    vtkErrorMacro("Ran out of disk space; deleting file: "
                  << this->FileName);
    unlink(this->FileName);
  }
}


char vtkAmeletHDFWriter::WriteAH5(
    vtkPoints *pts, vtkCellArray *polys, vtkCellArray *strips) {
  // Initialed AH5 mesh file and call VTK 2 AH5 converter.

  AH5_mesh_t mesh;
  hid_t fid;

  AH5_init_mesh(&mesh, /*nb_groups=*/1);
  AH5_init_msh_group(
      mesh.groups, "vtk2ah5",
      /*nb_meshs=*/1, /*nb_mesh_links=*/0);
  AH5_init_msh_instance(
      mesh.groups->msh_instances, GetMeshName().c_str(), MSH_UNSTRUCTURED);

  if (WriteAH5(
          &mesh.groups->msh_instances->data.unstructured,
          pts, polys, strips))
    return EXIT_FAILURE;

  if ((fid = AH5_create(
          this->FileName, H5F_ACC_TRUNC, /*entryPoint=*/"/mesh/vtk2ah5")) < 0) {
    vtkErrorMacro(<< "Couldn't open file: " << this->FileName);
    this->SetErrorCode(vtkErrorCode::CannotOpenFileError);
    return EXIT_FAILURE;
  }

  AH5_write_mesh(fid, &mesh);
  H5Fclose(fid);

  AH5_free_mesh(&mesh);

  return EXIT_SUCCESS;
}


char vtkAmeletHDFWriter::WriteAH5(
    AH5_umesh_t *mesh,
    vtkPoints *pts, vtkCellArray *polys, vtkCellArray *strips) {
  // Define umesh sizes.
  hsize_t nb_elementnodes = 0;
  hsize_t nb_elementtypes = 0;
  hsize_t nb_nodes = pts->GetNumberOfPoints();
  hsize_t nb_groups = 1;
  hsize_t nb_groupgroups = 0;
  hsize_t nb_som_tables = 0;

  std::vector<char> elementtypes;
  std::vector<int> elementnodes;

  vtkIdType npts = 0;
  vtkIdType *indx = 0;
  double point[3];

  // Decompose any triangle strips into triangles
  vtkSmartPointer<vtkCellArray> polyStrips =
      vtkSmartPointer<vtkCellArray>::New();
  if (strips->GetNumberOfCells() > 0) {
    vtkIdType *ptIds = 0;
    for (strips->InitTraversal(); strips->GetNextCell(npts, ptIds);) {
      vtkTriangleStrip::DecomposeStrip(npts, ptIds, polyStrips);
    }
  }

  // Convert out triangle strips
  for (polyStrips->InitTraversal(); polyStrips->GetNextCell(npts, indx); ) {
    elementtypes.push_back(UELE_TRI3);
    for (unsigned idx = 0; idx < 3; ++idx)
      elementnodes.push_back(indx[idx]);
  }

  for (polys->InitTraversal(); polys->GetNextCell(npts, indx); ) {
    if (npts > 3) {
      vtkErrorMacro(<<"STL file only supports triangles");
      this->SetErrorCode(vtkErrorCode::FileFormatError);
      return EXIT_FAILURE;
    }

    elementtypes.push_back(UELE_TRI3);
    for (unsigned idx = 0; idx < 3; ++idx)
      elementnodes.push_back(indx[idx]);
  }

  // Write the mesh
  nb_elementnodes = elementnodes.size();
  nb_elementtypes = elementtypes.size();

  AH5_init_umesh(
      mesh, nb_elementnodes, nb_elementtypes,
      nb_nodes, nb_groups, nb_groupgroups, nb_som_tables);

  for (vtkIdType iPoint = 0; iPoint < pts->GetNumberOfPoints(); ++iPoint) {
    pts->GetPoint(iPoint, point);
    for (unsigned idx = 0; idx < 3; ++idx)
      mesh->nodes[iPoint*3+idx] = static_cast<float>(point[idx]);
  }

  std::copy(elementnodes.begin(), elementnodes.end(), mesh->elementnodes);
  std::copy(elementtypes.begin(), elementtypes.end(), mesh->elementtypes);

  AH5_init_umsh_group(
      mesh->groups, GetMeshName().c_str(), nb_elementtypes,
      GROUP_ELEMENT, GROUP_FACE);
  for (unsigned i = 0; i < nb_elementtypes; ++i)
    mesh->groups[0].groupelts[i] = i;

  return EXIT_SUCCESS;
}


int vtkAmeletHDFWriter::FillInputPortInformation(
    int /*port*/, vtkInformation *info) {
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
  return 1;
}
