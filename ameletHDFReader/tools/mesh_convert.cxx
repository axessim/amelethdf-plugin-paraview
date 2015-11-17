/**
 * @file   mesh_convert.cxx
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Mon Nov  9 16:22:41 2015
 * 
 * @brief Convert mesh file
 *
 * from: (vtk, obj, ply, stl (stereo lithography), dem (elevation map))
 * to: (stl, ply, vtp, vtk)
 * 
 * 
 */


#include <ah5.h>

#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>
#include <vtkAmeletHDFExporter.h>


#include <string>


std::string GetMeshName(std::string fname) {
  return fname.substr(0, fname.find_last_of("."));
}


int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Required arguments: output.h5 input.vtk [*.vtk]" << std::endl;
    return EXIT_FAILURE;
  }

  vtkSmartPointer<vtkGenericDataObjectReader> reader
      = vtkSmartPointer<vtkGenericDataObjectReader>::New();
  AH5Umesh converter;
  for (int inputid = 2; inputid < argc; ++inputid) {
    std::string input_file_name = argv[inputid];
    reader->SetFileName(input_file_name.c_str());
    reader->Update();

    std::cout << "add inout: " << input_file_name << std::endl;
    converter.AddInput(GetMeshName(input_file_name), reader->GetOutput());
  }

  std::string output_file_name = argv[1];
  std::string mesh_name = GetMeshName(output_file_name);

  AH5_mesh_t mesh;
  AH5_init_mesh(&mesh, /*nb_groups=*/1);
  AH5_init_msh_group(
      mesh.groups, mesh_name.c_str(),
      /*nb_meshs=*/1, /*nb_mesh_links=*/0);
  AH5_init_msh_instance(
      mesh.groups->msh_instances, mesh_name.c_str(), MSH_UNSTRUCTURED);

  std::cout << "convert to AH5" << std::endl;
  converter.AH5Dump(&mesh.groups->msh_instances->data.unstructured);

  std::cout << "write output file: " << output_file_name << std::endl;
  hid_t fid;
  // Build AH5 entry point
  std::string entry_point = std::string("/mesh/") + mesh_name;
  if ((fid = AH5_create(
          output_file_name.c_str(), H5F_ACC_TRUNC, entry_point.c_str())) < 0) {
    return EXIT_FAILURE;
  }

  AH5_write_mesh(fid, &mesh);
  H5Fclose(fid);

  AH5_free_mesh(&mesh);

  return EXIT_SUCCESS;
}
