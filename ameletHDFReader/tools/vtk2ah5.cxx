/**
 * @file   ah52vtk.cxx
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Mon Nov  9 16:22:41 2015
 * 
 * @brief  Convert AH5 file to vtk
 * 
 * 
 */

#include <vtkSmartPointer.h>
#include <vtkGenericDataObjectReader.h>

#include "vtkAmeletHDFWriter.h"

#include <string>

int main(int argc, char *argv[]) {
  if (argc < 3) {
    std::cerr << "Required arguments: input.vtk output.h5" << std::endl;
    return EXIT_FAILURE;
  }

  std::string inputFileName = argv[1];
  std::string outputFileName = argv[2];

  vtkSmartPointer<vtkGenericDataObjectReader> reader
      = vtkSmartPointer<vtkGenericDataObjectReader>::New();
  reader->SetFileName(inputFileName.c_str());
  reader->Update();

  vtkSmartPointer<vtkAmeletHDFWriter> writer
      = vtkSmartPointer<vtkAmeletHDFWriter>::New();
  writer->SetFileName(outputFileName.c_str());
  writer->SetInputConnection(reader->GetOutputPort());
  writer->Update();

  return EXIT_SUCCESS;
}
