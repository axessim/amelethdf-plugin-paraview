/**
 * @file   test_vtkAmeletHDFWriter.cxx
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Thu Nov 12 12:56:03 2015
 * 
 * @brief  
 * 
 * 
 */

#include <ah5.h>

#include <vtkPolyData.h>
#include <vtkAmeletHDFWriter.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>

// Declares a test suite named "sanity_check"
#include <boost/test/unit_test.hpp>

#include <string>

// Friend class
class vtkAmeletHDFWriter_UnitTest: public vtkAmeletHDFWriter {
 public:
  static vtkAmeletHDFWriter_UnitTest *New() {
    return new vtkAmeletHDFWriter_UnitTest();
  }
  using vtkAmeletHDFWriter::WriteAH5;
  vtkAmeletHDFWriter_UnitTest() {}
};


BOOST_AUTO_TEST_SUITE(test_vtkAmeletHDFWriter)

// Simple test using VTK API
BOOST_AUTO_TEST_CASE(test_Write) {
  std::string filename = "test_vtk_ah5_writer.h5";

  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkAmeletHDFWriter> ah5Writer =
    vtkSmartPointer<vtkAmeletHDFWriter>::New();
  ah5Writer->SetFileName(filename.c_str());
  ah5Writer->SetInputConnection(sphereSource->GetOutputPort());
  ah5Writer->Write();
}


BOOST_AUTO_TEST_CASE(test_WriteAH5) {
  AH5_umesh_t mesh;

  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkAmeletHDFWriter_UnitTest> ah5Writer =
    vtkSmartPointer<vtkAmeletHDFWriter_UnitTest>::New();
  ah5Writer->SetInputConnection(sphereSource->GetOutputPort());

  vtkPoints *pts;
  vtkCellArray *polys;
  vtkCellArray *strips;
  vtkPolyData *input = ah5Writer->GetInput();
  polys = input->GetPolys();
  strips = input->GetStrips();
  pts = input->GetPoints();

  ah5Writer->WriteAH5(&mesh, pts, polys, strips);

  BOOST_CHECK_EQUAL(mesh.nb_elementnodes, 288);
  BOOST_CHECK_EQUAL(mesh.nb_elementtypes, 96);
  BOOST_CHECK_EQUAL(mesh.nb_nodes[0], 50);
  BOOST_CHECK_EQUAL(mesh.nb_nodes[1], 3);
}


BOOST_AUTO_TEST_SUITE_END()
