/**
 * @file   test_vtkAmeletHDFExporter.cxx
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Tue Nov 17 09:31:20 2015
 * 
 * @brief  
 * 
 * 
 */


#include <ah5.h>

#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>

#include <vtkAmeletHDFExporter.h>
#include <vtkAmeletHDFWriter.h>

// Declares a test suite named "sanity_check"
#include <boost/test/unit_test.hpp>

#include <string>
#include <map>
#include <vector>

std::vector<unsigned> range(unsigned start, unsigned stop) {
  std::vector<unsigned> xs(stop - start);
  std::vector<unsigned>::iterator x;
  unsigned val;
  for (x = xs.begin(), val = start; x != xs.end(); ++x, ++val)
    *x = val;
  return xs;
}

BOOST_AUTO_TEST_SUITE(test_vtkAmeletHDFExporter)

BOOST_AUTO_TEST_CASE(test_AddGgroup) {
  AH5Umesh converter;

  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkAmeletHDFWriter> ah5Writer =
    vtkSmartPointer<vtkAmeletHDFWriter>::New();
  ah5Writer->SetInputConnection(sphereSource->GetOutputPort());

  vtkPoints *pts;
  vtkCellArray *polys;
  vtkCellArray *strips;
  vtkPolyData *input = ah5Writer->GetInput();
  polys = input->GetPolys();
  strips = input->GetStrips();
  pts = input->GetPoints();

  converter.AddGroup("group1", pts, polys, strips);
  BOOST_CHECK_EQUAL(converter.nodes.size(), 150);
  BOOST_CHECK_EQUAL(converter.elementtypes.size(), 96);
  BOOST_CHECK_EQUAL(converter.elementnodes.size(), 288);
  BOOST_CHECK_EQUAL(converter.groups.size(), 1);
  std::map< std::string, std::vector<int> >::iterator grp
      = converter.groups.find("group1");
  BOOST_CHECK_EQUAL(grp->second.size(), 96);
  std::vector<unsigned> eles = range(0, 96);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      grp->second.begin(), grp->second.end(), eles.begin(), eles.end());
  BOOST_CHECK_EQUAL(converter.elementnodes.back(), 6);

  converter.AddGroup("group2", pts, polys, strips);
  BOOST_CHECK_EQUAL(converter.nodes.size(), 300);
  BOOST_CHECK_EQUAL(converter.elementtypes.size(), 192);
  BOOST_CHECK_EQUAL(converter.elementnodes.size(), 576);
  BOOST_CHECK_EQUAL(converter.groups.size(), 2);
  grp = converter.groups.find("group2");
  BOOST_CHECK_EQUAL(grp->second.size(), 96);
  eles = range(96, 192);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      grp->second.begin(), grp->second.end(), eles.begin(), eles.end());
  BOOST_CHECK_EQUAL(converter.elementnodes.back(), 56);
}


BOOST_AUTO_TEST_CASE(test_AddInput) {
  AH5Umesh converter;

  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkAmeletHDFWriter> ah5Writer =
    vtkSmartPointer<vtkAmeletHDFWriter>::New();
  ah5Writer->SetInputConnection(sphereSource->GetOutputPort());

  converter.AddInput("group1", ah5Writer->GetInput());
  BOOST_CHECK_EQUAL(converter.nodes.size(), 150);
  BOOST_CHECK_EQUAL(converter.elementtypes.size(), 96);
  BOOST_CHECK_EQUAL(converter.elementnodes.size(), 288);
  BOOST_CHECK_EQUAL(converter.groups.size(), 1);
  std::map< std::string, std::vector<int> >::iterator grp
      = converter.groups.find("group1");
  BOOST_CHECK_EQUAL(grp->second.size(), 96);
  std::vector<unsigned> eles = range(0, 96);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      grp->second.begin(), grp->second.end(), eles.begin(), eles.end());
  BOOST_CHECK_EQUAL(converter.elementnodes.back(), 6);

  converter.AddInput("group2", ah5Writer->GetInput());
  BOOST_CHECK_EQUAL(converter.nodes.size(), 300);
  BOOST_CHECK_EQUAL(converter.elementtypes.size(), 192);
  BOOST_CHECK_EQUAL(converter.elementnodes.size(), 576);
  BOOST_CHECK_EQUAL(converter.groups.size(), 2);
  grp = converter.groups.find("group2");
  BOOST_CHECK_EQUAL(grp->second.size(), 96);
  eles = range(96, 192);
  BOOST_CHECK_EQUAL_COLLECTIONS(
      grp->second.begin(), grp->second.end(), eles.begin(), eles.end());
  BOOST_CHECK_EQUAL(converter.elementnodes.back(), 56);
}


BOOST_AUTO_TEST_CASE(test_AH5Dump) {
  AH5Umesh converter;

  vtkSmartPointer<vtkSphereSource> sphereSource =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->Update();

  vtkSmartPointer<vtkAmeletHDFWriter> ah5Writer =
    vtkSmartPointer<vtkAmeletHDFWriter>::New();
  ah5Writer->SetInputConnection(sphereSource->GetOutputPort());

  vtkPoints *pts;
  vtkCellArray *polys;
  vtkCellArray *strips;
  vtkPolyData *input = ah5Writer->GetInput();
  polys = input->GetPolys();
  strips = input->GetStrips();
  pts = input->GetPoints();

  converter.AddGroup("group1", pts, polys, strips);
  converter.AddGroup("group2", pts, polys, strips);

  AH5_umesh_t mesh;
  converter.AH5Dump(&mesh);

  BOOST_CHECK_EQUAL(mesh.nb_elementnodes, 576);
  BOOST_CHECK_EQUAL(mesh.nb_elementtypes, 192);
  BOOST_CHECK_EQUAL(mesh.nb_nodes[0], 100);
  BOOST_CHECK_EQUAL(mesh.nb_nodes[1], 3);
  BOOST_CHECK_EQUAL(mesh.nb_groups, 2);
  BOOST_CHECK_EQUAL(mesh.groups[0].path, "group1");
  BOOST_CHECK_EQUAL(mesh.groups[1].path, "group2");
}

BOOST_AUTO_TEST_SUITE_END()
