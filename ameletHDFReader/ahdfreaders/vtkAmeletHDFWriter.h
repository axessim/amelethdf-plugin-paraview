/**
 * @file   vtkAmeletHDFWriter.h
 * @author Nathanaël MUOT <nathanael.muot@axessim.fr>
 * @date   Mon Nov  9 16:30:26 2015
 * 
 * @brief  
 * 
 * 
 */

#ifndef VTKAMELETHDFWRITER_H
#define VTKAMELETHDFWRITER_H

#include <ah5.h>

#include <vtkIOGeometryModule.h>  // For export macro
#include <vtkWriter.h>

#include <string>

class vtkCellArray;
class vtkPoints;
class vtkPolyData;


class vtkAmeletHDFWriter: public vtkWriter {
 public:
  // vtkWriter implementation
  static vtkAmeletHDFWriter *New();
  vtkTypeMacro(vtkAmeletHDFWriter, vtkWriter);
  virtual void PrintSelf(ostream& os, vtkIndent indent);  // NOLINT

  // Description:
  // Get the input to this writer.
  vtkPolyData* GetInput();
  vtkPolyData* GetInput(int port);

  // Description:
  // Specify file name of vtk polygon data file to write.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

 protected:
  vtkAmeletHDFWriter();
  ~vtkAmeletHDFWriter() {}
  void WriteData();

  virtual int FillInputPortInformation(int port, vtkInformation *info);

  /** 
   * Write VTK mesh into AH5 mesh file
   * 
   * @param pts 
   * @param polys 
   * @param strips 
   * 
   * @return 
   */
  char WriteAH5(vtkPoints *pts, vtkCellArray *polys, vtkCellArray *strips);

  /** 
   * Write VTK mesh into given unstructured mesh
   * 
   * @param mesh 
   * @param pts 
   * @param polys 
   * @param strips 
   * 
   * @return 
   */
  char WriteAH5(
      AH5_umesh_t *mesh, vtkPoints *pts, vtkCellArray *polys, vtkCellArray *strips);

  /** 
   * Return the mesh name
   */
  std::string GetMeshName();

  char *FileName;  // The output file name (.h5 extension)

 private:
  vtkAmeletHDFWriter(const vtkAmeletHDFWriter&);  // Not implemented.
  void operator=(const vtkAmeletHDFWriter&);  // Not implemented.
};

#endif /* VTKAMELETHDFWRITER_H */
