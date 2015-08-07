/*
 * vtkAmeletHDFDataReader.h
 *
 *  Created on: 8 janv. 2010
 *      Author: didier
 */

#ifndef VTKAMELETHDFDATAREADER_H_
#define VTKAMELETHDFDATAREADER_H_

#include <hdf5.h>
#include <vtkTable.h>
#include <vtkStructuredGrid.h>
#include <vtkAbstractArray.h>
#include <vtkFloatArray.h>
#include <vtkStringArray.h>
#include <vtkIntArray.h>

#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

#include "ah5.h"
#include "common/tools.h"


class vtkAmeletHDFDataReader
{
public:
	virtual int createMeshFromDimsData(hid_t file_id, vtkStructuredGrid *sgrid);
	virtual int readData(hid_t loc_id, vtkTable *table);
};
#endif /* VTKAMELETHDFDATAREADER_H_ */
