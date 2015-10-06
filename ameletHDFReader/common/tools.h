/*
 * tools.h
 *
 *  Created on: 5 août 2015
 *      Author: didier
 */

#ifndef COMMON_TOOLS_H_
#define COMMON_TOOLS_H_



#include <hdf5.h>
#include "ah5.h"
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include "vtkObjectFactory.h"

class commonTools
{
public:
	int getEntryPoint(hid_t file_id, std::string* entryPoint);
	float abs_complex(AH5_complex_t complex);
	char readFltDatasetSlice(hid_t file_id, const char *path,  hsize_t *offset, hsize_t *count, int rank, float *rdata);
	char readCpxDatasetSlice(hid_t file_id, const char *path,  hsize_t *offset, hsize_t *count, int rank, float *rdata);
	int readDims(hid_t file_id, int *dims_param, AH5_vector_t *dims);
	int readNbDims(hid_t file_id);
};

#endif /* COMMON_TOOLS_H_ */
