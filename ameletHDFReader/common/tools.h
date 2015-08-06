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
class commonTools
{
public:
	float abs_complex(AH5_complex_t complex);
	char AH5_read_flt_dataset_slice(hid_t file_id, const char *path,  hsize_t *offset, hsize_t *count, int rank, float *rdata);
	char AH5_read_cpx_dataset_slice(hid_t file_id, const char *path,  hsize_t *offset, hsize_t *count, int rank, float *rdata);

};

#endif /* COMMON_TOOLS_H_ */
