/*
 * tools.cxx

 *
 *  Created on: 5 août 2015
 *      Author: didier
 */
#include "tools.h"


float commonTools::abs_complex(AH5_complex_t complex)
{
	float module;
	module = (complex.re*complex.re)+(complex.im*complex.im);
	module = sqrt(module);
	return module;
}
char commonTools::AH5_read_flt_dataset_slice(hid_t file_id, const char *path,  hsize_t *offset,  hsize_t *count,
		                                       int rank,  float *rdata)
{
    char success = AH5_FALSE;
    hid_t dataset;
    int i, mn;
    hid_t       dataspace;
    hid_t       memspace;
    herr_t      status;
    hsize_t     *dimsm;
    hsize_t     *offset_out;

    mn=1;
    for (i=0;i<rank;i++)
    	mn=mn*count[i];

    dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    dataspace = H5Dget_space (dataset);
    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL,
                                       count, NULL);
    dimsm = (hsize_t *) malloc(rank * sizeof(hsize_t));
    for (i=0;i<rank;i++)
        dimsm[i] = count[i];
    memspace = H5Screate_simple (rank, dimsm, NULL);
    offset_out = (hsize_t *) malloc(rank * sizeof(hsize_t));

    for (i=0;i<rank;i++){
            offset_out[i]=0;
        }
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL,
                                      count, NULL);

    if( H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace,
                          H5P_DEFAULT, rdata)>=0)
        success = AH5_TRUE;
    H5Dclose(dataset);
    H5Sclose (dataspace);
    H5Sclose (memspace);
    return success;
}


// Read 1D complex float dataset
char commonTools::AH5_read_cpx_dataset_slice(hid_t file_id, const char *path,  hsize_t *offset,  hsize_t *count,
        int rank,  float *rdata)
{
    char success = AH5_FALSE;
    hid_t dataset;
    int i, mn;
    hid_t       dataspace;
    hid_t       memspace;
    hid_t       type_id;
    herr_t      status;
    hsize_t     *dimsm;
    hsize_t     *offset_out;
    float       *cplx_data;

    mn=1;
    for (i=0;i<rank;i++)
       	mn=mn*count[i];


    dataset = H5Dopen2(file_id, path, H5P_DEFAULT);
    dataspace = H5Dget_space (dataset);
    status = H5Sselect_hyperslab (dataspace, H5S_SELECT_SET, offset, NULL,
                                       count, NULL);
    dimsm = (hsize_t *) malloc(rank * sizeof(hsize_t));
    for (i=0;i<rank;i++)
        dimsm[i] = count[i];
    memspace = H5Screate_simple (rank, dimsm, NULL);
    offset_out = (hsize_t *) malloc(rank * sizeof(hsize_t));

    for (i=0;i<rank;i++){
            offset_out[i]=0;
        }
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_out, NULL,
                                      count, NULL);
    type_id = AH5_H5Tcreate_cpx_memtype();
    cplx_data = (float *)malloc(mn*2*sizeof(float));
    if( H5Dread (dataset, type_id, memspace, dataspace,
                          H5P_DEFAULT, cplx_data)>=0)
        success = AH5_TRUE;
    AH5_complex_t complex;
    for (i=0;i<mn;i++){
    	complex.re=cplx_data[2*i];
    	complex.im=cplx_data[2*i+1];
    	rdata[i]=abs_complex(complex);
    }

    H5Dclose(dataset);
    H5Sclose (dataspace);
    H5Sclose (memspace);
    free(cplx_data);
    return success;
}



