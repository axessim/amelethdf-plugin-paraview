/*
 * tools.cxx

 *
 *  Created on: 5 août 2015
 *      Author: didier
 */
#include "tools.h"



using namespace vtkstd;

int commonTools::getEntryPoint(hid_t file_id, std::string* entryPoint)
{

	int success = AH5_FALSE;
	char *entryPt = NULL;
	vtkstd::string path2;
	AH5_children_t children;
	hsize_t i, invalid_nb = -1;
	char invalid = AH5_FALSE;
	int nb_dims=0;




	if(!AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPt))
	{
			// test if floatingType group is present

			if(H5Lexists(file_id,"/floatingType",H5P_DEFAULT)!=AH5_FALSE)
			{
	    		path2="/floatingType";
				children = AH5_read_children_name(file_id, path2.c_str());
				path2=path2+ children.childnames[0];
				success = AH5_TRUE;

			}
			else
			{
				if (H5Lexists(file_id,"/mesh",H5P_DEFAULT)!=AH5_FALSE)
				{
				    path2="/mesh";
				    children = AH5_read_children_name(file_id, path2.c_str());
				    path2=path2+children.childnames[0];
				    success = AH5_TRUE;
				}
				//else path2=NULL;
			}

			*entryPoint=std::string (path2);
			for (i = 0; i < children.nb_children; i++)
				{
					free(children.childnames[i]);
				}
			free(children.childnames);

	}
	else
	{
		*entryPoint=std::string (entryPt);
		if	(strncmp(entryPt,"/floatingType", strlen("/floatingType"))==0)
    		success = AH5_TRUE;
        else if (strncmp(entryPt, "/mesh", strlen("/mesh"))==0)
    		success = AH5_TRUE;
        else
        	success = AH5_FALSE;
	}


	return success;
}
float commonTools::abs_complex(AH5_complex_t complex)
{
	float module;
	module = (complex.re*complex.re)+(complex.im*complex.im);
	module = sqrt(module);
	return module;
}
char commonTools::readFltDatasetSlice(hid_t file_id, const char *path,  hsize_t *offset,  hsize_t *count,
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



char commonTools::readCpxDatasetSlice(hid_t file_id, const char *path,  hsize_t *offset,  hsize_t *count,
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

int commonTools::readDims(hid_t file_id,  int *dims_param, AH5_vector_t *dims)
{

	char success = AH5_FALSE;
	std::string entryPoint;
	vtkstd::string path2;
	AH5_children_t children;
	hsize_t i, invalid_nb = -1;
	char invalid = AH5_FALSE;
	int nb_axis=0;

	int nb_dims=0;
	int timedim=-1;
	int componentdim=-1;
	int meshdim=-1;
	int xdim=-1;
	int ydim=-1;
	int zdim=-1;


	//AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPoint);
	getEntryPoint(file_id, &entryPoint);
	path2=entryPoint;
	path2=path2+AH5_G_DS;
	children = AH5_read_children_name(file_id, path2.c_str());
	nb_dims = children.nb_children;
	//dims = (AH5_vector_t *) malloc((size_t) children.nb_children * sizeof(AH5_vector_t));
	for (i = 0; i < children.nb_children; i++)
	{
		if (!invalid)
		{
			path2=entryPoint;
			path2=path2+ AH5_G_DS;
			path2=path2+children.childnames[i];
			if(!AH5_read_ft_vector(file_id, path2.c_str(), dims + i))
			{
				invalid_nb = i;
				invalid = AH5_TRUE;
			}
		}
		free(children.childnames[i]);
	}
	free(children.childnames);
	if (invalid)
	{
		for (i = 0; i < invalid_nb; i++)
			  AH5_free_ft_vector(dims + i);
		free(dims);
	}
	else
		success = AH5_TRUE;


	for ( i=0;i<nb_dims;i++)
	{
		for (int j=0;j<dims[i].opt_attrs.nb_instances;j++)
		{
			if(strcmp(dims[i].opt_attrs.instances[j].name,"physicalNature")==0)
			{
				if(strcmp(dims[i].opt_attrs.instances[j].value.s,"length")==0)
					for(int k=0;k<dims[i].opt_attrs.nb_instances;k++){
						if(strcmp(dims[i].opt_attrs.instances[k].name,"label")==0)
							if(strcmp(dims[i].opt_attrs.instances[k].value.s,"Xaxis")==0){
								nb_axis++;
								xdim=i;
							}
							else if(strcmp(dims[i].opt_attrs.instances[k].value.s,"Yaxis")==0){
								nb_axis++;
								ydim=i;
							}
							else if(strcmp(dims[i].opt_attrs.instances[k].value.s,"Zaxis")==0){
								nb_axis++;
								zdim=i;
							}
					}
				else if(strcmp(dims[i].opt_attrs.instances[j].value.s,"time")==0)
					timedim=i;
				else if(strcmp(dims[i].opt_attrs.instances[j].value.s,"frequency")==0)
					timedim=i;
				else if(strcmp(dims[i].opt_attrs.instances[j].value.s,"component")==0)
					componentdim=i;
				else if(strcmp(dims[i].opt_attrs.instances[j].value.s,"meshEntity")==0)
									meshdim=i;
			}
		}
	}
    dims_param[0]= timedim;
	dims_param[1]= componentdim;
	dims_param[2]= meshdim;
	dims_param[3]= xdim;
	dims_param[4]= ydim;
	dims_param[5]= zdim;

    return 1;
}
int commonTools::readNbDims(hid_t file_id)
{
	std::string entryPoint;
	vtkstd::string path2;
	AH5_children_t children;
	hsize_t i, invalid_nb = -1;
	char invalid = AH5_FALSE;
	int nb_dims=0;

	//AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPoint);
	getEntryPoint(file_id, &entryPoint);


	path2=entryPoint;
	path2=path2+ AH5_G_DS;
	children = AH5_read_children_name(file_id, path2.c_str());
	nb_dims = children.nb_children;
	for (i = 0; i < children.nb_children; i++)
	{
		free(children.childnames[i]);
	}
	free(children.childnames);
    return nb_dims;
}

