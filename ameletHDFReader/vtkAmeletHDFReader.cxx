#include <vtkAmeletHDFReader.h>





#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()


#define TRUE            1
#define FALSE           0

//using namespace vtkstd;

// Standard VTK Macros for vtkObject derived Classes
vtkCxxRevisionMacro(vtkAmeletHDFReader, "$Revision: 0.2.0$");
vtkStandardNewMacro(vtkAmeletHDFReader);


vtkInformationKeyMacro(vtkAmeletHDFReader, IS_UMESH, Integer);
vtkInformationKeyMacro(vtkAmeletHDFReader, IS_SMESH, Integer);
vtkInformationKeyMacro(vtkAmeletHDFReader, IS_DATAONMESH, Integer);
vtkInformationKeyMacro(vtkAmeletHDFReader, IS_DATA, Integer);
vtkInformationKeyMacro(vtkAmeletHDFReader, POINTS, ObjectBase);
vtkInformationKeyMacro(vtkAmeletHDFReader, POINT_DATA, ObjectBase);

//=============================================================================
// A convenience function that gets a block from a multiblock data set,
// performing allocation if necessary.
static vtkUnstructuredGrid *AllocateGetBlock(vtkMultiBlockDataSet *blocks,
                                             unsigned int blockno,
                                             vtkInformationIntegerKey *typeKey)
{
    if (blockno > 1000)
    {
        vtkGenericWarningMacro(<< "Unexpected block number: " << blockno);
        blockno = 0;
    }

  if (blocks->GetNumberOfBlocks() <= blockno)
  {
      blocks->SetNumberOfBlocks(blockno+1);
  }

  vtkUnstructuredGrid *grid
    = vtkUnstructuredGrid::SafeDownCast(blocks->GetBlock(blockno));
  if (!grid)
  {
      grid = vtkUnstructuredGrid::New();
      blocks->SetBlock(blockno, grid);
      blocks->GetMetaData(blockno)->Set(typeKey, 1);
      grid->Delete();     // Not really deleted.
  }

  return grid;
}
static vtkStructuredGrid *SAllocateGetBlock(vtkMultiBlockDataSet *blocks,
                                             unsigned int blockno,
                                             vtkInformationIntegerKey *typeKey)
{
    if (blockno > 1000)
    {
        vtkGenericWarningMacro(<< "Unexpected block number: " << blockno);
        blockno = 0;
    }

  if (blocks->GetNumberOfBlocks() <= blockno)
  {
      blocks->SetNumberOfBlocks(blockno+1);
  }

  vtkStructuredGrid *grid
    = vtkStructuredGrid::SafeDownCast(blocks->GetBlock(blockno));
  if (!grid)
  {
      grid = vtkStructuredGrid::New();
      blocks->SetBlock(blockno, grid);
      blocks->GetMetaData(blockno)->Set(typeKey, 1);
      grid->Delete();     // Not really deleted.
  }

  return grid;
}


// -----------------------------------------------------------------------------
// test if the file is readable

bool is_readable( const std:: string & file)
{
  std::ifstream ameletfile( file.c_str());
  return !ameletfile.fail();
}

// -----------------------------------------------------------------------------
//
int getAmeletHDFAttribut(hid_t loc_id, const char *attr_name, const char *attr_value)
{
    herr_t   status;
    hid_t    attr_id, attr_type, space, memtype;
    size_t sdim;
    hsize_t     dims[1]={1} ;
    int ndims, i, ret_val;
    char **rdata;

    ret_val = 0;
    attr_id = H5Aopen(loc_id,attr_name,H5P_DEFAULT);

    if (int(attr_id) >= 0)
    {
        attr_type = H5Aget_type(attr_id);
        sdim = H5Tget_size(attr_type);
        sdim++; /* Make room for null terminator */
        space = H5Aget_space (attr_id);
        ndims = H5Sget_simple_extent_dims (space, dims, NULL);
        rdata = (char **) malloc (dims[0] * sizeof (char *));
        rdata[0] = (char *) malloc (dims[0] * sdim * sizeof (char));
        /*
         * Set the rest of the pointers to rows to the correct addresses.
         */
        for (i=1; i<dims[0]; i++)
	    rdata[i] = rdata[0] + i * sdim;
        memtype = H5Tcopy (H5T_C_S1);
        status = H5Tset_size (memtype, sdim);
        status = H5Aread(attr_id, memtype, rdata[0]);
        // Check the value of FORMAT attribute.
        if(strcmp(rdata[0],attr_value)==0) ret_val=1;
        H5Aclose(attr_id);
        free(rdata);
    }
    return ret_val;

}
// -----------------------------------------------------------------------------
//

int getAmeletHDFDataType(char* fileName)
{
	int fileType=0;
	hid_t    file_id;

	file_id = H5Fopen (fileName, H5F_ACC_RDONLY, H5P_DEFAULT);

	// test if arrayset group is present


	if(H5Lexists(file_id,"/floatingType",H5P_DEFAULT)!=FALSE)
	{
		// Is a data on mesh ?
		if (H5Lexists(file_id,"/mesh",H5P_DEFAULT)!=FALSE) fileType = 1; // data on mesh
		else fileType = 2; // data
	}
	else
	{
		if (H5Lexists(file_id,"/mesh",H5P_DEFAULT)!=FALSE) fileType = 3; // mesh

	}
	H5Fclose(file_id);
	
    return fileType;

}
//----------------------------------------------------------------------------
class WithinTolerance: public std::binary_function<double, double, bool>
{
public:
  WithinTolerance(double tol) { this->tolerance = tol; }
  double tolerance;
  //
    result_type operator()(first_argument_type a, second_argument_type b) const
    {
      bool result = (fabs(a-b)<=(a*this->tolerance));
      return (result_type)result;
    }
};
// -----------------------------------------------------------------------------
//
int vtkAmeletHDFReader::ConvertDataToDataOnMesh(hid_t file_id, int nb_dims, int timedim, int componentdim,
		                                        int xdim, int ydim, int zdim, AH5_vector_t *dims, vtkMultiBlockDataSet *output)
{

    char path2[AH5_ABSOLUTE_PATH_LENGTH];
    AH5_children_t children;
    hsize_t i, invalid_nb = -1;
    char invalid = AH5_FALSE;
    int nb_axis=0;
    vtkAmeletHDFDataReader ahdfdata;
    char* entryPoint = NULL;
    commonTools tools;

    AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPoint);

	vtkStructuredGrid *grid = SAllocateGetBlock(output, 0,IS_DATAONMESH());
	ahdfdata.createMeshFromDimsData(file_id, grid);
	//create dataname
	int nbdataarray=1;
	for (i=0;i<nb_dims;i++)
		{
		if(i != xdim)
			if(i != ydim)
				if(i != zdim)
					if(i!=componentdim)
						if(i!=timedim)
							nbdataarray = nbdataarray*dims[i].nb_values;
		}

	int   datanameoffset[nbdataarray][nb_dims];
	vtkstd::string dataname[nbdataarray];


	for (i=0;i<nbdataarray;i++){
			dataname[i]=AH5_get_name_from_path(entryPoint);
			for(int j=0;j<nb_dims;j++)
				datanameoffset[i][j]=0;
	}
	int temp = nbdataarray;

	for( i=0;i<nb_dims;i++)
	{
		if((nb_dims-i-1) != xdim)
			if((nb_dims-i-1) != ydim)
				if((nb_dims-i-1) != zdim)
					if((nb_dims-i-1) != componentdim)
						if((nb_dims-i-1) != timedim)
						{
							temp = (int)temp/dims[nb_dims-i-1].nb_values;
							int j = 0 ;
							if(temp==0) j=nbdataarray;
							while(j<nbdataarray)
							{
								if(dims[nb_dims-i-1].type_class==H5T_COMPOUND)
									for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
									{
										std::ostringstream buf;
										buf<<"_"<<dims[nb_dims-i-1].values.c[l].re
												<<"_j"<<dims[nb_dims-i-1].values.c[l].im;
										for (int k=0;k<temp;k++)
										{
											vtkstd::string label="";
											for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
												if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
													label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
											dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
											datanameoffset[j+k][i]=l;
										}
										j=j+temp;

									}
								else if(dims[nb_dims-i-1].type_class==H5T_FLOAT)
									for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
									{
										std::ostringstream buf;
										buf<<"_"<<dims[nb_dims-i-1].values.f[l];
										for (int k=0;k<temp;k++)
										{
											vtkstd::string label="";
											for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
												if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
													label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
											dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
											datanameoffset[j+k][i]=l;
										}
										j=j+temp;
									}
								else if(dims[nb_dims-i-1].type_class==H5T_INTEGER)
									for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
									{
										std::ostringstream buf;
										buf<<"_"<<dims[nb_dims-i-1].values.i[l];
										for (int k=0;k<temp;k++)
										{
											vtkstd::string label="";
											for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
												if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
													label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
											dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
											datanameoffset[j+k][i]=l;
										}
										j=j+temp;
									}
								else if(dims[nb_dims-i-1].type_class==H5T_STRING)
									for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
									{
										std::ostringstream buf;
										buf<<"_"<<dims[nb_dims-i-1].values.s[l];
										for (int k=0;k<temp;k++)
										{
											vtkstd::string label="";
											for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
												if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
													label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
											dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
											datanameoffset[j+k][i]=l;
										}
										j=j+temp;
									}
							}
						}
	}
	hsize_t count[nb_dims];
	int nbelt=1;
	for (i=0;i<nb_dims;i++){
		if(i==xdim){
			count[nb_dims-i-1]=dims[xdim].nb_values;
			nbelt=nbelt*dims[xdim].nb_values;
		}
		else if(i==ydim){
			count[nb_dims-i-1]=dims[ydim].nb_values;
			nbelt=nbelt*dims[ydim].nb_values;
		}
		else if(i==zdim){
			count[nb_dims-i-1]=dims[zdim].nb_values;
			nbelt=nbelt*dims[zdim].nb_values;
		}
		else
			count[nb_dims-i-1]=1;
	}
	float data_tmp[nbelt];
	hsize_t *offset_tmp;
	offset_tmp = (hsize_t *)malloc(nb_dims * sizeof(hsize_t));

	strcpy(path2, entryPoint);
	strcat(path2, AH5_G_DATA);
	int nb_dim_data;
	H5LTget_dataset_ndims(file_id, path2, &nb_dim_data);
	hsize_t         *data_dims;
	size_t length;
	H5T_class_t     type_class;
	data_dims = (hsize_t *) malloc((nb_dim_data * sizeof(hsize_t)));
	H5LTget_dataset_info(file_id, path2, data_dims, &type_class, &length);

	for(i=0;i<nbdataarray;i++)
	{
		vtkFloatArray *floatscalar = vtkFloatArray::New();
		floatscalar->SetName(dataname[i].c_str());

		if(nb_dims>1)
		{
			if(timedim>-1)
			{
				int actualtimestep = this->ActualTimeStep;

				this->TimeStepMode = true;

				if(componentdim>-1)
				{
					if(dims[componentdim].nb_values<3)
						floatscalar->SetNumberOfComponents(dims[componentdim].nb_values);
					else
						floatscalar->SetNumberOfComponents(3);
					for(int j2=0;j2<dims[componentdim].nb_values;j2++)
					{
						if (j2<3){
						for (int ii=0;ii<nb_dims;ii++){
							if(ii!=nb_dims-timedim-1)
								offset_tmp[ii]=datanameoffset[i][ii];
							else
								offset_tmp[ii]=actualtimestep;
							if(ii==nb_dims-componentdim-1)
								offset_tmp[ii]=j2;
						}
						if(type_class == H5T_COMPOUND)
							tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
						else if(type_class == H5T_FLOAT)
							tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
						for(int k=0;k<nbelt;k++)
							floatscalar->InsertComponent(k,j2,data_tmp[k]);
						}
					}

				}
				else
				{
					for (int ii=0;ii<nb_dims;ii++){
						if(ii!=nb_dims-timedim-1)
							offset_tmp[ii]=datanameoffset[i][ii];
						else
							offset_tmp[ii]=actualtimestep;
					}
					if(type_class == H5T_COMPOUND)
						tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
					else if(type_class == H5T_FLOAT)
						tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
					for(int k=0;k<nbelt;k++)
						floatscalar->InsertTuple1(k,data_tmp[k]);
				}
				grid->GetPointData()->AddArray(floatscalar);

				double timevalue = (double)dims[timedim].values.f[actualtimestep];
				floatscalar->Delete();
				output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),timevalue);

			}
			else if(componentdim>-1)
			{
				if(dims[componentdim].nb_values<3)
					floatscalar->SetNumberOfComponents(dims[componentdim].nb_values);
				else
					floatscalar->SetNumberOfComponents(3);
				for(int j=0;j<dims[componentdim].nb_values;j++){
					if(j<3)
					{
						for (int ii=0;ii<nb_dims;ii++){
							if(ii!=nb_dims-componentdim-1)
								offset_tmp[ii]=datanameoffset[i][ii];
							else
								offset_tmp[ii]=j;

						}
						if(type_class == H5T_COMPOUND)
							tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
						else if(type_class == H5T_FLOAT)
							tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
						for (int k=0;k<nbelt;k++)
							floatscalar->InsertComponent(k,j,data_tmp[k]);
					}
				}

				grid->GetPointData()->AddArray(floatscalar);
				floatscalar->Delete();
			}
			else
			{
				for (int ii=0;ii<nb_dims;ii++) offset_tmp[ii]=datanameoffset[i][ii];
				if(type_class == H5T_COMPOUND)
					tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
				else if(type_class == H5T_FLOAT)
					tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
				for (int k=0;k<nbelt;k++)
					floatscalar->InsertTuple1(k,data_tmp[k]);
				grid->GetPointData()->AddArray(floatscalar);
				floatscalar->Delete();
			}

		}
		else
		{
			for (int ii=0;ii<nb_dims;ii++) offset_tmp[ii]=datanameoffset[i][ii];
			if(type_class == H5T_COMPOUND)
				tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
			else if(type_class == H5T_FLOAT)
				tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
			for (int k=0;k<nbelt;k++)
				floatscalar->InsertTuple1(k,data_tmp[k]);

			grid->GetPointData()->AddArray(floatscalar);
			floatscalar->Delete();
		}
	}
}
// -----------------------------------------------------------------------------
//
int vtkAmeletHDFReader::ReadDataOnMesh(hid_t file_id, vtkMultiBlockDataSet *output)
{
	vtkUnstructuredGrid *grid = AllocateGetBlock(output, 0,IS_DATAONMESH());
	vtkAmeletHDFMeshReader ahdfmesh;
	commonTools tools;

	char* entryPoint;
	char path2[AH5_ABSOLUTE_PATH_LENGTH];
    int nb_dims;
    int timedim, componentdim, meshdim, xdim, ydim, zdim;
	int nbdataarray=1;
	//hsize_t nb_dims=0;
	hsize_t nb_dims_data=0;
	hsize_t i, invalid_nb = -1;
    int dims_param[6];
	char invalid = AH5_FALSE;




	AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPoint);



    nb_dims = tools.readNbDims(file_id);

    AH5_vector_t *dims;
    dims = (AH5_vector_t *) malloc((size_t) nb_dims * sizeof(AH5_vector_t));
    tools.readDims(file_id, dims_param, dims);
    timedim = dims_param[0];
    componentdim = dims_param[1];
    meshdim = dims_param[2];
    xdim = dims_param[3];
    ydim = dims_param[4];
    zdim = dims_param[5];



	for (int i=0;i<nb_dims;i++)
	{
		if(i!=meshdim)
			if(i!=componentdim)
				if(i!=timedim)
					nbdataarray = nbdataarray*dims[i].nb_values;
	}
	int nbelt = 0;

    nbelt = ahdfmesh.readMeshGroup(file_id,dims[meshdim].values.s[0],grid);
    char *type;
    AH5_read_str_attr(file_id,dims[meshdim].values.s[0] , AH5_A_TYPE, &type);

	// set data name
	int   datanameoffset[nbdataarray][nb_dims];
	vtkstd::string dataname[nbdataarray];


	for (int i=0;i<nbdataarray;i++){
		dataname[i]=AH5_get_name_from_path(entryPoint);
		for(int j=0;j<nb_dims;j++)
			datanameoffset[i][j]=0;
	}

	int temp = nbdataarray;

	for(int i=0;i<nb_dims;i++)
	{
		if((nb_dims-i-1) != meshdim)
			if((nb_dims-i-1) != componentdim)
				if((nb_dims-i-1) != timedim)
				{
					temp = (int)temp/dims[nb_dims-i-1].nb_values;
					int j = 0 ;
					if(temp==0) j=nbdataarray;
					while(j<nbdataarray)
					{
						if(dims[nb_dims-i-1].type_class==H5T_COMPOUND)
							for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
							{
								std::ostringstream buf;
								buf<<"_"<<dims[nb_dims-i-1].values.c[l].re
										<<"_j"<<dims[nb_dims-i-1].values.c[l].im;
								for (int k=0;k<temp;k++)
								{
									vtkstd::string label="";
									for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
										if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
											label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
									dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
									datanameoffset[j+k][i]=l;
								}
								j=j+temp;

							}
						else if(dims[nb_dims-i-1].type_class==H5T_FLOAT)
							for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
							{
								std::ostringstream buf;
								buf<<"_"<<dims[nb_dims-i-1].values.f[l];
								for (int k=0;k<temp;k++)
								{
									vtkstd::string label="";
									for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
										if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
											label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
									dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
									datanameoffset[j+k][i]=l;
								}
								j=j+temp;
							}
						else if(dims[nb_dims-i-1].type_class==H5T_INTEGER)
							for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
							{
								std::ostringstream buf;
								buf<<"_"<<dims[nb_dims-i-1].values.i[l];
								for (int k=0;k<temp;k++)
								{
									vtkstd::string label="";
									for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
										if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
											label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
									dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
									datanameoffset[j+k][i]=l;
								}
								j=j+temp;
							}
						else if(dims[nb_dims-i-1].type_class==H5T_STRING)
							for (int l=0;l<dims[nb_dims-i-1].nb_values;l++)
							{
								std::ostringstream buf;
						        buf<<"_"<<dims[nb_dims-i-1].values.s[l];
								for (int k=0;k<temp;k++)
								{
									vtkstd::string label="";
									for(int ii=0;ii<dims[nb_dims-i-1].opt_attrs.nb_instances;ii++)
										if(strcmp(dims[nb_dims-i-1].opt_attrs.instances[ii].name,"label")==0)
											label=label+dims[nb_dims-i-1].opt_attrs.instances[ii].value.s;
									dataname[j+k]=dataname[j+k]+"_"+label+buf.str();
									datanameoffset[j+k][i]=l;
								}
								j=j+temp;
							}
					}
				}
	}

	// get type class of data
	strcpy(path2, entryPoint);
	strcat(path2, AH5_G_DATA);
    int nb_dim_data;
	H5LTget_dataset_ndims(file_id, path2, &nb_dim_data);
	hsize_t         *data_dims;
	size_t length;
	H5T_class_t     type_class;
	data_dims = (hsize_t *) malloc((nb_dim_data * sizeof(hsize_t)));
	H5LTget_dataset_info(file_id, path2, data_dims, &type_class, &length);


    hsize_t count[nb_dims];
    for (int i=0;i<nb_dims;i++){
    	if(i!=meshdim)
    		count[nb_dims-i-1]=1;
    	else
    		count[nb_dims-i-1]=nbelt;
    }
    float data_tmp[nbelt];
    hsize_t *offset_tmp;
    offset_tmp = (hsize_t *)malloc(nb_dims * sizeof(hsize_t));


    for(int i=0;i<nbdataarray;i++)
	{
		vtkFloatArray *floatscalar = vtkFloatArray::New();
		floatscalar->SetName(dataname[i].c_str());

		if(nb_dims>1)
		{
			if(timedim>-1)
			{
				int actualtimestep = this->ActualTimeStep;

				this->TimeStepMode = true;

				if(componentdim>-1)
				{
				    if(dims[componentdim].nb_values<3)
					    floatscalar->SetNumberOfComponents(dims[componentdim].nb_values);
					else
					    floatscalar->SetNumberOfComponents(3);
                    for(int j2=0;j2<dims[componentdim].nb_values;j2++)
					{
                    	if (j2<3){
                    	for (int ii=0;ii<nb_dims;ii++){
                            if(ii!=nb_dims-timedim-1)
                    		    offset_tmp[ii]=datanameoffset[i][ii];
                    		else
                    		    offset_tmp[ii]=actualtimestep;
                    		if(ii==nb_dims-componentdim-1)
                    	        offset_tmp[ii]=j2;
                    	}
                    	if(type_class == H5T_COMPOUND)
                    	    tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
                    	else if(type_class == H5T_FLOAT)
                    	    tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
					    for(int k=0;k<nbelt;k++)
					    	floatscalar->InsertComponent(k,j2,data_tmp[k]);
                    	}
                    }

				}
				else
				{
                	for (int ii=0;ii<nb_dims;ii++){
                        if(ii!=nb_dims-timedim-1)
                		    offset_tmp[ii]=datanameoffset[i][ii];
                		else
                		    offset_tmp[ii]=actualtimestep;
                	}
                	if(type_class == H5T_COMPOUND)
                	    tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
                	else if(type_class == H5T_FLOAT)
                	    tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
				    for(int k=0;k<nbelt;k++)
				    	floatscalar->InsertTuple1(k,data_tmp[k]);
				}
				if(grid->GetCell(0)->GetCellType()==VTK_VERTEX)
				    grid->GetPointData()->AddArray(floatscalar);
				else
					grid->GetCellData()->AddArray(floatscalar);

				double timevalue = (double)dims[timedim].values.f[actualtimestep];
				floatscalar->Delete();
				output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(),timevalue);

			}
			else if(componentdim>-1)
			{
				if(dims[componentdim].nb_values<3)
				    floatscalar->SetNumberOfComponents(dims[componentdim].nb_values);
				else
				    floatscalar->SetNumberOfComponents(3);
				for(int j=0;j<dims[componentdim].nb_values;j++){
					if(j<3)
					{
						for (int ii=0;ii<nb_dims;ii++){
							if(ii!=nb_dims-componentdim-1)
							    offset_tmp[ii]=datanameoffset[i][ii];
							else
								offset_tmp[ii]=j;

						}
						if(type_class == H5T_COMPOUND)
						    tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
						else if(type_class == H5T_FLOAT)
							tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
						for (int k=0;k<nbelt;k++)
						    floatscalar->InsertComponent(k,j,data_tmp[k]);
					}
				}



				if(strcmp(type,"node")==0)//grid->GetCell(0)->GetCellType()==VTK_VERTEX)
				{

					grid->GetPointData()->AddArray(floatscalar);
				}
				else
				{

					grid->GetCellData()->AddArray(floatscalar);
				}
				floatscalar->Delete();
			}
			else
			{
				for (int ii=0;ii<nb_dims;ii++) offset_tmp[ii]=datanameoffset[i][ii];
				if(type_class == H5T_COMPOUND)
				    tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
				else if(type_class == H5T_FLOAT)
					tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
				for (int k=0;k<nbelt;k++)
				    floatscalar->InsertTuple1(k,data_tmp[k]);
				if(strcmp(type,"node")==0)//if(grid->GetCell(0)->GetCellType()==VTK_VERTEX)
				{
					grid->GetPointData()->AddArray(floatscalar);
				}
				else
				{
					grid->GetCellData()->AddArray(floatscalar);
				}
				floatscalar->Delete();
			}

		}
		else
		{
			for (int ii=0;ii<nb_dims;ii++) offset_tmp[ii]=datanameoffset[i][ii];
			if(type_class == H5T_COMPOUND)
			    tools.readCpxDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
			else if(type_class == H5T_FLOAT)
				tools.readFltDatasetSlice( file_id, path2, offset_tmp, count, nb_dims, data_tmp);
			for (int k=0;k<nbelt;k++)
			    floatscalar->InsertTuple1(k,data_tmp[k]);

			if(strcmp(type,"node")==0)//grid->GetCell(0)->GetCellType()==VTK_VERTEX)
			{

				grid->GetPointData()->AddArray(floatscalar);
			}
			else
			{

				grid->GetCellData()->AddArray(floatscalar);
			}
			floatscalar->Delete();
		}
	}

    for (i = 0; i < nb_dims; i++)
        AH5_free_ft_vector(dims + i);
    free(dims);
    free(entryPoint);
    free(offset_tmp);


    return 1;
}

// -----------------------------------------------------------------------------
//
vtkAmeletHDFReader::vtkAmeletHDFReader()
{
	this->SetNumberOfInputPorts(0);

	this->FileName = NULL;
	this->TimeStep = 0;
	this->ActualTimeStep = 0;
	this->NumberOfTimeSteps = 0;
	this->TimeStepRange[0] = 0;
	this->TimeStepRange[1] = 0;
	this->TimeStepMode = false;

}

// -----------------------------------------------------------------------------
//
vtkAmeletHDFReader::~vtkAmeletHDFReader()
{
	this->SetFileName(NULL);

}

// -----------------------------------------------------------------------------
//
void vtkAmeletHDFReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
  os << indent << "FileName: " << this->FileName << endl;
}

// -----------------------------------------------------------------------------
//
int vtkAmeletHDFReader::FillInputPortInformation(int, vtkInformation *info)
{
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSetAlgorithm");
  return 1;
}

// -----------------------------------------------------------------------------
//
int vtkAmeletHDFReader::CanReadFile(const char *filename)
{

  hid_t    file_id;
  int ret_val;

  char *amelet;
  const char *path=".";
  ret_val = 0;
  if ( is_readable(filename))
  {
     
      file_id = H5Fopen (filename, H5F_ACC_RDONLY, H5P_DEFAULT);
      // Check FORMAT and AMELETHDF_FORMAT_VERSION attributes values
      AH5_read_str_attr(file_id, path, AH5_FILE_A_FORMAT, &amelet);
      if(strcmp(amelet,"AMELETHDF")==0)
    	  ret_val=1;
      H5Fclose(file_id);
      free(amelet);

      return ret_val;
  }
  else 
  {
    std::cout<<"file doesn't exist !!!"<<std::endl;
    return 0;
  }


}

// -----------------------------------------------------------------------------
//
int vtkAmeletHDFReader::RequestData( vtkInformation *request,
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector)
{
  int dataType = 0;
  hid_t file_id, loc_id;
  herr_t status;
  char* entryPoint = NULL;

  // get the info objects
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::GetData(outInfo);
  this->UpdateProgress(0.);

  if (!this->FileName)
     {
     vtkErrorMacro("No filename specified.");
     return 0;
     }
  if( !is_readable(this->FileName))
  {
      vtkErrorMacro("The file does not exist.");
      return 0;
  }
  // Test data type (unstructured mesh, structured mesh, data on mesh, or data ?)
  dataType = getAmeletHDFDataType(this->FileName);
  if (dataType==0)
      {
	  vtkErrorMacro("This is not a mesh, data or data on mesh ameletHDF file.");
	  return 0;
	  }

  file_id = H5Fopen (this->FileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPoint);

  if(dataType==1)
  {
	  //cout<<"data on mesh conversion"<<endl;
      int tsLength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
      double* steps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

	  if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) && this->TimeStepMode )
	  {
		  double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
		  int nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
		  double* steps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
          int cnt = 0;
          while(cnt < tsLength-1 && steps[cnt] < requestedTimeValue)
          {
        	  cnt++;
          }
          this->ActualTimeStep = cnt;
		  output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), TimeStepValues[this->ActualTimeStep]);
	  }
	  this->ReadDataOnMesh(file_id,output);
      this->UpdateProgress(1);
   	  //loc_id = H5Gopen(file_id,"/floatingType",H5P_DEFAULT);


  }
  else if(dataType==2)
  {
      //cout<<"data conversion"<<endl;
	  commonTools tools;
	  vtkAmeletHDFDataReader ahdfdata;
	  char path2[AH5_ABSOLUTE_PATH_LENGTH];
	  AH5_children_t children;
	  hsize_t i, invalid_nb = -1;
	  char invalid = AH5_FALSE;
	  int nb_axis=0;
      int dims_param[6];
      int nb_dims;
	  int timedim, componentdim, meshdim, xdim, ydim, zdim;
	  nb_dims = tools.readNbDims(file_id);
	  AH5_vector_t *dims;
	  dims = (AH5_vector_t *) malloc((size_t) nb_dims * sizeof(AH5_vector_t));
      tools.readDims(file_id, dims_param, dims);
      timedim = dims_param[0];
      componentdim = dims_param[1];
      meshdim = dims_param[2];
      xdim = dims_param[3];
      ydim = dims_param[4];
      zdim = dims_param[5];
      if(xdim>-1) nb_axis=nb_axis+1;
      if(ydim>-1) nb_axis=nb_axis+1;
      if(zdim>-1) nb_axis=nb_axis+1;

	  if(nb_axis>2)
	  {
	  		if(timedim>-1){
				int tsLength = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
				double* steps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());

				if(outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()) && this->TimeStepMode )
				{
				  double requestedTimeValue = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
				  int nSteps = outInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
				  double* steps = outInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
					int cnt = 0;
					while(cnt < tsLength-1 && steps[cnt] < requestedTimeValue)
					{
					  cnt++;
					}
					this->ActualTimeStep = cnt;
				  output->GetInformation()->Set(vtkDataObject::DATA_TIME_STEP(), TimeStepValues[this->ActualTimeStep]);
				}
	  		}
	  	    this->ConvertDataToDataOnMesh(file_id, nb_dims, timedim, componentdim, xdim, ydim, zdim, dims, output);

	  }
	  else
	  {
	      cout<<"Data read"<<endl;
          vtkTable *table = vtkTable::New();
          output->SetBlock(0,table);
          ahdfdata.readData(file_id,table);
   	  }
	  for (i = 0; i < nb_dims; i++)
	      AH5_free_ft_vector(dims + i);
	  free(dims);
      
  }
  else if(dataType==3) //mesh
  {

	  vtkAmeletHDFMeshReader ahdfmesh;
	  AH5_mesh_t mesh;
	  //read mesh "/mesh"
	  AH5_read_mesh(file_id, &mesh);

      int block_id=0;
      for(int i=0;i < mesh.nb_groups ; i++)
      {
    	  for (int j=0; j<mesh.groups[i].nb_msh_instances;j++)
    	  {

    		  if(mesh.groups[i].msh_instances[j].type == MSH_UNSTRUCTURED)
    		  {

    			  vtkUnstructuredGrid *ugrid = AllocateGetBlock(output, block_id, IS_UMESH());
    		      ahdfmesh.readUmesh(mesh.groups[i].msh_instances[j].data.unstructured,ugrid);
    		  }
    		  else if(mesh.groups[i].msh_instances[j].type == MSH_STRUCTURED)
    		  {

    		      vtkUnstructuredGrid *sgrid = AllocateGetBlock(output, block_id ,IS_SMESH());
    		      ahdfmesh.readSmesh(mesh.groups[i].msh_instances[j].data.structured,sgrid);

    		  }
    		  double prog = this->GetProgress()+(i+1)*j/(mesh.nb_groups);
    		  this->UpdateProgress(prog);
    		  block_id=block_id+1;
    	  }

      }

      AH5_free_mesh(&mesh);

  }
  else if(dataType>3 || dataType<1)
  {
	  vtkErrorMacro("This is not an ameletHDF data or mesh file .");
	  return 0;
  }


  H5Fclose(file_id);
  return 1;
}
// -----------------------------------------------------------------------------
//
int vtkAmeletHDFReader::RequestInformation(vtkInformation *vtkNotUsed(request),
                                           vtkInformationVector **vtkNotUsed(inputVector),
                                           vtkInformationVector *outputVector)
{

	vtkInformation* outInfo = outputVector->GetInformationObject(0);
	int dataType=0;
	char* entryPoint = NULL;
	commonTools tools;

	vtkDebugMacro("RequestInformation");
	if (!this->FileName)
	    {
	    vtkErrorMacro("No filename specified.");
	    return 0;
	    }
	// Test data type (unstructured mesh, structured mesh, data on mesh, or data ?)
        if( !is_readable(this->FileName))
        {
          vtkErrorMacro("The file does not exist.");
          return 0;
        }
    dataType = getAmeletHDFDataType(this->FileName);
    if(dataType==1 or dataType==2)
    {
    	//time series ?
    	hid_t file_id;
    	file_id=H5Fopen(this->FileName,H5F_ACC_RDONLY, H5P_DEFAULT);


        //loop on dims
    	TimeStepMode = false;
    	const char * path=".";
    	char path2[AH5_ABSOLUTE_PATH_LENGTH];
	    //AH5_children_t children;


	    hsize_t i = -1;
	    int dims_param[6];
        int nb_dims;
	    int timedim, componentdim, meshdim, xdim, ydim, zdim;
	    nb_dims = tools.readNbDims(file_id);
	    AH5_vector_t *dims;
	    dims = (AH5_vector_t *) malloc((size_t) nb_dims * sizeof(AH5_vector_t));
        tools.readDims(file_id, dims_param, dims);
        timedim = dims_param[0];
        componentdim = dims_param[1];
        meshdim = dims_param[2];
        xdim = dims_param[3];
        ydim = dims_param[4];
        zdim = dims_param[5];

        for(int i=0;i<nb_dims;i++)
        {

            for(int j=0;j<dims[i].opt_attrs.nb_instances;j++)
            	if(dims[i].opt_attrs.instances[j].type==H5T_STRING)
            	{

            		if((strcmp(dims[i].opt_attrs.instances[j].value.s,"time")==0) ||
            				(strcmp(dims[i].opt_attrs.instances[j].value.s,"frequency")==0))
            		{
            			outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
            			outInfo->Remove(vtkStreamingDemandDrivenPipeline::TIME_RANGE());
                    	this->TimeStepMode = true;
                        this->TimeStepValues = new double[dims[i].nb_values];
                    	for(int k=0;k<dims[i].nb_values;k++)
                    		this->TimeStepValues[k] = dims[i].values.f[k];
                    	this->TimeStepRange[0] = 0;
                    	this->TimeStepRange[1] = dims[i].nb_values - 1;
    		            outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_STEPS(),
    		            		&this->TimeStepValues[0], dims[i].nb_values);
    		            double timeRange[2];
    		            timeRange[0]=this->TimeStepValues[0];
    		            timeRange[1]=this->TimeStepValues[dims[i].nb_values-1];
    		            outInfo->Set(vtkStreamingDemandDrivenPipeline::TIME_RANGE(), timeRange,2);
            		}
            	}
        }

        H5Fclose(file_id);

        for (i = 0; i < nb_dims; i++)
            AH5_free_ft_vector(dims + i);
        free(dims);
    }


    if (dataType==0) return 0;

	return 1;
}


