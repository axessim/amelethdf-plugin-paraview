/*
 * vtkAmeletHDFDataReader.cxx
 *
 *  Created on: 8 janv. 2010
 *      Author: didier
 */

#include "vtkAmeletHDFDataReader.h"


#define TRUE            1
#define FALSE           0

int vtkAmeletHDFDataReader::createMeshFromDimsData(hid_t file_id, vtkStructuredGrid *sgrid)
{
	commonTools tools;
    char path2[AH5_ABSOLUTE_PATH_LENGTH];
    AH5_children_t children;
    AH5_vector_t    *dims;
    int nb_dims;
    hsize_t i, invalid_nb = -1;
    char invalid = AH5_FALSE;
    int nb_axis=0;
    int x_axis_dim=-1;
    int y_axis_dim=-1;
    int z_axis_dim=-1;
    int timedim=-1;
    int componentdim=-1;
    int nbx=1;
    int nby=1;
    int nbz=1;
    vtkPoints *points = vtkPoints::New();
    std::string entryPoint;


    tools.getEntryPoint(file_id, &entryPoint);
    strcpy(path2, entryPoint.c_str());
  	strcat(path2, AH5_G_DS);
  	children = AH5_read_children_name(file_id, path2);
  	nb_dims = children.nb_children;
  	dims = (AH5_vector_t *) malloc((size_t) children.nb_children * sizeof(AH5_vector_t));
  	for (i = 0; i < children.nb_children; i++)
  	{
  	    if (!invalid)
  	    {
  	        strcpy(path2, entryPoint.c_str());
  	        strcat(path2, AH5_G_DS);
  	        strcat(path2, children.childnames[i]);
  	        if(!AH5_read_ft_vector(file_id, path2, dims + i))
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
  					    		x_axis_dim=i;
  					    		nbx=dims[i].nb_values;
  					    	}
  					    	else if(strcmp(dims[i].opt_attrs.instances[k].value.s,"Yaxis")==0)
  					    	{
                                nb_axis++;
                                y_axis_dim=i;
                                nby=dims[i].nb_values;
  					    	}
  					    	else if(strcmp(dims[i].opt_attrs.instances[k].value.s,"Zaxis")==0)
  					    	{
                                nb_axis++;
                                z_axis_dim=i;
                                nbz=dims[i].nb_values;
  					    	}
  					}

  			}
  		}
  	}
    sgrid->SetDimensions(nbx,nby,nbz);

    for (int k=0 ;k<nbz;k++)
       	for (int j=0;j<nby;j++)
       		for(int i=0;i<nbx;i++)
       			if(z_axis_dim==-1)
                    points->InsertNextPoint(dims[x_axis_dim].values.f[i],
                		                    dims[y_axis_dim].values.f[j],
						                    0);
       			else if(y_axis_dim==-1)
                    points->InsertNextPoint(dims[x_axis_dim].values.f[i],
    		                                0,
											dims[z_axis_dim].values.f[k]);
       			else if(x_axis_dim==-1)
                    points->InsertNextPoint(0,
                    		                dims[y_axis_dim].values.f[j],
											dims[z_axis_dim].values.f[k]);
       			else
                    points->InsertNextPoint(dims[x_axis_dim].values.f[i],
                    		                dims[y_axis_dim].values.f[j],
											dims[z_axis_dim].values.f[k]);


    sgrid->SetPoints(points);
    points->Delete();


    						;
    for (i = 0; i < nb_dims; i++)
        AH5_free_ft_vector(dims + i);
    free(dims);
    return 1;
}

int vtkAmeletHDFDataReader::readData(hid_t file_id, vtkTable *table)
{
    vtkFloatArray *array;
    vtkIntArray *iarray;
    commonTools tools;
    std::string entryPoint;
    AH5_ft_t ft;

    tools.getEntryPoint(file_id,&entryPoint);

    AH5_read_floatingtype(file_id, entryPoint.c_str(), &ft);

    char path[AH5_ABSOLUTE_PATH_LENGTH];
    char strtemp[AH5_ABSOLUTE_PATH_LENGTH];
    int nbdataarray=1;
    int max=1;
    

    for(int i=0;i<ft.data.arrayset.nb_dims;i++)
            nbdataarray = nbdataarray*ft.data.arrayset.dims[i].nb_values;
    vtkstd::string dataname;
    dataname= AH5_get_name_from_path(entryPoint.c_str());
    int temp=nbdataarray;

    int offset = 0;

    for(int i=0;i<ft.data.arrayset.nb_dims;i++)
    {
      
        if(ft.data.arrayset.dims[i].type_class==H5T_FLOAT)
        {
            strcpy(strtemp,"dim");
			std::ostringstream buf;
			buf << i;
			std::string buffer = buf.str();
			char *buf2=new char[buffer.size()+1];
			buf2[buffer.size()]=0;
			memcpy(buf2,buffer.c_str(),buffer.size());
			strcat(strtemp,buf2);
			strcat(strtemp,"_");
			for(int ii=0;ii<ft.data.arrayset.dims[i].opt_attrs.nb_instances;ii++)
				if(strcmp(ft.data.arrayset.dims[i].opt_attrs.instances[ii].name,"label")==0)
					strcat(strtemp,ft.data.arrayset.dims[i].opt_attrs.instances[ii].value.s);
			array = vtkFloatArray::New();
			array->SetName(strtemp);
			table->AddColumn(array);
			array->Delete();
		}
        else if(ft.data.arrayset.dims[i].type_class==H5T_STRING)
        {
			strcpy(strtemp,"dim");
			std::ostringstream buf;
			buf << i;
			std::string buffer = buf.str();
			char *buf2=new char[buffer.size()+1];
			buf2[buffer.size()]=0;
			memcpy(buf2,buffer.c_str(),buffer.size());
			strcat(strtemp,buf2);
			strcat(strtemp,"_");

			for(int ii=0;ii<ft.data.arrayset.dims[i].opt_attrs.nb_instances;ii++)
				if(strcmp(ft.data.arrayset.dims[i].opt_attrs.instances[ii].name,"label")==0)
					strcat(strtemp,ft.data.arrayset.dims[i].opt_attrs.instances[ii].value.s);
			array = vtkFloatArray::New();
			array->SetName(strtemp);
			table->AddColumn(array);
			array->Delete();
        }
        else if(ft.data.arrayset.dims[i].type_class==H5T_INTEGER)
        {
			strcpy(strtemp,"dim");
			std::ostringstream buf;
			buf << i;
			std::string buffer = buf.str();
			char *buf2=new char[buffer.size()+1];
			buf2[buffer.size()]=0;
			memcpy(buf2,buffer.c_str(),buffer.size());
			strcat(strtemp,buf2);
			strcat(strtemp,"_");
			for(int ii=0;ii<ft.data.arrayset.dims[i].opt_attrs.nb_instances;ii++)
				if(strcmp(ft.data.arrayset.dims[i].opt_attrs.instances[ii].name,"label")==0)
					strcat(strtemp,ft.data.arrayset.dims[i].opt_attrs.instances[ii].value.s);
			iarray = vtkIntArray::New();
			iarray->SetName(strtemp);
			table->AddColumn(iarray);
			iarray->Delete();
        }
    }
    array = vtkFloatArray::New();
    char *buf2=new char[dataname.size()+1];


    array->SetName(dataname.c_str());
    table->AddColumn(array);
    array->Delete();
    table->SetNumberOfRows(nbdataarray);
    
    int offsettemp=1;
    for(int j=0;j<ft.data.arrayset.nb_dims;j++)
    {
      offset=0;
      for(int k=0;k<ft.data.arrayset.dims[j].nb_values;k++)
      {
        for(int itemp=0;itemp<offsettemp;itemp++)
        {

          if(ft.data.arrayset.dims[j].type_class==H5T_FLOAT)
          {
            table->SetValue(offset,j,ft.data.arrayset.dims[j].values.f[k]);
          }
          else if(ft.data.arrayset.dims[j].type_class==H5T_STRING)
          {
            table->SetValue(offset,j,float(k));
          }
          else if(ft.data.arrayset.dims[j].type_class==H5T_INTEGER)
          {
            table->SetValue(offset,j,ft.data.arrayset.dims[j].values.i[k]);
          }
          offset=offset+1;
        }
        if(offset<nbdataarray){
          if(k==(ft.data.arrayset.dims[j].nb_values-1)) k=-1;}
      }
      offsettemp=offsettemp*(ft.data.arrayset.dims[j].nb_values);
    }
    
    for(int i=0;i<nbdataarray;i++)
    {
      if(ft.data.arrayset.data.type_class==H5T_FLOAT)
        table->SetValue(i,ft.data.arrayset.nb_dims,ft.data.arrayset.data.values.f[i]);
      else if(ft.data.arrayset.data.type_class==H5T_COMPOUND)
      {
        float module;
        module=tools.abs_complex(ft.data.arrayset.data.values.c[i]);
        table->SetValue(i,ft.data.arrayset.nb_dims,module);
      }  
    }
    AH5_free_floatingtype(&ft);
    return 1;
}
