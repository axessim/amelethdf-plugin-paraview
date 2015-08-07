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
    char* entryPoint = NULL;

    AH5_read_str_attr(file_id, ".", AH5_A_ENTRY_POINT, &entryPoint);
    strcpy(path2, entryPoint);
  	strcat(path2, AH5_G_DS);
  	children = AH5_read_children_name(file_id, path2);
  	nb_dims = children.nb_children;
  	dims = (AH5_vector_t *) malloc((size_t) children.nb_children * sizeof(AH5_vector_t));
  	for (i = 0; i < children.nb_children; i++)
  	{
  	    if (!invalid)
  	    {
  	        strcpy(path2, entryPoint);
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
/*    vtkFloatArray *array;
    vtkIntArray *iarray;
    arrayset_t ars;
    children_t child;
    char  path[ABSOLUTE_PATH_NAME_LENGTH];
    char  attr_value[ELEMENT_NAME_LENGTH];
    char  attr[ELEMENT_NAME_LENGTH];
    char  strtemp[ELEMENT_NAME_LENGTH];
    char  strtemp2[ELEMENT_NAME_LENGTH];
    int xdim = -1; int ydim = -1; int zdim = -1;
    int nbdataarray=1;
    int max=1;
    
    child = read_children_name(file_id,"/floatingType");
    if(child.nbchild>1)
    {
    	cout<<"This is more than one arrayset in the ameletHDF file."<<endl;
        return 0;           
    }
    strcpy(attr,"floatingType");
    strcpy(path,"/floatingType");
    strcat(path,"/");
    strcat(path,child.childnames[0]);
    strcpy(attr_value,read_string_attribute(file_id,path,attr));
    if(strcmp(attr_value,"arraySet")!=0)
    {
      	cout<<("This is not an arrayset in the ameletHDF file.")<<endl;
      	return 0;
    }

    ars = read_arrayset(file_id,path);

    for(int i=0;i<ars.nbdims;i++)
            nbdataarray = nbdataarray*ars.dims[i].nbvalues;
    vtkstd::string dataname;
    dataname= child.childnames[0];
    int temp=nbdataarray;
    
    int offset = 0;

    for(int i=0;i<ars.nbdims;i++)
    {
      
      if(ars.dims[i].rvalue!=NULL)
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
        strcat(strtemp,ars.dims[i].single.label);
        array = vtkFloatArray::New();
        array->SetName(strtemp);
        table->AddColumn(array);
        array->Delete();
      }
      else if(ars.dims[i].svalue!=NULL)
      {
        array = vtkFloatArray::New();
        std::cout<<ars.dims[i].single.label<<std::endl;
        array->SetName(ars.dims[i].single.label);
        table->AddColumn(array);
        array->Delete();
      }
      else if(ars.dims[i].ivalue!=NULL)
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
        strcat(strtemp,ars.dims[i].single.label);
        iarray = vtkIntArray::New();
        iarray->SetName(strtemp);
        table->AddColumn(iarray);
        iarray->Delete();
      }
    }
    array = vtkFloatArray::New();
    array->SetName(child.childnames[0]);
    table->AddColumn(array);
    array->Delete();
    table->SetNumberOfRows(nbdataarray);
    
    int offsettemp=1;
    for(int j=0;j<ars.nbdims;j++)
    {
      offset=0;
      for(int k=0;k<ars.dims[j].nbvalues;k++)
      {
        for(int itemp=0;itemp<offsettemp;itemp++)
        {
          //std::cout<<"offset="<<offset<<std::endl;
          if(ars.dims[j].rvalue!=NULL)
          {
            table->SetValue(offset,j,ars.dims[j].rvalue[k]);
          }
          else if(ars.dims[j].svalue!=NULL)
          {
            table->SetValue(offset,j,float(k));
          }
          else if(ars.dims[j].ivalue!=NULL)
          {
            table->SetValue(offset,j,ars.dims[j].ivalue[k]);
          }
          offset=offset+1;
        }
        if(offset<nbdataarray){
          if(k==(ars.dims[j].nbvalues-1)) k=-1;}
      }
      offsettemp=offsettemp*(ars.dims[j].nbvalues);
    }
    
    for(int i=0;i<nbdataarray;i++)
    {
      if(ars.data.rvalue!=NULL)
        table->SetValue(i,ars.nbdims,ars.data.rvalue[i]);
      else if(ars.data.cvalue!=NULL)
      {
        float module;
        module=abs_complex(ars.data.cvalue[i]);
        table->SetValue(i,ars.nbdims,module);
      }  
    }
   
    if(child.nbchild>1) free(child.childnames[0]);
        free(child.childnames);*/
    return 1;
}
