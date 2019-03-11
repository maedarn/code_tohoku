#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include <hdf5.h>

// The number of cells in the X, Y dimensions
#define NX 512
#define NY 512
#define NZ 512
#define last 1



void
write_hdf5_data()
{
  int time;
  int ren;
  //char filename[160];
  //char filename1[160];
  char filename[75];
  char filename1[74];
  char pwd[40];
  //char pwd[100];
  char psps[30];
  char DAT[5];
  char HFILE[4];


  sprintf(psps,"/Users/maeda/Desktop/kaiseki/" );
  sprintf(DAT,".dat" );
  sprintf(HFILE,".h5" );
  //printf("文字列を入力してください\n");
  sprintf(pwd,"multigrid-plane-512(512core)-gosa/final" );
  //Trim(pwd);
  //multigrid-plane-512(512core)
  //printf("文字列を入力してください\n");
  //scanf("%s",pwd);
  //Trim(pwd);
  for(time=0;time<last;time++){
    printf("%d",time);
    hid_t     file_id;
    //char s1[6] = {'\0'};
    //sprintf(s1, "%6d", time);
    //strcat(str1, str2);
    //snprintf(filename, 160, "%s%s%s", psps, pwd , DAT);
    snprintf(filename, 75, "%s%s%s", psps, pwd , DAT);
    //Trim(filename);
    //snprintf(filename1, 160, "%s%s%s", psps, pwd , HFILE);
    snprintf(filename1, 74, "%s%s%s", psps, pwd , HFILE);
    //Trim(filename1);
    //sprintf(filename, strcat(psps,pwd,".DAT"));
    //sprintf(filename1, "/Users/maeda/Desktop/kaiseki/testcode25/final.h5");
    file_id = H5Fcreate(filename1, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create the coordinate data.


    //mesh
    FILE *fp; // FILE型構造体
    int i;
    //float *c = (float *) malloc(1*sizeof(float));



    float *Phir = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *Phil = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *Phim = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *a    = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *b    = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *c    = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *d    = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    float *e    = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    //char filename[50];
    int j,k,n,tot;

    //sprintf(filename, "final.DAT");
    fp = fopen(filename  , "r");
    printf("%s\n", filename);
    if(fp == NULL) {
      printf("ファイルを開くことが出来ませんでした．\n");
      return;
    }
        for(k=0;k<NZ;k++){
            for(j=0;j<NY;j++){
       //n = 0;
                for(n=0;n<NX;n++){
                    tot= n+(NX)*j+(NX)*(NZ)*k;
       //  ファイルが終わりでない「かつ」配列を飛び出さないうちは，読み込みを続ける
              //while ( ! feof(fp) && n < NX) {
              // fscanf(fp,"%f,%f,%f,%f,%f,%f,%f,%f", &(rho[i][j][k]),&(velocityx[i][j][k]),\
              //&(velocityy[n][j][k]),&(velocityz[n][j][k]),&(pressure[n][j][k]),&(Bfieldx[n][j][k]),&(Bfieldy[n][j][k]),&(Bfieldz[n][j][k])) //,&(rho[i][j][k]));  // 読み込む個数
              //        n++; }


		    //fscanf(fp,"%f %f %f",&(Phir[tot]),&(Phil[tot]),&(Phim[tot]));
                    fscanf(fp,"%f %f %f %f %f %f %f %f",&(a[tot]),&(b[tot]),&(Phir[tot]),&(d[tot]),&(Phil[tot]),&(e[tot]),&(Phim[tot]),&(c[tot])); //%f,%f,%f,%f,%f,%f,%f",\  点で区切らない
                            //&(rho[tot]),&(velocityx[tot]),&(velocityy[tot]),&(velocityz[tot]),&(pressure[tot]),&(Bfieldx[tot]),&(Bfieldy[tot]),&(Bfieldz[tot]));
                    //printf("%f,%f,%f,%f,%f\n",a[i], b[i], Phir[i], Phil[i], Phim[i]);
                    //printf("%d,%d,%d\n",n,j,k);
                }
            }
        }

    fclose(fp);


    // Write the data file.
    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[3];
    herr_t    status;
    const char *coordNames[] = {"/X_1D", "/Y_1D","/Z_1D"}; //name

    /* Write separate coordinate arrays for the x and y coordinates. */
    //for(int did = 0; did < 3; ++did)
    //{
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      //dataspace_id = H5Screate_simple(1, dims, NULL);
      // Write the scalar data.

    //Phi
        dims[0] = (NZ);
        dims[1] = (NY);
        dims[2] = (NX);                /**/
        dataspace_id = H5Screate_simple(3, dims, NULL);

        dataset_id = H5Dcreate(file_id, "/Phimns", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phir);

        status = H5Dclose(dataset_id);

        status = H5Sclose(dataspace_id);
      
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/fxexa", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phil);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);
      
      
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/Phical", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, a);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);
      
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/Phiexa", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, b);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);
      
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/fxcal", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, d);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);
      
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/fxmns", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, e);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);
      
      
      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/fxratio", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phim);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);

      dims[0] = (NZ);
      dims[1] = (NY);
      dims[2] = (NX);
      dataspace_id = H5Screate_simple(3, dims, NULL);
      
      dataset_id = H5Dcreate(file_id, "/fparsent", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, c);
      
      status = H5Dclose(dataset_id);
      
      status = H5Sclose(dataspace_id);
      
      free(Phir);
      free(Phil);
      free(Phim);
      free(a);
      free(b);
      free(c);
      free(d);
      free(e);
    status = H5Fclose(file_id);
    printf("point15 \n");
  }
}
void
write_xdmf_xml()
{
    FILE *xmf = 0;

    /*
     * Open the file and write the XML description of the mesh..
     */
    xmf = fopen("xdmf3d.xmf", "w");
    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n");
    //fprintf(xmf, "     <Topology TopologyType=\"2DSMesh\" NumberOfElements=\"%d %d\"/>\n", NY+1, NX+1); //Curvilinear
    fprintf(xmf, "     <Topology TopologyType=\"3DCoRectMesh\" NumberOfElements=\"%d %d %d\"/>\n", NY, NX, NZ); //Axis are perpendicular and spacing is constant
    //fprintf(xmf, "     <Geometry GeometryType=\"X_Y\">\n");
    //
    //typeで決める
    fprintf(xmf, "     <Geometry GeometryType=\"VXVYVZ\">\n");

    fprintf(xmf, "     <Attribute Name=\"Phir\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/Phir\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "     <Attribute Name=\"Phil\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/Phil\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");

    fprintf(xmf, "     <Attribute Name=\"Phim\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/Phim\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
}

int
main(int argc, char *argv[])
{
    printf("point1 \n");
    write_hdf5_data();
    printf("point2 \n");
 //   write_xdmf_xml();
    printf("point3 \n");
    return 0;
}
