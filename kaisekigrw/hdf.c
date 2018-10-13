#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <hdf5.h>

// The number of cells in the X, Y dimensions
#define NX 128
#define NY 128
#define NZ 128
#define M_PI 3.1415926535897932


void
write_hdf5_data()
{
    hid_t     file_id;
    file_id = H5Fcreate("xdmf3d.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create the coordinate data.
    float *x = (float *) malloc((NX) * sizeof(float)); //メモリ領域の確保 //double? float?
    float *y = (float *) malloc((NY) * sizeof(float));
    float *z = (float *) malloc((NZ) * sizeof(float));


    //mesh
    FILE *fp; // FILE型構造体
    char fname[] = "cdnt.DAT";
    int i;

    fp = fopen( fname , "r"); // ファイルを開く。失敗するとNULLを返す。
    if(fp == NULL) {
      printf("%s file not open!\n", fname);
      return -1;
    } else {
      //printf("%s file opened!\n", fname);
       for(i=0; i<NX; i++){
        fscanf(file,"%f",&(x[i]));
      }
    }
    fclose(fp);

    for(i=0; i<NX; i++){
        y[i] = x[i]
        z[i] = x[i]
       }
    //mesh



    /* //mesh structure
       float *x = (float *) malloc((NX+1)*(NY+1) * sizeof(float));
       float *y = (float *) malloc((NX+1)*(NY+1) * sizeof(float));
    int ndx = 0;
    for (int j = 0; j < NY+1; j++)
    {
        float yt = j*1.0 / NY;
        float angle = yt * M_PI;
        for (int i = 0; i < NX+1; i++)
        {
            float xt = i*1.0 / NX;
            float R = (1.-xt)*2. + xt*5.;

            x[ndx] = R * cos(angle);
            y[ndx] = R * sin(angle);
            ndx++;
        }
    }
    */


    //while ver.
    /*
    fp = fopen("sincos.csv", "r");
    if(fp == NULL) {
      printf("ファイルを開くことが出来ませんでした．¥n");
      return;
    }
    n = 0;
    //  ファイルが終わりでない「かつ」配列を飛び出さないうちは，読み込みを続ける
    while ( ! feof(fp) && n < 512) {
      fscanf(fp, "%f", &(siny[n]));
      n++;
    }
    fclose(fp);
    */


    // Create the scalar data.
    float *rho = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *velocityx = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *velocityy = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *velocityx = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *pressure = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Bfieldx = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Bfieldy = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Bfieldz = (float *) malloc(NX*NY*NZ * sizeof(float));
    /*
    float *nH = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *np = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *nH2 = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *nHe = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *nHep = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *C = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *CO = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Cp = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Phi = (float *) malloc(NX*NY*NZ * sizeof(float));
    */
    char filename[50];
    int j,k,n;
    /*
      FILE *file;
      file = fopen("test.dat","wb"); //binari
      fclose(file);
    */
    //sprintf(filename, "2Dall%03d.DAT", i);
    sprintf(filename, "Dall800.DAT");
    fp = fopen(filename , "r");
    if(fp == NULL) {
      printf("ファイルを開くことが出来ませんでした．\n");
      return;
    }
    for(k=0;k<NZ;k++){
      for(j=0;j<NZ;j++){
       n = 0;
       //  ファイルが終わりでない「かつ」配列を飛び出さないうちは，読み込みを続ける
       while ( ! feof(fp) && n < NX) {
         fscanf(fp,"%f,%f,%f,%f,%f,%f,%f,%f", &(rho[i][j][k]),&(velocityx[i][j][k]),\
&(velocityy[n][j][k]),&(velocityz[n][j][k]),&(pressure[n][j][k]),&(Bfieldx[n][j][k]),&(Bfieldy[n][j][k]),&(Bfieldz[n][j][k])) //,&(rho[i][j][k]));  // 読み込む個数
         n++;
       }
      }
    }
    fclose(fp);


    /*
    for (int j = 0; j < NY; j++)
    {
        for (int i = 0; i < NX; i++)
        {
            int ndx = j * NX + i;
            pressure[ndx] = (float) j;
        }
    }

    //float *velocityx = (float *) malloc((NX+1)*(NY+1) * sizeof(float));

    for (int j = 0; j < NY+1; j++)
    {
        for (int i = 0; i < NX+1; i++)
        {
            int ndx = j * (NX+1) + i;
            velocityx[ndx] = (float) i;
        }
    }
    */


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
      dataspace_id = H5Screate_simple(1, NX, NULL);

      dataset_id = H5Dcreate(file_id, coordNames[0], H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,\
H5P_DEFAULT, x);

      status = H5Dclose(dataset_id);

      status = H5Sclose(dataspace_id);
      //
      dataspace_id = H5Screate_simple(1, NY, NULL);

      dataset_id = H5Dcreate(file_id, coordNames[1], H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,\
H5P_DEFAULT, y);

      status = H5Dclose(dataset_id);

      status = H5Sclose(dataspace_id);
      //
      dataspace_id = H5Screate_simple(1, NZ, NULL);

      dataset_id = H5Dcreate(file_id, coordNames[2], H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,\
H5P_DEFAULT, z);

      status = H5Dclose(dataset_id);

      status = H5Sclose(dataspace_id);
      //}

    // Write the scalar data.
    //rho
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/rho", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rho);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    //vx
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/velocityx", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocityx);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    //vy
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/velocityy", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocityy);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    //vz
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/velocityz", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocityz);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);

    //p
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/Pressure", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, pressure);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    //bx
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/Bfieldx", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Bfieldx);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    //by
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/Bfieldy", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Bfieldy);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    //bz
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/Bfieldz", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Bfieldz);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);


    /*
    dims[0] = NY + 1;
    dims[1] = NX + 1;
    dataspace_id = H5Screate_simple(2, dims, NULL);

    dataset_id = H5Dcreate(file_id, "/VelocityX", H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, velocityx);

    status = H5Dclose(dataset_id);

    status = H5Sclose(dataspace_id);
    */

    // Free the data.
    free(x);
    free(y);
    free(z);
    free(rho);
    free(velocityx);
    free(velocityy);
    free(velocityz);
    free(pressure);
    free(Bfieldx);
    free(Bfieldy);
    free(Bfieldz);

    status = H5Fclose(file_id);
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
    /*
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1) , (NZ+1));
    fprintf(xmf, "        xdmf3d.h5:/X\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", (NY+1), (NX+1));
    fprintf(xmf, "        xdmf2d.h5:/Y\n");
    fprintf(xmf, "       </DataItem>\n");
    */
    fprintf(xmf, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",(NX));
    fprintf(xmf, "           xdmf3d.h5:/X_1D\n");
    fprintf(xmf, "        </DataItem>\n");
    fprintf(xmf, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",(NY));
    fprintf(xmf, "           xdmf3d.h5:/Y_1D\n");
    fprintf(xmf, "        </DataItem>\n");
    fprintf(xmf, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",(NZ));
    fprintf(xmf, "           xdmf3d.h5:/Z_1D\n");
    fprintf(xmf, "     </Geometry>\n");
    //
    //rho
    fprintf(xmf, "     <Attribute Name=\"rho\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX, NZ);
    fprintf(xmf, "        xdmf3d.h5:/rho\n");
    fprintf(xmf, "       </DataItem>\n");
    //v
    fprintf(xmf, "     <Attribute Name=\"vectors\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"Function\" Function=\"join($0, $1, $2)\" Dimensions=\"%d %d %d 3\">\n", NY, NX, NZ); //3=dimension
    fprintf(xmf, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX*NY*NZ);
    fprintf(xmf, "           xdmf3d.h5:/velocityx\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX*NY*NZ);
    fprintf(xmf, "           xdmf3d.h5:/velocityy\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX*NY*NZ);
    fprintf(xmf, "           xdmf3d.h5:/velocityz\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    //p
    fprintf(xmf, "     <Attribute Name=\"Pressure\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY, NX, NZ);
    fprintf(xmf, "        xdmf3d.h5:/Pressure\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    //B
    fprintf(xmf, "     <Attribute Name=\"vectors\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem ItemType=\"Function\" Function=\"join($0, $1, $2)\" Dimensions=\"%d %d %d 3\">\n", NY, NX, NZ); //3=dimension
    fprintf(xmf, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX*NY*NZ);
    fprintf(xmf, "           xdmf3d.h5:/Bfieldx\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX*NY*NZ);
    fprintf(xmf, "           xdmf3d.h5:/Bfieldy\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "         <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",NX*NY*NZ);
    fprintf(xmf, "           xdmf3d.h5:/Bfieldz\n");
    fprintf(xmf, "         </DataItem>\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    /*
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1, NZ+1);
    fprintf(xmf, "        xdmf3d.h5:/VelocityX\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);
    */
}

int
main(int argc, char *argv[])
{
    write_hdf5_data();
    write_xdmf_xml();

    return 0;
}
