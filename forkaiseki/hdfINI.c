#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <hdf5.h>

// The number of cells in the X, Y dimensions
#define NX 512
#define NY 512
#define NZ 512
//#define M_PI 3.1415926535897932


void
write_hdf5_data()
{
    hid_t     file_id;
    file_id = H5Fcreate("xdmf3d.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // Create the coordinate data.
    
    /*
    float *x = (float *) malloc((NX) * sizeof(float)); //メモリ領域の確保 //double? float?
    float *y = (float *) malloc((NY) * sizeof(float));
    float *z = (float *) malloc((NZ) * sizeof(float));


    //mesh
    FILE *fp; // FILE型構造体
    char fname[] = "cdnt.DAT";
    int i;

    
    printf("point11 \n");
    fp =  fopen( fname , "r"); // ファイルを開く。失敗するとNULLを返す。
    if(fp == NULL) {
      printf("%s file not open!\n", fname);
    // return -1;
    } else {
      //printf("%s file opened!\n", fname);
       for(i=0; i<NX; i++){
        fscanf(fp,"%f",&(x[i]));
      }
    }
    fclose(fp);

    for(i=0; i<NX; i++){
        y[i] = x[i];
        z[i] = x[i];
       }
    //mesh

    printf("point12 \n");

    */
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
    /*
    float *rho = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *velocityx = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *velocityy = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *velocityz = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *pressure = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Bfieldx = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Bfieldy = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Bfieldz = (float *) malloc(NX*NY*NZ * sizeof(float));
    */
    /*
    float *nH = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *np = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *nH2 = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *nHe = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *nHep = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *C = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *CO = (float *) malloc(NX*NY*NZ * sizeof(float));
    float *Cp = (float *) malloc(NX*NY*NZ * sizeof(float));
    */
    float *Phi = (float *) malloc((NX)*(NY)*(NZ) * sizeof(float));
    
    char filename[50];
    int j,k,n,tot;
    /*
      FILE *file;
      file = fopen("test.dat","wb"); //binari
      fclose(file);
    */
    //sprintf(filename, "2Dall%03d.DAT", i);
    sprintf(filename, "final.DAT");
    fp = fopen(filename , "r");
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

                    fscanf(fp,"%f",&(Phi[tot])); //%f,%f,%f,%f,%f,%f,%f",\
                            //&(rho[tot]),&(velocityx[tot]),&(velocityy[tot]),&(velocityz[tot]),&(pressure[tot]),&(Bfieldx[tot]),&(Bfieldy[tot]),&(Bfieldz[tot]));
                }
            }
        }

    fclose(fp);
    
  
    
   
    //Phi
    dims[0] = (NZ);
    dims[1] = (NY);
    dims[2] = (NX);
    dataspace_id = H5Screate_simple(3, dims, NULL);
    
    dataset_id = H5Dcreate(file_id, "/Phi", H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, Phi);
    
    status = H5Dclose(dataset_id);
    
    status = H5Sclose(dataspace_id);

 
    free(Phi);
 //   free(Phiexa);
 //   free(dPhi);
 //   free(dprPhi);
    //free(rhoini);
    status = H5Fclose(file_id);
    printf("point15 \n");
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
    /*
    fprintf(xmf, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",(NX));
    fprintf(xmf, "           xdmf3d.h5:/X_1D\n");
    fprintf(xmf, "        </DataItem>\n");
    fprintf(xmf, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",(NY));
    fprintf(xmf, "           xdmf3d.h5:/Y_1D\n");
    fprintf(xmf, "        </DataItem>\n");
    fprintf(xmf, "        <DataItem Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",(NZ));
    fprintf(xmf, "           xdmf3d.h5:/Z_1D\n");
    fprintf(xmf, "     </Geometry>\n");
     */
    //
    /*
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
     */
    fprintf(xmf, "     <Attribute Name=\"Phi\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/Phi\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    
    /*
    fprintf(xmf, "     <Attribute Name=\"Phiexa\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/Phiexa\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    
    
    
    fprintf(xmf, "     <Attribute Name=\"dPhi\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/dPhi\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    
    
    fprintf(xmf, "     <Attribute Name=\"dprPhi\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/dprPhi\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    
    
    fprintf(xmf, "     <Attribute Name=\"rhoini\" AttributeType=\"Scalar\" Center=\"Node\">\n"); //node=結び目 , center , ....
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NX, NY, NZ);
    fprintf(xmf, "        xdmf3d.h5:/rhoini\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
     */
    
    /*
    fprintf(xmf, "     <Attribute Name=\"VelocityX\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", NY+1, NX+1, NZ+1);
    fprintf(xmf, "        xdmf3d.h5:/VelocityX\n");
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
     */
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
    write_xdmf_xml();
    printf("point3 \n");
    return 0;
}
