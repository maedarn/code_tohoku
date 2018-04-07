
#include<stdio.h>
#include<math.h>

void Force(double v[3],double B[3],double E[3],double F[3])
{
    F[0] = v[1] * B[2] - v[2] * B[1] + E[0];
    F[1] = v[2] * B[0] - v[0] * B[2] + E[1];
    F[2] = v[0] * B[1] - v[1] * B[0] + E[2];
}

void RK4(double x[3],double v[3],double B[3],double E[3],double dt)
{
    int i;
    double kx1[3],kx2[3],kx3[3],kx4[3];
    double kv1[3],kv2[3],kv3[3],kv4[3];
    double u[3];
    double F[3];
    //1st step
    Force(v,B,E,F);
    for(i=0;i<3;i++)
    {
        kv1[i]=F[i]*dt;
        kx1[i]=v[i]*dt;
    }
    //2nd step
    for(i=0;i<3;i++)
    {
        u[i] = v[i] + 0.5*kv1[i];
    }
    Force(u,B,E,F);
    for(i=0;i<3;i++)
    {
        kv2[i]=F[i]*dt;
        kx2[i]=u[i]*dt;
    }
    //3rd step
    for(i=0;i<3;i++)
    {
        u[i] = v[i] + 0.5*kv2[i];
    }
    Force(u,B,E,F);
    for(i=0;i<3;i++)
    {
        kv3[i]=F[i]*dt;
        kx3[i]=u[i]*dt;
    }
    //4th step
    for(i=0;i<3;i++)
    {
        u[i] = v[i] + kv3[i];
    }
    Force(u,B,E,F);
    for(i=0;i<3;i++)
    {
        kv4[i]=F[i]*dt;
        kx4[i]=u[i]*dt;
    }
    //result
    for(i=0;i<3;i++)
    {
        v[i] = v[i] + (kv1[i] + 2.0 * kv2[i] + kv3[i] + kv4[i])/6.0;
        x[i] = x[i] + (kx1[i] + 2.0 * kx2[i] + kx3[i] + kx4[i])/6.0;
    }
}

int main(void)
{
    double x[3] = {0,0,0};
    double v[3] = {1,0,0};
    double E[3] = {1.00,0.00,0.00};
    double B[3] = {0.00,0.00,1.00};
    FILE *fp;
    fp = fopen("../plasma-cal/ChargedParticle.dat","w");
    double t,dt=0.01;
    for(t=dt;t<=20*M_PI;t=t+dt)
    {
        RK4(x,v,B,E,dt);
        printf("%f %f %f %f\n",t,x[0],x[1],x[2]);
        fprintf(fp,"%f %f %f %f\n",t,x[0],x[1],x[2]);
    }
    fclose(fp);
    return 0;
}
