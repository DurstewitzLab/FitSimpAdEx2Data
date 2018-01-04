/* computes the voltage trace of a single cell modeled by the simpAdEx with abitrary input */

/********************/
/* Include libaries */
/********************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stddef.h>
#include "mex.h"

/******************************************/
/* define global parameters and variables */
/******************************************/

struct NeuronPara
{
    double Cm;
    double gL;
    double EL;
    double sf;
    double Vup;
    double tcw;
    double a;
    double b;
    double Vr;
    double Vth;
    double Iinj;
    
    double v[2];
    double dv[2];
};

double wV;
double D0;

/***************************/
/* computational functions */
/***************************/

/* differential equations defining the simpAdEx model */
void DiffEq (struct NeuronPara np, double *v, double dt, double *dv)
{
    double dD0, I_ex;
    
    /* Exponential term */   
    I_ex=np.gL*np.sf*exp((v[0]-np.Vth)/np.sf);
    
    /* V-nullcline at v[0] */
    wV = np.Iinj - np.gL*(v[0]-np.EL) + I_ex;
    
    /* Calculation of D_0 */
    D0=(np.Cm/np.gL) * wV;
    
    /* Compute membrane potential derivative from all currents */
    dv[0] = (np.Iinj - np.gL*(v[0]-np.EL) - v[1] + I_ex)/np.Cm;    
    
    /* derivative of D_0 with respect to V */
    dD0=np.Cm*(exp((v[0]-np.Vth)/np.sf)-1);
    
    /* second differential equation */
    if((v[1]>wV-D0/np.tcw) && (v[1]<wV+D0/np.tcw) && v[0]<=np.Vth)
    {
        dv[1]=-(np.gL*(1-exp((v[0]-np.Vth)/np.sf)) + dD0/np.tcw)*dv[0];}
    else
    {
        dv[1]=0;}   
}

/* numerical solution of the differential equations defined in DiffEq */
void update (struct NeuronPara *np, double dt)
{
    double dv1[2],dv2[2],v[2];
    int i;
    
    for (i=0;i<2;i++) 
        v[i]=(*np).v[i];
    
    DiffEq(*np,v,0.0,dv1);
    
    for (i=0;i<2;i++)
        v[i]+=dt*dv1[i];
    
    DiffEq(*np,v,dt,dv2); 
    for (i=0;i<2;i++)
    {
        (*np).v[i]+=dt/2.0*(dv1[i]+dv2[i]);
        (*np).dv[i]=dt/2.0*(dv1[i]+dv2[i]);
    }
    
    /* Make a jump in w when it approaches the nullcline wV from the right to avoid singularities */
    if(((*np).v[1]>wV-D0/(*np).tcw) && ((*np).v[1]<wV+D0/(*np).tcw) && (*np).v[0]<=(*np).Vth) 
    {
        (*np).v[1]=wV-(D0/(*np).tcw);
    }
}

/* detect externally defined events */
double decide_event(double t, double *event ,int N, int *x)
{
    int l;
    
    if (t<event[1] || t>=event[3*(N-1)+2])
    {
        return(0);
    }
    
    for(l=*x;l<3*N;l+=3)
    {   
        if(t>=event[l+1] && t<event[l+2])
        {
            *x = l;
            return(event[l]);
            break;
        } 
    }
}    


/****************/
/* Mex function */
/****************/

void mexFunction(int nlhs, mxArray *plhs[], 
                 int nrhs, const mxArray *prhs[])
{
    double *NeuPar, *Events, *InitialPar, *SimPar;
    double T0,Tstop,dt0,ti,vp,wp,tsp,w_end;
    int N;
    int l0=0;
    FILE *fp;
    struct NeuronPara *NPtr;
    
    NeuPar=mxGetPr(prhs[0]);
    NPtr=(struct NeuronPara *) mxMalloc(sizeof(struct NeuronPara));
    (*NPtr).Cm=NeuPar[0];
    (*NPtr).gL=NeuPar[1];
    (*NPtr).EL=NeuPar[2];
    (*NPtr).sf=NeuPar[3];
    (*NPtr).Vup=NeuPar[4];
    (*NPtr).tcw=NeuPar[5];
    (*NPtr).a=NeuPar[6];
    (*NPtr).b=NeuPar[7];
    (*NPtr).Vr=NeuPar[8];
    (*NPtr).Vth=NeuPar[9];
    
    Events=mxGetPr(prhs[1]);
    (*NPtr).Iinj=Events[0];
    N=mxGetN(prhs[1]);  
    
    InitialPar=mxGetPr(prhs[2]);
    (*NPtr).v[0]=InitialPar[0];
    (*NPtr).v[1]=InitialPar[1];
    
    SimPar=mxGetPr(prhs[3]);
    T0=SimPar[0];
    Tstop=SimPar[1];
    dt0=SimPar[2];  
    
    
    
    /* actual simulation */
    fp = fopen("voltage_data.dat", "w");
    
    fprintf(fp,"%f %f %f",T0,(*NPtr).v[0],(*NPtr).v[1]);
    fprintf(fp,"\n");

     for (ti=T0;ti<=Tstop;) 
     {      
        while((*NPtr).v[0]<=(*NPtr).Vup && ti<=Tstop)
        {           
            (*NPtr).Iinj=decide_event(ti,Events,N,&l0);
            vp=(*NPtr).v[0];
            wp=(*NPtr).v[1];
            update(NPtr,dt0);            
            
            if((*NPtr).v[0]<=(*NPtr).Vup)
            {
                fprintf(fp,"%f %f %f",ti+dt0,(*NPtr).v[0],(*NPtr).v[1]);
                fprintf(fp,"\n");
            }
            ti=ti+dt0; 
        }
        
        if (ti<=Tstop)
        {
            /* stop upswing when Vup is reached and compute spike time*/
            tsp=((*NPtr).Vup-vp)*((dt0)/((*NPtr).v[0]-vp)) + (ti-dt0);
            w_end=(((*NPtr).v[1]-wp)/(dt0))*(tsp-(ti-dt0)) + wp;            
            fprintf(fp,"%f %f %f",tsp,(*NPtr).Vup,w_end);
            fprintf(fp,"\n");
            
            /* artificial reset */
            ti=tsp;
            (*NPtr).v[0]=(*NPtr).Vr;
            (*NPtr).v[1]=w_end+(*NPtr).b;
            fprintf(fp,"%f %f %f",ti,(*NPtr).v[0],(*NPtr).v[1]);
            fprintf(fp,"\n");
        }
     }
    
    fclose(fp);
    
    return;
}


/* (c) 2012 L. Hertaeg, J. Hass and D. Durstewitz,
  Central Institute of Mental Health, Mannheim University of Heidelberg 
  and BCCN Heidelberg-Mannheim */
