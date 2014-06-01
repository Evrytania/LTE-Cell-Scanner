/*Writed by jxj.All rights reserved*/
/*Revised by jxj. Jan. 2009. for INF; NaN; etc... */
#include "mex.h"
#define MAXSTATES 64
#define MAXINPUTS 2
#define MAXDATALEN 4096
#define APZERO (1e-20)

static double a[MAXSTATES*(MAXDATALEN+1)] = {0.0};
static double b[MAXSTATES*MAXDATALEN] = {0.0};
static double y[MAXSTATES*MAXDATALEN*MAXINPUTS] = {0.0};

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
/*-------------------------output---------------------------input-----------*/
{
int N,M,S,O,L1,L2,sizeofa,sizeofb,k,m1,m,tk,d,o;
double *pci,*pdi,reciofnumStates,*pco,*pdo,*nextStates,*outputs,*priorStates,*m1ivo,sumtmp,invM, invO;

pci=mxGetPr(prhs[1]);
pdi=mxGetPr(prhs[2]);
N=mxGetN(prhs[2]);

M=mxGetScalar(mxGetField(prhs[0],0,"numInputSymbols"));
S=mxGetScalar(mxGetField(prhs[0],0,"numStates"));
O=mxGetScalar(mxGetField(prhs[0],0,"numOutputSymbols"));
L1=mxGetScalar(mxGetField(prhs[0],0,"L1"));
L2=mxGetScalar(mxGetField(prhs[0],0,"L2"));

sizeofa=S*(N+1);
memset(a, 0, sizeofa*sizeof(double));
//a=mxCalloc(S*(N+1),sizeof(double));
//little=0;
//(*a)=1-(S-1)*little;
//for(k=1;k<S;k++) (*(a+k))=little;
(*a)=1;
        
sizeofb=sizeofa-S;
memset(b, 0, sizeofb*sizeof(double));
//b=mxCalloc(sizeofb,sizeof(double));
reciofnumStates=(1.00/(double)S);
invM=(1.00/(double)M);
invO=(1.00/(double)O);
for(k=sizeofb-1;k>sizeofb-S-1;k--) 
{
    (*(b+k))=reciofnumStates;
    //(*(b+k))=1;
}

plhs[0]=mxCreateDoubleMatrix(O,N,mxREAL);
plhs[1]=mxCreateDoubleMatrix(M,N,mxREAL);
pco=mxGetPr(plhs[0]);
pdo=mxGetPr(plhs[1]);

nextStates=mxGetPr(mxGetField(prhs[0],0,"nextStates"));
outputs=mxGetPr(mxGetField(prhs[0],0,"outputs"));
priorStates=mxGetPr(mxGetField(prhs[0],0,"priorStates"));
m1ivo=mxGetPr(mxGetField(prhs[0],0,"m1ivo"));

//y=mxCalloc(sizeofb*M,sizeof(double));

for(k=0;k<N;k++)
{
    for(m1=0;m1<S;m1++)
    {
        for(d=0;d<M;d++)
        {
            (*(y+m1+k*S+d*sizeofb))=(*(pci+(int)(*(outputs+m1+d*S))+k*O)) * (*(pdi+d+k*M));
        }
    }
   sumtmp=0;
   for(m=0;m<S;m++)
    {
        for(tk=0;tk<L1;tk++)
        {
            (*(a+m+(k+1)*S))=(*(a+m+(k+1)*S))+(*(a+(int)(*(priorStates+m+2*tk*S))+k*S))*(*(y+(int)(*(priorStates+m+2*tk*S))+k*S+(int)(*(priorStates+m+(2*tk+1)*S))*sizeofb));
        }
        sumtmp=sumtmp+(*(a+m+(k+1)*S));
    }
   if(sumtmp == 0.0)
    //if(sumtmp < APZERO)
    {
        for(m1=0;m1<S;m1++)
        {
           (*(a+m1+(k+1)*S))=reciofnumStates;
           //(*(a+m1+(k+1)*S))=APZERO;
        }
    }
    else
    {
        for(m1=0;m1<S;m1++)
        {
           (*(a+m1+(k+1)*S))=(*(a+m1+(k+1)*S))/sumtmp;
        }
    }
}

for(k=N-2;k>-1;k--)
{
    sumtmp=0;
    for(m1=0;m1<S;m1++)
    {
        for(tk=0;tk<M;tk++){(*(b+m1+k*S))=(*(b+m1+k*S))+(*(b+(int)(*(nextStates+m1+tk*S))+(k+1)*S))*(*(y+m1+(k+1)*S+tk*sizeofb));}
        sumtmp=sumtmp+(*(b+m1+k*S));
    }
    if(sumtmp == 0.0)
    //if(sumtmp < APZERO)
    {
        for(m=0;m<S;m++)
        {
            (*(b+m+k*S))=reciofnumStates;
            //(*(b+m+k*S))=APZERO;
        }
    }
    else
    {
        for(m=0;m<S;m++)
        {
            (*(b+m+k*S))=(*(b+m+k*S))/sumtmp;
        }
    }
     
}

for(k=0;k<N;k++)
{
    sumtmp=0;
    for(d=0;d<M;d++)
    {
        for(tk=0;tk<S;tk++){(*(pdo+d+k*M))=(*(pdo+d+k*M))+(*(a+tk+k*S))*(*(y+tk+k*S+d*sizeofb))*(*(b+(int)(*(nextStates+tk+d*S))+k*S));}
        sumtmp=sumtmp+(*(pdo+d+k*M));
    }
    if(sumtmp == 0.0)
    //if(sumtmp < APZERO)
    {
        for(d=0;d<M;d++)
        {
            (*(pdo+d+k*M))=invM;
        }
    }
    else
    {
        for(d=0;d<M;d++)
        {
            (*(pdo+d+k*M))=(*(pdo+d+k*M))/sumtmp;
        }
    }
    sumtmp=0;
    for(o=0;o<O;o++)
    {
        for(tk=0;tk<L2;tk++){(*(pco+o+k*O))=(*(pco+o+k*O))+(*(a+(int)(*(m1ivo+o+2*tk*O))+k*S))*(*(y+(int)(*(m1ivo+o+2*tk*O))+k*S+(int)(*(m1ivo+o+(2*tk+1)*O))*sizeofb))*(*(b+(int)(*(nextStates+(int)(*(m1ivo+o+2*tk*O))+(int)(*(m1ivo+o+(2*tk+1)*O))*S))+k*S));}
        sumtmp=sumtmp+(*(pco+o+k*O));
    }
    if(sumtmp == 0.0)
    //if(sumtmp < APZERO)
    {
        for(o=0;o<O;o++)
        {
            (*(pco+o+k*O))=invO;
        }
    }
    else
    {
        for(o=0;o<O;o++)
        {
            (*(pco+o+k*O))=(*(pco+o+k*O))/sumtmp;
        }
    }
}
}
