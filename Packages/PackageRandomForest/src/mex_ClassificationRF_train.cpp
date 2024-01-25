#include <math.h>
#include "memory.h"
#include "mex.h"

#ifdef OCTAVE
    #include "stdio.h"
#endif

#define DEBUG_ON 0
void classRF(double *x, int *dimx, int *cl, int *ncl, int *cat, int *maxcat,
	     int *sampsize, int *strata, int *Options, int *ntree, int *nvar,
	     int *ipi, double *classwt, double *cut, int *nodesize,
	     int *outcl, int *counttr, double *prox,
	     double *imprt, double *impsd, double *impmat, int *nrnodes,
	     int *ndbigtree, int *nodestatus, int *bestvar, int *treemap,
	     int *nodeclass, double *xbestsplit, double *errtr,
	     int *testdat, double *xts, int *clts, int *nts, double *countts,
	     int *outclts, int labelts, double *proxts, double *errts,
         int *inbag, int print_verbose_tree_progression);

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    if(nrhs==23);
    else{
		mexPrintf("Too less/many parameters: You supplied %d",nrhs);
		return;
	}
    
    double *_tmp_d;
        
    int i;
    int p_size = (int)mxGetScalar(prhs[16]);
    int n_size = (int)mxGetScalar(prhs[15]);
    double *x = mxGetPr(prhs[0]);
    int *y = (int*)mxGetData(prhs[1]);
    
    if (DEBUG_ON){
        //print few of the values
        //for(i=0;i<10;i++)
        //    mexPrintf("%d,",y[i]);
    }

    int nclass = (int)((double)mxGetScalar(prhs[2]));
	double nclass_d = ((double)mxGetScalar(prhs[2]));
    int* cat = (int*)mxGetData(prhs[5]);//calloc(p_size,sizeof(int));
                                        //for(i=0;i<p_size;i++) cat[i]=1;
    
    
	if (DEBUG_ON){ mexPrintf("n_size %d, p_size %d, nclass %d (%f)",n_size,p_size,nclass,(double)mxGetScalar(prhs[2]));}
	fflush(stdout);
    int maxcat=*((int*)mxGetData(prhs[6]));
    int* sampsize=(int*)mxGetData(prhs[7]);
    int nsum = *((int*)mxGetData(prhs[14]));
    int* strata = (int*)mxGetData(prhs[8]);
    //int Options[]={addclass,importance,localImp,proximity,oob_prox,do_trace,keep_forest,replace,stratify,keep_inbag};
    int* Options = (int*)mxGetData(prhs[9]);
    
    // now get individual values from the options so they can be decomposed and appropriate
    // array sizes can be set.
    int addclass = Options[0];
    int importance=Options[1];
    int localImp  =Options[2];
    int proximity =Options[3];
    int oob_prox  =Options[4];
    int do_trace  =Options[5];
    if (DEBUG_ON){ do_trace=1;} else {do_trace=0;}    
    int keep_forest=Options[6];
    int replace   =Options[7];
    int stratify  =Options[8];
    int keep_inbag=Options[9];
    
    int nsample=(int)mxGetScalar(prhs[17]);
    int dimx[]={p_size, n_size};
    
    int tst_available=(int)mxGetScalar(prhs[18]);
    
    double *xts;
    int *yts;
    int nts;
    int *outclts;
    int labelts=0;
    double* proxts;
    double* errts;

    
    
    int ntree;
    int mtry;
    ntree = (int)mxGetScalar(prhs[3]);
    mtry =  (int)mxGetScalar(prhs[4]);
        
    if (DEBUG_ON) mexPrintf("\nntree %d, mtry=%d, nsample=%d\n",ntree,mtry, nsample );
    
    int nt=ntree;
    int ipi=*((int*)mxGetData(prhs[10])); // ipi:      0=use class proportion as prob.; 1=use supplied priors
    plhs[10] = mxCreateDoubleScalar(mtry);
    plhs[3] = mxCreateNumericMatrix(nclass, 1, mxDOUBLE_CLASS, mxREAL);
    double *classwt = (double*) mxGetData(plhs[3]);
    _tmp_d = (double*) mxGetData(prhs[11]);
    
    //NOW COPY THE CLASSWT'S
    memcpy(classwt,_tmp_d,nclass*sizeof(double));
    
    
    plhs[4] = mxCreateNumericMatrix(nclass, 1, mxDOUBLE_CLASS, mxREAL);
    double *cutoff= (double*) mxGetData(plhs[4]);
    _tmp_d = (double*) mxGetData(prhs[12]);

    //NOW COPY THE CUTOFF's
    memcpy(cutoff,_tmp_d,nclass*sizeof(double));
    
    for(i=0;i<nclass;i++){
        //classwt[i]=1;
        //cutoff[i]=1/((double)nclass);
		if (DEBUG_ON){printf("%f,",cutoff[i]);}
    }
    int nodesize=*((int*) mxGetData(prhs[13]));
    
    plhs[11] = mxCreateNumericMatrix(nsample, 1, mxINT32_CLASS, mxREAL);
    int* outcl=(int*) mxGetData(plhs[11]); //calloc(nsample,sizeof(int));
    
    plhs[12] = mxCreateNumericMatrix(nclass, nsample, mxINT32_CLASS, mxREAL);
    int* counttr=(int*) mxGetData(plhs[12]); //calloc(nclass*nsample,sizeof(int));
    
    double* prox;
    if (proximity){
        plhs[13] = mxCreateNumericMatrix(n_size, n_size, mxDOUBLE_CLASS, mxREAL);
        prox = (double*) mxGetData(plhs[13]); //calloc(nsample*nsample,sizeof(double));
    }else{
        plhs[13] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        prox = (double*) mxGetData(plhs[13]); //calloc(1,sizeof(double));
        prox[0]=1;
    }
    double* impout;
    double* impmat; 
    double* impSD;
    
    if (localImp){
        if (addclass) {
            plhs[14] = mxCreateNumericMatrix(n_size * 2, p_size, mxDOUBLE_CLASS, mxREAL);
            impmat = (double*) mxGetData(plhs[14]); //calloc(n_size*p_size,sizeof(double));
        } else {
            plhs[14] = mxCreateNumericMatrix(n_size, p_size, mxDOUBLE_CLASS, mxREAL);
            impmat = (double*) mxGetData(plhs[14]); //calloc(n_size*p_size,sizeof(double));
        }
    }else{
        plhs[14] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
        impmat = (double*) mxGetData(plhs[14]); //calloc(1,sizeof(double));
        impmat[0]=1;
    }
        
    if (importance){    
        plhs[15] = mxCreateNumericMatrix(p_size,(nclass+2), mxDOUBLE_CLASS, mxREAL);
        plhs[16] = mxCreateNumericMatrix(p_size,(nclass+1), mxDOUBLE_CLASS, mxREAL);
        
        impout=(double*) mxGetData(plhs[15]); //calloc(p_size*(nclass+2),sizeof(double));
        impSD =(double*) mxGetData(plhs[16]); //calloc(p_size*(nclass+1),sizeof(double));
    }else{
        plhs[15] = mxCreateNumericMatrix(p_size,1, mxDOUBLE_CLASS, mxREAL);
        plhs[16] = mxCreateNumericMatrix(1,1, mxDOUBLE_CLASS, mxREAL);
        
        impout=(double*) mxGetData(plhs[15]); //calloc(p_size,sizeof(double));
        impSD =(double*) mxGetData(plhs[16]); //calloc(1,sizeof(double));
    }
    int nrnodes = 2 * (int)(nsum / nodesize) + 1;
    
    if (DEBUG_ON) { mexPrintf("\nnrnodes=%d, nsum=%d, nodesize=%d, mtry=%d\n",nrnodes,nsum,nodesize,mtry);}
    
    //int* ndbigtree = (int*) calloc(ntree,sizeof(int)); 
    plhs[9] = mxCreateNumericMatrix(1, nt, mxINT32_CLASS, mxREAL);
    int *ndbigtree = (int*) mxGetData(plhs[9]);
    
    //int* nodestatus = (int*) calloc(nt*nrnodes,sizeof(int));
    plhs[6] = mxCreateNumericMatrix(nrnodes, nt, mxINT32_CLASS, mxREAL);
    int *nodestatus = (int*) mxGetData(plhs[6]);

    //int* bestvar = (int*) calloc(nt*nrnodes,sizeof(int));
    plhs[8] = mxCreateNumericMatrix(nrnodes, nt, mxINT32_CLASS, mxREAL);
    int *bestvar = (int*) mxGetData(plhs[8]);
    
    //int* treemap = (int*) calloc(nt * 2 * nrnodes,sizeof(int));
    plhs[5] = mxCreateNumericMatrix(nrnodes, 2*nt, mxINT32_CLASS, mxREAL);
    int *treemap = (int*) mxGetData(plhs[5]);
    
    
    //int* nodepred = (int*) calloc(nt * nrnodes,sizeof(int));
    plhs[7] = mxCreateNumericMatrix(nrnodes, nt, mxINT32_CLASS, mxREAL);
    int *nodepred = (int*) mxGetData(plhs[7]);
    
    
    //double* xbestsplit = (double*) calloc(nt * nrnodes,sizeof(double));
    plhs[2] = mxCreateNumericMatrix(nrnodes, nt, mxDOUBLE_CLASS, mxREAL);
    double *xbestsplit = (double*) mxGetData(plhs[2]);

    plhs[17] = mxCreateNumericMatrix((nclass+1), ntree, mxDOUBLE_CLASS, mxREAL);
    double* errtr = (double*) mxGetData(plhs[17]); //calloc((nclass+1) * ntree,sizeof(double));
    
    int testdat=0;
    double* countts;
    int* inbag;
    if (keep_inbag){
        plhs[18] = mxCreateNumericMatrix(n_size, ntree, mxINT32_CLASS, mxREAL);
        inbag = (int*) mxGetData(plhs[18]);//calloc(n_size,sizeof(int));
    }else{
        plhs[18] = mxCreateNumericMatrix(n_size, 1, mxINT32_CLASS, mxREAL);
        inbag = (int*) mxGetData(plhs[18]);//calloc(n_size*ntree,sizeof(int));
    }
    if (DEBUG_ON){
        //printf few of the values
        //for(i=0;i<10;i++)
        //    mexPrintf("%f,",x[i]);
    }
    plhs[0] = mxCreateDoubleScalar(nrnodes);
    plhs[1] = mxCreateDoubleScalar(ntree);

    
    if (tst_available){
        xts = mxGetPr(prhs[19]);
        yts = (int*)mxGetData(prhs[20]);
        nts = (int)mxGetScalar(prhs[21]);
        plhs[19] = mxCreateNumericMatrix(nts, 1, mxINT32_CLASS, mxREAL);
        outclts = (int*)mxGetData(plhs[19]);
        countts = (double*) mxCalloc(nclass * nts,sizeof(double));
        if (proximity){
            plhs[20] = mxCreateNumericMatrix(nts, nts + n_size, mxDOUBLE_CLASS, mxREAL);
            proxts = (double*) mxGetData(plhs[20]); //calloc(nsample*nsample,sizeof(double));
        }else{
            plhs[20] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
            proxts = (double*) mxGetData(plhs[20]); //calloc(1,sizeof(double));
            proxts[0]=1;
        }
        plhs[21] = mxCreateNumericMatrix((nclass+1), ntree, mxDOUBLE_CLASS, mxREAL);
        errts = (double*) mxGetData(plhs[21]); //calloc((nclass+1) * ntree,sizeof(double));
        labelts=1;
        testdat=1;
    }else{
        xts = (double*)malloc(sizeof(double)*1);
        yts = (int*)malloc(sizeof(int)*1);
        nts = 1;
        plhs[19] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
        outclts = (int*)mxGetData(plhs[19]);
        countts = (double*) mxCalloc(nclass * nts,sizeof(double));
        if (proximity){
            plhs[20] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
            proxts = (double*) mxGetData(plhs[20]); //calloc(nsample*nsample,sizeof(double));
        }else{
            plhs[20] = mxCreateNumericMatrix(1, 1, mxDOUBLE_CLASS, mxREAL);
            proxts = (double*) mxGetData(plhs[20]); //calloc(1,sizeof(double));
            proxts[0]=1;
        }
        plhs[21] = mxCreateNumericMatrix((nclass+1), ntree, mxDOUBLE_CLASS, mxREAL);
        errts = (double*) mxGetData(plhs[21]); //calloc((nclass+1) * ntree,sizeof(double));
        labelts=0;
        testdat=0;
    }
     
	int print_verbose_tree_progression = (int)mxGetScalar(prhs[22]);

    //mexPrintf("tst available %d  %d\n",tst_available,nts);
    classRF(x, dimx, y, &nclass, cat, &maxcat,
	     sampsize, strata, Options, &ntree, &mtry,&ipi, 
         classwt, cutoff, &nodesize,outcl, counttr, prox,
	     impout, impSD, impmat, &nrnodes,ndbigtree, nodestatus, 
         bestvar, treemap,nodepred, xbestsplit, errtr,&testdat, 
         xts, yts, &nts, countts,outclts, labelts, 
         proxts, errts,inbag,print_verbose_tree_progression);
    
            
    mxFree(countts);
    
    //free(outcl);
    //free(counttr);
    //free(errtr);
    //free(countts);
    //free(xts);
    //free(yts);
    //free(inbag);
    //free(impout);
    //free(impmat);
    //free(impSD);
    //free(prox);
    
    // Below are allocated via matlab and will be needed for prediction
    //free(classwt);
    //free(cutoff);
    //free(ndbigtree);
    //free(nodestatus);
    //free(bestvar);
    //free(treemap);
    //free(nodepred);
    //free(xbestsplit);
    
    
}
