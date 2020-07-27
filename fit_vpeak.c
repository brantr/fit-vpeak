#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fit_vpeak.h"
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <float.h>
#include "multinest.h"


//these are the flat priors
#define A_MIN -15.0
#define A_MAX 1.0
#define V_CUT_MIN 0
#define V_CUT_MAX 10.0
#define ALPHA_MIN 0.0
#define ALPHA_MAX 6.0
#define BETA_MIN -2.0
#define BETA_MAX  2.0



int      n_bins;			//number of magnitude bins
double  *vp_binned;			//magnitude of each bin
double  *vp_fnc_binned;		//number density of galaxies in each bin
double  *vp_fnc_err_binned;	//uncertainty on number density of galaxies in each bin


/*! \fn double vpeak_function(double vpeak, void *p)
 *  \brief Peak velocity function. */
double v_peak_function(double v_peak, double log10_A, double v_cut, double alpha, double beta)
{
	double x = pow(10,v_peak)/pow(10,v_cut);

    //x = 10**log10vpeak_arr/10**log10vcut
    //output = log10(10**log10A * (1+x)**(-beta) * exp(- x**alpha))
	return log10_A+log10(pow(1.+x,-1*beta)*exp(-1*pow(x,alpha)));
}

/*! \fn void Read_V_Peak_Data(char fname[])
 *  \brief Read in the vpeak function info */
void Read_V_Peak_Data(char fname[])
{
	//This subroutine would read in the galaxy
	//data, the binned LF vs magnitude.

	//instead, we use it to make some fake data
	//here for testing.

	FILE *fp = fopen(fname,"r");

	int i;	//dummy index

	n_bins = 20;

	//print the number of luminosity bins
	printf("n_bins %d\n",n_bins);

	//allocate arrays to hold the binned lf
	vp_binned = (double *) malloc(n_bins*sizeof(double));
	vp_fnc_binned = (double *) malloc(n_bins*sizeof(double));
	vp_fnc_err_binned = (double *) malloc(n_bins*sizeof(double));

	double zb;
	fscanf(fp,"# %f\n",&zb);
	printf("zb = %e\n",zb);
	for(i=0;i<n_bins;i++)
	{
		fscanf(fp,"%lf\t%lf\t%lf\n",&vp_binned[i],&vp_fnc_binned[i],&vp_fnc_err_binned[i]);
		printf("%e\t%e\t%e\n",vp_binned[i],vp_fnc_binned[i],vp_fnc_err_binned[i]);
		vp_fnc_err_binned[i] = 0.1;
	}
	fclose(fp);

}




void LogLike(double *Cube, int *ndim, int *npars, double *lnew, void *context)
{

	double l = 0;	//the log likelihood

	//mapping of Cube[0] to A prior range
	double A_min  = A_MIN;
	double A_max  = A_MAX;
	double log10_A = Cube[0]*(A_max - A_min) + A_min;

	//mapping of Cube[1] to v_cut
	double v_cut_min = V_CUT_MIN;
	double v_cut_max = V_CUT_MAX;
	double v_cut     = Cube[1]*(v_cut_max-v_cut_min) + v_cut_min;

	//mapping of Cube[2] to alpha prior range
	double alpha_min = ALPHA_MIN;
	double alpha_max = ALPHA_MAX;
	double alpha = Cube[2]*(alpha_max-alpha_min) + alpha_min;

	//mapping of Cube[3] to alpha prior range
	double beta_min = BETA_MIN;
	double beta_max = BETA_MAX;
	double beta = Cube[3]*(beta_max-beta_min) + beta_min;

	int i;
	double y;	//measured value
	double ym;	//model value

	l = 0.0;


	//find the likelihood 
	for(i=0;i<n_bins;i++)
	{
	    //find difference between model and data
			ym = v_peak_function(vp_binned[i],log10_A, v_cut, alpha, beta);
			y = (vp_fnc_binned[i] - ym)/vp_fnc_err_binned[i];

	    //add chi^2 values
	    l += y*y;
	}

	//set the cube values,
	//which tells the fitting routine
	//what values the likelihood were
	//computed for.
	Cube[0] = log10_A;
	Cube[1] = v_cut;
	Cube[2] = alpha;
	Cube[3] = beta;

	//printf("log10_A %e v_cut %e alpha %e beta %e\n",log10_A,v_cut,alpha,beta);

	//record log likelihood
	*lnew = -0.5 * l;
}



/***********************************************************************************************************************/




/************************************************* dumper routine ******************************************************/

// The dumper routine will be called every updInt*10 iterations
// MultiNest doesn not need to the user to do anything. User can use the arguments in whichever way he/she wants
//
//
// Arguments:
//
// nSamples 						= total number of samples in posterior distribution
// nlive 						= total number of live points
// nPar 						= total number of parameters (free + derived)
// physLive[1][nlive * (nPar + 1)] 			= 2D array containing the last set of live points (physical parameters plus derived parameters) along with their loglikelihood values
// posterior[1][nSamples * (nPar + 2)] 			= posterior distribution containing nSamples points. Each sample has nPar parameters (physical + derived) along with the their loglike value & posterior probability
// paramConstr[1][4*nPar]:
// paramConstr[0][0] to paramConstr[0][nPar - 1] 	= mean values of the parameters
// paramConstr[0][nPar] to paramConstr[0][2*nPar - 1] 	= standard deviation of the parameters
// paramConstr[0][nPar*2] to paramConstr[0][3*nPar - 1] = best-fit (maxlike) parameters
// paramConstr[0][nPar*4] to paramConstr[0][4*nPar - 1] = MAP (maximum-a-posteriori) parameters
// maxLogLike						= maximum loglikelihood value
// logZ							= log evidence value
// logZerr						= error on log evidence value
// context						void pointer, any additional information

void dumper(int *nSamples, int *nlive, int *nPar, double **physLive, double **posterior, double **paramConstr, double *maxLogLike, double *logZ, double *logZerr, void *context)
{
	// convert the 2D Fortran arrays to C arrays
	
	
	// the posterior distribution
	// postdist will have nPar parameters in the first nPar columns & loglike value & the posterior probability in the last two columns
	
	int i, j;
	
	double postdist[*nSamples][*nPar + 2];
	for( i = 0; i < *nPar + 2; i++ )
		for( j = 0; j < *nSamples; j++ )
			postdist[j][i] = posterior[0][i * (*nSamples) + j];
	
	
	
	// last set of live points
	// pLivePts will have nPar parameters in the first nPar columns & loglike value in the last column
	
	double pLivePts[*nlive][*nPar + 1];
	for( i = 0; i < *nPar + 1; i++ )
		for( j = 0; j < *nlive; j++ )
			pLivePts[j][i] = physLive[0][i * (*nlive) + j];
}

/***********************************************************************************************************************/




/************************************************** Main program *******************************************************/



int main(int argc, char *argv[])
{
	int i;
	
	// set the MultiNest sampling parameters
	
	
	//int mmodal = 1;					// do mode separation?
	int mmodal = 0;					// do mode separation?
	
	int ceff = 0;					// run in constant efficiency mode?
	
	int nlive = 2000;				// number of live points
	
	double efr = 1.0;				// set the required efficiency
	
	double tol = 0.5;				// tol, defines the stopping criteria
	
	int ndims = 4;					// dimensionality (no. of free parameters)
	int nPar = 4;					// total no. of parameters including free & derived parameters
	int nClsPar = 4;				// no. of parameters to do mode separation on


	
	int updInt = 100;				// after how many iterations feedback is required & the output files should be updated
							// note: posterior files are updated & dumper routine is called after every updInt*10 iterations
	
	double Ztol = -1E90;				// all the modes with logZ < Ztol are ignored
	
	int maxModes = 100;				// expected max no. of modes (used only for memory allocation)
	
	int pWrap[ndims];				// which parameters to have periodic boundary conditions?
	for(i = 0; i < ndims; i++) pWrap[i] = 0;
	
	char root[100] = "chains/vpeak.";		// root for output files
	
	int seed = -1;					// random no. generator seed, if < 0 then take the seed from system clock
	
	int fb = 1;					// need feedback on standard output?
	
	int resume = 0;					// resume from a previous job?
	
	int outfile = 1;				// write output files?
	
	int initMPI = 1;				// initialize MPI routines?, relevant only if compiling with MPI
							// set it to F if you want your main program to handle MPI initialization
	
	double logZero = -DBL_MAX;			// points with loglike < logZero will be ignored by MultiNest
	
	int maxiter = 0;				// max no. of iterations, a non-positive value means infinity. MultiNest will terminate if either it 
							// has done max no. of iterations or convergence criterion (defined through tol) has been satisfied
	
	void *context = 0;				// not required by MultiNest, any additional information user wants to pass


	char fname[200];	//file name for the galaxy population
	
	sprintf(fname,"HVF_500.txt");


	if(argc>1)
		sprintf(fname,argv[1]);

	printf("Reading data from %s....\n",fname);

	//read in the v_peak function data
	Read_V_Peak_Data(fname);


	//print v_peak data to screen
	for(i=0;i<n_bins;i++)
	    printf("i %d v %e fnc %e fnc_err %e\n",i,vp_binned[i],vp_fnc_binned[i],vp_fnc_err_binned[i]);
	fflush(stdout);

	//exit(0);
	// calling MultiNest
	run(mmodal, ceff, nlive, tol, efr, ndims, nPar, nClsPar, maxModes, updInt, Ztol, root, seed, pWrap, fb, resume, outfile, initMPI, 
	logZero, maxiter, LogLike, dumper, context);

}


/***********************************************************************************************************************/
