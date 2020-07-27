#ifndef FIT_VPEAK
#define FIT_VPEAK


int      n_bins;			//number of vpeak bins
double  *vp_binned;			//vpeak of each bin
double  *vp_fnc_binned;		//galaxies in each vpeak bin
double  *vp_fnc_err_binned;	//uncertainty on vpeak function


/*! \fn double vpeak_function(double vpeak, void *p)
 *  \brief Peak velocity function. */
double v_peak_function(double v_peak, double log10_A, double v_cut, double alpha, double beta);

/*! \fn void Read_V_Peak_Data(char fname[])
 *  \brief Read in the vpeak function info */
void Read_V_Peak_Data(char fname[]);

#endif //FIT_VPEAK
