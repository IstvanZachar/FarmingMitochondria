// Mitochondrial origins
// individual based model
// István Zachar
// 2017


#define _USE_MATH_DEFINES // Force the usage of math.h defined values

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <unistd.h> // using: `access`
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include "randomGenerator.c"


#define VERSION   (18)
#define CONFIG    ("parameters.cfg")
#define OUTPUT    ("data.dat")
#define HARVEST   (linearSwitch)
#define SKIP      (fscanf(file, "%*[^\n]\n")) // skip rest of line to next line

// If M_PI is not defined in math.h, define it here
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif


// Globals
int 		ps;
//int			periodLen;
int     *size    = NULL; // cell size vector
int     *farm    = NULL; // farm size vector
double  *harvest = NULL; // culling ("harvest") rate vector
double  *farming = NULL; // farm allocation ("farming") rate vector;
double  *cost    = NULL; // cost or benefit of farming
double  *periodRandoms = NULL;
gsl_rng *gauss;


// Functions

int randomSplit(int n) { // random variable from a binomial distribution (between 0 and `n`)
	int i, k = 0;
	for(i = 1; i <= n; i++) if(randd() < 0.5) k++;
	return(k);
}

double mean(int *v) {
	int i = 0;
	double sum = 0;
	for(i = 0; i < ps; i++) sum += v[i];
	return(sum/(double)ps);
}

double sd(int *v) {
	int i = 0;
	double sum = 0, m = mean(v);
	for(i = 0; i < ps; i++) sum += pow(v[i]-m, 2);
	return(sqrt(sum/(double)(ps-1)));
}

double meanf(double *v) {
	int i = 0;
	double sum = 0;
	for(i = 0; i < ps; i++) sum += v[i];
	return(sum/(double)ps);
}

double sdf(double *v) {
	int i = 0;
	double sum = 0, m = meanf(v);
	for(i = 0; i < ps; i++) sum += pow(v[i]-m, 2);
	return(sqrt(sum/(double)(ps-1)));
}

int clip(int x, int min, int max) {
	if (x > max) return(max); else if (x < min) return(min); else return(x);
}

double clipf(double x, double min, double max) {
	if (x > max) return(max); else if (x < min) return(min); else return(x);
}

double linearSwitch(double r, double s) {
	if (s <= 0.5) return(2*s*r); else if ((s >= 1) || (r/2 > (1 - s))) return(1); else return(r/(2 - (2*s)));
}

double trapezoidWave(int t, int v, double bl, double gl) {
	double i, b = v*bl, g = v*gl;
	int m = t % v;
	i = (v-b-g)/2;
	if (m < b) return(0); else if (m < (b+i)) return((m-b)/i); else if (m < (b+i+g)) return(1); else return(1 - ((m-b-i-g)/i));
}

double randomTrapezoidWave(int t, int v, double bl, double gl) { 
	double i, b, g;
	int m = t % v, p = (t/v); // integer index of period
	b = v*bl*periodRandoms[p*2];
	g = v*gl*periodRandoms[(p*2)+1];
	i = (v-b-g)/2;
	if (m < b) return(0); else if (m < (b+i)) return((m-b)/i); else if (m < (b+i+g)) return(1); else return(1-((m-b-i-g)/i));
}

double logisticWave(int t, int v, int s, int a) {
	if (v == 0) return(1); else return(2/(1. + exp(-(s*pow(cos(M_PI*(t/(double)v)), 2*a))/M_PI)) - 1);
}

double boxWave(int t, int v, double bl) {
	int m = t % v, b, p = (t/v); // integer index of period
	b = v*bl*periodRandoms[p];
	if (m < b) return(0); else return(1);
}




// Main

int main(int argc, char** argv) {
	
	long int randomseed;
	int    tmax, smax, amax, csmax, fsmax;
	int    initSize, initFarm;
	double initSizeSD, initFarmSD, initHarvest, initHarvestSD, initFarming, initFarmingSD, initCost, initCostSD;
	double basecost, maxbenefit, farmgrowth, farmingmut, harvestmut, costmut, farmingSD, harvestSD, costSD;
	double rmin, rmax, rrange, poorLen, richLen;
  double periodSD, r = 0, ra, rb, h;
	int    i, k, t, a = 0, saveAt, periodLen, pn, rfun, lfSteep, lfSharp;
	int    split = 0, growA = 0, growB = 0, farmA = 0, pay = 0, benefit = 0, hSum = 0; // event counters
	FILE   *file = NULL;
	
	
	{ // Read parameters
		if (access(CONFIG, F_OK) != -1) {
			file = fopen(CONFIG, "r");
			
			fscanf(file, "%lu",  &randomseed); SKIP; // `signed long`: `seed` only accepts positive `long`, does not accept 0
			fscanf(file, "%d",   &tmax);       SKIP; // maximal number of time steps
			fscanf(file, "%d",   &smax);       SKIP; // maximal number of stored steps
			fscanf(file, "%d",   &ps);         SKIP; // population size; minimum is 2.
			fscanf(file, "%d",   &amax);       SKIP; // maximum of percieved external prey
			fscanf(file, "%d",   &csmax);      SKIP; // maximal cell size
			fscanf(file, "%d",   &fsmax);      SKIP; // maximal farm size
			
			fscanf(file, "%lf",  &basecost);   SKIP; // base cost of maintaining an apparatus for farming (0, 1)
			fscanf(file, "%lf",  &maxbenefit); SKIP; // maximum value of benefit of the farm
			fscanf(file, "%lf",  &farmgrowth); SKIP; // at each step, farm f is increased (0, ...)
			fscanf(file, "%lf",  &farmingmut); SKIP; // mutation rate of farm allocation trait (0, 1)
			fscanf(file, "%lf",  &harvestmut); SKIP; // mutation rate of harvest probability trait (0, 1)
			fscanf(file, "%lf",  &costmut);    SKIP; // mutation rate of farm cost trait (-1, 1)
			fscanf(file, "%lf",  &farmingSD);  SKIP; // standard deviation of gaussian mutation of farming allocation trait (0, 1)
			fscanf(file, "%lf",  &harvestSD);  SKIP; // standard deviation of gaussian mutation of harvest probability trait (0, 1)
			fscanf(file, "%lf",  &costSD);     SKIP; // standard deviation of gaussian mutation of farm cost trait (0, 1)
			
			fscanf(file, "%d",   &rfun);  		 SKIP; // resource function: 0 = randomTrapezoidWave; 1 = trapezoidWave; 2 = logisticWave
			fscanf(file, "%d",   &periodLen);  SKIP; // integer period length of poor + rich subperiods + intermediate states
			fscanf(file, "%lf",  &periodSD);   SKIP; // period length SD
			fscanf(file, "%lf",  &rmin);       SKIP; // resource function minimum, function is scaled between `rmin` and `rmax`
			fscanf(file, "%lf",  &rmax);       SKIP; // resource function maximum
			fscanf(file, "%lf",  &poorLen);    SKIP; // length of poor period as a fraction of `periodLen`
			fscanf(file, "%lf",  &richLen);    SKIP; // length of rich period as a fraction of `periodLen`
			fscanf(file, "%d",   &lfSteep);    SKIP; // steepness of logistic wave function
			fscanf(file, "%d",   &lfSharp);    SKIP; // sharpness of logistic wave function
			
			fscanf(file, "%d",   &initSize);      SKIP; // initial cell size mean
			fscanf(file, "%lf",  &initSizeSD);    SKIP; // initial cell size SD
			fscanf(file, "%d",   &initFarm);      SKIP; // initial farm size mean
			fscanf(file, "%lf",  &initFarmSD);    SKIP; // initial farm size SD
			fscanf(file, "%lf",  &initHarvest);   SKIP; // initial harvest rate mean
			fscanf(file, "%lf",  &initHarvestSD); SKIP; // initial harvest rate SD
			fscanf(file, "%lf",  &initFarming);   SKIP; // initial farming rate mean
			fscanf(file, "%lf",  &initFarmingSD); SKIP; // initial farming rate SD
			fscanf(file, "%lf",  &initCost);      SKIP; // initial cost mean
			fscanf(file, "%lf",  &initCostSD);    SKIP; // initial cost SD
			fclose(file);
			
			if (0) {
				printf("randomseed: %lu\n",  randomseed);
				printf("tmax:       %d\n",   tmax);
				printf("smax:       %d\n",   smax);
				printf("ps:         %d\n",   ps);
				printf("amax:       %d\n",   amax);
				printf("csmax:      %d\n",   csmax);
				printf("fsmax:      %d\n",   fsmax);
				printf("basecost:   %lf\n",  basecost);
				printf("maxbenefit: %lf\n",  maxbenefit);
				printf("farmgrowth: %lf\n",  farmgrowth);
				printf("farmingmut: %lf\n",  farmingmut);
				printf("harvestmut: %lf\n",  harvestmut);
				printf("costmut:    %lf\n",  costmut);
				printf("farmingSD:  %lf\n",  farmingSD);
				printf("harvestSD:  %lf\n",  harvestSD);
				printf("costSD:     %lf\n",  costSD);
				printf("rfun:       %d\n",   rfun);
				printf("periodLen:  %d\n",   periodLen);
				printf("periodSD:   %lf\n",  periodSD);
				printf("rmin:       %lf\n",  rmin);
				printf("rmax:       %lf\n",  rmax);
				printf("poorLen:    %lf\n",  poorLen);
				printf("richLen:    %lf\n",  richLen);
				printf("lfSteep:    %d\n",   lfSteep);
				printf("lfSharp:    %d\n",   lfSharp);
				
			  printf("initSize:      %d\n",  initSize);
			  printf("initSizeSD:    %lf\n", initSizeSD);
			  printf("initFarm:      %d\n",  initFarm);
			  printf("initFarmSD:    %lf\n", initFarmSD);
			  printf("initHarvest:   %lf\n", initHarvest);
			  printf("initHarvestSD: %lf\n", initHarvestSD);
		  	printf("initFarming:   %lf\n", initFarming);
		  	printf("initFarmingSD: %lf\n", initFarmingSD);
		  	printf("initCost:      %lf\n", initCost);
		  	printf("initCostSD:    %lf\n", initCostSD);
			}
			
		} else { printf("ERROR. Parameter file %s does not exist. Aborting.\n", CONFIG); exit(0); }

	}
		
	
	{ // Initialization
		seed(randomseed);
		gauss = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(gauss, randomseed);
		
		if(smax <= 2) saveAt = tmax; else saveAt = clip(tmax/(smax-1), 1, tmax); // save at every `saveAt` timestep; integer
		pn = (tmax/periodLen)+1; // number of periods
		
		rrange = rmax-rmin;
		
		// initial vectors
		if (ps < 2) {	printf("ERROR. Population size %d is too small. Aborting.\n", ps); exit(0); }
		
		size    = (int*)calloc(ps, sizeof(int));
		farm    = (int*)calloc(ps, sizeof(int));
		harvest = (double*)calloc(ps, sizeof(double));
		farming = (double*)calloc(ps, sizeof(double));
		cost    = (double*)calloc(ps, sizeof(double));
		for(i = 0; i < ps; i++) {
			size[i]    = clip(initSize     + gsl_ran_gaussian(gauss, initSizeSD),    0, csmax);
			harvest[i] = clipf(initHarvest + gsl_ran_gaussian(gauss, initHarvestSD), 0, 1);
			farming[i] = clipf(initFarming + gsl_ran_gaussian(gauss, initFarmingSD), 0, 1);
			cost[i]    = clipf(initCost    + gsl_ran_gaussian(gauss, initCostSD),    -maxbenefit, maxbenefit);
			if ((farming[i] == 0) || (fsmax == 0) || (initFarmSD == 0)) {
				farm[i] = 0;
			} else {
				farm[i] = clip(initFarm + gsl_ran_gaussian(gauss, initFarmSD), 0, fsmax);
			}
		}
		
		if (0) {
			for(i = 0; i < ps; i++) printf("%d\t",  size[i]);    printf("\n");
			for(i = 0; i < ps; i++) printf("%d\t",  farm[i]);    printf("\n");
			for(i = 0; i < ps; i++) printf("%lf\t", harvest[i]); printf("\n");
			for(i = 0; i < ps; i++) printf("%lf\t", farming[i]); printf("\n");
			for(i = 0; i < ps; i++) printf("%lf\t", cost[i]);    printf("\n");
		}
		
		
		//random gaussian period length reals
		if (periodSD == 0) {
			periodRandoms = (double*)calloc(pn*2, sizeof(double));
			for(i = 0; i < (pn*2); i++) periodRandoms[i] = 1;
		} else {
			if (rfun == 0 || rfun == 1) { // 2 for each period
				periodRandoms = (double*)calloc(pn*2, sizeof(double));
				for(i = 0; i < (pn*2); i++) periodRandoms[i] = clipf(1 + gsl_ran_gaussian(gauss, periodSD), 0, 2);
			};
			if (rfun == 3) { // 1 for each period
				periodRandoms = (double*)calloc(pn, sizeof(double));
				for(i = 0; i < pn; i++) periodRandoms[i] = clipf(1 + gsl_ran_gaussian(gauss, periodSD), 0, 2);
			};
		}

		
		file = fopen(OUTPUT, "w");
	}
	
	
	
	
	{ // Main loop
		
		for(t = 0; t <= tmax; t++) {
			
			switch (rfun) {
				case 0: { r = randomTrapezoidWave(t, periodLen, poorLen, richLen); break; }
				case 1: { r = trapezoidWave(t, periodLen, poorLen, richLen); break; }
				case 2: { r = logisticWave(t, periodLen, lfSteep, lfSharp); break; }
				case 3: { r = boxWave(t, periodLen, poorLen); break; }
			};
			r = rmin + (rrange*r);
			a = round(amax*r);
			
			//FUSION

			if (((t % saveAt) == 0) || (t == tmax))  { // sava data
				int cP = 0;
			
				for(i = 0; i < ps; i++) if (farming[i] == 0) cP++; // count non-farmers

				//printf("%d \r", t);
				
				/*
				if (0) {
					printf("%d %d %d %d ", t, a, ps-cP, cP);
					printf("%d %d %d %d %d %d %d ", split, growA, growB, farmA, pay, benefit, hSum);
					printf("%f %f %f %f ", mean(size), sd(size), mean(farm), sd(farm));
					printf("%f %f %f %f %f %f\n", meanf(harvest), sdf(harvest), meanf(farming), sdf(farming), meanf(cost), sdf(cost));
				}
				*/
				
				fprintf(file, "%d %d %d %d ",          t, a, ps-cP, cP);
				fprintf(file, "%d %d %d %d %d %d %d ", split, growA, growB, farmA, pay, benefit, hSum);
				fprintf(file, "%f %f %f %f ",          mean(size), sd(size), mean(farm), sd(farm));
				fprintf(file, "%f %f %f %f %f %f\n",   meanf(harvest), sdf(harvest), meanf(farming), sdf(farming), meanf(cost), sdf(cost));
				fflush(file);
				
				split = 0;
				farmA = 0;
				growA = 0;
				growB = 0;
				pay   = 0;
				hSum  = 0.;
				benefit = 0;
			}
			
			
			for(k = 0; k < ps; k++) {
				i = randl(ps);
				if ((fsmax == 0) || (farm[i] == 0) || (farming[i] == 0)) {
					farm[i] = 0;
				} else {
					farm[i] = clip(farm[i] + (r * farm[i] * (farmgrowth - 1)), 0, fsmax);
				}
				ra = (double)a/(double)amax;
				rb = (double)farm[i]/(double)fsmax;
				if ((rb == 0) || (fsmax == 0) || ((ra+rb) == 0)) {
					h = 0;
				} else {
					h = HARVEST(rb/(ra+rb), harvest[i]);
				}
				hSum = hSum + h;
				
				if ((farming[i] > 0) && (randd() < abs(basecost))) { // explicit cost/benefit of apparatus.
					size[i] = clip(size[i] - 1, 0, csmax);
				}
				
				if ((farming[i] > 0) && (rb > 0)) { // farmsize-dependent cost/benefit of farm.
				// deterministically add/subtract as many units to `size` as the integer part of `cost[i]` is
				// probabilistically add/subtract 1 unit to `size` according to the absolute fractional part of `cost[i]`
					int base = (int)cost[i], frac; // `base` is signed
					if (cost[i] < 0) { 
						frac = -(randd() < -(rb*(cost[i] - base)));
						pay++;
					} else { 
						frac = (randd() < (rb*(cost[i] - base)));
						benefit++;
					}
					size[i] = clip(size[i] + base + frac, 0, csmax);
				}
				
				
				if (randd() < ra) { // gets A
					a = clip(a - 1, 0, amax);
					if ((fsmax > 0) && (randd() < farming[i]) && (randd() < (1 - rb))) { // allocate A to farm
						farm[i] = clip(farm[i] + 1, 0, fsmax);
						farmA++;
					} else { // eats A and grows
						size[i] = clip(size[i] + 1, 0, csmax);
						growA++;
					}
				} else { // doesn't get A
					// NOTE: That the previous harvest value (instead of actual) is used.
					if ((farm[i] > 0) && (randd() < h)) { // harvest farm and grow
						farm[i] = clip(farm[i] - 1, 0, fsmax);
						size[i] = clip(size[i] + 1, 0, csmax);
						growB++;
					}
				} // eat & grow
				
				// Split
				if (size[i] == csmax) {
					int j = i;
					split++;
					while(i == j) j = randl(ps); // select another individual than `i`
					size[i] = 0;
					size[j] = 0;
					if (farm[i] > 0) {
						int z = randomSplit(farm[i]);
						farm[j] = z;
						farm[i] = farm[i] - z;
					} else {
						farm[j] = 0;
					}
					farming[j] = farming[i];
					harvest[j] = harvest[i];
					cost[j]    = cost[i];
					if (randd() < farmingmut) farming[j] = clipf(farming[i] + gsl_ran_gaussian(gauss, farmingSD),  0, 1);
					if (randd() < harvestmut) harvest[j] = clipf(harvest[i] + gsl_ran_gaussian(gauss, harvestSD),  0, 1);
					if (randd() < costmut   ) cost[j]    = clipf(cost[i]    + gsl_ran_gaussian(gauss, costSD   ), -maxbenefit, maxbenefit);					
				} // split
				
			} // population
						
		} // time
		
	}
	
	
	fclose(file);
	free(periodRandoms);
	free(size);
	free(farm);
	free(harvest);
	free(farming);
	free(cost);			
	gsl_rng_free(gauss);
	
	//printf("Ready.\n");
	
	return(0);
}