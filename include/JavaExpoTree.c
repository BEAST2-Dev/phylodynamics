#include <stdlib.h>
#include <math.h>
#include <jni.h>
#include "beast_phylodynamics_ExactBDSIS.h"
#include "expoTree.h"



JNIEXPORT void JNICALL Java_beast_evolution_speciation_DummySIS_SayHello(JNIEnv *env, jobject thisObj) {

	printf("Hello from C!\n");
} 


/* 

	TODO: - change printf to System.out.println in Java, or make it print to Java error stream
		  - what about INFINITY?  What is INIFINITY in Java?
		  - check mapping between j<...> and C-types, should be platform independent


	INPUT: *env, thisObj - Java environment and object
		   K_in			 - Population sizes
		   beta_in 		 - Infection rates
		   mu_in		 - Death rates
		   psi_in		 - Sampling rates through time
		   rho 			 - Sampling ratio
		   times_in		 - Branching event times
		   ttypes_in	 - Different types of nodes in the tree
		   survival		 - Survival probability or tree probability
		   extant 		 - Number of extant lineages at the root
		   rescaling     - Rescaling on/off
		   verbose       - Verbose model

*/
JNIEXPORT jdoubleArray JNICALL Java_beast_phylodynamics_ExactBDSIS_DirtyExpoTree(JNIEnv      *env, 
																				   jobject      thisObj, 
																			 	   jdoubleArray K_in, 
																			 	   jdoubleArray beta_in, 
																			 	   jdoubleArray mu_in,
																			 	   jdoubleArray psi_in, 
																			 	   jdouble      rho,
																				   jdoubleArray times_in,
																				   jintArray    ttypes_in,
																				   jboolean     survival,
																				   jint         extant, 
																				   jboolean	    rescaling, 
																				   jint	        verbose) 
{


	int i;
	jdouble zeroTol = 1e-12;
	jdouble inf     = INFINITY;

//	int parVecLen = 1;
//	int parVecCol = 1;
//  parVecLen = INTEGER(parDim)[0];
//  parVecCol = INTEGER(parDim)[1];
//  if (vf) Rprintf("nrow = %d, ncol = %d\n",parVecLen,parVecCol);

//  double* pars = NUMERIC_POINTER(parameters);
//  double* N    = pars;

  	jdouble* N       = (*env)->GetDoubleArrayElements(env, K_in, NULL);
	jdouble* beta    = (*env)->GetDoubleArrayElements(env, beta_in, NULL);
	jdouble* mu      = (*env)->GetDoubleArrayElements(env, mu_in, NULL);
	jdouble* psi     = (*env)->GetDoubleArrayElements(env, psi_in, NULL);

	jint Nlength     = (*env)->GetArrayLength(env, K_in);
	jint betalength  = (*env)->GetArrayLength(env, beta_in);
	jint mulength    = (*env)->GetArrayLength(env, mu_in);
	jint psilength   = (*env)->GetArrayLength(env, psi_in);

	jdouble* ptimes  = (*env)->GetDoubleArrayElements(env, times_in, NULL);
	jint*    pttypes = (*env)->GetIntArrayElements(env, ttypes_in, NULL);
	jint     nt      = (*env)->GetArrayLength(env, times_in);

	if (verbose) {
		printf("Verbose mode is on.\n");
		if (rescaling) printf("Rescaling is on.\n");
		printf("%d lineage(s) extant at the root.\n",extant);
		printf("%d switche(s)\n", Nlength);
	}

	/* convert N reals to integers */
	int Nmax = 0;
	for (i = 0; i < Nlength; ++i) {
		if (Nmax < ceil(N[i])) Nmax = (int) ceil(N[i]);
	}
	if (verbose) printf("Maximum dimension: N = %d\n",Nmax);



//	SEXP p;
//	PROTECT(p = NEW_NUMERIC(Nmax+1));  /* +1 for p(0) */
//	double* p0 = NUMERIC_POINTER(p);


	int SI = 1;
	jint     ki = 0;
	jdouble* p0 = (double*) malloc((Nmax+1)*sizeof(double));
	jdouble  t0 = 0.0;
	jdouble  scale = 0.0;


	int root = (ptimes[nt-1] == ptimes[nt-2]) ? 1 : 0;
	if (root && verbose) printf("Root correction required.\n");


	int maxExtant = 0;
	for (i = nt-1; i >= 0; --i) {
		switch (pttypes[i]) {
		  case 1:
		    ++extant;
		    break;
		  case 0:
		  case 2:
		  case 4:
		    --extant;
		    break;
		  default:
		    break;
		}
		if (maxExtant < extant) maxExtant = extant;
	}
	if (verbose) {
		printf("%d extant lineages at the present; %d maximum.\n",
	        	extant,maxExtant);
	}




	int goodParams = 1;
	if (betalength != Nlength  ||
		betalength != mulength || 
		betalength != psilength) goodParams = 0;

	for (i = 0; i < betalength && goodParams; ++i) {
		if (beta[i] <= 0.0 || mu[i] < 0.0 || psi[i] < 0.0) goodParams = 0;
	}

	if (! goodParams || rho < 0.0 || rho > 1.0) {
		if (verbose) {
		  printf("Illegal parameters:\n");
		  printf("N\tbeta\tmu\tpsi\trho\n");
		  for (i = 0; i < betalength; ++i) {
		    printf("%f %f %f %f %f\n",N[i],beta[i],mu[i],psi[i],rho);
		  }
		}
		for (i = 0; i <= Nmax; ++i) p0[i] = -inf;
	} else {

		int vf = (int) verbose;
		int rs = (int) rescaling;

		if (survival) {
			p0[0] = 1.0;
			for (i = 1; i <= N[0]; ++i) {
				p0[i] = (rho <= 0.0 || rho >= 1.0) ? 0.0 : pow(1.-rho,i);
			}
			for (i = N[0]+1; i <= Nmax; ++i) p0[i] = 0.0;

			// Actually run rExpotree	
	  		rExpoTree(N,&ki,beta,mu,psi,&nt,&betalength,ptimes,pttypes,p0,&t0,&SI,&vf,&rs);
		  
			for (i = 0; i <= N[Nlength-1]; ++i) {
				if (p0[i] >= zeroTol) {
				  printf("Log survival probability is non-negative! p(%d) = %g\n",i,p0[i]);
				}
				p0[i] = (p0[i] < 0.0) ? log(1.-exp(p0[i])) : -inf;
			}

		} else {

			// set initial value of p
			if (extant == 0) {
				p0[0] = 0.0;
				for (i = 1; i <= N[0]; ++i) p0[i] = psi[0];
				ki = 1;
				t0 = ptimes[0];

				// what does this do? -- cant have a rho, move times to present, set prob = psi
				ptimes = ptimes+1;
				pttypes = pttypes+1;
				--nt;
  			} else {
				ki = extant;
				p0[0] = 0.0;
				scale = extant*log(rho);
				for (i = 1; i <= N[0]; ++i) {
				  	if (i < extant) 
				  		p0[i] = 0.0;
				  	else 
				  		p0[i] = pow(1.-rho,i-extant);
				}
			}
			for (i = N[0]+1; i <= Nmax; ++i) p0[i] = 0.0;

			// Actually run rExpotree	
			rExpoTree(N,&ki,beta,mu,psi,&nt,&betalength,ptimes,pttypes,p0,&t0,&SI,&vf,&rs);

			for (i = 0; i <= N[Nlength-1]; ++i) {
				p0[i] += scale;
				if (root) p0[i] -= M_LN2 + log(beta[betalength-1]) + log(1.-1./N[Nlength-1]);
			}
    	}
  	}

  	// What exactly should be returned?  p0[1]?
	jdoubleArray p_out = (*env)->NewDoubleArray(env, Nmax+1);  
	if (p_out == NULL) return NULL;
	(*env)->SetDoubleArrayRegion(env, p_out, 0 , Nmax+1, p0); 

	return p_out;
}











