#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "inhom_v3.h"
#include "constants.h"
#include "pfr_cufft1.h"
#include "global_variables.h"

//#define DISLOCATION_ACTIVITY_ON
/*
#define Cijkla (400)
#define Cijklb (272)
//#define Cijklc (6700)
#define Cijklc (98)
*/
/*
#define Cijkla (200000)
#define Cijklb (136000)
//#define Cijklc (6700)
#define Cijklc (49000)
*/
#define Cijkla (18300)
#define Cijklb (14600)
//#define Cijklc (6700)
#define Cijklc (4600)
/*
#define Cijkla (1.83e6)
#define Cijklb (1.46e6)
//#define Cijklc (6700)
#define Cijklc (4.6e5)
*/
#define MISFIT_STRAIN (0.01) // Misfit strain

#define CURR_LINE() {printf("%s at %d:\n",__FILE__,__LINE__);}

double temp4_r[LMN];//temporary variable for fourier transform
double temp_r[LMN];//temporary variable for fourier transform
complx temp4_k[HLMN];//temporary variable for fourier transform
void cufft_fr(voigt *fr, voigtk *fk);
void cufft_fk(voigtk *fk, voigt *fr);
voigt te_p[V];
voigt T_tr[LMN];


voigt Tau_0[LMN],Tau[LMN];
voigt iEps0_last[LMN];

//voigt total_strain[LMN];//Output total strain
ten4th GAMMA[HLMN];//Green's operator
voigtk kSig[HLMN];//To store temporary variables
voigt rSig[LMN];//to store temporary variables
voigt vEps[LMN];//virtual strain

ten4th Sijkl0 = {{{{ 0.0 }}}};
voigt66 Sij0 = {{ 0.0 }};

void LU_dcmp(real **a, int n, int *indx, real *d);
void LU_bksb(real **a, int n, int *indx, real b[]);
void LU_inv_66(voigt66 c);
static inline void M3Inverse(real a[3][3], real ia[3][3]);
void CalcGAMMA(ten4th cijkl0, ten4th *gamma);
static inline double UpdateVirtualStrain(int opt);
void Set_rot(double rott[3][3]);

real normepsavg = 0.0;
int ElastBC=0;//stress controlled boundary condition,BC=1 stress control, bc=0 strain control
ten4th Cijkl0 = {{{{ 0.0 }}}};//Modulus of the matrix (B2 phase) C0
ten4th Cijkl1 = {{{{ 0.0 }}}};//Modulus of the ppt    (disordered phase) C1
ten4th Cijkl2 = {{{{ 0.0 }}}};    // the value of C0-C1
voigt66 Cij0 = {{ 0.0 }};
voigt66 Cij2 = {{ 0.0 }};
voigt EpsAvg,SigAvg;// Average strain and stress (external)
voigt Sig[LMN];//Output stress
voigt E_oo={MISFIT_STRAIN,MISFIT_STRAIN,MISFIT_STRAIN,0,0,0};
double rot[3][3];

void inhom_init()
{
		double f=0.0;
		double temp[LMN];
		char filename[500];
	// Macroscopic stress or strain
	V6_loop{
        EpsAvg[mi] = 0.0;
        SigAvg[mi] = 0.0;
    	}//----------------------------
/*        SigAvg[0] = 38.0936;
        SigAvg[1] = 59.5212;
        SigAvg[2] = 2.3808;
        SigAvg[3] = -11.9042/RSQ2;
        SigAvg[4] = 9.5234/RSQ2;
        SigAvg[5] = -47.6170/RSQ2;
*/
	//normepsavg used in convergence criterion
	V6_loop{
        normepsavg += EpsAvg[mi]*EpsAvg[mi];
    	}
    	normepsavg = sqrt(normepsavg);
	//--------------------------------------
real c11,c12,c44,c33,c13;
c11=c33=Cijkla;
c12=c13=Cijklb;
c44=Cijklc;
      Cijkl0[0][0][0][0] =c11;	Cijkl1[0][0][0][0] =f*c11;
      Cijkl0[1][1][1][1] =c11;	Cijkl1[1][1][1][1] =f*c11;
      Cijkl0[2][2][2][2] =c33;	Cijkl1[2][2][2][2] =f*c33;
      Cijkl0[0][0][1][1] =c12;	Cijkl1[0][0][1][1] =f*c12;
      Cijkl0[1][1][2][2] =c13;	Cijkl1[1][1][2][2] =f*c13;
      Cijkl0[2][2][0][0] =c13;	Cijkl1[2][2][0][0] =f*c13;
      Cijkl0[0][1][0][1] =c44;	Cijkl1[0][1][0][1] =f*c44;
      Cijkl0[1][2][1][2] =c44;	Cijkl1[1][2][1][2] =f*c44;
      Cijkl0[2][0][2][0] =c44;	Cijkl1[2][0][2][0] =f*c44;
      Cijkl0[1][1][0][0] =c12;	Cijkl1[1][1][0][0] =f*c12;
      Cijkl0[2][2][1][1] =c13;	Cijkl1[2][2][1][1] =f*c13;
      Cijkl0[0][0][2][2] =c13;	Cijkl1[0][0][2][2] =f*c13;
      Cijkl0[1][0][1][0] =c44;	Cijkl1[1][0][1][0] =f*c44;
      Cijkl0[0][1][1][0] =c44;	Cijkl1[0][1][1][0] =f*c44;
      Cijkl0[1][0][0][1] =c44;	Cijkl1[1][0][0][1] =f*c44;
      Cijkl0[2][1][2][1] =c44;	Cijkl1[2][1][2][1] =f*c44;
      Cijkl0[2][1][1][2] =c44;	Cijkl1[2][1][1][2] =f*c44;
      Cijkl0[1][2][2][1] =c44;	Cijkl1[1][2][2][1] =f*c44;
      Cijkl0[0][2][0][2] =c44;	Cijkl1[0][2][0][2] =f*c44;
      Cijkl0[0][2][2][0] =c44;	Cijkl1[0][2][2][0] =f*c44;
      Cijkl0[2][0][0][2] =c44;	Cijkl1[2][0][0][2] =f*c44;
//v1
/*
te_p[0][0]=-0.0065;
te_p[0][1]=0.0263;
te_p[0][2]=0.0263;
te_p[0][3]=0.1340;
te_p[0][4]=-0.0971;
te_p[0][5]=-0.0971;
//v2
te_p[1][0]=-0.0065;
te_p[1][1]=0.0263;
te_p[1][2]=0.0263;
te_p[1][3]=0.1340;
te_p[1][4]=0.0971;
te_p[1][5]=0.0971;
//v3
te_p[2][0]=-0.0065;
te_p[2][1]=0.0263;
te_p[2][2]=0.0263;
te_p[2][3]=-0.1340;
te_p[2][4]=-0.0971;
te_p[2][5]=0.0971;
//v4
te_p[3][0]=-0.0065;
te_p[3][1]=0.0263;
te_p[3][2]=0.0263;
te_p[3][3]=-0.1340;
te_p[3][4]=0.0971;
te_p[3][5]=-0.0971;
*/
//v1
te_p[0][0]=-0.0381;
te_p[0][1]=0.0289;
te_p[0][2]=0.0289;
te_p[0][3]=0.0916;
te_p[0][4]=-0.0709;
te_p[0][5]=-0.0709;
//v2
te_p[1][0]=-0.0381;
te_p[1][1]=0.0289;
te_p[1][2]=0.0289;
te_p[1][3]=0.0916;
te_p[1][4]=0.0709;
te_p[1][5]=0.0709;
//v3
te_p[2][0]=-0.0381;
te_p[2][1]=0.0289;
te_p[2][2]=0.0289;
te_p[2][3]=-0.0916;
te_p[2][4]=-0.0709;
te_p[2][5]=0.0709;
//v4
te_p[3][0]=-0.0381;
te_p[3][1]=0.0289;
te_p[3][2]=0.0289;
te_p[3][3]=-0.0916;
te_p[3][4]=0.0709;
te_p[3][5]=-0.0709;

//v5
/*
te_p[2][0]=0.0263;
te_p[2][1]=-0.0065;
te_p[2][2]=0.0263;
te_p[2][3]=-0.0971;
te_p[2][4]=0.1340;
//v5

te_p[2][0]=0.0263;
te_p[2][1]=-0.0065;
te_p[2][2]=0.0263;
te_p[2][3]=-0.0971;
te_p[2][4]=0.1340;
te_p[2][5]=-0.0971;
*/
/*
te_p[2][0]=0.0263;
te_p[2][1]=0.0263;
te_p[2][2]=-0.0065;
te_p[2][3]=-0.0971;
te_p[2][4]=-0.0971;
te_p[2][5]=0.1340;
*/
//v7
/*
te_p[3][0]=0.0263;
te_p[3][1]=-0.0065;
te_p[3][2]=0.0263;
te_p[3][3]=-0.0971;
te_p[3][4]=-0.1340;
te_p[3][5]=0.0971;
*/
//v11
/*
te_p[2][0]=0.0263;
te_p[2][1]=0.0263;
te_p[2][2]=-0.0065;
te_p[2][3]=0.0971;
te_p[2][4]=-0.0971;
te_p[2][5]=-0.1340;
*/
//v12
/*
te_p[3][0]=0.0263;
te_p[3][1]=0.0263;
te_p[3][2]=-0.0065;
te_p[3][3]=-0.0971;
te_p[3][4]=0.0971;
te_p[3][5]=-0.1340;
*/

te_p[4][0]=0.02291;
te_p[4][1]=0.00195;
te_p[4][2]=0.00195;
te_p[4][3]=0.0103;
te_p[4][4]=0;
te_p[4][5]=0;

te_p[5][0]=0.02291;
te_p[5][1]=0.00195;
te_p[5][2]=0.00195;
te_p[5][3]=-0.0103;
te_p[5][4]=0;
te_p[5][5]=0;

te_p[6][0]=0.00195;
te_p[6][1]=0.02291;
te_p[6][2]=0.00195;
te_p[6][3]=0.0;
te_p[6][4]=0.0103;
te_p[6][5]=0;

te_p[7][0]=0.00195;
te_p[7][1]=0.02291;
te_p[7][2]=0.00195;
te_p[7][3]=0.0;
te_p[7][4]=-0.0103;
te_p[7][5]=0;

te_p[8][0]=0.00195;
te_p[8][1]=0.00195;
te_p[8][2]=0.02291;
te_p[8][3]=0.0;
te_p[8][4]=0.0;
te_p[8][5]=0.0103;

te_p[9][0]=0.00195;
te_p[9][1]=0.00195;
te_p[9][2]=0.02291;
te_p[9][3]=0.0;
te_p[9][4]=0.0;
te_p[9][5]=-0.0103;

/*
		rot[0][0]=-0.6172;
		rot[0][1]=-0.5345;
		rot[0][2]=0.5774;
		rot[1][0]=0.7715;
		rot[1][1]=-0.2673;
		rot[1][2]=0.5774;
		rot[2][0]=-0.1543;
		rot[2][1]=0.8018;
		rot[2][2]=0.5774;
*/

                rot[0][0]=1.0;
		rot[0][1]=0.0;
		rot[0][2]=0.0;
		rot[1][0]=0.0;
		rot[1][1]=1.0;
		rot[1][2]=0.0;
		rot[2][0]=-0.0;
		rot[2][1]=0.0;
		rot[2][2]=1.0;

Set_rot(rot);
	C4_loop{Cijkl2[mi][mj][mk][ml]=Cijkl0[mi][mj][mk][ml]-Cijkl1[mi][mj][mk][ml];}
	chg_basis_Kelvin_4(Cijkl0,Cij0);
	chg_basis_Kelvin_4(Cijkl2,Cij2);
//Computing Sij0
C6_loop{
Sij0[mi][mj] = Cij0[mi][mj];
}
LU_inv_66(Sij0);

CURR_LINE();
/*
te_p[0][0]=-0.0065;
te_p[0][1]=0.0263;
te_p[0][2]=0.0263;
te_p[0][3]=0.1339;
te_p[0][4]=-0.0972;
te_p[0][5]=-0.0972;

te_p[1][0]=0.0263;
te_p[1][1]=-0.0065;
te_p[1][2]=0.0263;
te_p[1][3]=-0.0972;
te_p[1][4]=0.1339;
te_p[1][5]=-0.0972;

te_p[2][0]=0.0263;
te_p[2][1]=0.0263;
te_p[2][2]=-0.0065;
te_p[2][3]=-0.0972;
te_p[2][4]=-0.0972;
te_p[2][5]=0.1339;

te_p[3][0]=-0.0065;
te_p[3][1]=0.0263;
te_p[3][2]=0.0263;
te_p[3][3]=0.1339;
te_p[3][4]=0.0972;
te_p[3][5]=0.0972;
*/
/*
te_p[4][0]=-0.0065;
te_p[4][1]=0.0263;
te_p[4][2]=0.0263;
te_p[4][3]=0.1339;
te_p[4][4]=0.0972;
te_p[4][5]=0.0972;

te_p[5][0]=-0.0065;
te_p[5][1]=0.0263;
te_p[5][2]=0.0263;
te_p[5][3]=0.1339;
te_p[5][4]=0.0972;
te_p[5][5]=0.0972;
*/
/*
te_p[0][0]=0.0224;
te_p[0][1]=0.0024;
te_p[0][2]=0.0024;
te_p[0][3]=0.0051;
te_p[0][4]=0.0;
te_p[0][5]=0.0;

te_p[1][0]=0.0224;
te_p[1][1]=0.0024;
te_p[1][2]=0.0024;
te_p[1][3]=-0.0051;
te_p[1][4]=0.0;
te_p[1][5]=0.0;

te_p[2][0]=0.0024;
te_p[2][1]=0.0224;
te_p[2][2]=0.0024;
te_p[2][3]=0.0;
te_p[2][4]=0.0051;
te_p[2][5]=0.0;

te_p[3][0]=0.0024;
te_p[3][1]=0.0224;
te_p[3][2]=0.0024;
te_p[3][3]=0.0;
te_p[3][4]=-0.0051;
te_p[3][5]=0.0;
*/
#pragma acc enter data create(GAMMA)
#pragma acc enter data copyin(Sij0,Cij0,Cijkl0,EpsAvg,SigAvg)
#pragma acc enter data create(Sig)
#pragma acc enter data create(iEps0_last,vEps,kSig,Tau_0,Tau,rSig)
#pragma acc enter data copyin(Cij2,E_oo,te_p)
	CalcGAMMA(Cijkl0, GAMMA);	// calculate Green operator in k-space
#pragma acc data present(vEps,te_p,eta_r)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
           // vEps[pIDX][mi] = heta_r[pxyz]*E_oo[mi];       // Initial guess of virtual strain field
	for(int i=0;i<V;i++)
	{
            vEps[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // Initial guess of virtual strain field
	}
	}
}


#pragma acc enter data create(T_tr)
#pragma acc data present(T_tr,te_p,eta_r)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
           // vEps[pIDX][mi] = heta_r[pxyz]*E_oo[mi];       // Initial guess of virtual strain field
        T_tr[pIDX][mi]=0;
	for(int i=0;i<V;i++)
	{
            T_tr[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // Initial guess of virtual strain field
	}
	}
}

//#pragma acc data present(T_tr,vEps)
/*
#pragma acc enter data create(temp_r)
#pragma acc data present(T_tr)
#pragma acc parallel loop collapse(3) 
x_loop{ 
  temp_r[pIDX]=T_tr[pIDX][0];
}
#pragma acc update host(temp_r)
			sprintf(filename,"output_files/T_tr11.vtk");
			write_vtk(temp_r,filename,"vEps11");
*/
CURR_LINE();

}//end of inhom_init()

void inhom_finish()
{
#pragma acc exit data delete(GAMMA)
#pragma acc exit data delete(Sij0,Cij0,EpsAvg,SigAvg)
#pragma acc exit data copyout(Sig)
#pragma acc exit data delete(iEps0_last,vEps,kSig,Tau_0,Tau,rSig,T_tr)
#pragma acc exit data delete(Cij2,E_oo)
}//end of inhom_finish()

void inhom_solver()
{
	int step = 0;// to store the number of iterations
	  real convg=0.0;
	  char filename[500];
	  static int num=0;
	  num++;
 	//starting the iterations
 	 // printf("entering inhom_solver");
    while(step<MAXITER){
        step++;
        convg = UpdateVirtualStrain(0);
#ifdef DEBUG_ON
			sprintf(temp,"vEps_t%d.vtk",step);
			write_vtk_voigt(vEps,2,temp,"vEps");
#endif
		if(fabs(convg) <(ERROR_REL_ELAST*normepsavg)||fabs(convg)<ERROR_ABS_ELAST){
			break;
		}
    }
	//Calculating the stress field
	UpdateVirtualStrain(1);
      if(num<=total_step/2) 
       { SigAvg[0] += s_app_inc;}
	else
       { SigAvg[0] -= s_app_inc;}
/*	if(s_app_inc!=0){
	printf("Applied stress = %f ,number =%d \n",SigAvg[0],num);
	}
*/
//			sprintf(temp,"Stress_t%d.vtk",step);
//			write_vtk(Sig_out,temp,"vEps");
//update transformatin strain
/*
#pragma acc data present(te_p,eta_r,T_tr)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
	for(int i=0;i<V;i++)
	{
            T_tr[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // update transformation strain
	}
	}
}
*/
/*
#pragma acc enter data create(temp_r)
#pragma acc data present(T_tr)
#pragma acc parallel loop collapse(3) 
x_loop{ 
  temp_r[pIDX]=T_tr[pIDX][0];
}
#pragma acc update host(temp_r)
			sprintf(filename,"output_files/T_tr11_%d.vtk",num);
			write_vtk(temp_r,filename,"vEps11");

*/
}//end of inhom_solver()


void LU_dcmp(real **a, int n, int *indx, real *d)
{
	/* Adopted from "Numerical Recipeis in C" by W Press et. al.

	   Given a matrix a[1..n][1..n], this function REPLACE it by
	   the LU decomposition of a rowwise permutation of itself.
INPUT:
	a -- original matrix
	n -- dimension
OUTPUT:
	a -- matrix containing L and U
	indx -- row permuation effected by the partial pivoting
	d -- +1/-1, depending on the number of row interchanges was even or odd, respectively. */

	int i,imax,j,k;
	real big,dum,sum,temp;
	real *vv; //vv stores the implicit scaling of each row.

	vv=(real*)malloc(n*sizeof(real));
	*d=1.0;		//No row interchanges yet.
	for (i=1;i<=n;i++) { //Loop over rows to get the implicit scaling information
		big=0.0; 
		for (j=1;j<=n;j++)
		if ((temp=fabs(a[i-1][j-1])) > big) big=temp;
		if (fabs(big)<1E-10){
			printf("Singular matrix in function LU_dcmp()!!\n");
			for(j=1;j<=n;j++){
					printf("LU_dcmp: a[%d][%d]=%lf\n",i,j,a[i-1][j-1]);
					fflush(stdout);
			}
			exit(35);
		}
		//No nonzero largest element.
		vv[i-1]=1.0/big; //Save the scaling.
	}
	for (j=1;j<=n;j++) {	//This is the loop over columns of Crout's method.
		for (i=1;i<j;i++) {		//This is equation (2.3.12) except for i = j.
			sum=a[i-1][j-1];
			for (k=1;k<i;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
		}
		big=0.0;	//Initialize for the search for largest pivot element.
		for (i=j;i<=n;i++) {	//This is i = j of equation (2.3.12) and i = j+1...N
			sum=a[i-1][j-1];		//of equation (2.3.13).
			for (k=1;k<j;k++) sum -= a[i-1][k-1]*a[k-1][j-1];
			a[i-1][j-1]=sum;
			if ( (dum=vv[i-1]*fabs(sum)) >= big) {
				// Is the figure of merit for the pivot better than the best so far?
				big=dum;
				imax=i;
			}
		}
		//Do we need to interchange rows?
		if (j != imax) {	//Yes, do so...
			for (k=1;k<=n;k++) {	
				dum=a[imax-1][k-1];
				a[imax-1][k-1]=a[j-1][k-1];
				a[j-1][k-1]=dum;
			}
			*d = -(*d);		//...and change the parity of d.
			vv[imax-1]=vv[j-1];	//Also interchange the scale factor.
		}
		indx[j-1]=imax;
		if (fabs(a[j-1][j-1]) < 1E-5) a[j-1][j-1]=TINY;
		/* If the pivot element is zero the matrix is singular (at least to the precision of the
		algorithm). For some applications on singular matrices, it is desirable to substitute
		TINY for zero.*/
		if (j != n) {	//Now, finally, divide by the pivot element.
			dum=1.0/(a[j-1][j-1]);
			for (i=j+1;i<=n;i++) a[i-1][j-1] *= dum;
		}
	}	//Go back for the next column in the reduction.
	
	free(vv);

	return;
}/*end LU_dcmp()*/

void LU_bksb(real **a, int n, int *indx, real b[])
{
	/* Adopted from "Numerical Recipeis in C" by W Press et. al.

	   Solve the set of n linear equations A*X = B.
INPUT:
	a -- matrix A in the LU decomposition form obtained through LU_dcmp()
	n -- dimension
	indx -- output of LU_dcmp()
	b -- vextor B
OUTPUT:
	b -- solution X.							*/

	int i,ii=0,ip,j;
	real sum;

	/* When ii is set to a positive value, it will become the
		index of the first nonvanishing element of b. We now
		do the forward substitution, equation (2.3.6). The
		only new wrinkle is to unscramble the permutation
		as we go.*/
	for (i=1;i<=n;i++) {
		ip=indx[i-1];
		sum=b[ip-1];
		b[ip-1]=b[i-1];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i-1][j-1]*b[j-1];
		else if (fabs(sum)>1E-5) ii=i;		//A nonzero element was encountered, so from now on we
								//will have to b[i]=sum; do the sums in the loop above.
		b[i-1] = sum;
	}
	for (i=n;i>=1;i--) {	//Now we do the backsubstitution, equation (2.3.7).
		sum=b[i-1];
		for (j=i+1;j<=n;j++) sum -= a[i-1][j-1]*b[j-1];
		b[i-1]=sum/a[i-1][i-1];	//Store a component of the solution vector X.
	}	// All done!
}/*end LU_bksb()*/

void LU_inv_66(voigt66 c)
{
	/* Inverse the matrix using LU decomposition.
	   The original matrix will be replaced with
	   its inverse.
		*A special case for 6x6 matrix with voigt66 type,
		modified from LU_inverse(real**, int) */

	int indx[6];
	int n;
	real d;
	int i, j;
	voigt col;
	voigt66 y;
	real *a[6] = {c[0],c[1],c[2],c[3],c[4],c[5]};

	n = 6;
	LU_dcmp(a,n,indx,&d);
	for(j=0;j<n;j++){
		for(i=0;i<n;i++) col[i] = 0.0;
		col[j] = 1.0;
		LU_bksb(a,n,indx,col);
		for(i=0;i<n;i++) y[i][j] = col[i];
	}

	for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			c[i][j] = y[i][j];

	return;
}/*end LU_inv_66()*/

static inline void M3Inverse(real a[3][3], real ia[3][3])
{
	int i, j;
	real deta = 0.0;

	deta = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1]) \
		   -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0]) \
		   +a[0][2]*(a[1][0]*a[2][1]-a[2][0]*a[1][1]);
	ia[0][0] = a[1][1]*a[2][2] - a[1][2]*a[2][1];
	ia[0][1] = a[0][2]*a[2][1] - a[0][1]*a[2][2];
	ia[0][2] = a[0][1]*a[1][2] - a[0][2]*a[1][1];
	ia[1][0] = a[1][2]*a[2][0] - a[1][0]*a[2][2];
	ia[1][1] = a[0][0]*a[2][2] - a[0][2]*a[2][0];
	ia[1][2] = a[0][2]*a[1][0] - a[0][0]*a[1][2];
	ia[2][0] = a[1][0]*a[2][1] - a[1][1]*a[2][0];
	ia[2][1] = a[0][1]*a[2][0] - a[0][0]*a[2][1];
	ia[2][2] = a[0][0]*a[1][1] - a[0][1]*a[1][0];
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			ia[i][j] /=deta;

	return;
}/*end M3Inverse()*/

void CalcGAMMA(ten4th cijkl0, ten4th *gamma)
{

#pragma acc data present(gamma)
#pragma acc parallel loop collapse(3)
	k_loop{
		C4_loop{
			gamma[kIDX][mi][mj][mk][ml] = 0.0;
		}
	}

#pragma acc data present(gamma,g_v,cijkl0)
#pragma acc parallel loop collapse(3)
	k_loop{
	real nvector[3];
	real norm_nvector;
	real iomega[3][3], omega[3][3];
		nvector[0] = g_v[kxyz][0];
		nvector[1] = g_v[kxyz][1];
		nvector[2] = g_v[kxyz][2];
		norm_nvector=VSecNorm(nvector);
		VScale(nvector,1.0/norm_nvector);
		if(kIDX==0){		
			T2_loop{
				T2p_loop{
						gamma[kIDX][mi][mj][mip][mjp] = 0.0;
				}
			}
		}
		else{
			// Calculate Omega
			for(int i=0;i<3;i++){
				for(int k=0;k<3;k++){
					iomega[i][k] = 0.0;
					for(int j=0;j<3;j++)
						for(int l=0;l<3;l++)
							iomega[i][k] += cijkl0[i][j][k][l]*nvector[j]*nvector[l];
				}
			}

			M3Inverse(iomega, omega);

			T2_loop{
				T2p_loop{
					gamma[kIDX][mi][mj][mip][mjp] = omega[mi][mip]*nvector[mj]*nvector[mjp];
				}
			}
		}
	}

	return;
}/*end CalcGAMMA()*/

static inline double UpdateVirtualStrain(int opt)
{
    /*  
     *  opt = 
     *  0: update virtual strain only
     *  1: after 0, calculate stress field
     *  */

    /*  record the virtual SFTS of last time */
		char temp[200];

#pragma acc data present(iEps0_last,vEps)
#pragma acc parallel loop collapse(3)
    x_loop{
        V6_loop{
            iEps0_last[pIDX][mi] = vEps[pIDX][mi];
        }
    }
    // Stress of SFTS (current virtual strain) in k-space
	cufft_fr(vEps,kSig);

#pragma acc data present(Cij0,kSig)
#pragma acc parallel loop collapse(3)
    k_loop{
	    voigtk sym_du;
		V6_loop{
			sym_du[mi].x = 0.0;
			sym_du[mi].y = 0.0;
			V6p_loop{
				sym_du[mi].x += Cij0[mi][mip]*kSig[kIDX][mip].x;
				sym_du[mi].y += Cij0[mi][mip]*kSig[kIDX][mip].y;
					}
			}
		V6_loop{
			kSig[kIDX][mi].x = sym_du[mi].x;
			kSig[kIDX][mi].y = sym_du[mi].y;
			}
		}

	// Calculate the homogeneous elastic part and store to Tau_0
#pragma acc data present(kSig,GAMMA)
#pragma acc parallel loop collapse(3)
    k_loop{
    ten2ndk sig0;
    ten2ndk tmpeps;
        chg_basis_Kelvin_1k(kSig[kIDX],sig0);
        T2_loop{
		    tmpeps[mi][mj].x = 0.0;
		    tmpeps[mi][mj].y = 0.0;
		    T2p_loop{
                tmpeps[mi][mj].x += (GAMMA[kIDX][mi][mj][mip][mjp]+
                                  GAMMA[kIDX][mj][mi][mip][mjp])*sig0[mip][mjp].x;
                tmpeps[mi][mj].y += (GAMMA[kIDX][mi][mj][mip][mjp]+
                                    GAMMA[kIDX][mj][mi][mip][mjp])*sig0[mip][mjp].y;
            }
        }
        T2_loop{
            tmpeps[mi][mj].x *= 0.5;
            tmpeps[mi][mj].y *= 0.5;
        }
        chg_basis_Kelvin_2k(tmpeps,kSig[kIDX]);
    }//end of k loop

	cufft_fk(kSig,rSig);
#pragma acc data present(Tau_0,Sij0,Cij0,Cij2,hc_r,rSig)
#pragma acc parallel loop collapse(3)
    x_loop{
        V6_loop{
            Tau_0[pIDX][mi] = 0.0;
            V6p_loop{
                for(int ix=0; ix<6; ++ix){
                    Tau_0[pIDX][mi] += Sij0[mi][mip]*(hc_r[pxyz]*Cij2[mip][ix])*rSig[pIDX][ix];
                }
            }
        }
    }
    // Add the isolated strain terms to Tau_0
#pragma acc data present(te_p,eta_r,T_tr)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
        T_tr[pIDX][mi]=0;
	for(int i=0;i<V;i++)
	{
            T_tr[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // update transformation strain
	}
	}
}
#pragma acc data present(Tau_0,vEps,E_oo,T_tr,eta_r)
#pragma acc parallel loop collapse(3)
    x_loop{
        V6_loop{
           // Tau_0[pIDX][mi] -= vEps[pIDX][mi] - heta_r[pxyz]*E_oo[mi];
            Tau_0[pIDX][mi] -= vEps[pIDX][mi] - T_tr[pIDX][mi];
           // Tau_0[pIDX][mi] -= vEps[pIDX][mi] - eta_r[px][py][pz]*E_oo[mi];
        }
    }

    // The term related to the B.C., stored in Tau.
    if(ElastBC==0){     // strain-controlled
#pragma acc data present(te_p,eta_r,T_tr)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
        T_tr[pIDX][mi]=0;
	for(int i=0;i<V;i++)
	{
            T_tr[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // update transformation strain
	}
	}
}
#pragma acc data present(Tau,Sij0,Cij0,Cij2,hc_r,EpsAvg,T_tr)
#pragma acc parallel loop collapse(3)
	    x_loop{
	        V6_loop{
	            Tau[pIDX][mi] = 0.;
	            V6p_loop{
	                for(int ix=0; ix<6; ++ix){
	                    Tau[pIDX][mi] += Sij0[mi][mip]*(hc_r[pxyz]*Cij2[mip][ix])*
	                        (EpsAvg[ix]-T_tr[pIDX][ix]);
	                  //  Tau[pIDX][mi] += Sij0[mi][mip]*(hc_r[pxyz]*Cij2[mip][ix])*
	                  //      (EpsAvg[ix]-eta_r[px][py][pz]*E_oo[ix]);
	                }
	            }
	        }
	    }
#ifdef DEBUG_ON
	{int rID=32*64*64+32*64+32;
			CURR_LINE();
	printf("Tau=%e\n",Tau[rID][2]);
	}
#endif
    }

    else if(ElastBC==1){    // stress-controlled
        voigt veps_avg = {0.};  // average strain include virtual and macro parts
		V6_loop{
	        real global_strain = 0.0;
#pragma acc data present(vEps)
#pragma acc parallel loop collapse(3) reduction(+:global_strain)
	        x_loop{
	            global_strain += vEps[pIDX][mi];
	        }
			global_strain /= LMN;
			veps_avg[mi] = global_strain;
	    }
        V6_loop{
            V6p_loop{
                veps_avg[mi] += Sij0[mi][mip]*SigAvg[mip];
            }
        }
#pragma acc data present(te_p,eta_r,T_tr)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
        T_tr[pIDX][mi]=0;
	for(int i=0;i<V;i++)
	{
            T_tr[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // update transformation strain
	}
	}
}

#pragma acc data present(Tau,Sij0,Cij0,Cij2,hc_r,E_oo,eta_r,T_tr)
#pragma acc data copyin(veps_avg)
#pragma acc parallel loop collapse(3)
	    x_loop{
	        V6_loop{
	            Tau[pIDX][mi] = 0.;
	            V6p_loop{
	                for(int ix=0; ix<6; ++ix){
	                    Tau[pIDX][mi] += Sij0[mi][mip]*(hc_r[pxyz]*Cij2[mip][ix])*
	                        (veps_avg[ix]-T_tr[pIDX][ix]);
	               //     Tau[pIDX][mi] += Sij0[mi][mip]*(hc_r[pxyz]*Cij2[mip][ix])*
	               //         (veps_avg[ix]-eta_r[px][py][pz]*E_oo[ix]);
	                }
	            }
	        }
	    }
    }
    else{
        printf("Invalid elasticity boundary condition!\n");
    }

    // Update virtual strain
#pragma acc data present(vEps,Tau_0,Tau)
#pragma acc parallel loop collapse(3)
    x_loop{
        V6_loop{
            vEps[pIDX][mi] += (Tau_0[pIDX][mi] + Tau[pIDX][mi]);
        }
    }

	if(opt==1){

    // Stress of SFTS (current virtual strain) in k-space
	cufft_fr(vEps,kSig);
#pragma acc data present(kSig,Cij0)
#pragma acc parallel loop collapse(3)
	        k_loop{
	    voigtk sym_du;
		V6_loop{
			sym_du[mi].x = 0.0;;
			sym_du[mi].y = 0.0;;
			V6p_loop{
				sym_du[mi].x += Cij0[mi][mip]*kSig[kIDX][mip].x;
				sym_du[mi].y += Cij0[mi][mip]*kSig[kIDX][mip].y;
					}
			}
		V6_loop{
			kSig[kIDX][mi].x = sym_du[mi].x;
			kSig[kIDX][mi].y = sym_du[mi].y;
			}
		}

	    // Calculate the total strain and store in kSig_r temporarily
#pragma acc data present(kSig,GAMMA)
#pragma acc parallel loop collapse(3)
    k_loop{
    ten2ndk sig0;
    ten2ndk tmpeps;
        chg_basis_Kelvin_1k(kSig[kIDX],sig0);
        T2_loop{
		    tmpeps[mi][mj].x = 0.0;
		    tmpeps[mi][mj].y = 0.0;
		    T2p_loop{
                tmpeps[mi][mj].x += (GAMMA[kIDX][mi][mj][mip][mjp]+
                                  GAMMA[kIDX][mj][mi][mip][mjp])*sig0[mip][mjp].x;
                tmpeps[mi][mj].y += (GAMMA[kIDX][mi][mj][mip][mjp]+
                                    GAMMA[kIDX][mj][mi][mip][mjp])*sig0[mip][mjp].y;
            }
        }
        T2_loop{
            tmpeps[mi][mj].x *= 0.5;
            tmpeps[mi][mj].y *= 0.5;
        }
        chg_basis_Kelvin_2k(tmpeps,kSig[kIDX]);
    }//end of k loop

	cufft_fk(kSig,rSig);
 
	    // Update stress
	    if(ElastBC==0){
#pragma acc data present(Sig,Cij0,EpsAvg,rSig,vEps)
#pragma acc parallel loop collapse(3)
		    x_loop{
		        V6_loop{
		            Sig[pIDX][mi] = 0.0;
		            V6p_loop{
		                Sig[pIDX][mi] += Cij0[mi][mip]*(EpsAvg[mip]+rSig[pIDX][mip]-vEps[pIDX][mip]);
		            }
		        }
		    }
        }
        else if(ElastBC==1){
	        voigt veps_avg = {0.};  // average virtual strain
			V6_loop{
		        real global_strain = 0.;
#pragma acc data present(vEps)
#pragma acc parallel loop collapse(3) reduction(+:global_strain)
		        x_loop{
		            global_strain += vEps[pIDX][mi];
		        }
				global_strain /= LMN;
				veps_avg[mi] = global_strain;
		    }

#pragma acc data present(Sig,SigAvg,Cij0,rSig,vEps)
#pragma acc data copyin(veps_avg)
#pragma acc parallel loop collapse(3)
		    x_loop{
		        V6_loop{
		            Sig[pIDX][mi] = SigAvg[mi];
		            V6p_loop{
		                Sig[pIDX][mi] += Cij0[mi][mip]*(veps_avg[mip]+rSig[pIDX][mip]-vEps[pIDX][mip]);
		            }
		        }
		    }
        }
        else{
            fprintf(stderr,"Invalid elasticity boundary condition!\n");
        }

	}   // end of opt==1

    //  check convergence 
    real convg = 0.;
#pragma acc data present(vEps,iEps0_last)
#pragma acc parallel loop collapse(3) reduction(+:convg)
	x_loop{
		V6_loop{
			convg += (vEps[pIDX][mi]-iEps0_last[pIDX][mi])*(vEps[pIDX][mi]-iEps0_last[pIDX][mi]);
		}
	}
	convg /= LMN;
	convg = sqrt(convg);
    return convg;
}/*  end of UpdateVirtualStrain() */

void cufft_fr(voigt *fr, voigtk *fk)
{
#pragma acc data create(temp4_r,temp4_k)
V6_loop{
#pragma acc data present(fr)
#pragma acc parallel loop collapse(3)
	x_loop{temp4_r[pIDX]=fr[pIDX][mi];}
	cufftrc3k(temp4_r,temp4_k);
#pragma acc data present(fk)
#pragma acc parallel loop collapse(3)
	k_loop{fk[kIDX][mi].x=temp4_k[kIDX].x;
	       fk[kIDX][mi].y=temp4_k[kIDX].y;}
	}
}// end of cufft_fr()

void cufft_fk(voigtk *fk, voigt *fr)
{
#pragma acc data create(temp4_r,temp4_k)
V6_loop{
#pragma acc data present(fk)
#pragma acc parallel loop collapse(3)
	k_loop{temp4_k[kIDX].x=fk[kIDX][mi].x;
	       temp4_k[kIDX].y=fk[kIDX][mi].y;}
	cufftcr3k(temp4_k,temp4_r);
#pragma acc data present(fr)
#pragma acc parallel loop collapse(3)
	x_loop
		{fr[pIDX][mi]=temp4_r[pIDX];}
	}
}//end of cufft_fk()

double cal_Felas()
{
double sum=0;
#pragma acc data present(te_p,eta_r,T_tr)
#pragma acc parallel loop collapse(3) 
	x_loop{ 
	V6_loop{
        T_tr[pIDX][mi]=0;
	for(int i=0;i<V;i++)
	{
            T_tr[pIDX][mi] += eta_r[i][px][py][pz]*te_p[i][mi];       // update transformation strain
	}
	}
}
#pragma acc data present(vEps,T_tr,E_oo,Sij0,Sig)
#pragma acc parallel loop collapse(3) reduction(+:sum)
x_loop{
//	V6_loop{sum+=Sig[pIDX][mi]*(vEps[pIDX][mi]-heta_r[pxyz]*E_oo[mi]);
	V6_loop{sum+=Sig[pIDX][mi]*(vEps[pIDX][mi]-T_tr[pIDX][mi]);
		V6p_loop{sum+=Sig[pIDX][mi]*Sij0[mi][mip]*Sig[pIDX][mip];}
	       }
      }
sum*=LMN*0.5;
return sum;
}//end of calFelas()

void Set_rot(double rott[3][3])
{
int ii,jj,kk,ll;
int ip,jp,kp,lp;
real te_p_temp[6]={0.0};
real te_p_temp33[3][3]={{0.0}};
real te_p_rot[3][3]={{0.0}};
//real te_p[3][3];
ten4th Cijkl0_temp = {{{{ 0.0 }}}};

field_loop{
V6_loop{
te_p_temp[mi]=te_p[ip][mi];
}
chg_basis_Kelvin_1(te_p_temp,te_p_temp33);
T2_loop{
te_p_rot[mi][mj]=0.0;
   T2p_loop{
   te_p_rot[mi][mj]+=te_p_temp33[mip][mjp]*rott[mi][mip]*rott[mj][mjp];
   }
}
chg_basis_Kelvin_2(te_p_rot,te_p_temp);
V6_loop{
te_p[ip][mi]=te_p_temp[mi];
printf("tep_V%d_%d=%lf\n",ip+1,mi+1,te_p[ip][mi]);
}
}

C4_loop{
Cijkl0_temp[mi][mj][mk][ml]=0.0;
		for(ii=0;ii<3;ii++)
		for(jj=0;jj<3;jj++)
		for(kk=0;kk<3;kk++)
		for(ll=0;ll<3;ll++)
		{
			Cijkl0_temp[mi][mj][mk][ml]+=rott[mi][ii]*rott[mj][jj]*rott[mk][kk]*rott[ml][ll]*Cijkl0[ii][jj][kk][ll];
		}

}
C4_loop{
Cijkl0[mi][mj][mk][ml]=Cijkl0_temp[mi][mj][mk][ml];
printf("C0_%d_%d_%d_%d=%lf\n",mi+1,mj+1,mk+1,ml+1,Cijkl0[mi][mj][mk][ml]);
}
C4_loop{
Cijkl0_temp[mi][mj][mk][ml]=0.0;
		for(ii=0;ii<3;ii++)
		for(jj=0;jj<3;jj++)
		for(kk=0;kk<3;kk++)
		for(ll=0;ll<3;ll++)
		{
			Cijkl0_temp[mi][mj][mk][ml]+=rott[mi][ii]*rott[mj][jj]*rott[mk][kk]*rott[ml][ll]*Cijkl1[ii][jj][kk][ll];
		}

}
C4_loop{
Cijkl1[mi][mj][mk][ml]=Cijkl0_temp[mi][mj][mk][ml];
printf("C1_%d_%d_%d_%d=%lf\n",mi+1,mj+1,mk+1,ml+1,Cijkl1[mi][mj][mk][ml]);
}
}
