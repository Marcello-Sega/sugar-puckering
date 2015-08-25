#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <getopt.h>



#define real double

#define arccos acos
#define Sqrt sqrt
#define X 0
#define Y 1 
#define Z 2

#define TA 109.471220634491
#define TD 0.154
double ccc=TA ,coc=TA ,occ=TA, cc=TD ,co=TD;

real
SQNorm (real * v)
{
  return (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

real
Norm (real * v)
{
  return Sqrt (SQNorm (v));
}


void rotmatr(double phi,double *R){
	int i,j,dim=3;
	for(i=0;i<dim;i++) for(j=0;j<dim;j++) R[j+dim*i]=0.;;
	R[0+0*dim]=cos(phi); R[1+0*dim]=-sin(phi); R[2+0*dim]=0.;
	R[0+1*dim]=sin(phi); R[1+1*dim]=cos(phi); R[2+1*dim]=0.;
	R[0+2*dim]=0.; R[1+2*dim]=0.; R[2+2*dim]=1.;
}

real
IProd (real * a, real * b)
{
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

real *
VProd (real * a, real * b, real * w)
{
  w[X] = a[Y] * b[Z] - a[Z] * b[Y];
  w[Y] = -a[X] * b[Z] + a[Z] * b[X];
  w[Z] = a[X] * b[Y] - a[Y] * b[X];
  return w;
}

void
DirectionSelf (real * v)
{
  real d;
  d = 1. / Norm (v);
  v[X] *= d;
  v[Y] *= d;
  v[Z] *= d;
}

enum dihedraltypes {DIHEDRAL_LINEAR, DIHEDRAL_BRANCHED};
#define SMALL 1e-5
real
Dihedral (real * a, real * b, real * c, real * d)
{
  real r[3], s[3], t[3], v[3], w[3], sign=1.0, scalarprod=0;
  r[X] = b[X] - a[X];
  r[Y] = b[Y] - a[Y];
  r[Z] = b[Z] - a[Z];
  s[X] = c[X] - b[X];
  s[Y] = c[Y] - b[Y];
  s[Z] = c[Z] - b[Z];
  t[X] = d[X] - c[X];
  t[Y] = d[Y] - c[Y];
  t[Z] = d[Z] - c[Z];
  VProd (s, r, v);
  sign = ((IProd (v, t) >= 0.0) ? -1.0 : 1.0); 
  VProd (t, s, w);
  DirectionSelf (v);
  DirectionSelf (w);
  scalarprod = IProd (v, w);
  if (scalarprod < -(1.-SMALL)){
    return M_PI;
  }else if (scalarprod > (1.-SMALL)) {
          return 0.0;
  }else
  /* sign is positive if (s x r) points in the same direction as t */
  return (sign * arccos (scalarprod));	
}



int pickett(double Q,double theta, double phi, double * coords, double *alphas,int rad_out_flag, int angstrom_out_flag)
{
	int N=6,dim=3;
	int i,j,k;
	int counter=0;
	double q2,phi2,q3;
	double d1[dim],d2[dim],d3[dim],d4[dim],d5[dim],d6[dim];
	double z[6];
	double rij[6],betaijk[6];
	double rpij[6],betapijk[6],COSbetapijk[6],SINbetapijk[6];
	double OP,QP,OQ;
	double S11[dim],S12[dim],S13[dim];
	double S23[dim],S24[dim],S25[dim];
	double S35[dim],S36[dim],S31[dim];
	double rho1,rho2,rho3;
	double pO[dim],pP[dim],pQ[dim];
	double rhoPS1,rhoPS2,sigma1,sigma2,sigma3;
	double p1[dim],p2[dim],p3[dim],p4[dim],p5[dim],p6[dim];
	double R1[3*dim];
	double R2[3*dim];
	double R3[3*dim];
	double G[dim],rhoG;
	double R4[3*dim];

	//for(i=0;i<N;i++)  rij[i]=b;
	rij[0]=co*10.;
	rij[5]=co*10.;
	rij[1]=rij[2]=rij[3]=rij[4]=cc*10.;
	betaijk[0]=betaijk[4]=occ*M_PI/180.;
	betaijk[1]=betaijk[2]=betaijk[3]=ccc*M_PI/180.;
	betaijk[5]=coc*M_PI/180.;
	for(j=0;j<N;j++)
		z[j]= Q * sin(theta) * cos( phi + (2.*M_PI*j)/3. )/sqrt(3.) + Q * cos(theta) * pow(-1,j)/sqrt(6);

	d1[2]=z[0]; d2[2]=z[1]; d3[2]=z[2]; d4[2]=z[3]; d5[2]=z[4]; d6[2]=z[5];

	//PROJECTED Bond lenght
	for(i=0;i<N;i++)
		rpij[i]=sqrt( rij[i]*rij[i] - (z[i]-z[(i+1)%N])*(z[i]-z[(i+1)%N]) );
	//PROJECTED Bond angles 
	for(i=0;i<N;i++){
		COSbetapijk[i]= (  pow((z[(i+2)%N]-z[i]),2)-pow((z[(i+1)%N]-z[i]),2)-pow((z[(i+2)%N]-z[(i+1)%N]),2) + 2.*rij[i]*rij[(i+1)%N]*cos(betaijk[i]) )/( 2.*rpij[i]*rpij[(i+1)%N] ) ;
		SINbetapijk[i]=sqrt(1-COSbetapijk[i]*COSbetapijk[i]);
	}	// betaijk[0]=beta123; betaijk[1]=beta234; ... betaijk[5]=beta612;
	OP=sqrt(rpij[0]*rpij[0] + rpij[1]*rpij[1] - 2.*rpij[0]*rpij[1]*COSbetapijk[0]);
	QP=sqrt(rpij[2]*rpij[2] + rpij[3]*rpij[3] - 2.*rpij[2]*rpij[3]*COSbetapijk[2]);
	OQ=sqrt(rpij[4]*rpij[4] + rpij[5]*rpij[5] - 2.*rpij[4]*rpij[5]*COSbetapijk[4]);
	S11[0] = 0.;		S11[1] = 0.;
	S12[0] = -rpij[0];	S12[1] = 0.;
	S13[0] = -rpij[0]+rpij[1]*COSbetapijk[0];	S13[1] = rpij[1]*SINbetapijk[0];
	S23[0] = OQ+rpij[3]-rpij[2]*COSbetapijk[2];	S23[1] = rpij[2]*SINbetapijk[2];
	S24[0] = OQ+rpij[3];	S24[1] = 0.;
	S25[0] = OQ;		S25[1] = 0.;
	S35[0] = rpij[5]-rpij[4]*COSbetapijk[4];	S35[1] = rpij[4]*SINbetapijk[4];
	S36[0] = rpij[5];	S36[1] = 0.;
	S31[0] = 0.;	S31[1] = 0.;
	S11[2]=S12[2]=S13[2]=S23[2]=S24[2]=S25[2]=S35[2]=S36[2]=S31[2]=0.;
	rho1 = atan2(S13[1],S13[0]);
	rho2 = atan2(S23[1],S23[0]-OQ);
	rho3 = atan2(S35[1],S35[0]);
	//printf("Rotation angles (first): %lf %lf %lf\n",rho1,rho2,rho3);
	pO[0] = 0.;	pO[1] = 0.;
	pP[0] = (OP*OP+OQ*OQ-QP*QP)/(2.*OQ);
	pP[1] = sqrt( OP*OP - ( (OP*OP+OQ*OQ-QP*QP)*(OP*OP+OQ*OQ-QP*QP) )/(4.*OP*OP) );
	pQ[0] = OQ;	pQ[1] = 0.;
	pO[2]=pP[2]=pQ[2]=0.;
	rhoPS1 = atan2(pP[1],pP[0]);
	rhoPS2 = atan2(pP[1],pP[0]-OQ);
	sigma1 = rho1 - rhoPS1 ;
	sigma2 = rhoPS2 - rho2;
	sigma3 = rho3;
	//printf("Rotation angles (second): %lf %lf %lf %lf %lf\n",rhoPS1,rhoPS2,sigma1,sigma2,sigma3);
	for(i=0;i<dim;i++){
		p1[i] = pO[i];
		p3[i] = pP[i];
		p5[i] = pQ[i];
	}

	rotmatr(-sigma1,R1);
	for(i=0;i<dim;i++){
		p2[i]=0.;
		for(j=0;j<dim;j++){
			p2[i]+=R1[j+i*dim]*S12[j];
		}
	}


	rotmatr(sigma2,R2);
	for(i=0;i<dim;i++){
		p4[i]=0.;
		for(j=0;j<dim;j++){
			p4[i]+=R2[j+i*dim]*(S24[j]-pQ[j]);
		}
		p4[i]+=pQ[i];
	}

	rotmatr(-sigma3,R3);
	for(i=0;i<dim;i++){
		p6[i]=0.;
		for(j=0;j<dim;j++){
			p6[i]+=R3[j+i*dim]*S36[j];
		}
	}

	G[0]=G[1]=G[2]=0.;
	for(j=0;j<dim;j++){
		G[j]=p1[j]+p2[j]+p3[j]+p4[j]+p5[j]+p6[j];
		G[j]/=6.;
	}
	rhoG=M_PI/2.+atan2(G[1],G[0]);

	rotmatr(-rhoG,R4);
	for(i=0;i<2;i++){
		d1[i]=0.;
		for(j=0;j<dim;j++){
			d1[i]+=R4[j+i*dim]*(p1[j]-G[j]);
		}
		d2[i]=0.;
		for(j=0;j<dim;j++){
			d2[i]+=R4[j+i*dim]*(p2[j]-G[j]);
		}
		d3[i]=0.;
		for(j=0;j<dim;j++){
			d3[i]+=R4[j+i*dim]*(p3[j]-G[j]);
		}
		d4[i]=0.;
		for(j=0;j<dim;j++){
			d4[i]+=R4[j+i*dim]*(p4[j]-G[j]);
		}
		d5[i]=0.;
		for(j=0;j<dim;j++){
			d5[i]+=R4[j+i*dim]*(p5[j]-G[j]);
		}
		d6[i]=0.;
		for(j=0;j<dim;j++){
			d6[i]+=R4[j+i*dim]*(p6[j]-G[j]);
		}
	}
	for(j=0;j<3;j++) { coords[counter++]=d1[j]; } 
	for(j=0;j<3;j++) { coords[counter++]=d2[j]; } 
	for(j=0;j<3;j++) { coords[counter++]=d3[j]; } 
	for(j=0;j<3;j++) { coords[counter++]=d4[j]; } 
	for(j=0;j<3;j++) { coords[counter++]=d5[j]; } 
	for(j=0;j<3;j++) { coords[counter++]=d6[j]; } 
	if(angstrom_out_flag==0) for(j=0;j<18;j++) coords[j]/=10.0; 
        alphas[0] = 180./M_PI * Dihedral (d5,d1,d3,d2) - 180. ; 
        alphas[1] = 180./M_PI * Dihedral (d1,d3,d5,d4) - 180. ;
        alphas[2] = 180./M_PI * Dihedral (d3,d5,d1,d6) - 180. ;
	for(j=0;j<3;j++) { 
               while (alphas[j]>180.0) alphas[j]-=360.0; 
               while (alphas[j]<-180.0) alphas[j]+=360.0; 
	}
	if(rad_out_flag==1) for(j=0;j<3;j++) alphas[j]*=M_PI/180.; 
        	
	return 0;
}

void print_coords(double * coords,int verbose_flag){
	int i,j,counter=0;
	char name[6][3]={"O5\0","C1\0","C2\0","C3\0","C4\0","C5\0"};
	for(i=0;i<6;i++){
           if(verbose_flag) printf ("%s : ",name[i]);
	   for(j=0;j<3;j++){
		printf("%f ",coords[counter++]);
           }
	   printf("\n");
        }
}
void print_pickett(double * alphas,int verbose_flag){
	int i;
	for(i=0;i<3;i++){
           if(verbose_flag) printf ("alpha_%d : ",i);
	   printf("%f ",alphas[i]);
        }
	printf("\n");
}

void usage(double coc,double ccc,double occ, double cc, double co,char ** argv) { 
	printf("\nReconstructs atomic coordinates (in nm) in 6-membered rings, from the Q,theta,phi Cremer-Pople coordinates\n");
	printf("-h | --help\t\t This help\n");
	printf("-v | --verbose\t\t Print additional information\n");
//	printf("-t | --tetrahydropyran\t\t Changes default angles and bond lengths from those of cyclohexane to those of tetrahydropyran\n");
	printf("-r | --rad_in\t\t Input Cremer-Pople angles in rad\n");
	printf("-R | --rad_out\t\t Output Strauss-Pickett angles in rad\n");
	printf("-a | --angstrom_in\t Input Q in angstrom\n");
	printf("-A | --angstrom_out\t Output atomic positions in angstrom \n");
	printf("-p | --pickett\t\t Dump the values of Strauss-Pickett out-of-plane dihedral angles (in deg) \n");
	printf("-n | --no_coords\t Do not dump atomic coordinates\n");
	printf("-Q | --coc\t\t Set the value of the c-o-c angle (in deg). Default: %f\n",coc);
	printf("-W | --ccc\t\t Set the value of the c-c-c angle (in deg). Default: %f\n",ccc);
	printf("-E | --occ\t\t Set the value of the o-c-c angle (in deg). Default: %f\n",occ);
	printf("-q | --cc\t\t Set the value of the c-c bond (in nm). Default: %f\n",cc);
	printf("-w | --co\t\t Set the value of the c-o bond (in nm). Default: %f\n",co);
	printf("\nexample (generate an inverted chair): echo \"0.06 180.0 0.0\" |  %s -v -p\n\n",argv[0]);
        exit(0);
}
int main (int argc, char ** argv){ 

    double alphas[3],coords[18];
    int pickett_flag=0,verbose_flag=0,no_coords_flag=0,angstrom_in_flag=0,angstrom_out_flag=0,rad_in_flag=0,rad_out_flag=0,THP_flag=0;
    double Q,theta,phi;
    char str[1024],ch,*rets;

    extern char *optarg;
    extern int optind;
    extern int optopt;
    extern int opterr;
    extern int optreset;

    static struct option longopts[] = {
             { "help",       no_argument,            NULL,           'h' },
             { "pickett",    no_argument,            NULL,           'p' },
             { "tetrahydropyran",    no_argument,    NULL,           't' },
             { "verbose",    no_argument,            NULL,           'v' },
             { "no_coords",  no_argument,            NULL,           'n' },
             { "rad_in",     no_argument,            NULL,           'r' },
             { "rad_out",    no_argument,            NULL,           'R' },
             { "angstrom_in",no_argument,            NULL,           'a' },
             { "angstrom_out",no_argument,           NULL,           'A' },
             { "coc",        required_argument,      NULL,           'Q' },
             { "ccc",        required_argument,      NULL,           'W' },
             { "occ",        required_argument,      NULL,           'E' },
             { "cc",         required_argument,      NULL,           'q' },
             { "co",         required_argument,      NULL,           'w' },
             { NULL,         0,                      NULL,           0 }
     };

     while ((ch = getopt_long(argc, argv, "pntrRaAvhQ:W:E:q:w:", longopts, NULL)) != -1)
             switch (ch) {
             case 'p':   pickett_flag = 1;   break;
             case 'n':   no_coords_flag = 1; break;
	     case 't':   THP_flag = 1;    break;
	     case 'r':   rad_in_flag = 1;    break;
	     case 'R':   rad_out_flag = 1;   break;
	     case 'a':   angstrom_in_flag = 1;    break;
	     case 'A':   angstrom_out_flag = 1;   break;
	     case 'Q':   coc=atof(optarg);   break;
             case 'W':   ccc=atof(optarg);   break;
             case 'E':   occ=atof(optarg);   break;
             case 'q':   cc=atof(optarg);    break;
             case 'w':   co=atof(optarg);    break;
             case 'v':   verbose_flag=1;    break;
	     case 'h':  default: usage(coc,ccc,occ,cc,co,argv);
     }
     argc -= optind;
     argv += optind;
	 
    while(1){ 
         rets=fgets(str,1024,stdin);
	 if (rets==NULL) exit(0);
         sscanf(str,"%lf %lf %lf",&Q,&theta,&phi);
         if(Q<0) { printf("Warning: negative Q\n") ; } else { 
  	   if(theta > 180.  && rad_in_flag==0 || theta > M_PI && rad_in_flag==1){ printf("Warning: theta > 180 deg\n"); } else { 
	   if (theta < 0) {  printf("Warning: negative theta\n") ; } else { 
                if(angstrom_in_flag==0) { Q*=10;  } 
                if(rad_in_flag==0) { theta*=M_PI/180; phi*=M_PI/180;  } 
                pickett(Q,theta,phi,coords,alphas,rad_out_flag,angstrom_out_flag);
	        if(!no_coords_flag) print_coords(coords,verbose_flag);
	        if(pickett_flag)    print_pickett(alphas,verbose_flag);
         }
      }
     }
   } 

}
