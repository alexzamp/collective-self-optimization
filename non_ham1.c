/* LANGEVIN DYNAMICS FOR PARTICLES COMMUNICATING VIA OPTIMAL-FIELD COMMUNICATION RULES.
 Corresponding article: "Collective self-optimization of communicating active particles" by A. V. Zampetaki et al.
 Returns a folder of ".dat" files with names representative of the time step.
 Each ".dat" file has 5 columns  representing: time step, x-coordinate, y-coordinate, total force projection in x, total force projection in y
 in a different row for each particle.
 The simulation is performed in simulation box fo size Lx,Ly and is subjected to periodic 
boundary conditions (PB). More information is provided in the  "Simulations" section.
Author: Alexandra Zampetaki
 */


#define RAND48_SEED_0   (0x330e)
#define RAND48_SEED_1 (0xabcd)
#define RAND48_SEED_2 (0x1234)
#define RAND48_MULT_0 (0xe66d)
#define RAND48_MULT_1 (0xdeec)
#define RAND48_MULT_2 (0x0005)
#define RAND48_ADD (0x000b)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include<sys/stat.h> 
#include<sys/types.h>
#include <time.h>

/*useful constants*/
#define Pi 3.141592653589793
#define TwoPi 6.283185307179586
#define EMPTY -1
#define Lx 10// box dimension in x-axis
#define Ly 10//box dimension in y-axis
#define nparticles  14//number of particles
#define number_step 500001//number of time steps
#define write_after 4000 //write every this ammount of number of steps
#define dt 0.0001 


/*particles */
double particle_x[nparticles]; // x position of particles
double particle_y[nparticles]; // y position of particles 
double fptyol[nparticles]; // force projection in y axis
double fptxol[nparticles];// force projection in x axis

/* Parameters */
double kappa=1; // Yukawa screening parameter (see Eq. [2])
double Tcom=6; //optimal  field value
double lamb=1; //force strength value (see Eq. [3])
double Dt=sqrt(3.0)*0.001*1;//translational diffusion coefficient



// Initialization
int step=0;
int seir=53; // to reproduce different noise samples change this value
FILE *output;

/* Routines for random number generation*/
unsigned short _rand48_seed[3] = {
        RAND48_SEED_0,
         RAND48_SEED_1,
         RAND48_SEED_2
};
unsigned short _rand48_mult[3] = {
         RAND48_MULT_0,
         RAND48_MULT_1,
         RAND48_MULT_2
 };
unsigned short _rand48_add = RAND48_ADD;

void
 _dorand48(unsigned short xseed[3])
 {
	         unsigned long accu;
	         unsigned short temp[2];
	
	         accu = (unsigned long)_rand48_mult[0] * (unsigned long)xseed[0] +
	          (unsigned long)_rand48_add;
	         temp[0] = (unsigned short)accu;        /* lower 16 bits */
	         accu >>= sizeof(unsigned short)* 8;
	         accu += (unsigned long)_rand48_mult[0] * (unsigned long)xseed[1] +
	          (unsigned long)_rand48_mult[1] * (unsigned long)xseed[0];
	         temp[1] = (unsigned short)accu;        /* middle 16 bits */
	         accu >>= sizeof(unsigned short)* 8;
	         accu += _rand48_mult[0] * xseed[2] + _rand48_mult[1] * xseed[1] + _rand48_mult[2] * xseed[0];
	         xseed[0] = temp[0];
	         xseed[1] = temp[1];
	         xseed[2] = (unsigned short)accu;
}

double erand48(unsigned short xseed[3])
{
         _dorand48(xseed);
         return ldexp((double) xseed[0], -48) +
                ldexp((double) xseed[1], -32) +
                ldexp((double) xseed[2], -16);
}

double drand48(){
	return erand48(_rand48_seed);
}

void srand48(long seed){
	_rand48_seed[0] = RAND48_SEED_0;
	_rand48_seed[1] = (unsigned short)seed;
	_rand48_seed[2] = (unsigned short)(seed >> 16);
	_rand48_mult[0] = RAND48_MULT_0;
	_rand48_mult[1] = RAND48_MULT_1;
	_rand48_mult[2] = RAND48_MULT_2;
	_rand48_add = RAND48_ADD;
}

void streamfile(int step);

/* Main Program */
int main(int argc, char *argv[])
{
	
	//useful initializations
  int n,i1,j1,nsqr, rat1, mod1;
  double sid1;
  double fp1[nparticles];
  double fp2x[nparticles];
  double fptx[nparticles];
  double fp2y[nparticles];
  double fpty[nparticles];
  
 
 
  /*Initialize particles*/
    int start = 0; 
    srand48(seir); //noise sample
    
    nsqr=ceil(sqrt(nparticles));
    printf( "ns:%d  \n",nsqr);
    
    
    //initialize particle positions in a square grid of side length "sid1"
    for(n=0;n<nparticles;n++){
        sid1=1.0; // side of the intialization square box 
        rat1=n/nsqr;
        mod1=n%nsqr;
        particle_x[n] =3.5+((double)(sid1))*drand48();
        particle_y[n] = 3.5+((double)(sid1))*drand48();
        fptx[n]=0; //force initialization
        fpty[n]=0;   //force initialization
    }
    streamfile(1);
    
    
  /*Main for loop for dynamics*/
  for(step = start; step < number_step; step++){ 
        
        for (i1=0; i1<nparticles; i1++){
		fp1[i1]=0.0; //initialize fp1
		fp2x[i1]=0.0; //initialize fp2x
		fp2y[i1]=0.0; //initialize fp2y
	}
        
        /*Force calculation*/
         for(i1=0; i1< nparticles; i1++){ 
			 for(j1=0; j1< nparticles; j1++){ 
				 if(j1!=i1){
				 double dist=sqrt(pow(particle_x[i1]-particle_x[j1],2.0)+pow(particle_y[i1]-particle_y[j1],2.0)); // interparticle distance rij
				 double xij=particle_x[i1]-particle_x[j1]; //xij
				 double yij=particle_y[i1]-particle_y[j1]; //yij
	             if(dist>0.0001){
	             fp1[i1] +=exp(-kappa*dist)/(dist);
	             fp2x[i1] -=exp(-kappa*dist)*(kappa*dist+1)*(xij)/pow(dist,3);
	             fp2y[i1] -=exp(-kappa*dist)*(kappa*dist+1)*(yij)/pow(dist,3);
	             }
				}
			 }
			 
			 fptx[i1]=-lamb*(fp1[i1]-Tcom)*fp2x[i1]; // Force according to Eq. [3]
			 fpty[i1]=-lamb*(fp1[i1]-Tcom)*fp2y[i1];
			 fptyol[i1]=fpty[i1]; //total force acting on i1 particle in x
			 fptxol[i1]=fptx[i1]; //total force acting on i1 particle in x
			// printf( "fa1:%16.8f  \n",fp1[1]);
	     }		      
	     
      
  
  
    /*Propagate particles (Langevin Dynamics)*/
    for(n = 0; n < nparticles; n++){
      
     
      particle_x[n] += dt*fptx[n]+ sqrt(2.0*Dt*dt)*sqrt(3.0)*(1.0-2.0*drand48());
      particle_y[n] += dt*fpty[n]+ sqrt(2.0*Dt*dt)*sqrt(3.0)*(1.0-2.0*drand48());
      
      
      /*Periodic boundary conditions*/
      if(particle_x[n]>Lx) particle_x[n] -= Lx;
      if(particle_x[n]<0) particle_x[n] += Lx;
      if(particle_y[n]>Ly) particle_y[n] -= Ly;
      if(particle_y[n]<0) particle_y[n] += Ly;

    }


    //output
    if(step%write_after == 0){
      streamfile(step);
      printf("%d\n",step);
    }

  }

}

/*Write  data*/
void streamfile(int step)
{
  int n;

  char p_dat[81];
  char folder[81];	
  char file[81];
  sprintf(p_dat, "A%d_T%.3f_Lx%d_Np%d_Dt%.5f/",seir,Tcom,Lx,nparticles,Dt);  //Name of folder
  strcpy(folder, p_dat);
  sprintf(file, "%d.dat", step); //Name of file
  strcat(p_dat, file);  
  if((output = fopen(p_dat, "w")) == NULL){
      mkdir(folder);
      output = fopen(p_dat, "w");
  }  
   if(step!=0){  
    for(n=0; n<nparticles; n++) {
    fprintf(output, "%d %16.8f %16.8f %16.8f %16.8f\n",step,particle_x[n],particle_y[n],fptxol[n],fptyol[n]);
    }
    fclose(output);
  }
}






