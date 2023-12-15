// QUESTO CODICE E' UGUALE ALL'ALTRO, MA USA UNA FUNZIONE PER IL METROPOLIS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define N 80
#define sigma 2  //nm
#define A 0.04 //eV
#define beta 1.0 // è il beta adimensionale, cioè A/(kT)
#define rho_adimensionale 0.75 // è la densità adimensionale, cioè rho*sigma^3
#define KB 8.617333262145*0.00001 // eV/K
#define AVO 6.02214076*pow(10.,23.) //mol^-1

#define therm_check 0
#define compute_g 0

double dist(double ri[3], double rj[3], double);
double interaction(double);
double V_tr_sh(double, double);
double V_trunc(double, double);
double virial_trunc(double, double);
double energy_tail(double);
double pressure_tail(double);
double shell_volume(double, double);
int Metropolis(double, double);


int main(){
    int Nit=500, M;
    double therm_time;
    double L, rc, Delta;
    double position[N][3], rnew[3], d_old, d_new, ddd; // le 3 sono x,y,z
    double Uold, Unew, deltaU, rate;
    double energy[Nit], av_energy=0, pressure_id, pressure_exc[Nit], av_pressure_exc=0, av_pressure;
    double epsilon, h[300];
    double sigma_metri, sigma_centimetri, A_Joule;
    int i, j, k, t, n;
    FILE *fp, *fp1;
    fp = fopen("file.dat", "w");
    fp1 = fopen("g2.dat", "w");

    L = pow( (double)N/rho_adimensionale, 1./3.); // è adimensionale
    // in units of sigma, if N = 80, L = 4.25 sigma = 8.5 nm
    // in units of sigma, if N = 200, L = 4.25 sigma = 8.5 nm
    rc = L/2. ; // in units of sigma
    Delta = L/10.; // in units of sigma

    if(compute_g!=0){
        epsilon = 1./10. ; //in units of sigma
        M = (int) (L/epsilon);
        //printf("\n M = %d \n\n", M);
        for(n=0; n<M; n++){
            h[n] = 0.;
        }
    }


    if(therm_check != 0){
        therm_time = 0;
    }
    else if(therm_check==0){
        if(Delta==L/10.){
            therm_time = 20; // se delta è L/10
        }
        else if(Delta==L/10.){
            therm_time = 10;
        }
        else if(Delta==L/10.){
            therm_time = 8;
        }
        else if(Delta==L/10.){
            therm_time = 6;
        }
        else {
            therm_time = 4;
        }
    }

    pressure_id = rho_adimensionale/beta ; // è rho*sigma^3*KT/A
    pressure_id = pressure_id * beta ; // è la pressione in unità di kT/sigma^3

    srand48(time(0)); //uso drand48() che genera in [0, 1)

    // STARTING CONFIGURATION
    for(i=0; i<N; i++){
        for(k=0; k<3; k++){
            position[i][k] = L*drand48(); // in units of sigma
        }
    }

    rate=0;
    // ITERATION
    for(t=0; t<Nit; t++){

        //SEQUENTIAL UPDATE
        for(i=0; i<N; i++){

            // PROPOSE THE MOVE
            for(k=0; k<3; k++){
                rnew[k] = position[i][k] + Delta * (drand48() - 0.5);
            }

            // IS THE MOVE ACCEPTED?
            Uold=0;
            Unew=0;
            for(j=0; j<N; j++){
                if(j != i){
                    d_old = dist(position[i], position[j], L);
                    d_new = dist(rnew, position[j], L);

                    Uold += V_tr_sh(d_old, rc);
                    Unew += V_tr_sh(d_new, rc);
                    deltaU = Unew-Uold;
                }
            }

            // IF IT IS ACCEPTED, UPDATE
            if(Metropolis(deltaU, beta)==1){
              for(k=0; k<3; k++){
                  position[i][k] = rnew[k];
              }
              if(t>=therm_time){
                  rate++;
              }
            }
            // ELSE DO NOTHING

        }

        //Once the system has been updated, MEASURE in equilibrium
        if(t>=therm_time){ // così scarto le iterazioni 0....(therm_time-1) che sono in tutto therm_time

            //energy per particle and pressure
            energy[t] = 0. ;
            pressure_exc[t] = 0. ;
            for(i=1; i<N; i++){
                for(j=0; j<i; j++){ // così la distanza tra i e j la conto una volta sola
                    ddd = dist(position[i], position[j], L) ;
                    energy[t] += V_trunc(ddd, rc)/N ; //energy per particle
                    pressure_exc[t] += -1./(3.*pow(L, 3)) * virial_trunc(ddd, rc);
                    if(compute_g!=0){
                        n = (int)(ddd/epsilon);
                        if(n>M){
                            printf("Errore");
                        }
                        h[n]++;
                    }
                }
            }
            fprintf(fp, "%d %f %f \n", t-(int)(therm_time-1), energy[t], pressure_exc[t]);
            av_energy += energy[t] ;
            av_pressure_exc += pressure_exc[t] ;


        }

    }

    fclose(fp);

    // ricordo chi era Delta
    printf("Nell'update Delta è L/%i\n", (int) (L/Delta)) ;

    Nit = Nit - therm_time;

    // what acceptance rate?
    rate = rate/(N*Nit);
    printf("Acceptance rate: %f \n", rate) ;

    // what energy per particle?
    av_energy = av_energy/(Nit) ; // in units of A
    av_energy = av_energy * beta ; // in units of KT
    printf("\nENERGY\n");
    printf("Average energy per particle without tail: %f in units of kT\n", av_energy) ;
    printf("Energy per particle, tail: %f in units of kT\n", energy_tail(rc));
    printf("Average energy per particle with tail: %f in units of kT\n", av_energy + energy_tail(rc)) ;

    // what pressure?
    printf("\nPRESSURE\n");
    av_pressure_exc = av_pressure_exc/Nit * beta ;
    printf("Average excess pressure without tail: %f in units of kT/s^3\n", av_pressure_exc);
    printf("Excess pressure, tail: %f in units of kT/s^3\n", pressure_tail(rc));
    printf("Average excess pressure with tail: %f in units of kT/s^3\n", av_pressure_exc + pressure_tail(rc));
    // printf("Ideal Pressure: %f in units of kT/s^3\n", pressure_id);
    printf("Average pressure = ideal + excess without tail: %f in units of kT/s^3\n", av_pressure_exc + pressure_id) ;
    av_pressure = av_pressure_exc + pressure_id + pressure_tail(rc) ;
    printf("Average pressure = ideal + excess with tail: %f in units of kT/s^3\n", av_pressure) ;
    printf("\n");


    if(N==200){
        // physical units
        printf("\n");
        printf("Temperature: %f K \n", A/(KB*beta));
        sigma_centimetri = sigma*pow(10.,-7.);
        printf("Molar density: %f mol/cm^3 \n", rho_adimensionale/(pow(sigma_centimetri,3.)*AVO));
        printf("Average energy per particle: %f eV \n", (av_energy + energy_tail(rc))/beta*A ) ;
        sigma_metri = sigma*pow(10.,-9.);
        A_Joule = A* (1.602176634 * pow(10.,-19.));
        printf("Average pressure: %f kPa\n", av_pressure*A_Joule/(beta*pow(sigma_metri,3.))/1000. ) ;
        //printf("Average pressure: %f atm\n", (av_pressure*A_Joule/(beta*pow(sigma_metri,3.)))/101325. ) ;
        printf("\n");

    }

    if(compute_g!=0){
        for(n=0; n<M; n++){
            h[n] = h[n]/Nit;
            h[n] = h[n] * pow(L,3)/shell_volume(n*epsilon, epsilon) * 2. / (N*(N-1.)) ;
            fprintf(fp1, "%d %f %f %f\n", n, n*epsilon, n*epsilon*sigma, h[n]);
        }
    }
    fclose(fp1);

}





// funzione che calcola la distanza, tutto adimensionale
double dist(double ri[3], double rj[3], double L){
    double xij, dij[3], distanza;
    int k;
    distanza=0;
    for(k=0; k<3; k++){
        xij = ri[k]-rj[k];
        dij[k] = xij - L*floor(0.5 + xij/L) ;
        distanza += dij[k]*dij[k] ;
    //  if(dij[k]>L/2){
    //     printf("errore\n");
    //  }
    }
    return sqrt(distanza);
}



// funzione che calcola il pair potential
double interaction(double r){
    double V;
    V = exp(-r)/(r*r); //adimensionale
    return V;
}



double V_tr_sh(double r, double rc){
    double V;
    if(r<rc){
        V = interaction(r) - interaction(rc);
    }
    else {
        V=0.;
    }
    return V ;
}



double V_trunc(double r, double rc){
    double V;
    if(r<rc){
        V = interaction(r);
    }
    else {
        V=0.;
    }
    return V ;
}



int Metropolis(double deltaE, double b){
  int accept;
  double Axy, ran;
  // dà 1 se la mossa è accettata, zero altrimenti
  if(deltaE<0){
    accept=1;
  }
  else {
    Axy = exp(-b*deltaE);
    if(drand48()<Axy){
      accept=1;
    }
    else{
      accept=0;
    }
  }
  return(accept);
}


double virial_trunc(double x, double xc){ // x e xc sono adimensionali
    double V;
    if(x<xc){
        V = - exp(-x)/x * (2./x + 1.) ; //per avere il viriale devo moltiplicare x*V'
    }
    else {
        V=0.;
    }
    return V ;
}


// energy per particle
double energy_tail(double rc){
    double e;
    e = 2*M_PI*rho_adimensionale*exp(-rc); // in units of A, quindi adimensionale
    e = e*beta; // in units of kT
    return e ;
}

double pressure_tail(double rc){
    double p;
    p = 2.*M_PI/3. * pow(rho_adimensionale,2) * exp(-rc) * (rc + 3.) ; //in units of sigma^3/A
    p = p * beta ; //in units of sigma^3/kT
    return p ;
}

double shell_volume(double R, double e){ // alla fine è la shell tra (n+1)*epsilon ed  n*epsilon
    double v;
    v = 4./3. * M_PI * pow((R+e), 3) - 4./3. * M_PI * pow(R, 3) ;
    return v ;
}
