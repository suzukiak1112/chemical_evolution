#include  <math.h>
#include  <stdlib.h>
//#include  <time.h>
#include  <sys/time.h>
#include  <iostream>
#include  <fstream>
#include  <sstream>
#include "physc.h"
#include "parameter.h"
#include "chemical.h"
using namespace std;

//
// a function returning the lifetime of the star with the initial mass Mstar
//
double lifetime(double Mstar)
{
    double tau;
    if(Mstar<0.6){
        tau=160.0;
    }
    else if(Mstar<=6.6) {
        tau=1.790-0.2232*(7.764-log10(Mstar));
        tau=pow(10.0,(0.334-sqrt(tau))/0.1116);
    }
    else {
        tau=1.2*pow(Mstar,-1.85)+0.003;
    }
    return tau;// return lifetime tau(M_star)
}

//
// a function returning the delay time of SNe Ia
//
double delay_time_distribution(chemical & gce,
                               double td)
{
  //return gce.K_snia*exp(-(td-gce.T_snia_min)/gce.T_snia)/gce.T_snia;
  if(td<gce.T_snia_min) return 0.0;
  else return gce.K_snia*pow(td,-1.0);
}

//
// a function setting parameters and initial condition for galactic chemical evolution (gce)
//
void set_params(chemical &gce)
{
    gce.n_cur=0;                 // n_cur: current step number. initially zero
    //
    // time scales
    //
    gce.T_sf=1.0;                // star formation timescale, T_sf=1.0 means the timescale is set to to 1 Gyr
    gce.T_cool=3.0;              // gas cooling timescale, T_cool=3.0 means the timescale is set to 3 Gyr
    //
    // initial conditions
    //
    gce.M_hot_total[0]=4.0e11;   // initial total hot gas mass is set to 4e11 Msun
    gce.M_cold_total[0]=1.0e7;   // initial total cold gas mass is set to 1e7 Msun
    for(int is=0;is<gce.Nsp;is++){ // this loop set initial element masses (is=0 to Nsp) 
      gce.X_cold[is][0]=gce.Xi[is];   // gce.Xi[is] is the initial mass fraction of is'th element
      gce.M_hot[is][0]=gce.Xi[is]*gce.M_hot_total[0];
      gce.M_cold[is][0]=gce.Xi[is]*gce.M_cold_total[0];
    }
    //
    // Initial Mass Fucntion
    //
#ifdef SALPETER    // Salpeter IMF
    gce.idx_a=1.35;
    for(int im=0;im<gce.Nm;im++){
      gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a);
    }
#endif
#ifdef KROUPA01      // Kroupa (2001) IMF
    for(int im=0;im<gce.Nm;im++){
        if(gce.M_g[im]<0.1){
            gce.idx_a=-0.7; 
            gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a)/pow(0.5,-1.0-0.3)*pow(0.5,-1.0-1.3)/pow(0.1,-1.0+0.7)*pow(0.1,-1.0-0.3);
        }
        else if(gce.M_g[im]<0.5){
            gce.idx_a=0.3;
            gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a)/pow(0.5,-1.0-gce.idx_a)*pow(0.5,-1.0-1.3);
        }
        else {
            gce.idx_a=1.3;
            gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a);
        }
    }
#endif
#ifdef SCALO86       // Scalo (1986)
    for(int im=0;im<gce.Nm;im++){
        if(gce.M_g[im]<2.0){
            gce.idx_a=1.35;
            gce.dNdM[im]=pow(gce.M_g[im]/2.0,-1.0-gce.idx_a)/pow(2.0,-1.0-gce.idx_a)*pow(2.0,-1.0-1.7);
        }
        else {
            gce.idx_a=1.7;
            gce.dNdM[im]=pow(gce.M_g[im]/2.0,-1.0-gce.idx_a);
        }
    }
#endif
#ifdef KROUPA93       // Kroupa (1993)
    for(int im=0;im<gce.Nm;im++){
        if(gce.M_g[im]<0.5){
            gce.idx_a=-0.3;
            gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a)/pow(0.5,-1.0-gce.idx_a)*pow(0.5,-1.0-1.2);
        }
        else if(gce.M_g[im]<1.0){
            gce.idx_a=1.2;
            gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a);
        }
        else {
            gce.idx_a=1.7;
            gce.dNdM[im]=pow(gce.M_g[im],-1.0-gce.idx_a);
        }
    }
#endif
    //
    // normalizing the initial mass function
    //  so that the integration of M*dN/M gives 1Msol. 
    //
    for(int im=0;im<gce.Nm;im++){  // gce.Nm is the maximum number of the mass grid
      if(im==0)gce.dM_g[im]=0.5*(gce.M_g[im+1]-gce.M_g[im]);
      else if(im<gce.Nm-1)gce.dM_g[im]=0.5*(gce.M_g[im+1]-gce.M_g[im-1]);
      else gce.dM_g[im]=0.5*(gce.M_g[im]-gce.M_g[im-1]);
    }
    gce.norm=0.0;
    for(int im=0;im<gce.Nm;im++){
      gce.norm+=gce.M_g[im]*gce.dNdM[im]*gce.dM_g[im];
    }
    for(int im=0;im<gce.Nm;im++)gce.dNdM[im]/=gce.norm;
    //
    // integration of the initial mass function in different mass ranges
    //
    double sum0=0.0;
    double sum1=0.0;
    double sum2=0.0;
    double sum3=0.0;
    double sum4=0.0;
    for(int im=0;im<gce.Nm;im++){
      sum0+=gce.dNdM[im]*gce.dM_g[im];
      if((9.0<=gce.M_g[im])&&(gce.M_g[im]<=17.0))sum1+=gce.dNdM[im]*gce.dM_g[im];
      if((25.0<=gce.M_g[im]))sum2+=gce.dNdM[im]*gce.dM_g[im];
      if((30.0<=gce.M_g[im]))sum3+=gce.dNdM[im]*gce.dM_g[im];
      if((60.0<=gce.M_g[im]))sum4+=gce.dNdM[im]*gce.dM_g[im];
    }
    cout<<sum0<<" "<<sum1<<" "<<sum2<<" "<<sum3<<" "<<sum2/sum1<<" "<<sum3/sum1<<" "<<sum4/sum1<<endl;
    //
    // SNIa parameters
    //
    gce.T_snia=1.0;
    gce.T_snia_min=0.05;
    gce.N_snia=1.5e-3;
    gce.K_snia=gce.N_snia/log(gce.Te/gce.T_snia_min);
    //
    // outflow parameter
    //
    gce.eta_out=2.5;
}


void derivative(chemical &gce)
{
    int i_s,i_e,iz;
    double dt,t_s,t_e,dt_sub,dtd;
    double limitter;
    //
    // initialize the dataset SFT and derivative is set to zero
    // the system is updated from gce.n_cur to gce.n_cur+1 separated by timestep dt
    dt=gce.age[gce.n_cur+1]-gce.age[gce.n_cur];
    gce.SFR[gce.n_cur]=0.0;
    gce.dM_star[gce.n_cur]=0.0;
    for(int is=0;is<gce.Nsp;is++){
        gce.dM_hot[is][gce.n_cur]=0.0;
        gce.dM_cold[is][gce.n_cur]=0.0;
    }
    //
    // star formation rate
    //
    gce.SFR[gce.n_cur]=gce.M_cold_total[gce.n_cur]/gce.T_sf;
    //
    // inflow
    //
    for(int is=0;is<gce.Nsp;is++){
        gce.dM_hot[is][gce.n_cur]-=gce.M_hot[is][gce.n_cur]/gce.T_cool*dt;
        gce.dM_cold[is][gce.n_cur]+=gce.M_hot[is][gce.n_cur]/gce.T_cool*dt;
        //gce.dM_cold[is][gce.n_cur]+=gce.age[gce.n_cur]*gce.M_hot[is][gce.n_cur]/(gce.T_cool*gce.T_cool)*exp(-gce.age[gce.n_cur]/gce.T_cool)*dt;
    }
    //
    // star formation associated cold mass decrease and stellar mass increase
    //
    for(int is=0;is<gce.Nsp;is++){
        gce.dM_cold[is][gce.n_cur]-=gce.M_cold[is][gce.n_cur]/gce.T_sf*dt;
        gce.dM_star[gce.n_cur]+=gce.M_cold[is][gce.n_cur]/gce.T_sf*dt;
    }
    //
    // outflow reducing the cold gass mass
    //
    for(int is=0;is<gce.Nsp;is++){
        gce.dM_cold[is][gce.n_cur]-=gce.M_cold[is][gce.n_cur]/gce.M_cold_total[gce.n_cur]*gce.SFR[gce.n_cur]*gce.eta_out*dt;
        //gce.dM_hot[is][gce.n_cur]+=gce.M_cold[is][gce.n_cur]/gce.M_cold_total[gce.n_cur]*gce.SFR[gce.n_cur]*gce.eta_out*dt;
    }
    //
    // CCSNe
    //
    gce.R_ccsn[gce.n_cur]=0.0;
    for(int im=0;im<gce.Nm-1;im++){  // thi loop run throug hthe stellar mass grid
        // limitter for limiting the explosion mass range
        limitter=1.0;
        //if((8.0<gce.M_g[im])&&(gce.M_g[im]<20.0))limitter=0.0;
        //if(17.0<gce.M_g[im])limitter=0.0;  // model B
        //if((17.0<gce.M_g[im])&&(gce.M_g[im]<20.0))limitter=0.0;  // model C
        //if((17.0<gce.M_g[im])&&(gce.M_g[im]<30.0))limitter=0.0;  // model D
        //if((17.0<gce.M_g[im])&&(gce.M_g[im]<40.0))limitter=0.0;  // model E
        //if(50.0<gce.M_g[im])limitter=0.0;  // dependence on mass range
        if(gce.age[gce.n_cur+1]<lifetime(gce.M_g[im])){
            // if the age of the galaxy is less than the life time of stars with mass gce.M_g[im]
            // no contribution is expected
            t_s=0.0;
            i_s=0;
            t_e=0.0;
            i_e=0;
        }
        else {
            if(gce.age[gce.n_cur]<lifetime(gce.M_g[im])){
                t_s=0.0;
                i_s=0;
                t_e=gce.age[gce.n_cur+1]-lifetime(gce.M_g[im]);
                i_e=gce.find_nt(t_e)+1;
            }
            else {
                t_s=gce.age[gce.n_cur]-lifetime(gce.M_g[im]);
                i_s=gce.find_nt(t_s);
                t_e=gce.age[gce.n_cur+1]-lifetime(gce.M_g[im]);
                i_e=gce.find_nt(t_e)+1;
            }
        }
        for(int it=i_s;it<i_e;it++){
            // set model index
            gce.find_nz(gce.Z_g[it]);
            if(gce.M_g[im]<gce.M_l)iz=gce.nz_ag;
            else {
                iz=gce.nz_cc;
                if(iz==0)iz=1;
            }
            //
            if(gce.age[it]<t_s)dt_sub=t_s;
            else dt_sub=gce.age[it];
            if(t_e<gce.age[it+1])dt_sub=t_e-dt_sub;
            else dt_sub=gce.age[it+1]-dt_sub;
            for(int is=0;is<gce.Nsp;is++){
                gce.dM_cold[is][gce.n_cur]+=gce.X_e[im][iz][is]*gce.M_e[im][iz]*gce.SFR[it]*dt_sub*gce.dNdM[im]*gce.dM_g[im]*limitter;
                gce.dM_cold[is][gce.n_cur]+=gce.X_cold[is][it]*gce.M_w[im][iz]*gce.SFR[it]*dt_sub*gce.dNdM[im]*gce.dM_g[im]*limitter;
            }
            gce.dM_star[gce.n_cur]-=gce.M_e[im][iz]*gce.SFR[it]*dt_sub*gce.dNdM[im]*gce.dM_g[im]*limitter;
            gce.dM_star[gce.n_cur]-=gce.M_w[im][iz]*gce.SFR[it]*dt_sub*gce.dNdM[im]*gce.dM_g[im]*limitter;
            gce.R_ccsn[gce.n_cur]+=gce.SFR[it]*dt_sub*gce.dNdM[im]*gce.dM_g[im]*limitter;
        }
    }
    // SNIa
    gce.R_snia[gce.n_cur]=0.0;
    if(gce.age[gce.n_cur+1]<gce.T_snia_min){
        // no contribution
    }
    else {
        i_e=gce.find_nt(gce.age[gce.n_cur+1]-gce.T_snia_min);
        for(int it=0;it<i_e+1;it++){
            dtd=delay_time_distribution(gce,gce.age[gce.n_cur]-gce.age[it]);
            if(it==i_e)dt_sub=gce.age[gce.n_cur+1]-gce.age[it];
            else dt_sub=gce.age[it+1]-gce.age[it];
            for(int is=0;is<gce.Nsp;is++){
                gce.dM_cold[is][gce.n_cur]+=gce.X_Ia[is]*gce.Mej_Ia*dtd*gce.SFR[it]*dt_sub*dt;
            }
            gce.R_snia[gce.n_cur]+=dtd*gce.SFR[it]*dt_sub;
        }
    }
}

//
// fuction advance to push the system from t -> t+dt
//
void advance(chemical &gce)
{
    //
    // stellar mass
    //
    gce.M_star[gce.n_cur+1]=gce.M_star[gce.n_cur]+gce.dM_star[gce.n_cur];
    //
    // hot and cold gas mass for species is
    //
    for(int is=0;is<gce.Nsp;is++){
        gce.M_hot[is][gce.n_cur+1]=gce.M_hot[is][gce.n_cur]+gce.dM_hot[is][gce.n_cur];
        gce.M_cold[is][gce.n_cur+1]=gce.M_cold[is][gce.n_cur]+gce.dM_cold[is][gce.n_cur];
    }
    gce.n_cur++;
    //
    // computing the total hot and cold gass mass
    //
    gce.M_hot_total[gce.n_cur]=0.0;
    gce.M_cold_total[gce.n_cur]=0.0;
    gce.M_cold_metal[gce.n_cur]=0.0;
    for(int is=0;is<gce.Nsp;is++){
        gce.M_cold_total[gce.n_cur]+=gce.M_cold[is][gce.n_cur];
        gce.M_hot_total[gce.n_cur]+=gce.M_hot[is][gce.n_cur];
        if(is>1){
            gce.M_cold_metal[gce.n_cur]+=gce.M_cold[is][gce.n_cur];
        }
    }
    //
    // mass fraction of species is in the cold gas phase
    //
    for(int is=0;is<gce.Nsp;is++){
        gce.X_cold[is][gce.n_cur]=gce.M_cold[is][gce.n_cur]/gce.M_cold_total[gce.n_cur];
    }
    //
    // metallicity in the cold gas phase
    //
    gce.Z_g[gce.n_cur]=gce.M_cold_metal[gce.n_cur]/gce.M_cold_total[gce.n_cur];
}


int main()
{
    //
    // output filenames
    //
    ofstream abundance("abundance.txt");
    ofstream evolution("evolution.txt");
    //
    // define the object gce with specified mass range [mass_min,mass_max]
    // loading data from  "../data/list_species.txt" "../data/list_models_cl04.txt"
    //
    chemical gce(mass_min,mass_max);
    gce.read_data("../data/","list_species.txt","list_models_cl04.txt");
    //gce.read_data("./data/","list_species.txt","list_models.txt");
    //gce.read_data("./data/","list_species.txt","list_models_lc06.txt");
    //
    // initialize the object gce by function init
    // you need to specifi the following parameters
    // time_s   : initial age in Gyr
    // time_e   : final age in Gyr
    // num_time : number of integration in [time_s,time_e]
    // mass_up  : upper mass limit
    // mass_lw  : lower mass limit
    // num_mass : number of mass grid
    // lifetime : another class for the lifetimes of massive stars 
    //
    gce.init(time_s,time_e,num_time,mass_up,mass_lw,num_mass,lifetime);
    set_params(gce);
    //
    // before the integration, we check the number of expected number of CCSNe per unit mass Nccsn
    //
    double Nccsn=0.0;
    for(int im=0;im<gce.Nm;im++){
        if(im!=0){
            if(gce.M_g[im-1]>9.0){  // the threshold mass for CCSN is set to 9Msun here
                Nccsn+=(gce.M_g[im]-gce.M_g[im-1])*0.5*(gce.dNdM[im]+gce.dNdM[im-1]);
            }
        }
    }
    //
    // start integrating the system from nt=0 to nt=num_time-2
    //
    for(int nt=0;nt<num_time-1;nt++){
        // dervative(gce) gives dM/dt, ... , and so on        
        derivative(gce);
        // advance(gce) push the system forward from n_curr to n_curr+1
        advance(gce);
        // output in every 100 step
        if(nt%(100)==0){
            // gce.n_cur         :   current integration number. gce.qu[]
            // gce.age           :   current age of the system [Gyr]
            // gce.M_cold_total  :   current mass of the cold gas
            // gce.Z_g           :   current metaliccity
            cout<<gce.n_cur<<" "<<gce.age[gce.n_cur]<<" "<<gce.M_cold_total[gce.n_cur]<<" "<<gce.Z_g[gce.n_cur];
            cout<<" "<<log10(gce.M_cold[25][gce.n_cur]/gce.M_cold[0][gce.n_cur]);
            cout<<" "<<log10(gce.M_cold[7][gce.n_cur]/gce.M_cold[0][gce.n_cur]);
            cout<<" "<<gce.R_ccsn[nt]/gce.R_snia[nt];
            cout<<endl;
        }
        evolution<<gce.age[nt];                    // age in Gyr
        evolution<<" "<<gce.SFR[nt];               // star-formation rate in [Msun/yr]
        evolution<<" "<<gce.M_cold_total[nt];      // current total mass of the cold gas
        evolution<<" "<<gce.M_cold_metal[nt];      // current total metal mass in the cold gas
        evolution<<" "<<gce.M_star[nt];            // current stellar mass
        evolution<<" "<<gce.M_hot_total[nt];
        evolution<<" "<<gce.R_ccsn[nt]/1.0e9;
        evolution<<" "<<gce.R_snia[nt]/1.0e9;
        evolution<<endl;
        for(int is=0;is<gce.Nsp;is++){
            abundance<<gce.M_cold[is][gce.n_cur]<<" ";
        }
        abundance<<endl;
        // 
        // some output for results: number ratio of O,Fe and so on
        //
        if((gce.age[nt]<time_sol)&&(time_sol<=gce.age[nt+1])){  
            // when the age exceeds the formation time of the sun, the code output the elemental abundance at the time
            double N_hd,N_ox,N_fe;
            N_hd=gce.M_cold[0][gce.n_cur];           // M_cold[0][*] is the array for hydrogen in cold gas
            N_ox=gce.M_cold[7][gce.n_cur]/15.999;    // M_cold[7][*] is the array for oxygen in cold gas
            N_fe=gce.M_cold[25][gce.n_cur]/55.845;   // M_cold[25][*] is the array for iron in cold gas
            cout<<gce.age[gce.n_cur]<<" "<<N_hd<<" "<<N_ox<<" "<<N_fe<<endl;
            cout<<"log10(e(O))="<<log10(N_ox/N_hd)+12.0<<endl;
            cout<<"log10(e(Fe))="<<log10(N_fe/N_hd)+12.0<<endl;
            cout<<"[O/Fe]="<<log10(N_ox)-log10(N_fe)-8.69+7.50<<endl;
        }
    }
    //
    // some output for checking assume parameters
    // 
    cout<<gce.n_cur<<endl;
    cout<<gce.K_snia<<" "<<gce.T_snia<<endl;
    cout<<" "<<gce.R_ccsn[num_time-2]/1.0e9*100.0;
    cout<<" "<<gce.R_snia[num_time-2]/1.0e9*100.0;
    cout<<" "<<gce.R_ccsn[num_time-2]/gce.R_snia[num_time-2]<<endl;
    return 0;
}

