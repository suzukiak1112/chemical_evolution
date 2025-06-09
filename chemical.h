//****************************************************************//
//  class chemical for galactic chemical evolution                //
//                                   written by Akihiro Suzuki    //
//****************************************************************//
#ifndef _CHEMICAL_H
#define _CHEMICAL_H
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class chemical{
public: 
    int n_cur;
    int nz;
    int nz_cc;
    int nz_ag;
    int nt;
    double idx_a;
    double M_min;
    double M_max;
    double Mc_min;
    double Mwd_min;
    double Ts;
    double Te;
    double age_sf;
    int Nt;
    int Nsp;
    double dtmp;
    int itmp,itr;
    std::string stmp;
    std::string* spec;
    // I/O
    int num_ia,num_cc,num_hn,num_ag,num_files;
    int grid_num[100];
    std::string files[100];
    std::string ext_data;
    std::stringstream ss;
    std::ifstream in;
    // parameters
    double N_snia;         // SNe Ia parameter
    double K_snia;         // SNe Ia parameter, normlization for delay time distribution
    double T_snia;         // SNe Ia parameter, delay time constant
    double T_snia_min;     // SNe Ia parameter, minimum delay time
    double T_cool;         // gas cooling timescale
    double T_sf;           // star formation timescale
    double eta_out;        // outflow parameter
    //
    // variables for galaxy evolution
    //
    double* M_star;
    double** M_cold;
    double** M_hot;
    double* dM_star;
    double** dM_cold;
    double** dM_hot;
    double** X_cold;
    double* M_hot_total;
    double* M_cold_total;
    double* M_cold_metal;
    double* age;
    double* Z_g;
    double* SFR;
    double* R_snia;
    double* R_ccsn;
    double* Xi;
    double** M_z;
    //
    // variables storing yield data
    //
    double* X_Ia;
    double Mej_Ia;
    double*** X_cc;
    double** E_cc;
    double** M_cc;
    double** Mej_cc;
    double** Mc_cc;
    double* Z_cc;
    int* N_cc;
    double*** X_hn;
    double** E_hn;
    double** M_hn;
    double** Mej_hn;
    double** Mc_hn;
    double* Z_hn;
    int* N_hn;
    double*** X_ag;
    double** M_ag;
    double** Mc_ag;
    double* Z_ag;
    int* N_ag;
    //
    // variables and parameters for stellar mass grid
    //
    int Nm;
    double norm;
    double M_u;
    double M_l;
    double*** X;
    double* M_g;
    double* dM_g;
    double* T_g;
    double** M_c;
    double** M_e;
    double** M_w;
    double*** X_e;
    double*** X_w;
    double* Z;
    double* dNdM;
    //
    // initialization of chemical class
    //
    chemical(double M_min_in,double M_max_in){
        // M_min_in = input minimum stellar mass
        // M_max_in = input maximum stellar mass
        n_cur=0;
        ext_data=".dat";
        Mc_min=1.4;
        Mwd_min=M_min_in;
        M_min=M_min_in;
        M_max=M_max_in;
    }
    void read_data(const char*,const char*,const char*);
    void init(double,double,int,double,double,int,double(*)(double));
    void find_nz(double);
    int find_nt(double);
};

//
// reading yield data
//
void chemical::read_data(const char* folder,
             const char* list_spec,
             const char* list)
{
    // read the list of species
    ss<<folder<<list_spec;
    std::cout<<ss.str()<<std::endl;
    in.open(ss.str());
    ss.str("");
    in>>stmp>>Nsp;
    spec=new std::string[Nsp];
    Xi=new double[Nsp];
    for(int is=0;is<Nsp;is++){
        in>>spec[is]>>Xi[is];
    }
    in.close();
    // read the list of files
    ss<<folder<<list;
    std::cout<<ss.str()<<std::endl;
    in.open(ss.str());
    ss.str("");
    num_files=0;
    num_ia=num_cc=num_hn=num_ag=0;
    do{
      in>>files[num_files];
      if(files[num_files].find("SNIa")==0)num_ia++;
      if(files[num_files].find("CCSN")==0)num_cc++;
      if(files[num_files].find("HNe")==0)num_hn++;
      if(files[num_files].find("AGB")==0)num_ag++;
      num_files++;
    }while(!in.eof());
    num_files--;
    // define variables
    X_Ia=new double[Nsp];  
    X_cc=new double**[num_cc];
    E_cc=new double*[num_cc];
    M_cc=new double*[num_cc];
    Mej_cc=new double*[num_cc];
    Mc_cc=new double*[num_cc];
    Z_cc=new double[num_cc];
    N_cc=new int[num_cc];
    X_hn=new double**[num_hn];
    E_hn=new double*[num_hn];
    M_hn=new double*[num_hn];
    Mej_hn=new double*[num_hn];
    Mc_hn=new double*[num_hn];
    Z_hn=new double[num_hn];
    N_hn=new int[num_hn];
    X_ag=new double**[num_ag];
    M_ag=new double*[num_ag];
    Mc_ag=new double*[num_ag];
    Z_ag=new double[num_ag];
    N_ag=new int[num_ag];
    // read yield data
    in.close();
    num_ia=num_cc=num_hn=num_ag=0;
    for(int i=0;i<num_files;i++){
      ss<<folder<<files[i]<<ext_data;
      in.open(ss.str());
      ss.str("");
      if(files[i].find("SNIa")==0){
        // read SNIa data
        std::cout<<"reading SNIa:"<<files[i]<<std::endl;
        in>>stmp>>Mej_Ia;
        for(int is=0;is<Nsp;is++){
        in>>stmp>>X_Ia[is];
        }
        num_ia++;
      }
      else if(files[i].find("CCSN")==0){
        // read CCSNe data
        std::cout<<"reading CCSNe: "<<files[i]<<std::endl;
        in>>stmp>>Z_cc[num_cc];
        in>>stmp>>N_cc[num_cc];
        std::cout<<N_cc[num_cc]<<std::endl;
        E_cc[num_cc]=new double[N_cc[num_cc]];
        M_cc[num_cc]=new double[N_cc[num_cc]];
        Mej_cc[num_cc]=new double[N_cc[num_cc]];
        Mc_cc[num_cc]=new double[N_cc[num_cc]];
        X_cc[num_cc]=new double*[N_cc[num_cc]];
        for(int im=0;im<N_cc[num_cc];im++)X_cc[num_cc][im]=new double[Nsp];
        in>>stmp;
        for(int im=0;im<N_cc[num_cc];im++)in>>E_cc[num_cc][im];
        in>>stmp;
        for(int im=0;im<N_cc[num_cc];im++)in>>M_cc[num_cc][im];
        in>>stmp;
        for(int im=0;im<N_cc[num_cc];im++)in>>Mej_cc[num_cc][im];
        in>>stmp;
        for(int im=0;im<N_cc[num_cc];im++)in>>Mc_cc[num_cc][im];
        for(int is=0;is<Nsp;is++){
        in>>stmp;
        for(int im=0;im<N_cc[num_cc];im++)in>>X_cc[num_cc][im][is];
        }
        num_cc++;
      }
      else if(files[i].find("HNe")==0){
        // read HNe data
        std::cout<<"reading HNe: "<<files[i]<<std::endl;
        in>>stmp>>Z_hn[num_hn];
        in>>stmp>>N_hn[num_hn];
        E_hn[num_hn]=new double[N_hn[num_hn]];
        M_hn[num_hn]=new double[N_hn[num_hn]];
        Mej_hn[num_hn]=new double[N_hn[num_hn]];
        Mc_hn[num_hn]=new double[N_hn[num_hn]];
        X_hn[num_hn]=new double*[N_hn[num_hn]];
        for(int im=0;im<N_hn[num_hn];im++)X_hn[num_hn][im]=new double[Nsp];
        in>>stmp;
        for(int im=0;im<N_hn[num_hn];im++)in>>E_hn[num_hn][im];
        in>>stmp;
        for(int im=0;im<N_hn[num_hn];im++)in>>M_hn[num_hn][im];
        in>>stmp;
        for(int im=0;im<N_hn[num_hn];im++)in>>Mej_hn[num_hn][im];
        in>>stmp;
        for(int im=0;im<N_hn[num_hn];im++)in>>Mc_hn[num_hn][im];
        for(int is=0;is<Nsp;is++){
        in>>stmp;
        for(int im=0;im<N_hn[num_hn];im++)in>>X_hn[num_hn][im][is];
        }
        num_hn++;
      }
      else if(files[i].find("AGB")==0){
        // read AGB wind data
        std::cout<<"reading AGB: "<<files[i]<<std::endl;
        in>>stmp>>Z_ag[num_ag];
        in>>stmp>>N_ag[num_ag];
        M_ag[num_ag]=new double[N_ag[num_ag]];
        Mc_ag[num_ag]=new double[N_ag[num_ag]];
        X_ag[num_ag]=new double*[N_ag[num_ag]];
        for(int im=0;im<N_ag[num_ag];im++)X_ag[num_ag][im]=new double[Nsp];
        in>>stmp;
        for(int im=0;im<N_ag[num_ag];im++)in>>M_ag[num_ag][im];
        in>>stmp;
        for(int im=0;im<N_ag[num_ag];im++)in>>Mc_ag[num_ag][im];
        for(int is=0;is<Nsp;is++){
        in>>stmp;
        for(int im=0;im<N_ag[num_ag];im++)in>>X_ag[num_ag][im][is];
        }
        num_ag++;
      }
      in.close();
    }
}

void chemical::init(double Ts_in,             // initial time
                    double Te_in,             // final time
                    int Nt_in,                // number of time grid
                    double Mu_in,             // upper stellar mass
                    double Ml_in,             // lower stellar mass
                    int Nm_in,                // number of stellar mass grid
                    double (*life)(double))   // function giveing stellar lifetime
{
    // define variables
    Nt=Nt_in;                            // number of time grid
    M_cold=new double*[Nsp];             // cold gas mass for Nsp species 
    M_hot=new double*[Nsp];              // hot gas mass for Nsp species
    dM_cold=new double*[Nsp];            // derivative dM/dt for cold gas mass
    dM_hot=new double*[Nsp];             // derivative dM/dt for hot gas mass
    X_cold=new double*[Nsp];             // cold gas mass fraction for Nsp species
    for(int is=0;is<Nsp;is++){           // For variables for each element is, we define Nt time grid
        M_cold[is]=new double[Nt];
        M_hot[is]=new double[Nt];
        dM_cold[is]=new double[Nt];
        dM_hot[is]=new double[Nt];
        X_cold[is]=new double[Nt];
    }
    M_cold_total=new double[Nt];         // total cold gas mass
    M_cold_metal=new double[Nt];         // total cold gas metal mass
    M_hot_total=new double[Nt];          // total hot gas mass
    M_star=new double[Nt];               // stellar mass
    dM_star=new double[Nt];              // dM/dt for stellar mass
    age=new double[Nt];                  // age of the galaxy
    Z_g=new double[Nt];                  // metal mass fractioin in the gas phase
    SFR=new double[Nt];                  // star formation rate
    R_snia=new double[Nt];               // SN Ia rate
    R_ccsn=new double[Nt];               // CCSN rate
    // mass grid
    Nm=Nm_in;
    M_u=Mu_in;
    M_l=Ml_in;
    M_g=new double[Nm];
    dM_g=new double[Nm];
    T_g=new double[Nm];
    for(int im=0;im<Nm;im++){
      M_g[im]=M_min*pow(M_max/M_min,im/(Nm-1.0));
      T_g[im]=life(M_g[im]);
    }
    dNdM=new double[Nm];
    // stellar yields
    M_c=new double*[Nm];
    M_e=new double*[Nm];
    M_w=new double*[Nm];
    X_e=new double**[Nm];
    for(int im=0;im<Nm;im++){
        if(M_g[im]<M_l){
            // LM or IM stars prodicing white dwarfs
            M_c[im]=new double[num_ag];
            M_e[im]=new double[num_ag];
            M_w[im]=new double[num_ag];
            X_e[im]=new double*[Nm];
            for(int iz=0;iz<num_ag;iz++){
                X_e[im][iz]=new double[Nsp];
                if(M_g[im]<M_ag[iz][0]){
                    M_c[im][iz]=((M_g[im]-M_l)*Mc_ag[iz][0]+(M_ag[iz][0]-M_g[im])*Mwd_min)/(M_ag[iz][0]-M_l);
                    M_e[im][iz]=M_g[im]-M_c[im][iz];
                    if(M_e[im][iz]<0.0){
                        M_e[im][iz]=0.0;
                        M_c[im][iz]=M_g[im];
                    }
                    for(int is=0;is<Nsp;is++){
                        X_e[im][iz][is]=X_ag[iz][0][is];
                    }
                }
                else if(M_g[im]<M_ag[iz][N_ag[iz]-1]){
                    for(int j=0;j<N_ag[iz];j++){
                        if((M_ag[iz][j]<=M_g[im])&&(M_g[im]<M_ag[iz][j+1])){
                            M_c[im][iz]=((M_g[im]-M_ag[iz][j])*Mc_ag[iz][j+1]+(M_ag[iz][j+1]-M_g[im])*Mc_ag[iz][j])/(M_ag[iz][j+1]-M_ag[iz][j]);
                            M_e[im][iz]=M_g[im]-M_c[im][iz];
                            if((M_g[im]-M_ag[iz][j])<(M_ag[iz][j+1]-M_g[im])){
                                for(int is=0;is<Nsp;is++){
                                    X_e[im][iz][is]=X_ag[iz][j][is];
                                }
                            }
                            else {
                                for(int is=0;is<Nsp;is++){
                                    X_e[im][iz][is]=X_ag[iz][j+1][is];
                                }
                            }
                            break;
                        }
                    }
                }
                else {
                    M_c[im][iz]=Mc_ag[iz][N_ag[iz]-1];
                    M_e[im][iz]=M_g[im]-M_c[im][iz];
                    for(int is=0;is<Nsp;is++){
                        X_e[im][iz][is]=X_ag[iz][N_ag[iz]-1][is];
                    }
                }
            }
        }
        else if(M_g[im]<M_u){
            // HM stars producing neutron stars and black holes through SNe
            M_e[im]=new double[num_cc];
            M_c[im]=new double[num_cc];
            M_w[im]=new double[num_cc];
            X_e[im]=new double*[num_cc];
            for(int iz=0;iz<num_cc;iz++){
                X_e[im][iz]=new double[Nsp];
                if(M_g[im]<M_cc[iz][0]){
                    M_c[im][iz]=((M_g[im]-M_l)*Mc_cc[iz][0]+(M_cc[iz][0]-M_g[im])*Mc_min)/(M_cc[iz][0]-M_l);
                    M_e[im][iz]=((M_g[im]-M_l)*Mej_cc[iz][0]+(M_cc[iz][0]-M_g[im])*(M_l-Mc_min))/(M_cc[iz][0]-M_l);
                    M_w[im][iz]=M_g[im]-M_e[im][iz]-M_c[im][iz];
                    if(M_w[im][iz]<0.0)M_w[im][iz]=0.0;
                    for(int is=0;is<Nsp;is++){
                        X_e[im][iz][is]=X_cc[iz][0][is];
                    }
                }
                else {
                    for(int j=0;j<N_cc[iz];j++){
                        if((M_cc[iz][j]<=M_g[im])&&(M_g[im]<M_cc[iz][j+1])){
                            M_c[im][iz]=((M_g[im]-M_cc[iz][j])*Mc_cc[iz][j+1]+(M_cc[iz][j+1]-M_g[im])*Mc_cc[iz][j])/(M_cc[iz][j+1]-M_cc[iz][j]);
                            M_e[im][iz]=((M_g[im]-M_cc[iz][j])*Mej_cc[iz][j+1]+(M_cc[iz][j+1]-M_g[im])*Mej_cc[iz][j])/(M_cc[iz][j+1]-M_cc[iz][j]);
                            M_w[im][iz]=M_g[im]-M_e[im][iz]-M_c[im][iz];
                            if(M_w[im][iz]<0.0)M_w[im][iz]=0.0;
                            if((M_g[im]-M_cc[iz][j])<(M_cc[iz][j+1]-M_g[im])){
                                for(int is=0;is<Nsp;is++){
                                    X_e[im][iz][is]=X_cc[iz][j][is];
                                }
                            }
                            else {
                                for(int is=0;is<Nsp;is++){
                                    X_e[im][iz][is]=X_cc[iz][j+1][is];
                                }
                            }
                            break;
                        }
                        else if(j==N_cc[iz]-1){
                            M_c[im][iz]=Mc_cc[iz][N_cc[iz]-1];
                            M_e[im][iz]=Mej_cc[iz][N_cc[iz]-1];
                            M_w[im][iz]=M_g[im]-M_e[im][iz]-M_c[im][iz];
                            if(M_w[im][iz]<0.0)M_w[im][iz]=0.0;
                            for(int is=0;is<Nsp;is++){
                                X_e[im][iz][is]=X_cc[iz][N_cc[iz]-1][is];
                            }
                        }
                    }
                }
            }
        }
        else {
            // direct BH formation
            M_e[im]=new double[Nm];
            M_c[im]=new double[Nm];
            M_w[im]=new double[Nm];
            X_e[im]=new double*[Nm];
            for(int iz=0;iz<num_cc;iz++){
                M_c[im][iz]=M_g[im];
                M_e[im][iz]=0.0;
                M_w[im][iz]=0.0;
                X_e[im][iz]=new double[Nsp];
                for(int is;is<Nsp;is++){
                    X_e[im][iz][is]=0.0;
                }    
            }
        }
    }
    //
    // initialize time coordinate
    //
    Ts=Ts_in;
    Te=Te_in;
    for(int n=0;n<Nt;n++){
        //age[n]=Ts+(Te-Ts)/(Nt-1)*n;
        age[n]=Ts*pow(Te/Ts,n/(Nt-1.0));
    }
}

void chemical::find_nz(double Zin)
{
    int ip,im;
    nz_cc=num_cc-1;
    for(int iz=0;iz<num_cc-1;iz++){
        if(Zin<Z_cc[iz+1]){
            nz_cc=iz;
            break;
        }
    }
    nz_ag=num_ag-1;
    for(int iz=0;iz<num_ag-1;iz++){
        if(Zin<Z_cc[iz+1]){
            nz_ag=iz;
            break;
        }
    }  
}

int chemical::find_nt(double age_in)
{
    int ip,im;
    ip=Nt;
    im=0;
    do{
        nt=(ip+im)/2;
        if(age_in<age[nt]){
            ip=nt;
        }
        else {
            im=nt;
        }
    }while(ip>im+1);
    nt=im;
    return nt;
}

#endif
