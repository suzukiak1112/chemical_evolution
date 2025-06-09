#ifndef  _PARAMETER_H_INCLUDE
#define  _PARAMETER_H_INCLUDE

// parameters
#define time_s 1.0e-7         // initial time in Gyr
#define time_e 13.8           // final time in Gyr
#define time_sol (time_e-4.6) // solar formation in Gyr
#define num_time 10000        // time steps
#define mass_max 100.0        // maximum stellar mass in M_solar
#define mass_min 0.08         // minimum stellar mass in M_solar
#define mass_up 100.0         // maximum stellar mass for CCSN in M_solar
#define mass_lw 9.0           // minimum stellar mass for CCSN in M_solar
#define num_mass (256)  

// IMF switcher, uncomment one of them
//#define SALPETER
#define KROUPA01
//#define KROUPA93
//#define SCALO86

// some parameter for numerical integration
#define TOL 1.0e-10
#define Ng (1024)
#define Npar 10

#endif
