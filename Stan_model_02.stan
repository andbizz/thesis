//--- Stan function to obtain the distribution of the parameters in the model

//--- Define the ODE's system, whRDh a function called sir

functions{
  real[] sir(real t,       // the time t is a real variable
             real[] y,     // the vector field is made with real funtions
             real[] theta, // the parameters vector
             real[] x_r,   // these are vabiables that depend from the parameters
             int[] x_i)    // integer variables, it is used to pass the populations from parameter to input 
{ 
//--- Definition of the variables from the vector input (y) to the various componentime_sequence (S,E,I,R,H,RD)
//                                                      (x_i) to the various componentime_sequence (N_0_19,...)
//--- Group 0-19
             real S_0_19 = y[1];
             real E_0_19 = y[2];
             real I_0_19 = y[3];
             real H_0_19 = y[4];
             
             real N_0_19 = x_i[1];
//--- Group 20-39
             real S_20_39 = y[5];
             real E_20_39 = y[6];
             real I_20_39 = y[7];
             real H_20_39 = y[8];
             
             real N_20_39 = x_i[2];
//--- Group 40-59
             real S_40_59 = y[9];
             real E_40_59 = y[10];
             real I_40_59 = y[11];
             real H_40_59 = y[12];
             
             real N_40_59 = x_i[3];
//--- Group 60-79             
             real S_60_79 = y[13];
             real E_60_79 = y[14];
             real I_60_79 = y[15];
             real H_60_79 = y[16];

             real N_60_79 = x_i[4];
//--- Group over80             
             real S_over80 = y[17];
             real E_over80 = y[18];
             real I_over80 = y[19];
             real H_over80 = y[20];

             real N_over80 = x_i[5];
//--- Definition of the parameters (R0, sigma, gamma, omega_*, delta)             
// R0:
             real beta = theta[1];
// sigma:             
             real sigma = 0.2;
// Gamma:             
             real gamma = 0.5;

// Omega:             
             real omega_0_19 = theta[2];
             real omega_20_39 = theta[3];
             real omega_40_59 = theta[4];
             real omega_60_79 = theta[5];
             real omega_over80 = theta[6];

// Delta:             


//--- Definition of the componentime_sequence of the contact matrix, they are used to compute the force of infection             
// pre capital contactime_sequence of the the group 0-19 with the others groups:
             real c_0_19_0_19  = x_r[1];
             real c_0_19_20_39 = x_r[2];
             real c_0_19_40_59 = x_r[3];
             real c_0_19_60_79 = x_r[4];
             real c_0_19_over80 = x_r[5];
// pre capital contactime_sequence of the the group 20-39 with the others groups:           
             real c_20_39_0_19 = x_r[6];
             real c_20_39_20_39 = x_r[7];
             real c_20_39_40_59 = x_r[8];
             real c_20_39_60_79 = x_r[9];
             real c_20_39_over80 = x_r[10];
// pre capital contactime_sequence of the the group 40-59 with the others groups:             
             real c_40_59_0_19 = x_r[11];
             real c_40_59_20_39 = x_r[12];
             real c_40_59_40_59 = x_r[13];
             real c_40_59_60_79 = x_r[14];
             real c_40_59_over80 = x_r[15];
// pre capital contactime_sequence of the the group 60-79 with the others groups:             
             real c_60_79_0_19 = x_r[16];
             real c_60_79_20_39 = x_r[17];
             real c_60_79_40_59 = x_r[18];
             real c_60_79_60_79 = x_r[19];
             real c_60_79_over80 = x_r[20];
// pre capital contactime_sequence of the the group over80 with the others groups:             
             real c_over80_0_19 = x_r[21];
             real c_over80_20_39 = x_r[22];
             real c_over80_40_59 = x_r[23];
             real c_over80_60_79 = x_r[24];
             real c_over80_over80 = x_r[25];
             
//--- Definiction of beta

//--- Definition of the force of infection
// Force of infection for group 0-19:
             real lambda_0_19 = beta * (c_0_19_0_19 .* I_0_19 ./N_0_19  + c_0_19_20_39 .* I_20_39 ./N_20_39
                                      + c_0_19_40_59 * I_40_59 ./N_40_59 + c_0_19_60_79 * I_60_79 ./N_60_79
                                      + c_0_19_over80 * I_over80 ./N_over80) ;
// Force of infection for group 20-39:             
             real lambda_20_39 = beta * (c_20_39_0_19 * I_0_19/N_0_19  + c_20_39_20_39 * I_20_39/N_20_39
                                       + c_20_39_40_59 * I_40_59/N_40_59 + c_20_39_60_79 * I_60_79/N_60_79
                                       + c_20_39_over80 * I_over80/N_over80) ;
// Force of infection for group 40-59:             
             real lambda_40_59 = beta * (c_40_59_0_19 * I_0_19/N_0_19 + c_40_59_20_39 * I_20_39/N_20_39
                                       + c_40_59_40_59 * I_40_59/N_40_59 + c_40_59_60_79 * I_60_79/N_60_79  
                                       + c_40_59_over80 * I_over80/N_over80) ;
// Force of infection for group 60-79:             
             real lambda_60_79 = beta * (c_60_79_0_19 * I_0_19/N_0_19 + c_60_79_20_39 * I_20_39/N_20_39
                                       + c_60_79_40_59 * I_40_59/N_40_59 + c_60_79_60_79 * I_60_79/N_60_79 
                                       + c_60_79_over80 * I_over80/N_over80) ;
// Force of infection for group over80:             
             real lambda_over80 = beta * (c_over80_0_19 * I_0_19/N_0_19 + c_over80_20_39 * I_20_39/N_20_39
                                       + c_over80_40_59 * I_40_59/N_40_59 + c_over80_60_79 * I_60_79/N_60_79 
                                       + c_over80_over80 * I_over80/N_over80) ;


//--- Definition of the differential equations
// Group 0-19:
real dS_0_19_dt = - lambda_0_19 * S_0_19              ;
real dE_0_19_dt = + lambda_0_19 * S_0_19 - sigma * E_0_19;
real dI_0_19_dt =                        + sigma * E_0_19 - gamma * I_0_19;
real dH_0_19_dt = + gamma * omega_0_19 * I_0_19;


// Group 20-39:
real dS_20_39_dt = - lambda_20_39 * S_20_39              ;
real dE_20_39_dt = + lambda_20_39 * S_20_39 - sigma * E_20_39;
real dI_20_39_dt =                        + sigma * E_20_39 - gamma * I_20_39;
real dH_20_39_dt = + gamma * omega_20_39 * I_20_39;

// Group 40-59:
real dS_40_59_dt = - lambda_40_59 * S_40_59              ;
real dE_40_59_dt = + lambda_40_59 * S_40_59 - sigma * E_40_59;
real dI_40_59_dt =                        + sigma * E_40_59 - gamma * I_40_59;
real dH_40_59_dt = + gamma * omega_40_59 * I_40_59;

// Group 60-79:
real dS_60_79_dt = - lambda_60_79 * S_60_79              ;
real dE_60_79_dt = + lambda_60_79 * S_60_79 - sigma * E_60_79;
real dI_60_79_dt =                        + sigma * E_60_79 - gamma * I_60_79;
real dH_60_79_dt = + gamma * omega_60_79 * I_60_79;

// Group over80:
real dS_over80_dt = - lambda_over80 * S_over80              ;
real dE_over80_dt = + lambda_over80 * S_over80 - sigma * E_over80;
real dI_over80_dt =                        + sigma * E_over80 - gamma * I_over80;
real dH_over80_dt = + gamma * omega_over80 * I_over80;

// print("lambda_0_19: ",lambda_0_19); 
// print("lambda_40_59: ",lambda_40_59); 
// print("lambda_60_79: ",lambda_60_79); 


return{dS_0_19_dt, dE_0_19_dt, dI_0_19_dt, dH_0_19_dt,
       dS_20_39_dt, dE_20_39_dt, dI_20_39_dt,dH_20_39_dt,
       dS_40_59_dt, dE_40_59_dt, dI_40_59_dt,dH_40_59_dt,  
       dS_60_79_dt, dE_60_79_dt, dI_60_79_dt, dH_60_79_dt, 
       dS_over80_dt, dE_over80_dt, dI_over80_dt,dH_over80_dt};
}
}
//-----------------------------------
//--- Data block --------------------
//-----------------------------------
// Inside this block we define the data, it is executed per chain, the variables are global and 
// it doesn't affect the posterior distribution
data { int<lower = 1> n_simulation; // n_simulation is an integer data greater than 1
       int<lower = 1> n_fitting;

          // the initial conditions are inside a vector y0 with length 30
                              // an real component, in this assigment rstan defines it as an array
       
       real initial_time;    
       
       real time_sequence[n_simulation];       // the discretize time period is # of days long
       
       int N_data[5];         // the size of population in the 5 groups
       
       int H_0_19_data[n_fitting];  // the data of each group, they are used in the likelihood
       int H_20_39_data[n_fitting];
       int H_40_59_data[n_fitting];
       int H_60_79_data[n_fitting];
       int H_over80_data[n_fitting];
       
       real contact_matrix[5,5]; // the contact matrix is 5x5 size
}

//-------------------------------
//---Transformed Data block------
//-------------------------------
// Inside this block we can define some variables real and integer that are a trasformation of the data.
// It is executed per chain, the variables are global and in doesn't modify the posterior distribution
transformed data {
// I pass the contact matrix (form the data) as a real vector  x_r that is a variable that depends from the data
   real x_r[25] =  {contact_matrix[1,1],contact_matrix[1,2],contact_matrix[1,3],contact_matrix[1,4],contact_matrix[1,5],
                     contact_matrix[2,1],contact_matrix[2,2],contact_matrix[2,3],contact_matrix[2,4],contact_matrix[2,5],
                     contact_matrix[3,1],contact_matrix[3,2],contact_matrix[3,3],contact_matrix[3,4],contact_matrix[3,5],
                     contact_matrix[4,1],contact_matrix[4,2],contact_matrix[4,3],contact_matrix[4,4],contact_matrix[4,5],
                     contact_matrix[5,1],contact_matrix[5,2],contact_matrix[5,3],contact_matrix[5,4],contact_matrix[5,5]
                     };
// I pass the size population                  
   int x_i[5] = {N_data[1],N_data[2],N_data[3],N_data[4],N_data[5]};
}

//-------------------------------
//---Parameters block------------
//-------------------------------
// We define the parameters in this block. The parameters defined are global and they are saved.
parameters {
// phi, the dispersion parameter in the negative binomial distribution is positive
// since I want to fit 5 trajectories I assume 5 different phi
  real<lower=0> phi_inv;
// I assume all the parameters are positive   
  real<lower=0> beta;
  
  real<lower=0> omega_0_19;
  real<lower=0> omega_20_39;
  real<lower=0> omega_40_59;
  real<lower=0> omega_60_79;
  real<lower=0> omega_over80;
  
  real<lower=0> i0;
}

//------------------------------------
//---Trasformed Parameters block------
//------------------------------------
// In this block we can define variables depending from the parameters 
// It is executed many times: each leapfrog step, the variables are global and it
// doesn't affect the posterior distribution
transformed parameters{

//  the solution of the ode's system is a function that depends from the parameters 
 real y_ode[n_simulation, 20]; 
  real phi = 1. / phi_inv;
  real y0[20];
  real theta[6];
  real incidence_0_19[n_simulation];
  real incidence_20_39[n_simulation];
  real incidence_40_59[n_simulation];
  real incidence_60_79[n_simulation];
  real incidence_over80[n_simulation];
  
    y0 = {N_data[1],0,0,0,
          N_data[2],0,0,0,
          N_data[3]-i0 ./3,0,i0 ./3,0,
          N_data[4]-i0 ./3,0,i0 ./3,0,
          N_data[5]-i0 ./3,0,i0 ./3,0};
    
      
    theta[1] = beta;
    
    theta[2] = omega_0_19;
    theta[3] = omega_20_39;
    theta[4] = omega_40_59;
    theta[5] = omega_60_79;
    theta[6] = omega_over80;

    
// Compute the solution of the ode's system after the setting of the parameters in theta
 y_ode = integrate_ode_rk45(sir, y0, initial_time, time_sequence, theta, x_r, x_i);

  incidence_0_19[1] = 0;
  incidence_20_39[1] = 0;
  incidence_40_59[1] = 0;
  incidence_60_79[1] = 0;
  incidence_over80[1] = 0;

 for (i in 2:n_simulation){
  incidence_0_19[i] = y_ode[i,4]-y_ode[i-1,4];
  incidence_20_39[i] = y_ode[i,8]-y_ode[i-1,8];
  incidence_40_59[i] = y_ode[i,12]-y_ode[i-1,12];
  incidence_60_79[i] = y_ode[i,16]-y_ode[i-1,16];
 incidence_over80[i] = y_ode[i,20]-y_ode[i-1,20];
 }
 
}


//------------------------------------
//---Model block----------------------
//------------------------------------
// In this block we set the priors and the likelihood. It is executed each leapfrog step, 
// the variables defined in this block are local and it is the unique block that affectime_sequence
// the posterior distribution
model {
//Priors:

  beta ~ uniform(0.015,0.02);

  omega_0_19 ~ uniform(0.00002,0.0001);
  omega_20_39 ~ uniform(0.0005,0.0008);
  omega_40_59 ~ uniform(0.0025,0.0035);
  omega_60_79 ~ uniform(0.0075,0.009);
  omega_over80 ~ uniform(0.012,0.015);
  
  i0 ~ uniform(500, 1000);   // ( prob of to need intensive care )

  phi_inv ~ uniform(0.001,100);


// Likelihood:
// I suppose the data distributed following a negative binomial with mean the solution of the ode's system
// and variance the parameter phi

// for (i in 1:n_fitting){
//   H_0_19_data[i]  ~  neg_binomial_2((y_ode[16+i,3] +y_ode[17+i,3] +y_ode[18+i,3]
//                                  +y_ode[19+i,3] +y_ode[20+i,3] +y_ode[21+i,3] +y_ode[22+i,3]) ./7 .* omega_0_19 ,phi);
// 
//   H_20_39_data[i] ~  neg_binomial_2((y_ode[16+i,6] +y_ode[17+i,6] +y_ode[18+i,6]
//                                  +y_ode[19+i,6] +y_ode[20+i,6] +y_ode[21+i,6] +y_ode[22+i,6]) ./7 .* omega_20_39 ,phi);
// 
//   H_40_59_data[i] ~  neg_binomial_2((y_ode[16+i,9] +y_ode[17+i,9] +y_ode[18+i,9]
//                                  +y_ode[19+i,9] +y_ode[20+i,9] +y_ode[21+i,9] +y_ode[22+i,9]) ./7 .* omega_40_59 ,phi);
// 
//   H_60_79_data[i] ~  neg_binomial_2((y_ode[16+i,12] +y_ode[17+i,12] +y_ode[18+i,12]
//                                  +y_ode[19+i,12] +y_ode[20+i,12] +y_ode[21+i,12] +y_ode[22+i,12]) ./7 .* omega_60_79 ,phi);
// 
//   H_over80_data[i] ~ neg_binomial_2((y_ode[16+i,15] +y_ode[17+i,15] +y_ode[18+i,15]
//                                  +y_ode[19+i,15] +y_ode[20+i,15] +y_ode[21+i,15] +y_ode[22+i,15]) ./7 .* omega_over80  ,phi);
//                                  }

for (i in 1:n_fitting){
  H_0_19_data[i]  ~  neg_binomial_2((incidence_0_19[16+i] +incidence_0_19[17+i] +incidence_0_19[18+i]
                                 +incidence_0_19[19+i] +incidence_0_19[20+i] +incidence_0_19[21+i] +incidence_0_19[22+i]) ./7 ,phi);

  H_20_39_data[i] ~  neg_binomial_2((incidence_20_39[16+i] +incidence_20_39[17+i] +incidence_20_39[18+i]
                                 +incidence_20_39[19+i] +incidence_20_39[20+i] +incidence_20_39[21+i] +incidence_20_39[22+i]) ./7 ,phi);

  H_40_59_data[i] ~  neg_binomial_2((incidence_40_59[16+i] +incidence_40_59[17+i] +incidence_40_59[18+i]
                                 +incidence_40_59[19+i] +incidence_40_59[20+i] +incidence_40_59[21+i] +incidence_40_59[22+i]) ./7 ,phi);

  H_60_79_data[i] ~  neg_binomial_2((incidence_60_79[16+i] +incidence_60_79[17+i] +incidence_60_79[18+i]
                                 +incidence_60_79[19+i] +incidence_60_79[20+i] +incidence_60_79[21+i] +incidence_60_79[22+i]) ./7 ,phi);

  H_over80_data[i] ~ neg_binomial_2((incidence_over80[16+i] +incidence_over80[17+i] +incidence_over80[18+i]
                                 +incidence_over80[19+i] +incidence_over80[20+i] +incidence_over80[21+i] +incidence_over80[22+i]) ./7 ,phi);
                                 }
}


//------------------------------------
//---Generated Quantities block-------
//------------------------------------
// this block is used to check our work, it is executed per sample and the variables are local.
// it does not modify the posterior distribution.
generated quantities {
  

  real pred_H_0_19[n_simulation];
  real pred_H_20_39[n_simulation];
  real pred_H_40_59[n_simulation];
  real pred_H_60_79[n_simulation];
  real pred_H_over80[n_simulation];

// // I would like that the following trajectories fit the data in a good way

// 1


// print("y_ode: ",y_ode);

// pred_H_0_19[1]  =  neg_binomial_2_rng((y_ode[1,3] +y_ode[2,3] +y_ode[3,3] +y_ode[4,3]) ./4 .* omega_0_19 ,phi);
// 
// pred_H_20_39[1] =  neg_binomial_2_rng((y_ode[1,6] +y_ode[2,6] +y_ode[3,6] +y_ode[4,6]) ./4 .* omega_20_39 ,phi);
//                                  
// pred_H_40_59[1] =  neg_binomial_2_rng((y_ode[1,9] +y_ode[2,9] +y_ode[3,9] +y_ode[4,9]) ./4 .* omega_40_59 ,phi);
// 
// pred_H_60_79[1] =  neg_binomial_2_rng((y_ode[1,12] +y_ode[2,12] +y_ode[3,12] +y_ode[4,12]) ./4 .* omega_60_79 ,phi);
// 
// pred_H_over80[1] = neg_binomial_2_rng((y_ode[1,15] +y_ode[2,15] +y_ode[3,15] +y_ode[4,15]) ./4 .* omega_over80  ,phi);


 //  incidence_0_19[1] = 0;
 //  incidence_20_39[1] = 0;
 //  incidence_40_59[1] = 0;
 //  incidence_60_79[1] = 0;
 //  incidence_over80[1] = 0;
 // 
 // for (i in 2:n_simulation){
 //  incidence_0_19[i] = y_ode[i,4]-y_ode[i-1,4];
 //  incidence_20_39[i] = y_ode[i,8]-y_ode[i-1,8];
 //  incidence_40_59[i] = y_ode[i,12]-y_ode[i-1,12];
 //  incidence_60_79[i] = y_ode[i,16]-y_ode[i-1,16];
 //  incidence_over80[i] = y_ode[i,20]-y_ode[i-1,20];
 //  }

pred_H_0_19[1]  =  neg_binomial_2_rng((incidence_0_19[1] +incidence_0_19[2] +incidence_0_19[3] +incidence_0_19[4]) ./4 ,phi);

pred_H_20_39[1]  =  neg_binomial_2_rng((incidence_20_39[1] +incidence_20_39[2] +incidence_20_39[3] +incidence_20_39[4]) ./4 ,phi);
                                 
pred_H_40_59[1]  =  neg_binomial_2_rng((incidence_40_59[1] +incidence_40_59[2] +incidence_40_59[3] +incidence_40_59[4]) ./4 ,phi);

pred_H_60_79[1]  =  neg_binomial_2_rng((incidence_60_79[1] +incidence_60_79[2] +incidence_60_79[3] +incidence_60_79[4]) ./4 ,phi);

pred_H_over80[1]  =  neg_binomial_2_rng((incidence_over80[1] +incidence_over80[2] +incidence_over80[3] +incidence_over80[4]) ./4 ,phi);


// 2

// pred_H_0_19[2]  =  neg_binomial_2_rng((y_ode[1,3]+y_ode[2,3]+y_ode[3,3]+y_ode[4,3] +y_ode[5,3]) ./5 .* omega_0_19 ,phi);
// 
// pred_H_20_39[2] =  neg_binomial_2_rng((y_ode[1,6]+y_ode[2,6]+y_ode[3,6]+y_ode[4,6] +y_ode[5,6]) ./5 .* omega_20_39 ,phi);
//                                  
// pred_H_40_59[2] =  neg_binomial_2_rng((y_ode[1,9]+y_ode[2,9]+y_ode[3,9]+y_ode[4,9] +y_ode[5,9]) ./5 .* omega_40_59 ,phi);
// 
// pred_H_60_79[2] =  neg_binomial_2_rng((y_ode[1,12]+y_ode[2,12]+y_ode[3,12]+y_ode[4,12] +y_ode[5,12]) ./5 .* omega_60_79 ,phi);
// 
// pred_H_over80[2] = neg_binomial_2_rng((y_ode[1,15]+y_ode[2,15]+y_ode[3,15]+y_ode[4,15] +y_ode[5,15]) ./5 .* omega_over80  ,phi);

pred_H_0_19[2]  =  neg_binomial_2_rng((incidence_0_19[1] +incidence_0_19[2] +incidence_0_19[3]
                                      +incidence_0_19[4] +incidence_0_19[5]) ./5 ,phi);

pred_H_20_39[2]  =  neg_binomial_2_rng((incidence_20_39[1] +incidence_20_39[2] +incidence_20_39[3]
                                       +incidence_20_39[4]+incidence_20_39[5]) ./5 ,phi);
                                 
pred_H_40_59[2]  =  neg_binomial_2_rng((incidence_40_59[1] +incidence_40_59[2] +incidence_40_59[3]
                                       +incidence_40_59[4]+incidence_40_59[5]) ./5 ,phi);

pred_H_60_79[2]  =  neg_binomial_2_rng((incidence_60_79[1] +incidence_60_79[2] +incidence_60_79[3]
                                       +incidence_60_79[4]+incidence_60_79[5]) ./5 ,phi);

pred_H_over80[2]  =  neg_binomial_2_rng((incidence_over80[1] +incidence_over80[2] +incidence_over80[3]
                                        +incidence_over80[4]+incidence_over80[5]) ./5 ,phi);

// 3 

// pred_H_0_19[3]  =  neg_binomial_2_rng((y_ode[1,3] +y_ode[2,3] +y_ode[3,3]
//                                  +y_ode[4,3] +y_ode[5,3] +y_ode[6,3]) ./6 .* omega_0_19 ,phi);
// 
// pred_H_20_39[3] =  neg_binomial_2_rng((y_ode[1,6] +y_ode[2,6] +y_ode[3,6]
//                                  +y_ode[4,6] +y_ode[5,6] +y_ode[6,6]) ./6 .* omega_20_39 ,phi);
//                                  
// pred_H_40_59[3] =  neg_binomial_2_rng((y_ode[1,9] +y_ode[2,9] +y_ode[3,9]
//                                  +y_ode[4,9] +y_ode[5,9] +y_ode[6,9]) ./6 .* omega_40_59 ,phi);
// 
// pred_H_60_79[3] =  neg_binomial_2_rng((y_ode[1,12] +y_ode[2,12] +y_ode[3,12]
//                                  +y_ode[4,12] +y_ode[5,12] +y_ode[6,12]) ./6 .* omega_60_79 ,phi);
// 
// pred_H_over80[3] = neg_binomial_2_rng((y_ode[1,15] +y_ode[2,15] +y_ode[3,15]
//                                  +y_ode[4,15] +y_ode[5,15] +y_ode[6,15]) ./6 .* omega_over80  ,phi);

pred_H_0_19[3]  =  neg_binomial_2_rng((incidence_0_19[1] +incidence_0_19[2] +incidence_0_19[3]
                                      +incidence_0_19[4] +incidence_0_19[5]+incidence_0_19[6]) ./6 ,phi);

pred_H_20_39[3]  =  neg_binomial_2_rng((incidence_20_39[1] +incidence_20_39[2] +incidence_20_39[3]
                                       +incidence_20_39[4]+incidence_20_39[5]+incidence_20_39[6]) ./6 ,phi);
                                 
pred_H_40_59[3]  =  neg_binomial_2_rng((incidence_40_59[1] +incidence_40_59[2] +incidence_40_59[3]
                                       +incidence_40_59[4]+incidence_40_59[5]+incidence_40_59[6]) ./6 ,phi);

pred_H_60_79[3]  =  neg_binomial_2_rng((incidence_60_79[1] +incidence_60_79[2] +incidence_60_79[3]
                                       +incidence_60_79[4]+incidence_60_79[5]+incidence_60_79[6]) ./6 ,phi);

pred_H_over80[3]  =  neg_binomial_2_rng((incidence_over80[1] +incidence_over80[2] +incidence_over80[3]
                                        +incidence_over80[4]+incidence_over80[5]+incidence_over80[6]) ./6 ,phi);




                                 
// middle

// for (i in 4:(n_simulation-3)){
//  pred_H_0_19[i]  =  neg_binomial_2_rng((y_ode[i-1,3] +y_ode[i-2,3] +y_ode[i-3,3]
//                                  +y_ode[i,3] +y_ode[i+1,3] +y_ode[i+2,3] +y_ode[i+3,3]) ./7 .* omega_0_19 ,phi);
// 
//  pred_H_20_39[i] =  neg_binomial_2_rng((y_ode[i-1,6] +y_ode[i-2,6] +y_ode[i-3,6]
//                                  +y_ode[i,6] +y_ode[i+1,6] +y_ode[i+2,6] +y_ode[i+3,6]) ./7 .* omega_20_39 ,phi);
//                                  
//   pred_H_40_59[i] =  neg_binomial_2_rng((y_ode[i-1,9] +y_ode[i-2,9] +y_ode[i-3,9]
//                                  +y_ode[i,9] +y_ode[i+1] +y_ode[i+2,9] +y_ode[i+3,9]) ./7 .* omega_40_59 ,phi);
// 
//   pred_H_60_79[i] =  neg_binomial_2_rng((y_ode[i-1,12] +y_ode[i-2,12] +y_ode[i-3,12]
//                                  +y_ode[i,12] +y_ode[i+1,12] +y_ode[i+2,12] +y_ode[i+3,12]) ./7 .* omega_60_79 ,phi);
// 
//   pred_H_over80[i] = neg_binomial_2_rng((y_ode[i-1,15] +y_ode[i-2,15] +y_ode[i-3,15]
//                                  +y_ode[i,15] +y_ode[i+1,15] +y_ode[i+2,15] +y_ode[i+3,15]) ./7 .* omega_over80  ,phi);
//                                  }

for (i in 4:(n_simulation-3)){
 pred_H_0_19[i]  =  neg_binomial_2_rng((incidence_0_19[i-1] +incidence_0_19[i-2] +incidence_0_19[i-3]
                                 +incidence_0_19[i] +incidence_0_19[i+1] +incidence_0_19[i+2] +incidence_0_19[i+3]) ./7  ,phi);

 pred_H_20_39[i] =  neg_binomial_2_rng((incidence_20_39[i-1] +incidence_20_39[i-2] +incidence_20_39[i-3]
                                 +incidence_20_39[i] +incidence_20_39[i+1] +incidence_20_39[i+2] +incidence_20_39[i+3]) ./7  ,phi);
                                 
  pred_H_40_59[i] =  neg_binomial_2_rng((incidence_40_59[i-1] +incidence_40_59[i-2] +incidence_40_59[i-3]
                                 +incidence_40_59[i] +incidence_40_59[i+1] +incidence_40_59[i+2] +incidence_40_59[i+3]) ./7  ,phi);

  pred_H_60_79[i] =  neg_binomial_2_rng((incidence_60_79[i-1] +incidence_60_79[i-2] +incidence_60_79[i-3]
                                 +incidence_60_79[i] +incidence_60_79[i+1] +incidence_60_79[i+2] +incidence_60_79[i+3]) ./7 ,phi);

  pred_H_over80[i] = neg_binomial_2_rng((incidence_over80[i-1] +incidence_over80[i-2] +incidence_over80[i-3]
                                 +incidence_over80[i] +incidence_over80[i+1] +incidence_over80[i+2] +incidence_over80[i+3]) ./7 ,phi);
                                 }
                                 
                                 
// n_simulation - 2

// pred_H_0_19[n_simulation-2]  =  neg_binomial_2_rng((y_ode[n_simulation-5,3] +y_ode[n_simulation-4,3] +y_ode[n_simulation-3,3]
//                                  +y_ode[n_simulation-2,3] +y_ode[n_simulation-1,3] +y_ode[n_simulation,3]) ./6 .* omega_0_19 ,phi);
// 
// pred_H_20_39[n_simulation-2] =  neg_binomial_2_rng((y_ode[n_simulation-5,6] +y_ode[n_simulation-4,6] +y_ode[n_simulation-3,6]
//                                  +y_ode[n_simulation-2,6] +y_ode[n_simulation-1,6] +y_ode[n_simulation,6]) ./6 .* omega_20_39 ,phi);
//                                  
// pred_H_40_59[n_simulation-2] =  neg_binomial_2_rng((y_ode[n_simulation-5,9] +y_ode[n_simulation-4,9] +y_ode[n_simulation-3,9]
//                                  +y_ode[n_simulation-2,9] +y_ode[n_simulation-1,9] +y_ode[n_simulation,9]) ./6 .* omega_40_59 ,phi);
// 
// pred_H_60_79[n_simulation-2] =  neg_binomial_2_rng((y_ode[n_simulation-5,12] +y_ode[n_simulation-4,12] +y_ode[n_simulation-3,12]
//                                  +y_ode[n_simulation-2,12] +y_ode[n_simulation-1,12] +y_ode[n_simulation,12]) ./6 .* omega_60_79 ,phi);
// 
// pred_H_over80[n_simulation-2] = neg_binomial_2_rng((y_ode[n_simulation-5,15] +y_ode[n_simulation-4,15] +y_ode[n_simulation-3,15]
//                                  +y_ode[n_simulation-2,15] +y_ode[n_simulation-1,15] +y_ode[n_simulation,15]) ./6 .* omega_over80  ,phi);
                                 
pred_H_0_19[n_simulation-2]  =  neg_binomial_2_rng((incidence_0_19[n_simulation-5] +incidence_0_19[n_simulation-4]
                                                  +incidence_0_19[n_simulation-3]  +incidence_0_19[n_simulation-2]
                                                  +incidence_0_19[n_simulation-1] +incidence_0_19[n_simulation]) ./6 ,phi);

pred_H_20_39[n_simulation-2] =  neg_binomial_2_rng((incidence_20_39[n_simulation-5] +incidence_20_39[n_simulation-4]
                                                  +incidence_20_39[n_simulation-3]  +incidence_20_39[n_simulation-2]
                                                  +incidence_20_39[n_simulation-1] +incidence_20_39[n_simulation]) ./6 ,phi);
                                 
pred_H_40_59[n_simulation-2] =  neg_binomial_2_rng((incidence_40_59[n_simulation-5] +incidence_40_59[n_simulation-4]
                                                  +incidence_40_59[n_simulation-3]  +incidence_40_59[n_simulation-2]
                                                  +incidence_40_59[n_simulation-1] +incidence_40_59[n_simulation]) ./6 ,phi);

pred_H_60_79[n_simulation-2] =  neg_binomial_2_rng((incidence_60_79[n_simulation-5] +incidence_60_79[n_simulation-4]
                                                  +incidence_60_79[n_simulation-3]  +incidence_60_79[n_simulation-2]
                                                  +incidence_60_79[n_simulation-1] +incidence_60_79[n_simulation]) ./6 ,phi);

pred_H_over80[n_simulation-2] = neg_binomial_2_rng((incidence_over80[n_simulation-5] +incidence_over80[n_simulation-4]
                                                  +incidence_over80[n_simulation-3]  +incidence_over80[n_simulation-2]
                                                  +incidence_over80[n_simulation-1] +incidence_over80[n_simulation]) ./6 ,phi);                                
                                 

// n_simulation - 1
                                 
// pred_H_0_19[n_simulation-1]  =  neg_binomial_2_rng((y_ode[n_simulation-4,3] +y_ode[n_simulation-3,3] 
//                                  +y_ode[n_simulation-2,3] +y_ode[n_simulation-1,3] +y_ode[n_simulation,3] ) ./5 .* omega_0_19 ,phi);
//                                  
// pred_H_20_39[n_simulation-1] =  neg_binomial_2_rng((y_ode[n_simulation-4,6] +y_ode[n_simulation-3,6]
//                                  +y_ode[n_simulation-2,6] +y_ode[n_simulation-1,6] +y_ode[n_simulation,6] ) ./5 .* omega_20_39 ,phi);
//                                  
// pred_H_40_59[n_simulation-1] =  neg_binomial_2_rng((y_ode[n_simulation-4,9] +y_ode[n_simulation-3,9] 
//                                  +y_ode[n_simulation-2,9] +y_ode[n_simulation-1,9] +y_ode[n_simulation,9] ) ./5 .* omega_40_59 ,phi);
// 
// pred_H_60_79[n_simulation-1] =  neg_binomial_2_rng((y_ode[n_simulation-4,12] +y_ode[n_simulation-3,12] 
//                                  +y_ode[n_simulation-2,12] +y_ode[n_simulation-1,12] +y_ode[n_simulation,12] ) ./5 .* omega_60_79 ,phi);
// 
// pred_H_over80[n_simulation-1] = neg_binomial_2_rng((y_ode[n_simulation-4,15] +y_ode[n_simulation-3,15] 
//                                  +y_ode[n_simulation-2,15] +y_ode[n_simulation-1,15] +y_ode[n_simulation,15] ) ./5 .* omega_over80  ,phi);                

pred_H_0_19[n_simulation-1]  =  neg_binomial_2_rng((incidence_0_19[n_simulation-4]
                                                  +incidence_0_19[n_simulation-3]  +incidence_0_19[n_simulation-2]
                                                  +incidence_0_19[n_simulation-1] +incidence_0_19[n_simulation]) ./5 ,phi);

pred_H_20_39[n_simulation-1] =  neg_binomial_2_rng((incidence_20_39[n_simulation-4]
                                                  +incidence_20_39[n_simulation-3]  +incidence_20_39[n_simulation-2]
                                                  +incidence_20_39[n_simulation-1] +incidence_20_39[n_simulation]) ./5 ,phi);
                                 
pred_H_40_59[n_simulation-1] =  neg_binomial_2_rng((incidence_40_59[n_simulation-4]
                                                  +incidence_40_59[n_simulation-3]  +incidence_40_59[n_simulation-2]
                                                  +incidence_40_59[n_simulation-1] +incidence_40_59[n_simulation]) ./5 ,phi);

pred_H_60_79[n_simulation-1] =  neg_binomial_2_rng((incidence_60_79[n_simulation-4]
                                                  +incidence_60_79[n_simulation-3]  +incidence_60_79[n_simulation-2]
                                                  +incidence_60_79[n_simulation-1] +incidence_60_79[n_simulation]) ./5 ,phi);

pred_H_over80[n_simulation-1] = neg_binomial_2_rng((incidence_over80[n_simulation-4]
                                                  +incidence_over80[n_simulation-3]  +incidence_over80[n_simulation-2]
                                                  +incidence_over80[n_simulation-1] +incidence_over80[n_simulation]) ./5 ,phi);

                                 
// n_simulation                                  
                                 
// pred_H_0_19[n_simulation]  =  neg_binomial_2_rng((y_ode[n_simulation-3,3] +y_ode[n_simulation-2,3] 
//                                  +y_ode[n_simulation-1,3] +y_ode[n_simulation,3] ) ./4 .* omega_0_19 ,phi);
// 
// pred_H_20_39[n_simulation] =  neg_binomial_2_rng((y_ode[n_simulation-3,6] +y_ode[n_simulation-2,6]
//                                   +y_ode[n_simulation-1,6] +y_ode[n_simulation,6] ) ./4 .* omega_20_39 ,phi);
//                                  
// pred_H_40_59[n_simulation] =  neg_binomial_2_rng((y_ode[n_simulation-3,9] +y_ode[n_simulation-2,9] 
//                                   +y_ode[n_simulation-1,9] +y_ode[n_simulation,9] ) ./4 .* omega_40_59 ,phi);
// 
// pred_H_60_79[n_simulation] =  neg_binomial_2_rng((y_ode[n_simulation-3,12] +y_ode[n_simulation-2,12] 
//                                   +y_ode[n_simulation-1,12] +y_ode[n_simulation,12] ) ./4 .* omega_60_79 ,phi);
// 
// pred_H_over80[n_simulation] = neg_binomial_2_rng((y_ode[n_simulation-3,15] +y_ode[n_simulation-2,15] 
//                                   +y_ode[n_simulation-1,15] +y_ode[n_simulation,15] ) ./4 .* omega_over80  ,phi);                                 

pred_H_0_19[n_simulation]  =  neg_binomial_2_rng((incidence_0_19[n_simulation-3]  +incidence_0_19[n_simulation-2]
                                                  +incidence_0_19[n_simulation-1] +incidence_0_19[n_simulation]) ./4 ,phi);

pred_H_20_39[n_simulation] =  neg_binomial_2_rng((incidence_20_39[n_simulation-3]  +incidence_20_39[n_simulation-2]
                                                  +incidence_20_39[n_simulation-1] +incidence_20_39[n_simulation]) ./4 ,phi);
                                 
pred_H_40_59[n_simulation] =  neg_binomial_2_rng((incidence_40_59[n_simulation-3]  +incidence_40_59[n_simulation-2]
                                                  +incidence_40_59[n_simulation-1] +incidence_40_59[n_simulation]) ./4 ,phi);

pred_H_60_79[n_simulation] =  neg_binomial_2_rng((incidence_60_79[n_simulation-3]  +incidence_60_79[n_simulation-2]
                                                  +incidence_60_79[n_simulation-1] +incidence_60_79[n_simulation]) ./4 ,phi);

pred_H_over80[n_simulation] = neg_binomial_2_rng((incidence_over80[n_simulation-3]  +incidence_over80[n_simulation-2]
                                                  +incidence_over80[n_simulation-1] +incidence_over80[n_simulation]) ./4 ,phi);
                                 
}


//-------- Considerations:
// (1) If we want to generate some trajectories using priors parameters draws,
//     we can move the ODE integration from model block (executed per leapfrog step)
//     to generated quantities (executed per sample) 
// (2) Whether priors are sensible does not affect the time it takes to do one gradient evaluation.
//     On the other hand the size of the ODE's system affectime_sequence this time. 
