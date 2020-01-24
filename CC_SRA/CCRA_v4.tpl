

GLOBALS_SECTION

// =====================================================================

TOP_OF_MAIN_SECTION
  //gradient_structure::set_NUM_DEPENDENT_VARIABLES(20000);
  gradient_structure::set_MAX_NVAR_OFFSET(20000);
  arrmblsize = 400000000;
  //gradient_structure::set_GRADSTACK_BUFFER_SIZE(35000000); 
  //gradient_structure::set_CMPDIF_BUFFER_SIZE(35000000);   

// =====================================================================

DATA_SECTION

  !!USER_CODE ad_comm::change_datafile_name("CCRA.dat");

  init_int Nyears;
  init_int AgeMax;

  int RecDev_First;
  int RecDev_Last;
 LOCAL_CALCS  
    RecDev_First = 1 - AgeMax;
    RecDev_Last = Nyears;
 END_CALCS

  init_int F_method; //# F_method=1: Explicit F; F_method=2: Hybrid method (not implemented)
  init_vector ln_R0_prior(1,6);   //# 
  init_vector M_prior(1,6);
  init_vector h_prior(1,6); 
  init_vector S50_prior(1,6);
  init_vector Sslope_prior(1,6);
  init_vector F_t_prior(1,7);
  init_vector D_prior(1,3);
  init_vector SigmaR_prior(1,6);
  init_vector RecDev_prior(1,6);
  init_vector RecDev_biasadj(RecDev_First,RecDev_Last);
  init_vector Cw_t(1,Nyears);
  init_vector W_a(0,AgeMax);
  init_vector M_a(0,AgeMax);
  init_matrix AgeComp_at(0,AgeMax,1,Nyears);

  !! cout << "# F_method " << endl; cout << F_method << endl;
  !! cout << "# ln_R0_prior " << endl; cout << ln_R0_prior << endl;
  !! cout << "# M_prior " << endl; cout << M_prior << endl;
  !! cout << "# h_prior " << endl; cout << h_prior << endl;
  !! cout << "# S50_prior " << endl; cout << S50_prior << endl;
  !! cout << "# Sslope_prior " << endl; cout << Sslope_prior << endl;
  !! cout << "# F_t_prior " << endl; cout << F_t_prior << endl;

  !! cout << "# D_prior " << endl; cout << D_prior << endl;
  !! cout << "# SigmaR_prior " << endl; cout << SigmaR_prior << endl;
  !! cout << "# RecDev_prior " << endl; cout << RecDev_prior << endl;
  !! cout << "# RecDev_biasadj " << endl; cout << RecDev_biasadj << endl;
  !! cout << "# Cw_t " << endl; cout << Cw_t << endl;
  !! cout << "# W_a " << endl; cout << W_a << endl;
  !! cout << "# M_a " << endl; cout << M_a << endl;
  
  init_int TestVal;
  !! if (TestVal != 123456) { cout << "Test Number is not 123456" << endl; exit(1); }
  !! if(TestVal == 123456) { cout << "Data loaded" << endl; }

  int ln_R0_phase;
  number ln_R0_LO;
  number ln_R0_HI;
  int M_phase;
  number M_LO;
  number M_HI;
  int h_phase;
  number h_LO;
  number h_HI;
  int F_t_length; 
  int F_t_phase;
  number F_HI;
  int S50_phase;
  int Sslope_phase;
  int SigmaR_phase;
  number SigmaR_LO;
  number SigmaR_HI;
  int RecDev_phase;
  number RecDev_LO;
  number RecDev_HI;

 LOCAL_CALCS  
    //# ln_R0
    ln_R0_LO = ln_R0_prior(1);
    ln_R0_HI = ln_R0_prior(2);
    ln_R0_phase = ln_R0_prior(6);
    //# M
    M_LO = M_prior(1);
    M_HI = M_prior(2);
    M_phase = M_prior(6);
    //# h
    h_LO = h_prior(1);
    h_HI = h_prior(2);
    h_phase = h_prior(6);
    //# S50_phase
    S50_phase = S50_prior(6);
    //# Sslope    
    Sslope_phase = Sslope_prior(6);
    //# F_t
    F_t_phase = F_t_prior(6);
    F_HI = F_t_prior(2);
    F_t_length = F_t_prior(7);
    //# SigmaR
    SigmaR_phase = SigmaR_prior(6);
    SigmaR_LO = SigmaR_prior(1);
    SigmaR_HI = SigmaR_prior(2);
    //# RecDev
    RecDev_phase = RecDev_prior(6);
    RecDev_LO = RecDev_prior(1);
    RecDev_HI = RecDev_prior(2);
    RecDev_First = 1 - AgeMax;
    RecDev_Last = Nyears;
 END_CALCS
   
  !!CLASS ofstream post("post.txt");        //# C++ code for outputting to a file
  
// =====================================================================

INITIALIZATION_SECTION
  //# This needs to be before the PARAMETER_SECTION for ADMB-RE

// =====================================================================

PARAMETER_SECTION
  init_bounded_number ln_R0(ln_R0_LO,ln_R0_HI,ln_R0_phase);
  init_bounded_number M(M_LO,M_HI,M_phase);
  init_bounded_number h(h_LO,h_HI,h_phase);
  init_bounded_number S50(0,AgeMax,S50_phase);
  init_bounded_number Sslope(0.1,10,Sslope_phase);
  init_bounded_number SigmaR(SigmaR_LO,SigmaR_HI,SigmaR_phase);
  init_bounded_vector F_t_input(1,F_t_length,0.001,F_HI,F_t_phase);
  init_bounded_vector RecDev(RecDev_First,RecDev_Last,RecDev_LO,RecDev_HI,RecDev_phase);
  init_number Dummy(4);

  sdreport_number R0;
  number CatchCV;
  sdreport_number SB0;
  number pi;
  number Joint; 
  sdreport_number SBPR0;
  
  sdreport_vector Param_hat(1,6);
  sdreport_vector RecDev_hat(RecDev_First,RecDev_Last);
  vector F_t(1,Nyears);
  sdreport_vector ln_F_t(1,Nyears);
  vector SB_t(1,Nyears);
  sdreport_vector ln_SB_t(1,Nyears);
  vector D_t(1,Nyears);
  sdreport_vector ln_D_t(1,Nyears);
  vector Cw_t_hat(1,Nyears);
  vector S_a(0,AgeMax);
  sdreport_vector Rprop_t(1,Nyears);
  
  matrix N_at(0,AgeMax,1,Nyears);
  matrix Zn_at(0,AgeMax,1,Nyears);
  matrix Dn_at(0,AgeMax,1,Nyears);
  matrix Cn_at(0,AgeMax,1,Nyears);
  
  objective_function_value nll

// =====================================================================

PRELIMINARY_CALCS_SECTION
  ln_R0 = ln_R0_prior(3);
  if(F_method==-1 | F_method==1){ F_t_input = F_t_prior(3); }
  S50 = S50_prior(3);
  Sslope = Sslope_prior(3);
  h = h_prior(3);
  M = M_prior(3);
  SigmaR = SigmaR_prior(3);
  RecDev = RecDev_prior(3);
  
  cout << setprecision(4);
  
// =====================================================================

PROCEDURE_SECTION

  //cout << "ln_R0 = " << ln_R0 << endl;
  //cout << "F_t_input = " << F_t_input << endl;
  //cout << "S50 = " << S50 << endl;
  //cout << "Sslope = " << Sslope << endl;
  
  //# Initialize global variables
  pi = 3.141592;
  nll = 0;
  
  //# Settings
  if(F_method==-2 | F_method==-1){
    for(int YearI=1; YearI<=Nyears; YearI++){
      F_t(YearI) = F_t_input(1);
    }
  }
  if(F_method==1){
    F_t = F_t_input;
    if( current_phase()==1 ){ CatchCV = 0.1; }
    if( current_phase()==2 ){ CatchCV = 0.01; }
    if( current_phase()>=3 ){ CatchCV = 0.001; }
  }
  if(F_method==2){
    CatchCV = 100;
  }
  
  //# Transform parameters
  SBPR0 = 0;
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){ SBPR0 += exp(-M*AgeI) * W_a(AgeI) * M_a(AgeI); }
  R0 = mfexp(ln_R0);
  SB0 = R0 * SBPR0; 
  
  //# Initialization
  //# Abundance
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){
    S_a(AgeI) = 1 / (1 + exp( -Sslope * (AgeI - S50) )); 
    N_at(AgeI,1) = R0 * exp(-M * AgeI) * mfexp(RecDev(1-AgeI) - RecDev_biasadj(1-AgeI)*square(SigmaR)/2);
  }
  //# Calculate F for Pope's approximation
  if(F_method==2){ 
    F_t(1) = Cw_t(1) / (mfexp(-M/2) * sum( elem_prod(elem_prod( S_a, W_a ), extract_column(N_at,1)) )); 
    //cout <<"F_t(1) = " << F_t(1) << endl;
    Joint = 1 / (1 + exp(30*(F_t(1)-0.95)));
    //cout <<"Joint = " << Joint << endl;
    F_t(1) =  Joint*F_t(1) + (1-Joint)*0.95;
    //cout <<"F_t(1) = " << F_t(1) << endl;
  }
  //cout << "Initialization (A) done" << endl;  
  //# Deaths and removals
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){
    if(F_method==-1 | F_method==1){   //# Continuous F
      Zn_at(AgeI,1) = N_at(AgeI,1) * (1 - exp( -M - F_t(1)*S_a(AgeI) ));
      Dn_at(AgeI,1) = Zn_at(AgeI,1) * (M) / (M + F_t(1)*S_a(AgeI));
      Cn_at(AgeI,1) = Zn_at(AgeI,1) * (F_t(1)*S_a(AgeI)) / (M + F_t(1)*S_a(AgeI));
    }
    if(F_method==-2 | F_method==2){   //# Instantaneous F at mid-point of year 
      Dn_at(AgeI,1) = N_at(AgeI,1) * (1 - exp(-M/2));
      Cn_at(AgeI,1) = N_at(AgeI,1) * exp(-M/2) * S_a(AgeI)*F_t(1);
      Dn_at(AgeI,1) += N_at(AgeI,1) * exp(-M/2) * (1 - S_a(AgeI)*F_t(1)) * (1 - exp(-M/2));
      Zn_at(AgeI,1) = Cn_at(AgeI,1) + Dn_at(AgeI,1);
    }
  }
  //# Summaries
  SB_t(1) = sum( elem_prod(elem_prod(W_a,M_a), extract_column(N_at,1)) );
  Cw_t_hat(1) = sum( elem_prod( extract_column(Cn_at,1), W_a) );
  //cout << "Initialization (B) done" << endl;
  
  //# Projection
  for(int YearI=2; YearI<=Nyears; YearI++){
    //# Survival
    for(int AgeI=1; AgeI<=AgeMax; AgeI++){
      N_at(AgeI,YearI) = N_at(AgeI-1,YearI-1) - Zn_at(AgeI-1,YearI-1);
    }
    //# Spawning biomass
    SB_t(YearI) = sum( elem_prod(elem_prod(W_a,M_a), extract_column(N_at,YearI)) );
    //# Recruitment
    N_at(0,YearI) = 4 * h * R0 * SB_t(YearI) / ( SB0*(1-h) + SB_t(YearI)*(5*h-1) )  * mfexp(RecDev(YearI) - RecDev_biasadj(YearI)*square(SigmaR)/2); 
    //# Calculate F for Pope's approximation
    if(F_method==2){ 
      F_t(YearI) = Cw_t(YearI) / (mfexp(-M/2) * sum( elem_prod(elem_prod( S_a, W_a ), extract_column(N_at,YearI)) )); 
      Joint = 1 / (1 + exp(30*(F_t(YearI)-0.95)));
      F_t(YearI) =  Joint*F_t(YearI) + (1-Joint)*0.95;
    }
    //# Deaths and removals
    for(int AgeI=0; AgeI<=AgeMax; AgeI++){
      //# Removals
      if(F_method==-1 | F_method==1){   //# Continuous F
        Zn_at(AgeI,YearI) = N_at(AgeI,YearI) * (1 - exp( -M - F_t(YearI)*S_a(AgeI) ));
        Dn_at(AgeI,YearI) = Zn_at(AgeI,YearI) * (M) / (M + F_t(YearI)*S_a(AgeI));
        Cn_at(AgeI,YearI) = Zn_at(AgeI,YearI) * (F_t(YearI)*S_a(AgeI)) / (M + F_t(YearI)*S_a(AgeI));
      }
      if(F_method==-2 | F_method==2){    //# Instantaneous F at mid-point of year
        Dn_at(AgeI,YearI) = N_at(AgeI,YearI) * (1 - exp(-M/2));
        Cn_at(AgeI,YearI) = N_at(AgeI,YearI) * exp(-M/2) * S_a(AgeI)*F_t(YearI);
        Dn_at(AgeI,YearI) += N_at(AgeI,YearI) * exp(-M/2) * (1 - S_a(AgeI)*F_t(YearI)) * (1 - exp(-M/2));
        Zn_at(AgeI,YearI) = Cn_at(AgeI,YearI) + Dn_at(AgeI,YearI);
      }
    }
    //# Catch
    Cw_t_hat(YearI) = sum( elem_prod( extract_column(Cn_at,YearI), W_a) );
  }
  if(isnan(value(nll)) | isinf(value(nll)) | isinf(-1*value(nll))){ WriteCrash(); }
  //cout << "Projection done" << endl;
    
  //# Reporting 
  ln_F_t = log(F_t);
  ln_SB_t = log(SB_t);
  D_t = SB_t / (SBPR0 * R0);
  ln_D_t = log(D_t);
  Rprop_t = extract_row(N_at,0) / R0;
  Param_hat(1) = ln_R0;
  Param_hat(2) = M;
  Param_hat(3) = h;
  Param_hat(4) = S50;
  Param_hat(5) = Sslope;
  Param_hat(6) = SigmaR;
  RecDev_hat = RecDev;
  //cout << "F_t = " << F_t << endl;
  //cout << "N_at = " << N_at << endl;
  //cout << "Cw_t_hat = " << Cw_t_hat << endl;

  //# Objective function -- catches
  if(F_method==1 | F_method==2){ 
    for(int YearI=1; YearI<=Nyears; YearI++){  
      if(Cw_t(YearI)>0){ nll -= ( -log(2*pi)/2 - log(Cw_t(YearI)) - log(CatchCV)  - square(log(Cw_t(YearI))-log(Cw_t_hat(YearI)))/(2*square(CatchCV)) ); }
    }
  }
  if(isnan(value(nll)) | isinf(value(nll)) | isinf(-1*value(nll))){ WriteCrash(); }
  //cout << "Objective function -- catches done" << endl;

  //# Objective function -- compositional data
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){  
  for(int YearI=1; YearI<=Nyears; YearI++){
    //////cout << " Prob[AgeComp_at(AgeI,Year)] = " << Cn_at(AgeI,YearI)/sum(extract_column(Cn_at,YearI))*0.9999 + 0.0001/(AgeMax+1) << endl;
    nll -= ( AgeComp_at(AgeI,YearI) * log( Cn_at(AgeI,YearI)/sum(extract_column(Cn_at,YearI))*0.9999 + 0.0001/(AgeMax+1) ) );
  }}
  if(isnan(value(nll)) | isinf(value(nll)) | isinf(-1*value(nll))){ WriteCrash(); }
  //cout << "Objective function -- compositional data done" << endl;
  
  //# Objective function -- priors
    //# M
    nll -= ( -log(2*pi)/2 - log(M_prior(4)) - log(M_prior(5))  - square(log(M)-log(M_prior(4)))/(2*square(M_prior(5))) ); //# M
    //# h
    nll -= ( (h_prior(4)-1)*log((h-h_prior(1))/(h_prior(2)-h_prior(1))) + (h_prior(5)-1)*log(1-(h-h_prior(1))/(h_prior(2)-h_prior(1))) ); //# h
    //# Depletion
    if(D_prior(3)==1){ nll -= ( -log(2*pi)/2 - log(D_t(Nyears)) - log(D_prior(2))  - square(log(D_t(Nyears))-log(D_prior(1)))/(2*square(D_prior(2))) ); } //# Final depletion
    //# SigmaR
    nll -= ( -log(2*pi)/2 - log(SigmaR_prior(5))  - square(SigmaR-SigmaR_prior(4))/(2*square(SigmaR_prior(5))) );
    //# RecDevs
    for(int Index=1-AgeMax;Index<=Nyears;Index++){ nll -= ( -log(2*pi)/2 - log(SigmaR)  - square(RecDev(Index)-0)/(2*square(SigmaR)) ); }

  //# Objective function -- dummy
  nll += square(Dummy);

// ===========================================================================

FUNCTION  WriteCrash  

  post << "# Crash restart file" << endl;
  post << "# R0 " << endl; post << R0 << endl;
  post << "# F_t " << endl; post << F_t << endl;
  post << "# S50 " << endl; post << S50 << endl;
  post << "# Sslope " << endl; post << Sslope << endl;
  post << "# SB_t " << endl; post << SB_t << endl;
  post << "# Cw_t " << endl; post << Cw_t << endl;
  post << "# Cn_at " << endl; post << Cn_at << endl;

  exit(1); 
    
// =====================================================================

REPORT_SECTION
  report << "# ln_F_t_hat " << endl; report << ln_F_t << endl;
  report << "# ln_SB_t_hat " << endl; report << ln_SB_t << endl;
  report << "# ln_D_t_hat " << endl; report << ln_D_t << endl;
  report << "# RecDev_hat " << endl; report << RecDev_hat << endl;
  report << "# Param_hat " << endl; report << Param_hat << endl;

  report << "# Cw_t " << endl; report << Cw_t << endl;
  report << "# Cw_t_hat " << endl; report << Cw_t_hat << endl;
  report << "# N_at " << endl; report << N_at << endl;
  report << "# Cn_at " << endl; report << Cn_at << endl;

// =====================================================================

FINAL_SECTION

// ===========================================================================

RUNTIME_SECTION
  maximum_function_evaluations 10000
  convergence_criteria 1e-2, 1e-3, 1e-6

