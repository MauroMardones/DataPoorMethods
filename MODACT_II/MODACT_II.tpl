DATA_SECTION

  init_vector parbiol(1,8)
  init_number nedades  
  init_vector edades(1,nedades) 
  init_number ntallas
  init_vector Tallas(1,ntallas) 
  init_vector Fr_tallas(1,ntallas)
  
  init_number A50prior // 
  init_number rangoprior // 
  init_number Fcrprior // 
  init_number Loprior //
  init_number s1prior //
  init_number s2prior // 


 // Coeficientes de variación prior
  init_number cv4 // A50
  init_number cv99 // rango selectividad
  init_number cv100 // Fcr
  init_number cv1 // Lo
  init_number cv2 // aedad
  init_number cv3 // bedad

 // Fases de estimacion
  init_int  f3 // A50
  init_int  f4 // rango selectividad
  init_int  f2 // Fcr
  init_int  f7 // Lo
  init_int  f5 // s1
  init_int  f6 // s2

  number logA50ini
  number lograngoini
  number logFcrini
  number logLoini
  number logs1ini
  number logs2ini

 !! logA50ini=log(A50prior);
 !! lograngoini=log(rangoprior);
 !! logFcrini=log(Fcrprior);
 !! logLoini=log(Loprior);
 !! logs1ini=log(s1prior);
 !! logs2ini=log(s2prior);


  init_number dts
  init_int npbr
  init_vector ratio(1,npbr)


  init_number nm // Ctallaa
  init_number  phi

  

INITIALIZATION_SECTION

  log_Fcr    logFcrini
  log_A50    logA50ini
  log_rango  lograngoini
  log_s1     logs1ini
  log_s2     logs2ini
  log_Lo     logLoini



PARAMETER_SECTION

 init_number log_Ro(1) 
 init_number log_Fcr(f2) 
 init_number log_A50(f3) 
 init_number log_rango(f4) 
 init_number log_s1(f5)
 init_number log_s2(f6)
 init_number log_Lo(f7)


 vector N0(1,nedades)
 vector N40(1,nedades)
 
 vector N(1,nedades)
 vector Sel(1,nedades)
 vector F(1,nedades)
 vector Z(1,nedades)
 vector S(1,nedades)
 vector mu_edad(1,nedades)
 vector wmed(1,nedades)
 vector msex(1,nedades)
 vector sigma_edad(1,nedades)
 vector pred_Ctot_a(1,nedades)
 vector pred_Ctot(1,ntallas)
 vector likeval(1,10)

 vector prop_obs(1,ntallas)
 vector prop_pred(1,ntallas)


 number Linf
 number k
 number Lo
 number M
 number Nest
 number Nobs
 number LMLC

 matrix Prob_talla(1,nedades,1,ntallas) 
 matrix FrecL(1,nedades,1,ntallas)

 vector Fref(1,500)
 vector YPR(1,500)
 vector BPR(1,500)
 vector Fpbr(1,npbr)
 number SPRFcr
 number SPRF01
 number F01

 number alfa
 number beta
 number Bo
 number h


 objective_function_value f



PRELIMINARY_CALCS_SECTION


 Linf=parbiol(1);
 k=parbiol(2);
 M=parbiol(3);
 h=parbiol(8);



PROCEDURE_SECTION


 Eval_prob_talla_edad();
 Eval_Dinamica_equilibrio();
 Eval_logverosim();
  

FUNCTION Eval_prob_talla_edad

  int i, j;


// genero una clave edad-talla para otros calculos. Se modela desde L(1)
 mu_edad(1)=exp(log_Lo);
 for (i=2;i<=nedades;i++)
  {
  mu_edad(i)=Linf*(1-exp(-k))+exp(-k)*mu_edad(i-1);
  }

  sigma_edad=exp(log_s1)+exp(log_s2)*mu_edad;
  
  Prob_talla = ALK( mu_edad, sigma_edad, Tallas);

  wmed=exp(parbiol(4))*pow(mu_edad,parbiol(5));


  msex=Sel=1./(1+exp(-log(19)*(mu_edad-parbiol(6))/parbiol(7)));



//----------------------------------------------------------------------
FUNCTION dvar_matrix ALK(dvar_vector& mu, dvar_vector& sig, dvector& x)
	//RETURN_ARRAYS_INCREMENT();
	int i, j;
	dvariable z1;
	dvariable z2;
	int si,ni; si=mu.indexmin(); ni=mu.indexmax();
	int sj,nj; sj=x.indexmin(); nj=x.indexmax();
	dvar_matrix pdf(si,ni,sj,nj);
	pdf.initialize();
	double xs=0.5*(x[sj+1]-x[sj]);
	for(i=si;i<=ni;i++) //loop over ages
	{
		 for(j=sj;j<=nj;j++) //loop over length bins
		{
			z1=((x(j)-xs)-mu(i))/sig(i);
			z2=((x(j)+xs)-mu(i))/sig(i);
			pdf(i,j)=cumd_norm(z2)-cumd_norm(z1);
		}//end nbins
		pdf(i)/=sum(pdf(i));
	}//end nage
	//RETURN_ARRAYS_DECREMENT();
	return(pdf);
//----------------------------------------------------------------------

FUNCTION Eval_Dinamica_equilibrio

// SELECTIVIDAD


  Sel=Sel=1./(1+exp(-log(19)*(edades-exp(log_A50))/exp(log_rango)));
  F=exp(log_Fcr)*Sel;

  Z=F+M;
  S=exp(-1.*Z);


// condición inicial Bo
  N0(1)=1.;
  for (int j=2;j<=nedades;j++)
  { N0(j)=N0(j-1)*exp(-1.*M);
    N0(nedades)=N0(nedades)/(1-exp(-1.*M));
  }
  Bo=sum(elem_prod(N0*exp(-dts*M),elem_prod(wmed,msex)));

// parametros S/R
  alfa=4*h/(5*h-1);
  beta=(1-h)/(5*h-1)*Bo;

  // se estima la sobrevivencia por edad y año
  N(1)=exp(log_Ro);
  for (int i=2;i<=nedades;i++){
  N(i)=N(i-1)*exp(-Z(i-1));
  }

  N(nedades)=N(nedades)/(1-exp(-Z(nedades)));
  pred_Ctot_a=elem_prod(elem_div(F,Z),elem_prod(1.-S,N));

  pred_Ctot=pred_Ctot_a*Prob_talla;


  prop_obs=Fr_tallas/sum(Fr_tallas);
  prop_pred=pred_Ctot/sum(pred_Ctot);
  
FUNCTION Eval_logverosim


  likeval(1)=-nm*sum(elem_prod(prop_obs,log(prop_pred)));// Fr_tallas
  likeval(2)=0.5*square((log_Lo-logLoini)/cv1);
  likeval(3)=0.5*square((log_s1-logs1ini)/cv2);
  likeval(4)=0.5*square((log_s2-logs2ini)/cv3);
  likeval(5)=0.5*square((log_A50-logA50ini)/cv4);
  likeval(6)=0.5*square((log_rango-lograngoini)/cv99);
  likeval(7)=0.5*square((log_Fcr-logFcrini)/cv100);


  likeval(8)=0.5*square((log(sum(Fr_tallas))-log(sum(pred_Ctot)))/0.05);

  f=phi*sum(likeval);




 /*

 //------------- GRABA DATOS PARA PBR-----------------------------------------------
 ofstream escribe_pbr("PBRmodel.dat");


 escribe_pbr<<"#M     steepness   dt_desove"<< endl;
 escribe_pbr<< M <<" "<< h <<" "<< dt << endl;
 escribe_pbr<<"#Numero de PBR "<<endl;
 escribe_pbr<< npbr << endl;
 escribe_pbr<<"#SSB/SSBo to evaluate"<<endl;
 escribe_pbr<<ratio<<endl;
 escribe_pbr<<"#Numero de edades "<<endl;
 escribe_pbr<<nedades<<endl;
 escribe_pbr<<"#proporcion de madurez a la edad "<<endl;
 escribe_pbr<<msex<<endl;
 escribe_pbr<<"#peso medio a la edad "<<endl;
 escribe_pbr<<wmed<<endl;
 escribe_pbr<<"#Selectividad a la edad "<<endl;
 escribe_pbr<<Sel<<endl;
 
 */ 


REPORT_SECTION

  for (int i=1;i<=nedades;i++){
    FrecL(i)=Prob_talla(i)*pred_Ctot_a(i);
  }


  report << "Modelo de Analisis Captura a la Talla MODACT_II " << endl;
  report << "PUCV 2017 (cristian.canales.r@pucv.cl) " << endl;
  report << "-----------------------------------------------------------------------------" << endl;

  report << "Frecuencias de tallas observadas y predichas" << endl;
  report <<  Tallas << endl;
  report << "---------------------------------------------------------------------------- " << endl;
  report <<  Fr_tallas << endl;
  report <<  pred_Ctot << endl;
  report << " " << endl;
  report << "Frecuencia de tallas de la captura (columnas) por grupo de edad(filas)" << endl;
  report <<  FrecL  << endl;
  report << " " << endl;
  report << "Probabilidad de la talla (columnas) a la edad (filas)" << endl;
  report <<  Prob_talla  << endl;
  report << " " << endl;
  report << " " << endl;

 //-----------------------BPR----------------------------------------

  Fref.fill_seqadd(1e-10,0.01);


  for (int j=1;j<=500;j++){ // l

    F=Fref(j)*Sel;
    Z=F+M;
    S=exp(-1.*Z);

  // se estima la sobrevivencia por edad y año
    N(1)=exp(log_Ro);
      for (int i=2;i<=nedades;i++){
    N(i)=N(i-1)*exp(-Z(i-1));
    }

  N(nedades)=N(nedades)/(1-exp(-Z(nedades)));

  BPR(j)=alfa*sum(elem_prod(elem_prod(N,exp(-dts*Z)),elem_prod(wmed,msex)))-beta;
  YPR(j)=(alfa*BPR(j)/(beta+BPR(j)))*sum(elem_prod(wmed,elem_prod(elem_div(F,Z),elem_prod(1.-S,N))));


 if(j==1){N0=N;} // guardo la comps edad virginal
 if(BPR(j)/BPR(1)>0.4){N40=N;} // guardo la comps edad 40%

   }




 for (int i=1;i<=npbr;i++){
    for (int j=1;j<=500;j++){
        if(BPR(j)/BPR(1)-ratio(i)<0){Fpbr(i)=0.5*(Fref(j)+Fref(j-1));
        j=500;}
    }
 }


    for (int j=1;j<=500;j++){
        if(Fref(j)>exp(log_Fcr)){SPRFcr=BPR(j)/BPR(1);
        j=500;}
    }


// rutina para el F0.1 aprox
    for (int j=1;j<=500;j++){
    
        if(YPR(j)>0.1*(YPR(2)-YPR(1))/(Fref(2)-Fref(1)))
        {F01=0.5*(Fref(j)+Fref(j-1));
        j=500;}
     }


    for (int j=1;j<=500;j++){
        if(Fref(j)>F01){SPRF01=BPR(j)/BPR(1);
        j=500;}
    }


  report << "Edad  Talla  Dev.st  Npobl  Capt  Mad  Peso  Sel  F " << endl;
  report << "----------------------------------------------------- " << endl;
  for (int j=1;j<=nedades;j++){ // l
  report << edades(j) <<" "<<mu_edad(j)<<" "<<sigma_edad(j)<<" "<<N(j)<<" "<<pred_Ctot_a(j)<<" "<<msex(j)<<" "<<wmed(j)<<" "<<Sel(j)<<" "<<Sel(j)*exp(log_Fcr)<<endl; 
  }


 report << " "<<endl;
 report << " "<<endl;

 report << "Comps de edades de la población explotable virginal y al 40%Bo" << endl;
 report << "----------------------------------------------------- " <<endl;
 report << elem_prod(N0,Sel)<<endl;
 report << elem_prod(N40,Sel)<<endl;
 report << " "<<endl;
 report << "Comps tallas de la población explotable virginal y al 40%B0" << endl;
 report << "----------------------------------------------------- " <<endl;
 report << elem_prod(N40,Sel)*Prob_talla<<endl;
 report << elem_prod(N0,Sel)*Prob_talla<<endl;

 report << " "<<endl;
 report << " "<<endl;
 
 report << "Parámetros del modelo " << endl;
 report << "----------------------------------------------------- " <<endl;
 report<<"Ro	Fcr	A50	rango	s1_edad	  s2_edad     Lo"<<endl;
 report<<exp(log_Ro)<<" "<<exp(log_Fcr)<<" "<<exp(log_A50)<<" "<<exp(log_rango)<<" "<<exp(log_s1)<<" "<<exp(log_s2)<<" "<<exp(log_Lo)<<endl;

 report << " "<<endl;
 report << " "<<endl;


 report << "Puntos Biológicos de Referencia - Fpbr" << endl;
 report << "----------------------------------------------------- " <<endl;
 report << "% Biomasa desovante virginal (%Bo)"<<endl;
 report << ratio<<endl;
 report << " "<<endl;
 report << "Mortalidad por pesca para alcanzar %Bo"<<endl;
 report << Fpbr<<endl;
 report << " "<<endl;
 report << "Mortalidad por pesca Fo.1 y %Bo(F=Fo.1)"<<endl;
 report << F01<<" "<<SPRF01<<endl;
 report << " "<<endl;
 report << "Mortalidad por pesca actual Fcr y %Bo(F=Fcr)"<<endl;
 report << exp(log_Fcr)<<" "<<SPRFcr<<endl;
 report << " "<<endl;


 report << "Curvas de BPR y YPR" << endl;
 report << "F      BPR      YPR      %B0"<<endl;
 report << "----------------------------------------------------- " << endl;
 for (int i=1;i<=500;i++){ // l
 report << Fref(i)<<" "<<BPR(i)<<" "<<YPR(i)<<" "<<BPR(i)/BPR(1)<<endl;
 }


 
  report << " " << endl;
  report << "Componentes de log-verosimilitud (verosimilitud + priors)" << endl;
  report << "  " << endl;
  report << "pC         Lo        s1        s2        A50        rango          Fcr         Ntot" << endl;
  report << likeval<<endl;
  report << "Total  " << endl;
  report << sum(likeval)<<endl;



