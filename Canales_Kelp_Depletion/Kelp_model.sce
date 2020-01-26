clear

//-------------------------------------------------------------------
ydp=2;    //numero de años promedio representan el agotamiento percibido...
nsim=300; // numero de simulaciones...
opE=1;    //opciones de error de proceso, 1: error en el reclutamiento, 2: error en M, 3: error en produccion...


// error_std  
sigmaR=0.3;
sigmaM=0.2;
a_h=15.105; // parametros de la distr_ de h tipo beta determinados previamente IC95% 0.6-0.9
b_h=4.3464;
M=0.3;

// limites de proporción máxima de varado en los desembarques
rho1=0.1;
rho2=0.3;


//-------------------------------------------------------------------
// FUNCION OBJETIVO
//------------------------------------------------------------------------------
function [fun,B,P,alfa,betta,indice]=rutina(Xo,hv,dep,n,Captura,varado,Mv,errR,opE,sigmaR)

indice=0;

K=exp(Xo);
B(1)=K;
m=1-exp(-Mv);
alfa=4*hv*m*K/(5*hv-1);
betta=(1-hv)*K/(5*hv-1);
Capt=(1-varado).*Captura;

for i=2:n
    
    // variabilidad solo en el reclutamiento
    P(i)=max([(alfa*B(i-1)/(betta+B(i-1)))*exp(errR(i-1)+0.5*sigmaR^2)-m*B(i-1),0]); 
    

    if opE==2; // variabilidad en M
    P(i)=alfa*B(i-1)/(betta+B(i-1))-m*B(i-1)*exp(errR(i-1)+0.5*sigmaR^2);
    end


    if opE==3; // variabilidad en la fnc de producción
    P(i)=(alfa*B(i-1)/(betta+B(i-1))-m*B(i-1))*exp(errR(i-1)+0.5*sigmaR^2);
    end
    
    B(i)=B(i-1)+P(i)-Capt(i-1); 
end

if(min(m*B-varado.*Captura)<0);
    indice=-1;
end


fun=(B(n)-K*dep).^2;


endfunction
//----------------------------------------------------------------------------

    Captura=fscanfMat("desembarques.txt"); // capturas de Caldera
    Deple=fscanfMat("deple.txt"); // importa valores empiricos
   
   Ko=sum(Captura);
   largo=length(Captura);


  
for j=1:nsim // loop de simulaciones
    disp([j])
    errR(:,j)= grand(largo-1, 1, "nor", 0, sigmaR);// error proceso produccion latente en R
    errM= grand(1, 1, "nor", 0, sigmaM);// error proceso produccion latente en M. Varía entre escenarios
    hv(j)= grand(1, 1, "bet", a_h, b_h);// error proceso produccion latente en h. Varía entre escenarios.
    rho(j)= grand(1, 1, "unf", rho1, rho2);// error proceso varado anual no informativo. Varía entre escenarios.
    Valor=Deple(round(grand(1, 1, "unf", 1, length(Deple))));// muestra aleatoria empírico del agotamiento
    
// los valores de deplecion son agrupados en tres categorias 20-40; 40-60 y 60-80. Se realiza un remuestreo aleatorio proporcional a la distribucion empirica resultante y dentro del rango seleccionado, se elige un valor desde una distribución uniforme acotada entre los valores del rango. La recionalidad es que el pescador puede entregar con mayor precision, rangos redondeadosdentro de los cuales cualquier valor puede resultar igualmente probable


    if Valor<0.40 then
    deple(j)=grand(1, 1, "unf", 0.2, 0.4);
    end
    
    if Valor>=0.40 & Valor<0.60 then
    deple(j)=grand(1, 1, "unf", 0.4, 0.6);
    end
//
    if Valor>=0.60 then
    deple(j)=grand(1, 1, "unf", 0.6, 0.8);
    end
    
    Mv(j)=M*exp(errM); // valor aleatorio de M

    varado=1-(1-rho(j)).*Captura/max(Captura);
    Capt=(1-varado).*Captura;
  
         
   [log_K]=fminsearch(list(rutina,hv(j),deple(j),length(Captura),Captura,varado,Mv(j),errR(:,j), opE, sigmaR),log(Ko));
   [fun(j),Bio,P,alfa(j),betta(j),indice(j)]=rutina(log_K,hv(j),deple(j),length(Captura),Captura,varado,Mv(j),errR(:,j), opE,sigmaR);

   Pobla(:,j)=Bio;
   Produ(:,j)=P;
   B0(j)=Bio(1);
    

   Bmsy(j)=-betta(j)+sqrt(alfa(j).*betta(j)./(1-exp(-Mv(j))));
   MSY(j)=alfa(j)+(1-exp(-Mv(j))).*betta(j)-2*sqrt((1-exp(-Mv(j))).*alfa(j).*betta(j));


end // fin simulaciones


  subplot(2,2,1), histplot(10,Bmsy); title('Bmsy','FontSize',5)
  subplot(2,2,2), histplot(10,MSY); title('MSY','FontSize',5)
  subplot(2,2,3), histplot(10,deple); title('Agotamiento','FontSize',5)
  subplot(2,2,4), plot(Pobla,'g'); plot(median(Pobla,'c'),'k','Linewidth',3); title('Biomasa','FontSize',5)

  save("salidas.sod","Pobla","Captura","Mv","hv","deple","rho","Deple","B0","Bmsy","MSY");




  
