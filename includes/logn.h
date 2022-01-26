double logn2(double x, double eta, double sigma, double Epeak, double norm) {

  Double_t Aterm;
  Double_t logterms0,s0;
  Double_t logn,logterm;
  Double_t expterm;
  Double_t pigreco=3.14159265;
  

 ////////////////////////////// results from fit 
 //Double_t eta = -03914;
 //Double_t sigma = 0.0004587;
 //Double_t Epeak = 0.0007011;
 //Double_t norm = 0.1997;
 /////////////////////////////////////////////
  
  logterms0 = eta*2.35/2+sqrt(1+pow((eta*2.35/2),2));
  s0 = (2/2.35)*log(logterms0);

  Aterm = eta/(sqrt(2*pigreco)*sigma*s0);

  logterm = 1-(eta/sigma)*(x-Epeak);  

  
  if(logterm<0){
      logterm = 0.0001;
      }
  expterm = log(logterm)/s0;
  expterm = -0.5*pow(expterm,2);

  logn = norm*Aterm *exp(expterm);      
  return logn;

}


Double_t logn(Double_t *x, Double_t *par){
      Double_t Epeak, sigma, eta, norm;
      Double_t Aterm;
      Double_t logterms0,s0;
      Double_t logn,logterm;
      Double_t expterm;
      Double_t pigreco=3.14159265;
      
      eta = par[0];
      sigma = par[1];
      Epeak = par[2];
      norm = par[3];
      
      
     ////////////////////////////// results from fit 
     //Double_t eta = -03914;
     //Double_t sigma = 0.0004587;
     //Double_t Epeak = 0.0007011;
     //Double_t norm = 0.1997;
     /////////////////////////////////////////////
      
      logterms0 = eta*2.35/2+sqrt(1+pow((eta*2.35/2),2));
      s0 = (2/2.35)*log(logterms0);
    
      Aterm = eta/(sqrt(2*pigreco)*sigma*s0);

      logterm = 1-(eta/sigma)*(x[0]-Epeak);  

      
      if(logterm<0){
          logterm = 0.0001;
          }
      expterm = log(logterm)/s0;
      expterm = -0.5*pow(expterm,2);

      logn = norm*Aterm *exp(expterm);      
      return logn;
}  // fine function

