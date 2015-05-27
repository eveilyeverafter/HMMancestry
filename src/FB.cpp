#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

/*

Function Declarations

*/

// This is an internal function for n choose k
unsigned nChoosek( unsigned n, unsigned k );

// This is the haldane function
double haldane(double d);

// This is the main function to export (haploid version)
//'@export
// [[Rcpp::export]]
DataFrame c_est_fwd_back(NumericVector snp_locations, NumericVector k0, NumericVector k1, double p_assign, double p_trans) {
//    if (!dat.inherits("data.frame")) stop("Input must be a data.frame with 5 columns\nc(\"Tetrad\", \"Spore\", \"Chr\", \"Snp\", \"p0\", \"p1\")");
    
    /*

    INITIALIZE

    */

    // n_snps is the number of snps
    int n_snps(snp_locations.size());
    //cout << "There are " << n_snps << " snps present." << endl;
    
    // displace holds the physical distance (in bps) between each snp
    NumericVector displace((n_snps-1));
    //cout << "The size of the displacement vector is: " << displace.size() << endl;
    
    // Populate the displacement vector with the distance between snps
    for (int i = 0; i <displace.size(); i++)
    {
        displace[i] = snp_locations[i+1] - snp_locations[i];
        //cout << displace[i] << "\t";
        //if((i+1) % 10 == 0)
        //{
        //    cout << endl;
        //}
    }

    // Initial probability (at first snp position)
    double pi_initial(0.5);

    /*

    SETUP & RUN

    */

    // 1) establish matrices/vectors:
        NumericMatrix emissions(n_snps, 2), forward(n_snps,2), backward(n_snps,2), posterior(n_snps,2);
        NumericVector scale(n_snps), scaleb(n_snps), states_inferred(n_snps);

    // 2) calculate emission probabilities for all states at all positions
        for(int i=0; i<n_snps; i++)
        {
            emissions(i,0) = nChoosek((k0[i]+k1[i]),k0[i])*pow(p_assign,k0[i])*pow((1 - p_assign),((k0[i]+k1[i])-k0[i]));
            emissions(i,1) = nChoosek((k0[i]+k1[i]),k1[i])*pow(p_assign,k1[i])*pow((1 - p_assign),((k0[i]+k1[i])-k1[i]));
            // cout << emissions(i,1) << "\t" << emissions(i,2) << endl;
        }
        
    // 3) calculate forward probabilities and scaling factor:
    
        // For the first position 
        for(int i=0;i<2;i++)
        {
            forward(0,i) = pi_initial*emissions(0,i);
            // cout << pi_initial << " times " << emissions(0,i) << " equals " << forward(0,i) << "\t";
        }
        // cout << endl;

        // Calculate the scale factor for each snp (used to avoid underflow)
        scale[0] = (forward(0,0) + forward(0,1));
        // cout << "The first scale factor is " << scale[0] << endl;
        
        // Rescale the first forward probabilites for the first snp
        for(int i=0;i<2;i++)
        {
            forward(0,i) = forward(0,i)/scale[0];
            // cout << "Rescaled forward is " << forward(0,i) << "\t";
        }

        // For the remainding positions 

        for(int i=1; i<n_snps; i++)
        {
            // int i=1;
            double a(1-(haldane(displace[i-1])*p_trans)), b((haldane(displace[i-1])*p_trans)), c((haldane(displace[i-1])*p_trans)), d(1-(haldane(displace[i-1])*p_trans));
            double e(emissions(i,0)), f(0.0), g(0.0), h(emissions(i,1));
           
           // Calculate forward probs
           forward(i,0) = (((e*a+f*c)*forward(i-1,0))+((e*b+f*d)*forward(i-1,1)));
           forward(i,1) = (((g*a+h*c)*forward(i-1,0))+((g*b+h*d)*forward(i-1,1)));           
            //cout << "T = " << a << " " << b << " " << c << " " << d << endl;
            //cout << "E = " << e << " " << f << " " << g << " " << h << endl;
            //cout << "F before scaling = " << forward(i,0) << " " << forward(i,1) << endl;
           // Resscale
           scale[i] = (forward(i,0) + forward(i,1));
            //cout << "The scaling is: " << scale[i] << endl;
           forward(i,0) = forward(i,0)/scale[i];
           forward(i,1) = forward(i,1)/scale[i];
            //cout << "F after scaling = " << forward(i,0) << " " << forward(i,1) << endl;
        }

    // 4) Calculate the backward probabilities and scaling factors
    
        // For the last position:
        backward((n_snps-1),0) = 1;
        backward((n_snps-1),1) = 1;

        // For the rest:
        for(int i=(n_snps-2); i>=0; i--)
        {
            double a(1-(haldane(displace[i]))*p_trans), b((haldane(displace[i]))*p_trans), c((haldane(displace[i]))*p_trans), d(1-(haldane(displace[i]))*p_trans);  
            double e(emissions(i+1,0)), f(0.0), g(0.0), h(emissions(i+1,1));
           backward(i,0) = (((a*e+b*g)*backward(i+1,0))+((a*f+b*h)*backward(i+1,1)));
           backward(i,1) = (((c*e+d*g)*backward(i+1,0))+((c*f+d*h)*backward(i+1,1)));
           scaleb[i] = (backward(i,0) + backward(i,1));
           backward(i,0) = backward(i,0)/scaleb[i];
           backward(i,1) = backward(i,1)/scaleb[i];
           
        }
        // for(int i=0; i<n_snps; i++)
        // {
        //     cout << "B after scaling = " << backward(i,0) << " " << backward(i,1) << endl;
        // }

    // 5) Calculate posteriors
        for(int i=0; i<n_snps; i++)
        {
            /*  In rare cases the expression:
             *  forward(i,0)*backward(i,0) + forward(i,1)*backward(i,1)
             *  returns 0*0=0, which introduces NaNs. Avoid that by creating
             *  a condition that checks if it exists and, if so, return a missing
             *  data flag (-99).
             */
            if((forward(i,0)*backward(i,0) + forward(i,1)*backward(i,1))==0)
            {
                posterior(i,0) = -99.0;
                posterior(i,1) = -99.0;  
            } else 
            {
                for(int j=0; j<2; j++)
                {
                    posterior(i,j) = (forward(i,j)*backward(i,j))/(forward(i,0)*backward(i,0) + forward(i,1)*backward(i,1));
                }
            }
            //cout << posterior(i,0) << " " << posterior(i,1) << endl;
        }

    // 6) total likelihood:
        double lnL(0.0);  
        for(int i=0; i<n_snps; i++)
        {
            lnL+=log(scale[i]);
        }  
        //cout << "lnL = " << lnL << endl;
        
        for(int i=0; i<n_snps; i++)
        {
            states_inferred[i] = posterior(i,0) > posterior(i,1) ? 0.0 : 1.0;
            //cout << states_inferred[i] << " ";
        }


    // Output data

    DataFrame out = DataFrame::create(
        Named("Snp")=snp_locations,
        Named("p0")=k0,
        Named("p1")=k1,
        Named("emiss")=emissions,
        Named("forward")=forward,
        Named("backward")=backward,
        Named("Fscale")=scale,
        Named("Bscale")=scaleb,
        Named("posterior")=posterior,
        Named("states_inferred")=states_inferred,
        Named("lnL")=lnL
        );
    
    return out;
}

//'@export
// [[Rcpp::export]]
DataFrame c_est_fwd_back_diploid(NumericVector snp_locations, NumericVector k0, NumericVector k1, double p_assign, double p_trans) {
//    if (!dat.inherits("data.frame")) stop("Input must be a data.frame with 5 columns\nc(\"Tetrad\", \"Spore\", \"Chr\", \"Snp\", \"p0\", \"p1\")");
    
    /*

    INITIALIZE

    */

    // n_snps is the number of snps
    int n_snps(snp_locations.size());
    //cout << "There are " << n_snps << " snps present." << endl;
    
    // displace holds the physical distance (in bps) between each snp
    NumericVector displace((n_snps-1));
    //cout << "The size of the displacement vector is: " << displace.size() << endl;
    
    // Populate the displacement vector with the distance between snps
    for (int i = 0; i <displace.size(); i++)
    {
        displace[i] = snp_locations[i+1] - snp_locations[i];
        //cout << displace[i] << "\t";
        //if((i+1) % 10 == 0)
        //{
        //    cout << endl;
        //}
    }

    // Initial probability (at first snp position)
    double pi_initial(0.3333);

    /*

    SETUP & RUN

    */

    // 1) establish matrices/vectors:
        NumericMatrix emissions(n_snps, 3), forward(n_snps,3), backward(n_snps,3), posterior(n_snps,3);
        NumericVector scale(n_snps), scaleb(n_snps), states_inferred(n_snps);

    // 2) calculate emission probabilities for all states at all positions
        for(int i=0; i<n_snps; i++)
        {
            emissions(i,0) = nChoosek((k0[i]+k1[i]),k0[i])*pow(p_assign,k0[i])*pow((1 - p_assign),((k0[i]+k1[i])-k0[i]));
            emissions(i,1) = nChoosek((k0[i]+k1[i]),k0[i])*pow(0.5, k1[i])*pow(0.5, k0[i]);
            emissions(i,2) = nChoosek((k0[i]+k1[i]),k1[i])*pow(p_assign,k1[i])*pow((1 - p_assign),((k0[i]+k1[i])-k1[i]));
            // cout << emissions(i,0) << "\t" << emissions(i,1) << "\t" << emissions(i,2) << endl;
        }
        
    // 3) calculate forward probabilities and scaling factor:
    
        // For the first position 
        for(int i=0;i<3;i++)
        {
            forward(0,i) = pi_initial*emissions(0,i);
            // cout << pi_initial << " times " << emissions(0,i) << " equals " << forward(0,i) << "\t";
        }
        // cout << endl;

        // Calculate the scale factor for each snp (used to avoid underflow)
        scale[0] = (forward(0,0) + forward(0,1) + forward(0,2));
        // cout << "The first scale factor is " << scale[0] << endl;
        
        // Rescale the first forward probabilites for the first snp
        for(int i=0;i<3;i++)
        {
            forward(0,i) = forward(0,i)/scale[0];
            // cout << "Rescaled forward is " << forward(0,i) << "\t";
        }

        // For the remainding positions 

        for(int i=1; i<n_snps; i++)
        {
            double D(displace[i-1]*p_trans); // mapping distance used for haldane function
            double R(haldane(D)); // R = recombination rate between the two snps
            // cout << D << endl;

            // double a(1-R-pow(R,2));
            // double b(R);
            // double c(pow(R,2));
            // double d(R);
            // double e(1-(2*R));
            // double f(R);
            // double g(pow(R,2));
            // double h(R);
            // double ii(1-R-pow(R,2));
            double a(pow(1-R,2));
            double b(2*R*(1-R));
            double c(pow(R,2));
            double d(R*(1-R));
            double e(pow(1-R, 2) + pow(R,2));
            double f(R*(1-R));
            double g(pow(R,2));
            double h(2*R*(1-R));
            double ii(pow(1-R,2));

            double j(emissions(i,0));
            double k(0.0);
            double l(0.0);
            double m(0.0);
            double n(emissions(i,1));
            double o(0.0);
            double p(0.0);
            double q(0.0);
            double r(emissions(i,2));


            
           // Calculate forward probs
           forward(i,0) = (forward(i-1,0)*(j*a+k*d+l*g) + forward(i-1, 1)*(j*b+k*e+l*h) + forward(i-1, 2)*(j*c+k*f+l*ii));
           forward(i,1) = (forward(i-1,0)*(m*a+n*d+o*g) + forward(i-1, 1)*(m*b+n*e+o*h) + forward(i-1, 2)*(m*c+n*f+o*ii));
           forward(i,2) = (forward(i-1,0)*(p*a+q*d+r*g) + forward(i-1, 1)*(p*b+q*e+r*h) + forward(i-1, 2)*(p*c+q*f+r*ii));          
            //cout << "T = " << a << " " << b << " " << c << " " << d << endl;
            //cout << "E = " << e << " " << f << " " << g << " " << h << endl;
            //cout << "F before scaling = " << forward(i,0) << " " << forward(i,1) << endl;
           // Resscale
           scale[i] = (forward(i,0) + forward(i,1) + forward(i,2));
            //cout << "The scaling is: " << scale[i] << endl;
           forward(i,0) = forward(i,0)/scale[i];
           forward(i,1) = forward(i,1)/scale[i];
           forward(i,2) = forward(i,2)/scale[i];           
            // cout << "F after scaling = " << forward(i,0) << " " << forward(i,1) << " " << forward(i,2) << endl;
        }

    // 4) Calculate the backward probabilities and scaling factors
    
        // For the last position:
        backward((n_snps-1),0) = 1;
        backward((n_snps-1),1) = 1;
        backward((n_snps-1),2) = 1;

        // For the rest:
        for(int i=(n_snps-2); i>=0; i--)
        {
            double D(displace[i]*p_trans); // mapping distance used for haldane function
            double R(haldane(D)); // R = recombination rate between the two snps
            // cout << D << endl;

            // double a(1-R-pow(R,2));
            // double b(R);
            // double c(pow(R,2));
            // double d(R);
            // double e(1-(2*R));
            // double f(R);
            // double g(pow(R,2));
            // double h(R);
            // double ii(1-R-pow(R,2));
            
            double a(pow(1-R,2));
            double b(2*R*(1-R));
            double c(pow(R,2));
            double d(R*(1-R));
            double e(pow(1-R, 2) + pow(R,2));
            double f(R*(1-R));
            double g(pow(R,2));
            double h(2*R*(1-R));
            double ii(pow(1-R,2));

            double j(emissions(i+1,0));
            double k(0.0);
            double l(0.0);
            double m(0.0);
            double n(emissions(i+1,1));
            double o(0.0);
            double p(0.0);
            double q(0.0);
            double r(emissions(i+1,2));
          
           // T %*% E %*% B+1
           backward(i,0) = (backward(i+1,0)*(a*j+b*m+c*p) + backward(i+1,1)*(a*k+b*n+c*q) + backward(i+1,2)*(a*l+b*o+c*r));
           backward(i,1) = (backward(i+1,0)*(d*j+e*m+f*p) + backward(i+1,1)*(d*k+e*n+f*q) + backward(i+1,2)*(d*l+e*o+f*r));
           backward(i,2) = (backward(i+1,0)*(g*j+h*m+ii*p) + backward(i+1,1)*(g*k+h*n+ii*q) + backward(i+1,2)*(g*l+h*o+ii*r));
           
           scaleb[i] = (backward(i,0) + backward(i,1) + backward(i,2));
           backward(i,0) = backward(i,0)/scaleb[i];
           backward(i,1) = backward(i,1)/scaleb[i];
           backward(i,2) = backward(i,2)/scaleb[i];           
        }
        // for(int i=0; i<n_snps; i++)
        // {
        //     cout << "B after scaling = " << backward(i,0) << " " << backward(i,1) << endl;
        // }

    // 5) Calculate posteriors
        for(int i=0; i<n_snps; i++)
        {
            /*  In rare cases the expression:
             *  forward(i,0)*backward(i,0) + forward(i,1)*backward(i,1) + forward(i,2)*backward(i,2)
             *  returns 0*0*0=0, which introduces NaNs. Avoid that by creating
             *  a condition that checks if it exists and, if so, return a missing
             *  data flag (-99).
             */
            if((forward(i,0)*backward(i,0) + forward(i,1)*backward(i,1) + forward(i,2)*backward(i,2))==0)
            {
                posterior(i,0) = -99.0;
                posterior(i,1) = -99.0;  
                posterior(i,2) = -99.0; 
            } else 
            {
                for(int j=0; j<3; j++)
                {
                    posterior(i,j) = (forward(i,j)*backward(i,j))/(forward(i,0)*backward(i,0) + forward(i,1)*backward(i,1) + forward(i,2)*backward(i,2));
                }
            }
            // cout << posterior(i,0) << " " << posterior(i,1) << " " << posterior(i,2) << endl;
        }

    // 6) total likelihood:
        double lnL(0.0);  
        for(int i=0; i<n_snps; i++)
        {
            lnL+=log(scale[i]);
            
        }  
        //cout << "lnL = " << lnL << endl;
        
        for(int i=0; i<n_snps; i++)
        {
            double maxj(0.0);
            int counter(-1);
            for(int j=0; j<3; j++)
            {
                if(posterior(i,j)>maxj){
                    maxj = posterior(i,j);
                    counter++;
                }
            }
            states_inferred[i] = counter;
            // cout << states_inferred[i] << " ";
        }

    // Output data

    DataFrame out = DataFrame::create(
        Named("Snp")=snp_locations,
        Named("p0")=k0,
        Named("p1")=k1,
        Named("emiss")=emissions,
        Named("forward")=forward,
        Named("backward")=backward,
        Named("Fscale")=scale,
        Named("Bscale")=scaleb,
        Named("posterior")=posterior,
        Named("states_inferred")=states_inferred,
        Named("lnL")=lnL
        );
    
    return out;
}

/* 

MINOR FUNCTIONS

*/

// To define n choose k
unsigned nChoosek( unsigned n, unsigned k )
{
if (k > n) return 0;
if (k * 2 > n) k = n-k;
if (k == 0) return 1;

int result = n;
for( int i = 2; i <= k; ++i ) {
    result *= (n-i+1);
    result /= i;
}
return result;
}

// To calculate the probability of an odd number of crossover events given the 
// distance (bps) between two snps
double haldane(double d)
{
    double out = (0.5)*(1.0 - exp(-2.0*d));
    return out;
}





