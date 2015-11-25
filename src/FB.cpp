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

//' @title Inferring hidden ancestry states from haploids
//'
//' @description Uses the forward-backward algorithm to estimate ancestry 
//' along a given chromosome for a given genotyped tetrad or simulated data.
//'
//' @param snp_locations a numeric vector specifying the locations of each snp (in bps). This
//' vector is assumed to be ordered (sorted from smallest to largest snp).
//' 
//' @param p0 a vector specifying the number of reads that mapped to parent 0. 
//' p0 is assumed to be in the smae order as snp_locations. 
//'
//' @param p1 a vector specifying the number of reads that mapped to parent 1. 
//' p1 is assumed to be in the smae order as snp_locations.
//'
//' @param p_assign a value specifying the assignment probabilty (see details).
//'
//' @param scale a numeric specifying the the genome-wide recombination rate (Morgans / bp). \code{scale} is assumed to be 
//' between 0 and 1 but in practice it is usually quite small. 
//'
//' @details The function \code{fb_haploid} attempts to estimate 
//' parental genotypic 'states' along a chromosome given empirical or 
//' simulated F2 cross data. 
//' Next-generation data inherits both sequencing error and missing data 
//' -- especially when sequencing coverage is low. Also, parental populations
//' can be polymorphic with respect to gene frequencies before admixture. In these cases 
//' the parental state is ambiguous or unknown. \code{fb_haploid} takes 
//' these uncertainties into account and estimates the most likely sequence of
//' ancestry. 
//' 
//' The three steps taken in \code{fb_haploid} are: 
//' \enumerate{
//'  \item Calculate forward probabilities for each state (5' to 3')
//'  \item Calculate backward probabilities for each state (3' to 5')
//'  \item From these two probabilities, calculate the posterior probability
//'  that a snp location is of a given parental state.
//' } 
//' 
//' The two element vector of forward probabilities, \eqn{f_{i}}, for each parental state 
//' at each position \eqn{ i } is calculated as:
//' \deqn{ f_{i} = \mathbf{e_{i}}~\mathbf{T_{i}}~\mathbf{f_{i-1}} }
//' 
//' where \eqn{e_{i}} is the emission probabilities for each state (see below), \eqn{T_{i}}
//' is a 2-by-2 matrix describing the transition (recombination) probabilities
//' between two states, and \eqn{f_{i-1}} is the forward probability at the 
//' previous position along the chromosome (5' of position \eqn{i}). It is assumed 
//' that each state is equally likely to occur at the first position. 
//' 
//' 
//' The emmission probabilities are calculated independently for 
//' each snp and depends on the sequence reads assigned to parent "0" 
//' and parent "1" and \code{p_assign}. The emission probabilities are
//' calculated using the binomial equation. For example, the emission
//' probability for parental state "0" at snp position i is: 
//' 
//' \deqn{e_{i}^{(0)} = {n \choose{k_{0}}} p^{k_{0}} (1-p)^{n-k_{0}}}
//' 
//' where \eqn{n} is the sum of the number of reads from both parents at snp \eqn{i} and \eqn{p} 
//' is the assignment probability, \code{p_assign}. 
//' 
//' Like the emission probabilities, the transition probabilities are also
//' calculated for each snp position. The transition probabilities in matrix \eqn{\mathbf{T_{i}}}
//' are calcuated from the genome-wide recombination rate, \code{scale}, and the physical
//' distance (in bps) between the two snps. Specifically, Haldane's function is used
//' to estimate the probabilty of getting an odd number of crossovers between snp \eqn{i}
//' and snp \eqn{i-1}. For more details, see Hether et al. (in prep).
//' 
//' To avoid underflow, we rescale the forward (and backward) probabilities
//' each iteration to that they sum to unity. 
//' 
//' The backward probabilities are calculated similarly but in the 3' to 5' 
//' direction. It is assumed that the backward probability is equally likely 
//' in each state. Following Durbin et al. (1998) the backward probability at
//' snp position \eqn{i} is: 
//' 
//' \deqn{b_{i} = \mathbf{T_{i}}~\mathbf{e_{i+1}}~\mathbf{b_{i+1}}}
//' 
//' Again, these probabilities are rescaled to avoid underflow. 
//' 
//' The posterior probability that the state at position \eqn{i} is \eqn{j}
//' given the observed sequence read counts, \eqn{x_{i}}
//' is calculated by:  
//' 
//' \deqn{P(\pi_{i} = j~|~x_{i}) = \frac{\mathbf{f}_{i}^{(j)} \mathbf{b}_{i}^{(j)}}{\mathbf{f}_{i}~\mathbf{b}_{i}}}
//' 
//' Ties in posterior probabilities are rare with sufficient coverage and/or signal. However, in the event of a tie state 1 is picked.
//'
//' Finally, we can determine the log likelihood of the whole sequence 
//' of observations by summing up the log of all the scale factors in the 
//' forward probability calculation.
//'
//' @return a dataframe with the following columns: 
//' \itemize{
//'     \item \code{snp_locations}
//'     \item \code{p0}
//'     \item \code{p1}
//'     \item The emission matrix (2 columns) specifying the emmision probability of belong
//' to the parent 0 (emiss.1) and parent 2 (emiss.2) 
//'     \item the forward probability matrix (2 columns) giving the scaled
//'     forward probabilities of belong
//' to the parent 0 (forward.1) and parent 2 (forward.2)
//'     \item the backward probability matrix (2 columns) giving the scaled
//'     backward probabilities of belong
//' to the parent 0 (backward.1) and parent 2 (backward.2)
//'     \item \code{Fscale} The forward scaling factor for each snp
//'     \item \code{Bscale} The backward scaling factor for each snp
//'     \item posterior probability matrix (2 columns) posterior
//'     probabilities for belong to each parental state
//'     \item states_inferred a vector of length snp_locations giving the state with the
//'     highest posterior. 
//'     \item lnL a numeric giving the total likelihood of the data (see details).
//' }
//'
//' @references Hether, T.D., C. G. Wiench1, and P.A. Hohenlohe (in review). 2015. Novel molecular and analytical tools 
//' for efficient estimation of rates of meiotic crossover, non-crossover and gene conversion
//'
//' Drubin, R. S. Eddy, A. Krogh, and G. Mitchison. 1998. Biological Sequence Analysis: Probabilistic Models
//' of proteins and nucleic acids. Cambridge University Press, Cambridge CB2 8RU, UK.
//'
//' @seealso \code{\link{fb_diploid}}
//'
//' @author Tyler D. Hether
//'
//' @examples
//' set.seed(1234567)        # For reproducibility
//' n_spores <- 1            # number of spores
//' l <- 75                  # number of snps to simulate
//' c <- 3.5e-05             # recombination rate between snps (Morgan/bp)
//' snps <- c(1:l)*1.3e4     # snps are evenly spaced 20kbp apart
//' p_a <- 0.95              # assignment probability
//' coverage <- 2.1          # mean coverage
//' # Now simulate
//' sim1 <- sim_en_masse(n.spores=n_spores, scale=c, snps=snps, 
//'  p.assign=p_a, mu.rate=0, f.cross=0.8, f.convert=0.3, 
//'  length.conversion=2e3, coverage=coverage)
//' fb_haploid(snp_locations=sim1$Snp, p0=sim1$p0, p1=sim1$p1, p_assign=p_a, scale=c)
//'
//' @export
// [[Rcpp::export]]
DataFrame fb_haploid(NumericVector snp_locations, NumericVector p0, NumericVector p1, double p_assign, double scale) {
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
        // cout << displace[i] << "\t";
        // if((i+1) % 10 == 0)
        // {
        //    cout << endl;
        // }
    }
    
    // Initial probability (at first snp position)
    double pi_initial(0.5);
    
    /*
     
     SETUP & RUN
     
     */
    
    // 1) establish matrices/vectors:
    NumericMatrix emissions(n_snps, 2), forward(n_snps,2), backward(n_snps,2), posterior(n_snps,2);
    NumericVector scalef(n_snps), scaleb(n_snps), states_inferred(n_snps);
    
    // 2) calculate emission probabilities for all states at all positions
    for(int i=0; i<n_snps; i++)
    {
        emissions(i,0) = nChoosek((p0[i]+p1[i]),p0[i])*pow(p_assign,p0[i])*pow((1 - p_assign),((p0[i]+p1[i])-p0[i]));
        emissions(i,1) = nChoosek((p0[i]+p1[i]),p1[i])*pow(p_assign,p1[i])*pow((1 - p_assign),((p0[i]+p1[i])-p1[i]));
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
    scalef[0] = (forward(0,0) + forward(0,1));
    // cout << "The first scale factor is " << scalef[0] << endl;
    
    // Rescale the first forward probabilites for the first snp
    for(int i=0;i<2;i++)
    {
        forward(0,i) = forward(0,i)/scalef[0];
        // cout << "Rescaled forward is " << forward(0,i) << "\t";
    }
    
    // For the remainding positions
    
    for(int i=1; i<n_snps; i++)
    {
        // int i=1;
        // std::cout << displace[i-1] << "\t";
        double a(1-(haldane(displace[i-1]*scale)*0.5));
        double b((haldane(displace[i-1]*scale)*0.5));
        double c(b);
        double d(1-(haldane(displace[i-1]*scale)*0.5));
        
        double e(emissions(i,0)), f(0.0), g(0.0), h(emissions(i,1));
        
        // Calculate forward probs
        forward(i,0) = (((e*a+f*c)*forward(i-1,0))+((e*b+f*d)*forward(i-1,1)));
        forward(i,1) = (((g*a+h*c)*forward(i-1,0))+((g*b+h*d)*forward(i-1,1)));
        //cout << "T = " << a << " " << b << " " << c << " " << d << endl;
        //cout << "E = " << e << " " << f << " " << g << " " << h << endl;
        //cout << "F before scaling = " << forward(i,0) << " " << forward(i,1) << endl;
        // Resscale
        scalef[i] = (forward(i,0) + forward(i,1));
        //cout << "The scaling is: " << scalef[i] << endl;
        forward(i,0) = forward(i,0)/scalef[i];
        forward(i,1) = forward(i,1)/scalef[i];
        //cout << "F after scaling = " << forward(i,0) << " " << forward(i,1) << endl;
    }
    
    // 4) Calculate the backward probabilities and scaling factors
    
    // For the last position:
    backward((n_snps-1),0) = 1;
    backward((n_snps-1),1) = 1;
    
    // For the rest:
    for(int i=(n_snps-2); i>=0; i--)
    {
        // std::cout << displace[i] << "\t";
        // double a(1-(haldane(displace[i]*scale))), b((haldane(displace[i]*scale)), c(b), d(1-(haldane(displace[i]*scale)));
        
        double a(1-(haldane(displace[i]*scale)*0.5));
        double b(haldane(displace[i]*scale)*0.5);
        double c(b);
        double d(1-(haldane(displace[i]*scale)*0.5));
        
        
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
        lnL+=log(scalef[i]);
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
                                      Named("p0")=p0,
                                      Named("p1")=p1,
                                      Named("emiss")=emissions,
                                      Named("forward")=forward,
                                      Named("backward")=backward,
                                      Named("Fscale")=scalef,
                                      Named("Bscale")=scaleb,
                                      Named("posterior")=posterior,
                                      Named("states_inferred")=states_inferred,
                                      Named("lnL")=lnL
                                      );
    
    return out;
}

//' @title Inferring hidden ancestry states from diploids
//'
//' @description Uses the forward-backward algorithm to estimate ancestry 
//' along a given chromosome for a given genotyped diploid.
//'
//' @param snp_locations a numeric vector specifying the locations of each snp (in bps). This
//' vector is assumed to be ordered (sorted from smallest to largest snp).
//' 
//' @param p0 a vector specifying the number of reads that mapped to parent 0. 
//' p0 is assumed to be in the smae order as snp_locations. 
//'
//' @param p1 a vector specifying the number of reads that mapped to parent 1. 
//' p1 is assumed to be in the smae order as snp_locations.
//'
//' @param p_assign a value specifying the assignment probabilty (see details).
//'
//' @param scale a numeric specifying the the genome-wide recombination rate (Morgans / bp). \code{scale} is assumed to be 
//' between 0 and 1 but in practice it is usually quite small. 
//'
//' @details This is an extension of \code{\link{fb_haploid}} with three hidden states:
//' two homozygous states and a heterozygous state. 
//'
//' @return a dataframe with the following columns: 
//' \itemize{
//'     \item \code{snp_locations}
//'     \item \code{p0}
//'     \item \code{p1}
//'     \item The emission matrix (3 columns) specifying the emmision probability of belong
//' to the only parent 0 (emiss.1), both parents (emiss.1), and only parent 2 (emiss.3) 
//'     \item the forward probability matrix (3 columns) giving the scaled
//'     forward probabilities of belong
//' to the only parent 0 (forward.1), both parents (forward.2), and only parent 2 (forward.3)
//'     \item the backward probability matrix (3 columns) giving the scaled
//'     backward probabilities of belong
//' to the only parent 0 (backward.1), both parents (backward.2), and only parent 2 (backward.3)
//'     \item \code{Fscale} The forward scaling factor for each snp
//'     \item \code{Bscale} The backward scaling factor for each snp
//'     \item posterior probability matrix (3 columns) posterior
//'     probabilities for belong to each parental state
//'     \item states_inferred a vector of length snp_locations giving the state with the
//'     highest posterior. 
//'     \item lnL a numeric giving the total likelihood of the data (see \code{\link{fb_haploid}} ).
//' }
//'
//' @references Hether, T.D., C. G. Wiench1, and P.A. Hohenlohe (in review). 2015. Novel molecular and analytical tools 
//' for efficient estimation of rates of meiotic crossover, non-crossover and gene conversion
//'
//' Drubin, R. S. Eddy, A. Krogh, and G. Mitchison. 1998. Biological Sequence Analysis: Probabilistic Models
//' of proteins and nucleic acids. Cambridge University Press, Cambridge CB2 8RU, UK.
//'
//' @seealso \code{\link{fb_haploid}}
//'
//' @author Tyler D. Hether
//'
//' @examples
//' set.seed(1234567)        # For reproducibility
//' n_spores <- 1            # number of spores
//' l <- 75                  # number of snps to simulate
//' c <- 3.5e-05             # recombination rate between snps (Morgan/bp)
//' snps <- c(1:l)*1.3e4     # snps are evenly spaced 20kbp apart
//' p_a <- 0.95              # assignment probability
//' coverage <- 2.1          # mean coverage
//' # Now simulate two haploids
//' sim1 <- sim_en_masse(n.spores=n_spores, scale=c, snps=snps, 
//'  p.assign=p_a, mu.rate=0, f.cross=0.8, f.convert=0.3, 
//'  length.conversion=2e3, coverage=coverage)
//' sim2 <- sim_en_masse(n.spores=n_spores, scale=c, snps=snps, 
//'  p.assign=p_a, mu.rate=0, f.cross=0.8, f.convert=0.3, 
//'  length.conversion=2e3, coverage=coverage)
//' # Now merge the two haploids to make a diploid
//' p0 <- sim1$p0+sim2$p0
//' p1 <- sim1$p1+sim2$p1
//' res <- fb_diploid(snp_locations=sim1$Snp, p0=p0, p1=p1, p_assign=p_a, scale=c)
//' res
//'@export
// [[Rcpp::export]]
DataFrame fb_diploid(NumericVector snp_locations, NumericVector p0, NumericVector p1, double p_assign, double scale) {
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
    NumericVector scalef(n_snps), scaleb(n_snps), states_inferred(n_snps);
    
    // 2) calculate emission probabilities for all states at all positions
    for(int i=0; i<n_snps; i++)
    {
        emissions(i,0) = nChoosek((p0[i]+p1[i]),p0[i])*pow(p_assign,p0[i])*pow((1 - p_assign),((p0[i]+p1[i])-p0[i]));
        emissions(i,1) = nChoosek((p0[i]+p1[i]),p0[i])*pow(0.5, p1[i])*pow(0.5, p0[i]);
        emissions(i,2) = nChoosek((p0[i]+p1[i]),p1[i])*pow(p_assign,p1[i])*pow((1 - p_assign),((p0[i]+p1[i])-p1[i]));
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
    scalef[0] = (forward(0,0) + forward(0,1) + forward(0,2));
    // cout << "The first scale factor is " << scalef[0] << endl;
    
    // Rescale the first forward probabilites for the first snp
    for(int i=0;i<3;i++)
    {
        forward(0,i) = forward(0,i)/scalef[0];
        // cout << "Rescaled forward is " << forward(0,i) << "\t";
    }
    
    // For the remainding positions
    
    for(int i=1; i<n_snps; i++)
    {
        double D(displace[i-1]*scale); // mapping distance used for haldane function
        double R(haldane(D)*0.5); 
        // cout << D << endl;
        
        
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
        scalef[i] = (forward(i,0) + forward(i,1) + forward(i,2));
        //cout << "The scaling is: " << scalef[i] << endl;
        forward(i,0) = forward(i,0)/scalef[i];
        forward(i,1) = forward(i,1)/scalef[i];
        forward(i,2) = forward(i,2)/scalef[i];
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
        double D(displace[i]*scale); // mapping distance used for haldane function
        // R = recombination rate between ancestry i to j (or j to i)
        // cout << D << endl;
        double R(haldane(D)*0.5); 
        
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
        lnL+=log(scalef[i]);
        
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
                                      Named("p0")=p0,
                                      Named("p1")=p1,
                                      Named("emiss")=emissions,
                                      Named("forward")=forward,
                                      Named("backward")=backward,
                                      Named("Fscale")=scalef,
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





