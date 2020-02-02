#ifndef CTestComboNormUnbalance_h
#define CTestComboNormUnbalance_h

#include <ConvergenceTest.h>
#include <bool.h>
class EquiSolnAlgo;
class LinearSOE;


class CTestComboNormUnbalance: public ConvergenceTest
{
public:
    // constructors
    CTestComboNormUnbalance();	    	
    CTestComboNormUnbalance(double tol, int maxNumIter, int printFlag, int normType=2);

    // destructor
    ~CTestComboNormUnbalance();
    
    ConvergenceTest  *getCopy(int interations);
    
    void setTolerance(double newTol);
    int setEquiSolnAlgo(EquiSolnAlgo &theAlgo);
    
    int test(void);
    int start(void);
    
    int getNumTests(void);
    int getMaxNumTests(void);        
    double getRatioNumToMax(void);                
    const Vector &getNorms(void);    
    
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    
    
protected:
    
private:
    LinearSOE *theSOE;
    double tol;         // the tol on the norm used to test for convergence
    
    int maxNumIter;     // max number of iterations
    int currentIter;    // number of times test() has been invokes since last start()
    int printFlag;      // a flag indicating if to print on test
    int nType;          // type of norm to use (1-norm, 2-norm, p-norm, max-norm)
    
    Vector norms;       // vector to hold the norms
};

#endif
