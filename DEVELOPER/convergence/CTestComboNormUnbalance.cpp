#include <CTestComboNormUnbalance.h>
#include <Vector.h>
#include <Channel.h>
#include <EquiSolnAlgo.h>
#include <LinearSOE.h>
//#include <iostream>



CTestComboNormUnbalance::CTestComboNormUnbalance()
	: ConvergenceTest(CONVERGENCE_TEST_CTestComboNormUnbalance),
	theSOE(0), tol(0.0), maxNumIter(0), currentIter(0), printFlag(0),
	norms(1), nType(2)
{

}


CTestComboNormUnbalance::CTestComboNormUnbalance(double theTol, int maxIter, int printIt, int normType)
	: ConvergenceTest(CONVERGENCE_TEST_CTestComboNormUnbalance),
	theSOE(0), tol(theTol), maxNumIter(maxIter), currentIter(0), printFlag(printIt),
	norms(maxNumIter), nType(normType)
{

}


CTestComboNormUnbalance::~CTestComboNormUnbalance()
{

}


ConvergenceTest* CTestComboNormUnbalance::getCopy(int iterations)
{
	CTestComboNormUnbalance *theCopy;
	theCopy = new CTestComboNormUnbalance(this->tol, iterations, this->printFlag, this->nType);

	theCopy->theSOE = this->theSOE;

	return theCopy;
}


void CTestComboNormUnbalance::setTolerance(double newTol)
{
	tol = newTol;
}


int CTestComboNormUnbalance::setEquiSolnAlgo(EquiSolnAlgo &theAlgo)
{
	theSOE = theAlgo.getLinearSOEptr();
	if (theSOE == 0) {
		opserr << "WARNING: CTestComboNormUnbalance::setEquiSolnAlgo - no SOE\n";
		return -1;
	}
	else
		return 0;
}


int CTestComboNormUnbalance::test(void)
{
	// check to ensure the SOE has been set - this should not happen if the 
	// return from start() is checked
	if (theSOE == 0)
		return -2;

	// check to ensure the algo does invoke start() - this is needed otherwise
	// may never get convergence later on in analysis!
	if (currentIter == 0) {
		opserr << "WARNING: CTestComboNormUnbalance::test() - start() was never invoked.\n";
		return -2;
	}

	// get the B vector & determine it's norm & save the value in norms vector
	const Vector &x = theSOE->getB();
	const Vector &y = theSOE->getC();

	double normX = x.pNorm(nType);
	double normY = y.pNorm(nType);

	if (currentIter <= maxNumIter)
		norms(currentIter - 1) = normX;

	// print the data if required
	if (printFlag == 1) {
		opserr << "CTestComboNormUnbalance::test() - iteration: " << currentIter;
		opserr << " current Norm: " << normX << " (max: " << tol;
		opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
	}
	if (printFlag == 4) {
		opserr << "CTestComboNormUnbalance::test() - iteration: " << currentIter;
		opserr << " current Norm: " << normX << " (max: " << tol << ")\n";
		opserr << "\tNorm deltaX: " << theSOE->getX().pNorm(nType) << ", Norm deltaR: " << normX << endln;
		opserr << "\tdeltaX: " << theSOE->getX() << "\tdeltaR: " << x;
	}

	//
	// check if the algorithm converged
	//

	// if converged - print & return ok
	if (normX <= fmax(100. * tol, tol * normY)) {

			// do some printing first
			if (printFlag != 0) {
				if (printFlag == 1 || printFlag == 4)
					opserr << endln;
				else if (printFlag == 2 || printFlag == 6) {
					opserr << "CTestComboNormUnbalance::test() - iteration: " << currentIter;
					opserr << " current Norm: " << normX << " (max: " << fmax(100.0 * tol, tol * normY);
					opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
				}
			}

			// return the number of times test has been called
			return currentIter;
	}

	// algo failed to converged after specified number of iterations - but RETURN OK
	else if ((printFlag == 5 || printFlag == 6) && currentIter >= maxNumIter) {
		opserr << "WARNING: CTestComboNormUnbalance::test() - failed to converge but going on -";
		opserr << " current Norm: " << normX << " (max: " << tol;
		opserr << ", Norm deltaX: " << theSOE->getX().pNorm(nType) << ")\n";
		return currentIter;
	}

	// algo failed to converged after specified number of iterations - return FAILURE -2
	else if (currentIter >= maxNumIter) { // the algorithm failed to converge
		opserr << "WARNING: CTestComboNormUnbalance::test() - failed to converge \n";
		opserr << "after: " << currentIter << " iterations\n";
		currentIter++;  // we increment in case analysis does not check for convergence
		return -2;
	}

	// algorithm not yet converged - increment counter and return -1
	else {
		currentIter++;
		return -1;
	}
}


int CTestComboNormUnbalance::start(void)
{
	if (theSOE == 0) {
		opserr << "WARNING: CTestComboNormUnbalance::test() - no SOE returning true\n";
		return -1;
	}

	// set iteration count = 1
	norms.Zero();
	currentIter = 1;
	return 0;
}


int CTestComboNormUnbalance::getNumTests()
{
	return currentIter;
}


int CTestComboNormUnbalance::getMaxNumTests(void)
{
	return maxNumIter;
}


double CTestComboNormUnbalance::getRatioNumToMax(void)
{
	double div = maxNumIter;
	return currentIter / div;
}


const Vector& CTestComboNormUnbalance::getNorms()
{
	return norms;
}


int CTestComboNormUnbalance::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector x(4);
	x(0) = tol;
	x(1) = maxNumIter;
	x(2) = printFlag;
	x(3) = nType;
	res = theChannel.sendVector(this->getDbTag(), cTag, x);
	if (res < 0)
		opserr << "CTestComboNormUnbalance::sendSelf() - failed to send data\n";

	return res;
}


int CTestComboNormUnbalance::recvSelf(int cTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector x(4);
	res = theChannel.recvVector(this->getDbTag(), cTag, x);

	if (res < 0) {
		opserr << "CTestComboNormUnbalance::sendSelf() - failed to send data\n";
		tol = 1.0e-8;
		maxNumIter = 25;
		printFlag = 0;
		nType = 2;
	}
	else {
		tol = x(0);
		maxNumIter = (int)x(1);
		printFlag = (int)x(2);
		nType = (int)x(3);
		norms.resize(maxNumIter);
	}
	return res;
}
