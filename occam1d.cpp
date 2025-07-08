#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <complex>

#include <Eigen/Dense>
#include "/home/m104x/programs/matplotlib-cpp/matplotlibcpp.h"

namespace plt = matplotlibcpp;


const double PI = 3.1415926535897932384626433832795028841971;
const double MU = 4 * PI * 1e-7;  // Magnetic permeability [H/m]
const std::complex<double> imaginaryUnit(0.0, 1.0);

std::string dataFile;
int numStation, stationNumber, numData;
double distance;

int maxIteration = 100;
double differentialH = 1e-3;
double initialDelta = 1;
double deltaThreshold = 1;
double ndiv, vStart, vEnd, initialResistivity, vGamma, hGamma;
double lambda = 2e-4;

double factor;
int numParameter;
int sumNumData(0);

// Observed data
std::vector<int> numDataVec;
std::vector<double> distanceVec;
std::vector<std::vector<double>> freqVec;
std::vector<std::vector<double>> appResWeightVec;
std::vector<std::vector<double>> appResVec;
std::vector<std::vector<double>> phaseWeightVec;
std::vector<std::vector<double>> phaseVec;

// Model data
std::vector<double> depthVector;
std::vector<double> thicknessVector;

// Calculated data
std::vector<double> resistivityVector;
std::vector<std::vector<double>> resistivityVectorAll;
std::vector<double> apparentResistivityVector;
std::vector<double> phaseVector;
std::vector<std::complex<double>> impedanceVector;

std::vector<double> rmsVec;

double delta(1);

double xlimMin, xlimMax, ylimMin, ylimMax;


double res2logRes(double resistivity) {
	return std::log10(resistivity);
}


double logRes2res(double logResistivity) {
	return std::pow(10, logResistivity);
}


double deg2rad(double degree) {
	return degree * PI / 180;
}


double rad2deg(double radian) {
	return radian * 180 / PI;
}


std::vector<double> unwrapPhase(std::vector<double> phaseVector) {
	if (phaseVector.empty()) {
		return {};
	}

	std::vector<double> unwrappedPhase = phaseVector;

	for (size_t i = 1 ; i < unwrappedPhase.size() ; i++) {
		double diff = unwrappedPhase[i] - unwrappedPhase[i - 1];

		if (diff > 0.5 * PI && diff < 0.8 * 2 * PI) {
			unwrappedPhase[i] -= PI;
		} else if (diff < -0.5 * PI && diff > -0.8 * 2 * PI) {
			unwrappedPhase[i] += PI;
		}
	}

	return unwrappedPhase;
}


std::tuple<std::string, double, double, double, double, double, double, double, double, double, double> readParameter(std::string paramFile) {

	std::cout << "Read parameters from '" << paramFile << "' ." << std::endl;

	std::ifstream ifs(paramFile.c_str(), std::ios::in);
	ifs >> dataFile >> maxIteration >> differentialH >> deltaThreshold >> ndiv >> vStart >> vEnd >> initialResistivity >> vGamma >> hGamma >> lambda;;

	std::cout << " Filename            : " << dataFile << std::endl << " Max Iteration       : " << maxIteration << std::endl << " Differential h      : " << differentialH << std::endl << " Delta Threshold     : " << deltaThreshold << std::endl << " ndiv                : " << ndiv << std::endl << " Range of v          : " << vStart << " - " << vEnd << std::endl << " Initial Resistivity : " << initialResistivity << std::endl << " Vertical Gamma      : " << vGamma << std::endl << " Horizontal Gamma    : " << hGamma << std::endl << " Initial Lambda          : " << lambda << std::endl << std::endl;

	return std::make_tuple(dataFile, maxIteration, differentialH, deltaThreshold, ndiv, vStart, vEnd, initialResistivity, vGamma, hGamma, lambda);

}


std::tuple<int, std::vector<int>, std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> readData (std::string dataFile) {

	std::cout << "Read data from '" << dataFile << "' ." << std::endl;

	std::ifstream ifs(dataFile.c_str(), std::ios::in);
	ifs >> numStation;
	std::cout << " Number of Stations : " << numStation << std::endl;

	numDataVec     .resize(numStation);
	distanceVec    .resize(numStation);
	freqVec        .resize(numStation);
	appResWeightVec.resize(numStation);
	appResVec      .resize(numStation);
	phaseWeightVec .resize(numStation);
	phaseVec       .resize(numStation);

	for (int i = 0 ; i < numStation ; i++) {

		ifs >> stationNumber >> numData >> distance;

		numDataVec[i]  = numData;
		freqVec[i]        .resize(numData);
		appResWeightVec[i].resize(numData);
		appResVec[i]      .resize(numData);
		phaseWeightVec[i] .resize(numData);
		phaseVec[i]       .resize(numData);

		if (i == 0) {
			distanceVec[i] = distance;
		} else {
			distanceVec[i] += distance;
		}

		std::cout << " Read station : " << stationNumber << std::endl;

		for (int j = 0 ; j < numData ; j++) {

			double freq, appResWeight, appRes, phaseWeight, phase;
			ifs >> freq >> appRes >> appResWeight >> phase >> phaseWeight;
			freqVec[i][j]         = freq;
			appResVec[i][j]       = res2logRes(appRes);
			appResWeightVec[i][j] = appResWeight;
			phaseVec[i][j]        = deg2rad(phase);
			phaseWeightVec[i][j]  = phaseWeight;

		}

		phaseVec[i] = unwrapPhase(phaseVec[i]);

	}

	std::cout << std::endl;

	for (int i = 0 ; i < numStation ; i++) {

		std::cout << " Station : " << i+1 << std::endl;
		std::cout << "  Freq                App.Res     Weight(App.Res)              Phase       Weight(Phase)" << std::endl;

		for (int j = 0 ; j < numDataVec[i] ; j++) {

			std::cout << "  " << std::setw(20) << std::left << freqVec[i][j] << std::setw(20) << std::left << appResVec[i][j] << std::setw(20) << std::left << appResWeightVec[i][j] << std::setw(20) << std::left << phaseVec[i][j] << std::setw(20) << std::left << phaseWeightVec[i][j] << std::endl; 

		}

		std::cout << std::endl;

	}

	return std::make_tuple(numStation, numDataVec, distanceVec, freqVec, appResWeightVec, appResVec, phaseWeightVec, phaseVec);

}


std::tuple<std::vector<double>, std::vector<double>> forwardCalc(std::vector<double> freqVec, std::vector<double> resistivityVector, std::vector<double> thicknessVector) {

//	std::cout << "Forward calculation" << std::endl;

	int numLayer = resistivityVector.size();

	double resistivityBuf, thicknessBuf;
	double apparentResistivityBuf, phaseBuf;
	
	apparentResistivityVector.resize(freqVec.size());
	phaseVector              .resize(freqVec.size());
	impedanceVector          .resize(numLayer);

	for (int i = 0 ; i < static_cast<int>(freqVec.size()) ; i++) {

//		std::cout << " FREQUENCY = " << freqVec[i] << std::endl;
		double OMEGA = 2 * PI * freqVec[i];

		impedanceVector[numLayer - 1] = std::sqrt(OMEGA * MU * resistivityVector[numLayer - 1] * imaginaryUnit);

		for (int j = numLayer - 2 ; j > -1 ; j--) {

			resistivityBuf = resistivityVector[j];
			thicknessBuf   = thicknessVector[j];

			std::complex<double> dj, wj, ej, rj, re, zj, belowImpedance;

			dj = std::sqrt((OMEGA * MU * (1.0 / resistivityBuf)) * imaginaryUnit);
			wj = dj * resistivityBuf;
			ej = std::exp(-2 * thicknessBuf * dj);

			belowImpedance = impedanceVector[j + 1];

			rj = (wj - belowImpedance) / (wj + belowImpedance);
			re = rj * ej;
			zj = wj * ((std::complex(1.0, 0.0) - re) / (std::complex(1.0, 0.0) + re));
			impedanceVector[j] = zj;

		}

//		std::cout << "impedanceVector[0] = " << impedanceVector[0] << std::endl;
//		std::cout << "abs(impedanceVector[0]) = " << std::abs(impedanceVector[0]) << std::endl;
//		std::cout << "abs(impedanceVector[0])^2 = " << std::abs(impedanceVector[0]) * std::abs(impedanceVector[0]) << std::endl;
//		std::cout << "OMEGA * MU = " << OMEGA * MU << std::endl;
//		std::cout << "(abs(impedanceVector[0])^2) / (OMEGA * MU) = " << (std::abs(impedanceVector[0]) * std::abs(impedanceVector[0])) / (OMEGA * MU) << std::endl;

		apparentResistivityBuf = std::log10((std::abs(impedanceVector[0]) * std::abs(impedanceVector[0])) / (OMEGA * MU));
		phaseBuf               = std::atan2(impedanceVector[0].imag(), impedanceVector[0].real());

		apparentResistivityVector[i] = apparentResistivityBuf;
		phaseVector[i]               = phaseBuf;

//		std::cout << "apparentResistivityVector[" << i << "] = " << apparentResistivityVector[i] << std::endl;
//		std::cout << "phaseVector[" << i << "] = " << phaseVector[i] << std::endl;
		
//		std::cout << std::endl;

	}

	return std::make_tuple(apparentResistivityVector, phaseVector);

}


void drawSoundingCurve(std::string dataType, std::string dataUnit, std::vector<double> freq, std::vector<double> calculatedData, std::vector<double> observedData, double xlimMin, double xlimMax, double ylimMin, double ylimMax, int stationNumber, std::string outputDirectoryName){

	std::string figTitle = outputDirectoryName + dataType + "_" + std::to_string(stationNumber+1)+".png";
	std::cout << "  Draw " << figTitle << std::endl;
	std::string calculatedDataLabel = dataType+"-Cal";
	std::string observedDataLabel = dataType+"-Obs";
	std::string optionX = "bo-";
	std::string optionY = "ro";

	if( dataType == "AppResis" ){
		plt::named_loglog(calculatedDataLabel, freq, calculatedData, optionX);
		plt::named_loglog(observedDataLabel, freq, observedData, optionY);
	}else if( dataType == "Phase" ){
		plt::named_semilogx(calculatedDataLabel, freq, calculatedData, optionX);
		plt::named_semilogx(observedDataLabel, freq, observedData, optionY);
	}

	plt::xlim(xlimMin, xlimMax);
	plt::ylim(ylimMin, ylimMax);
	plt::xlabel("Frequency[Hz]");
	plt::ylabel(dataType+"["+dataUnit+"]");
	plt::legend();
	plt::title("ST"+std::to_string(stationNumber+1)+" "+dataType);
	plt::save(figTitle);
	plt::clf();
	plt::close();

}

 
void draw2DSection(std::vector<std::vector<double>> resistivityVector, std::vector<double> depthVector, std::vector<double> distanceVector, int numStation, int numParameter, std::string outputDirectoryName) {

	std::string outfileTitle = "outputDataForGnuplot.txt";
	std::ofstream outfile(outfileTitle);

	if (!outfile) {
		std::cerr << "Error: Cannot open " << outfileTitle << " for writing." << std::endl;
		return;
	}

//	outfile << "# Station_ID  Depth(km)  log10_Resistivity" << std::endl;
//
//	for (int i = 0 ; i < numStation ; i++) {
//		for (int j = 0 ; j < numParameter ; j++) {
//			outfile << (i + 1) << " " << depthVector[j] << " " << res2logRes(resistivityVector[i][j]) << std::endl;
//		}
//	}

	outfile << "# Distance(km)  Depth(km)  log10_Resistivity" << std::endl;

	for (int i = 0 ; i < numStation ; i++) {
		for (int j = 0 ; j < numParameter ; j++) {
			outfile << distanceVector[i] << " " << depthVector[j] << " " << res2logRes(resistivityVector[i][j]) << std::endl;
		}
	}

	outfile.close();

 }


int main () {

	// Read parameters
	auto [dataFile, maxIteration, differentialH, deltaThreshold, ndiv, vStart, vEnd, initialResistivity, vGamma, hGamma, lambda] = readParameter("param.dat");

	// Read data
	auto [numStation, numDataVec, distanceVec, freqVec, observedAppResWeightVec, observedAppResVec, observedPhaseWeightVec, observedPhaseVec] = readData(dataFile);

	// Set initial parameters
	std::cout << "Initial parameters" << std::endl;
	factor = std::pow(10, 1/ndiv);
	std::cout << " factor = " << factor << std::endl;
	numParameter = (std::log10(vEnd) - std::log10(vStart)) * ndiv + 3;
	std::cout << " numParameter = " << numParameter << std::endl;
	
	depthVector         .resize(numParameter+1);
	thicknessVector     .resize(numParameter);
	resistivityVector   .resize(numParameter);
	resistivityVectorAll.resize(numStation);

	depthVector[0]     = vStart;
	thicknessVector[0] = vStart;

	for (int i = 1 ; i < numParameter ; i++) {
		depthVector[i]      = depthVector[i-1] * factor;
		thicknessVector[i]  = depthVector[i] - depthVector[i-1];
	}

	for (int i = 0 ; i < numParameter ; i++) {
		resistivityVector[i] = initialResistivity;
	}

	for (int i = 0 ; i < numStation ; i++) {
		resistivityVectorAll[i] = resistivityVector;
	}

	depthVector[numParameter] = depthVector[numParameter-1] * 1.1;

	// Calculate data
	std::vector<std::vector<double>> apparentResistivityVector(numStation), phaseVector(numStation);
	std::vector<double> apparentResistivityVectorBufVector, phaseVectorBufVector;
	std::vector<std::vector<double>> residualApparentResistivityVector(numStation), residualPhaseVector(numStation);

	double objectiveFunction(0);
//	double residualBuf(0);

	for (int i = 0 ; i < numStation ; i++) {

		std::cout << " Station No. " << i+1 << std::endl;

		auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resistivityVector, thicknessVector);

//		phaseVectorBufVector = unwrapPhase(phaseVectorBufVector);

		sumNumData += numDataVec[i];
		std::cout << "  sumNumData = " << sumNumData << std::endl;
		
		std::vector<double> residualApparentResistivityVectorBufVector(numDataVec[i]), residualPhaseVectorBufVector(numDataVec[i]);

		for (int j = 0 ; j < numDataVec[i] ; j++) {
			
			std::cout << "  Data No. " << j << std::endl;
//			std::cout << "   observedAppResVec[" << i << "][" << j << "] = " << observedAppResVec[i][j] << std::endl;
//			std::cout << "   observedAppResWeightVec[" << i << "][" << j << "] = " << observedAppResWeightVec[i][j] << std::endl;
//			std::cout << "   apparentResistivityVectorBufVector[" << j << "] = " << apparentResistivityVectorBufVector[j] << std::endl;
			
			residualApparentResistivityVectorBufVector[j] = (observedAppResVec[i][j] - apparentResistivityVectorBufVector[j]) * std::sqrt(observedAppResWeightVec[i][j]);
//			std::cout << "   residualApparentResistivityVectorBufVector[" << j << "] = " << residualApparentResistivityVectorBufVector[j] << std::endl;
			residualPhaseVectorBufVector[j]               = (observedPhaseVec[i][j]  - phaseVectorBufVector[j]              ) * std::sqrt(observedPhaseWeightVec[i][j]) ;
//			std::cout << "   residualPhaseVectorBufVector[" << j << "] = " << residualPhaseVectorBufVector[j] << std::endl;
//			residualBuf += std::pow(residualApparentResistivityVectorBufVector[j], 2) + std::pow(residualPhaseVectorBufVector[j], 2);
//			std::cout << "   residualBuf = " << residualBuf << std::endl;
//			objectiveFunction += std::pow(residualApparentResistivityVectorBufVector[j], 2) + residualApparentResistivityVectorBufVector[j] * residualPhaseVectorBufVector[j];
			objectiveFunction += std::pow(residualApparentResistivityVectorBufVector[j], 2) + std::pow(residualPhaseVectorBufVector[j], 2);
//			std::cout << "   objectiveFunction = " << objectiveFunction << std::endl;
		}

		apparentResistivityVector[i]         = apparentResistivityVectorBufVector;
		phaseVector[i]                       = phaseVectorBufVector;
		residualApparentResistivityVector[i] = residualApparentResistivityVectorBufVector;
		residualPhaseVector[i]               = residualPhaseVectorBufVector;

	}

//	rmsVec.push_back(std::sqrt(residualBuf / (2 * sumNumData)));
	rmsVec.push_back(std::sqrt(objectiveFunction / (2 * sumNumData)));

	std::cout << " Initial RMS                = " << rmsVec[0]         << std::endl;
	std::cout << " Initial Objective Function = " << objectiveFunction << std::endl;

	std::cout << std::endl;

	int successfulIter(0);
	while (successfulIter < maxIteration) {
		std::cout << "\n<-- Iteration: " << successfulIter + 1 << " -->" << std::endl;
		std::cout << "  Current lambda: " << lambda << std::endl;

		// Calculate Jacobian
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(sumNumData * 2, numParameter * numStation);

		int numDataBuf(0), numParameterBuf(0);

		for (int i = 0 ; i < numStation ; i++) {
			
			for (int j = 0 ; j < numParameter ; j++) {

				double resistivityBuf = resistivityVectorAll[i][j];
				resistivityVectorAll[i][j] = resistivityBuf + differentialH;
				auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resistivityVectorAll[i], thicknessVector);

				for (int k = 0 ; k < numDataVec[i] ; k++) {

					A(k + numDataBuf                , j + numParameterBuf) = (apparentResistivityVectorBufVector[k] - apparentResistivityVector[i][k]) / differentialH * std::sqrt(observedAppResWeightVec[i][k]);
					A(k + numDataBuf + numDataVec[i], j + numParameterBuf) = (phaseVectorBufVector[k]               - phaseVector[i][k]              ) / differentialH * std::sqrt(observedPhaseWeightVec[i][k]);

				}

				resistivityVectorAll[i][j] = resistivityBuf;

			}

			numDataBuf += 2 * numDataVec[i];
			numParameterBuf += numParameter;

		}

		// Construct Roughness Matrix C
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(numParameter * numStation, numParameter * numStation);

		numDataBuf = 0;
		numParameterBuf = 0;

		for (int i = 0 ; i < numStation ; i++) {

			for (int j = 0 ; j < numParameter ; j++) {

				// R1 (Vertical Smoothing)
				if (j != 0) {
					C(j + numParameterBuf, j + numParameterBuf    ) +=      std::sqrt(vGamma);
					C(j + numParameterBuf, j + numParameterBuf - 1) += -1 * std::sqrt(vGamma);
				}

				// R2 (Horizontal Smoothing)
				if (i == 0) {
					C(j + numParameterBuf, j + numParameterBuf               ) +=      std::sqrt(hGamma);
					C(j + numParameterBuf, j + numParameterBuf + numParameter) += -1 * std::sqrt(hGamma);
				} else if (i == numStation - 1) {
					C(j + numParameterBuf, j + numParameterBuf               ) +=      std::sqrt(hGamma);
					C(j + numParameterBuf, j + numParameterBuf - numParameter) += -1 * std::sqrt(hGamma);
				} else {
					C(j + numParameterBuf, j + numParameterBuf               ) +=  2 * std::sqrt(hGamma);
					C(j + numParameterBuf, j + numParameterBuf - numParameter) += -1 * std::sqrt(hGamma);
					C(j + numParameterBuf, j + numParameterBuf + numParameter) += -1 * std::sqrt(hGamma);
				}

			}

			numParameterBuf += numParameter;

		}

		numParameterBuf = 0;

		// Construct Residual Vector B
		Eigen::MatrixXd B = Eigen::VectorXd::Zero(sumNumData * 2);

		std::vector<double> numDataVecBuf(numDataVec.begin(), numDataVec.end());
		numDataVecBuf.insert(numDataVecBuf.begin(), 0.0);

		for (int i = 0 ; i < numStation ; i++) {

			numDataBuf += 2 * numDataVecBuf[i];

			for (int j = 0 ; j < numDataVec[i] ; j++) {

//				std::cout << "i = " << i << "  j = " << j << std::endl;
				B(j + numDataBuf                ) = residualApparentResistivityVector[i][j];
				B(j + numDataBuf + numDataVec[i]) = residualPhaseVector[i][j];

			}

		}

		numDataBuf = 0;

		// Nonlilnlear Least Squares Method
		Eigen::MatrixXd C5 = lambda * C;
		Eigen::MatrixXd AandC5(A.rows() + C5.rows(), A.cols());
		AandC5 << A, C5;

		Eigen::VectorXd H = Eigen::VectorXd::Zero(numParameter * numStation);
		Eigen::VectorXd BandH(B.size() + H.size());
		BandH << B, H;

		// Cholesky分解
//		Eigen::MatrixXd AtA = AandC5.transpose() * AandC5;
//		Eigen::VectorXd Atb = AandC5.transpose() * BandH;
//		Eigen::VectorXd dxCholesky = AtA.ldlt().solve(Atb);

		// QR分解
//		Eigen::VectorXd dxQr = AandC5.colPivHouseholderQr().solve(BandH);

		// SVD分解
		Eigen::VectorXd dxSvd = AandC5.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(BandH);

		// Calculate new data
		std::vector<std::vector<double>> resistivityVectorAllBuf(numStation);

		for (int i = 0 ; i < numStation ; i++) {

			auto segment = dxSvd.segment(i * numParameter, numParameter);
			std::vector<double> updatedResistivity(numParameter);
			
			for (int j = 0 ; j < numParameter ; j++) {

				updatedResistivity[j] = resistivityVectorAll[i][j] + segment(j);

			}

			resistivityVectorAllBuf[i] = updatedResistivity;

		}

		// Calculate new objective function
		std::vector<std::vector<double>> newApparentResistivityVector(numStation), newPhaseVector(numStation), newResidualApparentResistivityVector(numStation), newResidualPhaseVector(numStation);
		double beforeObjectiveFunction = objectiveFunction;
		objectiveFunction = 0;

		for (int i = 0 ; i < numStation ; i++) {

			auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resistivityVectorAllBuf[i], thicknessVector);

//			phaseVectorBufVector = unwrapPhase(phaseVectorBufVector);

			std::vector<double> residualApparentResistivityVectorBufVector(numDataVec[i]), residualPhaseVectorBufVector(numDataVec[i]);

			for (int j = 0 ; j < numDataVec[i] ; j++) {

				residualApparentResistivityVectorBufVector[j] = (observedAppResVec[i][j] - apparentResistivityVectorBufVector[j]) * std::sqrt(observedAppResWeightVec[i][j]);
				residualPhaseVectorBufVector[j]               = (observedPhaseVec[i][j]  - phaseVectorBufVector[j]              ) * std::sqrt(observedPhaseWeightVec[i][j]) ;
//				residualBuf += residualApparentResistivityVectorBufVector[j] + residualPhaseVectorBufVector[j];
				objectiveFunction += std::pow(residualApparentResistivityVectorBufVector[j], 2) + std::pow(residualPhaseVectorBufVector[j], 2);

			}

			newApparentResistivityVector[i]         = apparentResistivityVectorBufVector;
			newPhaseVector[i]                       = phaseVectorBufVector;
			newResidualApparentResistivityVector[i] = residualApparentResistivityVectorBufVector;
			newResidualPhaseVector[i]               = residualPhaseVectorBufVector;

		}

		// Check convergence
		if (objectiveFunction < beforeObjectiveFunction) {

			std::cout << "  STATUS: SUCCESS - Objective function improved - " << beforeObjectiveFunction << " --> " << objectiveFunction << std::endl;
			successfulIter++;
			lambda *= 0.5; // Decrease lambda

			// Update model
			resistivityVectorAll              = resistivityVectorAllBuf;

			std::vector<std::vector<double>> beforeApparentResistivityVector, beforePhaseVector;

			beforeApparentResistivityVector   = apparentResistivityVector;
			beforePhaseVector                 = phaseVector;

			apparentResistivityVector         = newApparentResistivityVector;
			phaseVector                       = newPhaseVector;
			residualApparentResistivityVector = newResidualApparentResistivityVector;
			residualPhaseVector               = newResidualPhaseVector;

			// Calculate RMS
			rmsVec.push_back(std::sqrt(objectiveFunction / (2 * sumNumData)));
			std::cout << "  RMS = " << std::sqrt(objectiveFunction / (2 * sumNumData)) << std::endl;
			std::cout << "  RMS = " << rmsVec[successfulIter] << std::endl;

			// Check delta for final convergence
			double abf(0), abdf(0);
			for (int i = 0 ; i < numStation ; i++) {
				for (int j = 0 ; j < numDataVec[i] ; j++) {
					abf  += std::abs(newApparentResistivityVector[i][j]) + std::abs(newPhaseVector[i][j]);
					abdf += std::abs(beforeApparentResistivityVector[i][j] - newApparentResistivityVector[i][j]) + std::abs(beforePhaseVector[i][j] - newPhaseVector[i][j]);
				}
			}
			delta = abdf / abf;
			std::cout << "  Delta = " << delta << std::endl;
			std::cout << "  Threshold of Delta = " << deltaThreshold << std::endl;

			if (delta <= deltaThreshold) {
				std::cout << "  STATUS: CONVERGED by delta threshold." << std::endl;
				break;
			}

		} else {
			std::cout << "  STATUS: FAILED - Objective function did not improve (New objective function = " << objectiveFunction << ")" << std::endl;
			lambda *= 10;  // Increase lambda

			if (lambda > 1e20) {
				std::cout << "  ERROR: Lambda is too large, stopping inversion" << std::endl;
				break;
			}
		}
	}


//	// Iteration start
//	for (int iter = 0 ; iter <= maxIteration ; iter++) {
//
//		std::cout << "\nIteration: " << iter << std::endl;
//
//		// Calculate Jacobian
////		std::cout << "Calculate Jacobian" << std::endl;
//		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(sumNumData * 2, numParameter * numStation);
//		Eigen::MatrixXd B = Eigen::VectorXd::Zero(sumNumData * 2);
//
//		int numDataBuf(0), numParameterBuf(0);
//		double resistivityBuf;
//
//		for (int i = 0 ; i < numStation ; i++) {
//
//			for (int j = 0 ; j < numParameter ; j++) {
//
//				resistivityBuf = resistivityVectorAll[i][j];
//				resistivityVectorAll[i][j] = resistivityBuf + differentialH;
//				auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resistivityVectorAll[i], thicknessVector);
//
//				for (int k = 0 ; k < numDataVec[i] ; k++) {
//
//					A(k + 2 * numDataBuf, j + numParameterBuf) = (apparentResistivityVectorBufVector[k] - apparentResistivityVector[i][k]) / differentialH * std::sqrt(observedAppResWeightVec[i][k]);
//					A(k + numDataVec[i] + 2 * numDataBuf, j + numParameterBuf) = (phaseVectorBufVector[k] - phaseVector[i][k]) / differentialH * std::sqrt(observedPhaseWeightVec[i][k]);
//
//				}
//			
//				resistivityVectorAll[i][j] = resistivityBuf;
//
//			}
//
//			numDataBuf += numDataVec[i];
//			numParameterBuf += numParameter;
//
//		}
//
//		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(numParameter * numStation, numParameter * numStation);
////		std::cout << "numParameter * numStation = " << numParameter * numStation << std::endl;
//
//		numDataBuf = 0;
//		numParameterBuf = 0;
//
//		for (int i = 0 ; i < numStation ; i++) {
//
//			for (int j = 0 ; j < numParameter ; j++) {
//
////				std::cout << "i = " << i << "  j = " << j << std::endl;
//
//
////				// A
////				if (i == 0 || i == numStation-1) {
////					if (j == 0) {
////						C(j + numParameterBuf, j + numParameterBuf) = std::sqrt(hGamma);
////					}
////					if (j != 0) {
////						C(j + numParameterBuf, j + numParameterBuf) = std::sqrt(vGamma) + std::sqrt(hGamma);
////						C(j + numParameterBuf, j + numParameterBuf - 1) = -1 * std::sqrt(vGamma);
////					}
////				}
////
////				// B
////				if (i != 0 && i != numStation-1) {
////					if (j == 0) {
////						C(j + numParameterBuf, j + numParameterBuf) = std::sqrt(hGamma);
////					}
////					if (j != 0) {
////						C(j + numParameterBuf, j + numParameterBuf) = std::sqrt(vGamma) + 2 * std::sqrt(hGamma);
////						C(j + numParameterBuf, j + numParameterBuf - 1) = -1 * std::sqrt(vGamma);
////					}
////				}
////
////				//C
////				if (i != 0 && i != numStation-1) {
////					C(j + numParameterBuf, j + numParameterBuf - numParameter) = -1 * std::sqrt(hGamma);
////					C(j + numParameterBuf - numParameter, j + numParameterBuf) = -1 * std::sqrt(hGamma);
////				}
//
//				// R1
//				if (j != 0) {
//					C(j + numParameterBuf, j + numParameterBuf    ) =      std::sqrt(vGamma);
//					C(j + numParameterBuf, j + numParameterBuf - 1) = -1 * std::sqrt(vGamma);
//				}
//
//				// R2
//				if (i == 0) {
//					C(j + numParameterBuf, j + numParameterBuf               ) =      std::sqrt(hGamma);
//					C(j + numParameterBuf, j + numParameterBuf + numParameter) = -1 * std::sqrt(hGamma);
//				} else if (i == numStation - 1) {
//					C(j + numParameterBuf, j + numParameterBuf               ) =      std::sqrt(hGamma);
//					C(j + numParameterBuf, j + numParameterBuf - numParameter) = -1 * std::sqrt(hGamma);
//				} else {
//					C(j + numParameterBuf, j + numParameterBuf               ) =  2 * std::sqrt(hGamma);
//					C(j + numParameterBuf, j + numParameterBuf - numParameter) = -1 * std::sqrt(hGamma);
//					C(j + numParameterBuf, j + numParameterBuf + numParameter) = -1 * std::sqrt(hGamma);
//				}
//
//			}
//
//			numParameterBuf += numParameter;
//
//		}
//
//		numParameterBuf = 0;
//		std::vector<double> numDataVecBuf(numDataVec.begin(), numDataVec.end());
//		numDataVecBuf.insert(numDataVecBuf.begin(), 0.0);
//
//		for (int i = 0 ; i < numStation ; i++) {
//
//			numDataBuf += numDataVecBuf[i] * 2;
//
//			for (int j = 0 ; j < numDataVec[i] ; j++) {
//
////				std::cout << "i = " << i << "  j = " << j << std::endl;
//				B(j + numDataBuf) = residualApparentResistivityVector[i][j];
//				B(j + numDataBuf + numDataVec[i]) = residualPhaseVector[i][j];
//
//			}
//
//		}
//
//		numDataBuf = 0;
//
//	// Nonlilnlear least squares method
//		Eigen::MatrixXd C5 = lambda * C;
//		Eigen::MatrixXd AandC5(A.rows() + C5.rows(), A.cols());
//		AandC5 << A, C5;
//
//		Eigen::VectorXd H = Eigen::VectorXd::Zero(numParameter * numStation);
//		Eigen::VectorXd BandH(B.size() + H.size());
//		BandH << B, H;
//
//		// Cholesky分解
////		Eigen::MatrixXd AtA = AandC5.transpose() * AandC5;
////		Eigen::VectorXd Atb = AandC5.transpose() * BandH;
////		Eigen::VectorXd dxCholesky = AtA.ldlt().solve(Atb);
//
//		// QR分解
////		Eigen::VectorXd dxQr = AandC5.colPivHouseholderQr().solve(BandH);
//
//		// SVD分解
//		Eigen::VectorXd dxSvd = AandC5.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(BandH);
//
//		// Calculate new data
//		std::vector<std::vector<double>> beforeApparentResistivityVector = apparentResistivityVector;
//		std::vector<std::vector<double>> beforePhaseVector = phaseVector;
//
//		std::vector<std::vector<double>> resistivityVectorAllBuf(numStation);
//
//		for (int i = 0 ; i < numStation ; i++) {
//
//			auto segment = dxSvd.segment(i * numParameter, numParameter);
//			std::vector<double> updatedResistivity(numParameter);
//			
//			for (int j = 0 ; j < numParameter ; j++) {
//
//				updatedResistivity[j] = resistivityVectorAll[i][j] + segment(j);
//
//			}
//
//			resistivityVectorAllBuf[i] = updatedResistivity;
//
//		}
//
//		double beforeObjectiveFunction = objectiveFunction;
//		objectiveFunction = 0;
//
//		for (int i = 0 ; i < numStation ; i++) {
//
//			auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resistivityVectorAllBuf[i], thicknessVector);
//
//			std::vector<double> residualApparentResistivityVectorBufVector(numDataVec[i]), residualPhaseVectorBufVector(numDataVec[i]);
//
//			for (int j = 0 ; j < numDataVec[i] ; j++) {
//
//				residualApparentResistivityVectorBufVector[j] = (observedAppResVec[i][j] - apparentResistivityVectorBufVector[j]) * std::sqrt(observedAppResWeightVec[i][j]);
//				residualPhaseVectorBufVector[j]               = (observedPhaseVec[i][j]  - phaseVectorBufVector[j]              ) * std::sqrt(observedPhaseWeightVec[i][j]) ;
////				residualBuf += residualApparentResistivityVectorBufVector[j] + residualPhaseVectorBufVector[j];
//				objectiveFunction += std::pow(residualApparentResistivityVectorBufVector[j], 2) + std::pow(residualPhaseVectorBufVector[j], 2);
//
//			}
//
//		apparentResistivityVector[i]         = apparentResistivityVectorBufVector;
//		phaseVector[i]                       = phaseVectorBufVector;
//		residualApparentResistivityVector[i] = residualApparentResistivityVectorBufVector;
//		residualPhaseVector[i]               = residualPhaseVectorBufVector;
//
//		}
//
//		std::cout << "  Objective Function = " << objectiveFunction << std::endl;
//
//		}
//
//		// Calculate RMS
//		double res2 = 0;
//		for (int i = 0 ; i < numStation ; i++) {
//
//			std::vector<double> resrho(numDataVec[i]), resphi(numDataVec[i]);
//
//			for (int j = 0 ; j < numDataVec[i] ; j++) {
//
//				resrho[j] = (observedAppResVec[i][j] - apparentResistivityVector[i][j]) * std::sqrt(observedAppResWeightVec[i][j]);
//				resphi[j] = (observedPhaseVec[i][j] - phaseVector[i][j]) * std::sqrt(observedPhaseWeightVec[i][j]);
//				res2 += std::pow(resrho[j], 2) + std::pow(resphi[j], 2);
//
//			}
//
//		}
//
//		rmsVec.push_back(std::sqrt(res2 / (2 * sumNumData)));
//		std::cout << "  RMS = " << rmsVec[iter] << std::endl;
//
//		// Check convergence
//		if (objectiveFunction < beforeObjectiveFunction) {
//			double abf(0), abdf(0);
//			for (int i = 0 ; i < numStation ; i++) {
//				for (int j = 0 ; j < numDataVec[i] ; j++) {
//					abf  += std::abs(apparentResistivityVector[i][j]) + std::abs(phaseVector[i][j]);
//					abdf += std::abs(beforeApparentResistivityVector[i][j] - apparentResistivityVector[i][j]) + std::abs(beforePhaseVector[i][j] - phaseVector[i][j]);
//				}
//			}
//
//			delta = abdf / abf;
//			std::cout << "  Delta = " << delta << std::endl;
//			std::cout << "  Threshold of Delta = " << deltaThreshold << std::endl;
//
//			if (delta <= deltaThreshold) {
//				for (int i = 0 ; i < numStation ; i++) {
//
//					auto segment = dxSvd.segment(i * numParameter, numParameter);
//					std::vector<double> updatedResistivity(numParameter);
//
//					for (int j = 0 ; j < numParameter ; j++) {
//
//						updatedResistivity[j] = resistivityVectorAll[i][j] + segment(j);
//
//						}
//
//						resistivityVectorAll[i] = updatedResistivity;
//
//					}
//					break;
//
//				} else {
//					for (int i = 0 ; i < numStation ; i++) {
//
//						auto segment = dxSvd.segment(i * numParameter, numParameter);
//						std::vector<double> updatedResistivity(numParameter);
//
//						for (int j = 0 ; j < numParameter ; j++) {
//
//							updatedResistivity[j] = resistivityVectorAll[i][j] + segment(j);
//
//						}
//
//						resistivityVectorAll[i] = updatedResistivity;
//
//					}
//
//				}
//
//			} else {
//					for (int i = 0 ; i < numStation ; i++) {
//
//						for (int j = 0 ; j < numParameter ; j++) {
//
//							dxSvd(j + i * numParameter) = 0;
//
//						}
//
//					}
//				delta = 1;
//			}
//
//		}

	// Output final result
	std::cout << "\nRESULT" << std::endl;
	for (int i = 0 ; i < numStation ; i++) {

		for (int j = 0 ; j < numDataVec[i] ; j++) {

			observedAppResVec[i][j]          = logRes2res(observedAppResVec[i][j]);
			apparentResistivityVector[i][j]  = logRes2res(apparentResistivityVector[i][j]);
			observedPhaseVec[i][j]           = rad2deg(observedPhaseVec[i][j]);
			phaseVector[i][j]                = rad2deg(phaseVector[i][j]);

		}

	}

	for (int i = 0 ; i < numStation ; i++) {

		std::cout << "\nPARAM. THICKNESS   DEPTH   RESISTIVITY" << std::endl;

		for ( int j = 0 ; j < numParameter ; j++) {

			std::cout << std::right << std::setw(2) << j << " "
                      << std::fixed << std::setprecision(4) 
                      << std::setw(10) << thicknessVector[j] << "  "
                      << std::setw(10) << depthVector[j] << "  "
                      << std::setw(8) << resistivityVectorAll[i][j] << std::endl;

		}

	}

	for (int i = 0 ; i < numStation ; i++) {

		std::ofstream ofs("result_resistivity_st" + std::to_string(i) + ".csv");
		ofs << "param,thickness,depth,resistivity" << std::endl;

		for ( int j = 0 ; j < numParameter ; j++) {

			ofs << j << ","
                << thicknessVector[j] << ","
                << depthVector[j] << ","
                << resistivityVectorAll[i][j] << std::endl;

		}

		ofs.close();

	}

	// Draw sounding curves
	std::string outputDirectoryName = "./";

	xlimMin = 100000;
	xlimMax = 0.1;

	for (int i = 0 ; i < numStation ; i++) {
		drawSoundingCurve("AppResis", "ohm.m", freqVec[i], apparentResistivityVector[i], observedAppResVec[i], xlimMin, xlimMax, ylimMin=0.1, ylimMax=10000, i, outputDirectoryName);
		drawSoundingCurve("Phase"   , "deg." , freqVec[i], phaseVector[i], observedPhaseVec[i], xlimMin, xlimMax, ylimMin=-180, ylimMax=180, i, outputDirectoryName);
	}

	// Contour
	draw2DSection(resistivityVectorAll, depthVector, distanceVec, numStation, numParameter, outputDirectoryName);

	return 0;

}
