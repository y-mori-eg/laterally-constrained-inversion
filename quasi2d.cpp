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
double ndiv, vStart, vEnd, initialResistivityValue, vGamma, hGamma;
double lambda = 2e-4;
std::string initialResistivity;

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

std::vector<std::vector<double>> observedAppResVec;
std::vector<std::vector<double>> observedPhaseVec;
std::vector<std::vector<double>> observedAppResWeightVec;
std::vector<std::vector<double>> observedPhaseWeightVec;

// Model data
std::vector<double> depthVector;
std::vector<double> thicknessVector;

// Calculated data
std::vector<std::vector<double>> resistivityVectorAll;
std::vector<std::vector<double>> apparentResistivityVector;
std::vector<std::vector<double>> phaseVector;
std::vector<std::vector<double>> residualApparentResistivityVector;
std::vector<std::vector<double>> residualPhaseVector;
std::vector<std::complex<double>> impedanceVector;

std::vector<double> rmsVec;

double delta(1);
int errCount(0);

double xlimMin, xlimMax, ylimMin, ylimMax;

std::string strBuf;


// Resistivity to Log Resistivity
double res2logRes(double resistivity) {
	return std::log10(resistivity);
}


// Log Resistivity to Resistivity
double logRes2res(double logResistivity) {
	return std::pow(10, logResistivity);
}


// Degree to Radian
double deg2rad(double degree) {
	return degree * PI / 180;
}


// Radian to Degree
double rad2deg(double radian) {
	return radian * 180 / PI;
}


// Check decimal
bool checkDecimal (std::string str) {
	try {
		std::stod(str);
		return true;
	} catch (const std::invalid_argument& e) {
		return false;
	}
}


// Phase unwrapping
std::vector<double> unwrapPhase(std::vector<double> phaseVector) {
	if (phaseVector.empty()) {
		return {};
	}

	std::vector<double> unwrappedPhase = phaseVector;

	for (size_t i = 1 ; i < unwrappedPhase.size() ; i++) {
		double diff = unwrappedPhase[i] - unwrappedPhase[i - 1];

		if (diff > 0.7 * PI && diff < 0.8 * 2 * PI) {
			unwrappedPhase[i] -= PI;
		} else if (diff < -0.7 * PI && diff > -0.8 * 2 * PI) {
			unwrappedPhase[i] += PI;
		}
	}

	return unwrappedPhase;
}


// Read parameters
std::tuple<std::string, double, double, double, double, double, double, std::string, double, double, double> readParameter(std::string paramFile) {

	std::cout << "Read parameters from '" << paramFile << "' ." << std::endl;

	std::ifstream ifs(paramFile.c_str(), std::ios::in);
	ifs >> dataFile >> maxIteration >> differentialH >> deltaThreshold >> ndiv >> vStart >> vEnd >> initialResistivity >> vGamma >> hGamma >> lambda;;

	std::cout << " Filename            : " << dataFile << std::endl << " Max Iteration       : " << maxIteration << std::endl << " Differential h      : " << differentialH << std::endl << " Delta Threshold     : " << deltaThreshold << std::endl << " ndiv                : " << ndiv << std::endl << " Range of v          : " << vStart << " - " << vEnd << std::endl << " Initial Resistivity : " << initialResistivity << std::endl << " Vertical Gamma      : " << vGamma << std::endl << " Horizontal Gamma    : " << hGamma << std::endl << " Initial Lambda          : " << lambda << std::endl << std::endl;

	return std::make_tuple(dataFile, maxIteration, differentialH, deltaThreshold, ndiv, vStart, vEnd, initialResistivity, vGamma, hGamma, lambda);

}


// Read observed data
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

//		phaseVec[i] = unwrapPhase(phaseVec[i]);

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


std::vector<std::vector<double>> setInitialResistivity(std::string initialResistivity) {

	std::cout << "Set initial resistivity: ";

	if (checkDecimal(initialResistivity)) {

		initialResistivityValue = stoi(initialResistivity);

		std::cout << initialResistivityValue << " Ohm-m." << std::endl << std::endl;

		for (int i = 0 ; i < numStation ; i++) {
			resistivityVectorAll[i].resize(numParameter);
			for (int j = 0 ; j < numParameter ; j++) {
				resistivityVectorAll[i][j] = res2logRes(initialResistivityValue);
			}
		}

	} else {
		std::ifstream ifs(initialResistivity.c_str(), std::ios::in);
		int numResisLayer;
		ifs >> numResisLayer;
		std::vector<double> resisValue(numResisLayer), resisThickness(numResisLayer-1);
		if (numResisLayer == 1) {
			// Read resistivity settings
			ifs >> strBuf;
			resisValue[0] = stod(strBuf);
			// Apply initial resistivity
			for (int i = 0 ; i < numStation ; i++) {
				resistivityVectorAll[i].resize(numParameter);
				for (int j = 0 ; j < numParameter ; j++) {
					resistivityVectorAll[i][j] = res2logRes(resisValue[0]);
				}
			}
		} else if (numResisLayer > 1) {
			// Read resistivity settings
			for (int i = 0 ; i < numResisLayer - 1 ; i++) {
				ifs >> resisValue[i] >> resisThickness[i];
			}
			ifs >> resisValue[numResisLayer - 1];
			// Apply initial resistivity
			for (int i = 0 ; i < numStation ; i++) {
				resistivityVectorAll[i].resize(numParameter);
				for (int j = 0 ; j < numParameter ; j++) {
					double resisDepthBuf(0);
					for (int k = 0 ; k < numResisLayer - 1 ; k++) {
						resisDepthBuf += resisThickness[k];
						if (depthVector[j] <= resisDepthBuf && depthVector[j] > resisDepthBuf - resisThickness[k]) {
							resistivityVectorAll[i][j] = res2logRes(resisValue[k]);
						} else if (depthVector[j] > resisDepthBuf){
							resistivityVectorAll[i][j] = res2logRes(resisValue[numResisLayer - 1]);
						}
					}
				}
			}
		} else {
			std::cerr << "Error: Number of layers not appropriate ." << std::endl;
		}
	}

	return resistivityVectorAll;

}


// Forward calculation
std::tuple<std::vector<double>, std::vector<double>> forwardCalc(std::vector<double> freqVec, std::vector<double> logResistivityVector, std::vector<double> thicknessVector) {

	std::vector<double> calculatedApparentResistivityVector(freqVec.size());
	std::vector<double> calculatedPhaseVector(freqVec.size());

	std::vector<double> resistivityVector(logResistivityVector.size());

	for (size_t i = 0 ; i < logResistivityVector.size() ; i++) {
		resistivityVector[i] = logRes2res(logResistivityVector[i]);
	}

	int numLayer = resistivityVector.size();

	double resistivityBuf, thicknessBuf;
	double apparentResistivityBuf, phaseBuf;
	
	apparentResistivityVector.resize(freqVec.size());
	phaseVector              .resize(freqVec.size());
	impedanceVector          .resize(numLayer);

	for (int i = 0 ; i < static_cast<int>(freqVec.size()) ; i++) {

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

		apparentResistivityBuf = std::log10((std::abs(impedanceVector[0]) * std::abs(impedanceVector[0])) / (OMEGA * MU));
		phaseBuf               = std::atan2(impedanceVector[0].imag(), impedanceVector[0].real());

		calculatedApparentResistivityVector[i] = apparentResistivityBuf;
		calculatedPhaseVector[i]               = phaseBuf;

	}

	return std::make_tuple(calculatedApparentResistivityVector, calculatedPhaseVector);

}


// Apparent Resistivity and Phase calculation
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> calcApparentResistivityAndPhase(std::vector<std::vector<double>> appResVector, std::vector<std::vector<double>> phsVector, std::vector<std::vector<double>> residualAppResVector, std::vector<std::vector<double>> residualPhsVector, std::vector<std::vector<double>> resisVector) {
	
	std::vector<double> apparentResistivityVectorBufVector, phaseVectorBufVector;

	sumNumData = 0;

	for (int i = 0 ; i < numStation ; i++) {

		std::cout << " Station No. " << i+1 << std::endl;

		auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resisVector[i], thicknessVector);

		sumNumData += numDataVec[i];
		
		std::vector<double> residualApparentResistivityVectorBufVector(numDataVec[i]), residualPhaseVectorBufVector(numDataVec[i]);

		for (int j = 0 ; j < numDataVec[i] ; j++) {
			
			residualApparentResistivityVectorBufVector[j] = (observedAppResVec[i][j] - apparentResistivityVectorBufVector[j]) * std::sqrt(observedAppResWeightVec[i][j]);
			residualPhaseVectorBufVector[j]               = (observedPhaseVec[i][j]  - phaseVectorBufVector[j]              ) * std::sqrt(observedPhaseWeightVec[i][j]) ;
		}

		appResVector[i]         = apparentResistivityVectorBufVector;
		phsVector[i]            = phaseVectorBufVector;
		residualAppResVector[i] = residualApparentResistivityVectorBufVector;
		residualPhsVector[i]    = residualPhaseVectorBufVector;

	}

	return std::make_tuple(appResVector, phsVector, residualAppResVector, residualPhsVector);

}


// Objective Function culculation
double calcObjectiveFunction(std::vector<std::vector<double>> residualAppResVector, std::vector<std::vector<double>> residualPhsVector) {

	double objectiveFunctionBuf(0);

	for (int i = 0 ; i < numStation ; i++) {
		for (int j = 0 ; j < numDataVec[i] ; j++) {
			objectiveFunctionBuf += std::pow(residualAppResVector[i][j], 2) + std::pow(residualPhsVector[i][j], 2);
		}
	}

	std::cout << "  Objective Function = " << objectiveFunctionBuf << std::endl;

	return objectiveFunctionBuf;

}


// RMS calculation
std::vector<double> calcRMS(std::vector<double> rmsVec, double objectiveFunction) {

	rmsVec.push_back(std::sqrt(objectiveFunction / (2 * sumNumData)));

	std::cout << "  RMS                = " << std::sqrt(objectiveFunction / (2 * sumNumData)) << std::endl;

	return rmsVec;

}


// Construct Jacobian A
Eigen::MatrixXd constructJacobian(Eigen::MatrixXd A, std::vector<std::vector<double>> apparentResistivityVector, std::vector<std::vector<double>> phaseVector) {

	std::cout << "  Construct Jacobian" << std::endl;

	int numDataBuf(0), numParameterBuf(0);

	for (int i = 0 ; i < numStation ; i++) {

		for (int j = 0 ; j < numParameter ; j++) {

			double resistivityBuf = resistivityVectorAll[i][j];
			resistivityVectorAll[i][j] = resistivityBuf + differentialH;
			auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resistivityVectorAll[i], thicknessVector);

			for (int k = 0 ; k < numDataVec[i] ; k++) {
				A(k + numDataBuf                , j + numParameterBuf) = (apparentResistivityVectorBufVector[k] - apparentResistivityVector[i][k]) / differentialH * std::sqrt(observedAppResWeightVec[i][k]);
				A(k + numDataBuf + numDataVec[i], j + numParameterBuf) = (phaseVectorBufVector[k]               - phaseVector[i][k]              ) / differentialH * std::sqrt(observedPhaseWeightVec[i][k] );
			}

			resistivityVectorAll[i][j] = resistivityBuf;

		}
		numDataBuf += 2 * numDataVec[i];
		numParameterBuf += numParameter;
	}

	return A;

}


// Construct Roughning Matrix C
Eigen::MatrixXd constructRoughningMatrix(Eigen::MatrixXd C) {

	std::cout << "  Construct Roughning Matrix" << std::endl;

	int numParameterBuf(0);

	for (int i = 0 ; i < numStation ; i++) {
		for (int j = 0 ; j < numParameter ; j++) {

			// R1 (Vertical Smoothing)
			if (j != 0) {
				C(j + numParameterBuf, j + numParameterBuf    ) +=      std::sqrt(vGamma);
				C(j + numParameterBuf, j + numParameterBuf - 1) += -1 * std::sqrt(vGamma);
			}

			// R2 (Horizontal Smoothing)
			if (numStation > 1) {
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

		}
		numParameterBuf += numParameter;
	}

	return C;

}


// Construct Residual Vector B
Eigen::VectorXd constructResidualVector(Eigen::VectorXd B, std::vector<std::vector<double>> residualApparentResistivityVector, std::vector<std::vector<double>> residualPhaseVector) {

	std::cout << "  Construct Residual Vector" << std::endl;

	int numDataBuf(0);

	std::vector<double> numDataVecBuf(numDataVec.begin(), numDataVec.end());
	numDataVecBuf.insert(numDataVecBuf.begin(), 0.0);

	for (int i = 0 ; i < numStation ; i++) {

		numDataBuf += 2 * numDataVecBuf[i];

		for (int j = 0 ; j < numDataVec[i] ; j++) {
			B(j + numDataBuf                ) = residualApparentResistivityVector[i][j];
			B(j + numDataBuf + numDataVec[i]) = residualPhaseVector[i][j];
		}

	}

	return B;

}


// Nonlinear Least Squares Method
std::vector<std::vector<double>> calcNewResistivity(Eigen::MatrixXd A, Eigen::MatrixXd B, Eigen::MatrixXd C, std::vector<std::vector<double>> resistivityVectorAllBuf) {

	Eigen::MatrixXd C5 = lambda * C;
	Eigen::MatrixXd AandC5(A.rows() + C5.rows(), A.cols());
	AandC5 << A, C5;

	Eigen::VectorXd H = Eigen::VectorXd::Zero(numParameter * numStation);
	Eigen::VectorXd BandH(B.size() + H.size());
	BandH << B, H;

	// Cholesky分解
//	Eigen::MatrixXd AtA = AandC5.transpose() * AandC5;
//	Eigen::VectorXd Atb = AandC5.transpose() * BandH;
//	Eigen::VectorXd dxCholesky = AtA.ldlt().solve(Atb);

	// QR分解
//	Eigen::VectorXd dxQr = AandC5.colPivHouseholderQr().solve(BandH);

	// SVD
	Eigen::VectorXd dxSvd = AandC5.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(BandH);

	for (int i = 0 ; i < numStation ; i++) {

		auto segment = dxSvd.segment(i * numParameter, numParameter);
		std::vector<double> updatedResistivity(numParameter);

		for (int j = 0 ; j < numParameter ; j++) {
			updatedResistivity[j] = resistivityVectorAll[i][j] + segment(j);
		}

		resistivityVectorAllBuf[i] = updatedResistivity;

	}

	return resistivityVectorAllBuf;

}


// Calculating Delta
double calcDelta(std::vector<std::vector<double>> newAppResVector, std::vector<std::vector<double>> newPhsVector, std::vector<std::vector<double>> beforeAppResVector, std::vector<std::vector<double>> beforePhsVector) {
	double abf(0), abdf(0);

	for (int i = 0 ; i < numStation ; i++) {
		for (int j = 0 ; j < numDataVec[i] ; j++) {
			abf  += std::abs(newAppResVector[i][j]                           ) + std::abs(newPhsVector[i][j]);
			abdf += std::abs(beforeAppResVector[i][j] - newAppResVector[i][j]) + std::abs(beforePhsVector[i][j] - newPhsVector[i][j]);
		}
	}

	delta = abdf / abf;

	std::cout << "  Delta              = " << delta << std::endl;
	std::cout << "  (ThresholdDelta    = " << deltaThreshold << ")" << std::endl;

	return delta;

}


// Draw sounding curves
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

 
// Draw 2D sections
void draw2DSection(std::vector<std::vector<double>> resistivityVector, std::vector<double> depthVector, std::vector<double> distanceVector, int numStation, int numParameter, std::string outputDirectoryName) {

	std::string outfileTitle = "outputDataForGnuplot.txt";
	std::ofstream outfile(outfileTitle);

	if (!outfile) {
		std::cerr << "Error: Cannot open " << outfileTitle << " for writing." << std::endl;
		return;
	}

	outfile << "# Distance(km)  Depth(km)  log10_Resistivity" << std::endl;

	for (int i = 0 ; i < numStation ; i++) {
		for (int j = 0 ; j < numParameter ; j++) {
			outfile << distanceVector[i] << " " << depthVector[j] << " " << resistivityVector[i][j] << std::endl;
		}
	}

	outfile.close();

 }


int main () {

	// Read parameters
	auto [dataFile, maxIteration, differentialH, deltaThreshold, ndiv, vStart, vEnd, initialResistivity, vGamma, hGamma, lambda] = readParameter("param.dat");

	// Read data
	std::tie (numStation, numDataVec, distanceVec, freqVec, observedAppResWeightVec, observedAppResVec, observedPhaseWeightVec, observedPhaseVec) = readData(dataFile);

	// Set initial parameters/options
	std::cout << "Initial parameters" << std::endl;
	factor = std::pow(10, 1.0/ndiv);
	std::cout << " factor = " << factor << std::endl;
	numParameter = (std::log10(vEnd) - std::log10(vStart)) * ndiv + 3;
	std::cout << " numParameter = " << numParameter << std::endl << std::endl;
	
	depthVector         .resize(numParameter+1);
	thicknessVector     .resize(numParameter);
	resistivityVectorAll.resize(numStation);

	depthVector[0]     = vStart;
	thicknessVector[0] = vStart;

	for (int i = 1 ; i < numParameter ; i++) {
		depthVector[i]      = depthVector[i-1] * factor;
		thicknessVector[i]  = depthVector[i] - depthVector[i-1];
	}

	// Set initial resistivity
	resistivityVectorAll = setInitialResistivity(initialResistivity);

	depthVector[numParameter] = depthVector[numParameter-1] * 1.1;

	// Calculate data
	std::vector<std::vector<double>> apparentResistivityVector(numStation), phaseVector(numStation);
	std::vector<std::vector<double>> residualApparentResistivityVector(numStation), residualPhaseVector(numStation);

	std::tie(apparentResistivityVector, phaseVector, residualApparentResistivityVector, residualPhaseVector) = calcApparentResistivityAndPhase(apparentResistivityVector, phaseVector, residualApparentResistivityVector, residualPhaseVector, resistivityVectorAll);

//	double residualBuf(0);

	double objectiveFunction(0);
	objectiveFunction = calcObjectiveFunction(residualApparentResistivityVector, residualPhaseVector);

	rmsVec = calcRMS(rmsVec, objectiveFunction);

	std::cout << std::endl;

	int iter(0);
	while (delta >= deltaThreshold) {

		if (iter >= maxIteration) {
			std::cout << " STATUS: Iteration number reached the threshold." << std::endl;
			break;
		}

		std::cout << "\n<-- Iteration: " << iter + 1 << " -->" << std::endl;
		std::cout << "  Current lambda: " << lambda << std::endl;

		// Calculate Jacobian
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(sumNumData * 2, numParameter * numStation);
		A = constructJacobian(A, apparentResistivityVector, phaseVector);

		// Construct Roughness Matrix C
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(numParameter * numStation, numParameter * numStation);
		C = constructRoughningMatrix(C);

		// Construct Residual Vector B
		Eigen::VectorXd B = Eigen::VectorXd::Zero(sumNumData * 2);
		B = constructResidualVector(B, residualApparentResistivityVector, residualPhaseVector);

		// Nonlilnlear Least Squares Method
		std::vector<std::vector<double>> resistivityVectorAllBuf(numStation);
		resistivityVectorAllBuf = calcNewResistivity(A, B, C, resistivityVectorAllBuf);

		// Calculate new objective function
		std::vector<std::vector<double>> newApparentResistivityVector(numStation), newPhaseVector(numStation), newResidualApparentResistivityVector(numStation), newResidualPhaseVector(numStation);

		std::tie(newApparentResistivityVector, newPhaseVector, newResidualApparentResistivityVector, newResidualPhaseVector) = calcApparentResistivityAndPhase(newApparentResistivityVector, newPhaseVector, newResidualApparentResistivityVector, newResidualPhaseVector, resistivityVectorAllBuf);

		double beforeObjectiveFunction = objectiveFunction;
		objectiveFunction = calcObjectiveFunction(newResidualApparentResistivityVector, newResidualPhaseVector);

		// Check convergence
		if (objectiveFunction < beforeObjectiveFunction) {

			std::cout << " STATUS: SUCCESS - Objective function improved - " << beforeObjectiveFunction << " --> " << objectiveFunction << std::endl;
			iter++;

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
			rmsVec = calcRMS(rmsVec, objectiveFunction);

			// Check delta for final convergence
			delta = calcDelta(apparentResistivityVector, phaseVector, beforeApparentResistivityVector, beforePhaseVector);

			if (delta <= deltaThreshold) {
				std::cout << " STATUS: CONVERGED by delta threshold." << std::endl;
				break;
			}

		} else {
			std::cout << " STATUS: FAILED - Objective function did not improve (New objective function = " << objectiveFunction << ")" << std::endl;
			errCount += 1;
			if (errCount > 100) {
				break;
			}
		}
	}

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
                      << std::setw(8) << logRes2res(resistivityVectorAll[i][j]) << std::endl;

		}

	}

	std::ofstream ofs("result_app_resis_and_phase.dat");
	ofs << numStation << std::endl;
	for (int i = 0 ; i < numStation ; i++) {

		ofs << i+1 << "    " << numDataVec[i] << "    " << distanceVec[i] << std::endl;

		for (int j = 0 ; j < numDataVec[i] ; j++) {
			ofs << std::scientific << std::setprecision(10) << freqVec[i][j] << "    "
				<< std::scientific << std::setprecision(10) << apparentResistivityVector[i][j] << "    "
				<< std::scientific << std::setprecision(10) << appResWeightVec[i][j] << "    "
				<< std::scientific << std::setprecision(10) << phaseVector[i][j] << "    "
				<< std::scientific << std::setprecision(10) << phaseWeightVec[i][j] << std::endl;
		}

	}

	ofs.close();

	for (int i = 0 ; i < numStation ; i++) {

		std::ofstream ofs("result_resistivity_st" + std::to_string(i) + ".csv");
		ofs << "param,thickness,depth,resistivity" << std::endl;

		for ( int j = 0 ; j < numParameter ; j++) {

			ofs << j << ","
                << thicknessVector[j] << ","
                << depthVector[j] << ","
                << logRes2res(resistivityVectorAll[i][j]) << std::endl;

		}

		ofs.close();

	}

	// Draw sounding curves
	std::string outputDirectoryName = "./";

	xlimMin = 100000;
	xlimMax = 0.001;

	for (int i = 0 ; i < numStation ; i++) {
		drawSoundingCurve("AppResis", "ohm.m", freqVec[i], apparentResistivityVector[i], observedAppResVec[i], xlimMin, xlimMax, ylimMin=0.1, ylimMax=10000, i, outputDirectoryName);
		drawSoundingCurve("Phase"   , "deg." , freqVec[i], phaseVector[i], observedPhaseVec[i], xlimMin, xlimMax, ylimMin=-180, ylimMax=180, i, outputDirectoryName);
	}

	// Contour
	draw2DSection(resistivityVectorAll, depthVector, distanceVec, numStation, numParameter, outputDirectoryName);

	return 0;

}
