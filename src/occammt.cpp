#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>
#include <cmath>
#include <complex>

#include <Eigen/Dense>

#include <chrono>
#include <omp.h>


const double PI = 3.1415926535897932384626433832795028841971;
const double MU = 4 * PI * 1e-7;  // Magnetic permeability [H/m]
const std::complex<double> imaginaryUnit(0.0, 1.0);

std::chrono::system_clock::time_point clockStart;
std::ofstream ofsLog("occammt.log");

std::string dataFile;
int numStation, stationNumber, numData;
int numStationRow, numStationColumn;
double coordinateX;
double coordinateY;

int maxIteration = 100;
double differentialH = 1e-3;
double initialDelta = 1;
double deltaThreshold = 1;
std::string initialResistivity, decompositionMethod;
double ndiv, vStart, vEnd, initialResistivityValue, vGamma, h1Gamma, h2Gamma;
double lambda = 2e-4;
double increaseFactor = 2;
double decreaseFactor = 2;
int openMpTF;

double factor;
int numParameter;
int sumNumData(0);

// Observed data
std::vector<int> numDataVec;
std::vector<double> coordinateXVec;
std::vector<double> coordinateYVec;
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
double rmsBuf(0);

double delta(1e+20);
int errCount(0);

std::string outputFileConvergence = "result_convergence.dat";

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


std::chrono::system_clock::time_point clockNow() {
	return std::chrono::system_clock::now();
}


auto elapsedTime(std::chrono::system_clock::time_point start, std::chrono::system_clock::time_point end) {
	return std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
}


// Read parameters
void readParameter(std::string paramFile) {

	ofsLog << "Read parameters from '" << paramFile << "' . (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

	std::ifstream ifs(paramFile.c_str(), std::ios::in);
	ifs >> dataFile >> maxIteration >> differentialH >> deltaThreshold >> ndiv >> vStart >> vEnd >> initialResistivity >> vGamma >> h1Gamma >> h2Gamma >> lambda >> increaseFactor >> decreaseFactor >> decompositionMethod >> openMpTF;

	ofsLog << " Filename               : " << dataFile << std::endl << " Max Iteration          : " << maxIteration << std::endl << " Differential h         : " << differentialH << std::endl << " Delta Threshold        : " << deltaThreshold << std::endl << " ndiv                   : " << ndiv << std::endl << " Range of v             : " << vStart << " - " << vEnd << std::endl << " Initial Resistivity    : " << initialResistivity << std::endl << " Vertical Gamma         : " << vGamma << std::endl << " Horizontal Gamma       : " << h1Gamma << ", " << h2Gamma << std::endl << " Initial Lambda          : " << lambda << std::endl << " Increase Factor         : " << increaseFactor << std::endl << " Decrease Factor         : " << decreaseFactor << std::endl << " Decomposition Method    : " << decompositionMethod << std::endl << " OpenMP (0:OFF, 1:ON)    : " << openMpTF << std::endl;

}


// Read observed data
std::tuple<int, std::vector<int>, std::vector<double>, std::vector<double>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> readData (const std::string& dataFile) {

	ofsLog << "Read data from '" << dataFile << "' . (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

	std::ifstream ifs(dataFile.c_str(), std::ios::in);
	ifs >> numStation >> numStationRow >> numStationColumn;
//	ifs >> numStation;
	ofsLog << " Number of Stations : " << numStation << std::endl;

	numDataVec     .resize(numStation);
	coordinateXVec .resize(numStation);
	coordinateYVec .resize(numStation);
	freqVec        .resize(numStation);
	appResWeightVec.resize(numStation);
	appResVec      .resize(numStation);
	phaseWeightVec .resize(numStation);
	phaseVec       .resize(numStation);

	for (int i = 0 ; i < numStation ; i++) {

		ifs >> stationNumber >> numData >> coordinateX >> coordinateY;

		numDataVec[i]  = numData;
		freqVec[i]        .resize(numData);
		appResWeightVec[i].resize(numData);
		appResVec[i]      .resize(numData);
		phaseWeightVec[i] .resize(numData);
		phaseVec[i]       .resize(numData);

		coordinateXVec[i] = coordinateX;
		coordinateYVec[i] = coordinateY;

		ofsLog << " Read station : " << stationNumber << std::endl;

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

	ofsLog << std::endl;

	for (int i = 0 ; i < numStation ; i++) {

		ofsLog << " Station : " << i+1 << std::endl;
		ofsLog << "  Freq                App.Res     Weight(App.Res)              Phase       Weight(Phase)" << std::endl;

		for (int j = 0 ; j < numDataVec[i] ; j++) {

			ofsLog << "  " << std::setw(20) << std::left << freqVec[i][j]
                           << std::setw(20) << std::left << appResVec[i][j]
                           << std::setw(20) << std::left << appResWeightVec[i][j]
                           << std::setw(20) << std::left << phaseVec[i][j]
                           << std::setw(20) << std::left << phaseWeightVec[i][j] << std::endl; 

		}

		ofsLog << std::endl;

	}

	return std::make_tuple(numStation, numDataVec, coordinateXVec, coordinateYVec, freqVec, appResWeightVec, appResVec, phaseWeightVec, phaseVec);

}


std::vector<std::vector<double>> setInitialResistivity(std::string initialResistivity) {
// 現在、水平多層でしか初期値を与えられない。
// ファイルとして比抵抗構造を与えることで、測点ごとに異なる初期値を設定できるように修正を加える。

	ofsLog << "Set initial resistivity: ";

	if (checkDecimal(initialResistivity)) {

		initialResistivityValue = stoi(initialResistivity);

		ofsLog << initialResistivityValue << " Ohm-m." << std::endl << std::endl;

		for (int i = 0 ; i < numStation ; i++) {
			resistivityVectorAll[i].resize(numParameter);
			for (int j = 0 ; j < numParameter ; j++) {
				resistivityVectorAll[i][j] = res2logRes(initialResistivityValue);
			}
		}

	} else {
		std::ifstream ifs(initialResistivity.c_str(), std::ios::in);
		int numResisLayer;
		for (int i = 0 ; i < numStation ; i++) {
			ifs >> stationNumber >> numResisLayer;
			std::vector<double> resisValue(numResisLayer), resisThickness(numResisLayer-1);
			if (numResisLayer == 1) {
				// Read resistivity settings
				ifs >> strBuf;
				resisValue[0] = stod(strBuf);
				// Apply initial resistivity
				resistivityVectorAll[i].resize(numParameter);
				for (int j = 0 ; j < numParameter ; j++) {
					resistivityVectorAll[i][j] = res2logRes(resisValue[0]);
				}
			} else if (numResisLayer > 1) {
				// Read resistivity settings
				for (int k = 0 ; k < numResisLayer - 1 ; k++) {
					ifs >> resisValue[k] >> resisThickness[k];
				}
				ifs >> resisValue[numResisLayer - 1];
				// Apply initial resistivity
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
			} else {
				std::cerr << "Error: Number of layers not appropriate ." << std::endl;
			}
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

	for (int i = 0 ; i < numStation ; i++) {

//		ofsLog << " Station No. " << i+1 << std::endl;

		auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], resisVector[i], thicknessVector);

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


// Apparent Resistivity and Phase calculation
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>> calcApparentResistivityAndPhaseParallel(std::vector<std::vector<double>> appResVector, std::vector<std::vector<double>> phsVector, std::vector<std::vector<double>> residualAppResVector, std::vector<std::vector<double>> residualPhsVector, std::vector<std::vector<double>> resisVector) {
	
	std::vector<double> apparentResistivityVectorBufVector, phaseVectorBufVector;

	sumNumData = 0;

	#pragma omp parallel for reduction(+:sumNumData)
	for (int i = 0 ; i < numStation ; i++) {

//		ofsLog << " Station No. " << i+1 << std::endl;

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

	ofsLog << "  Objective Function = " << objectiveFunctionBuf << std::endl;

	return objectiveFunctionBuf;

}


// RMS calculation
std::vector<double> calcRMS(std::vector<double> rmsVec, double objectiveFunction) {

	rmsBuf = std::sqrt(objectiveFunction / (2 * sumNumData));
	rmsVec.push_back(rmsBuf);

	ofsLog << "  RMS                = " << rmsBuf << std::endl;

	return rmsVec;

}


// Construct Jacobian A
void constructJacobian(Eigen::MatrixXd& A, std::vector<std::vector<double>> apparentResistivityVector, std::vector<std::vector<double>> phaseVector) {

	ofsLog << "  Construct Jacobian (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

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

}


// Construct Jacobian A (OpenMP)
void constructJacobianParallel(Eigen::MatrixXd& A, std::vector<std::vector<double>>& apparentResistivityVector, std::vector<std::vector<double>>& phaseVector) {

	ofsLog << "  Construct Jacobian (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

	#pragma omp parallel for
	for (int i = 0 ; i < numStation ; i++) {

		int numDataBuf(0);

		for(int k=0; k < i; ++k) {
			numDataBuf += 2 * numDataVec[k];
		}

		int numParameterBuf = i * numParameter;

		auto localResistivityVector = resistivityVectorAll[i];

		for (int j = 0 ; j < numParameter ; j++) {

			double resistivityBuf = localResistivityVector[j];
			localResistivityVector[j] = resistivityBuf + differentialH;
			auto [apparentResistivityVectorBufVector, phaseVectorBufVector] = forwardCalc(freqVec[i], localResistivityVector, thicknessVector);

			for (int k = 0 ; k < numDataVec[i] ; k++) {
				A(k + numDataBuf                , j + numParameterBuf) = (apparentResistivityVectorBufVector[k] - apparentResistivityVector[i][k]) / differentialH * std::sqrt(observedAppResWeightVec[i][k]);
				A(k + numDataBuf + numDataVec[i], j + numParameterBuf) = (phaseVectorBufVector[k]               - phaseVector[i][k]              ) / differentialH * std::sqrt(observedPhaseWeightVec[i][k] );
			}

			localResistivityVector[j] = resistivityBuf;

		}
	}

}


// Construct Roughning Matrix C
void constructRoughningMatrix(Eigen::MatrixXd& C) {

	ofsLog << "  Construct Roughning Matrix (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

	int numParameterBuf(0);

	for (int i = 0 ; i < numStation ; i++) {
		for (int j = 0 ; j < numParameter ; j++) {

			// R1 (Vertical Smoothing)
			if (j != 0) {
				C(j + numParameterBuf, j + numParameterBuf    ) +=      std::sqrt(vGamma);
				C(j + numParameterBuf, j + numParameterBuf - 1) += -1 * std::sqrt(vGamma);
			}

			// R2 (Horizontal Smoothing - 1)
			if (numStationRow > 1) {
				if (i % numStationRow == 0) {
					C(j + numParameterBuf, j + numParameterBuf               ) +=      std::sqrt(h1Gamma);
					C(j + numParameterBuf, j + numParameterBuf + numParameter) += -1 * std::sqrt(h1Gamma);
				} else if (i % numStationRow == numStationRow - 1) {
					C(j + numParameterBuf, j + numParameterBuf               ) +=      std::sqrt(h1Gamma);
					C(j + numParameterBuf, j + numParameterBuf - numParameter) += -1 * std::sqrt(h1Gamma);
				} else {
					C(j + numParameterBuf, j + numParameterBuf               ) +=  2 * std::sqrt(h1Gamma);
					C(j + numParameterBuf, j + numParameterBuf - numParameter) += -1 * std::sqrt(h1Gamma);
					C(j + numParameterBuf, j + numParameterBuf + numParameter) += -1 * std::sqrt(h1Gamma);
				}
			}

			// R3 (Horizontal Smoothing - 2)
			//ofsLog << "numStationColumn = " << numStationColumn << std::endl;
			if (numStationColumn > 1) {
				if (i < numStationRow) {
					C(j + numParameterBuf, j + numParameterBuf                               ) +=      std::sqrt(h2Gamma);
					C(j + numParameterBuf, j + numParameterBuf + numParameter * numStationRow) += -1 * std::sqrt(h2Gamma);
				} else if (i > numStation - numStationRow - 1) {
					C(j + numParameterBuf, j + numParameterBuf                               ) +=      std::sqrt(h2Gamma);
					C(j + numParameterBuf, j + numParameterBuf - numParameter * numStationRow) += -1 * std::sqrt(h2Gamma);
				} else {
					C(j + numParameterBuf, j + numParameterBuf                               ) +=  2 * std::sqrt(h2Gamma);
					C(j + numParameterBuf, j + numParameterBuf - numParameter * numStationRow) += -1 * std::sqrt(h2Gamma);
					C(j + numParameterBuf, j + numParameterBuf + numParameter * numStationRow) += -1 * std::sqrt(h2Gamma);
				}
			}

		}
		numParameterBuf += numParameter;
	}

}


// Construct Residual Vector B
void constructResidualVector(Eigen::VectorXd& B, std::vector<std::vector<double>> residualApparentResistivityVector, std::vector<std::vector<double>> residualPhaseVector) {

	ofsLog << "  Construct Residual Vector (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

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

}


// Nonlinear Least Squares Method
void calcNewResistivity(const Eigen::MatrixXd& A, const Eigen::VectorXd& B, const Eigen::MatrixXd& C, std::vector<std::vector<double>>& resistivityVectorAllBuf) {

//	ofsLog << "Checkpoint 4-1" << std::endl;
//	Eigen::MatrixXd C5 = lambda * C;
//	Eigen::MatrixXd AandC5(A.rows() + C5.rows(), A.cols());
//	AandC5 << A, C5;


//	ofsLog << "Checkpoint 4-2" << std::endl;
//	Eigen::VectorXd H = Eigen::VectorXd::Zero(numParameter * numStation);
//	Eigen::VectorXd BandH(B.size() + H.size());
//	BandH << B, H;

//	Eigen::VectorXd dxQr = AandC5.colPivHouseholderQr().solve(BandH);

	Eigen::VectorXd dx;

	if (decompositionMethod == "ldlt") {
	// QR分解

//		ofsLog << "Checkpoint 4-3" << std::endl;
	
		Eigen::MatrixXd AtA = A.transpose() * A;
		Eigen::MatrixXd CtC = C.transpose() * C;
		Eigen::VectorXd AtB = A.transpose() * B;
		Eigen::MatrixXd KtK = AtA + (lambda * lambda) * CtC;
		dx = KtK.ldlt().solve(AtB);

		ofsLog << "  Solve determinants - ldlt (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;
	
	} else if (decompositionMethod == "svd") {
	// SVD

//		ofsLog << "Checkpoint 4-3" << std::endl;

		Eigen::MatrixXd lambdaC = lambda * C;
		Eigen::MatrixXd AandlambdaC(A.rows() + lambdaC.rows(), A.cols());
		AandlambdaC << A, lambdaC;
	
		Eigen::VectorXd H = Eigen::VectorXd::Zero(numParameter * numStation);
		Eigen::VectorXd BandH(B.size() + H.size());
		BandH << B, H;
	
		dx = AandlambdaC.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(BandH);

		ofsLog << "  Solve determinants - svd (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

	} else {
		std::cerr << "Error: Inappropriate word (decompositionMethod) ." << std::endl;
		return;
	}

//	ofsLog << "Checkpoint 4-4" << std::endl;
	for (int i = 0 ; i < numStation ; i++) {

		auto segment = dx.segment(i * numParameter, numParameter);
		std::vector<double> updatedResistivity(numParameter);

		for (int j = 0 ; j < numParameter ; j++) {
			updatedResistivity[j] = resistivityVectorAll[i][j] + segment(j);
		}

		resistivityVectorAllBuf[i] = updatedResistivity;

	}

	ofsLog << "  Update resistivity vector (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

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

	ofsLog << "  Delta              = " << delta << std::endl;
	ofsLog << "  (ThresholdDelta    = " << deltaThreshold << ")" << std::endl;

	return delta;

}


void outputForwardResult(std::string dataType, int stationNumber, std::vector<double> freq, std::vector<double> observedData, std::vector<double> calculatedData) {

	std::ofstream ofsForwardResult("result_forward_" + dataType + "_st" + std::to_string(stationNumber) + ".dat");
	ofsForwardResult << "                  freq         observed data       calculated data" << std::endl;
	for ( int i = 0 ; i < numDataVec[stationNumber] ; i++ ) {
		ofsForwardResult << std::right << std::setw(22) << std::scientific << std::setprecision(10) << freq[i]
                         << std::right << std::setw(22) << std::scientific << std::setprecision(10) << observedData[i]
                         << std::right << std::setw(22) << std::scientific << std::setprecision(10) << calculatedData[i] << std::endl;
	}

}

 
// Draw 2D sections
void draw2DSection(std::vector<std::vector<double>> resistivityVector, std::vector<double> depthVector, std::vector<double> coordinateXVec, std::vector<double> coordinateYVec, int numStation, int numParameter, std::string outputDirectoryName) {

	for (int i = 0 ; i < numStationColumn ; i++) {

		std::string outfileTitle = "outputDataForGnuplot_column" + std::to_string(i) + ".txt";
		std::ofstream outfile(outfileTitle);

		if (!outfile) {
			std::cerr << "Error: Cannot open " << outfileTitle << " for writing." << std::endl;
			return;
		}

		outfile << "# Distance(km)   Depth(km)   log10_Resistivity" << std::endl;

		for (int j = 0 ; j < numStationRow ; j++) {
			for (int k = 0 ; k < numParameter ; k++) {
				outfile << coordinateXVec[i * numStationRow + j] << "  " << depthVector[k] << "  " << resistivityVector[i * numStationRow + j][k] << std::endl;
			}
		}

		outfile.close();

	}

	for (int i = 0 ; i < numStationRow ; i++) {

		std::string outfileTitle = "outputDataForGnuplot_row" + std::to_string(i) + ".txt";
		std::ofstream outfile(outfileTitle);

		if (!outfile) {
			std::cerr << "Error: Cannot open " << outfileTitle << " for writing." << std::endl;
			return;
		}

		outfile << "# Distance(km)   Depth(km)   log10_Resistivity" << std::endl;

		for (int j = 0 ; j < numStation ; j += numStationRow) {
			for (int k = 0 ; k < numParameter ; k++) {
				outfile << coordinateYVec[i + j] << "  " << depthVector[k] << "  " << resistivityVector[i + j][k] << std::endl;
			}
		}

		outfile.close();

	}

}


int main () {

	clockStart = clockNow();

	// Read parameters
	readParameter("param.dat");

	ofsLog << "Maximum threads: " << omp_get_max_threads() << std::endl;

	// Read data
	std::tie (numStation, numDataVec, coordinateXVec, coordinateYVec, freqVec, observedAppResWeightVec, observedAppResVec, observedPhaseWeightVec, observedPhaseVec) = readData(dataFile);

	sumNumData = 0;
	for (int i = 0 ; i < numStation ; i++) {
		sumNumData += numDataVec[i];
	}

	// Set initial parameters/options
	ofsLog << "Initial parameters" << std::endl;
	factor = std::pow(10, 1.0/ndiv);
	ofsLog << " factor = " << factor << std::endl;
	numParameter = (std::log10(vEnd) - std::log10(vStart)) * ndiv + 3;
	ofsLog << " numParameter = " << numParameter << std::endl << std::endl;
	
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

	if (openMpTF == 0) {
		std::tie(apparentResistivityVector, phaseVector, residualApparentResistivityVector, residualPhaseVector) = calcApparentResistivityAndPhase        (apparentResistivityVector, phaseVector, residualApparentResistivityVector, residualPhaseVector, resistivityVectorAll);
	} else if (openMpTF == 1) {
		std::tie(apparentResistivityVector, phaseVector, residualApparentResistivityVector, residualPhaseVector) = calcApparentResistivityAndPhaseParallel(apparentResistivityVector, phaseVector, residualApparentResistivityVector, residualPhaseVector, resistivityVectorAll);
	} else {
		std::cerr << "Error: Inappropriate number (openMpTF) ." << std::endl;
		return 1;
	}

//	double residualBuf(0);

	double objectiveFunction(0);
	objectiveFunction = calcObjectiveFunction(residualApparentResistivityVector, residualPhaseVector);

	rmsVec = calcRMS(rmsVec, objectiveFunction);

	std::ofstream ofsConv(outputFileConvergence);
	ofsConv << "  Iter.                  lambda         Objective Func.                     RMS                   Delta" << std::endl;

	int iter(0);
	while (delta >= deltaThreshold) {

		if (iter >= maxIteration) {
			ofsLog << " STATUS: Iteration number reached the threshold." << std::endl;
			break;
		}

		ofsLog << "\n<-- Iteration: " << iter + 1 << " -->" << std::endl;
//		ofsLog << "Checkpoint 1" << std::endl;
//		ofsLog << "  Current lambda: " << lambda << std::endl;

		// Calculate Jacobian
		Eigen::MatrixXd A = Eigen::MatrixXd::Zero(sumNumData * 2, numParameter * numStation);
		if (openMpTF == 0) {
			constructJacobian        (A, apparentResistivityVector, phaseVector);
		} else if (openMpTF == 1) {
			constructJacobianParallel(A, apparentResistivityVector, phaseVector);
		} else {
			std::cerr << "Error: Inappropriate number (openMpTF) ." << std::endl;
			return 1;
		}
//		ofsLog << "Checkpoint 2" << std::endl;
		ofsLog << "  Construct Jacobian (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

		// Construct Roughness Matrix C
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(numParameter * numStation, numParameter * numStation);
		constructRoughningMatrix(C);
//		ofsLog << "Checkpoint 3" << std::endl;
		ofsLog << "  Construct roughness matrix (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

		// Construct Residual Vector B
		Eigen::VectorXd B = Eigen::VectorXd::Zero(sumNumData * 2);
		constructResidualVector(B, residualApparentResistivityVector, residualPhaseVector);
//		ofsLog << "Checkpoint 4" << std::endl;
		ofsLog << "  Construct residual vector (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

		// Nonlilnlear Least Squares Method
		std::vector<std::vector<double>> resistivityVectorAllBuf(numStation);
		calcNewResistivity(A, B, C, resistivityVectorAllBuf);
//		ofsLog << "Checkpoint 5" << std::endl;
		ofsLog << "  Construct residual vector (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

		// Calculate new objective function
		std::vector<std::vector<double>> newApparentResistivityVector(numStation), newPhaseVector(numStation), newResidualApparentResistivityVector(numStation), newResidualPhaseVector(numStation);
//		ofsLog << "Checkpoint 6" << std::endl;
		ofsLog << "  Calculate new objective function (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

		if (openMpTF == 0) {
			std::tie(newApparentResistivityVector, newPhaseVector, newResidualApparentResistivityVector, newResidualPhaseVector) = calcApparentResistivityAndPhase        (newApparentResistivityVector, newPhaseVector, newResidualApparentResistivityVector, newResidualPhaseVector, resistivityVectorAllBuf);
		} else if (openMpTF == 1) {
			std::tie(newApparentResistivityVector, newPhaseVector, newResidualApparentResistivityVector, newResidualPhaseVector) = calcApparentResistivityAndPhaseParallel(newApparentResistivityVector, newPhaseVector, newResidualApparentResistivityVector, newResidualPhaseVector, resistivityVectorAllBuf);
		} else {
			std::cerr << "Error: Inappropriate number (openMpTF) ." << std::endl;
			return 1;
		}

//		ofsLog << "Checkpoint 7" << std::endl;
		ofsLog << "  Calculate new apparent resistivity and phase data (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

		double beforeObjectiveFunction = objectiveFunction;
		objectiveFunction = calcObjectiveFunction(newResidualApparentResistivityVector, newResidualPhaseVector);

		// Check convergence
		if (objectiveFunction < beforeObjectiveFunction) {

			ofsLog << " STATUS: SUCCESS - Objective function improved - " << beforeObjectiveFunction << " --> " << objectiveFunction << std::endl;
			iter++;
			lambda /= decreaseFactor;

			// Update model
			std::vector<std::vector<double>> beforeApparentResistivityVector, beforePhaseVector;
			resistivityVectorAll              = resistivityVectorAllBuf;

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
				ofsLog << " STATUS: CONVERGED by delta threshold." << std::endl;
				ofsConv << std::right << std::setw(7) << iter
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << lambda
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << objectiveFunction
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << rmsBuf
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << delta << std::endl;
				break;
			}

		} else {
			ofsLog << " STATUS: FAILED - Objective function did not improve (New objective function = " << objectiveFunction << ")" << std::endl;
			lambda *= increaseFactor;
			errCount += 1;
			if (errCount > 100) {
				ofsConv << std::right << std::setw(7) << iter
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << lambda
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << objectiveFunction
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << rmsBuf
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << delta << std::endl;
				break;
			}
		}

		ofsConv << std::right << std::setw(7) << iter
                        << std::right << std::setw(24) << std::scientific << std::setprecision(10) << lambda
                << std::right << std::setw(24) << std::scientific << std::setprecision(10) << objectiveFunction
                << std::right << std::setw(24) << std::scientific << std::setprecision(10) << rmsBuf
                << std::right << std::setw(24) << std::scientific << std::setprecision(10) << delta << std::endl;

		ofsLog << "(" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;
	}
	ofsConv.close();

	// Output final result
	ofsLog << "\nRESULT" << std::endl;
	for (int i = 0 ; i < numStation ; i++) {

		for (int j = 0 ; j < numDataVec[i] ; j++) {

			observedAppResVec[i][j]          = logRes2res(observedAppResVec[i][j]);
			apparentResistivityVector[i][j]  = logRes2res(apparentResistivityVector[i][j]);
			observedPhaseVec[i][j]           = rad2deg(observedPhaseVec[i][j]);
			phaseVector[i][j]                = rad2deg(phaseVector[i][j]);

		}

	}

	for (int i = 0 ; i < numStation ; i++) {

		ofsLog << std::endl;
		ofsLog << "Station No." << i + 1 << std::endl;
		ofsLog << "PARAM. THICKNESS   DEPTH   RESISTIVITY" << std::endl;

		for ( int j = 0 ; j < numParameter ; j++) {

			ofsLog << std::right << std::setw(2) << j << " "
                      << std::fixed << std::setprecision(4) 
                      << std::setw(10) << thicknessVector[j] << "  "
                      << std::setw(10) << depthVector[j] << "  "
                      << std::setw(8) << logRes2res(resistivityVectorAll[i][j]) << std::endl;

		}

	}

	std::ofstream ofs("result_app_resis_and_phase.dat");
	ofs << numStation << "    " << numStationRow << "    " << numStationColumn << std::endl;
	for (int i = 0 ; i < numStation ; i++) {

		ofs << i+1 << "    " << numDataVec[i] << "    " << coordinateXVec[i] << "    " << coordinateYVec[i] << std::endl;

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
//
//	xlimMin = 100000;
//	xlimMax = 0.001;
//
//	for (int i = 0 ; i < numStation ; i++) {
//		drawSoundingCurve("AppResis", "ohm.m", freqVec[i], apparentResistivityVector[i], observedAppResVec[i], xlimMin, xlimMax, ylimMin=0.1, ylimMax=10000, i, outputDirectoryName);
//		drawSoundingCurve("Phase"   , "deg." , freqVec[i], phaseVector[i], observedPhaseVec[i], xlimMin, xlimMax, ylimMin=-180, ylimMax=180, i, outputDirectoryName);
//	}

	// Output sounding curve data
	for (int i = 0 ; i < numStation ; i++) {
		outputForwardResult("app_resis", i, freqVec[i], observedAppResVec[i], apparentResistivityVector[i]);
		outputForwardResult("phase",     i, freqVec[i], observedPhaseVec[i] , phaseVector[i]              );
	}

	// Contour
	draw2DSection(resistivityVectorAll, depthVector, coordinateXVec, coordinateYVec, numStation, numParameter, outputDirectoryName);

	ofsLog << "Finish inversion (" << elapsedTime(clockStart, clockNow()) << " msec)" << std::endl;

	ofsLog.close();

	return 0;

}
