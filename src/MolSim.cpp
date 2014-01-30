#define _USE_MATH_DEFINES

#include "outputWriter/VTKWriter.h"
#include "FileReader.h"

#include "Particle.h"
#include "Cuboid.h"
#include "Sphere.h"

#include "utils/Vector.h"
#include "utils/ParticleContainer.h"
#include "utils/ParticleIterator.h"
#include "utils/ParticleGenerator.h"
#include "utils/LCParticleContainer.h"
#include "utils/LCInnerParticleIterator.h"
#include "utils/LCOuterParticleIterator.h"
#include "utils/Thermostat.h"
#include "utils/BoundaryHandler.h"

#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/xml/domconfigurator.h>

#include <cppunit/TestCase.h>
#include <cppunit/TestSuite.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestRunner.h>

#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/CompilerOutputter.h>

#include "tests/ParticleIteratorTest.h"
#include "tests/ParticleContainerTest.h"
#include "tests/ParticleGeneratorTest.h"
#include "tests/LCParticleContainerTest.h"
#include "tests/LCInnerParticleIteratorTest.h"
#include "tests/LCOuterParticleIteratorTest.h"
#include "tests/ThermostatTest.h"
#include "tests/FileReaderTest.h"
#include "tests/MembraneTest.h"

#include <list>
#include <cassert>
#include <cstring>
#include <string>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <sys/time.h>
#include <math.h>
#include <omp.h>

using namespace std;
using namespace log4cxx;
using namespace log4cxx::xml;
using namespace log4cxx::helpers;

/**
 *  forward declaration of the simulation functions
 */
void simulate();
void LCsimulate();
void LCsimulateArgon();

/**
 * force calculation methods(with and without linked-cell)
 */
void calculateFLJ();
void LCcalculateF();

/**
 * calculate the position for all particles
 */
void calculateX();
void LCcalculateX();

/**
 * calculate the velocity for all particles
 */
void calculateV();
void LCcalculateV();

/**
 * function pointer to the used force calculation method
 */
void (*computeForce)(Particle&, Particle&);

/**
 * lennard jones force calculation
 */
void computeLennardJones(Particle& p1, Particle& p2);
void computeSmoothedLennardJones(Particle& p1, Particle& p2);
void computeMembraneLennardJones(Particle& p1, Particle& p2);
void calculateDiffusion(int para);
void calculateRDF(int para);
void getIntegerInput(string &str, int &input);
void resizeEpsSig(int inSize);
void fillEpsSig(int inSize);

// plot the particles to VTKWriter-File
void plotVTK(int iteration);
void LCplotVTK(int iteration);
void writeOutputFile(list<Particle *> parList);

void runAllTests();
void initializeSimulation(string inpFile);

double start_time = 0;
double end_time = 1000;
double delta_t = 0.014;

/**
 * contains the inputNames for simulation
 * possibilities: InputCuboids, InputParticles, InputSpheres
 */
list<string> inputNames;

/**
 * contains the inputTypes for the simulation
 * possibilities: cuboids, particles, spheres
 */
list<string> inputTypes;

string outputMask = "VTU_Archive/MD_vtk";
int freq = 10;

/**
 * For Linked Cell Algorithm
 */
double R_CUTOFF = 3.0;
utils::Vector<double, 3> domainSize;
bool outflow_flag = false;
const double h = 1.2;
int depth = 0;

/**
 * For gravity
 */
bool gravity = false;
vector<double> mass;
double G_CONST = 0; //-12.44;
/** only for cuboids and spheres, not for particles alone */
/** type functions as the index here */
double gDirMass[] = { 0.0, 1.0, 0.0 }; //will store mass (without G_CONST*)
vector<utils::Vector<double, 3> > gravForce;
vector<vector<double> > EPS;
vector<vector<double> > SIG;

/**	variables for membrane **/

double gDirMassMem[] = { 0.0, 0.0, 1.0 };
bool membraneRunning = false;
double current_time = 0.0;
double k = 300;
double rDirect = 2.2;
double rDiag = 2.2 * sqrt(2);
double rFLJ = pow(2, 1 / 6);
double FZUpArr[] = { 0.0, 0.0, 0.8 };
utils::Vector<double, 3> FZUp(FZUpArr);
int id1, id2, id3, id4;

/** other variables */

double r = 2.2;
double deltaR = 0.05;
int inputSize = 0;
list<Particle> particleList;
utils::ParticleContainer container;
utils::ParticleGenerator pgen;
utils::LCParticleContainer lcContainer;

std::vector<std::vector<double> > rdfSave;
std::vector<utils::Vector<double, 3> > diffSave;

vector<int> domainCondition;
/**
 * domainCondition contains information about
 * the conditions prevailing at the domain borders
 * it saves either 4 int values in 2D or 6 int values in 3D
 * 0: left border	1: right border
 * 2: lower border	3:upper border
 * 4: front border	5: back border
 *
 * the values stand for the following:
 * 1 stands for outflow condition
 * 2 stands for reflecting condition
 * 3 stands for periodic condition
 */

/** stores the name of the setting file */
string fileName;

/**
 * can be used to check the
 * length of the simulation
 * */
struct timeval tim;
double t1;

/** Thermostat option: only called when thermo.getEnabled() == true */
Thermostat thermo;

/** loggers */
log4cxx::LoggerPtr molsimlogger(log4cxx::Logger::getLogger("molsim"));
log4cxx::LoggerPtr argon1rdflogger(log4cxx::Logger::getLogger("argon1rdf"));
log4cxx::LoggerPtr argon1difflogger(log4cxx::Logger::getLogger("argon1diff"));
log4cxx::LoggerPtr argon2rdflogger(log4cxx::Logger::getLogger("argon2rdf"));
log4cxx::LoggerPtr argon2difflogger(log4cxx::Logger::getLogger("argon2diff"));

/**
 * @param argsv the first parameter is the choice
 * if it is empty, one can choose an experiment
 * out of a limited selection
 */
int main(int argc, char* argsv[]) {

	PropertyConfigurator::configure("Log4cxxConfig.cfg");
	LOG4CXX_INFO(molsimlogger, "Arrived @ main.");

	/* Format input command line:
	 *
	 * ./MolSim								: for manual choice
	 *
	 * ./MolSim --test						: for running tests
	 *
	 * ./MolSim --falling-drop-init			: for initialization of falling drop experiment
	 *
	 * ./MolSim --falling-drop-end			: for the falling drop experiment
	 *
	 * ./MolSim --rayleigh-taylor-small		: for the small Rayleigh-Taylor experiment
	 *
	 * ./MolSim --rayleigh-taylor-big		: for the big Rayleigh-Taylor experiment
	 *
	 * ./MolSim --rayleigh-taylor-3D		: for the Rayleigh-Taylor experiment in 3D
	 *
	 * ./MolSim --membrane			 		: for the membrane experiment
	 *
	 * ./MolSim --cooling-argon-init		: for the initialization of cooling argon experiment
	 *
	 * ./MolSim --cooling-argon-1			: for the first cooling argon experiment
	 *
	 * ./MolSim --cooling-argon-2			: for the second cooling argon experiment
	 */

	gettimeofday(&tim, NULL);
	t1 = tim.tv_sec + (tim.tv_usec / 1000000.0);
	if (argc > 2) {
		cout << "Maximal 2 arguments." << endl;
		return -1;
	}

	if (argc == 2) {
		string arg1 = argsv[1];

		if (arg1 == "--test") {
			LOG4CXX_INFO(molsimlogger, "Started testing.");

			string str;
			int option = 0;

			cout << "Here are the testing options: " << endl;
			cout
					<< "Hint: If you enter '1' or '2' you will return to this menu";
			cout << endl << endl;
			while (option != 3) {
				cout << "Enter '1' to test all unit tests." << endl;
				cout << "Enter '2' to test a specific unit test." << endl;
				cout << "Enter '3' to exit the 'test suite'." << endl;

				/** Check for correct input */
				getIntegerInput(str, option);
				if (option == 1) {
					LOG4CXX_DEBUG(molsimlogger, "Running all Tests.");
					runAllTests();
				} else if (option == 2) {
					cout << "Enter '1' to test the ParticleContainer." << endl;
					cout << "Enter '2' to test the ParticleIterator." << endl;
					cout << "Enter '3' to test the ParticleGenerator." << endl;
					cout << "Enter '4' to test the Thermostat." << endl;
					cout << "Enter '5' to test the FileReader." << endl;
					cout << "Enter '6' to test the Membrane." << endl;
					cout << "Enter '7' to test the LCParticleContainer."
							<< endl;
					cout << "Enter '8' to test the LCInnerParticleIterator."
							<< endl;
					cout << "Enter '9' to test the LCOuterParticleIterator."
							<< endl;

					/** Check for correct input */
					int opT;
					getIntegerInput(str, opT);

					CppUnit::TextUi::TestRunner runner;
					switch (opT) {
					case 1:
						LOG4CXX_DEBUG(molsimlogger,
								"Testing ParticleContainer alone.")
						;
						runner.addTest(ParticleContainerTest::suite());
						break;
					case 2:
						LOG4CXX_DEBUG(molsimlogger,
								"Testing ParticleIterator alone.")
						;
						runner.addTest(ParticleIteratorTest::suite());
						break;
					case 3:
						LOG4CXX_DEBUG(molsimlogger,
								"Testing ParticleGenerator alone.")
						;
						runner.addTest(ParticleGeneratorTest::suite());
						break;
					case 4:
						LOG4CXX_DEBUG(molsimlogger, "Testing Thermostat alone.")
						;
						runner.addTest(ThermostatTest::suite());
						break;
					case 5:
						LOG4CXX_DEBUG(molsimlogger, "Testing FileReader alone.")
						;
						runner.addTest(FileReaderTest::suite());
						break;
					case 6:
						LOG4CXX_DEBUG(molsimlogger, "Testing Membrane alone.")
						;
						runner.addTest(MembraneTest::suite());
						break;

					case 7:
						LOG4CXX_DEBUG(molsimlogger,
								"Testing LCParticleContainer alone.")
						;
						runner.addTest(LCParticleContainerTest::suite());
						break;
					case 8:
						LOG4CXX_DEBUG(molsimlogger,
								"Testing LCInnerParticleIterator alone.")
						;
						runner.addTest(LCInnerParticleIteratorTest::suite());
						break;
					case 9:
						LOG4CXX_DEBUG(molsimlogger,
								"Testing LCOuterParticleIterator alone.")
						;
						runner.addTest(LCOuterParticleIteratorTest::suite());
						break;
					}
					runner.run();
				} else if (option == 3) {
					LOG4CXX_INFO(molsimlogger, "Stopped testing.");
					cout << "Alright! Goodbye then!" << endl;
				} else {
					LOG4CXX_INFO(molsimlogger, "Aborted testing.");
					cout << "Invalid number, please try again.\n" << endl;
				}
			}
		} else if (arg1 == "--falling-drop-init") {
			LOG4CXX_INFO(molsimlogger, "Doing FallingDropInit.");
			string input = "input/FallingDropInitSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--falling-drop-end") {
			LOG4CXX_INFO(molsimlogger, "Doing FallingDropEnd.");
			string input = "input/FallingDropEndSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--rayleigh-taylor-small") {
			LOG4CXX_INFO(molsimlogger, "Doing RayleighTaylorSmall.");
			string input = "input/RayleighTaylorSmallSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--rayleigh-taylor-big") {
			LOG4CXX_INFO(molsimlogger, "Doing RayleighTaylorBig.");
			string input = "input/RayleighTaylorBigSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--rayleigh-taylor-3D") {
			LOG4CXX_INFO(molsimlogger, "Doing RayleighTaylor3D.");
			string input = "input/RayleighTaylor3DSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--membrane") {
			string input = "input/MembraneSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--cooling-argon-init") {
			LOG4CXX_INFO(molsimlogger, "Doing CoolingArgonInit.");
			string input = "input/CoolingArgonInitSetting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--cooling-argon-1") {
			LOG4CXX_INFO(molsimlogger, "Doing CoolingArgon1.");
			string input = "input/CoolingArgon1Setting.xml";
			initializeSimulation(input);
		} else if (arg1 == "--cooling-argon-2") {
			LOG4CXX_INFO(molsimlogger, "Doing CoolingArgon2.");
			string input = "input/CoolingArgon2Setting.xml";
			initializeSimulation(input);
		} else {
			LOG4CXX_ERROR(molsimlogger, "Invalid choice input.");
			cout << arg1 << " invalid!" << endl;
			return -1;
		}
	} else {
		/**
		 * argc==1
		 * ./MolSim is being executed
		 */
		LOG4CXX_INFO(molsimlogger, "Manual file decision.");
		string str;
		FileReader fileReader;
		cout << "Hello from Andreas, Matthias and Son!" << endl;
		cout << endl;
		/** variable for choice */
		int option1;
		cout << "Enter '1' to run the program with a particle file." << endl;
		cout << "Enter '2' to run program with a cuboid file." << endl;
		cout << "Enter '3' to run program with xml input files." << endl;

		getIntegerInput(str, option1);

		switch (option1) {
		case 1:
			fileName = "eingabe-sonne.txt";
			outputMask = "VTU_Archive/MolSimSonne_vtk";
			cout << fileName << " will be used.\n" << endl;
			delta_t = 0.014;
			end_time = 1000;
			break;
		case 2:
			fileName = "eingabe-brownian.txt";
			outputMask = "VTU_Archive/MolSimBrownian_vtk";
			cout << fileName << " will be used.\n" << endl;
			delta_t = 0.0002;
			end_time = 5;
			break;
		case 3:
			cout << "XML input files are stored in MolSim." << endl;
			cout << "There are 3 sources of input files:" << endl;
			cout
					<< "\tInputSetting: contains start_time, end_time, delta_t,\
					\n\t\t inputfile name, inputfile type, output mask\
					\n\t\t and output frequency."
					<< endl;
			cout
					<< "\tInputParticles: contains all information needed for particles."
					<< endl;
			cout
					<< "\tInputCuboids: contains all information needed for cuboids."
					<< endl;
			cout
					<< "\tInputSpheres: contains all information needed for spheres.\n"
					<< endl;
			outputMask = "VTU_Archive/MolSimXML_vtk";
			break;

		default:
			cout << "Invalid number!" << endl;
			return -1;
			break;
		}
		if (option1 == 3) {
			/** getting information from InputSetting first */
			fileName = "input/InputSetting.xml";
			pgen.extractSetting(fileName, start_time, end_time, delta_t,
					inputNames, inputTypes, outputMask, freq, domainSize,
					R_CUTOFF, domainCondition, G_CONST, inputSize);
			particleList.clear();
			list<string>::iterator itT = inputTypes.begin();
			int i = 1;

			/** initialize the size of gravForce */
			gravForce.resize(inputSize);
			resizeEpsSig(inputSize);

			for (list<string>::iterator itN = inputNames.begin();
					itN != inputNames.end(); itN++) {
				if (*itT == "particles") {
					cout << "Input file #" << i << ": " << "[particles]."
							<< endl;
					pgen.extractParticles(*itN);
					for (list<Particle>::iterator it =
							pgen.getParticleList().begin();
							it != pgen.getParticleList().end(); it++) {
						gDirMass[1] = (*it).getM();
						gravForce[(*it).getType()] = utils::Vector<double, 3>(
								gDirMass);
					}
					pgen.mergeWithParticleList(particleList);
				} else if (*itT == "cuboids") {
					cout << "Input file #" << i << ": " << "[cuboids]" << endl;
					pgen.extractCuboids(*itN);
					/** For each type */
					for (list<Cuboid>::iterator it =
							pgen.getCuboidList().begin();
							it != pgen.getCuboidList().end(); it++) {
						gDirMass[1] = (*it).getMass();
						gravForce[(*it).getType()] = utils::Vector<double, 3>(
								gDirMass);
						/** fill the diagonal first */
						EPS[(*it).getType()][(*it).getType()] =
								(*it).getEpsilon();
						SIG[(*it).getType()][(*it).getType()] =
								(*it).getSigma();
					}

					pgen.cuboidsToList();
					pgen.mergeWithParticleList(particleList);

				} else {
					cout << "Input file #" << i << ": " << "[spheres]" << endl;
					pgen.extractSpheres(*itN);
					/** For each type */
					for (list<Sphere>::iterator it =
							pgen.getSphereList().begin();
							it != pgen.getSphereList().end(); it++) {
						gDirMass[1] = (*it).getM();
						gravForce[(*it).getType()] = utils::Vector<double, 3>(
								gDirMass);
						/** fill the diagonal first */
						EPS[(*it).getType()][(*it).getType()] =
								(*it).getEpsilon();
						SIG[(*it).getType()][(*it).getType()] =
								(*it).getSigma();
					}
					pgen.spheresToList();
					pgen.mergeWithParticleList(particleList);
				}
				itT++;
				i++;
			}
			//TODO: cout << "\n";
			/**======================GRAVITY=====================*/
			int ag;
			cout
					<< "Do you want to add gravity?\nEnter '1' to confirm, any other number to decline."
					<< endl;
			getIntegerInput(str, ag);
			if (ag == 1) {
				cout << "Gravity enabled.\n" << "G = " << G_CONST << "."
						<< endl;
			} else {
				cout << "Gravity disabled." << endl;
				G_CONST = 0;
			}
			/** set the gravity for each type */
			for (int i = 0; i < gravForce.size(); i++) {
				gravForce[i] = G_CONST * gravForce[i];
			}
			/**======================GRAVITY=====================*/

			fillEpsSig(inputSize);
		}
		/**======================SIMULATION=====================*/
		//TODO: cout << "Running simulation..." << endl;
		int wo;
		if (option1 == 3) {
			lcContainer.initialize(particleList, domainSize, R_CUTOFF);
			LCsimulate();
			cout << "\nDo you want to write particles out?" << endl;
			cout << "Enter '1' to confirm, any other number to decline."
					<< endl;
			getIntegerInput(str, wo);
			if (wo == 1) {
				writeOutputFile((lcContainer.getList()));
				cout << "ParListStatus.txt written." << endl;
			}
		} else {
			char *cstr = new char[fileName.length() + 1];
			strcpy(cstr, fileName.c_str());
			pgen.readCuboids(cstr);
			pgen.cuboidsToList();
			particleList = pgen.getParticleList();

			resizeEpsSig(1);
			EPS[0][0] = 5;
			SIG[0][0] = 1;
			container.initialize(particleList);

			simulate();
		}
		/**======================SIMULATION=====================*/
	}

	LOG4CXX_INFO(molsimlogger, "Simulation ended successfully.");

	return 0;
}

void runAllTests() {
	CppUnit::TextUi::TestRunner runner;
	runner.addTest(ParticleIteratorTest::suite());
	runner.addTest(ParticleContainerTest::suite());
	runner.addTest(ParticleGeneratorTest::suite());
	runner.addTest(MembraneTest::suite());
	runner.addTest(LCParticleContainerTest::suite());
	runner.addTest(LCInnerParticleIteratorTest::suite());
	runner.addTest(LCOuterParticleIteratorTest::suite());
	runner.addTest(ThermostatTest::suite());
	runner.addTest(FileReaderTest::suite());
	runner.run();
}

void initializeSimulation(string inpFile) {
	LOG4CXX_INFO(molsimlogger, "Arrived @ initialization call.");
	fileName = inpFile;
	if (fileName == "input/FallingDropEndSetting.xml"
			|| fileName == "input/CoolingArgon1Setting.xml"
			|| fileName == "input/CoolingArgon2Setting.xml") {
		pgen.getParticleList().clear();
		double eps1 = 1.0;
		double sig1 = 1.0;
		FileReader fileReader;
		string inTextS = "input/ParListArgon.txt";
		if (fileName == "input/FallingDropEndSetting.xml") {
			inTextS = "input/ParListStatus.txt";
			sig1 = 1.2;
		}
		char *inText = new char[inTextS.length() + 1];
		strcpy(inText, inTextS.c_str());
		/** extract particles from the particle file */
		fileReader.readStatus(pgen.getParticleList(), eps1, sig1, inText);
		particleList = pgen.getParticleList();
		particleList.size();
	}
	membraneRunning = fileName == "input/MembraneSetting.xml";
	pgen.extractSetting(fileName, start_time, end_time, delta_t, inputNames,
			inputTypes, outputMask, freq, domainSize, R_CUTOFF, domainCondition,
			G_CONST, inputSize);

	if (fileName == "input/FallingDropEndSetting.xml") {
		/** Gravity + EpsSig */
		gravForce.resize(2);
		resizeEpsSig(2);
		gDirMass[1] = G_CONST * ((*particleList.begin()).getM());
		gravForce[0] = utils::Vector<double, 3>(gDirMass);
		EPS[0][0] = 1.0;
		SIG[0][0] = 1.2;
		pgen.extractSpheres(*inputNames.begin());
		list<Sphere>::iterator itS = pgen.getSphereList().begin();
		gDirMass[1] = G_CONST * ((*itS).getM());
		gravForce[1] = utils::Vector<double, 3>(gDirMass);
		EPS[1][1] = (*itS).getEpsilon();
		SIG[1][1] = (*itS).getSigma();
		fillEpsSig(2);

		pgen.spheresToList();
		pgen.mergeWithParticleList(particleList);

	} else if (fileName == "input/FallingDropInitSetting.xml"
			|| fileName == "input/CoolingArgonInitSetting.xml") {
		particleList.clear();
		/** Gravity + EpsSig */
		gravForce.resize(1);
		resizeEpsSig(1);
		pgen.extractCuboids(*inputNames.begin());
		list<Cuboid>::iterator itC = pgen.getCuboidList().begin();
		gDirMass[1] = G_CONST * ((*itC).getMass());
		gravForce[0] = utils::Vector<double, 3>(gDirMass);
		EPS[0][0] = (*itC).getEpsilon();
		SIG[0][0] = (*itC).getSigma();
		if (fileName == "input/FallingDropInitSetting.xml"
				|| fileName == "input/CoolingArgonInitSetting.xml") {
			pgen.cuboidsToList();
			particleList = pgen.getParticleList();
		}

	} else if (fileName == "input/CoolingArgon1Setting.xml"
			|| fileName == "input/CoolingArgon2Setting.xml") {
		std::list<Particle>::iterator iterator = particleList.begin();
		/** Gravity + EpsSig */
		gravForce.resize(1);
		resizeEpsSig(1);
		pgen.extractCuboids(*inputNames.begin());
		list<Cuboid>::iterator itC = pgen.getCuboidList().begin();
		gDirMass[1] = G_CONST * ((*itC).getMass());
		gravForce[0] = utils::Vector<double, 3>(gDirMass);
		EPS[0][0] = (*itC).getEpsilon();
		SIG[0][0] = (*itC).getSigma();

	} else if (fileName == "input/MembraneSetting.xml") {
		membraneRunning = true;
		pgen.extractSetting(fileName, start_time, end_time, delta_t, inputNames,
				inputTypes, outputMask, freq, domainSize, R_CUTOFF,
				domainCondition, G_CONST, inputSize);
		particleList.clear();
		gravForce.resize(1);
		resizeEpsSig(1);
		pgen.extractCuboids(*inputNames.begin());
		list<Cuboid>::iterator itC = pgen.getCuboidList().begin();
		(*itC).initNeighbours();

		/** Gravity + EpsSig */
		gDirMassMem[0] = 0;
		gDirMassMem[1] = 0;
		gDirMassMem[2] = -0.001;
		gravForce[0] = utils::Vector<double, 3>(gDirMassMem);
		EPS[0][0] = (*itC).getEpsilon();
		SIG[0][0] = (*itC).getSigma();

		/** set the IDs of 4 special particles */
		id1 = 24 * ((*itC).getWidth()) + 17;        //(17, 24)
		id2 = 25 * ((*itC).getWidth()) + 17;        //(17, 25)
		id3 = 24 * ((*itC).getWidth()) + 18;        //(18, 24)
		id4 = 25 * ((*itC).getWidth()) + 18;        //(18, 25)

		pgen.cuboidsToList();
		particleList = pgen.getParticleList();
	} else {

		particleList.clear();
		/** Gravity + EpsSig */
		gravForce.resize(inputSize);
		resizeEpsSig(inputSize);
		pgen.extractCuboids(*inputNames.begin());
		int t;
		for (list<Cuboid>::iterator itC = pgen.getCuboidList().begin();
				itC != pgen.getCuboidList().end(); itC++) {
			gDirMass[1] = G_CONST * ((*itC).getMass());
			t = (*itC).getType();
			gravForce[t] = utils::Vector<double, 3>(gDirMass);
			EPS[t][t] = (*itC).getEpsilon();
			SIG[t][t] = (*itC).getSigma();
		}
		fillEpsSig(inputSize);
		pgen.cuboidsToList();
		particleList = pgen.getParticleList();
	}
	thermo = Thermostat(fileName);

	lcContainer.initialize(particleList, domainSize, R_CUTOFF);
	LOG4CXX_INFO(molsimlogger, "Arrived @ simulation call.");
	if (fileName == "input/CoolingArgonInitSetting.xml"
			|| fileName == "input/CoolingArgon1Setting.xml"
			|| fileName == "input/CoolingArgon2Setting.xml") {
		LCsimulateArgon();
	} else {
		LCsimulate();
	}
	if (fileName == "input/FallingDropInitSetting.xml"
			|| fileName == "input/CoolingArgonInitSetting.xml") {
		writeOutputFile((lcContainer.getList()));
	}
}

void simulate() {
	LOG4CXX_INFO(molsimlogger,
			"Size of container: " << container.size() << " particles.");
	LOG4CXX_DEBUG(molsimlogger,
			"Starting force calculation for the first time...");
	calculateFLJ();

	double current_time = start_time;

	int iteration = 0;

	double temperature = thermo.getT_init();
	bool target_temp_reached = false;

	/** for this loop, we assume: current x, current f and current v are set */
	while (current_time < end_time) {
		/** const clock_t beginTime = clock(); */

		calculateX();

		calculateFLJ();

		calculateV();

		iteration++;
		cout << "\r" << "Iteration " << iteration << " completed." << flush;

		if (iteration % freq == 0) {
			plotVTK(iteration);
		}
		LOG4CXX_TRACE(molsimlogger, "Iteration " << iteration << " finished.");

		current_time += delta_t;

		/**
		 * cout << float( clock () - beginTime ) /  CLOCKS_PER_SEC;
		 * cin.ignore();
		 */
	}
	cout << "\nOutput written. Terminating..." << endl;
}

/**
 * This method calculates the forces between the particles.
 * The calculation obeys the Lennard-Jones force between two molecules.
 */
void calculateFLJ() {

	utils::Vector<double, 3> zero((double) 0);
	utils::ParticleIterator iterator;
	utils::Vector<double, 3> sumF[container.size()];
	for (int i = 0; i < container.size(); i++) {
		sumF[i] = zero;
	}
	iterator = container.begin();
	int i = 0;
	while (iterator != container.end()) {
		utils::ParticleIterator innerIterator;
		innerIterator = iterator;
		++innerIterator;
		int j = i + 1;

		while (innerIterator != container.end()) {

			Particle& p1 = *iterator;
			Particle& p2 = *innerIterator;

			//calculations
			utils::Vector<double, 3> tempD = p2.getX() - p1.getX();
			double tempDNorm = tempD.L2Norm();

			double tempDSigDivNormPowSix = pow(SIG[0][0] / tempDNorm, 6);
			utils::Vector<double, 3> tempF =
					24 * EPS[0][0] * pow(1 / tempDNorm, 2)
							* (tempDSigDivNormPowSix
									- 2 * pow(tempDSigDivNormPowSix, 2))
							* tempD;

			sumF[i] += tempF;
			sumF[j] += (-1) * tempF;

			++innerIterator;
			++j;
		}
		++iterator;
		++i;
	}
}

/**
 *  This method calculates the position of the particles.
 *  It obeys the Velocity-Stoermer-Verlet-Algorithm.
 */
void calculateX() {
	utils::ParticleIterator iterator;
	iterator = container.begin();
	while (iterator != container.end()) {
		Particle& p = *iterator;
		utils::Vector<double, 3> tempX = p.getX() + delta_t * p.getV()
				+ ((delta_t) * (delta_t) / (2 * p.getM())) * p.getOldF();
		p.getX() = tempX;
		++iterator;
	}
}

/**
 *  This method calculates the position of the particles.
 *  It obeys the Velocity-Stoermer-Verlet-Algorithm.
 */
void calculateV() {
	utils::ParticleIterator iterator;
	iterator = container.begin();
	while (iterator != container.end()) {
		Particle& p = *iterator;
		utils::Vector<double, 3> tempV = p.getV()
				+ (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
		p.getV() = tempV;
		++iterator;
	}
}

/**
 * this method catches user input for experiment choice
 */
void getIntegerInput(string &str, int &input) {
	while (true) {
		getline(cin, str);
		stringstream myStream(str);
		if (myStream >> input)
			break;
		cout << "Invalid number, please try again" << endl;
	}
}

/**
 * this method writes the output VTK files
 */
void plotVTK(int iteration) {

	LOG4CXX_TRACE(molsimlogger, "Arrived @ plotVTK.");

	outputWriter::VTKWriter writer;
	utils::ParticleIterator iterator;
	iterator = container.begin();
	writer.initializeOutput(container.size());
	while (iterator != container.end()) {
		Particle& p = *iterator;
		writer.plotParticle(p);
		++iterator;
	}
	string out_name(outputMask);
	writer.writeFile(out_name, iteration);
}

/** LINKED CELL ALGORITHMS */

/*! To compare the execution time between naive and LC Algorithm:
 *  \image html GraphLCWithoutPlot.jpg
 *  \image html GraphLCWithPlot.jpg
 */

void LCsimulate() {
	LOG4CXX_INFO(molsimlogger, "lcContainer.size: " << lcContainer.size());
	LOG4CXX_INFO(molsimlogger, "particleList.size: " << particleList.size());

	/**
	 * set the function pointer to the normal
	 * LJ formula or to the smoothed version
	 */
	computeForce = computeLennardJones;
	if (membraneRunning) {
		computeForce = computeMembraneLennardJones;
	}
	utils::BoundaryHandler boundHandler(domainCondition, lcContainer, h,
			computeLennardJones);
	LCcalculateF();
	current_time = start_time;
	double temperature = thermo.getT_init();
	bool target_temp_reached = false;

	int iteration = 0;
	/** for this loop, we assume: current x, current f and current v are known */
	while (current_time < end_time) {
		/** const clock_t beginTime = clock(); */
		LCcalculateX();
		boundHandler.applyOutflow();
		boundHandler.applyPeriodicMoving();
		lcContainer.updateCells();
		boundHandler.applyReflecting();
		boundHandler.applyPeriodic();
		LCcalculateF();
		LCcalculateV();
		iteration++;
		cout << "\r" << "Iteration " << iteration << " completed." << flush;
		if (thermo.getEnabled()) {
			if (!target_temp_reached) {
				if (iteration % thermo.getn_delta() == 0) {
					temperature += thermo.getDelta_T();
					if (temperature > thermo.getT_target()) {
						temperature -= thermo.getDelta_T();
						target_temp_reached = true;
					}
				}
			}
			if (iteration % thermo.getn_thermo() == 0) {
				thermo.setThermo(lcContainer.getList(), 2, temperature);
			}
		}
		if (iteration % freq == 0) {
			LCplotVTK(iteration);
		}

		LOG4CXX_TRACE(molsimlogger, "Iteration " << iteration << " finished.");

		/**
		 * if(iteration == 1000){
		 * 	gettimeofday(&tim, NULL);
		 * 	double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
		 * 	cout << (t2-t1) << endl;
		 * }
		 */
		current_time += delta_t;

		/**
		 * cout << float( clock () - beginTime ) /  CLOCKS_PER_SEC;
		 * cin.ignore();
		 */
	}

	cout << "\nOutput written. Terminating..." << endl;
}

void LCsimulateArgon() {
	computeForce = computeSmoothedLennardJones;
	utils::BoundaryHandler boundHandler(domainCondition, lcContainer, h,
			computeLennardJones);
	LCcalculateF();
	current_time = start_time;
	rdfSave.resize((end_time / delta_t) + 1);
	for (int tem = 0; tem < rdfSave.size(); tem++) {
		int z = r / deltaR;
		rdfSave[tem].resize(z);
		for (int temm = 0; temm < z; temm++) {
			rdfSave[tem][temm] = 0;
		}
	}
	diffSave.resize((end_time / delta_t) + 3);
	utils::Vector<double, 3> zero((double) 0);
	for (int tem = 0; tem < (end_time / delta_t) + 2; tem++) {
		diffSave[tem] = zero;
	}

	double temperature = thermo.getT_init();
	bool target_temp_reached = false;

	int iteration = 0;
	while (current_time < end_time) {
		LCcalculateX();
		boundHandler.applyOutflow();
		boundHandler.applyPeriodicMoving();
		lcContainer.updateCells();
		boundHandler.applyReflecting();
		boundHandler.applyPeriodic();
		LCcalculateF();
		LCcalculateV();

		/** calculating RDF and Diffusion */
		if ((iteration % 1000) < 5 && iteration > 10) {
			int i = iteration / 1000;

			calculateRDF(i);
			calculateDiffusion(i);
			if ((iteration % 1000) == 4) {
				double ri = 0;
				diffSave[i] = diffSave[i] * 0.2;
				double diffSum = 0;
				for (int o = 0; o < 3; o++) {
					diffSum += diffSave[i][o];
				}
				diffSum = diffSum * (1.0 / 3.0);
				if (fileName == "input/CoolingArgon1Setting.xml") {
					LOG4CXX_INFO(argon1rdflogger,
							"Iteration: " << iteration << "----------");
					LOG4CXX_INFO(argon1difflogger,
							"Iteration " << iteration << ": " << diffSum);
					for (int o = 0; o < r / deltaR; o++) {
						ri = o * deltaR;
						LOG4CXX_INFO(argon1rdflogger,
								"]" << ri << ", " << (ri + deltaR) << "]: " << (rdfSave[i][o] * 0.2));
					}
				} else if (fileName == "input/CoolingArgon2Setting.xml") {
					LOG4CXX_INFO(argon2rdflogger,
							"Iteration: " << iteration << "----------");
					LOG4CXX_INFO(argon2difflogger,
							"Iteration " << iteration << ": " << diffSum);
					for (int o = 0; o < r / deltaR; o++) {
						ri = o * deltaR;
						LOG4CXX_INFO(argon2rdflogger,
								"]" << ri << ", " << (ri + deltaR) << "]: " << (rdfSave[i][o] * 0.2));
					}
				}
			}
		}

		if (iteration % freq == 0) {
			LCplotVTK(iteration);
		}

		iteration++;
		cout << "\r" << "Iteration " << iteration << " completed." << flush;

		if (thermo.getEnabled()) {
			if (!target_temp_reached) {
				if (iteration % thermo.getn_delta() == 0) {
					temperature += thermo.getDelta_T();
					if (temperature > thermo.getT_target()) {
						temperature -= thermo.getDelta_T();
						target_temp_reached = true;
					}
				}
			}
			if (iteration % thermo.getn_thermo() == 0) {
				thermo.setThermo(lcContainer.getList(), 2, temperature);
			}
		}
		LOG4CXX_TRACE(molsimlogger, "Iteration " << iteration << " finished.");
		current_time += delta_t;
	}
	cout << "\nOutput written. Terminating..." << endl;
}

void LCcalculateF() {

	utils::Vector<double, 3> zero((double) 0);
	utils::Vector<double, 3> sumF((double) 0);

#ifndef _OPENMP
	utils::LCOuterParticleIterator iterator = lcContainer.beginOuter();
	while (iterator != lcContainer.endOuter()) {
		sumF = zero;
		utils::LCInnerParticleIterator innerIterator = lcContainer.beginInner(
				iterator);
		while (innerIterator != lcContainer.endInner(iterator.getCellNumber())) {
			Particle& p1 = *iterator;
			Particle& p2 = *innerIterator;
			if (p1 == p2) {
				++innerIterator;
				continue;
			} else {
				assert(!(p1 == p2));
				computeForce(p1, p2);
			}
			++innerIterator;
		}
		++iterator;
	}
#endif
#ifdef _OPENMP
	utils::LCOuterParticleIterator iterator;
	utils::LCInnerParticleIterator innerIterator;
	int i = 0;
#pragma omp parallel private(i, iterator, innerIterator, zero, sumF)
	{
#pragma omp for schedule(dynamic)
		for (i = 0; i < omp_get_num_threads(); i++) {
			iterator = lcContainer.beginOuter(i);
			while (iterator != lcContainer.endOuter(i)) {
				assert(iterator.getCellNumber() <= lcContainer.endOuter(i).getCellNumber());
				sumF = zero;
				innerIterator = lcContainer.beginInner(iterator);
				while (innerIterator
						!= lcContainer.endInner(iterator.getCellNumber())) {
					Particle& p1 = *iterator;
					Particle& p2 = *innerIterator;
					if (p1 == p2) {
						assert( p1 == p2);
						++innerIterator;
						continue;
					} else {
						computeForce(p1, p2);
					}
					++innerIterator;
				}
				++iterator;
			}

		}
	}
#endif

	/** setF */
	utils::LCOuterParticleIterator iterator2 = lcContainer.beginOuter();
	while (iterator2 != lcContainer.endOuter()) {
		bool set = false;
		if (membraneRunning) {
			if (current_time <= 150) {
				if (((*iterator2).getID() == id1)
						|| ((*iterator2).getID() == id2)
						|| ((*iterator2).getID() == id3)
						|| ((*iterator2).getID() == id4)) {
					(*iterator2).setF(
							gravForce[(*iterator2).getType()]
									+ (*iterator2).getTempF() + FZUp);
					(*iterator2).deleteTempF();
					set = true;
				}
			}
		}
		if (set == false) {
			(*iterator2).setF(
					gravForce[(*iterator2).getType()]
							+ (*iterator2).getTempF());
			(*iterator2).deleteTempF();
		}
		++iterator2;
	}
}

void LCcalculateX() {
#ifndef _OPENMP
	utils::LCOuterParticleIterator iterator = lcContainer.beginOuter();
	while (iterator != lcContainer.endOuter()) {
		Particle& p = *iterator;
		utils::Vector<double, 3> tempX = p.getX() + delta_t * p.getV()
				+ ((delta_t) * (delta_t) / (2 * p.getM())) * p.getOldF();
		p.getX() = tempX;
		++iterator;
	}
#else
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for(int i = 0; i < omp_get_num_threads(); i++) {
			utils::LCOuterParticleIterator iterator = lcContainer.beginOuter(i);
			while (iterator != lcContainer.endOuter(i)) {
				Particle& p = *iterator;
				utils::Vector<double, 3> tempX = p.getX() + delta_t * p.getV()
				+ ((delta_t) * (delta_t) / (2 * p.getM())) * p.getOldF();
				p.getX() = tempX;
				++iterator;
			}
		}
	}
#endif
}

void LCcalculateV() {
#ifndef _OPENMP
	utils::LCOuterParticleIterator iterator = lcContainer.beginOuter();
	while (iterator != lcContainer.endOuter()) {
		Particle& p = *iterator;
		utils::Vector<double, 3> tempV = p.getV()
				+ (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
		p.getV() = tempV;
		++iterator;
	}
#else
#pragma omp parallel
	{
#pragma omp for schedule(dynamic)
		for (int i = 0; i < omp_get_num_threads(); i++) {
			utils::LCOuterParticleIterator iterator = lcContainer.beginOuter(i);
			while (iterator != lcContainer.endOuter(i)) {
				Particle& p = *iterator;
				utils::Vector<double, 3> tempV = p.getV()
				+ (delta_t / (2 * p.getM())) * (p.getF() + p.getOldF());
				p.getV() = tempV;
				++iterator;
			}
		}
	}
#endif
}

/**
 * @param p1 the particle of the outer iterator
 * @param p2 the particle of the inner iterator
 */
void computeLennardJones(Particle& p1, Particle& p2) {
	utils::Vector<double, 3> tempD = p1.getX() - p2.getX();
	double tempDNorm = tempD.L2Norm();

	if (tempDNorm < R_CUTOFF) {

		double tempDSigDivNorm = pow(
				SIG[p1.getType()][p2.getType()] / tempDNorm, 6);

		utils::Vector<double, 3> tempF = 24 * EPS[p1.getType()][p2.getType()]
				* pow(1 / tempDNorm, 2)
				* (tempDSigDivNorm - 2 * pow(tempDSigDivNorm, 2))
				* (p2.getX() - p1.getX());
		//TODO: tempD * (-1)

		p2.updateTempF((-1) * tempF);
		p1.updateTempF(tempF);
	}
}

/**
 *
 * @param p1 the particle of the outer iterator
 * @param p2 the particle of the inner iterator
 */
void computeMembraneLennardJones(Particle& p1, Particle& p2) {
	utils::Vector<double, 3> tempD = p2.getX() - p1.getX();
	double tempDNorm = tempD.L2Norm();
	utils::Vector<double, 3> tempF((double) 0);

	if (tempDNorm <= rFLJ) {
		double tempDSigDivNorm = pow(
				SIG[p1.getType()][p2.getType()] / tempDNorm, 6);
		tempF = 24 * EPS[p1.getType()][p2.getType()] * pow(1 / tempDNorm, 2)
				* (tempDSigDivNorm - 2 * pow(tempDSigDivNorm, 2)) * tempD;
	}

	if (p1.isDirectNeighbour(p2)) {
		tempF = tempF + k * (tempDNorm - rDirect) * (1 / (tempDNorm)) * tempD;
	} else if (p1.isDiagNeighbour(p2)) {
		tempF = tempF + k * (tempDNorm - rDiag) * (1 / (tempDNorm)) * tempD;
	}

	p2.updateTempF((-1) * tempF);
	p1.updateTempF(tempF);

}

void computeSmoothedLennardJones(Particle& p1, Particle& p2) {
	utils::Vector<double, 3> tempD = p1.getX() - p2.getX();
	double tempDNorm = tempD.L2Norm();
	assert(tempDNorm > 0);

	double rl = 1.9;
	utils::Vector<double, 3> tempS;
	if (tempDNorm <= rl) {
		tempS = 1.0;
	} else if (tempDNorm <= R_CUTOFF) {
		tempS = 1
				- (tempDNorm - rl) * (tempDNorm - rl)
						* ((3 * R_CUTOFF - rl - (2 * tempDNorm))
								/ pow(R_CUTOFF - rl, 3));
	} else {
		tempS = 0.0;

	}

	double tempDSigDivNorm = pow(SIG[p1.getType()][p2.getType()] / tempDNorm,
			6);

	utils::Vector<double, 3> tempF = 24 * EPS[p1.getType()][p2.getType()]
			* pow(1 / tempDNorm, 2)
			* (tempDSigDivNorm - 2 * pow(tempDSigDivNorm, 2)) * tempS
			* (p2.getX() - p1.getX());
	//TODO: tempD * (-1)

	p2.updateTempF((-1) * tempF);
	p1.updateTempF(tempF);
}

void LCplotVTK(int iteration) {

	outputWriter::VTKWriter writer;
	utils::LCOuterParticleIterator iterator = lcContainer.beginOuter();
	writer.initializeOutput(lcContainer.size());
	while (iterator != lcContainer.endOuter()) {
		Particle& p = *iterator;
		writer.plotParticle(p);
		++iterator;
	}
	string out_name(outputMask);
	writer.writeFile(out_name, iteration);
}

void writeOutputFile(list<Particle *> parList) {
	ofstream file;
	if (fileName == "input/CoolingArgonInitSetting.xml") {
		file.open("input/ParListArgon.txt", ios::trunc);
	} else {
		file.open("input/ParListStatus.txt", ios::trunc);
	}
	file << "# file format:\n"
			<< "# Lines of comment start with '#' and are only allowed at the beginning of the file\n"
			<< "# Empty lines are not allowed.\n"
			<< "# The first line not being a comment has to be "
			<< "# <int: number of particles> <double: epsilon> <double: sigma>\n"
			<< "# molecule data sets.\n" << "#\n"
			<< "# Molecule data consists of\n"
			<< "# * xyz-coordinates (3 double values)\n"
			<< "# * velocities (3 double values)\n"
			<< "# * force (3 double values)\n"
			<< "# * old force (3 double values)\n"
			<< "# * mass (1 double value)\n" << "# * type (1 int value)\n"
			<< "#\n" << "# " << setw(45) << "xyz-coord" << setw(45)
			<< "velocity" << setw(45) << "force" << setw(45) << "old force"
			<< setw(15) << "mass" << setw(10) << "type\n" << setw(10)
			<< parList.size() << setw(10)
			<< EPS[(*parList.begin())->getType()][(*parList.begin())->getType()]
			<< setw(10)
			<< SIG[(*parList.begin())->getType()][(*parList.begin())->getType()]
			<< endl;
	for (list<Particle *>::iterator it = parList.begin(); it != parList.end();
			it++) {
		file << setw(15) << (*it)->getX()[0] << setw(15) << (*it)->getX()[1]
				<< setw(15) << (*it)->getX()[2]

				<< setw(15) << (*it)->getV()[0] << setw(15) << (*it)->getV()[1]
				<< setw(15) << (*it)->getV()[2]

				<< setw(15) << (*it)->getF()[0] << setw(15) << (*it)->getF()[1]
				<< setw(15) << (*it)->getF()[2]

				<< setw(15) << (*it)->getOldF()[0] << setw(15)
				<< (*it)->getOldF()[1] << setw(15) << (*it)->getOldF()[2]

				<< setw(15) << (*it)->getM() << setw(10) << (*it)->getType()
				<< endl;
	}
	file.close();
}

void resizeEpsSig(int inSize) {
	EPS.resize(inSize);
	SIG.resize(inSize);
	for (int i = 0; i < inSize; i++) {
		EPS[i].resize(inSize);
		SIG[i].resize(inSize);
	}
}

/** the diagonal of the 2D matrix must have been filled before */
void fillEpsSig(int inSize) {
	for (int i = 0; i < inSize; i++) {
		for (int j = 0; j < inSize; j++) {
			/** Lorentz-Berthelot Mixing rule */
			EPS[i][j] = (EPS[i][i] + EPS[j][j]) / 2;
			SIG[i][j] = sqrt(SIG[i][i] * SIG[j][j]);
		}
	}
}

void calculateDiffusion(int para) {
	utils::Vector<double, 3> diffusion = 0.0;
	utils::LCOuterParticleIterator iterator = lcContainer.beginOuter();
	int i = 0; /** number of particles */
	while (iterator != lcContainer.endOuter()) {
		Particle& p = *iterator;
		diffusion += (p.getX() - p.getX0()) * (p.getX() - p.getX0());
		++i;
		++iterator;
	}
	diffusion = (1 / (double) i) * diffusion;
	diffSave[para] += diffusion;
}

void calculateRDF(int para) {
	std::vector<double> number;
	number.resize((r / deltaR));
	for (int i = 0; i < number.size(); i++) {
		number[i] = 0;
	}
	utils::LCOuterParticleIterator outerIterator = lcContainer.beginOuter();
	while (outerIterator != lcContainer.endOuter()) {
		utils::LCInnerParticleIterator innerIterator = lcContainer.beginInner(
				outerIterator);
		while (innerIterator
				!= lcContainer.endInner(outerIterator.getCellNumber())) {
			utils::Vector<double, 3> tempD = (*innerIterator).getX()
					- (*outerIterator).getX();
			double tempDNorm = tempD.L2Norm();
			if (tempDNorm <= r) {
				number[(tempDNorm / deltaR)] += 2;
			}
			++innerIterator;
		}
		++outerIterator;
	}
	double ri = 0;
	for (int i = 0; i < number.size(); i++) {
		ri = i * deltaR;
		number[i] =
				number[i]
						/ ((double) lcContainer.size()
								* ((4 * M_PI) / 3
										* (pow(ri + deltaR, 3) - pow(ri, 3))));
		rdfSave[para][i] += number[i];
	}
}
