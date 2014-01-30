/*
 * MembraneTest.cpp
 *
 *  Created on: Jan 12, 2014
 *      Author: son
 */

#include <string>
#include <cstring>
#include <algorithm>
#include <iostream>

#include "MembraneTest.h"
#include "utils/Vector.h"

MembraneTest::MembraneTest() {
	// TODO Auto-generated constructor stub

}

MembraneTest::~MembraneTest() {
	// TODO Auto-generated destructor stub
}

void MembraneTest::setUp(){
	std::string fileName = "src/tests/testFiles/TestMembrane.xml";
	char *cstr = new char[fileName.length() + 1];
	strcpy(cstr, fileName.c_str());

	pgen.extractCuboids(cstr);
	cub = *pgen.getCuboidList().begin();
	cub.initNeighbours();
	parList = cub.getCuboid();
}

void MembraneTest::tearDown(){

}

void MembraneTest::testInitNeighbours(){
	testGetDirectNeighbours();
	testGetDiagNeighbours();
}

void MembraneTest::testGetParticleAtID(){
	for (int i = 0; i<100; i++){
		Particle * p = cub.getParticleAtID(i);
		//CPPUNIT_ASSERT(p->getID() != -1);
		CPPUNIT_ASSERT((*p).getID() == i);
	}

	Particle * p = cub.getParticleAtID(-1);
	CPPUNIT_ASSERT(p == NULL);
	p = cub.getParticleAtID(-2);
	CPPUNIT_ASSERT(p == NULL);
	p = cub.getParticleAtID(-20);
	CPPUNIT_ASSERT(p == NULL);

	p = cub.getParticleAtID(101);
	CPPUNIT_ASSERT(p == NULL);
	p = cub.getParticleAtID(102);
	CPPUNIT_ASSERT(p == NULL);
	p = cub.getParticleAtID(200);
	CPPUNIT_ASSERT(p == NULL);

}

void MembraneTest::testGetDirectNeighbours(){
	for (std::list<Particle>::iterator it = parList.begin();
			it != parList.end(); it++){
		CPPUNIT_ASSERT((*it).getDirectNeighbours().size() >= 2
						&& (*it).getDirectNeighbours().size() <= 4);
		for(std::list<int>::iterator its = (*it).getDirectNeighbours().begin();
				its != (*it).getDirectNeighbours().end(); its++){
			//no more pNull in the list
			//CPPUNIT_ASSERT(*its != NULL);

			double distance = ((*it).getX() - (*cub.getParticleAtID(*its)).getX()).L2Norm();
			CPPUNIT_ASSERT(distance <= 1.0);
		}
	}
}

void MembraneTest::testGetDiagNeighbours(){
	for (std::list<Particle>::iterator it = parList.begin();
			it != parList.end(); it++){
		CPPUNIT_ASSERT((*it).getDiagNeighbours().size() >= 1
						&& (*it).getDiagNeighbours().size() <= 4);
		for(std::list<int>::iterator its = (*it).getDiagNeighbours().begin();
				its != (*it).getDiagNeighbours().end(); its++){
			//no more pNull in the list
			//CPPUNIT_ASSERT(*its != NULL);

			double distance = ((*it).getX() - (*cub.getParticleAtID(*its)).getX()).L2Norm();
			CPPUNIT_ASSERT(distance <= 1.0*sqrt(2.0));
		}
	}
}

void MembraneTest::testGetID(){
	std::list<int> idList;
	idList.clear();
	for (std::list<Particle>::iterator it = parList.begin();
			it != parList.end(); it++){
		idList.push_back((*it).getID());
	}
	//test size
	CPPUNIT_ASSERT(idList.size() == parList.size());

	//test uniqueness of each ID (each ID must exist only once)
	idList.sort();

	int i = 0;
	for (std::list<int>::iterator iti = idList.begin();
			iti != idList.end(); iti++){
		CPPUNIT_ASSERT(*iti == i);
		i++;
	}

	//test coordinates
	//input: h=1, size=10x10, origin={0,0,0}
	for (std::list<Particle>::iterator it = parList.begin();
			it != parList.end(); it++){
		int id = (*it).getID();
		int pos_of_line = id/10;
		int pos_of_col = id % 10;
		CPPUNIT_ASSERT((*it).getX()[0] == pos_of_col*1.0);
		CPPUNIT_ASSERT((*it).getX()[1] == pos_of_line*1.0);
	}
}

void MembraneTest::testIsDirectNeighbour(){
	/*
	for (std::list<Particle>::iterator it1 = parList.begin();
			it1 != parList.end(); it1++){
		std::list<Particle> listS = (*it1).getDirectNeighbours();
		CPPUNIT_ASSERT(listS.size() >= 2
						&& listS.size() <= 4);
		for (std::list<Particle>::iterator it2 = parList.begin();
				it2 != parList.end(); it2++){
			for(std::list<Particle>::iterator its = listS.begin();
					its != listS.end(); its++){
				if ((*its).getID() == (*it2).getID()){
					CPPUNIT_ASSERT((*it1).isDirectNeighbour(*it2) == true);
				}
				else
					CPPUNIT_ASSERT((*it1).isDirectNeighbour(*it2) == false);
			}
		}

	}
	*/

	for (std::list<Particle>::iterator it1 = parList.begin();
			it1 != parList.end(); it1++){
		CPPUNIT_ASSERT((*it1).getDirectNeighbours().size() >= 2
						&& (*it1).getDirectNeighbours().size() <= 4);
		for (std::list<Particle>::iterator it2 = parList.begin();
				it2 != parList.end(); it2++){
			double distance = ((*it1).getX() - (*it2).getX()).L2Norm();
			if ((*it1).isDirectNeighbour(*it2)){
				CPPUNIT_ASSERT(distance <= 1.0);
			}
		}
	}
}

void MembraneTest::testIsDiagNeighbour(){
	/*
	for (std::list<Particle>::iterator it1 = parList.begin();
			it1 != parList.end(); it1++){
		std::list<Particle> listS = (*it1).getDiagNeighbours();
		CPPUNIT_ASSERT(listS.size() >= 1
						&& listS.size() <= 4);
		for (std::list<Particle>::iterator it2 = parList.begin();
				it2 != parList.end(); it2++){
			for(std::list<Particle>::iterator its = listS.begin();
					its != listS.end(); its++){
				if ((*its).getID() == (*it2).getID()){
					CPPUNIT_ASSERT((*it1).isDiagNeighbour(*it2) == true);
				}
				else
					CPPUNIT_ASSERT((*it1).isDiagNeighbour(*it2) == false);
			}
		}
	}
	*/

	for (std::list<Particle>::iterator it1 = parList.begin();
			it1 != parList.end(); it1++){
		CPPUNIT_ASSERT((*it1).getDiagNeighbours().size() >= 1
						&& (*it1).getDiagNeighbours().size() <= 4);
		for (std::list<Particle>::iterator it2 = parList.begin();
				it2 != parList.end(); it2++){
			double distance = ((*it1).getX() - (*it2).getX()).L2Norm();
			if ((*it1).isDiagNeighbour(*it2)){
				CPPUNIT_ASSERT(distance <= 1.0*sqrt(2.0));
			}
		}
	}
}

CppUnit::Test *MembraneTest::suite() {
	CppUnit::TestSuite *testSuite = new CppUnit::TestSuite(
			"ParticleContainerTest");

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testGetDiagNeighbours",
					&MembraneTest::testGetDiagNeighbours));

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testGetDirectNeighbours",
					&MembraneTest::testGetDirectNeighbours));

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testGetID",
					&MembraneTest::testGetID));

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testGetParticleAtID",
					&MembraneTest::testGetParticleAtID));

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testInitNeighbours",
					&MembraneTest::testInitNeighbours));

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testIsDiagNeighbour",
					&MembraneTest::testIsDiagNeighbour));

	testSuite->addTest(
			new CppUnit::TestCaller<MembraneTest>("testIsDirectNeighbour",
					&MembraneTest::testIsDirectNeighbour));

	return testSuite;
}
