/*
 * Cuboid.h
 *
 *  Created on: 09.11.2013
 *      Author: Son
 */

#ifndef CUBOID_H_
#define CUBOID_H_
#include <vector>
#include <Particle.h>
#include "utils/Vector.h"
#include <list>

/** \class Cuboid
 *  \brief This is a class representing a 3D cuboid of particles.
 * 
 *  One cuboid can save a 3D set of particles, which are located next to each other.
 *  There is no isolated particle in one cuboid.
 */
class Cuboid {
private:

	utils::Vector<double, 3> origin; 		/*!< A 3D vector indicating the location of the lower left corner of a cuboid. */

	std::list<Particle> cub;			/*!< The particles of a cuboid are saved in one single list. */

	int cHeight;					/*!< Cuboid's height (measured in particles). */
	int cWidth;					/*!< Cuboid's width (measured in particles). */
	int cDepth;					/*!< Cuboid's depth (measured in particles). */

	double h;					/*!< The mesh width between every 2 particles. */

	double m;					/*!< The mass of each particle in cuboid. */

	utils::Vector<double, 3> startV;		/*!< The velocity at the beginning of each particle. */

	double meanV;					/*!< The mean velocity (aka. brownian factor) of each particle. */

	int parType;				/*< Type of particles. */

	/**
	 * For Lennard Jones Force.
	 */
	double EPSILON;

	/**
	 * For Lennard Jones Force.
	 */
	double SIGMA;

public:
	/** Default constructor.
	 */
	Cuboid();

	// Constructor: Height first, then width and depth
	// Height --> Width --> Depth
	/** /brief The main constructor.
	 * 
	 * Constructs a new cuboid with all information needed.
	 * @param[in] height Cuboid's height in particles.
	 * @param[in] width Cuboid's width in particles.
	 * @param[in] depth Cuboid's depth in particles.
	 * @param[in] distance Cuboid's mesh width.
	 * @param[in] mass Mass of each particle in cuboid.
	 * @param[in] ori Cuboid's lower left corner's 3D coordinate.
	 * @param[in] startVelocity Velocity at the beginning of each particle in cuboid.
	 * @param[in] meanVelocity Mean velocity (aka. brownian factor) of each partcile in cuboid.
	 * @param[in] parType Particle type of cuboid.
	 * @param[in] EPSILON For Lennard Jones Force.
	 * @param[in] SIGMA For Lennard Jones Force.
	 */
	Cuboid(int height, int width, int depth, double distance, double mass,
			utils::Vector<double, 3> ori, utils::Vector<double, 3> startVelocity,
			double meanVelocity, int parType, double EPSILON, double SIGMA);

	/** \returns Cuboid's lower left corner's location.
	 */
	utils::Vector<double, 3>& getOrigin();

	/** \returns Velocity at the beginning of each particle in cuboid.
	 */
	utils::Vector<double, 3>& getStartV();

	/** \returns Cuboid's list of particles.
	 */
	std::list<Particle>& getCuboid();

	/** \returns Cuboid's height in particles.
	 */
	int getHeight();
	
	/** Setter for height.
	 * @param[in] newH The new value of height to be set.
	 */
	void setHeight(double newH);

	/** \returns Cuboid's width in particles.
	 */
	int getWidth();
	
	/** Setter for width.
	 * @param[in] newW The new value of width to be set.
	 */
	void setWidth(double newW);

	/** \returns Cuboid's depth in particles.
	 */
	int getDepth();
	
	/** Setter for depth.
	 * @param[in] newD The new value of width to be set.
	 */
	void setDepth(double newD);

	/** \returns Cuboid's mesh width between every 2 particles.
	 */
	double getDistance();
	
	/** Setter for mesh width between every 2 particles.
	 * @param[in] newD The new value of mesh width to be set.
	 */
	void setDistance(double newD);

	/** \returns Cuboid's mas of every particle.
	 */
	double getMass();
	
	/** Setter for mass every particle.
	 * @param[in] newM The new value of mass to be set.
	 */
	void setMass(double newM);

	/** \returns Cuboid's mean velocity of every particle.
	 */
	double getMeanV();
	
	/** Setter for mean velocity of every particle.
	 * @param[in] newV The new value of mean velocity to be set.
	 */
	void setMeanV(double newV);

    /**
     * initializes the neighbor lists of particles.
     */
    void initNeighbours();

    Particle * getParticleAtID(int id);

	/** \returns Cuboid's size in particles.
	 */
	int getSize();

	/** \returns Cuboid's particle type.
	 */
	int& getType();

	/** \returns Cuboid's epsilon for Lennard Jones Force.
	 */
	double& getEpsilon();

	/** \returns Cuboid's sigma for Lennard Jones Force.
	 */
	double& getSigma();

	virtual ~Cuboid();
};

#endif /* CUBOID_H_ */
