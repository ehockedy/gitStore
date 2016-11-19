// Translate this file with
//
// g++ -O3 spaceboddies.c -o spaceboddies
//
// Run it with
//
// ./spaceboddies
//
// Open Paraview (www.paraview.org) and do the following:
// - Select File/open and select all the results files. Press the Apply button.
// - Click into the left visualisation window (usually titled Layout #1).
// - Click the result-* item in the window Pipeline Browser. Ensure that your Layout #1 and the item result-* is marked.
// - Select Filters/Alphabetical/TableToPoints. Press Apply button.
// - Switch the representation (on top) from Surface into Points.
// - Press the play button and adopt colours and point sizes.
// - For some Paraview versions, you have to mark your TableToPoints item (usually called TableToPoints1) and explicitly select that X Column is x, Y Column is y, and Z Column is z.
// - What is pretty cool is the Filter TemporalParticlesToPathlines. If you set Mask Points to 1, you see a part of the trajactory.
//
// (C) 2015 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <omp.h>

using namespace std;

string questionNum = "21"; 
int numParticles = 2;

double (*particlePositions)[3];//[numParticles][3]; //Holds the positions of the particles
double (*particleVelocities)[3];//[numParticles][3]; 

double (*forceStore)[3]; //Holds the current calculated forces. Ths prevents calculating the same force twice.

double s = pow(10, -5);//pow(3.405, -10);
double a = pow(10, -5);//pow(1.654, -21);

double timeStepSize = 1e6;//0.0000001;

bool closest = false;
bool first = true; //TESTING

double trackDist = 0.2;

double maxVelocity = 0;

double minDist = 1.0;

double distOfParticles = 1;// 1.122 * s; //1.122 from 1st paragraph of https://en.wikipedia.org/wiki/Lennard-Jones_potential

double randPos() //from http://stackoverflow.com/questions/6218399/how-to-generate-a-random-number-between-0-and-1
{
    return (double)rand() / (double)RAND_MAX; //Generates a random number between -1 and 1 exclusive, then scale to particle distances
}

void setUp(string q) {
	srand(1); //Set the seed so can reproduce results.
	
	if(q == "1" || q == "2")
	{//Position
		for(int p = 0; p < numParticles ; p++) //For each particle
		{
			particlePositions[p][0] = randPos(); //Generate a random x position.
			particlePositions[p][1] = randPos(); //Generate a random y position.
			particlePositions[p][2] = randPos(); //Generate a random z position.
			
			particleVelocities[p][0] = 0.0;//randPos(); //Generate a random x Velocity.
			particleVelocities[p][1] = 0.0;//randPos(); //Generate a random y Velocity.
			particleVelocities[p][2] = 0.0;//randPos(); //Generate a random z Velocity.
			
			forceStore[p][0] = 0.0;
			forceStore[p][1] = 0.0;
			forceStore[p][2] = 0.0;
		}
	}
	else if(q == "21")
	{
		particlePositions[0][0] = 0.4; //Generate a random x position.
		particlePositions[0][1] = 0.5; //Generate a random y position.
		particlePositions[0][2] = 0.5;
		particlePositions[1][0] = 0.6; //Generate a random x position.
		particlePositions[1][1] = 0.5; //Generate a random y position.
		particlePositions[1][2] = 0.5;
		
		particleVelocities[0][0] = 0.0;
		particleVelocities[0][1] = 0.0;
		particleVelocities[0][2] = 0.0;
		particleVelocities[1][0] = 0.0;
		particleVelocities[1][1] = 0.0;
		particleVelocities[1][2] = 0.0;
		
		forceStore[0][0] = 0.0;
		forceStore[0][1] = 0.0;
		forceStore[0][2] = 0.0;
		forceStore[1][0] = 0.0;
		forceStore[1][1] = 0.0;
		forceStore[1][2] = 0.0;
	}
	else if(q == "22")
	{
		particlePositions[0][0] = 0.1; //Generate a random x position.
		particlePositions[0][1] = 0.5; //Generate a random y position.
		particlePositions[0][2] = 0.5;
		particlePositions[1][0] = 0.9; //Generate a random x position.
		particlePositions[1][1] = 0.5; //Generate a random y position.
		particlePositions[1][2] = 0.5;
		
		particleVelocities[0][0] = 0.0;
		particleVelocities[0][1] = 0.0;
		particleVelocities[0][2] = 0.0;
		particleVelocities[1][0] = 0.0;
		particleVelocities[1][1] = 0.0;
		particleVelocities[1][2] = 0.0;
		
		forceStore[0][0] = 0.0;
		forceStore[0][1] = 0.0;
		forceStore[0][2] = 0.0;
		forceStore[1][0] = 0.0;
		forceStore[1][1] = 0.0;
		forceStore[1][2] = 0.0;
	}
}

void printCSVFile(double counter) {
  std::stringstream filename;
  filename << "data/result-" << counter <<  ".csv";
  std::ofstream out( filename.str().c_str() );

  out << "x, y, z" << std::endl;
  for (int i=0; i<numParticles; i++) {
    out << particlePositions[i][0]
        << ","
        << particlePositions[i][1]
        << ","
        << particlePositions[i][2]
        << std::endl;
  }
}

void updateBody(int t, string q)
{
	//cout << particlePositions[0][0] << " efwefefef " << endl;
    //double force[3];
	
	//The problem with a large time step is that because the velocity is based on the previous velocity and the force times the time step, then one large time step will make the new velocity very big, and so in the next interation this velocity is also very big. SO it gets to the point where the velocity keeps increasing as the forces are so tiny it doesn't affect the velocity once it is huge. However, if there is a small time step, the velocity can be influenced by the forces, so doesn't increase so dramatically.	
  
	/*if(q == 2)
	{
		if(t >= 7550)
		{
			timeStepSize = 0.00005;
		}
		else if(t >= 6620)
		{
			timeStepSize = 0.0002;
		}
		else if(t >= 1300)
		{
			timeStepSize = 0.01;
		}
		else if(t >= 900)
		{
			timeStepSize = 1;
		}
		else if(t >= 300)
		{
			timeStepSize = 10;
		}
	}*/
	
	for(int p = 0; p < numParticles ; p++) //For each particle
	{
		forceStore[p][0] = 0.0;
		forceStore[p][1] = 0.0;
		forceStore[p][2] = 0.0;
	}
    double distTest = 0;
	double distTest2 = 0;
    //force and velocity
	#pragma opm parallel for
	{
    for(int i=0; i<numParticles; i++){ //For each particle
        /*force[0] = 0.0; //Force in each direction
        force[1] = 0.0;
        force[2] = 0.0;
		
		#pragma opm parallel for
		{
		
        for(int j=0; j<numParticles; j++){ //For each particle. j = i + 1 as don't need to do cach calculation twice (e.g. for particle 1 and 2 then 2 and 1, just do both in one go.)
            if(i == j) continue; //If it is not the same particle. continue skips rest of current iteration if condition is true
            const double distance = sqrt( //Calculate the absolute distance between the particles
			(particlePositions[i][0]-particlePositions[j][0]) * (particlePositions[i][0]-particlePositions[j][0]) +
            (particlePositions[i][1]-particlePositions[j][1]) * (particlePositions[i][1]-particlePositions[j][1]) +
            (particlePositions[i][2]-particlePositions[j][2]) * (particlePositions[i][2]-particlePositions[j][2])
            );
			//printf("DIST: %lf\n", distance);
            double rx = (particlePositions[i][0]-particlePositions[j][0]); //Proportion of force in the x direction
			double ry = (particlePositions[i][1]-particlePositions[j][1]);  //j - i because it is the force the other particle exerts on this particle
			double rz = (particlePositions[i][2]-particlePositions[j][2]); 
			
            force[0] += 4 * a * (12 * (pow(s, 12)/pow(distance, 13)) - 6 * pow(s, 6)/pow(distance, 7)) * (rx/distance); //Calculate force in the x direction and update the overall force on the particle in theis direction. Overall since being forced by every particle
            force[1] += 4 * a * (12 * (pow(s, 12)/pow(distance, 13)) - 6 * pow(s, 6)/pow(distance, 7)) * (ry/distance); //Force in y dir
            force[2] += 4 * a * (12 * (pow(s, 12)/pow(distance, 13)) - 6 * pow(s, 6)/pow(distance, 7)) * (rz/distance); //Force in z dir
		}
		
		} //END PRAGMA*/
		
		
		
		//#pragma opm parallel for
		//{
		
        for(int j=i+1; j<numParticles; j++){ //For each particle. j = i + 1 as don't need to do cach calculation twice (e.g. for particle 1 and 2 then 2 and 1, just do both in one go.)
            const double distance = sqrt( //Calculate the absolute distance between the particles
			(particlePositions[i][0]-particlePositions[j][0]) * (particlePositions[i][0]-particlePositions[j][0]) +
            (particlePositions[i][1]-particlePositions[j][1]) * (particlePositions[i][1]-particlePositions[j][1]) +
            (particlePositions[i][2]-particlePositions[j][2]) * (particlePositions[i][2]-particlePositions[j][2])
            );
			
			/*if(distance < minDist)
			{
				minDist = distance;
			}
			else
			{*/
			distTest = distance;
			
			//trackDist = distance;
				//cout << "mindist: " << minDist << "  distance: " << distance << " time " << t << endl;
			//}
			
            double rx = (particlePositions[i][0]-particlePositions[j][0]); //Proportion of force in the x direction
			double ry = (particlePositions[i][1]-particlePositions[j][1]);  //j - i because it is the force the other particle exerts on this particle
			double rz = (particlePositions[i][2]-particlePositions[j][2]); 
			
			double f = 4 * a * (12 * (pow(s, 12)/pow(distance, 13)) - 6 * pow(s, 6)/pow(distance, 7));
			
			double fx = f * (rx/distance);
			double fy = f * (ry/distance);
			double fz = f * (rz/distance);
			
            forceStore[i][0] += fx; //Calculate force in the x direction and update the overall force on the particle in this direction. Overall since being forced by every particle
            forceStore[i][1] += fy; //Force in y dir
            forceStore[i][2] += fz; //Force in z dir
			
			forceStore[j][0] -= fx; //Calculate force in the x direction and update the overall force on the particle in this direction. Overall since being forced by every particle
            forceStore[j][1] -= fy; //Force in y dir
            forceStore[j][2] -= fz;//Force in z dir
			
			//cout << distance << "    " << (rx/distance) << "    " << (ry/distance) << "     " << (rz/distance) << endl;
			
			//Wrap around distances
			//cout << distance << endl;
			//if(t<3){cout << "Dist: " << distance << endl;}
			double altDistance = fmod(1.0 - distance, 1.0); //do fmod incase distance is 0 (even though it never is)
			
			distTest2 = altDistance;
			
			//double altDistance = 1-distance;
			//cout << altDistance;
			/*if(distance > 0)
			{
				altDistance = 0
			}*/
			/*if(t < 3)
			{
				cout << "1: " << (1 - fabs(ry)) << "   " << (ry/fabs(ry)) << endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			if(rx != 0)
			{
				rx = -1 * (1 - fabs(rx)) * (rx/fabs(rx));//1-rx;//-1*rx; //
			}
			else
			{
				rx = 0; //because of divide by zero from rx/fabs(rx)
			}
			
			if(ry != 0)
			{
				ry = -1 * (1 - fabs(ry)) * (ry/fabs(ry));//1-rx;//-1*rx; //
			}
			else
			{
				ry = 0;
			}
			
			if(rz != 0)
			{
				rz = -1 * (1 - fabs(rz)) * (rx/fabs(rz));//1-rx;//-1*rx; //
			}
			else
			{
				rz = 0;
			}
			
			//ry = -1 * (1 - fabs(ry)) * (ry/fabs(ry));//1-ry;//1 - ry;//-1*ry; //
			//rz = -1 * (1 - fabs(rz)) * (ry/fabs(rz));//1-rz;//(1 - fabs(rz)) * (rz/fabs(rz));//1 - rz;//-1*rz; //
			/*if(t < 3)
			{
				cout << "2: " << ry << endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			
			
			/*if(t < 3)
			{
				cout << "AD: " << altDistance << endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			
			/*if(distance != 0) //WILL NEED TO CHANGE FOR ALL 3 DIRS
			{
				rx = -1*fabs(1-rx);
				ry = -1*fabs(1-ry);
				rz = -1*fabs(1-rz);
			}
			else
			{
				rx = 0;
				ry = 0;
				rz = 0;
			}*/
			
			double f2 = 4 * a * (12 * (pow(s, 12)/pow(altDistance, 13)) - 6 * pow(s, 6)/pow(altDistance, 7));
			/*if(t < 3)
			{
				cout << "F2:  " << f2<< endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			double f2x = f2 * (rx/altDistance);
			/*if(t < 3)
			{
				cout << "ry " << ry << " altDist: " << altDistance << endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			double f2y = f2 * (ry/altDistance);
			double f2z = f2 * (rz/altDistance);
			
			/*if(t < 3)
			{
				cout << "FORCE f2y:  " << f2y<< endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			
			forceStore[i][0] += f2x; //Calculate force in the x direction and update the overall force on the particle in this direction. Overall since being forced by every particle
            forceStore[i][1] += f2y; //Force in y dir
            forceStore[i][2] += f2z; //Force in z dir
			
			forceStore[j][0] -= f2x; //Calculate force in the x direction and update the overall force on the particle in this direction. Overall since being forced by every particle
            forceStore[j][1] -= f2y; //Force in y dir
            forceStore[j][2] -= f2z; //Force in z dir
			
			/*if(t < 3)
			{
				cout << "Force Store: [j][1] " << forceStore[j][1]<< endl;//"re: " << rx << "fabs" << 1- fabs(rx) << endl;
			}*/
			
			//cout << altDistance << "    " << (rx/altDistance) << "    " << (ry/altDistance) << "     " << (rz/altDistance) << endl;
		
			/*if(forceStore[i][0] + forceStore[j][0] == 0)
			{
				cout << "FORCES SAME  i=" << i << "  j=" << j << "  t=" << t << endl; 
			}*/
		}
		
		
		//} //END PRAGMA
        
        /*particleVelocities[i][0] = particleVelocities[i][0] + timeStepSize * force[0]; //Update velocities
        particleVelocities[i][1] = particleVelocities[i][1] + timeStepSize * force[1]; 
        particleVelocities[i][2] = particleVelocities[i][2] + timeStepSize * force[2];*/
		/*if(t == 100000 || t == 200000)
		{
			cout << i << " velocity: " << particleVelocities[i][0] << endl;
			cout << i << " force: " << forceStore[i][0] << endl;
		}*/
		particleVelocities[i][0] = particleVelocities[i][0] + timeStepSize * forceStore[i][0]; //Update velocities
        particleVelocities[i][1] = particleVelocities[i][1] + timeStepSize * forceStore[i][1]; 
        particleVelocities[i][2] = particleVelocities[i][2] + timeStepSize * forceStore[i][2];
		
		/*if(t < 2)
		{
			cout << "INLOOPv1: " << particleVelocities[0][0] << endl;
			cout << "INLOOPv2: " << particleVelocities[1][0] << endl;
			cout << "INLOOPv after:" << (particleVelocities[0][0] + particleVelocities[1][0]) << endl;
			//cout << (particlePositions[0][1] - particlePositions[1][1]) << endl;
			//cout << (particlePositions[0][2] - particlePositions[1][2]) << endl;
		}*/
    }
	}//end pragma
	
	if(distTest < trackDist*0.75)
	{
		timeStepSize = timeStepSize/2;
		trackDist = distTest;
		cout << "NORM t: " << t << " d: " << distTest << " v: " << particleVelocities[0][0] << " s: " << timeStepSize << endl;
	}
	else if(distTest2 < trackDist*0.75)
	{
		timeStepSize = timeStepSize/2;
		trackDist = distTest2;
		cout << "WRAP t: " << t << " d: " << distTest2 << " v: " << particleVelocities[0][0] << " s: " << timeStepSize << endl;
	}
	/*if(distTest2 < trackDist2*0.75)
	{
		timeStepSize = timeStepSize/2;
		trackDist2 = distTest2;
		cout << "time: " << t << " distance: " << distTest << " time: " << t << " velocity: " << particleVelocities[0][0] << "step: " << timeStepSize << endl;
	}*/
	
	//bool doneThisRound = false;
	/*for (int i = 0; i < numParticles; i++)
	{
		//if(doneThisRound == false)
		//{
			double overallVelocity = sqrt(particleVelocities[i][0] * particleVelocities[i][0]
											+ particleVelocities[i][1] * particleVelocities[i][1]
											+ particleVelocities[i][2] * particleVelocities[i][2]);
												
			//if(distTest >= minDist and fmod(t, 100000) == 0)
			//{
				//cout << "mindist: " << minDist << "  distance: " << distTest << " time " << t << " velocity: " << overallVelocity << endl;
				//timeStepSize = timeStepSize*2;
			//}
											
			int mult = 3;
			int mult2 = 2;
			
			
			
			if(distTest < 0.0005)//overallVelocity <)
			{
				//cout << "DISTANCE " << distTest << endl;
				//timeStepSize = 0.00000005;
				mult2 = 10;
			}
			if(distTest < 0.00005)//overallVelocity <)
			{
				//cout << "DISTANCE " << distTest << endl;
				timeStepSize = 0.00000005;
				mult2 = 2;
			}
			
			if(maxVelocity == 0)
			{
				maxVelocity = overallVelocity;
			}
			else if(overallVelocity >= mult*maxVelocity) //This should hopefully mean that whenever the velocity gets too large, the timestep changes accordingly
			{
				cout << "CHANGE: " << t << "  step size: " << timeStepSize << "  velocity: " << particleVelocities[i][0] << "  distance: " << distTest << endl;
				maxVelocity = overallVelocity;
				timeStepSize = timeStepSize/mult2;
				i = numParticles;
			}
			
			
			
			//if(timeStepSize < 1000)
			//{
			//	mult = 5;
			//}
		//}
	}*/
	
	/*if(t < 2)
	{
		cout << "v1: " << particleVelocities[0][0] << endl;
		cout << "v2: " << particleVelocities[1][0] << endl;
		cout << "v after:" << (particleVelocities[0][0] + particleVelocities[1][0]) << endl;
		//cout << (particlePositions[0][1] - particlePositions[1][1]) << endl;
		//cout << (particlePositions[0][2] - particlePositions[1][2]) << endl;
	}*/
		
	/*if(minDist <= 0.00001)
	{
		timeStepSize = 0.00001;
	}
	else */
	
	/*if(closest == false)
	{
		if(minDist <= 0.000005)
		{
			timeStepSize = 0.00001;
			closest = true;
		}
		else if(minDist <= 0.00005)
		{
			timeStepSize = 0.0001;
		}
		else if(minDist <= 0.0001)
		{
			timeStepSize = 0.001;//0.005
		}
		else if(minDist <= 0.0005)
		{
			timeStepSize = 0.01;//0.5;
		}
		else if(minDist <= 0.001)
		{
			timeStepSize = 0.1;//5
		}
		else if(minDist <= 0.005)
		{
			timeStepSize = 10;//5000;
		}
		else if(minDist <= 0.05)
		{
			timeStepSize = 1000;//1e6;
		}
		else if(minDist <= 0.5)
		{
			timeStepSize = 1000000;//1e8;
		}
	}*/
		
    //printf("%lf\n", force[0]);
    //position
	
	/*if(((particlePositions[0][0] - (1- particlePositions[1][0])) > 4.69e-013) && first == true)
	{
		cout << "POSITIONS DIFFERENT BEFORE t=" << t << "  difference=" << (particlePositions[0][0] - (1- particlePositions[1][0])) <<endl; 
		first = false;
	}
	
	if(((particleVelocities[0][0] + particleVelocities[1][0]) > 1e-15) && first == true)
	{
		cout << "VELOCITY DIFFERENT BEFORE t=" << t << "  difference=" << (particleVelocities[0][0] - (-1*particleVelocities[1][0])) << endl; 
		first = false;
	}
	
	if(t < 2)
	{
		//cout << "x pos before: " << (particlePositions[0][0] - particlePositions[1][0]) << endl;
		//cout << (particlePositions[0][1] - particlePositions[1][1]) << endl;
		//cout << (particlePositions[0][2] - particlePositions[1][2]) << endl << endl;
	}*/	
	
	//Put pragma out here, 
	#pragma opm parallel for
	{
    for(int i=0; i<numParticles; i++){ //Update positions of all particles 
		//#pragma opm parallel for
		//{
		
		for(int n = 0; n <3; n++)
		{
			double updatedPos = particlePositions[i][n] + timeStepSize * particleVelocities[i][n];
			/*if(t > 100000 && t < 100010)
			{
				cout << "INLOOP " << updatedPos << endl;
			}*/
			/*if(t == 10237 && i == 2)
			{
				cout << "pP: " << particlePositions[i][n] << endl << "pV: " << particleVelocities[i][n] << endl << "uP: " << updatedPos << endl <<"fa(uP):"  << fabs(updatedPos) << endl <<"fl(fa(uP)):"  << floor(fabs(updatedPos)) <<endl << "1-fa-fl(fa): " << 1 - (fabs(updatedPos) - floor(fabs(updatedPos))) << endl << endl; 
			}*/
			/*if((t == 10200 || t == 10300 || t == 10400) && i == 2)
			{
				cout << "i = " << i << ", n =" << n << endl << "t = " << t << ": particlePosition = " << particlePositions[i][n] << endl; 
			}*/
			if(updatedPos >= 1) //Could be something to do with whether 1 is inclusive or not as to why 0.4/0.6 and 0.1/0.9 differ
			{
				updatedPos = updatedPos - floor(updatedPos);
				//cout << "1" << endl;
			}
			else if(updatedPos < 0)
			{
				updatedPos = 1 - (fabs(updatedPos) - floor(fabs(updatedPos)));
				//cout << "2" << endl;
			}
			/*if((t == 10200 || t == 10300 || t == 10400) && i == 2)
			{
				cout << "t = " << t << ": updatedPos = " << updatedPos << endl;
			}*/
			particlePositions[i][n] = updatedPos;
			//if(t<3){cout << particlePositions[i][n] << " " << i << "  "<< n << endl;}
		}
		
		//cout << particlePositions[i][0] << i << endl;
		//}
		/*double updatedPosx = particlePositions[i][0] + timeStepSize * particleVelocities[i][0];
		double updatedPosy = particlePositions[i][1] + timeStepSize * particleVelocities[i][1];
		double updatedPosz = particlePositions[i][2] + timeStepSize * particleVelocities[i][2];
	
        if()
		particlePositions[i][0] =  //Since fmd wasn't working with negatives, this does it just on the absolute value of the new value, then makes it pos or neg
        particlePositions[i][1] = fmod(fabs(updatedPosy), distOfParticles) * updatedPosy/fabs(updatedPosy); //fabs is for doubles, fabsf for floats, abs is for ints
        particlePositions[i][2] = fmod(fabs(updatedPosz), distOfParticles) * updatedPosz/fabs(updatedPosz); 
		*/
	}
	}//end pragma
	
	
	/*if(t > 100 && t < 120)
	{
		cout << "x1: " << particlePositions[0][0] << endl;
		cout << "x2: " << particlePositions[1][0] << endl;
		cout << "xdist:" << (particlePositions[1][0] - particlePositions[0][0]) << endl << endl;
		//cout << (particlePositions[0][1] - particlePositions[1][1]) << endl;
		//cout << (particlePositions[0][2] - particlePositions[1][2]) << endl;
	}*/
	
	if(((1-(particlePositions[0][1] + particlePositions[1][1])) > 4.69e-013) && first == true)
	{
		cout << "POSITIONS DIFFERENT AFTER t=" << t << "  difference=" << (particlePositions[0][0] - (1- particlePositions[1][0])) <<endl; 
		first = false;
	}
	
	if(((particleVelocities[0][0] - (-1*particleVelocities[1][0])) > 1e-10) && first == true)
	{
		cout << "VELOCITY DIFFERENT AFTER t=" << t << endl; 
		first = false;
	}
	
	

}
/*
stepsize:1e-007, v:5.44e-8
ss: 1.95e-10, v=0.00107


*/


int main(int argc, char *argv[]) {
	//cout << fmod(0.9123,1) << "  " << fmod(1.0,1.0);
	//cout << argv[1] << "   " << atoi(argv[2]);
	questionNum = argv[1]; //THE QUESTION NUMBER. OPTIONS ARE: 1, 2, 21, 22
	numParticles = atoi(argv[2]);
	if(questionNum == "21" || questionNum == "22"){numParticles = 2;} //Make sure only have 2 particles for q2 parts 1 and 2 
	particlePositions = new double[numParticles][3];
	particleVelocities = new double[numParticles][3];
	forceStore = new double[numParticles][3];
  setUp(questionNum);
  printCSVFile(0);
  cout.precision(10);
  const double timeSteps = 1600000000;//5e10;//30000;
  double plotEveryKthStep = 800000;//1e6;//100;
  for (double i=0; i<timeSteps; i++) { //2145000000
    updateBody(i, questionNum);
	if(fmod(i,plotEveryKthStep) == 0)
	{
		cout << "i: " << i << "   pos: " << particlePositions[0][0] << endl;
		//cout << i/plotEveryKthStep+1 << " of " << timeSteps/plotEveryKthStep <<  endl;
		printCSVFile(i/plotEveryKthStep+1); // Please switch off all IO if you do performance tests.
	}

  }
  delete [] particlePositions;
  delete [] particleVelocities;

  return 0;
}
