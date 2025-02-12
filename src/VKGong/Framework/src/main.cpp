/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * Top level code, command line parsing, main loop
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <iostream>
#include <fstream>
using namespace std;

#include <unistd.h>

#ifdef WIN32
#include <windows.h>
#include <shlobj.h>
#else
#include <sys/stat.h>
#include <sys/types.h>
#endif

#include "PhysicalSystem.h"
#include "PhysicalSystemParser.h"
#include "ScoreParser.h"
#include "Input.h"

#include "GlobalSettings.h"
#include "SettingsManager.h"
#include "Logger.h"
#include "Profiler.h"

static PhysicalSystem *physicalSystem;

static string physicalSystemName = "physicalSystem.txt";
static string scoreName = "score.txt";

static string outputName = "output";
static bool saveRawOutputs = false;
static bool saveIndividualOutputs = true;

// name of executable, not including any path components
static string exeName;

static void getExecutableLocation(char *buf, int bufsize)
{
#ifndef WIN32
    size_t sz = readlink("/proc/self/exe", buf, bufsize);
    buf[sz] = 0;
    char *slash = strrchr(buf, '/');
    if (slash) *slash = 0;
#else
    GetModuleFileName(NULL, buf, bufsize);
    char *slash = strrchr(buf, '\\');
    if (slash) *slash = 0;
#endif
}


/*
 * Load physicalSystem and score
 *
 * Returns 0 on failure
 */
static int initialise()
{
    
    // get simulation length from score file first!
    double duration = ScoreParser::getDuration(scoreName);

    if (duration < 0.0) {
	logMessage(5, "Unable to get duration from score file!");
	return 0;
    }    
    GlobalSettings::getInstance()->setDuration(duration);
     
    physicalSystem = PhysicalSystemParser::parsePhysicalSystem(physicalSystemName);
   
    if (!physicalSystem) {
	logMessage(5, "Unable to parse physicalSystem file");
	return 0;
    }

    if (!ScoreParser::parseScore(scoreName, physicalSystem)) {
	logMessage(5, "Unable to parse score file");
	delete physicalSystem;
	return 0;
    }


    return 1;
}

/*
 * Run the main simulation loop
 */
static void mainLoop()
{
    int NF = GlobalSettings::getInstance()->getNumTimesteps();
    int n;
    int start;

    if (GlobalSettings::getInstance()->getEstimate()) {
	if (NF > 1000) {
	    NF = 1000;
	}
	start = 0;
    }
    else {
	// start from when first input happens. everything before
	// that is all zeros
	start = Input::getFirstInputTimestep();
	logMessage(1, "Starting from timestep %d", start);
    }
    
    if (GlobalSettings::getInstance()->getStormerVerletScheme()){
		for (n = start; n < NF; n++) {
		if ((n % 1000) == 0) {
			logMessage(1, "Iteration %d", n);
	
			// write progress to a file every 1000 iterations, allowing the frontend
			// to monitor long running simulations
			ofstream of("progress.txt");
			of << n << endl;
	
			of.close();
		}
		physicalSystem->runTimestepVerlet(n);
		physicalSystem->endTimestepVerlet(n);
		}
    }
    else {
		for (n = start; n < NF; n++) {
		if ((n % 1000) == 0) {
			logMessage(1, "Iteration %d", n);
	
			// write progress to a file every 1000 iterations, allowing the frontend
			// to monitor long running simulations
			ofstream of("progress.txt");
			of << n << endl;
	
			of.close();
		}
		physicalSystem->runTimestepECS(n);
		physicalSystem->endTimestepECS(n);
		}
    }
   
}

/*
 * Save outputs and free resources
 */
static void finalise()
{
    if (!GlobalSettings::getInstance()->getEstimate()) {  	
	physicalSystem->saveOutputs(outputName, saveIndividualOutputs, saveRawOutputs);
    }
    delete physicalSystem;
}


/*
 * Parse the command line and configure the code (mostly via GlobalSettings)
 *
 * Returns 0 on failure (invalid command line etc.)
 */
static int parseCommandLine(int argc, char *argv[])
{
    int i;
    GlobalSettings *gs = GlobalSettings::getInstance();
    SettingsManager *sm = SettingsManager::getInstance();

    for (i = 1; i < argc; i++) {
	if (strchr(argv[i], ':')) {
	    // component-specific setting
	    if (i < (argc-1)) {
		if (sm->putSetting(argv[i], argv[i+1])) i++;
	    }
	    else {
		sm->putSetting(argv[i]);
	    }
	    continue;
	}

	if (!strcmp(argv[i], "-i")) {
	    // physicalSystem filename
	    if (i == (argc-1)) return 0;
	    i++;
	    physicalSystemName = argv[i];
	}
	else if (!strcmp(argv[i], "-s")) {
	    // score filename
	    if (i == (argc-1)) return 0;
	    i++;
	    scoreName = argv[i];
	}
	
	else if (!strcmp(argv[i], "-r")) {
		    // save raw outputs
		    saveRawOutputs = true;
	}

	else if (!strcmp(argv[i], "-e")) {
	    // estimate
	    gs->setEstimate(true);
	}
	else if (!strcmp(argv[i], "-o")) {
	    // output base name
	    if (i == (argc-1)) return 0;
	    i++;
	    outputName = argv[i];	    
	}


	else if (!strcmp(argv[i], "-c")) {
	    // save only stereo mix or individual outputs as well
	    if (i == (argc-1)) return 0;
	    i++;
	    if (!strcmp(argv[i], "stereo")) {
		saveIndividualOutputs = false;
	    }
	    else if (!strcmp(argv[i], "all")) {
		saveIndividualOutputs = true;
	    }
	    else return 0;
	}
	
	else if (!strcmp(argv[i], "-t")) {
	    // Choose time integration scheme
	    if (i == (argc-1)) return 0;
	    i++;
	    if (!strcmp(argv[i], "verlet")) {
	    	gs->setStormerVerletScheme(true);
	    	logMessage(1, "Time stepping scheme: Stormer-Verlet");
	    }
	    else if (!strcmp(argv[i], "ecs")) {
	    	gs->setStormerVerletScheme(false);
	    	logMessage(1, "Time stepping scheme: Energy conserving scheme");
	    }
	    else return 0;
	}

	else {
	    return 0;
	}
    }
    return 1;
}

static void usage(char *progname)
{
    fprintf(stderr, "Usage: %s [-i <physicalSystem file>] [-s <score file>] [-o <output filename base>] [-c stereo|all] [-t verlet|ecs] [-r] [-e]\n", progname);

}

static char *profileNames[] = {
    "startup", "main loop", "finalise"
};

int main(int argc, char *argv[])
{    
    if (!parseCommandLine(argc, argv)) {
	usage(argv[0]);
	return 1;
    }

    // get executable name
    exeName = argv[0];
    int lastslash = exeName.rfind("/");
    if (lastslash != string::npos) {
	exeName = exeName.substr(lastslash + 1);
    }
    logMessage(1, "Executable name: %s", exeName.c_str());

    Profiler profiler(3, profileNames);
    
    profiler.start(0);

    if (!initialise()) {
	return 1;    
    }

    
    profiler.end(0);

    profiler.start(1);
    mainLoop();
    profiler.end(1);

    

    if (GlobalSettings::getInstance()->getEstimate()) {
	double simthousandtime = profiler.get(1);
	double simestimatetime = profiler.get(0) + simthousandtime * GlobalSettings::getInstance()->getNumTimesteps() / 1000.0;
	logMessage(5, "Simulation time for 1000 timesteps : %fs", simthousandtime);
	logMessage(5, "Estimated time for full run : %fs", simestimatetime);	
    }

    profiler.start(2);
    finalise();
    profiler.end(2);

    printf("Profile: %s\n", profiler.print().c_str());
    


    return 0;
}
