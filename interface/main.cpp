#include "Surface.h"
#include "SurfPoint.h"
#include "SurfData.h"
#include "PolynomialSurface.h"
#include "KrigingSurface.h"
//#include "ANNSurface.h"
//#include "ParticleSwarmMLP.h"
//#include "BackpropMLP.h"
#include "SurfDataIterator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <string>

typedef struct {
    double min;
    double max;
    int pts;
    double interval;
} Axis;


using namespace std;

vector<Axis> axes;
vector<int> point;
vector<double> surfptx;
int ndims;
int npts;


void readInput(istream& is = cin)
{
    npts = 1;
    const int MAX_CHAR = 200;
    char line[MAX_CHAR];
    is.getline(line,MAX_CHAR);
    string sline = line;
    istringstream streamline(sline);
    streamline >> ndims;
    axes.resize(ndims);
    for (int i = 0; i < ndims; i++) {
        is.getline(line,MAX_CHAR);
        sline = line;
	streamline.str(sline);
        char c;
        streamline >> c;
	if (c == 'v') {
	    streamline >> axes[i].min >> axes[i].max >> axes[i].pts;
            axes[i].interval = (axes[i].max - axes[i].min) / (axes[i].pts - 1);
	    npts *= axes[i].pts;
        } else if (c == 'f') {
            streamline >> axes[i].min;
	    axes[i].max = axes[i].min;
	    axes[i].pts = 1;
	    axes[i].interval = 0.0;
	} else {
    	    cerr << "Expected 'f' or 'v'" << endl;
	}
    }
}    
     
void initialize()
{
    point.resize(ndims);
    for (int i = 0; i < ndims; i++) {
	point[i] = 0;
    }
}

void nextPoint()
{
    int curDim = ndims-1;
    while (axes[curDim].pts == 1 || point[curDim] == (axes[curDim].pts - 1)) {
	curDim--;
    }
    if (curDim < ndims) {
        point[curDim]++;
	curDim++;
	while(curDim < ndims) {
	    point[curDim] = 0;
	    curDim++;
	}
    }
}
void populateSurfData(SurfData* surfData, string filename, Surface* surface=0)
{
    ifstream infile(filename.c_str(),ios::in);
    if (!infile) {
	    cerr << "File " << filename << " not found." << endl;
	    return;
    }
    readInput(infile);
    infile.close();
    initialize();
    surfptx.resize(ndims);
    for (int i = 0; i < npts; i++) {
        for (int j = 0; j < ndims; j++) {
	    surfptx[j] = axes[j].min + axes[j].interval*point[j];
    	    //cout << setw(10) << surfptx[j];
	}
	//cout << endl;
	//double fx = surface ? surface->evaluate(surfptx) : rastrigin(surfptx)*.5;
	surfData->shallowAddPoint(new SurfPoint(surfptx));
	nextPoint();
    }
    //cout << "Number of points: " << npts << endl;
}

void gridPoints(vector< string >& args) 
{
    SurfData sd;
    populateSurfData(&sd, args[1]);

    ofstream outfile(args[2].c_str(), ios::out);
    if (!outfile) {
	    cerr << "Error: unable to open " << args[2] << " for output." << endl;
	    return;
    }

    sd.write(outfile);
    outfile.close();
}

void createSurface(vector< string >& args)
{
    Surface* s = 0;
    SurfData sd;
    ifstream infile(args[1].c_str(), ios::in);
    if (!infile) {
	    cerr << "Error: unable to open" << args[1] << "." << endl;
	    return;
    }
    sd.read(infile);
    
    //ofstream outfile(args[2].c_str(), ios::out);
    //if (!outfile) {
    //        cerr << "Error: unable to open " << args[2] << " for output." << endl;
    //        return;
    //}

    if (args[3] == "kriging") {
	    s = new KrigingSurface(&sd);
	    s->build();
	    s->test(&sd);
    } else if (args[3] == "polynomial") {
	    int order = atoi(args[4].c_str());
	    s = new PolynomialSurface(&sd,order);
	    s->build();
    }
    if (s) {
	    s->saveBinary(args[2]);
    }
    infile.close();
    delete s;
}

void evaluateSurface(vector< string >& args)
{
    Surface* s = 0;
    SurfData sd;
    ifstream infile(args[1].c_str(), ios::in);
    if (!infile) {
	    cerr << "Error: unable to open" << args[1] << " for input." << endl;
	    return;
    }
    sd.read(infile);

    ifstream insurface(args[3].c_str(), ios::in | ios::binary);
    if (!insurface) {
	    cerr << "Error: unable to open" << args[3] << " for input." << endl;
	    return;
    }
    
    ofstream outfile(args[2].c_str(), ios::out);
    if (!outfile) {
	    cerr << "Error: unable to open " << args[2] << " for output." << endl;
	    return;
    }
    int surfaceID;
    insurface.read((char*)&surfaceID, sizeof(int));
    insurface.close();
    switch (surfaceID) {
	    case Surface::polynomialSurfaceID: 
	        s=new PolynomialSurface(args[3]);
		cout << "Polynomial" << endl;
	        break;
	    case Surface::krigingSurfaceID: 
		s = new KrigingSurface(args[3]);
	        cout << "Kriging" << endl; 
		break;
	    default:
		cout << "Unknown Surface" << endl;
    }
    s->evaluate(&sd);
    
    sd.write(outfile);
    infile.close();
    outfile.close();
    delete s;
}

void computeErrorMetric(vector< string>& args)
{
    ifstream insurface(args[1].c_str(), ios::in | ios::binary);
    if (!insurface) {
	    cerr << "Error: unable to open" << args[1] << " for input." << endl;
	    return;
    }

    int surfaceID;
    insurface.read((char*)&surfaceID, sizeof(int));
    insurface.close();
    Surface* s = 0;
    switch (surfaceID) {
	    case Surface::polynomialSurfaceID: 
	        s=new PolynomialSurface(args[1]);
		cout << "Polynomial" << endl;
	        break;
	    case Surface::krigingSurfaceID: 
		s = new KrigingSurface(args[1]);
	        cout << "Kriging" << endl; 
		break;
	    default:
		cout << "Unknown Surface" << endl;
    }
    double errorValue = s->errorMetric(args[2]);
    cout << args[2] << ": " << errorValue << endl;
    delete s;


}
void executeCommand(vector< string >& args) 
{
	cout << "Executing command: ";
	for (unsigned i = 0; i < args.size(); i++) {
		cout << args[i] << " ";
	}
	cout << endl;

	if (args[0] == "gridpoints") {
		gridPoints(args);
	} else if (args[0] == "create") {
		createSurface(args);
	} else if (args[0] == "evaluate") {
		evaluateSurface(args);
	} else if (args[0] == "error") {
		computeErrorMetric(args);
	} else {
		cout << "Unrecognized command" << endl;
	}
	
}

void executeScript(string filename) 
{
	cout << "Executing script: " << filename << endl;
	ifstream inscript(filename.c_str(), ios::in);
	if (!inscript) {
		cerr << "Error: file not found (" << filename << "). " << endl;
	} else {
		const int MAX_LINE = 200;
		char command[MAX_LINE];
		vector< string > args;
		while(!inscript.eof()) {
		    args.resize(0);
		    inscript.getline(command, MAX_LINE);
		    string commandstring = command;
		    if (commandstring != "\0") { // else it's an empty line; do nothing
		        istringstream incommand(commandstring);
		        while (!incommand.eof()) {
		            string nextarg;
		            incommand >> nextarg;
		            //cout << "Next arg: " << nextarg << endl;
		            args.push_back(nextarg);
		        }
		        //cout << "End of arguments " << endl;
		        if (args.size() > 0 && args[0][0] != '#') {
		                executeCommand(args);
		        }
		    }
		}
	}
}

void readEvalPrint()
{
	const int MAX_LINE = 200;
	//cout << "Read eval print" << endl;
	char command[MAX_LINE];
        vector< string > args;
	do {
	    args.resize(0);
	    cout << "surfpack> " << flush;
	    cin.getline(command, MAX_LINE);
	    string commandstring = command;
	    istringstream incommand(commandstring);
	    while (!incommand.eof()) {
                string nextarg;
		incommand >> nextarg;
		//cout << "Next arg: " << nextarg << endl;
		args.push_back(nextarg);
	    }
	    //cout << "End of arguments " << endl;
	    if (args.size() > 0 && args[0] != "quit") {
		    executeCommand(args);
	    }
	} while (! (args.size() > 0 && args[0] == "quit") );

	    
}

int main(int argc, char* argv[]) 
{
    if (argc == 1) {
	    readEvalPrint();
    } else if (argc == 2) {
	    executeScript(argv[1]);
    } else {
	    vector< string > args;
	    for (int i = 1; i < argc; i++) {
		    args.push_back(argv[i]);
	    }
	    executeCommand(args);
    }
    return 0;
}    
