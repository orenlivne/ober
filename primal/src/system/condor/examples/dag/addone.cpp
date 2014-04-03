//=============================================================
// An example C++ program that reads numbers from an input
// file and writes an output file with all numbers, incremented
// by 1.
//=============================================================
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <climits>
#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include  <time.h>
using namespace std;

//==============================================================
// Input arguments
//==============================================================
struct globalArgs_t
{
  int failure;                /* -f option - induce intermittent failures */
  char *inputFile;            /* input file name */
  char *outputFile;           /* output file name */
} globalArgs;

static const char *optString = "f";

static const struct option longOpts[] =
  {
    { "failure", no_argument, NULL, 'f' },
    { "help", no_argument, NULL, 'h' },
    { NULL, no_argument, NULL, 0 }
  };

void displayUsage()
{
  cerr << "Usage: addone [options] <infile> <outfile>" << endl << endl;
  cerr << "Reads numbers from the file infile, increments each by 1," << endl;
  cerr << "and writes the result to the file outfile." << endl << endl;
  cerr << "Options:" << endl;
  cerr << "\t-f\tInduce intermittent failures" << endl;
  exit(-1);
}

void parseArgs(int argc, char* const argv[])
{
  if (argc < 3)
    displayUsage();

  /* Initialize globalArgs before we get to work. */
  int opt = 0;
  globalArgs.failure = 0;     /* false */
  globalArgs.inputFile = NULL;
  globalArgs.outputFile = NULL;
  
  int longIndex;
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while ( opt != -1 )
    {
      switch ( opt )
	{
	case 'f':
	  globalArgs.failure = 1; /* true */
	  break;
	case 'h':   /* fall-through is intentional */
	case '?':
	  displayUsage();
	  break;
	default:
	  /* You won't actually get here. */
	  break;
	}
      opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
    }
  globalArgs.inputFile  = *(argv + optind);
  globalArgs.outputFile = *(argv + optind + 1);
}

//==============================================================
// Helper functions
//==============================================================

int str2int(const char* const s)
// Convert a string to an integer.
{
  int i;
  int correctly_processed = sscanf(s, "%d", &i);
  if (correctly_processed < 1)
    throw "Failed to parse integer " + (*s);
  return i;
}

//==============================================================
// Main function
//==============================================================
int main(int argc, char* const argv[])
{
  // Read and validate arguments
  parseArgs(argc, argv);

  char* const inFileName = globalArgs.inputFile;
  ifstream inFile(inFileName);
  if (!inFile.is_open())
    {
      cerr << "Unable to open input file " << inFileName << endl;
      exit(-1);
    }

  char* const outFileName = globalArgs.outputFile;
  ofstream outFile(outFileName);
  if (!outFile.is_open())
    {
      cerr << "Unable to open output file " << outFileName << endl;
      exit(-1);
    }
  
  srand ( time(NULL) );
  int fail = ((rand() % 3) == 0);
  if (globalArgs.failure && fail)
    {
      cerr << "Intermittent failure" << endl;
      inFile.close();
      outFile.close();
      exit(-1);
    }

  // Read lines from input file
  while (inFile.good())
    {
      // Read line. Don't grab the last line in the file if it's empty
      string line;
      getline (inFile, line);
      if (inFile.eof())
	break;
      
      // Convert to integer
      int i = str2int(line.c_str());
      
      // Write to output file
      outFile << (i+1) << endl;
    }
  inFile.close();
  outFile.close();
  
  return 0;
}
