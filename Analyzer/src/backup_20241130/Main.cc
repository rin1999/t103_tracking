/*
  Main.cc

  2024/04 K.Shirotori
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#include <string>
#include <vector>
#include <signal.h>

#include "ConfMan.hh"
#include "VEvent.hh"

#include <TROOT.h>
#include <TFile.h>

#include <Compression.h>

extern TFile* oFile;

//if this is 0, all events will be analyzed.
unsigned int Nanalysis=-1;

enum ArgName{
  exec, 
  conffile, 
  inputfile, 
  outputfile,
  nanalysis
};

bool find_opt(int argc, char* argv[], const std::string& item)
{
  bool ret=false;
  std::vector<std::string> v(argv, argv+argc);
  for (unsigned int i=0; i<v.size(); ++i) {
    if (v[i] == item) {
      ret = true;
      break;
    }
  }
  return ret;
}

const int MaxChar = 400;

void closeFile(int sig)
{
  if(oFile){
    oFile->Write();
    oFile->Close();
  }
}

TROOT theROOT("sample-ana", "test");
TFile* oFile;

int main(int argc, char* argv[])
{
  if(argc<=3){
    std::cout <<" #D Usage: "
	      <<" [analyzer config file]"
	      <<" [input root file]"
	      <<" [output root file]"
	      <<" [number of analyzed events (optional)]"
	      << std::endl;
    return 0;
  } 

  std::vector<std::string> arg(argv, argv + argc);
  const std::string& confFile = arg[conffile];
  const std::string& dataFile = arg[inputfile];
  const std::string& rootFile = arg[outputfile];

  if(argc == 5){
    Nanalysis = atoi(argv[nanalysis]);
    std::cout << Nanalysis  << " events will be analyzed" << std::endl;
  }

  ConfMan* gconfManager =new ConfMan( confFile );

  signal(SIGINT,closeFile);

  TFile* iFile = new TFile(dataFile.c_str(), "READ");
  if( !iFile ) return -1;
  std::cout << "#D read root file : " << iFile << std::endl;

  oFile = new TFile(rootFile.c_str(), "recreate");
  std::cout << "#D recreate root file : " << oFile << std::endl;
  
  if (!gconfManager->Initialize())
     return -1;

  VEvent* event = gconfManager->EventAllocator();
  if ( event->ProcessingNormal( iFile, Nanalysis ) )
     delete event;

  oFile->Write();
  oFile->Close();

  //std::cout << "# of analyzed events: " << std::dec << evNum << std::endl;

  return 0;
}

