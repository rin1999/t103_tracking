/*
  Main.cc

  2012/5  K.Shirotori
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

#include "Main.hh"
#include "ConfMan.hh"
#include "VEvent.hh"

#include <TROOT.h>
#include <TFile.h>

const int MaxChar = 400;

void closeFile(int sig)
{
  if(oFile)
    {
      oFile->Write();
      oFile->Close();
    }
}

TROOT theROOT("sample-ana", "test");
TFile* oFile;

int main(int argc, char* argv[])
{
  if (argc!=kArgc)
    {
      std::cout <<" #D Usage:"
		<<" [analyzer config file]"
		<<" [input root file]"
		<<" [output root file]"
		<< std::endl;
      return 0;
    }

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  // std::ifstream InputData(argv[kArgInFile]);
  const std::string& iFile_name = arg[kArgInRootFile];
  const std::string& oFile_name = arg[kArgOutRootFile];

  ConfMan* gconfManager =new ConfMan( confFile );

  signal(SIGINT,closeFile);

  TFile* iFile = new TFile(iFile_name.c_str(), "READ");
  if( !iFile ) return -1;
  std::cout << "#D read root file : " << iFile << std::endl;


  oFile = new TFile(oFile_name.c_str(), "recreate");
  std::cout << "#D recreate root file : " << oFile << std::endl;

  if (!gconfManager->Initialize())
    return -1;

  VEvent* event = gconfManager->EventAllocator();
  if ( event->ProcessingNormal( iFile ) )
      delete event;

  oFile->Write();
  oFile->Close();
  
  return 0;
}

