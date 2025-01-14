/*
  Main.cc
*/

#include <iostream>

#include <string>
#include <vector>
#include <signal.h>

#include "ConfMan.hh"
#include "UnpackerManager.hh"
#include "VEvent.hh"

#include <TROOT.h>
#include <TFile.h>

enum eArg
  {
    kArgProcessName,
    kArgConfFile,
    kArgInFile,
    kArgOutRootFile,
    kArgc
  };


void closeFile(int sig)
{
  if(gFile)
    {
      gFile->Write();
      gFile->Close();
    }
}


TROOT theROOT("sample-ana", "test");

int main(int argc, char* argv[])
{
  if (argc!=kArgc)
    {
      std::cout <<" #D Usage: "
		<<" [analyzer config file]"
		<<" [data input stream] "
		<<" [output root file]"
		<< std::endl;
      return 0;
    }

  std::vector<std::string> arg(argv, argv + argc);

  const std::string& confFile = arg[kArgConfFile];
  const std::string& inFile   = arg[kArgInFile];
  const std::string& rootFile = arg[kArgOutRootFile];


  //  ConfMan& gConfManager = ConfMan::GetInstance();
  ConfMan* gconfManager =new ConfMan(confFile);

  signal(SIGINT,closeFile);

  TFile f(rootFile.c_str(), "recreate");
  std::cout << "#D recreate root file : "
	    << rootFile << std::endl;

  if (!gconfManager->Initialize())
    return 0;
  hddaq::unpacker::UnpackerManager& gUnpacker 
    = hddaq::unpacker::GUnpacker::get_instance();

  if (gconfManager->GetEvDispFlag()) {
    std::cout << "Create Event Display" << std::endl;
    gconfManager->InitializeEvDisp();
  }

  gUnpacker.set_istream(inFile);
  gUnpacker.initialize();

  for(;!gUnpacker.eof();++gUnpacker){
    VEvent* event = gconfManager->EventAllocator();
    event->ProcessingBegin();
    event->ProcessingNormal();
    event->ProcessingEnd();
    delete event;
  }
  gFile->Write();
  gFile->Close();
  
  return 0;
}

