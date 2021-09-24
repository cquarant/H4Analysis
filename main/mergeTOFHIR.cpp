#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <experimental/filesystem>

#include "TFile.h"
#include "TChain.h"

#include "interface/PositionTree.h"

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"


using namespace std;

//**********MAIN**************************************************************************
int main(int argc, char* argv[])
{
    if(argc < 2)
    {
        cout << argv[0] << " cfg file " << "[run] " << "[spill] " <<endl; 
        return -1;
    }

    //---load options---    
    CfgManager opts;
    opts.ParseConfigFile(argv[1]);

    //-----input setup-----    
    int spill=-1;
    if(argc > 2)
    {
        vector<string> run(1, argv[2]);
        opts.SetOpt("merge.run", run);
    }

    if(argc > 3)
        spill = atoi(argv[3]);

    auto tofhir_path = opts.GetOpt<string>("merge.tofhirPath");
    auto h4daq_path = opts.GetOpt<string>("merge.h4daqPath");
    auto output_path = opts.GetOpt<string>("merge.outputPath");

    string run = opts.GetOpt<string>("merge.run");
    
    //TOFHIR input
    TFile* inputTOFHIR=TFile::Open(Form("%s/%s/%d_ped_e.root",tofhir_path.c_str(),run.c_str(),spill),"UPDATE");
    TTree* data=(TTree*)inputTOFHIR->Get("data");
    int channelIdx[128];
    vector<long int>* times = 0;
    data->SetBranchAddress("channelIdx",channelIdx);
    data->SetBranchAddress("time",&times);

    TTree* outTree=new TTree("h4","h4");
    int iev_TOFHIR;
    int iev_H4DAQ;
    unsigned long t_TOFHIR;
    unsigned long t_H4DAQ;
    float x_WC;
    float y_WC;
    int nhits_WC;
    int nclusters_WC;

    outTree->Branch("iev_TOFHIR",&iev_TOFHIR,"iev_TOFHIR/I");
    outTree->Branch("iev_H4DAQ",&iev_H4DAQ,"iev_H4DAQ/I");
    outTree->Branch("t_TOFHIR",&t_TOFHIR,"t_TOFHIR/l");
    outTree->Branch("t_H4DAQ",&t_H4DAQ,"t_H4DAQ/l");
    outTree->Branch("nclusters_WC",&nclusters_WC,"nclusters_WC/I");
    outTree->Branch("nhits_WC",&nhits_WC,"nhits_WC/I");
    outTree->Branch("x_WC",&x_WC,"x_WC/F");
    outTree->Branch("y_WC",&y_WC,"y_WC/F");
    
    //H4DAQ input
    TFile* inputH4DAQ=TFile::Open(Form("%s/%s.root",h4daq_path.c_str(),run.c_str()));
    TTree* h4=(TTree*)inputH4DAQ->Get("h4");
    vector<unsigned long int>* timesH4DAQ = 0;
    vector<PositionMeasurement>* wc_measurements = 0;
    unsigned int runH4DAQ;
    unsigned int spillH4DAQ;
    unsigned int eventH4DAQ;

    h4->SetBranchAddress("time_stamps",&timesH4DAQ);
    h4->SetBranchAddress("run",&runH4DAQ);
    h4->SetBranchAddress("spill",&spillH4DAQ);
    h4->SetBranchAddress("event",&eventH4DAQ);
    h4->SetBranchAddress("clusters",&wc_measurements);

    unsigned long int nEntriesTOFHIR=data->GetEntries();
    unsigned long int iEventTOFHIR=0;

    unsigned long int nEntriesH4DAQ=h4->GetEntries();
    unsigned long int iEventH4DAQ=0;

    unsigned long int firstTriggerTimeTOFHIR=0;
    unsigned long int firstTriggerTimeH4DAQ=0;

    long int timeTOFHIR=0;
    long int timeH4DAQ=0;

    bool matchFound=0;
    bool spillFound=0;
    int maxSearchWindow=10;
    float matchTimeWindow=150;

    while(iEventTOFHIR<nEntriesTOFHIR)
    //while(iEventTOFHIR<11000)
      {
	matchFound=0;
	iev_TOFHIR=iEventTOFHIR;
	iev_H4DAQ=-1;
	x_WC=-999;
	y_WC=-999;
	nclusters_WC=-1;
	nhits_WC=-1;
	data->GetEntry(iEventTOFHIR++);

	//skip events without trigger
	if (! (channelIdx[96]>-1) )
	  {
	    outTree->Fill();
	    continue;
	  }

	//first event with trigger for TOFHIR
	if (firstTriggerTimeTOFHIR==0)
	  {
	    firstTriggerTimeTOFHIR=times->at(channelIdx[96]);
	    std::cout << "First trigger TOFHIR event " << iEventTOFHIR << " time " << firstTriggerTimeTOFHIR << std::endl; 
	  }
	
	timeTOFHIR=(times->at(channelIdx[96])-firstTriggerTimeTOFHIR)/1E6; //tofhir time in ps, convert to mus
	t_TOFHIR=timeTOFHIR;
	t_H4DAQ=0;
	int searchIndex=0;

	//matching loop
	while(!matchFound && ( (spillFound&&searchIndex<maxSearchWindow) || !spillFound) )
	  {
	    if ((iEventH4DAQ+searchIndex)<nEntriesH4DAQ)
	      {
		h4->GetEntry(iEventH4DAQ+searchIndex);
		searchIndex++;
	      }
	    else
	      {
		std::cout << "Spill not found  " << spill << std::endl;
		goto theEnd;
	      }

	    //skip events not in this spill
	    if (spill!=spillH4DAQ)
	      continue;

	    //spillFound
	    if (!spillFound)
	      {
		std::cout << "Spill found " << spill << std::endl;
		spillFound=true;
	      }

	    //first event in spill
	    if (firstTriggerTimeH4DAQ==0)
	      {
		firstTriggerTimeH4DAQ=timesH4DAQ->at(0);
		std::cout << "First trigger H4DAQ event " << iEventH4DAQ+searchIndex << " time " << firstTriggerTimeH4DAQ << std::endl; 
	      }

	    timeH4DAQ=timesH4DAQ->at(0) - firstTriggerTimeH4DAQ; //h4daq time already in mus 
	    std::cout << "H4DAQ event " << iEventH4DAQ+searchIndex << " time " << timeH4DAQ << std::endl;
	    if (abs(timeH4DAQ-timeTOFHIR+70*(timeH4DAQ/4300000)+50)<matchTimeWindow)
	      {
		std::cout << "Found match TOFHIR event " << iEventTOFHIR << " H4DAQ event " << iEventH4DAQ+searchIndex << " deltaT " << timeH4DAQ-timeTOFHIR+70*(timeH4DAQ/4300000)+50 << " x " << wc_measurements->at(0).X() << " y " << wc_measurements->at(0).Y() << std::endl;
		matchFound=true;
		iev_H4DAQ=iEventH4DAQ+searchIndex-1;
		t_H4DAQ=timeH4DAQ;
		nclusters_WC=wc_measurements->size();
		if (nclusters_WC>0)
		  {
		    nhits_WC=wc_measurements->at(0).nHits();
		    x_WC=wc_measurements->at(0).X();
		    y_WC=wc_measurements->at(0).Y();
		  }
	      }
	  } //end-loop

	outTree->Fill();
	//resetting search position
	if (matchFound)
	  iEventH4DAQ+=searchIndex;
	
      }

 theEnd:
    //---close
    // outROOT->cd();
    inputTOFHIR->cd();
    outTree->Write();
    data->AddFriend(outTree);
    data->Write();
    inputTOFHIR->Close();
    // outROOT->Close();
    exit(0);
}    


