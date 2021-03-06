#include "TLeaf.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "Mylib.h"
#include <string>
#include <fstream>
#include <math.h>

void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats, std::string runNum, std::string MCP);
float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOut, std::string RunStats, std::string MCP);
float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOut, std::string RunStats, std::string MCP);

void AveragePulseShape(std::string FileIn, std::string detector, Float_t bound, std::string runNum, std::string MCP)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, hodoCfg;
	Float_t X[2], Y[2], Ampl[6];
	Float_t X1true, Y1true, Xshift, Yshift, Xavg, Yavg;
	Float_t Gain_val, Energy_val;
	float XMax, YMax;

	float fitRange=5;
	//setStyle();
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/";

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();

	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string RunStats = "E"+Energy+"_G"+Gain;

	AmplitudeProfilesFit(FileIn, detector, bound, &XMax, &YMax);
	
	if(fabs(XMax)<100 && fabs(YMax)<100)
	{
		PulseShapes(h4, detector, 0, XMax, YMax, 2, pathToOutput, RunStats, runNum, MCP);
	}
	else cout<< "SBALLATO!" << endl;
	
}


void PulseShapes(TTree* h4, std::string detector, int plane, float XMax, float YMax, float range, std::string pathToOutput, std::string RunStats, std::string runNum, std::string MCP)
{
	int i;
	float AmpMean, AmpSigma;

	std::string TimeShift;
	std::string TimeMCP;
	std::string AmpMean_str;
	std::string AmpSigma_str;

	TProfile2D* p2D_amp_vs_time = new TProfile2D("p2D_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	TH2F* h2_amp_vs_time = new TH2F("h2_amp_vs_time", "", 300, -10, 40, 300, -0.5, 1.5);

	std::string Selection = "fabs(X[" + std::to_string(plane) + "]-(" + std::to_string(XMax) + "))<" + std::to_string(range) + " && fabs(Y[" + std::to_string(plane) + "]-(" + std::to_string(YMax) + "))<" + std::to_string(range) + " && amp_max["+MCP+"]>100";	

	TimeMCP = std::to_string(MeanTimeMCP(h4, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = Selection + " && fabs(time["+MCP+"]-("+TimeMCP+"))<1";

	Selection = Selection + " && fabs(time["+detector+"]-time["+MCP+"])<7";

	AmplitudeHist(h4, detector, Selection, pathToOutput, RunStats, &AmpMean, &AmpSigma);
	AmpMean_str = std::to_string(AmpMean);
	AmpSigma_str = std::to_string(AmpSigma);	
	Selection = Selection + " && fabs(amp_max["+detector+"]-("+AmpMean_str+"))<5*"+AmpSigma_str;
	TimeShift = std::to_string(MeanTimeShift(h4, detector, Selection, pathToOutput+"PulseShapes/", RunStats, MCP));
	Selection = "WF_ch == " + detector + " && " + Selection;	
	

	cout << Selection << endl;
        h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+")>>p2D_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 1" << endl;
	
	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time["+MCP+"]-("+TimeShift+") >> h2_amp_vs_time").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

	TObjArray aSlices;
	h2_amp_vs_time->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices);
	
	TProfile *waveForm = new TProfile(("XTAL_"+detector+"_"+RunStats+"_prof_").c_str(), "", 300, -10, 40);
	waveForm = (TProfile*)aSlices[1];

	p2D_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
    	h2_amp_vs_time->GetXaxis()->SetTitle((std::string("WF_time-time["+MCP+"] (ns)")).c_str());
    	waveForm->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	p2D_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time->GetZaxis()->SetTitle("amp_max");	
    	    	
	gStyle->SetOptStat(0);
	
    	TCanvas* c1 = new TCanvas();
    	c1->cd();
    	h2_amp_vs_time->Draw("COLZ");
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_h2.png").c_str());
    	c1 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+detector+"_h2.pdf").c_str());
	
    	TCanvas* c2 = new TCanvas();
    	c2->cd();
    	p2D_amp_vs_time->Draw("COLZ");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.png").c_str(),"png");
    	c2 -> Print(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_Amp.pdf").c_str(),"pdf");
	
	TCanvas* c0 = new TCanvas();
    	c0->cd();
	waveForm->GetYaxis()->SetRangeUser(0, 1.05);
	waveForm->SetName(("XTAL_"+detector+"_"+RunStats+"_prof_").c_str());
    	waveForm->Draw();
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.png").c_str());
    	c0 -> SaveAs(std::string(pathToOutput+"PulseShapes/AllCuts/PS_"+detector+"_"+RunStats+"_"+runNum+"_profile.pdf").c_str());    	

    	TFile* output_Waveform = new TFile(std::string("WaveForms/"+detector+"_"+RunStats+"_Waveform.root").c_str(),"RECREATE");
    	output_Waveform->cd();
	
	TTree *TRunStats = new TTree("TRunStats", "TRunStats");

	TBranch *RUN_NUM = new TBranch(TRunStats, "RUNSTATS", &RunStats, "RUNSTATS/C");

	p2D_amp_vs_time->SetName("Profile2DWaveforms_");
	h2_amp_vs_time->SetName("H2Waveforms_");

	p2D_amp_vs_time->Write();
    	h2_amp_vs_time->Write();
    	waveForm->Write();
	TRunStats->Write();

    	output_Waveform->Close();
	cout << waveForm->GetName() << endl; 
}

//returns the mean time shift between detector and MCP selected
float MeanTimeShift(TTree* h4, std::string detector, std::string Selection, std::string pathToOut, std::string RunStats, std::string MCP)
{
	TH1F* raw_time_dist = new TH1F("raw_time_dist", "", 400, -6, 6);

	h4->Draw((std::string("time["+detector+"]-time["+MCP+"]>>raw_time_dist")).c_str(), Selection.c_str());

	raw_time_dist->GetXaxis()->SetTitle((std::string("time["+detector+"]-time["+MCP+"] (ns)")).c_str());
	raw_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* c0 = new TCanvas();
    	c0->cd();
	raw_time_dist->Fit("gaus", "", "", -6, 6);
    	raw_time_dist->Draw();
    	c0 -> SaveAs(std::string(pathToOut+"RawTimeDistribution/RawTimeDist"+MCP+"_"+detector+"_"+RunStats+".png").c_str());
    	c0 -> SaveAs(std::string(pathToOut+"RawTimeDistribution/RawTimeDist"+MCP+"_"+detector+"_"+RunStats+".pdf").c_str());

	cout << raw_time_dist->GetMean() << endl;
	cout << raw_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return raw_time_dist->GetFunction("gaus")->GetMaximumX();
}

//returns the mean time measured by MCP selected
float MeanTimeMCP(TTree* h4, std::string Selection, std::string pathToOut, std::string RunStats, std::string MCP)
{
	TH1F* MCP_time_dist = new TH1F("MCP_time_dist", "", 200, 0, 50);

	h4->Draw((std::string("time["+MCP+"]>>MCP_time_dist")).c_str(), Selection.c_str());
	
	MCP_time_dist->GetXaxis()->SetTitle(("time["+MCP+"] (ns)").c_str());
	MCP_time_dist->GetYaxis()->SetTitle("events");

	TCanvas* ca = new TCanvas();
    	ca->cd();
	MCP_time_dist->Fit("gaus", "", "", 0, 50);
    	MCP_time_dist->Draw();
    	ca -> SaveAs(std::string(pathToOut+"MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".png").c_str());
    	ca -> SaveAs(std::string(pathToOut+"MCPTimeDistribution/"+MCP+"_TimeDist_"+RunStats+".pdf").c_str());

	cout << MCP_time_dist->GetMean() << endl;
	cout << MCP_time_dist->GetFunction("gaus")->GetMaximumX() << endl;

	return MCP_time_dist->GetFunction("gaus")->GetMaximumX();
}


//TProfile2D* p2D_amp_vs_time_no_cut = new TProfile2D("p2D_amp_vs_time_no_cut", "", 300, -10, 40, 300, -0.5, 1.5, 0., 10000.);
	//TH2F* h2_amp_vs_time_no_cut = new TH2F("h2_amp_vs_time_no_cut", "", 300, -10, 40, 300, -0.5, 1.5);

/*
	//Drawing WITHOUT CUTS on HODOSCOPE
	Selection = "amp_max[MCP1]>200 && fabs(time[MCP1]-("+TimeMCP+"))<8 && X[0]>-800 && Y[0]>-800";

	h4->Draw(("amp_max["+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+")>>p2D_amp_vs_time_no_cut").c_str(),Selection.c_str());
        //std::cout << (std::string("amp_max[")+detector+"]:WF_val/amp_max["+detector+"]:WF_time-time[MCP1]-("+TimeShift+") >> p2D_amp_vs_time") << std::endl;	
	//cout << Selection << endl;
	cout << "Draw 1" << endl;

	h4->Draw(("WF_val/amp_max["+detector+"]:WF_time-time[MCP1]"+TimeShift+" >> h2_amp_vs_time_no_cut").c_str(),Selection.c_str());
	cout << "Draw 2" << endl;

    	TProfile* waveForm_no_cut = h2_amp_vs_time_no_cut->ProfileX(); 
    	waveForm_no_cut->SetName(std::string(detector+std::string("_waveform_prof_no_cut")).c_str());

	p2D_amp_vs_time_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	h2_amp_vs_time_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	waveForm_no_cut->GetXaxis()->SetTitle((std::string("WF_time-time[MCP1] (ns)")).c_str());
    	p2D_amp_vs_time_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	h2_amp_vs_time_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	waveForm_no_cut->GetYaxis()->SetTitle((std::string("WF_val/amp_max[")+detector+std::string("]")).c_str());
    	p2D_amp_vs_time_no_cut->GetZaxis()->SetTitle("amp_max");
   	h2_amp_vs_time_no_cut->GetZaxis()->SetTitle("amp_max");	
	*/

/*
	TCanvas* c4 = new TCanvas();
    	c5->cd();
    	h2_amp_vs_time_no_cut->Draw("COLZ");
    	c5 -> SaveAs(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_h2_no_cut.png").c_str());
    	c5 -> SaveAs(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_h2_no_cut.pdf").c_str());

    	TCanvas* c5 = new TCanvas();
    	c6->cd();
    	p2D_amp_vs_time_no_cut->Draw("COLZ");
    	c6 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_Amp_no_cut.png").c_str(),"png");
    	c6 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_Amp_no_cut.pdf").c_str(),"pdf");

    	TCanvas* c6 = new TCanvas();
    	c7->cd();
    	waveForm_no_cut->Draw("P");
    	c7 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_profile_no_cut.png").c_str(),"png");
    	c7 -> Print(std::string(pathToOutput+"PulseShapes/NoCutAmplitude/PS_"+detector+"_"+RunStats+"_profile_no_cut.pdf").c_str(),"pdf");
	*/



