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
#include "TGaxis.h"
#include "TPaletteAxis.h"
#include "stdio.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include <string>
#include <fstream>
#include <math.h>

float HodoPlaneShift(TTree* h4, std::string detector, std::string pathToOut, std::string RunStats, std::string axis);

void AmplitudeProfilesFit(std::string FileIn, std::string detector, Float_t bound, float* XMax, float* YMax)
{
	Int_t Nentries, i, CH, C3=0, C0APD1=0, C0APD2=0, runNum, hodoCfg;
	Float_t X[2], Y[2], Ampl[6];
	Float_t X1true, Y1true, Xshift, Yshift, Xavg, Yavg;
	Float_t Gain_val, Energy_val;

	float fitRange=5;
	gStyle->SetOptStat();
	gStyle->SetOptFit();

	std::string pathToOutput = "/afs/cern.ch/user/c/cquarant/www/Amplitude_profiles/";
	std::string pathToLinFitOutput = "/afs/cern.ch/user/c/cquarant/www/";

	TFile *f = TFile::Open(FileIn.c_str());
	TTree *h4 = (TTree*)f->Get("h4");
	
	h4->GetEntry(0);
	CH=h4->GetLeaf(detector.c_str())->GetValue(0);
	Nentries = h4->GetEntries();
	runNum = h4->GetLeaf("run")->GetValue(0);

	auto *AmpX0 = new TProfile("AmpX0", "", 32, -16, 16, 0, 10000);
	auto *AmpY0 = new TProfile("AmpY0", "", 32, -16, 16, 0, 10000);
	auto *AmpX1 = new TProfile("AmpX1", "", 32, -16, 16, 0, 10000);
	auto *AmpY1 = new TProfile("AmpY1", "", 32, -16, 16, 0, 10000);
	auto *AmpXavg = new TProfile("AmpXavg", "", 32, -16, 16, 0, 10000);
	auto *AmpYavg = new TProfile("AmpYavg", "", 32, -16, 16, 0, 10000);
	
	std::string Gain = std::to_string((int)h4->GetLeaf("CHGain")->GetValue(0));
	std::string Energy = std::to_string((int)h4->GetLeaf("Energy")->GetValue(0));
	std::string fileOutpdf;
	std::string fileOutpng;
	std::string RunStats = Energy+"Gev_G"+Gain+"_"+std::to_string(runNum);
	
	//Calculating shift of hodo plane 1 respect to hodo plane 0 faro X and Y axis
	Xshift = HodoPlaneShift(h4, detector, pathToLinFitOutput, RunStats, "X");
	Yshift = HodoPlaneShift(h4, detector, pathToLinFitOutput, RunStats, "Y");

	//Filling Amplitude profile histograms Draw method
	std::string bound_str = std::to_string(bound);
	std::string Xshift_str = std::to_string(Xshift);
	std::string Yshift_str = std::to_string(Yshift);

	std::string varexp = "amp_max[" + detector + "]:X[0]>>AmpX0";
	std::string selection = "Y[0]>-" + bound_str + " && Y[0]<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());
	
	varexp = "amp_max[" + detector + "]:Y[0]>>AmpY0"; 
	selection = "X[0]>-" + bound_str + " && X[0]<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:(X[1]" + Xshift_str + ")>>AmpX1";
	selection = "Y[1]" + Yshift_str + ">-" + bound_str + " && Y[1]" + Yshift_str + "<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:(Y[1]" + Yshift_str + ")>>AmpY1";
	selection = "X[1]" + Xshift_str + ">-" + bound_str + " && X[1]" + Xshift_str + "<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:0.5*(X[0]+X[1]" + Xshift_str + ")>>AmpXavg";
	selection = "0.5*(Y[0]+Y[1]" + Yshift_str + ")>-" + bound_str + " && 0.5*(Y[0]+Y[1]" + Yshift_str + ")<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());

	varexp = "amp_max[" + detector + "]:0.5*(Y[0]+Y[1]" + Yshift_str + ")>>AmpYavg";
	selection = "0.5*(X[0]+X[1]" + Xshift_str + ")>-" + bound_str + " && 0.5*(X[0]+X[1]" + Xshift_str + ")<" + bound_str;
	h4->Draw(varexp.c_str(), selection.c_str());
	

	//Drawing and fitting X0 histogram
	TCanvas* c1 = new TCanvas("c1","c1");
 	TH2F* H1 = new TH2F("H1","", 32, -16, 16, 50, 0, (AmpX0->GetMaximum())*1.25);
	H1->GetXaxis()->SetTitle("X[0]");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw("AXIS");
	AmpX0->Fit("pol2", "", "", -fitRange, fitRange);
	AmpX0->Draw("SAME");

	fileOutpdf = pathToOutput + "p0/" + Energy + "Gev/" + "AmpX0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p0/" + Energy + "Gev/" + "AmpX0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  	c1->SaveAs(fileOutpdf.c_str());
  	c1->SaveAs(fileOutpng.c_str());
	
	
	//Drawing and fitting Y0 histogram
	TCanvas* c2 = new TCanvas("c2","c2");
 	H1->GetXaxis()->SetTitle("Y[0]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw("AXIS");	
	AmpY0->Fit("pol2", "", "", -fitRange, fitRange);	
	AmpY0->Draw("SAME");

	fileOutpdf = pathToOutput + "p0/" + Energy + "Gev/" + "AmpY0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p0/" + Energy + "Gev/" + "AmpY0_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c2->SaveAs(fileOutpdf.c_str());
  	c2->SaveAs(fileOutpng.c_str());


	//Drawing X1 histogram
	TCanvas* c3 = new TCanvas("c3","c3");
 	H1->GetXaxis()->SetTitle("X[1]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpX1->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpX1->Draw("SAME");
  
	fileOutpdf = pathToOutput + "p1/" + Energy + "Gev/" + "AmpX1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p1/" + Energy + "Gev/" + "AmpX1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c3->SaveAs(fileOutpdf.c_str());
  	c3->SaveAs(fileOutpng.c_str());

	//Drawing Y1 histogram
	TCanvas* c4 = new TCanvas("c4","c4");
 	H1->GetXaxis()->SetTitle("Y[1]");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpY1->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpY1->Draw("SAME");

	fileOutpdf = pathToOutput + "p1/" + Energy + "Gev/" + "AmpY1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "p1/" + Energy + "Gev/" + "AmpY1_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c4->SaveAs(fileOutpdf.c_str());
  	c4->SaveAs(fileOutpng.c_str());

	//Drawing Xavg histogram
	TCanvas* c5 = new TCanvas("c5","c5");
 	H1->GetXaxis()->SetTitle("Xavg");    
 	H1->GetYaxis()->SetTitle("amp_max");
  
	H1->Draw();
	AmpXavg->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpXavg->Draw("SAME");

	TF1 *fitResX = AmpXavg->GetFunction("pol2");
  
	fileOutpdf = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpXAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpXAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";
  
  
  	c5->SaveAs(fileOutpdf.c_str());
  	c5->SaveAs(fileOutpng.c_str());
	
	//Drawing Yavg histogram
	TCanvas* c6 = new TCanvas("c6","c6");
 	H1->GetXaxis()->SetTitle("Yavg");    
 	H1->GetYaxis()->SetTitle("amp_max");

	H1->Draw();  
	AmpYavg->Fit("pol2", "", "", -fitRange, fitRange);		
	AmpYavg->Draw("SAME");
	
	TF1* fitResY = AmpYavg->GetFunction("pol2");
	
  	fileOutpdf = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpYAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".pdf";
  	fileOutpng = pathToOutput + "pAVG/" + Energy + "Gev/" + "AmpYAVG_" + std::to_string(runNum) + "_" + detector + "_G" + Gain + ".png";

    
  	c6->SaveAs(fileOutpdf.c_str());
  	c6->SaveAs(fileOutpng.c_str());

	*XMax = fitResX->GetMaximumX();
	*YMax = fitResY->GetMaximumX();
}

float HodoPlaneShift(TTree* h4, std::string detector, std::string pathToOut, std::string RunStats, std::string axis)
{
	auto *DXvsX = new TProfile(("D"+axis+"vs"+axis+"").c_str(), "", 128, -16, 16, -10, 10);

	//Filling DeltaX histogram
	h4->Draw(("("+axis+"[0]-"+axis+"[1]):"+axis+"[0]>>D"+axis+"vs"+axis+"").c_str(), (axis+"[0]>-800 && "+axis+"[1]>-800").c_str());

	//Fitting and Drawing DeltaX
	TCanvas *cDeltaX = new TCanvas(("cDelta"+axis+"").c_str(), ("cDelta"+axis+"").c_str());
	TH2F* Hset = new TH2F(("Hset"+axis).c_str(),"", 128, -16, 16, 100, -5, 5);     
	Hset->GetXaxis()->SetTitle((axis+"[0]").c_str());    
 	Hset->GetYaxis()->SetTitle((axis+"[0]-"+axis+"[1]").c_str());
  
	Hset->Draw();	
	DXvsX->Fit("pol1", "Q", "", -4, 4);
	DXvsX->Draw("SAME");
	
	std::string fileOutpdf = pathToOut + "Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".pdf";
  	std::string fileOutpng = pathToOut + "Delta"+axis+"/D"+axis+"vs"+axis+"_" + detector + "_" + RunStats + ".png";
    
  	cDeltaX->SaveAs(fileOutpdf.c_str());
  	cDeltaX->SaveAs(fileOutpng.c_str());

	return DXvsX->GetFunction("pol1")->GetParameter(0);
}
