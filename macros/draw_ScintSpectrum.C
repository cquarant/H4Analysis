#include "TFile.h"
#include "TStyle.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h" 
#include "TCanvas.h"
#include "TLegend.h"
#include "TColor.h"
#include "TLatex.h"
#include "TGraphAsymmErrors.h"
#include <iostream>
#include "FPCanvasStyle.C"
#include "setStyle.C"

#include<iostream>
#include<string>
#include<fstream>

void draw_ScintSpectrum()
{
    gStyle->SetOptTitle(0); 
    //gStyle->SetOptStat(1110); 
    gStyle->SetOptStat(0000); 
    //gStyle->SetOptFit(1); 
    gStyle->SetOptFit(0); 
    gStyle->SetErrorX(0);

    
    TFile* inputFile = TFile::Open("root://eoscms.cern.ch//store/group/dpg_ecal/comm_ecal/upgrade/testbeam/TimingTB_BTF_Jun2016/ntuples/v5/ntuples_v5/btf2016_RU5_2378.root");
     
    TTree* h4 = (TTree*)inputFile->Get("h4");
    TH1F* h_amp = new TH1F("h_amp","",225,0.,4500.);
    h4->Draw("adc_data[scint] >> h_amp");  
    h_amp->SetLineColor(kCyan+1);
    h_amp->SetFillColor(kCyan+1);
    
    TH1F* h_amp_sel = new TH1F("h_amp_sel","",225,0.,4500.);
    h4->Draw("adc_data[scint] >> h_amp_sel","adc_data[scint]>200. && adc_data[scint] <700.");  
    h_amp_sel->SetLineColor(kYellow+1);
    h_amp_sel->SetFillColor(kYellow+1);
    //h_amp_sel->SetFillStyle(3003);

    float norm = 1./h_amp->Integral();

    h_amp->Scale(norm);
    h_amp_sel->Scale(norm);

    setStyle();  

    TH2F* H2 = new TH2F("H2","",225,0.,4500.,100,0.00001,1.);
    H2->GetXaxis()->SetTitle("amplitude (ADC counts)");
    H2->GetYaxis()->SetTitle("a.u.");
    H2->GetXaxis()->SetLabelSize(0.04);
    //H2->GetYaxis()->SetLabelSize(0.04);

    TLegend* legend = new TLegend(0.55, 0.65, 0.65, 0.82);
    legend -> SetFillColor(kWhite);
    legend -> SetFillStyle(1000);
    legend -> SetLineWidth(0);
    legend -> SetLineColor(kWhite);
    legend -> SetTextFont(42);  
    legend -> SetTextSize(0.04);
    
    legend -> AddEntry(h_amp,"Scintillator","F");
    legend -> AddEntry(h_amp_sel,"Selected range","F");
    
    TCanvas* c1 = new TCanvas();
    FPCanvasStyle(c1);
    c1->SetLogy();
    H2->Draw();
    h_amp->Draw("H,same");
    h_amp_sel->Draw("HF,same");
    c1->RedrawAxis("sameaxis");
    legend->Draw("same");
    TLatex latex2(0.65, 0.94,"#bf{#bf{Electrons at 491 MeV}}");;
    latex2.SetTextSize(0.04);
    latex2.SetNDC(kTRUE);
    latex2.Draw(); 
    c1 -> Print("Scint_spectrum.png","png");
    c1 -> Print("Scint_spectrum.pdf","pdf");  
}
