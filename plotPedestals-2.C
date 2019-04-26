#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "stdio.h"
#include "stdlib.h"

#include "TAxis.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLine.h"
#include "TLegend.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TObject.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TTree.h"

using namespace std;

static const int NSAMPLE = 32;
static const int NCHAN = 64;
static const int NMOD = 4;

// int eventNumber;
// int bclk[NMOD];
// int chanADC[NCHAN*NMOD][NSAMPLE];
// float mean[NCHAN*NMOD];
// float rms[NCHAN*NMOD];
// int NSamples;

TFile *file[10];
TFile *file1, *file2;

float mean[NCHAN*NMOD];
float rms[NCHAN*NMOD];
int chanADC[NCHAN*NMOD][NSAMPLE];
int NSamples;


TFile *openTree(TString rootFileName) { 
  return (new TFile(rootFileName)); 
}

void comparePedestals(int module, int channel) { 

  int j;

  file[0] = openTree("root/emcaltest_LEDBackground-23June17-off.root");
  file[1] = openTree("root/emcaltest_LEDBackground-23June17-1.root");
  file[2] = openTree("root/emcaltest_LEDBackground-23June17-2.root");
  file[3] = openTree("root/emcaltest_LEDBackground-23June17-3.root");
  file[4] = openTree("root/emcaltest_LEDBackground-23June17-4.root");
  file[5] = openTree("root/emcaltest_LEDBackground-23June17-5.root");
  file[6] = openTree("root/emcaltest_LEDBackground-23June17-6.root");
  file[7] = openTree("root/emcaltest_LEDBackground-23June17-8.root");
  file[8] = openTree("root/emcaltest_LEDBackground-23June17-9.root");
  int nfile = 9;

  TString LED[9] = {"LED Off", "LED-1","LED-2","LED-3","LED-4","LED-5",
		"LED-6","LED-7","LED-8"};
  
  TTree *adc;
  int nenteries = 0;
  
  TH1 *mean_hist[9], *rms_hist[9];
  TString hid, title;
  
  TGraph *pulse[9];
  int xpnt[31], ypnt[32];

  for (int ifile = 0; ifile < nfile; ifile++) {
    hid = "ped_";
    hid += ifile;
    title = "Channel Pedestal: ";
    title  += module; title += "-";title += channel;
    title += ": "; title += LED[ifile];
    mean_hist[ifile] = new TH1F(hid, title, 1000, 500., 2500.);

    hid = "rms_";
    hid += ifile;
    title = "Channel RMS: ";
    title  += module; title += "-";title += channel;
    title += ": "; title += LED[ifile];
    rms_hist[ifile] = new TH1F(hid, title, 100, 0., 100.);
  }


  for (int ifile = 0; ifile < nfile; ifile++) {

    std::cout << "Processing File " << ifile << std::endl; 
    file[ifile]->GetObject("ADC", adc);
    nenteries = adc->GetEntries();

    adc->SetBranchAddress("mean", mean);
    adc->SetBranchAddress("rms", rms);
    adc->SetBranchAddress("chanADC", chanADC);
    adc->SetBranchAddress("nsamples", &NSamples);
    
    for ( int i = 0; i < nenteries; i++) { 
      if (i%1000 == 0) std::cout << "processing event: " << i << std::endl;
      adc->GetEntry(i);
      if (rms[module*NCHAN+channel] > 0.)
	mean_hist[ifile]->Fill(mean[module*NCHAN + channel], 1.);
      if (rms[module*NCHAN+channel] > 0.)
	rms_hist[ifile]->Fill(rms[module*NCHAN + channel], 1.);	

      if (i == 100) { 
	for (j = 0; j < NSamples; j++) { 
	  xpnt[j] = j;
	  ypnt[j] = chanADC[module*NCHAN + channel][j];
	}
	pulse[ifile] = new TGraph(j, xpnt, ypnt);
      }
    } 
  }

  TCanvas *C1, *C2, *C3, *C4;
  TString cid;

  cid = "ped_1"; 
  title = "Channel Pedestal/RMS ";
  C1 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  C1->Divide(2,4);

  for (int i = 0; i < 4; i++) {
    C1->cd(i*2+1);
    mean_hist[i]->Draw();
    C1->cd(i*2+2);
    rms_hist[i]->Draw();
  }
    
  cid = "ped_2"; 
  title = "Channel Pedestal/RMS ";
  C2 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  C2->Divide(2,4);
  for (int i = 0; i < 4; i++) {
    C2->cd(i*2+1);
    mean_hist[i+4]->Draw();
    C2->cd(i*2+2);
    rms_hist[i+4]->Draw();
  }

  cid = "pulse-1"; 
  title = "Channel 10 Pulse ";
  C3 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  C3->Divide(2,4);
  for (int i = 0; i < 8; i++) { 
    C3->cd(i+1);
    pulse[i]->SetMinimum(0.0);
    pulse[i]->SetMaximum(3500.0);
    pulse[i]->GetXaxis()->SetTitle("time bin");
    pulse[i]->GetYaxis()->SetTitle("ADC Counts");
    pulse[i]->SetTitle(LED[i]);
    pulse[i]->Draw("AC*");
  }

  cid = "pulse-2"; 
  title = "Channel 10 Pulse ";
  C3 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  C3->Divide(2,4);
  for (int i = 0; i < 8; i++) { 
    C3->cd(i+1);
    if ( i < 4) {
      pulse[i]->SetMinimum(1000.0);
      pulse[i]->SetMaximum(2000.0);
    } else {
      pulse[i]->SetMinimum(0.0);
      pulse[i]->SetMaximum(3500.0);
    }
    pulse[i]->GetXaxis()->SetTitle("time bin");
    pulse[i]->GetYaxis()->SetTitle("ADC Counts");
    pulse[i]->SetTitle(LED[i]);
    pulse[i]->Draw("AC*");
  }
  
  TGraph *grms, *gmean;
  int i;
  for (i = 0; i < 8; i++) { 
    xpnt[i] = i;
    ypnt[i] = mean_hist[i]->GetMean(1);
  }
  gmean = new TGraph(i, xpnt, ypnt);
  for (i = 0; i < 8; i++) { 
    xpnt[i] = i;
    ypnt[i] = rms_hist[i]->GetMean(1);
  }
  grms = new TGraph(i, xpnt, ypnt);

  cid = "C4"; 
  title = "LED vs Mean Ped/RMS";
  C4 = new TCanvas(cid, title, 0., 0., 600., 300.);
  C4->Divide(2,1); 
  C4->cd(1);
  gmean->SetMinimum(1400.0);
  gmean->SetMaximum(1600.0);
  gmean->GetXaxis()->SetTitle("LED Setting");
  gmean->GetYaxis()->SetTitle("Mean Pedestal");
  gmean->SetTitle("LED vs Mean Pedestal");
  gmean->Draw("A*");
  C4->cd(2);
  grms->SetMinimum(0.0);
  grms->SetMaximum(50.0);
  grms->GetXaxis()->SetTitle("LED Setting");
  grms->GetYaxis()->SetTitle("Mean Pedestal RMS");
  grms->SetTitle("LED vs Mean Pedestal RMS");
  grms->Draw("A*");
 

}
  
void test (TString rootFileName1, TString rootFileName2) { 


  TH2D *mean, *rms;

  file1 = openTree(rootFileName1);
  // file1->ls();
  TTree *adc1, *adc2;
  int nentries1, nentries2;
  
  file1->GetObject("ADC", adc1);
 
  nentries1 =  adc1->GetEntries();
  std::cout << nentries1 << std::endl;
  

  file2 = openTree(rootFileName2);
  // file2->ls();
  file2->GetObject("ADC", adc2);
 
  nentries2 =  adc2->GetEntries();
  std::cout << nentries2 << std::endl;


  TString hid;
  TH1 *hist1, *hist2;

  mean = new TH2D("mean", "Mean Pedestal keithly vs Wiener", 1000, 1000., 2000.,
	       1000, 1000., 2000.);

  rms = new TH2D("rms", "Mean RMS keithly vs Wiener", 100, 0., 10.,
	       100, 0., 10.);

  float mean1, mean2, rmd1, rms2;
  
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 64; j++) {

      hid = "mean_";
      hid += i;
      hid += "-";
      hid += j;

      hist1 = (TH1 *)file1->Get(hid);
      hist2 = (TH1 *)file2->Get(hid);

      mean1 = hist1->GetMean(1);
      mean2 = hist2->GetMean(1);

      mean->Fill(mean1, mean2, 1);

      hid = "rms_";
      hid += i;
      hid += "-";
      hid += j;

      hist1 = (TH1 *)file1->Get(hid);
      hist2 = (TH1 *)file2->Get(hid);

      mean1 = hist1->GetMean(1);
      mean2 = hist2->GetMean(1);

      rms->Fill(mean1, mean2, 1);
    }
  }
  
  
  TCanvas *C1 = new TCanvas("c1","Mean Pedestal" );
  mean->Draw();
  TCanvas *C2 = new TCanvas("c2","Mean RMS" );
  rms->Draw();

  
}
