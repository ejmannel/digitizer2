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

int eventNumber;
int chanADC[NCHAN*NMOD][NSAMPLE];
float mean[NCHAN*NMOD];
float rms[NCHAN*NMOD];
int NSamples;

void plot_mean(TString rootFileName) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  TH1 *mean_hist[NMOD][NCHAN];
  TString hid, title;
  float min,max;
  min = 0.;
  max = 5.;
  int nbin = 100;
  
  for (int j = 0; j < NMOD; j++) { 
    for (int i = 0; i < NCHAN; i++) { 

      hid = "mean_";
      hid += j;  hid += "-"; hid += i;
      title = "Channel Mean ";
      title  += j; title += "-";title += i;
        
      mean_hist[j][i] = new TH1F(hid, title, 1000, 1000., 2000.);
    }
  }

  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("mean", mean);
  adc->SetBranchAddress("rms", rms);
  adc->SetBranchAddress("nsamples", &NSamples);
  
  cout << "Number of samples/event: " << NSamples << endl;
  
  int nenteries = adc->GetEntries();
  
  for ( int i = 0; i < nenteries; i++) { 
    adc->GetEntry(i);
    for (int mod = 0; mod < NMOD; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
	  
	mean_hist[mod][chan]->Fill((float)mean[mod*NCHAN + chan], 1.);
	
      }
    }
  }

  TCanvas *C2 = new TCanvas("c2"," Temp " );
  mean_hist[1][0]->Draw();

  TCanvas *C1[NMOD*4];
  TString cid;
  for (int i = 0; i < NMOD; i++) {
    for (int j = 0; j < 4; j++) { 
      
      cid = "mean_"; cid += i; cid += "-"; cid += j; 
      title = "Mean- Module: "; title += i;
      title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
      C1[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 600., 400.);
      C1[i*NMOD+j]->Divide(4,4);
      
      for ( int k = 0; k < 16; k++){
  	C1[i*NMOD+j]->cd((k%16)+1);
  	mean_hist[i][j*16+k]->Draw();
      }
    }
  }
}

void plot_rms(TString rootFileName) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);


  TH1 *rms_hist[NMOD][NCHAN];
  TString hid, title;
  
  float min,max;
  min = 0.;
  max = 5.;
  int nbin = 100;
    
  for (int j = 0; j < NMOD; j++) { 
    for (int i = 0; i < NCHAN; i++) { 

      hid = "rms_";
      hid += j;  hid += "-"; hid += i;
      title = "Channel RMS ";
      title  += j; title += "-";title += i;
        
      rms_hist[j][i] = new TH1F(hid, title, 1000, 0., 10.);
    }
  }

  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("mean", mean);
  adc->SetBranchAddress("rms", rms);
  adc->SetBranchAddress("nsamples", &NSamples);

  cout << "Number of samples/event: " << NSamples << endl;

  int nenteries = adc->GetEntries();

  for ( int i = 0; i < nenteries; i++) { 
    adc->GetEntry(i);
    for (int mod = 0; mod < NMOD; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
	  
	  rms_hist[mod][chan]->Fill((float)rms[mod*NCHAN + chan], 1.);
	  
      }
    }
  }

  TCanvas *RMS[NMOD*4];
  TString cid;
  for (int i = 0; i < NMOD; i++) {
    for (int j = 0; j < 4; j++) { 
      
      cid = "rms_"; cid += i; cid += "-"; cid += j; 
      title = "RMS- Module: "; title += i;
      title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
      RMS[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 600., 400.);
      RMS[i*NMOD+j]->Divide(4,4);
      
      for ( int k = 0; k < 16; k++){
  	RMS[i*NMOD+j]->cd((k%16)+1);
  	rms_hist[i][j*16+k]->Draw();
      }
    }
  }
}


void plot_adc(TString rootFileName) { 
  
  TFile *rootFile = new TFile(rootFileName);
  
  TTree *adc;
  
  rootFile->GetObject("ADC", adc);
  
  TH1 *adc_hist[NMOD][NCHAN];
  TString hid, title;
  for (int j = 0; j < NMOD; j++) { 
    for (int i = 0; i < NCHAN; i++) { 
      hid = "adc_";
      hid += j;  hid += "-"; hid += i;
      
      title = "Channel ADC Distribution ";
      title  += j; title += "-";title += i;
      adc_hist[j][i] = new TH1F(hid, title, 300, 1300., 1600.);
    }  
  }

  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("nsamples", &NSamples);

  int nenteries = adc->GetEntries();

  cout << nenteries << endl;;

  int i;
  for ( i = 0; i < nenteries; i++) {

    adc->GetEntry(i);
    
    for (int mod = 0; mod < NMOD; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
  	for (int sample = 0; sample < NSamples; sample++) { 
	  
  	  adc_hist[mod][chan]->Fill((float)chanADC[mod*NCHAN + chan][sample], 1.);
	  
  	}
      }
    }
  }
  
  TCanvas *C1[NMOD*4];
  TString cid;
  for (int i = 0; i < NMOD; i++) {
    for (int j = 0; j < 4; j++) { 
      
      cid = "adc_"; cid += i; cid += "-"; cid += j; 
      title = "ADC- Module: "; title += i;
      title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
      C1[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 600., 400.);
      C1[i*NMOD+j]->Divide(4,4);
      
      for ( int k = 0; k < 16; k++){
  	C1[i*NMOD+j]->cd((k%16)+1);
  	adc_hist[i][j*16+k]->Draw();
      }
    }
  }
  
}
  

void plot_area(TString rootFileName, int max_evt = 0) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  TH1 *area_hist[NMOD][NCHAN];
  TH1 *ped_hist[NMOD][NCHAN];
  TH1 *rms_hist[NMOD][NCHAN];
  TString hid, title;
  float ped, area;

  float min = 0;
  float max = 0;
  float scale = 16.;

  TGraph *pulse[64];
  int pulseCnt = 0;
  int xpnt[32];
  int ypnt[32];
  
  for (int j = 0; j < NMOD; j++) { 
    for (int i = 0; i < NCHAN; i++) { 
      hid = "adc_";
      hid += j;  hid += "-"; hid += i;
      
      title = "Channel ADC Distribution ";
      title  += j; title += "-";title += i;

      area_hist[j][i] = new TH1F(hid, title, 500, 0., 40000.);

      hid = "ped_";
      hid += j;  hid += "-"; hid += i;
      
      title = "Channel PED Distribution ";
      title  += j; title += "-";title += i;	
      ped_hist[j][i] = new TH1F(hid, title, 350, 1250., 1600);

      hid = "rms_";
      hid += j;  hid += "-"; hid += i;
      
      title = "Channel RMS Distribution ";
      title  += j; title += "-";title += i;

      rms_hist[j][i] = new TH1F(hid, title, 100, 0., 10);
    }
  }

  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("nsamples", &NSamples);

  cout << "Number of samples/event: " << NSamples << endl;

  int nenteries = adc->GetEntries();

  if ( max_evt <= 0) max_evt = nenteries;

  int i;
  for ( i = 0; i < max_evt; i++) { 
    if ((i%100) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc->GetEntry(i);
    int nsample;
    int first_sample = 0;
    int last_sample = 5;

    float rms;

    for (int mod = 0; mod < NMOD; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
	ped = area = rms = 0.;
	nsample = 0;
	for (int sample = first_sample; sample <= last_sample; sample++) { 
	  ped += (float)chanADC[mod*NCHAN + chan][sample];
	
	  rms += ((float)chanADC[mod*NCHAN + chan][sample] * 
		  (float)chanADC[mod*NCHAN + chan][sample]);
	  
	  nsample++;
	}

	ped /= (float)(nsample);
	ped_hist[mod][chan]->Fill(ped, 1.);

	rms = rms/nsample - ped*ped;
	if (rms > 0.0) 
	  rms = sqrt(rms);
	else
	  rms = 0;

	rms_hist[mod][chan]->Fill(rms, 1.);
	
	int sample = 0;

	for (sample = 0; sample < NSamples; sample++) { 
	  area += ((float)chanADC[mod*NCHAN + chan][sample] - ped);	  
	  
	  if (pulseCnt < 64 && mod == 0 && chan == 18) { 
	    xpnt[sample] = sample;
	    ypnt[sample] = chanADC[mod*NCHAN+chan][sample];
	  }

	}	

	if (pulseCnt < 64 && mod == 0 && chan == 18) { 
	  pulse[pulseCnt] = new TGraph(sample-1, xpnt, ypnt);
	  pulse[pulseCnt]->SetMinimum(0.0);
	  pulse[pulseCnt]->SetMaximum(7500.);
	  pulseCnt++;
	}

	area_hist[mod][chan]->Fill(area, 1.);
     }
    }
  }

  std::cout << "Events processed: " << i << std::endl;

  int verbose = 0;
  if (verbose) { 
    TCanvas *PULSE[NMOD*4];
    TString cid;
    for (int i = 0; i < 1; i++) {
      for (int j = 0; j < 2; j++) { 
      
	cid = "pulse_"; cid += i; cid += "-"; cid += j; 
	title = "Pulse Area- Module: "; title += i;
	title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
	PULSE[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 600., 400.);
	PULSE[i*NMOD+j]->Divide(4,4);
	
	for ( int k = 0; k < 16; k++){
	  PULSE[i*NMOD+j]->cd((k%16)+1);
	  area_hist[i][j*16+k]->Draw();
	}
      }
    }
    
    TCanvas *PED[NMOD*4];
    for (int i = 0; i < 1; i++) {
      for (int j = 0; j < 2; j++) { 
	
	cid = "ped_"; cid += i; cid += "-"; cid += j; 
	title = "Pedestal- Module: "; title += i;
	title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
	PED[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 600., 400.);
	PED[i*NMOD+j]->Divide(4,4);
	
	for ( int k = 0; k < 16; k++){
	  PED[i*NMOD+j]->cd((k%16)+1);
	  ped_hist[i][j*16+k]->Draw();
	}
      }
    }
    
    TCanvas *RMS[NMOD*4];
    for (int i = 0; i < 1; i++) {
      for (int j = 0; j < 2; j++) { 
	
	cid = "rms_"; cid += i; cid += "-"; cid += j; 
	title = "Pedestal RMS- Module: "; title += i;
	title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
	RMS[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 600., 400.);
	RMS[i*NMOD+j]->Divide(4,4);
	
	for ( int k = 0; k < 16; k++){
	  RMS[i*NMOD+j]->cd((k%16)+1);
	  rms_hist[i][j*16+k]->Draw();
	}
      }
    }
    
    TCanvas *T1 = new TCanvas("t1", title, 0., 0., 600., 400.);
    T1->Divide(4,4);
    for (int i = 0; i < 16; i++) {
      T1->cd(i+1);
      pulse[i]->Draw("AC*");
    }
    
    TCanvas *ejm = new TCanvas("fig1", "Sample Channel", 0., 0., 600., 200.);
    ejm->Divide(3,1);
    ejm->cd(1); 
    pulse[1]->Draw("AC*");
    ejm->cd(2); 
    ped_hist[0][18]->Draw();
    ejm->cd(3); 
    area_hist[0][18]->Draw();
  }

  gStyle->SetOptStat(0000000);
  gStyle->SetOptFit();
  ped_hist[0][16]->Fit("gaus");
  ped_hist[0][17]->Fit("gaus");
  TCanvas *ejm1 = new TCanvas("fig2", "Pedestal", 0., 0., 600., 200.);
  ejm1->Divide(2,1);
  ejm1->cd(1); 
  ped_hist[0][16]->Draw("");
  ejm1->cd(2); 
  ped_hist[0][17]->SetAxisRange(1460., 1520., "X");
  ped_hist[0][17]->Draw("");

}


void plot_pulse(TString rootFileName, int module=0, int channel=0, 
	   int nevent=1, int fevent=0) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  if (nevent > 64) { nevent = 64; };

  TGraph *pulse[64];
  int xpnt[32];
  int ypnt[32];

  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("mean", mean);
  adc->SetBranchAddress("rms", rms);
  adc->SetBranchAddress("nsamples", &NSamples);

  int nenteries = adc->GetEntries();

  std::cout << "Number of events: " << nenteries << std::endl;

  if (fevent + nevent > nenteries) { 
    nevent = nenteries - fevent;
  } 

  int k = 0;
  int j = 0;

  std::cout << "Processing event range: " << fevent << " " 
	    << fevent+nevent << std::endl;  
  for ( int i = fevent; i < fevent+nevent; i++) { 
    adc->GetEntry(i);

    std::cout << "Number of samples/event: " << NSamples << std::endl;

    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j;
      ypnt[j] = chanADC[module*NCHAN+channel][j];
    }

    pulse[i-fevent] = new TGraph(j, xpnt, ypnt);
    pulse[i-fevent]->SetMinimum(0.0);
    pulse[i-fevent]->SetMaximum(7500. + (channel%2) * 0000);

  }

  TString title = "Channel ";
  title += channel;

  TCanvas *T1 = new TCanvas("t1", title, 0., 0., 600., 400.);
  T1->Divide(4,4);
  if (nevent > 16) { nevent = 16; }
  for (int i = 0; i < nevent; i++) {
    T1->cd(i+1);
    pulse[i]->Draw("AC*");
  }

}

void plot_pulse2(TString rootFileName, int module=0, int fevent=0) { 


  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);


  TGraph *pulse[64];
  int xpnt[32];
  int ypnt[32];

  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("mean", mean);
  adc->SetBranchAddress("rms", rms);
  adc->SetBranchAddress("nsamples", &NSamples);

  int nenteries = adc->GetEntries();

  if (fevent  > nenteries) { 
    fevent = nenteries-1;
  } 

  adc->GetEntry(fevent);

  cout << "Entry " << fevent << endl;
  cout << "Number of samples/event: " << NSamples << endl;

  int i, j;

  for (i = 0; i < 64; i++) { 

    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j;
      ypnt[j] = chanADC[module*NCHAN+i][j];

    }

    pulse[i] = new TGraph(j, xpnt, ypnt);
    pulse[i]->SetMinimum(0000.);
    pulse[i]->SetMaximum(5000. + (i%2)*7500.); /* low gain/High gain */

  }

  TCanvas *T1 = new TCanvas("t1", "Channels 0-15", 0., 0., 600., 400.);
  T1->Divide(4,4);
  for (int i = 0; i < 16; i++) {
    T1->cd(i+1);
    pulse[i]->Draw("AC*");
  }

  TCanvas *T2 = new TCanvas("t2", "Channels 16-31", 0., 0., 600., 400.);
  T2->Divide(4,4);
  for (int i = 16; i < 32; i++) {
    T2->cd(i-15);
    pulse[i]->Draw("AC*");
  }

  TCanvas *T3 = new TCanvas("t3", "Channels 32-47", 0., 0., 600., 400.);
  T3->Divide(4,4);
  for (int i = 32; i < 48; i++) {
    T3->cd(i-31);
    pulse[i]->Draw("AC*");
  }

  TCanvas *T4 = new TCanvas("t4", "Channels 48-63", 0., 0., 600., 400.);
  T4->Divide(4,4);
  for (int i = 48; i < 64; i++) {
    T4->cd(i-47);
    pulse[i]->Draw("AC*");
  }

  TCanvas *T5 = new TCanvas("t5", "Channels 0-3", 0., 0., 600., 400.);
  pulse[0]->Draw("AC*");
  for (int i = 1; i < 4; i++) {
    pulse[i]->Draw("C*");
  }

  TCanvas *T6 = new TCanvas("t6", "Channels 4-63", 0., 0., 600., 400.);
  pulse[4]->Draw("AC*");
  for (int i = 5; i < 64; i++) {
    pulse[i]->Draw("C*");
  }

  
}

void plot_FEMheader(TString rootFileName) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  int modulo = 1000;

  rootFile->GetObject("ADC", adc);

  int xpnt0[10000];
  int ypnt0[10000];
  int xpnt1[10000];
  int ypnt1[10000];

  int eventNumber, detId, module, flag, bco;
  int FEMSlot[NMOD], FEMEvtNum[NMOD], FEMBCO[NMOD];
  
  adc->SetBranchAddress("eventNum", &eventNumber);
  adc->SetBranchAddress("detId", &detId);
  adc->SetBranchAddress("module", &module);
  adc->SetBranchAddress("flag", &flag);
  adc->SetBranchAddress("bco", &bco);
  adc->SetBranchAddress("nsamples", &NSamples);
  adc->SetBranchAddress("FEMSlot", FEMSlot);
  adc->SetBranchAddress("FEMEvtNum", FEMEvtNum);
  adc->SetBranchAddress("FEMBCO", FEMBCO);
  adc->SetBranchAddress("chanADC", chanADC);
  adc->SetBranchAddress("mean", mean);
  adc->SetBranchAddress("rms", rms);

  TH1 *h_bco1D[20];
  TH2 *h_bco2D[20];

  TH1D *h_bco4, *h_bco5, *h_bco6;
  TH1D *h_bco7, *h_bco8, *h_bco9, *h_bco10;

  h_bco2D[0] = new TH2D("bco1", "Event BCO vs Module 0 BCO", 1000, 0., 66000.,
		   1000, 0., 66000.);
  h_bco2D[1] = new TH2D("bco2", "Event BCO vs Module 1 BCO", 1000, 0., 66000.,
		   1000, 0., 66000);
  h_bco2D[2] = new TH2D("bco2", "Event BCO vs Module 2 BCO", 1000, 0., 66000.,
		   1000, 0., 66000);

  TString title = "Global BCO Modulo "; 
  title += modulo;
  h_bco1D[0] = new TH1D("bco_0", title, modulo, 0., (float)modulo);

  h_bco4 = new TH1D("bco4", "Delta BCO: Event-Module0", 401, -20.5, 200.5);
  h_bco5 = new TH1D("bco5", "Delta BCO: Event-Module1", 401, -200.5, 200.5);
  h_bco6 = new TH1D("bco6", "Delta BCO: Event-Module2", 401, -200.5, 200.5);

  h_bco7 = new TH1D("bco7", "Delta BCO Module0", 65536, 0., 65536.);
  h_bco8 = new TH1D("bco8", "Delta BCO Module1", 65536, 0., 65536.);
  h_bco9 = new TH1D("bco9", "Delta BCO Module2", 65536, 0., 65536.);
  h_bco10 = new TH1D("bco10", "Delta BCO Global", 65536, 0., 65536.);

  int nenteries = adc->GetEntries();

  std::cout << "Number of events: " << nenteries << std::endl;

  int i = 0;
  float delta = 0.;

  int lastFEMBCO[3];
  int lastBCO;
  int first = 1;
  
  for (i = 0; i < nenteries; i++) { 

    if (i+1 < 1000) { 
      if ((i+1)%100 == 0)
	std::cout << "Processing event " << i+1 << std::endl;
    } else {
      if ((i+1)%1000 == 0)
	std::cout << "Processing event " << i+1 << std::endl;
    }

    adc->GetEntry(i);
   
    if ( i < 10000) {
      xpnt0[i] = eventNumber;
      ypnt0[i] = FEMEvtNum[0];
    
      xpnt1[i] = eventNumber;
      ypnt1[i] = FEMEvtNum[1];
    }

    if (first) {
      lastBCO = bco;
      lastFEMBCO[0] = FEMBCO[0];
      lastFEMBCO[1] = FEMBCO[1];
      lastFEMBCO[2] = FEMBCO[2];
      first = 0;
    } else { 
      delta = FEMBCO[0]-lastFEMBCO[0];
      if (delta < 0) delta += 32767.;
      h_bco7->Fill(delta, 1.);
      delta = FEMBCO[1]-lastFEMBCO[1];
      if (delta < 0) delta += 32767.;
      h_bco8->Fill(delta, 1.);
      delta = FEMBCO[2]-lastFEMBCO[2];
      if (delta < 0) delta += 32767.;
      h_bco9->Fill(delta, 1.);
      
      delta = fabs(bco-lastBCO);
      if (delta < 0) delta += 32767.;
      h_bco10->Fill(delta, 1.);
    }

    //update last BCO for next event
    lastBCO = bco;
    lastFEMBCO[0] = FEMBCO[0];
    lastFEMBCO[1] = FEMBCO[1];
    lastFEMBCO[2] = FEMBCO[2];
    
    h_bco2D[0]->Fill(bco, FEMBCO[0], 1.);
    h_bco2D[1]->Fill(bco, FEMBCO[1], 1.);
    h_bco2D[2]->Fill(bco, FEMBCO[2], 1.);

    h_bco1D[0]->Fill(bco%modulo, 1.);

    delta = bco-FEMBCO[0];
    if (delta < -1000) delta = -200;
    if (delta > 1000) delta = 200;
    h_bco4->Fill(delta, 1.);
 
    delta = bco-FEMBCO[1];
    if (delta < -1000) delta = -200;
    if (delta > 1000) delta = 200;
    h_bco5->Fill(delta, 1.);

    delta = bco-FEMBCO[2];
    if (delta < -1000) delta = -200;
    if (delta > 1000) delta = 200;
    h_bco6->Fill(delta, 1.);    

  }


  int ii = 0;
  if (i >= 10000)
    ii = 10000;
  else
    ii = i;

  TGraph *gr[2];
  gr[0] = new TGraph(ii, xpnt0, ypnt0);
  gr[0]->SetMinimum(0.0);
  gr[0]->SetMaximum(nenteries*30);
  gr[0]->GetXaxis()->SetTitle("Event Run Number");
  gr[0]->GetXaxis()->SetTitleOffset(0.3);
  gr[0]->GetXaxis()->SetLabelSize(0.025);
  gr[0]->GetXaxis()->SetLabelOffset(0.05);
  gr[0]->GetYaxis()->SetTitle("FEM Run Number");
  gr[0]->GetYaxis()->SetTitleOffset(0.3);
  gr[0]->GetYaxis()->SetLabelSize(0.025);
  gr[0]->GetYaxis()->SetLabelOffset(0.05);
  gr[0]->SetLineColor(1);
  
  gr[1] = new TGraph(ii, xpnt1, ypnt1);
  gr[1]->SetMinimum(0.0);
  gr[1]->SetMaximum(nenteries*30);
  gr[1]->GetXaxis()->SetTitle("Event Run Number");
  gr[1]->GetXaxis()->SetTitleOffset(0.3);
  gr[1]->GetXaxis()->SetLabelSize(0.025);
  gr[1]->GetXaxis()->SetLabelOffset(0.05);
  gr[1]->GetYaxis()->SetTitle("FEM Run Number");
  gr[1]->GetYaxis()->SetTitleOffset(0.3);
  gr[1]->GetYaxis()->SetLabelSize(0.025);
  gr[1]->GetYaxis()->SetLabelOffset(0.05);
  gr[1]->SetLineColor(2);
    
  title = "Event Number";
  TCanvas *T1 = new TCanvas("t1", title, 0., 0., 600., 400.);
  gr[0]->Draw("A*");
  gr[1]->Draw("L");

  title = "BCO Comparison";
  TCanvas *T2 = new TCanvas("t2", title, 0., 0., 600., 400.);
  T2->Divide(2,2);
  for (int i = 0; i < 3; i++) { 
    T2->cd(i+1);
    h_bco2D[i]->Draw();
  }

  title = "Delta BCO Comparison";
  TCanvas *T3 = new TCanvas("t3", title, 0., 0., 600., 400.);
  T3->Divide(2,2);
  T3->cd(1);
  h_bco4->Draw();
  T3->cd(2);
  h_bco5->Draw();
  T3->cd(3);
  h_bco6->Draw();

  title = "Delta BCO ";
  TCanvas *T4 = new TCanvas("t4", title, 0., 0., 600., 400.);
  T4->Divide(2,2);
  T4->cd(1);
  h_bco7->Draw();
  T4->cd(2);
  h_bco8->Draw();
  T4->cd(3);
  h_bco9->Draw();
  T4->cd(4);
  h_bco10->Draw();

  title = "Delta Global BCO ";
  TCanvas *T5 = new TCanvas("t5", title, 0., 0., 600., 400.);
  h_bco10->SetAxisRange(780., 830., "X");
  h_bco10->Draw();

  title = "Delta Module 0 BCO ";
  TCanvas *T6 = new TCanvas("t6", title, 0., 0., 600., 400.);
  h_bco7->SetAxisRange(780., 830., "X");
  h_bco7->Draw();

  title = "Global vs Module 1 BCO ";
  TCanvas *T7 = new TCanvas("t7", title, 0., 0., 600., 400.);
  h_bco2D[0]->Draw();

  title = "Global BCO ";
  TCanvas *T8 = new TCanvas("t8", title, 0., 0., 600., 400.);  
  TH1D *h_bco1x = h_bco2D[0]->ProjectionX();
  h_bco1x->Draw();

  title = "Global BCO Module 10000 ";
  TCanvas *T9 = new TCanvas("t9", title, 0., 0., 600., 400.);
  h_bco1D[0]->Draw();

}

