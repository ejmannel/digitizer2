#include "../inc/std_headers.h"
#include "../inc/root.h" 
#include "../inc/utility.h" 

using namespace std;

static const int NSAMPLE = 32;
static const int NCHAN = 64;
static const int NMOD = 4;

// Function prototypes
// void plot_mean(TString, int, int);
// void plot_rms(TString, int, int);
// void plot_adc(TString);
// void plot_area(TString, int);
// void plot_pulsne(TString, int, int, int, int);
// void plot_pulse2(TString, int, int);
// void plot_FEMheader(TString);

int isLED(int *, int *);
void setBranch(TTree *);

//global variables for the branch structure
int eventNumber, detId, module, flag, bco;
int FEMSlot[NMOD], FEMEvtNum[NMOD], FEMBCO[NMOD];
int bclk[NMOD];
int chanADC[NCHAN*NMOD][NSAMPLE];
float pmean[NCHAN*NMOD];
float prms[NCHAN*NMOD];
int NSamples;

// map ADC channel number to EMCal position
// int tower_map[64] = 
//   {0, 1, 16, 17, 2, 3, 18, 19, 4, 5, 20, 21, 6, 7, 22, 23,
//    32, 33, 48, 49, 34, 35, 50, 51, 36, 37, 52, 53, 38, 39,54, 55,
//    8, 9, 24, 25, 10, 11, 26, 27, 12, 13, 28, 29, 14, 15, 30, 31,
//    40, 41, 56, 57, 42, 43, 58, 59, 44, 45, 60, 61, 46, 47, 62, 63};

// int tower_map[64] = 
//   {14, 15, 12, 13, 10, 11, 8, 9, 6, 7, 4, 5, 2, 3, 0, 1,
//    30, 31, 28, 29, 26, 27, 24, 25, 22, 23, 20, 21, 18, 19, 16, 17,
//    46, 47, 44, 45, 42, 43, 40, 41, 38, 39, 36, 37, 34, 35, 32, 33,
//    62, 63, 60, 61, 58, 59, 56, 57, 54, 55, 52, 53, 50, 51, 48, 49};


int tower_map[64] = 
  {14, 12, 10, 8, 15, 13, 11, 9, 30, 28, 26, 24, 31, 29, 27, 25,
   6, 4, 2, 0, 7, 5, 3, 1, 22, 20, 18, 16, 23, 21, 19, 17,
   46, 44, 42, 40, 47, 45, 43, 41, 62, 60, 58, 56, 63, 61, 59, 57,
   38, 36, 34, 32, 39, 37, 35, 33, 54, 52, 50, 48, 55, 53, 51, 49};

int oHCal_map[32] = { 5, 5, 1, 1, 13,13, 9, 9, 
		      6, 6, 2, 2, 14, 14, 10, 10, 
		      7, 7, 3, 3, 15, 15, 11, 11, 
		      8, 8, 4, 4, 16, 16, 12, 12};

int iHCal_map[16] = { 64, 68, 72, 76, 65, 69, 73, 77, 
		      66, 70, 74, 78, 67, 71, 75, 79};

int isLED(int *start, int *end) {
  int LEDMod = 2;
  int LEDChan = 17;
  int maxADC = 0;
  
  *start = *end = 0;
  int ped = 0;

  for (int i = 0; i < NSamples; i++)  {
    if ((chanADC[LEDMod*NCHAN + LEDChan][i] & 0x3FFF) > maxADC) {
      maxADC = chanADC[LEDMod*NCHAN + LEDChan][i] & 0x3FFF;
    }
    if ( i == 0) 
      ped = (chanADC[LEDMod*NCHAN + LEDChan][i] & 0x3FFF);
    else if (*start == 0 && (chanADC[LEDMod*NCHAN + LEDChan][i] & 0x3FFF)
	     > ped+30)
      *start = i;

    if (*start > 0 && *end == 0 && 
	(chanADC[LEDMod*NCHAN + LEDChan][i] & 0x3FFF) < ped+30)
      *end = i -1;

  }
  
  if (maxADC > 6000) 
    return(1);
  else
    return(0);
}

void setBranch (TTree *adc) { 
 
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
  adc->SetBranchAddress("mean", pmean);
  adc->SetBranchAddress("rms", prms);

  return;
}

void plot_mean(TString rootFileName, int fmod=0, int nummod=3) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  if (fmod > NMOD-1) { 
    std::cout << "first module can not be greater then " << NMOD << std::endl;
    return;
  }
  if (nummod > NMOD-fmod) { 
    std::cout << "number of modules can not be greater then " << NMOD-fmod 
	      << std::endl;
    return;
  }

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
        
      mean_hist[j][i] = new TH1F(hid, title, 300, 1400., 1700.);
    }
  }

  setBranch(adc);

  int nenteries = adc->GetEntries();
  
  for ( int i = 0; i < nenteries; i++) { 
    adc->GetEntry(i);
    if ((i+1 < 1000 && (i+1)%100 == 0) ||
	(i+1)%1000 == 0) { 
      std::cout << "Processing event " << i+1 << std::endl;
    }

    for (int mod = fmod; mod < nummod; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
	  
	mean_hist[mod][chan]->Fill((float)pmean[mod*NCHAN + chan], 1.);
	
      }
    }
  }

  Double_t ejm_chan_num[200], ejm_chan_err[200]; 
  Double_t ejm_ped_mean[200], ejm_ped_rms[200];
  
  for (int mod = fmod; mod < nummod; mod++) { 
    for ( int chan = 0; chan < NCHAN; chan++) { 
      ejm_chan_num[mod*64+chan] = (mod*64+chan);
      ejm_chan_err[mod*64+chan] = 0.0;
      ejm_ped_mean[mod*64+chan] = mean_hist[mod][chan]->GetMean(1);
      ejm_ped_rms[mod*64+chan] = mean_hist[mod][chan]->GetRMS(1);
    }
  }

  int pntCnt;
  pntCnt = 200;
  TGraphErrors *gr_ped1 = 
    new TGraphErrors(pntCnt, ejm_chan_num, ejm_ped_mean,
		     ejm_chan_err, ejm_ped_rms);

  // gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  TCanvas *C4 = new TCanvas("c4", "Channel vs Mean Pedestal");

  gr_ped1->SetMarkerStyle(2);
  gr_ped1->SetMarkerSize(0.5);
  gr_ped1->SetMarkerColor(2);
  gr_ped1->SetMinimum(1200.);
  gr_ped1->SetMaximum(2000.);
  gr_ped1->GetXaxis()->SetLimits(0.0, 200.);
  gr_ped1->Fit("pol1");
  gr_ped1->Draw("AP");
  
  TCanvas *C2 = new TCanvas("c2"," Module 0" );
  C2->Divide(2, 2);
  C2->cd(1);
  mean_hist[0][8]->Draw();
  C2->cd(2);
  mean_hist[0][24]->Draw();
  C2->cd(3);
  mean_hist[0][32]->Draw();
  C2->cd(4);
  mean_hist[0][56]->Draw();

  TString plotFileName;
  //Strip off leading directories
  TString delimiter = "/";
  TObjArray *tempString = rootFileName.Tokenize(delimiter);
  TIter iString(tempString);
  TObjString *os = 0;
  while ((os = (TObjString *)iString())) {
    plotFileName = os->GetString().Data();
  }  
  plotFileName.ReplaceAll(".root", "_mean_module0.png");
  plotFileName = "plots/" + plotFileName;
  C2->Print(plotFileName);

  TCanvas *C3 = new TCanvas("c3"," Module 2" );
  C3->Divide(2, 2);
  C3->cd(1);
  mean_hist[2][8]->Draw();
  C3->cd(2);
  mean_hist[2][24]->Draw();
  C3->cd(3);
  mean_hist[2][32]->Draw();
  C3->cd(4);
  mean_hist[2][56]->Draw();

  plotFileName.ReplaceAll("module0", "module1");
  C3->Print(plotFileName);

  TCanvas *C1[NMOD*4];
  TString cid;
  for (int i = fmod; i < nummod; i++) {
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

void plot_rms(TString rootFileName, int fmod=0, int lmod=2) { 

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
        
      rms_hist[j][i] = new TH1F(hid, title, 100, 0., 10.);
    }
  }

  setBranch(adc);

  int nenteries = adc->GetEntries();

  for ( int i = 0; i < nenteries; i++) { 
    adc->GetEntry(i);
    for (int mod = 0; mod < NMOD; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
	  
	  rms_hist[mod][chan]->Fill(prms[mod*NCHAN + chan], 1.);
	  
      }
    }
  }

  TCanvas *R2 = new TCanvas("R2"," Module 0 " );
  R2->Divide(2, 2);
  R2->cd(1);
  rms_hist[0][8]->Draw();
  R2->cd(2);
  rms_hist[0][24]->Draw();
  R2->cd(3);
  rms_hist[0][32]->Draw();
  R2->cd(4);
  rms_hist[0][56]->Draw();

  TString plotFileName;
  //Strip off leading directories
  TString delimiter = "/";
  TObjArray *tempString = rootFileName.Tokenize(delimiter);
  TIter iString(tempString);
  TObjString *os = 0;
  while ((os = (TObjString *)iString())) {
    plotFileName = os->GetString().Data();
  }  
  plotFileName.ReplaceAll(".root", "_rms_module0.png");
  plotFileName = "plots/" + plotFileName;
  R2->Print(plotFileName);

  TCanvas *R3 = new TCanvas("R3"," Module 2" );
  R3->Divide(2, 2);
  R3->cd(1);
  rms_hist[2][8]->Draw();
  R3->cd(2);
  rms_hist[2][24]->Draw();
  R3->cd(3);
  rms_hist[2][32]->Draw();
  R3->cd(4);
  rms_hist[2][56]->Draw();

  plotFileName.ReplaceAll("module0", "_module2");
  R3->Print(plotFileName);

  TCanvas *RMS[NMOD*4];
  TString cid;
  for (int i = fmod; i < lmod+1; i++) {
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

  setBranch(adc);

  int nenteries = adc->GetEntries();

  cout << nenteries << endl;;

  int i;
  for ( i = 0; i < nenteries; i++) {

    adc->GetEntry(i);
    
    for (int mod = 0; mod < NMOD; mod++) { 
      for ( int chan = 0; chan < NCHAN; chan++) { 
  	for (int sample = 0; sample < NSamples; sample++) { 
	  
  	  adc_hist[mod][chan]->
	    Fill((float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF), 1.);
	  
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

  TString hid, title;
  float ped, area;

  float min = -100;
  float max = 2000;
  float scale = 16.;
  int pstart, pend;
  
  for (int j = 0; j < NMOD; j++) { 
    for (int i = 0; i < NCHAN; i++) { 
      hid = "adc_";
      hid += j;  hid += "-"; hid += i;
      
      title = "Channel ADC Distribution ";
      title  += j; title += "-";title += i;

      area_hist[j][i] = new TH1F(hid, title, 5000, 0., 500000.);

      hid = "ped_";
      hid += j;  hid += "-"; hid += i;
      
      title = "Channel Pedestal Distribution ";
      title  += j; title += "-";title += i;

      ped_hist[j][i] = new TH1F(hid, title, 500, 1000., 2000.);
    }
  }

  setBranch(adc);

  int nenteries = adc->GetEntries();

  if ( max_evt >= 0) max_evt = nenteries;

  for ( int i = 0; i < max_evt; i++) { 
    if ((i%1000) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc->GetEntry(i);

    // if (isLED(&pstart, &pend)) { 

      int nsample;
      int maxADC = 0; 
      int maxTime = 0;
      for (int mod = 0; mod < NMOD; mod++) { 
	for ( int chan = 0; chan < NCHAN; chan++) { 
	  // Use first 5 samples for pedestal
	  ped = area = 0.;
	  nsample = 0;
	  for (int sample = 0; sample < 5; sample++) { 
	    ped += (float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF);
	    nsample++;
	  }
	  
	  ped /= (float)(nsample);

	  maxADC = maxTime = 0;
	  for (int sample = 5; sample < NSamples; sample++) { 
	    if (maxADC < (chanADC[mod*NCHAN + chan][sample] & 0x3FFF)) { 
	      maxADC = chanADC[mod*NCHAN + chan][sample] & 0x3FFF;
	      maxTime = sample;
	    } 
	    area += (float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF) - ped;	  
	  }	
	  
	  area_hist[mod][chan]->Fill(area, 1.);
	  ped_hist[mod][chan]->Fill(ped, 1.);

	}
      }
    // }
  }

  TCanvas *EJM, *EJM2;
  TString cid;

  cid = "p_";
  title = "LEDs Pulse Area for  4 Towers";

  EJM = new TCanvas(cid, title, 0., 0., 1200., 800.);
  EJM->Divide(2,2);
  
  EJM->cd(1);
  area_hist[0][0]->Draw();
  EJM->cd(2);
  area_hist[0][19]->Draw();
  EJM->cd(3);
  area_hist[0][35]->Draw();
  EJM->cd(4);
  area_hist[0][58]->Draw();

  cid = "p1_";
  title = "Pedestals for 4 Towers";

  EJM2 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  EJM2->Divide(2,2);
  
  EJM2->cd(1);
  ped_hist[0][0]->Draw();
  EJM2->cd(2);
  ped_hist[0][19]->Draw();
  EJM2->cd(3);
  ped_hist[0][35]->Draw();
  EJM2->cd(4);
  ped_hist[0][58]->Draw();

  TCanvas *PULSE[NMOD*4];

  for (int i = 0; i < 1; i++) {
    for (int j = 0; j < 4; j++) { 
      
      cid = "pulse_"; cid += i; cid += "-"; cid += j; 
      if (j == 0) 
	title = "Top Left";
      else if (j == 1) 
	title = "Top Right";
      else if (j == 2) 
	title = "Lower Left";
      else 
	title = "Lower Right";

      PULSE[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 1200., 800.);
      PULSE[i*NMOD+j]->Divide(4,4);
      
      for ( int k = 0; k < 16; k++){
  	PULSE[i*NMOD+j]->cd((k%16)+1);
  	area_hist[i][tower_map[j*16+k]]->Draw();
      }
    }
  }

  TCanvas *PED[NMOD*4];
  for (int i = 0; i < 1; i++) {
    for (int j = 0; j < 4; j++) { 
      
      cid = "ped_"; cid += i; cid += "-"; cid += j; 
      title = "Ped- Module: "; title += i;
      title += " Channels: "; title += j*16; title += "-"; title += j*16+15; 
      PED[i*NMOD+j] = new TCanvas(cid, title, 0., 0., 1200., 800.);
      PED[i*NMOD+j]->Divide(4,4);
      
      for ( int k = 0; k < 16; k++){
  	PED[i*NMOD+j]->cd((k%16)+1);
  	ped_hist[i][tower_map[j*16+k]]->Draw();
      }
    }
  }

}

void check_duplicate(TString rootFileName, int mod, int max_evt = 0) { 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  int pstart, pend;
  int dup[NCHAN];

  memset(dup, 0, sizeof(int)*NCHAN);
  
  setBranch(adc);

  int nenteries = adc->GetEntries();

  std::cout << "Processing " << nenteries << " events" << std::endl;

  if ( max_evt <= 0) max_evt = nenteries;

  for ( int i = 0; i < max_evt; i++) { 
    if ((i%1000) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc->GetEntry(i);

    if (isLED(&pstart, &pend)) { 
      pstart -= 2;
      pend -=2;
      for ( int chan = 0; chan < NCHAN; chan++) { 
	int isDuplicate = 0;
	for (int sample = pstart; sample < pend; sample++) { 
	  if ((chanADC[mod*NCHAN + chan][sample] & 0x3FFF) < 16380 &&
	      ((chanADC[mod*NCHAN + chan][sample] & 0x3FFF) ==
	       (chanADC[mod*NCHAN + chan][sample+1]&0x3FFF))) {
	    dup[chan]++;
	    isDuplicate = 1;
	  }
	}

	if (isDuplicate) {
	  std::cout << "Event: " << i << " Channel: " << chan 
		    << " pulse: " << pstart << " : " << pend << std::endl;
	  for (int sample = pstart; sample < pend+1; sample++) { 
	    std::cout << sample << ": " 
		      << (chanADC[mod*NCHAN + chan][sample]&0x3FFF) 
		      << " ";
	  }
	  std::cout << std::endl;
	}
      }
    }
  }

  for (int i = 0; i < NCHAN; i++) { 
    if (i%8 == 0) 
      std::cout << std::endl << i << ": ";
    std::cout << dup[i] << " ";
  }
  std::cout << std::endl;

}

void check_duplicate2(TString rootFileName, int module, int max_evt = 0) { 

  TFile *rootFile = new TFile(rootFileName);

  TH2D *hist_dupCnt = new TH2D("hdupCnt", "Number of Duplicate Hits/Event/Channel", 
			       32, -0.5, 31.5, 32, -0.5, 31.5);

  std::vector<float> evtNum;
  std::vector<float> duplicateCnt[6];
  std::vector<float> sequence;
  
  int lastTime;
  int delta[6][32];
  memset(delta, 0, 6*32*sizeof(int));

  TTree *adc;
  
  //  int ramp_chan[6] = { 0, 2, 4, 6, 9, 11};
  // int ramp_chan[6] = { 16, 18, 20, 22, 25, 27};
  int ramp_chan[6] = { 0, 1, 2, 3, 4, 5};
  // int ramp_chan[6] = { 32, 34, 36, 38, 41, 43};
  // int ramp_chan[6] = { 48, 50, 52, 54, 57, 59};
  
  rootFile->GetObject("ADC", adc);
  
  setBranch(adc);
  
  int nenteries = adc->GetEntries();
  
  std::cout << "Processing " << nenteries << " events" << std::endl;
  
  if ( max_evt <= 0) max_evt = nenteries;
  int seq = 0; 
  int lastSequence = 0;
  
  for ( int i = 0; i < max_evt; i++) { 
    if ((i%1000) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc->GetEntry(i);
    
    // std::cout << "Event : " << eventNumber << std::endl; 
    
    evtNum.push_back(eventNumber);
    int dupCnt = 0;
    
    for ( int channel = 0; channel < NCHAN; channel++) { 
      dupCnt = 0;
      
      // std::cout << "Chan: " << channel << " ";
      for (int sample = 0; sample < NSamples-1; sample++) { 
	seq = (chanADC[module*NCHAN + channel][sample] & 0xC000) >> 14;
	// std:: cout << seq << " ";
	lastSequence = seq;
      }
      // std::cout << endl;

      for (int ii = 0; ii < 4; ii++) { 
	
	if (channel == ramp_chan[ii]) {
	  // std::cout << "Ramp Channel: " << ii << " : " << chan << " : " 
	  // 	    << ramp_chan[ii] << std::endl;
	  
	  dupCnt = 0;
	  lastTime = -1;
	  int adc_1, adc_2;
	  for (int sample = 0; sample < NSamples-1; sample++) { 
	    adc_1 = chanADC[module*NCHAN + channel][sample] & 0x3FFF;
	    adc_2 = chanADC[module*NCHAN + channel][sample+1] & 0x3FFF;
	    if (adc_1 == adc_2) {
	      dupCnt++;
	      if (lastTime == -1)
		lastTime = sample;
	      else {
		delta[ii][sample-lastTime]++;
		lastTime = sample;
	      }
	    }
	  }
	}
	
	hist_dupCnt->Fill(channel, dupCnt);
      }	
    }
  }


  
  TCanvas *CDuplicate = new TCanvas("Cdup","Duplicate Hit Plots", 
  				    0., 0., 400., 300.);
  CDuplicate->Divide(2,2);
  CDuplicate->cd(1);
  hist_dupCnt->Draw("p");
  CDuplicate->cd(2);
  
  for (int i = 0; i < 6; i++) {
    for (int j = 0; j < 32; j++) 
      std::cout << delta[i][j] << "  ";
    std::cout << std::endl;
  }
}

void check_sequence(TString rootFileName, int module, int max_evt = 0) { 

  TFile *rootFile = new TFile(rootFileName);

  std::vector<float> evtNum;
  std::vector<float> duplicateCnt[6];
  std::vector<float> sequence;
  
  int lastTime;
  int delta[6][32];
  memset(delta, 0, 6*32*sizeof(int));

  TTree *adc;
  
  rootFile->GetObject("ADC", adc);
  
  setBranch(adc);
  
  int nenteries = adc->GetEntries();
  
  std::cout << "Processing " << nenteries << " events" << std::endl;
  
  if ( max_evt <= 0) max_evt = nenteries;
  int seq = 0; 
  int lastSequence = 0;
  
  for ( int i = 0; i < max_evt; i++) {
    if ((i%1000) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc->GetEntry(i);

   
    //std::cout << "Event : " << eventNumber << std::endl; 
    
    for ( int channel = 0; channel < NCHAN; channel++) { 

      std::cout << "Chan: " << channel << " ";
      for (int sample = 0; sample < NSamples-1; sample++) { 
	seq = (chanADC[module*NCHAN + channel][sample] & 0xC000) >> 14;
	printf("%X ", chanADC[module*NCHAN + channel][sample]);
	// std:: cout << seq << " ";
	lastSequence = seq;
      }
       std::cout << endl;
    }
  }
}

void plot_pulse(TString rootFileName, int module=0, int channel=0, 
	   int nevent=1, int fevent=0) { 


  int evtNumLast = 0;
  int evtCount = 0;
  
  TFile *rootFile = 0; 
  rootFile = new TFile(rootFileName);
  if (!rootFile) { 
    std::cout << "File " << rootFile << " not found, exiting" << std::endl;
    return;
  }

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  if (nevent > 128) { nevent = 128; };

  TGraph *pulse[128];
  TGraph *pulse2[128];
  int xpnt[32];
  int ypnt[32];

  setBranch(adc);

  int nenteries = adc->GetEntries();

  std::cout << "Number of events: " << nenteries << std::endl;

  if (fevent + nevent > nenteries) { 
    nevent = nenteries - fevent;
  } 

  int k = 0;
  int j = 0;
  
  int pstart, pend;

  std::cout << "Processing event range: " << fevent << " " 
	    << fevent+nevent << std::endl;  

  TString title = "Channel ";
  title += channel;

  for ( int i = fevent; i < fevent+nevent; i++) { 
    adc->GetEntry(i);

    pstart = pend = 0;

    if (eventNumber == evtNumLast) {
      evtNumLast = eventNumber;
      continue;
    }
    
    evtNumLast = eventNumber;

    int maxADC = 0;
    int minADC = 16000;
    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j;
      ypnt[j] = chanADC[module*NCHAN+channel][j] & 0x3FFF;
      if (ypnt[j] < minADC) minADC = ypnt[j];
      if (ypnt[j] > maxADC) maxADC = ypnt[j];
    }

    pulse[evtCount] = new TGraph(j, xpnt, ypnt);
    pulse[evtCount]->SetMinimum(1200.0);
    pulse[evtCount]->SetName(title);
    pulse[evtCount]->SetTitle(title);
    pulse[evtCount]->SetMarkerColor((i-fevent)%8+1);
    pulse[evtCount]->SetLineColor((i-fevent)%8+1);
		 
    pulse[evtCount]->GetXaxis()->SetTitle("ADC Time Bins");
    pulse[evtCount]->GetYaxis()->SetTitle("ADC Counts");

    pulse[evtCount]->GetXaxis()->SetTitleOffset(1.0);
    pulse[evtCount]->GetYaxis()->SetTitleOffset(1.3);

    pulse[evtCount]->GetXaxis()->SetLabelSize(0.03);
    pulse[evtCount]->GetYaxis()->SetLabelSize(0.03);
    
    pulse[evtCount]->SetMinimum(1200.0);
    pulse[evtCount]->SetMaximum(16000.0);

    pulse2[evtCount] = new TGraph(j, xpnt, ypnt);
    pulse2[evtCount]->SetMinimum(1200.0);
    pulse2[evtCount]->SetName(title);

    TString title2 = "Event ";
    title2 += eventNumber;
    pulse2[evtCount]->SetTitle(title);
    // pulse2[evtCount]->SetTextSize(0.5);
    pulse2[evtCount]->SetMarkerColor((i-fevent)%8+1);
    pulse2[evtCount]->SetLineColor((i-fevent)%8+1);

    pulse2[evtCount]->GetXaxis()->SetTitle("ADC Time Bins");
    pulse2[evtCount]->GetYaxis()->SetTitle("ADC Counts");

    pulse2[evtCount]->GetXaxis()->SetTitleOffset(1.0);
    pulse2[evtCount]->GetYaxis()->SetTitleOffset(1.3);

    pulse2[evtCount]->GetXaxis()->SetLabelSize(0.03);
    pulse2[evtCount]->GetYaxis()->SetLabelSize(0.03);
    
    pulse2[evtCount]->SetMinimum(minADC-100.);
    pulse2[evtCount]->SetMaximum(maxADC/0.8);
    // pulse2[evtCount]->SetMinimum(0.);
    // pulse2[evtCount]->SetMaximum(2500.);

    evtCount++;
    //    pulse[evtCount]->SetMaximum(1500. + (channel%2) * 7500);
				       
  }


  TCanvas *T1 = new TCanvas("t1", title, 0., 0., 1200., 800.);
  // int devent = evtCount;
  // if (devent < 5) {
  //   T1->Divide(2,2);
  // } else if (devent < 17) {
  //   devent = 17;
  //   T1->Divide(4,4);
  // } else if (devent < 37) {
  //   devent = 36;
  //   T1->Divide(6,6);
  // } else {
  //   devent = 64;
  //   T1->Divide(8,8);
  // }
    
  pulse2[0]->Draw("AL*");
  for (int i = 1; i < evtCount; i++) {
    pulse2[i]->Draw("L*");
  }

  TCanvas *T2 = new TCanvas("t2", title, 0., 0., 1200., 800.);
  pulse[0]->Draw("AL*");
  for (int i = 1; i < evtCount; i++) {
    pulse[i]->Draw("L*");
  }
  

  std::cout << "Event Count: " << nevent << " - " << evtCount << std::endl;
  
}

void plot_pulse2(TString rootFileName, int module=0, int fevent=0, 
		 int hilow = 0) { 


  TFile *rootFile; 
  rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  TGraph *pulse[64];
  int xpnt[32];
  int ypnt[32];

  setBranch(adc);

  int nenteries = adc->GetEntries();

  if (fevent  > nenteries) { 
    fevent = nenteries-1;
  } 

  adc->GetEntry(fevent);

  cout << "Entry " << fevent << endl;
  cout << "Number of samples/event: " << NSamples << endl;

  int i, j;
  TString title;

  int max_chan[64];
  memset (max_chan, 0, 64*sizeof(int));
  
  for (i = 0; i < 64; i++) { 

    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j;
      ypnt[j] = (chanADC[module*NCHAN+i][j] & 0x3FFF);

      if (max_chan[i] < ypnt[j])
	max_chan[i] = ypnt[j];

    }

    pulse[i] = new TGraph(j, xpnt, ypnt);

    pulse[i]->SetMinimum(0000.);

    if (hilow == 0)  
      pulse[i]->SetMaximum(5000.);
    else if (hilow == 1) 
      pulse[i]->SetMaximum(15000.);
    else if (hilow == 2)
      pulse[i]->SetMaximum(1800.); 
    else if (hilow == 3)  
      pulse[i]->SetMaximum(max_chan[i]*1.2);
    else
      pulse[i]->SetMaximum(20000.);

    title = "Chan "; title += i;
    pulse[i]->SetName(title);
    pulse[i]->SetTitle(title);
    pulse[i]->GetXaxis()->SetTitle("Time Bin (16 nSec bins)");
    pulse[i]->GetYaxis()->SetTitle("Amplitude (ADC counts)");
    pulse[i]->GetYaxis()->SetTitleOffset(1.3);
    pulse[i]->SetMarkerColor(i%8+1);
    pulse[i]->SetLineColor(i%8+1);
  }

  TCanvas *T1, *T2, *T3, *T4;
  if (module == 0) {
    T1 = new TCanvas("t1", "Top Left", 0., 0., 900., 600.);
    T1->Divide(4,4);
    for (int i = 0; i < 16; i++) {
      T1->cd(i+1);
      if (module == 0) 
	pulse[tower_map[i]]->Draw("AC*");
      else
	pulse[i]->Draw("AC*");
    }
  
    T2 = new TCanvas("t2", "Top Right", 0., 0., 900., 600.);
    T2->Divide(4,4);
    for (int i = 16; i < 32; i++) {
      T2->cd(i-15);
      if (module == 0) 
	pulse[tower_map[i]]->Draw("AC*");
      else
	pulse[i]->Draw("AC*");
    }
    
    T3 = new TCanvas("t3", "Bottom Left", 0., 0., 900., 600.);
    T3->Divide(4,4);
    for (int i = 32; i < 48; i++) {
      T3->cd(i-31);
      if (module == 0) 
	pulse[tower_map[i]]->Draw("AC*");
      else
	pulse[i]->Draw("AC*");
    }
    
    T4 = new TCanvas("t4", "Bottom Right", 0., 0., 900., 600.);
    T4->Divide(4,4);
    for (int i = 48; i < 64; i++) {
      T4->cd(i-47);
      if (module == 0) 
	pulse[tower_map[i]]->Draw("AC*");
      else
	pulse[i]->Draw("AC*");
    }

    TMultiGraph *mg_pulse[16];
    int flag = 0;

    Float_t min_scale = 1000.;
    Float_t max_scale = 20000.;

    if (hilow == 0)
      max_scale = 5000.;
    else if (hilow == 1) 
      max_scale = 20000.;
    else if (hilow == 2) 
      max_scale = 1800.;
    else if (hilow == 3) {
      max_scale = 0.;
      min_scale = 0.;
    }
      
    ULong_t tStart = 0;
    ULong_t tStop = 32;
    TString xTitle, yTitle;
    xTitle = "Clock Ticks";
    yTitle = "ADC Counts";
    
    TCanvas *pCanvas[16];
    
    for (int mg_id = 0; mg_id < 16; mg_id++) { 
      mg_pulse[mg_id] = new TMultiGraph();
      for ( int pulse_id = 4*mg_id; pulse_id < 4*mg_id + 4; pulse_id++) { 
	mg_pulse[mg_id]->Add(pulse[pulse_id]);
      }
      
      TString cId = "pC_";
      cId += mg_id;
      TString cTitle = "Channels ";
      cTitle += 4*mg_id;
      cTitle += "-";
      cTitle += 4*mg_id+3;
      
      mg_pulse[mg_id]->SetTitle(cTitle);
      pCanvas[mg_id] = new TCanvas(cId, cTitle, 0., 0., 600., 400.);
      mg_pulse[mg_id]->Draw("AC*");

      if (hilow == 3) {
	max_scale = 0;
	for (i = mg_id*4; i < mg_id*4+4; i ++){
	  if (max_scale < max_chan[i]*1.2)
	    max_scale = max_chan[i]*1.2;
	}
      }
 
      setupMGraph(mg_pulse[mg_id], flag, 
		  min_scale, max_scale, 
		  tStart, tStop, 
		  xTitle, yTitle);    
      
    }
    
    int map[16] = { 4, 3, 2, 1, 8, 7, 6, 5, 12, 11, 10, 9, 16, 15, 14, 13};
    
    TCanvas *allPulse = new TCanvas("allPulse", "Channels 0-63", 900., 600.);
    allPulse->Divide(4,4);
    for (int mg_id = 0;  mg_id < 16; mg_id++) { 
      allPulse->cd(map[mg_id]);
      mg_pulse[mg_id]->Draw("AC*");
    }
    
    TString plotFileName;
    //Strip off leading directories
    TString delimiter = "/";
    TObjArray *tempString = rootFileName.Tokenize(delimiter);
    TIter iString(tempString);
    TObjString *os = 0;
    while ((os = (TObjString *)iString())) {
      plotFileName = os->GetString().Data();
    }  
    plotFileName.ReplaceAll(".root", "_LED_1-16.png");
    plotFileName = "plots/" + plotFileName;
    allPulse->Print(plotFileName);
    
  } else if (module == 1) { 
    TCanvas *iHCal = new TCanvas("iHCal", "iHCal Towers", 900., 600.);
    iHCal->Divide(4,4);
    for (int i = 0; i < 16; i++) { 
      iHCal->cd(i+1);
      pulse[iHCal_map[i]-64]->Draw("AC");
    }

  } 

}

void plot_FEMheader(TString rootFileName) { 

  static const int Max_BCO = 10000; 

  TFile *rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  int xpnt0[Max_BCO];
  int ypnt0[Max_BCO];
  int xpnt1[Max_BCO];
  int ypnt1[Max_BCO];

  setBranch(adc);
  
  TH2D *h_bco1, *h_bco2, *h_bco3;
  TH1D *h_bco4, *h_bco5, *h_bco6;
  TH1D *h_bco7, *h_bco8, *h_bco9;
  TH1D *h_evtBCO1[4];
  TH1D *h_evtBCO2[4];


  h_bco1 = new TH2D("bco1", "XMIT BCO vs Module 2 BCO", 6555, 0., 65550.,
		   6555, 0., 65550.);
  h_bco2 = new TH2D("bco2", "XMIT BCO vs Module 2 BCO", 6555, 0., 65550,
		   6555, 0., 65550.);
  h_bco3 = new TH2D("bco3", "Module 1 BCO vs Module 2 BCO", 6555, 0., 65550,
		   6555, 0., 65550.);

  h_bco4 = new TH1D("bco4", "Delta BCO: XMIT-Module0", 401, -20.5, 200.5);
  h_bco5 = new TH1D("bco5", "Delta BCO: XMIT-Module1", 401, -200.5, 200.5);
  h_bco6 = new TH1D("bco6", "Delta BCO: Module0-Module1", 401, -200.5, 200.5);

  TString htitle, hid; 
  for (int i = 0; i < 4; i++) { 
    hid = "h_evtbco_"; 
    hid += i; 
    if ( i == 0) 
      htitle = "Event XMIT Delta BCO";
    else { 
      htitle = "Event Module ";
      htitle += (i-1); 
      htitle += " Delta BCO";
    } 
    h_evtBCO1[i] = new TH1D(hid, htitle, 1000, 0., 1000.);
    h_evtBCO2[i] = new TH1D(hid, htitle, 30000, 0., 30000.);
  } 

  int nenteries = adc->GetEntries();

  std::cout << "Number of events: " << nenteries << std::endl;

  int i = 0;
  int delta = 0;

  int lastBCO[4];
  int first = 1;
  
  for (i = 0; i < nenteries; i++) { 
    adc->GetEntry(i);
    
    if (i < Max_BCO) { 
      xpnt0[i] = eventNumber;
      ypnt0[i] = FEMEvtNum[0];

      xpnt1[i] = eventNumber;
      ypnt1[i] = FEMEvtNum[1];
    }
    
    h_bco1->Fill(bco, FEMBCO[0], 1.);
    h_bco2->Fill(bco, FEMBCO[1], 1.);
    h_bco3->Fill(FEMBCO[0], FEMBCO[1], 1.);

    delta = bco-FEMBCO[0];
    if (delta < -1000) delta = -200;
    if (delta > 1000) delta = 200;
    h_bco4->Fill(delta, 1.);
 
    delta = bco-FEMBCO[1];
    if (delta < -1000) delta = -200;
    if (delta > 1000) delta = 200;
    h_bco5->Fill(delta, 1.);

    delta = FEMBCO[0]-FEMBCO[1];
    if (delta < -1000) delta = -200;
    if (delta > 1000) delta = 200;
    h_bco6->Fill(delta, 1.);
   
    if (first) {
      first = 0;
    } else {
      h_evtBCO1[0]->Fill(bco-lastBCO[0], 1.);
      h_evtBCO2[0]->Fill(bco-lastBCO[0], 1.);
      for (int i = 0; i < 3; i++) {
	h_evtBCO1[i+1]->Fill(FEMBCO[i] - lastBCO[i+1], 1.);
	h_evtBCO2[i+1]->Fill(FEMBCO[i] - lastBCO[i+1], 1.);
      }
    }      
    lastBCO[0] = bco;
    for (int i = 0; i < 3; i++) {
      lastBCO[i+1] = FEMBCO[i];
    }
        
  }

  TGraph *gr[2];
  gr[0] = new TGraph(i, xpnt0, ypnt0);
  gr[0]->SetMinimum(0.0);
  gr[0]->SetMaximum(nenteries*30);
  gr[0]->GetXaxis()->SetTitle("FEM 0 Event Number");
  gr[0]->GetXaxis()->SetTitleOffset(0.3);
  gr[0]->GetXaxis()->SetLabelSize(0.025);
  gr[0]->GetXaxis()->SetLabelOffset(0.05);
  gr[0]->GetYaxis()->SetTitle("XMIT Run Number");
  gr[0]->GetYaxis()->SetTitleOffset(0.3);
  gr[0]->GetYaxis()->SetLabelSize(0.025);
  gr[0]->GetYaxis()->SetLabelOffset(0.05);
  gr[0]->SetLineColor(1);
  
  gr[1] = new TGraph(i, xpnt1, ypnt1);
  gr[1]->SetMinimum(0.0);
  gr[1]->SetMaximum(nenteries*30);
  gr[1]->GetXaxis()->SetTitle("FEM 1 Event Number");
  gr[1]->GetXaxis()->SetTitleOffset(0.3);
  gr[1]->GetXaxis()->SetLabelSize(0.025);
  gr[1]->GetXaxis()->SetLabelOffset(0.05);
  gr[1]->GetYaxis()->SetTitle("XMIT Run Number");
  gr[1]->GetYaxis()->SetTitleOffset(0.3);
  gr[1]->GetYaxis()->SetLabelSize(0.025);
  gr[1]->GetYaxis()->SetLabelOffset(0.05);
  gr[1]->SetLineColor(2);
    
  TString title = "Event Number";
  TCanvas *T1 = new TCanvas("t1", title, 0., 0., 600., 400.);
  gr[0]->Draw("A*");
  gr[1]->Draw("L");

  title = "BCO Comparison";
  TCanvas *T2 = new TCanvas("t2", title, 0., 0., 600., 400.);
  T2->Divide(2,2);
  T2->cd(1);
  h_bco1->Draw();
  T2->cd(2);
  h_bco2->Draw();
  T2->cd(3);
  h_bco3->Draw();

  title = "Delta BCO Comparison";
  TCanvas *T3 = new TCanvas("t3", title, 0., 0., 600., 400.);
  T3->Divide(2,2);
  T3->cd(1);
  h_bco4->Draw();
  T3->cd(2);
  h_bco5->Draw();
  T3->cd(3);
  h_bco6->Draw();

  title = "EventBCO Delta Comparison";
  TCanvas *T4 = new TCanvas("t4", title, 0., 0., 600., 400.);
  T4->Divide(2,2);
  for ( int i = 0; i < 4; i++) { 
    T4->cd(i+1);
    h_evtBCO1[i]->Draw();
  }  

  TCanvas *T5 = new TCanvas("t5", title, 0., 0., 600., 400.);
  T5->Divide(2,2);
  for ( int i = 0; i < 4; i++) { 
    T5->cd(i+1);
    h_evtBCO2[i]->Draw();
  }  
  
}

void plot_pulse3(TString rootFileName, int module=0, int channel=0, 
	   int nevent=1, int fevent=0) { 


  int evtNumLast = 0;
  int evtCount = 0;
  
  TFile *rootFile = 0; 
  rootFile = new TFile(rootFileName);
  if (!rootFile) { 
    std::cout << "File " << rootFile << " not found, exiting" << std::endl;
    return;
  }

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  if (nevent > 512) { nevent = 512; };

  TGraph *pulse[512];
  TGraph *pulse2[512];
  TGraph *gr1;

  int xpnt[64];
  int ypnt[64];

  setBranch(adc);

  int nenteries = adc->GetEntries();

  std::cout << "Number of events: " << nenteries << std::endl;

  if (fevent + nevent > nenteries) { 
    nevent = nenteries - fevent;
  } 

  int k = 0;
  int j = 0;

  for ( j = 0; j < 64; j++) { 
    xpnt[j] = j -32;
    ypnt[j] = 0;
  }
  gr1 = new TGraph(j, xpnt, ypnt);
  
  int pstart, pend;

  std::cout << "Processing event range: " << fevent << " " 
	    << fevent+nevent << std::endl;  

  TString title = "Channel ";
  title += channel;

  for ( int i = fevent; i < fevent+nevent; i++) { 
    adc->GetEntry(i);

    pstart = pend = 0;

    if (eventNumber == evtNumLast) {
      evtNumLast = eventNumber;
      continue;
    }
    
    evtNumLast = eventNumber;

    int ped = 0;
    for ( j = 0; j < 5; j++) { 
      ped += chanADC[module*NCHAN+channel][j] & 0x3FFF;
    }
    ped /= 5;

    int maxADC = 0;
    int minADC = 16000;
    int maxTime = 0;
    for ( j = 0; j < NSamples; j++) { 
      ypnt[j] = chanADC[module*NCHAN+channel][j] & 0x3FFF;
      if (ypnt[j] < (minADC - ped)) minADC = ypnt[j];
      if (ypnt[j] > (maxADC-ped)) { 
	maxADC = ypnt[j];
	maxTime = j;
      }  
    }
    
    cout << "maxTime: " << maxTime << endl;

    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j + (5 - maxTime);
      ypnt[j] = chanADC[module*NCHAN+channel][j] & 0x3FFF - ped;
    }

    pulse[evtCount] = new TGraph(j, xpnt, ypnt);
    pulse[evtCount]->SetMinimum(0.0);
    pulse[evtCount]->SetName(title);
    pulse[evtCount]->SetTitle(title);
    pulse[evtCount]->SetMarkerColor((i-fevent)%8+1);
    pulse[evtCount]->SetLineColor((i-fevent)%8+1);
		 
    pulse[evtCount]->GetXaxis()->SetTitle("ADC Time Bins");
    pulse[evtCount]->GetYaxis()->SetTitle("ADC Counts");

    pulse[evtCount]->GetXaxis()->SetTitleOffset(1.0);
    pulse[evtCount]->GetYaxis()->SetTitleOffset(1.3);

    pulse[evtCount]->GetXaxis()->SetLabelSize(0.03);
    pulse[evtCount]->GetYaxis()->SetLabelSize(0.03);
    
    pulse[evtCount]->SetMinimum(1200.0);
    pulse[evtCount]->SetMaximum(1700.0);

    pulse2[evtCount] = new TGraph(j, xpnt, ypnt);
    pulse2[evtCount]->SetMinimum(0.0);
    pulse2[evtCount]->SetName(title);
    TString title2 = "Event ";
    title2 += eventNumber;
    pulse2[evtCount]->SetTitle(title2);
    // pulse2[evtCount]->SetTextSize(0.5);

    pulse2[evtCount]->GetXaxis()->SetTitle("ADC Time Bins");
    pulse2[evtCount]->GetYaxis()->SetTitle("ADC Counts");

    pulse2[evtCount]->GetXaxis()->SetTitleOffset(1.0);
    pulse2[evtCount]->GetYaxis()->SetTitleOffset(1.3);

    pulse2[evtCount]->GetXaxis()->SetLabelSize(0.03);
    pulse2[evtCount]->GetYaxis()->SetLabelSize(0.03);
    
    // pulse2[evtCount]->SetMinimum(minADC-100.);
    // pulse2[evtCount]->SetMaximum(maxADC/0.8);
    pulse2[evtCount]->SetMinimum(0.);
    pulse2[evtCount]->SetMaximum(500.);

    evtCount++;
    //    pulse[evtCount]->SetMaximum(1500. + (channel%2) * 7500);
				       
  }

  TCanvas *T1 = new TCanvas("t1", title, 0., 0., 1200., 800.);
  int devent = evtCount;
  if (devent < 5) {
    T1->Divide(2,2);
  } else if (devent < 17) {
    devent = 17;
    T1->Divide(4,4);
  } else if (devent < 37) {
    devent = 36;
    T1->Divide(6,6);
  } else {
    devent = 64;
    T1->Divide(8,8);
  }
    
  for (int i = 0; i < devent; i++) {
    if (devent > 1) 
      T1->cd(i+1);
    pulse2[i]->Draw("A*");
  }


  TCanvas *T2 = new TCanvas("t2", title, 0., 0., 1200., 800.);
  gr1->SetMinimum(0.);
  gr1->SetMaximum(500.);
  gr1->Draw("AL*");
  for (int i = 0; i < evtCount; i++) {
    pulse[i]->Draw("L*");
  }
  

  std::cout << "Event COUnt: " << nevent << " - " << evtCount << std::endl;
  
}

void plot_area2(TString rootFileName1, TString rootFileName2, int max_evt = 0) { 

  TFile *rootFile1 = new TFile(rootFileName1);
  TFile *rootFile2 = new TFile(rootFileName2);

  TTree *adc1, *adc2;
 
  rootFile1->GetObject("ADC", adc1);
  rootFile2->GetObject("ADC", adc2);

  TH1 *area_hist1[NMOD][NCHAN];
  TH1 *ped_hist1[NMOD][NCHAN];
  TH1 *area_hist2[NMOD][NCHAN];
  TH1 *ped_hist2[NMOD][NCHAN];

  TString hid, title;
  float ped, area;

  float min = -100;
  float max = 2000;
  float scale = 16.;
  int pstart, pend;
  
  for (int j = 0; j < NMOD; j++) { 
    for (int i = 0; i < NCHAN; i++) { 

      title = "Channel ADC Distribution ";
      title  += j; title += "-";title += i;

      hid = "adc1_";
      hid += j;  hid += "-"; hid += i;
      area_hist1[j][i] = new TH1F(hid, title, 1000, 0., 100000.);
      area_hist1[j][i]->SetFillColor(1);

      hid = "adc2_";
      hid += j;  hid += "-"; hid += i;
      area_hist2[j][i] = new TH1F(hid, title, 1000, 0., 100000.);
      area_hist2[j][i]->SetFillColor(2);
      
      title = "Channel Pedestal Distribution ";
      title  += j; title += "-";title += i;

      hid = "ped1_";
      hid += j;  hid += "-"; hid += i;
      ped_hist1[j][i] = new TH1F(hid, title, 400, 1300., 1700.);
      ped_hist1[j][i]->SetFillColor(1);

      hid = "ped2_";
      hid += j;  hid += "-"; hid += i;
      ped_hist2[j][i] = new TH1F(hid, title, 400, 1300., 1700.);
      ped_hist2[j][i]->SetFillColor(2);
    }
  }

  setBranch(adc1);

  int nenteries = adc1->GetEntries();

  if ( max_evt <= 0) max_evt = nenteries;

  for ( int i = 0; i < max_evt; i++) { 
    if ((i%1000) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc1->GetEntry(i);

    // if (isLED(&pstart, &pend)) { 

      int nsample;
      int maxADC = 0; 
      int maxTime = 0;
      for (int mod = 0; mod < NMOD; mod++) { 
	for ( int chan = 0; chan < NCHAN; chan++) { 
	  // Use first 5 samples for pedestal
	  ped = area = 0.;
	  nsample = 0;
	  for (int sample = 0; sample < 5; sample++) { 
	    ped += (float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF);
	    nsample++;
	  }
	  
	  ped /= (float)(nsample);
	  
	  maxADC = maxTime = 0;

	  for (int sample = 0; sample < NSamples; sample++) { 
	    if (maxADC < (chanADC[mod*NCHAN + chan][sample] & 0x3FFF)) { 
	      maxADC = chanADC[mod*NCHAN + chan][sample] & 0x3FFF;
	      maxTime = sample;
	    }
	    area += (float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF) - ped;	  
	  }	
	  if (maxTime >= 15) { 
	    area_hist1[mod][chan]->Fill(area, 1.);
	    ped_hist1[mod][chan]->Fill(ped, 1.);
	  }
	}
      }
    // }
  }

  setBranch(adc2);

  nenteries = adc2->GetEntries();

  if ( max_evt <= 0) max_evt = nenteries;

  for ( int i = 0; i < max_evt; i++) { 
    if ((i%1000) == 0) 
      std::cout << "Processing event " << i << std::endl;
    adc2->GetEntry(i);

    // if (isLED(&pstart, &pend)) { 

      int nsample;
      int maxADC = 0;
      int maxTime = 0;

      for (int mod = 0; mod < NMOD; mod++) { 
	for ( int chan = 0; chan < NCHAN; chan++) { 
	  // Use first 5 samples for pedestal
	  ped = area = 0.;
	  nsample = 0;
	  for (int sample = 0; sample < 5; sample++) { 
	    ped += (float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF);
	    nsample++;
	  }
	  
	  ped /= (float)(nsample);
	  
	  maxADC = maxTime = 0;
	  for (int sample = 0; sample < NSamples; sample++) { 
	    area += (float)(chanADC[mod*NCHAN + chan][sample] & 0x3FFF) - ped;	  
	    if (maxADC < (chanADC[mod*NCHAN + chan][sample] & 0x3FFF)) { 
	      maxADC = chanADC[mod*NCHAN + chan][sample] & 0x3FFF;
	      maxTime = sample;
	    }
	  }	
	  
	  // if (chan == 0) { 
	  //   cout << maxTime << "  " << maxADC << endl;
	  // }
	  if (maxTime >= 15) { 
	    area_hist2[mod][chan]->Fill(area, 1.);
	    ped_hist2[mod][chan]->Fill(ped, 1.);
	  }
	}
      }
    // }
  }

  TCanvas *EJM, *EJM2, *EJM3, *EJM4;
  TString cid;

  cid = "p_";
  title = "LEDs Pulse Area for  4 Towers";

  EJM = new TCanvas(cid, title, 0., 0., 1200., 800.);
  EJM->Divide(2,2);
  
  int dchan[4] = {13, 25, 37, 59};

  EJM->cd(1);
  area_hist1[0][dchan[0]]->Draw();
  area_hist2[0][dchan[0]]->Draw("same");
  EJM->cd(2);
  area_hist1[0][dchan[1]]->Draw();
  area_hist2[0][dchan[1]]->Draw("same");
  EJM->cd(3);
  area_hist1[0][dchan[2]]->Draw();
  area_hist2[0][dchan[2]]->Draw("same");
  EJM->cd(4);
  area_hist1[0][dchan[3]]->Draw();
  area_hist2[0][dchan[3]]->Draw("same");

  cid = "p1_";
  title = "Pedestals for 4 Towers";

  EJM2 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  EJM2->Divide(2,2);
  
  EJM2->cd(1);
  ped_hist1[0][dchan[0]]->Draw();
  ped_hist2[0][dchan[0]]->Draw("same");
  EJM2->cd(2);
  ped_hist1[0][dchan[1]]->Draw();
  ped_hist2[0][dchan[1]]->Draw("same");
  EJM2->cd(3);
  ped_hist1[0][dchan[2]]->Draw();
  ped_hist2[0][dchan[2]]->Draw("same");
  EJM2->cd(4);
  ped_hist1[0][dchan[3]]->Draw();
  ped_hist2[0][dchan[3]]->Draw("same");

  cid = "p3_";
  title = "LEDs Pulse Area for  4 Towers";

  EJM3 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  EJM3->Divide(2,2);
  
  EJM3->cd(1);
  area_hist2[0][dchan[0]]->Draw();
  EJM3->cd(2);
  area_hist2[0][dchan[1]]->Draw();
  EJM3->cd(3);
  area_hist2[0][dchan[2]]->Draw();
  EJM3->cd(4);
  area_hist2[0][dchan[3]]->Draw();

  cid = "p4_";
  title = "Pedestals for 4 Towers";

  EJM4 = new TCanvas(cid, title, 0., 0., 1200., 800.);
  EJM4->Divide(2,2);
  
  EJM4->cd(1);
  ped_hist2[0][dchan[0]]->Draw();
  EJM4->cd(2);
  ped_hist2[0][dchan[1]]->Draw();
  EJM4->cd(3);
  ped_hist2[0][dchan[2]]->Draw();
  EJM4->cd(4);
  ped_hist2[0][dchan[3]]->Draw();


}

void fit_pulse(TString rootFileName, int module=0, int channel=0, 
	   int nevent=1, int fevent=0) { 


  int evtNumLast = 0;
  int evtCount = 0;

  TFile *rootFile = 0; 
  rootFile = new TFile(rootFileName);
  if (!rootFile) { 
    std::cout << "File " << rootFile << " not found, exiting" << std::endl;
    return;
  }
  
  TTree *adc;
  
  rootFile->GetObject("ADC", adc);
  
  if (nevent > 128) { nevent = 128; };
  
  TGraph *pulse[128];
  TGraph *pulse2[128];
  int xpnt[32];
  int ypnt[32];
  
  setBranch(adc);
  
  int nenteries = adc->GetEntries();
  
  std::cout << "Number of events: " << nenteries << std::endl;
  
  if (fevent + nevent > nenteries) { 
    nevent = nenteries - fevent;
  } 
  
  int k = 0;
  int j = 0;
  
  int pstart, pend;
  
  std::cout << "Processing event range: " << fevent << " " 
	    << fevent+nevent << std::endl;  
  
  TString title = "Channel ";
  title += channel;
  
  for ( int i = fevent; i < fevent+nevent; i++) { 
    adc->GetEntry(i);
    
    pstart = pend = 0;
    
    if (eventNumber == evtNumLast) {
      evtNumLast = eventNumber;
      continue;
    }
    
    evtNumLast = eventNumber;
    
    int maxADC = 0;
    int minADC = 16000;
    int maxTime = 0;

    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j;
      ypnt[j] = chanADC[module*NCHAN+channel][j] & 0x3FFF;
      if (ypnt[j] < minADC) minADC = ypnt[j];
      if (ypnt[j] > maxADC) { 
	maxADC = ypnt[j];
	maxTime = j;
      }
    }
    
    pulse[evtCount] = new TGraph(j, xpnt, ypnt);
    pulse[evtCount]->SetMinimum(0.0);
    pulse[evtCount]->SetName(title);
    pulse[evtCount]->SetTitle(title);
    pulse[evtCount]->SetMarkerColor((i-fevent)%8+1);
    pulse[evtCount]->SetLineColor((i-fevent)%8+1);
    
    pulse[evtCount]->GetXaxis()->SetTitle("ADC Time Bins");
    pulse[evtCount]->GetYaxis()->SetTitle("ADC Counts");
    
    pulse[evtCount]->GetXaxis()->SetTitleOffset(1.0);
    pulse[evtCount]->GetYaxis()->SetTitleOffset(1.3);
    
    pulse[evtCount]->GetXaxis()->SetLabelSize(0.03);
    pulse[evtCount]->GetYaxis()->SetLabelSize(0.03);
    
    pulse[evtCount]->SetMinimum(0.0);
    pulse[evtCount]->SetMaximum(16000.0);
    pulse[evtCount]->Fit("gaus", "R", "", maxTime-3, maxTime+5);
    
    pulse2[evtCount] = new TGraph(j, xpnt, ypnt);
    pulse2[evtCount]->SetMinimum(0.0);
    pulse2[evtCount]->SetName(title);
    
    TString title2 = "Event ";
    title2 += eventNumber;
    pulse2[evtCount]->SetTitle(title);
    // pulse2[evtCount]->SetTextSize(0.5);
    pulse2[evtCount]->SetMarkerColor((i-fevent)%8+1);
    pulse2[evtCount]->SetLineColor((i-fevent)%8+1);
    
    pulse2[evtCount]->GetXaxis()->SetTitle("ADC Time Bins");
    pulse2[evtCount]->GetYaxis()->SetTitle("ADC Counts");
    
    pulse2[evtCount]->GetXaxis()->SetTitleOffset(1.0);
    pulse2[evtCount]->GetYaxis()->SetTitleOffset(1.3);
    
    pulse2[evtCount]->GetXaxis()->SetLabelSize(0.03);
    pulse2[evtCount]->GetYaxis()->SetLabelSize(0.03);
    
    pulse2[evtCount]->SetMinimum(minADC-100.);
    pulse2[evtCount]->SetMaximum(maxADC/0.8);
    pulse2[evtCount]->SetMinimum(0.);
    pulse2[evtCount]->SetMaximum(5000.);
    
    pulse2[evtCount]->Fit("gaus", "R", "", maxTime-3, maxTime+5);
    
    evtCount++;
    
  }
  TCanvas *T1 = new TCanvas("t1", title, 0., 0., 1200., 800.);
    
  pulse2[0]->Draw("AL*");
  for (int i = 1; i < evtCount; i++) {
    pulse2[i]->Draw("L*");
  }

  TCanvas *T2 = new TCanvas("t2", title, 0., 0., 1200., 800.);
  pulse[0]->Draw("AL*");
  for (int i = 1; i < evtCount; i++) {
    pulse[i]->Draw("L*");
  }
  

  std::cout << "Event Count: " << nevent << " - " << evtCount << std::endl;

}

void plot_oHCal(TString rootFileName,int nevent=1, int fevent=0) { 
  TFile *rootFile; 
  rootFile = new TFile(rootFileName);

  TTree *adc;
 
  rootFile->GetObject("ADC", adc);

  TGraph *pulse[32];
  int xpnt[32];
  int ypnt[32];

  setBranch(adc);

  int nenteries = adc->GetEntries();

  if (fevent  > nenteries) { 
    fevent = nenteries-1;
  } 

  adc->GetEntry(fevent);

  cout << "Entry " << fevent << endl;
  cout << "Number of samples/event: " << NSamples << endl;

  int i, j;
  TString title;
  int max;

  for (i = 0; i < 32; i++) { 
    max = 0;
    for ( j = 0; j < NSamples; j++) { 
      xpnt[j] = j;
      ypnt[j] = (chanADC[i+112][j] & 0x3FFF);
      if (ypnt[j] > max) max = ypnt[j];
    }

    pulse[i] = new TGraph(j, xpnt, ypnt);
    pulse[i]->SetMinimum(0000.);
    pulse[i]->SetMaximum(max/0.8); /* low gain/High gain */

    title = "Chan "; title += (i+112);
    pulse[i]->SetName(title);
    pulse[i]->SetTitle(title);
    pulse[i]->GetXaxis()->SetTitle("Time Bin (16 nSec bins)");
    pulse[i]->GetYaxis()->SetTitle("Amplitude (ADC counts)");
    pulse[i]->GetYaxis()->SetTitleOffset(1.3);
    pulse[i]->SetMarkerColor(i%8+1);
    pulse[i]->SetLineColor(i%8+1);
  }

  TCanvas *normal_Canvas = new TCanvas("norm", "oHCal Normal Gain", 800, 600);
  normal_Canvas->Divide(4,4);
  for ( int i = 0; i < 32; i+=2) { 
    normal_Canvas->cd(oHCal_map[i]);
    pulse[i]->Draw("AL");
  }

  TCanvas *higain_Canvas = new TCanvas("higain", "oHCal High Gain", 800, 600);
  higain_Canvas->Divide(4,4);
  for ( int i = 1; i < 32; i+=2) { 
    higain_Canvas->cd(oHCal_map[i]);
    pulse[i]->Draw("AL");
  }
  
}
