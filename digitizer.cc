#include <iostream>
#include <cmath>
#include <cstdlib>
#include <map>

#include "Event.h"
#include "Eventiterator.h"
#include "fileEventiterator.h"
#include <packet_A.h>
#include <packetHeaders.h>

#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TDatime.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TObjString.h"

#include "digitizerrawdata.h"

using namespace std;

// program to look at digitzer raw data
// E.J. Mannel, BNL
// 16-Sept-2013

// Prototypes

void bookHist(void);
void plotHist(void);
void saveHist(void);
int digitizer(char *);

// Global variables

TH1 *h_ADC[MaxChannels];
TH1 *h_Mean[MaxChannels];
TH1 *h_RMS[MaxChannels];

int digitizer(int pkt,  char *file ) {

  int buffer[32767];
  int nw;
  int iword;
  int n_events = 0;
  int n_recs = 0;
  int status;

  Packet *p;
  Event *evt;
  Eventiterator *it;

  TString filename;

  adcdata *adc;
  adc = new adcdata();
  
  // first make a nice filename by making use of sthe string class

  TString rootfilename = file;
  TString delimiter = "/";
  
  // strip off the directory info
  TObjArray *tempString = rootfilename.Tokenize(delimiter);
  TIter iString(tempString);
  TObjString *os = 0;
  while ((os = (TObjString *)iString())) { 
    rootfilename = os->GetString().Data();
  }
  

  rootfilename.ReplaceAll(".PRDFF",".root");
  rootfilename.ReplaceAll(".prdf",".root");
  rootfilename.ReplaceAll(".evt",".root");

  cout << "creating root file: " << rootfilename << endl;


  // now make the root file
  TFile f(rootfilename,"recreate");

  bookHist();

  TTree ADC("ADC", "Digitizer Raw Data");
  ADC.SetMaxTreeSize(50000000000LL);

  ADC.Branch("runNumber", &adc->runNumber, "runNumber/I");
  ADC.Branch("eventNum", &adc->eventNumber, "eventNumber/I");
  ADC.Branch("detId", &adc->detId, "detId/I");
  ADC.Branch("module", &adc->module, "module/I");
  ADC.Branch("flag", &adc->flag, "flag/I");

  TString branchString = "bco[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("bco", adc->bco, branchString);

  branchString = "FEMHeader[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMHeader", adc->FEMHeader, branchString);

  branchString = "FEMSlot[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMSot", adc->FEMSlot, branchString);

  branchString = "FEMBCO[";
  branchString += NMOD;
  branchString += "]/I";
  ADC.Branch("FEMBCO", adc->FEMBCO, branchString);

  ADC.Branch("nsamples", &adc->nsamples, "nsamples/I");
  
  ADC.Branch("parity", &adc->parity, "parity[2][2]/I");

  branchString = "chanADC[";
  branchString += MaxChannels;
  branchString += "][";
  branchString += MaxSamples;
  branchString += "]/I";
  ADC.Branch("chanADC", adc->chanADC, branchString);

  branchString = "mean[";
  branchString += MaxChannels;
  branchString += "]/F";
  ADC.Branch("mean", adc->mean, branchString);

  branchString = "rms[";
  branchString += MaxChannels;
  branchString += "]/F";
  ADC.Branch("rms", adc->rms, branchString);
  
  filename = file;

  std::cout << "file: " << filename << std::endl;
  
  it =  new fileEventiterator (filename, status);

  if (status) {
    std::cout  << "Couldn't open input file " << filename << std::endl;
    delete it;
    exit(1);
  }
  
  std::cout << "Processing file" << std::endl;

  while ( (evt = it->getNextEvent())) {

    n_events++;

    for (Int_t raw_pkt = pkt; raw_pkt < pkt+1; raw_pkt++) {
      
      p = evt->getPacket( raw_pkt );
      
      if ( p ) {
	
      	p->fillIntArray(buffer, 32767, &nw, "DATA");
	
	std::cout << nw << std::endl;

	// there are a lot of kludges in this since there
      	// are still issues wit the packet format 9-Dec-2016
	
      	// this extracts the number of samples based on packet format 
      	// as of 9-Dec-2016
      	adc->nsamples = (nw - 11)/64;

      	adc->parity[0][0] = 0;
      	adc->parity[0][1] = 0;
      	adc->parity[1][0] = 0;
      	adc->parity[1][1] = 0;

      	for (int i = 0; i < MaxChannels; i++) { 
      	  adc->mean[i] = 0.;
      	  adc->rms[i] = 0.;
      	}

      	iword = 0;

      	adc->eventNumber =  buffer[iword] & 0xFFFF;
      	// adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 1
      	adc->flag  = buffer[iword] & 0xFFFF;
      	// adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 2
      	adc->detId  = buffer[iword] & 0xFFFF;
      	// adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 3
      	adc->module  = buffer[iword] & 0xFFFF;
      	// adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;            
      	iword++;

	// iword == 4
      	adc->bco[0]  = buffer[iword] & 0xFFFF;
      	// adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 5
      	adc->FEMHeader[0] = buffer[iword] & 0xFFFF;
      	adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 6
      	adc->FEMSlot[0] = buffer[iword] & 0xFFFF;
      	adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 7
      	adc->FEMEvtNum[0] = buffer[iword] & 0xFFFF;
      	adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;

	// iword == 8
      	adc->FEMBCO[0] = buffer[iword] & 0xFFFF;
      	adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;
      	iword++;
	
      	int ichan = 0;
      	int isample = 0;

	// iword == 9...
      	while ( iword < nw-2 ) {
      	  adc->chanADC[ichan][isample] = buffer[iword] & 0x00003FFF;
      	  adc->parity[1][iword%2] ^= buffer[iword] & 0x0000FFFF;

      	  adc->chanADC[ichan+1][isample] = buffer[iword+1] & 0x00003FFF;
      	  adc->parity[1][(iword+1)%2] ^= buffer[iword+1] & 0x0000FFFF;

      	  h_ADC[ichan]->Fill(adc->chanADC[ichan][isample], 1.);
      	  h_ADC[ichan+1]->Fill(adc->chanADC[ichan+1][isample], 1.);

      	  adc->mean[ichan] += (float) adc->chanADC[ichan][isample];
      	  adc->mean[ichan+1] += (float) adc->chanADC[ichan+1][isample];

      	  adc->rms[ichan] += (float) adc->chanADC[ichan][isample] *
      	    (float) adc->chanADC[ichan][isample];

      	  adc->rms[ichan+1] += (float) adc->chanADC[ichan+1][isample] *
      	    (float) adc->chanADC[ichan+1][isample];

      	  isample++;
      	  if (isample == adc->nsamples){

      	    adc->mean[ichan] /= (float)isample;
      	    adc->mean[ichan+1] /= (float) isample;

      	    h_Mean[ichan]->Fill(adc->mean[ichan], 1.);
      	    h_Mean[ichan+1]->Fill(adc->mean[ichan+1], 1.);

      	    adc->rms[ichan] = sqrt(adc->rms[ichan]/adc->nsamples -
      				   adc->mean[ichan]*adc->mean[ichan]);
      	    adc->rms[ichan+1] = sqrt(adc->rms[ichan+1]/adc->nsamples -
      				   adc->mean[ichan+1]*adc->mean[ichan+1]);


      	    h_RMS[ichan]->Fill(adc->rms[ichan], 1.);
      	    h_RMS[ichan+1]->Fill(adc->rms[ichan+1], 1.);

      	    isample = 0;
      	    ichan += 2;
      	  }

      	  iword+=2;

      	}

      	adc->parity[0][1] = buffer[iword] & 0x0000FFFF;
      	iword++;
      	adc->parity[0][0] = buffer[iword] & 0x0000FFFF;
      	iword++;

      	if ((adc->parity[0][0] != adc->parity[1][0]) ||
      	    (adc->parity[0][1] != adc->parity[1][1])) { 
      	  std::cout << "!!! Parity Error !!! " << std::endl;
      	  std::cout << "Event: "  << adc->eventNumber << std::endl;
      	  std::cout << "Packet size: " << nw << std::endl;
      	  printf ("Even Parity Words: %4.4X %4.4X\n", adc->parity[0][0],
      		  adc->parity[1][0]);
      	  printf ("Odd Parity Words: %4.4X %4.4X\n", adc->parity[0][1],
      		  adc->parity[1][1]);
      	  std::cout << std::endl;
      	}

      	n_recs++;

      	ADC.Fill();
      } else { 
	std::cout<< "Packet: " << raw_pkt << " not found" << std::endl;
      }
      
      delete p;
    }
    delete evt;

  }
  
  if (n_events > 0) { 

    std::cout << "Building Tree" << std::endl;

    ADC.BuildIndex("runNumber", "eventNumber");
    ADC.Write();
  
    // saveHist();
    f.Close();

  }

  delete adc;
  delete it;
  
  cout << "Digitizer raw: Final n_recs: " << n_recs << endl;
  cout << "Digitizer raw: Total events: " << n_events << endl;

  return 0;  
}



void saveHist(void) { 
  
  for ( int i = 0; i < MaxChannels; i++) { 
    h_ADC[i]->Write();
    h_Mean[i]->Write();
    h_RMS[i]->Write();
  }

  return;
}


void bookHist(void) {

  TString name, title;

  for (int i = 0; i < MaxChannels; i++) { 
    name = "adc_";
    name += i;
    title = "ADC ";
    title += i;
    h_ADC[i] = new TH1D(name, title, 16384, 0., 1684.);

    name = "mean_";
    name += i;
    title = "mean ";
    title += i;
    h_Mean[i] = new TH1D(name, title, 16384, 0., 1684.);

    name = "rms_";
    name += i;
    title = "rms ";
    title += i;
    h_RMS[i] = new TH1D(name, title, 100, 0., 10.);

  }

  return;

}

