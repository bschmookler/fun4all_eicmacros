//Code to compute EmCal E/p distribution for CORE detector configuration. 
//Single electrons and pions are generated separately from vertex=(0,0,0) 
//in the backwards eta region with several fixed momenta.
//The momentum is reconstructed using the track and saved in the standard tracking 
//evaluator, which the calorimeter energy is saved to the standard EEMC evaluator.

#include "Riostream.h"
#include "TApplication.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TChain.h"
#include "TVector3.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TSystem.h"
#include "TROOT.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

//---------------------
int main(int argc, char **argv){
	
#ifdef WITHRINT
	TRint *myapp = new TRint("RootSession",&argc,argv,NULL,0);
#else
	TApplication *myapp = new TApplication("myapp",0,0);
#endif

	//Set Style (seem to have to do it before histograms if compiling)
	gStyle->SetPadBorderMode(0);
    	gStyle->SetFrameBorderMode(0);
   	gStyle->SetFrameLineWidth(2);
    	gStyle->SetOptStat(0);
    	gStyle->SetOptFit(0);
    	gStyle->SetLabelSize(0.035,"X");
    	gStyle->SetLabelSize(0.035,"Y");
    	//gStyle->SetLabelOffset(0.01,"X");
    	//gStyle->SetLabelOffset(0.01,"Y");
    	gStyle->SetTitleXSize(0.055);
    	gStyle->SetTitleXOffset(0.85);
    	gStyle->SetTitleYSize(0.055);
    	gStyle->SetTitleYOffset(0.85);

	//Define momentum values
	const int num_mom = 5;
	double mom_val[num_mom] = {0.5,1.0,5.0,10.0,20.0};

	//Define Histograms, graphs, and function
	TH1 *he[num_mom];
	TH1 *hpi[num_mom];

	TH1 *he_eta[num_mom];
	TH1 *hpi_eta[num_mom];

	TH1 *he_lim[num_mom]; //generated eta from [-2.5,2]
        TH1 *hpi_lim[num_mom]; //generated eta from [-2.5,2]
	TH1 *he_limcut[num_mom]; //generated eta from [-2.5,2], E/p>0.8
        TH1 *hpi_limcut[num_mom]; //generated eta from [-2.5,2], E/p>0.8

	for(int ihist=0;ihist<num_mom;ihist++){
		he[ihist] = new TH1D(Form("he[%d]",ihist),"",100,0,1.25);
		he[ihist]->SetLineColor(kBlue);he[ihist]->SetLineWidth(3);
		he[ihist]->GetXaxis()->CenterTitle();he[ihist]->GetXaxis()->SetTitle("E/p");
		he[ihist]->GetYaxis()->CenterTitle();he[ihist]->GetYaxis()->SetTitle("Yield (arb.)");

		hpi[ihist] = new TH1D(Form("hpi[%d]",ihist),"",100,0,1.25);
		hpi[ihist]->SetLineColor(kRed);hpi[ihist]->SetLineWidth(3);
		hpi[ihist]->GetXaxis()->CenterTitle();hpi[ihist]->GetXaxis()->SetTitle("E/p");
                hpi[ihist]->GetYaxis()->CenterTitle();hpi[ihist]->GetYaxis()->SetTitle("Yield (arb.)");	
	
		he_eta[ihist] = new TH2D(Form("he_eta[%d]",ihist),Form("Momentum = %.1f GeV/c",mom_val[ihist]),100,-4,0,100,0,1.25);
		he_eta[ihist]->GetXaxis()->CenterTitle();he_eta[ihist]->GetXaxis()->SetTitle("True Electron #eta");
                he_eta[ihist]->GetYaxis()->CenterTitle();he_eta[ihist]->GetYaxis()->SetTitle("Electron E/p");
	
		hpi_eta[ihist] = new TH2D(Form("hpi_eta[%d]",ihist),Form("Momentum = %.1f GeV/c",mom_val[ihist]),100,-4,0,100,0,1.25);
                hpi_eta[ihist]->GetXaxis()->CenterTitle();hpi_eta[ihist]->GetXaxis()->SetTitle("True Neg. Pion #eta");
                hpi_eta[ihist]->GetYaxis()->CenterTitle();hpi_eta[ihist]->GetYaxis()->SetTitle("Neg. Pion E/p");

		he_lim[ihist] = new TH1D(Form("he_lim[%d]",ihist),"",100,0,1.25);
                he_lim[ihist]->SetLineColor(kBlue);he_lim[ihist]->SetLineWidth(3);
                he_lim[ihist]->GetXaxis()->CenterTitle();he_lim[ihist]->GetXaxis()->SetTitle("E/p");
                he_lim[ihist]->GetYaxis()->CenterTitle();he_lim[ihist]->GetYaxis()->SetTitle("Yield (arb.)");

                hpi_lim[ihist] = new TH1D(Form("hpi_lim[%d]",ihist),"",100,0,1.25);
                hpi_lim[ihist]->SetLineColor(kRed);hpi_lim[ihist]->SetLineWidth(3);
                hpi_lim[ihist]->GetXaxis()->CenterTitle();hpi_lim[ihist]->GetXaxis()->SetTitle("E/p");
                hpi_lim[ihist]->GetYaxis()->CenterTitle();hpi_lim[ihist]->GetYaxis()->SetTitle("Yield (arb.)");

		he_limcut[ihist] = new TH1D(Form("he_limcut[%d]",ihist),"",100,0,1.25);
                he_limcut[ihist]->SetLineColor(kBlue);he_limcut[ihist]->SetLineWidth(3);
                he_limcut[ihist]->GetXaxis()->CenterTitle();he_limcut[ihist]->GetXaxis()->SetTitle("E/p");
                he_limcut[ihist]->GetYaxis()->CenterTitle();he_limcut[ihist]->GetYaxis()->SetTitle("Yield (arb.)");

                hpi_limcut[ihist] = new TH1D(Form("hpi_limcut[%d]",ihist),"",100,0,1.25);
                hpi_limcut[ihist]->SetLineColor(kRed);hpi_limcut[ihist]->SetLineWidth(3);
                hpi_limcut[ihist]->GetXaxis()->CenterTitle();hpi_limcut[ihist]->GetXaxis()->SetTitle("E/p");
                hpi_limcut[ihist]->GetYaxis()->CenterTitle();hpi_limcut[ihist]->GetYaxis()->SetTitle("Yield (arb.)");

	}

	//Define variables
	float rec_px(0),rec_py(0),rec_pz(0);
	float gen_px(0),gen_py(0),gen_pz(0);
	float rec_e(0);
	int charge(0);
	TVector3 rec_p,gen_p;

	//Loop over ROOT files for each momentum
	for(int iFile=0;iFile<num_mom;iFile++)
	for(int itype=0;itype<2;itype++){
	
		//Load ROOT files for single particles
		TChain *tracks = new TChain("tracks");
		TChain *emcal = new TChain("ntp_gshower");
		for(int iRun=0;iRun<20;iRun++){
			if(itype==0){ //electrons
				tracks->Add(Form("../output/electron_single_emcal/G4COREDetector_g4tracking_eval_%.1f_%d.root",
					mom_val[iFile],iRun));
				emcal->Add(Form("../output/electron_single_emcal/G4COREDetector_g4eemc_eval_%.1f_%d.root",
					mom_val[iFile],iRun));
			}
			else if(itype==1){ //pions
				tracks->Add(Form("../output/pion_single_emcal/G4COREDetector_g4tracking_eval_%.1f_%d.root",
                                        mom_val[iFile],iRun));
                                emcal->Add(Form("../output/pion_single_emcal/G4COREDetector_g4eemc_eval_%.1f_%d.root",
                                        mom_val[iFile],iRun));
			}
		}

		tracks->AddFriend(emcal,"ntp_gshower");
		
		tracks->SetBranchAddress("gpx",&gen_px);
		tracks->SetBranchAddress("gpy",&gen_py);
		tracks->SetBranchAddress("gpz",&gen_pz);
		tracks->SetBranchAddress("px",&rec_px);
		tracks->SetBranchAddress("py",&rec_py);
		tracks->SetBranchAddress("pz",&rec_pz);
		tracks->SetBranchAddress("charge",&charge);
		tracks->SetBranchAddress("e",&rec_e);

		int nevents = tracks->GetEntries();
		std::cout <<"Number of events = "<<nevents<<std::endl;

		//Loop over events
		for(int iEvent=0;iEvent<nevents;iEvent++){
			if(iEvent%10000==0) std::cout<<"Events Analysed = "<<iEvent<<"!"<<std::endl;
        		tracks->GetEntry(iEvent);

			if(charge==-1 && rec_e>0){ //Make sure track is reconstructed and hits the calorimeter
				gen_p.SetXYZ(gen_px,gen_py,gen_pz);
				rec_p.SetXYZ(rec_px,rec_py,rec_pz);
				if(itype==0){
					he[iFile]->Fill(rec_e/rec_p.Mag());
					he_eta[iFile]->Fill(gen_p.Eta(),rec_e/rec_p.Mag());
					if(gen_p.Eta()>-2.5 && gen_p.Eta()<-2.0){
						he_lim[iFile]->Fill(rec_e/rec_p.Mag());
						if(rec_e/rec_p.Mag()>0.8) he_limcut[iFile]->Fill(rec_e/rec_p.Mag());
					}
				}
				else if(itype==1){
					hpi[iFile]->Fill(rec_e/rec_p.Mag());
					hpi_eta[iFile]->Fill(gen_p.Eta(),rec_e/rec_p.Mag());
					if(gen_p.Eta()>-2.5 && gen_p.Eta()<-2.0){
                                                hpi_lim[iFile]->Fill(rec_e/rec_p.Mag());
						if(rec_e/rec_p.Mag()>0.8) hpi_limcut[iFile]->Fill(rec_e/rec_p.Mag());
					}
				}
			}
		}//End event loop
	}//End loop for this type+momentum

	//Scale down histograms
	for(int ihist=0;ihist<num_mom;ihist++){
		he[ihist]->Scale(1./he[ihist]->GetEntries());
		hpi[ihist]->Scale(1./hpi[ihist]->GetEntries());
		he_lim[ihist]->Scale(1./he_lim[ihist]->GetEntries());
                hpi_lim[ihist]->Scale(1./hpi_lim[ihist]->GetEntries());
	}

	//Draw Plots
	TCanvas *c1[num_mom];
	TPaveText *pt1[num_mom],*pt2[num_mom],*pt3[num_mom];	

	for(int iCan=0;iCan<num_mom;iCan++){
		c1[iCan] = new TCanvas(Form("c1[%d]",iCan));
		
		if(he[iCan]->GetMaximum()>hpi[iCan]->GetMaximum()){
			he[iCan]->Draw("hist");hpi[iCan]->Draw("hist same");
		}
		else{
			hpi[iCan]->Draw("hist");he[iCan]->Draw("hist same");
		}

		//Draw text
		pt1[iCan] = new TPaveText(0.2,0.8,0.6,0.9,"NDCNB");
		pt1[iCan]->SetFillStyle(4000);
       	 	pt1[iCan]->SetTextFont(63);pt1[iCan]->SetTextSize(20);
        	pt1[iCan]->SetTextColor(kBlack);
        	pt1[iCan]->AddText(Form("Momentum = %.1f GeV/c",mom_val[iCan]));
        	pt1[iCan]->Draw();

		pt2[iCan] = new TPaveText(0.2,0.75,0.6,0.85,"NDCNB");
                pt2[iCan]->SetFillStyle(4000);
                pt2[iCan]->SetTextFont(63);pt2[iCan]->SetTextSize(20);
                pt2[iCan]->SetTextColor(kBlue);
                pt2[iCan]->AddText("Electrons");
                pt2[iCan]->Draw();

		pt3[iCan] = new TPaveText(0.2,0.7,0.6,0.8,"NDCNB");
                pt3[iCan]->SetFillStyle(4000);
                pt3[iCan]->SetTextFont(63);pt3[iCan]->SetTextSize(20);
                pt3[iCan]->SetTextColor(kRed);
                pt3[iCan]->AddText("Negative Pions");
                pt3[iCan]->Draw();
	
		if(iCan==0) c1[iCan]->Print("plots/emcal_core.pdf[");
		c1[iCan]->Print("plots/emcal_core.pdf");
	}

	TCanvas *c2[num_mom];

	for(int iCan=0;iCan<num_mom;iCan++){
                c2[iCan] = new TCanvas(Form("c2[%d]",iCan));
		c2[iCan]->Divide(2,1);
		c2[iCan]->cd(1);he_eta[iCan]->Draw("colz");
		c2[iCan]->cd(2);hpi_eta[iCan]->Draw("colz");

		c2[iCan]->Print("plots/emcal_core.pdf");
	}

	TCanvas *c3[num_mom];
	TPaveText *pt1a[num_mom],*pt4[num_mom],*pt5[num_mom],*pt6[num_mom];

        for(int iCan=0;iCan<num_mom;iCan++){
                c3[iCan] = new TCanvas(Form("c3[%d]",iCan));
		gPad->SetLogy();
		
                if(he_lim[iCan]->GetMaximum()>hpi_lim[iCan]->GetMaximum()){
                        he_lim[iCan]->Draw("hist");hpi_lim[iCan]->Draw("hist same");
                }
                else{
                        hpi_lim[iCan]->Draw("hist");he_lim[iCan]->Draw("hist same");
                }
		
		pt1a[iCan] = new TPaveText(0.15,0.8,0.65,0.9,"NDCNB");
                pt1a[iCan]->SetFillStyle(4000);
                pt1a[iCan]->SetTextFont(63);pt1a[iCan]->SetTextSize(20);
                pt1a[iCan]->SetTextColor(kBlack);
                pt1a[iCan]->AddText(Form("Momentum = %.1f GeV/c, -2.5 < #eta < 2.0",mom_val[iCan]));
                pt1a[iCan]->Draw();
		pt2[iCan]->Draw();pt3[iCan]->Draw();

		pt4[iCan] = new TPaveText(0.75,0.8,0.85,0.9,"NDCNB");
                pt4[iCan]->SetFillStyle(4000);
                pt4[iCan]->SetTextFont(63);pt4[iCan]->SetTextSize(20);
                pt4[iCan]->SetTextColor(kBlack);
                pt4[iCan]->AddText("For E/p > 0.8:");
                pt4[iCan]->Draw();

		double e_eff = ((double)he_limcut[iCan]->GetEntries())/he_lim[iCan]->GetEntries();

		pt5[iCan] = new TPaveText(0.775,0.75,0.875,0.85,"NDCNB");
                pt5[iCan]->SetFillStyle(4000);
                pt5[iCan]->SetTextFont(63);pt5[iCan]->SetTextSize(20);
                pt5[iCan]->SetTextColor(kBlue);
                pt5[iCan]->AddText(Form("Eff._{e} = %.2f",e_eff));
                pt5[iCan]->Draw();

		double pi_sup = ((double)hpi_lim[iCan]->GetEntries())/hpi_limcut[iCan]->GetEntries();
		
		pt6[iCan] = new TPaveText(0.775,0.7,0.875,0.8,"NDCNB");
                pt6[iCan]->SetFillStyle(4000);
                pt6[iCan]->SetTextFont(63);pt6[iCan]->SetTextSize(20);
                pt6[iCan]->SetTextColor(kRed);
                pt6[iCan]->AddText(Form("Rej._{#pi} = %.0f",pi_sup));
                pt6[iCan]->Draw();
	
		c3[iCan]->Print("plots/emcal_core.pdf");
		if(iCan==(num_mom-1))c3[iCan]->Print("plots/emcal_core.pdf]");
	}
	
	myapp->Run();
	return 0;
}

