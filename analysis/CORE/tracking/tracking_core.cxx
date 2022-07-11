//Code to compute tracking resolution and acceptance
//for CORE detector configuration. Single muons are
//generated from vertex=(0,0,0) uniformly in eta=[-4.25,4.25]
//and momentum=[0.1,40]GeV/c and phi=[-PI,PI]. 
//We use the standard tracking evaluator to analyse the simultion.

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

#include "PadMxN.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>

//---------------------
//Detector matrix momentum resolution [in %]
double matrix_res(Double_t *x, Double_t *par){

    double mom = x[0];

    double sigp_p = sqrt( pow(par[0]*mom,2) + pow(par[1],2) );
    return sigp_p;
}
//---------------------

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

	//Define Histograms, graphs, and functions
	TH2 *hh1 = new TH2D("hh1","",200,-5,5,200,-5,45); //Generated momentum vs. eta
	hh1->GetXaxis()->SetTitle("Generated #eta");hh1->GetXaxis()->CenterTitle();
	hh1->GetYaxis()->SetTitle("Generated Momentum [GeV/c]");hh1->GetYaxis()->CenterTitle();

	TH2 *hh2 = new TH2D("hh2","",200,-5,5,400,0.01,45); //Generated transverse momentum vs. eta
	hh2->GetXaxis()->SetTitle("Generated #eta");hh2->GetXaxis()->CenterTitle();
	hh2->GetYaxis()->SetTitle("Generated Trans. Mom. [GeV/c]");hh2->GetYaxis()->CenterTitle();

	TH2 *hh3 = new TH2D("hh3","",200,-5,5,400,0.01,45); //Reconstructed transverse momentum vs. eta
	hh3->GetXaxis()->SetTitle("Reconstructed #eta");hh3->GetXaxis()->CenterTitle();
	hh3->GetYaxis()->SetTitle("Reconstructed Trans. Mom. [GeV/c]");hh3->GetYaxis()->CenterTitle();

	TH1 *h1 = new TH1D("h1","",200,-5,45); //Generated Momentum
	h1->GetXaxis()->SetTitle("Momentum [GeV/c]");h1->GetXaxis()->CenterTitle();
    h1->SetLineColor(kBlue);h1->SetLineWidth(3);
	
	TH1 *h2 = new TH1D("h2","",200,-5,5); //Generated Eta
	h2->GetXaxis()->SetTitle("#eta");h2->GetXaxis()->CenterTitle();
	h2->SetLineColor(kBlue);h2->SetLineWidth(3);
	
	TH1 *h3 = new TH1D("h3","",200,2,2); //Generated Phi
	h3->GetXaxis()->SetTitle("#phi [Rad]");h3->GetXaxis()->CenterTitle();
    h3->SetLineColor(kBlue);h3->SetLineWidth(3);	

	TH1 *h4 = new TH1D("h4","",200,-5,45);//Reconstructed Momentum
	h4->GetXaxis()->SetTitle("Momentum [GeV/c]");h4->GetXaxis()->CenterTitle();
    h4->SetLineColor(kRed);h4->SetLineWidth(3);

	TH1 *h5 = new TH1D("h5","",200,-5,5); //Reconstructed Eta
	h5->GetXaxis()->SetTitle("#eta");h5->GetXaxis()->CenterTitle();
	h5->SetLineColor(kRed);h5->SetLineWidth(3);

	TH1 *h6 = new TH1D("h6","",200,2,2); //Reconstructed Phi
	h6->GetXaxis()->SetTitle("#phi [Rad]");h6->GetXaxis()->CenterTitle();
    h6->SetLineColor(kRed);h6->SetLineWidth(3);

	//Define eta bins
	const int neta = 5;
	double eta_low[] = {-3.5, -2.5, -1.0, 1.0, 2.5};
	double eta_hi[] =  {-2.5, -1.0,  1.0, 2.5, 3.5};

	double par_mom[] = {0.2, 0.04, 0.04, 0.04, 0.2};
	double par_con[] = {5  , 2   , 1   , 2   , 5  };

	//Define detector matrix momentum resolution functions for B=1.5T
	//Based on https://github.com/eic/eicsmeardetectors/blob/master/SmearMatrixDetector_0_2_B1_5T.cxx
	//Not the same as detector matrix requirements in EIC Yellow Report Figure 10.6
	TF1 *fp[neta];
	
	for(int ieta=0;ieta<neta;ieta++){
		fp[ieta] = new TF1(Form("fp[%d]",ieta),matrix_res,0,40,2);
		fp[ieta]->SetParameter(0,par_mom[ieta]);fp[ieta]->SetParameter(1,par_con[ieta]);
		fp[ieta]->SetLineColor(kRed);fp[ieta]->SetLineWidth(3);
	}

	//Define momentum bins
	const int nmom = 40; //Number of momentum bins
	double mom_low[nmom];
	double mom_hi[nmom];
	double mom_cen[nmom];

	int mom_bin_low[neta] = {-50, -30, -25, -30, -50};
	int mom_bin_hi[neta]  = { 50,  30,  25,  30,  50};

	TH1 *hp[neta][nmom]; //for each eta bin, 40 bins in momentum (0-1GeV/c, 1-2GeV/c, etc...)

	for(int ieta=0;ieta<neta;ieta++){
		for(int imom=0;imom<nmom;imom++){

			mom_low[imom] = (double)imom;
			mom_hi[imom] = (double) (imom+1);
			mom_cen[imom] = (mom_low[imom]+mom_hi[imom]) / 2.;

			hp[ieta][imom] = new TH1F(Form("hp[%d][%d]",ieta,imom),Form("%.1f < #eta < %.1f, %.0f < Mom. < %.0f",eta_low[ieta],eta_hi[ieta],
									  mom_low[imom],mom_hi[imom]),400,mom_bin_low[ieta],mom_bin_hi[ieta]);
			hp[ieta][imom]->SetLineColor(kBlue);hp[ieta][imom]->SetLineWidth(3);
			hp[ieta][imom]->GetXaxis()->CenterTitle();hp[ieta][imom]->GetXaxis()->SetTitle("100 x #sigma_{p}/p");
		}
	}

	//Define variables
	float gen_px(0),gen_py(0),gen_pz(0);	
	float rec_px(0),rec_py(0),rec_pz(0);
	int charge(0);
	TVector3 gen_p,rec_p;

	//Load ROOT files
	TChain *t = new TChain("tracks");
	for(int i=0;i<10;i++){
		t->Add(Form("/sphenix/user/baschmoo/analysis/CORE/output/muon_single/G4COREDetector_g4tracking_eval_%d.root",i));
	}

	t->SetBranchAddress("gpx",&gen_px);
	t->SetBranchAddress("gpy",&gen_py);
	t->SetBranchAddress("gpz",&gen_pz);
	t->SetBranchAddress("px",&rec_px);
	t->SetBranchAddress("py",&rec_py);
	t->SetBranchAddress("pz",&rec_pz);
	t->SetBranchAddress("charge",&charge);

	int nevents = t->GetEntries();
	std::cout <<"Number of events = "<<nevents<<std::endl;

	//Loop over events
	for(int iEvent=0;iEvent<nevents;iEvent++){
		if(iEvent%10000==0) std::cout<<"Events Analysed = "<<iEvent<<"!"<<std::endl;
        t->GetEntry(iEvent);

		//Fill generated quantities/histograms
		gen_p.SetXYZ(gen_px,gen_py,gen_pz);
		
		hh1->Fill(gen_p.Eta(),gen_p.Mag());
		hh2->Fill(gen_p.Eta(),gen_p.Perp());
		h1->Fill(gen_p.Mag());
		h2->Fill(gen_p.Eta());
		h3->Fill(gen_p.Phi());

		if(charge==-1){ //Make sure track is reconstructed
			rec_p.SetXYZ(rec_px,rec_py,rec_pz);
			
			hh3->Fill(gen_p.Eta(),gen_p.Perp());
			h4->Fill(rec_p.Mag());
			h5->Fill(rec_p.Eta());
			h6->Fill(rec_p.Phi());

			for(int ieta=0;ieta<neta;ieta++){
				for(int imom=0;imom<nmom;imom++){
					if(gen_p.Eta()>eta_low[ieta] && gen_p.Eta()<eta_hi[ieta] && gen_p.Mag()>mom_low[imom] && gen_p.Mag()<mom_hi[imom])
						hp[ieta][imom]->Fill( 100. * (gen_p.Mag()-rec_p.Mag())/(mom_cen[imom]) );
				}
			}

		}
	}

	//Draw Plots
	TCanvas *c1 = new TCanvas("c1");
	hh1->Draw("colz");

	TLatex *tex1 = new TLatex(-4,47,"Single Muons generated at r = (0,0,0)");
	tex1->SetTextSize(0.05);
	tex1->Draw();

	c1->Print("plots/tracking_core.pdf[");
    c1->Print("plots/tracking_core.pdf");

	TCanvas *c2 = new TCanvas("c2");
	h1->Draw();h4->Draw("same");
	
	TPad* pad2 = new TPad("pad2","",0,0,1,1);
	pad2->SetFillStyle(4000);pad2->SetBorderMode(0);pad2->SetBorderSize(0);
	pad2->Draw();pad2->cd();
	
	TPaveText *pt1 = new TPaveText(0.15,0.925,0.45,0.975,"NDCNB");
	pt1->SetFillStyle(4000);	
	pt1->SetTextFont(63);pt1->SetTextSize(30);
	pt1->SetTextColor(kBlue);
	pt1->AddText("Generated");
	pt1->Draw();
	
	TPaveText *pt2 = new TPaveText(0.6,0.925,0.9,0.975,"NDCNB");
	pt2->SetFillStyle(4000);
    pt2->SetTextFont(63);pt2->SetTextSize(30);
    pt2->SetTextColor(kRed);
    pt2->AddText("Reconstructed");
    pt2->Draw();
	
	c2->Print("plots/tracking_core.pdf");

	TCanvas *c3 = new TCanvas("c3");
	h2->Draw();h5->Draw("same");
	pad2->Draw();pad2->cd();pt1->Draw();pt2->Draw();

	c3->Print("plots/tracking_core.pdf");

	TCanvas *c4 = new TCanvas("c4");
	h3->Draw();h6->Draw("same");
	pad2->Draw();pad2->cd();pt1->Draw();pt2->Draw();
	
	c4->Print("plots/tracking_core.pdf");

	TCanvas *c5 = new TCanvas("c5");
	c5->SetLogy();
	hh2->Draw("colz");

	c5->Print("plots/tracking_core.pdf");

	TCanvas *c6 = new TCanvas("c6");
	c6->SetLogy();
	hh3->Draw("colz");

	c6->Print("plots/tracking_core.pdf");

	//Draw and fit resolution histograms
	//Extract fit sigma to tgrapherrors
	TGraphErrors *gr[neta];

	for(int ieta=0;ieta<neta;ieta++){
		gr[ieta] = new TGraphErrors();
		gr[ieta]->SetMarkerColor(kBlue);gr[ieta]->SetLineColor(kBlue);
		gr[ieta]->SetMarkerStyle(20);gr[ieta]->SetMarkerSize(2);
		gr[ieta]->SetLineWidth(1);
	}


	PadMxN *pad4x4;
	PadMxN *pad4x2;
	TPad *mypad = {0};
	int counter = 0;

	for(int ieta=0;ieta<neta;ieta++){

		//Draw momentum from 0-16 GeV/c
		pad4x4 = new PadMxN(Form("c7a_%d",ieta),500,500,150,100,100,150,4,4);
    	pad4x4->Draw();
		mypad = {0};
		counter = 0;

		for(int imom=0;imom<16;imom++){
			mypad = pad4x4->GetPad(imom+1);
			hp[ieta][counter]->Draw();
			
			if(hp[ieta][counter]->GetEntries() > 10){
				hp[ieta][counter]->Fit("gaus");
		
				gr[ieta]->SetPoint(counter ,mom_cen[counter], hp[ieta][counter]->GetFunction("gaus")->GetParameter(2) );
				gr[ieta]->SetPointError(counter, 0, hp[ieta][counter]->GetFunction("gaus")->GetParError(2) );
			}
			else{
				gr[ieta]->SetPoint(counter ,mom_cen[counter], 0 );
				gr[ieta]->SetPointError(counter, 0, 0 );
			}

			counter++;
		}

		gROOT->ProcessLine( Form("c7a_%d->Print(\"plots/tracking_core.pdf\");",ieta) );

		//Draw momentum from 16-32 GeV/c
		pad4x4 = new PadMxN(Form("c7b_%d",ieta),500,500,150,100,100,150,4,4);
    	pad4x4->Draw();

		for(int imom=0;imom<16;imom++){
			mypad = pad4x4->GetPad(imom+1);
			hp[ieta][counter]->Draw();

			if(hp[ieta][counter]->GetEntries() > 10){
				hp[ieta][counter]->Fit("gaus");
		
				gr[ieta]->SetPoint(counter ,mom_cen[counter], hp[ieta][counter]->GetFunction("gaus")->GetParameter(2) );
				gr[ieta]->SetPointError(counter, 0, hp[ieta][counter]->GetFunction("gaus")->GetParError(2) );
			}
			else{
				gr[ieta]->SetPoint(counter ,mom_cen[counter], 0 );
				gr[ieta]->SetPointError(counter, 0, 0 );
			}
			

			counter++;
		}

		gROOT->ProcessLine( Form("c7b_%d->Print(\"plots/tracking_core.pdf\");",ieta) );

		//Draw momentum from 32-40 GeV/c
		pad4x2 = new PadMxN(Form("c7c_%d",ieta),500,500,150,100,100,150,4,2);
    	pad4x2->Draw();

		for(int imom=0;imom<8;imom++){
			mypad = pad4x2->GetPad(imom+1);
			hp[ieta][counter]->Draw();

			if(hp[ieta][counter]->GetEntries() > 10){
				hp[ieta][counter]->Fit("gaus");
		
				gr[ieta]->SetPoint(counter ,mom_cen[counter], hp[ieta][counter]->GetFunction("gaus")->GetParameter(2) );
				gr[ieta]->SetPointError(counter, 0, hp[ieta][counter]->GetFunction("gaus")->GetParError(2) );
			}
			else{
				gr[ieta]->SetPoint(counter ,mom_cen[counter], 0 );
				gr[ieta]->SetPointError(counter, 0, 0 );
			}
			
			counter++;
		}

		gROOT->ProcessLine( Form("c7c_%d->Print(\"plots/tracking_core.pdf\");",ieta) );
	}

	//Momentum Resolution
	PadMxN *pad3x2 = new PadMxN("c8",500,500,150,100,100,150,3,2);
    pad3x2->Draw();
	TH2 *htemp = new TH2F("htemp","",10,0,40,10,0,15);
        
	for(int ieta=0;ieta<neta+1;ieta++){
		mypad = pad3x2->GetPad(ieta+1);

		htemp->Draw();
		htemp->GetXaxis()->SetTitle("Momentum [GeV/c]");htemp->GetYaxis()->SetTitle("100 x #sigma_{p}/p [%]");
		htemp->GetXaxis()->SetLabelFont(63);htemp->GetYaxis()->SetLabelFont(63);
      	htemp->GetXaxis()->SetLabelSize(25);htemp->GetYaxis()->SetLabelSize(25);
      	htemp->GetXaxis()->SetLabelOffset(0.01);htemp->GetYaxis()->SetLabelOffset(0.01);
      	htemp->GetXaxis()->CenterTitle(1);htemp->GetYaxis()->CenterTitle(1);
      	htemp->GetXaxis()->SetTitleSize(35);htemp->GetXaxis()->SetTitleOffset(2.5); 
      	htemp->GetYaxis()->SetTitleSize(35);htemp->GetYaxis()->SetTitleOffset(3.0);

		if(ieta<5){
			fp[ieta]->Draw("same");
			gr[ieta]->Draw("P same");
		}
	}

	// get the text pad
    pad3x2->GetPad(7);

	TPaveText* pave7_1 = new TPaveText(0.1,0.825,0.25,0.875,"NDCNB");
    pave7_1->AddText("-3.5 < #eta < -2.5");
	pave7_1->SetFillStyle(4000);
    pave7_1->SetTextFont(63);pave7_1->SetTextSize(45);
    pave7_1->Draw();

	TPaveText* pave7_2 = new TPaveText(0.375,0.825,0.575,0.875,"NDCNB");
    pave7_2->AddText("-2.5 < #eta < -1.0");
	pave7_2->SetFillStyle(4000);
    pave7_2->SetTextFont(63);pave7_2->SetTextSize(45);
    pave7_2->Draw();

	TPaveText* pave7_3 = new TPaveText(0.65,0.825,0.85,0.875,"NDCNB");
    pave7_3->AddText("-1.0 < #eta < 1.0");
	pave7_3->SetFillStyle(4000);
    pave7_3->SetTextFont(63);pave7_3->SetTextSize(45);
    pave7_3->Draw();

	TPaveText* pave7_4 = new TPaveText(0.1,0.375,0.25,0.5,"NDCNB");
    pave7_4->AddText("1.0 < #eta < 2.5");
	pave7_4->SetFillStyle(4000);
    pave7_4->SetTextFont(63);pave7_4->SetTextSize(45);
    pave7_4->Draw();

	TPaveText* pave7_5 = new TPaveText(0.375,0.375,0.575,0.5,"NDCNB");
    pave7_5->AddText("2.5 < #eta < 3.5");
	pave7_5->SetFillStyle(4000);
    pave7_5->SetTextFont(63);pave7_5->SetTextSize(45);
    pave7_5->Draw();

	TPaveText* pave7_6 = new TPaveText(0.695,0.375,0.895,0.5,"NDCNB");
    pave7_6->AddText("Detector Matrix Resolutions");
	pave7_6->SetFillStyle(4000);
    pave7_6->SetTextFont(63);pave7_6->SetTextSize(35);
	pave7_6->SetTextColor(kRed);
    pave7_6->Draw();

	TPaveText* pave7_7 = new TPaveText(0.625,0.375,0.825,0.425,"NDCNB");
    pave7_7->AddText("for B = 1.5 T");
	pave7_7->SetFillStyle(4000);
    pave7_7->SetTextFont(63);pave7_7->SetTextSize(35);
	pave7_7->SetTextColor(kRed);
    pave7_7->Draw();

	TPaveText* pave7_8 = new TPaveText(0.65,0.3,0.85,0.35,"NDCNB");
    pave7_8->AddText("CORE Resolutions");
	pave7_8->SetFillStyle(4000);
    pave7_8->SetTextFont(63);pave7_8->SetTextSize(35);
	pave7_8->SetTextColor(kBlue);
    pave7_8->Draw();

	TPaveText* pave7_9 = new TPaveText(0.665,0.25,0.865,0.3,"NDCNB");
    pave7_9->AddText("for Uniform 1.5 T field");
	pave7_9->SetFillStyle(4000);
    pave7_9->SetTextFont(63);pave7_9->SetTextSize(35);
	pave7_9->SetTextColor(kBlue);
    pave7_9->Draw();

	gROOT->ProcessLine( "c8->Print(\"plots/tracking_core.pdf\");" );
	gROOT->ProcessLine( "c8->Print(\"plots/tracking_core.pdf]\");" );

	myapp->Run();
	return 0;
}

