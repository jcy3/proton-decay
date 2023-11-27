#include<TMath.h>
#include<TRandom3.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TFile.h>

#include<iostream>

void proton_decay() {
    const double pi = TMath::Pi();
     
    // Particle masses in MeV/c2
    const double m_p = 938.27;   // proton
    
    const double m_pi = 134.98;  // pion
    const double m_e = 0.511;  // positron

    // Calculate momentum
    const double m_p2 = m_p * m_p;
    const double m_pi2 = m_pi * m_pi;
    const double m_e2 = m_e * m_e;
    
    const double frac = (m_p2 + m_e2 - m_pi2) / (2 * m_p);
    const double pp = sqrt(frac * frac - m_e2);

    // std::cout << "Momentum: " << p << " MeV/c" << "\n";

    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed(624);
    std::cout << "Random Seed:" << rnd->GetSeed() << "\n";

    // create histograms
    TH1D *htheta = new TH1D("htheta", "Distribution of Theta;theta;count", 100, 0, 2 * pi);
    TH1D *hphi = new TH1D("hphi", "Distribution of Phi;phi;count", 100, 0, 2 * pi);
    TH2D *h2 = new TH2D("h2", "Theta vs Phi;theta;phi", 100, 0, 2 * pi, 100, 0, 2 * pi);

    // Generate theta and phi for 1000 events
    // theta, phi are uniformly distributed

    double theta = -999;
    double phi = -999;
    double momentum = -999;

    TTree *eventTree = new TTree("eventTree", "Events");
    eventTree->Branch("Theta", &theta);
    eventTree->Branch("Phi", &phi);
    eventTree->Branch("Momentum", &momentum);
    
    for (int ii = 0; ii < 1000; ii++) {

	 theta = rnd->Uniform(0, 2 * pi);
	 phi = rnd->Uniform(0, 2 * pi);
	 momentum = pp;
	 
	 htheta->Fill(theta);
	 hphi->Fill(phi);
	 h2->Fill(theta, phi);

	 eventTree->Fill();
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    TCanvas *c2 = new TCanvas("c2", "c2");
    TCanvas *c3 = new TCanvas("c3", "c3");
    c1->cd();
    htheta->Draw();
    c2->cd();
    hphi->Draw();
    c3->cd();
    h2->Draw();

    c1->Print("theta.png");
    c2->Print("phi.png");
    c3->Print("thetaphi.png");

    // Store pi0 events in tree
    TFile *file = TFile::Open("events.root", "RECREATE");
    
    

    eventTree->Write();
    file->Save();
    file->Close();
}
			
			
