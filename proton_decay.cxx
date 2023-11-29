#include<TMath.h>
#include<TRandom3.h>
#include<TH1D.h>
#include<TH2D.h>
#include<TH3D.h>
#include<TTree.h>
#include<TCanvas.h>
#include<TFile.h>

#include<iostream>

void proton_decay(const int nevents=1000) {
     
    printf("nevents: %d\n", nevents);
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

    TList *lout = new TList;
    // create histograms
    TH1D *htheta = new TH1D("htheta", "Distribution of Theta;theta;count", 100, 0, 2 * pi);
    lout->Add(htheta);
    TH1D *hphi = new TH1D("hphi", "Distribution of Phi;phi;count", 100, 0, 2 * pi);
    lout->Add(hphi);
    TH2D *h2 = new TH2D("h2", "Theta vs Phi;theta;phi", 100, 0, 2 * pi, 100, 0, 2 * pi);
    lout->Add(h2);

    TH3D *hxyz = new TH3D("hxyz", "Distribution in 3D;x;y;z", 100, -pp, pp, 100, -pp, pp, 100, -pp, pp); lout->Add(hxyz);

    // Generate theta and phi for 1000 events
    // theta, phi have uniform spherical distribution

    double theta = -999;
    double phi = -999;
    double momentum = -999;

    double xx = -999;
    double yy = -999;
    double zz = -999;

    TTree *eventTree = new TTree("eventTree", "Events");
    eventTree->Branch("Theta", &theta);
    eventTree->Branch("Phi", &phi);
    eventTree->Branch("Momentum", &momentum);
    
    eventTree->Branch("x", &xx);
    eventTree->Branch("y", &yy);
    eventTree->Branch("z", &zz);

    for (int ii = 0; ii < nevents; ii++) {

	 theta = acos(rnd->Uniform(-1, 1));
	 phi = rnd->Uniform(0, 2 * pi);
	 momentum = pp;

	 // convert spherical coordinates to cartesian
	 xx = pp * sin(theta) * cos(phi);
	 yy = pp * sin(theta) * sin(phi);
	 zz = pp * cos(theta);
	 
	 htheta->Fill(theta);
	 hphi->Fill(phi);
	 h2->Fill(theta, phi);

	 hxyz->Fill(xx, yy, zz);

	 eventTree->Fill();
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    TCanvas *c2 = new TCanvas("c2", "c2");
    TCanvas *c3 = new TCanvas("c3", "c3");
    TCanvas *c4 = new TCanvas("c4", "c4");
    c1->cd();
    htheta->Draw();
    c2->cd();
    hphi->Draw();
    c3->cd();
    h2->Draw();
    c4->cd();
    hxyz->Draw();

    c1->Print("theta.png");
    c2->Print("phi.png");
    c3->Print("thetaphi.png");
    c4->Print("xyz.png");

    // Store pi0 events in tree
    TFile *file = TFile::Open("events.root", "RECREATE");
    lout->Write();
    

    eventTree->Write();
    file->Save();
    file->Close();
}
			
			
