#include<TMath.h>
#include<TRandom3.h>

#include<TH1D.h>
#include<TH2D.h>
#include<TH3D.h>

#include<TTree.h>
#include<TCanvas.h>
#include<TFile.h>

#include "Math/Vector4D.h"

#include<iostream>

void proton_decay(const int nevents=1000) {

    printf("nevents: %d\n", nevents);
    const double pi = TMath::Pi();

    // Particle masses in MeV/c2
    const double m_p = 938.27;   // proton

    const double m_pi = 134.98;  // pion
    const double m_e = 0.511;  // positron

    const double p_p = 200;  // proton initial momentum

    const double EE2 = p_p * p_p + m_p * m_p;  // total energy squared

    const double gamma = sqrt(1 + (p_p/m_p) * (p_p/m_p));
    const double beta = sqrt(1 - (1 / (gamma * gamma)));

    printf("gamma: %lf, beta: %lf\n", gamma, beta);

    // Calculate momentum of decay products
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
    TH1D *htheta = new TH1D("htheta", "Distribution of Theta;theta;count", 100, 0, pi);
    lout->Add(htheta);
    TH1D *hphi = new TH1D("hphi", "Distribution of Phi;phi;count", 100, 0, 2 * pi);
    lout->Add(hphi);
    TH2D *h2 = new TH2D("h2", "Theta vs Phi;theta;phi", 100, 0, pi, 100, 0, 2 * pi);
    lout->Add(h2);

    TH3D *hxyz = new TH3D("hxyz", "Distribution in 3D;x;y;z", 100, 1, -0, 100, -pp, pp, 100, -pp, pp);
    lout->Add(hxyz);

    // Generate theta and phi for 1000 events
    // theta, phi have uniform spherical distribution in rest frame

    double theta = -999;
    double phi = -999;
    double momentum = -999;

    double xx = -999;
    double yy = -999;
    double zz = -999;

    double theta_p = -999;
    double phi_p = -999;
    double momentum_p = -999;

    double x_p = -999;
    double y_p = -999;
    double z_p = -999;

    double theta_pi = -999;
    double phi_pi = -999;
    double p_pi = -999;

    double theta_e = -999;
    double phi_e = -999;
    double p_e = -999;
 
    double x_pi = -999;
    double y_pi = -999;
    double z_pi = -999;

    double x_e = -999;
    double y_e = -999;
    double z_e = -999;

    TTree *eventTree = new TTree("eventTree", "Events");

    /*
    eventTree->Branch("Theta Pi", &theta_pi);
    eventTree->Branch("Phi Pi", &phi_pi);
    eventTree->Branch("Momentum Pi", &p_pi);

    eventTree->Branch("Theta Positron", &theta_e);
    eventTree->Branch("Phi Positron", &phi_e);
    eventTree->Branch("Momentum Positron", &p_e);
    */

    eventTree->Branch("x Proton", &x_p);
    eventTree->Branch("y Proton", &y_p);
    eventTree->Branch("z Proton", &z_p);

    eventTree->Branch("x Pi", &x_pi);
    eventTree->Branch("y Pi", &y_pi);
    eventTree->Branch("z Pi", &z_pi);

    eventTree->Branch("x Positron", &x_e);
    eventTree->Branch("y Positron", &y_e);
    eventTree->Branch("z Positron", &z_e);

    for (int ii = 0; ii < nevents; ii++) {

	// generate proton momentum direction
	theta_p = acos(rnd->Uniform(-1, 1));
	phi_p = rnd->Uniform(0, 2 * pi);
	momentum_p = p_p;

	x_p = p_p * sin(theta) * cos(phi);
	y_p = p_p * sin(theta) * sin(phi);
	z_p = p_p * cos(theta);

	// get boost vector of proton
	TLorentzVector protonVect(x_p, y_p, z_p, sqrt(EE2));
	// TLorentzVector protonVect(p_p, 0, 0, sqrt(EE2));
	auto pBoostVect = protonVect.BoostVector();

	// boost proton into rest frame
	protonVect.Boost(-pBoostVect);

	// calculate decay products in proton rest frame
	theta = acos(rnd->Uniform(-1, 1));
	phi = rnd->Uniform(0, 2 * pi);
	momentum = pp;

	// convert spherical coordinates to cartesian
	xx = pp * sin(theta) * cos(phi);
	yy = pp * sin(theta) * sin(phi);
	zz = pp * cos(theta);

	const double ER_pi2 = pp * pp + m_pi * m_pi;
	TLorentzVector piVect(xx, yy, zz, sqrt(ER_pi2));

	const double ER_e2 = pp * pp + m_e * m_e;
	TLorentzVector posVect(-xx, -yy, -zz, sqrt(ER_e2));

	// transform to lab frame
	piVect.Boost(pBoostVect);
	posVect.Boost(pBoostVect);

	x_pi = piVect.Px();
	y_pi = piVect.Py();
	z_pi = piVect.Pz();

	x_e = posVect.Px();
	y_e = posVect.Py();
	z_e = posVect.Pz();

	htheta->Fill(theta);
	hphi->Fill(phi);
	h2->Fill(theta, phi);

	hxyz->Fill(x_pi, y_pi, z_pi);

	eventTree->Fill();
    }

    /*
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
    */

    // Store pi0 events in tree
    TFile *file = TFile::Open("events.root", "RECREATE");
    lout->Write();
    eventTree->Write();
    file->Save();
    file->Close();
}
			
			
