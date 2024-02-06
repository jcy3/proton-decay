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

    // proton has initial momentum
    // gaussian magnitude
    // isotropic direction

    const double p_p_mean = 200;
    const double p_p_sigma = 50;

    // Calculate momentum of proton decay products
    const double m_p2 = m_p * m_p;
    const double m_pi2 = m_pi * m_pi;
    const double m_e2 = m_e * m_e;

    const double frac = (m_p2 + m_e2 - m_pi2) / (2 * m_p);
    const double pp = sqrt(frac * frac - m_e2);

    // Calculate momentum of pion decay products
    const double p_gamma = m_pi / 2;    // in pion rest frame

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

    TH3D *hxyzproton = new TH3D("hxyzproton", "Distribution of proton momentum in 3D;x;y;z", 100, 1, -0, 100, 1, 0, 100, 1, 0);
    lout->Add(hxyzproton);

    TH1D *hxproton = new TH1D("hxproton", "Distribution of x component of pion momentum", 100, 1, -0);
    lout->Add(hxproton);
    TH1D *hyproton = new TH1D("hyproton", "Distribution of y component of pion momentum", 100, 1, -0);
    lout->Add(hyproton);
    TH1D *hzproton = new TH1D("hzproton", "Distribution of z component of pion momentum", 100, 1, -0);
    lout->Add(hzproton);

    TH3D *hxyzpi = new TH3D("hxyzpi", "Distribution of pion momentum in 3D;x;y;z", 100, 1, -0, 100, 1, 0, 100, 1, 0);
    lout->Add(hxyzpi);

    TH1D *hxpi = new TH1D("hxpi", "Distribution of x component of pion momentum", 100, 1, -0);
    lout->Add(hxpi);
    TH1D *hypi = new TH1D("hypi", "Distribution of y component of pion momentum", 100, 1, -0);
    lout->Add(hypi);
    TH1D *hzpi = new TH1D("hzpi", "Distribution of z component of pion momentum", 100, 1, -0);
    lout->Add(hzpi);

    TH3D *hxyzgamma = new TH3D("hxyzgamma", "Distribution of photon momentum in 3D;x;y;z", 100, 1, 0, 100, 1, 0, 100, 1, 0);
    lout->Add(hxyzgamma);

    TH3D *hxyzgamma1 = new TH3D("hxyzgamma1", "Distribution of photon momentum in 3D;x;y;z", 100, 1, 0, 100, 1, 0, 100, 1, 0);
    lout->Add(hxyzgamma1);

    TH3D *hxyzgamma2 = new TH3D("hxyzgamma2", "Distribution of photon momentum in 3D;x;y;z", 100, 1, 0, 100, 1, 0, 100, 1, 0);
    lout->Add(hxyzgamma2);

    TH1D *hxgamma1 = new TH1D("hxgamma1", "Distribution of x component of photon momentum", 100, 1, -0);
    lout->Add(hxgamma1);
    TH1D *hygamma1 = new TH1D("hygamma1", "Distribution of y component of photon momentum", 100, 1, -0);
    lout->Add(hygamma1);
    TH1D *hzgamma1 = new TH1D("hzgamma1", "Distribution of z component of photon momentum", 100, 1, -0);
    lout->Add(hzgamma1);

    TH1D *hxgamma2 = new TH1D("hxgamma2", "Distribution of x component of photon momentum", 100, 1, -0);
    lout->Add(hxgamma2);
    TH1D *hygamma2 = new TH1D("hygamma2", "Distribution of y component of photon momentum", 100, 1, -0);
    lout->Add(hygamma2);
    TH1D *hzgamma2 = new TH1D("hzgamma2", "Distribution of z component of photon momentum", 100, 1, -0);
    lout->Add(hzgamma2);

    // initialise variables

    double theta = -999;
    double phi = -999;
    double momentum = -999;

    double xx = -999;
    double yy = -999;
    double zz = -999;

    double theta_p = -999;
    double phi_p = -999;
    double p_p = -999;
    
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

    double theta_gamma = -999;
    double phi_gamma = -999;

    double x_gamma = -999;
    double y_gamma = -999;
    double z_gamma = -999;

    double x_gamma1 = -999;
    double y_gamma1 = -999;
    double z_gamma1 = -999;

    double x_gamma2 = -999;
    double y_gamma2 = -999;
    double z_gamma2 = -999;

    // create event tree
    TTree *eventTree = new TTree("eventTree", "Events");

    eventTree->Branch("Proton Momentum", &p_p);
    eventTree->Branch("x Proton", &x_p);
    eventTree->Branch("y Proton", &y_p);
    eventTree->Branch("z Proton", &z_p);

    eventTree->Branch("x Pi", &x_pi);
    eventTree->Branch("y Pi", &y_pi);
    eventTree->Branch("z Pi", &z_pi);

    eventTree->Branch("x Positron", &x_e);
    eventTree->Branch("y Positron", &y_e);
    eventTree->Branch("z Positron", &z_e);

    eventTree->Branch("x Photon 1", &x_gamma1);
    eventTree->Branch("y Photon 1", &y_gamma1);
    eventTree->Branch("z Photon 1", &z_gamma1);

    eventTree->Branch("x Photon 2", &x_gamma2);
    eventTree->Branch("y Photon 2", &y_gamma2);
    eventTree->Branch("z Photon 2", &z_gamma2);


    for (int ii = 0; ii < nevents; ii++) {

	// generate proton momentum magnitude
	p_p = rnd->Gaus(p_p_mean, p_p_sigma);
	const double E_p2 = p_p * p_p + m_p * m_p;

	// generate proton momentum direction
	theta_p = acos(rnd->Uniform(-1, 1));
	phi_p = rnd->Uniform(0, 2 * pi);

	x_p = p_p * sin(theta_p) * cos(phi_p);
	y_p = p_p * sin(theta_p) * sin(phi_p);
	z_p = p_p * cos(theta_p);

	// get boost vector of proton
	TLorentzVector protonVect(x_p, y_p, z_p, sqrt(E_p2));
	auto pBoostVect = protonVect.BoostVector();

	// calculate decay products in proton rest frame
	theta_pi = acos(rnd->Uniform(-1, 1));
	phi_pi = rnd->Uniform(0, 2 * pi);
	momentum = pp;

	// convert spherical coordinates to cartesian
	xx = pp * sin(theta_pi) * cos(phi_pi);
	yy = pp * sin(theta_pi) * sin(phi_pi);
	zz = pp * cos(theta_pi);

	const double ER_pi2 = pp * pp + m_pi * m_pi;
	TLorentzVector piVect(xx, yy, zz, sqrt(ER_pi2));

	const double ER_e2 = pp * pp + m_e * m_e;
	TLorentzVector posVect(-xx, -yy, -zz, sqrt(ER_e2));

	// pion decay into two photons in pion rest frame

	// get boost vector of pion
	auto piBoostVect = piVect.BoostVector();
	// generate direction
	theta_gamma = acos(rnd->Uniform(-1, 1));
	phi_gamma = rnd->Uniform(0, 2 * pi);

	// calculate momentum in pion rest frame
	x_gamma = p_gamma * sin(theta_gamma) * cos(phi_gamma);
	y_gamma = p_gamma * sin(theta_gamma) * sin(phi_gamma);
	z_gamma = p_gamma * cos(theta_gamma);

	TLorentzVector gammaVect1(x_gamma, y_gamma, z_gamma, p_gamma);
	TLorentzVector gammaVect2(-x_gamma, -y_gamma, -z_gamma, p_gamma);

	// transform pion and positron to lab frame
	piVect.Boost(pBoostVect);
	posVect.Boost(pBoostVect);

	x_pi = piVect.Px();
	y_pi = piVect.Py();
	z_pi = piVect.Pz();

	x_e = posVect.Px();
	y_e = posVect.Py();
	z_e = posVect.Pz();

	// for photon, transform to proton rest frame, then lab frame
	gammaVect1.Boost(piBoostVect);
	gammaVect1.Boost(pBoostVect);

	gammaVect2.Boost(piBoostVect);
	gammaVect2.Boost(pBoostVect);

	// gamma1 should have higher energy, swap photons if needed
	if (gammaVect1.E() > gammaVect2.E()) {
	    x_gamma1 = gammaVect1.Px();
	    y_gamma1 = gammaVect1.Py();
	    z_gamma1 = gammaVect1.Pz();

	    x_gamma2 = gammaVect2.Px();
	    y_gamma2 = gammaVect2.Py();
	    z_gamma2 = gammaVect2.Pz();
	} else {
	    x_gamma2 = gammaVect1.Px();
	    y_gamma2 = gammaVect1.Py();
	    z_gamma2 = gammaVect1.Pz();

	    x_gamma1 = gammaVect2.Px();
	    y_gamma1 = gammaVect2.Py();
	    z_gamma1 = gammaVect2.Pz();
	}

	// fill histograms
	htheta->Fill(theta);
	hphi->Fill(phi);
	h2->Fill(theta, phi);

	hxyzproton->Fill(x_p, y_p, z_p);

	hxproton->Fill(x_p);
	hyproton->Fill(y_p);
	hzproton->Fill(z_p);

	hxyzpi->Fill(x_pi, y_pi, z_pi);

	hxpi->Fill(x_pi);
	hypi->Fill(y_pi);
	hzpi->Fill(z_pi);

	hxyzgamma->Fill(x_gamma1, y_gamma1, z_gamma1);
	hxyzgamma->Fill(x_gamma2, y_gamma2, z_gamma2);

	hxyzgamma1->Fill(x_gamma1, y_gamma1, z_gamma1);
	hxyzgamma2->Fill(x_gamma2, y_gamma2, z_gamma2);

	hxgamma1->Fill(x_gamma1);
	hygamma1->Fill(y_gamma1);
	hzgamma1->Fill(z_gamma1);

	hxgamma2->Fill(x_gamma2);
	hygamma2->Fill(y_gamma2);
	hzgamma2->Fill(z_gamma2);
	

	eventTree->Fill();
    }

    // Store event tree and histograms in root file
    TFile *file = TFile::Open("events.root", "RECREATE");
    lout->Write();
    eventTree->Write();
    file->Save();
    file->Close();
}
