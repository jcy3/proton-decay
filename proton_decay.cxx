
void proton_decay() {
    double pi = TMath::Pi();
     
    // Particle masses in MeV/c2
    double m_p = 938.27;   // proton
    
    double m_pi = 134.98;  // pion
    double m_e = 0.511;  // positron

    // Calculate momentum
    double p;
    {
	 double m_p2 = m_p * m_p;
	 double m_pi2 = m_pi * m_pi;
	 double m_e2 = m_e * m_e;

	 double frac = (m_p2 + m_e2 - m_pi2) / (2 * m_p);
	 p = sqrt(frac * frac - m_e2);
    }
    std::cout << "Momentum: " << p << " MeV/c" << "\n";

    TRandom3 *rnd = new TRandom3();
    rnd->SetSeed();

    std::vector<double> theta( 1000 );
    std::vector<double> phi( 1000 );

    // Generate theta and phi for 1000 events
    // theta, phi are uniformly distributed
    for (int i = 0; i < 1000; i++) {
	 theta[i] = rnd->Uniform(0, 2 * pi);
	 phi[i] = rnd->Uniform(0, 2 * pi);
    }

    TH1D *htheta = new TH1D("htheta", "theta histogram", 100, 0, 2 * pi);
    TH1D *hphi = new TH1D("hphi", "phi histogram", 100, 0, 2 * pi);
    TH2D *h2 = new TH2D("h2", "theta vs phi", 100, 0, 2 * pi, 100, 0, 2 * pi);

    for (int i = 0; i < 1000; i++) {
	 htheta->Fill(theta[i]);
	 hphi->Fill(phi[i]);
	 h2->Fill(theta[i], phi[i]);
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

    // Store pi0 events in tree
    TTree *events = new TTree("eventTree", "Events");
    auto branch = events->Branch("theta", &theta, 1000, 0);
    
}
			
			
