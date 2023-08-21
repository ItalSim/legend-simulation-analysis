void steplevelplots(void)
{

  //A ROOT/C++ script for the neutron kinematic study which tallies step-level information, without paying particular attention to only boundary crossings
  
  TFile input("NeutronKinematics.root","READ");
  TTree *fTree = (TTree*)input.Get("Steps");
  int entries = fTree->GetEntries();


  //OUTPUT1: Consecutive steps
  TH1D *consecutivestepsinlarhisto = new TH1D("consecutivestepsinlarhisto","Number of uninterrupted steps in material before absorption or exiting; #Steps; #Neutrons",500,0,500);
  TH1D *consecutivestepsinpmmahisto = new TH1D("consecutivestepsinpmmahisto","Number of uninterrupted steps in material before absorption or exiting; #Steps; #Neutrons",500,0,500);

  //OUTPUT2: DeltaE/E for neutrons in the LAr or in the shield, separately
  TH2D *energylossvselarhisto = new TH2D("energylossvselarhisto","Fraction of energy lost vs. kinetic energy;Energy lost (%);Kinetic energy upon entering the material (keV)",100,0,100,500,0,500);
  TH2D *energylossvsepmmahisto = new TH2D("energylossvsepmmahisto","Fraction of energy lost vs. kinetic energy;Energy lost (%);Kinetic energy upon entering the material (keV)",100,0,100,500,0,500);

  //OUTPUT3: Consecutive steps again, but this time with reflection and transmission
  TH1D *consecutivestepsreflectedhisto = new TH1D("consecutivestepsreflectedhisto","Number of uninterrupted steps in shield; #Steps; #Neutrons",500,0,500);
  TH1D *consecutivestepstransmittedhisto = new TH1D("consecutivestepstransmittedhisto","Number of uninterrupted steps in shield; #Steps; #Neutrons",500,0,500);


  
  int stepnumber = 0;
  int tracknumber = 0;
  int eventnumber = 0;
  double X = 0;
  double Y = 0;
  double Z = 0;
  double energy = 0;
  Char_t materialc[100];
  Char_t processc[100];
  
  fTree->SetBranchAddress("NeutronStepTrackID",&tracknumber);
  fTree->SetBranchAddress("NeutronStepID",&stepnumber);
  fTree->SetBranchAddress("NeutronStepX",&X);//mm
  fTree->SetBranchAddress("NeutronStepY",&Y);//mm
  fTree->SetBranchAddress("NeutronStepZ",&Z);//mm
  fTree->SetBranchAddress("NeutronStepMaterial",&materialc);
  fTree->SetBranchAddress("NeutronStepProcess",&processc);
  fTree->SetBranchAddress("NeutronStepKineticEnergy",&energy);//keV
  
  int firsttracknumber = 0;
  int secondtracknumber = 0;
  double radius = 0;
  int currentsteps = 0;//Number of steps in current material (changes for each boundary crossing)
  double initialenergy = 0;
  double deltaEoverE = 0;

  bool insideneutron = false;
  bool outsideneutron = false;
  
  //There's some weird trick with the multithread ntuple and the way it stores 'strings'
  //For some reason it stores all strings as character arrays
  //Obviously, C++ throws a tantrum about this, so we have to convert each entry to a string
  //after reading it if we want to operate on it like a string
  string materialone;
  string materialtwo;
  string process;

  //Prime the comparators that need priming
  fTree->GetEntry(0);
  lasttracknumber = tracknumber;
  initialenergy = energy;
  
  for (int j = 1; j<(entries-1); j++)
    {
      fTree->GetEntry(j);
      secondtracknumber = tracknumber;

      //Make sure the vertices before/after belong to the same neutron
      //This doesn't account for the vanishingly rare case where the last neutron
      //in one event shares a tracknumber with the first neutron in the next event
      fTree->GetEntry(j-1);	  
      firsttracknumber = tracknumber;


      if(firsttracknumber == secondtracknumber)
	{//The same neutron

	  materialone = materialc;

	  fTree->GetEntry(j);
	  process = processc;
	  materialtwo = materialc;
	  radius = sqrt(X*X+Y*Y);


	  if(materialone!=materialtwo)
	    {//Boundary crossing, but materialone is the material of THIS step
	      currentsteps++;
	      deltaEoverE = 100*(1 - energy/initialenergy);

	      //Fill any histograms that tally parameters in a certain volume
	      if(materialone == "G4_lAr")
		{
		  consecutivestepsinlarhisto->Fill(currentsteps);
		  energylossvselarhisto->Fill(deltaEoverE,initialenergy);

		  if(materialtwo == "PMMA")//More reflection/transmission stuff, unfortunately
		    {
		      if(radius < 2050 && Z > -1630 && Z < 470)//Inside neutron
			insideneutron = true;
		      else
			outsideneutron = true;
		    }
		  
		}
	      else if(materialone == "PMMA")
		{
		  consecutivestepsinpmmahisto->Fill(currentsteps);
		  energylossvsepmmahisto->Fill(deltaEoverE,initialenergy);

		  //Transmission/reflection
		  if(materialtwo == "G4_lAr")
		    {
		      if(radius < 2050 && Z > -1630 && Z < 470)//Inside neutron
			{
			  if(insideneutron)//Reflection
			    consecutivestepsreflectedhisto->Fill(currentsteps);
			  else if(outsideneutron)//Transmission
			    consecutivestepstransmittedhisto->Fill(currentsteps);
			}
		      else//Outside neutron
			{
			  if(insideneutron)//Transmission
			    consecutivestepstransmittedhisto->Fill(currentsteps);
			  else if(outsideneutron)//Reflection
			    consecutivestepsreflectedhisto->Fill(currentsteps);
			}
		    }
		  insideneutron = false;
		  outsideneutron = false;
		}

	      
	      //Reset tallies
	      currentsteps = 0;
	      initialenergy = energy;
	    }

	  else//Same neutron in same material, tally as usual
	      currentsteps++;
	  
	  
       	}//Same neutron


      
      else//Different neutron
	{
	  currentsteps++;
	  //Fill any histograms that need filling one last time

	  if(materialone == "G4_lAr")
	      consecutivestepsinlarhisto->Fill(currentsteps);	  
	  else if(materialone == "PMMA")
	    consecutivestepsinpmmahisto->Fill(currentsteps);

	  

	  //Make sure ALL tallies are reset
	  currentsteps = 0;
	  initialenergy = energy;
	  deltaEoverE = 0;
	}//Different neutron

      
	
    }//Loop over all tree entries

  
  gStyle->SetOptStat(0000);

  //OUTPUT1
  /*
  consecutivestepsinlarhisto->SetDirectory(gROOT);
  consecutivestepsinpmmahisto->SetDirectory(gROOT);  

  TCanvas *consecutivestepscanvas = new TCanvas;

  consecutivestepsinlarhisto->SetLineColor(kTeal+3);
  consecutivestepsinpmmahisto->SetLineColor(kGray+2);
  consecutivestepsinlarhisto->SetLineWidth(5);
  consecutivestepsinpmmahisto->SetLineWidth(5);
  
  consecutivestepsinlarhisto->Draw();
  consecutivestepsinpmmahisto->Draw("same");
  consecutivestepscanvas->SetLogy();

  TLegend *consecutivestepslegend = new TLegend(0.7,0.77,0.9,0.9);
  consecutivestepslegend->AddEntry(consecutivestepsinlarhisto,"LAr");
  consecutivestepslegend->AddEntry(consecutivestepsinpmmahisto,"PMMA");

  consecutivestepslegend->Draw("same");
  */
  //OUTPUT1

  //OUTPUT2
  /*
  energylossvselarhisto->SetDirectory(gROOT);
  energylossvsepmmahisto->SetDirectory(gROOT);

  TCanvas *energylossvsecanvas = new TCanvas;

  energylossvselarhisto->SetMarkerColor(kTeal+3);
  energylossvsepmmahisto->SetMarkerColor(kGray+2);
  energylossvselarhisto->SetMarkerStyle(21);
  energylossvsepmmahisto->SetMarkerStyle(21);
  energylossvselarhisto->SetMarkerSize(.2);
  energylossvsepmmahisto->SetMarkerSize(.2);

  
  //energylossvselarhisto->Draw();
  energylossvsepmmahisto->Draw("same");
  //energylossvsecanvas->SetLogy();

  TLegend *energylossvselegend = new TLegend(0.7,0.77,0.9,0.9);
  //energylossvselegend->AddEntry(energylossvselarhisto,"LAr");
  energylossvselegend->AddEntry(energylossvsepmmahisto,"PMMA");

  energylossvselegend->Draw("same");
  */
  //OUTPUT2

  //OUTPUT3
  /*  
  consecutivestepsreflectedhisto->SetDirectory(gROOT);
  consecutivestepstransmittedhisto->SetDirectory(gROOT);

  TCanvas *consecutivestepscanvas2 = new TCanvas;

  consecutivestepsreflectedhisto->SetLineColor(kBlue-4);
  consecutivestepstransmittedhisto->SetLineColor(kRed-3);
  consecutivestepsreflectedhisto->SetLineWidth(3);
  consecutivestepstransmittedhisto->SetLineWidth(3);

  consecutivestepsreflectedhisto->Draw();
  consecutivestepstransmittedhisto->Draw("same");
  consecutivestepscanvas2->SetLogy();

  TLegend *consecutivestepslegend2 = new TLegend(0.7,0.77,0.9,0.9);
  consecutivestepslegend2->AddEntry(consecutivestepsreflectedhisto,"Reflected");
  consecutivestepslegend2->AddEntry(consecutivestepstransmittedhisto,"Transmitted");

  consecutivestepslegend2->Draw("same");

  cout << "TOTAL REFLECTED NEUTRONS (new): " << consecutivestepsreflectedhisto->Integral() << endl;
  cout << "TOTAL TRANSMITTED NEUTRONS (new): " << consecutivestepstransmittedhisto->Integral() << endl;
  */
  //OUTPUT3
  
}//EOF
