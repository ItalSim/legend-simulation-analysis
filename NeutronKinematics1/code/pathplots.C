void pathplots(void)
{

  //A C++/ROOT script to study neutron behavior as a function of distance travelled, or distance between interactions


  TFile input("NeutronKinematics.root","READ");
  TTree *fTree = (TTree*)input.Get("Steps");
  int entries = fTree->GetEntries();


  //OUTPUT1: 1D plots of pathlengths (step level)
  TH1D *pathlengthstepslarhisto = new TH1D("pathlengthstepslarhisto","Path length until interaction (no boundary crossings); Path length (m); #Neutrons",600,0,6);
  TH1D *pathlengthstepspmmahisto = new TH1D("pathlengthstepspmmahisto","Path length until interaction (no boundary crossings); Path length (m); #Neutrons",600,0,6);

  //OUTPUT2: 1D plots of total pathlengths (cumulative track level)
  TH1D *pathlengthtotallarhisto = new TH1D("pathlengthtotallarhisto","Total path length of each neutron; Path length (m); #Neutrons",1000,10,50);
  TH1D *pathlengthtotalpmmahisto = new TH1D("pathlengthtotalpmmahisto","Total path length of each neutron; Path length (m); #Neutrons",1000,10,50);
  TH1D *pathlengthtotalbothhisto = new TH1D("pathlengthtotalbothhisto","Total path length of each neutron; Path length (m); #Neutrons",1000,10,50);

  //OUTPUT3: 2D plots of pathlengths vs kinetic energy
  TH2D *pathlengthstepsvselarhisto = new TH2D("pathlengthstepsvselarhisto","Free path length vs. kinetic energy;Path length (m);Kinetic energy (MeV)",600,0,6,90,10,100);
  TH2D *pathlengthstepsvsepmmahisto = new TH2D("pathlengthstepsvsepmmahisto","Free path length vs. kinetic energy;Path length (m);Kinetic energy (MeV)",600,0,6,90,10,100);

  
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
  double initialenergy = 0;//In this script, energy of neutron upon first entering any relevant volume

  double X1 = 0;
  double Y1 = 0;
  double Z1 = 0;
  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;//Used for path length calculation
  double pathlengthofstep = 0;

    
  double totallarpath = 0;//Cumulative path length of neutron in LAr
  double totalpmmapath = 0;//Cumulative path length of neutron in PMMA
  
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
  
  for (int j = 1; j<(entries-1); j++)
    {
      fTree->GetEntry(j);
      secondtracknumber = tracknumber;      
      X2 = X;
      Y2 = Y;
      Z2 = Z;
      
      //Make sure the vertices before/after belong to the same neutron
      //This doesn't account for the vanishingly rare case where the last neutron
      //in one event shares a tracknumber with the first neutron in the next event
      fTree->GetEntry(j-1);	  
      firsttracknumber = tracknumber;
      materialone = materialc;
      X1 = X;
      Y1 = Y;
      Z1 = Z;

      pathlengthofstep = sqrt((X2-X1)*(X2-X1)+(Y2-Y1)*(Y2-Y1)+(Z2-Z1)*(Z2-Z1));
      pathlengthofstep = pathlengthofstep/1000.;
      //If this is a neutron not previously seen before in the LAr or shield
      //(we know it hasn't been seen before if initialenergy is not assigned)
      if((materialone == "G4_lAr" || materialone == "PMMA") && !initialenergy)
	initialenergy = energy;


	 if(firsttracknumber == secondtracknumber)
	{//The same neutron

	  fTree->GetEntry(j);
	  process = processc;
	  materialtwo = materialc;
	  radius = sqrt(X*X+Y*Y);


	  if(materialone!=materialtwo)
	    {//Boundary crossing, but materialone is the material of THIS step
	      currentsteps++;

	      //Fill any histograms that tally parameters in a certain volume
	      
	      if(materialone == "G4_lAr")  //Add to LAr tallies
		  totallarpath += pathlengthofstep;
		  		
	      else if(materialone == "PMMA")//Add to PMMA tallies
		  totalpmmapath += pathlengthofstep;
	      
	      //Reset tallies that depend on current material remaining consistent
	      currentsteps = 0;

	    }//Different materials

	  
	  else//Same neutron in same material, tally as usual
	    {

	      if(materialone == "G4_lAr")
		{//Add to track-levelLAr tallies/fill step-level LAr histos
		  totallarpath += pathlengthofstep;
		  pathlengthstepslarhisto->Fill(pathlengthofstep);
		  pathlengthstepsvselarhisto->Fill(pathlengthofstep,energy);
		  //if(pathlengthofstep>4.5)//Diagnostics
		  //cout << pathlengthofstep << endl << X1 << "  " << Y1 << "  " << Z1 << "  " << X2 << "  " << Y2 << "  " << Z2 << endl << endl;
		}
		
	      else if(materialone == "PMMA")
		{//Add to track-level PMMA tallies/fill step-level PMMA histos
		  totalpmmapath += pathlengthofstep;
		  pathlengthstepspmmahisto->Fill(pathlengthofstep);
		  pathlengthstepsvsepmmahisto->Fill(pathlengthofstep,energy);
		  //if(pathlengthofstep>4.5)//Diagnostics
		  //cout << pathlengthofstep << endl << X1 << "  " << Y1 << "  " << Z1 << "  " << X2 << "  " << Y2 << "  " << Z2 << endl << endl;

		}
	      
	      currentsteps++;
	      
	    }
	  
	  
       	}//Same neutron

      
      else//Different neutron
	{
	  currentsteps++;
	  //Fill track-level histos

	  if(totallarpath)
	    pathlengthtotallarhisto->Fill(totallarpath);
	  if(totalpmmapath)
	    pathlengthtotalpmmahisto->Fill(totalpmmapath);	  
	  //Recycle one variable since it will be reset soon anyway
	  totalpmmapath+=totallarpath;
	  if(totalpmmapath)
	    pathlengthtotalbothhisto->Fill(totalpmmapath);	  

	  //Make sure ALL tallies are reset
	  currentsteps = 0;
	  initialenergy = 0;
	  totallarpath = 0;
	  totalpmmapath = 0;
	  
	}//Different neutron
      
	
    }//Loop over all tree entries

    
  gStyle->SetOptStat(0000);

  //OUTPUT1
  /*
  pathlengthstepslarhisto->SetDirectory(gROOT);
  pathlengthstepspmmahisto->SetDirectory(gROOT);
  
  TCanvas *pathlengthstepscanvas = new TCanvas;
  pathlengthstepscanvas->SetLogy();
  pathlengthstepslarhisto->SetLineColor(kTeal+3);
  pathlengthstepspmmahisto->SetLineColor(kGray+2);
  pathlengthstepslarhisto->SetLineWidth(3);
  pathlengthstepspmmahisto->SetLineWidth(3);
  pathlengthstepspmmahisto->Draw();
  pathlengthstepslarhisto->Draw("same");

  TLegend *pathlengthstepslegend = new TLegend(0.7,0.77,0.9,0.9);
  pathlengthstepslegend->AddEntry(pathlengthstepslarhisto,"LAr");
  pathlengthstepslegend->AddEntry(pathlengthstepspmmahisto,"PMMA");
  pathlengthstepslegend->Draw("same");
  */

  //OUTPUT2
  /*
  pathlengthtotallarhisto->SetDirectory(gROOT);
  pathlengthtotalpmmahisto->SetDirectory(gROOT);
  pathlengthtotalbothhisto->SetDirectory(gROOT);
  
  TCanvas *pathlengthtotalcanvas = new TCanvas;
  pathlengthtotalcanvas->SetLogy();
  pathlengthtotallarhisto->SetLineColor(kTeal+3);
  pathlengthtotalpmmahisto->SetLineColor(kGray+2);
  pathlengthtotalbothhisto->SetLineColorAlpha(kRed,0.3);
  pathlengthtotallarhisto->SetLineWidth(3);
  pathlengthtotalpmmahisto->SetLineWidth(3);
  Pathlengthtotalbothhisto->SetLineWidth(3);
  pathlengthtotalbothhisto->Draw();
  pathlengthtotalpmmahisto->Draw("same");
  pathlengthtotallarhisto->Draw("same");


  TLegend *pathlengthtotallegend = new TLegend(0.7,0.77,0.9,0.9);
  pathlengthtotallegend->AddEntry(pathlengthtotallarhisto,"LAr");
  pathlengthtotallegend->AddEntry(pathlengthtotalpmmahisto,"PMMA");
  pathlengthtotallegend->AddEntry(pathlengthtotalbothhisto,"LAr+PMMA");
  pathlengthtotallegend->Draw("same");
  */

  //OUTPUT3
  /*
  pathlengthstepsvselarhisto->SetDirectory(gROOT);
  pathlengthstepsvsepmmahisto->SetDirectory(gROOT);

  TCanvas *pathlengthstepsvsecanvas = new TCanvas;
  //pathlengthstepsvsecanvas->SetLogy();
  pathlengthstepsvselarhisto->SetMarkerColor(kTeal+3);
  pathlengthstepsvsepmmahisto->SetMarkerColor(kGray+2);
  pathlengthstepsvselarhisto->SetMarkerStyle(21);
  pathlengthstepsvsepmmahisto->SetMarkerStyle(21);
  pathlengthstepsvselarhisto->SetMarkerSize(0.5);
  pathlengthstepsvsepmmahisto->SetMarkerSize(0.5);
  pathlengthstepsvselarhisto->Draw();
  pathlengthstepsvsepmmahisto->Draw("same");

  TLegend *pathlengthstepsvselegend = new TLegend(0.7,0.77,0.9,0.9);
  pathlengthstepsvselegend->AddEntry(pathlengthstepsvselarhisto,"LAr");
  pathlengthstepsvselegend->AddEntry(pathlengthstepsvsepmmahisto,"PMMA");
  pathlengthstepsvselegend->Draw("same");
  */

  
}//EOF
