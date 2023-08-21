void transverseplots(void)
{

  //A C++/ROOT script to study neutron behavior into the moderator
  //The algorithm was a bit too complicated to include in the other scripts
  //Otherwise, consider this an extension of dedxplots.C


  TFile input("NeutronKinematics.root","READ");
  TTree *fTree = (TTree*)input.Get("Steps");
  int entries = fTree->GetEntries();


  //OUTPUT1:
  TH1D *transversereflectedhisto = new TH1D("transversereflectedhisto","Horizontal distance travelled in shield (reflected neutrons); Horizontal distance (mm); #Neutrons",101,0,101);
  TH1D *transversetransmittedhisto = new TH1D("transversetransmittedhisto","Horizontal distance travelled in shield (transmitted neutrons); Horizontal distance (mm); #Neutrons",101,0,101);
  //TH1D *pathlengthstepslarhisto = new TH1D("pathlengthstepslarhisto","Path length until interaction (no boundary crossings); Path length (m); #Neutrons",600,0,6);

  //OUTPUT2:
  //TH2D *pathlengthstepsvselarhisto = new TH2D("pathlengthstepsvselarhisto","Free path length vs. kinetic energy;Path length (m);Kinetic energy (MeV)",600,0,6,90,10,100);

  
  int stepnumber = 0;
  int tracknumber = 0;
  int eventnumber = 0;
  double X = 0;
  double Y = 0;
  double Z = 0;
  double energy = 0;
  Char_t materialc[100];
  Char_t processc[100];
  Char_t volumec[100];
  
  fTree->SetBranchAddress("NeutronStepTrackID",&tracknumber);
  fTree->SetBranchAddress("NeutronStepID",&stepnumber);
  fTree->SetBranchAddress("NeutronStepX",&X);//mm
  fTree->SetBranchAddress("NeutronStepY",&Y);//mm
  fTree->SetBranchAddress("NeutronStepZ",&Z);//mm
  fTree->SetBranchAddress("NeutronStepMaterial",&materialc);
  fTree->SetBranchAddress("NeutronStepProcess",&processc);
  fTree->SetBranchAddress("NeutronStepVolume",&volumec);
  fTree->SetBranchAddress("NeutronStepKineticEnergy",&energy);//keV
  
  int firsttracknumber = 0;
  int secondtracknumber = 0;
  double radiusone = 0;
  double radiustwo = 0;
  double initialenergy = 0;
  double initialradius = 0;
  double biggesttransverse = 0;
  
  bool insideneutron = false;
  bool outsideneutron = false;

  double X1 = 0;
  double Y1 = 0;
  double Z1 = 0;
  double X2 = 0;
  double Y2 = 0;
  double Z2 = 0;//Used for path length calculations
  double energyone = 0;
  double energytwo = 0;
  
  double radius = 0;  
  //There's some weird trick with the multithread ntuple and the way it stores 'strings'
  //For some reason it stores all strings as character arrays
  //Obviously, C++ throws a tantrum about this, so we have to convert each entry to a string
  //after reading it if we want to operate on it like a string
  string materialone;
  string materialtwo;
  string processone;
  string processtwo;
  string volumeone;
  string volumetwo;
  
  for (int j = 1; j<(entries-1); j++)
    {
      fTree->GetEntry(j);
      secondtracknumber = tracknumber;      
      materialtwo = materialc;
      radiustwo = sqrt(X*X+Y*Y);
      processtwo = processc;
      volumetwo = volumec;
      energytwo = energy;
      
      //Make sure the vertices before/after belong to the same neutron
      //This doesn't account for the vanishingly rare case where the last neutron
      //in one event shares a tracknumber with the first neutron in the next event
      fTree->GetEntry(j-1);	  
      firsttracknumber = tracknumber;
      materialone = materialc;
      radiusone = sqrt(X*X+Y*Y);
      processone = processc;
      volumeone = volumec;
      energyone = energy;
      
      if(firsttracknumber == secondtracknumber)
	{//The same neutron

	  if(materialone!=materialtwo)
	    {//Boundary crossing - materialone is the material before boundary change

	      if(materialone == "G4_lAr" && materialtwo == "PMMA" && Z < 300 && Z > -1500)
		{//Neutron enters the side of the shield (boundaries at 2000 & 2100 mm) from the LAr
		  //I'll call this an 'assigned neutron' in other comments		  

		  if(radiusone < 2050)//Inside neutron
		    insideneutron = true;
		  else if(radiusone > 2050)
		    outsideneutron = true;

		  initialenergy = energyone;
		  initialradius = radiusone;
		  //cout << materialone << "   " << materialtwo << endl;
		  //cout << volumeone << "   " << volumetwo << endl;
		  //cout << processone << "   " << processtwo << endl;
		  //cout << radiusone << "   " << radiustwo << endl;
		  
		}
	      
	      else if(materialone == "PMMA" && (insideneutron || outsideneutron) && Z < 300 && Z > -1500)
		{//An assigned neutron has completed its journey
		  //Keeping the Z requirement because neutrons entering the side and exiting the top mucks the numbers

		  //Tally the results, based on if the neutron was reflected or transmitted through the shield

		  //Check one last time for transverse calculations
		  if(initialradius && biggesttransverse < (abs(initialradius - radiusone)))
		    biggesttransverse = abs(initialradius - radiusone);
                     
		  if((insideneutron && radiusone < 2050) || (outsideneutron && radiusone > 2050))//Reflected
			{
			  transversereflectedhisto->Fill(biggesttransverse);
			}
		      else if((insideneutron && radiusone > 2050) || (outsideneutron && radiusone < 2050))//Transmitted
			{
			  transversetransmittedhisto->Fill(biggesttransverse);
			}
		      else
			{
			cerr << "Ambiguous neutron - check assignment algorithm" << endl;
			cout << materialone << "   " << materialtwo << endl;
			cout << volumeone << "   " << volumetwo << endl;
			cout << processone << "   " << processtwo << endl;
			cout << setprecision(16) << initialradius << "   " <<radiusone << "   " << radiustwo << endl;
			cout << insideneutron << "   " << outsideneutron << endl << endl;;
			}
		    
		  
		  insideneutron = false;
		  outsideneutron = false;
		  initialenergy = 0;
		  initialradius = 0;
		  biggesttransverse = 0;
		}

	      else
		{//For safety and edge cases, reset tallies one more time
		  insideneutron = false;
                  outsideneutron = false;
                  initialenergy = 0;
                  initialradius = 0;
                  biggesttransverse = 0;
		}
	      
	    }

	  
	  
	  else//Same neutron in same material, perform step-level checks
	    {
	      if(initialradius && biggesttransverse < (abs(initialradius - radiusone)))
		biggesttransverse = abs(initialradius - radiusone);
	      //Basically, if an assigned neutron has travelled further into the shield horizontally	      
	    }

	}//Same neutron


      
      else//Different neutron
	{

	  //Make sure ALL tallies are reset
	  initialenergy = 0;
	  initialradius = 0;
	  insideneutron = false;
	  outsideneutron = false;
	  biggesttransverse = 0;
	  
	}//Different neutron

      	
    }//Loop over all tree entries

  
  
  gStyle->SetOptStat(0000000);
  //h->SetDirectory(gROOT);  
  //TCanvas *c = new TCanvas;
  //c->SetLogy();
  //h->Draw();
  //h->->SetMarkerColor();
  //h->SetMarkerStyle(21);
  //h->SetMarkerSize(.2);
  //h->->SetLineColor();
  //h->SetLineWidth();
  // TLegend *l = new TLegend(0.7,0.77,0.9,0.9);
  //l->AddEntry(,"");
  //l->Draw("same");


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

  transversereflectedhisto->SetDirectory(gROOT);
  transversetransmittedhisto->SetDirectory(gROOT);
  
  TCanvas *transversecanvas = new TCanvas;
  transversecanvas->SetLogy();
  transversereflectedhisto->SetLineColor(kBlue);
  transversetransmittedhisto->SetLineColor(kRed);
  transversereflectedhisto->SetLineWidth(3);
  transversetransmittedhisto->SetLineWidth(3);
  transversereflectedhisto->Draw();
  transversetransmittedhisto->Draw("same");

  TLegend *transverselegend = new TLegend(0.7,0.77,0.9,0.9);
  transverselegend->AddEntry(transversereflectedhisto,"Reflected");
  transverselegend->AddEntry(transversetransmittedhisto,"Transmitted");
  transverselegend->Draw("same");

  
}//EOF
