void scatteringangleplots(void)
{

  //Given 3 points in space, calculate the scattering angle which occurred at the second vertice for neutrons
  //A simple concept, but there are many many things we have to consider

  //Firstly, not all the vertices we're looking at belong to neutrons
  //Secondly, we have to make sure that the vertices we look at are physical events
  //Thirdly, we have to take care about which material the vertice is in etc

  //For three points, the important thing to do is make sure the middle point is a scatter
  //Also, that the first and third point belong to the same neutron as the middle point
  //Material discrimination can be considered later, if necessary

  //Once three candidate vertices have been found, draw an imaginary triangle between them,
  //determine the length of each side of the triangle, and use the law of cosines to find
  //the angle associated with the middle vertex
  //This will be the complementary angle of the scattering angle

  TRandom3 *rng = new TRandom3(0);
  //cout << rng->Rndm();

  //Let's see what the individual steps in order look like!
  
  //TFile input("output11758804.root","READ");
  TFile input("NeutronKinematics.root","READ");
  TTree *fTree = (TTree*)input.Get("Steps");
  int entries = fTree->GetEntries();

  TH1D *scatterhistoallall = new TH1D("scatterhistoallall",              "Scattering angle for neutrons (all energies, all volumes);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistolowall = new TH1D("scatterhistolowall",              "Scattering angle for neutrons (energy < 10 keV, all volumes);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistomedall = new TH1D("scatterhistomedall",              "Scattering angle for neutrons (energy 10 to 1000 keV, all volumes);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistohighall = new TH1D("scatterhistohighall",            "Scattering angle for neutrons (energy > 1 MeV, all volumes);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistoalllar = new TH1D("scatterhistoalllar",              "Scattering angle for neutrons (all energies, LAr volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistolowlar = new TH1D("scatterhistolowlar",              "Scattering angle for neutrons (energy < 10 keV, LAr volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistomedlar = new TH1D("scatterhistomedlar",              "Scattering angle for neutrons (energy 10 to 1000 keV, LAr volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistohighlar = new TH1D("scatterhistohighlar",            "Scattering angle for neutrons (energy > 1 MeV, LAr volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistoallshield = new TH1D("scatterhistoallshield",        "Scattering angle for neutrons (all energies, shield volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistolowshield = new TH1D("scatterhistolowshield",        "Scattering angle for neutrons (energy < 10 keV, shield volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistomedshield = new TH1D("scatterhistomedshield",        "Scattering angle for neutrons (energy 10 to 1000 keV, shield volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistohighshield = new TH1D("scatterhistohighshield",      "Scattering angle for neutrons (energy > 1 MeV, shield volume);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistoalllarshield = new TH1D("scatterhistoalllarshield",  "Scattering angle for neutrons (all energies, LAr+shield);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistolowlarshield = new TH1D("scatterhistolowlarshield",  "Scattering angle for neutrons (energy < 10 keV, LAr+shield);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistomedlarshield = new TH1D("scatterhistomedlarshield",  "Scattering angle for neutrons (energy 10 to 1000 keV, LAr+shield);Scattering angle (degrees);Counts",180,0,180);
  TH1D *scatterhistohighlarshield = new TH1D("scatterhistohighlarshield","Scattering angle for neutrons (energy > 1 MeV, LAr+shield);Scattering angle (degrees);Counts",180,0,180);

  TH2D *scatteringevsanglelow = new TH2D("scatteringevsanglelow","Kineticenergy vs scattering angle, all cryostat volumes (0.5 to 10 keV);Scattering angle (degrees);Kinetic energy (keV)",18, 0, 180, 19,0.5,10);
  TH2D *scatteringevsanglemed = new TH2D("scatteringevsanglemed","Kineticenergy vs scattering angle, all cryostat volumes (10 to 1000 keV);Scattering angle (degrees);Kinetic energy (keV)",18, 0, 180, 20,10,1000);
  TH2D *scatteringevsanglehigh = new TH2D("scatteringevsanglehigh","Kineticenergy vs scattering angle, all cryostat volumes (1 to 10 MeV with energy overflow bin);Scattering angle (degrees);Kinetic energy (MeV)",18, 0, 180, 18,1,10); 

  TH2D *deltaevsanglelowlar = new TH2D("deltaevsanglelowlar","Fraction of energy lost vs scattering angle, shield only (0.5 to 10 keV);Scattering angle (degrees);Energy lost (%)",180, 0, 180, 100,0,100);
  TH2D *deltaevsanglelowpmma = new TH2D("deltaevsanglelowpmma","Fraction of energy lost vs scattering angle, shield only (0.5 to 10 keV);Scattering angle (degrees);Energy lost (%)",180, 0, 180, 100,0,100);
  TH2D *deltaevsanglemedlar = new TH2D("deltaevsanglemedlar","Fraction of energy lost vs scattering angle, shield only (10 to 1000 keV);Scattering angle (degrees);Energy lost (%)",180, 0, 180, 100,0,100);
  TH2D *deltaevsanglemedpmma = new TH2D("deltaevsanglemedpmma","Fraction of energy lost vs scattering angle, shield only (10 to 1000 keV);Scattering angle (degrees);Energy lost (%)",180, 0, 180, 100,0,100);
  TH2D *deltaevsanglehighlar = new TH2D("deltaevsanglehighlar","Fraction of energy lost vs scattering angle, shield only (1 to 10 MeV);Scattering angle (degrees);Energy lost (%)",180, 0, 180, 100,0,100);
  TH2D *deltaevsanglehighpmma = new TH2D("deltaevsanglehighpmma","Fraction of energy lost vs scattering angle, shield only (1 to 10 MeV);Scattering angle (degrees);Energy lost (%)",180, 0, 180, 100,0,100);
  
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
  fTree->SetBranchAddress("NeutronStepX",&X);
  fTree->SetBranchAddress("NeutronStepY",&Y);
  fTree->SetBranchAddress("NeutronStepZ",&Z);
  fTree->SetBranchAddress("NeutronStepMaterial",&materialc);
  fTree->SetBranchAddress("NeutronStepProcess",&processc);
  fTree->SetBranchAddress("NeutronStepKineticEnergy",&energy);
  
  int firsttracknumber = 0;
  int secondtracknumber = 0;
  int thirdtracknumber = 0;
  double firstenergy = 0;
  double secondenergy = 0;
  double thirdenergy = 0;
  double x[3];
  double y[3];
  double z[3];
  double AB = 0;
  double BC = 0;
  double AC = 0;
  double internalangle = 0;
  double scatteringangle = 0;
  double energyMeV = 0;
  double deltaE = 0;
  bool prescatter = true;//Set to 'false' if you want the energy post-scatter
  
  //There's some weird trick with the multithread ntuple and the way it stores 'strings'
  //For some reason it stores all strings as character arrays
  //Obviously, C++ throws a tantrum about this, so we have to convert each entry to a string
  //after reading it if we want to operate on it like a string
  string material;
  string process;
  
  for (int j = 1; j<(entries-1); j++)
    {
      fTree->GetEntry(j);
      secondtracknumber = tracknumber;
      secondenergy = energy;
      //if(PID==2112&&strstr(process->c_str(),"hadElastic"))//Scattered neutron vertex
      //{
	  //Make sure the vertices before/after belong to the same neutron
	  //Only one neutron can be in each event at these energies
	  fTree->GetEntry(j-1);
	  
	  firsttracknumber = tracknumber;
	  firstenergy = energy;
	  x[0] = X;
	  y[0] = Y;
	  z[0] = Z;
	  fTree->GetEntry(j+1);
	  thirdtracknumber = tracknumber;
	  

	  if(firsttracknumber == secondtracknumber && secondtracknumber == thirdtracknumber)
	    {//All the same neutron, no change in tracks
	      //Split up the variable assignments to save a little runtime
	      x[2] = X;
	      y[2] = Y;
	      z[2] = Z;
	      fTree->GetEntry(j);
	      process = processc;
	      material = materialc;
	      if(process!="hadElastic") continue;//We only care if the middle vertex is a scatter
	      x[1] = X;
	      y[1] = Y;
	      z[1] = Z;
	      
	      //Angle calculating algorithm
	      AB = TMath::Sqrt(((x[1]-x[0])*(x[1]-x[0]))+((y[1]-y[0])*(y[1]-y[0]))+((z[1]-z[0])*(z[1]-z[0])));
	      BC = TMath::Sqrt(((x[1]-x[2])*(x[1]-x[2]))+((y[1]-y[2])*(y[1]-y[2]))+((z[1]-z[2])*(z[1]-z[2])));
	      AC = TMath::Sqrt(((x[0]-x[2])*(x[0]-x[2]))+((y[0]-y[2])*(y[0]-y[2]))+((z[0]-z[2])*(z[0]-z[2])));

	      internalangle = TMath::ACos((AB*AB + BC*BC - AC*AC)/(2*AB*BC))*TMath::RadToDeg();//Law of cosines, solved for an angle

	      scatteringangle = 180 - internalangle;
	      //cout << scatteringangle << endl;

	      //Calculate change in energy due to the scatter
	      deltaE = 100 * (1 - secondenergy/firstenergy);
	      
	      //Fill the various histograms, based on neutron energy and the material of the second vertex
	      scatterhistoallall->Fill(scatteringangle);

	      if(prescatter)//We need the first energy, not the second, to know the energy before scattering
		fTree->GetEntry(j-1);
	      
	      if(energy < 10)
		{
		  scatterhistolowall->Fill(scatteringangle);
		  scatteringevsanglelow->Fill(scatteringangle, energy);
		  
		}
	      
	      else if(energy > 10 && energy < 1000)
		{
		scatterhistomedall->Fill(scatteringangle);
		scatteringevsanglemed->Fill(scatteringangle, energy);
		}	      
	      else
		{
		  energyMeV = energy/1000;
		scatterhistohighall->Fill(scatteringangle);
		scatteringevsanglehigh->Fill(scatteringangle, energyMeV);
		}
		  
	      if(material == "G4_lAr")
		{
		  scatterhistoalllar->Fill(scatteringangle);
		  scatterhistoalllarshield->Fill(scatteringangle);
		  
		  if(energy < 10)
		    {
		      scatterhistolowlar->Fill(scatteringangle);
		      scatterhistolowlarshield->Fill(scatteringangle);
		      deltaevsanglelowlar->Fill(scatteringangle,deltaE);
		    }
		  
		  else if(energy > 10 && energy < 1000)
		    {
		      scatterhistomedlar->Fill(scatteringangle);
		      scatterhistomedlarshield->Fill(scatteringangle);
		      deltaevsanglemedlar->Fill(scatteringangle,deltaE);
		    }
		  else
		    {
		      scatterhistohighlar->Fill(scatteringangle);
		      scatterhistohighlarshield->Fill(scatteringangle);
		      deltaevsanglehighlar->Fill(scatteringangle,deltaE);
		    }
		  
		}//LAr
	      
	      
	      else if(material == "PMMA")
		{
		  scatterhistoallshield->Fill(scatteringangle);
		  scatterhistoalllarshield->Fill(scatteringangle);

		  if(energy < 10)
		    {
		      scatterhistolowshield->Fill(scatteringangle);
		      scatterhistolowlarshield->Fill(scatteringangle);
		      deltaevsanglelowpmma->Fill(scatteringangle,deltaE);
		    }

		  else if(energy > 10 && energy < 1000)
		    {
		      scatterhistomedshield->Fill(scatteringangle);
		      scatterhistomedlarshield->Fill(scatteringangle);
		      deltaevsanglemedpmma->Fill(scatteringangle,deltaE);
		    }
		  else
		    {
		      scatterhistohighshield->Fill(scatteringangle);
		      scatterhistohighlarshield->Fill(scatteringangle);
		      deltaevsanglehighpmma->Fill(scatteringangle,deltaE);
		    }

		}//Shield

	      
	      
	    }//Same neutron

    }//int j

  scatterhistoallall->SetDirectory(gROOT);
  scatterhistolowall->SetDirectory(gROOT);
  scatterhistomedall->SetDirectory(gROOT);
  scatterhistohighall->SetDirectory(gROOT);
  scatterhistoalllar->SetDirectory(gROOT);
  scatterhistolowlar->SetDirectory(gROOT);
  scatterhistomedlar->SetDirectory(gROOT);
  scatterhistohighlar->SetDirectory(gROOT);
  scatterhistoallshield->SetDirectory(gROOT);
  scatterhistolowshield->SetDirectory(gROOT);
  scatterhistomedshield->SetDirectory(gROOT);
  scatterhistohighshield->SetDirectory(gROOT);
  scatterhistoalllarshield->SetDirectory(gROOT);
  scatterhistolowlarshield->SetDirectory(gROOT);
  scatterhistomedlarshield->SetDirectory(gROOT);
  scatterhistohighlarshield->SetDirectory(gROOT);
  deltaevsanglelowlar->SetDirectory(gROOT);
  deltaevsanglemedlar->SetDirectory(gROOT);
  deltaevsanglehighlar->SetDirectory(gROOT);
  deltaevsanglelowpmma->SetDirectory(gROOT);
  deltaevsanglemedpmma->SetDirectory(gROOT);
  deltaevsanglehighpmma->SetDirectory(gROOT);
  scatteringevsanglelow->SetDirectory(gROOT);
  scatteringevsanglemed->SetDirectory(gROOT);
  scatteringevsanglehigh->SetDirectory(gROOT);

  gStyle->SetOptStat(0000);

  /*  cout << 100*scatterhistolowlar->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistomedlar->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistohighlar->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistoalllar->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistolowshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistomedshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistohighshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistoallshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistolowlarshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistomedlarshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistohighlarshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistoalllarshield->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistolowall->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistomedall->GetEntries()/scatterhistoallall->GetEntries() << endl;
  cout << 100*scatterhistohighall->GetEntries()/scatterhistoallall->GetEntries() << endl;

  cout << scatterhistoallall->GetEntries() << endl;
  */  
  
  /*  TCanvas *c1 = new TCanvas;
  scatterhistoallall->Draw();
  TCanvas *c2 = new TCanvas;
  scatterhistolowall->Draw();
  TCanvas *c3 = new TCanvas;
  scatterhistomedall->Draw();
  TCanvas *c4 = new TCanvas;
  scatterhistohighall->Draw();
  TCanvas *c5 = new TCanvas;
  scatterhistoalllar->Draw();
  TCanvas *c6 = new TCanvas;
  scatterhistolowlar->Draw();
  TCanvas *c7 = new TCanvas;
  scatterhistomedlar->Draw();
  TCanvas *c8 = new TCanvas;
  scatterhistohighlar->Draw();
  TCanvas *c9 = new TCanvas;
  scatterhistoallshield->Draw();
  TCanvas *c10 = new TCanvas;
  scatterhistolowshield->Draw();
  TCanvas *c11 = new TCanvas;
  scatterhistomedshield->Draw();
  TCanvas *c12 = new TCanvas;
  scatterhistohighshield->Draw();
  TCanvas *c13 = new TCanvas;
  scatterhistoalllarshield->Draw();
  TCanvas *c14 = new TCanvas;
  scatterhistolowlarshield->Draw();
  TCanvas *c15 = new TCanvas;
  scatterhistomedlarshield->Draw();
  TCanvas *c16 = new TCanvas;
  scatterhistohighlarshield->Draw();
  

  scatteringevsanglehigh->GetYaxis()->SetRange(1, scatteringevsanglehigh->GetNbinsY() + 1);
  
  TCanvas *c17 = new TCanvas;
  scatteringevsanglelow->Draw("colz");
  TCanvas *c18 = new TCanvas;
  scatteringevsanglemed->Draw("colz");
  TCanvas *c19 = new TCanvas;
  scatteringevsanglehigh->Draw("colz");
  */

  deltaevsanglelowlar->SetMarkerColor(kTeal+3);
  deltaevsanglemedlar->SetMarkerColor(kTeal+3);
  deltaevsanglehighlar->SetMarkerColor(kTeal+3);
  deltaevsanglelowpmma->SetMarkerColor(kGray+2);
  deltaevsanglemedpmma->SetMarkerColor(kGray+2);
  deltaevsanglehighpmma->SetMarkerColor(kGray+2);


  //Use THStack to guarantee same range on colz
  //THStack *deltaelowstack = new THStack;
  //deltaelowstack->Add(deltaevsanglelowlar);
  //deltaelowstack->Add(deltaevsanglelowpmma);
  //THStack *deltaemedstack = new THStack;
  //THStack *deltaehighstack = new THStack;
  
  TCanvas *c20 = new TCanvas;
  //deltaevsanglelowlar->Draw("colz");
  deltaevsanglelowpmma->Draw("colz same");
  TLegend *deltaevsanglelowlegend = new TLegend(0.7,0.77,0.9,0.9);
  deltaevsanglelowlegend->AddEntry(deltaevsanglelowlar,"LAr");
  deltaevsanglelowlegend->AddEntry(deltaevsanglelowpmma,"PMMA");
  //deltaevsanglelowlegend->Draw("same");
  TCanvas *c21 = new TCanvas;
  //deltaevsanglemedlar->Draw("colz");
  deltaevsanglemedpmma->Draw("colz same");
  TLegend *deltaevsanglemedlegend = new TLegend(0.7,0.77,0.9,0.9);
  deltaevsanglemedlegend->AddEntry(deltaevsanglemedlar,"LAr");
  deltaevsanglemedlegend->AddEntry(deltaevsanglemedpmma,"PMMA");
  //deltaevsanglemedlegend->Draw("same");
  TCanvas *c22 = new TCanvas;
  //deltaevsanglehighlar->Draw("colz");
  deltaevsanglehighpmma->Draw("colz same");
  TLegend *deltaevsanglehighlegend = new TLegend(0.7,0.77,0.9,0.9);
  deltaevsanglehighlegend->AddEntry(deltaevsanglehighlar,"LAr");
  deltaevsanglehighlegend->AddEntry(deltaevsanglehighpmma,"PMMA");
  //deltaevsanglehighlegend->Draw("same");
  

  

}//EOF

  /*Function used to validate the scattering angle calculation
  double x[3];
  double y[3];
  double z[3];

  double AB;
  double BC;
  double AC;

  double internalangle = 0;
  double scatteringangle = 0;

  for (int i = 0; i < 3; i++)
    {
      x[i] = rng->Rndm()*10;
      y[i] = rng->Rndm()*10;
      z[i] = rng->Rndm()*10;

      cout << setprecision(3) << x[i] << "  " << y[i] << "  " << z[i] << endl;

    }

  AB = TMath::Sqrt(((x[1]-x[0])*(x[1]-x[0]))+((y[1]-y[0])*(y[1]-y[0]))+((z[1]-z[0])*(z[1]-z[0])));
  BC = TMath::Sqrt(((x[1]-x[2])*(x[1]-x[2]))+((y[1]-y[2])*(y[1]-y[2]))+((z[1]-z[2])*(z[1]-z[2])));
  AC = TMath::Sqrt(((x[0]-x[2])*(x[0]-x[2]))+((y[0]-y[2])*(y[0]-y[2]))+((z[0]-z[2])*(z[0]-z[2])));

  cout << endl << setprecision(3) << AB << "  " << BC << "  " << AC << endl;

  internalangle = TMath::ACos((AB*AB + BC*BC - AC*AC)/(2*AB*BC))*TMath::RadToDeg();//Law of cosines, solved for an angle



    scatteringangle = 180 - internalangle;

    cout << endl << setprecision(3) << internalangle << endl;
  */
