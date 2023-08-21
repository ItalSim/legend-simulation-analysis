void elossinshieldplots(void)
{

  //A C++/ROOT script to compile step-level information about individual neutron tracks into a set of meaningful plots to extract data about the neutron kinematics
  //Along with energy loss in the shield, special care is taken to track each individual neutron, to determine its behaviour around the shield
  //This includes information about how many times the neutron enters/exits the shield, and whether it passed through or was reflected
  //All of this information is a bit difficult to squeeze into an intuitive format, but we'll do our best!


  //Open file and tree
  TFile input("NeutronKinematics.root","READ");
  TTree *fTree = (TTree*)input.Get("Steps");
  int entries = fTree->GetEntries();

  //Open all of the potential output histograms
  //These really don't take any time or memory, so it's fine to leave them uncommented when not in use
  //Read the title of each histogram to determine what it does
  
  //OUTPUT1
  TH1D *sameneutroncrossingshisto = new TH1D("sameneutroncrossingshisto","Number of times neutron entered/exited neutron shield;Number of boundary crossings",20,0,20);

  //OUTPUT2
  TH1D *energyoutsideinwardallhisto = new TH1D("energyoutsideinwardallhisto","Neutron kinetic energy (outside of shield, aimed inwards);Energy (keV); # Neutrons",100,0,10000);
  TH1D *energyoutsideoutwardallhisto = new TH1D("energyoutsideoutwardallhisto","Neutron kinetic energy (outside of shield, aimed outwards);Energy (keV); # Neutrons",100,0,10000);
  TH1D *energyinsideinwardallhisto = new TH1D("energyinsideinwardallhisto","Neutron kinetic energy (inside of shield, aimed inwards);Energy (keV); # Neutrons",100,0,10000);
  TH1D *energyinsideoutwardallhisto = new TH1D("energyinsideoutwardallhisto","Neutron kinetic energy (inside of shield, aimed outwards);Energy (keV); # Neutrons",100,0,10000);

  //OUTPUT3
  TH1D *energyfirstcontactoutsidehisto = new TH1D("energyfirstcontactoutsidehisto","Kinetic energy of neutrons when first touching shield (outside);Energy (keV); # Neutrons",100,0,10000);
  TH1D *energyfirstcontactinsidehisto = new TH1D("energyfirstcontactinsidehisto","Kinetic energy of neutrons when first touching shield (inside);Energy (keV); # Neutrons",100,0,10000);

  //OUTPUT4
  TH1D *energylossreflectedoutsidehisto = new TH1D("energyreflectedoutsidehisto","Fraction of energy lost in the shield (outside neutrons, reflected);Energy lost (%); # Neutrons",20,0,100);
  TH1D *energylosstransmittedoutsidehisto = new TH1D("energytransmittedoutsidehisto","Fraction of energy lost in the shield (outside neutrons, through-going);Energy lost (%); # Neutrons",20,0,100);
  TH1D *energylossreflectedinsidehisto = new TH1D("energyreflectedinsidehisto","Fraction of energy lost in the shield (inside neutrons, reflected);Energy lost (%); # Neutrons",20,0,100);
  TH1D *energylosstransmittedinsidehisto = new TH1D("energytransmittedinsidehisto","Fraction of energy lost in the shield (inside neutrons, through-going);Energy lost (%); # Neutrons",20,0,100);

  //OUTPUT5
  TH2D *energylossvseoutsidehisto = new TH2D("energylossvseoutsidehisto","Fraction of energy lost vs. kinetic energy (outside neutrons);Energy lost (%);Kinetic energy upon entering (keV)",100,0,100,20,0,500);
  TH2D *energylossvseinsidehisto = new TH2D("energylossvseinsidehisto","Fraction of energy lost vs. kinetic energy (inside neutrons);Energy lost (%);Kinetic energy upon entering (keV)",100,0,100,20,0,500);

  //OUTPUT6
  TH1D *consecutivestepsreflectedhisto = new TH1D("consecutivestepsreflectedhisto","Number of uninterrupted steps in shield; #Steps; #Neutrons",500,0,500);
  TH1D *consecutivestepstransmittedhisto = new TH1D("consecutivestepstransmittedhisto","Number of uninterrupted steps in shield; #Steps; #Neutrons",500,0,500);

  
  //Initialize the variables to read from the tree
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


  //There's some weird trick with the multithread ntuple and the way it stores 'strings'
  //For some reason it stores all strings as character arrays
  //Obviously, C++ throws a tantrum about this, so we have to convert each entry to a string
  //each time we read it if we want to operate on it like a string
  string materialone;
  string materialtwo;
  string process;


  
  //Analysis variables

  //These first few variables change all the time
  int firsttracknumber = 0;
  int secondtracknumber = 0;//These two are used to determine when the processing loop switches from one neutron to another
  double radius = 0;//Used to determine if a neutron enters the shield from the inside or the outside

  //The following are persistent variables for all-track analysis and should never reset
  //------------------------------------------------------------------------------------
  int OLtoPcrossings = 0;
  int ILtoPcrossings = 0;
  int PtoOLcrossings = 0;
  int PtoILcrossings = 0;//Used to tally the type of boundary crossing

  int outsidereflected = 0;
  int insidereflected = 0;
  int outsidetransmitted = 0;
  int insidetransmitted = 0;//Used to tally if the neutron went through the shield (transmission) or not (reflection)
  
  int lasttracknumber = 0;
  int uniqueneutrons = 0;
  int sameneutroncrossings = 0;//Used to check parameters about neutrons
  //------------------------------------------------------------------------------------


  //These next variables should change any time a neutron enters the shield
  
  int enteredfromoutside = 0;//Remembers whether the neutron enters the shield from the inside or outside surface
  //Sadly I can't use a bool, because some neutrons are created INSIDE the shield, and would mess up the numbers
  //We'll define a value of 0 as 'hasnt entered', a value of one as 'entered from inside', and a value of 2 as 'entered from outside'
  double energyuponentering = 0;
  double deltaEoverE = 0;
  int tracknumberuponentering = 0;
  int stepnumberuponentering = 0;


  //Prime the comparators that need priming
  //Since the master loop needs information about the previous step, we cannot start on step 0 inside the loop
  fTree->GetEntry(0);
  lasttracknumber = tracknumber;

  
  for (int j = 1; j<(entries-1); j++)
    {//Master loop

      fTree->GetEntry(j);
      secondtracknumber = tracknumber;

      //Make sure the vertices before/after belong to the same neutron
      //This doesn't account for the vanishingly rare case where the last neutron
      //in one event shares a tracknumber with the first neutron in the next event

      fTree->GetEntry(j-1);	  
      firsttracknumber = tracknumber;
      materialone = materialc;

      if(firsttracknumber == secondtracknumber)
	{//The current step is the same neutron as the last step

	  fTree->GetEntry(j);
	  process = processc;
	  materialtwo = materialc;

	  	  
	  if(materialone=="G4_lAr"&&materialtwo=="PMMA")
	    {//Entering the shield, LAr to PMMA, should affect most variables

	      radius = sqrt(X*X+Y*Y);

	      if(lasttracknumber!=tracknumber)
		{//Different neutron than previously recorded, used by histograms which tally track-level information for the whole run
		  lasttracknumber = tracknumber;
		  uniqueneutrons++;
		  sameneutroncrossingshisto->Fill(sameneutroncrossings);
		  sameneutroncrossings = 1;
		  enteredfromoutside = 0;
		}

	      else
		sameneutroncrossings++;

	      
	      if(radius < 2050 && Z > -1630 && Z < 470)//Neutron is in the volume enclosed by the shield, touching an inner wall
		{
		ILtoPcrossings++;
		enteredfromoutside = 1;
		energyinsideoutwardallhisto->Fill(energy);
		if(sameneutroncrossings == 1)//First time touching shield
		  energyfirstcontactinsidehisto->Fill(energy);
		energyuponentering = energy;
		tracknumberuponentering = tracknumber;
		stepnumberuponentering = stepnumber;

		}
	      else//Neutron is outside of the volume enclosed by the shield
		{
		OLtoPcrossings++;
		enteredfromoutside = 2;
		energyoutsideinwardallhisto->Fill(energy);
		if(sameneutroncrossings == 1)//First time touching shield
		  energyfirstcontactoutsidehisto->Fill(energy);
		energyuponentering = energy;
		tracknumberuponentering = tracknumber;
		stepnumberuponentering = stepnumber;

		}
	    }//LAr to PMMA

	  
	  
	  else if(materialone=="PMMA"&&materialtwo=="G4_lAr")
	    {//Exiting the shield, PMMA to LAr

	      radius = sqrt(X*X+Y*Y);
	      deltaEoverE = 100*(energyuponentering-energy)/energyuponentering;


	      if(lasttracknumber!=tracknumber)//Again, persistent tallying comes first
		{
		  lasttracknumber = tracknumber;
		  uniqueneutrons++;
		  sameneutroncrossingshisto->Fill(sameneutroncrossings);
		  sameneutroncrossings = 1;
		  enteredfromoutside = 0;
		}
	      else
		sameneutroncrossings++;

	      
	      if(radius < 2050 && Z > -1630 && Z < 470)
		{//Neutron goes inside the volume enclosed by the shield
		  PtoILcrossings++;
		  energyinsideinwardallhisto->Fill(energy);

		  if(enteredfromoutside == 1)//Started inside the shield as well
		    {
		      insidereflected++;
		      if(tracknumberuponentering==tracknumber)
			{
			  energylossreflectedinsidehisto->Fill(deltaEoverE);
			  energylossvseinsidehisto->Fill(deltaEoverE,energy);
			  consecutivestepsreflectedhisto->Fill(stepnumber - stepnumberuponentering);
			}
		    }
		  else if(enteredfromoutside == 2)//Started outside the shield
		    {
		      outsidetransmitted++;
		      if(tracknumberuponentering==tracknumber)
			{
			  energylosstransmittedoutsidehisto->Fill(deltaEoverE);
			  energylossvseoutsidehisto->Fill(deltaEoverE,energy);
			  consecutivestepstransmittedhisto->Fill(stepnumber - stepnumberuponentering);
			}
		    }
		  energyuponentering = 0;
		  tracknumberuponentering = 0;
		  stepnumberuponentering = 0;

		}//Enclosed in the shield
	      

              else //Neutron goes outside away from the shield
		{
		  PtoOLcrossings++;
		  energyoutsideoutwardallhisto->Fill(energy);
		  if(enteredfromoutside == 1)//Started inside the shield
		    {
		      insidetransmitted++;
		      if(tracknumberuponentering==tracknumber)
			{
			  energylosstransmittedinsidehisto->Fill(energyuponentering);
			  energylossvseinsidehisto->Fill(energyuponentering,energy);
			  consecutivestepstransmittedhisto->Fill(stepnumber - stepnumberuponentering);
			}
		    }
                  else if(enteredfromoutside == 2)//Started outside the shield as well
		    {
		      outsidereflected++;
		      if(tracknumberuponentering==tracknumber)
                        {
			  energylossreflectedoutsidehisto->Fill(energyuponentering);
			  energylossvseoutsidehisto->Fill(energyuponentering,energy);
			  consecutivestepsreflectedhisto->Fill(stepnumber - stepnumberuponentering);
			}
		    }
		  energyuponentering = 0;
		  tracknumberuponentering = 0;
		  stepnumberuponentering = 0;
		}//Outside the shield

	      
	      enteredfromoutside = 0;

	    }//PMMA to LAr

	  
	  
	}//Same neutron

      //This analysis script does NOT do step-level analysis/tallying - that's taken care of in another script in this module
      //So, if the neutron is different between steps, there's nothing to be done
      
      
    }//Master loop


  
  //Print out some relevant numbers from the persistent tallies, for use in tables and such
  cout << "Total unique neutrons making at least one boundary crossing    " << uniqueneutrons << endl;
  cout << "Total neutron boundary crossings (PMMA or LAr to LAr or PMMA)   " << OLtoPcrossings+ILtoPcrossings+PtoOLcrossings+PtoILcrossings << endl;
  cout << "Total neutron boundary crossings (LAr to PMMA only)   " << OLtoPcrossings+ILtoPcrossings  << endl;
  cout << "Total neutron boundary crossings (PMMA to LAr only)   " << PtoOLcrossings+PtoILcrossings  << endl;
  cout << "Neutron boundary crossings (outer LAr to PMMA)   " <<OLtoPcrossings  << endl;
  cout << "Neutron boundary crossings (inner LAr to PMMA)   " <<ILtoPcrossings  << endl;
  cout << "Neutron boundary crossings (PMMA to outer LAr)   " <<PtoOLcrossings  << endl;
  cout << "Neutron boundary crossings (PMMA to inner LAr)   " <<PtoILcrossings  << endl;
  cout << "Outer neutrons reflected:   " << outsidereflected << endl;
  cout << "Outer neutrons transmitted:   " << outsidetransmitted << endl;
  cout << "Inner neutrons reflected:   " << insidereflected << endl;
  cout << "Inner neutrons transmitted:   " << insidetransmitted << endl;

  gStyle->SetOptStat(0000);

  sameneutroncrossingshisto->SetDirectory(gROOT);
  energyinsideinwardallhisto->SetDirectory(gROOT);
  energyinsideoutwardallhisto->SetDirectory(gROOT);
  energyoutsideinwardallhisto->SetDirectory(gROOT);
  energyoutsideoutwardallhisto->SetDirectory(gROOT);
  energyfirstcontactinsidehisto->SetDirectory(gROOT);
  energyfirstcontactoutsidehisto->SetDirectory(gROOT);
  energylossreflectedoutsidehisto->SetDirectory(gROOT);
  energylosstransmittedoutsidehisto->SetDirectory(gROOT);
  energylossreflectedinsidehisto->SetDirectory(gROOT);
  energylosstransmittedinsidehisto->SetDirectory(gROOT);
  energylossvseoutsidehisto->SetDirectory(gROOT);
  energylossvseinsidehisto->SetDirectory(gROOT);
  consecutivestepsreflectedhisto->SetDirectory(gROOT);
  consecutivestepstransmittedhisto->SetDirectory(gROOT);

  //OUTPUT1
  /*
  TCanvas *sameneutroncrossingscanvas = new TCanvas;
  sameneutroncrossingscanvas->SetLogy();
  sameneutroncrossingshisto->Draw();
  */

  //OUTPUT2
  /*
  TCanvas *energyinsideinwardallcanvas = new TCanvas;
  energyinsideinwardallcanvas->SetLogy();
  energyinsideinwardallhisto->Draw();
  TCanvas *energyinsideoutwardallcanvas = new TCanvas;
  energyinsideoutwardallcanvas->SetLogy();
  energyinsideoutwardallhisto->Draw();
  TCanvas *energyoutsideinwardallcanvas = new TCanvas;
  energyoutsideinwardallcanvas->SetLogy();
  energyoutsideinwardallhisto->Draw();
  TCanvas *energyoutsideoutwardallcanvas = new TCanvas;
  energyoutsideoutwardallcanvas->SetLogy();
  energyoutsideoutwardallhisto->Draw();
  */

  //OUTPUT3
  /*
  TCanvas *energyfirstcontactinsidecanvas = new TCanvas;
  energyfirstcontactinsidecanvas->SetLogy();
  energyfirstcontactinsidehisto->Draw();
  TCanvas *energyfirstcontactoutsidecanvas = new TCanvas;
  energyfirstcontactoutsidecanvas->SetLogy();
  energyfirstcontactoutsidehisto->Draw();
  */
  
  //OUTPUT4
  /*
  TCanvas *energylossreflectedoutsidecanvas = new TCanvas;
  energylossreflectedoutsidehisto->Draw();
  TCanvas *energylosstransmittedoutsidecanvas = new TCanvas;
  energylosstransmittedoutsidehisto->Draw();
  TCanvas *energylossreflectedinsidecanvas = new TCanvas;
  energylossreflectedinsidehisto->Draw();
  TCanvas *energylosstransmittedinsidecanvas = new TCanvas;
  energylosstransmittedinsidehisto->Draw();
  */

  //OUTPUT5
  /*
  energylossvseoutsidehisto->SetMarkerStyle(21);
  energylossvseinsidehisto->SetMarkerStyle(21);
  energylossvseoutsidehisto->SetMarkerSize(.2);
  energylossvseinsidehisto->SetMarkerSize(.2);  
  TCanvas *energylossvseoutsidecanvas = new TCanvas;
  energylossvseoutsidehisto->Draw();  
  TCanvas *energylossvseinsidecanvas = new TCanvas;
  energylossvseinsidehisto->Draw();  
  */

  //OUTPUT6
  /*
  consecutivestepsreflectedhisto->SetLineColor(kBlue-4);
  consecutivestepstransmittedhisto->SetLineColor(kRed-3);
  consecutivestepsreflectedhisto->SetLineWidth(3);
  consecutivestepstransmittedhisto->SetLineWidth(3);

  TCanvas *consecutivestepscanvas = new TCanvas;
  consecutivestepsreflectedhisto->Draw();
  consecutivestepstransmittedhisto->Draw("same");

  TLegend *consecutivestepslegend = new TLegend(0.7,0.77,0.9,0.9);
  consecutivestepslegend->AddEntry(consecutivestepsreflectedhisto,"Reflected");
  consecutivestepslegend->AddEntry(consecutivestepstransmittedhisto,"Transmitted");
  consecutivestepslegend->Draw("same");
  */
  
}//EOF
