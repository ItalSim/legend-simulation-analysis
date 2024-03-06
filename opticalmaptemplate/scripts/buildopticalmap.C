void buildopticalmap(int startnumber)
{
  //There are two general philosophies which may be followed when trying to construct an optical map for G4 simulations:
  //1) Save all of the simulation data in an unstructured output, then format the output into an optical map
  //2) Do some runtime processing inside of the simulation in order to have a structured output

  //The second involves fewer steps, in theory, and also can save some data storage space
  //However, this could have unpredictable effects on runtime/job stability
  //Also, the first is easier to implement and more flexible
  //This is an attempt at the first

  //Given an unstructured set of ordered (x,y,z) and (x',y',z') coordinate pairs, construct an optical map
  //In the final configuration, the optical map should take advantage of the 12-fold symmetry within the simulation,
  //reducing the structured output filesize by at least this factor. This can be done in a unique way, by using a
  //translation in radial coordinates to effectively move all coordinate pairs so that they overlap in one sector
  //This explanation may be a bit lacking, but the implementation should be more easy to understand

  //Find the .root file(s) to be processed

  const char* allfilesbeginningwith = Form("opmaprawdata%i.", startnumber);
  const char* indirectory = "/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap/4-innermap-Dec2023/output/raw/";

  char* directory = gSystem->ExpandPathName(indirectory);
  void* directoryP = gSystem->OpenDirectory(directory);

  const char* checkedname;
  const char* acceptedname[1000];
  int wan = 0;
  TString checkedstring;

  while((checkedname = (char*)gSystem->GetDirEntry(directoryP)))
    {
      checkedstring = checkedname;
      cout << "Checked string: " << checkedstring << endl;
      if(checkedstring.BeginsWith(allfilesbeginningwith))
	{
	  acceptedname[wan] = checkedname;
	  wan++;
	}
    }


  //Declare persistent variables, such as the output histogram(s)

  //For reference: current 'standard shield is a 12-sided polygon
  //The reentrance tube has a radius of 95 cm
  //The largest circle which can be drawn in the LAr without touching the shield (XY plane) is 200 cm radius
  //The smallest circle which contains the entire inner LAr volume (some overlap with shield) is 207.055 cm radius
  //The largest circle which is completely within the shielding material (inscribed) is 210 cm radius
  //The smallest circle which completely contains (circumscribes) the shield is 217.408 cm in radius
  //The largest circle which can be inscribed in the cryostat is 269 cm in radius
  //The limits of the histograms are subject to change in the future

  //More reference calculations
  //A 30 degree chord across a 269 cm circle is 139.2 cm in length - rounding up, 140 cm
  //Maximum straight-line distance from shield to cryostat is of course 269-210 = 59 cm
  //z height, as discussed in detail later, is 320 cm (or 300 cm for the active area)
  //Therefore, we have to fill a 140x320x59 cm rectangular prism, of total volume 2.6432 m^3
  //If each voxel is a cubic centimeter, this results in 2.64 million voxels

  //Define voxel size and boundaries

  //This is a personal preference - but I prefer to iterate with x as the left-right dimension and y as the up-down dimension,
  //if viewing the geometry along the z axis looking downwards
  //Given the above dimensions, our bounding box is therefore:

  //For outer map (small cryo):
  //-70 < x < 70
  //210 < y < 269
  //-160 < z < 160 (after transforms)
  //Expected volume = 2643200 cm^3

  //For outer map (new cryo ) Oct 2023):
  //-85 < x < 85
  //210 < y < 325
  //-160 < z < 160
  //Expected volume = 6256000 cm^3
  
  //For inner map:
  //Chord length: 0.491756
  //For the inner map, we cannot neglect the curvature of the reentrance tube, so y can be < 95 if it's along the outer sides
  //I calculated the minimum y value between 91 and 92, so of course we round down to 91 for safety
  //-54 < x < 54
  //91 < y < 200*
  //-151 < z < 151**
  //Expected volume = 3555144 cm^3

  //*In theory, the radius can be out to 207 cm for the inner photons, but the x and y individually cannot exceed 200 cm
  //**For safety/to catch all the photons at exactly z=+-150


  //Important note: the expected volume may not be equal to the final volume, depending on the voxel size. There may be extra
  

  const int voxelsize = 2;//In cm on each side

  //All dimensions in cm
  //Inner map settings
  const int xmin = -54;
  const int xmax = -1*xmin;
  int xbins = floor((xmax-xmin)/voxelsize);
  const int ymin = 91;
  const int ymax = 200;
  int ybins = floor((ymax-ymin)/voxelsize);

  
  //Currently, the shield is 3.2 m high (3m panels + 0.1 m thick top and bottom pieces)
  //The shield has 0.58 m offset from the absolute world center
  //We define the radius from the world center, and the z value from the shield center
  //This is because the shield is guaranteed to be centered radially, but not necessarily in z
  const int zoffset= 58;
  //const int totalz = 320; For outer map
  const int totalz = 302; //For inner map
  int zbins     = floor(totalz/voxelsize);
  const int zmin   = -1*totalz/2;
  const int zmax   =    totalz/2;

  //Correction for when the map dimensions and the voxel dimensions don't line up perfectly
  if((xmax - xmin) % voxelsize)
    xbins++;
  if((ymax - ymin) % voxelsize)
    ybins++;
  if((zmax - zmin) % voxelsize)
    zbins++;
  
  cout << "xmin : " << xmin << "xmax : " << xmax << "ymin : " << ymin << "ymax : " << ymax << "zmin : " << zmin << "zmax : " << zmax<< endl; 
  const int matrixsize = xbins*ybins*zbins;
  cout << xbins << "  "  << ybins << "  " << zbins << endl;
  cout << endl << "Matrix size (in voxels): " <<  matrixsize << endl << endl;
  
  //Define the output file
  //Each voxel will have the following values:
  //index (int), x' (float array), y' (float array), z' (float array)
  //The number of hits in the voxel will be the size of the array, and can be used to
  //calculate the probability implicitly, in the next processing script

  int index = 0;
  vector<Float_t> *Xprime = 0;
  vector<Float_t> *Yprime = 0;
  vector<Float_t> *Zprime = 0;
  
  TFile *output = new TFile(Form("/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap/4-innermap-Dec2023/output/built/builtmap%i.root",startnumber),"recreate");	
  TTree *opmaptree = new TTree("opmap","opmap");
  opmaptree->Branch("index",&index);
  opmaptree->Branch("Xprime",&Xprime);
  opmaptree->Branch("Yprime",&Yprime);
  opmaptree->Branch("Zprime",&Zprime);

  //File-by-file processing  
  for(int iwan = 0;iwan < wan;iwan++)
    {
      cout << acceptedname[iwan] << endl;

      
         //Initialize variables
      //I/O
      TFile *input = new TFile(Form("%s%s",indirectory,acceptedname[iwan]),"read");
      input->cd();

      TTree *OpMapData = (TTree*)input->Get("OpticalMapData");

      //The original idea was to make a vector of vectors of 3-vectors
      //But we can instead save a vector of vector of ints for the map
      //This reduces the RAM requirements by a factor of ~eight,
      //but doubles the runtime. This could be reverted if necessary.

      //The idea is to save the entry# in the tree for map building,
      //rather than save the prime position values at that entry#
      vector<vector<int>> map;
      
      map.resize(matrixsize, vector<int>(0));

      
      //setup to work with input tree

      int entries = OpMapData->GetEntries();

      float X1 = 0;
      float Y1 = 0;
      float Z1 = 0;
      float X2 = 0;
      float Y2 = 0;
      float Z2 = 0;

      OpMapData->SetBranchAddress("X1",&X1);
      OpMapData->SetBranchAddress("Y1",&Y1);
      OpMapData->SetBranchAddress("Z1",&Z1);
      OpMapData->SetBranchAddress("X2",&X2);
      OpMapData->SetBranchAddress("Y2",&Y2);
      OpMapData->SetBranchAddress("Z2",&Z2);

      OpMapData->GetEntry(0);
      cout << "Entries in tree for this file: " << entries << endl << endl;

      int onepercentofentries = entries/100.;//Used for the processing log
      
      //Analysis variables
      double radius = 0;
      double phi    = 0;
      double phidecimalvalue = 0;
      double phidifference = 0;
      double xtr;//Transforms
      double ytr;
      double ztr;
      int xindex;
      int yindex;
      int zindex;
      int finalindex;
      TVector3 fillvector;

      bool verbose = false;

      if(verbose)
	{
	  entries = 10;//For testing
	  cout << "Map size: " << map.size() << endl;
	}

      
      cout << "Processing loop 1 - generating virtual map structure" << endl << endl;
            
      for(int i = 0; i < entries; i++)
	{//Main processing loop

	  OpMapData->GetEntry(i);

	  if(i%onepercentofentries == 0)
	    cout << i*1/onepercentofentries << "% complete." << endl;
	  
	  if(verbose)
	    cout << "X1: " << X1 << endl
	         << "Y1: " << Y1 << endl
		 << "Z1: " << Z1 << endl << endl;

	  //Input is in mm, but calculations/map are in cm, so divide by 10
	  radius = TMath::Sqrt(X1*X1+Y1*Y1)/10.;
	  phi = TMath::RadToDeg()*TMath::ATan2(Y1,X1);
	  ztr = zoffset + Z1/10.;

	  if(verbose)
	    cout << "PHI1: " << phi << endl;
	  
	  //Incorporate the 12-fold symmetry
	  //We've chosen the coordinates of our map such that at the center, x=0, y=positive
	  //This suggests a central value of phi of 90 degrees, assuming phi=0 lies along x
	  //With 12-fold symmetry, we only need to cover a 30 degree range, or +-15 from 90
	  //Our valid range is (75,105)
	  //It will be easier to just transform our phi values by increments of 30 degrees
	  while(phi < 75 || phi > 105)
	    {
	      if(phi < 75)
		phi += 30;
	      else
		phi -= 30;
	    }
	      if(verbose)
		cout << "PHI2: " << phi << endl;
	  
	  xtr = radius*cos(phi*TMath::DegToRad());
	  ytr = radius*sin(phi*TMath::DegToRad());
	    
	  //Always round down (is this necessary?)
	  //I don't think it's necessary, but keeping it for legacy
	  //xtr = floor(xtr);
	  //ytr = floor(ytr);
	  //ztr = floor(ztr);

	  if(verbose)
	    cout << "XTR: " << xtr << endl
	         << "YTR: " << ytr << endl
		 << "ZTR: " << ztr << endl << endl;

	  //Calculate the indice of the matrix to fill

	  xindex = floor(abs((xtr - xmin))/voxelsize);//Should be from 0 to xbins-1
	  yindex = floor(abs((ytr - ymin))/voxelsize);//Should be from 0 to ybins-1
	  zindex = floor(abs((ztr - zmin))/voxelsize);//Should be from 0 to zbins-1

	  if(verbose)
	    cout << "XIND: " << xindex << endl 
	         << "YIND: " << yindex << endl
		 << "ZIND: " << zindex << endl << endl;

	  //The index is a very important value, implicitly encoding the voxel positions
	  //and keeping each voxel distinct in future processing.
	  //For a full explanation of how this works, review some of my talks - CJ
	  finalindex = xindex + (yindex*xbins) + (zindex*xbins*ybins);

	  if(verbose)
	    cout << "INDEX: " << finalindex << endl << endl;

	  map[finalindex].push_back(i);//Save tree entry number with this index

	  //It may not be obvious, but we cannot save the map values in this processing loop
	  
	}//Main processing loop 1

      
      if(verbose)
	{
	  for(int iter1 = 0; iter1 < map.size(); iter1++)
	    {
	      if(map[iter1].size())
		cout << iter1 << "   ";
	      for(int iter2 = 0; iter2 < map[iter1].size();iter2++)
		cout << map[iter1][iter2] << "   ";
	      if(map[iter1].size())
		cout << endl;	
	    }
	}
      
      //With the map generated, we can fill each leaf of the TTree sequentially

      //At this point, consider how the map will be used
      //For now, I assume the user will have access to this code and the index calculation
      //In that case, it is better to store the index value, and allow the user to unfold
      //the positional coordinates themselves

      cout << endl << "Processing loop 2 - adding to output file" << endl << endl;

      onepercentofentries = matrixsize/100;
      
      for(int j1 = 0; j1 < matrixsize; j1++)
	{//Main processing loop 2

	  if(j1%onepercentofentries == 0)
	    cout << j1*1/onepercentofentries << "% complete." << endl;
	  
	  for(int j2 = 0; j2 < map[j1].size(); j2++)
	    {//For all shield hits in each voxel:
	      
	                if(verbose)
			  cout << "X2: " << X2 << endl
			       << "Y2: " << Y2 << endl
			       << "Z2: " << Z2 << endl << endl;
			
	      OpMapData->GetEntry(map[j1][j2]);

	      radius = TMath::Sqrt(X2*X2+Y2*Y2)/10.;
	      ztr = zoffset + Z2/10.;//Unchanged
	      
	      //Use the 12-fold symmetry
	      //This case is a bit more nuanced, because the prime coordinates can be on any side of the shield
	      //which is to say, we want to preserve the relationship between X and Xprime, etc
	      //We want to transform the prime coordinates in EXACTLY the same way as the origin coordinates
	      //So, let's once again calculate the phi transform for the origin coordinates, and store it

	      phi = TMath::RadToDeg()*TMath::ATan2(Y1,X1);//Use origin coordinates first
	      phidifference = phi;//Store this value
	      
	  while(phi < 75 || phi > 105)
	    {//Copypasted from processing loop 1, for consistency
	      if(phi < 75)
		phi += 30;
	      else
		phi -= 30;
	    }	      	      

	      phidifference -= phi;//We now have the phi transform value and can apply it to the primes

	      //Now use destination coordinates
	      phi = TMath::RadToDeg()*TMath::ATan2(Y2,X2);

	      phi -= phidifference;//This does NOT have to be between 75 and 105

	      xtr = radius*cos(phi*TMath::DegToRad());
	      ytr = radius*sin(phi*TMath::DegToRad());
	      
	      Xprime->push_back(xtr);
	      Yprime->push_back(ytr);
	      Zprime->push_back(ztr);

	      if(verbose)
		{
		  cout << "XTR2: " << xtr << endl  << "YTR2: " << ytr << endl << "ZTR2: " << ztr  << endl << endl;
		cout << radius << endl;
		cout << phi << endl;
		cout << phidifference << endl;
		cout << index << endl;
		}
	    }

	  //Fill the output with the index and the filled prime coordinates array
	  index = j1;
	  output->cd();
	  opmaptree->Fill();
	  input->cd();

	  //Cleanup
	  Xprime->clear();
	  Xprime->shrink_to_fit();
	  Yprime->clear();
	  Yprime->shrink_to_fit();
	  Zprime->clear();
	  Zprime->shrink_to_fit();

	  map[j1].clear();
	  map[j1].shrink_to_fit();
	  
	}//Main processing loop 2

      
    }//File by file processing

  //Write the output
  output->cd();
  opmaptree->Write();

  
}//EOF
