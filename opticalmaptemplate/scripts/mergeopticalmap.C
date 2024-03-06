void mergeopticalmap(int startnumber)
{

  //Once the optical map is built, the files need to be combined. Typically, for
  //ROOT files, we would simply use hadd to throw all the files together. However,
  //this does not work due to the structure of the map. In the map, each index
  //value has exactly one entry, and the leaves of that entry contain our primes.
  //When we use hadd for, say, 10 files with 100 index numbers, instead of
  //getting 100 index numbers, one entry each, 10x as many leaves, we get
  //100 index numbers, 10 entries each, and the leaves are kept separate.

  //To illustrate:
  /*
    This is what we want:
                         [1][2]     [1][2]     [1][2]
                         [L][L]     [L][L]     [L][L]
                         [L][L]  +  [L][L]  =  [L][L]
                         [L][L]                [L][L]
                                               [L][L]
			         	       [L][L]


    But with hadd we get:
                         [1][2]     [1][2]     [1][2][1][2]
                         [L][L]     [L][L]     [L][L][L][L]
                         [L][L]  +  [L][L]  =  [L][L][L][L]
                         [L][L]                [L][L]
  */                                         
  //So now we have to merge the leaves ourselves.

  //While we're at it, this is a good time to calculate and add the probability
  //value for each voxel/index. We do this using the voxel volume compared to the
  //total volume of the region covered by the map, and the number of photons
  //sampled, to determine on average how many photons should have been sampled
  //in each voxel. Then, the probability is just the number of leaves in the
  //voxel divided by this expectation value of number of samples per voxel.
  //Does that make sense? If not, too bad, attend more of my talks.

  //Unfortunately, if we wanted to do this iteratively, probably our only
  //options in C++ would be to either load all the input files into virtual
  //memory, or to open/close them every time we filled one index, which would
  //be millions of calls to the disk. It's less of a headache to explicitly
  //open all the files at once, because I don't believe they'll all be loaded
  //into virtual memory, and we can flush the buffer after each index is filled.


  //This is a rather complicated program, so, sorry about that.

  
    //Find the .root file(s) to be processed

  const char* allfilesbeginningwith = "builtmap";
  const char* indirectory = "/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap/4-innermap-Dec2023/output/built/";

  char* directory = gSystem->ExpandPathName(indirectory);
  void* directoryP = gSystem->OpenDirectory(directory);

  const char* checkedname;
  const char* acceptedname[1000];
  int wan = 0;
  TString checkedstring;

  while((checkedname = (char*)gSystem->GetDirEntry(directoryP)))
    {
      checkedstring = checkedname;
      //cout << "Checked string: " << checkedstring << endl;
      if(checkedstring.BeginsWith(allfilesbeginningwith))
        {
          acceptedname[wan] = checkedname;
          wan++;
        }
    }

  
  //Declare variables etc
  
  TFile *input;// = new TFile("/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap1/output/built/builtmap1.root","read");
  TTree *inputtree; //= (TTree*)input1->Get("opmap");  

  vector<float> *Xprime1 = 0;
  vector<float> *Yprime1 = 0;
  vector<float> *Zprime1 = 0;
  
  //Input the map parameters yourself, because they're impossible to determine with the info given
  //Change some ints to doubles so that they act right in the floating point calculations
  
  //These values are for the inner optical map
  
  int voxelsize = 2;//In cm on each side

  //All dimensions in cm
  int xmin = -54;
  int xmax = -1*xmin;
  int xbins = floor((xmax-xmin)/voxelsize);
  int ymin = 91;
  int ymax = 200;
  int ybins = floor((ymax-ymin)/voxelsize);
    
  //The shield has -0.58 m offset from the absolute world center
  //We define the radius from the world center, and the z value from the shield center
  //This is because the shield is guaranteed to be centered radially, but not necessarily in z
  double zoffset= 58;
  int totalz = 302;
  int zbins     = floor(totalz/voxelsize);
  int zmin   = -1*totalz/2;
  int zmax   =    totalz/2;
    
  //Correction for when the map dimensions and the voxel dimensions don't line up perfectly
  if((xmax - xmin) % voxelsize)
    xbins++;
  if((ymax - ymin) % voxelsize)
    ybins++;
  if((zmax - zmin) % voxelsize)
    zbins++;

  int matrixsize = xbins*ybins*zbins;
    
  cout << "xbins: " << xbins  << " ybins: " << ybins  << " zbins: " << zbins << endl << "Matrix size: " << matrixsize << " voxels" << endl;

  
  //Calculate the expected number of initial photons per voxel (not the number which reach the shield)
  
  double numberofinitialphotons = 1000000000;//10^9

  //There is an intentional design "flaw" in the inner map simulation
  //The LAr between the bottom of the reentrance tube and the shield isn't sampled
  //This accounts for ~5% of the total AAr volume inside the shield, but including
  //it in the photon sampling  would increase the size of the map by a factor of ~2
  //So, we have a dodecagonal prism, with a cylinder missing from the center
        
  double cylindervolume = TMath::Pi() * 95 * 95 * 300;//Reentrance tube + unsampled LAr volume

  //Easier to just calculate elsewhere than to add the formula here
  double shieldvolume = 128616 * 300;//Total volume enclosed by the shield

  double samplingvolume = (shieldvolume - cylindervolume) / 12;//divide by 12 for the 12-fold symmetry

  cout <<"Sampling volume size (m^3): " << samplingvolume/1000000. << endl << endl;

  double voxelvolume = (double)voxelsize*(double)voxelsize*(double)voxelsize;

  double voxelfraction = voxelvolume / samplingvolume;

  double photonspervoxel = numberofinitialphotons * voxelfraction;

  cout <<"Expected INITIAL photons per voxel: " << photonspervoxel << endl << endl;


  
  TFile *output = new TFile(Form("/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap/4-innermap-Dec2023/output/merged/mergedmap%i.root",startnumber),"recreate");
  TTree *mergedopmap = new TTree("opmap","opmap");

  //Output variables
  int index = 0;
  vector<vector<float>> Xprimearray;
  vector<vector<float>> Yprimearray;
  vector<vector<float>> Zprimearray;
  double probability = 0;

  vector<float> *Xprime = 0;
  vector<float> *Yprime = 0;
  vector<float> *Zprime = 0;
    
    
  mergedopmap->Branch("index",&index);
  mergedopmap->Branch("Xprime",&Xprime);
  mergedopmap->Branch("Yprime",&Yprime);
  mergedopmap->Branch("Zprime",&Zprime);
  mergedopmap->Branch("probability",&probability);


  int numbertoprocess = 4485;//Should be number of matrices divided by the number of jobs submitted, rounded up
  //For inner map, number of voxels is 448470
    
  //For bit-by bit processing, set the range of i

  int starti = numbertoprocess*startnumber;
  int endi = starti+numbertoprocess;

  //Makes sure the index doesn't go out of range
    
  if(starti > matrixsize || starti==matrixsize)//Sent too many jobs, silly!
    {
      cout << "Fatal error: the index which this routine would begin processing is out of range." << endl;
      return;
    }
  if(endi > matrixsize)//The last file, most likely. Can be processed, but not over the whole range
    {
      cout << endl << "File only partially in range - Adjusting range" << endl;
      numbertoprocess = numbertoprocess - endi + matrixsize;
      cout << "New number to process: " << numbertoprocess << endl;
      endi = matrixsize;
    }

    
  //File loop start
  for(int iwan = 0;iwan < wan;iwan++)
    {
      cout << indirectory << acceptedname[iwan] << endl;	  

      //Declare inputtree and input, open files
      input = new TFile(Form("%s%s",indirectory,acceptedname[iwan]),"read");
      input->cd();
      inputtree = (TTree*)input->Get("opmap");	  
      inputtree->SetBranchAddress("Xprime",&Xprime1);
      inputtree->SetBranchAddress("Yprime",&Yprime1);
      inputtree->SetBranchAddress("Zprime",&Zprime1);

      Xprimearray.resize(numbertoprocess);
      Yprimearray.resize(numbertoprocess);
      Zprimearray.resize(numbertoprocess);

      //First main processing loop, to merge the files
      //Note that the vector needs to start from 0
      //but the primes can be loaded in from anywhere
      index = 0;

      for(int i = starti; i < endi; i++)
	{
	  inputtree->GetEntry(i);
	  Xprimearray[index].insert(Xprimearray[index].end(), Xprime1->begin(), Xprime1->end());
	  Yprimearray[index].insert(Yprimearray[index].end(), Yprime1->begin(), Yprime1->end());
	  Zprimearray[index].insert(Zprimearray[index].end(), Zprime1->begin(), Zprime1->end());
	  index++;

	}//First main processing loop end
		  
    }//File loop end

  
  //Second processing loop - only once
  //With the merging complete, we can calculate the probability of the voxel's photons
  //reaching the shield and fill the output
  
  for(int i = 0; i < numbertoprocess; i++)
    {
      index = starti+i;
      Xprime->insert(Xprime->end(),Xprimearray[i].begin(),Xprimearray[i].end());
      Yprime->insert(Yprime->end(),Yprimearray[i].begin(),Yprimearray[i].end());
      Zprime->insert(Zprime->end(),Zprimearray[i].begin(),Zprimearray[i].end());
      probability = (double)Zprime->size()/photonspervoxel;
      //cout << index << "   " <<  probability*100 << endl;
	
      mergedopmap->Fill();

      Xprime->clear();
      Yprime->clear();
      Zprime->clear();

      Xprime->shrink_to_fit();
      Yprime->shrink_to_fit();
      Zprime->shrink_to_fit();
	
    }//Processing loop 2

  cout << "Fin" << endl;    
  //Finally, write the output file
  output->cd();
  mergedopmap->Write();

    
}//EOF
