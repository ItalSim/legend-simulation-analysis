void opticalmappostproc(void)
{

  //Post-processing script for the optical maps generated in the LAr between
  //the LEGEND-1000  plastic neutron shield and the steel cryostat.

  //The map consists of a series of cubic voxels, but the liquid argon volume
  //is only flat and orthoganal on three of the six sides: the side in contact
  //with the plastic shield, the top (artificial boundary), and the bottom
  //(artificial boundary). The side in contact with the steel cryostat is
  //curved, and the two sides with edges touching both the shield and the
  //cryostat (artificial boundaries) are flat but not orthoganal/parallel.

  //As a result, some of the voxels will be 'incomplete', if these voxels
  //are intersected by one of the latter bounding surfaces. It's technically
  //possible to calculate the volume of each voxel which is filled/missing,
  //given the boundary surfaces, but this would be tedious and subject to
  //change each time the map is generated with a different number of voxels.

  //Instead, a heuristic approach will be pursued. The number of incomplete
  //voxels will be determined, and then each incomplete voxel will be filled
  //with the value of the adjacent complete voxel.

  //Small point of order: it is possible for two adjacent voxels to both be
  //incomplete, if the bounding surface cuts through in a certain way:


  //====================================================================
  //||                               ||  /                            ||
  //||                               || /                             ||
  //||                               ||/                              ||
  //||                               |/                               ||
  //||                               /|                               ||
  //||                              /||                               ||
  //||                             / ||                               ||
  //||                            /  ||                               ||
  //||                           /   ||                               ||
  //||                          /    ||                               ||
  //||                         /     ||                               ||
  //||                        /      ||                               ||
  //||                       /       ||                               ||
  //||                      /        ||                               ||
  //====================================================================



  
  //However, the bounding surfaces are 'close' to parallel, even though they
  //are not parallel. According to my calculations, in the absolute worst
  //case, where the bounding surface just barely cuts through one of the
  //adjacent voxels, and then maximally cuts through the other, the more
  //full voxel will still be 86.6% complete. So, we will treat this voxel
  //as though it were complete, and pass its value along to the less
  //complete voxel.

  //There is another essential question to answer: how do we determine which
  //voxels are incomplete? It would be easier to do 'by eye': we would just
  //locate the voxels with abnormally low efficiency values compared to the
  //neighboring voxels. However, there is another trick we can use, where
  //we take advantage of the fact that the map is 'too large'. Next to the
  //most incomplete voxel in each row, there will be an empty voxel (the
  //unused space at the map's edge). So, if we find an incomplete voxel,
  //and we know which 'side' it is closest to, we know which direction to
  //inherit a complete voxel value from.

  //In the event that the user has very large voxels (like in the test
  //files...), the voxels at the edges may not be empty. However, they will
  //at the very least be incomplete. For safety, we'll make all edge voxels
  //inherit probability values from adjacent voxels if they're not empty.

  //Finally, note the directionality dependence of the adjacency. Voxels
  //on the left side will have negative x values, and should look left for
  //emptiness, and voxels on the right side will have positive x values,
  //and should look right for empty voxels. For voxels near the cryostat,
  //it's the same, but with forward/backward. Because of the directionality
  //dependence, the order of operations may change.

  
  //--------------------------------------------------------------------------

      
  //Open the built, merged optical map
  TFile *input = new TFile("/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap/4-innermap-Dec2023/output/merged/mergedmap.root","read");
  TTree *maptree = (TTree*)input->Get("opmap");

  double prob = 0;

  TFile *output = new TFile("/lfs/l1/legend/users/cbarton/simulations/campaigns/opticalmap/4-innermap-Dec2023/output/map.root","recreate");
  output->cd();
  TTree *map = maptree->CloneTree();
  //Add a branch for the 'fixed' probability values
  TBranch *probbranch = map->Branch("prob",&prob,"prob/D");

  input->cd();  

  //We only need the index and the probability for what we want to accomplish
  
  int index = 0;
  int xindex = 0;
  int yindex = 0;
  int zindex = 0;
  double probability = 0;
  
  maptree->SetBranchAddress("index",&index);
  maptree->SetBranchAddress("probability",&probability);
  
  int entries = maptree->GetEntries();
  cout << "Number of entries in file: " << entries << endl;
  
  
  //Manually add the number of bins in each dimension
  int xbins = 54;
  int ybins = 55;
  int zbins = 151;

  //Master processing loop

  for(int i = 0; i < entries; i++)
    {

      if(i%xbins==0)
	cout << endl;

      maptree->GetEntry(i);
     
      if(probability == 0)
	{
	  prob = 0;
	  cout <<"0";
	  //cout << prob << endl;
	  probbranch->BackFill();
	  continue;	
	}
      
      else
	{
	  //If the map is built correctly, i and index should have the same value at all times	  

	  //Scenario 1: account for all the edge voxels
	  
	  //Edge voxels, first case: along the 'left' side of the optical map
	  if(i%xbins == 0)
	    {
	      maptree->GetEntry(i+1);
	      prob = probability;
	      cout <<"L";
	      //cout << prob << endl;
	      probbranch->BackFill();
	      continue;	
	    }
      
	  //Edge voxels, second case: along the 'right' side of the optical map
	  else if(i%xbins == (xbins-1))
	    {
	      maptree->GetEntry(i-1);
	      prob = probability;
	      cout <<"R";
	      //cout << prob << endl;
	      probbranch->BackFill();
	      continue;		  
	    }
      
	  //Edge voxels, third case: along the cryostat wall
	  else if((int)floor(i/xbins)%ybins == (ybins-1))
	    {
	      maptree->GetEntry(i - xbins);//Roll back one entire row
	      prob = probability;
	      cout <<"W";
	      //cout << prob << endl;
	      probbranch->BackFill();
	      continue;
	    }

      
	  //Scenario 2: not an edge voxel
	  //However, we still need to check if it's adjacent to an empty voxel
	  //Note that we don't need to worry as much about the index going out of range,
	  //because three of the edges are accounted for. Our last concern is with the edge
	  //which doesn't need to be fixed, so be wary of that.
	  else
	    {
	  
	      maptree->GetEntry(i-1);
	      if(probability == 0)
		{//Grab from the other side
		  maptree->GetEntry(i+1);
		  prob = probability;
		  cout <<"l";
		  //cout << prob << endl;
		  probbranch->BackFill();
		  continue;
		}
	      maptree->GetEntry(i+1);
	      if(probability == 0)
		{//Grab from the other side
		  maptree->GetEntry(i-1);
		  prob = probability;
		  cout <<"r";
		  //cout << prob << endl;
		  probbranch->BackFill();
		  continue;
		}
	      //if(i > (entries - xbins) || i < xbins)
	      //{
	      //  cout<<"A";
	      //}
	      maptree->GetEntry(i+xbins);
	      if(probability == 0)
		{//Grab from the other side
		  maptree->GetEntry(i-xbins);
		  prob = probability;
		  cout <<"w";
		  //cout << prob << endl;
		  probbranch->BackFill();
		  continue;
		}
	      //If absolutely nothing else, then probability stays as it is
	  
	      maptree->GetEntry(i);
	      prob = probability;	  
	      cout<<"A";
	      //cout << prob << endl;
	      probbranch->BackFill();
	      continue;
	    }//else
      
	}//If probability isn't 0

      cout <<"You shouldn't be seeing this" << endl;
			               
    }//Main processing loop

  
  cout << endl;

  //Make new output file and swap out the probability branches
  output->cd();  
  map->SetBranchStatus("probability",0);
  map->Print();
  map->Write();
  output->Write();
}//EOF
