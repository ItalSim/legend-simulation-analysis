void jobsubmit()
{
#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

  for(int i=1;i<11;i++)
    {

      gSystem->Exec(Form("sbatch --export=counter=%i runsims.sh",i));
      sleep(1);
    }
}
