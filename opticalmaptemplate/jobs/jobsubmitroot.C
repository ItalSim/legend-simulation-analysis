void jobsubmitroot()
{
#include<time.h>
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include <chrono>
#include <thread>

  for(int i=0;i<100;i++)
  {
    
    gSystem->Exec(Form("sbatch --export=counter=%i submitroot.sh",i));
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }
}
