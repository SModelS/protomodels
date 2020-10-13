int Test()
{
TFile *MyFile = new TFile("CMS-SUS-19-006_Figure_013-a.root","read");
TH2F * h1 = new TH2F("h1","h1 title" , 128, 600, 2500,160,0,2000);
h1 = (TH2F*)MyFile->Get("MassScan2D");



ofstream file;

file.open("result.txt");
for(int i=0 ; i<128 ; i++)
{
for( int j=0 ; j<160 ;j++)

{
	file<<i*(2500-600)/128+600<<" "<<j*(2000)/160<<" "<<h1->GetBinContent(i,j)<<std::endl;

}

file<<" "<<std::endl;
}

file.close();
return 0;
}
