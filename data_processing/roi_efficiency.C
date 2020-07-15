// Efficiency of ROI
// Run on root: protodune_results/protodune-data-check-roi.root
// ROI is +/-10 channels, +/-120 tbins
// compare gauss filtered results with roi and without
#include <string>
#include <fstream>
{
    TFile f("protodune_results/protodune-data-check-roi.root");
    TIter next(f.GetListOfKeys());
    TKey *key;
    while ((key=(TKey*)next())) 
    {
        string key_name = key->GetName();
        if (key_name.find("gauss") != string::npos)
        {
            printf("key: %s \n", key->GetName());
            TH2F *gauss = (TH2F*)f.Get(key->GetName()); // 2D histogram of floats

            // setup csv
            std::ofstream myfile;
            string outfile_name = "protodune_results/" + key_name + "_roi.csv";
            myfile.open (outfile_name);

            // loop through each bin in hist and save bin contents to csv
            int n_bins_x = gauss->GetNbinsX();
            int n_bins_y = gauss->GetNbinsY();
            for (int xbin=1; xbin<=n_bins_x; xbin++)
            {
                for (int ybin=1; ybin<=n_bins_y; ybin++)
                {
                    // get bin content
                    int bin_id = gauss->GetBin(xbin, ybin);
                    float bin_content = gauss->GetBinContent(bin_id);

                    // save to csv
                    if (bin_content != 0.0){
                        myfile << xbin << "," << ybin << "," << bin_content << "\n";
                    }
                }
            }

            myfile.close();

        }
    }
}
