
#include "WireCellSigProc/RegionOfInterestFilter.h"

#include "WireCellSigProc/Diagnostics.h"

#include "WireCellUtil/Response.h"

#include "WireCellIface/SimpleFrame.h"
#include "WireCellIface/SimpleTrace.h"

#include "WireCellUtil/NamedFactory.h"
//#include "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04d/Linux64bit+3.10-2.17-e19-prof/include/RConfigure.h"
//#include "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_18_04d/Linux64bit+3.10-2.17-e19-prof/include/ROOT/RConfig.hxx"
//#include "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_22_08d/Linux64bit+3.10-2.17-e20-p392-prof/include/TH2F.h"
//#include "/cvmfs/larsoft.opensciencegrid.org/products/root/v6_22_08d/Linux64bit+3.10-2.17-e20-p392-prof/include/TFile.h"
#include "TFile.h"
#include "TH2F.h"
// #include "WireCellUtil/ExecMon.h" // debugging

//#include "FrameUtils.h"          // fixme: needs to move to somewhere more useful.

#include <unordered_map>
#include <random>
#define PEAK 20
#define ROI_t 120
#define ROI_ch 10

// testing
#include <iostream>
using namespace std;

WIRECELL_FACTORY(RegionOfInterestFilter, WireCell::SigProc::RegionOfInterestFilter,
                 WireCell::IFrameFilter, WireCell::IConfigurable)

using namespace WireCell;

using namespace WireCell::SigProc;

RegionOfInterestFilter::RegionOfInterestFilter(const std::string& roi_tag, const std::string& old_tag)
: m_roi_tag(roi_tag)
, m_old_tag(old_tag)
, m_frame_tag("sigproc")
, log(Log::logger("sigproc"))
{
}
RegionOfInterestFilter::~RegionOfInterestFilter()
{
}

void RegionOfInterestFilter::configure(const WireCell::Configuration& cfg)
{
  m_roi_tag = get(cfg,"roi_tag",m_roi_tag);   
  m_old_tag = get(cfg,"old_tag",m_old_tag);   
  m_frame_tag = get(cfg,"frame_tag",m_frame_tag); 
}

WireCell::Configuration RegionOfInterestFilter::default_configuration() const
{
    Configuration cfg;    
    cfg["roi_tag"] = m_roi_tag;
    cfg["old_tag"] = m_old_tag;
    cfg["frame_tag"] = m_frame_tag;
    return cfg;
}

static bool ispeak(float x) { return (x <0 or x>0); }
static bool isZero(float x) { return x == 0.0; }

bool RegionOfInterestFilter::operator()(const input_pointer& inframe, output_pointer& outframe)
{
  log->debug("RegionOfInterestFilter: inside operator roi_tag {}, old_tag {}", m_roi_tag, m_old_tag);

  if (!inframe) {             // eos
    outframe = nullptr;
    return true;
  }
  // Convert to OSP cmm indexed by OSB sequential channels, NOT WCT channel ID.
  m_cmm.clear();
  for (auto cm : inframe->masks())
  {  //in->masks give back a map (ChannelMaskMap) done with a string (in name) and on other map (in m) ChannelMasks
  const std::string name = cm.first;
  for (auto m: cm.second)
  {   //m is the map ChannelMasks done with an int (wct_channel_ident) and a BinRangeList
    const int wct_channel_ident = m.first;
    const OspChan& och = m_channel_map[wct_channel_ident]; //from config
    if (och.plane < 0) continue;               // in case user gives us multi apa frame
      m_cmm[name][och.channel] = m.second;  //m.second è la BinRangeList
    }
  }

  IFrame::trace_list_t roi_traces, old_traces;
  ITrace::vector* newtraces = new ITrace::vector; // will become shared_ptr.

  // outframe = inframe;

  auto traces = inframe->traces();
  auto trace_tags = inframe->trace_tags();

  log->debug("RegionOfInterestFilter: num tags {}", trace_tags.size());

  for (auto ttag : trace_tags)
  log->debug("RegionOfInterestFilter: tag {}", ttag);

  int num_channels = traces->size(); // number of channels
  ITrace::ChargeSequence old_array[num_channels]; // charges of traces, stored as array with sequential channels
  ITrace::ChargeSequence ROI_array[num_channels]; // ROI charges of traces, stored as array

  int lowest_ch = -1;
  int num_bins = traces->at(0)->charge().size();

  // Get lowest channel number and lowest number of bins
  for (auto trace : *traces.get())
  {
    int channel = trace->channel();
    if (channel<lowest_ch or lowest_ch==-1) 
      lowest_ch = channel;
    if ((int)trace->charge().size() < num_bins) 
      num_bins = trace->charge().size();
  }

  // cout << "lowest_ch: " << lowest_ch << "\n";
  int num_store_ch[num_bins] = {}; // number of upper ch roi to store per bin
  int upper_time_roi_counter = 0; // counts time bins above a peak to skip rewriting bkg

  TFile bkgFile("/dune/data/users/mman/wire-cell-toolkit/roi_analysis/background.root", "READ"); // <--

  TH2F *h1 = (TH2F*)bkgFile.Get("hMerged"); // "hMerged"

  for (auto trace : *traces.get())
  {
    int channel = trace->channel();
    // int tbin = trace->tbin();
    auto const& charges = trace->charge();
    int bkg_channel = channel%15320; //TODO workaround for applying bkg based on different geometry

    ITrace::ChargeSequence newcharge(charges.size(), 0.0);
    ITrace::ChargeSequence oldcharge(charges.size(), 0.0);

    auto chargessort = charges;

    std::nth_element(chargessort.begin(), chargessort.begin() + chargessort.size()/2, chargessort.end());
    float median = chargessort[chargessort.size()/2];

    //log->debug("RegionOfInterestFilter: channel {}, initial time {}, size {}", channel, tbin, (int)charges.size());  


    // Define random generator with Gaussian distribution
    const float mean = 0.0;
    const float stddev = 10;
    std::default_random_engine generator;
    std::normal_distribution<float> dist(mean, stddev);

    // Time ROI
    for (int bin = 0; bin < (int)charges.size(); ++bin)
    {
      float central_value = charges[bin]- median;
      oldcharge.at(bin) = central_value;//charges[bin]; //central_value;
      if(upper_time_roi_counter > 0){
        upper_time_roi_counter = upper_time_roi_counter - 1;
      }
      else{
        newcharge.at(bin) = h1->GetBinContent(bkg_channel, bin);
      }

      if(central_value<-PEAK or central_value>PEAK)
      {
        //log->debug("RegionOfInterestFilter: peak in the bin {} = {}, median {}, Cvalue {}, ispeak {}", bin, charges[bin], median, central_value, ispeak(charges[bin]) );
        for(int delta = -ROI_t; delta < ROI_t; ++delta)
        {
          int newbin = bin+delta;
          if(newbin>-1 and newbin<(int)charges.size())
            newcharge.at(newbin) = charges[newbin]- median;//1000;
        }
        upper_time_roi_counter = ROI_t;
      }
    }

    // Write to charge arrays
    ROI_array[channel-lowest_ch] = newcharge;
    old_array[channel-lowest_ch] = oldcharge;
  }
  //cout << "Time ROI complete \n";

  // Channel ROI
  for (int ch_ind = 1; ch_ind < num_channels; ch_ind++)
  {
    for(int bin=0; bin<num_bins; bin++)
    {
      // Lower channel ROI
      // peak in channel and zero in channel-1 (start of track)
      if(ispeak(ROI_array[ch_ind].at(bin)) and isZero(ROI_array[ch_ind-1].at(bin)))
      {
        // fill lower channel ROI
        for(int j=-ROI_ch; j<0; j++)
        {
          int update_channel = ch_ind+j;
          if(update_channel>-1)
          {
            ROI_array[update_channel].at(bin) = old_array[update_channel].at(bin);
          }
        }
      }

      // Upper channel ROI
      // zero in channel and peak in channel-1 (end of track)
      if(num_store_ch[bin]==0 and isZero(ROI_array[ch_ind].at(bin)) and ispeak(ROI_array[ch_ind-1].at(bin)))
      {
        // fill upper channel ROI
        ROI_array[ch_ind].at(bin) = old_array[ch_ind].at(bin);
        num_store_ch[bin] = ROI_ch;
      }
      else if (num_store_ch[bin]>0)
        ROI_array[ch_ind].at(bin) = old_array[ch_ind].at(bin);

      // iterate num_store_ch
      if(num_store_ch[bin]>1)
        num_store_ch[bin] -= 1; // >0 means store that many more channels in ROI
      else if(num_store_ch[bin]==1)
        num_store_ch[bin] = -1; // -1 means end of upper ROI
      else if(num_store_ch[bin]==-1)
        num_store_ch[bin] = 0; // 0 means ready to find next end of track
    }
  }
  // cout << "Channel ROI complete \n";

  // Write ROI_array to newtraces
  for (auto trace : *traces.get())
  {
    int channel = trace->channel();
    int tbin = trace->tbin();
    ITrace::ChargeSequence newcharge = ROI_array[channel-lowest_ch];

    SimpleTrace *tracetemp = new SimpleTrace(channel, tbin, newcharge);
    const size_t roi_trace_index = newtraces->size();
    roi_traces.push_back(roi_trace_index);
    newtraces->push_back(ITrace::pointer(tracetemp));

    //const size_t old_trace_index = newtraces->size();
    //old_traces.push_back(old_trace_index);
    //newtraces->push_back(ITrace::pointer(trace));
  }

  SimpleFrame* sframe = new SimpleFrame(inframe->ident(), inframe->time(),
  ITrace::shared_vector(newtraces),
  inframe->tick(), inframe->masks());
  sframe->tag_frame(m_frame_tag);
  sframe->tag_traces(m_roi_tag, roi_traces);
  //sframe->tag_traces(m_old_tag, old_traces);

  outframe = IFrame::pointer(sframe);

  //auto checktraces = outframe->traces();

  // for (auto checktrace : *checktraces.get())
  // {
  //     int channel = checktrace->channel();
  //     int tbin = checktrace->tbin();
  //     auto const& charges = checktrace->charge();

  //     log->debug("RegionOfInterestFilter: newtraces channel {}, start time {}, size {}", channel, tbin,(int)charges.size());

  //     for (int bin = 0; bin < (int)charges.size(); ++bin)
  //     {
  //       log->debug("RegionOfInterestFilter: newtrace bin {} value = {}", bin, charges[bin] );
  //     }
  // }

  // outframe = IFrame::pointer(sframe);
  log->debug("RegionOfInterestFilter: end operator");
  //bkgFile->Close();
  return true;
}


// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
