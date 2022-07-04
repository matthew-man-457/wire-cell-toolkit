/** Remove all possible noise from a microboone-like detector.
 *
 * This filter is a kitchen sink class and is a candidate for
 * factoring.
 */

#ifndef WIRECELLSIGPROC_REGIONOFINTERESTFILTER
#define WIRECELLSIGPROC_REGIONOFINTERESTFILTER

#include "WireCellIface/IFrameFilter.h"
#include "WireCellIface/IConfigurable.h"
#include "WireCellIface/IChannelNoiseDatabase.h"
#include "WireCellIface/IChannelFilter.h"

#include "WireCellUtil/Waveform.h"
#include "WireCellUtil/Logging.h"

#include <vector>
#include <map>
#include <string>

namespace WireCell {
    namespace SigProc {

	/**
	   The Omnibus Noise Filter applies two series of IChannelFilter
	   objects. The first series is applied on a per-channel basis and
	   the second is applied on groups of channels as determined by
	   its channel grouping.
	*/
	class RegionOfInterestFilter : public WireCell::IFrameFilter, public WireCell::IConfigurable {
	public:
	    typedef std::vector< std::vector<int> > grouped_channels_t;

	    /// Create an RegionOfInterestFilter.
	    RegionOfInterestFilter(const std::string& roi_tag = "roi", const std::string& old_tag = "raw");
	    virtual ~RegionOfInterestFilter();

	    /// IFrameFilter interface.
	    virtual bool operator()(const input_pointer& in, output_pointer& out);

	    /// IConfigurable interface.
	    virtual void configure(const WireCell::Configuration& config);
	    virtual WireCell::Configuration default_configuration() const;



	private:

      struct OspChan {
        int channel;            // between 0 and nwire_u+nwire_v+nwire_w-1
        int wire;               // between 0 and nwire_{u,v,w,}-1 depending on plane
        int plane;              // 0,1,2
        int ident;              // wct ident, opaque non-negative number.  set in wires geom file
        OspChan(int c=-1, int w=-1, int p=-1, int id=-1) : channel(c), wire(w), plane(p), ident(id) {}
        std::string str() const;
      };
      
		std::string m_roi_tag;
		std::string m_old_tag;
		std::string m_frame_tag;
        Log::logptr_t log;
		Waveform::ChannelMaskMap m_cmm; 

		// Need to go from WCT channel ident to {OSP channel, wire and plane}
		std::map<int,OspChan> m_channel_map;
	};

    }
}
#endif

// Local Variables:
// mode: c++
// c-basic-offset: 4
// End:
