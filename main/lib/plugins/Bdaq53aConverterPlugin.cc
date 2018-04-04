/* Raw data conversion for the Bdaq53a producer
 *
 *
 * jorge.duarte.campderros@cern.ch
 * 
 */

#include "eudaq/DataConverterPlugin.hh"
#include "eudaq/StandardEvent.hh"
#include "eudaq/Utils.hh"
#include "eudaq/RawDataEvent.hh"
#include "eudaq/Timer.hh"
#include "eudaq/Logger.hh"
#include "eudaq/Bdaq53a.hh"

#include <stdlib.h>

// All LCIO-specific parts are put in conditional compilation blocks
// so that the other parts may still be used if LCIO is not available.
#if USE_LCIO
#  include "IMPL/LCEventImpl.h"
#  include "IMPL/TrackerRawDataImpl.h"
#  include "IMPL/LCCollectionVec.h"
#  include "lcio.h"
#  include "IMPL/TrackerDataImpl.h"
#  include "IMPL/LCGenericObjectImpl.h"
#  include "UTIL/CellIDEncoder.h"
#endif

#if USE_EUTELESCOPE
#  include "EUTELESCOPE.h"
#  include "EUTelSetupDescription.h"
#  include "EUTelEventImpl.h"
#  include "EUTelTrackerDataInterfacerImpl.h"
#  include "EUTelGenericSparsePixel.h"
#  include "EUTelAPIXMCDetector.h"
#  include "EUTelRunHeaderImpl.h"
using eutelescope::EUTELESCOPE;
#endif

namespace eudaq 
{
    // The event type for which this converter plugin will be registered
    // Modify this to match your actual event type (from the Producer)
    static const char* EVENT_TYPE = "bdaq53a";

#if USE_LCIO && USE_EUTELESCOPE
    static const int chip_id_offset = 20;
#endif

    // Declare a new class that inherits from DataConverterPlugin
    class Bdaq53aConverterPlugin : public DataConverterPlugin
    {
        public:
            // This is called once at the beginning of each run.
            // You may extract information from the BORE and/or configuration
            // and store it in member variables to use during the decoding later.
            virtual void Initialize(const Event & bore, const Configuration & cnf) 
            {
#ifndef WIN32
                (void)cnf; // just to suppress a warning about unused parameter cnf
#endif
            }
            
            // This should return the trigger ID (as provided by the TLU)
            // if it was read out, otherwise it can either return (unsigned)-1,
            // or be left undefined as there is already a default version.
            virtual unsigned GetTriggerID(const Event & ev) const 
            {
                // Make sure the event is of class RawDataEvent
                if(const RawDataEvent * rev = dynamic_cast<const RawDataEvent *> (&ev)) 
                {
                    if(rev->IsBORE() || rev->IsEORE())
                    {
                        return 0;
                    }
                    bool is_data_header = false;
                    for(unsigned int i = 0; i < rev->NumBlocks(); ++i)
                    {
                        const RawDataEvent::data_t & word = rev->GetBlock(i);
                        
                        // Obtain the data word (could be at the FE high or low word)
                        //const RawDataEvent::data_t & data_word = BDAQ53A_DATA_HEADER_MACRO(word))
                    }
                    if(rev->NumBlocks() > 0) 
                    {
                        return -1;
                        //return BDAQ53A_DATA_HEADER_TRIGGER_ID_MACRO(rev->GetBlock(0));
                    }
                }
                // Not proper event or not enable to extract Trigger ID
                return (unsigned) -1;
            }

            // Here, the data from the RawDataEvent is extracted into a StandardEvent.
            // The return value indicates whether the conversion was successful.
            // Again, this is just an example, adapted it for the actual data layout.
            virtual bool GetStandardSubEvent(StandardEvent & sev, const Event & ev) const 
            {
                if(ev.IsBORE() || ev.IsEORE())
                {
                    // nothing to do
                    return true;
                }

                // If we get here it must be a data event
                const RawDataEvent & ev_raw = dynamic_cast<const RawDataEvent &>(ev);
                std::cout << "Bdaq53aConverterPlugin NumBlocks="<<ev_raw.NumBlocks() << std::endl;
                for(size_t i = 0; i < ev_raw.NumBlocks(); ++i) 
                {
                    std::cout << "BDAQ53ACP EV_RAW.GetID: " << ev_raw.GetID(i) << " i:" << i 
                        << " size raw_data: " << ev_raw.GetBlock(i).size() << std::endl;
                    const RawDataEvent::data_t & raw_data = ev_raw.GetBlock(i);
                    std::cout << " Hoow " << ( raw_data & 0x00010000 ) << std::endl;
                    //--sev.AddPlane(ConvertPlane(ev_raw.GetBlock(i), ev_raw.GetID(i)));
                }
                return true;
            }

#if USE_LCIO && USE_EUTELESCOPE
            // This is where the conversion to LCIO is done
            virtual lcio::LCEvent * GetLCIOEvent(const Event * /*ev*/) const 
            {
                return 0;
            }
            
            virtual bool GetLCIOSubEvent(lcio::LCEvent & lcioEvent, const Event & eudaqEvent) const 
            {
                //std::cout << "getlciosubevent (I4) event " << eudaqEvent.GetEventNumber() << " | " << GetTriggerID(eudaqEvent) << std::endl;
                if (eudaqEvent.IsBORE()) 
                {
                    // shouldn't happen
                    return true;
                } 
                else if (eudaqEvent.IsEORE()) 
                {
                    // nothing to do
                    return true;
                }
                
                // set type of the resulting lcio event
                lcioEvent.parameters().setValue( eutelescope::EUTELESCOPE::EVENTTYPE, eutelescope::kDE );
                // pointer to collection which will store data
                LCCollectionVec * zsDataCollection;
                // it can be already in event or has to be created
                bool zsDataCollectionExists = false;
                try 
                {
                    zsDataCollection = static_cast< LCCollectionVec* > ( lcioEvent.getCollection( "zsdata_apix" ) );
                    zsDataCollectionExists = true;
                } 
                catch( lcio::DataNotAvailableException& e ) 
                {
                    zsDataCollection = new LCCollectionVec( lcio::LCIO::TRACKERDATA );
                }
                //	create cell encoders to set sensorID and pixel type
                CellIDEncoder< TrackerDataImpl > zsDataEncoder   ( eutelescope::EUTELESCOPE::ZSDATADEFAULTENCODING, zsDataCollection  );
                
                // this is an event as we sent from Producer
                // needs to be converted to concrete type RawDataEvent
                const RawDataEvent & ev_raw = dynamic_cast <const RawDataEvent &> (eudaqEvent);
                
                std::vector< eutelescope::EUTelSetupDescription * >  setupDescription;

        for (size_t chip = 0; chip < ev_raw.NumBlocks(); ++chip) {
          const std::vector <unsigned char> & buffer=dynamic_cast<const std::vector<unsigned char> &> (ev_raw.GetBlock(chip));

          if (lcioEvent.getEventNumber() == 0) {
            eutelescope::EUTelPixelDetector * currentDetector = new eutelescope::EUTelAPIXMCDetector(2);
            currentDetector->setMode( "ZS" );

            setupDescription.push_back( new eutelescope::EUTelSetupDescription( currentDetector )) ;
          }

          zsDataEncoder["sensorID"] = ev_raw.GetID(chip) + chip_id_offset + first_sensor_id; // formerly 14
          zsDataEncoder["sparsePixelType"] = eutelescope::kEUTelGenericSparsePixel;

          // prepare a new TrackerData object for the ZS data
          // it contains all the hits for a particular sensor in one event
          std::unique_ptr<lcio::TrackerDataImpl > zsFrame( new lcio::TrackerDataImpl );
          // set some values of "Cells" for this object
          zsDataEncoder.setCellID( zsFrame.get() );

          // this is the structure that will host the sparse pixel
          // it helps to decode (and later to decode) parameters of all hits (x, y, charge, ...) to
          // a single TrackerData object (zsFrame) that will correspond to a single sensor in one event
          std::unique_ptr< eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > >
            sparseFrame( new eutelescope::EUTelTrackerDataInterfacerImpl< eutelescope::EUTelGenericSparsePixel > ( zsFrame.get() ) );

          unsigned int ToT = 0;
          unsigned int Col = 0;
          unsigned int Row = 0;
          unsigned int lvl1 = 0;

          for (unsigned int i=4; i < buffer.size(); i += 4) {
            unsigned int Word = (((unsigned int)buffer[i]) << 24) | (((unsigned int)buffer[i+1]) << 16) | (((unsigned int)buffer[i+2]) << 8) | (unsigned int)buffer[i+3];

            if (BDAQ53A_DATA_HEADER_MACRO(Word)) {
              lvl1++;
            } else {
              // First Hit
              if (getHitData(Word, false, Col, Row, ToT)) {
                sparseFrame->emplace_back( Col, Row, ToT, lvl1-1 );
              }
              // Second Hit
              if (getHitData(Word, true, Col, Row, ToT)) {
                sparseFrame->emplace_back( Col, Row, ToT, lvl1-1 );
              }
            }
          }

          // write TrackerData object that contains info from one sensor to LCIO collection
          zsDataCollection->push_back( zsFrame.release() );
        }

        // add this collection to lcio event
        if ( ( !zsDataCollectionExists )  && ( zsDataCollection->size() != 0 ) ) lcioEvent.addCollection( zsDataCollection, "zsdata_apix" );

        if (lcioEvent.getEventNumber() == 0) {
          // do this only in the first event
          LCCollectionVec * apixSetupCollection = NULL;

          bool apixSetupExists = false;
          try {
            apixSetupCollection = static_cast< LCCollectionVec* > ( lcioEvent.getCollection( "apix_setup" ) ) ;
            apixSetupExists = true;
          } catch (...) {
            apixSetupCollection = new LCCollectionVec( lcio::LCIO::LCGENERICOBJECT );
          }

          for ( size_t iPlane = 0 ; iPlane < setupDescription.size() ; ++iPlane ) {
            apixSetupCollection->push_back( setupDescription.at( iPlane ) );
          }

          if (!apixSetupExists) lcioEvent.addCollection( apixSetupCollection, "apix_setup" );
        }
        return true;

      }
#endif

    private:

      // The constructor can be private, only one static instance is created
      // The DataConverterPlugin constructor must be passed the event type
      // in order to register this converter for the corresponding conversions
      // Member variables should also be initialized to default values here.
      Bdaq53aConverterPlugin() : DataConverterPlugin(EVENT_TYPE) {}

      // The single instance of this converter plugin
      static Bdaq53aConverterPlugin m_instance;
  };

  // Instantiate the converter plugin instance
  Bdaq53aConverterPlugin Bdaq53aConverterPlugin::m_instance;

} // namespace eudaq
