// Copyright 2012 Evrytania LLC (http://www.evrytania.com)
//
// Written by James Peroulas <james@evrytania.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>
#include <itpp/itbase.h>
#include <itpp/signal/transforms.h>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/thread.hpp>
#include <boost/thread/condition.hpp>
#include <list>
#include <sstream>
#include <signal.h>
#include <queue>
#include <curses.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "lte_lib.h"
#include "constants.h"
#include "capbuf.h"
#include "itpp_ext.h"
#include "searcher.h"
#include "dsp.h"
#include "rtl-sdr.h"
#include "LTE-Tracker.h"

using namespace itpp;
using namespace std;

// Process that takes samples and distributes them to the appropriate
// process.
void display_thread(
  sampbuf_sync_t & sampbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list
) {
  global_thread_data.display_thread_id=syscall(SYS_gettid);

  initscr();
  // What does this do?
  (void)echo();
  while (true) {
    clear();
    printw("LTE-Tracker v%i.%i.%i\n\n",MAJOR_VERSION,MINOR_VERSION,PATCH_LEVEL);
    //{
    //  boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
    double frequency_offset=global_thread_data.frequency_offset();
    //cout << "System frequency offset: " << frequency_offset << endl;
    printw("Global FO: %6.0lf",frequency_offset);
    //}

    {
      boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
      //cout << "Raw sample buffer size: " << sampbuf_sync.fifo.size() << "/" << sampbuf_sync.fifo_peak_size << endl;
      printw(" Raw buffer status: %6li/%6li\n",sampbuf_sync.fifo.size(),sampbuf_sync.fifo_peak_size);
    }
    printw("Searcher cycle time: %.1lf s\n\n",global_thread_data.searcher_cycle_time());
    printw("Searcher thread ID: %5i\n",global_thread_data.searcher_thread_id);
    printw("Producer thread ID: %5i\n",global_thread_data.producer_thread_id);
    printw("Main     thread ID: %5i\n",global_thread_data.main_thread_id);
    printw("Display  thread ID: %5i\n\n",global_thread_data.display_thread_id);

    printw("Cell seconds dropped: %i\n\n",global_thread_data.cell_seconds_dropped());

    {
      boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
      list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
      while (it!=tracked_cell_list.tracked_cells.end()) {
        tracked_cell_t & tracked_cell=(*(*it));

        // Deadlock possible???
        boost::mutex::scoped_lock lock1(tracked_cell.fifo_mutex);
        boost::mutex::scoped_lock lock2(tracked_cell.meas_mutex);

        //cout << "Cell ID " << tracked_cell.n_id_cell << " TO: " << setprecision(10) << (*(*it)).frame_timing << " peak input buffer size: " << tracked_cell.fifo.size() << "/" << tracked_cell.fifo_peak_size << endl;
        printw("Cell ID %3i TO: %8.2lf cell buffer status %5i/%5i MIB Failures: %3.0lf\n",
          tracked_cell.n_id_cell,
          tracked_cell.frame_timing(),
          tracked_cell.fifo.size(),
          tracked_cell.fifo_peak_size,
          tracked_cell.mib_decode_failures
        );

        for (uint8 t=0;t<tracked_cell.n_ports;t++) {
          printw("  P%i CRS SP/NP/SNR %5.1lf/%5.1lf/%5.1lf CRS(AVG) SP/NP/SNR %5.1lf/%5.1lf/%5.1lf",
            t,
            db10(tracked_cell.crs_sp(t)),
            db10(tracked_cell.crs_np(t)),
            db10(tracked_cell.crs_sp(t)/tracked_cell.crs_np(t)),
            db10(tracked_cell.crs_sp_av(t)),
            db10(tracked_cell.crs_np_av(t)),
            db10(tracked_cell.crs_sp_av(t)/tracked_cell.crs_np_av(t))
          );
          int8 CB=-1;
          for (uint8 k=1;k<12;k++) {
            if (abs(tracked_cell.ac_fd(k))<=abs(tracked_cell.ac_fd(0))/2) {
              CB=k;
              break;
            }
          }
          if (CB==-1) {
            printw(" >990 kHz\n");
          } else {
            printw(" %4i kHz\n",CB*90);
          }
        }

        /*
        for (uint8 t=0;t<tracked_cell.n_ports;t++) {
          stringstream ss;
          ss << abs(tracked_cell.ac_fd()) << endl;
          printw("%s",ss.str());
        }
        */

        ++it;
      }
    }
    refresh();

    sleep(1);
  }
}

