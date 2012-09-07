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
void producer_thread(
  sampbuf_sync_t & sampbuf_sync,
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list,
  double & fc
) {
  // Main loop which distributes data to the appropriate subthread.
  //Real_Timer tt;
  double sample_time=-1;
  bool searcher_capbuf_filling=false;
  uint32 searcher_capbuf_idx=0;
  uint32 n_samps_read=0;
  unsigned long long int sample_number=0;
  // Elevate privileges of the producer thread.
  //int retval=nice(-10);
  //if (retval==-1) {
  //  cerr << "Error: could not set elevated privileges" << endl;
  //  exit(-1);
  //}
  //int policy=SCHED_RR;
  //struct sched_param param;
  //param.sched_priority=50;
  ////pthread_getschedparam(pthread_self(), &policy, &param);
  //if (pthread_setschedparam(pthread_self(),policy,&param)) {
  //  cerr << "Error: could not elevate main thread priority" << endl;
  //  exit(-1);
  //}
  //tt.tic();
  while (true) {
    // Each iteration of this loop processes one sample.
    double k_factor;
    double frequency_offset;
    {
      boost::mutex::scoped_lock lock(global_thread_data.frequency_offset_mutex);
      frequency_offset=global_thread_data.frequency_offset;
      k_factor=(fc-frequency_offset)/fc;
      if (mod(sample_number++,19200*100)==0) {
        //tt.toc_print();
        //cout << sample_number << endl;
        // Status message displayed every second.
        //cout << "System frequency offset is currently: " << frequency_offset << endl;
        boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
        list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
        while (it!=tracked_cell_list.tracked_cells.end()) {
          //cout << "Cell ID " << (*(*it)).n_id_cell << " TO: " << setprecision(10) << (*(*it)).frame_timing << endl;
          ++it;
        }
      }
    }

    // Get the next sample
    complex <double> sample;
    {
      boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
      while (sampbuf_sync.fifo.size()<2) {
        sampbuf_sync.condition.wait(lock);
      }
      uint8 real=sampbuf_sync.fifo[0];
      uint8 imag=sampbuf_sync.fifo[1];
      sample=complex <double>((real-127.0)/128.0,(imag-127.0)/128.0);
      sampbuf_sync.fifo.pop_front();
      sampbuf_sync.fifo.pop_front();
    }
    n_samps_read++;
    sample_time+=1.0/k_factor;
    sample_time=WRAP(sample_time,0.0,19200.0);

    // Should we begin filling the capture buffer?
    {
      boost::mutex::scoped_lock lock(capbuf_sync.mutex);
      if ((capbuf_sync.request)&&(!searcher_capbuf_filling)&&(abs(WRAP(sample_time-0,-19200.0/2,19200.0/2))<0.5)) {
        //cout << "searcher data cap beginning" << endl;
        capbuf_sync.request=false;
        searcher_capbuf_filling=true;
        searcher_capbuf_idx=0;
        capbuf_sync.late=WRAP(sample_time-0,-19200.0/2,19200.0/2);
      }
    }

    // Populate the capture buffer
    if (searcher_capbuf_filling) {
      capbuf_sync.capbuf(searcher_capbuf_idx++)=sample;
      if (searcher_capbuf_idx==(unsigned)capbuf_sync.capbuf.size()) {
        // Buffer is full. Signal the searcher thread.
        searcher_capbuf_filling=false;
        boost::mutex::scoped_lock lock(capbuf_sync.mutex);
        capbuf_sync.condition.notify_one();
        //cout << "searcher data cap finished" << endl;
      }
    }

    // Loop for each tracked cell and save data, if necessary. Also delete
    // threads that may have lost lock.
    {
      boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
      list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
      while (it!=tracked_cell_list.tracked_cells.end()) {
        tracked_cell_t & tracked_cell=(*(*it));
        boost::mutex::scoped_lock lock2(tracked_cell.mutex);

        // Delete the tracker if lock has been lost.
        if (tracked_cell.kill_me) {
          // TODO: Fix memory leak here!
          //delete (*it);
          it=tracked_cell_list.tracked_cells.erase(it);
          continue;
        }

        // See if we should start filling the buffer.
        if (tracked_cell.tracker_thread_ready&&!tracked_cell.filling) {
          double tdiff=WRAP(sample_time-(tracked_cell.frame_timing+tracked_cell.target_cap_start_time),-19200.0/2,19200.0/2);
          if (
            // Ideal start time is 0.5 samples away from current time
            (abs(tdiff)<0.5) ||
            // It's possible for the frame timing to change between iterations
            // of the outer sample loop and because of this, it's possible that
            // we missed the best start. Start capturing anyways.
            ((tdiff>0)&&(tdiff<2))
          ) {
            // Configure parameters for this capture
            tracked_cell.filling=true;
            tracked_cell.late=tdiff;
            tracked_cell.buffer_offset=0;
          }
        }

        // Save this sample if our state indicates we are filling the
        // buffer.
        if (tracked_cell.filling) {
          tracked_cell.buffer(tracked_cell.buffer_offset++)=sample;
          if (tracked_cell.buffer_offset==128) {
            // Buffer is full!
            // Send PDU
            td_fifo_pdu_t p;
            p.data=tracked_cell.buffer;
            p.slot_num=tracked_cell.slot_num;
            p.sym_num=tracked_cell.sym_num;
            p.late=tracked_cell.late;
            // Record the frequency offset and frame timing as they were
            // during the capture.
            p.frequency_offset=frequency_offset;
            p.frame_timing=tracked_cell.frame_timing;
            tracked_cell.fifo.push(p);
            tracked_cell.fifo_peak_size=MAX(tracked_cell.fifo.size(),tracked_cell.fifo_peak_size);
            tracked_cell.condition.notify_one();
            //cout << "Sleeping..." << endl;
            //sleep(0.0005/8);
            //cout << "fifo size: " << tracked_cell.fifo.size() << endl;
            // Calculate trigger parameters of next capture
            tracked_cell.filling=false;
            if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
              tracked_cell.target_cap_start_time+=32+128;
            } else {
              tracked_cell.target_cap_start_time+=(tracked_cell.sym_num==6)?128+10:128+9;
            }
            tracked_cell.target_cap_start_time=mod(tracked_cell.target_cap_start_time,19200);
            slot_sym_inc(tracked_cell.n_symb_dl(),tracked_cell.slot_num,tracked_cell.sym_num);
          }
        }
        ++it;
      }
    }
  }
}

