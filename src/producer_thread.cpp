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

// This is local storage for each cell that is being tracked.
typedef struct {
  uint32 serial_num;
  uint32 target_cap_start_time;
  bool filling;
  uint16 buffer_offset;
  td_fifo_pdu_t pdu;
} cell_local_t;

// Process that takes samples and distributes them to the appropriate
// process.
void producer_thread(
  sampbuf_sync_t & sampbuf_sync,
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list,
  double & fc
) {
  global_thread_data.producer_thread_id=syscall(SYS_gettid);

  // Main loop which distributes data to the appropriate subthread.
  // Local storage for each cell.
  cell_local_t cell_local[504];
  for (uint16 t=0;t<504;t++) {
    cell_local[t].serial_num=0;
  }

  //Real_Timer tt;
  double sample_time=-1;
  bool searcher_capbuf_filling=false;
  uint32 searcher_capbuf_idx=0;
  //unsigned long long int sample_number=0;
  // Elevate privileges of the producer thread.
  //int retval=nice(-10);
  //if (retval==-1) {
  //  cerr << "Error: could not set elevated privileges" << endl;
  //  ABORT(-1);
  //}
  //int policy=SCHED_RR;
  //struct sched_param param;
  //param.sched_priority=50;
  ////pthread_getschedparam(pthread_self(), &policy, &param);
  //if (pthread_setschedparam(pthread_self(),policy,&param)) {
  //  cerr << "Error: could not elevate main thread priority" << endl;
  //  ABORT(-1);
  //}
  //tt.tic();
#define BLOCK_SIZE 10000
  while (true) {
    // Each iteration of this loop processes one block of data.
    const double frequency_offset=global_thread_data.frequency_offset();
    const double k_factor=(global_thread_data.fc_requested-frequency_offset)/global_thread_data.fc_programmed;
    //const double k_factor_inv=1/k_factor;
    const double & fs_programmed=global_thread_data.fs_programmed;

    // Get the next block
    //complex <double> sample;
    cvec samples(BLOCK_SIZE);
    vec samples_timestamp(BLOCK_SIZE);
    uint32 n_samples;
    {
      boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
      while (sampbuf_sync.fifo.size()<2) {
        sampbuf_sync.condition.wait(lock);
      }
      // Dump data if there is too much in the fifo
      while (sampbuf_sync.fifo.size()>2*FS_LTE/16*1.5) {
        for (uint32 t=0;t<(unsigned)round_i(fs_programmed*k_factor);t++) {
          sampbuf_sync.fifo.pop_front();
        }
        global_thread_data.raw_seconds_dropped_inc();
      }
      n_samples=BLOCK_SIZE;
      complex <double> sample_temp;
      for (uint16 t=0;t<BLOCK_SIZE;t++) {
        if (sampbuf_sync.fifo.size()<2) {
          n_samples=t;
          break;
        }
        sample_temp.real()=(sampbuf_sync.fifo.front()-127.0)/128.0;
        sampbuf_sync.fifo.pop_front();
        sample_temp.imag()=(sampbuf_sync.fifo.front()-127.0)/128.0;
        sampbuf_sync.fifo.pop_front();
        samples(t)=sample_temp;
        sample_time+=(FS_LTE/16)/(fs_programmed*k_factor);
        //sample_time=itpp_ext::matlab_mod(sample_time,19200.0);
        if (sample_time>19200.0)
          sample_time-=19200.0;
        samples_timestamp(t)=sample_time;
      }
    }

    // Handle the searcher capture buffer
    for (uint32 t=0;t<n_samples;t++) {
      if ((capbuf_sync.request)&&(abs(WRAP(samples_timestamp(t)-0,-19200.0/2,19200.0/2))<0.5)) {
        //cout << "searcher data cap beginning" << samples_timestamp(t) << endl;
        capbuf_sync.request=false;
        searcher_capbuf_filling=true;
        searcher_capbuf_idx=0;
        capbuf_sync.late=WRAP(samples_timestamp(t)-0,-19200.0/2,19200.0/2);
      }

      // Populate the capture buffer
      if (searcher_capbuf_filling) {
        capbuf_sync.capbuf(searcher_capbuf_idx++)=samples(t);
        if (searcher_capbuf_idx==(unsigned)capbuf_sync.capbuf.size()) {
          // Buffer is full. Signal the searcher thread.
          searcher_capbuf_filling=false;
          boost::mutex::scoped_lock lock(capbuf_sync.mutex);
          capbuf_sync.condition.notify_one();
          //cout << "searcher data cap finished" << endl;
        }
      }
    }

    // Loop for each tracked cell and save data, if necessary. Also delete
    // threads that may have lost lock.
    {
      boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
      list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
      while (it!=tracked_cell_list.tracked_cells.end()) {
        tracked_cell_t & tracked_cell=(*(*it));
        // See if this thread has been launched yet. If not, launch it.
        if (!tracked_cell.launched) {
          tracked_cell.thread=boost::thread(tracker_thread,boost::ref(tracked_cell),boost::ref(global_thread_data));
          tracked_cell.launched=true;
        }
        double frame_timing=tracked_cell.frame_timing();

        cell_local_t & cl=cell_local[tracked_cell.n_id_cell];

        // Initialize local storage if necessary
        if (tracked_cell.serial_num!=cl.serial_num) {
          cl.serial_num=tracked_cell.serial_num;
          cl.pdu.slot_num=0;
          cl.pdu.sym_num=0;
          cl.target_cap_start_time=(tracked_cell.cp_type==cp_type_t::NORMAL)?10:32;
          cl.filling=false;
          cl.buffer_offset=0;
          if (cl.serial_num==1)
            cl.pdu.data.set_size(128);
        }

        // Delete the tracker if lock has been lost.
        if (tracked_cell.kill_me) {
          tracked_cell_t * temp=(*it);
          it=tracked_cell_list.tracked_cells.erase(it);
          delete temp;
          continue;
        }

        // Loop for each sample in the buffer.
        for (uint32 t=0;t<n_samples;t++) {
          // See if we should start filling the buffer.
          if (tracked_cell.tracker_thread_ready&&!cl.filling) {
            double tdiff=WRAP(samples_timestamp(t)-(frame_timing+cl.target_cap_start_time),-19200.0/2,19200.0/2);
            if (
              // Ideal start time is 0.5 samples away from current time
              (abs(tdiff)<0.5) ||
              // It's possible for the frame timing to change between iterations
              // of the outer loop and because of this, it's possible that
              // we missed the best start. Start capturing anyways.
              ((tdiff>0)&&(tdiff<3))
            ) {
              // Configure parameters for this capture
              cl.filling=true;
              cl.pdu.late=tdiff;
              cl.buffer_offset=0;
              // Record the frequency offset and frame timing as they were
              // at the beginning of the capture.
              cl.pdu.frequency_offset=frequency_offset;
              cl.pdu.frame_timing=frame_timing;
            }
          }

          // Save this sample if our state indicates we are filling the
          // buffer.
          if (cl.filling) {
            cl.pdu.data(cl.buffer_offset++)=samples(t);
            if (cl.buffer_offset==128) {
              // Buffer is full! Send PDU
              {
                boost::mutex::scoped_lock lock2(tracked_cell.fifo_mutex);
                tracked_cell.fifo.push(cl.pdu);
                tracked_cell.fifo_peak_size=MAX(tracked_cell.fifo.size(),tracked_cell.fifo_peak_size);
                tracked_cell.fifo_condition.notify_one();
              }
              //cout << "fifo size: " << tracked_cell.fifo.size() << endl;
              // Calculate trigger parameters of next capture
              cl.filling=false;
              if (tracked_cell.cp_type==cp_type_t::EXTENDED) {
                cl.target_cap_start_time+=32+128;
              } else {
                cl.target_cap_start_time+=(cl.pdu.sym_num==6)?128+10:128+9;
              }
              cl.target_cap_start_time=mod(cl.target_cap_start_time,19200);
              slot_sym_inc(tracked_cell.n_symb_dl(),cl.pdu.slot_num,cl.pdu.sym_num);
            }
          }
        }
        ++it;
      }
    }
  }
}

