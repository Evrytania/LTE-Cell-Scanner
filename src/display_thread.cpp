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

#define RED 1
#define GREEN 2
#define YELLOW 3
#define BLUE 4
#define MAGENTA 5
#define CYAN 6
#define WHITE 7

// Minimum terminal screensize required to run this program.
#define LINES_MIN 24
#define COLS_MIN 80

using namespace itpp;
using namespace std;

// Struct used to keep trach of which cell and which port has been displayed
// on a particular line.
typedef struct {
  bool occupied;
  // -1 indicates this row is is not one of the ports
  int16 n_id_cell;
  // -1 indicates this is the 'synchronization' port.
  int8 port_num;
} row_desc_t;

// ncurses right justified print of a string
void print_right(
  string s
) {
  int r,c;
  getyx(stdscr,r,c);

  move(r,COLS-s.size());
  printw("%s",s.c_str());
}

// ncurses left justified print of a string
void print_left(
  string s
) {
  int r,c;
  getyx(stdscr,r,c);

  move(r,0);
  printw("%s",s.c_str());
}

// ncurses center justified print of a string
void print_center(
  string s
) {
  int r,c;
  getyx(stdscr,r,c);

  c=(COLS-s.size())/2;
  c=MAX(c,0);
  move(r,c);
  printw("%s",s.c_str());
}

// Display status information for a certain cell on a certain line on the
// screen.
void display_cell(
  tracked_cell_t & tracked_cell,
  uint8 row,
  const bool & fifo_status,
  const bool & avg_values,
  const bool & expert_mode
) {
  move(row,0);
  attron(COLOR_PAIR(BLUE));
  printw("[Cell ID %3i][TO: %7.1lf]",
    tracked_cell.n_id_cell,
    tracked_cell.frame_timing()
  );
  if (expert_mode) {
    printw("[UOS pwr: %5.1lf dB]",db10(tracked_cell.sync_np_blank_av));
  }
  const double health=100-tracked_cell.mib_decode_failures/CELL_DROP_THRESHOLD*100;
  attroff(COLOR_PAIR(BLUE));
  uint32 attrs=0;
  if (health<25) {
    attrs=COLOR_PAIR(RED)|A_REVERSE;
  } else if (health<75) {
    attrs=COLOR_PAIR(YELLOW)|A_REVERSE;
  } else {
    attrs=COLOR_PAIR(BLUE);
  }
  attron(attrs);
  printw("[Health:");
  printw("%3.0lf",health);
  printw("]");
  attroff(attrs);
  attron(COLOR_PAIR(BLUE));
  if (fifo_status) {
    printw("[buffer %5i/%5i]",
      tracked_cell.fifo.size(),
      tracked_cell.fifo_peak_size
    );
  }
  printw("\n");
  attroff(COLOR_PAIR(BLUE));

  attron(COLOR_PAIR(YELLOW));
  for (uint8 t=0;t<tracked_cell.n_ports;t++) {
    printw("  P%i ",t);
    if (expert_mode) {
      if (avg_values) {
        printw("%5.1lf/%5.1lf/%5.1lf",
          db10(tracked_cell.crs_sp_raw_av(t)),
          db10(tracked_cell.crs_np_av(t)),
          db10(tracked_cell.crs_sp_raw_av(t)/tracked_cell.crs_np_av(t))
        );
      } else {
        printw("%5.1lf/%5.1lf/%5.1lf",
          db10(tracked_cell.crs_sp_raw(t)),
          db10(tracked_cell.crs_np(t)),
          db10(tracked_cell.crs_sp_raw(t)/tracked_cell.crs_np(t))
        );
      }
      int8 CB=-1;
      for (uint8 k=1;k<12;k++) {
        if (abs(tracked_cell.ac_fd(k))<=0.5) {
          CB=k;
          break;
        }
      }
      if (CB==-1) {
        printw(" >990 kHz\n");
      } else {
        printw(" %4i kHz\n",CB*90);
      }
    } else {
      if (isfinite(db10(tracked_cell.crs_sp_raw_av(t)/tracked_cell.crs_np_av(t)))) {
        printw("%5.1lf dB SNR\n",
          db10(tracked_cell.crs_sp_raw_av(t)/tracked_cell.crs_np_av(t))
        );
      } else {
        printw(" -Inf dB SNR\n");
      }
    }
  }
  if (expert_mode) {
    if (avg_values) {
      printw("  S  %5.1lf/%5.1lf/%5.1lf\n",
        db10(tracked_cell.sync_sp_av),
        db10(tracked_cell.sync_np_av),
        db10(tracked_cell.sync_sp_av/tracked_cell.sync_np_av)
      );
    } else {
      printw("  S  %5.1lf/%5.1lf/%5.1lf\n",
        db10(tracked_cell.sync_sp),
        db10(tracked_cell.sync_np),
        db10(tracked_cell.sync_sp/tracked_cell.sync_np)
      );
    }
  } else {
    if (isfinite(db10(tracked_cell.sync_sp_av/tracked_cell.sync_np_av))) {
      printw("  S  %5.1lf dB SNR\n",
        db10(tracked_cell.sync_sp_av/tracked_cell.sync_np_av)
      );
    } else {
      printw("  S   -Inf dB SNR\n");
    }
  }
  attroff(COLOR_PAIR(YELLOW));
}

// Will a certain cell display fit starting at a certain row?
bool will_fit(
  const vector <row_desc_t> & row_desc,
  const uint16 & desired_row,
  const uint8 & n_rows_required
) {
  for (uint16 t=desired_row;t<desired_row+n_rows_required;t++) {
    if (row_desc[t].occupied)
      return false;
  }
  return true;
}

// Indicate that certain rows are occupied
void set_occupied(
  const tracked_cell_t & tracked_cell,
  vector <row_desc_t> & row_desc,
  const uint16 & print_row
) {
  row_desc[print_row].occupied=true;
  for (uint8 t=0;t<tracked_cell.n_ports;t++) {
    row_desc[print_row+t+1].occupied=true;
    row_desc[print_row+t+1].n_id_cell=tracked_cell.n_id_cell;
    row_desc[print_row+t+1].port_num=t;
  }
  row_desc[print_row+tracked_cell.n_ports+1].occupied=true;
  row_desc[print_row+tracked_cell.n_ports+1].n_id_cell=tracked_cell.n_id_cell;
  row_desc[print_row+tracked_cell.n_ports+1].port_num=-1;
}

// Helper function to plot a transfer function on the screen.
void plot_trace(
  // Trace to plot
  const vec & Y,
  const vec & X,
  // Range of trace that is to be plotted
  const double & x_min,
  const double & x_max,
  const double & x_tick,
  const double & x_supermark,
  const double & y_min,
  const double & y_max,
  const double & y_tick,
  // Coordinate of the upper left corner of the plot.
  const uint16 & ul_row,
  const uint16 & ul_col,
  // Coordinate of the lower right corner of the plot.
  const uint16 & lr_row,
  const uint16 & lr_col,
  // Whether to connect the dots in the y dimension.
  bool connect_the_dots
) {
  ASSERT(ul_row<lr_row);
  ASSERT(ul_col<lr_col);

  // This is the area including the axes
  const uint16 canvas_width=lr_col-ul_col+1;
  const uint16 canvas_height=lr_row-ul_row+1;
  // This is the actual plot area, not including axes
  const uint16 plot_width=canvas_width-5;
  const uint16 plot_height=canvas_height-1;
  const uint16 plot_ul_row=ul_row;
  const uint16 plot_ul_col=ul_col+5;
  const uint16 plot_lr_row=lr_row-1;
  const uint16 plot_lr_col=lr_col;

  // Plot the border
  attron(COLOR_PAIR(BLUE));
  for (uint16 t=plot_ul_row;t<=plot_lr_row;t++) {
    move(t,plot_ul_col-1);
    //addch('|');
    addch(ACS_VLINE);
  }
  move(plot_lr_row+1,plot_ul_col-1);
  //addch('+');
  addch(ACS_LLCORNER);
  for (uint16 t=plot_ul_col;t<=plot_lr_col;t++) {
    move(plot_lr_row+1,t);
    //addch('-');
    addch(ACS_HLINE);
  }

  // Add tick marks
  // Y axis
  const double first_y_tick=ceil(y_min/y_tick)*y_tick;
  vec y_ticks=itpp_ext::matlab_range(first_y_tick,y_tick,y_max);
  for (uint32 t=0;t<(unsigned)y_ticks.size();t++) {
    const uint16 tick_row=round_i(plot_lr_row-(y_ticks(t)-y_min)/((y_max-y_min)/(plot_height-1)));
    move(tick_row,plot_ul_col-1);
    //addch('-');
    addch(ACS_RTEE);
    move(tick_row,plot_ul_col-5);
    stringstream ss;
    ss << setw(4) << right << y_ticks(t) << endl;
    addnstr(ss.str().c_str(),4);
  }
  // X axis
  const double first_x_tick=ceil(x_min/x_tick)*x_tick;
  vec x_ticks=itpp_ext::matlab_range(first_x_tick,x_tick,x_max);
  for (uint32 t=0;t<(unsigned)x_ticks.size();t++) {
    const uint16 tick_col=round_i(plot_ul_col+(x_ticks(t)-x_min)/((x_max-x_min)/(plot_width-1)));
    move(plot_lr_row+1,tick_col);
    //addch('|');
    addch(ACS_BTEE);
  }
  // Add the supermark
  if (isfinite(x_supermark)) {
    const int16 sm_col=round_i((x_supermark-x_min)/((x_max-x_min)/(plot_width-1)));
    if ((sm_col>=0)&&(sm_col<=plot_width-1)) {
      move(plot_lr_row+1,sm_col+plot_ul_col);
      addch('*');
    }
  }
  attroff(COLOR_PAIR(BLUE));

  // Plot trace
  attron(COLOR_PAIR(YELLOW));
  vec x=linspace(x_min,x_max,plot_width);
  vec y=interp1(X,Y,x);
  int16 prev_row=-1;
  for (uint16 t=0;t<plot_width;t++) {
    if (!isfinite(y(t))) {
      prev_row=-1;
      continue;
    }
    int16 curr_row=round_i(plot_lr_row-(y(t)-y_min)/((y_max-y_min)/(plot_height-1)));
    char ch='*';
    if (curr_row>plot_lr_row) {
      //ch=ACS_DARROW;
      ch='-';
      curr_row=plot_lr_row;
    }
    if (curr_row<plot_ul_row) {
      //ch=ACS_UARROW;
      ch='^';
      curr_row=plot_ul_row;
    }
    move(curr_row,plot_ul_col+t);
    addch(ch);
    if (connect_the_dots&&(prev_row!=-1)&&(abs(curr_row-prev_row)>1)) {
      const int8 dir=(curr_row>prev_row)?1:-1;
      const uint16 mid=round_i((prev_row+curr_row)/2);
      for (uint16 k=prev_row+dir;k!=curr_row;k+=dir) {
        uint16 col;
        if (dir==1)
          col=(k<mid)?(t-1):t;
        else
          col=(k<mid)?t:t-1;
        move(k,plot_ul_col+col);
        addch(ACS_BULLET);
      }
    }
    prev_row=curr_row;
  }
  attron(COLOR_PAIR(YELLOW));

}

// Process that displays the status of all the tracker threads.
enum disp_mode_t {STD, DETAIL};
void display_thread(
  sampbuf_sync_t & sampbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list,
  bool & expert_mode
) {
  global_thread_data.display_thread_id=syscall(SYS_gettid);

  // Initialize the curses screen
  initscr();
  start_color();
  use_default_colors();
  if (LINES<LINES_MIN) {
    //endwin();
    cerr << "Error: resize window so it has at least " << LINES_MIN << " rows and " << COLS_MIN << " columns" << endl;
    ABORT(-1);
  }
  if (COLS<COLS_MIN) {
    //endwin();
    cerr << "Error: resize window so it has at least " << LINES_MIN << " rows and " << COLS_MIN << " columns" << endl;
    ABORT(-1);
  }
  // Do not echo input chars to screen.
  noecho();
  // Turn off display of cursor.
  curs_set(0);
  // Enable function keys, arrow keys, etc.
  keypad(stdscr,true);
  //init_color(COLOR_BLACK,0,0,0);
  //init_pair(3,COLOR_GREEN,-1);
  init_pair(RED,COLOR_RED,-1);
  init_pair(GREEN,COLOR_GREEN,-1);
  init_pair(YELLOW,COLOR_YELLOW,-1);
  init_pair(BLUE,COLOR_BLUE,-1);
  init_pair(MAGENTA,COLOR_MAGENTA,-1);
  init_pair(CYAN,COLOR_CYAN,-1);
  init_pair(WHITE,COLOR_WHITE,-1);
  // Reduce delay between pressing 'Esc' and having this program
  // respond.
  // In ms.
  set_escdelay(50);

  // Start, end, and size of cell info display area.
  const uint16 CELL_DISP_START_ROW=2;
  const uint16 CELL_DISP_END_ROW=LINES-6;
  const uint16 CELL_DISP_N_ROWS=CELL_DISP_END_ROW-CELL_DISP_START_ROW+1;

  // Record where a particular cell was most recently printed.
  ivec disp_history(504);
  disp_history=-1;
  vector <row_desc_t> row_desc(CELL_DISP_N_ROWS);
  //bvec row_occupied(CELL_DISP_N_ROWS);

  // Settings that the user can change in realtime
  bool auto_refresh=true;
  double refresh_delay_sec=1;
  bool fifo_status=false;
  bool avg_values=true;
  disp_mode_t disp_mode=STD;
  int8 detail_type=0;
#define N_DETAILS 2
  // Send control chars directly to program.
  //cbreak();
  halfdelay(round_i(refresh_delay_sec*10.0));
  int16 highlight_row=-1;
  while (true) {
    clear();
    //attron(COLOR_PAIR(3));

    // No matter which display mode we are in, always update row_disc
    // and highlight_row.
    bool all_cells_displayed=true;
    for (uint16 t=0;t<CELL_DISP_N_ROWS;t++) {
      row_desc[t].occupied=false;
      row_desc[t].n_id_cell=-1;
      row_desc[t].port_num=-2;
      //row_occupied(t)=false;
    }
    move(2,0);
    {
      boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
      list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
      vector <tracked_cell_t *> pass2_display;
      while (it!=tracked_cell_list.tracked_cells.end()) {
        tracked_cell_t & tracked_cell=(*(*it));

        // Deadlock possible???
        boost::mutex::scoped_lock lock1(tracked_cell.fifo_mutex);
        boost::mutex::scoped_lock lock2(tracked_cell.meas_mutex);

        uint8 n_rows_required=tracked_cell.n_ports+2;

        // If this cell has been displayed before, try to display it
        // in the same location.
        if (disp_history(tracked_cell.n_id_cell)!=-1) {
          uint16 row_desired=disp_history(tracked_cell.n_id_cell);
          if (will_fit(row_desc,row_desired,n_rows_required-1)) {
            if (disp_mode==STD)
              display_cell(tracked_cell,row_desired+CELL_DISP_START_ROW,fifo_status,avg_values,expert_mode);
            set_occupied(tracked_cell,row_desc,row_desired);
          } else {
            pass2_display.push_back(&tracked_cell);
          }
        } else {
          pass2_display.push_back(&tracked_cell);
        }
        ++it;
      }

      // Display cells that cannot be displayed where they have previously
      // been displayed.
      for (uint8 t=0;t<pass2_display.size();t++) {
        tracked_cell_t & tracked_cell=(*pass2_display[t]);

        // Deadlock possible???
        boost::mutex::scoped_lock lock1(tracked_cell.fifo_mutex);
        boost::mutex::scoped_lock lock2(tracked_cell.meas_mutex);

        uint8 n_rows_required=tracked_cell.n_ports+2;
        bool placed=false;
        for (uint16 k=0;k<CELL_DISP_N_ROWS-n_rows_required+1;k++) {
          if (will_fit(row_desc,k,n_rows_required)) {
            if (disp_mode==STD)
              display_cell(tracked_cell,k+CELL_DISP_START_ROW,fifo_status,avg_values,expert_mode);
            set_occupied(tracked_cell,row_desc,k);
            disp_history(tracked_cell.n_id_cell)=k;
            placed=true;
            break;
          }
        }
        if (!placed)
          all_cells_displayed=false;
      }
    }

    // Always highlight a row, if a cell exists.
    if (highlight_row==-1) {
      for (uint16 t=0;t<CELL_DISP_N_ROWS;t++) {
        if (row_desc[t].n_id_cell!=-1) {
          highlight_row=t;
          break;
        }
      }
    }

    if (disp_mode==STD) {
      // Header and footer
      {
        stringstream ss;
        ss << "LTE-Tracker v" << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_LEVEL << " -- www.evrytania.com";
        move(0,0);
        attron(COLOR_PAIR(CYAN));
        print_center(ss.str());
        attroff(COLOR_PAIR(CYAN));
      }

      move(CELL_DISP_END_ROW+1,0);
      if (!all_cells_displayed) {
        attron(COLOR_PAIR(YELLOW));
        printw("Warning: some tracked cells could not be displayed! (screen has too few rows)\n");
        attroff(COLOR_PAIR(YELLOW));
      } else {
        printw("\n");
      }
      if (global_thread_data.cell_seconds_dropped()||global_thread_data.raw_seconds_dropped()) {
        attron(COLOR_PAIR(RED));
        printw("[dropped cell/raw data: %i/%i s]\n",global_thread_data.cell_seconds_dropped(),global_thread_data.raw_seconds_dropped());
        attroff(COLOR_PAIR(RED));
      } else {
        printw("\n");
      }
      if (expert_mode) {
        if (avg_values) {
          printw("Showing average measurements\n");
        } else {
          printw("Showing instant measurements\n");
        }
      } else {
        printw("\n");
      }
      attron(COLOR_PAIR(MAGENTA));
      printw("[FO: %6.0lf Hz]",global_thread_data.frequency_offset());
      //attron(A_BOLD);
      printw("[searcher delay: %.1lf s]",global_thread_data.searcher_cycle_time());
      //attroff(A_BOLD);

      if (fifo_status) {
        boost::mutex::scoped_lock lock(sampbuf_sync.mutex);
        stringstream ss;
        uint8 w=ceil(log10(sampbuf_sync.fifo_peak_size));
        ss << "[inp buf: " << setw(w) << sampbuf_sync.fifo.size() << "/" << sampbuf_sync.fifo_peak_size << "]";
        //attron(A_BOLD);
        printw("%s\n",ss.str().c_str());
        //attroff(A_BOLD);
      } else {
        printw("\n");
      }
      attroff(COLOR_PAIR(MAGENTA));

      /*
      printw("Searcher thread ID: %5i\n",global_thread_data.searcher_thread_id);
      printw("Producer thread ID: %5i\n",global_thread_data.producer_thread_id);
      printw("Main     thread ID: %5i\n",global_thread_data.main_thread_id);
      printw("Display  thread ID: %5i\n\n",global_thread_data.display_thread_id);
      */

      attron(COLOR_PAIR(GREEN));
      //printw("Up/down = select, left/right = details, (h)elp, (q)uit\n");
      print_right("right arrow > TF");
      print_center("(h)elp, (q)uit");
      attroff(COLOR_PAIR(GREEN));
      //attroff(A_BOLD);

      // Highlight one of the rows
      if (highlight_row!=-1) {
        if (row_desc[highlight_row].n_id_cell==-1) {
          highlight_row=-1;
        } else {
          move(CELL_DISP_START_ROW+highlight_row,0);
          chgat(-1,A_REVERSE,0,NULL);
        }
      }

    } else {
      // Zoom into port details

      // Shortcuts
      const int16 & n_id_cell=row_desc[highlight_row].n_id_cell;
      const int8 & port_num=row_desc[highlight_row].port_num;

      {
        boost::mutex::scoped_lock lock(tracked_cell_list.mutex);
        list <tracked_cell_t *>::iterator it=tracked_cell_list.tracked_cells.begin();
        vector <tracked_cell_t *> pass2_display;
        bool cell_found=false;
        while (it!=tracked_cell_list.tracked_cells.end()) {
          tracked_cell_t & tracked_cell=(*(*it));
          if ((tracked_cell.n_id_cell==n_id_cell)&&((port_num==-1)||(port_num<tracked_cell.n_ports))) {
            cell_found=true;
            break;
          }
          ++it;
        }
        if (cell_found) {
          tracked_cell_t & tracked_cell=(*(*it));
          boost::mutex::scoped_lock lock2(tracked_cell.meas_mutex);
          if (detail_type==0) {
            // Plot transfer function mag
            vec trace;
            if (port_num==-1) {
              trace=db10(sqr(tracked_cell.sync_ce));
            } else {
              trace=db10(sqr(tracked_cell.ce.get_row(port_num)));
            }
            plot_trace(
              // Trace desc.
              trace,itpp_ext::matlab_range(0.0,71.0),
              // X axis
              0,71,12,NAN,
              // Y axis
              -50,0,10,
              // UL corner
              1,0,
              // LR corner
              LINES-2,0+72+4,
              true
            );
            move(0,0);
            stringstream ss;
            ss << "Cell " << n_id_cell;
            if (port_num==-1) {
              ss << " Sync channel magnitude\n";
            } else {
              ss << " port " << port_num << " magnitude\n";
            }
            attron(COLOR_PAIR(CYAN));
            print_center(ss.str());
            attroff(COLOR_PAIR(CYAN));
            move(LINES-1,0);
            attron(COLOR_PAIR(GREEN));
            print_center("(h)elp, (q)uit");
            print_right("right arrow > phase resp");
            print_left("top menu < left arrow");
            attroff(COLOR_PAIR(GREEN));
          } else if (detail_type==1) {
            // Plot transfer function phase
            vec trace;
            double mean_ang;
            if (port_num==-1) {
              trace=arg(tracked_cell.sync_ce);
              mean_ang=arg(sum(exp(J*trace(5,66))));
              trace=arg(tracked_cell.sync_ce*exp(J*-mean_ang));
              trace(0)=NAN;
              trace(1)=NAN;
              trace(2)=NAN;
              trace(3)=NAN;
              trace(4)=NAN;
              trace(67)=NAN;
              trace(68)=NAN;
              trace(69)=NAN;
              trace(70)=NAN;
              trace(71)=NAN;
            } else {
              trace=arg(tracked_cell.ce.get_row(port_num));
              mean_ang=arg(sum(exp(J*trace)));
              trace=arg(tracked_cell.ce.get_row(port_num)*exp(J*-mean_ang));
            }
            trace=trace/pi*180;
            plot_trace(
              // Trace desc.
              trace,itpp_ext::matlab_range(0.0,71.0),
              // X axis
              0,71,12,(mean_ang+pi)/(2*pi)*71,
              // Y axis
              -40,40,10,
              // UL corner
              1,0,
              // LR corner
              LINES-2,0+72+4,
              false
            );
            move(0,0);
            stringstream ss;
            ss << "Cell " << n_id_cell;
            if (port_num==-1) {
              ss << " Sync channel phase\n";
            } else {
              ss << " port " << port_num << " phase\n";
            }
            attron(COLOR_PAIR(CYAN));
            print_center(ss.str());
            attroff(COLOR_PAIR(CYAN));
            move(LINES-1,0);
            attron(COLOR_PAIR(GREEN));
            print_center("(h)elp, (q)uit");
            //print_right("right arrow > phase response");
            print_left("mag resp < left arrow");
            attroff(COLOR_PAIR(GREEN));
          } else if (detail_type==2) {
            // Frequency domain autocorrelation
            const vec trace=abs(tracked_cell.ac_fd);
            plot_trace(
              // Trace desc.
              trace,itpp_ext::matlab_range(0.0,11.0),
              // X axis
              0,11,2,NAN,
              // Y axis
              0,1.2,.5,
              // UL corner
              0,0,
              // LR corner
              LINES-3,0+72+4,
              true
            );
            move(LINES-2,0);
            printw("Cell ID: %i\n",n_id_cell);
            printw("Frequency domain channel autocorrelation function. x-axis spans 1.26MHz\n");
          } else if (detail_type==3) {
            // Time domain autocorrelation
            const vec trace=abs(tracked_cell.ac_td);
            plot_trace(
              // Trace desc.
              trace,itpp_ext::matlab_range(0.0,71.0)*.0005,
              // X axis
              0,71*.0005,.010,NAN,
              // Y axis
              0,3.2,.5,
              // UL corner
              0,0,
              // LR corner
              LINES-3,0+72+4,
              true
            );
            move(LINES-2,0);
            printw("Cell ID: %i\n",n_id_cell);
            printw("Time domain channel autocorrelation function. x-axis spans 35.5ms\n");
          }
        } else {
          move(1,0);
          printw("Cell is no longer being tracked. Press left arrow to go back!\n");
        }
      }

    }

    refresh();

    // Handle keyboard input.
    // Previous halfdelay() function ensures that this will not block.
    int ch=getch();
    switch (ch) {
      case 'q':
      case 'Q':
        //endwin();
        ABORT(-1);
        break;
      case 'r':
      case 'R':
        auto_refresh=!auto_refresh;
        if (auto_refresh) {
          halfdelay(1+round_i(refresh_delay_sec*10.0));
        } else {
          cbreak();
        }
        break;
      case '-':
      case '_':
        refresh_delay_sec=MIN(15,refresh_delay_sec*1.5);
        halfdelay(round_i(refresh_delay_sec*10.0));
        break;
      case '+':
      case '=':
        refresh_delay_sec=MAX(0.001,refresh_delay_sec/1.5);
        halfdelay(round_i(refresh_delay_sec*10.0));
        break;
      case 'f':
      case 'F':
        fifo_status=!fifo_status;
        break;
      case 'a':
      case 'A':
        avg_values=!avg_values;
        break;
      case 27:
        // Escape key
        disp_mode=STD;
        break;
      case 'k':
      case 'K':
      case KEY_UP:
        for (int16 t=highlight_row-1;t>=0;t--) {
          if (row_desc[t].n_id_cell!=-1) {
            highlight_row=t;
            break;
          }
        }
        break;
      case 'j':
      case 'J':
      case KEY_DOWN:
        for (uint16 t=highlight_row+1;t<CELL_DISP_N_ROWS;t++) {
          if (row_desc[t].n_id_cell!=-1) {
            highlight_row=t;
            break;
          }
        }
        break;
      case 'l':
      case 'L':
      case KEY_RIGHT:
      case '\n':
        if (disp_mode==STD) {
          disp_mode=DETAIL;
          detail_type=0;
        } else {
          detail_type=MIN(detail_type+1,N_DETAILS-1);
        }
        break;
      case KEY_LEFT:
        if (disp_mode==DETAIL) {
          if (detail_type==0) {
            disp_mode=STD;
          } else {
            detail_type-=1;
          }
        }
        break;
      case 'h':
      case 'H':
        clear();
        move(0,0);
        if (disp_mode==STD) {
          print_center("LTE-Tracker main display help menu\n");
          printw("\n");
          //printw("Display frequency response:\n");
          printw("up/down keys move selection bar\n");
          printw("'Enter' or right-arrow displays the frequency response\n");
          printw("\n");
          printw("P0 : SNR received from port 0\n");
          printw("S  : SNR of the PSS/SSS synchronization channel\n");
          printw("TO : frame timing referenced to the dongle's timescale of 19200 samples\n");
          printw("FO : dongle's current residual frequency offset\n");
          printw("\n");
          printw("'q' exits the program\n");
          printw("+/- increases or decreases screen refresh rate\n");
          printw("'r' turns automatic screen refresh on/off\n");
          printw("'f' displays the current status of the fifo's\n");
          printw("Esc always returns to the top menu\n");
          printw("\n");
          printw("Health is a measure of the amount of time since the last time the MIB\n");
          printw("was successfully decoded for that cell. Cells are dropped when health\n");
          printw("reaches 0.\n");
          printw("\n");
          printw("searcher delay is the amount of time (in seconds) that it takes for the\n");
          printw("searcher to complete a full search cycle\n");
          printw("\n");
          printw("Press any key to exit help screen!\n");
          cbreak();
          getch();
          if (auto_refresh) {
            halfdelay(round_i(refresh_delay_sec*10.0));
          }
        } else {
          print_center("LTE-Tracker mag/phase response help menu\n");
          printw("\n");
          printw("up/down keys change port and cell ID\n");
          printw("left/right keys change to phase response and main menu\n");
          printw("\n");
          printw("'q' exits the program\n");
          printw("+/- increases or decreases screen refresh rate\n");
          printw("'r' turns automatic screen refresh on/off\n");
          printw("Esc always returns to the top menu\n");
          printw("\n");
          printw("For the phase plot, the average phase offset is removed from the trace, but\n");
          printw("the value of this phase offset is indicated by an asterisk on the x axis.\n");
          printw("\n");
          printw("Press any key to exit help screen!\n");
          cbreak();
          getch();
          if (auto_refresh) {
            halfdelay(round_i(refresh_delay_sec*10.0));
          }
        }
        break;
    }
  }
}

