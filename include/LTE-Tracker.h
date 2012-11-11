#ifndef HAVE_LTE_TRACKER_H
#define HAVE_LTE_TRACKER_H

//
// Data structures used to communicate between threads.
//
// A packet of information that is sent from the main thread to each
// tracker thread.
typedef struct {
  itpp::cvec data;
  uint8 slot_num;
  uint8 sym_num;
  double late;
  double frequency_offset;
  double frame_timing;
} td_fifo_pdu_t;

// Structure to describe a cell which is currently being tracked.
class tracked_cell_t {
  public:
    // Initializer
    tracked_cell_t(
      const uint16 & n_id_cell,
      const int8 & n_ports,
      const cp_type_t::cp_type_t & cp_type,
      const int8 & n_rb_dl,
      const phich_duration_t::phich_duration_t & phich_duration,
      const phich_resource_t::phich_resource_t & phich_resource,
      const double & ft,
      const uint32 & serial_num
    ) :
      n_id_1(floor(n_id_cell/3.0)),
      n_id_2(n_id_cell-3*floor(n_id_cell/3.0)),
      n_id_cell(n_id_cell),
      n_ports(n_ports),
      cp_type(cp_type),
      n_rb_dl(n_rb_dl),
      phich_duration(phich_duration),
      phich_resource(phich_resource),
      serial_num(serial_num)
    {
      frame_timing_private=ft;
      fifo_peak_size=0;
      kill_me=false;
      ac_fd.set_size(12);
      ac_fd=std::complex <double> (0,0);
      ac_td.set_size(72);
      ac_td=std::complex <double> (0,0);
      tracker_thread_ready=false;
      mib_decode_failures=0;
      crs_tp=itpp::vec(4);
      crs_sp_raw=itpp::vec(4);
      crs_sp_raw=NAN;
      crs_np=itpp::vec(4);
      crs_np=NAN;
      crs_tp_av=itpp::vec(4);
      crs_sp_raw_av=itpp::vec(4);
      crs_sp_raw_av=NAN;
      crs_np_av=itpp::vec(4);
      crs_np_av=NAN;
      sync_tp=NAN;
      sync_sp=NAN;
      sync_np=NAN;
      sync_np_blank=NAN;
      sync_tp_av=NAN;
      sync_sp_av=NAN;
      sync_np_av=NAN;
      sync_np_blank_av=NAN;
      launched=false;
    }
    inline uint8 const n_symb_dl() const {
      return (cp_type==cp_type_t::NORMAL)?7:((cp_type==cp_type_t::EXTENDED)?6:-1);
    }
    // Constants that do not change and can be read freely.
    const uint8 n_id_1;
    const uint8 n_id_2;
    const uint16 n_id_cell;
    const int8 n_ports;
    const cp_type_t::cp_type_t cp_type;
    const int8 n_rb_dl;
    const phich_duration_t::phich_duration_t phich_duration;
    const phich_resource_t::phich_resource_t phich_resource;
    const uint32 serial_num;

    // Do we need this?
    boost::thread thread;

    // Mutex and data structures for the flow of information from the
    // producer thread to the tracker thread.
    boost::mutex fifo_mutex;
    boost::condition fifo_condition;
    std::queue <td_fifo_pdu_t> fifo;
    uint32 fifo_peak_size;

    // Indicates that the tracker process is ready to receive data.
    bool tracker_thread_ready;
    // Indicates that the thread has been launched
    bool launched;

    // Mutex and measurement data produced by the tracker thread and read by
    // the display thread.
    boost::mutex meas_mutex;
    double mib_decode_failures;
    itpp::vec crs_tp;
    itpp::vec crs_sp_raw;
    itpp::vec crs_np;
    itpp::vec crs_tp_av;
    itpp::vec crs_sp_raw_av;
    itpp::vec crs_np_av;
    itpp::cmat ce;
    double sync_tp;
    double sync_sp;
    double sync_np;
    double sync_np_blank;
    itpp::cvec sync_ce;
    double sync_tp_av;
    double sync_sp_av;
    double sync_np_av;
    double sync_np_blank_av;
    // Frequency domain channel autocorrelation.
    itpp::cvec ac_fd;
    // Time domain channel autocorrelation.
    itpp::cvec ac_td;

    // Read/write frame_timing (via mutex).
    // Only one thread (tracker_thread) can update frame_timing. We
    // only need to ensure that no thread reads a partial value before
    // the new value is completely written.
    inline double frame_timing() {
      boost::mutex::scoped_lock lock(frame_timing_mutex);
      double r=frame_timing_private;
      return r;
    }
    inline void frame_timing(const double & ft) {
      boost::mutex::scoped_lock lock(frame_timing_mutex);
      frame_timing_private=ft;
    }

    bool kill_me;

  private:
    // Frame timing info is needed by producer, tracker, and display
    // threads.
    boost::mutex frame_timing_mutex;
    double frame_timing_private;
};

// Structure that stores the list of all the tracked cells.
typedef struct {
  // List of cells which are currently being tracked.
  // Only the searcher can add elements to this list.
  // Only the main thread can remove elements from this list.
  boost::mutex mutex;
  std::list <tracked_cell_t *> tracked_cells;
} tracked_cell_list_t;

// Global data shared by all threads
class global_thread_data_t {
  public:
    // Constructor
    global_thread_data_t(
      const double & fc_requested,
      const double & fc_programmed,
      const double & fs_programmed
    ) :
      fc_requested(fc_requested),
      fc_programmed(fc_programmed),
      fs_programmed(fs_programmed)
    {
      searcher_cycle_time_private=0;
      cell_seconds_dropped_private=0;
      raw_seconds_dropped_private=0;
    }
    // These values will never change.
    const double fc_requested;
    const double fc_programmed;
    const double fs_programmed;
    // Read/write frequency offset (via mutex).
    // Mutex makes sure that no read or write is interrupted when
    // only part of the data has been read.
    inline double frequency_offset() {
      boost::mutex::scoped_lock lock(frequency_offset_mutex);
      double r=frequency_offset_private;
      return r;
    }
    inline void frequency_offset(const double & f) {
      boost::mutex::scoped_lock lock(frequency_offset_mutex);
      frequency_offset_private=f;
    }
    // Read/write searcher cycle time (via mutex).
    // Mutex makes sure that no read or write is interrupted when
    // only part of the data has been read.
    inline double searcher_cycle_time() {
      boost::mutex::scoped_lock lock(searcher_cycle_time_mutex);
      double r=searcher_cycle_time_private;
      return r;
    }
    inline void searcher_cycle_time(const double & t) {
      boost::mutex::scoped_lock lock(searcher_cycle_time_mutex);
      searcher_cycle_time_private=t;
    }
    inline uint32 cell_seconds_dropped() {
      boost::mutex::scoped_lock lock(cell_seconds_dropped_mutex);
      double r=cell_seconds_dropped_private;
      return r;
    }
    inline void cell_seconds_dropped_inc() {
      boost::mutex::scoped_lock lock(cell_seconds_dropped_mutex);
      cell_seconds_dropped_private+=1;
    }
    inline uint32 raw_seconds_dropped() {
      boost::mutex::scoped_lock lock(raw_seconds_dropped_mutex);
      double r=raw_seconds_dropped_private;
      return r;
    }
    inline void raw_seconds_dropped_inc() {
      boost::mutex::scoped_lock lock(raw_seconds_dropped_mutex);
      raw_seconds_dropped_private+=1;
    }
    uint32 searcher_thread_id;
    uint32 producer_thread_id;
    uint32 main_thread_id;
    uint32 display_thread_id;
  private:
    // The frequency offset of the dongle. This value will be updated
    // continuously.
    boost::mutex frequency_offset_mutex;
    double frequency_offset_private;
    boost::mutex searcher_cycle_time_mutex;
    double searcher_cycle_time_private;
    boost::mutex cell_seconds_dropped_mutex;
    uint32 cell_seconds_dropped_private;
    boost::mutex raw_seconds_dropped_mutex;
    uint32 raw_seconds_dropped_private;
};

// IPC between main thread and searcher thread covering data capture issues.
typedef struct {
  boost::mutex mutex;
  boost::condition condition;
  bool request;
  itpp::cvec capbuf;
  double late;
} capbuf_sync_t;

// IPC between main thread and producer thread.
typedef struct {
  boost::mutex mutex;
  boost::condition condition;
  std::deque <uint8> fifo;
  uint32 fifo_peak_size;
} sampbuf_sync_t;

// Small helper function to increment the slot number and the symbol number.
inline void slot_sym_inc(
  const uint8 n_symb_dl,
  uint8 & slot_num,
  uint8 & sym_num
) {
  sym_num=itpp::mod(sym_num+1,n_symb_dl);
  if (sym_num==0)
    slot_num=itpp::mod(slot_num+1,20);
}

// Prototypes for all the threads.
void producer_thread(
  sampbuf_sync_t & sampbuf_sync,
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list,
  double & fc
);
void tracker_thread(
  tracked_cell_t & tracked_cell,
  global_thread_data_t & global_thread_data
);
void searcher_thread(
  capbuf_sync_t & capbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list
);
void display_thread(
  sampbuf_sync_t & sampbuf_sync,
  global_thread_data_t & global_thread_data,
  tracked_cell_list_t & tracked_cell_list,
  bool & expert_mode
);

#endif

