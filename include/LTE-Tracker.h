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
      const double & frame_timing
    ) : n_id_cell(n_id_cell), n_ports(n_ports), cp_type(cp_type), frame_timing(frame_timing) {
      fifo_peak_size=0;
      kill_me=false;
      sym_num=0;
      slot_num=0;
      target_cap_start_time=(cp_type==cp_type_t::NORMAL)?10:32;
      filling=0;
      buffer.set_size(128);
      buffer_offset=0;
      ac_fd.set_size(12);
      ac_fd=std::complex <double> (0,0);
      bulk_phase_offset=0;
      tracker_thread_ready=false;
      mib_decode_failures=0;
      crs_sp=itpp::vec(4);
      crs_sp=NAN;
      crs_np=itpp::vec(4);
      crs_np=NAN;
      crs_sp_av=itpp::vec(4);
      crs_sp_av=NAN;
      crs_np_av=itpp::vec(4);
      crs_np_av=NAN;
    }
    uint8 const n_symb_dl() const {
      return (cp_type==cp_type_t::NORMAL)?7:((cp_type==cp_type_t::EXTENDED)?6:-1);
    }
    boost::mutex mutex;
    boost::condition condition;
    boost::thread thread;
    // Indicates that the tracker process is ready to receive data.
    bool tracker_thread_ready;
    // These are not allowed to change
    const uint16 n_id_cell;
    const int8 n_ports;
    const cp_type_t::cp_type_t cp_type;
    // These are constantly changing
    double frame_timing;
    std::queue <td_fifo_pdu_t> fifo;
    uint32 fifo_peak_size;
    bool kill_me;
    // Calculated values returned by the cell tracker process.
    // Changing constantly.
    itpp::cvec ac_fd;
    // Tracker process uses these as local variables.
    double bulk_phase_offset;
    // The producer process (main) will use these members as
    // local variables...
    uint8 sym_num;
    uint8 slot_num;
    uint32 target_cap_start_time;
    double late;
    bool filling;
    itpp::cvec buffer;
    uint16 buffer_offset;
    double mib_decode_failures;
    itpp::vec crs_sp;
    itpp::vec crs_np;
    itpp::vec crs_sp_av;
    itpp::vec crs_np_av;
  private:
};
// Structure that is used to record all the tracked cells.
typedef struct {
  // List of cells which are currently being tracked.
  // Only the searcher can add elements to this list.
  // Only the main thread can remove elements from this list.
  boost::mutex mutex;
  std::list <tracked_cell_t *> tracked_cells;
} tracked_cell_list_t;
// Global data shared by all threads
typedef struct {
  // The frequency offset of the dongle. This value will be updated
  // continuously.
  boost::mutex frequency_offset_mutex;
  double frequency_offset;
  // This value will never change.
  double fc;
} global_thread_data_t;
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
  tracked_cell_list_t & tracked_cell_list
);

#endif

