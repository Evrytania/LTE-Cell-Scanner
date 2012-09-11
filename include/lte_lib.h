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

#ifndef HAVE_LTE_LIB_H
#define HAVE_LTE_LIB_H

// LTE PN generator
itpp::bvec lte_pn(
  const uint32 & c_init,
  const uint32 & len
);

// PSS in the time and frequency domain
class PSS_fd {
  public:
    // Constructor
    PSS_fd(void);
    // Return vector
    const itpp::cvec & operator[](const uint8 & idx) const;
  private:
    std::vector <itpp::cvec> table;
};
class PSS_td {
  public:
    // Constructor
    PSS_td(void);
    // Return vector
    const itpp::cvec & operator[](const uint8 & idx) const;
  private:
    std::vector <itpp::cvec> table;
};

// SSS is only needed in the frequency domain
class SSS_fd {
  public:
    // Constructor
    SSS_fd(void);
    // Return vector
    const itpp::ivec & operator()(const uint8 & n_id_1,const uint8 & n_id_2,const uint8 & n_slot) const;
  private:
    std::vector < std::vector < std::vector <itpp::ivec> > > table;
};
class SSS_td {
  public:
    // Constructor
    SSS_td(void);
    // Return vector
    const itpp::cvec & operator()(const uint8 & n_id_1,const uint8 & n_id_2,const uint8 & n_slot) const;
  private:
    static std::vector < std::vector < std::vector <itpp::cvec> > > table;
};

// Class used to precompute all DL RS.
class RS_DL {
  public:
    RS_DL(
      const uint16 & n_id_cell,
      const uint8 & n_rb_dl,
      const cp_type_t::cp_type_t & cp_type
    );
    const itpp::cvec & get_rs(
      const uint8 & slot_num,
      const uint8 & sym_num
    ) const;
    double get_shift(
      const uint8 & slot_num,
      const uint8 & sym_num,
      const uint8 & port_num
    ) const;
  private:
    uint8 n_symb_dl;
    std::vector <itpp::cvec> table;
    itpp::mat shift_table;
};

// LTE convolutional ratematching and deratematching
itpp::cvec lte_conv_ratematch(
  const itpp::cmat & d,
  const uint32 & n_e
);
itpp::mat lte_conv_deratematch(
  const itpp::vec & e_est,
  const uint32 & n_c
);

// LTE convolutional encoding and decoding
itpp::bmat lte_conv_encode(
  const itpp::bvec & c
);
itpp::bvec lte_conv_decode(
  const itpp::mat & d_est
);

// Class to precompute all LTE modulation maps
class Mod_map {
  public:
    Mod_map(void);
    const itpp::cvec & operator()(
      const modulation_t::modulation_t & mod
    ) const;
  private:
    itpp::Array <itpp::cvec> table;
};

// Modulate and demodulate bits to/from symbols according to LTE specs.
itpp::cvec lte_modulate(
  const itpp::bvec & bits,
  const modulation_t::modulation_t modulation
);
itpp::vec lte_demodulate(
  const itpp::cvec & syms,
  const itpp::vec & np,
  const modulation_t::modulation_t & modulation
);

// Calculate one of the LTE CRC values.
itpp::bvec lte_calc_crc(
  const itpp::bvec & a,
  const crc_t crc
);

#endif

