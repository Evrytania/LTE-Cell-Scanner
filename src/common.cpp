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

#include <itpp/itbase.h>
#include <boost/math/special_functions/gamma.hpp>
#include "rtl-sdr.h"
#include "common.h"
#include "macros.h"
#include "itpp_ext.h"
#include "dsp.h"

using namespace std;

// Cell class member functions
int16 const Cell::n_id_cell() const {
  return ((n_id_1>=0)&&(n_id_2>=0))?(n_id_2+3*n_id_1):-1;
}
int8 const Cell::n_symb_dl() const {
  return (cp_type==cp_type_t::NORMAL)?7:((cp_type==cp_type_t::EXTENDED)?6:-1);
}
// Default constructor initializes structure to impossible values.
Cell::Cell() :
  fc_requested(NAN),
  fc_programmed(NAN),
  pss_pow(NAN),
  ind(-1),
  freq(NAN),
  n_id_2(-1),

  n_id_1(-1),
  cp_type(cp_type_t::UNKNOWN),
  frame_start(NAN),
  freq_fine(NAN),

  freq_superfine(NAN),

  n_ports(-1),
  n_rb_dl(-1),
  phich_duration(phich_duration_t::UNKNOWN),
  phich_resource(phich_resource_t::UNKNOWN),
  sfn(-1)
{}

// Overload << to allow easy printing of 'Cell'.
ostream & operator<< (
  ostream & os,
  const Cell & c
) {
  if (isnan(c.fc_requested)&&isnan(c.fc_programmed)&&isnan(c.pss_pow)&&(c.ind==-1)&&isnan(c.freq)&&(c.n_id_2==-1)) {
    os << "<EMPTY>";
    return os;
  }
  os << "fc_requested = " << c.fc_requested/1e6 << " MHz" << endl;
  os << "fc_programmed = " << c.fc_programmed/1e6 << " MHz" << endl;
  os << "pss_pow = " << db10(c.pss_pow) << " dB" << endl;
  os << "ind = " << c.ind << endl;
  os << "freq = " << c.freq << endl;
  os << "n_id_2 = " << c.n_id_2;

  if ((c.n_id_1==-1)&&(c.cp_type==cp_type_t::UNKNOWN)&&isnan(c.frame_start)&&isnan(c.freq_fine))
    return os;

  os << endl;
  os << "n_id_1 = " << c.n_id_1 << endl;
  os << "n_id_cell = " << c.n_id_cell() << endl;
  os << "cp_type = " << c.cp_type << endl;
  os << "frame_start = " << c.frame_start;

  if (isnan(c.freq_fine))
    return os;

  cout << endl;
  os << "freq_fine = " << c.freq_fine;

  if (isnan(c.freq_superfine))
    return os;

  os << endl;
  os << "freq_superfine = " << c.freq_superfine;

  if ((c.n_ports==-1)&&(c.n_rb_dl==-1)&&(c.phich_duration==phich_duration_t::UNKNOWN)&&(c.phich_resource==phich_resource_t::UNKNOWN)&&(c.sfn==-1))
    return os;

  os << endl;
  os << "n_ports = " << c.n_ports << endl;
  os << "n_rb_dl = " << c.n_rb_dl << endl;
  os << "phich_duration = " << c.phich_duration << endl;
  os << "phich_resource = " << c.phich_resource << endl;
  os << "sfn = " << c.sfn;

  return os;
}

