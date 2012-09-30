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

#ifndef HAVE_CONSTANTS_H
#define HAVE_CONSTANTS_H

// Place all the tables into one structure
class Rom_tables {
  public:
    PSS_fd pss_fd;
    PSS_td pss_td;
    SSS_fd sss_fd;
    Mod_map mod_map;
};

extern Rom_tables ROM_TABLES;

#define FS_LTE 30720000.0
#define N_RB_MAXDL 110

#define CELL_DROP_THRESHOLD 400.0

#endif

