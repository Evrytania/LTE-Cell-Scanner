// An OpenCL accelerated LTE Cell Scanner
//
// Written by Jiao Xianjun <putaoshu@gmail.com>
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

// Header file to store coefficients of FIR 6RB channel filter

#ifndef HAVE_FILTER_COEF_H
#define HAVE_FILTER_COEF_H

// filter coefficients
// This is for non-OpenCL version. You should change coefficients in OpenCL kernel files for OpenCL version
const float chn_6RB_filter_coef[47] = { \
8.193313185354206e-04,     3.535548569572820e-04,    -1.453429245341695e-03,     1.042805860697287e-03,     1.264224526451337e-03, \
  -3.219586065044259e-03,     1.423981657254563e-03,     3.859884310477692e-03,    -6.552708013395765e-03,     8.590509694961493e-04, \
  9.363722386299336e-03,    -1.120357391780316e-02,    -2.423088424232164e-03,     1.927528718829535e-02,    -1.646405738285926e-02, \
  -1.143040384534755e-02,     3.652830082843752e-02,    -2.132986170036144e-02,    -3.396829121834471e-02,     7.273086636811442e-02, \
  -2.476823886110626e-02,    -1.207789042999466e-01,     2.861583432079335e-01,     6.398255789896659e-01,     2.861583432079335e-01, \
  -1.207789042999466e-01,    -2.476823886110626e-02,     7.273086636811442e-02,    -3.396829121834471e-02,    -2.132986170036144e-02, \
  3.652830082843752e-02,    -1.143040384534755e-02,    -1.646405738285926e-02,     1.927528718829535e-02,    -2.423088424232164e-03, \
  -1.120357391780316e-02,     9.363722386299336e-03,     8.590509694961493e-04,    -6.552708013395765e-03,     3.859884310477692e-03, \
  1.423981657254563e-03,    -3.219586065044259e-03,     1.264224526451337e-03,     1.042805860697287e-03,    -1.453429245341695e-03, \
  3.535548569572820e-04,     8.193313185354206e-04
};

#endif
