// Code originally came from osmocom:
// http://sdr.osmocom.org/trac/wiki/rtl-sdr
// From tuner_e4k.c

/*
 * Elonics E4000 tuner driver
 *
 * (C) 2011-2012 by Harald Welte <laforge@gnumonks.org>
 * (C) 2012 by Sylvain Munaut <tnt@246tNt.com>
 * (C) 2012 by Hoernchen <la@tfc-server.de>
 *
 * All Rights Reserved
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
#include <queue>
#include <sys/stat.h>
#include <sys/syscall.h>
#include <sys/types.h>
#include <curses.h>
#include "rtl-sdr.h"
#include "common.h"

struct e4k_pll_params {
  uint32 fosc;
  uint32 intended_flo;
  uint32 flo;
  uint16 x;
  uint8 z;
  uint8 r;
  uint8 r_idx;
  uint8 threephase;
};

#define ARRAY_SIZE(arr) (sizeof(arr) / sizeof((arr)[0]))

#define KHZ(x)	((x)*1000)

#define E4K_PLL_Y		65536

struct pll_settings {
	uint32_t freq;
	uint8_t reg_synth7;
	uint8_t mult;
};

static const struct pll_settings pll_vars[] = {
	{KHZ(72400),	(1 << 3) | 7,	48},
	{KHZ(81200),	(1 << 3) | 6,	40},
	{KHZ(108300),	(1 << 3) | 5,	32},
	{KHZ(162500),	(1 << 3) | 4,	24},
	{KHZ(216600),	(1 << 3) | 3,	16},
	{KHZ(325000),	(1 << 3) | 2,	12},
	{KHZ(350000),	(1 << 3) | 1,	8},
	{KHZ(432000),	(0 << 3) | 3,	8},
	{KHZ(667000),	(0 << 3) | 2,	6},
	{KHZ(1200000),	(0 << 3) | 1,	4}
};

/* \brief compute Fvco based on Fosc, Z and X
 * \returns positive value (Fvco in Hz), 0 in case of error */
static uint64_t compute_fvco(uint32_t f_osc, uint8_t z, uint16_t x)
{
	uint64_t fvco_z, fvco_x, fvco;

	/* We use the following transformation in order to
	 * handle the fractional part with integer arithmetic:
	 *  Fvco = Fosc * (Z + X/Y) <=> Fvco = Fosc * Z + (Fosc * X)/Y
	 * This avoids X/Y = 0.  However, then we would overflow a 32bit
	 * integer, as we cannot hold e.g. 26 MHz * 65536 either.
	 */
	fvco_z = (uint64_t)f_osc * z;

	fvco_x = ((uint64_t)f_osc * x) / E4K_PLL_Y;

	fvco = fvco_z + fvco_x;

	return fvco;
}

static uint32_t compute_flo(uint32_t f_osc, uint8_t z, uint16_t x, uint8_t r)
{
	uint64_t fvco = compute_fvco(f_osc, z, x);
	if (fvco == 0)
		return -EINVAL;

	return fvco / r;
}

double compute_fc_programmed(const double & fosc,const double & intended_flo)
{
        e4k_pll_params oscp_actual;
        e4k_pll_params * oscp=&oscp_actual;

	uint32_t i;
	uint8_t r = 2;
	uint64_t intended_fvco, remainder;
	uint64_t z = 0;
	uint32_t x;
	int flo;
	int three_phase_mixing = 0;
	oscp->r_idx = 0;

	//if (!is_fosc_valid(fosc))
	//	return 0;

	for(i = 0; i < ARRAY_SIZE(pll_vars); ++i) {
		if(intended_flo < pll_vars[i].freq) {
			three_phase_mixing = (pll_vars[i].reg_synth7 & 0x08) ? 1 : 0;
			oscp->r_idx = pll_vars[i].reg_synth7;
			r = pll_vars[i].mult;
			break;
		}
	}

	//fprintf(stderr, "[E4K] Fint=%u, R=%u\n", intended_flo, r);

	/* flo(max) = 1700MHz, R(max) = 48, we need 64bit! */
	intended_fvco = ((uint64_t)intended_flo) * r;

	/* compute integral component of multiplier */
	z = intended_fvco / fosc;

	/* compute fractional part.  this will not overflow,
	* as fosc(max) = 30MHz and z(max) = 255 */
	remainder = intended_fvco - (fosc * z);
	/* remainder(max) = 30MHz, E4K_PLL_Y = 65536 -> 64bit! */
	x = (remainder * E4K_PLL_Y) / fosc;
	/* x(max) as result of this computation is 65536 */

	flo = compute_flo(fosc, z, x, r);

	oscp->fosc = fosc;
	oscp->flo = flo;
	oscp->intended_flo = intended_flo;
	oscp->r = r;
//	oscp->r_idx = pll_vars[i].reg_synth7 & 0x0;
	oscp->threephase = three_phase_mixing;
	oscp->x = x;
	oscp->z = z;

	return flo;
}

