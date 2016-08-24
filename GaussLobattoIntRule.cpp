// Copyright (C) 2016 William H. Greene
//
// This program is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.

#include <math.h>
#include <assert.h>

#include "GaussLobattoIntRule.h"


static const double xi5p2 = sqrt(21.)/7.;
static const double wt5p2 = 49./90.;

static GaussLobattoIntRule::IRule rule5Pt[] = {
    { -1, .1 }, { -xi5p2, wt5p2 }, { 0., 32. / 45. },
    { xi5p2, wt5p2 }, { 1, .1 }
};


static  GaussLobattoIntRule::IRule *rules[] = {
  0, 0, 0, 0, 0, rule5Pt };

GaussLobattoIntRule::GaussLobattoIntRule(int npts) : numPts(npts)
{
  assert(npts==5);
  rule = rules[npts];
}

const GaussLobattoIntRule::IRule &GaussLobattoIntRule::getPoint(int iPt) const
{
  return rule[iPt];
}