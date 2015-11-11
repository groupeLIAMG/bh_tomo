/*
 *  ttRaisCourbes.h
 *  ttRaisCourbes
 *
 *  Created by Bernard Giroux on 08-07-01.
 *  Copyright 2008 Bernard Giroux. All rights reserved.
 *
 */

/*
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
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

#ifndef __TTRAISCOURBES_T_H__
#define __TTRAISCOURBES_T_H__

const double small = 1.e-10;

template<typename T>
struct sxz {
    T x;
    T z;
};

template<typename T>
struct sxyz {
    T x;
    T y;
    T z;
};

template<typename T>
struct siv {
    size_t i;   // index
    T v;        // value
};

template<typename T>
struct siv2 {
    size_t i;    // index
    T v1;        // first value
    T v2;        // second value
};

template<typename T>
class CompareSiv {
public:
    bool operator()(const siv<T> n1, const siv<T> n2) const {
        return n1.i < n2.i;
    }
};

template<typename T>
class CompareSiv2 {
public:
    bool operator()(const siv2<T> n1, const siv2<T> n2) const {
        return n1.i < n2.i;
    }
};


template<typename T>
struct txPar {
    sxz<T> pt;
    T t0;
    T theta;
    T diam;
    bool inWater;
};

template<typename T>
struct rxPar {
    std::vector<sxz<T> > pts;
    std::vector<T> theta;
    std::vector<T> diam;
    std::vector<bool> inWater;
};



#endif
