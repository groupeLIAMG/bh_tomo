/*
 *  AntennaCorrection.h
 *  ttRaisCourbes
 *
 *  Created by Bernard Giroux on 08-05-11.
 *  Copyright 2008 École Polytechnique de Montréal. All rights reserved.
 *
 */

/*
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
 *
 */


#ifndef __ANTENNACORRECTION_H__
#define __ANTENNACORRECTION_H__

class AntennaCorrection {
public:
    virtual ~AntennaCorrection() = 0;
    virtual double getSlowness(double, double, bool) const = 0;
    virtual double getLength() const = 0;
	virtual const char *getName() const = 0;
};

class Ramac_250 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const;
    double getLength() const { return 1.0; }
	const char *getName() const { return "Ramac 250"; }
    
private:
    
    static const size_t nd = 12, nk = 18;
    
    static const double h_diam[];
    static const double k_form[];
    static const double slown_air[];
    static const double slown_water[];
};

class Fixed_05 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 5.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 05"; }
};    

class Fixed_06 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 6.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 06"; }
};    

class Fixed_07 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 7.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 07"; }
};    

class Fixed_08 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 8.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 08"; }
};    

class Fixed_09 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 9.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 09"; }
};    

class Fixed_10 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 10.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 10"; }
};    

class Fixed_11 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 11.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 11"; }
};    

class Fixed_12 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 12.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 12"; }
};    

class Fixed_13 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 13.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 13"; }
};    

class Fixed_14 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 14.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 14"; }
};    

class Fixed_15 : public AntennaCorrection {
public:
    double getSlowness(double k, double d, bool inWater) const { return 15.0; }
    double getLength() const { return 1.0; }
	const char *getName() const { return "Fixed 15"; }
};    

#endif
