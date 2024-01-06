/*******************************************************************************
* Copyright (c) 2015-2017
* School of Electrical, Computer and Energy Engineering, Arizona State University
* PI: Prof. Shimeng Yu
* All rights reserved.
* 
* This source code is part of NeuroSim - a device-circuit-algorithm framework to benchmark 
* neuro-inspired architectures with synaptic devices(e.g., SRAM and emerging non-volatile memory). 
* Copyright of the model is maintained by the developers, and the model is distributed under 
* the terms of the Creative Commons Attribution-NonCommercial 4.0 International Public License 
* http://creativecommons.org/licenses/by-nc/4.0/legalcode.
* The source code is free and you can redistribute and/or modify it
* by providing that the following conditions are met:
* 
*  1) Redistributions of source code must retain the above copyright notice,
*     this list of conditions and the following disclaimer.
* 
*  2) Redistributions in binary form must reproduce the above copyright notice,
*     this list of conditions and the following disclaimer in the documentation
*     and/or other materials provided with the distribution.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
* FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
* DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
* SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
* OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
* OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Developer list: 
*   Pai-Yu Chen	    Email: pchen72 at asu dot edu 
*                    
*   Xiaochen Peng   Email: xpeng15 at asu dot edu
********************************************************************************/

#include <cstdio>
#include <random>
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <sstream>
#include <chrono>
#include <algorithm>
#include "constant.h"
#include "formula.h"
#include "Param.h"
#include "InputParameter.h"
#include "Technology.h"
#include "MemCell.h"
#include "Definition.h"

#include "DFF.h"

using namespace std;


int main(int argc, char * argv[]) {   

	auto start = chrono::high_resolution_clock::now();
	gen.seed(0);
	{
		inputParameter.temperature = param->temp;
		inputParameter.processNode = param->technode;
		tech.Initialize(inputParameter.processNode, inputParameter.deviceRoadmap, inputParameter.transistorType);

		cell.resistanceOn = param->resistanceOn;
		cell.resistanceOff = param->resistanceOff;
		cell.resistanceAvg = (cell.resistanceOn + cell.resistanceOff)/2;
		cell.readVoltage = param->readVoltage;
		cell.readPulseWidth = param->readPulseWidth;
		cell.accessVoltage = param->accessVoltage;
		cell.resistanceAccess = param->resistanceAccess;
		cell.featureSize = param->featuresize;
		cell.writeVoltage = param->writeVoltage;

		if (cell.memCellType == Type::SRAM) {   // SRAM
			cell.heightInFeatureSize = param->heightInFeatureSizeSRAM;                   // Cell height in feature size
			cell.widthInFeatureSize = param->widthInFeatureSizeSRAM;                     // Cell width in feature size
			cell.widthSRAMCellNMOS = param->widthSRAMCellNMOS;
			cell.widthSRAMCellPMOS = param->widthSRAMCellPMOS;
			cell.widthAccessCMOS = param->widthAccessCMOS;
			cell.minSenseVoltage = param->minSenseVoltage;
		} else {
			cell.heightInFeatureSize = (cell.accessType==CMOS_access)? param->heightInFeatureSize1T1R : param->heightInFeatureSizeCrossbar;         // Cell height in feature size
			cell.widthInFeatureSize = (cell.accessType==CMOS_access)? param->widthInFeatureSize1T1R : param->widthInFeatureSizeCrossbar;            // Cell width in feature size
		} 

		if (cell.memCellType == Type::SRAM) {
			double capCellAccess = CalculateDrainCap(cell.widthAccessCMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, cell.widthInFeatureSize / ((tech.featureSize <= 14*1e-9)? (MAX_TRANSISTOR_HEIGHT_FINFET/MAX_TRANSISTOR_HEIGHT):1) * tech.featureSize, tech);
			cell.capSRAMCell = capCellAccess + CalculateDrainCap(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, NMOS, cell.widthInFeatureSize / ((tech.featureSize <= 14*1e-9)? (MAX_TRANSISTOR_HEIGHT_FINFET/MAX_TRANSISTOR_HEIGHT):1) * tech.featureSize, tech) + CalculateDrainCap(cell.widthSRAMCellPMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, PMOS, cell.widthInFeatureSize / ((tech.featureSize <= 14*1e-9)? (MAX_TRANSISTOR_HEIGHT_FINFET/MAX_TRANSISTOR_HEIGHT):1) * tech.featureSize, tech) + CalculateGateCap(cell.widthSRAMCellNMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech) + CalculateGateCap(cell.widthSRAMCellPMOS * ((tech.featureSize <= 14*1e-9)? 2:1) * tech.featureSize, tech);
		} else if (cell.memCellType == Type::RRAM || cell.memCellType == Type::FeFET) {
			if (cell.accessType == CMOS_access) {	// 1T1R
				cell.resCellAccess = cell.resistanceOn * IR_DROP_TOLERANCE;    //calculate access CMOS resistance
				cell.widthAccessCMOS = CalculateOnResistance(((tech.featureSize <= 14*1e-9)? 2:1)*tech.featureSize, NMOS, inputParameter.temperature, tech) * LINEAR_REGION_RATIO / cell.resCellAccess;   //get access CMOS width
				cell.resMemCellOn = cell.resCellAccess + cell.resistanceOn;        //calculate single memory cell resistance_ON
				cell.resMemCellOff = cell.resCellAccess + cell.resistanceOff;      //calculate single memory cell resistance_OFF
				cell.resMemCellAvg = cell.resCellAccess + cell.resistanceAvg;      //calculate single memory cell resistance_AVG
			} else {
				cell.resMemCellOn = cell.resistanceOn;
				cell.resMemCellOff = cell.resistanceOff;
				cell.resMemCellOnAtHalfVw = cell.resistanceOn;
				cell.resMemCellOffAtHalfVw = cell.resistanceOff;
				cell.resMemCellOnAtVw = cell.resistanceOn;
				cell.resMemCellOffAtVw = cell.resistanceOff;
				cell.resMemCellAvg = cell.resistanceAvg;
				cell.resMemCellAvgAtHalfVw = cell.resistanceAvg;
				cell.resMemCellAvgAtVw = cell.resistanceAvg;
			}

		}
	}


	DFF *dff1 = new DFF(inputParameter, tech, cell);
	dff1->Initialize(2, 1e9);
	dff1->CalculateArea(1e-5, 0, NONE);
	dff1->CalculateLatency(1e20, 1);
	dff1->CalculatePower(1, 2, true);

	printf("DFF area = %e (um^2)\n", dff1->area);
	printf("DFF read latency = %e (ns)\n", dff1->readLatency);
	printf("DFF read dynamic energy = %e (nJ)\n", dff1->readDynamicEnergy);
	printf("DFF leakage = %e (nA)\n", dff1->leakage);

	return 0;
}	



