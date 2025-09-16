/**
 *	Author	: Jaeuk Kim
 *	Email	: phy000.kim@gmail.com
 *	Date	: April 14th 2019 */

/** \file MonteCarlo.cpp 
 * \brief Function implementations for the MonteCarlo class and its derivatives. */

#include "MonteCarlo.h"
#include "Move.h"

virtual void MonteCarlo::Move(double Temperature, size_t Repeat, MCMove<System, AdditionalInfo> & move, AdditionalInfo add)
	{
		for(size_t i=0; i<Repeat; i++)
		{
			MoveCount++;
			double PreFactor=1.0;
			double dE=move.DeltaEnergy(*this->pSys, add, this->gen, PreFactor, LockStepSize, RecordEnergy);
			if(Temperature ==0)
			{
				if(dE>0)
					continue;
			}
			else if(gen.RandomDouble() > PreFactor*std::exp((-1)*dE/Temperature))
				continue;
				
			move.Accept(*this->pSys, add);
			this->RecordEnergy+=dE;
		}
	}

virtual void MonteCarlo::Anneal(AdditionalInfo add, CoolingSchedule & cool, MCMove<System, AdditionalInfo> & move, std::function<void(double & Energy, double Temperature, System & sys)> CallBack)
	{
		while(cool.Continue)
		{
			this->Move(cool.Temperature, cool.NumTry, move, add);
			CallBack(this->RecordEnergy, cool.Temperature, *pSys);

			cool.Report(this->RecordEnergy);
		}
	}