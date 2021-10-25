#include <iostream>
#include <stdio.h>
#include <string.h>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <math.h>
#include <float.h>
#include <vector>
#include <algorithm>

const int n = 20;
using namespace std;

int CataBM = 1;
double Rdeg = 0.00005;
double xEnv[n];
double popFreqMin = 0.0001;
int intT = 18000;
int nt = 1, Nenzyme = n/5;
double dt = 0.005;
double D_S = 1.0, Dout = 20.0; 
double Senv, Venv, T;
int pathNum = 2 * n;
int inter = 20000;
double kCoefficient = 1.0;


double getUniformRandom(void){
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double sign(double A){
	return (A > 0.0) - (A < 0.0);
}

void xEnvInitialize() {
	for(int i = 0; i < n; i++){
		if(i < nt) xEnv[i] = Senv;
		else xEnv[i] = 0.0;
	}
}

class Cell{
public:
	double mu, popFreq;
	double x[n];
	double Difs[n];
	vector< pair<int, int> > CataReact[n];
	double xEnvIso[n];
	int id;
	double prepopFreq, preLambda;
	int tau;
	int tauAfterTransition;
	double dD;
	int dDIndex;
	double muBefore;

	
	double k[n][4];
	double xTemp[n];
	double muTemp;


	Cell(double prop){
		popFreq = prop;
		mu = 0.0;
		tauAfterTransition = 0;
		id = -1;

		for(int i = 0; i < n; i++) x[i] = getUniformRandom();

		for(int i = 0; i < n; i++){
			if(i < nt) xEnvIso[i] = Senv;
			else xEnvIso[i] = 0.0;
		}

		for(int i = 0; i < nt; i++) Difs[i] = D_S;
		for(int i = nt; i < n; i++) Difs[i] = 0.0;
	}

	int makeRandomNet();
	int IsoDynamics();
};


inline int Cell::makeRandomNet(){
	for(int i = 0; i < n; i++) CataReact[i].clear();
	int pathNumTemp=0;
	while(pathNumTemp < pathNum){
		int sub = int( ( 1.0*(n-1) ) * getUniformRandom() );
		if(sub < nt){
			int pro = int( (1.0*(n-(nt+Nenzyme)))*getUniformRandom() ) + nt + Nenzyme; 
			int cata = int( (1.0*Nenzyme)*getUniformRandom() ) + nt;
			while(pro == cata) cata = int( (1.0*Nenzyme)*getUniformRandom() ) + nt;
			CataReact[sub].push_back( pair<int, int>(pro, cata) );
			pathNumTemp++;
		}
		if(sub >= nt+Nenzyme){
			int pro = int( (1.0*(n-nt))*getUniformRandom() ) + nt;
			while(sub == pro) pro = int( (1.0*(n-nt))*getUniformRandom() ) + nt;
			int cata = int( (1.0*Nenzyme)*getUniformRandom() ) + nt;
			while(pro == cata) cata = int((1.0*Nenzyme)*getUniformRandom() ) + nt;
			CataReact[sub].push_back( pair<int, int>(pro, cata) );
			pathNumTemp++;
		}
	}


	dDIndex = int( ( 1.0*(n - nt - Nenzyme) )*getUniformRandom() ) + nt+Nenzyme;

	return 0;
}

inline int Cell::IsoDynamics(){
	
	for(int i = 0; i < n; i++) for(int count = 0; count < 4; count++) k[i][count] = 0.0;

	
	for(int count = 0; count < 4; count++){
		for(int i = 0; i < n; i++){
			xTemp[i] = x[i];// 1st
			if(count == 1) if(!isnan(k[i][0])) xTemp[i] += k[i][0]*dt/2.0; // 2nd
			if(count == 2) if(!isnan(k[i][1])) xTemp[i] += k[i][1]*dt/2.0; // 3rd
			if(count == 3) if(!isnan(k[i][2])) xTemp[i] += k[i][2]*dt; // 4th
		}

		k[n-1][count] -= kCoefficient * xTemp[n-1] * xTemp[CataBM]; //Biomass Synthesis Reaction
  		for(int i = 0; i < n; i++){
  			for(int r = 0; r < CataReact[i].size(); r++){
      			int j = CataReact[i][r].first;
				int l = CataReact[i][r].second;
				double cDelta = kCoefficient * xTemp[i] * xTemp[l];
				k[j][count] += cDelta;
				k[i][count] -= cDelta;
			}
		}

		for(int i = 0; i < n; i++) if(Difs[i] > 0.0) k[i][count] += Difs[i] * (xEnvIso[i] - xTemp[i]);

		muTemp = kCoefficient * xTemp[n-1] * xTemp[CataBM]; //Growth Rate
		for(int i = 0; i < n; i++) k[i][count] -= muTemp * xTemp[i]; //Growth-dilution

		for(int i = 0; i < n; i++) if(isnan(k[i][count])) k[i][count] = 0.0;
	}

	for(int i = 0; i < n; i++) x[i] += (k[i][0] + 2.0*k[i][1] + 2.0*k[i][2] + k[i][3])*dt/6.0;
	for(int i = 0; i < n; i++) if(isnan(x[i])) x[i] = 0.0;
	
	mu = kCoefficient * x[n-1] * x[CataBM];

  return 0;
}


vector<Cell> Ecosystem;

int PopDynamics(){
	
	
	int m = Ecosystem.size();
	int varSize = m*n + n + m;

	// {x_i^α}:m Cell Species * n Chemical Species, {x_i^(env)}:n Chemical Species, {p_α}:m Cell Species
	vector< vector<double> > k( varSize, vector<double>(4, 0.0) );
	
	

	
	vector<double> xTemp( varSize, 0.0 );
	
	
	for(int count = 0; count < 4; count++){
		// 1st
		for(int alpha = 0; alpha < m; alpha++) for(int i = 0; i < n; i++) xTemp[alpha*n + i] = Ecosystem[alpha].x[i];
		for(int i = 0; i < n; i++) xTemp[m*n + i] = xEnv[i];
		for(int alpha = 0; alpha < m; alpha++) xTemp[m*n + n + alpha] = Ecosystem[alpha].popFreq;

		if(count == 1) for(int i = 0; i < varSize; i++) xTemp[i] += k[i][0]*dt/2.0; // 2nd
		if(count == 2) for(int i = 0; i < varSize; i++) xTemp[i] += k[i][1]*dt/2.0; // 3rd
		if(count == 3) for(int i = 0; i < varSize; i++) xTemp[i] += k[i][2]*dt; // 4th


		//Intracellular Chemical Reactions
		for(int alpha = 0; alpha < m; alpha++){
			k[alpha*n + n-1][count] -= kCoefficient * xTemp[alpha*n + n-1] * xTemp[alpha*n + CataBM]; // Biomass Synthesis Reaction
  			for(int i = 0; i < n; i++){
  				for(int r = 0; r < Ecosystem[alpha].CataReact[i].size(); r++){
      				int j = Ecosystem[alpha].CataReact[i][r].first;
					int l = Ecosystem[alpha].CataReact[i][r].second;
					double cDelta = kCoefficient*xTemp[alpha*n + i]*xTemp[alpha*n + l];
					k[alpha*n + j][count] += cDelta;
					k[alpha*n + i][count] -= cDelta;
				}
			}
		}
		//Intracellular Growth-dilution
		for(int alpha = 0; alpha < m; alpha++){
			Ecosystem[alpha].muTemp = kCoefficient * xTemp[alpha*n + n-1] * xTemp[alpha*n + CataBM];
			for(int i = 0; i < n; i++) k[alpha*n + i][count] -= Ecosystem[alpha].muTemp*xTemp[alpha*n + i];
		}

		//Chemical Exchange between Cells and Env.
		for(int alpha = 0; alpha < m; alpha++){
			for(int i = 0; i < n; i++){
				if(Ecosystem[alpha].Difs[i] > 0.0){
					double inflow = Ecosystem[alpha].Difs[i]*(xTemp[m*n + i] - xTemp[alpha*n + i]);
					k[alpha*n + i][count] += inflow;
					k[m*n + i][count] -= inflow*xTemp[m*n + n + alpha]/Venv;
				}
			}
		}
		// Environmental Chemical Dilution
		for(int i = 0; i < n; i++) k[m*n + i][count] -= Rdeg*xTemp[m*n + i];
		// Nutrient Supply into Env.
		for(int i = 0; i < nt; i++) k[m*n + i][count] += Dout*(Senv - xTemp[m*n + i]);

		//Population Dynamics
		double averaged_muTemp = 0.0;
		for(int alpha = 0; alpha < m; alpha++) averaged_muTemp += xTemp[m*n + n + alpha] * Ecosystem[alpha].muTemp;
		for(int alpha = 0; alpha < m; alpha++) k[m*n + n + alpha][count] += (Ecosystem[alpha].muTemp - averaged_muTemp) * xTemp[m*n + n + alpha];


		for(int i = 0; i < varSize; i++) if(isnan(k[i][count])) k[i][count] = 0.0;
	}


	for(int alpha = 0; alpha < m; alpha++){
		for(int i = 0; i < n; i++){
			Ecosystem[alpha].x[i] += (k[alpha*n + i][0] + 2.0*k[alpha*n + i][1] + 2.0*k[alpha*n + i][2] + k[alpha*n + i][3])*dt/6.0;
			if(isnan(Ecosystem[alpha].x[i])) Ecosystem[alpha].x[i] = 0.0; 
		}
		Ecosystem[alpha].mu = kCoefficient * Ecosystem[alpha].x[n-1] * Ecosystem[alpha].x[CataBM];
	}
	for(int i = 0; i < n; i++){
		xEnv[i] += (k[m*n + i][0] + 2.0*k[m*n + i][1] + 2.0*k[m*n + i][2] + k[m*n + i][3])*dt/6.0;
		if(isnan(xEnv[i])) xEnv[i] = 0.0; 
	}
	for(int alpha = 0; alpha < m; alpha++) Ecosystem[alpha].popFreq += (k[m*n + n + alpha][0] + 2.0*k[m*n + n + alpha][1] + 2.0*k[m*n + n + alpha][2] + k[m*n + n + alpha][3])*dt/6.0;
	for(int alpha = Ecosystem.size()-1; alpha >= 0; alpha--) if(isnan(Ecosystem[alpha].popFreq) || Ecosystem[alpha].popFreq < popFreqMin) Ecosystem.erase(Ecosystem.begin() + alpha);
	double popFreqSum=0.0;
	for(int alpha = 0; alpha < m; alpha++) popFreqSum += Ecosystem[alpha].popFreq;
	for(int alpha = 0; alpha < m; alpha++) Ecosystem[alpha].popFreq /= popFreqSum;

	return 0;
}


double Dleak;
vector<double> Dlevels; 

int selectOptimalDifs(int dNum){
	for(int i = 0; i < n; i++) Ecosystem[dNum].xEnvIso[i] = xEnv[i];
	for(int i = nt; i < n; i++) Ecosystem[dNum].Difs[i] = 0.0;
	for(double t = 0.0; t < T/30.0; t+=dt) Ecosystem[dNum].IsoDynamics();

	int intDLeakLevel[n];
	for(int i = 0; i < n; i++) intDLeakLevel[i] = -1;

	int LeakChemNum = 0;
	for(int z = nt+Nenzyme; z < n; z++){
		intDLeakLevel[z]=0;
		LeakChemNum++;
	}

	Cell cellTemp(1.0);
	cellTemp = Ecosystem[dNum];
	while(LeakChemNum > 0){
		for(int z = nt+Nenzyme; z < n; z++){
			if(intDLeakLevel[z] >= 0){
				cellTemp.Difs[z] = Dlevels[intDLeakLevel[z]];
				for(double t = 0.0; t < T/30.0 ; t+=dt) cellTemp.IsoDynamics();

				if(cellTemp.mu > Ecosystem[dNum].mu){
					Ecosystem[dNum].mu = cellTemp.mu;
					Ecosystem[dNum].Difs[z] = cellTemp.Difs[z];
					intDLeakLevel[z]++;
				}else{
					intDLeakLevel[z] = -1;
					LeakChemNum--;
					cellTemp.mu = Ecosystem[dNum].mu;
					cellTemp.Difs[z] = Ecosystem[dNum].Difs[z];
				}
				if(intDLeakLevel[z] == Dlevels.size()){
					intDLeakLevel[z] = -1;
					LeakChemNum--;
				}
			}
		}
	}

	return 0;
}

int selectOptimalDifsMut(){
	int LeakChemNum = 0;
	int intDLeakLevel[n];
	for(int y = 0; y < n; y++) intDLeakLevel[y] = -1;
	for(int y = nt+Nenzyme; y < n; y++){
		intDLeakLevel[y] = 0;
		LeakChemNum++;
	}

	Cell cellTemp(1.0);
	cellTemp = Ecosystem[0];
	while(LeakChemNum > 0){
		for(int z = nt+Nenzyme; z < n; z++){
			if(intDLeakLevel[z] >= 0){
				Ecosystem[0].Difs[z] = Dlevels[intDLeakLevel[z]];
				xEnvInitialize();
				for(double t = 0.0; t < T/30.0 ; t+=dt) PopDynamics();

				if(Ecosystem[0].mu > cellTemp.mu){
					cellTemp.mu = Ecosystem[0].mu;
					cellTemp.Difs[z] = Ecosystem[0].Difs[z];
					intDLeakLevel[z]++;
				}else{
					intDLeakLevel[z] = -1;
					LeakChemNum--;
					Ecosystem[0].mu = cellTemp.mu;
					Ecosystem[0].Difs[z] = cellTemp.Difs[z];
				}
				if(intDLeakLevel[z] == Dlevels.size()){
					intDLeakLevel[z] = -1;
					LeakChemNum--;
				}
			}
		}
	}

	return 0;
}

vector<int> invadingOrder;
void shuffle() {
    for(int i = 0; i < invadingOrder.size(); i++) {
        int j = rand()%invadingOrder.size();
        int t = invadingOrder[i];
        invadingOrder[i] = invadingOrder[j];
        invadingOrder[j] = t;
    }
}


int main(int argc, char *argv[]){
	int tane = atoi(argv[1]);
	Senv = atof(argv[2]);
	Venv = atof(argv[3]);

	int SpeNum = 50;
	int SpeNumPrep = 75;
	double Dl = 0.001;
	for(int u1 = 0; u1 < 15; u1++){
		Dlevels.push_back(Dl*pow(2.0,1.0*u1));
		if(Dl*pow(2.0,1.0*u1)>1.0) break;
	}
	Dleak = 0.001;
	double TauOrder = 10.0;
	double pChange = 0.1;
	double deltaP = 0.0025;
	double muMin = 5.0e-5;
	T = dt*(1.0*inter)*(1.0*intT);

	char output3B[256] = "Output3B.txt";
	char output3C[256] = "Output3C.txt";
	ofstream of;
	of.open(output3C, ios::trunc);
	of << "# n = "<< n <<", Senv = "<< Senv <<", Venv = "<< Venv;
	of.close();
	of.open(output3B, ios::trunc);
	of << "# n = "<< n <<", Senv = "<< Senv <<", Venv = "<< Venv;

	int sampleSize = 1;
	vector<int> CoExNums;
	vector<Cell> InvasionSpecies;
	Ecosystem.reserve(SpeNum); 

	vector< pair<int, double> > SymGrowthRates;
	vector<double> SymGrowthRatio;
	for(int index = 0; index < sampleSize; index++){
		srand( tane );
		tane+=5;

		InvasionSpecies.clear();
		while(InvasionSpecies.size() < SpeNumPrep){
			Cell cellTemp(1.0);
			cellTemp.makeRandomNet();
			Ecosystem.clear();
			Ecosystem.push_back(cellTemp);
			xEnvInitialize();

			for(double t = 0.0; t < T/10.0 ; t+=dt) PopDynamics();

			if(!isnan(Ecosystem[0].mu) && !isinf(Ecosystem[0].mu) && Ecosystem[0].mu > muMin){
				Ecosystem[0].tau = int(TauOrder/Ecosystem[0].mu);
				InvasionSpecies.push_back(Ecosystem[0]);
			}
		}

		for(int u1 = 0; u1 < InvasionSpecies.size(); u1++){
			for(int u2 = u1+1; u2 < InvasionSpecies.size(); u2++){
				if(InvasionSpecies[u1].mu > InvasionSpecies[u2].mu){
					Cell cSort(1.0);
					cSort = InvasionSpecies[u2];
					InvasionSpecies[u2] = InvasionSpecies[u1];
					InvasionSpecies[u1] = cSort;
				}
			}
		}
		for(int alpha = 0; alpha < SpeNum; alpha++) InvasionSpecies[alpha].id = alpha;

		Ecosystem.clear();
		Ecosystem.push_back(InvasionSpecies[SpeNum - 1]);
		xEnvInitialize();
		selectOptimalDifsMut();
		double moreFromStrongestCell=0.0;
		for(int i = nt; i < n; i++) moreFromStrongestCell += Ecosystem[0].Difs[i];
		int eraseSpeciesNum;
		while(moreFromStrongestCell < Dleak/10.0){
			eraseSpeciesNum = int( (1.0*SpeNum)*getUniformRandom() );
			InvasionSpecies.erase(InvasionSpecies.begin() + eraseSpeciesNum);

			for(int alpha = 0; alpha < SpeNum; alpha++) InvasionSpecies[alpha].id = alpha;

			Ecosystem.clear();
			Ecosystem.push_back(InvasionSpecies[SpeNum-1]);
			xEnvInitialize();
			selectOptimalDifsMut();
			moreFromStrongestCell=0.0;
			for(int alpha = nt; alpha < n; alpha++) moreFromStrongestCell += Ecosystem[0].Difs[alpha];
		}
		
		for(int alpha = InvasionSpecies.size()-1; alpha >= SpeNum; alpha--) InvasionSpecies.erase(InvasionSpecies.begin() + alpha);
		
		vector<Cell> EcosystemBefore;
		double xEnvBefore[n];
		CoExNums.clear();
		CoExNums.push_back(SpeNum-1);
		double Di, influx;
		int dDNext;
		invadingOrder.clear();
		for(int alpha = SpeNum-1; alpha >= 0; alpha--) invadingOrder.push_back(alpha);
	  	//shuffle();

		of << "# # of invasion, Growth rate, # of coexisting species, whether successful invaision or not (>0 or -1)";
		if(moreFromStrongestCell > Dleak/10.0){
			for(int y = 0; y < SpeNum; y++){
				int invadingSpecies = invadingOrder[y];

				int judgeCoEx = 1;
				for(int u1 = 0; u1 < Ecosystem.size(); u1++) if(invadingSpecies == Ecosystem[u1].id) judgeCoEx = -1;
				if(judgeCoEx > 0){
					EcosystemBefore.clear();
					EcosystemBefore = Ecosystem;
					for(int i = 0; i < n; i++) xEnvBefore[i]=xEnv[i];
					for(int u1 = 0; u1 < Ecosystem.size(); u1++) Ecosystem[u1].popFreq = Ecosystem[u1].popFreq*(1.0-deltaP);
					InvasionSpecies[invadingSpecies].popFreq = deltaP;
					Ecosystem.push_back(InvasionSpecies[invadingSpecies]);
					selectOptimalDifs(Ecosystem.size()-1);

					for(double t = 0.0; t < T/4.0 ; t+=dt) PopDynamics();

					for(int y2 = nt; y2 < n; y2++){
						judgeCoEx = -1;
						for(int u1 = 0; u1 < Ecosystem.size(); u1++) if(Ecosystem[u1].popFreq > 0.0 && Ecosystem[u1].Difs[y2] > 0.0) judgeCoEx = 1;
						if(judgeCoEx < 0) xEnv[y2] = 0.0;
					}

					int PositivepopFreqNum=0;
					for(int u1 = 0; u1 < Ecosystem.size(); u1++){
						if(Ecosystem[u1].popFreq > 0.0) PositivepopFreqNum++;
						if(isnan(Ecosystem[u1].popFreq)) PositivepopFreqNum = -2 * SpeNum;
					}
					if(PositivepopFreqNum <= 0){
						Ecosystem.clear();
						Ecosystem = EcosystemBefore;
						for(int i = 0; i <n; i++) xEnv[i] = xEnvBefore[i];
					}
					if(PositivepopFreqNum >= 1){
						if(Ecosystem[0].mu < muMin){
							Ecosystem.clear();
							Ecosystem = EcosystemBefore;
							for(int i = 0; i <n; i++) xEnv[i] = xEnvBefore[i];
						}
					}

					int AdaptationDynamicsJudge=-1;
					if(Ecosystem.size() > CoExNums.size()) AdaptationDynamicsJudge = 1;
					if(Ecosystem.size() == CoExNums.size()){
						if(Ecosystem[Ecosystem.size()-1].id == invadingSpecies) AdaptationDynamicsJudge = 1;
						else AdaptationDynamicsJudge = -1;
					}
					if(Ecosystem.size() < CoExNums.size()){
						if(Ecosystem[Ecosystem.size()-1].id == invadingSpecies) AdaptationDynamicsJudge = 1;
						else{
							Ecosystem.clear();
							Ecosystem = EcosystemBefore;
							AdaptationDynamicsJudge = -1;
						}
					}
					

					of << endl << y << " " << Ecosystem[0].mu << " " << Ecosystem.size() << " ";
					if(AdaptationDynamicsJudge > 0) of << Ecosystem.size(); 
					else of << -1;

					if(AdaptationDynamicsJudge > 0){ 
						int tauSym=0;
						for(int s = 0; s < Ecosystem.size(); s++) if(Ecosystem[s].tau > tauSym) tauSym = Ecosystem[s].tau;

						int it=0;
						int itSym=0;
						for(double t = 0.0; t < T ; t+=dt){
							PopDynamics();

							if(it % 200==0){
								for(int u1 = 0; u1 < Ecosystem.size(); u1++){
									Ecosystem[u1].prepopFreq = Ecosystem[u1].popFreq;
									if(Ecosystem[u1].popFreq == 1.0) t = T;
									if((itSym/200) % tauSym == 0){
										Ecosystem[u1].dD = sign(Ecosystem[u1].mu - Ecosystem[u1].muBefore)*sign(Ecosystem[u1].dD)*Dleak;
										Ecosystem[u1].Difs[Ecosystem[u1].dDIndex] += Ecosystem[u1].dD;

										Ecosystem[u1].muBefore = Ecosystem[u1].mu;
										Ecosystem[u1].tau = int(TauOrder/Ecosystem[u1].mu);
										Ecosystem[u1].tauAfterTransition++;

										if(Ecosystem[u1].Difs[Ecosystem[u1].dDIndex] < 0.0){
											Ecosystem[u1].Difs[Ecosystem[u1].dDIndex] = 0.0;
											dDNext = Ecosystem[u1].dDIndex;
											while(dDNext == Ecosystem[u1].dDIndex) dDNext = int((1.0*(n-(nt+Nenzyme)))*getUniformRandom())+nt+Nenzyme;
											Ecosystem[u1].dDIndex = dDNext;

											Ecosystem[u1].dD = Dleak;
											Ecosystem[u1].Difs[Ecosystem[u1].dDIndex] += Ecosystem[u1].dD;

											Ecosystem[u1].tauAfterTransition = 0;
										}

										if(getUniformRandom()<pChange && Ecosystem[u1].tauAfterTransition >= 2){
											Ecosystem[u1].dD = Dleak;
											dDNext = Ecosystem[u1].dDIndex;
											while(dDNext == Ecosystem[u1].dDIndex) dDNext = int((1.0*(n-(nt+Nenzyme)))*getUniformRandom())+nt+Nenzyme;
											Ecosystem[u1].dDIndex = dDNext;

											Ecosystem[u1].Difs[Ecosystem[u1].dDIndex] += Ecosystem[u1].dD;

											Ecosystem[u1].tauAfterTransition=0;
										}
									}
								}
								if((itSym/200) % tauSym == 0){
									itSym -= tauSym*200;
									for(int s = 0; s < Ecosystem.size(); s++) if(Ecosystem[s].popFreq > popFreqMin && Ecosystem[s].tau > tauSym) tauSym = Ecosystem[s].tau;
								}
							}
							it++;
							itSym++;
						}

						for(double t=0.0; t<T/4.0; t +=dt) PopDynamics();
						of << endl << (1.0*y)+0.5 << " " << Ecosystem[0].mu << " " << Ecosystem.size() << " " << -1; // after adaptation of diffusion coefficients
					}


					CoExNums.clear();
					for(int u1 = 0; u1 < Ecosystem.size(); u1++) CoExNums.push_back(Ecosystem[u1].id);
				}
			}
		}


		of.close();
		of.open(output3C, ios::app);
		of << endl << "# Vertical axis: Chemical species. Horizontal axis: Cell species.";
		for(int u1 = 0; u1 < Ecosystem.size(); u1++){
			for(int u2 = 0; u2 < n; u2++){
				influx = Ecosystem[u1].Difs[u2]*(xEnv[u2] - Ecosystem[u1].x[u2]);
				of << endl << u1 << " " << u2 << " " << influx;
				of << endl << u1 << " " << u2+1 << " " << influx;
			}
			of << endl;
			for(int u2 = 0; u2 < n; u2++){
				influx=Ecosystem[u1].Difs[u2]*(xEnv[u2] - Ecosystem[u1].x[u2]);
				of << endl << u1+1 << " " << u2 << " " << influx;
				of << endl << u1+1 << " " << u2+1 << " " << influx;
			}
			of << endl;
		}
	}

	of.close();

	return 0;
}
