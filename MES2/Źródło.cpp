#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <regex>
#include <iomanip>
#include <math.h>

using namespace std;

#define Number_Of_Dimensions 2
#define Number_Of_Points 4

double f1(double x) {
	return 5 * pow(x, 2) + 3 * x + 6;
}

double f2(double x, double y) {
	return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
}

double det(double matrix[2][2]) {
	return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
}

// Struktury


struct GlobalData {
	unsigned int SimulationTime;
	unsigned int SimulationStepTime;
	int Conductivity;
	int Alfa;
	int Tot;
	int InitialTemp;
	int Density;
	int SpecificHeat;
};

GlobalData GD;

struct node {
	double x;
	double y;
	short int BC = 0;
};

struct integration { 
	unsigned int NumberOfPoints; // Iloœæ punktów
	double* Tpoints = new double[NumberOfPoints];
	double* Tweights = new double[NumberOfPoints];
	integration(unsigned int Npoints) : NumberOfPoints(Npoints) {
		switch (NumberOfPoints) {
		case 1:
			Tpoints[0] = 0.;
			Tweights[0] = 2.;
			break;
		case 2:
			Tpoints[0] = -sqrt(1. / 3.);
			Tpoints[1] = sqrt(1. / 3.);
			Tweights[0] = 1.;
			Tweights[1] = 1.;
			break;
		case 3:
			Tpoints[0] = -sqrt(3. / 5.);
			Tpoints[1] = 0.;
			Tpoints[2] = sqrt(3. / 5.);
			Tweights[0] = 5. / 9.;
			Tweights[1] = 8. / 9.;
			Tweights[2] = 5. / 9.;
			break;
		case 4:
			Tpoints[0] = -sqrt(3. / 7. + (2. / 7.) * sqrt(6. / 5.));
			Tpoints[1] = -sqrt(3. / 7. - (2. / 7.) * sqrt(6. / 5.));
			Tpoints[2] = sqrt(3. / 7. - (2. / 7.) * sqrt(6. / 5.));
			Tpoints[3] = sqrt(3. / 7. + (2. / 7.) * sqrt(6. / 5.));
			Tweights[0] = (18. - sqrt(30.)) / 36.;
			Tweights[1] = (18. + sqrt(30.)) / 36.;
			Tweights[2] = (18. + sqrt(30.)) / 36.;
			Tweights[3] = (18. - sqrt(30.)) / 36.;
			break;
		default:
			break;
		}
	}
	~integration() {
		delete Tpoints;
		delete Tweights;
	}
	double Gauss() {
		double wynik = 0;
		if (Number_Of_Dimensions == 1)
			for (int i = 0; i < NumberOfPoints; i++)
				wynik += Tweights[i] * f1(Tpoints[i]);
		else if (Number_Of_Dimensions == 2)
			for (int i = 0; i < NumberOfPoints; i++)
				for (int j = 0; j < NumberOfPoints; j++)
					wynik += Tweights[i] * Tweights[j] * f2(Tpoints[i], Tpoints[j]);
		else
			return 0;
		return wynik;
	}
};

struct ElementUniw {
	unsigned int NumberOfPoints;
	vector<vector<double>> TKsi;
	vector<vector<double>> TEta;
	vector<vector<double>> TN;
	struct Surface{
		unsigned int npc;
		vector<vector<double>> N;
		double Hbc[4][4];
		Surface() {};
		Surface(unsigned int Npoints) : npc(Npoints) {
			N = vector<vector<double>>(npc, vector<double>(4));
		}
	};
	

	Surface surface[4];

	ElementUniw(unsigned int Npoints) : NumberOfPoints(Npoints) {
		for (int i = 0; i < 4; i++) {
			surface[i] = Surface(NumberOfPoints);
		}
		TKsi = vector<vector<double>>(NumberOfPoints * NumberOfPoints, vector<double>(4));
		TEta = vector<vector<double>>(NumberOfPoints * NumberOfPoints, vector<double>(4));
		TN = vector<vector<double>>(NumberOfPoints * NumberOfPoints, vector<double>(4));
		for (int i = 0; i < NumberOfPoints * NumberOfPoints; i++) {
			
			double Ksi, Eta;
			switch (NumberOfPoints) {
			case 1:
				Ksi = 0.;
				Eta = 0.;
				break;
			case 2:
				if (i % 2 == 0) {
					Ksi = -sqrt(1. / 3.);
				}
				else {
					Ksi = sqrt(1. / 3.);
				}
				if (i < 2) {
					Eta = -sqrt(1. / 3.);
				}
				else {
					Eta = sqrt(1. / 3.);
				}
				break;
			case 3:
				if (i % 3 == 0) {
					Ksi = -sqrt(3. / 5.);
				}
				else if (i % 3 == 1) {
					Ksi = 0.;
				}
				else {
					Ksi = sqrt(3. / 5.);
				}
				if (i < 3) {
					Eta = -sqrt(3. / 5.);
				}
				else if (i < 6) {
					Eta = 0.;
				}
				else {
					Eta = sqrt(3. / 5.);
				}
				break;
			case 4:
				if (i % 4 == 0) {
					Ksi = -sqrt(3. / 7. + (2. / 7.) * sqrt(6. / 5.));
				}
				else if (i % 4 == 1) {
					Ksi = -sqrt(3. / 7. - (2. / 7.) * sqrt(6. / 5.));
				}
				else if (i % 4 == 2) {
					Ksi = sqrt(3. / 7. - (2. / 7.) * sqrt(6. / 5.));
				}
				else {
					Ksi = sqrt(3. / 7. + (2. / 7.) * sqrt(6. / 5.));
				}
				if (i < 4) {
					Eta = -sqrt(3. / 7. + (2. / 7.) * sqrt(6. / 5.));
				}
				else if (i < 8) {
					Eta = -sqrt(3. / 7. - (2. / 7.) * sqrt(6. / 5.));
				}
				else if (i < 12) {
					Eta = sqrt(3. / 7. - (2. / 7.) * sqrt(6. / 5.));
				}
				else {
					Eta = sqrt(3. / 7. + (2. / 7.) * sqrt(6. / 5.));
				}
				break;
			default:
				break;
			}
			Deriv(Ksi, Eta, i);
		}
		Nsurface();
	}
	void Deriv(double Ksi, double Eta, int i) {
		for (int j = 0; j < 4; j++) {
			switch (j) {
			case 0:
				TKsi[i][j] = -0.25 * (1 - Eta);
				TEta[i][j] = -0.25 * (1 - Ksi);
				TN[i][j] = 0.25 * (1 - Ksi) * (1 - Eta);
				break;
			case 1:
				TKsi[i][j] = 0.25 * (1 - Eta);
				TEta[i][j] = -0.25 * (1 + Ksi);
				TN[i][j] = 0.25 * (1 + Ksi) * (1 - Eta);
				break;
			case 2:
				TKsi[i][j] = 0.25 * (1 + Eta);
				TEta[i][j] = 0.25 * (1 + Ksi);
				TN[i][j] = 0.25 * (1 + Ksi) * (1 + Eta);
				break;
			case 3:
				TKsi[i][j] = -0.25 * (1 + Eta);
				TEta[i][j] = 0.25 * (1 - Ksi);
				TN[i][j] = 0.25 * (1 - Ksi) * (1 + Eta);
				break;
			default:
				break;
			}
		}
	}
	void Nsurface() {
		integration integ(NumberOfPoints);
		for (int k = 0; k < 4; k++) {
			for (int i = 0; i < NumberOfPoints; i++) {
				double Ksi;
				double Eta;
				switch (k) {
				case 0:
					Ksi = integ.Tpoints[i];
					Eta = -1.;
					break;
				case 1:
					Ksi = 1.;
					Eta = integ.Tpoints[i];
					break;
				case 2:
					Ksi = integ.Tpoints[i];
					Eta = 1.;
					break;
				case 3:
					Ksi = -1.;
					Eta = integ.Tpoints[i];
					break;
				default:
					break;
				}
				for (int j = 0; j < 4; j++) {
					switch (j) {
					case 0:
						surface[k].N[i][j] = 0.25 * (1 - Ksi) * (1 - Eta);
						break;
					case 1:
						surface[k].N[i][j] = 0.25 * (1 + Ksi) * (1 - Eta);
						break;
					case 2:
						surface[k].N[i][j] = 0.25 * (1 + Ksi) * (1 + Eta);
						break;
					case 3:
						surface[k].N[i][j] = 0.25 * (1 - Ksi) * (1 + Eta);
						break;
					default:
						break;
					}
				}
			}
		}
	}
};

struct element {
	unsigned int NumberOfPoints;
	unsigned int ID[4];
	double H[4][4];
	double Hbc[4][4];
	double P[4];
	double C[4][4];
	element(unsigned int Npoints) : NumberOfPoints(Npoints) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				H[i][j] = 0.;
				Hbc[i][j] = 0.;
				C[i][j] = 0.;
			}
			P[i] = 0.;
		}
	}
	void fHbc(ElementUniw eu);
	void buildElement();
};


struct grid {
	unsigned int nN;
	unsigned int nE;
	vector<node> Tnode;
	vector<element> Telem;

	vector<vector<double>> HG;
	vector<double> PG;

	vector<vector<double>> CG;

	void buildGrid() {
		HG = vector<vector<double>>(nN, vector<double>(nN));
		PG = vector<double>(nN);
		CG = vector<vector<double>>(nN, vector<double>(nN));
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				HG[i][j] = 0.;
				CG[i][j] = 0.;
			}
			PG[i] = 0.;
		}

		for (int e = 0; e < nE; e++) {
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					HG[Telem[e].ID[i]][Telem[e].ID[j]] += Telem[e].H[i][j];
					CG[Telem[e].ID[i]][Telem[e].ID[j]] += Telem[e].C[i][j];
				}
				PG[Telem[e].ID[i]] += Telem[e].P[i];
			}
		}
	}
};

grid g;

void element::buildElement() {
	ElementUniw eu(NumberOfPoints);
	integration integ(NumberOfPoints);
	double TNx[4];
	double TNy[4];
	for (int p = 0; p < NumberOfPoints * NumberOfPoints; p++) {
		double pc[2][2] = { {0.,0.},{0.,0.} };

		// Punkty calkowania
		for (int i = 0; i < 4; i++) {
			pc[0][0] += eu.TKsi[p][i] * g.Tnode[ID[i]].x;
			pc[0][1] += eu.TKsi[p][i] * g.Tnode[ID[i]].y;
			pc[1][0] += eu.TEta[p][i] * g.Tnode[ID[i]].x;
			pc[1][1] += eu.TEta[p][i] * g.Tnode[ID[i]].y;
		}

		double detMat = pc[0][0] * pc[1][1] - pc[0][1]*pc[1][0];

		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 2; j++) {
				pc[i][j] = pc[i][j] / detMat;
				// Zamiana z | a b | na | d -b |
				//           | c d |    | -c a |
				if (i != j) {
					pc[i][j] = -pc[i][j];
				}
			}
		}
		double pc_switch = pc[0][0];
		pc[0][0] = pc[1][1];
		pc[1][1] = pc_switch;
		// Pochodne dNi/dx, dNi/dy
		for (int i = 0; i < 4; i++) {
			TNx[i] = pc[0][0] * eu.TKsi[p][i] + pc[0][1] * eu.TEta[p][i];
			TNy[i] = pc[1][0] * eu.TKsi[p][i] + pc[1][1] * eu.TEta[p][i];
		}	
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				// Macierz H
				H[i][j] += GD.Conductivity * (TNx[i] * TNx[j] + TNy[i] * TNy[j]) * detMat * integ.Tweights[p%NumberOfPoints] * integ.Tweights[p/NumberOfPoints]; 
				// Macierz C
				C[i][j] += GD.SpecificHeat * GD.Density * detMat * eu.TN[p][i] * eu.TN[p][j] * integ.Tweights[p % NumberOfPoints] * integ.Tweights[p / NumberOfPoints];
			}
		}
	}
	// Macierz Hbc
	fHbc(eu);
	// Dodanie macierzy Hbc do macierzy H
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			H[i][j] += Hbc[i][j];
		}
	}
}

void element::fHbc(ElementUniw eu) {
	integration integ(NumberOfPoints);
	for (int k = 0; k < 4; k++) {
		if (g.Tnode[ID[k]].BC == 1 && g.Tnode[ID[(k+1)%4]].BC == 1) {
			double L = pow((g.Tnode[ID[k]].x - g.Tnode[ID[(k + 1)%4]].x), 2) + pow((g.Tnode[ID[k]].y - g.Tnode[ID[(k + 1)%4]].y), 2);
			double det = sqrt(L) / 2.;
			for (int p = 0; p < NumberOfPoints; p++) {
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						Hbc[i][j] += GD.Alfa * det * eu.surface[k].N[p][i] * eu.surface[k].N[p][j] * integ.Tweights[p];
					}
					P[i] += GD.Alfa * det * GD.Tot * eu.surface[k].N[p][i] * integ.Tweights[p];
				}
			}
		}
	}
}

bool isStringEmpty(const std::string& str) {
	const std::regex pattern(R"(^s*$)");
	return std::regex_match(str, pattern);
}

void ReadFromFile(string FileName, int NumberOfPoints) {
	ifstream File(FileName);
	string line, str;
	vector<string> strArr;
	int x = 0;
	while (getline(File, line)) {
		stringstream ss(line);
		while (getline(ss, str, ' ')) {
			if (!isStringEmpty(str))
				strArr.push_back(str);
		}
	}
	int x1 = 0;

	GD.SimulationTime = stoi(strArr[1]);
	GD.SimulationStepTime = stoi(strArr[3]);
	GD.Conductivity = stoi(strArr[5]);
	GD.Alfa = stoi(strArr[7]);
	GD.Tot = stoi(strArr[9]);
	GD.InitialTemp = stoi(strArr[11]);
	GD.Density = stoi(strArr[13]);
	GD.SpecificHeat = stoi(strArr[15]);

	g.nN = stoi(strArr[18]);
	g.nE = stoi(strArr[21]);


	for (int i = 0; i < g.nN; i++) {
		int z = i * 2;
		strArr[i + 24 + z].pop_back();
		// Stworzenie obiektu wezla
		node n;
		n.x = stod(strArr[i + 24 + z]);
		n.y = stod(strArr[i + 25 + z]);
		// Wstawienie obiektu wezla do tablicy w siatce
		g.Tnode.push_back(n);
	}

	for (int i = 0; i < g.nE; i++) {
		int z = i * 4;
		int j = 26 + 3 * g.nN;
		strArr[i + j + z].pop_back();
		strArr[i + 1 + j + z].pop_back();
		strArr[i + 2 + j + z].pop_back();
		// Stworzenie obiektu elementu
		element e(NumberOfPoints);
		e.ID[0] = stoi(strArr[i + j + z]) - 1;
		e.ID[1] = stoi(strArr[i + j + z + 1]) - 1;
		e.ID[2] = stoi(strArr[i + j + z + 2]) - 1;
		e.ID[3] = stoi(strArr[i + j + z + 3]) - 1;
		// Wstawienie obiektu elementu do tablicy w siatce
		g.Telem.push_back(e);
	}

	int i = 0;
	int j = 26 + 3 * g.nN + 5 * g.nE;
	bool flag1 = false;
	while (!flag1) {
		if (strArr[i + j].find(',') != string::npos)
			strArr[i + j].pop_back();
		else
			flag1 = true;
		g.Tnode[stoi(strArr[i + j]) - 1].BC = 1;
		i++;
	}
	File.close();
}

void display(vector<double> matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		cout << "x" << i << " = " << matrix[i] << endl;
	}
}

void display(vector<vector<double>> matrix) {
	for (int i = 0; i < matrix.size(); i++) {
		for (int j = 0; j < matrix[0].size(); j++) {
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

void PosProste(vector<vector<double>>& matrix, double n1, double n2) {
	double mnoznik = 0;
	for (int k = 0; k < n1; k++) {
		for (int i = 1; i < n1; i++) {
			mnoznik = matrix[i][k] / matrix[k][k];
			for (int j = 0; j < n2; j++) {
				if (i > k)
					matrix[i][j] -= matrix[k][j] * mnoznik;
			}
		}
	}
}

vector<double> PosOdwrotne(vector<vector<double>> matrix, double n1, double n2) {
	vector<double> xk(n2-1);
	
	for (int i = n1 - 1; i >= 0; i--) {
		double suma = 0;
		for (int k = i + 1; k < n2 - 1; k++) {
			suma += matrix[i][k] * xk[k];
		}
		if (matrix[i][i] == 0)
			cout << "na przekatnej znajduje sie 0!" << endl;
		else
			xk[i] = (matrix[i][n2 - 1] - suma) / matrix[i][i];
	}
	return xk;
}

void ParaView(int fileNumber, vector<double> t0) {
	ofstream File("ParaView/Foo" + to_string(fileNumber) + ".vtk");
	File << "# vtk DataFile Version 2.0" << endl;
	File << "Unstructured Grid Example" << endl;
	File << "ASCII" << endl;
	File << "DATASET UNSTRUCTURED_GRID" << endl;
	File << endl;
	File << "POINTS "<< g.nN <<"  float" << endl;
	for (int i = 0; i < g.nN; i++) {
		File << g.Tnode[i].x << " " << g.Tnode[i].y << " 0" << endl;
	}
	File << endl;
	File << "CELLS " << g.nE << " " << 5 * g.nE << endl;
	for (int i = 0; i < g.nE; i++) {
		File << "4 " << g.Telem[i].ID[0] << " " << g.Telem[i].ID[1] << " " << g.Telem[i].ID[2] << " " << g.Telem[i].ID[3] << endl;
	}
	File << endl;
	File << "CELL_TYPES " << g.nE << endl;
	for (int i = 0; i < g.nE; i++) {
		File << 9 << endl;
	}
	File << endl;
	File << "POINT_DATA " << g.nN << endl;
	File << "SCALARS Temp float 1" << endl;
	File << "LOOKUP_TABLE default" << endl;
	for (int i = 0; i < g.nN; i++) {
		File << t0[i] << endl;
	}
	File.close();
}

void main() {
	int x = Number_Of_Points;
	//ReadFromFile("Test4_31_31_trapez.txt", x);
	//ReadFromFile("Test3_31_31_kwadrat.txt", x);
	ReadFromFile("Test2_4_4_MixGrid.txt", x);
	//ReadFromFile("Test1_4_4.txt", x);

	// Zbudowanie elementów i siatki, czyli obliczenie macierzy H, Hbc, C i wektora P
	for (int i = 0; i < g.Telem.size(); i++) {
		g.Telem[i].buildElement();
	}
	g.buildGrid();

	// Wypisanie obliczonych macierzy i wektora
	
	/*
	for (int k = 0; k < g.Telem.size(); k++) {
		cout << "Element " << k << ":" << endl;

		cout << "Macierz H:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << g.Telem[k].H[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;

		cout << "Macierz Hbc:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << g.Telem[k].Hbc[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;

		cout << "Macierz C:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << g.Telem[k].C[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;

		cout << "Wektor P:" << endl;
		for (int i = 0; i < 4; i++) {
			cout << g.Telem[k].P[i] << " ";
		}
		cout << endl << endl;
	}

	*/
	cout << "Macierz H globalna" << endl;
	//display(g.HG);
	cout << "Wektor P globalny" << endl;
	//display(g.PG);
	cout << "Macierz C globalna" << endl;
	//display(g.CG);
	
	// Obliczenie temperatury w wezlach

	for (int i = 0; i < g.nN; i++) {
		for (int j = 0; j < g.nN; j++) {
			g.HG[i][j] += g.CG[i][j] / GD.SimulationStepTime;
		}
	}
	
	vector<double> Ptemp = g.PG;
	for (int i = 0; i < g.nN; i++) {
		for (int j = 0; j < g.nN; j++) {
			g.PG[i] += (g.CG[i][j] / GD.SimulationStepTime) * GD.InitialTemp;
		}
	}
	
	cout << "Macierz globalna [H]+[C]_dT:" << endl;
	//display(g.HG);
	cout << "Wektor globalny {P}+[C]_dT*{t0}:" << endl;
	//display(g.PG);

	for (int i = 0; i < g.nN; i++) {
		g.HG[i].push_back(g.PG[i]);
	}
	vector<vector<double>> Htemp = g.HG;
	vector<vector<double>> Htemp2 = g.HG;
	PosProste(Htemp, g.nN, g.nN + 1);
	
	cout << endl << "Rozwiazanie ukladu rownañ dla czasu " << GD.SimulationStepTime << ", temperatura: " << endl;
	vector<double> t0 = PosOdwrotne(Htemp, g.nN, g.nN + 1);
	//display(t0);
	int FileNumber = 0;
	ParaView(FileNumber, t0);
	for (int timeNow = 2*GD.SimulationStepTime; timeNow <= GD.SimulationTime; timeNow+= GD.SimulationStepTime) {
		Htemp = Htemp2;
		for (int i = 0; i < g.nN; i++) {
			Htemp[i][g.nN] = Ptemp[i];
			for (int j = 0; j < g.nN; j++) {
				Htemp[i][g.nN] += (g.CG[i][j] / GD.SimulationStepTime) * t0[j];
			}	 
		}

		PosProste(Htemp, g.nN, g.nN + 1);
		t0 = PosOdwrotne(Htemp, g.nN, g.nN + 1);
		cout << endl << "Rozwiazanie ukladu rownañ dla czasu " << timeNow << ", temperatura: " << endl;
		FileNumber++;
		display(t0);
		ParaView(FileNumber, t0);
	}
}