#include <iostream>
#include <string>
#include <fstream>


class traject {public: std::string flmp; int ndata; std::string fxyz; std::string fener; double ww;};

int main()
{
  int num_file, ndata;
  std::string flmp, fxyz, fener;
  int i;
  double w;
  std::ifstream fin;

  fin.open("input_datas_all");

  fin >> num_file ;
  std::cout << num_file << '\n';

  traject inputs[num_file];
  
  for (int i = 0; i < num_file; i++)
    {
      fin >> inputs[i].flmp >> inputs[i].ndata >> inputs[i].fxyz >> inputs[i].fener >> inputs[i].ww;
    }
  fin.close();
  
  for (int i = 0; i < num_file; i++)
    {
      std::cout << inputs[i].flmp << '\t' << inputs[i].ndata << '\t' << inputs[i].fxyz << '\t' << inputs[i].fener << '\t' << inputs[i].ww << '\n';
    }
}
// END PROGRAM MAIN
