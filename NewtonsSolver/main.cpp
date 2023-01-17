#include <iostream>
#include <fstream>
#include <vector>


using namespace std;

/// <summary>
/// Получение функций F_i
/// </summary>
/// <param name="funcNum"> - номер функции для вызова;</param>
/// <param name="vars"> - параметры функции</param>
/// <returns>Значение функции в точке [vars]</returns>
double F(int funcNum, const vector<double>& vars) {
   switch (funcNum)
   {
      //case 0: return vars[0] + vars[1];
      //case 1: return vars[1];

      default: throw runtime_error("Неправильный номер функции F");
   }
}

/// <summary>
/// Получение производных функций F_i
/// </summary>
/// <param name="funcNum"> - номер функции для выбора производной;</param>
/// <param name="varNum"> - номер параметра, по которому берётся производная;</param>
/// <param name="vars"> - параметры функции</param>
/// <returns>Значение функции в точке [vars]</returns>
double div_F(int funcNum, int varNum, const vector<double>& vars) {
   switch (funcNum)
   {
      //case 0:
      //{
      //   switch (varNum)
      //   {
      //      case 0:
      //         return 1;
      //      case 1:
      //         return 1;

      //      default: throw runtime_error("Неправильный номер параметра функции");
      //   }
      //}


      //case 1:
      //{
      //   switch (varNum)
      //   {
      //      case 0:
      //         return 0;
      //      case 1:
      //         return 1;

      //      default: throw runtime_error("Неправильный номер параметра функции");
      //   }
      //}

      default: throw runtime_error("Неправильный номер функции F");
   }
}

int main() {


}