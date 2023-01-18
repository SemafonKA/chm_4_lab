#include <iostream>
#include <fstream>
#include <vector>
#include "NewtonsSolver.h"

using namespace std;
using NewtonsSolver = Newtons::NewtonsSolver;

/// <summary>
/// Получение функций F_i
/// </summary>
/// <param name="funcNum"> - номер функции для вызова;</param>
/// <param name="x"> - параметры функции</param>
/// <returns>Значение функции в точке [x]</returns>
double F(size_t funcNum, const std::vector<double>& x) {
   switch (funcNum)
   {
      case 0: return pow((x[0] - 2), 2) + x[1] * x[1] - 4;
      case 1: return pow((x[0] + 2), 2) + x[1] * x[1] - 4;
      //case 2: return x[0] * x[0] + pow((x[1] - 2), 2) - 4;

      default: throw std::runtime_error("Неправильный номер функции F");
   }
}

/// <summary>
/// Получение производных функций F_i
/// </summary>
/// <param name="funcNum"> - номер функции для выбора производной;</param>
/// <param name="varNum"> - номер параметра, по которому берётся производная;</param>
/// <param name="x"> - параметры функции</param>
/// <returns>Значение функции в точке [x]</returns>
double div_F(size_t funcNum, size_t varNum, const std::vector<double>& x) {
   switch (funcNum)
   {
      case 0:
      {
         switch (varNum)
         {
            case 0:
               return 2*x[0] - 4;
            case 1:
               return 2*x[1];

            default: throw std::runtime_error("Неправильный номер параметра функции");
         }
      }


      case 1:
      {
         switch (varNum)
         {
            case 0:
               return 2 * x[0] + 4;
            case 1:
               return 2 * x[1];

            default: throw std::runtime_error("Неправильный номер параметра функции");
         }
      }

      //case 2:
      //{
      //   switch (varNum)
      //   {
      //      case 0:
      //         return 2 * x[0];
      //      case 1:
      //         return 2 * x[1] - 4;

      //      default: throw std::runtime_error("Неправильный номер параметра функции");
      //   }
      //}

      default: throw std::runtime_error("Неправильный номер функции F");
   }
}

double div_F_numeric(size_t funcNum, size_t varNum, const std::vector<double>& vars) {
   static double step = 1e-6;
   static vector<double> r;
   static vector<double> l;
   r = vars;
   l = vars;
   r[varNum] += step;
   l[varNum] -= step;

   return (F(funcNum, r) - F(funcNum, l)) / (2.0 * step);
}

int main() {
   setlocale(LC_ALL, "ru-RU");

   constexpr size_t varCount = 2;
   constexpr size_t funcCount = 2;

   auto solver = NewtonsSolver(varCount, funcCount, F, div_F);

   NewtonsSolver::TraceVector traceVector;
   solver.EnableTracing(traceVector);

   double eps;
   vector<double> x = { 25.0, -4.0 };

   int a = solver.Solve(x, eps, true);

   cout << "Полученное решение:";
   for (auto& el : x)
   {
      cout << format("{:23.15e}", el);
   }
   cout << "\n\n";

   ofstream outFile;
   outFile.open("iofiles/jsonTrace.json");
   if (!outFile.is_open())
      throw runtime_error("OutFile can't be open!");
   outFile << traceVector.ToJson().dump(2);
   outFile.close();
}