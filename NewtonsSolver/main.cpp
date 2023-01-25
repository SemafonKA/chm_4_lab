#include <iostream>
#include <fstream>
#include <vector>
#include "NewtonsSolver.h"
#include "GraphicDrawer.h"

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
      case 0: return pow(x[0] - 2, 2) + pow(x[1], 2) - 4;
      case 1: return pow(x[0] + 2, 2) + pow(x[1], 2) - 4;
      //case 2: return x[0];
      // 
      //case 0: return sin(x[0])*sin(x[0]) - x[1];
      //case 1: return 2*exp(x[0]) - x[1] - 5;
      //case 2: return x[0] * x[0] + pow((x[1] - 2), 2) - 4;

      //case 0: return x[1] + x[0] - 4;
      //case 1: return x[1] - 2 * x[0];
      //case 2: return x[1] - x[0] + 4;

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
               return 2*sin(x[0])*cos(x[0]);
            case 1:
               return -1;

            default: throw std::runtime_error("Неправильный номер параметра функции");
         }
      }


      case 1:
      {
         switch (varNum)
         {
            case 0:
               return 2*exp(x[0]);
            case 1:
               return -1;

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
   static double step = 1e-7;
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

   auto solver = NewtonsSolver(varCount, funcCount, F, div_F_numeric);

   NewtonsSolver::TraceVector traceVector;
   solver.EnableTracing(traceVector);

   double eps;
   vector<double> x = { 4.0, 1.0 };

   int a = solver.Solve(x, eps, true);

   cout << "Полученное решение:";
   for (auto& el : x)
   {
      cout << format("{:23.15e}", el);
   }
   cout << "\n\n";


   if (varCount == 2)
   {
      constexpr size_t width = 800, height = 600;
      constexpr float scale = 0.035f;
      constexpr float nearToFunc = 0.1f;

      sf::Vector2f midPoint{ width / 2.0, height / 2.0 };

      Newtons::GraphicDrawer drawer(width, height, L"Графики функций и движения метода", scale);

      // Рисует функции, тепловое поле невязки и оси координат
      drawer.DrawAll(18, 15, nearToFunc, funcCount, F);

      // Рисует движение точки х в процессе решения метода
      for (size_t i = 0; i < traceVector.Size(); i++)
      {
         auto beginPoint = Newtons::GetCoordinatesIJ(sf::Vector2f(traceVector[i].prevX[0], traceVector[i].prevX[1]), midPoint, scale);
         auto endPoint = Newtons::GetCoordinatesIJ(sf::Vector2f(traceVector[i].X[0], traceVector[i].X[1]), midPoint, scale);

         sf::Vertex line[] = {
            sf::Vertex(beginPoint, sf::Color::Blue),
            sf::Vertex(endPoint, sf::Color::Blue)
         };
         drawer.window.draw(line, 2, sf::Lines);
      }

      // Рисует точку начала движения метода
      sf::CircleShape startPoint(4);
      auto startCirclePosition = Newtons::GetCoordinatesIJ(sf::Vector2f(traceVector[0].prevX[0], traceVector[0].prevX[1]), midPoint, scale);
      startCirclePosition.x -= 4;
      startCirclePosition.y -= 4;
      startPoint.setPosition(startCirclePosition);
      startPoint.setFillColor(sf::Color::Blue);
      drawer.window.draw(startPoint);

      // Рисует точку конца движения метода
      sf::CircleShape endPoint(4);
      auto endCirclePosition = Newtons::GetCoordinatesIJ(sf::Vector2f(traceVector[traceVector.Size()-1].X[0], traceVector[traceVector.Size() - 1].X[1]), midPoint, scale);
      endCirclePosition.x -= 4;
      endCirclePosition.y -= 4;
      endPoint.setPosition(endCirclePosition);
      endPoint.setFillColor(sf::Color::Blue);
      drawer.window.draw(endPoint);

      drawer.window.display();
      drawer.AwaitCloseSync();
   }
}