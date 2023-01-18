#pragma once
#include "LU solver/headers/ProfileLU.h"
#include "json.hpp"
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include <format>

using json = nlohmann::json;

namespace Newtons {

   
   namespace Vec {
      inline double Scalar(const std::vector<double>& l, const std::vector<double>& r);

      inline double Norm(const std::vector<double>& vec);

      // ans = left + coef * right
      inline void AddVec(
         const std::vector<double>& left,
         double coef,
         const std::vector<double>& right,
         std::vector<double>& ans) 
      {
         for (size_t i = 0; i < right.size(); i++)
         {
            ans[i] = left[i] + coef * right[i];
         }
      }
   }


   class NewtonsSolver {
   public:

      struct TraceElement {
         int iterationNum{};
         std::vector<double> prevX;
         std::vector<double> X;
         std::vector<double> dX;
         double eps{};
         double prevEps{};

         TraceElement() noexcept {}

         TraceElement(const TraceElement& elem) noexcept {
            iterationNum = elem.iterationNum;
            prevX = (elem.prevX);
            X = (elem.X);
            dX = (elem.dX);
            eps = elem.eps;
            prevEps = elem.prevEps;
         }

         TraceElement(TraceElement&& elem) noexcept {
            iterationNum = elem.iterationNum;
            prevX = std::move(elem.prevX);
            X = std::move(elem.X);
            dX = std::move(elem.dX);
            eps = elem.eps;
            prevEps = elem.prevEps;
         }
      };

      class TraceVector {
      private:
         std::vector<TraceElement> _traceVec;

      public:
         TraceElement& operator[](size_t ind) {
            return _traceVec[ind];
         }

         void Push(TraceElement&& elem) {
            _traceVec.push_back(elem);
         }
         void Push(TraceElement& elem) {
            _traceVec.push_back(elem);
         }

         void Clear() {
            _traceVec.clear();
         }

         json ToJson() {
            json js;
            std::vector<json> semi_jsons;
            semi_jsons.resize(_traceVec.size());
            for (size_t i = 0; i < _traceVec.size(); i++)
            {
               semi_jsons[i]["iterNum"] = _traceVec[i].iterationNum;
               semi_jsons[i]["prevX"] = _traceVec[i].prevX;
               semi_jsons[i]["X"] = _traceVec[i].X;
               semi_jsons[i]["dX"] = _traceVec[i].dX;
               semi_jsons[i]["eps"] = _traceVec[i].eps;
               semi_jsons[i]["prevEps"] = _traceVec[i].prevEps;
            }
            js = semi_jsons;

            return js;
         }
      };

   private:

      // ���������� ��� ������ � ��������, ������������ � �������
      ProfileMatrix _profMat;
      Matrix _mat;
      std::vector<double> _F;
      std::vector<double> _x;
      std::vector<double> _dx;
      std::vector<double> _dx_trim;

      // ���������� ��� �������� ������� � �� ��������������
      std::function<double(std::size_t, const std::vector<double>&)> _functions;
      std::function<double(std::size_t, const std::size_t, std::vector<double>&)> _differentials;

      // ���������� ��� �������� ���������� ������� � ����������
      size_t _funcCount;
      size_t _varCount;

      // ������ ��� ������������ �����
      std::vector<std::pair<size_t, double>> _pairVec;

      // ��� ����������� ������� �����
      enum class MaskType {
         None,
         MoreVars,
         MoreFuncs
      };
      MaskType _maskType = MaskType::None;

      // ����� ��� ���������� ������ ���������� ��� �������
      std::vector<bool> _mask;

      // ��������� �� ������ ��� ����������� ������ (��������� ���������� ���������� �� ������ ����)
      TraceVector* _traceVector;

   public:

      /// <summary>
      /// ������������� ������� ������� �������
      /// </summary>
      /// <param name="variableCount"> - ���������� ���������� � �������</param>
      /// <param name="funcCount"> - ���������� ������� � �������</param>
      /// <param name="functions"> - �������, � ������� ������ ��� ������� F_j �������</param>
      /// <param name="differentials"> - �������, � ������� ������ ��� ������������� �������. </param>
      NewtonsSolver(
         size_t variableCount,
         size_t funcCount,
         std::function<double(size_t, const std::vector<double>&)> functions,
         std::function<double(size_t, size_t, const std::vector<double>&)> differentials)
      {
         _varCount = variableCount;
         _funcCount = funcCount;

         _x.resize(variableCount);
         _dx.resize(variableCount);

         _functions = functions;
         _differentials = differentials;

         size_t minSize = std::min(variableCount, funcCount);
         _mat.resize(minSize, minSize);
         _F.resize(minSize);

         if (variableCount != funcCount)
         {
            if (variableCount > funcCount)
            {
               _dx_trim.resize(funcCount);
            }

            _mask.resize(std::max(variableCount, funcCount));
            _pairVec.resize(std::max(variableCount, funcCount));
         }
      }


   private:

      // ���������� ��� ����� � �������� �����, ������ _mask � _maskType
      inline void _GetMask();

      // ������� ������� ����� � ������ ����� � ���������� � � _mat
      void _GetJacobi();

      // ������� ������ ������� F ��� ������� �������
      void _GetF();

      // ������� ����� ������� �� �������� ������ F � ����� x
      double _GetNormF(const std::vector<double>& x) {
         double res = 0;
         double t;
         for (size_t i = 0; i < _funcCount; i++)
         {
            t = _functions(i, x);
            res += t * t;
         }
         return std::sqrt(res);
      }

   public:

      // ����������� �������� ������� ������� �������
      double minEps = 1e-5;

      // ������������ ����� ��������
      int maxIter = 100;

      // ����������� �����������, ����� �������� ����� ����������� � ������� ����������
      // �� ��������� ����� �����, ���������������� 6 ���������� ������������ �� 2
      double criticalCoef = 1 / (1 << 6);

      // ����� ��� ������� ������� ���������� ���������
      // - init_x - ��������� �����������, � ��� ����� �������� �������
      // - eps - ���������� ������� �������
      // ��������� �������:
      // -1 - ������ ���������� (��� ����� ��������� ��������� ���� ������� ����������)
      // -2 - ����� �� ���������� ����� ��������
      // ������������� ����� - ����� �������� ���������� ������
      int Solve(std::vector<double>& init_x, double& eps, const bool debugOutput = false);

      class TraceVector;

      void EnableTracing(TraceVector& traceVector) {
         _traceVector = &traceVector;
         traceVector.Clear();
      }
   };


}