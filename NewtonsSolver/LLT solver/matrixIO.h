#pragma once
#include <fstream>
#include <iostream>

/// <summary>
/// Функция получения числа типа size_t из файла
/// </summary>
/// <param name="filePath"> - путь до файла </param>
/// <returns> число из файла </returns>
inline size_t GetSizetFromFile(const char* filePath)
{
   auto fin = std::ifstream(filePath);
   size_t matrix_size = 0;
   fin >> matrix_size;
   fin.close();
   return matrix_size;
}

/// <summary>
/// Функция считывания массива данных с файла
/// </summary>
/// <typeparam name="T"> - тип данных, который необходимо считать </typeparam>
/// <param name="filePath"> - путь до файла</param>
/// <param name="count"> - число элементов, необходимое считать </param>
/// <returns> указатель на начало считанного массива </returns>
template<typename T>
T* GetArrayFromFile(const char* filePath, const size_t count)
{
   T* arr = new T[count];
   auto fin = std::ifstream(filePath);
   for (int i = 0; i < count; i++)
   {
      fin >> arr[i];
   }
   fin.close();

   return arr;
}

/// <summary>
/// Вывод массива в произвольный поток
/// </summary>
/// <typeparam name="T"> - тип данных для вывода</typeparam>
/// <param name="_arr"> - указатель на начало массива </param>
/// <param name="_size"> - длина массива </param>
/// <param name="_out">  - поток вывода (cout по умолчанию)</param>
template <typename T>
void PrintArray(const T* _arr, const size_t _size, const int precision = 7, std::ostream& _out = std::cout)
{
   auto prevPrec = _out.precision();
   auto prevFlags = _out.flags();
   _out.precision(precision);
   _out.setf(std::ios_base::fixed);

   for (int i = 0; i < _size; i++) {
      _out << _arr[i] << std::endl;
   }

   _out.precision(prevPrec);
   _out.flags(prevFlags);
}
