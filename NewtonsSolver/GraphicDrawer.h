#pragma once
#include <SFML/Graphics.hpp>


namespace Newtons {


   using namespace sf;


   std::vector<sf::Color> functionColors = {
      sf::Color::Blue,
      sf::Color::Green,
      sf::Color::Magenta,
      sf::Color::Cyan,
      sf::Color::Yellow,
   };

   inline sf::Vector2f GetCoordinatesXY(size_t i, size_t j, sf::Vector2f midPoint, float scale) {
      sf::Vector2f point;
      point.x = (i - midPoint.x) * scale;
      point.y = -(j - midPoint.y) * scale;

      return point;
   }

   inline sf::Vector2f GetCoordinatesIJ(sf::Vector2f xy, sf::Vector2f midPoint, float scale) {
      sf::Vector2f point;
      point.x = (xy.x) / scale + midPoint.x;
      point.y = -(xy.y) / scale + midPoint.y;

      return point;
   }

   inline double GetNormF(double x, double y, size_t funcCount, std::function<double(size_t, const std::vector<double>&)> functions) {
      double res = 0;
      double t;
      for (size_t i = 0; i < funcCount; i++)
      {
         t = functions(i, { x,y });
         res += t * t;
      }
      return std::sqrt(res);
   }



   class GraphicDrawer {

   private:
      std::size_t _width;
      std::size_t _height;
      Vector2f _midPoint;
      std::wstring _title;
      float _scale = 1;

   public:

      GraphicDrawer(std::size_t width, std::size_t height, std::wstring title, float scale) {
         _scale = scale;
         _width = width;
         _height = height;
         _title = title;
         _midPoint = Vector2f(width / 2.0f, height / 2.0f);

         window.create(VideoMode(_width, _height), _title);
         //window.setVerticalSyncEnabled(true);
         window.clear(Color::White);
      }

   public:
      RenderWindow window;

      void DrawCoodrLines(size_t xDivisorsCount, size_t yDivisorsCount) {
         // Отрисовка главных осей
         Vertex xLine[] = {
            Vertex(Vector2f(_midPoint.x, 0.f), Color::Black),
            Vertex(Vector2f(_midPoint.x, window.getSize().y * 1.f), Color::Black)
         };
         Vertex yLine[] = {
            Vertex(Vector2f(0.f, _midPoint.y), Color::Black),
            Vertex(Vector2f(window.getSize().x * 1.f, _midPoint.y), Color::Black)
         };
         
         window.draw(xLine, 2, sf::Lines);
         window.draw(yLine, 2, sf::Lines);

         // крупные разделители, каждые 2 значения X и Y
         if (xDivisorsCount % 2 == 0) xDivisorsCount++;
         if (yDivisorsCount % 2 == 0) yDivisorsCount++;

         Font font;
         if (!font.loadFromFile("Kelvinch-Roman.otf"))
         {
            throw std::runtime_error("Error when loading font");
         }
         Text text;
         text.setFont(font);
         text.setFillColor(Color::Black);
         text.setCharacterSize(15);

         int xStep = (_width / xDivisorsCount);
         for (int i = -(static_cast<int>(xDivisorsCount) / 2); i <= static_cast<int>(xDivisorsCount) / 2; i++)
         {
            Vertex line[] = {
               Vertex(Vector2f(_midPoint.x + i * xStep, _midPoint.y - 6), Color::Black),
               Vertex(Vector2f(_midPoint.x + i * xStep, _midPoint.y + 6), Color::Black)
            };
            window.draw(line, 2, sf::Lines);
            text.setPosition(Vector2f(_midPoint.x + i * xStep - 15, _midPoint.y + 10));
            text.setString(std::format("{: 4.1f}", GetCoordinatesXY(_midPoint.x + i * xStep, 0, _midPoint, _scale).x));
            window.draw(text);
         }

         int yStep = (_height / yDivisorsCount);
         for (int i = -(static_cast<int>(yDivisorsCount) / 2); i <= static_cast<int>(yDivisorsCount) / 2; i++)
         {
            if (i == 0) continue;

            Vertex line[] = {
               Vertex(Vector2f(_midPoint.x - 6, _midPoint.y + i * yStep), Color::Black),
               Vertex(Vector2f(_midPoint.x + 6, _midPoint.y + i * yStep), Color::Black)
            };
            window.draw(line, 2, sf::Lines);
            text.setPosition(Vector2f(_midPoint.x + 10, _midPoint.y + i * yStep - 10));
            text.setString(std::format("{: 4.1f}", GetCoordinatesXY(0, _midPoint.y + i * yStep, _midPoint, _scale).y));
            window.draw(text);
         }
      }

      void DrawHeatMap(VertexArray& points, std::size_t funcCount, std::function<double(std::size_t, const std::vector<double>&)> F) {
         for (std::size_t i = 0; i < _width; i++)
         {
            for (std::size_t j = 0; j < _height; j++)
            {
               auto point = GetCoordinatesXY(i, j, _midPoint, _scale);
               double norm = GetNormF(point.x, point.y, funcCount, F);
               if (norm >= 25)
               {
                  points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), sf::Color::White);
               }
               else if (norm > 15.0)
               {
                  points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), sf::Color(255, 0, 0, 50));
               }
               else if (norm > 5.0)
               {
                  points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), sf::Color(255, 0, 0, 100));
               }
               else if (norm > 1)
               {
                  points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), sf::Color(255, 0, 0, 150));
               }
               else if (norm > 0.1)
               {
                  points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), sf::Color(255, 0, 0, 200));
               }
               else
               {
                  points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), sf::Color(255, 0, 0, 230));
               }
            }
         }
      }

      void DrawGraphics(VertexArray& points, float nearToFunc, std::size_t funcCount, std::function<double(std::size_t, const std::vector<double>&)> F) {
         for (std::size_t i = 0; i < _width; i++)
         {
            for (std::size_t j = 0; j < _height; j++)
            {
               auto point = Newtons::GetCoordinatesXY(i, j, _midPoint, _scale);
               for (int fn = 0; fn < funcCount; fn++)
               {
                  if (abs(F(fn, { point.x, point.y })) < nearToFunc)
                  {
                     points[i * _height + j] = sf::Vertex(sf::Vector2f(i, j), functionColors[fn]);
                     break;
                  }
               }
            }
         }
      }

      void DrawAll(size_t xDivCount, size_t yDivCount, float nearToFunc, std::size_t funcCount, std::function<double(std::size_t, const std::vector<double>&)> F) {
         VertexArray points(Points, _width * _height);

         DrawHeatMap(points, funcCount, F);
         DrawGraphics(points, nearToFunc, funcCount, F);
         window.draw(points);

         DrawCoodrLines(xDivCount, yDivCount);
      }

      void AwaitCloseSync() {
         while (window.isOpen())
         {
            Event event;
            while (window.waitEvent(event))
            {
               if (event.type == Event::Closed)
               {
                  window.close();
               }
            }

         }
      }

   };
}