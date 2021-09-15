#include "Algorithm/ConvexHull.h"

#include "BuildingBlocks/Api.h"
#include "BuildingBlocks/Structures.h"
#include "BuildingBlocks/Enums.h"

#include <algorithm>
#include <vector>

#define _USE_MATH_DEFINES
#include <math.h>
#include <cstring>

namespace BB
{
   namespace Dim2
   {
      namespace GrahamScan
      {
         /// <summary>
         /// <algorithm name="Graham's Scan Convex Hull"></algorithm> 
         /// <runtime complexity="n*logn"></runtime>
         /// </summary>
         /// <param name="vertices"></param>
         /// <param name="convexHull"></param>
         void ConvexHull(const std::vector<Vertex>& vertices, std::vector<Vertex>& convexHull)
         {
            if (vertices.size() > 3)
            {
               // Sort helper expression
               static auto sortByXThenY = [](const Vertex& lhs, const Vertex& rhs)
               {
                  return (lhs.x < rhs.x) || (lhs.x == rhs.x && lhs.y < rhs.y);
               };

               // Sort Vertices Left to Right
               std::vector<Vertex> sortedVertices = vertices;
               std::sort(sortedVertices.begin(), sortedVertices.end(), sortByXThenY);

               convexHull.push_back(sortedVertices[0]);
               convexHull.push_back(sortedVertices[1]);

               for (int i = 2; i < sortedVertices.size(); i++)
               {
                  auto& vertex = sortedVertices[i];
                  while (convexHull.size() >= 2 && Orientation(vertex, convexHull[convexHull.size() - 1], convexHull[convexHull.size() - 2]) != ORIENTATION_CW)
                  {
                     convexHull.pop_back();
                  }
                  convexHull.push_back(vertex);
               }

               int numVerticesInLowerHull = static_cast<int>(convexHull.size());
               convexHull.push_back(sortedVertices[sortedVertices.size() - 2]);
               for (int i = sortedVertices.size() - 3; i >= 0; i--)
               {
                  auto& vertex = sortedVertices[i];
                  while (static_cast<int>(convexHull.size()) >= static_cast<int>(numVerticesInLowerHull) + 1 && Orientation(vertex, convexHull[convexHull.size() - 1], convexHull[convexHull.size() - 2]) != ORIENTATION_CW)
                  {
                     convexHull.pop_back();
                  }

                  if (i > 0)
                  {
                     convexHull.push_back(vertex);
                  }
               }
            }
            else
            {
               convexHull = vertices;
            }
         }
      }

      namespace GiftWrap
      {
         /// <summary>
         /// <algorithm name="Gift Wrap Convex Hull"></algorithm> 
         /// <runtime complexity="n²"></runtime>
         /// </summary>
         /// <param name="vertices"></param>
         /// <param name="convexHull"></param>
         void ConvexHull(const std::vector<Vertex>& vertices, std::vector<Vertex>& convexHull)
         {
            if (vertices.size() > 3)
            {
               // Sort helper expression
               static auto sortByXThenY = [](const Vertex& lhs, const Vertex& rhs)
               {
                  return (lhs.x < rhs.x) || (lhs.x == rhs.x && lhs.y < rhs.y);
               };

               // Find left most point
               Vertex leftMostPoint = vertices[0];
               for (auto& vertex : vertices)
               {
                  if (sortByXThenY(vertex, leftMostPoint))
                  {
                     leftMostPoint = vertex;
                  }
               }

               Vertex endPoint = leftMostPoint;
               do
               {
                  convexHull.push_back(endPoint);
                  endPoint = (endPoint == vertices[0]) ? vertices[1] : vertices[0];
                  for (int i = 0; i < static_cast<int>(vertices.size()); i++)
                  {
                     if ((vertices[i] != endPoint) && Orientation(convexHull.back(), endPoint, vertices[i]) == ORIENTATION_CW)
                     {
                        endPoint = vertices[i];
                     }
                  }
               } while (convexHull[0] != endPoint);
            }
            else
            {
               convexHull = vertices;
            }
         }
      }

      namespace MergeHull
      {
         /// <summary>
         /// <algorithm name="Merge Hull"></algorithm> 
         /// <runtime complexity="n*logn"></runtime>
         /// </summary>
         /// <param name="vertices"></param>
         /// <param name="convexHull"></param>
         void MergeHull(const std::vector<Vertex>& lhs, const std::vector<Vertex>& rhs, std::vector<Vertex>& convexHull)
         {
            // Sort helper expression
            static auto sortByXThenY = [](const Vertex& lhs, const Vertex& rhs)
            {
               return (lhs.x < rhs.x) || (lhs.x == rhs.x && lhs.y < rhs.y);
            };

            // Find left most point of rhs
            int indexLeftMostPointRhs = 0;
            for (int i = 0; i < static_cast<int>(rhs.size()); i++)
            {
               if (sortByXThenY(rhs[i], rhs[indexLeftMostPointRhs]))
               {
                  indexLeftMostPointRhs = i;
               }
            }

            // Find right most point of lhs
            int indexRightMostPointLhs = 0;
            for (int i = 0; i < static_cast<int>(lhs.size()); i++)
            {
               if (!sortByXThenY(lhs[i], lhs[indexRightMostPointLhs]))
               {
                  indexRightMostPointLhs = i;
               }
            }

            // Find upper tangent
            int lhsUpperTangentPoint = indexRightMostPointLhs;
            int rhsUpperTangentPoint = indexLeftMostPointRhs;
            int prevLhsUpperTangentPoint = lhsUpperTangentPoint;
            int prevRhsUpperTangentPoint = rhsUpperTangentPoint;
            do
            {
               prevLhsUpperTangentPoint = lhsUpperTangentPoint;
               prevRhsUpperTangentPoint = rhsUpperTangentPoint;
               while (Orientation(lhs[lhsUpperTangentPoint], rhs[rhsUpperTangentPoint], rhs[(rhs.size() + rhsUpperTangentPoint - 1) % rhs.size()]) == ORIENTATION_CCW)
               {
                  rhsUpperTangentPoint = (rhsUpperTangentPoint <= 0) ? static_cast<int>(rhs.size() - 1) : rhsUpperTangentPoint - 1;
               }
               while (Orientation(rhs[rhsUpperTangentPoint], lhs[lhsUpperTangentPoint], lhs[(static_cast<long long>(lhsUpperTangentPoint) + 1) % lhs.size()]) == ORIENTATION_CW)
               {
                  lhsUpperTangentPoint = (lhsUpperTangentPoint >= lhs.size()) ? 0 : lhsUpperTangentPoint + 1;
               }
            } while(!(lhsUpperTangentPoint == prevLhsUpperTangentPoint && rhsUpperTangentPoint == prevRhsUpperTangentPoint));

            // Find lower tangent
            int lhsLowerTangentPoint = indexRightMostPointLhs;
            int rhsLowerTangentPoint = indexLeftMostPointRhs;
            int prevLhsLowerTangentPoint = lhsLowerTangentPoint;
            int prevRhsLowerTangentPoint = rhsLowerTangentPoint;
            do
            {
               prevLhsLowerTangentPoint = lhsLowerTangentPoint;
               prevRhsLowerTangentPoint = rhsLowerTangentPoint;
               while (Orientation(lhs[lhsLowerTangentPoint], rhs[rhsLowerTangentPoint], rhs[(static_cast<long long>(rhsLowerTangentPoint) + 1) % rhs.size()]) == ORIENTATION_CW)
               {
                  rhsLowerTangentPoint = (rhsLowerTangentPoint >= rhs.size()) ? 0 : rhsLowerTangentPoint + 1;
               }
               while (Orientation(rhs[rhsLowerTangentPoint], lhs[lhsLowerTangentPoint], lhs[(lhs.size() + lhsLowerTangentPoint - 1) % lhs.size()]) == ORIENTATION_CCW)
               {
                  lhsLowerTangentPoint = (lhsLowerTangentPoint <= 0) ? static_cast<int>(lhs.size() - 1) : lhsLowerTangentPoint - 1;
               }
            } while (!(lhsLowerTangentPoint == prevLhsLowerTangentPoint && rhsLowerTangentPoint == prevRhsLowerTangentPoint));

            // Find left most point of lhs
            int indexLeftMostPointLhs = 0;
            for (int i = 0; i < static_cast<int>(lhs.size()); i++)
            {
               if (sortByXThenY(lhs[i], lhs[indexLeftMostPointLhs]))
               {
                  indexLeftMostPointLhs = i;
               }
            }

            // Combine Convex Hulls
            for (int i = indexLeftMostPointLhs; i != lhsLowerTangentPoint;)
            {
               convexHull.push_back(lhs[i]);
               i = (i == lhs.size() - 1) ? 0 : i + 1;
            }
            convexHull.push_back(lhs[lhsLowerTangentPoint]);

            for (int i = rhsLowerTangentPoint; i != rhsUpperTangentPoint; )
            {
               convexHull.push_back(rhs[i]);
               i = (i == rhs.size() - 1) ? 0 : i + 1;
            }
            convexHull.push_back(rhs[rhsUpperTangentPoint]);
            
            for (int i = lhsUpperTangentPoint; i != indexLeftMostPointLhs; )
            {
               convexHull.push_back(lhs[i]);
               i = (i == lhs.size() - 1) ? 0 : i + 1;
            }
         }

         /// <summary>
         /// <algorithm name="Merge Convex Hull"></algorithm> 
         /// <runtime complexity="n*logn"></runtime>
         /// </summary>
         /// <param name="vertices"></param>
         /// <param name="convexHull"></param>
         void ConvexHull(const std::vector<Vertex>& vertices, std::vector<Vertex>& convexHull)
         {
            if (vertices.size() < 6)
            {
               GiftWrap::ConvexHull(vertices, convexHull);
            }
            else
            {
               std::size_t const half_size = vertices.size() / 2;
               std::vector<Vertex> split_loConvexHull;
               std::vector<Vertex> split_hiConvexHull;
               std::vector<Vertex> split_lo(vertices.begin(), vertices.begin() + half_size);
               std::vector<Vertex> split_hi(vertices.begin() + half_size, vertices.end());

               ConvexHull(split_lo, split_loConvexHull);
               ConvexHull(split_hi, split_hiConvexHull);

               MergeHull(split_loConvexHull, split_hiConvexHull, convexHull);
            }
         }
      }

      namespace QuickHull
      {
         enum TriHullEdge
         {
            EDGE_AP = 0,
            EDGE_PB,

            LAST_TRI_HULL_EDGE
         };

         enum QuickHullEdge
         {
            EDGE_AB = 0,
            EDGE_BC,
            EDGE_CD,
            EDGE_DA,

            LAST_QUICK_HULL_EDGE
         };

         enum QuickHullVertex
         {
            VERTEX_A = 0,
            VERTEX_B,
            VERTEX_C,
            VERTEX_D,

            LAST_QUICK_HULL_VERTEX
         };

         void TriHull(const std::vector<Vertex>& vertices, const Edge& edge, std::vector<Vertex>& convexHull)
         {
            if (vertices.size() == 1)
            {
               convexHull.push_back(vertices[0]);
            }
            else if(vertices.size() > 1)
            {
               // Find point furthest from edge
               float heuristic = 0.0f;
               int hullIndex = 0;
               for (int i = 0; i < vertices.size(); i++)
               {
                  Edge vectorFromStart(Vertex(vertices[i].x - edge.start.x, vertices[i].y - edge.start.y));
                  Edge orthoVector(Vertex(edge.start.y - edge.end.y, edge.end.x - edge.start.x));
                  float dotProduct = vectorFromStart * orthoVector;
                  float lenSquared = orthoVector * orthoVector;
                  float currHeuristic = dotProduct * dotProduct / lenSquared;
                  if (currHeuristic > heuristic)
                  {
                     heuristic = currHeuristic;
                     hullIndex = i;
                  }
               }

               // Get Point Subsets
               std::vector<Vertex> verticesS[LAST_QUICK_HULL_EDGE];
               for (int i = 0; i < static_cast<int>(vertices.size()); i++)
               {
                  if (Orientation(edge.start, vertices[hullIndex], vertices[i]) == ORIENTATION_CW)
                  {
                     verticesS[EDGE_AP].push_back(vertices[i]);
                  }
                  else if (Orientation(vertices[hullIndex], edge.end, vertices[i]) == ORIENTATION_CW)
                  {
                     verticesS[EDGE_PB].push_back(vertices[i]);
                  }
               }

               TriHull(verticesS[EDGE_AP], Edge(edge.start, vertices[hullIndex]), convexHull);
               convexHull.push_back(vertices[hullIndex]);
               TriHull(verticesS[EDGE_PB], Edge(vertices[hullIndex], edge.end), convexHull);
            }
         }

         /// <summary>
         /// <algorithm name="Quick Convex Hull"></algorithm> 
         /// <runtime complexity="n²"></runtime>
         /// </summary>
         /// <param name="vertices"></param>
         /// <param name="convexHull"></param>
         void ConvexHull(const std::vector<Vertex>& vertices, std::vector<Vertex>& convexHull)
         {
            if (vertices.size() < 6)
            {
               GiftWrap::ConvexHull(vertices, convexHull);
            }
            else
            {
               int hullPointIndices[LAST_QUICK_HULL_VERTEX];
               ::memset(hullPointIndices, 0, static_cast<int>(LAST_QUICK_HULL_VERTEX) * sizeof(int));

               // Get Extreme Points
               for (int i = 1; i < static_cast<int>(vertices.size()); i++)
               {
                  if (vertices[i].y > vertices[hullPointIndices[VERTEX_A]].y)
                  {
                     hullPointIndices[VERTEX_A] = i;
                  }
                  else if (vertices[i].y < vertices[hullPointIndices[VERTEX_C]].y)
                  {
                     hullPointIndices[VERTEX_C] = i;
                  }

                  if (vertices[i].x < vertices[hullPointIndices[VERTEX_B]].x)
                  {
                     hullPointIndices[VERTEX_B] = i;
                  }
                  else if (vertices[i].x > vertices[hullPointIndices[VERTEX_D]].x)
                  {
                     hullPointIndices[VERTEX_D] = i;
                  }
               }

               // Get Point Subsets
               std::vector<Vertex> verticesS[LAST_QUICK_HULL_EDGE];
               for (int i = 0; i < static_cast<int>(vertices.size()); i++)
               {
                  if (Orientation(vertices[hullPointIndices[VERTEX_B]], vertices[hullPointIndices[VERTEX_C]], vertices[i]) == ORIENTATION_CW)
                  {
                     verticesS[EDGE_BC].push_back(vertices[i]);
                  }
                  else if (Orientation(vertices[hullPointIndices[VERTEX_C]], vertices[hullPointIndices[VERTEX_D]], vertices[i]) == ORIENTATION_CW)
                  {
                     verticesS[EDGE_CD].push_back(vertices[i]);
                  }
                  else if (Orientation(vertices[hullPointIndices[VERTEX_D]], vertices[hullPointIndices[VERTEX_A]], vertices[i]) == ORIENTATION_CW)
                  {
                     verticesS[EDGE_DA].push_back(vertices[i]);
                  }
                  else if (Orientation(vertices[hullPointIndices[VERTEX_A]], vertices[hullPointIndices[VERTEX_B]], vertices[i]) == ORIENTATION_CW)
                  {
                     verticesS[EDGE_AB].push_back(vertices[i]);
                  }
               }

               // Add leftmost point first
               convexHull.push_back(vertices[hullPointIndices[VERTEX_B]]);

               if (vertices[hullPointIndices[VERTEX_B]] != vertices[hullPointIndices[VERTEX_C]])
               {
                  TriHull(verticesS[EDGE_BC], Edge(vertices[hullPointIndices[VERTEX_B]], vertices[hullPointIndices[VERTEX_C]]), convexHull);
                  convexHull.push_back(vertices[hullPointIndices[VERTEX_C]]);
               }

               if (vertices[hullPointIndices[VERTEX_C]] != vertices[hullPointIndices[VERTEX_D]])
               {
                  TriHull(verticesS[EDGE_CD], Edge(vertices[hullPointIndices[VERTEX_C]], vertices[hullPointIndices[VERTEX_D]]), convexHull);
                  convexHull.push_back(vertices[hullPointIndices[VERTEX_D]]);
               }

               if (vertices[hullPointIndices[VERTEX_D]] != vertices[hullPointIndices[VERTEX_A]])
               {
                  TriHull(verticesS[EDGE_DA], Edge(vertices[hullPointIndices[VERTEX_D]], vertices[hullPointIndices[VERTEX_A]]), convexHull);
                  if (vertices[hullPointIndices[VERTEX_B]] != vertices[hullPointIndices[VERTEX_A]])
                  {
                     convexHull.push_back(vertices[hullPointIndices[VERTEX_A]]);
                  }
               }

               if (vertices[hullPointIndices[VERTEX_B]] != vertices[hullPointIndices[VERTEX_A]])
               {
                  TriHull(verticesS[EDGE_AB], Edge(vertices[hullPointIndices[VERTEX_A]], vertices[hullPointIndices[VERTEX_B]]), convexHull);
               }
            }
         }
      }
   }
}