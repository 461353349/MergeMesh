# MergeMesh
将两张有重叠区域的mesh捏合成一张mesh，三维建筑模型
算法流程：
(1)找重叠区域
(2)找两个mesh的边界，并将边界投影到经纬度坐标，最终得到了两个二维的多边形
(3)求两个多边形的交点
(4)求两个交点在一个mesh上的最短路，最短路的距离是当前mesh上的点与另一个mesh上点距离的最小值（可以用kdtree快速求得）
(5)两个mesh分别向之前最短路路径扩散，这样会形成细缝
(6)缝补细缝