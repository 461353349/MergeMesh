#ifndef MESHMERGEZIPPER_H
#define MESHMERGEZIPPER_H

#include <set>
#include <list>
#include <vector>
#include <queue>
#include <map>

#include "PUMBase\GeomBase\Mesh.h"
#include "BaseClass.h"
#include "triangulate.h"
#include "PlugIn_Alg_Holefilling\EasyHoleFilling.h"
#include "PUMBase\Utility\KDTree.h"
#include "PlugIn_Alg_MeshRepair_Prereqs.h"
#include <unordered_map>
//#include "Bounds.h"

class MergeBoundryEdge;

class OverlapPatch;

class CubeBounds{
public:
	PUM::Point3D center;
	PUM::Vector3D length;
	CubeBounds(){}
	CubeBounds(PUM::Point3D& c,PUM::Vector3D& v):center(c),length(v){}
	void setBounds(PUM::Point3D& c,PUM::Vector3D& v){
		center = c;
		length = v;
	}
};

class MergeBoundryTriangle{
public:
	int							triangleIndex;
	vector<MergeBoundryEdge>	boundryEdge;
	vector<PUM::Point3D>		cutPointVec;
};

class Vertexs{
public:
	set<int> vertexs;
	void addVertex(int vertexId){vertexs.insert(vertexId);}
	int getSize(){return vertexs.size();}
};

class EdgeLoop{
public:
	vector<OctreeBase::Edge> edgeVec;
public:
	int addEdge(OctreeBase::Edge& edge){
		edgeVec.push_back(edge);
		return 0;
	}
	int getEdgeNum(){return edgeVec.size();}
	int getDirection(PUM::Mesh *mesh,PUM::Vector3D &dirSum);
	//给定两个端点的编号将edgeloop分成两半
	int partLoop(int e1,int e2,EdgeLoop &el1,EdgeLoop &el2);
	void generateBounds(CubeBounds& bound,PUM::Mesh *mesh);
	void generateVertexSet(set<PUM::Vertex_Index>& ptset);
	void deleteNotInOverlapArea(set<PUM::Vertex_Index>& ptset,OverlapPatch &curr);
	int clear();
};

class MergeBoundryEdge{
public:
	int		vertex1;
	int		vertex2;
	int		adjTrianglePatchId;//指向它所在的triangle
	SimpleTriagle triWall[4];
public:
	MergeBoundryEdge(int v1,int v2,int triId):vertex1(v1),vertex2(v2),adjTrianglePatchId(triId){}
	int generateWallVertex(PUM::Mesh *mesh,OverlapPatch &op,double epsilon,int v,PUM::Point3D &v1,PUM::Point3D &v2);
	int generateTriWall(PUM::Mesh *mesh,OverlapPatch &op,double epsilon);
	int getVertexAdjBoundryEdge(int v1,PUM::Mesh *mesh);//返回该边界边相邻的三角形id
};

class HierarchicalNode{
public:
	PUM::Vertex_Index vertId;
	int depth;
	HierarchicalNode():vertId(-1),depth(0){}
	HierarchicalNode(PUM::Vertex_Index v,int dep):vertId(v),depth(dep){}
};

class PLUGIN_ALG_MESHREPAIR_CLASS OverlapPatch{
public:
	OverlapPatch();
	int modifyPatchVertex(int vertex1,int vertex2);
	int covertSet2Vector();
	int deleteListEle(int vertexId);
	int generatetBoundryTriangle(PUM::Mesh *mesh);
	int generateVertexVec(PUM::Mesh *mesh);
	int generateHierarchicalStructure(PUM::Mesh *mesh,EdgeLoop &el);
	void gengerateBounds(CubeBounds& bound,PUM::Mesh *mesh);
	int testAdjVec(vector<PUM::Vertex_Index> &vertexVec,map<int,int> &flag);
	int clear();
public:
	set<int> triangleSet;
	vector<int> triangleVec;
	vector<MergeBoundryTriangle> boundryTriangleVec;
	bool isConvert;
	vector<int> idVec;
	set<PUM::Vertex_Index> vertexSet;
	vector<HierarchicalNode> hiNodeVec;
	vector<PUM::Face_Index> partBoundryFaceVec;
};

class PLUGIN_ALG_MESHREPAIR_CLASS Zipper{

public:
	Zipper();
	Zipper(PUM::Mesh* mesh1,PUM::Mesh* mesh2);
	int				setCepsilon(double cp);
	int				setDestMesh(PUM::Mesh* mesh);
	int				setSourceMesh(PUM::Mesh* mesh);
	int				getDestPatchSize();
	int				getSourcePatchSize();
	EdgeLoop		&getSourceEdgeLoop();
	EdgeLoop		&getDestEdgeLoop();
	PUM::Mesh		*getSourceMesh();
	PUM::Mesh		*getDestMesh();
	OverlapPatch	&getSourcePatch();
	OverlapPatch	&getDestPatch();

	int drawPatch(OverlapPatch &op,PUM::Mesh *mesh,int color);
	int drawEdgeLoop(EdgeLoop &edgeLoop,PUM::Mesh* mesh,PUM::Color& color);
	int drawHoleBoundry();
	int drawPath();
	int draw();

	double findNearestEleKdTree(PUM::Point3D &queryPoint,PUM::KDTree *queryKdTree,int &resultId);

	int findHole();
	void findEdgeLoopHolePart(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &srcLoop,EdgeLoop &destLoop,EdgeLoop &resultLoop1,EdgeLoop &resultLoop2);
	int findNearLoop(PUM::Mesh *mesh,EdgeLoop &loop1,EdgeLoop &loop2,PUM::KDTree *loopKdtree);
	int getBoundryEdgeNextPoint(PUM::Mesh *mesh,PUM::Vertex_Index vi,vector<PUM::Vertex_Index>& adjVertexs,unordered_map<int,int> &flag);
	int getBoundryEdgeNextPoint(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,vector<PUM::Vertex_Index>& adjVertexs,vector<PUM::Vertex_Index>& nextPoints);
	//生成一个完整的edgeLoop
	int generateOneEdgeLoop(PUM::Mesh *mesh,PUM::Vertex_Index vi,unordered_map<int,int> &flag,EdgeLoop &loop);
	//生成mesh的边界edgeloop
	int generateEdgeLoop(PUM::Mesh *mesh,vector<EdgeLoop> &edgeLoopVec);
	bool dfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteFaces);
	bool bfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,unordered_map<PUM::Vertex_Index,int>& pres,set<PUM::Vertex_Index>& deleteFaces);
	void deleteOneLoopCircle(PUM::Mesh *mesh,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteFaces);
	bool deleteLoopCircle(PUM::Mesh *mesh);
	
	//
	bool deleteNonMainfoldVerticle(PUM::Mesh *mesh);
	//找到边数目最多的edgeloop
	int findMaxLoop(vector<EdgeLoop> &loopVec);
	//生成overlapPatch的edgeloop
	int generateOverLapPatchEdgeLoop(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &edgeLoop);
	//将重合部分边界分成三部分
	int generateOverLapPatchLoops(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop1,EdgeLoop &loop2);
	//生成距离与另一个edgeloopPart相近的部分(实际三维距离)
	int findNearOtherEdgeLoopPart(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &loop1,EdgeLoop &loop2,EdgeLoop &loopPart1);
	//同上，用的是经纬度数据
	int findNearOtherEdgeLoopPartUsingLola(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &loop1,EdgeLoop &loop2,EdgeLoop &loopPart1);
	int closeHole(PUM::Mesh *mesh);
	//监测hole绕向和周围面片法向是否一致
	bool isHoleDirectionRight(PUM::Mesh *mesh,int index0,int index1);
	//生成边界区域的边界点（用户指定的两个端点）
	int generateOverLapPatchBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop);
	int generateOVerLapPatchBoundry();
	int generateOverLapSingleAreaBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop);
	//生成非相同深度的边界
	int generateNotSameDepthBoundry();
	int generateNotSameDepthBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop1,EdgeLoop &loop2);
	//生成重叠区域的顶点
	int generateOverlapPatchVertex();
	//找到重叠区域（三角面片）
	int findOverlappArea();
	int findOverlappArea(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	void findOverlapAreaUsing2d(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	void findOverlapAreaLolat(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	//找到mesh相交区域的两个端点
	bool findMesh2dIntersectionPoint(PUM::Mesh* mesh1,PUM::Mesh* mesh2);
	int findPath(PUM::Mesh *mesh1,PUM::Mesh *mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	//两块mesh简单的侵蚀（丢掉重叠的部分）
	int eatBack(PUM::Mesh* currMesh,PUM::Mesh* oppositeMesh,OverlapPatch &currPatch,OverlapPatch &opposite);
	int eatBack();
	//mesh由边界向最短路扩散
	void sourceMeshEatBack(PUM::Mesh* currMesh,OverlapPatch &currPatch,EdgeLoop &loop);
	//vertexset是边界点，边界点和path组成一个封闭的区域
	void sourceMeshEatBack(PUM::Mesh* currMesh,set<PUM::Vertex_Index> &vertexset);
	//利用mesh的经纬度信息进行侵蚀
	void eatBackUsing2d(PUM::Mesh* currMesh,PUM::Mesh* oppositeMesh,OverlapPatch &currPatch,OverlapPatch &oppositePatch);
	//边界被startIndex,endIndex分成两部分，判断那一部分是在重叠区域
	int edgeLoopInoverlapArea(PUM::Mesh* mesh,set<PUM::Vertex_Index> &vertexs1,set<PUM::Vertex_Index> &vertexs2,OverlapPatch &patch);
	
	//获得sourcepatch和destPatch的边界边以及生成他们的wall
	int cippingInit(PUM::Mesh *mesh,OverlapPatch &op);
	int clipping();
	//将destMesh合并到sourcemesh上，同时修改destmesh上的edgeLoop(边顶点序号的修改)
	int mergeMesh(PUM::Mesh *sourceMesh,PUM::Mesh *destMesh,EdgeLoop &el);
	int mergeMesh();
	//产生kdtree
	int generateKDtree(set<PUM::Vertex_Index> &vertexSet,PUM::Mesh *mesh,PUM::KDTree *newKdtree);
	int generateLolaKDtree(set<PUM::Vertex_Index> &vertexSet,PUM::Mesh *mesh,PUM::KDTree *newKdtree);
	//寻找合适的扩展范围
	int findSuitExpandDepth(PUM::Mesh *scMesh,OverlapPatch &scPatch,EdgeLoop &el,PUM::Mesh *oppMesh,OverlapPatch &oppPatch);
	//重合区域的扩散
	int overLapPatchExtend(PUM::Mesh *scMesh,OverlapPatch &scPatch,int depth);
	//去除mesh中的孤立点
	int deleteIsolatePoints(PUM::Mesh *scMesh);
	//删除面
	void deleteFaces(PUM::Mesh *scMesh,vector<PUM::Face_Index>& faces);
	void deleteFaces(PUM::Mesh *scMesh,set<PUM::Face_Index>& faces);
	void deleteVertexs(PUM::Mesh *scMesh,set<PUM::Face_Index>& vertexs);
	//
	int clear();
	//产生扩散边界点
	int generateZipperPoint(PUM::Mesh *srcMesh,PUM::Mesh *destMesh,OverlapPatch &srcPatch,OverlapPatch& destPatch);
	void findMaxLoop(PUM::Mesh* mesh,EdgeLoop& maxLoop);
	void covert3dPoint2LoLadata(PUM::Point3D &threedPoint,PUM::Point3D &lolaPoint);
	void meshMergeManager();
	PUM::Mesh* meshMerge(PUM::Mesh* meshs[],int left,int right,const char *path);
	PUM::Mesh* meshMerge(PUM::Mesh* meshs[], int num,const char *path);
	void setMatrixPath(string &ma);
	void setMatrixPath(const char *path);
private:
	int					smallPartNum;
	double				lltitudeDiffer;		//经纬度误差
	double				cepsilon;			//三维空间距离误差

	OverlapPatch		destPatch;			//目的重叠区域
	OverlapPatch		sourcePatch;		//源重叠区域
	PUM::Mesh			*destMesh;			//目的mesh
	PUM::Mesh			*sourceMesh;		//源mesh
	EdgeLoop			sourceEdgeLoop;		//源mesh边界
	EdgeLoop			destEdgeLoop;		//目的mesh边界
	vector<PUM::Vertex_Index> holeVertexVec; //洞顶点
	
	PUM::Vertex_Index	startIndex;			//sourceMesh最短路的起始顶点
	PUM::Vertex_Index	endIndex;			//sourceMesh最短路的结束顶点
	PUM::Point3D		startPoint;
	PUM::Point3D		endPoint;
	vector<PUM::Vertex_Index> shortestPath;
	string				matrixPath;
};

#endif