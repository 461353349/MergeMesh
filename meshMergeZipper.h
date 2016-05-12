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
	//���������˵�ı�Ž�edgeloop�ֳ�����
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
	int		adjTrianglePatchId;//ָ�������ڵ�triangle
	SimpleTriagle triWall[4];
public:
	MergeBoundryEdge(int v1,int v2,int triId):vertex1(v1),vertex2(v2),adjTrianglePatchId(triId){}
	int generateWallVertex(PUM::Mesh *mesh,OverlapPatch &op,double epsilon,int v,PUM::Point3D &v1,PUM::Point3D &v2);
	int generateTriWall(PUM::Mesh *mesh,OverlapPatch &op,double epsilon);
	int getVertexAdjBoundryEdge(int v1,PUM::Mesh *mesh);//���ظñ߽�����ڵ�������id
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
	//����һ��������edgeLoop
	int generateOneEdgeLoop(PUM::Mesh *mesh,PUM::Vertex_Index vi,unordered_map<int,int> &flag,EdgeLoop &loop);
	//����mesh�ı߽�edgeloop
	int generateEdgeLoop(PUM::Mesh *mesh,vector<EdgeLoop> &edgeLoopVec);
	bool dfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteFaces);
	bool bfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,unordered_map<PUM::Vertex_Index,int>& pres,set<PUM::Vertex_Index>& deleteFaces);
	void deleteOneLoopCircle(PUM::Mesh *mesh,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteFaces);
	bool deleteLoopCircle(PUM::Mesh *mesh);
	
	//
	bool deleteNonMainfoldVerticle(PUM::Mesh *mesh);
	//�ҵ�����Ŀ����edgeloop
	int findMaxLoop(vector<EdgeLoop> &loopVec);
	//����overlapPatch��edgeloop
	int generateOverLapPatchEdgeLoop(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &edgeLoop);
	//���غϲ��ֱ߽�ֳ�������
	int generateOverLapPatchLoops(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop1,EdgeLoop &loop2);
	//���ɾ�������һ��edgeloopPart����Ĳ���(ʵ����ά����)
	int findNearOtherEdgeLoopPart(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &loop1,EdgeLoop &loop2,EdgeLoop &loopPart1);
	//ͬ�ϣ��õ��Ǿ�γ������
	int findNearOtherEdgeLoopPartUsingLola(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &loop1,EdgeLoop &loop2,EdgeLoop &loopPart1);
	int closeHole(PUM::Mesh *mesh);
	//���hole�������Χ��Ƭ�����Ƿ�һ��
	bool isHoleDirectionRight(PUM::Mesh *mesh,int index0,int index1);
	//���ɱ߽�����ı߽�㣨�û�ָ���������˵㣩
	int generateOverLapPatchBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop);
	int generateOVerLapPatchBoundry();
	int generateOverLapSingleAreaBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop);
	//���ɷ���ͬ��ȵı߽�
	int generateNotSameDepthBoundry();
	int generateNotSameDepthBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop1,EdgeLoop &loop2);
	//�����ص�����Ķ���
	int generateOverlapPatchVertex();
	//�ҵ��ص�����������Ƭ��
	int findOverlappArea();
	int findOverlappArea(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	void findOverlapAreaUsing2d(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	void findOverlapAreaLolat(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	//�ҵ�mesh�ཻ����������˵�
	bool findMesh2dIntersectionPoint(PUM::Mesh* mesh1,PUM::Mesh* mesh2);
	int findPath(PUM::Mesh *mesh1,PUM::Mesh *mesh2,OverlapPatch &patch1,OverlapPatch &patch2);
	//����mesh�򵥵���ʴ�������ص��Ĳ��֣�
	int eatBack(PUM::Mesh* currMesh,PUM::Mesh* oppositeMesh,OverlapPatch &currPatch,OverlapPatch &opposite);
	int eatBack();
	//mesh�ɱ߽������·��ɢ
	void sourceMeshEatBack(PUM::Mesh* currMesh,OverlapPatch &currPatch,EdgeLoop &loop);
	//vertexset�Ǳ߽�㣬�߽���path���һ����յ�����
	void sourceMeshEatBack(PUM::Mesh* currMesh,set<PUM::Vertex_Index> &vertexset);
	//����mesh�ľ�γ����Ϣ������ʴ
	void eatBackUsing2d(PUM::Mesh* currMesh,PUM::Mesh* oppositeMesh,OverlapPatch &currPatch,OverlapPatch &oppositePatch);
	//�߽类startIndex,endIndex�ֳ������֣��ж���һ���������ص�����
	int edgeLoopInoverlapArea(PUM::Mesh* mesh,set<PUM::Vertex_Index> &vertexs1,set<PUM::Vertex_Index> &vertexs2,OverlapPatch &patch);
	
	//���sourcepatch��destPatch�ı߽���Լ��������ǵ�wall
	int cippingInit(PUM::Mesh *mesh,OverlapPatch &op);
	int clipping();
	//��destMesh�ϲ���sourcemesh�ϣ�ͬʱ�޸�destmesh�ϵ�edgeLoop(�߶�����ŵ��޸�)
	int mergeMesh(PUM::Mesh *sourceMesh,PUM::Mesh *destMesh,EdgeLoop &el);
	int mergeMesh();
	//����kdtree
	int generateKDtree(set<PUM::Vertex_Index> &vertexSet,PUM::Mesh *mesh,PUM::KDTree *newKdtree);
	int generateLolaKDtree(set<PUM::Vertex_Index> &vertexSet,PUM::Mesh *mesh,PUM::KDTree *newKdtree);
	//Ѱ�Һ��ʵ���չ��Χ
	int findSuitExpandDepth(PUM::Mesh *scMesh,OverlapPatch &scPatch,EdgeLoop &el,PUM::Mesh *oppMesh,OverlapPatch &oppPatch);
	//�غ��������ɢ
	int overLapPatchExtend(PUM::Mesh *scMesh,OverlapPatch &scPatch,int depth);
	//ȥ��mesh�еĹ�����
	int deleteIsolatePoints(PUM::Mesh *scMesh);
	//ɾ����
	void deleteFaces(PUM::Mesh *scMesh,vector<PUM::Face_Index>& faces);
	void deleteFaces(PUM::Mesh *scMesh,set<PUM::Face_Index>& faces);
	void deleteVertexs(PUM::Mesh *scMesh,set<PUM::Face_Index>& vertexs);
	//
	int clear();
	//������ɢ�߽��
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
	double				lltitudeDiffer;		//��γ�����
	double				cepsilon;			//��ά�ռ�������

	OverlapPatch		destPatch;			//Ŀ���ص�����
	OverlapPatch		sourcePatch;		//Դ�ص�����
	PUM::Mesh			*destMesh;			//Ŀ��mesh
	PUM::Mesh			*sourceMesh;		//Դmesh
	EdgeLoop			sourceEdgeLoop;		//Դmesh�߽�
	EdgeLoop			destEdgeLoop;		//Ŀ��mesh�߽�
	vector<PUM::Vertex_Index> holeVertexVec; //������
	
	PUM::Vertex_Index	startIndex;			//sourceMesh���·����ʼ����
	PUM::Vertex_Index	endIndex;			//sourceMesh���·�Ľ�������
	PUM::Point3D		startPoint;
	PUM::Point3D		endPoint;
	vector<PUM::Vertex_Index> shortestPath;
	string				matrixPath;
};

#endif