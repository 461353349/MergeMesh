#include "stdafx.h"

#include "PointPolygonLocation.h"
#include "meshMergeZipper.h"
#include "GLDraw\OpenGLDraw.h"
#include "PUMBase\GeomBase\Mesh.h"
#include "extraOperation.h"
#include "CoordTransform.h"
#include "SegmentIntersection.h"


#include <queue>
#include <stack>
#include <unordered_set>

Zipper::Zipper():lltitudeDiffer(2e-5),smallPartNum(50),startIndex(0),endIndex(0){};

Zipper::Zipper(PUM::Mesh* mesh1,PUM::Mesh* mesh2):sourceMesh(mesh1),destMesh(mesh2),lltitudeDiffer(2e-5),smallPartNum(50),startIndex(0),endIndex(0){}

int Zipper::setCepsilon(double cp){
	cepsilon = cp;
	return 0;
}

void Zipper:: setMatrixPath(string &ma){
	matrixPath = ma;
}

void Zipper::setMatrixPath(const char *path){
	matrixPath.assign(path);
}

int Zipper::setDestMesh(PUM::Mesh* mesh){
	destMesh = mesh;
	return 0;
}
int Zipper::setSourceMesh(PUM::Mesh* mesh){
	sourceMesh = mesh;
	return 0;
}
PUM::Mesh* Zipper::getSourceMesh(){
	return sourceMesh;
}
PUM::Mesh* Zipper::getDestMesh(){
	return destMesh;
}

OverlapPatch& Zipper::getSourcePatch(){
	return sourcePatch;
}
OverlapPatch& Zipper::getDestPatch(){
	return destPatch;
}

EdgeLoop& Zipper::getSourceEdgeLoop(){
	return sourceEdgeLoop;
}
EdgeLoop& Zipper::getDestEdgeLoop(){
	return destEdgeLoop;
}

int convertLinkedlist2Vector(PUM::LinkedList<PUM::Face_Index> &ll,vector<PUM::Face_Index> &dest ){
	while( !ll.EndOfList() )                  // 导入一邻域种子面
	{
		dest.push_back( ll.Data() );
		ll.Next();
	}
	return 0;
}

int Zipper::findOverlappArea(){
	//findOverlappArea(sourceMesh,destMesh,sourcePatch,destPatch);
	//findOverlappArea(destMesh,sourceMesh,destPatch,sourcePatch);
	findOverlapAreaLolat(sourceMesh,destMesh,sourcePatch,destPatch);
	findOverlapAreaLolat(destMesh,sourceMesh,destPatch,sourcePatch);
	//findOverlapAreaUsing2d(sourceMesh,destMesh,sourcePatch,destPatch);
	//findOverlapAreaUsing2d(destMesh,sourceMesh,destPatch,sourcePatch);
	OUTPUT_TRACE_T_2("size:%d",sourcePatch.triangleSet.size());
	OUTPUT_TRACE_T_2("size:%d",destPatch.triangleSet.size());
	return 0;
}

void Zipper::findOverlapAreaUsing2d(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2){
	vector<EdgeLoop> loopVec;
	generateEdgeLoop(mesh2,loopVec);
	int loopIndex = findMaxLoop(loopVec);
	CCoordTransform::Init(matrixPath);
	if(loopIndex != -1){
		//寻找相对mesh的最外层边界（由点集构成的多边形）
		EdgeLoop& oppLoop = loopVec[loopIndex];
		const int currSize = oppLoop.edgeVec.size()+1;
		vector<Point> points(currSize+1);
		set<PUM::Vertex_Index> pointset;
		oppLoop.generateVertexSet(pointset);
		for(auto it = pointset.begin(); it != pointset.end(); it++){
			PUM::Point3D &pt = mesh2->Point(*it);
			//转换成经纬度
			auto lonlat = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
			points.push_back(Point(lonlat.first,lonlat.second));
			//WriteString2File_Debug("D:\\twoloops.obj","v " +toString(lonlat.first)+ " "+ toString(lonlat.second) + " 0.0");
		}
		points.push_back(points[0]);
		//寻找mesh上在上边所求多边形内部的面片(至少有一个顶点在多边形内部)，并加入到patch上
		for(int i = 0;i<mesh1->GetFaceNum();i++){
			PUM::Face& currFace = mesh1->GetFace(i);
			int vertexNum = currFace.GetVertexNum();
			for(int j = 0;j<vertexNum;j++){
				int pointIndex = currFace.GetVertex(j);
				PUM::Point3D& currPoint = mesh1->Point(pointIndex);
				auto lonlat = CCoordTransform::TransformLonlat(currPoint.x(),currPoint.y(),currPoint.z());
				Point curr(lonlat.first,lonlat.second);
				//WriteString2File_Debug("D:\\jingweidu_mesh.obj","v " +toString(lonlat.first)+ " "+ toString(lonlat.second) + " 0.0");
				//if(point_in_polygon_check_edge(curr,points)){
				if(cn_PnPoly(curr,points) == 1){
				//if(){
					patch1.triangleSet.insert(i);
					break;
				}
			}
		}
	}
}

void Zipper::findOverlapAreaLolat(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2){
	CCoordTransform::Init(matrixPath);
	PUM::KDTree *myKdTree = new PUM::KDTree(3);
	//const PUM::Array<PUM::Point3D> &destArray = mesh2->GetPtList();
	int pointNum = mesh2->GetVertexNum();
	for(int i = 0;i< pointNum;i++){
		PUM::Point3D &pt = mesh2->Point(i);
		auto lonlat = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
		myKdTree->insert3(lonlat.first,lonlat.second,0,i);
		//WriteString2File_Debug("D:\\mesh_lola.obj","v " +toString(lonlat.first)+ " "+ toString(lonlat.second) +" "+ toString(0));
	}
	const PUM::Array<PUM::Point3D> &sourceArray = mesh1->GetPtList();
	const double angelCos = 0.9;
	for(int i = 0;i<mesh1->GetVertexNum();i++){
		PUM::Point3D &currPoint = mesh1->Point(i);
		auto lonlat1 = CCoordTransform::TransformLonlat(currPoint.x(),currPoint.y(),currPoint.z());
		vector<PUM::Face_Index> vertexAdjFaceVec;
		//sourceMesh->GetAdjFacessOfVertex_kneighbor(i,vertexAdjFaceVec);
		PUM::LinkedList<PUM::Face_Index> adjFaceList;
		mesh1->GetFacesOfVertex(i, adjFaceList );
		convertLinkedlist2Vector(adjFaceList,vertexAdjFaceVec);
		PUM::kdres* kd_result = myKdTree->nearest3(lonlat1.first, lonlat1.second,0);
		//if(!kd_result)
		//	continue;
		double pos[3];
		int temp_id = myKdTree->res_item(kd_result, pos);
		myKdTree->kd_res_free(kd_result);//检测内存泄漏
		PUM::Vector3D& normal1 = mesh1->GetPointNormal(i);
		PUM::Vector3D& normal2 = mesh2->GetPointNormal(temp_id);

		PUM::Point2D sourcelola(lonlat1.first,lonlat1.second);
		PUM::Point2D resultlola(pos[0],pos[1]);
		double dis = Distance(sourcelola,resultlola);
		if( dis < lltitudeDiffer ){
			vector<PUM::Face_Index> destVertexAdjFaceVec;
			//destMesh->GetAdjFacessOfVertex_kneighbor(temp_id,destVertexAdjFaceVec);
			PUM::LinkedList<PUM::Face_Index> destAdjFaceList;
			mesh2->GetFacesOfVertex(temp_id, destAdjFaceList );
			convertLinkedlist2Vector(destAdjFaceList,destVertexAdjFaceVec);
			patch1.triangleSet.insert(vertexAdjFaceVec.begin(),vertexAdjFaceVec.end());
			patch2.triangleSet.insert(destVertexAdjFaceVec.begin(),destVertexAdjFaceVec.end());
		}

	}
	delete myKdTree;
}

int Zipper::findOverlappArea(PUM::Mesh* mesh1,PUM::Mesh* mesh2,OverlapPatch &patch1,OverlapPatch &patch2){
	PUM::KDTree *myKdTree = new PUM::KDTree(3);
	//const PUM::Array<PUM::Point3D> &destArray = mesh2->GetPtList();
	int pointNum = mesh2->GetVertexNum();
	for(int i = 0;i< pointNum;i++){
		PUM::Point3D &pt = mesh2->Point(i);
		myKdTree->insert3(pt.x(),pt.y(),pt.z(),i);
	}
	const PUM::Array<PUM::Point3D> &sourceArray = mesh1->GetPtList();
	CCoordTransform::Init(matrixPath);
	const double lltitudeDiffer = 3e-5 ;
	const double angelCos = 0.9;
	for(int i = 0;i<mesh1->GetVertexNum();i++){
		vector<PUM::Face_Index> vertexAdjFaceVec;
		//sourceMesh->GetAdjFacessOfVertex_kneighbor(i,vertexAdjFaceVec);
		PUM::LinkedList<PUM::Face_Index> adjFaceList;
		mesh1->GetFacesOfVertex(i, adjFaceList );
		convertLinkedlist2Vector(adjFaceList,vertexAdjFaceVec);
		const PUM::Point3D& currPoint = sourceArray[i];
		PUM::kdres* kd_result = myKdTree->nearest3(currPoint.x(), currPoint.y(),currPoint.z());
		//if(!kd_result)
		//	continue;
		double pos[3];
		int temp_id = myKdTree->res_item(kd_result, pos);
		myKdTree->kd_res_free(kd_result);//检测内存泄漏
		PUM::Point3D resultPoint(pos[0],pos[1],pos[2]);
		PUM::Vector3D& normal1 = mesh1->GetPointNormal(i);
		PUM::Vector3D& normal2 = mesh2->GetPointNormal(temp_id);
		PUM::Point3D &pt1 = mesh1->Point(i);
		auto lonlat1 = CCoordTransform::TransformLonlat(pt1.x(),pt1.y(),pt1.z());
		auto lonlat2 = CCoordTransform::TransformLonlat(pos[0],pos[1],pos[2]);
		double d = Distance(currPoint,resultPoint);	
		double lld = Distance(Point2D(lonlat1.first,lonlat1.second),Point2D(lonlat2.first,lonlat2.second));
		double dis = normal1 * normal2;
		if( d<cepsilon || lld < lltitudeDiffer /*|| dis > angelCos*/){
			vector<PUM::Face_Index> destVertexAdjFaceVec;
			//destMesh->GetAdjFacessOfVertex_kneighbor(temp_id,destVertexAdjFaceVec);
			PUM::LinkedList<PUM::Face_Index> destAdjFaceList;
			mesh2->GetFacesOfVertex(temp_id, destAdjFaceList );
			convertLinkedlist2Vector(destAdjFaceList,destVertexAdjFaceVec);
			patch1.triangleSet.insert(vertexAdjFaceVec.begin(),vertexAdjFaceVec.end());
			patch2.triangleSet.insert(destVertexAdjFaceVec.begin(),destVertexAdjFaceVec.end());
		}

	}
	delete myKdTree;
	return 0;
}

int Zipper::generateLolaKDtree(set<PUM::Vertex_Index> &vertexSet,PUM::Mesh *mesh,PUM::KDTree *kdtree){
	CCoordTransform::Init(matrixPath);
	int pointNum = mesh->GetVertexNum();
	for(int i = 0;i< pointNum;i++){
		PUM::Point3D &pt = mesh->Point(i);
		auto lonlat = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
		kdtree->insert3(lonlat.first,lonlat.second,0,i);
	}
	return 0;
}

bool Zipper::findMesh2dIntersectionPoint(PUM::Mesh* mesh1,PUM::Mesh* mesh2){
	//WriteString2File_Debug("test.txt","fuck");
	CCoordTransform::Init(matrixPath);
	EdgeLoop loop1,loop2;
	findMaxLoop(mesh1,loop1);
	sourceEdgeLoop = loop1;
	findMaxLoop(mesh2,loop2);
	destEdgeLoop = loop2;
	/*set<PUM::Vertex_Index> vertexs1,vertexs2;
	loop1.generateVertexSet(vertexs1);
	loop2.generateVertexSet(vertexs2);
	PUM::KDTree *kdtree1 = new PUM::KDTree(3);
	PUM::KDTree *kdtree2 = new PUM::KDTree(3);
	generateLolaKDtree(vertexs1,mesh1,kdtree1);
	generateLolaKDtree(vertexs2,mesh2,kdtree2);*/
	vector<int> edgeIndexs;
	for(int i =0;i<loop1.edgeVec.size();i++){
		for(int j = 0;j<loop2.edgeVec.size();j++){
			Edge edge1 = loop1.edgeVec[i];
			Edge edge2 = loop2.edgeVec[j];
			PUM::Point3D lola1,lola2;
			covert3dPoint2LoLadata(mesh1->Point(edge1.index[0]),lola1);
			covert3dPoint2LoLadata(mesh1->Point(edge1.index[1]),lola2);
			Segment s1(PUM::Vector2f(lola1.x(),lola1.y()),PUM::Vector2f(lola2.x(),lola2.y()));
			covert3dPoint2LoLadata(mesh2->Point(edge2.index[0]),lola1);
			covert3dPoint2LoLadata(mesh2->Point(edge2.index[1]),lola2);
			Segment s2(PUM::Vector2f(lola1.x(),lola1.y()),PUM::Vector2f(lola2.x(),lola2.y()));
			//PUM::Vector2f inter1,inter2;
			//if(intersect2D_2Segments(s1,s2,&inter1,&inter2)!= 0){
			if(FasterLineSegmentIntersection(s1.P0,s1.P1,s2.P0,s2.P1)){
				edgeIndexs.push_back(i);
				//sourceEdgeLoop.addEdge(edge1);
			}
			//WriteString2File_Debug("D:\\twoloops_02.obj","v " +toString(lola1.x())+ " "+ toString(lola1.y()) +" "+ toString(lola1.z()));
			//WriteString2File_Debug("D:\\twoloops_02.obj","v " +toString(lola2.x())+ " "+ toString(lola2.y()) +" "+ toString(lola2.z()));
		}
	}
	//求距离最远的点对
	if(edgeIndexs.size() < 2) return false;
	int startEdgeIndex = edgeIndexs[0];
	int v1,v2;
	v1 = loop1.edgeVec[startEdgeIndex].index[0];
	v2 = loop1.edgeVec[startEdgeIndex].index[1];
	double dis = Distance(mesh1->Point(v1),mesh1->Point(v2));
	//OUTPUT_TRACE_T_2("dis:：%lf",dis);
	for(int i = 0;i<edgeIndexs.size();i++){
		int edgeIndex1 = edgeIndexs[i];
		int vertexIndex1 = loop1.edgeVec[edgeIndex1].index[0];
		for(int j = i+1;j<edgeIndexs.size();j++){
			int edgeIndex2 = edgeIndexs[j];			
			int vertexIndex2 = loop1.edgeVec[edgeIndex2].index[0];
			double currDis = Distance(mesh1->Point(vertexIndex1),mesh1->Point(vertexIndex2));
			//OUTPUT_TRACE_T_2("dis:：%lf",currDis);
			if(currDis>dis){
				dis = currDis;
				v1 = vertexIndex1;
				v2 = vertexIndex2;
			}
		}
	}
	startIndex = v1;
	endIndex = v2;
	startPoint = mesh1->Point(startIndex);
	endPoint = mesh1->Point(endIndex);
	OUTPUT_TRACE_T_2("start vertex：%d",startIndex);
	OUTPUT_TRACE_T_2("end vertex：%d",endIndex);
	OUTPUT_TRACE_T_2("inter edge size：%d",edgeIndexs.size());
	/*std::priority_queue<pair<double,PUM::Vertex_Index>> heaps;
	for(auto it = vertexs2.begin(); it != vertexs2.end(); it++){
		PUM::Point3D &pt = mesh2->Point(*it);
		auto lonlat = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
		PUM::Point3D queryPoint;
		int resultId;
		double d = findNearestEleKdTree(queryPoint,kdtree1,resultId);
		heaps.push(make_pair(d,*it));
	}*/
	return true;
}


int Zipper::findPath(PUM::Mesh *mesh1,PUM::Mesh *mesh2,OverlapPatch &patch1,OverlapPatch &patch2){
	unordered_set<PUM::Vertex_Index> visit;
	unordered_map<PUM::Vertex_Index,PUM::Vertex_Index> pres;
	unordered_map<PUM::Vertex_Index,double> ds;
	patch1.generateVertexVec(mesh1);
	patch2.generateVertexVec(mesh2);
	PUM::KDTree *kdtree = new PUM::KDTree(3);
	generateKDtree(patch2.vertexSet,mesh2,kdtree);
	priority_queue<ShortestNode> heaps;
	set<PUM::Vertex_Index> &overlapVertexSets = patch1.vertexSet;
	for(int i = 0;i<mesh1->GetVertexNum();i++){
		ds[i] = INT_MAX;
	}
	int resultId;
	double d = findNearestEleKdTree(mesh1->Point(startIndex),kdtree,resultId);
	ds[startIndex] = d;
	heaps.push(ShortestNode(startIndex,d));
	pres[startIndex] = -1;
	while(!heaps.empty()){
		ShortestNode top = heaps.top();
		heaps.pop();
		if(visit.count(top.vertexIndex)) continue;
		if(top.vertexIndex == endIndex) break;
		visit.insert(top.vertexIndex);
		vector<PUM::Vertex_Index> adjvertexs;
		mesh1->Get_Adj_Vertex(top.vertexIndex,adjvertexs);
		
		for(auto it = adjvertexs.begin(); it != adjvertexs.end(); it++){
			if(/*overlapVertexSets.count(*it) &&*/ !visit.count(*it)){
				d = findNearestEleKdTree(mesh1->Point(*it),kdtree,resultId);
				if(ds[*it] > top.dis + d){
					ds[*it] = top.dis + d;
					heaps.push(ShortestNode(*it,ds[*it]));
					pres[*it] = top.vertexIndex;
				}
			}
		}
	}
	OUTPUT_TRACE_T_2("shortest path visit size：%d",pres.size());
	PUM::Vertex_Index index = endIndex;
	while(index != -1){
		shortestPath.push_back(index);
		index = pres[index];
		//WriteString2File_Debug("D:\\path.txt",toString(index));
	}
	delete kdtree;
	return 0;
}

int Zipper::generateOVerLapPatchBoundry(){
	//debug
	//sourcePatch.idVec.push_back(2559);
	//sourcePatch.idVec.push_back(25338);
	//
	generateZipperPoint(sourceMesh,destMesh,sourcePatch,destPatch);
	
	generateOverLapPatchBoundry(sourceMesh,sourcePatch,sourceEdgeLoop);
	generateOverLapPatchBoundry(destMesh,destPatch,destEdgeLoop);
	//generateOverLapSingleAreaBoundry(sourceMesh,sourcePatch,sourceEdgeLoop);
	//generateOverLapSingleAreaBoundry(destMesh,destPatch,destEdgeLoop);
	return 0;
}

void findMaxDis(CubeBounds& bound1,CubeBounds& bound2,int& minDisIndex,int& maxDisIndex,int& ignoreIndex){
	ignoreIndex = 0;
	double minLength = bound1.length[0];
	double dis[3];
	for(int i = 0;i<3;i++){
		if(bound1.length[i] < minLength){
			minLength = bound1.length[i];
			ignoreIndex = i;
		}
		double max1 = bound1.center[i]+bound1.length[i]/2;
		double min1 = bound1.center[i]-bound1.length[i]/2;
		double max2 = bound2.center[i]+bound2.length[i]/2;
		double min2 = bound2.center[i]-bound2.length[i]/2;
		dis[i] = fabs(max1-max2) + fabs(min1-min2);
	}
	switch(ignoreIndex){
		case 0: if(dis[1]<dis[2]) {
					minDisIndex = 1;
					maxDisIndex = 2;
				}else {
					minDisIndex = 2;
					maxDisIndex = 1;
				}break;
		case 1: if(dis[0]<dis[2]){
					minDisIndex = 0;
					maxDisIndex = 2;
				}else{
					minDisIndex = 0;
					maxDisIndex = 2;
				}break;
		case 2: if(dis[0]<dis[1]){
					minDisIndex = 0;
					maxDisIndex = 1;
				}else{
					minDisIndex = 1;
					maxDisIndex = 0;
				}break;
	}
}

int Zipper::generateZipperPoint(PUM::Mesh *srcMesh,PUM::Mesh *destMesh,OverlapPatch &srcPatch,OverlapPatch& destPatch){
	findMaxLoop(srcMesh,sourceEdgeLoop);
	findMaxLoop(destMesh,destEdgeLoop);
	CubeBounds sourceBounds,destBounds;
	//为每一个mesh的边界点产生一个包围盒，并且求出包围盒相交部分顶点距离mesh最近的点
	sourceEdgeLoop.generateBounds(sourceBounds,srcMesh);
	destEdgeLoop.generateBounds(destBounds,destMesh);
	int minDisIndex,maxDisIndex,ignorIndex;
	//确定包围盒相交扩散的主轴和次轴
	findMaxDis(sourceBounds,destBounds,minDisIndex,maxDisIndex,ignorIndex);
	double bound1max = sourceBounds.center[maxDisIndex]+sourceBounds.length[maxDisIndex]/2;
	double bound1min = sourceBounds.center[maxDisIndex]-sourceBounds.length[maxDisIndex]/2;
	double bound2max = destBounds.center[maxDisIndex]+destBounds.length[maxDisIndex]/2;
	double bound2min = destBounds.center[maxDisIndex]+destBounds.length[maxDisIndex]/2;
	PUM::Point3D bound1p1,bound1p2,bound2p1,bound2p2;
	bound1p2[ignorIndex] = bound1p1[ignorIndex] = sourceBounds.center[ignorIndex];
	bound1p1[minDisIndex] = sourceBounds.center[minDisIndex]+sourceBounds.length[minDisIndex]/2;
	bound1p2[minDisIndex] = sourceBounds.center[minDisIndex]-sourceBounds.length[minDisIndex]/2;

	bound2p2[ignorIndex] = bound2p1[ignorIndex] = destBounds.center[ignorIndex];
	bound2p1[minDisIndex] = destBounds.center[minDisIndex]+destBounds.length[minDisIndex]/2;
	bound2p2[minDisIndex] = destBounds.center[minDisIndex]-destBounds.length[minDisIndex]/2;
	if(bound1min<bound2min){
		bound1p2[maxDisIndex] = bound1p1[maxDisIndex] = bound1max;
		bound2p2[maxDisIndex] = bound2p1[maxDisIndex] = bound2min;
	}else{
		bound1p2[maxDisIndex] = bound1p1[maxDisIndex] = bound1min;
		bound2p2[maxDisIndex] = bound2p1[maxDisIndex] = bound2max;
	}
	//将edgeloop上的点压入kdtree，求距离包围盒相交顶点最近的点
	PUM::KDTree *kdtree1 = new PUM::KDTree(3);
	PUM::KDTree *kdtree2 = new PUM::KDTree(3);
	set<PUM::Vertex_Index> vertexs1,vertexs2;
	int sourceNearestIndex1,sourceNearestIndex2,destNearestIndex1,destNearestIndex2;
	sourceEdgeLoop.generateVertexSet(vertexs1);
	sourceEdgeLoop.deleteNotInOverlapArea(vertexs1,srcPatch);
	destEdgeLoop.generateVertexSet(vertexs2);
	destEdgeLoop.deleteNotInOverlapArea(vertexs2,destPatch);
	generateKDtree(vertexs1,srcMesh,kdtree1);
	findNearestEleKdTree(bound1p1,kdtree1,sourceNearestIndex1);
	findNearestEleKdTree(bound1p2,kdtree1,sourceNearestIndex2);
	generateKDtree(vertexs2,destMesh,kdtree2);
	findNearestEleKdTree(bound2p1,kdtree2,destNearestIndex1);
	findNearestEleKdTree(bound2p2,kdtree2,destNearestIndex2);
	srcPatch.idVec.push_back(sourceNearestIndex1);
	srcPatch.idVec.push_back(sourceNearestIndex2);
	destPatch.idVec.push_back(destNearestIndex1);
	destPatch.idVec.push_back(destNearestIndex2);
	delete kdtree1;delete kdtree2;
	return 0;
}



int Zipper::generateNotSameDepthBoundry(){
	return 0;
}

/*
int Zipper::generateNotSameDepthBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop1,EdgeLoop &loop2){
	for()
}*/



int Zipper::generateOverlapPatchVertex(){
	sourcePatch.generateVertexVec(sourceMesh);
	destPatch.generateVertexVec(destMesh);
	return 0;
}

int Zipper::generateKDtree(set<PUM::Vertex_Index> &vertexSet,PUM::Mesh *mesh,PUM::KDTree *newKdtree){
	for(auto it = vertexSet.begin(); it != vertexSet.end(); it++){
		PUM::Point3D currPoint = mesh->Point(*it);
		newKdtree->insert3(currPoint.x(),currPoint.y(),currPoint.z(),*it);
	}
	return 0;
}


int Zipper::findSuitExpandDepth(PUM::Mesh *scMesh,OverlapPatch &scPatch,EdgeLoop &el,PUM::Mesh *oppMesh,OverlapPatch &oppPatch){
	scPatch.generateHierarchicalStructure(scMesh,el);
	PUM::KDTree *oppositeKdTree = new PUM::KDTree(3);
	generateKDtree(oppPatch.vertexSet,oppMesh,oppositeKdTree);
	double errorAccumulative = 0.0;
	int currDepth = 0;
	int count = 0;
	vector<double> errorVec;
	for(auto it = scPatch.hiNodeVec.begin(); it != scPatch.hiNodeVec.end(); it++){
		PUM::Point3D pt = scMesh->Point(it->vertId);
		PUM::kdres* kd_result = oppositeKdTree->nearest3(pt.x(), pt.y(),pt.z());
		//if(!kd_result)
		//	continue;
		double pos[3];
		int temp_id = oppositeKdTree->res_item(kd_result, pos);
		oppositeKdTree->kd_res_free(kd_result);//释放资源
		PUM::Point3D resultPoint(pos[0],pos[1],pos[2]);
		double d = Distance(pt,resultPoint);
		if(it->depth == currDepth){
			errorAccumulative += d;
			count ++;
		}else{
			errorVec.push_back(errorAccumulative/count);
			count = 1;
			errorAccumulative  = d;
			currDepth = it->depth;
		}
	}
	OUTPUT_TRACE_T_2("depth size：%d",errorVec.size());
	delete oppositeKdTree;
	return min_element(errorVec.begin(),errorVec.end()) - errorVec.begin();
}

int Zipper::overLapPatchExtend(PUM::Mesh *scMesh,OverlapPatch &scPatch,int depth){
	set<PUM::Face_Index> deleteFaceSet;
	//vector<PUM::Face_Index> deleteFaceVec;
	for(auto it = scPatch.hiNodeVec.begin(); it != scPatch.hiNodeVec.end(); it++){
		if(it->depth < depth){
			PUM::Face_List faceList;
			vector<PUM::Face_Index> adjFace;
			scMesh->GetFacesOfVertex(it->vertId,faceList);
			convertLinkedlist2Vector(faceList,adjFace);
			deleteFaceSet.insert(adjFace.begin(),adjFace.end());
			//deleteFaceSet.insert(it->vertId);
		}
	}
	//deleteFaceVec.assign(deleteFaceSet.begin(),deleteFaceSet.end());
	//sort(deleteFaceVec.begin(),deleteFaceVec.end());
	OUTPUT_TRACE_T_2("source expand deleteFacevecSize：%d",deleteFaceSet.size());
	for(auto it = deleteFaceSet.rbegin(); it != deleteFaceSet.rend(); it++){
		//OUTPUT_TRACE_T_2("deleteFaceSet:%d",*it);
		scMesh->Delete_a_Face(*it);
	}
	//找出孤立点，并删除之
	deleteIsolatePoints(scMesh);
	deleteLoopCircle(scMesh);
	//deleteNonMainfoldVerticle(scMesh);
	return 0;
}

int Zipper::generateOverLapSingleAreaBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop){
	std::set<int>& triangleSet = patch.triangleSet;
	WriteString2File_Debug("D:\\test.txt",std::to_string((long long)triangleSet.size()));
	
	for(auto tri= triangleSet.begin();tri != triangleSet.end(); tri++){
		PUM::Face& face = mesh->GetFace(*tri);
		int vertexNum = face.GetVertexNum();
		for(int i = 0;i<vertexNum;i++){
			int index0 = face.GetVertex(i);
			int index1 = face.GetVertex((i+1)%vertexNum);
			vector<PUM::Face_Index> faces;
			mesh->GetCommonFacesByVertexes(index0,index1,faces);
			if(faces.size() == 1) continue;
			if((patch.triangleSet.count(faces[0]) && !patch.triangleSet.count(faces[1]))
			|| (patch.triangleSet.count(faces[1]) && !patch.triangleSet.count(faces[0])))
			loop.addEdge(OctreeBase::Edge(index0,index1));
		}
	}
	WriteString2File_Debug("D:\\test.txt",std::to_string((long long)loop.edgeVec.size()));
	return 0;
}

void Zipper::findMaxLoop(PUM::Mesh* mesh,EdgeLoop& maxLoop){
	vector<EdgeLoop> loopVec;
	generateEdgeLoop(mesh,loopVec);
	int loopIndex = findMaxLoop(loopVec);
	if(loopIndex != -1)
		maxLoop = loopVec[loopIndex];
}

int Zipper::generateOverLapPatchBoundry(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop){
	//generateEdgeLoop(sourceMesh,)
	vector<EdgeLoop> loopVec;
	EdgeLoop loopMax,part1,part2;
	generateEdgeLoop(mesh,loopVec);
	int loopIndex = findMaxLoop(loopVec);
	if(loopIndex != -1)
		loopMax = loopVec[loopIndex];
	if(!patch.idVec.empty()){
		OUTPUT_TRACE_T_2("patch.idVec[0]:%d",patch.idVec[0]);
		OUTPUT_TRACE_T_2("patch.idVec[1]:%d",patch.idVec[1]);
		loopMax.partLoop(patch.idVec[0],patch.idVec[1],part1,part2);
	}
	for(auto it = part1.edgeVec.begin(); it != part1.edgeVec.end(); it++){
		int index1 = it->index[0];
		int index2 = it->index[1];
		if(find(patch.vertexSet.begin(), patch.vertexSet.end(),index1) == patch.vertexSet.end()){
			loop = part2;
			return 0;
		}
	}
	loop = part1;
	//loop = loopMax;
	return 0;
}

bool isEdgeBoundry(PUM::Face_Index vertexId1,int vertexId2,PUM::Mesh *mesh){
	PUM::LinkedList<PUM::Face_Index> adjFaceList1,adjFaceList2;
	mesh->GetFacesOfVertex(vertexId1, adjFaceList1);
	mesh->GetFacesOfVertex(vertexId2, adjFaceList2);
	set<PUM::Face_Index> adjFaceSet;
	while(!adjFaceList1.EndOfList()){
		adjFaceSet.insert(adjFaceList1.Data());
		adjFaceList1.Next();
	}
	while(!adjFaceList2.EndOfList()){
		adjFaceSet.insert(adjFaceList2.Data());
		adjFaceList2.Next();
	}
	return (adjFaceList1.Size() + adjFaceList2.Size()- adjFaceSet.size() == 1);
}

bool isVertexBoundry(int vertexId,PUM::Mesh *mesh){
	PUM::LinkedList<PUM::Vertex_Index> adjFaceList;
	mesh->GetAdjVertsOfVertex(vertexId,adjFaceList);
	while( !adjFaceList.EndOfList() )
	{
		if(isEdgeBoundry(adjFaceList.Data(),vertexId,mesh))
			return true;
		adjFaceList.Next();
	}
	return false;
}



int Zipper::eatBack(){
	//eatBack(sourceMesh,destMesh,sourcePatch,destPatch);
	//eatBack(destMesh,sourceMesh,destPatch,sourcePatch);
	sourceMeshEatBack(sourceMesh,sourcePatch,sourceEdgeLoop);
	clear();
	deleteLoopCircle(sourceMesh);
	findOverlappArea();
	eatBackUsing2d(destMesh,sourceMesh,destPatch,sourcePatch);
	return 0;
}

int Zipper::edgeLoopInoverlapArea(PUM::Mesh* mesh,set<PUM::Vertex_Index> &vertexs1,set<PUM::Vertex_Index> &vertexs2,OverlapPatch &patch){
	PUM::KDTree* kdtree = new PUM::KDTree(3);
	for(auto it = patch.vertexSet.begin(); it != patch.vertexSet.end(); it ++){
		PUM::Point3D &pt = mesh->Point(*it);
		kdtree->insert3(pt.x(),pt.y(),pt.z(),*it);
	}
	double d1 = 0,d2 = 0;
	for(auto it = vertexs1.begin(); it != vertexs1.end(); it++){
		int resultId;
		d1 += findNearestEleKdTree(mesh->Point(*it),kdtree,resultId);
	}
	for(auto it = vertexs2.begin(); it != vertexs2.end(); it++){
		int resultId;
		d2 += findNearestEleKdTree(mesh->Point(*it),kdtree,resultId);
	}
	delete kdtree;
	return d1<d2?1:2;
}


void Zipper::sourceMeshEatBack(PUM::Mesh* currMesh,set<PUM::Vertex_Index> &vertexset){
	queue<PUM::Vertex_Index> queues;
	unordered_set<PUM::Vertex_Index> visit;
	unordered_set<PUM::Vertex_Index> shortPath;
	set<PUM::Vertex_Index> vertexs;

	shortPath.insert(shortestPath.begin(),shortestPath.end());

	for(auto it = vertexset.begin(); it != vertexset.end(); it++){
		//if(*it != startIndex && *it != endIndex){
		if(!shortPath.count(*it)){
			queues.push(*it);
			visit.insert(*it);
			vertexs.insert(*it);
		}
	}
	
	
	//vertexs.insert(vertexset.begin(),vertexset.end());
	while(!queues.empty()){
		PUM::Vertex_Index front = queues.front();
		queues.pop();
		vector<PUM::Vertex_Index> adjVertexs;
		currMesh->Get_Adj_Vertex(front,adjVertexs);
		for(auto it = adjVertexs.begin(); it != adjVertexs.end(); it ++){
			if(!shortPath.count(*it) && !visit.count(*it)){
				vertexs.insert(*it);
				queues.push(*it);
				visit.insert(*it);
			}
		}
	}
	//deleteFaces(currMesh,faces);
	deleteVertexs(currMesh,vertexs);
}

void Zipper::deleteVertexs(PUM::Mesh *scMesh,set<PUM::Face_Index>& vertexs){
	for(auto it = vertexs.rbegin(); it != vertexs.rend(); it++){
		scMesh->Delete_a_Vertex(*it);
	}
}

void Zipper::sourceMeshEatBack(PUM::Mesh* currMesh,OverlapPatch &currPatch,EdgeLoop &loop){
	EdgeLoop loop1,loop2,innerLoop;
	loop.partLoop(startIndex,endIndex,loop1,loop2);
	set<PUM::Vertex_Index> vertexset1,vertexset2;
	loop1.generateVertexSet(vertexset1);
	loop2.generateVertexSet(vertexset2);
	
	if(edgeLoopInoverlapArea(currMesh,vertexset1,vertexset2,currPatch) == 1){
		//sourceEdgeLoop = loop1;
		sourceMeshEatBack(currMesh,vertexset1);
	}else {
		//sourceEdgeLoop = loop2;
		sourceMeshEatBack(currMesh,vertexset2);
	}
	OUTPUT_TRACE_T_2("eat edge loop size:%d",sourceEdgeLoop.edgeVec.size());
}

void Zipper::eatBackUsing2d(PUM::Mesh* currMesh,PUM::Mesh* oppositeMesh,OverlapPatch &currPatch,OverlapPatch &oppositePatch){
	vector<PUM::Face_Index> faces;
	faces.assign(currPatch.triangleSet.begin(),currPatch.triangleSet.end());
	//删除面
	deleteFaces(currMesh,faces);
	//找出孤立点，并删除之
	deleteIsolatePoints(currMesh);
	/*vector<EdgeLoop> loopVec;
	generateEdgeLoop(oppositeMesh,loopVec);
	int loopIndex = findMaxLoop(loopVec);
	CCoordTransform::Init("local.viwo.local");
	auto lonlat = CCoordTransform::TransformLonlat(-37.5528, -83.7114, 22.028);
	OUTPUT_TRACE_T_2("jingweidu:%lf%lf",lonlat.first,lonlat.second);
	if(loopIndex != -1){
		//寻找相对patch的最外层边界（由点集构成的多边形）
		EdgeLoop& oppLoop = loopVec[loopIndex];
		const int currSize = oppLoop.edgeVec.size()+1;
		vector<Point> points(currSize+1);
		set<PUM::Vertex_Index> pointset;
		oppLoop.generateVertexSet(pointset);
		for(auto it = pointset.begin(); it != pointset.end(); it++){
			PUM::Point3D &pt = oppositeMesh->Point(*it);
			//转换成经纬度
			auto lonlat = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
			points.push_back(Point(lonlat.first,lonlat.second));
			//WriteString2File_Debug("D:\\loop_00.obj","v " +toString(lonlat.first)+ " "+ toString(lonlat.second) + " 0.0");
			//WriteString2File_Debug("D:\\loop3d_01.obj","v " +toString(pt.x())+ " "+ toString(pt.y()) +" "+ toString(pt.z()));
		}
		points.push_back(points[0]);
		//寻找eate patch上在上边所求多边形内部的面片(至少有一个顶点在多边形内部)，并删除
		currPatch.covertSet2Vector();
		vector<PUM::Face_Index> deletefaceVec;
		for(auto triIt = currPatch.triangleVec.begin();triIt != currPatch.triangleVec.end(); ){
			auto currIt = triIt;
			int currFaceId = *triIt;
			++triIt;
			PUM::Face currFace = currMesh->GetFace(currFaceId);
			int vertexNum = currFace.GetVertexNum();
			bool isDelete = false;
			for(int i = 0;i<3;i++){
				int pointIndex = currFace.GetVertex(i);
				PUM::Point3D& currPoint = currMesh->Point(pointIndex);
				auto lonlat = CCoordTransform::TransformLonlat(currPoint.x(),currPoint.y(),currPoint.z());
				Point curr(lonlat.first,lonlat.second);
				if(point_in_polygon_check_edge(curr,points)){
					isDelete = true;
					break;
				}
			}
			if(isDelete)
				deletefaceVec.push_back(currFaceId);
		}
		//删除面
		deleteFaces(currMesh,deletefaceVec);
		//找出孤立点，并删除之
		deleteIsolatePoints(currMesh);
	}*/

}

int Zipper::eatBack(PUM::Mesh* currMesh,PUM::Mesh* oppositeMesh,OverlapPatch &currPatch,OverlapPatch &oppositePatch){
	PUM::KDTree *oppositeKdTree = new PUM::KDTree(3);
	bool isFirst = true;
	//为oppositePatch的顶点建立kdtree
	oppositePatch.covertSet2Vector();
	for(auto oppIt = oppositePatch.triangleVec.begin(); oppIt != oppositePatch.triangleVec.end(); oppIt ++){
		PUM::Face& currFace = oppositeMesh->GetFace(*oppIt);
		int vertexNum = currFace.GetVertexNum();
		for(int i = 0;i<vertexNum;i++){
			int pointIndex = currFace.GetVertex(i);
			PUM::Point3D& currPoint = oppositeMesh->Point(pointIndex);
			if(isFirst){
				oppositeKdTree->insert3(currPoint.x(),currPoint.y(),currPoint.z(),pointIndex);
				isFirst = false;
			}else{
				PUM::kdres* kd_result = oppositeKdTree->nearest3(currPoint.x(), currPoint.y(),currPoint.z());
				double pos[3];
				int temp_id = oppositeKdTree->res_item(kd_result, pos);
				oppositeKdTree->kd_res_free(kd_result);//释放查找结果
				PUM::Point3D pt(pos[0],pos[1],pos[2]);
				if(Distance(currPoint,pt) > SMALL_NUM){
					oppositeKdTree->insert3(currPoint.x(),currPoint.y(),currPoint.z(),pointIndex);
				}
			}
		}
	}
	
	currPatch.covertSet2Vector();
	vector<PUM::Face_Index> deletefaceVec;
	for(auto triIt = currPatch.triangleVec.begin();triIt != currPatch.triangleVec.end(); /*++triIt*/){
		auto currIt = triIt;
		int currFaceId = *triIt;
		++triIt;
		PUM::Face& currFace = currMesh->GetFace(currFaceId);
		int vertexNum = currFace.GetVertexNum();
		bool isDelete = false;
		for(int i = 0;i<3;i++){
			int pointIndex = currFace.GetVertex(i);
			PUM::Point3D& currPoint = currMesh->Point(pointIndex);
			PUM::kdres* kd_result = oppositeKdTree->nearest3(currPoint.x(), currPoint.y(),currPoint.z());
			double pos[3];
			if(kd_result == nullptr) continue;
			int temp_id = oppositeKdTree->res_item(kd_result, pos);
			oppositeKdTree->kd_res_free(kd_result);
			PUM::Point3D pt(pos[0],pos[1],pos[2]);
			double d = Distance(currPoint,pt);
			if( d < cepsilon){
				isDelete = true;
				break;
			}
		}
		if(isDelete)
			deletefaceVec.push_back(currFaceId);
	}
	delete oppositeKdTree;
	OUTPUT_TRACE_T_2("eat back deleteFaceVec size:%d",deletefaceVec.size());
	//删除面
	deleteFaces(currMesh,deletefaceVec);
	//找出孤立点，并删除之
	deleteIsolatePoints(currMesh);
	return 0;
}

int Zipper::deleteIsolatePoints(PUM::Mesh *scMesh){
	for(int i = scMesh->GetVertexNum()-1;i >= 0;i--){
		PUM::Face_List fl;
		scMesh->GetFacesOfVertex(i,fl);
		if(fl.Size() == 0)
			scMesh->Delete_a_Vertex(i);
	}
	return 0;
}

void Zipper::deleteFaces(PUM::Mesh *scMesh,vector<PUM::Face_Index>& faces){
	sort(faces.begin(),faces.end());
	for(auto it = faces.rbegin(); it != faces.rend(); it++)
		scMesh->Delete_a_Face(*it);
}

void Zipper::deleteFaces(PUM::Mesh *scMesh,set<PUM::Face_Index>& faces){
	for(auto it = faces.rbegin(); it != faces.rend();it++){
		scMesh->Delete_a_Face(*it);
	}
}

int Zipper::cippingInit(PUM::Mesh *mesh,OverlapPatch &op){
	op.generatetBoundryTriangle(mesh);
	for(auto it = op.boundryTriangleVec.begin(); it != op.boundryTriangleVec.end(); it++){
		for(auto edgeIt = it->boundryEdge.begin();edgeIt != it->boundryEdge.end(); edgeIt ++){
			edgeIt->generateTriWall(mesh,op,cepsilon);
		}
	}
	return 0;
}

int Zipper::clipping(){
	
	cippingInit(sourceMesh,sourcePatch);
	cippingInit(destMesh,destPatch);

	return 0;
}

int Zipper::drawEdgeLoop(EdgeLoop &edgeLoop,PUM::Mesh* mesh,PUM::Color& color){
	//OpenGLColor(COLOR_BLUE);
	glPointSize(5.0);
	for(auto it = edgeLoop.edgeVec.begin(); it != edgeLoop.edgeVec.end(); it++){
		int vertexId1 = it->index[0];
		int vertexId2 = it->index[1];
		PUM::Point3D pt1,pt2;
		pt1 = mesh->Point(vertexId1);
		pt2 = mesh->Point(vertexId2);
		DrawPoint(pt1);
		DrawPoint(pt2);
		OpenGLColor(color);
		glLineWidth(5);
		DrawLine(pt1,pt2);
	}
	return 0;
}

int Zipper::drawPatch(OverlapPatch &op,PUM::Mesh *mesh,int color){
	if(color)
		OpenGLColor(COLOR_GREEN);
	else OpenGLColor(COLOR_RED);
	glPointSize(5.0);
	if(!op.isConvert) {
		op.covertSet2Vector();
		op.isConvert = true;
	}
	for(auto it = op.triangleVec.begin();it != op.triangleVec.end(); it++){
		PUM::Face currFace = mesh->GetFace(*it);
		int vertexNum = currFace.GetVertexNum();
		for(int j = 0;j<vertexNum;j++){
			DrawPoint(mesh->Point(currFace.GetVertex(j)));
		}
	}
	return 0;
}

int Zipper::drawHoleBoundry(){
	if(!holeVertexVec.empty()){
		int pre = holeVertexVec.front();
		PUM::Point3D pt = sourceMesh->Point(pre);
		PUM::Point3D endP = sourceMesh->Point(holeVertexVec.back());
		glPointSize(5.0);
		OpenGLColor(COLOR_BLACK);
		DrawPoint(pt);
		OpenGLColor(COLOR_YELLOW);
		DrawPoint(endP);
		for(int i = 1; i < holeVertexVec.size(); i++){
			int curr = holeVertexVec[i];
			PUM::Point3D pt1 = sourceMesh->Point(pre);
			PUM::Point3D pt2 = sourceMesh->Point(curr);
			glLineWidth(3);
			OpenGLColor(COLOR_BLUE);
			DrawLine(pt1,pt2);
			pre = curr;
		}
	}
	return 0;
}

int Zipper::drawPath(){
	if(shortestPath.empty()) return 0;
	glPointSize(5.0);
	OpenGLColor(COLOR_RED);
	for(int i = 0;i<shortestPath.size();i ++){
		DrawPoint(sourceMesh->Point(shortestPath[i]));
	}
	return 0;
}

int Zipper::draw(){
	//drawPatch(sourcePatch,sourceMesh,1);
	//drawPatch(destPatch,destMesh,0);
	//drawEdgeLoop(sourceEdgeLoop,sourceMesh,COLOR_BLUE);
	//drawEdgeLoop(destEdgeLoop,destMesh,COLOR_RED);
	//drawHoleBoundry();
	drawPath();
	OpenGLColor(COLOR_YELLOW);
	glPointSize(5.0);
	DrawPoint(sourceMesh->Point(startIndex));
	DrawPoint(sourceMesh->Point(endIndex));
	
	return 0;
}

int Zipper::getDestPatchSize(){
	return destPatch.triangleVec.size();
}
int Zipper::getSourcePatchSize(){
	return sourcePatch.triangleVec.size();
}

int Zipper::findMaxLoop(vector<EdgeLoop> &loopVec){
	int maxValue = 0;
	int maxIndex = -1;
	for(int i = 0;i< loopVec.size();i++){
		if(loopVec[i].getEdgeNum() >maxValue){
			maxValue = loopVec[i].getEdgeNum();
			maxIndex = i;
		}
	}
	return maxIndex;
}

int Zipper::generateOverLapPatchEdgeLoop(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &edgeLoop){
	vector<EdgeLoop> loopVec;
	EdgeLoop loopMax;
	generateEdgeLoop(mesh,loopVec);
	int loopIndex = findMaxLoop(loopVec);
	if(loopIndex != -1)
		loopMax = loopVec[loopIndex];
	//将loopmax的edgeVec复制一份
	loopMax.edgeVec.insert(loopMax.edgeVec.end(),loopMax.edgeVec.begin(),loopMax.edgeVec.end());
	bool isStart = false;
	for(auto it = loopMax.edgeVec.begin(); it != loopMax.edgeVec.end();it++){
		int vertex1 = it->index[0];
		int vertex2 = it->index[1];
		if(find(patch.vertexSet.begin(),patch.vertexSet.end(),vertex1) == patch.vertexSet.end() ||
			find(patch.vertexSet.begin(),patch.vertexSet.end(),vertex2) == patch.vertexSet.end() ){
				isStart = true;
		}
		if(isStart){
			if(find(patch.vertexSet.begin(),patch.vertexSet.end(),vertex1) != patch.vertexSet.end() ||
				find(patch.vertexSet.begin(),patch.vertexSet.end(),vertex2) != patch.vertexSet.end()){
					edgeLoop.addEdge(*it);
			}else
				break;
		}
	}
	return 0;
}

int Zipper::generateOverLapPatchLoops(PUM::Mesh *mesh,OverlapPatch &patch,EdgeLoop &loop1,EdgeLoop &loop2){
	EdgeLoop boundryLoop;
	generateOverLapPatchBoundry(mesh,patch,sourceEdgeLoop);
	generateOverLapPatchEdgeLoop(mesh,patch,boundryLoop);
	int id1 = patch.idVec[0];
	int id2 = patch.idVec[1];
	int count = 0;

	for(auto it = boundryLoop.edgeVec.begin(); it != boundryLoop.edgeVec.end(); it++){
		int vertex1 = it->index[0];
		int vertex2 = it->index[1];
		if(vertex1 == id1 || vertex2 == id2 || vertex1 == id2 || vertex2 == id1)
			count ++;
		if(count == 0){
			loop1.addEdge(*it);
		}else if(count == 2){
			loop2.addEdge(*it);
		}
	}
	return 0;
}

int Zipper::findHole(){
	EdgeLoop loop1,loop2;
	findMaxLoop(sourceMesh,loop1);
	findMaxLoop(destMesh,loop2);
	//findNearOtherEdgeLoopPart(sourceMesh,destMesh,loop1,loop2,sourceEdgeLoop);
	//findNearOtherEdgeLoopPart(destMesh,sourceMesh,loop2,loop1,destEdgeLoop);
	//findNearOtherEdgeLoopPartUsingLola(sourceMesh,destMesh,loop1,loop2,sourceEdgeLoop);
	//findNearOtherEdgeLoopPartUsingLola(destMesh,sourceMesh,loop2,loop1,destEdgeLoop);
	//sourceEdgeLoop = loop1;
	//destEdgeLoop = loop2;
	findEdgeLoopHolePart(sourceMesh,destMesh,loop1,loop2,sourceEdgeLoop,destEdgeLoop);
	return 0;
}

void Zipper::findEdgeLoopHolePart(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &srcLoop,EdgeLoop &destLoop,EdgeLoop &resultLoop1,EdgeLoop &resultLoop2){
	set<PUM::Vertex_Index> vertexs1,vertexs2;
	srcLoop.generateVertexSet(vertexs1);
	destLoop.generateVertexSet(vertexs2);
	PUM::KDTree *loopKdtree1 = new PUM::KDTree(3);
	PUM::KDTree *loopKdtree2 = new PUM::KDTree(3);
	generateKDtree(vertexs1,mesh1,loopKdtree1);
	generateKDtree(vertexs2,mesh2,loopKdtree2);
	int currStartIndex,currEndIndex;
	EdgeLoop loop1,loop2;
	//找到距离最短路起始点的index,根据index分割edge loop,然后找到距离洞口最近的一部分（由于在删除的时候可能使得最短路起始顶点index发生变化，所以这里增加步骤来重新求）
	findNearestEleKdTree(startPoint,loopKdtree1,currStartIndex);
	findNearestEleKdTree(endPoint,loopKdtree1,currEndIndex);
	srcLoop.partLoop(currStartIndex,currEndIndex,loop1,loop2);
	if(findNearLoop(mesh1,loop1,loop2,loopKdtree2) == 1) resultLoop1 = loop1;
	else resultLoop1 = loop2;
	//同上，求洞的另一侧
	findNearestEleKdTree(startPoint,loopKdtree2,currStartIndex);
	findNearestEleKdTree(endPoint,loopKdtree2,currEndIndex);
	loop1.clear();loop2.clear();
	destLoop.partLoop(currStartIndex,currEndIndex,loop1,loop2);
	if(findNearLoop(mesh2,loop1,loop2,loopKdtree1) == 1) resultLoop2 = loop1;
	else resultLoop2 = loop2;
	delete loopKdtree1;
	delete loopKdtree2;
}

int Zipper::findNearLoop(PUM::Mesh *mesh,EdgeLoop &loop1,EdgeLoop &loop2,PUM::KDTree *loopKdtree){
	double d1 = 0,d2 = 0;
	set<PUM::Vertex_Index> vertexs1,vertexs2;
	loop1.generateVertexSet(vertexs1);
	loop2.generateVertexSet(vertexs2);
	for(auto it = vertexs1.begin();it!= vertexs1.end();it++){
		int id;
		d1 += findNearestEleKdTree(mesh->Point(*it),loopKdtree,id);
	}
	d1 /= vertexs1.size();
	for(auto it = vertexs2.begin();it!= vertexs2.end();it++){
		int id;
		d2 += findNearestEleKdTree(mesh->Point(*it),loopKdtree,id);
	}
	d2 = d2 /= vertexs2.size();
	return d1<d2?1:2;
}

int Zipper::findNearOtherEdgeLoopPartUsingLola(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &loop1,EdgeLoop &loop2,EdgeLoop &loopPart1){
	CCoordTransform::Init(matrixPath);
	//将loop2的顶点压入kdtree
	PUM::KDTree *loopKdtree = new PUM::KDTree(3);
	OctreeBase::Edge edge = loop2.edgeVec[0];
	PUM::Point3D& pt = mesh2->Point(edge.index[0]);
	auto lonlat = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
	loopKdtree->insert3(lonlat.first, lonlat.second, 0, edge.index[0]);
	for(auto it = loop2.edgeVec.begin(); it != loop2.edgeVec.end(); it++){
		PUM::Point3D& pt = mesh2->Point(it->index[1]);
		auto lonlattemp = CCoordTransform::TransformLonlat(pt.x(),pt.y(),pt.z());
		loopKdtree->insert3(lonlattemp.first, lonlattemp.second, 0, it->index[1]);
	}
	//
	bool isStart = false;
	for(auto it = loop1.edgeVec.begin(); it != loop1.edgeVec.end(); it++){
		PUM::Point3D& vertex1 = mesh1->Point(it->index[0]);
		PUM::Point3D& vertex2 = mesh1->Point(it->index[1]);
		auto lonlat1 = CCoordTransform::TransformLonlat(vertex1.x(),vertex1.y(),vertex1.z());
		auto lonlat2 = CCoordTransform::TransformLonlat(vertex2.x(),vertex2.y(),vertex2.z());
		PUM::Point3D query1(lonlat1.first,lonlat1.second,0);
		PUM::Point3D query2(lonlat2.first,lonlat2.second,0);
		int vertexId1,vertexId2;
		double d1 = findNearestEleKdTree(query1,loopKdtree,vertexId1);
		double d2 = findNearestEleKdTree(query2,loopKdtree,vertexId2);
		if(d1 > 10*lltitudeDiffer && d2 > 10*lltitudeDiffer){
			//seedIndex = vertex1;
			isStart = true;
			//break;
		}else{
			if(isStart){
				loopPart1.addEdge(*it);
			}
		}
	}
	for(auto it = loop1.edgeVec.begin(); it != loop1.edgeVec.end(); it++){
		int vertexIndex = it->index[0];
		PUM::Point3D &vertex1 = mesh1->Point(vertexIndex);
		PUM::Point3D lolaPt;
		covert3dPoint2LoLadata(vertex1,lolaPt);
		int resultId;
		double d = findNearestEleKdTree(lolaPt,loopKdtree,resultId);
		if(d <= 10*lltitudeDiffer){
			loopPart1.addEdge(*it);
		}else 
			break;
	}
	delete loopKdtree;
	return 0;
}

void Zipper::covert3dPoint2LoLadata(PUM::Point3D &threedPoint,PUM::Point3D &lolaPoint){
	auto lonlat = CCoordTransform::TransformLonlat(threedPoint.x(),threedPoint.y(),threedPoint.z());
	lolaPoint.x() = lonlat.first;
	lolaPoint.y() = lonlat.second;
	lolaPoint.z() = 0;
}

int Zipper::findNearOtherEdgeLoopPart(PUM::Mesh *mesh1,PUM::Mesh *mesh2,EdgeLoop &loop1,EdgeLoop &loop2,EdgeLoop &loopPart1){
	//将loop2的顶点压入kdtree
	PUM::KDTree *loopKdtree = new PUM::KDTree(3);
	OctreeBase::Edge edge = loop2.edgeVec[0];
	PUM::Point3D pt = mesh2->Point(edge.index[0]);
	loopKdtree->insert3(pt.x(), pt.y(), pt.z(), edge.index[0]);
	for(auto it = loop2.edgeVec.begin(); it != loop2.edgeVec.end(); it++){
		PUM::Point3D pt = mesh2->Point(it->index[1]);
		loopKdtree->insert3(pt.x(), pt.y(), pt.z(), it->index[1]);
	}
	//
	bool isStart = false;
	for(auto it = loop1.edgeVec.begin(); it != loop1.edgeVec.end(); it++){
		PUM::Point3D vertex1 = mesh1->Point(it->index[0]);
		PUM::Point3D vertex2 = mesh1->Point(it->index[1]);
		int vertexId;
		double d1 = findNearestEleKdTree(vertex1,loopKdtree,vertexId);
		double d2 = findNearestEleKdTree(vertex2,loopKdtree,vertexId);
		if(d1 > 2*cepsilon && d2 > 2*cepsilon){
			//seedIndex = vertex1;
			isStart = true;
			//break;
		}else{
			if(isStart){
				loopPart1.addEdge(*it);
			}
		}
	}
	for(auto it = loop1.edgeVec.begin(); it != loop1.edgeVec.end(); it++){
		int vertexIndex = it->index[0];
		PUM::Point3D vertex1 = mesh1->Point(vertexIndex);
		PUM::kdres* kd_result = loopKdtree->nearest3(vertex1.x(), vertex1.y(),vertex1.z());
		//if(!kd_result)
		//	continue;
		double pos[3];
		int temp_id = loopKdtree->res_item(kd_result, pos);
		loopKdtree->kd_res_free(kd_result);//释放资源
		PUM::Point3D resultPoint(pos[0],pos[1],pos[2]);
		double d = Distance(vertex1,resultPoint);
		if(d <= 2*cepsilon){
			loopPart1.addEdge(*it);
		}else 
			break;
	}
	delete loopKdtree;
	return 0;
}


int Zipper::getBoundryEdgeNextPoint(PUM::Mesh *mesh,PUM::Vertex_Index vi,vector<PUM::Vertex_Index>& adjVertexs,unordered_map<int,int> &flag){
	for(auto it = adjVertexs.begin(); it != adjVertexs.end(); it++){
		if(flag[*it] == 0){
			vector<PUM::Face_Index> neiFaceVec;
			mesh->GetCommonFacesByVertexes(*it,vi,neiFaceVec);
			if(/*mesh->isEdgeBoundary(vi,*it)*/neiFaceVec.size() == 1){
				//flag[*it] = 1;
				return *it;
			}
		}
	}
	return -1;
}

int getAdjVertexsOfVertex(PUM::Mesh *mesh,PUM::Vertex_Index vi,vector<PUM::Vertex_Index> &adjVertexsVec){
	set<PUM::Vertex_Index>  adjVertexSet;
	vector<PUM::Face_Index> adjFaces;
	PUM::Face_List adjFaceList;
	mesh->GetFacesOfVertex(vi,adjFaceList);
	while(!adjFaceList.EndOfList()){
		adjFaces.push_back(adjFaceList.Data());
		adjFaceList.Next();
	}

	for(int i = 0; i<adjFaces.size();i++){
		PUM::Face& face = mesh->GetFace(adjFaces[i]);
		int faceNum = face.GetVertexNum();
		for(int j =0;j<faceNum;j++){
			if(vi != face.Vertex(j))
				adjVertexSet.insert(face.Vertex(j));
		}
	}
	adjVertexsVec.assign(adjVertexSet.begin(),adjVertexSet.end());
	//mesh->getvertex
	return 0;
}


int Zipper::generateOneEdgeLoop(PUM::Mesh *mesh,PUM::Vertex_Index vi,unordered_map<int,int> &flag,EdgeLoop &loop){
	vector<PUM::Vertex_Index> adjVertexsVec;
	getAdjVertexsOfVertex(mesh,vi,adjVertexsVec);
	int pre = vi;
	int result = getBoundryEdgeNextPoint(mesh,pre,adjVertexsVec,flag);
	while(result != -1){
		OctreeBase::Edge edge(pre,result);
		loop.addEdge(edge);
		flag[result] = 1;
		pre = result;
		adjVertexsVec.clear();
		//mesh->Get_Adj_Vertex(pre,adjVertexsVec);
		getAdjVertexsOfVertex(mesh,pre,adjVertexsVec);
		//mesh->GetAdjVertsOfVertex_kneighbor(pre,adjVertexsVec);
		result = getBoundryEdgeNextPoint(mesh,pre,adjVertexsVec,flag);
	}
	return 0;
}

int Zipper::generateEdgeLoop(PUM::Mesh *mesh,vector<EdgeLoop> &edgeLoopVec){
	//const PUM::Array<PUM::Point3D> &ptArray = mesh->GetPtList();
	while(deleteLoopCircle(mesh));
	//vector<int> visitFlag(mesh->GetVertexNum(),0);
	unordered_map<int,int> visitFlag;
	for(int i = 0;i<mesh->GetVertexNum();i++){
		if(visitFlag[i] == 0){
			if(mesh->isVertexBoundary(i)){
				EdgeLoop loop;
				visitFlag[i] = 1;
				generateOneEdgeLoop(mesh,i,visitFlag,loop);
				if(loop.getEdgeNum() != 0)
					edgeLoopVec.push_back(loop);
			}
		}
	}
	return 0;
}

int Zipper::getBoundryEdgeNextPoint(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,vector<PUM::Vertex_Index>& adjVertexs,vector<PUM::Vertex_Index>& nextPoints){
	for(auto it = adjVertexs.begin(); it != adjVertexs.end(); it++){
		if(*it == pre) continue;
		vector<PUM::Face_Index> neiFaceVec;
		mesh->GetCommonFacesByVertexes(*it,vi,neiFaceVec);
		if(/*mesh->isEdgeBoundary(vi,*it)*/neiFaceVec.size() == 1){
			nextPoints.push_back(*it);
		}
	}
	return 0;
}

bool Zipper::dfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteVertexs){
	if(flag[vi]==-1){
		if(vi == path[0]) return true;
		int i = path.size() -1;
		while( i>=0 && path[i] != vi){
			deleteVertexs.insert(path[i]);
			i--;
		}
		WriteString2File_Debug("D:\\test.txt",std::to_string((long long)mesh->GetVertexNum()));
		return false;
	}
	path.push_back(vi);
	flag[vi] = -1;
	vector<PUM::Vertex_Index> adjVertexsVec,nextPoints;
	//mesh->GetAdjVertsOfVertex_kneighbor(vi,adjVertexsVec);
	//mesh->Get_Adj_Vertex(vi,adjVertexsVec);
	getAdjVertexsOfVertex(mesh,vi,adjVertexsVec);
	getBoundryEdgeNextPoint(mesh,pre,vi,adjVertexsVec,nextPoints);
	for(int i = 0;i<nextPoints.size();i++){
		if(flag[nextPoints[i]] == 1) continue;
		if(!dfsDeleteCircle(mesh,vi,nextPoints[i],flag,path,deleteVertexs)) {
			flag[vi] = 1;
			path.pop_back();
			return false;
		}
	}
	flag[vi] = 1;
	path.pop_back();
	return true;
}

//bool Zipper::dfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteVertexs){
//	stack<pair<PUM::Vertex_Index,PUM::Vertex_Index>> stacks;
//	stacks.push(make_pair(pre,vi));
//	while(!stacks.empty()){
//		if(flag[vi]==-1){
//			if(vi == path[0]) return true;
//			int i = path.size() -1;
//			while( i>=0 && path[i] != vi){
//				deleteVertexs.insert(path[i]);
//				i--;
//			}
//			//WriteString2File_Debug("D:\\test.txt",std::to_string((long long)mesh->GetVertexNum()));
//			return false;
//		}
//		path.push_back(vi);
//		flag[vi] = -1;
//		vector<PUM::Vertex_Index> adjVertexsVec,nextPoints;
//		//mesh->GetAdjVertsOfVertex_kneighbor(vi,adjVertexsVec);
//		//mesh->Get_Adj_Vertex(vi,adjVertexsVec);
//		getAdjVertexsOfVertex(mesh,vi,adjVertexsVec);
//		getBoundryEdgeNextPoint(mesh,pre,vi,adjVertexsVec,nextPoints);
//		for(int i = 0;i<nextPoints.size();i++){
//			if(flag[nextPoints[i]] == 1) continue;
//			stacks.push(make_pair(vi,nextPoints[i]));
//		}
//		flag[vi] = 1;
//		path.pop_back();
//	}
//	return true;
//}

PUM::Vertex_Index findSameParent(PUM::Vertex_Index v1,PUM::Vertex_Index v2,unordered_map<PUM::Vertex_Index,int>& pres){
	unordered_set<int> visitVertex;
	int temp1 = v1,temp2 = v2;
	while(temp1 != -1){
		visitVertex.insert(temp1);
		temp1 = pres[temp1];
	}
	while(temp2 != -1){
		if(visitVertex.count(temp2))
			return temp2;
		temp2 = pres[temp2];
	}
	return -1;
}

bool Zipper::bfsDeleteCircle(PUM::Mesh *mesh,PUM::Vertex_Index pre,PUM::Vertex_Index vi,unordered_map<int,int> &flag,unordered_map<PUM::Vertex_Index,int>& pres,set<PUM::Vertex_Index>& deleteVertexs){
	queue<PUM::Vertex_Index> vertexs;
	//unordered_map<PUM::Vertex_Index,int> visitFlag;
	
	vertexs.push(vi);
	flag[vi] = -1;
	pres[vi] = pre;
	while(!vertexs.empty()){
		PUM::Vertex_Index top = vertexs.front();
		vertexs.pop();
		flag[top] = 1;
		vector<PUM::Vertex_Index> adjVertexsVec,nextPoints;
		//mesh->GetAdjVertsOfVertex_kneighbor(vi,adjVertexsVec);
		//mesh->Get_Adj_Vertex(vi,adjVertexsVec);
		getAdjVertexsOfVertex(mesh,top,adjVertexsVec);
		getBoundryEdgeNextPoint(mesh,pres[top],top,adjVertexsVec,nextPoints);
		for(int i = 0;i<nextPoints.size();i++){
			if(flag[nextPoints[i]] == 1||flag[nextPoints[i]] == -1) {
				PUM::Vertex_Index index1 = top,index2 = nextPoints[i];
				PUM::Vertex_Index pa = findSameParent(index1,index2,pres);
				vector<PUM::Vertex_Index> tempVertexs;
				//if(pres[pa] != -1){
					while(index1 != pa){
						//deleteVertexs.insert(index1);
						tempVertexs.push_back(index1);
						index1 = pres[index1];
					}
					while(index2 != pa){
						tempVertexs.push_back(index2);
						//deleteVertexs.insert(index2);
						index2 = pres[index2];
					}
				//}
				if(tempVertexs.size() < smallPartNum && pres[pa] == -1)
					deleteVertexs.insert(tempVertexs.begin(),tempVertexs.end());
				if(pres[pa] != -1)
					deleteVertexs.insert(tempVertexs.begin(),tempVertexs.end());
				continue;
			}
			vertexs.push(nextPoints[i]);
			flag[nextPoints[i]] = -1;
			pres[nextPoints[i]] = top;
		}
	}
	return true;
}

void Zipper::deleteOneLoopCircle(PUM::Mesh *mesh,PUM::Vertex_Index vi,unordered_map<int,int> &flag,vector<PUM::Vertex_Index>& path,set<PUM::Vertex_Index>& deleteVertexs){
	flag[vi] = -1;
	path.push_back(vi);
	vector<PUM::Vertex_Index> adjVertexsVec;
	getAdjVertexsOfVertex(mesh,vi,adjVertexsVec);
	//mesh->GetAdjVertsOfVertex_kneighbor(vi,adjVertexsVec);
	//mesh->Get_Adj_Vertex(vi,adjVertexsVec);
	//int curr = getBoundryEdgeNextPoint(mesh,vi,adjVertexsVec,flag);
	for(int i = 0;i<adjVertexsVec.size();i++){
		if(mesh->isEdgeBoundary(vi,adjVertexsVec[i]) && flag[adjVertexsVec[i]] == 0){
			dfsDeleteCircle(mesh,vi,adjVertexsVec[i],flag,path,deleteVertexs);
		}
	}
	//if(curr != -1 && flag[curr] == 0)
	//	dfsDeleteCircle(mesh,vi,curr,flag,path,deleteVertexs);
	flag[vi] = 1;
	path.pop_back();
}



bool Zipper::deleteLoopCircle(PUM::Mesh *mesh){
	unordered_map<int,int> visitFlag;
	unordered_map<PUM::Vertex_Index,int> pres;
	set<PUM::Vertex_Index> deleteVertexs;
	for(int i = 0;i<mesh->GetVertexNum();i++){
		if(visitFlag[i] == 0){
			if(mesh->isVertexBoundary(i)){
				vector<PUM::Vertex_Index> path;
				//deleteOneLoopCircle(mesh,i,visitFlag,path,deleteVertexs);
				bfsDeleteCircle(mesh,-1,i,visitFlag,pres,deleteVertexs);
				/*if(pres.count(i))
					bfsDeleteCircle(mesh,pres[i],i,visitFlag,pres,deleteVertexs);
				else
					bfsDeleteCircle(mesh,-1,i,visitFlag,pres,deleteVertexs);*/
			}
		}
	}
	OUTPUT_TRACE_T_2("CircleSize：%d",deleteVertexs.size());
	for(auto it = deleteVertexs.rbegin(); it != deleteVertexs.rend(); it++){
		mesh->Delete_a_Vertex(*it);
	}
	//找出孤立点，并删除之
	deleteIsolatePoints(mesh);
	return !deleteVertexs.empty();
}

bool Zipper::deleteNonMainfoldVerticle(PUM::Mesh *mesh){
	int vertexNum = mesh->GetVertexNum();
	set<int> deleteVerticles;
	for(int i = 0;i<vertexNum;i++){
		if(mesh->isVertexNonMainfold(i)){
			deleteVerticles.insert(i);
		}
	}
	for(auto it = deleteVerticles.rbegin(); it!=deleteVerticles.rend(); it++){
		vector<PUM::Vertex_Index> adjVertexs;
		//mesh->GetAdjVertsOfVertex_kneighbor(*it,adjVertexs);
		//mesh->Get_Adj_Vertex(*it,adjVertexs);
		getAdjVertexsOfVertex(mesh,*it,adjVertexs);
	}
	return !deleteVerticles.empty();
}


double Zipper::findNearestEleKdTree(PUM::Point3D &queryPoint,PUM::KDTree *queryKdTree,int &resultId){
	PUM::kdres* kd_result = queryKdTree->nearest3(queryPoint.x(), queryPoint.y(),queryPoint.z());
	//if(!kd_result)
	//	continue;
	double pos[3];
	resultId = queryKdTree->res_item(kd_result, pos);
	queryKdTree->kd_res_free(kd_result);//消除内存泄漏
	PUM::Point3D resultPoint(pos[0],pos[1],pos[2]);
	return Distance(queryPoint,resultPoint);
}

int Zipper::mergeMesh(){
	OUTPUT_TRACE_T_2("before mesh merge src face size：%d",sourceMesh->GetFaceNum());
	mergeMesh(sourceMesh,destMesh,destEdgeLoop);
	//destMesh = sourceMesh;
	OUTPUT_TRACE_T_2("after mesh merge src face size：%d",sourceMesh->GetFaceNum());
	return 0;
}

int Zipper::mergeMesh(PUM::Mesh *scMesh,PUM::Mesh *deMesh,EdgeLoop &el){
	PUM::KDTree *sourceKdTree = new PUM::KDTree(3);
	const PUM::Array<PUM::Point3D> &scArray = scMesh->GetPtList();
	PUM::Size_Type sourcePtSize = scMesh->GetVertexNum();
	for(int i = 0;i < sourcePtSize;i++){
		sourceKdTree->insert3(scArray[i].x(),scArray[i].y(),scArray[i].z(),i);
	}
	PUM::Size_Type deFaceSize = deMesh->GetFaceNum();
	for(int faceId = 0;faceId < deFaceSize;faceId++){
		const PUM::Face& face = deMesh->GetFace(faceId);
		PUM::Face insertFace;
		int vertexNum = face.GetVertexNum();
		insertFace.Create(vertexNum);
		for(int faceVertexId = 0; faceVertexId < vertexNum; faceVertexId++){
			int pointId = face.GetVertex(faceVertexId);
			PUM::Point3D currPoint = deMesh->Point(pointId);
			int temp_id;
			double d = findNearestEleKdTree(currPoint,sourceKdTree,temp_id);
			if( d > SMALL_NUM){
				int id = scMesh->Add_a_Vertex(currPoint);
				sourceKdTree->insert3(currPoint.x(),currPoint.y(),currPoint.z(),id);
				insertFace.Vertex(faceVertexId) = id;
			}else
				insertFace.Vertex(faceVertexId) = temp_id;
		}
		scMesh->Add_a_Face(insertFace);
	}
	//修改edgeLoop的编号
	for(auto it = el.edgeVec.begin(); it != el.edgeVec.end(); it++){
		int id1,id2;
		int index1 = it->index[0];
		PUM::Point3D pt1 = deMesh->Point(index1);
		double d = findNearestEleKdTree(pt1,sourceKdTree,id1);
		if(d < SMALL_NUM) it->index[0] = id1;
		int index2 = it->index[1];
		PUM::Point3D pt2 = deMesh->Point(index2);
		findNearestEleKdTree(pt2,sourceKdTree,id2);
		d = findNearestEleKdTree(pt2,sourceKdTree,id2);
		if(d < SMALL_NUM) it->index[1] = id2;
	}

	delete sourceKdTree;
	return 0;
}

bool Zipper::isHoleDirectionRight(PUM::Mesh *mesh,int index0,int index1){
	vector<PUM::Face_Index> faces;
	mesh->GetCommonFacesByVertexes(index0,index1,faces);
	PUM::Face &face = mesh->GetFace(faces[0]);
	int vertexNum = face.GetVertexNum();
	for(int i = 0;i<vertexNum;i++){
		if(face.GetVertex(i) == index0 && face.GetVertex((i+1)%vertexNum) == index1)
			return true;
	}
	return false;
}

int Zipper::closeHole(PUM::Mesh *mesh){
	//求洞的边界
	int vertexIndex = sourceEdgeLoop.edgeVec[0].index[0];
	holeVertexVec.push_back(vertexIndex);
	for(auto it = sourceEdgeLoop.edgeVec.begin(); it != sourceEdgeLoop.edgeVec.end(); it++){
		int vertexIndex = it->index[1];
		holeVertexVec.push_back(vertexIndex);
	}
	PUM::Vector3D v1,v2;
	sourceEdgeLoop.getDirection(mesh,v1);
	destEdgeLoop.getDirection(mesh,v2);
	double re = v1 *v2;
	if(re < 0){        
		holeVertexVec.push_back(destEdgeLoop.edgeVec[0].index[0]);
		for(auto it = destEdgeLoop.edgeVec.begin(); it != destEdgeLoop.edgeVec.end(); it++){
			holeVertexVec.push_back(it->index[1]);
		}
	}else{
		holeVertexVec.push_back(destEdgeLoop.edgeVec.back().index[1]);
		for(auto it = destEdgeLoop.edgeVec.rbegin(); it != destEdgeLoop.edgeVec.rend(); it++){
			holeVertexVec.push_back(it->index[0]);
		}
	}
	//PUM::Vector3D pointNormal;
	//mesh->ComputeNormal();
	//getPointsNormal(mesh,pointNormal);
	//double rr = pointNormal *(v1^v2);
	int vertexIndex0 = holeVertexVec[0];
	int vertexIndex1 = holeVertexVec[1];
	if(isHoleDirectionRight(mesh,vertexIndex0,vertexIndex1)) 
		reverse(holeVertexVec.begin(),holeVertexVec.end());
	//补洞
	EasyHoleFilling ehf;
	ehf.Fill_A_Hole(mesh,holeVertexVec,true);
	//mesh->ComputeFaceNormal();
	mesh->ComputeNormal();
	return 0;
}

int Zipper::clear(){
	sourcePatch.clear();
	destPatch.clear();
	sourceEdgeLoop.clear();
	destEdgeLoop.clear();
	return 0;
}


void Zipper::meshMergeManager(){
	findMesh2dIntersectionPoint(sourceMesh,destMesh);
	findOverlappArea();
	findPath(sourceMesh,destMesh,sourcePatch,destPatch);
	eatBack();
	findHole();
	mergeMesh();
	closeHole(sourceMesh);
}

PUM::Mesh* Zipper::meshMerge(PUM::Mesh* meshs[],int n,const char *path){
	return meshMerge(meshs,0,n-1,path);
}

PUM::Mesh* Zipper::meshMerge(PUM::Mesh* meshs[],int left,int right,const char *path){
	if(right < left) return nullptr;
	if(right - left==1){
		Zipper zipper(meshs[left],meshs[right]);
		zipper.setMatrixPath(path);
		zipper.meshMergeManager();	
		return zipper.getSourceMesh();
	}else if(left == right)
		return meshs[left];
	int middle = (left + right)/2;
	PUM::Mesh *leftMesh = meshMerge(meshs,left,middle,path);
	PUM::Mesh *rightMesh = meshMerge(meshs,middle+1,right,path);
	Zipper zipper(leftMesh,rightMesh);
	zipper.setMatrixPath(path);
	zipper.meshMergeManager();	
	return zipper.getSourceMesh();
}

OverlapPatch::OverlapPatch(){
	isConvert = false;
}

int OverlapPatch::modifyPatchVertex(int vertex1,int vertex2){
	if(!isConvert)
		covertSet2Vector();
	vector<int>::iterator it = std::find(triangleVec.begin(),triangleVec.end(),vertex1);
	if(it != triangleVec.end())
		*it = vertex2;
	return 0;
}

int OverlapPatch::covertSet2Vector(){
	if(!isConvert){
		triangleVec.assign(triangleSet.begin(),triangleSet.end());
		sort(triangleVec.begin(),triangleVec.end());
		triangleSet.clear();
		isConvert = true;
	}
	return 0;
}

int OverlapPatch::deleteListEle(int vertexId){
	if(!isConvert)
		covertSet2Vector();
	vector<int>::iterator it = std::find(triangleVec.begin(),triangleVec.end(),vertexId);
	if(it != triangleVec.end()){
		triangleVec.erase(it); 
		return 1;
	}
	return 0;
}

int OverlapPatch::generatetBoundryTriangle(PUM::Mesh *mesh){
	for(set<int>::iterator it = triangleSet.begin(); it != triangleSet.end(); it++){
		PUM::Face face = mesh->GetFace(*it);
		int vertexNum = face.GetVertexNum();
		MergeBoundryTriangle mbt;
		mbt.triangleIndex = *it;
		for(int i = 0; i < vertexNum; i++){
			int j = (i + 1) % vertexNum;
			int vertexId1 = face.GetVertex(i);
			int vertexId2 = face.GetVertex(j);
			if(mesh->isEdgeBoundary(vertexId1,vertexId2)){
				MergeBoundryEdge mbe(vertexId1,vertexId2,boundryTriangleVec.size());
				mbt.boundryEdge.push_back(mbe);
			}
		}
	}
	return 0;
}

int OverlapPatch::generateVertexVec(PUM::Mesh *mesh){
	//set<int> tempVertexSet;
	if(isConvert){
		for(auto it = triangleVec.begin(); it != triangleVec.end(); it++){
			PUM::Face currFace = mesh->GetFace(*it);
			int vertexNum = currFace.GetVertexNum();
			for(int i = 0;i<vertexNum;i++){
				int pointIndex = currFace.GetVertex(i);
				vertexSet.insert(pointIndex);
			}
		}
	}else{
		for(auto it = triangleSet.begin(); it != triangleSet.end(); it++){
			PUM::Face currFace = mesh->GetFace(*it);
			int vertexNum = currFace.GetVertexNum();
			for(int i = 0;i<vertexNum;i++){
				int pointIndex = currFace.GetVertex(i);
				vertexSet.insert(pointIndex);
			}
		}
	}
	//vertexVec.assign(tempVertexSet.begin(),tempVertexSet.end());
	return 0;
}

void OverlapPatch::gengerateBounds(CubeBounds& bound,PUM::Mesh* mesh){
	if(vertexSet.empty()) return;
	int index1 = *vertexSet.begin();
	PUM::Point3D &pt = mesh->Point(index1);
	double minx = pt.x(),miny = pt.y() ,minz = pt.z();
	double maxx = pt.x(),maxy = pt.y() ,maxz = pt.z();
	for(auto it = vertexSet.begin(); it!=vertexSet.end(); it++){
		PUM::Point3D &currPoint = mesh->Point(*it);
		minx = min(minx,currPoint.x());
		maxx = max(maxx,currPoint.x());
		miny = min(miny,currPoint.y());
		maxy = max(maxy,currPoint.y());
		minz = min(minz,currPoint.z());
		maxz = max(maxz,currPoint.z());
	}
	PUM::Point3D center((minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2);
	PUM::Vector3D v(maxx-minx,maxy-miny,maxz-minz);
	bound.setBounds(center,v);
}

int OverlapPatch::testAdjVec(vector<PUM::Vertex_Index> &vertexVec,map<int,int> &flag){
	vector<PUM::Vertex_Index> tempVec;
	for(auto it = vertexVec.begin(); it != vertexVec.end(); it++){
		if(!flag[*it] && vertexSet.count(*it)){
			tempVec.push_back(*it);
		}
	}
	swap(vertexVec,tempVec);
	return 0;
}

int OverlapPatch::clear(){
	this->boundryTriangleVec.clear();
	this->vertexSet.clear();
	this->triangleSet.clear();
	this->triangleVec.clear();
	this->hiNodeVec.clear();
	//this->idVec.clear();
	this->isConvert = false;
	return 0;
}

int OverlapPatch::generateHierarchicalStructure(PUM::Mesh *mesh,EdgeLoop &el){

	/*vector<int> flag(vertexSet.size(),0);*/
	map<int,int> flag;
	//for(auto it = vertexSet.begin(); it != vertexSet.end(); it++){
	//	flag[*it] = 0;
	//}
	//int index0 = el.edgeVec.front().index[0];
	//HierarchicalNode hn(index0,0);
	//hiNodeVec.push_back(hn);
	queue<HierarchicalNode>	tempHiQue;
	for(auto it = el.edgeVec.begin(); it != el.edgeVec.end(); it++){
		int index0 = it->index[0];
		int index1 = it->index[1];
		if(!flag[index0]){
			tempHiQue.push(HierarchicalNode(index0,0));
			flag[index0] = 2;
		}
		if(!flag[index1]){
			tempHiQue.push(HierarchicalNode(index1,0));
			flag[index1] = 2;
		}
	}
	while(!tempHiQue.empty()){
		HierarchicalNode hn = tempHiQue.front();
		tempHiQue.pop();
		//if(flag[hn.vertId]) continue;
		hiNodeVec.push_back(hn);
		flag[hn.vertId] = 1;
		vector<PUM::Vertex_Index> vertexVec;
		mesh->GetAdjVertsOfVertex(hn.vertId,vertexVec);
		testAdjVec(vertexVec,flag);
		for(auto it = vertexVec.begin(); it != vertexVec.end(); it++){
			flag[*it] = 1;
			tempHiQue.push(HierarchicalNode(*it,hn.depth + 1));
		}
	}
	return 0;
}

int MergeBoundryEdge::generateWallVertex(PUM::Mesh *mesh,OverlapPatch &op,double epsilon,int vertexId,PUM::Point3D &v1,PUM::Point3D &v2){
	int triangleId = getVertexAdjBoundryEdge(vertexId,mesh);
	PUM::Vector3D vt1 = mesh->GetFaceNormal(triangleId);
	int triangleId2 = op.boundryTriangleVec[adjTrianglePatchId].triangleIndex;
	PUM::Vector3D vt2 = mesh->GetFaceNormal(triangleId2);
	PUM::Point3D pt = mesh->Point(vertexId);
	PUM::Vector3D dir = vt1 + vt2;
	v1 = pt - dir*epsilon;
	v2 = pt + dir*epsilon;
	return 0;
}

int MergeBoundryEdge::generateTriWall(PUM::Mesh *mesh,OverlapPatch &op,double epsilon){
	PUM::Point3D c1,c2,c3,c4,v1,v2;
	v1 = mesh->Point(vertex1);
	v2 = mesh->Point(vertex2);
	generateWallVertex(mesh,op,epsilon,vertex1,c1,c2);
	generateWallVertex(mesh,op,epsilon,vertex2,c3,c4);
	triWall[0].setAllVertexs(v1,c1,c3);
	triWall[1].setAllVertexs(v1,c3,v2);
	triWall[2].setAllVertexs(c2,v1,c4);
	triWall[3].setAllVertexs(v1,v2,c4);
	return 0;
}

int MergeBoundryEdge::getVertexAdjBoundryEdge(int v,PUM::Mesh *mesh){
	PUM::LinkedList<PUM::Vertex_Index> adjVertexList;
	mesh->GetAdjVertsOfVertex(v,adjVertexList);
	while(!adjVertexList.EndOfList()){
		//PUM::Face_Index fi[2];
		//int fiNum;
		vector<PUM::Face_Index> neighFaces;
		if(adjVertexList.Data() != vertex2){
			//mesh->GetFaceByEdge(vertex1,adjVertexList.Data(),fi,fiNum);
			mesh->GetCommonFacesByVertexes(vertex1,adjVertexList.Data(),neighFaces);
			if(neighFaces.size() == 1){
				return neighFaces[0];
			}
		}
	}
	return -1;
}

int EdgeLoop::getDirection(PUM::Mesh *mesh,PUM::Vector3D &dirSum){
	//PUM::Vector3D dirSum(0,0,0);
	for(auto it = edgeVec.begin(); it != edgeVec.end(); it++){
		PUM::Point3D pt1,pt2;
		pt1 = mesh->Point(it->index[0]);
		pt2 = mesh->Point(it->index[1]);
		PUM::Vector3D dir = pt2-pt1;
		dir.SetUnit();
		dirSum += dir;
	}
	dirSum.SetUnit();
	return 0;
}



int EdgeLoop::partLoop(int e1,int e2,EdgeLoop &el1,EdgeLoop &el2){
	bool isFirst = false,isSecond = false;
	for(auto it = edgeVec.begin(); it != edgeVec.end(); it++){
		int id1 = it->index[0];
		int id2 = it->index[1];
		if(isFirst && (id1 == e1 || id1 == e2)){
			isFirst = false;
			isSecond = true;
		}
		if(!isFirst && !isSecond && (id1 == e1 || id1 == e2)){
			isFirst = true;
		}
		if(isFirst){
			el1.addEdge(*it);
		}
		if(isSecond){
			el2.addEdge(*it);
		}
		if(id2 == e1 || id2 == e2){
			if(it == --edgeVec.end()) return 0;
		}
		
	}
	for(auto it = edgeVec.begin(); it != edgeVec.end();it ++){
		int id1 = it->index[0];
		int id2 = it->index[1];
		el2.addEdge(*it);
		if (id2 == e1 || id2 == e2 )
			break;
	}
	return 0;
}

int EdgeLoop::clear(){
	this->edgeVec.clear();
	return 0;
}

void EdgeLoop::deleteNotInOverlapArea(set<PUM::Vertex_Index>& ptSet,OverlapPatch &curr){
	vector<PUM::Vertex_Index> deleteSets;
	for(auto it = ptSet.begin();it != ptSet.end();it++){
		if(!curr.vertexSet.count(*it))
			deleteSets.push_back(*it);
	}
	for(auto it = deleteSets.begin(); it != deleteSets.end(); it++){
		ptSet.erase(*it);
	}
}


void EdgeLoop::generateBounds(CubeBounds& bound,PUM::Mesh *mesh){
	PUM::Point3D &pt = mesh->Point(edgeVec[0].index[0]);
	double minx = pt.x(),miny = pt.y() ,minz = pt.z();
	double maxx = pt.x(),maxy = pt.y() ,maxz = pt.z();
	for(int i = 0;i<edgeVec.size();i++){
		PUM::Point3D &currPoint = mesh->Point(edgeVec[i].index[1]);
		minx = min(minx,currPoint.x());
		maxx = max(maxx,currPoint.x());
		miny = min(miny,currPoint.y());
		maxy = max(maxy,currPoint.y());
		minz = min(minz,currPoint.z());
		maxz = max(maxz,currPoint.z());
	}
	PUM::Point3D center((minx+maxx)/2,(miny+maxy)/2,(minz+maxz)/2);
	PUM::Vector3D v(maxx-minx,maxy-miny,maxz-minz);
	bound.setBounds(center,v);
}

void EdgeLoop::generateVertexSet(set<PUM::Vertex_Index>& ptset){
	for(auto it = edgeVec.begin();it!=edgeVec.end();it++){
		ptset.insert(it->index[0]);
		ptset.insert(it->index[1]);
	}
}
