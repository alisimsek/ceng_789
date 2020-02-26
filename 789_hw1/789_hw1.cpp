#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include <queue> 
#include <utility>
#include <map>

#include "Painter.h"
#include "fibonacci.hpp"

void drawThickEdges(SoSeparator*, Mesh*, std::vector<Edge*>);
void dijkstra(Mesh*, int ,int);
void calculateThickEdges(Mesh *,int*, int, int);
void printDistToFile(FILE*, float* , int);

std::vector< Edge* > sedges;
bool fprintfEnabled = false;
const char * outputFile = "M for man0.txt";
int ifArray = 0, ifMinHeap =0;

int main(int, char ** argv)
{
	HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);

	SoSeparator * root = new SoSeparator;
	root->ref();

	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();

	const char * horsePath = "meshes/geodesic/fprint matrix/horse0.off";
	const char * manPath = "meshes/geodesic/fprint matrix/man0.off";
	const char * centaurPath = "meshes/geodesic/timing/centaur.off";
	const char * man2Path = "meshes/geodesic/timing/man.off";
	const char * spherePath = "meshes/geodesic/timing/weirdSphere.off";



	mesh->loadOff(manPath);

	root->addChild(painter->getShapeSep(mesh));


	dijkstra(mesh, 150, 501);

	drawThickEdges(root, mesh, sedges);

	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();

	return 0;
}

void drawThickEdges(SoSeparator* root, Mesh* mesh, std::vector< Edge* > sedges) {

	if (sedges.size() == 0)
		return;
	SoSeparator* thickEdgeSep = new SoSeparator;
	//material
	SoMaterial* ma = new SoMaterial;
	ma->diffuseColor.set1Value(0, 1.0f, 0.0f, 0.0f);
	thickEdgeSep->addChild(ma);
	SoDrawStyle* sty = new SoDrawStyle;	sty->lineWidth = 5.0f;	thickEdgeSep->addChild(sty);

	//shape
	SoIndexedLineSet* ils = new SoIndexedLineSet;
	SoCoordinate3* co = new SoCoordinate3;

	//assumes no edge in sedges is removed
	for (unsigned int se = 0; se < sedges.size(); se++)
	{
		SbVec3f end1 = mesh->verts[sedges[se]->v1i]->coords + SbVec3f(0, 0.0f, 0.0f),
			end2 = mesh->verts[sedges[se]->v2i]->coords + SbVec3f(0, 0.0f, 0.0f);
		co->point.set1Value(2 * se, end1);
		co->point.set1Value(2 * se + 1, end2);
	}

	for (unsigned int ci = 0; ci < sedges.size(); ci++)
	{
		ils->coordIndex.set1Value(3 * ci, 2 * ci);	ils->coordIndex.set1Value(3 * ci + 1, 2 * ci + 1);
		ils->coordIndex.set1Value(3 * ci + 2, -1); //end this edge with -1
	}
	thickEdgeSep->addChild(co);	thickEdgeSep->addChild(ils);
	root->addChild(thickEdgeSep);
}

int minDistance(float* dist,	bool* sptSet, int size)
{
	float min = std::numeric_limits<float>::infinity();
	int min_index;

	for (int v = 0; v < size; v++)
		if (sptSet[v] == false &&
			dist[v] <= min)
			min = dist[v], min_index = v;

	return min_index;
}

void calculateThickEdges(Mesh * mesh, int* parent, int srcVertex, int destVertex) {

	if (parent[destVertex] == -1 || srcVertex == destVertex)
		return;

	calculateThickEdges(mesh, parent, srcVertex, parent[destVertex]);
	auto edgeList = mesh->verts[destVertex]->edgeList;
	for (int edge : edgeList) {
		auto e = mesh->edges.at(edge);
		if ((e->v1i == destVertex && e->v2i == parent[destVertex]) || (e->v2i == destVertex && e->v1i == parent[destVertex])) {
			sedges.push_back(e);
			break;
		}
	}
}

void printDistToFile(FILE * fp, float * dist, int size) {

	for (int j = 0; j < size; j++) {
		fprintf(fp, "%g ", dist[j] );
	}
	fprintf(fp,"\n");

}

float calculateDistance(Vertex* v1, Vertex* v2) {
	float result = sqrt(pow((v1->coords[0] - v2->coords[0]), 2) + pow((v1->coords[1] - v2->coords[1]), 2) + pow((v1->coords[2] - v2->coords[2]), 2));
	return result;
}

void dijkstra(Mesh* mesh, int v1, int v2)
{
	int size = mesh->verts.size();
	int * parent = new int[size];
	std::map<int, std::map<int, float>> myMap;

	// map structure used to create graph
	for (int i = 0; i < size; i++) {

		Vertex * v = mesh->verts[i];
		std::vector<int> vertList = v->vertList;
		std::map<int, float> lengthMap;
		for (int j = 0; j < vertList.size(); j++) {
			lengthMap[vertList[j]] = calculateDistance(v, mesh->verts[vertList[j]]);
		}
		myMap[i] = lengthMap;
	}


	if (fprintfEnabled) {
		bool *sptSet = new bool[size];
		float * dist = new float[size];

		FILE * fp;
		fp = fopen(outputFile, "w");

		for (int src = 0; src < size; src++) {

			for (int i = 0; i < size; i++) {
				dist[i] = INT_MAX;
				sptSet[i] = false;
			}

			dist[src] = 0;
			parent[src] = -1;

			for (int count = 0; count < size - 1; count++)
			{
				int u = minDistance(dist, sptSet, size);
				sptSet[u] = true;

				int mapSize = myMap[u].size();
				auto it = myMap[u].begin();
				while (it != myMap[u].end()) {

					if (!sptSet[it->first] && dist[u] + it->second < dist[it->first])
					{
						parent[it->first] = u;
						dist[it->first] = dist[u] + it->second;
					}

					it++;
				}
			}

			if (src == v1) {
				calculateThickEdges(mesh, parent, v1, v2);
			}

			printDistToFile(fp, dist, size);

		}

		fclose(fp);

		delete[] dist;
		delete[] sptSet;
	}
	else {
		parent[v1] = -1;

		bool *sptSet = new bool[size];
		for (int i = 0; i < size; i++) {
			sptSet[i] = false;
		}

		if (ifArray) {
			float * dist = new float[size];
			
			for (int i = 0; i < size; i++) {
				dist[i] = std::numeric_limits<float>::infinity();
			}

			dist[v1] = 0;

			for (int count = 0; count < size - 1; count++)
			{
				int u = minDistance(dist, sptSet, size);
				sptSet[u] = true;

				auto it = myMap[u].begin();
				while (it != myMap[u].end()) {

					if (!sptSet[it->first] && dist[u] + it->second < dist[it->first])
					{
						parent[it->first] = u;
						dist[it->first] = dist[u] + it->second;
					}

					it++;
				}
			}

			calculateThickEdges(mesh, parent, v1, v2);
			printf("%g \n", dist[v2]);

			delete[] dist;
		} 
		else if (ifMinHeap) {
			typedef std::pair<float, int> myPair;
			std::priority_queue<myPair, std::vector<myPair>, std::greater<myPair> > pq;
			std::vector<float> dist(size, std::numeric_limits<float>::infinity());
			std::pair<float, int> p = std::make_pair(0, v1);
			pq.push(p);
			dist[v1] = 0;

			while (!pq.empty())
			{
				int u = pq.top().second;
				pq.pop();
				sptSet[u] = true;

				auto it = myMap[u].begin();
				while (it != myMap[u].end()) {
					
					if (!sptSet[it->first] && dist[u] + it->second < dist[it->first])
					{						
						parent[it->first] = u;
						dist[it->first] = dist[u] + it->second;
						pq.push(std::make_pair(dist[it->first], it->first));
					}

					it++;
				}
			}

			calculateThickEdges(mesh, parent, v1, v2);
			printf("%g \n", dist[v2]);
		}
		else { // fibonacci heap
			typedef std::pair<float, int> myPair;
			FibonacciHeap<myPair> h;
			std::vector<float> dist(size, std::numeric_limits<float>::infinity());
			h.insert(std::make_pair(0, v1));
			dist[v1] = 0;

			while (!h.isEmpty())
			{
				int u = h.removeMinimum().second;
				sptSet[u] = true;

				auto it = myMap[u].begin();
				while (it != myMap[u].end()) {

					if (!sptSet[it->first] &&  dist[u] + it->second < dist[it->first])
					{
						parent[it->first] = u;
						dist[it->first] = dist[u] + it->second;
						h.insert(std::make_pair(dist[it->first], it->first));
					}

					it++;
				}
			}

			calculateThickEdges(mesh, parent, v1, v2);
			printf("%g \n", dist[v2]);
		}

		delete[] sptSet;
	}
	
	delete[] parent;
	
}
