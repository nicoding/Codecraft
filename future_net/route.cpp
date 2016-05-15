#include "route.h"
#include "lib_record.h"
//////////////////////////////////全局变量声明###############################################
unsigned short edgeArray[MAX_EDGE][3];	//边集数组s,t,w
unsigned short nodeEdge[MAX_NODE_NUM][MAX_OUT_DEGREE + 1];	//第一列用来放出度边的条数
unsigned short nodeNode[MAX_NODE_NUM][MAX_OUT_DEGREE + 1];	//第一列用来放出度边的条数
unsigned short weightNodeNode[MAX_NODE_NUM][MAX_NODE_NUM];	//存放节点间边的权重
unsigned short nodetoedge[MAX_NODE_NUM][MAX_NODE_NUM];

bool inDemand[MAX_NODE_NUM] = { false };   //是否在vp中
vector<int> DemandVector;    //vp中的节点编号

unsigned short totalNodeNum = 0;	//总结点数
unsigned short inDemandNum = 0;	//总vp中的点数
int RealVertexNum, RealEdgeNum;  //考虑节点不连续或两点之间多条路径
int *WeightsofEdge;      //边的权重，可更新

int source, target;     //src 和 tar
int newsource, newtarget; // 重新编号的 src 和 tar
int ExtendCap;   //设置代码可扩展性，用于bfs
int JinnodeSize;  //设置bfs搜索的path数
unsigned long timestart;//计时函数

int VeryImportantParameter; //记录一个很重要的参数，用来调整

int dijiVpToVp[MAX_NODE_NUM][MAX_NODE_NUM];
unsigned short vpIndex[MAX_DEMAND_NUM];

unsigned short dchBegin[MAX_TOTAL_WEIGHT];
unsigned short dchEnd[MAX_TOTAL_WEIGHT];
unsigned short dchLength[MAX_TOTAL_WEIGHT];

unsigned short MAX_DCH_QUEUE_SIZE = 300;
unsigned short MAX_DCH_ARRAY_SIZE = 300;
unsigned short MAX_VP_SON_NUM = MAX_DEMAND_NUM;

unsigned short bigGraph=560;
unsigned short bigInDemad=25;


map<int, int> RemapVertex; 
map<int, int> RemapEdge; 

bool vphasdone = false;
bool ifItrate = false;

struct EdgeNode{
	int no;
	EdgeNode *nedge;
};

struct VexHuawei{
	EdgeNode *rsting;
};

VexHuawei EdgeNodeAdj[MAX_NODE_NUM+1];

struct BFSeModel {
    int disorder;
    int inorder;
    int froms;
    int tos;
    int weight;
};
BFSeModel TDedgeQueue[MAX_EDGE]; 
struct BFSvModel {
    int disorder;
    int inorder;
    bool isAvp;
    vector<int> nextEdge; 
    BFSvModel(){
        nextEdge.reserve(8);}
};
BFSvModel TDvertexQueue[MAX_NODE_NUM]; 


struct Node{
	unsigned short no;
	Node * father;
	bool * visit;
	Node * sons[MAX_DEMAND_NUM];
	unsigned short *path;
	unsigned short *sonsInDemand; //记录访问过的vp中的节点
	int *dijiDist;
	int f, g, h;
	unsigned short *dijiSonToFather;
	unsigned short level;
	bool hasBeenDiji;
	short notVisitedVpNum;
	unsigned short notVisitedVps[MAX_DEMAND_NUM];
	Node(unsigned short myNo , Node *hisFather = NULL){
		visit = new bool[totalNodeNum];
		if (hisFather != NULL){
			level = hisFather->level + 1;
			memcpy(visit, hisFather->visit, totalNodeNum*sizeof(bool));
			notVisitedVpNum = hisFather->notVisitedVpNum - 1;
			memcpy(notVisitedVps, hisFather->notVisitedVps, sizeof(unsigned short)*MAX_DEMAND_NUM);
			for (int i = 0; i < notVisitedVpNum; i++){
				if (notVisitedVps[i] == myNo){
					notVisitedVps[i] = notVisitedVps[notVisitedVpNum];
					break;
				}
			}
		}
		else{
			level = 0;
			memset(visit, false, totalNodeNum*sizeof(bool));
			notVisitedVpNum = inDemandNum;
			memcpy(notVisitedVps, vpIndex, sizeof(unsigned short)*MAX_DEMAND_NUM);
		}
		no = myNo;
		visit[no] = true;
		father = hisFather;
		path = new unsigned short[totalNodeNum];
		path[0] = 0;	//第一个位置，记录path中元素的个数
		f = 0;
		g = 0;
		h = 0;
		sonsInDemand = new unsigned short[inDemandNum + 1];
		sonsInDemand[0] = 0;	//第一个位置，记录找到的sonInDemand的个数
		dijiDist = new int[totalNodeNum];
		memset(dijiDist, 128, totalNodeNum*sizeof(int));
		dijiDist[no] = 0;
		dijiSonToFather = new unsigned short[totalNodeNum];
		hasBeenDiji = false;
	}
	void delNode(){}
};
Node* dchArray[MAX_TOTAL_WEIGHT][300];
multimap<int, Node*> dchQueue;

class NodePathInfo {
public:
    int UpperBound;
    int LowerBound;
    std::vector<int> pathviaNode;

    NodePathInfo(int i = 10)
    {
        UpperBound = INT_MAX;
        LowerBound = INT_MAX;
        pathviaNode.reserve(i);
    }
};
NodePathInfo *PathsInformation;

struct Jinnode
{
    int lastPath;
    int Bian;
    int Dian;
    int countSumWeight;
    int pathnodeSum;
    std::bitset<MAX_NODE_NUM> HaveVisited;
};
Jinnode *ZHSNodeQueue;


///////////////////////////////全局变量声明结束//////////////////////////////////////////
//////////////////////////////用于测试的打印函数///////////////////////////////////////
fstream filetest;
void printBoolArray(bool *arr) {
    for (size_t i = 0; i < totalNodeNum; i++) {
        filetest << arr[i] << " ";
    }
    filetest << endl;
    for (size_t i = 0; i < totalNodeNum; i++) {
        if(arr[i])
            filetest << i<<'\t';
    }
    filetest << endl;
}

void printNodeArr(Node **arr, unsigned short totalNum) {
    for (size_t i = 0; i < totalNum; i++){
        filetest << arr[i]->no << " ";
    }
    
    filetest << endl;
}

void printUnsignedShortArr(unsigned short *arr, unsigned short totalNum) {
    for (size_t i = 0; i < totalNum; i++) {
        filetest << arr[i] << " ";
    }
    filetest << endl;
}

void printShort(int *arr, unsigned short totalNum) {
    for (size_t i = 0; i < totalNum; i++) {
        filetest << arr[i] << " ";
    }
    filetest << endl;
}

void printNode(Node * node){
    if (node == NULL)
    {
        filetest << "Node is NUll" << endl;
        filetest << "=================================================================" << endl;
        return;
    }

    filetest <<"Node Num: " << node->no;
    if (node->father != NULL) {
        filetest << "\nFather Num: " << node->father->no << endl;
    }
    else{
        filetest << endl;
    }
    filetest << "Node Visit: ";
    printBoolArray(node->visit);

    filetest << "Node Sons: ";
    printNodeArr(node->sons, node->sonsInDemand[0]);
    
    filetest << "Path: ";
    printUnsignedShortArr(node->path, node->path[0]+1);
    
    filetest << "Sons in Demand: ";
    printUnsignedShortArr(node->sonsInDemand, node->sonsInDemand[0]+1);
    
    filetest << "DijiDist: ";
    printShort(node->dijiDist, totalNodeNum);
    
    filetest << "f = " << node->f << " g = " << node->g << " h = " << node->h << endl;
    
    filetest << "Diji Son to Father: ";
    printUnsignedShortArr(node->dijiSonToFather, totalNodeNum);
    
    filetest << "Node Level: " << node->level << "\nHas been Diji: " << node->hasBeenDiji << endl;
    filetest << "=======================================================================" << endl;
}
//////////////////////////////测试打印结束/////////////////////////////////////////////////////
//////////////////////////////一些保留函数，将来可能有用///////////////////////////////////

unsigned short GetSonNodeEdgeWeight(unsigned short father, unsigned short son){
    return weightNodeNode[father][son];
}


unsigned short GetAllWeight(unsigned short result[]){
    int w = 0;
    for (int i = result[0]; i>1; i--)
        w += GetSonNodeEdgeWeight(result[i-1], result[i]);
    return w;
}


void OutputResult(unsigned short result[], int tmpg){// windows 下面的打印函数
    ofstream f;
    f.open("C://Users//vip//Desktop//ooooxxxx.txt", ios::app);
    //cout<<"---------------start output-----------------"<<endl;
    f << tmpg << endl;
    for (int i = result[0]; i>0; i--)
        f << result[i] << '\t';
    f << endl;
    f.close();
}


/*void initNextT_zhs(){
    bool hasBeenVisited[MAX_NODE_NUM] = { false };
    hasBeenVisited[target] = true;
    int j;
    //把图翻转一下，从t开始往回走，检查weightNodeNode即可，扫一列
    queue<unsigned short> calArea;
    calArea.push(target);
    while (!calArea.empty()){
        //对待算区的front的子节点进行扩充
        for (j = 0; j < MAX_NODE_NUM; j++){
            if (weightNodeNode[j][calArea.front()]!=0 && !hasBeenVisited[j]){
                hasBeenVisited[j] = true;
                if (inDemand[j]){
                    nextT[j] = true;
                    nextTNum++;
                }
                else
                    calArea.push(j);
            }
        }
        calArea.pop();
    }
}


void initlevelVStart_zhs(){
    bool tmpNextT[MAX_NODE_NUM] = { false };
    vector<vector<unsigned short> > tmpNextVp(MAX_NODE_NUM);//存储当前vp的下一层，有哪些vp
    queue<unsigned short> calArea;
    calArea.push(target);
    bool hasBeenVisited[MAX_NODE_NUM] = { false };
    unsigned short level = inDemandNum;
    while (!calArea.empty()){
        for (int j = 0; j < MAX_NODE_NUM; j++){
            if (weightNodeNode[j][calArea.front()] != 0 && !hasBeenVisited[j]){
                hasBeenVisited[j] = true;
                if (inDemand[j]){
                    tmpNextT[j] = true;
                    nextTNum++;
                }
                else
                    calArea.push(j);
            }
        }
        calArea.pop();
    }
    unsigned short vpIndex = 0;
    for(int i=0; i<inDemandNum; i++){
        memset(hasBeenVisited, false, MAX_NODE_NUM*sizeof(bool));
        while (!inDemand[vpIndex]){
            vpIndex++;
            if (vpIndex>MAX_NODE_NUM)
                cout << "error!" << endl;
        }
        calArea.push(vpIndex);
        while (!calArea.empty()){
            for (int j = 0; j < MAX_NODE_NUM; j++){
                if (weightNodeNode[j][calArea.front()] != 0 && !hasBeenVisited[j]){
                    hasBeenVisited[j] = true;
                    if (inDemand[j]){
                        tmpNextVp[vpIndex].push_back(j);
                    }
                    else
                        calArea.push(j);
                }
            }
            calArea.pop();
        }
        vpIndex++;
    }
    level = inDemandNum;
    memcpy(levelVStart[level], tmpNextT, MAX_NODE_NUM*sizeof(bool));
    for (int i = 0; i < MAX_NODE_NUM; i++){
        if (tmpNextT[i]){
            tmpNextVp[target].push_back(i);
            levelVStart[level - 1][i] = true;
        }
    }
    unsigned short curLevelNo = nextTNum;
    while (--level){
        if (curLevelNo == inDemandNum){
            memcpy(levelVStart[level - 1], levelVStart[level], MAX_NODE_NUM*sizeof(bool));
            continue;
        }
        curLevelNo = 0;
        //进入这个循环，level至少为1
        for (int i = 0; i < MAX_NODE_NUM; i++){
            if (levelVStart[level][i]){
                for (int j = 0; j < tmpNextVp[i].size(); j++){
                    if (!levelVStart[level - 1][tmpNextVp[i][j]]){
                        levelVStart[level - 1][tmpNextVp[i][j]] = true;
                        curLevelNo++;
                    }
                }
            }
        }
    }
}


void initFloyed() {
    for (int i = 0; i < totalNodeNum; i++){
        for (int j = 0; j < totalNodeNum; j++){
            if (i == j) {
                weightFloyed[i][j] = 0;
            
            }
            else {
                if (weightNodeNode[i][j] != 0)
                    weightFloyed[i][j] = weightNodeNode[i][j];
            }
        }
    }
    
    for (int k = 0; k < totalNodeNum; k++){
        for (int i = 0; i < totalNodeNum; i++){
            for (int j = 0; j < totalNodeNum; j++){
                if (weightFloyed[i][k] + weightFloyed[k][j] < weightFloyed[i][j]) {
                    weightFloyed[i][j] = weightFloyed[i][k] + weightFloyed[k][j];

                    //floyedPath[i][j] = floyedPath[k][j];
                }
            }
        }
    }
}

void initReverseDijkstra(){
    // check whether the node get into set
    bool visitedFromTarget[MAX_NODE_NUM];

    //initialize the distToTarget
    for (int i = 0; i < totalNodeNum; ++i) {
        if (weightNodeNode[i][target] > 0 && i != target) {
            distToTarget[i] = weightNodeNode[i][target];
        }
        else {
            distToTarget[i] = MAX_DIST;
        }
        visitedFromTarget[i] = false;
    }
    distToTarget[target] = 0;
    visitedFromTarget[target] = true;

    for (int i = 1; i < totalNodeNum; ++i) {
        int minDist = MAX_DIST;
        int minDistNode;

        // look for the min weight node
        for (int j = 0; j < totalNodeNum; ++j) {
            if (visitedFromTarget[j] == false && distToTarget[j] < minDist) {
                minDist = distToTarget[j];
                minDistNode = j;
            }
        }

        visitedFromTarget[minDistNode] = true;

        //update distToTarget[]
        for (int k = 0; k < totalNodeNum; ++k) {
            if (visitedFromTarget[k] == false && weightNodeNode[k][minDistNode] > 0
                && minDist + weightNodeNode[k][minDistNode] < distToTarget[k])
            {
                distToTarget[k] = minDist + weightNodeNode[k][minDistNode];
            }
        }
    }
}

int randint(int l, int u)
{
    int temp;
    temp = floor(l + (1.0*rand() / RAND_MAX)*(u - l));
    return temp;
}*/

///////////////////////////////初始化数据区/////////////////////////////////////////////

void SplitCharStar(char *str, unsigned short re[4]){
	int a, b, c, d;
	sscanf(str, "%d,%d,%d,%d", &a, &b, &c, &d);
	re[0]=a;
	re[1]=b;
	re[2]=c;
	re[3]=d;
}


void SplitCharStarToSet(char *c){
	const char delimiter = '|';
	string s(c);

	int midval;
	string::size_type end;
	string::size_type start = 0;
	end = s.find(delimiter);
	while (end != s.npos){
		midval = atoi(s.substr(start, end - start).c_str());
		inDemand[midval] = true;
		vpIndex[inDemandNum] = midval;
		inDemandNum++;
		TDvertexQueue[RemapVertex[midval]].isAvp = 1;
        DemandVector.push_back(RemapVertex[midval]);
		start = end + 1;
		end = s.find(delimiter, start);

	}
	end = s.length();
	midval = atoi(s.substr(start, end - start).c_str());
	inDemand[midval] = true;
	vpIndex[inDemandNum] = midval;
	inDemandNum++;
	TDvertexQueue[RemapVertex[midval]].isAvp = 1;
    DemandVector.push_back(RemapVertex[midval]);

}


void zhsinit_graph(char *graph[5000], int edge_num, char *condition){
	RealVertexNum = 0;
    RealEdgeNum = 0;


	unsigned short edgeInfo[4];
	int i;
	EdgeNode *q;
	int  para2, para3;

	for (i = 0;i<MAX_NODE_NUM+1;i++) EdgeNodeAdj[i].rsting = NULL;
	//routeDataInit();
	bool count[MAX_NODE_NUM] = { false };
	for (i = 0; i<edge_num; i++){
		SplitCharStar(graph[i], edgeInfo);
		//cout << "edgeInfo:\t" << edgeInfo[0] << '\t' << edgeInfo[1] << '\t' << edgeInfo[2] << '\t' << edgeInfo[3] << endl;
		nodetoedge[edgeInfo[1]][edgeInfo[2]] = edgeInfo[0];
		edgeArray[edgeInfo[0]][0] = edgeInfo[1];
		edgeArray[edgeInfo[0]][1] = edgeInfo[2];
		edgeArray[edgeInfo[0]][2] = edgeInfo[3];
		if (count[edgeInfo[1]] == false){
			totalNodeNum++;
			count[edgeInfo[1]] = true;
		}
		if (count[edgeInfo[2]] == false){
			totalNodeNum++;
			count[edgeInfo[2]] = true;
		}
		nodeEdge[edgeInfo[1]][0]++;
		nodeEdge[edgeInfo[1]][nodeEdge[edgeInfo[1]][0]] = edgeInfo[0];
		nodeNode[edgeInfo[1]][0]++;
		nodeNode[edgeInfo[1]][nodeNode[edgeInfo[1]][0]] = edgeInfo[2];
		if (weightNodeNode[edgeInfo[1]][edgeInfo[2]] == 0)
			weightNodeNode[edgeInfo[1]][edgeInfo[2]] = edgeInfo[3];
		else
			weightNodeNode[edgeInfo[1]][edgeInfo[2]] = min(weightNodeNode[edgeInfo[1]][edgeInfo[2]], edgeInfo[3]);

		if(RemapVertex.count(int(edgeInfo[1]))) 
            para2 = RemapVertex[int(edgeInfo[1])];
        else{            
            RemapVertex.insert(pair<int, int>(int(edgeInfo[1]), RealVertexNum));
		    TDvertexQueue[RealVertexNum].disorder = int(edgeInfo[1]);
		    TDvertexQueue[RealVertexNum].inorder = RealVertexNum;
		    
		    para2 = RealVertexNum;
		    RealVertexNum++;
        }

        if(RemapVertex.count(int(edgeInfo[2])))
            para3 = RemapVertex[int(edgeInfo[2])];
        else{
            RemapVertex.insert(pair<int, int>(int(edgeInfo[2]), RealVertexNum));
		    TDvertexQueue[RealVertexNum].disorder = int(edgeInfo[2]);
		    TDvertexQueue[RealVertexNum].inorder = RealVertexNum;		    
		    para3 = RealVertexNum;
		    RealVertexNum++;
        }

        int tmpID;
        bool tmpflag = false;
        if(para2 >= RealVertexNum || para3 >= RealVertexNum)
	        tmpflag = false;
	    else{
		    vector<int> *nextEdges = &(TDvertexQueue[para2].nextEdge);
		    for(vector<int>::iterator it=nextEdges->begin(); it!=nextEdges->end(); it++) {
		        if(TDedgeQueue[*it].tos == para3) {
		            tmpID = *it;
		            tmpflag = true;
		            break;
		        }
		    }
		}

	    if(tmpflag) {
	        if(int(edgeInfo[3]) < TDedgeQueue[tmpID].weight) {
	            RemapEdge.erase(TDedgeQueue[tmpID].disorder);
	            RemapEdge.insert(pair<int, int>(int(edgeInfo[0]), tmpID));
	            TDedgeQueue[tmpID].disorder = int(edgeInfo[0]);
	            TDedgeQueue[tmpID].weight = int(edgeInfo[3]);
	        }
	    } else {
	        RemapEdge.insert(pair<int, int>(int(edgeInfo[0]), RealEdgeNum));
	        TDedgeQueue[RealEdgeNum].disorder = int(edgeInfo[0]);
	        TDedgeQueue[RealEdgeNum].inorder = RealEdgeNum;
	        TDedgeQueue[RealEdgeNum].weight = int(edgeInfo[3]);
	        TDedgeQueue[RealEdgeNum].froms = para2;
	        TDedgeQueue[RealEdgeNum].tos = para3;
	        TDvertexQueue[para2].nextEdge.push_back(RealEdgeNum);
	        RealEdgeNum++;
	    }

		q = new EdgeNode;
        q->no = edgeInfo[2];
        q->nedge = EdgeNodeAdj[edgeInfo[1]].rsting;
        EdgeNodeAdj[edgeInfo[1]].rsting = q;
	}


	//unsigned short will be wrong in sscanf
	char vDemandChar[200];
	sscanf(condition, "%d,%d,%s", &source, &target, vDemandChar);
	newsource = RemapVertex[source];
	newtarget = RemapVertex[target];
	SplitCharStarToSet(vDemandChar);
}
///////////////////////////////初始化数据区结束////////////////////////////////////
///////////////////////////////核心算法区开始：三种算法//////////////////////////////////////////
//暴力求解dfs：适合小规模
void search_route_dfs(){
	int stackfront;
	int vextmp;
	int backtrack = 1;
	int i;
	int nodeindemand;
	int pathWeightsum;

	int visit[MAX_NODE_NUM+1];
	int dfsstack[MAX_NODE_NUM+1];
	
	struct EdgeNode *tmp = NULL;
	short MinimumPathWeight = MAX_SHORT_NUM;
	short result_record[MAX_NODE_NUM+1] = { 0 };

	for (i = 0;i<MAX_NODE_NUM+1;i++) 
		visit[i] = 0;

	vextmp = source;
	visit[source] = 1;
	stackfront = 1;
	dfsstack[stackfront] = vextmp;

	do {
	if (backtrack == 1) {
		tmp = EdgeNodeAdj[vextmp].rsting;
		backtrack = 0;
	}
	else 
		tmp = tmp->nedge;
	if (tmp){
		if (visit[tmp->no] == 0){
			visit[tmp->no] = 1;
			stackfront++;
			dfsstack[stackfront] = tmp->no;
			if (tmp->no == target){

				if (stackfront - 2 >= inDemandNum) {
					nodeindemand = 0;
					for (i = 1;i <= stackfront;i++) {
						if (inDemand[dfsstack[i]] != false)
							nodeindemand++;							
					}
					if (inDemandNum == nodeindemand) {
						pathWeightsum = 0;

						for (i = 1;i < stackfront;i++) 
							pathWeightsum += weightNodeNode[dfsstack[i]][dfsstack[i + 1]];

						if (pathWeightsum < MinimumPathWeight) {

							MinimumPathWeight = pathWeightsum;

							for (i = 1;i < stackfront;i++) 
								result_record[i] = nodetoedge[dfsstack[i]][dfsstack[i + 1]];

							result_record[stackfront] = A_NEGATIVE_NUM;
						}
					}
				}

				visit[target] = 0;
				stackfront--;
				vextmp = dfsstack[stackfront];
				backtrack = 0;
			}
			else {
				vextmp = dfsstack[stackfront];
				backtrack = 1;
			}
		} 
	}
	else {
		visit[dfsstack[stackfront--]] = 0; 
		if (stackfront){

			tmp = EdgeNodeAdj[dfsstack[stackfront]].rsting;
			while (tmp->no != vextmp) 
				tmp = tmp->nedge;

			vextmp = dfsstack[stackfront];
			backtrack = 0;
		}
	}
	} while (stackfront);

	if (MinimumPathWeight != MAX_SHORT_NUM) {
		i = 1;
		while (result_record[i] != A_NEGATIVE_NUM) {
			record_result(result_record[i]);
			i++;
		}
	}	
}


////////////////////////
//广度搜索bfs，适合规模居中
void search_route_bfs()
{
	//更新graph的一些信息
    WeightsofEdge = new int[RealEdgeNum];
    for (int i = 0; i < RealEdgeNum; ++i)
    {
        if (TDedgeQueue[i].tos == newsource || TDedgeQueue[i].tos == newtarget)
        {
            WeightsofEdge[i] = TDedgeQueue[i].weight;
        }
        else if (TDvertexQueue[TDedgeQueue[i].tos].isAvp)
            WeightsofEdge[i] = TDedgeQueue[i].weight - PWW;
        else
            WeightsofEdge[i] = TDedgeQueue[i].weight;
    }

    PathsInformation = new NodePathInfo[RealVertexNum];  


    if(RealVertexNum <= 80){
        VeryImportantParameter = 50;
    } else if(RealVertexNum <= 100){
        VeryImportantParameter = RealVertexNum*6;  
    } else if(RealVertexNum <= 200) {
        VeryImportantParameter = RealVertexNum*4;
    }else if(RealVertexNum <= 300) {
        VeryImportantParameter = 250;  
    } else {
        VeryImportantParameter = 200;
    }
    //算法程序开始
    ExtendCap = MAX_BFS_SIZE;
    ZHSNodeQueue = new Jinnode[ExtendCap];
    JinnodeSize = 0;

    NodePathInfo *newPathsInformation = new NodePathInfo[RealVertexNum];

    //初始化
    PathsInformation[newsource].UpperBound = 0;
    PathsInformation[newsource].LowerBound = 0;

    if (JinnodeSize == ExtendCap)
    {
        Jinnode *newList = new Jinnode[2 * ExtendCap];
        for (int i = 0; i < JinnodeSize; ++i)
        {
            newList[i] = ZHSNodeQueue[i];
        }
        delete[] ZHSNodeQueue;
        ZHSNodeQueue = newList;
        ExtendCap *= 2;
    }
    ZHSNodeQueue[JinnodeSize].Dian = newsource;
    ZHSNodeQueue[JinnodeSize].Bian = -1;
    ZHSNodeQueue[JinnodeSize].lastPath = -1;
    
    ZHSNodeQueue[JinnodeSize].countSumWeight = 0;
    ZHSNodeQueue[JinnodeSize].pathnodeSum = 0;
    ZHSNodeQueue[JinnodeSize].HaveVisited.set(newsource);

    ++JinnodeSize;

    PathsInformation[newsource].pathviaNode.push_back(JinnodeSize - 1);

    newPathsInformation[newsource].UpperBound = 0;
    newPathsInformation[newsource].LowerBound = 0;
    newPathsInformation[newsource].pathviaNode.push_back(PathsInformation[newsource].pathviaNode[0]);

    bool beenmodified = false;
    for (int i = 0; i < RealVertexNum - 1; ++i)
    {
        beenmodified = false;
        for (int j = 0; j < RealEdgeNum; ++j)
        {
            int tmpedgesou = TDedgeQueue[j].froms;
            int tmpedgedes = TDedgeQueue[j].tos;
            if (PathsInformation[tmpedgesou].UpperBound == INT_MAX)
            {
                continue;
            }
            else if (PathsInformation[tmpedgesou].LowerBound + WeightsofEdge[j]< PathsInformation[tmpedgedes].UpperBound)
            {
                int tmpdstID = TDedgeQueue[j].tos;

			    NodePathInfo *src = &PathsInformation[tmpedgesou];
			    NodePathInfo *dst = &newPathsInformation[tmpdstID];

				unsigned int mm = 0;
			    unsigned int nn = 0;
			    int pathcalc = 0;
				
			    static vector<int> constantPath;

			    bool beenchanged = false;
			    while (pathcalc < VeryImportantParameter)
			    {
			        int newPathID;

			        if (mm == src->pathviaNode.size() && nn == dst->pathviaNode.size())
			        {
			            break;
			        }
			        else if (mm == src->pathviaNode.size())
			        {
			            constantPath.push_back(dst->pathviaNode[nn]);
			            ++nn;
			            ++pathcalc;
			        }
			        else if (ZHSNodeQueue[src->pathviaNode[mm]].HaveVisited[tmpdstID])
			        {
			            ++mm;
			            continue;
			        }
			        else if (nn == dst->pathviaNode.size())
			        {
						if (JinnodeSize == ExtendCap)
						{
						    Jinnode *newList = new Jinnode[2 * ExtendCap];
						    for (int mm = 0; mm < JinnodeSize; ++mm)
						    {
						        newList[mm] = ZHSNodeQueue[mm];
						    }
						    delete[] ZHSNodeQueue;
						    ZHSNodeQueue = newList;
						    ExtendCap *= 2;
						}
						ZHSNodeQueue[JinnodeSize].Bian = j;
						ZHSNodeQueue[JinnodeSize].Dian = tmpdstID;
						ZHSNodeQueue[JinnodeSize].lastPath = src->pathviaNode[mm];

						ZHSNodeQueue[JinnodeSize].countSumWeight = ZHSNodeQueue[src->pathviaNode[mm]].countSumWeight + WeightsofEdge[j];
						ZHSNodeQueue[JinnodeSize].pathnodeSum = ZHSNodeQueue[src->pathviaNode[mm]].pathnodeSum + 1;

						ZHSNodeQueue[JinnodeSize].HaveVisited = ZHSNodeQueue[src->pathviaNode[mm]].HaveVisited;
						ZHSNodeQueue[JinnodeSize].HaveVisited.set(tmpdstID);

						++JinnodeSize;
						newPathID = JinnodeSize - 1;

			            constantPath.push_back(newPathID);
			            ++mm;
			            ++pathcalc;
			            beenchanged = true;
			        }
			        else if (ZHSNodeQueue[dst->pathviaNode[nn]].lastPath == src->pathviaNode[mm])
			        {
			            ++mm;
			            continue;
			        }
			        else if (ZHSNodeQueue[src->pathviaNode[mm]].countSumWeight + WeightsofEdge[j] < ZHSNodeQueue[dst->pathviaNode[nn]].countSumWeight)
			        {
			            if (JinnodeSize == ExtendCap)
						{
						    Jinnode *newList = new Jinnode[2 * ExtendCap];
						    for (int mm = 0; mm < JinnodeSize; ++mm)
						    {
						        newList[mm] = ZHSNodeQueue[mm];
						    }
						    delete[] ZHSNodeQueue;
						    ZHSNodeQueue = newList;
						    ExtendCap *= 2;
						}
						ZHSNodeQueue[JinnodeSize].Bian = j;
						ZHSNodeQueue[JinnodeSize].Dian = tmpdstID;

						ZHSNodeQueue[JinnodeSize].lastPath = src->pathviaNode[mm];
						ZHSNodeQueue[JinnodeSize].countSumWeight = ZHSNodeQueue[src->pathviaNode[mm]].countSumWeight + WeightsofEdge[j];
						
                        ZHSNodeQueue[JinnodeSize].pathnodeSum = ZHSNodeQueue[src->pathviaNode[mm]].pathnodeSum + 1;
						ZHSNodeQueue[JinnodeSize].HaveVisited = ZHSNodeQueue[src->pathviaNode[mm]].HaveVisited;
						ZHSNodeQueue[JinnodeSize].HaveVisited.set(tmpdstID);

						++JinnodeSize;
						newPathID = JinnodeSize - 1;

			            constantPath.push_back(newPathID);
			            ++mm;
			            ++pathcalc;
			            beenchanged = true;
			        }
			        else
			        {
			            constantPath.push_back(dst->pathviaNode[nn]);
			            ++nn;
			            ++pathcalc;
			        }
			    }
			    if (beenchanged)
			    {
			        dst->pathviaNode = constantPath;
			    }
			    dst->LowerBound = ZHSNodeQueue[constantPath[0]].countSumWeight;
			    dst->UpperBound = ZHSNodeQueue[constantPath.back()].countSumWeight;

			    constantPath.clear();
			    
                beenmodified = true;
            }
        }
        if (!beenmodified)
        {
            break;
        }
        for (int j = 0; j < RealVertexNum; ++j)
        {
            PathsInformation[j] = newPathsInformation[j];
        }
    }

    delete[]newPathsInformation;

    if (PathsInformation[newtarget].LowerBound > -PWW * inDemandNum && PathsInformation[newtarget].LowerBound < -PWW * (inDemandNum - 1))
    {
        int *tmppath = new int[MAX_NODE_NUM];
	    int tmplength;

	    int currentNode = PathsInformation[newtarget].pathviaNode[0];
	    tmplength = ZHSNodeQueue[currentNode].pathnodeSum;
	    for (int i = tmplength - 1; i >= 0; --i)
	    {
	        tmppath[i] = ZHSNodeQueue[currentNode].Bian;
	        currentNode = ZHSNodeQueue[currentNode].lastPath;
	    }

	    for (int i = 0; i < tmplength; ++i)
	    {
	        record_result(tmppath[i]);
	    }
	    delete[] tmppath;
    }
    else
    {
        return;
    }
}

//////////////////////////////////
//分支限界，启发式搜索，适合大规模

int randint(int l, int u)
{
    int temp;
    temp = floor(l + (1.0*rand() / RAND_MAX)*(u - l));
    return temp;
}

unsigned long GetTickCount()
{
  struct timeval tv;
  if( gettimeofday(&tv, NULL) != 0 )
    return 0;
 
  return (tv.tv_sec * 1000) + (tv.tv_usec / 1000);
}


void setMyHFunction_0(Node * currNode) {
	currNode->h = 0;
	return;
}

void setMyHFunction(Node * currNode,int curNodeDCH) {
	if (currNode->no == target){
		currNode->h = 0;
		return;
	}
	if (currNode->notVisitedVpNum == 0){
		currNode->h = dijiVpToVp[currNode->no][target];
		return;
	}

	int tmpH = 0;
	int tmp_min=MAX_DIST;
	for (int i = 0; i < currNode->notVisitedVpNum; i++){
		tmp_min = min(tmp_min, dijiVpToVp[currNode->no][currNode->notVisitedVps[i]]);
	}
	if (tmp_min == MAX_DIST){
		currNode->h = MAX_DIST;
		return;
	}
	tmpH = tmp_min;
	unsigned short noPath = MAX_NODE_NUM+1;
	if (currNode->notVisitedVpNum>1){
		int tmp_max = 0;
		for (int i = 0; i < currNode->notVisitedVpNum; i++){
			tmp_min = MAX_DIST;
			for (int j = 0; j < currNode->notVisitedVpNum; j++){
				if (i != j)
					tmp_min = min(tmp_min, dijiVpToVp[currNode->notVisitedVps[i]][currNode->notVisitedVps[j]]);
			}
			if (tmp_min == MAX_DIST){
				if (noPath == MAX_NODE_NUM + 1){
					noPath = currNode->notVisitedVps[i];
				}
				else{
					currNode->h = MAX_DIST;
					return;
				}
			}
			tmp_max = max(tmp_max, tmp_min);
			tmpH += tmp_min;
		}
		tmpH -= tmp_max;
	}

	if (noPath == MAX_NODE_NUM + 1){
		tmp_min = MAX_DIST;
		for (int i = 0; i < currNode->notVisitedVpNum; i++){
			tmp_min = min(tmp_min, dijiVpToVp[currNode->notVisitedVps[i]][target]);
		}
		if (tmp_min == MAX_DIST){
			currNode->h = MAX_DIST;
			return;
		}
		tmpH += tmp_min;
		currNode->h = tmpH;
	}
	else{
		if (dijiVpToVp[noPath][target] == MAX_DIST){
			currNode->h = MAX_DIST;
			return;
		}
		else
			currNode->h = dijiVpToVp[noPath][target] + tmpH;
	}
    int myPath = (MAX_TOTAL_WEIGHT/(BESTROAD_NUM+POP_TOP_NUM));
    if(totalNodeNum>=bigGraph&&inDemandNum>=bigInDemad)
        currNode->h=(myPath/(currNode->h+curNodeDCH))?currNode->h:MAX_DIST;
}

void initDijiVpToVp(){
	//dijiVpToVp
	bool visited[MAX_NODE_NUM];
	int dist[MAX_NODE_NUM];
	int bestNode, bestDist;
	int i, j;
	for (i = 0; i < inDemandNum; i++){
		memset(visited, false, sizeof(bool)*MAX_NODE_NUM);
		memset(dist, 1, sizeof(int)*MAX_NODE_NUM);
		//visited[vpIndex[i]] = false;
		dist[vpIndex[i]] = 0;
		bestNode = vpIndex[i];
		while (1){
			bestDist = MAX_DIST;
			for (j = 0; j < MAX_NODE_NUM; j++){
				if (!visited[j]){
					if (dist[j] < bestDist){
						bestNode = j;
						bestDist = dist[j];
					}
				}
			}
			if (bestDist == MAX_DIST)
				break;
			visited[bestNode] = true;
			if (!inDemand[bestNode] || bestNode == vpIndex[i]){
				for (j = 1; j < nodeNode[bestNode][0] + 1; j++){
					if (!visited[nodeNode[bestNode][j]])
						dist[nodeNode[bestNode][j]] = min(dist[nodeNode[bestNode][j]], dist[bestNode] + weightNodeNode[bestNode][nodeNode[bestNode][j]]);
				}
			}
		}
		memcpy(dijiVpToVp[vpIndex[i]], dist,sizeof(int)*MAX_NODE_NUM);
	}
}
// =========================dhcArray opreate by zhs 4_8=================
bool empty_DchArray(unsigned short index){
	return dchLength[index] == 0;
}

bool full_DchArray(unsigned short index){
	return dchLength[index] == MAX_DCH_ARRAY_SIZE;
}

unsigned short findNextMinIndex(unsigned short min_index){
	while (empty_DchArray(min_index))
		min_index++;
	return min_index;
}

int findNextMaxIndex(int max_index){
	while (empty_DchArray(max_index))
		max_index--;
	return max_index;
}

void popBack_DchArray(unsigned short index){
	if (dchEnd[index] != dchBegin[index])
		dchEnd[index] = (dchEnd[index] + MAX_DCH_ARRAY_SIZE - 1) % MAX_DCH_ARRAY_SIZE;
	dchLength[index]--;
}

void popFront_DchArray(unsigned short index){
	if (dchEnd[index] != dchBegin[index])
		dchBegin[index] = (dchBegin[index] + MAX_DCH_ARRAY_SIZE + 1) % MAX_DCH_ARRAY_SIZE;
	dchLength[index]--;
}


unsigned short delCor(unsigned short index){
	unsigned short tmp = dchLength[index];
	dchLength[index] = 0;
	dchBegin[index] = 0;
	dchEnd[index] = 0;
	return tmp;
}

//#define hdti
//#define ti
//#define hd
unsigned short pushBack_DchArray(unsigned short index, Node* node){
	if (full_DchArray(index)){
		//从头删除，从尾巴插入
#ifdef hdti
		dchBegin[index] = (dchBegin[index] + 1) % MAX_DCH_ARRAY_SIZE;
		dchEnd[index] = (dchEnd[index] + 1) % MAX_DCH_ARRAY_SIZE;
		dchArray[index][dchEnd[index]] = node;
		//从尾巴删除，从头插入
#else
		dchBegin[index] = (dchBegin[index] - 1 + MAX_DCH_ARRAY_SIZE) % MAX_DCH_ARRAY_SIZE;
		dchEnd[index] = (dchEnd[index] - 1 + MAX_DCH_ARRAY_SIZE) % MAX_DCH_ARRAY_SIZE;
		dchArray[index][dchBegin[index]] = node;
#endif
		return 0;
	}
	else{
		if (empty_DchArray(index))
			dchArray[index][dchEnd[index]] = node;
		else{
			//从尾插插入
#ifdef ti
			dchEnd[index] = (dchEnd[index] + 1) % MAX_DCH_ARRAY_SIZE;
			dchArray[index][dchEnd[index]] = node;
			//从头插入
#else
			dchBegin[index] = (dchBegin[index] - 1 + MAX_DCH_ARRAY_SIZE) % MAX_DCH_ARRAY_SIZE;
			dchArray[index][dchBegin[index]] = node;
#endif
		}
		dchLength[index]++;
		return 1;
	}
}


void clearDCHArray(){
    memset(dchEnd,0,sizeof(unsigned short)*MAX_TOTAL_WEIGHT);
    memset(dchBegin,0,sizeof(unsigned short)*MAX_TOTAL_WEIGHT);
    memset(dchLength,0,sizeof(unsigned short)*MAX_TOTAL_WEIGHT);
    //menset(dchArray,0,sizeof()
}


void OutputResult_linux(unsigned short result[]){
	for (int i = result[0]; i>1; i--)
		record_result( nodetoedge[result[i]][result[i-1]] );
}


//dijistra终止条件：没加进来的集合为空集
//返回值，0：说明失败，1：说明搜到了vp中的点，2：说明遍历了所有vp中的点且到达t
Node* dijistra_zhs(Node * node){
	unsigned short bestNodeNo = 0;
	int bestNodeDist;
	unsigned short sonNo, newWeight;
	while (1){
		if (!node->hasBeenDiji){
			bestNodeNo = node->no;
			bestNodeDist = 0;
		}
		else{
			bestNodeDist = -MAX_DIST - 1;
			for (int i = 0; i < totalNodeNum; i++){
				if (!node->visit[i] && node->dijiDist[i]<0 && node->dijiDist[i]>bestNodeDist){
					bestNodeNo = i;
					bestNodeDist = node->dijiDist[i];
				}
			}
			if (bestNodeDist == -MAX_DIST - 1)
				//说明没有找到待加入节点
				return NULL;
		}
		bestNodeDist = -node->dijiDist[bestNodeNo];
		node->dijiDist[bestNodeNo] = bestNodeDist;

		if (bestNodeNo == target&&node->level != inDemandNum){
			continue;
		}
		if ((inDemand[bestNodeNo] || bestNodeNo == target) && node->hasBeenDiji){
			//找到了vp中的点或tagert
			//更新father
			node->sonsInDemand[0]++;
			node->sonsInDemand[node->sonsInDemand[0]] = bestNodeNo;
			Node *vpSonNode = new Node(bestNodeNo,node);
			node->sons[node->sonsInDemand[0] - 1] = vpSonNode;
			//设置子节点
			while (node->dijiSonToFather[bestNodeNo] != node->no){	//s和v1都没有压进path，它之间的节点压进去
				bestNodeNo = node->dijiSonToFather[bestNodeNo];
				vpSonNode->path[0]++;
				vpSonNode->path[vpSonNode->path[0]] = bestNodeNo;
				vpSonNode->visit[bestNodeNo] = true;
			}
			vpSonNode->g = node->g + bestNodeDist;
			setMyHFunction(vpSonNode,vpSonNode->g);
			vpSonNode->f = vpSonNode->g + vpSonNode->h;
			//vpSonNode->f = (unsigned short)(vpSonNode->g*1.0*(inDemandNum + 1) / vpSonNode->level);
			return vpSonNode;
		}
		node->hasBeenDiji = true;
		for (int i = 1; i < nodeNode[bestNodeNo][0] + 1; i++){
			sonNo = nodeNode[bestNodeNo][i];
			if (!node->visit[sonNo] && node->dijiDist[sonNo]<0){
				newWeight = bestNodeDist + weightNodeNode[bestNodeNo][sonNo];
				if (abs(node->dijiDist[sonNo])>newWeight){
					//更新路径
					node->dijiDist[sonNo] = -newWeight;
					node->dijiSonToFather[sonNo] = bestNodeNo;
				}
			}
		}
	}
}






int search_route_array(int maxf=MAX_DIST){

    if(!vphasdone){
	   initDijiVpToVp();
       vphasdone = true;
    }

	Node * head = new Node(source);
	head->f = 0;
	dchArray[0][0]=head;
	dchLength[0] = 1;
	unsigned short dch_size = 1;
	unsigned short dch_min_index = 0;
	int dch_max_index = 0;
	Node *findSource,* findResult;
	unsigned short successPath[MAX_NODE_NUM + 1];
	successPath[0] = 0;

	Node * bestNode=NULL;
	int bestF = 2147483647;
	unsigned short successNum = 0;

	unsigned short maxSonVpNum = MAX_VP_SON_NUM;

	while (dch_size != 0){
		
		
		dch_min_index = findNextMinIndex(dch_min_index);

		if (dch_min_index>bestF){
			dch_size = dch_size	-	delCor(dch_min_index);
			break;
		}

#ifdef hd
		findSource = dchArray[dch_min_index][dchBegin[dch_min_index]];
#else
		findSource = dchArray[dch_min_index][dchEnd[dch_min_index]];
#endif


		
		findResult = dijistra_zhs(findSource);
		if (findResult == NULL){
			//Source节点找不出子节点了
#ifdef hd
			popFront_DchArray(dch_min_index);
#else
			popBack_DchArray(dch_min_index);
#endif
			dch_size--;
			if (findSource->sonsInDemand[0] == 0){
				//说明它没有儿子节点，可以删除它
				//TODO
				//是不是只要做这一个操作？
				//dchQueue.erase(mapit);
				delete findSource;
			}
		}
		else if (findResult->no == target){
			//找到完整路径
			
			int tmpf = findResult->g;
			if (bestF>tmpf){
				bestF = tmpf;
				bestNode = findResult;
				successNum++;
			}
			
			//if (successNum < BESTROAD_NUM)
				continue;
			
			while (bestNode->no != source){
				successPath[0]++;
				successPath[successPath[0]] = bestNode->no;
				for (int i = 1; i <= bestNode->path[0]; i++){
					successPath[0]++;
					successPath[successPath[0]] = bestNode->path[i];
				}
				bestNode = bestNode->father;
			}
			successPath[0]++;
			successPath[successPath[0]] = source;
			//OutputResult(successPath, tmpf);
			OutputResult_linux(successPath);
			return MAX_DIST;
		}
		else if(findResult->f<maxf){
			if (findSource->sonsInDemand[0] == maxSonVpNum){
				//达到它的搜索上限，那么删除这个节点
#ifdef hd
				popFront_DchArray(dch_min_index);
#else
				popBack_DchArray(dch_min_index);
#endif
				dch_size--;

			}
			if (dch_size < MAX_DCH_QUEUE_SIZE){
				dch_size += pushBack_DchArray(findResult->f, findResult);
				dch_max_index = max(dch_max_index, findResult->f);
			}
			else if (findResult->f < dch_max_index){
				if (pushBack_DchArray(findResult->f, findResult)){
					//如果成功加入了，那就得再删除一个节点
#ifdef hd
					popFront_DchArray(dch_max_index);
#else
					popBack_DchArray(dch_max_index);
#endif
					dch_max_index = findNextMaxIndex(dch_max_index);
				}
			}
			
		}


	}
	if (bestF != 2147483647){
		while (bestNode->no != source){
			successPath[0]++;
			successPath[successPath[0]] = bestNode->no;
			for (int i = 1; i <= bestNode->path[0]; i++){
				successPath[0]++;
				successPath[successPath[0]] = bestNode->path[i];
			}
			bestNode = bestNode->father;
		}
		successPath[0]++;
		successPath[successPath[0]] = source;
        if(!ifItrate){
		  OutputResult_linux(successPath);
          return MAX_DIST;
        }
        else
		  return bestF;
	}
	return bestF;
}

void search_route_map(){
	Node * head = new Node(source);
	head->f = 0;
	dchQueue.insert(std::pair<int, Node*>(head->f, head));
	Node *findSource,* findResult;
	unsigned short successPath[MAX_NODE_NUM + 1];
	successPath[0] = 0;

	Node * bestNode = NULL;
	int bestF = 2147483647;
	unsigned short successNum = 0;

	int tmpcount = 0;

	tmpcount = 0;

	while (dchQueue.size()!=0){

		tmpcount = (tmpcount + 1) % POP_TOP_NUM;
		tmpcount = tmpcount % dchQueue.size();

		multimap<int, Node*>::iterator mapit = dchQueue.begin();
		unsigned short tmpcount_ = tmpcount;
		while (tmpcount_--)
			mapit++;
		findSource = mapit->second;
		findResult = dijistra_zhs(findSource);

		if (findResult == NULL){
			//Source节点找不出子节点了
			if (findSource->sonsInDemand[0] == 0){
				//说明它没有儿子节点，可以删除它
				//TODO
				//是不是只要做这一个操作？
				dchQueue.erase(mapit);
				delete findSource;
			}
			else{
				//有儿子，但是已经算不出其他儿子了，那么节点不删除，只从队列中删除
				dchQueue.erase(mapit);
			}
		}
		else if (findResult->no == target){
			//找到完整路径
			
			int tmpf = findResult->g;
			if (bestF>tmpf){
				bestF = tmpf;
				bestNode = findResult;
				successNum++;
			}
			
			// if (successNum < BESTROAD_NUM)
			continue;
			
			while (bestNode->no != source){
				successPath[0]++;
				successPath[successPath[0]] = bestNode->no;
				for (int i = 1; i <= bestNode->path[0]; i++){
					successPath[0]++;
					successPath[successPath[0]] = bestNode->path[i];
				}
				bestNode = bestNode->father;
			}
			successPath[0]++;
			successPath[successPath[0]] = source;
			//OutputResult(successPath, tmpf);
			OutputResult_linux(successPath);
			return;
		}
		else{
			dchQueue.insert(std::pair<int, Node*>(findResult->f, findResult));
		}

		if (dchQueue.size() > MAX_DCH_QUEUE_SIZE)
			dchQueue.erase(--dchQueue.end());
	}

	if (bestF != 2147483647){

		while (bestNode->no != source){
			successPath[0]++;
			successPath[successPath[0]] = bestNode->no;
			for (int i = 1; i <= bestNode->path[0]; i++){
				successPath[0]++;
				successPath[successPath[0]] = bestNode->path[i];
			}
			bestNode = bestNode->father;
		}
		successPath[0]++;
		successPath[successPath[0]] = source;
		//OutputResult(successPath, tmpf);
		OutputResult_linux(successPath);
		return;
	}
}


bool search_route_iteration_zhs(int lastBestf){

    if(!vphasdone){
       initDijiVpToVp();
       vphasdone = true;
    }


    unsigned short qsize = 400;
    Node * head = new Node(source);
    head->f = 0;
    multimap<int, Node*> zhsQueue[MAX_DEMAND_NUM + 1];
    unsigned short qLevel = 0;
    zhsQueue[0].insert(std::pair<int, Node*>(head->f, head));
    Node *findSource, *findResult;
    unsigned short successPath[MAX_NODE_NUM + 1];
    successPath[0] = 0;

    Node * bestNode = NULL;
    int bestF = 2147483647;
    unsigned short successNum = 0;

    multimap<int, Node*>::iterator mapit;
    while (1){
        while (zhsQueue[qLevel].empty()){
            qLevel++;
            if (qLevel == MAX_DEMAND_NUM+1)break;
        }
        if (qLevel == MAX_DEMAND_NUM+1)break;
        

        mapit = zhsQueue[qLevel].begin();
        findSource = mapit->second;
        findResult = dijistra_zhs(findSource);

        if (findResult == NULL){
            zhsQueue[qLevel].erase(mapit);
            if (findSource->sonsInDemand[0] == 0){
                delete findSource;
            }
        }
        else if (findResult->no == target){
            int tmpf = findResult->g;
            if (bestF > tmpf){
                bestF = tmpf;
                bestNode = findResult;
                successNum++;
            }
            //continue;
            break;
        }
        else if (findResult->f < lastBestf){
            //qLevel++;
            zhsQueue[qLevel + 1].insert(std::pair<int, Node*>(findResult->f, findResult));
        }
        if (zhsQueue[qLevel + 1].size() == qsize){
            //mapit = zhsQueue[qLevel + 1].begin();
            mapit = --zhsQueue[qLevel + 1].end();
            if (mapit->second->sonsInDemand[0] == 0)
                delete mapit->second;
            zhsQueue[qLevel + 1].erase(mapit);
        }
    }

    if (bestF != 2147483647){

        while (bestNode->no != source){
            successPath[0]++;
            successPath[successPath[0]] = bestNode->no;
            for (int i = 1; i <= bestNode->path[0]; i++){
                successPath[0]++;
                successPath[successPath[0]] = bestNode->path[i];
            }
            bestNode = bestNode->father;
        }
        successPath[0]++;
        successPath[successPath[0]] = source;
        OutputResult_linux(successPath);
        return true;
    }
    return false;
}

void search_route_big_array(){

	initDijiVpToVp();
	Node * head = new Node(source);
	head->f = 0;
	dchArray[0][0]=head;
	dchLength[0] = 1;
	unsigned short dch_size = 1;
	unsigned short dch_min_index = 0;
	int dch_max_index = 0;
	Node *findSource,* findResult;
	unsigned short successPath[MAX_NODE_NUM + 1];
	successPath[0] = 0;

	Node * bestNode=NULL;
	int bestF = 2147483647;
	unsigned short successNum = 0;

	
	int count_while = 0;

	unsigned short maxSonVpNum = MAX_VP_SON_NUM;

	while (dch_size != 0){
		
		++count_while;
		if(count_while == 500000){
			count_while = 0;
			if(GetTickCount()-timestart > 9000)
			
			    break;

		}
		
		dch_min_index = findNextMinIndex(dch_min_index);

		if (dch_min_index>bestF){
			dch_size = dch_size	-	delCor(dch_min_index);
			break;
		}

#ifdef hd
		findSource = dchArray[dch_min_index][dchBegin[dch_min_index]];
#else
		findSource = dchArray[dch_min_index][dchEnd[dch_min_index]];
#endif


		
		findResult = dijistra_zhs(findSource);
		if (findResult == NULL){
			//Source节点找不出子节点了
#ifdef hd
			popFront_DchArray(dch_min_index);
#else
			popBack_DchArray(dch_min_index);
#endif
			dch_size--;
			if (findSource->sonsInDemand[0] == 0){
				//说明它没有儿子节点，可以删除它
				//TODO
				//是不是只要做这一个操作？
				//dchQueue.erase(mapit);
				delete findSource;
			}
		}
		else if (findResult->no == target){
			//找到完整路径
			
			int tmpf = findResult->g;
			if (bestF>tmpf){
				bestF = tmpf;
				bestNode = findResult;
				successNum++;
			}
			
			//if (successNum < BESTROAD_NUM)
				continue;
			
			while (bestNode->no != source){
				successPath[0]++;
				successPath[successPath[0]] = bestNode->no;
				for (int i = 1; i <= bestNode->path[0]; i++){
					successPath[0]++;
					successPath[successPath[0]] = bestNode->path[i];
				}
				bestNode = bestNode->father;
			}
			successPath[0]++;
			successPath[successPath[0]] = source;
			//OutputResult(successPath, tmpf);
			OutputResult_linux(successPath);
			return;
		}
		else {
			if (findSource->sonsInDemand[0] == maxSonVpNum){
				//达到它的搜索上限，那么删除这个节点
#ifdef hd
				popFront_DchArray(dch_min_index);
#else
				popBack_DchArray(dch_min_index);
#endif
				dch_size--;

			}
			if (dch_size < MAX_DCH_QUEUE_SIZE){
				dch_size += pushBack_DchArray(findResult->f, findResult);
				dch_max_index = max(dch_max_index, findResult->f);
			}
			else if (findResult->f < dch_max_index){
				if (pushBack_DchArray(findResult->f, findResult)){
					//如果成功加入了，那就得再删除一个节点
#ifdef hd
					popFront_DchArray(dch_max_index);
#else
					popBack_DchArray(dch_max_index);
#endif
					dch_max_index = findNextMaxIndex(dch_max_index);
				}
			}
			
		}
	}
	if (bestF != 2147483647){
		while (bestNode->no != source){
			successPath[0]++;
			successPath[successPath[0]] = bestNode->no;
			for (int i = 1; i <= bestNode->path[0]; i++){
				successPath[0]++;
				successPath[successPath[0]] = bestNode->path[i];
			}
			bestNode = bestNode->father;
		}
		successPath[0]++;
		successPath[successPath[0]] = source;
		OutputResult_linux(successPath);
		return;
	}
	return;
}





bool search_route_iteration_zhs_big_graph(int lastBestf=MAX_DIST){

    if(!vphasdone){
       initDijiVpToVp();
       vphasdone = true;
    }


    unsigned short qsize = 600;
    Node * head = new Node(source);
    head->f = 0;
    multimap<int, Node*> zhsQueue[MAX_DEMAND_NUM + 1];
    unsigned short qLevel = 0;
    zhsQueue[0].insert(std::pair<int, Node*>(head->f, head));
    Node *findSource, *findResult;
    unsigned short successPath[MAX_NODE_NUM + 1];
    successPath[0] = 0;

    Node * bestNode = NULL;
    int bestF = 2147483647;
    unsigned short successNum = 0;

    multimap<int, Node*>::iterator mapit;
    while (1){
        while (zhsQueue[qLevel].empty()){
            qLevel++;
            if (qLevel == MAX_DEMAND_NUM+1)break;
        }
        if (qLevel == MAX_DEMAND_NUM+1)break;
        

        mapit = zhsQueue[qLevel].begin();
        findSource = mapit->second;
        findResult = dijistra_zhs(findSource);

        if (findResult == NULL){
            zhsQueue[qLevel].erase(mapit);
            if (findSource->sonsInDemand[0] == 0){
                delete findSource;
            }
        }
        else if (findResult->no == target){
            int tmpf = findResult->g;
            if (bestF > tmpf){
                bestF = tmpf;
                bestNode = findResult;
                successNum++;
            }
            //continue;
            break;
        }
        else if (findResult->f < lastBestf){
            //qLevel++;
            zhsQueue[qLevel + 1].insert(std::pair<int, Node*>(findResult->f, findResult));
        }
        if (zhsQueue[qLevel + 1].size() == qsize){
            //mapit = zhsQueue[qLevel + 1].begin();
            mapit = --zhsQueue[qLevel + 1].end();
            if (mapit->second->sonsInDemand[0] == 0)
                delete mapit->second;
            zhsQueue[qLevel + 1].erase(mapit);
        }
    }

    if (bestF != 2147483647){

        while (bestNode->no != source){
            successPath[0]++;
            successPath[successPath[0]] = bestNode->no;
            for (int i = 1; i <= bestNode->path[0]; i++){
                successPath[0]++;
                successPath[successPath[0]] = bestNode->path[i];
            }
            bestNode = bestNode->father;
        }
        successPath[0]++;
        successPath[successPath[0]] = source;
        OutputResult_linux(successPath);
        return true;
    }
    return false;
}


void search_route_double_array(){
        int tmp=search_route_array();
        ifItrate=false;
        clearDCHArray();
        MAX_DCH_QUEUE_SIZE = 2000;
        MAX_DCH_ARRAY_SIZE = 2000;
        search_route_array(tmp);
}



void search_route(char *graph[5000], int edge_num, char *condition){
	
	//timestart = GetTickCount(); 

	zhsinit_graph(graph, edge_num, condition);


	if(edge_num < 100){
		search_route_dfs();
        return;
    } else if(totalNodeNum < 150){
        search_route_bfs();
        return;
    } else if(totalNodeNum < 200){
        MAX_DCH_QUEUE_SIZE = 550;
        MAX_DCH_ARRAY_SIZE = 300;
        ifItrate=true;
        search_route_double_array();
        return;
    } else if(totalNodeNum < 250){
        search_route_bfs();
        return;
    } else if(totalNodeNum < 300){
 
		MAX_DCH_QUEUE_SIZE = 550;
		MAX_DCH_ARRAY_SIZE = 300;
        search_route_array();
        return;
    } else if(totalNodeNum < 350){

		MAX_DCH_QUEUE_SIZE = 550;
		MAX_DCH_ARRAY_SIZE = 300;
        ifItrate=true;
        MAX_VP_SON_NUM = 4;
        search_route_iteration_zhs_big_graph(search_route_array());
        return;
    } else if(totalNodeNum < 550){
        if(inDemandNum >= 25) {
    
            search_route_bfs();
        } else if(inDemandNum >= 23) {
  
			MAX_DCH_QUEUE_SIZE = 300;
			MAX_DCH_ARRAY_SIZE = 300;
            ifItrate=true;
            search_route_iteration_zhs(search_route_array());
        } else {
 
			MAX_DCH_QUEUE_SIZE = MAX_TOTAL_WEIGHT;
			MAX_DCH_ARRAY_SIZE = 50;
            ifItrate=true;
            search_route_iteration_zhs(search_route_array());
        }
    } else {
        if(inDemandNum >= 25) {
      
			MAX_DCH_QUEUE_SIZE = 300;
			MAX_DCH_ARRAY_SIZE = 300;
            search_route_iteration_zhs_big_graph(MAX_DIST);
        } else {

		MAX_DCH_QUEUE_SIZE = 550;
		MAX_DCH_ARRAY_SIZE = 300;
        MAX_VP_SON_NUM = 3;
        search_route_big_array();

        }
    }	
	return;
}