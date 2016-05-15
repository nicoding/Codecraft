#ifndef __ROUTE_H__
#define __ROUTE_H__
#include <stdlib.h>
#include <stdio.h>
#include<iostream>
#include<map>
#include <vector>
#include<stack>
//#include<unordered_set>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <string.h>
#include <time.h> 
#include <queue>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <cstring>
#include <climits>
#include <set>
#include <string>
#include <sstream>
#include <bitset>
using namespace std;

#define PWW 100000
#define MAX_OUT_DEGREE 8
#define MAX_NODE_NUM 600
#define MAX_SHORT_NUM 32766
#define A_NEGATIVE_NUM -666
#define MAX_BFS_SIZE 4096
#define MAX_TOTAL_WEIGHT 12000
#define MAX_DEMAND_NUM 50
#define MAX_EDGE 4800
#define MAX_DIST 16843009

#define BESTROAD_NUM 30  //每次找出x条满足要求的路径，返回最小值

#define POP_TOP_NUM 10   //对dch队列的头上x条路进行扩展

void search_route(char *graph[5000], int edge_num, char *condition);
void SplitCharStar(char *c,unsigned short re[]);
void SplitCharStarToSet(char *c);



#endif
