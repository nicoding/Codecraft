
#include "route.h"
#include "lib_io.h"
#include "lib_time.h"

using namespace std;

/*void init_org(char** &graph,int &edge_num){
graph[0]="0,0,13,15";
graph[1]="1,0,8,17";
graph[2]="2,0,19,1";
graph[3]="3,0,4,8";
graph[4]="4,1,0,4";
graph[5]="5,2,9,19";
graph[6]="6,2,15,8";
graph[7]="7,3,0,14";
graph[8]="8,3,11,12";
graph[9]="9,4,1,15";
graph[10]="10,4,5,17";
graph[11]="11,5,8,18";
graph[12]="12,5,9,14";
graph[13]="13,5,6,2";
graph[14]="14,6,17,4";
graph[15]="15,7,13,1";
graph[16]="16,7,16,19";
graph[17]="17,8,6,1";
graph[18]="18,8,12,17";
graph[19]="19,9,14,11";
graph[20]="20,10,12,1";
graph[21]="21,11,7,12";
graph[22]="22,11,4,7";
graph[23]="23,12,14,5";
graph[24]="24,13,17,12";
graph[25]="25,13,4,2";
graph[26]="26,14,19,9";
graph[27]="27,15,10,14";
graph[28]="28,15,18,2";
graph[29]="29,16,8,1";
graph[30]="30,17,9,14";
graph[31]="31,17,19,3";
graph[32]="32,17,18,10";
graph[33]="33,18,15,8";
graph[34]="34,18,3,8";
graph[35]="35,19,18,12";
graph[36]="36,2,3,20";
graph[37]="37,3,5,20";
graph[38]="38,5,7,20";
graph[39]="39,7,11,20";
graph[40]="40,11,13,20";
graph[41]="41,17,11,20";
graph[42]="42,11,19,20";
graph[43]="43,17,5,20";
graph[44]="44,5,19,20";
edge_num=45;
}*/

/*#include <windows.h>
void main(){
	char **graph;
	graph=new char*[5000];
	char * condition;

	int edge_num;
	
	//mini-scale test case
	graph[0] = "0,0,1,1";
	graph[1] = "1,0,2,2";
	graph[2] = "2,0,3,1";
	graph[3] = "3,2,1,3";
	graph[4] = "4,3,1,1";
	graph[5] = "5,2,3,1";
	graph[6] = "6,3,2,1";
	edge_num = 7;
	condition = "0,1,2|3";
	

	
	//small-scale test case
	init_org(graph,edge_num);
	condition="2,19,3|5|7|11|13|17";
	
	
	
	//middle-scale test case
	init(graph,edge_num);
	condition="0,1,4|7|8|14|24|44|33|32|25|48";
	//condition="19,87,55|22|33|15|198|123|134|156|255|258|236|77|27|233|85|20|66|222|238|79";

	cout<<"----------start----------------"<<endl;
	DWORD dwStart,dwEnd;
	dwStart = GetTickCount();
	//bool a[5][5] = { false };
	//cout << a[2][2] << endl;
	
	search_route(graph, edge_num, condition);
	
	dwEnd = GetTickCount();
	cout<<'\n'<<dwEnd-dwStart<<endl;
	system("pause");
}*/

int main(int argc, char *argv[]){
    print_time("Begin");
	char *graph[5000];
    int edge_num;
    char *condition;
    int condition_num;

    char *graph_file = argv[1];
    edge_num = read_file(graph, 5000, graph_file);
    char *condition_file = argv[2];
    condition_num = read_file(&condition, 1, condition_file);

	//DWORD dwStart,dwEnd;
	//dwStart = GetTickCount();
	//bool a[5][5] = { false };
	//cout << a[2][2] << endl;
	
	search_route(graph, edge_num, condition);


    char *result_file = argv[3];
    write_result(result_file);
    release_buff(graph, edge_num);
    release_buff(&condition, condition_num);
	print_time("End");

	//dwEnd = GetTickCount();
	//cout<<'\n'<<dwEnd-dwStart<<endl;
	return 0;
}