/*
The GAlib-based genetic algorithm code for the Travelling Salesman Problem (TSP) Berlin52.
作者：wying
单位：华南理工大学软件学院
*/

#include <math.h>
#include <ctime>
#include "ga/GASStateGA.h"
#include "ga/GASimpleGA.h"
#include "ga/GA1DArrayGenome.h"
#include "ga/garandom.h"
#include "ga/std_stream.h"

#define cout STD_COUT
#define cerr STD_CERR
#define endl STD_ENDL
#define ostream STD_OSTREAM
#define ifstream STD_IFSTREAM

// Set this up for your favorite TSP.  The sample one is a contrived problem
// with the towns laid out in a grid (so it is easy to figure out what the 
// shortest distance is, and there are many different paths with the same
// shortest path).  File format is that used by the TSPLIB problems.  You can 
// grab more problems from TSPLIB.
// 
#define MAX_TOWNS 100
#define TSP_FILE "berlin52.txt"

int x[MAX_TOWNS],y[MAX_TOWNS];//每个城市的x坐标和y坐标
int DISTANCE[MAX_TOWNS][MAX_TOWNS];//每两个城市之间的旅行成本，是对称的


float TSPObjective(GAGenome&);//计算染色体的旅行总费用的目标函数
void  TSPInitializer(GAGenome&);//TSP问题的染色体初始化算子
int   TSPMutator(GAGenome&, float);//针对TSP问题的染色体变异算子
int   TSPCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);//针对TSP问题的染色体交叉算子

void writeTSPPath(ostream & os, GAGenome& g);//将指定染色体的旅行路线输出到指定文件

int main() {
  cout << "The GAlib program for the Travelling Salesman Problem (TSP) Berlin52.\n" << endl;

  //从Berlin52.txt文件读出各城市坐标
  double CityID;
  ifstream in(TSP_FILE); 
  if(!in) {
    cerr << "could not read data file " << TSP_FILE << "\n";
    exit(1);
  }
  int ntowns=0;
  do {
    in >> CityID;
    in >> x[ntowns];
    in >> y[ntowns];
    ntowns++;
  } while(!in.eof() && ntowns < MAX_TOWNS);
  in.close();
  if(ntowns >= MAX_TOWNS) {
    cerr << "data file contains more towns than allowed for in the fixed\n";
    cerr << "arrays.  Recompile the program with larger arrays or try a\n";
    cerr << "smaller problem.\n";
    exit(1);
  }

  //计算任意两个城市间的旅行成本
  double dx,dy;
  for(int i=0;i<ntowns;i++) {
    for(int j=i; j<ntowns;j++) {
      dx=x[i]-x[j]; dy=y[i]-y[j];
      DISTANCE[j][i]=DISTANCE[i][j]=floor(0.5 + sqrt(dx*dx+dy*dy) );//注意取四舍五入之后的整数值
    }
  }

  //定义TSP问题的编码方案为一维的整数数组，其固定长度为城市个数
  GA1DArrayGenome<int> genome(ntowns);

  genome.evaluator(::TSPObjective);//为染色体指定计算目标值的函数
  genome.initializer(::TSPInitializer);//为染色体指定自定义的初始化算子
  genome.crossover(::TSPCrossover);//为染色体指定自定义的交叉算子
  genome.mutator(::TSPMutator);//为染色体指定自定义的变异算子

  GASteadyStateGA ga(genome); ga.nReplacement(2); ga.nGenerations(500000);//选用稳态遗传算法进行TSP问题求解，指定其染色体编码方式、每一代要替换的个体数=2、总的运行代数500000，那么搜索的总个体数=2*500000=1000000
  //GASimpleGA ga(genome); ga.elitist(gaTrue); ga.nGenerations(5000);//选用简单遗传算法进行TSP问题求解，采用精英保留策略，指定其染色体编码方式、总的运行代数10000，那么搜索的总个体数=200（种群大小）*5000=1000000
  ga.minimize();//为遗传算法指定优化目的是将目标函数值最小化
  ga.populationSize(200);//为遗传算法指定种群大小为200
  ga.pMutation(0.02);//为遗传算法指定变异概率
  ga.pCrossover(0.8);//为遗传算法指定交叉概率

  cout << "initializing..."<<"\n"; cout.flush();
  unsigned int seed = clock();
  ga.initialize(seed);//使用从时钟得到的随机种子初始化遗传算法

  cout << "evolving..."<<"\n"; cout.flush();
  std::fstream fgacurve;
  fgacurve.open("tspgacurve.txt",std::ios::out);

  //遗传算法开始迭代进化，直到达到指定的代数
  while(!ga.done()) {
    ga.step();//进化一代
    if(ga.generation() % (ga.nGenerations()/100) == 0) 
	{//进化过程中取100个采样点，记录进化过程中的最优目标值收敛信息到文件
		int bestscore = ga.statistics().bestIndividual().score();
		cout << ga.generation() << "    " << bestscore << "\n"; cout.flush();
		fgacurve << ga.generation() << "    " << bestscore << "\n";
    }
  }
  fgacurve.close();

  //遗传算法迭代终止后输出找到的最优旅行路线到文件
  genome = ga.statistics().bestIndividual();
  //cout << "\n" << "the shortest path found is "  << "\n";
  //writeTSPPath(cout, genome);
  std::fstream fbestpath;
  fbestpath.open("tsppath.txt",std::ios::out);
  writeTSPPath(fbestpath, genome);
  fbestpath.close();
  cout << "the distance of the shortest path found: " << genome.score() << "\n";

  return 0;
}


// Here are the genome operators that we want to use for this problem.
//计算染色体的旅行总费用的目标函数
float TSPObjective(GAGenome& g) {
  GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
  int genomelength = genome.size();//genome.size()获取染色体的长度
  float dist = 0;
  int xx;
  int yy;

    for(int i=0; i<genomelength; i++) {
      xx = genome.gene(i);
      yy = genome.gene( (i+1)%genomelength );
      dist += DISTANCE[xx-1][yy-1];
    }

  return dist;
}

//TSP问题的染色体初始化算子
void TSPInitializer(GAGenome& g) {
  GA1DArrayGenome<int> &genome=(GA1DArrayGenome<int> &)g;

  int genomelength = genome.size();
  int i,town;
  static bool visit[MAX_TOWNS];

  memset(visit, false, MAX_TOWNS*sizeof(bool));
  town=GARandomInt(1,genomelength);//GARandomInt(1,genomelength)生成1到genomelength之间的一个均匀随机整数
  visit[town-1]=true;
  genome.gene(0, town);//genome.gene(0, town)设置该染色体第0个基因位上的基因值为town
 
  for( i=1; i<genomelength; i++) {
    do {
      town=GARandomInt(1,genomelength);
    } while (visit[town-1]);
    visit[town-1]=true;
    genome.gene(i, town);
  }	
}

//针对TSP问题的染色体变异算子，pmut为变异概率
int TSPMutator(GAGenome& g, float pmut) {
  GA1DArrayGenome<int> &genome=(GA1DArrayGenome<int> &)g;
  int i;

  int genomelength = genome.size();
  float nmutator = pmut*genomelength;//要改变的边的数量

  int imutator=0;
  while( imutator<nmutator){
    if (GARandomFloat()<0.5) {//GARandomFloat()生成0到1之间的一个均匀随机浮点数
	  //以0.5概率使用相互交换变异
	  int swapIndex1 = GARandomInt(0,genomelength-1);
	  int swapIndex2 = GARandomInt(0,genomelength-1);
	  int tmp;
	  tmp = genome.gene(swapIndex2);
	  genome.gene(swapIndex2, genome.gene(swapIndex1) );
	  genome.gene(swapIndex1, tmp );// swap only one time
	  imutator+=4;
    }else
	  {
	  //以0.5概率使用反转变异
	  int inversion_start, inversion_end, tmp;
	  inversion_start = GARandomInt(0,genomelength-1);
	  inversion_end = GARandomInt(0,genomelength-1);
	  if(inversion_start > inversion_end)
	  {
		  tmp = inversion_start;
		  inversion_start = inversion_end;
		  inversion_end = tmp;
	  }
	   
	  for(i=inversion_start;inversion_start<inversion_end; inversion_start++, inversion_end-- )
	  {
		  tmp = genome.gene(inversion_start);
		  genome.gene(inversion_start, genome.gene(inversion_end) );
		  genome.gene(inversion_end, tmp ); 
	  }  
	  imutator+=2;
    }
  }

  return (1);
}


//针对TSP问题的染色体交叉算子
// 针对TSP问题的染色体交叉算子
int TSPCrossover(const GAGenome& g1, const GAGenome& g2, GAGenome* c1,
	GAGenome* c2) {
	GA1DArrayGenome<int> parent1 = (GA1DArrayGenome<int> &)g1;
	GA1DArrayGenome<int> parent2 = (GA1DArrayGenome<int> &)g2;

	int genomelength = parent1.size();
	// maxto

	int nc = 0;

	GA1DArrayGenome<int>& child1 = (GA1DArrayGenome<int> &) * c1;
	GA1DArrayGenome<int>& child2 = (GA1DArrayGenome<int> &) * c2;

	if (c1) {
		child1 = parent2;
		nc++;
	}
	if (c2) {
		child2 = parent1;
		nc++;
	}

	/////
	// 此处添加代码实现自己的交叉算子
	/////
	int start, end;
	start = GARandomInt(0, genomelength - 1);
	end = GARandomInt(0, genomelength - 1);
	if (start > end) {
		int tmp = start;
		start = end;
		end = tmp;
	}
	// visit数组用于记录基因是否已经访问过, 0表示不存在, 1表示存在
	int visit1[MAX_TOWNS] = { 0 };
	int visit2[MAX_TOWNS] = { 0 };
	// 交换child中的基因, 并将交换过的基因分别标记为已访问
	for (int i = start; i <= end; i++) {
		int tmp = child1.gene(i);
		child1.gene(i, child2.gene(i));
		child2.gene(i, tmp);
		visit1[child1.gene(i) - 1] = 1;
		visit2[child2.gene(i) - 1] = 1;
	}
	// 将parent中visit数组中为0的基因插入到child中
	int index1 = 0, index2 = 0;
	for (int i = 0; i < genomelength; i++) {
		if (i >= start && i <= end) continue;
		while (visit1[parent1.gene(index1) - 1] == 1) index1++;
		while (visit2[parent2.gene(index2) - 1] == 1) index2++;
		child1.gene(i, parent1.gene(index1));
		child2.gene(i, parent2.gene(index2));
		visit1[parent1.gene(index1) - 1] = 1;
		visit2[parent2.gene(index2) - 1] = 1;
	}

	return nc;
}

//将指定染色体的旅行路线输出到指定文件
void writeTSPPath(ostream & os, GAGenome& g) 
{
	GA1DArrayGenome<int> & genome = (GA1DArrayGenome<int> &)g;
	int genomelength = genome.size();
	for(int i=0; i<genomelength; i++)
	{
		int xx = genome.gene(i);
		os << xx <<"    " <<x[xx-1] << "      "<<y[xx-1] << "\n";
	}
}
