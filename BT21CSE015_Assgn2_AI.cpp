#include <bits/stdc++.h>
using namespace std;
class Graph
{
public:
	vector<vector<int>> adj;
	int n;
	Graph()
	{
		adj = vector<vector<int>>();
		n = 0;
	}
	Graph(int n)
	{
		this->n = n;
		adj.resize(n, vector<int>(n, -1));
	}
	Graph(const vector<vector<int>> &adj)
	{
		this->adj = adj;
		this->n = adj.size();
	}
	void addEdge(int u, int v, int weight)
	{
		adj[u][v] = weight;
		adj[v][u] = weight;
	}
};

class State
{
public:
	vector<int> path;
	int cost;
	int currNode; // Current node
	int gCost;	  // Cost taken till now to reach the current node
	int hCost;	  // Heuristic cost of the current node
	vector<bool> visited;
	set<int> unVisited;
	// Constructor with a given number of nodes, n
	State(int n)
	{
		path = vector<int>(1, 0);
		cost = 0;
		currNode = 0;
		gCost = 0;
		hCost = 0;
		visited = vector<bool>(n, false);
		visited[0] = true;
		for (int i = 1; i < n; i++)
		{
			unVisited.insert(i);
		}
	}

	// Constructor with given data members
	State(const vector<int> &path, int cost, int currNode, int gCost, int hCost, const set<int> &unVisited)
	{
		this->path = path;
		this->cost = cost;
		this->currNode = currNode;
		this->gCost = gCost;  
		this->hCost = hCost;
		this->visited = visited;
		this->unVisited = unVisited;
	}

	// Operator overloading for copying a State object
	State &operator=(const State &other)
	{
		path = other.path;
		cost = other.cost;
		currNode = other.currNode;
		gCost = other.gCost;
		hCost = other.hCost;
		visited = other.visited;
		unVisited = other.unVisited;
		return *this;
	}
};
// Class used as a comparator for the priority queue to sort states in ascending order of f = g + h
class CompareA
{
public:
	bool operator()(const State &a, const State &b)
	{
		return a.gCost + a.hCost > b.gCost + b.hCost; // Sort in  order of f = g + h
	}
};

// Class to solve the Traveling Salesperson Problem (TSP) using Minimum Spanning Tree (MST) heuristic
class TSP
{
public:
	Graph g;
	int n;
	priority_queue<State, vector<State>, CompareA> pq; //Min heap
	// Constructor with a given graph object
	TSP(const Graph &graph)
	{
		this->g = graph;
		this->n = graph.adj.size();
	}
	// Method to find the cost of the given Minimum Spanning Tree (MST)
	int costOfMST(vector<int> &parent, const set<int> &unVisited)
	{
		int cost = 0;
		for (auto &node : unVisited)
		{
			if (parent[node] == -1)
				continue;
			if (g.adj[node][parent[node]] == -1)
				continue;
			cost += g.adj[node][parent[node]];
		}
		return cost;
	}
	// Method to find the minimum key value in Prim's MST Algorithm
	int minKey(vector<int> &key, set<int> &mstSet, const set<int> &unVisited)
	{
		int mn = INT_MAX, minIndex;
		for (auto &x : unVisited)
		{
			if (mstSet.find(x) == mstSet.end() && key[x] < mn)
			{
				mn = key[x];
				minIndex = x;
			}
		}
		return minIndex;
	}

	vector<int> primsMinimumSpanningTree(const set<int> &unv)
	{
		set<int> unVisited = unv;
		if (unVisited.size() == 0){ 
			return vector<int>(n, -1);
		}
		vector<int> parent(n), key(n);// keys to keep track of min weights
		set<int> mstSet;// To store what is included in mst
		for (auto &v : unVisited)
		{
			key[v] = INT_MAX;
		}
		int start = *unVisited.begin();
		key[start] = start;
		parent[start] = -1;
		int V = unVisited.size();
		for (int i = 0; i < V - 1; i++)
		{
			int u = minKey(key, mstSet, unVisited);
			mstSet.insert(u);
			// Update key values and parent indices of adjacent vertices of u
			for (auto &x : unVisited)
			{
				if (mstSet.find(x) != mstSet.end())
					continue;
				if (g.adj[u][x] != -1 && g.adj[u][x] < key[x])
				{
					parent[x] = u;
					key[x] = g.adj[u][x];
				}
			}
		}
		return parent;
	}
	// Method to print the solution path
	void printSolution(const State &currState, int destination, int expandedNodes)
	{
		cout << "\nPath: ";
		for (int i = 0; i < currState.path.size(); i++)
		{
			cout << currState.path[i] << " -> ";
		}
		cout << destination << endl;
		cout << "\nNumber of nodes expanded: " << expandedNodes << endl;
	}
	int AStar()
	{
		int start = 0, destination = 0;
		int expandedNodes = 0;
		vector<int> pth = {start}; 
		set<int> unv;		
		for (int i = 0; i < n; i++)
		{
			if (i != start)
			{
				unv.insert(i);
			}
		}
		State startState(pth, 0, start, 0, 0, unv);
		// Find the heuristic cost of the start state
		vector<int> parent = primsMinimumSpanningTree(startState.unVisited);
		startState.hCost = costOfMST(parent, startState.unVisited);
		pq.push(startState);
		while (!pq.empty())
		{
			State currState = pq.top();
			pq.pop();
			expandedNodes++;
			// Destination node reached, and all nodes are visited
			if (currState.path.size() == n && g.adj[currState.currNode][destination] != -1)
			{
				// Print the path from the start node to the destination node
				printSolution(currState, destination, expandedNodes);
				return (currState.gCost + g.adj[currState.currNode][destination]);
			}
			// Expand over neighbors
			for (int i = 0; i < n; i++)
			{
				// If the node is visited then continue 
				if (currState.unVisited.find(i) == currState.unVisited.end())
					continue;
				if (g.adj[currState.currNode][i] == -1) //there is no an edge between the current node and the node so continue 
					continue;
				vector<int> newPath = currState.path;
				newPath.push_back(i);//adding currunt node in path 
				set<int> newUnv = currState.unVisited;
				newUnv.erase(i);// remove it from unvisited 
				State newState(newPath, 0, i, currState.gCost + g.adj[currState.currNode][i], 0, newUnv);
				vector<int> parent = primsMinimumSpanningTree(newState.unVisited);
				newState.hCost = costOfMST(parent, newState.unVisited);
				pq.push(newState);
			}
		}

		return -1;
	}
};

Graph readInput(string inputFile)
{
	ifstream fin(inputFile);
	int n;
	fin >> n;
	Graph g(n);
	for (int i = 0; i < g.adj.size(); i++)
	{
		for (int j = 0; j < g.adj.size(); j++)
		{
			int x;
			fin >> x;
			if (j >= i || x == -1)
				continue;
			if (x > 0)
				g.addEdge(i, j, x);
		}
	}
	fin.close();
	return g;
}

int main()
{
	Graph graph = readInput("input50_backup.txt");
	TSP tsp(graph);
	int res = tsp.AStar();
	if (res == -1)
		cout << "\nNo path found." << endl;
	else
		cout << "\nCost of the path: " << res << endl;
	return 0;
}