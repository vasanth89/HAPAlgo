#include <iostream>
using namespace std;

//Structure to store time, probability and cost  
struct Entry{
	int time;
    double prob;
    int cost;
};

//Struct that links each node with entries for a certain manufacturer
struct node_mfr_link{
    int mfrId;
    int nodeId;
    int no_of_entries;
    Entry *entries;
};

//Structure for input table
struct node_collection{
  int no_of_mfr;
  int no_of_nodes;
  node_mfr_link *set;
};

//Structure that links each node with time constraint
struct time_node_link{
	int nodeId;
	int time;
	int no_of_pairs;
	Entry *entries;
};

//Structure for final table
struct time_collection{
	int no_of_nodes;
	int maxTime;
	time_node_link *timeset;
};

//Structure declared to get various combinations for a give time constraint
struct combination{
    int x;
    int y;
};

//Initialization of various methods used
void displayTable(node_collection table);
void calculate_time_scale(int i);
int calculate_local_table_entries(int offset, int time, int localTableOffset);
void removeAt(int location, int no_of_pairs);
void selectionSort(int no_of_pairs);
int minLocation(int first, int last);
void swap(int first, int last);
bool compare(int item,int no_of_pairs);
void calculate_final_table_timescale(int i);
int get_combinations(int i, combination *combo, int node);
int remove_redundant_pairs(int no_of_pairs);
void cummulative_calculation(time_node_link t1, time_node_link t2, int &no_of_combination_pairs);

//Global Variable Declaration
node_collection modifiedTable;
node_collection table;
time_collection localTable;
time_collection finalTable;
Entry pairs[100];
int maxTime = 0;
int finalTableTime = 0;

//Starting point of program
int main(){
	
	int i,j,k;
	double cummProb = 0.0;
	int offset;
	int tempTime= 0;
	
	//Getting number of nodes and types to decide size of table
	cout<<"Enter number of nodes :"<<endl;
	cin>>table.no_of_nodes;
	cout<<"Enter number of types :"<<endl;
	cin>>table.no_of_mfr;

	//Setting the number of entries in input table
	table.set = new node_mfr_link[table.no_of_mfr*table.no_of_nodes];

	//Getting input data and setting values
	offset = 0;
	for(i=0; i<table.no_of_nodes;i++){
		for(j=0; j<table.no_of_mfr; j++){
			node_mfr_link t1;
			t1.nodeId = i;
			t1.mfrId = j;
			cout<<"Enter number of entries for node["<<i<<"]"<<"and Manufacturer["<<j<<"]"<<endl;
			cin>>t1.no_of_entries;
			if(t1.no_of_entries>0){
				t1.entries = new Entry[t1.no_of_entries];
				for(k=0;k<t1.no_of_entries;k++){
					Entry e1;
					cout<<"Enter time:"<<endl;\
					cin>>e1.time;
					cout<<"Enter probability:"<<endl;
					cin>>e1.prob;
					cout<<"Enter cost:"<<endl;
					cin>>e1.cost;
					t1.entries[k] = e1;
				}
			}
			table.set[offset+j] = t1;
			t1.entries=NULL;
			t1.no_of_entries = 0;
		}
		offset = offset+table.no_of_mfr;
	}

	//Display the new table
	displayTable(table);
	
	offset = 0;
	//Finding the maximum time for constructing local table
	for(i=0; i<table.no_of_nodes;i++){
		for(j=0; j<table.no_of_mfr; j++){
			node_mfr_link t1;
			t1 = table.set[offset+j];
			for(k=0;k<t1.no_of_entries;k++){
				Entry e1 = t1.entries[k];
				if(maxTime<e1.time)
					maxTime = e1.time;
			}
			t1.entries=NULL;
			t1.no_of_entries = 0;
		}
		offset = offset+table.no_of_mfr;
	}
	
	cout<<"Max Time :"<<maxTime<<endl;

	//Calculating values for modified table with cummulative probabilities
	modifiedTable = table;
	offset = 0;
	for(i=0; i<table.no_of_nodes;i++){
		for(j=0; j<table.no_of_mfr; j++){
			node_mfr_link t1;
			t1 = table.set[offset+j];
			if(t1.no_of_entries>1){
				cummProb = t1.entries[0].prob;
				for(k=1;k<=t1.no_of_entries;k++){
						cummProb = cummProb + t1.entries[k].prob;
						t1.entries[k].prob = cummProb;
				}
				modifiedTable.set[offset+j] = t1;
			}
			t1.entries=NULL;
			t1.no_of_entries = 0;
		}
		offset = offset+table.no_of_mfr;
	}

	//Displaying modified table
	displayTable(modifiedTable);

	//Initializing local table
	localTable.no_of_nodes = table.no_of_nodes;
	localTable.maxTime = maxTime;
	localTable.timeset = new time_node_link[localTable.no_of_nodes*localTable.maxTime];

	//Generating local table with probability, cost pairs 
	for(i=0;i<localTable.no_of_nodes;i++){
	calculate_time_scale(i);
	}

	offset = 0;
	//Displaying local table
	for(i=0; i<localTable.no_of_nodes;i++){
		for(j=0; j<localTable.maxTime; j++){
			time_node_link t1;
			t1 = localTable.timeset[offset+j];
			cout<<"Node["<<t1.nodeId<<"]"<<endl;
			cout<<"Time["<<t1.time<<"]"<<endl;
			for(int k=0;k<t1.no_of_pairs;k++){
				Entry e1 = t1.entries[k];
				cout<<"Entry["<<k<<"]\tProbability:"<<e1.prob<<"\tCost:"<<e1.cost<<endl;
			}
			t1.entries=NULL;
			t1.no_of_pairs = 0;
		}
		offset = offset+localTable.maxTime;
	}
	
	offset = 0;
	maxTime = 0;
	//Finding the maximum time for final table
	for(i=0; i<table.no_of_nodes;i++){
		for(j=0; j<table.no_of_mfr; j++){
			node_mfr_link t1;
			t1 = table.set[offset+j];
			for(k=0;k<t1.no_of_entries;k++){
				Entry e1 = t1.entries[k];
				if(maxTime<e1.time)
					maxTime = e1.time;
			}
			t1.entries=NULL;
			t1.no_of_entries = 0;
		}
		finalTableTime = finalTableTime + maxTime;
		offset = offset+table.no_of_mfr;
	}

	cout<<"Final Table Time:"<<finalTableTime<<endl;

	//Initializing Final Table
	finalTable.maxTime = finalTableTime;
	finalTable.no_of_nodes = table.no_of_nodes;
	finalTable.timeset = new time_node_link[finalTable.no_of_nodes*finalTable.maxTime];

	//Setting first row of local table to first row of modified table
	for(j=0;j<localTable.maxTime;j++){
		finalTable.timeset[j].nodeId = localTable.timeset[j].nodeId;
		finalTable.timeset[j].time = localTable.timeset[j].time;
		finalTable.timeset[j].no_of_pairs = localTable.timeset[j].no_of_pairs;
		finalTable.timeset[j].entries = localTable.timeset[j].entries;
	}
	for(j=localTable.maxTime;j<finalTable.maxTime;j++){
		finalTable.timeset[j].nodeId = finalTable.timeset[j-1].nodeId;
		finalTable.timeset[j].time = j+1;
		finalTable.timeset[j].no_of_pairs = finalTable.timeset[j-1].no_of_pairs;
		finalTable.timeset[j].entries = finalTable.timeset[j-1].entries;
	}

	//Calculating probability, cost pairs of final table
	for(i=1;i<finalTable.no_of_nodes;i++){
	calculate_final_table_timescale(i);
	}

	offset = 0;
	//Displaying final table
	for(i=0; i<finalTable.no_of_nodes;i++){
		for(j=0; j<finalTable.maxTime; j++){
			time_node_link t1;
			t1 = finalTable.timeset[offset+j];
			cout<<"Node["<<t1.nodeId<<"]"<<endl;
			cout<<"Time["<<t1.time<<"]"<<endl;
			for(int k=0;k<t1.no_of_pairs;k++){
				Entry e1 = t1.entries[k];
				cout<<"Entry["<<k<<"]\tProbability:"<<e1.prob<<"\tCost:"<<e1.cost<<endl;
			}
			t1.entries=NULL;
			t1.no_of_pairs = 0;
		}
		offset = offset+finalTable.maxTime;
	}

	cin>>i;
	return 0;
}

//This method is used to calculate the probability,cost pairs of final table
void calculate_final_table_timescale(int node){
	int offset = node * finalTable.maxTime;
	combination combo[100];
	int no_of_combination_pairs=0;
	for(int i=0; i<finalTable.maxTime; i++){
		finalTable.timeset[offset+i].nodeId = node;
		finalTable.timeset[offset+i].time = i + 1;
		finalTable.timeset[offset+i].no_of_pairs = 0;
		if((i+1)>node){
			//This step calls the method get_combinations to get the possible combinations for a given time constraint
			int no_of_combinations = get_combinations(finalTable.timeset[offset+i].time, combo, node);
			for(int j=0;j<no_of_combinations;j++){
				int x = combo[j].x-1 + (finalTable.maxTime*(node-1)) ;
				//If time is greater that local table's maximum time
				if(combo[j].y>localTable.maxTime)
					continue;
				int y = combo[j].y-1 + node * localTable.maxTime;
				//This step calls a method to get the entire list of all possible combinations of probability and cost
				cummulative_calculation(finalTable.timeset[x],localTable.timeset[y],no_of_combination_pairs);
			}
			//This step will remove the redundant pairs from the list using the respective algorithms
			finalTable.timeset[offset+i].no_of_pairs = remove_redundant_pairs(no_of_combination_pairs);
			no_of_combination_pairs=0;
			//Adding entries to the final table after filtering the redundant pairs
			if(finalTable.timeset[offset+i].no_of_pairs>0){
				finalTable.timeset[offset+i].entries = new Entry[finalTable.timeset[offset+i].no_of_pairs];
				for(int k=0; k<finalTable.timeset[offset+i].no_of_pairs;k++){
					finalTable.timeset[offset+i].entries[k] = pairs[k];
				}
			}
		}
	}
}

//This method is used to calculate the probability cost pairs for all entries of each node with previous entry of final table
void cummulative_calculation(time_node_link t1, time_node_link t2, int &no_of_combination_pairs){
	Entry entry;
	if(t1.no_of_pairs > 0 && t2.no_of_pairs > 0){
		for(int i=0;i<t1.no_of_pairs;i++){
			for(int j=0;j<t2.no_of_pairs;j++){
				entry.prob = t1.entries[i].prob * t2.entries[j].prob;
				entry.cost = t1.entries[i].cost + t2.entries[j].cost;
				pairs[no_of_combination_pairs] = entry;
				no_of_combination_pairs++;
			}
		}
	}
}

//This method contains the algorithm to remove redundant pairs
int remove_redundant_pairs(int no_of_pairs){
	selectionSort(no_of_pairs);

	bool del = false;
	for(int k=0;k<no_of_pairs;k++){
		//This step calls compare() method which removes redundant pairs
		del = compare(k,no_of_pairs);
		if(del){
			k=-1;
			no_of_pairs--;
		}
	}
	return no_of_pairs;
}

//This method is used to get all possible combinations for a given time constraint
int get_combinations(int time, combination *combo, int node){
	int n=0;
	for(int i=1;i<time;i++){
		for(int j=1;j<=finalTable.maxTime;j++){
			if((i+j) > time)
				break;
			else if ((i+j) != time)
				continue;
			combo[n].x = i;
			combo[n].y = j;
			n++;
		}
	}
	return n;
}

//This method is used to generate probability,cost pairs for local table for each node
void calculate_time_scale(int i){
	int offset = i * localTable.maxTime;
	for(int j=0; j<localTable.maxTime; j++){
		localTable.timeset[offset+j].time = j+1;
		localTable.timeset[offset+j].nodeId = i;
		localTable.timeset[offset+j].no_of_pairs = calculate_local_table_entries(i,localTable.timeset[offset+j].time,offset+j);
		if(localTable.timeset[offset+j].no_of_pairs>0){
			localTable.timeset[offset+j].entries = new Entry[localTable.timeset[offset+j].no_of_pairs];
			for(int k=0; k<localTable.timeset[offset+j].no_of_pairs;k++){
				localTable.timeset[offset+j].entries[k] = pairs[k];
			}
		}
	}
}

//This table is used to calculate entries for local table
int calculate_local_table_entries(int i, int time, int localTableOffset){
	int offset = i * modifiedTable.no_of_mfr;
	int no_of_pairs = 0;
	for(int k=0; k<modifiedTable.no_of_mfr;k++){
		node_mfr_link t1 = modifiedTable.set[offset+k];
		for(int l=0; l<t1.no_of_entries;l++){
			if(t1.entries[l].time<=time){
				pairs[no_of_pairs].prob = t1.entries[l].prob;
				pairs[no_of_pairs].cost = t1.entries[l].cost;
				no_of_pairs++;
			}
		}
	}
	//This step is called to sort the list of probability,cost pairs in ascending order of probability
	selectionSort(no_of_pairs);

	//This step will take care of redundant pairs
	bool del = false;
	for(int k=0;k<no_of_pairs;k++){
		del = compare(k,no_of_pairs);
		if(del){
			k=-1;
			no_of_pairs--;
		}
	}
	return no_of_pairs;
}

//This method is important as it contains the algorithm to remove redundant pairs
bool compare(int item,int no_of_pairs){
	bool del = false;
	int minIndex;
	minIndex = item;
	for (int loc = item + 1; loc < no_of_pairs; loc++){
		if(pairs[minIndex].prob == pairs[loc].prob){
				if(pairs[minIndex].cost >= pairs[loc].cost){
					removeAt(minIndex, no_of_pairs);
					del = true;
					break;
				} else{
					removeAt(loc, no_of_pairs);
					del = true;
					break;
				}
		}
		else{
			if(pairs[minIndex].cost >= pairs[loc].cost){
				removeAt(minIndex, no_of_pairs);
				del = true;
				break;
			}
		}
	}
	return del;
}

//Method to remove an entry
void removeAt(int location, int no_of_pairs)
{
	for (int i = location; i < no_of_pairs; i++){
		pairs[i] = pairs[i+1];
	}
}

//Selection sort for arranging list in ascending order of probability
void selectionSort(int no_of_pairs)
{
	int minIndex;
	for (int loc = 0; loc < no_of_pairs; loc++)
	{
		minIndex = minLocation(loc, no_of_pairs-1);
		swap(loc, minIndex);
	}
}

//Method used in selection sort
int minLocation(int first, int last)
{
	int minIndex;
	minIndex = first;
	for (int loc = first + 1; loc <= last; loc++)
		if(pairs[loc].prob < pairs[minIndex].prob)
		minIndex = loc;
	return minIndex;
}

//Method to swap two entries
void swap(int first, int second)
{
	Entry temp;
	temp = pairs[first];
	pairs[first] = pairs[second];
	pairs[second] = temp;
}

//Function to display the table
void displayTable(node_collection table){
	int offset = 0;
	for(int i=0; i<table.no_of_nodes;i++){
		for(int j=0; j<table.no_of_mfr; j++){
			node_mfr_link t1;
			t1 = table.set[offset+j];
			cout<<"Node["<<t1.nodeId<<"]"<<endl;
			cout<<"Manufacturer["<<t1.mfrId<<"]"<<endl;
			for(int k=0;k<t1.no_of_entries;k++){
				Entry e1 = t1.entries[k];
				cout<<"Entry["<<k<<"]\tTime:"<<e1.time<<"\tProbability:"<<e1.prob<<"\tCost:"<<e1.cost<<endl;
			}
			t1.entries=NULL;
			t1.entries = 0;
		}
		offset = offset+table.no_of_mfr;
	}
}