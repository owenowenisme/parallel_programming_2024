#pragma GCC optimize(3)

#include<mpi.h>
#include<iostream>
#include<vector>
#include<utility>
#include<queue>
#include<string.h>

#define MAX_N 5000
#define MAX_N_2 2500000000L
#define F first
#define S second
#define INF 1e9

#define min(a,b) ((a)<(b)?(a):(b))

short adjacencyMatrix[MAX_N_2];

inline long long calculateIndex(int sourceVertex, int targetVertex, int totalVertices) {
    long long source = sourceVertex;
    long long target = targetVertex;
    long long total = totalVertices;
    return source * total + target;
}

int main(int argc, char *argv[]){
    int vertexCount;
    int processCount;
    int processId;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);

    std::vector<int> shortestDistances;

    char inputFilename[100];
    if(processId == 0){
        scanf("%s",inputFilename);
    }
    MPI_Bcast(inputFilename,100,MPI_CHAR,0,MPI_COMM_WORLD);
    

    FILE *input_file = fopen(inputFilename,"r");
    if (input_file == NULL){
        printf("File not found\n");
        return 0;
    }
    // read vertexCount
    fscanf(input_file,"%d",&vertexCount);
    while(!feof(input_file)){
        int sourceVertex, targetVertex, edgeWeight;
        fscanf(input_file,"%d %d %d\n",&sourceVertex,&targetVertex,&edgeWeight);
        adjacencyMatrix[calculateIndex(sourceVertex, targetVertex, vertexCount)] = edgeWeight;
    }
    fclose(input_file);

    if(processId == 0){
        shortestDistances.resize(vertexCount,INF);
        shortestDistances[0] = 0;
    }

    // divide the load
    int workerCount = processCount - 1;
    int verticesPerWorker = vertexCount / workerCount;


    // coordinator
    if( processId == 0){
        std::vector<int> receivedDistances(vertexCount); // recv distances from workers
        /* Dijkstra Start */
        std::priority_queue<std::pair<int,int>,std::vector<std::pair<int,int> >,std::greater<std::pair<int,int> > > pq;
        pq.push(std::make_pair(0,0));
        while(!pq.empty()){
            std::pair<int,int> currentNode = pq.top();
            pq.pop();
            int currentVertex = currentNode.S;
            int currentDistance = currentNode.F;
            if(currentDistance > shortestDistances[currentVertex]) continue;

            // no need to distribute the work ( the graph is small )
            if( vertexCount <= workerCount ){
                for(int v=0;v<vertexCount;v++){
                    // printf("graph[%d][%d] = %d\n",u,v,graph[u][v]);
                    if(adjacencyMatrix[calculateIndex(currentVertex, v, vertexCount)] == 0) continue;
                    // printf("u = %d, v = %d, dis[v] = %d, dis[u] = %d, wt = %d\n",u,v,distances[v],distances[u],graph[u][v]);
                    if(shortestDistances[v] > shortestDistances[currentVertex] + (int)adjacencyMatrix[calculateIndex(currentVertex, v, vertexCount)]){
                        shortestDistances[v] = shortestDistances[currentVertex] + (int)adjacencyMatrix[calculateIndex(currentVertex, v, vertexCount)];
                        pq.push(std::make_pair(shortestDistances[v],v));
                    }
                }
                continue;
            }

            // distribute the work
            for(int ith_worker=0;ith_worker<workerCount;ith_worker++){
                MPI_Send(&currentVertex,1,MPI_INT,ith_worker+1,0,MPI_COMM_WORLD);
                // send the distances[u] to all workers, instead of 
                MPI_Send(&shortestDistances[currentVertex],1,MPI_INT,ith_worker+1,0,MPI_COMM_WORLD);
                int start_node = ith_worker * verticesPerWorker;
                int end_node = ((ith_worker) == (workerCount-1) ? vertexCount: start_node + verticesPerWorker);
            }

            for(int ith_worker=0;ith_worker<workerCount;ith_worker++){
                int start_node = ith_worker * verticesPerWorker;
                int end_node = ((ith_worker) == (workerCount-1) ? vertexCount: start_node + verticesPerWorker);
                MPI_Recv(receivedDistances.data()+start_node,end_node-start_node,MPI_INT,ith_worker+1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            }
            for(int v=0;v<vertexCount;v++){
                if(shortestDistances[v] > receivedDistances[v]){
                    shortestDistances[v] = receivedDistances[v];
                    pq.push(std::make_pair(shortestDistances[v],v));
                }
            }
        }
        /* Dijkstra End */

        // send exit message after the dijkstra is done
        for(int ith_worker=0;ith_worker<workerCount;ith_worker++){
            int from = -1;
            MPI_Send(&from,1,MPI_INT,ith_worker+1,0,MPI_COMM_WORLD);
        }
    }
    else{ // workers
        int currentVertex;
        int workerStartVertex = (processId-1) * verticesPerWorker;
        int workerEndVertex = ((processId-1) == (workerCount-1) ? vertexCount: workerStartVertex + verticesPerWorker);
        int sourceDistance;
        std::vector<int> workerDistances(workerEndVertex - workerStartVertex,INF);
        while(true){

            MPI_Recv(&currentVertex,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            if(currentVertex == -1){ // exit message received
                break;
            }

            MPI_Recv(&sourceDistance,1,MPI_INT,0,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            
            for(int i=0,v=workerStartVertex;v<workerEndVertex;i++,v++){
                if(adjacencyMatrix[calculateIndex(currentVertex, v, vertexCount)] == 0) continue;
                if(workerDistances[i] > sourceDistance + adjacencyMatrix[calculateIndex(currentVertex, v, vertexCount)]){
                    workerDistances[i] = sourceDistance + adjacencyMatrix[calculateIndex(currentVertex, v, vertexCount)];
                }
            }
            
            MPI_Send(workerDistances.data(),workerEndVertex - workerStartVertex,MPI_INT,0,0,MPI_COMM_WORLD);
        }
    }

    // print the result
    MPI_Barrier(MPI_COMM_WORLD);
    if(processId == 0){
        for(int i=0;i<vertexCount;i++){
            printf("%d ",shortestDistances[i]);
        }
    }

    MPI_Finalize();

    return 0;
}