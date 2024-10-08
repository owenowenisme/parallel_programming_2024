#include <stdio.h>
#include <mpi.h>
# define ll long long int
int coverage_num[50];
int cost[50];
ll coverage_arr[50];
int n,m;

int main(int argc, char** argv) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ll global_ans = 0;
    char filename[50];
    if (rank == 0){
        scanf("%s", filename);
        FILE *fp = fopen(filename, "r");
        if (fp == NULL){
            printf("Error: Failed to open file %s\n", filename);
            return 1;
        }
        fscanf(fp, "%d %d", &n, &m);
        for (int i=0;i<m;i++){
            fscanf(fp, "%d %d",coverage_num+i, cost+i);
            for (int j=0;j<coverage_num[i];j++){
                int tmp;
                fscanf(fp, "%d", &tmp);
                coverage_arr[i] |= (1<<tmp-1);
            }
        }

        
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(coverage_num, 50, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(coverage_arr, 50, MPI_LONG_LONG, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ll local_ans = 0;
    ll permutations = 1LL<<m;
    ll tmp = 0;
    ll chunk_size = permutations/size;
    for (ll i=rank*chunk_size;i<(rank+1)*chunk_size;i++){
        tmp = i;
        ll sum = 0;
        for(int j=0;j<m;j++){
            if (tmp & 1){
                sum |= coverage_arr[j];
            }
            tmp>>=1;
        }
        if (sum == (1<<n)-1){
            local_ans++;
        }
    }
    MPI_Reduce(&local_ans, &global_ans, 1, MPI_LONG_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0){
        printf("%lld", global_ans);
    }
    MPI_Finalize();
    return 0;
}cd ..

    
