#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;
#define pb push_back

typedef struct Point {
    int id, x, y;
} Point;
bool operator<(const Point &a, const Point &b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}
bool operator>(const Point &a, const Point &b) {
    return a.x > b.x || (a.x == b.x && a.y > b.y);
}
int sortPoint(Point a, Point b) {
    if (a.x == b.x) return a.y < b.y;
    return a.x < b.x;
}

int crossProduct(Point root, Point a,
                 Point b) {  // lazy to define another type "vector" so use
                             // point, it both contain x,y anyway
    Point va, vb;
    va.x = a.x - root.x;
    va.y = a.y - root.y;
    vb.x = b.x - root.x;
    vb.y = b.y - root.y;
    int product = va.x * vb.y - va.y * vb.x;
    if (product == 0) return 0;
    return va.x * vb.y - va.y * vb.x;
}
vector<Point> findConvexHull(Point points[], int n) {
    vector<Point> v, ans;
    v.pb(points[0]);
    for (int i = 1; i < n; i++) {
        if (v.size() <= 1)
            v.pb(points[i]);
        else {
            while (v.size() > 1 && crossProduct(*(v.end() - 2), *(v.end() - 1),
                                                points[i]) <= 0) {
                v.pop_back();
            }
            v.pb(points[i]);
        }
    }
    for (auto a : v) {
        ans.pb(a);
    }
    v.clear();
    v.pb(points[n - 1]);
    for (int i = n - 2; i >= 0; i--) {
        if (v.size() < 2)
            v.pb(points[i]);
        else {
            while (v.size() > 1 && crossProduct(*(v.end() - 2), *(v.end() - 1),
                                                points[i]) <= 0) {
                v.pop_back();
            }
            v.pb(points[i]);
        }
    }
    for (int i = 1; i < v.size() - 1; i++) {
        ans.pb(v[i]);
    }
    return ans;
}

vector<Point> merge_convex_hull(vector<Point> &a, vector<Point> &b) {
    int n = a.size();
    int m = b.size();
    int a_rightmost = max_element(a.begin(), a.end()) - a.begin(),
        b_leftmost = min_element(b.begin(), b.end()) - b.begin();

    // Find upper tagent
    int a_upper = a_rightmost, b_upper = b_leftmost;
    bool done = false;
    while (!done) {
        done = true;
        while (crossProduct(b[b_upper], a[a_upper],
                            a[(a_upper - 1 + a.size()) % a.size()]) > 0) {
            a_upper = (a_upper - 1 + a.size()) % a.size();
        }
        while (crossProduct(a[a_upper], b[b_upper],
                            b[(b_upper + 1) % b.size()]) < 0) {
            b_upper = (b_upper + 1) % b.size();
            done = false;
        }
    }

    // Find lower tagent
    int a_lower = a_rightmost, b_lower = b_leftmost;
    done = false;
    while (!done) {
        done = true;
        while (crossProduct(a[a_lower], b[b_lower],
                            b[(b_lower - 1 + b.size()) % b.size()]) > 0) {
            b_lower = (b_lower - 1 + b.size()) % b.size();
        }
        while (crossProduct(b[b_lower], a[a_lower],
                            a[(a_lower + 1) % a.size()]) < 0) {
            a_lower = (a_lower + 1) % a.size();
            done = false;
        }
    }

    // Merge convex hull
    vector<Point> convex_hull;
    for (int i = 0; i != a_upper; i = (i + 1) % a.size()) {
        convex_hull.push_back(a[i]);
    }
    convex_hull.push_back(a[a_upper]);
    for (int i = b_upper; i != b_lower; i = (i + 1) % b.size()) {
        convex_hull.push_back(b[i]);
    }
    convex_hull.push_back(b[b_lower]);
    for (int i = a_lower; i != 0; i = (i + 1) % a.size()) {
        convex_hull.push_back(a[i]);
    }

    return convex_hull;
}
int main(int argc, char **argv) {
    int rank, size, n, local_chunk_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Datatype MPI_Point;
    MPI_Type_contiguous(3, MPI_INT, &MPI_Point);
    MPI_Type_commit(&MPI_Point);
    Point local_points[20000];
    vector<vector<Point>> receive_local_hulls;
    char filename[50];
    Point points[20000];
    vector<Point> v;
    if (rank == 0) {
        scanf("%s", filename);
        FILE *fp = fopen(filename, "r");
        if (fp == NULL) {
            printf("Error: Failed to open file %s\n", filename);
            return 1;
        }
        fscanf(fp, "%d", &n);
        for (int i = 0; i < n; i++) {
            fscanf(fp, "%d %d", &points[i].x, &points[i].y);
            points[i].id = i + 1;
        }
        sort(points, points + n, sortPoint);
        if (size == 1 || n <= 40) {
            v = findConvexHull(points, n);
            reverse(v.begin() + 1, v.end());
            for (auto a : v) {
                printf("%d ", a.id);
            }
            MPI_Type_free(&MPI_Point);
            MPI_Finalize();
            return 0;
        } else {
            int chunk_size = n / size;
            for (int i = 0; i < size; i++) {
                int tmp_start = i * chunk_size;
                int tmp_end = tmp_start + chunk_size;
                int local_n = chunk_size;
                if (i == size - 1) {
                    local_n = n - tmp_start;
                }
                MPI_Send(&local_n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(points + tmp_start, local_n, MPI_Point, i, 0,
                         MPI_COMM_WORLD);
            }
        }
    }

    MPI_Recv(&local_chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    MPI_Recv(local_points, local_chunk_size, MPI_Point, 0, 0, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    vector<Point> local_hull = findConvexHull(local_points, local_chunk_size);

    for (int i = 0; i < size; i++) {
        // send local hull to root
        int local_hull_size = local_hull.size();
        MPI_Send(&local_hull_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(local_hull.data(), local_hull_size, MPI_Point, 0, 0,
                 MPI_COMM_WORLD);
    }

    if (rank == 0) {
        for (int i = 0; i < size; i++) {
            receive_local_hulls.pb(vector<Point>());
        }
        for (int i = 0; i < size; i++) {
            int receive_local_hull_size;
            MPI_Recv(&receive_local_hull_size, 1, MPI_INT, i, 0, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);
            receive_local_hulls[i].resize(receive_local_hull_size);
            MPI_Recv(receive_local_hulls[i].data(), receive_local_hull_size,
                     MPI_Point, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (rank == 0) {
        vector<vector<Point>> ans = receive_local_hulls;
        while (ans.size() > 1) {
            vector<vector<Point>> tmp;
            for (int i = 0; i < ans.size() / 2; i++) {
                tmp.pb(merge_convex_hull(ans[i * 2], ans[i * 2 + 1]));
            }

            ans = tmp;
        }
        reverse(ans[0].begin() + 1, ans[0].end());
        for (auto p : ans[0]) printf("%d ", p.id);
    }

    MPI_Type_free(&MPI_Point);
    MPI_Finalize();
    return 0;
}
