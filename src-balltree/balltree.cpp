/*************************************************
*
* Author: Nazneen Rajani
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <vector>
#include <ctype.h>
#include <assert.h>
#include <sys/time.h>
#include <string>
#include <iostream>
#include <time.h>
#include <omp.h>

using namespace std;

int read_X_from_file(double **X, int n, int D, char *filename);


struct datapoint {
    double *dim; //list of dimensions
    int idx;  //its own index
      datapoint(double *dim, int idx){
        this -> dim = dim;
        this -> idx = idx;
    }
};

struct centroid {
    double *dim; //list of dimensions
};

struct ballnode {
    vector<struct datapoint> *data; //list of nodes it owns
    double radius;
    struct centroid *pivot;
    vector<struct datapoint*> child1;
    vector<struct datapoint*> child2;
};


struct balltree {
    int num_points;
    struct ballnode *root;
};

struct balltree *balltree_create(int n)
{
    struct balltree *tree;
    tree->num_points = n;
    tree->root = 0;
    return tree;
}

unsigned int get_msec(void)
{
static struct timeval timeval, first_timeval;

gettimeofday(&timeval, 0);

if(first_timeval.tv_sec == 0) {
first_timeval = timeval;
return 0;
}
return (timeval.tv_sec - first_timeval.tv_sec) * 1000 + (timeval.tv_usec - first_timeval.tv_usec) / 1000;
}

pair<double,struct datapoint*> getRadius(const struct centroid& target, vector<datapoint>& pts, int D) {
    double radius =0.0;
    struct datapoint *child1;
    for (int k=0; k<pts.size(); k++) {
         double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target.dim[i] - pts[k].dim[i]) * (target.dim[i] - pts[k].dim[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            child1 = &pts[k];
        }
    }
    return make_pair(radius,child1);
}

struct datapoint* getMaxDist(const struct datapoint& target, vector<datapoint>& pts, int D) {
    double radius =0.0;
    struct datapoint *child2;
    for (int k=0; k<pts.size(); k++) {
        double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target.dim[i] - pts[k].dim[i]) * (target.dim[i] - pts[k].dim[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            child2 = &pts[k];
        }
    }
    return child2;
}

double getDistance( datapoint& key, datapoint& curr, int D) {
    double dist = 0.0;
    for (int i = 0; i < D; ++i) {
        dist += (key.dim[i] - curr.dim[i]) * (key.dim[i] - curr.dim[i]);
    }
    return sqrt(dist);
}


int balltree_insert_all(struct balltree *tree, double **dim_all, int n, int D)
{   vector<datapoint> points;
    int j;
    points.reserve(n);
    for (j=0; j<n; j++) {
        points.emplace_back(dim_all[j],j+1);
    }
    struct ballnode *root;
    root-> data = &points;
    struct centroid *centroid;
    double dim[D] = {};
    for (int k=0; k<D; k++) {
        for (int i=0; i<n; i++) {
            dim[k] += points.at(i).dim[k];
        }
    }
    centroid -> dim = dim;
    root->pivot = centroid;
    pair<double, struct datapoint*> answer = getRadius(*root->pivot,*root->data,D);
    root->radius = answer.first;
    struct datapoint *child1 = answer.second;
    struct datapoint *child2 = getMaxDist(*child1,*root->data,D);
    vector<datapoint*> child1_points;
    vector<datapoint*> child2_points;
    for (int k=0; k<points.size(); k++) {
        datapoint dt =points.at(k);
        if(getDistance(*child1,dt,D)<getDistance(*child2,dt,D))
            child1_points.push_back(&dt);
        else
            child2_points.push_back(&dt);
    }
    root->child1 = child1_points;
    root->child2 = child2_points;
    return 0;
}


int main(int argc, char **argv)
{
int i, j;
unsigned int msec, start;
double end, start1;
int n = 10000;                 // # data points
int D = 28*28;                  // data dimension
char fname[] = "/Users/nrajani/kdtree-knn/data/mnist-test.dat"; // data file
double **X;
X=(double**)malloc(sizeof(double*)*n);
for (i=0; i<n; i++) {
X[i] = (double*)malloc(sizeof(double)*D);
}

if (read_X_from_file(X, n, D, fname)==0)
{
printf("error open file.\n");
exit(1);
}
    
struct balltree *btree = balltree_create(n);
balltree_insert_all(btree,X,n,D);
    cout << "done "<<endl;
    return 0;
}

int read_X_from_file(double **X, int n, int D, char *filename)
{
FILE *fp = NULL;
if (!(fp = fopen(filename, "rb")))
return 0;
int i;
int num_in, num_total;
for (i = 0; i < n; i++)
{
num_total = 0;
while (num_total < D)
{
num_in = fread(X[i]+num_total, sizeof(double), D, fp);
num_total += num_in;
}
}

fclose(fp);
return 1;
}
