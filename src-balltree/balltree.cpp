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

vector<datapoint> pts;

struct centroid {
    double *dim; //list of dimensions
};

struct ballnode {
    vector<struct datapoint*> data; //list of nodes it owns
    double radius;
    struct centroid *pivot;
    struct ballnode* child1;
    struct ballnode* child2;
    ballnode(vector<struct datapoint*> data){
        this -> data = data;
    }
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

pair<double,struct datapoint*> getRadius(const struct centroid *target, vector<datapoint*>& pts, int D) {
    double radius =0.0;
    struct datapoint *child1;
    for (int k=0; k<pts.size(); k++) {
         double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target->dim[i] - pts[k]->dim[i]) * (target->dim[i] - pts[k]->dim[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            child1 = pts[k];
        }
    }
    return make_pair(radius,child1);
}

struct datapoint* getMaxDist(const struct datapoint *target, vector<datapoint*>& pts, int D) {
    double radius =0.0;
    struct datapoint *child2;
    for (int k=0; k<pts.size(); k++) {
        double dist = 0.0;
        for (int i = 0; i < D; ++i) {
            dist += (target->dim[i] - pts[k]->dim[i]) * (target->dim[i] - pts[k]->dim[i]);
        }
        if(sqrt(dist)>radius){
            radius = sqrt(dist);
            child2 = pts[k];
        }
    }
    return child2;
}

double getDistance(struct datapoint* key, struct datapoint* curr, int D) {
    double dist = 0.0;
    for (int i = 0; i < D; ++i) {
        dist += (key->dim[i] - curr->dim[i]) * (key->dim[i] - curr->dim[i]);
    }
    return sqrt(dist);
}

void recursive_insert(struct ballnode *root, int items, int D){
    struct centroid centroid;
    double dim[D] = {};
    for (int k=0; k<D; k++) {
        for (int i=0; i<items; i++) {
            dim[k] += root->data.at(i)->dim[k];
        }
    }
    centroid.dim = dim;
    root->pivot = &centroid;
    //cout<<"done"<<endl;
    pair<double, struct datapoint*> answer = getRadius(root->pivot,root->data,D);
    root->radius = answer.first;
    //cout<<"done"<<endl;
    struct datapoint *child1_point = answer.second;
    struct datapoint *child2_point = getMaxDist(child1_point,root->data,D);
    
    //cout<<"done"<<endl;
    
    vector<datapoint*> child1_points;
    vector<datapoint*> child2_points;
    //cout<<"done"<<endl;
    for (int k=0; k<items; k++) {
        datapoint *dt = root->data.at(k);
        if(getDistance(child1_point,dt,D)<getDistance(child2_point,dt,D))
            child1_points.push_back(dt);
        else
            child2_points.push_back(dt);
    }
    struct ballnode *child1 = new ballnode(child1_points);
    struct ballnode *child2 = new ballnode(child2_points);
    
    root->child1 = child1;
    root->child2 = child2;

    if(child1_points.size()>2){ //I might change
        recursive_insert(root->child1, child1_points.size(),D); //recursive insert on both childs
    }
    if(child2_points.size()>2) //I might change
        recursive_insert(root->child2, child2_points.size(),D);
}

int balltree_insert_all(struct balltree *tree, double **dim_all, int n, int D)
{
    int j;
    pts.reserve(n);
    for (j=0; j<n; j++) {
        pts.emplace_back(dim_all[j],j);
    }

    vector<datapoint*> points;
    points.reserve(n);
    for (j=0; j<n; j++) {
        points.push_back(&pts.at(j));
    }

    struct ballnode *root = new ballnode(points);
    
    struct centroid centroid;
    double dim[D] = {};
    
    for (int k=0; k<D; k++) {
        for (int i=0; i<n; i++) {
            dim[k] += points.at(i)->dim[k];
        }
    }
    centroid.dim = dim;
    root->pivot = &centroid;
    pair<double, struct datapoint*> answer = getRadius(root->pivot,root->data,D);
    root->radius = answer.first;

    struct datapoint* child1_point = answer.second;
    struct datapoint* child2_point = getMaxDist(child1_point,root->data,D);
    
    vector<datapoint*> child1_points;
    vector<datapoint*> child2_points;
    
    for (int k=0; k<points.size(); k++) {
        datapoint *dt = points.at(k);
        if(getDistance(child1_point,dt,D)<getDistance(child2_point,dt,D))
            child1_points.push_back(dt);
        else
            child2_points.push_back(dt);
    }
    
    struct ballnode *child1 = new ballnode(child1_points);
    struct ballnode *child2 = new ballnode(child2_points);
    
    root->child1 = child1;
    root->child2 = child2;
    
    if(child1_points.size()>2) // might change
        recursive_insert(root->child1, child1_points.size(),D); //recursive insert on both childs
    
    if(child2_points.size()>2) //might change
        recursive_insert(root->child2, child2_points.size(),D);
    
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

double start_build = omp_get_wtime();
struct balltree *btree = balltree_create(n);
balltree_insert_all(btree,X,n,D);
    double end_build = omp_get_wtime()-start_build;
    cout<<"Time to create ball tree: "<<end_build<<endl;
    
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
