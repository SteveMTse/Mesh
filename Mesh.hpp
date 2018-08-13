#ifndef _MESH_HPP
#define _MESH_HPP

#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

struct Point {
    int index;
    double x;
    double y;
};

struct Polygon {
    int index;
    vector<int> vertices;
};

template <typename A, typename B = A, typename C = B>
struct triplet {
    A first;
    B second;
    C third;
};

class Mesh {
    private:
        //Domain
        double xleft = -1;
        double xright = -1;
        double ybottom = -1;
        double ytop = -1;

        //Nodes and Elements and Coordinates
        int number_of_nodes = -1;
        int number_of_elements = -1;
        Point* nodes = NULL;

        Polygon* elements = NULL;

        triplet<int>* dirichlet = NULL;
        triplet<int>* neumann = NULL;

    public:
        Mesh(double xl, double xr, double yb, double yt, int NNode, int NElt) {
            this -> xleft = xl;
            this -> xright = xr;
            this -> ybottom = yb;
            this -> ytop = yb;
            this -> number_of_nodes = NNode;
            this -> number_of_elements = NElt;
            this -> nodes = new Point[NNode];
            this -> elements = new Polygon[NElt];
        }

        Mesh(string filename) {
            ifstream data;
            data.open(filename, ios::in);
            if(data.is_open()) {
                //TODO;
            }
        }

        void set_node(int index, double x, double y) {
            (this -> nodes[index]).index = index + 1;
            (this -> nodes[index]).x = x;
            (this -> nodes[index]).y = y;
        }

        void init_polygon(int index, int NVertices) {
            (this -> elements[index]).index = index + 1;
            (this -> elements[index]).vertices.reserve(NVertices);
        }

        void set_polygon(int index, int elt) {
            (this -> elements[index]).vertices.push_back(elt);
        }

        void init_dirichlet(int size) {
            this -> dirichlet = new triplet<int>[size];
        }

        void set_dirichlet(int index, int eltID, int faceID) {
            (this -> dirichlet[index]).first = index + 1;
            (this -> dirichlet[index]).second = eltID;
            (this -> dirichlet[index]).third = faceID;
        }

        void init_neumann(int size) {
            this -> neumann = new triplet<int>[size];
        }

        void set_neumann(int index, int eltID, int faceID) {
            (this -> neumann[index]).first = index + 1;
            (this -> neumann[index]).second = eltID;
            (this -> neumann[index]).third = faceID;
        }

        ~Mesh() {
            if(this -> nodes != NULL) delete[] this -> nodes;
            if(this -> elements != NULL) delete[] this -> elements;
            if(this -> dirichlet != NULL) delete[] this -> dirichlet;
            if(this -> neumann != NULL) delete[] this -> neumann;
        }
};

#endif