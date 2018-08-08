#ifndef _MESH_HPP
#define _MESH_HPP

#include <iostream>
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

        triplet<int>* boundary_condition = NULL;

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

        void init_boundary_condition(int size) {
            this -> boundary_condition = new triplet<int>[size];
        }

        void set_boundary_condition(int index, int eltID, int faceID) {
            (this -> boundary_condition[index]).first = index + 1;
            (this -> boundary_condition[index]).second = eltID;
            (this -> boundary_condition[index]).third = faceID;
        }

        ~Mesh() {
            delete[] this -> nodes;
            delete[] this -> elements;
            delete[] this -> boundary_condition;
        }
};

#endif