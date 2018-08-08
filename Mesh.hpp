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

        ~Mesh() {
            delete[] this -> nodes;
            delete[] this -> elements;
        }
};

#endif