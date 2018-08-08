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

    public:
        Mesh(double xl, double xr, double yb, double yt, int NNode, int NElt) {
            this -> xleft = xl;
            this -> xright = xr;
            this -> ybottom = yb;
            this -> ytop = yb;
            this -> number_of_nodes = NNode;
            this -> number_of_elements = NElt;
            this -> nodes = new Point[NNode];
        }

        ~Mesh() {
            delete[] this -> nodes;
        }
};

#endif