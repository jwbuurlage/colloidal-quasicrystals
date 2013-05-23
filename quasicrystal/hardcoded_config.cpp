#include <iostream>
using std::cout;
using std::endl;
#include <fstream>
using std::ofstream;
#include <sstream>
using std::ostringstream;

#include <vector>
using std::vector;
#include <string>
using std::string;

#include <math.h>

struct pos
{
    pos() { x = 0.0; y = 0.0; } 
    pos(double _x, double _y) { x = _x; y = _y; }
    double x;
    double y;
};

vector<pos> points;
double boxSize = 10.0;

void write_snapshot_to_file(double boxSize, int MCS, double pressure)
{
    // open and save file snapshot_MCS.dat
    ostringstream filenamestream;
    filenamestream << "dodecagonalpattern.dat";
    string s = filenamestream.str();
    ofstream fout(s.c_str());

    // number of particles
    fout << points.size() << endl;

    // origin of box
    fout << -boxSize/2.0 << " " << -boxSize/2.0 << " " << -0.5 << endl;

    // span(box)
    fout << boxSize << " " << 0.0 << " " << 0.0 << endl;
    fout << 0.0 << " " << boxSize << " " << 0.0 << endl;
    fout << 0.0 << " " << 0.0 << " " << 1.0 << endl;

    // particle coordinates
    for(unsigned int i = 0; i < points.size(); ++i)
        fout << points[i].x << " " << points[i].y << " " << 0.0 << " " << 1 << " " << 0 << endl;

    // close file
    fout.close();
}

double image_distance_2(pos a, pos b, double boxSize)
{
    double dx, dy;
    dx = fabs(a.x - b.x);
    dy = fabs(a.y - b.y);
    if(dx > boxSize / 2.0) dx -= boxSize;
    if(dy > boxSize / 2.0) dy -= boxSize;

    return dx*dx + dy*dy;
}

int main()
{
    // box size????

    vector<pos> A; // particle of type A
    vector<pos> B; // particle of type B

    A.push_back(pos(0, 0));
    // inner hexagon
    for(double theta = 0.0; theta < 2.0 * M_PI; theta += M_PI/3.0)
        A.push_back(pos(cos(theta), sin(theta)));
    // outer dodecagon
    for(double theta = M_PI / 12.0; theta < 2.0 * M_PI; theta += M_PI/6.0)
        A.push_back(pos(2.0*cos(theta), 2.0*sin(theta)));

    B.push_back(pos(0, 0));
    // inner hexagon
    for(double theta = M_PI/6.0; theta < 2.0 * M_PI; theta += M_PI/3.0)
        B.push_back(pos(cos(theta), sin(theta)));
    // outer dodecagon
    for(double theta = 0; theta < 2.0 * M_PI; theta += M_PI/6.0)
        B.push_back(pos(2.0*cos(theta), 2.0*sin(theta)));

    // this is the pattern
    // A B A ..
    // B A B ..
    // A B A ..
    // . . .
    // . . .

    double currentX = 0.0;
    double currentY = 0.0;
    
    for(int i = 0; i < A.size(); ++i)
    {
        points.push_back(pos(currentX + A[i].x,
                                currentY + A[i].y));
    }

    // remove duplicates
    for(int i = 0; i < points.size(); ++i)
        for(int j = i; i < points.size(); ++i)
        {
            // do i and j overlap then remove
        }

    // write snapshot
    write_snapshot_to_file(boxSize, 0, 0);

    return 0;
}
