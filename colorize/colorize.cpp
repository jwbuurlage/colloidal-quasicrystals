// takes a configuration
// and outputs a colored one
// it is okay if this is slow (it will be)

#include <iostream>
using std::cout;
using std::cin;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;
#include <sstream>
using std::ostringstream;
#include <string>
using std::string;

#include <vector>
using std::vector;

// for fabs
#include <cmath>

struct pos
{
    pos() { x = 0.0; y = 0.0; }
    pos(double _x, double _y) { x = _x; y = _y; }
    double x;
    double y;
};

double boxSize;

double image_distance_2(pos a, pos b)
{
    double dx, dy;
    dx = fabs(a.x - b.x);
    dy = fabs(a.y - b.y);
    if(dx > boxSize / 2.0) dx -= boxSize;
    if(dy > boxSize / 2.0) dy -= boxSize;

    return dx*dx + dy*dy;
}

void colorize_file(string f)
{
    vector<pos> particles;

    ifstream fin(f.c_str());
    ostringstream of;
    of << "colored_" << f;
    ofstream fout(of.str().c_str());
    cout << "outputting to: " << of.str() << endl;

	int n_particles;
	fin >> n_particles;

	double dummy;
	double origin_x, origin_y, origin_z;
	fin >> origin_x >> origin_y >> origin_z;

	double  span1x, span1y,
			span2x, span2y;
	fin >> span1x >> span1y;	
	fin >> dummy;
	fin >> span2x >> span2y;
	fin >> dummy >> dummy >> dummy >> dummy;

	double x, y;
	for(int i = 0; i < n_particles; ++i)
	{
		fin >> x >> y >> dummy >> dummy >> dummy;
		particles.push_back(pos(x, y));
	}

    // output header
	fout << n_particles << endl;
	fout << origin_x << " " << origin_y << "  " << origin_z << endl;
	fout << span1x << " " << 0.0 << "  " << 0.0 << endl;
	fout << 0.0 << " " << span2y << "  " << 0.0 << endl;
	fout << 0.0 << " " << 0.0 << "  " << 1.0 << endl;
	
	boxSize = span1x;

	// find neighbours of particle i, 
	for(int i = 0; i < particles.size(); ++i)
	{
		int particle_type = 0;
		// count first are particles that are touching
		int count_first = 0;
		// count secnod are particles that are slightly further away
		int count_second = 0;
		// look at the near particles and count them by either distance = SIGMA or distance ~= 1.4*SIGMA
		for(int j = 0; j < particles.size(); ++j)
		{
			if(i == j) continue;
			double d = image_distance_2(particles[i], particles[j]);
			if(d < 1.4)
				count_first++;
			else if(d < 2.4)
				count_second++;

			if(count_first == 6 && count_second == 0) // hexagonal
				particle_type = 1;
			else if(count_first == 5 && count_second == 2) // edge
				particle_type = 2;
			else if(count_first == 4 && count_second == 4) // square
				particle_type = 3;

		}

		fout << particles[i].x << " " << particles[i].y << " " << 0.0 << "  " << 1.0 << " " << particle_type << endl;
	}
}

int main(int argc, char* argv[])
{
    if(argc < 2) 
    {
        cout << "usage: colorize [files]" << endl;
        return 0;
    }
    
    for(int i = 1; i < argc; ++i)
    {
        string s = string(argv[i]);
        colorize_file(s);
    }

	return 0;
}
