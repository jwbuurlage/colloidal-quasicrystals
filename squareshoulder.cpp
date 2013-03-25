#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <sstream>
using std::ostringstream;
#include <string>
using std::string;
#include <vector>
using std::vector;

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>

// TODO: seed this properly
boost::mt19937 rng;
boost::random::uniform_real_distribution<> ur(0.0, 1.0);
boost::random::uniform_int_distribution<> ui(0, 100000);

// for some benchmarking
#include <sys/time.h> 

#include <cmath>
// for atof
#include <stdlib.h>

/////////////////////////////////
// CONSTANTS
/////////////////////////////////

#define fccsize 30
#define N_POINTS (fccsize * fccsize * 2)
#define SIGMA 1.0
double lambda = 0.0;
double epsilon = 0.0;
double pstart = 0.0;
double pend = 0.0;
double dp = 0.0;
int num_cycles = 0;

// for negative modulo VERY INEFFICIENT FOR LARGE NUMBERS
inline int small_mod(int a, int b)
{
    if(b < 0) return b + a;
    if(b >= a) return b - a;
    else return b;
}

/////////////////////////////////
// HELPERS
/////////////////////////////////

double rand(double min, double plus)
{
    return min + ur(rng)*(plus - min);
}

int rand(int max)
{
    return (ui(rng) % max);
}

struct pos
{
    pos() { x = 0.0; y = 0.0; } 
    pos(double _x, double _y) { x = _x; y = _y; }
    double x;
    double y;
};

void restrict_to_volume(pos* p, double boxSize)
{
    while(p->x > boxSize/2.0) p->x -= boxSize;
    while(p->x < -boxSize/2.0) p->x += boxSize;
    while(p->y > boxSize/2.0) p->y -= boxSize;
    while(p->y < -boxSize/2.0) p->y += boxSize;
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

struct cell
{
    // doubly linked list will be faster
    vector<int> cell_points;

    void add(int i)
    {
        cell_points.push_back(i);
    }

    void remove(int i)
    {
        for (std::vector<int>::iterator it = cell_points.begin(); it != cell_points.end(); ++it)
        {
            if((*it) == i)
            {
                cell_points.erase(it);
                return;
            }
        }
        cerr << "NO SUCH POINT: " << i << endl;
    }

    void clear()
    {
        cell_points.clear();
    }
};

struct particle
{
    particle(pos _p, int _cx, int _cy)
    {
        p = _p;
        cx = _cx;
        cy = _cy;
    }

    pos p;
    int cx;
    int cy;

    // maybe..
    double e;
};

vector<particle> points;
cell cells[N_POINTS][N_POINTS];
double cellSize = 0.0;
int cellCount = 10; 

vector<double> energy_contribution_of_particle;
vector<double> new_energy_contribution_of_particle;
double total_energy = 0.0;

// VERY SLOW ONLY DO SPARINGLY
void cell_indices_for_pos(int& cx, int& cy, pos p, double boxSize)
{
    // the bias should not be necesarry but does fix some bugs
    cx = std::max((int)((p.x + boxSize/2.0)/boxSize * cellCount - 1e6), 0);
    cy = std::max((int)((p.y + boxSize/2.0)/boxSize * cellCount - 1e6), 0);
}

/////////////////////////////////
// SIMULATION
/////////////////////////////////

// initialization
void initialize_points()
{
    double outerOffset = 2.0/fccsize;
    double innerOffset = (2.0 * (1.0 - outerOffset))/(fccsize - 2);

    // set points and enter cell
    for(int i = 0; i < fccsize; ++i)
        for(int j = 0; j < fccsize; ++j)
            points.push_back(particle(pos(-1.0 + outerOffset*(i + 0.5),
                            -1.0 + outerOffset*(j + 0.5)), 0, 0));

    for(int i = 0; i < fccsize; ++i)
        for(int j = 0; j < fccsize; ++j)
            points.push_back(particle(pos(-1.0 + outerOffset + innerOffset*i,
                            -1.0 + outerOffset + innerOffset*j), 0, 0));
}


void set_density(double density, double& boxSize)
{
    // set the appropriate boxSize
    boxSize = pow(points.size() / density, 0.5);
    for(unsigned int i = 0; i < points.size(); ++i)
    {
        points[i].p.x *= boxSize * 0.5;
        points[i].p.y *= boxSize * 0.5;
    }
}

void put_in_cells(double boxSize)
{
    int cx, cy;
    for(unsigned int i = 0; i < points.size(); ++i)
    {
        cell_indices_for_pos(cx, cy, points[i].p, boxSize);
        cells[cx][cy].cell_points.push_back(i);
        points[i].cx = cx;
        points[i].cy = cx;
    }
}

void optimize_cells(double boxSize)
{
    for(int i = 0; i < cellCount; ++i)
        for(int j = 0; j < cellCount; ++j)
            cells[i][j].clear();
    cellCount =  (int)(2.0 * boxSize / (2.0*(lambda * SIGMA))) + 1;
    put_in_cells(boxSize);
}

// physical quantities
inline double potential(pos p1, pos p2, double boxSize)
{
    return image_distance_2(p1, p2, boxSize) < lambda*lambda ? epsilon : 0.0;
}

// energy of a single particle
double energy_at_point(double boxSize, int index, pos newPos)
{
    double energy = 0.0;
    int cx = points[index].cx;
    int cy = points[index].cy;
    cell* c;
    // for cell and neighbours
    for(int dx = -1; dx <= 1; ++dx) 
        for(int dy = -1; dy <= 1; ++dy) 
        {
            c = &cells[small_mod(cellCount, cx + dx)][small_mod(cellCount, cy + dy)];
            for(unsigned int i = 0; i < c->cell_points.size(); ++i)
            {
                if(c->cell_points[i] == index) continue;
                energy += potential(points[c->cell_points[i]].p, newPos, boxSize);
            }
        }

    return energy;
}

// time: about 30 ms
double recalculate_system_energy(double boxSize, vector<double> &container)
{
    double energy = 0.0;

    for(unsigned int i = 0; i < points.size(); ++i)
    {
        container[i] = energy_at_point(boxSize, i, points[i].p);
        energy += container[i];
    }

    // we counted interactions twice
    return energy * 0.5;
}

// output
double average_rho = 0.0;
int n_samples = 0;
void sample(bool output, int MCS, double pressure, double boxSize, bool update = true)
{
    double rho = points.size() / (boxSize * boxSize);

    if(update)
    {
        ++n_samples;
        // update average
        average_rho = average_rho + (rho - average_rho)/n_samples;
    }

    if(output)
        cout << pressure << "  " << MCS << "  " << average_rho << " " << boxSize << "  " << total_energy << endl;
}

void write_snapshot_to_file(double boxSize, int MCS, double pressure)
{
    // open and save file snapshot_MCS.dat
    ostringstream filenamestream;
    filenamestream << "configs/configuration_" << points.size() << "_" << lambda << "_" << epsilon << "_" << pressure << "_" << MCS << ".dat"; // c, n_particles, lambda, pressure, MCS
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
        fout << points[i].p.x << " " << points[i].p.y << " " << 0.0 << " " << SIGMA << " " << 0 << endl;

    // close file
    fout.close();
}

bool particle_move(int index, double jumpSize, double boxSize)
{
    pos newPos = pos(points[index].p.x + rand(-1.0, 1.0) * jumpSize,
            points[index].p.y + rand(-1.0, 1.0) * jumpSize);
    restrict_to_volume(&newPos, boxSize);

    double new_total_energy = total_energy;
    double new_contribution = energy_at_point(boxSize, index, newPos);
    new_total_energy -= energy_contribution_of_particle[index];
    new_total_energy += new_contribution;

    double acc = exp(total_energy - new_total_energy);
    if(rand(0.0, 1.0) > acc)
        return false;

    int cx, cy;
    cell_indices_for_pos(cx, cy, newPos, boxSize);
    cell* c;
    // for cell and neighbours
    for(int dx = -1; dx <= 1; ++dx) 
        for(int dy = -1; dy <= 1; ++dy) 
        {
            c = &cells[small_mod(cellCount, cx + dx)][small_mod(cellCount, cy + dy)];
            for(unsigned int i = 0; i < c->cell_points.size(); ++i)
            {
                if(c->cell_points[i] == index) continue;
                if(image_distance_2(newPos, points[c->cell_points[i]].p, boxSize) < SIGMA*SIGMA)
                    return false;
            }
        }

    int cxp = points[index].cx;
    int cyp = points[index].cy;

    if(points[index].cx != cx || points[index].cy != cy)
    {
        cells[points[index].cx][points[index].cy].remove(index);
        cells[cx][cy].add(index);
        points[index].cx = cx;
        points[index].cy = cy;
    }

    // update energy contributions of nearby points
    for(int dx = -1; dx <= 1; ++dx) 
        for(int dy = -1; dy <= 1; ++dy) 
        {
            c = &cells[small_mod(cellCount, cxp + dx)][small_mod(cellCount, cyp + dy)];
            for(unsigned int i = 0; i < c->cell_points.size(); ++i)
            {
                if(c->cell_points[i] == index) continue;
                energy_contribution_of_particle[c->cell_points[i]] -= potential(points[index].p, points[c->cell_points[i]].p, boxSize);
            }
        }

    for(int dx = -1; dx <= 1; ++dx) 
        for(int dy = -1; dy <= 1; ++dy) 
        {
            c = &cells[small_mod(cellCount, cx + dx)][small_mod(cellCount, cy + dy)];
            for(unsigned int i = 0; i < c->cell_points.size(); ++i)
            {
                if(c->cell_points[i] == index) continue;
                energy_contribution_of_particle[c->cell_points[i]] += potential(newPos, points[c->cell_points[i]].p, boxSize);
            }
        }

    // update the particle itself
    energy_contribution_of_particle[index] = new_contribution;
    total_energy = new_total_energy;
    points[index].p = newPos;

    return true;
}

bool volume_move(double pressure, double vJumpSize, double &boxSize)
{
    double oldVolume = boxSize * boxSize;
    double newVolume = oldVolume + rand(-1.0, 1.0) * vJumpSize;
    if(newVolume <= 0.0) return false;
    double newBoxSize = pow(newVolume, 0.5);


    for(unsigned int i = 0; i < points.size(); ++i) {
        points[i].p.x *= newBoxSize / boxSize;
        points[i].p.y *= newBoxSize / boxSize;
    }

    double new_total_energy = recalculate_system_energy(newBoxSize, new_energy_contribution_of_particle);

    // calculate acc 
    double acc = exp((total_energy - new_total_energy) 
            - pressure * (newVolume - oldVolume) 
            + points.size() * log(newVolume/oldVolume));

    if(rand(0.0, 1.0) > acc)
    {
        for(unsigned int i = 0; i < points.size(); ++i) {
            points[i].p.x *= boxSize / newBoxSize;
            points[i].p.y *= boxSize / newBoxSize;
        }

        return false;
    }

    // check for overlap then still say NO even though initial acceptance rule is ok
    // because infinite energy means acc is zero (we dont have this case above)
    // need new positions here

    for(unsigned int i = 0; i < points.size(); ++i)
    {
        int cx = points[i].cx;
        int	cy = points[i].cy;
        cell* c;
        // for cell and neighbour
        for(int dx = -1; dx <= 1; ++dx) 
            for(int dy = -1; dy <= 1; ++dy) 
            {
                c = &cells[small_mod(cellCount, cx + dx)][small_mod(cellCount, cy + dy)];
                for(unsigned int j = 0; j < c->cell_points.size(); ++j)
                {
                    if(c->cell_points[j] == (int)i) continue;
                    if(image_distance_2(points[i].p, points[c->cell_points[j]].p, newBoxSize) < SIGMA*SIGMA)
                    {
                        // reset points and nothing else happens 
                        for(unsigned int i = 0; i < points.size(); ++i)
                        {
                            points[i].p.x *= boxSize / newBoxSize;
                            points[i].p.y *= boxSize / newBoxSize;
                        }

                        return false;
                    }
                }
            }
    }

    // it is fine
    boxSize = newBoxSize;
    total_energy = new_total_energy;
    energy_contribution_of_particle.swap(new_energy_contribution_of_particle);

    return true;
}

void step(int& accepted, int& vAccepted, double jumpSize, double vJumpSize, double &boxSize, double pressure)
{
    /* timeval t1, t2;
       double elapsedTime;

    // start timer
    gettimeofday(&t1, NULL); */

    for(unsigned int i = 0; i < points.size(); ++i)
    {
        if(particle_move(rand(points.size()), jumpSize, boxSize))
            accepted++;
    }

    /* // stop timer
       gettimeofday(&t2, NULL);

       elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
       elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
       cout << "particle: " << elapsedTime << " ms.\n";


       gettimeofday(&t1, NULL); */

    if(volume_move(pressure, vJumpSize, boxSize))
        vAccepted++;

    /* gettimeofday(&t2, NULL);

       elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
       elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
       cout << "volume: " << elapsedTime << " ms.\n"; */
}

void run()
{
    initialize_points();
    double boxSize = 10.0;
    double jumpSize = 1.0;
    double vJumpSize = 10.0;
    double initialDensity = 0.06;
    set_density(initialDensity, boxSize);
    optimize_cells(boxSize);

    energy_contribution_of_particle.resize(points.size(), 0.0);
    new_energy_contribution_of_particle.resize(points.size(), 0.0);
    // calculate initial energy
    total_energy = recalculate_system_energy(boxSize, energy_contribution_of_particle);

    write_snapshot_to_file(boxSize, 0, 0.0);

    int accepted = 0;
    int vAccepted = 0;
    int adjustInterval = 500;
    int heatCycles = 200000;
    int measureCycles = 50000;

    int MCS_total = 0;

    double pressure = pstart;
    for(; pressure <= pend; pressure += dp)
    {
        int MCS = 1;
        for(; MCS <= heatCycles; ++MCS)
        {
            step(accepted, vAccepted, jumpSize, vJumpSize, boxSize, pressure);

            if(MCS % adjustInterval == 0)
            {
                //adjust jumpsize according to rate
                double rate = ((double)accepted)/(points.size() * adjustInterval);
                if(rate > 0.32)
                    jumpSize *= 1.05;
                else if(rate < 0.28)
                    jumpSize *= 0.95;

                if(jumpSize > boxSize) jumpSize = boxSize;

                double vRate = ((double)vAccepted)/adjustInterval;
                if(vRate > 0.16)
                    vJumpSize *= 1.05;
                else if(vRate < 0.14)
                    vJumpSize *= 0.95;

                accepted = 0;
                vAccepted = 0;

                sample(true, MCS_total, pressure, boxSize, false);

                optimize_cells(boxSize);
            }

            if(MCS % 10000 == 0)
                write_snapshot_to_file(boxSize, MCS, pressure);
        }

        average_rho = 0.0;
        n_samples = 0;

        for(; MCS <= heatCycles + measureCycles; ++MCS)
        {
            step(accepted, vAccepted, jumpSize, vJumpSize, boxSize, pressure);

            if(MCS % 10000 == 0)
                write_snapshot_to_file(boxSize, MCS, pressure);

            sample((MCS % 500 == 0) ? true : false, MCS, pressure, boxSize);
        }


        MCS_total += MCS - 1;
        sample(true, MCS_total, pressure, boxSize);

        heatCycles = 50000;
    }
}

int main(int argc, char* argv[])
{
    if(argc != 7)
    {
        cout << "Usage: ./sqsh MCS lambda epsilon pstart pend dp" << endl;
        cout << "This will heat at pstart for 3*MCS and then spend MCS for every pressure in (pstart, pstart+dp, ..., pend)" << endl;
        return 0;
    }

    num_cycles = atoi(argv[1]);
    lambda  = atof(argv[2]);
    epsilon = atof(argv[3]);
    pstart  = atof(argv[4]);
    pend    = atof(argv[5]);
    dp      = atof(argv[6]);

    run();

    return 0;
}
