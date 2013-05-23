#include <iostream>
using std::cout;
using std::cin;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <sstream>
using std::ostringstream;
using std::stringstream;

#include <string.h>

#include <string>
using std::string;

struct Record
{
    Record();
    Record(string _pressure, 
            string _energy, 
            int _samples,
            float ca,
            float cb,
            float cc,
            float cd)
    {
        pressure = _pressure;
        energy = _energy;
        samples = _samples;
        count_averages[0] = ca;
        count_averages[1] = cb;
        count_averages[2] = cc;
        count_averages[3] = cd;
    }

    ~Record() { }

    string pressure;
    string energy;
    int samples;
    float count_averages[4]; // in [0, 1]
};

Record* currentRecord;

void writeCurrentRecordToFile()
{
    stringstream energy_file;
    energy_file << "e" << currentRecord->energy << "_counter.txt";

    // append to 'energy_file'
    ofstream fout(energy_file.str().c_str(), ofstream::out | ofstream::app);

    // APPEND to file....
    // --------------------------------------------------------------
    // | pressure | count 0 | count 1 | count 2 | count 3 | samples |
    // --------------------------------------------------------------

    fout << currentRecord->pressure << " " << currentRecord->count_averages[0] 
                                    << " " << currentRecord->count_averages[1] 
                                    << " " << currentRecord->count_averages[2] 
                                    << " " << currentRecord->count_averages[3] 
                                    << " " << currentRecord->samples << endl;
}

void count(string f)
{
    ifstream fin(f.c_str());

    // grab the name of the energy
    // file name is: configuration_n_l_e_p_MCS.dat
    
    // find position of last two underscores
    int b, a;
    int currentPosition = f.length() - 1;
    while(f[currentPosition] != '_')
        currentPosition--;
    b = currentPosition;
    currentPosition--;
    while(f[currentPosition] != '_')
        currentPosition--;
    a = currentPosition + 1;

    string pressure = f.substr(a, (b - a));

    b = currentPosition;
    currentPosition--;
    while(f[currentPosition] != '_')
        currentPosition--;
    a = currentPosition + 1;

    string energy = f.substr(a, (b - a));

    if(currentRecord->pressure != pressure 
            || currentRecord->energy != energy)
    {
        if(currentRecord->energy != "")
            writeCurrentRecordToFile();
        // make new record
        currentRecord = new Record(pressure, energy, 0, 0.0f, 0.0f, 0.0f, 0.0f);
    }


    // ignore header
    // yes this is very ugly should just skip the first x bytes YOLO
    // read only last 'particle' type and increase counter
    int count[4];
    memset(count, 0, sizeof(count));
	
    int npart;
    double ddummy;
    fin >> npart; // nparticles
    fin >> ddummy >> ddummy >> ddummy; // origins
    fin >> ddummy >> ddummy >> ddummy; // span1
    fin >> ddummy >> ddummy >> ddummy; // span2
    fin >> ddummy >> ddummy >> ddummy; // span3

    int type;
    for(int i = 0; i < npart; ++i)
    {
        fin >> ddummy >> ddummy >> ddummy >> ddummy >> type;
        count[type]++;
    }

    currentRecord->samples++;
    for(int t = 0; t < 4; ++t)
    {
        float fraction = count[t] / (float)npart;
        currentRecord->count_averages[t] += (fraction - currentRecord->count_averages[t])/currentRecord->samples;
    }
}

// can we input it sorted by pressure? would make things easier
// yes... we can

int main(int argc, char* argv[])
{
    if(argc < 2) 
    {
        cout << "usage: counter [files]" << endl;
        return 0;
    }
    
    currentRecord = new Record("", "", 0, 0.0f, 0.0f, 0.0f, 0.0f);

    for(int i = 1; i < argc; ++i)
    {
        string s = string(argv[i]);
        count(s);
    }

    writeCurrentRecordToFile();
}
