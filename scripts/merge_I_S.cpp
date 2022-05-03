#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

struct node{
    string loopName;
    string distance;
    string ra, rb, rab;
    string ES, FDR, poisson_p_value;
    long long int sta, ena;
    long long int stb, enb;
    string sig;
    friend ostream &operator<<(ostream &output, node& rhs){
        output << rhs.loopName << "\t" << rhs.sta << "\t" << rhs.ena << "\t" \
        << rhs.loopName << "\t" << rhs.stb << "\t" << rhs.enb << "\t" \
        << rhs.ES << "\t" << rhs.FDR << "\t" << rhs.poisson_p_value ;
        return output;
    }
    bool operator<(const node& rhs)const{
        if (loopName != rhs.loopName)
            return loopName < rhs.loopName;
        if (sta != rhs.sta)
            return sta < rhs.sta;
        if (ena != rhs.ena)
            return ena < rhs.ena;
        if (stb != rhs.stb)
            return stb < rhs.stb;
        return enb < rhs.enb;
    }
};
// loopId	distance	ra	rb	rab	ES	FDR	poisson_p-value	iva	ivb	significant
int main(int argc, char* argv[])
{
    if(argc != 3){
        cout << "argc != 3 , arc is "<< argc << endl;
        cout << "Two files in format bedpe is required !" << endl;
        return 0;
    }
    ifstream fin_I(argv[1]);
    ifstream fin_S(argv[2]);

    vector<node> vec;

    auto loadFileData = [&](ifstream &fin){
        string line;
        string ts;
        getline(fin, line);
        while (getline(fin, line)){
            stringstream ss(line);
            node t;
            if(ss){
                // loopId	distance	ra	rb	rab	ES	FDR	poisson_p-value	iva	ivb	significant
                ss>>t.loopName;
                t.loopName = t.loopName.substr(0, t.loopName.find_first_of('-'));
                ss>>t.distance;
                ss>>t.ra;
                ss>>t.rb;
                ss>>t.rab;
                ss>>t.ES;
                ss>>t.FDR;
                ss>>t.poisson_p_value;
                // iva
                ss>>ts;
                t.sta = atoll(ts.substr(ts.find_first_of(':')+1, ts.find_first_of('-')).c_str() );
                t.ena = atoll(ts.substr(ts.find_first_of('-')+1).c_str());
                // ivb
                ss>>ts;
                t.stb = atoll(ts.substr(ts.find_first_of(':')+1, ts.find_first_of('-')).c_str() );
                t.enb = atoll(ts.substr(ts.find_first_of('-')+1).c_str());

                vec.emplace_back(t);
            }
        }

        return;
    };
    loadFileData(fin_I);
    loadFileData(fin_S);

    sort(vec.begin(), vec.end());
    cout << "chrA\tstartA\tendA\tchrB\tstartB\tendB\tES\tFDR\tPoisson_p-value" << endl;

    for (auto &t: vec) {
        cout << t << endl;
    }

    return 0;
}
