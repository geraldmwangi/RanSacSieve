#include <iostream>
#include <multiRansac.hxx>
#include <vector>
#include <map>
#include <algorithm>
using namespace std;


class FrequencySieve: public Ransac<pair<float,float>, float>
{
public:
    FrequencySieve(vector<pair<float,float> >data,double thresh):
         Ransac(data,1,0.99,thresh)
    {
        computeRansac();
    }


    virtual double errorfunction(pair<float,float> data,const vector<float> &currF0)
    {
        //Are the two frequencies of the duplet multiples of the given f0
        int mult1=round(data.first/currF0[0]);
        int mult2=round(data.second/currF0[0]);
        double err=fabs(data.first-mult1*currF0[0])+fabs(data.second-mult2*currF0[0]);
        err+=fabs(fabs(data.first-data.second)-currF0[0]);
        err/=3.0;
        return fabs(err);///norm;
    }

    virtual vector<float> computeParameters(const vector<pair<float,float> > &data,int numunknowns)
    {
        vector<float> newF0;
        //We get a list of duplets and compute the f0 as the median of the difference data[i].first-data[i].second
        vector<float>freqs;
        for(int i=0;i<data.size();i++)
            freqs.push_back(fabs(data[i].first-data[i].second));//freq+=fabs(data[i].first-data[i].second);//min(fabs(data[i].first-data[i].second),freq);
        sort(freqs.begin(),freqs.end());
        double freq=freqs[freqs.size()/2];
        if(freqs.size()>1)
        {
            if(freqs.size()%2==0)
                freq=freqs[(freqs.size())/2.0];
            else
                freq=freqs[(freqs.size()+1)/2.0];
        }
        newF0.push_back(freq);

        return newF0;
    }
};
int main(int argc, char *argv[])
{
    //Generate harmonic frequencies
    double freq=82.41;
    double freq2=110.0;
    double freq3=146.8;
    int numharmonics=20;
    int numharmonics2=20;
    int numharmonics3=20;

    //Generate arbitrary peaks
    int numArbitraryPeaks=5;

    //The error the noisy harmonic combs can have compared to the (unknown) true combs
    double thresh=5;

    //The frequencies of the harmonic peaks are distorted +-noise Hz
    double noise=5;
    vector<float> frequencies;
    srand(time(0));
    for(int i=0;i<numharmonics;i++)
    {
        frequencies.push_back((i+1)*freq+(((double)rand())/RAND_MAX*2.0-1.0)*noise);
        cout<<freq<<": "<<frequencies.back()<<endl;
    }
    for(int i=0;i<numharmonics2;i++)
    {
        frequencies.push_back((i+1)*freq2+(((double)rand())/RAND_MAX*2.0-1.0)*noise);
        cout<<freq2<<": "<<frequencies.back()<<endl;
    }
    for(int i=0;i<numharmonics3;i++)
    {
        frequencies.push_back((i+1)*freq3+(((double)rand())/RAND_MAX*2.0-1.0)*noise);
        cout<<freq3<<": "<<frequencies.back()<<endl;
    }

    for(int i=0;i<numArbitraryPeaks;i++)
    {
        double freq=((double)rand())/RAND_MAX*(1000-10)+10;
        frequencies.push_back(freq);
    }
    vector<float> mixfrequencies;
    srand(time(0));
    //Mix up the frequencies
    while(frequencies.size())
    {
        int id=((double)rand())/RAND_MAX*frequencies.size();
        mixfrequencies.push_back(frequencies[id]);
        cout<<"Mix: "<<mixfrequencies.back()<<endl;
        vector<float> trunk;
        for(int i=0;i<id;i++)
            trunk.push_back(frequencies[i]);
        for(int i=(id+1);i<frequencies.size();i++)
            trunk.push_back(frequencies[i]);
        frequencies=trunk;
    }
    frequencies=mixfrequencies;

    //Our model depends on pairs of harmonics
    vector<pair<float,float> > data;
    for(int i=0;i<(frequencies.size()-1);i++)
        for(int j=i+1;j<frequencies.size();j++)
            data.push_back(pair<float,float>(frequencies[i],frequencies[j]));
    vector<double> unknowns;
    unknowns.push_back(0.0);
    srand(time(0));
    cout<<"Data size"<<data.size()<<endl;

    //Seed the pseudo-random number generator (PRNG)
    srand(time(0));
    FrequencySieve estim(data,thresh);
    cout << "Ransac Test!" << endl;
    //The estimated F0
    cout<<"Found F0: "<<estim.unknowns()[0]<<endl;
    cout<<estim.inliers().size()<<endl;

    //On the outlier data of the first run we run another sieve
    vector<pair<float,float> > outlierdata;
    vector<bool> inlier=estim.inliermask();
    for(int i=0;i<inlier.size();i++)
        if(!inlier[i])
            outlierdata.push_back(data[i]);
    vector<double> unknowns2;
    unknowns2.push_back(0.0);
    cout<<"outlierdata size"<<outlierdata.size()<<endl;

    //Re-seed the PRNG
    srand(time(0));
    FrequencySieve estim2(outlierdata,thresh);
    cout << "2Ransac Test!" << endl;
    cout<<"Found F0: "<<estim2.unknowns()[0]<<endl;
    cout<<estim2.inliers().size()<<endl;

    //And a third sieve
    vector<pair<float,float> > outlierdata2;
    inlier=estim2.inliermask();
    for(int i=0;i<inlier.size();i++)
        if(!inlier[i])
            outlierdata2.push_back(outlierdata[i]);
    vector<double> unknowns3;
    unknowns3.push_back(0.0);
    cout<<"outlierdata size"<<outlierdata2.size()<<endl;

    srand(time(0));
    FrequencySieve estim3(outlierdata2,thresh);
    cout << "2Ransac Test!" << endl;
    cout<<"Found F0: "<<estim3.unknowns()[0]<<endl;
    cout<<estim3.inliers().size()<<endl;

    return 0;
}
