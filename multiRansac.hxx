#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <cmath>
#include <cfloat>
#include <iostream>
#include <assert.h>
//#include <opencv2/opencv.hpp>
using namespace std;


template<class T,class U>
class Ransac
{
public:
    virtual vector<U> computeParameters(const vector<T> &data,int numunknowns)=0;
    virtual double errorfunction(T data,const vector<U> &unknowns)=0;
    void computeRansac()
    {
        cout<<"Computing RanSaC with global inlier probability: "<<m_prob<<" and inlier threshold: "<<m_inlierthresh<<endl;
        //Initialize algorithm
        if(m_data.size()==0)
        {
            cout<<"No input data!"<<endl;
            return;
        }

        vector<T> m_sampleset;
        vector<bool> samplechosen,best_samplechosen;
        m_sampleset=getSample(m_numunknowns,samplechosen);
        if(!m_useinitialparams)
            m_unknowns=computeParameters(m_sampleset,m_numunknowns);




        //Initialize
        int maxIters=1e5;
        double numsample=maxIters;

        int samplecount=0;
        double epsilon=1.0; //outlierfraction
        double minerror=100000000;

        ///Estimate inlier and model parameters
        /// Basically the algorithm goes like this:
        ///
        /// a) randomly sample the minimum amount of data points
        /// needed to compute a hypothesis of the model parameters with computeParameters().
        /// The current hypothesis is stored in m_unknowns
        ///
        /// b) Now iterate through all data point (datapt) in m_data.
        /// For each datapt compute it model error, aka how well  datapt fits to the current hypothesis
        /// in m_unknowns.
        /// If errorfunction(datapt,m_unknowns)<m_inlierthresh then datapt is an inlier for the current hypothesis.
        /// m_inlierthresh is model specific and provided by the user
        ///
        /// c) In the end the hypothesis with the largest set of inliers is the winner.
        ///
        /// d) How many times should we repeat a)-c) ?
        ///
        /// m_prob is a given parameter that specifies the fraction of m_data  that
        /// are inliers to the best hypothesis found after
        /// numsamples trials of running a)-c).
        /// 1-m_prob conversely determines the outlier fraction after
        /// numsamples trials of running a)-c).
        ///
        /// Since every run of a) is an independent sampling process, 1-m_prob
        /// is a product of the individual  outlierfractions per run:
        /// 1-m_prob=(1-m_inliers.size()/m_data.size())^numsamples
        /// This means we can dynamically compute the maximum number of runs (numsamples)
        /// by numsamples=log(1-m_prob)/log(1-m_inliers.size()/m_data.size())

        while(numsample>samplecount)
        {
            int numinlier=0;
            double toterr=0;
            vector<int> inlier;
            vector<double> inlierweight;
            int samplesize=m_numunknowns;//round(((double)rand())/RAND_MAX*(m_data.size()-1))+1;

            //Sample MINIMAL set of data needed to compute model parameters (current hypothesis)
            vector<T> sampleset=getSample(samplesize,samplechosen);

            //Compute parameters
            if(!m_useinitialparams)
                m_unknowns=computeParameters(sampleset,m_numunknowns);
            for(int i=0;i<m_data.size();i++)
            {
                //Test all m_data (except for those points used to compute current hypothesis)
                if(!samplechosen[i])
                {
                    double err=errorfunction(m_data[i],m_unknowns);
                    //Does m_data[i] fit the current hypothesis?
                    if(pow(err,2.0)<=pow(m_inlierthresh,2.0))
                    {
                        //Yes? Then add the error to the total error.
                        numinlier++;
                        toterr+=err;
                        inlier.push_back(i);
                        double w=1.0-fabs(err)/fabs(m_inlierthresh);
                        inlierweight.push_back(w);
                    }
                    else
                        toterr+=m_inlierthresh;///No? All outliers should contribute EQUALLY to the total error (we dont want to model the outliers)
                }
                else
                {
                    //The sampled data are also inliers
                    inlier.push_back(i);
                    inlierweight.push_back(1.0);
                }
            }
            numinlier+=m_numunknowns;
            if(sampleset.size())
            {
                double neweps=1.0-((double)numinlier)/m_data.size();
                //Is the outlier fraction and the total error of the current hypothesis
                //lower then the current best hypothesis
                if(neweps<epsilon&&toterr<minerror)//&&inlier.size()>m_inliers.size())
                {
                    //Update outlierfraction and best total error
                    epsilon=neweps;
                    minerror=toterr;
                    double num=max(1.0-m_prob,DBL_MIN);

                    double denom=1.0-pow(1.0-neweps,((double)samplesize));

                    //Update numsamples
                    if(denom<DBL_MIN)
                    {
                        numsample=0;
                        epsilon=0.0;
                        m_inliers=inlier;
                        m_inlierweight=inlierweight;
                    }
                    else
                    {
                        num = std::log(num);
                        denom = std::log(denom);
                        numsample=denom >= 0 || -num >= maxIters*(-denom) ? maxIters : num/denom;
                        cout<<"Numsample: "<<numsample<<endl;
                        m_sampleset=sampleset;
                        best_samplechosen=samplechosen;
                        m_inliers=inlier;
                        m_inlierweight=inlierweight;

                    }

                }
            }

            samplecount++;



        }

        samplechosen=best_samplechosen;
        //We probably always have outliers
        if(epsilon==0.0)
        {
            cout<<"Warning: unrealistic probability or threshold setting"<<endl;
            //m_inliers=m_data;
        }

        //
        int acceptablesamplesize=(1.0-epsilon)*m_data.size();
        cout<<"Found acceptable samplesize: "<<acceptablesamplesize<<endl;
//        numsample+=sqrt(1-pow((1-epsilon),m_numunknowns))/pow((1-epsilon),m_numunknowns);
//        numsample=round(numsample);
        cout<<"Number of iterations: "<<round(numsample)<<endl;
        cout<<"Computed outlierfraction: "<<epsilon<<endl;

        vector<T> datain;
        for(int din=0;din<m_inliers.size();din++)
            datain.push_back(m_data[m_inliers[din]]);
        m_unknowns=computeParameters(datain,m_numunknowns);











    }

private:
    vector<T> m_data;
    vector<U> m_unknowns;
    int m_numunknowns;
    bool m_useinitialparams;


    vector<int> m_inliers;
    vector<double> m_inlierweight;

    double m_prob; //Probability that sampleset s is free from outliers
    double m_inlierthresh; //Threshhold error that a point in sampleset s is an inlier
    vector<T> getSample(int samplesize,vector<bool> &samplechosenret)
    {
       vector<T> res;
       vector<bool> samplechosen(m_data.size(),false);

       for(int s=0;s<samplesize;s++)
       {

//           int sampleidx=rng.uniform(0,m_data.size()-1);//((double)rand())/RAND_MAX*((double)m_data.size()-1);
           int sampleidx=((double)rand())/RAND_MAX*((double)m_data.size()-1);
           while(samplechosen[sampleidx])
                sampleidx=((double)rand())/RAND_MAX*((double)m_data.size()-1);
           res.push_back(m_data[sampleidx]);
           samplechosen[sampleidx]=true;
       }
       samplechosenret=samplechosen;
       return res;
    }




public:
    Ransac(vector<T> &data,int numunknowns,double prob=0.99,double inlierthresh=3,bool useinitialparam=false)
    {
        m_data=data;


        m_numunknowns=numunknowns;

        m_prob=prob;
        m_inlierthresh=inlierthresh;
        m_useinitialparams=useinitialparam;
        srand(time(0));



    }
    vector<U> unknowns()
    {
        return m_unknowns;
    }
    vector<int> inliers()
    {
        return m_inliers;
    }
    vector<bool> inliermask()
    {
        vector<bool> res(m_data.size(),false);
        for(int i=0;i<m_inliers.size();i++)
            res[m_inliers[i]]=true;
        return res;

    }
    vector<double> inlierweights()
    {
        return m_inlierweight;
    }

    vector<double> inlierweightmask()
    {
        vector<double> res(m_data.size(),0.0);
        for(int i=0;i<m_inliers.size();i++)
            res[m_inliers[i]]=m_inlierweight[i];
        return res;

    }




};
