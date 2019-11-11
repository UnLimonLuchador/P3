/// @file



#include <iostream>
#include <fstream>
#include <math.h>
#include "pitch_analyzer.h"

using namespace std;

std::ofstream potencia;
std::ofstream ro;
std::ofstream romax;
std::ofstream zero;
std::ofstream Amdf;
std::ofstream preclip;
std::ofstream clip;
/// Name space of UPC
namespace upc {
  void PitchAnalyzer::autocorrelation(const vector<float> &x, vector<float> &r) const {

    for (unsigned int l = 0; l < r.size(); ++l) {
  		/// \HECHO Compute the autocorrelation r[l]
      float sum = 0;
      for(unsigned int i = 0; i < (r.size()-l); i++){
        sum = sum + x[i]*x[l+i];
        
      }
      r[l]=sum;
     

    }
    //cout << r[0] << '\t' << r[1] << '\t' << r[2] << '\t' << r[3] << '\t' << r[4] << '\t' << r[5] << endl;
    //cout << x[0] << '\t' << x[1] << '\t' << x[2] << '\t' << x[3] << '\t' << x[4] << '\t' << x[5] << endl;

 

    if (r[0] == 0.0F) //to avoid log() and divide zero 
      r[0] = 1e-10; 

    
  }
  void PitchAnalyzer::AMDF(const vector<float> &x, vector<float> &amdf) const {
      for (unsigned int l = 0; l < amdf.size(); ++l) {
          /// \HECHO Compute the autocorrelation r[l]
          float sum = 0;
          for(unsigned int i = 0; i < (amdf.size()-l); i++){
            sum += abs(x[i]-x[l+i]);            
          }
          amdf[l]=sum;
        

      }
  }
  float PitchAnalyzer::zerocrossing(const vector<float> &x) const {
      float zcr=0;
      int i=0;
      
      for(i=1; i<x.size(); i++) {
        if(x[i-1]<0 && x[i]>0 || x[i-1]>0 && x[i]<0){
          zcr++;
        }
      }
      zcr = (samplingFreq*zcr)/(2*(i-1));
      return zcr;
  }

  void PitchAnalyzer::set_window(Window win_type) {
    if (frameLen == 0)
      return;

    this->window.resize(this->frameLen);

    switch (win_type) {
    case HAMMING:
    {
      /// \HECHO Implement the Hamming window
      float a0 = 0.53836;
      float a1 = 0.46164;
      for(unsigned int i = 0; i<frameLen; i++){
        this->window[i] = a0 -a1*cos((2*M_PI*i)/(this->frameLen - 1));
      }
      break;
    }
    case RECT:
      window.assign(frameLen, 1);
      break;
    
    
    default:
      window.assign(frameLen, 1);
      break;
    }
  }

  void PitchAnalyzer::set_f0_range(float min_F0, float max_F0) {
    npitch_min = (unsigned int) samplingFreq/max_F0;
    if (npitch_min < 2)
      npitch_min = 2;  // samplingFreq/2

    npitch_max = 1 + (unsigned int) samplingFreq/min_F0;

    //frameLen should include at least 2*T0
    if (npitch_max > frameLen/2)
      npitch_max = frameLen/2;
  }

  bool PitchAnalyzer::unvoiced(float pot, float r1norm, float rmaxnorm, float zcr) const {
    /// \HECHO Implement a rule to decide whether the sound is voiced or not.
    /// * You can use the standard features (pot, r1norm, rmaxnorm),
    ///   or compute and use other ones.
  /*
    potencia.open("pot.txt", std::ios_base::app);
    ro.open("ro.txt", std::ios_base::app);
    romax.open("romax.txt", std::ios_base::app);
    zero.open("zcr.txt", std::ios_base::app);

    potencia << pot << '\n';   
    ro << r1norm << '\n'; 
    romax << rmaxnorm << '\n';    
    zero << zcr << '\n';         

    potencia.close();
    ro.close();
    romax.close();
    zero.close();
    // Cerrar el fichero, 
    // para luego poder abrirlo para lectura:
    */
       
    if(r1norm > 0 && rmaxnorm > 0.1 && zcr < 2000 && pot > -25){      
      return false;      
    }
    return true;
    
  }

  float PitchAnalyzer::compute_pitch(vector<float> & x) const {
    if (x.size() != frameLen)
      return -1.0F;

    //Window input frame
    for (unsigned int i=0; i<x.size(); ++i){
      x[i] *= window[i];
      
        //if( i >= 1 || i==x.size()) x[i]=(1/3)*(x[i-1]+x[i+1]+x[i]); 
          
    }

    //center clipping
    /*
    vector<float>::const_iterator iX = x.begin(), iXMax = iX;
    for(;iX<x.end();iX++){    
      if(*iX>*iXMax) iXMax = iX;
    }    
    
    for(unsigned int i = 0; i<x.size();i++){
      preclip.open("pre-clipping.txt", std::ios_base::app);
        preclip << x[i] << '\n';         
        preclip.close();
    }   
    
    float xth = *iXMax*0.3;

    for(unsigned int i = 0; i<x.size();i++){    
      
      //if(x[i]<-xth) x[i] += xth;       
      //else if(x[i]>xth) x[i] -= xth;
      //else x[i]=0;
      
    }
/*
     for(unsigned int i = 0; i<x.size();i++){
      clip.open("center-clipping.txt", std::ios_base::app);
        clip << x[i] << '\n';           
        clip.close();
    }
    */



    vector<float> r(npitch_max);
    vector<float> amdf(npitch_max);

    //Compute correlation
    autocorrelation(x, r);
    AMDF(x,amdf);
    float z = zerocrossing(x);

    vector<float>::const_iterator iR = r.begin(), iRMax = iR;

    /// \TODO 
	/// Find the lag of the maximum value of the autocorrelation away from the origin.<br>
	/// Choices to set the minimum value of the lag are:
	///    - The first negative value of the autocorrelation.
	///    - The lag corresponding to the maximum value of the pitch.
    ///	   .
	/// In either case, the lag should not exceed that of the minimum value of the pitch.
    
    
    vector<float>::const_iterator iA = amdf.begin(), iAMin = iA;
    //for(iA = amdf.begin(); *iA<1; ++iA);
    iAMin = iA + npitch_min;
    for(iA = amdf.begin() + npitch_min + 1;iA<amdf.end();iA++){    
      if(*iA<*iAMin) iAMin = iA;
    }
    

    for(iR = r.begin(); *iR>0; ++iR);
    iRMax = iR;
    for(;iR<r.end();iR++){
      if(*iR>*iRMax) iRMax = iR;
    }

    unsigned int lagIR = iRMax - r.begin();   
    unsigned int lagAM = iAMin - amdf.begin(); 

    float pot = 10 * log10(r[0]);

    //You can print these (and other) features, look at them using wavesurfer
    //Based on that, implement a rule for unvoiced
    //change to #if 1 and compile
#if 1
    if (r[0] > 0.0F)
      //cout << pot << '\t' << r[1]/r[0] << '\t' << r[lag]/r[0] << endl;
#endif
    
    if (unvoiced(pot, r[1]/r[0], r[lagIR]/r[0], z )){
      
      ofstream fs("corr_unvoiced.txt"); 

    // Enviamos una cadena al fichero de salida:
      for(unsigned int i = 0; i<r.size(); i++){
        fs << r[i] << endl;      
      }
      
      
    // Cerrar el fichero, 
    // para luego poder abrirlo para lectura:
      fs.close();
      ofstream fa("amdf_unvoiced.txt"); 

   // Enviamos una cadena al fichero de salida:
    for(unsigned int i = 0; i<r.size(); i++){
      fa << amdf[i] << endl;      
    }
    
    fa.close();

      return 0;
    }
    else{
         // Crea un fichero de salida
    ofstream fs("corr_voiced.txt"); 

   // Enviamos una cadena al fichero de salida:
    for(unsigned int i = 0; i<r.size(); i++){
      fs << r[i] << endl;      
    }
    
    fs.close();

    ofstream fa("amdf_voiced.txt"); 

   // Enviamos una cadena al fichero de salida:
    for(unsigned int i = 0; i<r.size(); i++){
      fa << amdf[i] << endl;      
    }
    
    fa.close();

      if(lagIR == 0){
        cout << "fallo" << endl;
        return 0;
      }
      return (float) samplingFreq/(float) lagIR;
    }
     
  }
  
}
